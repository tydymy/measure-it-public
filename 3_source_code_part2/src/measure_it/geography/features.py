"""Geographic burden, context and opportunity features (SPEC section E, Phase 4; Tests 4 and 6).

Pre-registered plan: docs/ANALYSIS_PLAN_GEOGRAPHY.md. Burden definitions / proxy map:
docs/BURDEN_DEFINITIONS.md (mirrored below as BURDEN_DEFINITIONS).

Outputs (data/processed):
  geo_condition_burden     long table: every burden estimate, harmonised, with evidence level,
                           source resolution, inherited flag and measure role (primary / proxy / alternate)
  geo_burden_definitions   the pre-specified primary measure per condition x geography level
  geo_context              one wide row per state / county / ZCTA (ACS, SVI, RUCC, HPSA/MUA, PLACES)
  geo_condition_features   one row per (state|county, condition): burden, vulnerability, relevant
                           provider density, HRSA, HPSA, trials, NIH, diagnostic_desert + components
Results: results/tables/test4_*.csv, test6_geography_*.csv, geo_*.csv; results/figures/drafts/geo_*.png;
results/GEOGRAPHY_RESULTS.md (generated from those files).

Guardrails enforced here: ecological joins only; a state value attached to a county keeps
source_geographic_resolution='state' and inherited=True; proxies are level C and never silently primary;
conditions without a measure get level D and a null burden with a reason.

Run: uv run python -m measure_it.geography.features
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

from .. import config
from ..config import FIGURES, TABLES, UNKNOWN, load_config, utc_now_iso
from ..provenance import PROVENANCE_COLUMNS
from ..store import read_table, table_exists, write_table
from .crosswalk import STATE_FIPS

PRODUCER = "measure_it.geography.features"
DRAFTS = FIGURES / "drafts"

# 50 states + DC: the analysis universe (PR ranks SVI separately; island areas lack ACS/SVI/PLACES/HPS).
TERRITORY_FIPS = {"60", "66", "69", "72", "78"}
UNIVERSE_STATES = tuple(sorted(set(STATE_FIPS) - TERRITORY_FIPS))

PRIMARY_CMS_YEAR = 2022
ALT_CMS_PRELIM_YEAR = 2023
PRIMARY_LYME_YEAR = 2023
ALT_LYME_PRE2022_YEAR = 2019
HPS_POOLED_FROM = "2024-01-01"
DESERT_TRIAL_RADIUS_KM = 50
PROVIDER_RADIUS_KM = 50
RADII_KM = (50, 100)
TOP_N_STATE = 10
N_BOOT = 2000
N_PERM = 1000
N_PERM_MORAN = 999
MORAN_BAND_KM = 100
CONTIGUITY_BUFFER_M = 100

ACTIVE_STATUSES = ("RECRUITING", "NOT_YET_RECRUITING", "ENROLLING_BY_INVITATION", "ACTIVE_NOT_RECRUITING")
LEVEL_MULT = load_config("scoring")["evidence_level_uncertainty_multiplier"]

# PLACES proxy map (docs/BURDEN_DEFINITIONS.md). county_proxy = the separate county C-proxy column.
PLACES_PROXY = {"long_covid": "PHLTH", "me_cfs": "PHLTH", "fibromyalgia": "ARTHRITIS"}
PLACES_ALTERNATES = {"long_covid": ["COGNITION", "DISABILITY"]}
PLACES_PROXY_WHY = {
    "PHLTH": "frequent physical distress (>=14 poor-physical-health days of 30); symptom/function proxy",
    "ARTHRITIS": "BRFSS arthritis item names fibromyalgia among its conditions; dominated by osteoarthritis",
    "COGNITION": "cognitive disability; alternate long-COVID symptom proxy",
    "DISABILITY": "any disability; alternate long-COVID function proxy",
}


# ============================================================================ burden definitions
def _def(condition, level, measure_id, ev, resolution, *, inherited=False, derivation="as_published",
         rationale="", alternates="", count_kind=None, proxy=None):
    return {"condition_id": condition, "geo_level": level, "primary_measure_id": measure_id,
            "burden_evidence_level": ev, "source_geographic_resolution": resolution, "inherited": inherited,
            "derivation": derivation, "rationale": rationale, "alternates": alternates,
            "burden_count_kind": count_kind, "county_proxy_places_measure": proxy}


_LC_WHY = ("HPS 'currently experiencing long COVID', all adults: the only direct public long-COVID measure "
           "(NCHS gsea-w83j; series ends 2024-09-16).")
_CMS51_WHY = ("CCW 'Fibromyalgia, Chronic Pain and Fatigue' (MMD 51), Medicare FFS, 2022 final claims, "
              "unsmoothed actual, all ages.")
_LYME_WHY = "CDC reported Lyme cases (confirmed+probable) by county of residence, 2023 (latest; 2022+ case definition)."
_D_WHY = {
    "pots": "No public population measure: no CCW algorithm contains G90.A; PLACES has no orthostatic item.",
    "dysautonomia": "No public population measure: G90.x absent from MMD-exposed CCW algorithms; no PLACES item.",
    "eds_hsd": "No public population measure: Q79.6 / M35.7 absent from MMD-exposed CCW algorithms.",
    "mcas": "No public population measure: D89.4 absent from MMD-exposed CCW algorithms.",
    "post_infectious_syndrome": "Grouping node; member measures cannot be unioned without overlap information.",
    "ibs": "No public population measure: K58 absent from MMD-exposed CCW algorithms; no PLACES item.",
    "gastroparesis": "No public population measure: K31.84 absent from MMD-exposed CCW algorithms; no PLACES item.",
    "endometriosis": "No public population measure: N80 absent from MMD-exposed CCW algorithms; no PLACES item.",
}
_CMS_ALTS = (f"cms_mmd_<code>_prev_agestd_all {PRIMARY_CMS_YEAR}; cms_mmd_<code>_prev_actual_all "
             f"{ALT_CMS_PRELIM_YEAR} (preliminary); cms_mmd_<code>_prev_actual_lt65 {PRIMARY_CMS_YEAR}")
_LYME_ALTS = f"lyme_incidence_per_100k_2022_2023 (pooled); lyme incidence {ALT_LYME_PRE2022_YEAR} (pre-2022 definition)"


def _build_definitions() -> list[dict]:
    d = []
    for lvl in ("national", "state"):
        d.append(_def("long_covid", lvl, "lc_current_pct_all_adults", "A", lvl, rationale=_LC_WHY
                      + (" Latest non-suppressed period per state." if lvl == "state" else " Latest period."),
                      alternates="lc_current_pct_all_adults_2024_mean; strict latest period (suppressed -> null)",
                      count_kind="survey_prevalence_x_acs_adults" if lvl == "state" else None))
    d.append(_def("long_covid", "county", "lc_current_pct_all_adults", "A", "state", inherited=True,
                  derivation="inherited_from_state",
                  rationale="No county long-COVID measure exists (HPS is national/state; PLACES 2025 has no COVID "
                            "item). The state value is carried as inherited context, never as county prevalence.",
                  alternates="places_PHLTH_crude (county C-proxy column); places_COGNITION/DISABILITY crude+ageadj",
                  proxy="PHLTH"))
    for lvl in ("national", "state", "county"):
        d.append(_def("me_cfs", lvl, "cms_mmd_51_prev_actual_all", "C", lvl,
                      rationale=_CMS51_WHY + " Contains R53.82 (chronic fatigue) but not G93.3x: symptom proxy.",
                      alternates=_CMS_ALTS.replace("<code>", "51") + ("; places_PHLTH_crude" if lvl == "county" else ""),
                      proxy="PHLTH" if lvl == "county" else None))
        d.append(_def("fibromyalgia", lvl, "cms_mmd_51_prev_actual_all", "B", lvl,
                      rationale=_CMS51_WHY + " Contains M79.7 (fibromyalgia); broad pain/fatigue composite.",
                      alternates=_CMS_ALTS.replace("<code>", "51") + ("; places_ARTHRITIS_crude" if lvl == "county" else ""),
                      proxy="ARTHRITIS" if lvl == "county" else None))
        d.append(_def("migraine", lvl, "cms_mmd_59_prev_actual_all", "B", lvl,
                      rationale="CCW 'Migraine and Other Chronic Headache' (MMD 59; G43 + all G44), Medicare FFS, "
                                "2022 final, unsmoothed actual, all ages.",
                      alternates=_CMS_ALTS.replace("<code>", "59")))
    d.append(_def("lyme_disease", "county", "lyme_incidence_per_100k", "A", "county", derivation="computed_rate",
                  rationale=_LYME_WHY + " Rate = cases / Census PEP population (ingestion's computation).",
                  alternates=_LYME_ALTS, count_kind="surveillance_cases"))
    for lvl in ("state", "national"):
        d.append(_def("lyme_disease", lvl, "lyme_incidence_per_100k_agg", "A", "county",
                      derivation="aggregated_from_county",
                      rationale=_LYME_WHY + " Sum of county cases / sum of county populations; excludes cases "
                                            "with unknown county.",
                      alternates=_LYME_ALTS.replace("lyme_incidence_per_100k_2022_2023",
                                                    "lyme_incidence_per_100k_agg_2022_2023"),
                      count_kind="surveillance_cases"))
    for lvl, mid in (("county", "lyme_incidence_per_100k"), ("state", "lyme_incidence_per_100k_agg"),
                     ("national", "lyme_incidence_per_100k_agg")):
        d.append(_def("ptlds", lvl, mid, "C", "county", derivation="copied_as_proxy",
                      rationale="Antecedent-exposure proxy: PTLDS follows treated Lyme disease; no PTLDS surveillance "
                                "or survey exists. Lyme incidence is not PTLDS prevalence.",
                      alternates=_LYME_ALTS))
    for cond, why in _D_WHY.items():
        for lvl in ("national", "state", "county"):
            d.append(_def(cond, lvl, None, "D", None, derivation=None, rationale=why))
    return d


BURDEN_DEFINITIONS: list[dict] = _build_definitions()


def definitions_frame() -> pd.DataFrame:
    df = pd.DataFrame(BURDEN_DEFINITIONS)
    df["uncertainty_multiplier"] = df["burden_evidence_level"].map(LEVEL_MULT)
    return df


def target_conditions() -> list[str]:
    return [c["id"] for c in load_config("conditions")["conditions"]]


# ============================================================================ pure helpers
def percentile_rank(s: pd.Series) -> pd.Series:
    """Percentile rank in (0, 1] with average ties; NaN stays NaN."""
    return s.rank(pct=True, method="average")


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(a))


def radius_matrix(src_lat, src_lon, dst_lat, dst_lon, radius_km: float) -> sparse.csr_matrix:
    """Boolean sparse matrix (n_src x n_dst): dst point within radius_km of src point (haversine)."""
    from sklearn.neighbors import BallTree
    dst = np.radians(np.column_stack([dst_lat, dst_lon]))
    tree = BallTree(dst, metric="haversine")
    src = np.radians(np.column_stack([src_lat, src_lon]))
    ind = tree.query_radius(src, r=radius_km / 6371.0088)
    rows = np.repeat(np.arange(len(ind)), [len(i) for i in ind])
    cols = np.concatenate(ind) if len(ind) else np.array([], dtype=int)
    return sparse.csr_matrix((np.ones(len(cols), dtype=np.int32), (rows, cols)), shape=(len(src), len(dst)))


def spearman(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    if ok.sum() < 3:
        return np.nan
    return float(stats.spearmanr(x[ok], y[ok]).statistic)


def bootstrap_spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = config.SEED, alpha: float = 0.05):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 5:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for b in range(n_boot):
        i = rng.integers(0, n, n)
        rx, ry = stats.rankdata(x[i]), stats.rankdata(y[i])
        out[b] = np.corrcoef(rx, ry)[0, 1] if rx.std() > 0 and ry.std() > 0 else np.nan
    return float(np.nanpercentile(out, 100 * alpha / 2)), float(np.nanpercentile(out, 100 * (1 - alpha / 2)))


def partial_spearman(x, y, z) -> float:
    """Rank-based partial correlation of x and y controlling for z (residuals of ranks on rank(z))."""
    df = pd.DataFrame({"x": x, "y": y, "z": z}).dropna()
    if len(df) < 5:
        return np.nan
    r = df.rank()
    zc = np.column_stack([np.ones(len(r)), r["z"].to_numpy()])
    res = []
    for c in ("x", "y"):
        beta, *_ = np.linalg.lstsq(zc, r[c].to_numpy(), rcond=None)
        res.append(r[c].to_numpy() - zc @ beta)
    if res[0].std() == 0 or res[1].std() == 0:
        return np.nan
    return float(np.corrcoef(res[0], res[1])[0, 1])


def top_n_ids(values: pd.Series, n: int) -> set:
    """Ids (index) of the n largest non-null values; deterministic tie-break by index."""
    v = values.dropna()
    s = pd.DataFrame({"v": v.to_numpy(), "k": v.index.astype(str)}, index=v.index)
    s = s.sort_values(["v", "k"], ascending=[False, True])
    return set(s.index[:n])


def topn_inclusion_prob(values: pd.Series, n: int) -> pd.Series:
    """Probability that each id is in the top-n when ties at the boundary are broken at random.

    1 for values strictly above the n-th value, (slots left)/(tied count) for values equal to it, else 0.
    Reviewer addition: CMS MMD prevalences are published as integers, so many counties tie at the top-N boundary
    and a deterministic tie-break (by FIPS) decides membership.
    """
    v = values.dropna()
    p = pd.Series(0.0, index=v.index)
    if len(v) == 0 or n <= 0:
        return p
    if n >= len(v):
        return p + 1.0
    thr = v.sort_values(ascending=False).iloc[n - 1]
    above, tied = v > thr, v == thr
    p[above] = 1.0
    p[tied] = (n - int(above.sum())) / int(tied.sum())
    return p


def topn_overlap(a: pd.Series, b: pd.Series, n: int) -> dict:
    """Overlap of the top-n of a and of b among ids where both are non-null.

    `overlap` uses the deterministic tie-break of top_n_ids (pre-specified); `overlap_tie_averaged` is the expected
    overlap under random tie-breaking and `ties_at_boundary_a/_b` count values tied at each boundary.
    """
    common = a.dropna().index.intersection(b.dropna().index)
    m = len(common)
    if m == 0 or n <= 0:
        return {"M": m, "N": n, "overlap": np.nan, "jaccard": np.nan, "expected_overlap": np.nan, "p_hypergeom": np.nan,
                "overlap_tie_averaged": np.nan, "ties_at_boundary_a": np.nan, "ties_at_boundary_b": np.nan}
    n_eff = min(n, m)
    ta, tb = top_n_ids(a.loc[common], n_eff), top_n_ids(b.loc[common], n_eff)
    k = len(ta & tb)
    pa, pb = topn_inclusion_prob(a.loc[common], n_eff), topn_inclusion_prob(b.loc[common], n_eff)

    def _ties(p):
        return int(((p > 0) & (p < 1)).sum())
    return {"M": m, "N": n_eff, "overlap": k, "jaccard": k / len(ta | tb),
            "expected_overlap": n_eff * n_eff / m,
            "p_hypergeom": float(stats.hypergeom.sf(k - 1, m, n_eff, n_eff)),
            "overlap_tie_averaged": float((pa * pb.reindex(pa.index)).sum()),
            "ties_at_boundary_a": _ties(pa), "ties_at_boundary_b": _ties(pb)}


def permute_within_groups(values: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle values among positions that share a group label (e.g. counties within a state)."""
    keys = rng.random(len(values))
    order = np.lexsort((keys, groups))           # positions sorted by group then random key
    base = np.lexsort((np.arange(len(values)), groups))  # positions sorted by group (stable)
    out = np.empty_like(values)
    out[base] = values[order]
    return out


def row_standardize(w: sparse.csr_matrix) -> sparse.csr_matrix:
    w = w.tocsr().astype(float)
    rs = np.asarray(w.sum(axis=1)).ravel()
    inv = np.divide(1.0, rs, out=np.zeros_like(rs), where=rs > 0)
    return sparse.diags(inv) @ w


def morans_i(x: np.ndarray, w: sparse.csr_matrix) -> float:
    """Moran's I for a row-standardised weight matrix (rows without neighbours must be removed first)."""
    z = np.asarray(x, float) - np.mean(x)
    den = z @ z
    if den == 0:
        return np.nan
    n, s0 = len(z), w.sum()
    return float((n / s0) * (z @ (w @ z)) / den)


def permutation_summary(observed: float, null: np.ndarray, no_structure: float) -> dict:
    """Summary of a permutation null.

    p_perm (pre-specified) is two-sided around the no-structure value. That is the right test for the national
    shuffle, whose null is centred there, but NOT for the within-state shuffle, whose null is centred on the
    between-state component. p_perm_vs_null_centre (reviewer addition) is two-sided around the null's own mean:
    for the within-state scheme it tests whether the observed value differs from what between-state
    differences alone produce (i.e. whether a within-state association exists).
    """
    null = null[~np.isnan(null)]
    lo, hi = np.percentile(null, [2.5, 97.5])
    dev_obs = abs(observed - no_structure)
    p = (1 + np.sum(np.abs(null - no_structure) >= dev_obs - 1e-12)) / (1 + len(null))
    c = null.mean()
    p_c = (1 + np.sum(np.abs(null - c) >= abs(observed - c) - 1e-12)) / (1 + len(null))
    return {"null_mean": float(null.mean()), "null_p2_5": float(lo), "null_p97_5": float(hi), "p_perm": float(p),
            "p_perm_vs_null_centre": float(p_c),
            "observed_outside_null95": bool(observed < lo or observed > hi),
            "retained_fraction": float((null.mean() - no_structure) / (observed - no_structure))
            if observed != no_structure else np.nan}


def holm(p) -> np.ndarray:
    """Holm step-down adjusted p-values (family = the values passed; NaN ignored)."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return out
    order = np.argsort(pv)
    adj = np.maximum.accumulate((m - np.arange(m)) * pv[order])
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def min_detectable_rho(n: float, alpha: float = 0.05, power: float = 0.8) -> float:
    """Smallest |rho| detectable with the given power (two-sided, Fisher z approximation)."""
    if n is None or not np.isfinite(n) or n <= 4:
        return np.nan
    z = (stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)) / np.sqrt(n - 3)
    return float(np.tanh(z))


def crh_effective_n(x, y, lat, lon, band_km: float = 100.0, max_km: float = 2000.0) -> dict:
    """Clifford, Richardson & Hemon (1989) effective sample size for a correlation between two spatially
    autocorrelated variables, applied to ranks (Spearman).

    Var(r) ~ (1/n^2) * [n + sum_k n_k rho_x(k) rho_y(k)], with n_k ordered pairs in distance class k
    (classes of `band_km` up to `max_km`, one class beyond) and rho(k) the Moran-type autocorrelation of each
    variable in class k. n_eff = 1 + 1/Var(r) (clipped to [3, n]); p from t = r sqrt((n_eff-2)/(1-r^2)).
    Reviewer addition: a national label shuffle ignores spatial autocorrelation, so its null is too narrow when
    both variables are spatially structured.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    ok = ~(np.isnan(x) | np.isnan(y) | np.isnan(lat) | np.isnan(lon))
    x, y, lat, lon = x[ok], y[ok], lat[ok], lon[ok]
    n = len(x)
    if n < 5:
        return {"n_eff_spatial": np.nan, "p_spatial_crh": np.nan}
    zx, zy = stats.rankdata(x), stats.rankdata(y)
    zx, zy = zx - zx.mean(), zy - zy.mean()
    d = haversine_km(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
    cls = np.minimum((d // band_km).astype(np.int64), int(max_km // band_km))
    np.fill_diagonal(cls, -1)
    m = cls >= 0
    c = cls[m]
    nk = np.bincount(c)
    sx = np.bincount(c, weights=np.outer(zx, zx)[m], minlength=len(nk))
    sy = np.bincount(c, weights=np.outer(zy, zy)[m], minlength=len(nk))
    with np.errstate(invalid="ignore", divide="ignore"):
        rx = np.where(nk > 0, (sx / nk) / (zx @ zx / n), 0.0)
        ry = np.where(nk > 0, (sy / nk) / (zy @ zy / n), 0.0)
    var_r = (n + float(np.sum(nk * rx * ry))) / n ** 2
    n_eff = float(np.clip(1 + 1 / var_r, 3, n)) if var_r > 0 else float(n)
    r = float(np.corrcoef(zx, zy)[0, 1])
    t = r * np.sqrt((n_eff - 2) / max(1e-12, 1 - r ** 2))
    return {"n_eff_spatial": n_eff, "p_spatial_crh": float(2 * stats.t.sf(abs(t), n_eff - 2))}


def cluster_bootstrap_spearman_ci(x, y, groups, n_boot: int = N_BOOT, seed: int = config.SEED, alpha: float = 0.05):
    """Percentile CI of Spearman rho resampling whole clusters (states) with replacement (reviewer addition:
    counties are spatially dependent, so an i.i.d. county bootstrap understates uncertainty)."""
    x, y, groups = np.asarray(x, float), np.asarray(y, float), np.asarray(groups)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y, groups = x[ok], y[ok], groups[ok]
    ug = np.unique(groups)
    if len(ug) < 5:
        return np.nan, np.nan
    members = [np.flatnonzero(groups == g) for g in ug]
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for b in range(n_boot):
        i = np.concatenate([members[j] for j in rng.integers(0, len(ug), len(ug))])
        rx, ry = stats.rankdata(x[i]), stats.rankdata(y[i])
        out[b] = np.corrcoef(rx, ry)[0, 1] if rx.std() > 0 and ry.std() > 0 else np.nan
    return float(np.nanpercentile(out, 100 * alpha / 2)), float(np.nanpercentile(out, 100 * (1 - alpha / 2)))


# ============================================================================ burden harmonisation
BURDEN_COLS = [
    "burden_row_id", "geo_id", "geo_level", "geo_name", "state_fips", "county_fips", "geo_vintage",
    "in_canonical_geographies", "condition_id", "measure_id", "measure_label", "metric_type", "value", "value_unit",
    "ci_low", "ci_high", "ci_level", "numerator", "denominator", "denominator_description", "period", "period_start",
    "period_end", "year", "stratum", "adjustment", "burden_evidence_level", "measure_role", "is_primary_measure",
    "inherited", "derivation", "proxy_for_condition", "suppressed", "estimate_kind", "geo_match_status",
    "source_table", "source_measure_id",
] + PROVENANCE_COLUMNS


def _geo_names() -> pd.DataFrame:
    g = read_table("geographies")
    return g


def _finish_burden(df: pd.DataFrame) -> pd.DataFrame:
    for c in BURDEN_COLS:
        if c not in df.columns:
            df[c] = None
    df["is_primary_measure"] = df["measure_role"].eq("primary")
    df["evidence_level"] = df["burden_evidence_level"]
    df["burden_row_id"] = (df["condition_id"] + "|" + df["derivation"].astype(str) + "|" + df["geo_level"] + ":"
                           + df["geo_id"].astype(str) + "|" + df["source_record_id"].astype(str))
    return df[BURDEN_COLS]


def _hps_rows() -> pd.DataFrame:
    lc = read_table("geo_condition_burden__cdc_long_covid")
    st_mask = (lc["measure_id"] == "lc_current_pct_all_adults") & (lc["subgroup_type"] == "By State")
    nat_mask = (lc["measure_id"] == "lc_current_pct_all_adults") & (lc["subgroup_type"] == "National Estimate")
    role = pd.Series("other", index=lc.index)
    # state primary: latest non-suppressed period per state
    st = lc[st_mask & lc["value"].notna()]
    prim_idx = st.sort_values("period_end").groupby("geo_id").tail(1).index
    role.loc[prim_idx] = "primary"
    latest = lc.loc[st_mask, "period_end"].max()
    strict = lc.index[st_mask & (lc["period_end"] == latest)]
    role.loc[[i for i in strict if role.loc[i] != "primary"]] = "alternate"
    nat = lc[nat_mask & lc["value"].notna()]
    role.loc[nat.sort_values("period_end").tail(1).index] = "primary"
    strat = np.where(lc["subgroup_type"].isin(["By State", "National Estimate"]), "all",
                     lc["subgroup_type"].astype(str) + "=" + lc["subgroup"].astype(str))
    out = pd.DataFrame({
        "geo_id": lc["geo_id"], "geo_level": lc["geo_level"], "geo_name": lc["geo_name"],
        "state_fips": lc["state_fips"], "county_fips": None, "geo_vintage": lc["geo_vintage"],
        "in_canonical_geographies": lc["in_canonical_geographies"], "condition_id": "long_covid",
        "measure_id": lc["measure_id"], "measure_label": lc["measure_label"], "metric_type": lc["metric_type"],
        "value": lc["value"], "value_unit": lc["value_unit"], "ci_low": lc["ci_low"], "ci_high": lc["ci_high"],
        "ci_level": lc["ci_level"], "numerator": np.nan, "denominator": np.nan,
        "denominator_description": lc["denominator_population"].astype(str) + "; " + lc["denominator_source"].astype(str),
        "period": lc["period_label"], "period_start": lc["period_start"], "period_end": lc["period_end"],
        "year": lc["year"].astype("float"), "stratum": strat, "adjustment": "survey_weighted_crude",
        "burden_evidence_level": lc["burden_evidence_level"], "measure_role": role, "inherited": False,
        "derivation": "as_published", "proxy_for_condition": None, "suppressed": lc["suppressed"],
        "estimate_kind": "survey",
        "geo_match_status": np.where(lc["geo_level"] == "national", "national",
                                     np.where(lc["in_canonical_geographies"], "matched_2024", "unmatched")),
        "source_table": "geo_condition_burden__cdc_long_covid", "source_measure_id": lc["measure_id"],
    })
    for c in PROVENANCE_COLUMNS:
        out[c] = lc[c]
    # alternate: mean of 2024 non-suppressed periods (state + national)
    pool = lc[(st_mask | nat_mask) & lc["value"].notna() & (lc["period_start"] >= HPS_POOLED_FROM)]
    agg = pool.groupby("geo_id").agg(value=("value", "mean"), n_periods=("value", "size"),
                                     period_start=("period_start", "min"), period_end=("period_end", "max"))
    first = pool.drop_duplicates("geo_id").set_index("geo_id")
    alt = pd.DataFrame({
        "geo_id": agg.index, "geo_level": first.loc[agg.index, "geo_level"].to_numpy(),
        "geo_name": first.loc[agg.index, "geo_name"].to_numpy(),
        "state_fips": first.loc[agg.index, "state_fips"].to_numpy(), "county_fips": None,
        "geo_vintage": first.loc[agg.index, "geo_vintage"].to_numpy(),
        "in_canonical_geographies": first.loc[agg.index, "in_canonical_geographies"].to_numpy(),
        "condition_id": "long_covid", "measure_id": "lc_current_pct_all_adults_2024_mean",
        "measure_label": "Currently experiencing long COVID, % of all adults: unweighted mean of 2024 HPS periods",
        "metric_type": "prevalence_pct", "value": agg["value"].to_numpy(), "value_unit": "percent",
        "ci_low": np.nan, "ci_high": np.nan, "ci_level": None, "numerator": np.nan, "denominator": np.nan,
        "denominator_description": "all adults 18+ (mean of " + agg["n_periods"].astype(str).to_numpy() + " periods)",
        "period": "2024 periods (" + agg["period_start"].to_numpy() + " .. " + agg["period_end"].to_numpy() + ")",
        "period_start": agg["period_start"].to_numpy(), "period_end": agg["period_end"].to_numpy(), "year": 2024.0,
        "stratum": "all", "adjustment": "survey_weighted_crude", "burden_evidence_level": "A",
        "measure_role": "alternate", "inherited": False, "derivation": "mean_of_periods",
        "proxy_for_condition": None, "suppressed": False, "estimate_kind": "survey",
        "geo_match_status": np.where(first.loc[agg.index, "geo_level"] == "national", "national", "matched_2024"),
        "source_table": "geo_condition_burden__cdc_long_covid", "source_measure_id": "lc_current_pct_all_adults",
    })
    for c in PROVENANCE_COLUMNS:
        alt[c] = first.loc[agg.index, c].to_numpy()
    alt["source_record_id"] = "hps2024mean:" + alt["geo_id"]
    alt["provenance_notes"] = ("Derived alternate: unweighted mean of the non-suppressed 2024 HPS period estimates "
                               "(no CI; periods are not independent samples of equal size).")
    return pd.concat([out, alt], ignore_index=True)


def _lyme_rows(geos: pd.DataFrame) -> pd.DataFrame:
    ly = read_table("geo_condition_burden__cdc_lyme")
    legacy = set(geos.loc[geos["ct_legacy"], "geo_id"])
    status = np.where(ly["in_canonical_geographies"], "matched_2024",
                      np.where(ly["geo_id"].isin(legacy), "matched_ct_legacy", "unmatched"))
    role = pd.Series("other", index=ly.index)
    canon = ly["in_canonical_geographies"].astype(bool)
    inc = ly["measure_id"] == "lyme_incidence_per_100k"
    role[canon & inc & (ly["year"] == PRIMARY_LYME_YEAR)] = "primary"
    role[canon & ~inc & (ly["year"] == PRIMARY_LYME_YEAR)] = "count_companion"
    role[canon & inc & (ly["year"] == ALT_LYME_PRE2022_YEAR)] = "alternate"
    out = pd.DataFrame({
        "geo_id": ly["geo_id"], "geo_level": ly["geo_level"], "geo_name": ly["geo_name"] + ", " + ly["state_name"],
        "state_fips": ly["state_fips"], "county_fips": ly["county_fips"], "geo_vintage": ly["geo_vintage"],
        "in_canonical_geographies": canon, "condition_id": "lyme_disease", "measure_id": ly["measure_id"],
        "measure_label": ly["measure_label"], "metric_type": ly["metric_type"], "value": ly["value"],
        "value_unit": ly["value_unit"], "ci_low": ly["ci_low"], "ci_high": ly["ci_high"], "ci_level": None,
        "numerator": ly["numerator"].astype("float"), "denominator": ly["denominator"].astype("float"),
        "denominator_description": ly["denominator_source"], "period": ly["period_label"],
        "period_start": ly["period_start"], "period_end": ly["period_end"], "year": ly["year"].astype(float),
        "stratum": "all", "adjustment": "crude", "burden_evidence_level": ly["burden_evidence_level"],
        "measure_role": role, "inherited": False,
        "derivation": np.where(inc, "computed_rate", "as_published"), "proxy_for_condition": None,
        "suppressed": ly["suppressed"], "estimate_kind": "surveillance", "geo_match_status": status,
        "source_table": "geo_condition_burden__cdc_lyme", "source_measure_id": ly["measure_id"],
    })
    for c in PROVENANCE_COLUMNS:
        out[c] = ly[c]
    prov = ly.iloc[0]

    # ---- derived: county pooled 2022-2023; state + national aggregates (2023, 2019, pooled 2022-2023)
    cases = ly[ly["measure_id"] == "lyme_reported_cases"][["geo_id", "state_fips", "year", "numerator",
                                                            "denominator", "in_canonical_geographies",
                                                            "state_name"]].copy()
    cases["numerator"] = cases["numerator"].astype(float)
    cases["denominator"] = cases["denominator"].astype(float)
    derived = []

    def _row(geo_id, level, name, state_fips, county_fips, mid, label, value, num, den, period, year, role_, deriv,
             note, metric="incidence_per_100k", unit="per 100,000 population"):
        return {"geo_id": geo_id, "geo_level": level, "geo_name": name, "state_fips": state_fips,
                "county_fips": county_fips, "geo_vintage": "2024" if level != "national" else "national",
                "in_canonical_geographies": level != "national", "condition_id": "lyme_disease",
                "measure_id": mid, "measure_label": label, "metric_type": metric, "value": value,
                "value_unit": unit, "ci_low": np.nan, "ci_high": np.nan, "ci_level": None, "numerator": num,
                "denominator": den, "denominator_description": "Census PEP July-1 population (ingestion)",
                "period": period, "period_start": f"{period[:4]}-01-01", "period_end": f"{period[-4:]}-12-31",
                "year": float(year), "stratum": "all", "adjustment": "crude", "burden_evidence_level": "A",
                "measure_role": role_, "inherited": False, "derivation": deriv, "proxy_for_condition": None,
                "suppressed": False, "estimate_kind": "surveillance",
                "geo_match_status": "national" if level == "national" else "matched_2024",
                "source_table": "geo_condition_burden__cdc_lyme", "source_measure_id": "lyme_reported_cases",
                "data_layer": "geographic", "source_name": prov["source_name"],
                "source_record_id": f"derived:lyme:{mid}:{level}:{geo_id}:{period}",
                "source_version": prov["source_version"], "retrieved_at": prov["retrieved_at"],
                "source_geographic_resolution": "county", "evidence_type": "surveillance_case_count",
                "provenance_notes": note}

    c = cases[cases["in_canonical_geographies"].astype(bool)]
    p = c[c["year"].isin([2022, 2023])].pivot_table(index="geo_id", columns="year",
                                                     values=["numerator", "denominator"], aggfunc="sum")
    p = p.dropna()
    names = out.drop_duplicates("geo_id").set_index("geo_id")
    for gid, r in p.iterrows():
        num = r[("numerator", 2022)] + r[("numerator", 2023)]
        den = r[("denominator", 2022)] + r[("denominator", 2023)]
        derived.append(_row(gid, "county", names.loc[gid, "geo_name"], gid[:2], gid,
                            "lyme_incidence_per_100k_2022_2023",
                            "Reported Lyme incidence per 100,000, pooled 2022-2023", num / den * 1e5, num, den,
                            "2022-2023", 2023, "alternate", "pooled_years",
                            "Derived alternate: (cases 2022 + 2023) / (population 2022 + 2023) x 1e5; only counties "
                            "with both years on the same 2024 FIPS (CT planning regions have 2023 only)."))

    def _agg(df, year_set):
        g = df[df["year"].isin(year_set)]
        withden = g[g["denominator"].notna()]
        num = withden.groupby("state_fips")["numerator"].sum()
        den = withden.groupby("state_fips")["denominator"].sum()
        miss = g[g["denominator"].isna()].groupby("state_fips")["numerator"].sum()
        years_present = g.groupby("state_fips")["year"].nunique()
        return num, den, miss, years_present

    snames = cases.drop_duplicates("state_fips").set_index("state_fips")["state_name"]
    specs = [((PRIMARY_LYME_YEAR,), str(PRIMARY_LYME_YEAR), "primary"),
             ((ALT_LYME_PRE2022_YEAR,), str(ALT_LYME_PRE2022_YEAR), "alternate"),
             ((2022, 2023), "2022-2023", "alternate")]
    for years, period, role_ in specs:
        num, den, miss, yp = _agg(cases, years)
        pooled = len(years) > 1
        mid = "lyme_incidence_per_100k_agg" + ("_2022_2023" if pooled else "")
        for sf in num.index:
            if sf in TERRITORY_FIPS or sf is None:
                continue
            if pooled and (yp.get(sf, 0) < 2 or miss.get(sf, 0) > 0 and sf == "09"):
                continue  # CT 2022 is on legacy counties without denominators
            note = (f"Aggregated from county case counts with a population denominator; {int(miss.get(sf, 0))} "
                    "cases in county-years without a denominator are excluded; cases with unknown county are "
                    "not in the county file (see state_cases_missing_from_county_file_pct in the source).")
            derived.append(_row(sf, "state", snames.get(sf), sf, None, mid,
                                f"Reported Lyme incidence per 100,000 (state aggregate of counties), {period}",
                                num[sf] / den[sf] * 1e5 if den[sf] > 0 else np.nan, num[sf], den[sf], period,
                                years[-1], role_, "aggregated_from_county", note))
            if not pooled and role_ == "primary":
                derived.append(_row(sf, "state", snames.get(sf), sf, None, "lyme_reported_cases_agg",
                                    f"Reported Lyme cases (state sum of counties), {period}", num[sf] + miss.get(sf, 0),
                                    num[sf] + miss.get(sf, 0), np.nan, period, years[-1], "count_companion",
                                    "aggregated_from_county", note, metric="case_count", unit="cases"))
        keep = [s for s in num.index if s in UNIVERSE_STATES]
        if pooled:
            keep = [s for s in keep if yp.get(s, 0) == 2 and s != "09"]
        nn, dd = num[keep].sum(), den[keep].sum()
        derived.append(_row("US", "national", "United States", None, None, mid,
                            f"Reported Lyme incidence per 100,000 (sum of counties, 50 states + DC), {period}",
                            nn / dd * 1e5, nn, dd, period, years[-1], role_, "aggregated_from_county",
                            "Aggregated from county case counts" + (" (CT excluded: 2022 on legacy counties)"
                                                                    if pooled else "")))
    return pd.concat([out, pd.DataFrame(derived)], ignore_index=True)


def _cms_rows(geos: pd.DataFrame) -> pd.DataFrame:
    cm = read_table("geo_condition_burden__cms_mmd")
    adj = np.where(cm["adjustment"] == "unsmoothed_actual", "actual", "agestd")
    age = np.where(cm["stratum_age"] == "all", "all", "lt65")
    mid = "cms_mmd_" + cm["source_condition_code"].astype(str) + "_prev_" + adj + "_" + age
    matched = cm["geo_match_status"].isin(["matched_2024", "national"])
    role = pd.Series("other", index=cm.index)
    is_actual_all = (adj == "actual") & (age == "all")
    role[matched & is_actual_all & (cm["year"] == PRIMARY_CMS_YEAR)] = "primary"
    role[matched & (adj == "agestd") & (age == "all") & (cm["year"] == PRIMARY_CMS_YEAR)] = "alternate"
    role[matched & is_actual_all & (cm["year"] == ALT_CMS_PRELIM_YEAR)] = "alternate"
    role[matched & (adj == "actual") & (age == "lt65") & (cm["year"] == PRIMARY_CMS_YEAR)] = "alternate"
    nm = geos.set_index("geo_id")
    gname = cm["geo_id"].map(nm["name"] + ", " + nm["state_name"]).where(cm["geo_level"] == "county",
                                                                        cm["geo_id"].map(nm["name"]))
    gname = gname.where(cm["geo_level"] != "national", "United States")
    vint = cm["geo_match_status"].map({"matched_2024": "2024", "matched_ct_legacy": "2020_ct_legacy",
                                       "unmatched": "cms_fips_not_in_2024", "national": "national"})
    out = pd.DataFrame({
        "geo_id": cm["geo_id"], "geo_level": cm["geo_level"], "geo_name": gname, "state_fips": cm["state_fips"],
        "county_fips": cm["county_fips"], "geo_vintage": vint, "in_canonical_geographies": cm["geo_match_status"] == "matched_2024",
        "condition_id": cm["condition_id"], "measure_id": mid,
        "measure_label": ("CMS MMD " + cm["source_condition_name"] + " prevalence, "
                          + np.where(adj == "actual", "unsmoothed actual", "unsmoothed age-standardized")
                          + np.where(age == "all", ", all ages", ", age < 65")),
        "metric_type": "prevalence_pct", "value": cm["prevalence_pct"].astype(float),
        "value_unit": "percent of Medicare FFS beneficiaries (integer as published)", "ci_low": np.nan,
        "ci_high": np.nan, "ci_level": None, "numerator": np.nan, "denominator": np.nan,
        "denominator_description": "Medicare FFS beneficiaries; size class " + cm["denominator_range"].astype(str),
        "period": "claims year " + cm["year"].astype(str) + " (" + cm["claims_data_status"].astype(str) + ")",
        "period_start": cm["year"].astype(str) + "-01-01", "period_end": cm["year"].astype(str) + "-12-31",
        "year": cm["year"].astype(float), "stratum": "age=" + cm["stratum_age"].astype(str),
        "adjustment": np.where(adj == "actual", "crude", "age_standardized"),
        "burden_evidence_level": cm["burden_evidence_level"], "measure_role": role, "inherited": False,
        "derivation": "as_published", "proxy_for_condition": np.where(cm["burden_evidence_level"] == "C",
                                                                      cm["condition_id"], None),
        "suppressed": cm["value_zero_possibly_suppressed"], "estimate_kind": "claims",
        "geo_match_status": cm["geo_match_status"], "source_table": "geo_condition_burden__cms_mmd",
        "source_measure_id": "mmd_condition_" + cm["source_condition_code"].astype(str) + "_prevalence",
    })
    for c in PROVENANCE_COLUMNS:
        out[c] = cm[c]
    return out


def _places_rows() -> pd.DataFrame:
    pl = read_table("geo_context__cdc_places")
    meas = sorted({*PLACES_PROXY.values(), *[m for v in PLACES_ALTERNATES.values() for m in v]})
    pl = pl[(pl["geo_level"] == "county") & pl["measure_id"].isin(meas)]
    frames = []
    pairs = [(c, m, "county_proxy") for c, m in PLACES_PROXY.items()]
    pairs += [(c, m, "alternate") for c, ms in PLACES_ALTERNATES.items() for m in ms]
    for cond, m, role in pairs:
        d = pl[pl["measure_id"] == m]
        crude = d["data_value_type_id"] == "CrdPrv"
        r = np.where(crude, role, "alternate")
        f = pd.DataFrame({
            "geo_id": d["geo_id"], "geo_level": "county", "geo_name": d["geo_name"] + ", " + d["state_name"],
            "state_fips": d["state_fips"], "county_fips": d["county_fips"], "geo_vintage": d["geo_vintage"],
            "in_canonical_geographies": d["in_canonical_geographies"], "condition_id": cond,
            "measure_id": "places_" + m + np.where(crude, "_crude", "_ageadj"),
            "measure_label": "PLACES " + d["measure_label"] + np.where(crude, " (crude)", " (age-adjusted)")
                             + f" - level C proxy for {cond}",
            "metric_type": "modeled_prevalence_pct", "value": d["value"], "value_unit": d["value_unit"],
            "ci_low": d["ci_low"], "ci_high": d["ci_high"], "ci_level": "95% (model)", "numerator": np.nan,
            "denominator": d["total_pop_18plus"].astype(float),
            "denominator_description": "PLACES total_pop_18plus (Census population used by PLACES)",
            "period": "BRFSS " + d["brfss_year"].astype(str) + " (PLACES " + d["places_release"].astype(str) + ")",
            "period_start": d["brfss_year"].astype(str) + "-01-01", "period_end": d["brfss_year"].astype(str) + "-12-31",
            "year": d["brfss_year"].astype(float), "stratum": "adults 18+",
            "adjustment": np.where(crude, "crude", "age_adjusted"), "burden_evidence_level": "C",
            "measure_role": r, "inherited": False, "derivation": "copied_as_proxy", "proxy_for_condition": cond,
            "suppressed": d["suppressed"], "estimate_kind": "modeled",
            "geo_match_status": np.where(d["in_canonical_geographies"], "matched_2024", "unmatched"),
            "source_table": "geo_context__cdc_places", "source_measure_id": d["measure_id"],
        })
        for c in PROVENANCE_COLUMNS:
            f[c] = d[c].to_numpy()
        f["evidence_level"] = "C"
        f["provenance_notes"] = (f"Level C proxy for {cond}: {PLACES_PROXY_WHY[m]}. Modeled small-area estimate "
                                 "(not observed prevalence; not a measure of the condition).")
        frames.append(f)
    out = pd.concat(frames, ignore_index=True)
    # state aggregates of the crude county proxies (alternate; Test 6 R5)
    cr = out[(out["measure_role"] == "county_proxy") & out["value"].notna()].copy()
    cr["w"] = cr["denominator"]
    cr["vw"] = cr["value"] * cr["w"]
    agg = cr.groupby(["condition_id", "measure_id", "state_fips"]).agg(
        vw=("vw", "sum"), w=("w", "sum"), n=("value", "size"), year=("year", "max"),
        source_record_id=("source_record_id", "first"), source_version=("source_version", "first"),
        retrieved_at=("retrieved_at", "first"), source_name=("source_name", "first")).reset_index()
    geos = read_table("geographies").query("geo_level == 'state'").set_index("geo_id")
    st = pd.DataFrame({
        "geo_id": agg["state_fips"], "geo_level": "state", "geo_name": agg["state_fips"].map(geos["name"]),
        "state_fips": agg["state_fips"], "county_fips": None, "geo_vintage": "2024",
        "in_canonical_geographies": True, "condition_id": agg["condition_id"],
        "measure_id": agg["measure_id"] + "_state_agg",
        "measure_label": "PLACES " + agg["measure_id"].str.replace("places_", "") + " state aggregate "
                         "(pop 18+-weighted mean of county crude estimates) - level C proxy",
        "metric_type": "modeled_prevalence_pct", "value": agg["vw"] / agg["w"], "value_unit": "%",
        "ci_low": np.nan, "ci_high": np.nan, "ci_level": None, "numerator": np.nan, "denominator": agg["w"],
        "denominator_description": "sum of county PLACES total_pop_18plus (" + agg["n"].astype(str) + " counties)",
        "period": "BRFSS " + agg["year"].astype(int).astype(str), "period_start": None, "period_end": None,
        "year": agg["year"], "stratum": "adults 18+", "adjustment": "crude", "burden_evidence_level": "C",
        "measure_role": "alternate", "inherited": False, "derivation": "aggregated_from_county",
        "proxy_for_condition": agg["condition_id"], "suppressed": False, "estimate_kind": "modeled",
        "geo_match_status": "matched_2024", "source_table": "geo_context__cdc_places",
        "source_measure_id": agg["measure_id"], "data_layer": "geographic", "source_name": agg["source_name"],
        "source_record_id": "stateagg:" + agg["measure_id"] + ":" + agg["state_fips"],
        "source_version": agg["source_version"], "retrieved_at": agg["retrieved_at"],
        "source_geographic_resolution": "county", "evidence_type": "modeled_small_area_estimate",
        "evidence_level": "C",
        "provenance_notes": "Derived: population(18+)-weighted mean of PLACES county crude estimates; used only for the "
                            "state-level proxy check. Not an official PLACES state estimate.",
    })
    return pd.concat([out, st], ignore_index=True)


def _ptlds_rows(lyme: pd.DataFrame) -> pd.DataFrame:
    src = lyme[lyme["measure_role"].isin(["primary", "alternate"])
               & lyme["metric_type"].eq("incidence_per_100k")].copy()
    src["source_measure_id"] = src["measure_id"]
    src["condition_id"] = "ptlds"
    src["burden_evidence_level"] = "C"
    src["evidence_level"] = "C"
    src["derivation"] = "copied_as_proxy"
    src["proxy_for_condition"] = "ptlds"
    src["measure_label"] = src["measure_label"] + " - level C antecedent-exposure proxy for PTLDS"
    src["provenance_notes"] = ("Level C proxy for PTLDS: Lyme incidence (antecedent exposure). PTLDS has no "
                               "surveillance or survey estimate; this is not PTLDS prevalence.")
    return src


def _inherited_rows(hps: pd.DataFrame, geos: pd.DataFrame) -> pd.DataFrame:
    prim = hps[(hps["geo_level"] == "state") & (hps["measure_role"] == "primary")].set_index("geo_id")
    cty = geos[(geos["geo_level"] == "county") & ~geos["ct_legacy"] & geos["state_fips"].isin(prim.index)]
    rows = prim.loc[cty["state_fips"]].reset_index(drop=True)
    rows["source_record_id"] = rows["source_record_id"].to_numpy()
    rows["geo_id"] = cty["geo_id"].to_numpy()
    rows["geo_level"] = "county"
    rows["geo_name"] = (cty["name"] + ", " + cty["state_name"]).to_numpy()
    rows["county_fips"] = cty["geo_id"].to_numpy()
    rows["geo_vintage"] = "2024"
    rows["inherited"] = True
    rows["derivation"] = "inherited_from_state"
    rows["source_geographic_resolution"] = "state"
    rows["geo_match_status"] = "matched_2024"
    rows["measure_label"] = rows["measure_label"] + " [STATE value inherited as county context; not county prevalence]"
    rows["provenance_notes"] = ("Inherited state context (CONVENTIONS 3.5): the value is the state's HPS estimate; "
                                "there is no county long-COVID measure. Never read as county prevalence.")
    return rows


SAE_TABLE = "geo_condition_burden__brfss_long_covid_sae"


def _long_covid_source() -> tuple[str, str]:
    import yaml
    from ..config import CONFIGS
    cfg = yaml.safe_load((CONFIGS / "scoring.yaml").read_text())
    return cfg.get("long_covid_burden", "hps_inherited"), cfg.get("long_covid_sae_measure",
                                                                  "lc_current_pct_all_adults_mmwr_def")


def _sae_rows(measure: str) -> pd.DataFrame:
    """BRFSS 2023 MRP small-area estimates (geography.sae) as Long COVID burden rows; the configured measure is
    primary at county, state and national level, the others 'other'. Modelled, not observed: estimate_kind and the
    provenance notes say so (results/BRFSS_SAE_RESULTS.md)."""
    sa = read_table(SAE_TABLE)
    out = pd.DataFrame({
        "geo_id": sa["geo_id"], "geo_level": sa["geo_level"],
        "geo_name": np.where(sa["geo_level"] == "county", sa["geo_name"] + ", " + sa["state_name"].astype(str),
                             sa["geo_name"]),
        "state_fips": sa["state_fips"], "county_fips": sa["county_fips"], "geo_vintage": sa["geo_vintage"],
        "in_canonical_geographies": sa["in_canonical_geographies"], "condition_id": "long_covid",
        "measure_id": sa["measure_id"], "measure_label": sa["measure_label"] + ", % of all adults",
        "metric_type": sa["metric_type"],
        "value": sa["value"], "value_unit": sa["value_unit"], "ci_low": sa["ci_low"], "ci_high": sa["ci_high"],
        "ci_level": sa["ci_level"], "numerator": sa["numerator"], "denominator": sa["denominator"],
        "denominator_description": sa["denominator_population"].astype(str), "period": sa["period_label"],
        "period_start": sa["period_start"], "period_end": sa["period_end"], "year": sa["year"].astype(float),
        "stratum": "all", "adjustment": "mrp_post_stratified", "burden_evidence_level": sa["burden_evidence_level"],
        "measure_role": np.where(sa["measure_id"] == measure, "primary", "other"), "inherited": False,
        "derivation": sa["derivation"], "proxy_for_condition": None, "suppressed": sa["suppressed"],
        "estimate_kind": "modeled_small_area",
        "geo_match_status": np.where(sa["geo_level"] == "national", "national",
                                     np.where(sa["in_canonical_geographies"], "matched_2024", "unmatched")),
        "source_table": SAE_TABLE, "source_measure_id": sa["measure_id"],
    })
    for c in PROVENANCE_COLUMNS:
        out[c] = sa[c]
    return out


def build_burden() -> pd.DataFrame:
    geos = _geo_names()
    hps = _hps_rows()
    lyme = _lyme_rows(geos)
    inherited = _inherited_rows(hps, geos)
    source, measure = _long_covid_source()
    lc_parts = []
    if source == "brfss_sae" and table_exists(SAE_TABLE):
        # the SAE replaces the HPS value as primary; HPS state rows and the inherited county rows stay as alternates
        hps = hps.assign(measure_role=hps["measure_role"].replace({"primary": "alternate"}))
        inherited = inherited.assign(measure_role="alternate")
        lc_parts = [_sae_rows(measure)]
    elif source == "brfss_sae":
        print(f"[geography] {SAE_TABLE} not built: Long COVID burden falls back to the inherited HPS state value")
    parts = [hps, lyme, _ptlds_rows(lyme), _cms_rows(geos), _places_rows(), inherited] + lc_parts
    df = pd.concat(parts, ignore_index=True)
    df = _finish_burden(df)
    for c in ("year", "numerator", "denominator", "value", "ci_low", "ci_high"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("inherited", "suppressed", "in_canonical_geographies"):
        df[c] = df[c].fillna(False).astype(bool)
    for c in ("state_fips", "county_fips", "geo_name", "proxy_for_condition", "ci_level", "period_start",
              "period_end", "geo_vintage"):
        df[c] = df[c].astype("string")
    return df.reset_index(drop=True)


# ============================================================================ context
ACS_COLS = ["total_population", "median_age", "pct_age_18_64", "pct_age_65_plus", "pct_with_disability",
            "pct_uninsured", "pct_below_poverty", "median_household_income", "pct_households_no_vehicle",
            "pct_households_broadband", "pop_density_per_km2"]
HPSA_COLS = ["pc_any_designation", "pc_whole_county", "pc_partial_county", "pc_coverage", "pc_max_score",
             "mh_any_designation", "mh_coverage", "dh_any_designation", "dh_coverage", "mua_whole_county",
             "mua_partial_county", "mua_coverage", "mua_min_imu_score"]
SVI_MAP = {"rpl_themes": "svi_overall", "rpl_theme1": "svi_socioeconomic", "rpl_theme2": "svi_household",
           "rpl_theme3": "svi_minority", "rpl_theme4": "svi_housing_transport"}


def _svi_county() -> pd.DataFrame:
    us = read_table("geo_vulnerability")
    pr = read_table("geo_vulnerability_puerto_rico")
    frames = []
    for df, uni in ((us, "US counties (50 states + DC)"), (pr, "Puerto Rico municipios only")):
        f = df[["geo_id", "e_totpop", *SVI_MAP]].rename(columns=SVI_MAP)
        f["svi_ranking_universe"] = uni
        f["svi_release"] = df["svi_release"].iloc[0] if "svi_release" in df else "SVI 2022"
        frames.append(f)
    out = pd.concat(frames, ignore_index=True)
    out["svi_derivation"] = "published_county_percentile"
    return out


def _svi_state(svi: pd.DataFrame) -> pd.DataFrame:
    s = svi.copy()
    s["state_fips"] = s["geo_id"].str[:2]
    cols = list(SVI_MAP.values())
    w = s["e_totpop"].astype(float)
    agg = {c: (s[c] * w).groupby(s["state_fips"]).sum() / w.groupby(s["state_fips"]).sum() for c in cols}
    out = pd.DataFrame(agg).reset_index().rename(columns={"state_fips": "geo_id"})
    out["svi_ranking_universe"] = np.where(out["geo_id"] == "72", "Puerto Rico municipios only",
                                           "US counties (50 states + DC)")
    out["svi_release"] = svi["svi_release"].iloc[0]
    out["svi_derivation"] = "state_popweighted_mean_of_county_percentiles (derived; no official state SVI)"
    return out


def _acs() -> tuple[pd.DataFrame, pd.DataFrame]:
    acs = read_table("geo_context__acs")
    acs["n_adults_18plus"] = acs["n_age_18_64"] + acs["n_age_65_plus"]
    moe = [c + "_moe" for c in ACS_COLS if c + "_moe" in acs.columns]
    keep = ["geo_id", "acs_vintage", "n_adults_18plus", *ACS_COLS, *moe]
    z = read_table("geo_context_zcta__acs")
    zkeep = ["zcta", "acs_vintage"] + [c for c in ACS_COLS + moe if c in z.columns]
    return acs[keep], z[zkeep]


def _places_wide() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pl = read_table("geo_context__cdc_places")
    pl = pl[pl["relevant_measure"].astype(bool) & (pl["geo_level"] == "county")]
    pl = pl.assign(col="places_" + pl["measure_id"] + np.where(pl["data_value_type_id"] == "CrdPrv", "_crude", "_ageadj"))
    wide = pl.pivot_table(index="geo_id", columns="col", values="value", aggfunc="first").reset_index()
    # state aggregates of crude values (pop 18+-weighted), labelled
    cr = pl[pl["data_value_type_id"] == "CrdPrv"].copy()
    cr["w"] = cr["total_pop_18plus"].astype(float)
    cr = cr[cr["value"].notna()]
    num = (cr["value"] * cr["w"]).groupby([cr["state_fips"], cr["col"]]).sum()
    den = cr["w"].groupby([cr["state_fips"], cr["col"]]).sum()
    st = (num / den).unstack().reset_index().rename(columns={"state_fips": "geo_id"})
    zp = read_table("geo_context_zcta__cdc_places")
    zp = zp[zp["relevant_measure"].astype(bool)]
    zw = zp.assign(col="places_" + zp["measure_id"] + "_crude").pivot_table(
        index="zcta", columns="col", values="value", aggfunc="first").reset_index()
    return wide, st, zw


def build_context() -> pd.DataFrame:
    geos = _geo_names()
    base = geos[~geos["ct_legacy"]][["geo_id", "geo_level", "name", "state_fips", "state_abbr", "state_name",
                                     "county_fips", "lat", "lon", "aland_m2", "vintage"]].copy()
    acs, zacs = _acs()
    svi_c = _svi_county()
    svi_s = _svi_state(svi_c)
    rucc = read_table("geo_context__rucc")[["geo_id", "rucc_2023", "rucc_description", "metro_status", "is_metro",
                                            "nonmetro_adjacent_to_metro"]]
    hp = read_table("geo_context__hpsa")[["geo_id", *HPSA_COLS]].rename(columns={c: "hpsa_" + c for c in HPSA_COLS})
    pw, pst, pz = _places_wide()
    cty = base[base["geo_level"] == "county"].merge(acs, on="geo_id", how="left") \
        .merge(svi_c.drop(columns="e_totpop"), on="geo_id", how="left").merge(rucc, on="geo_id", how="left") \
        .merge(hp, on="geo_id", how="left").merge(pw, on="geo_id", how="left")
    cty["places_derivation"] = np.where(cty["geo_id"].isin(pw["geo_id"]), "PLACES 2025 county estimates", None)
    sta = base[base["geo_level"] == "state"].merge(acs, on="geo_id", how="left") \
        .merge(svi_s, on="geo_id", how="left").merge(pst, on="geo_id", how="left")
    sta["places_derivation"] = np.where(sta["geo_id"].isin(pst["geo_id"]),
                                        "pop18+-weighted mean of PLACES county crude estimates (derived)", None)
    zc = read_table("zcta_centroids")[["zcta", "lat", "lon", "aland_m2", "county_fips", "state_fips", "state_abbr"]]
    zz = zc.merge(zacs, on="zcta", how="left").merge(pz, on="zcta", how="left")
    zz = zz.rename(columns={"zcta": "geo_id"})
    zz["geo_level"] = "zcta"
    zz["name"] = "ZCTA " + zz["geo_id"]
    zz["vintage"] = "2020 ZCTA (2024 Gazetteer)"
    zz["places_derivation"] = np.where(zz["places_PHLTH_crude"].notna() if "places_PHLTH_crude" in zz else False,
                                       "PLACES 2025 ZCTA crude estimates", None)
    ctx = pd.concat([sta, cty, zz], ignore_index=True, sort=False)
    ctx["object_id"] = np.where(ctx["geo_level"] == "zcta", "geo:zcta:" + ctx["geo_id"], "geo:" + ctx["geo_id"])
    versions = {t: read_table(t, columns=["source_version", "retrieved_at"]).iloc[0]
                for t in ("geo_context__acs", "geo_vulnerability", "geo_context__rucc", "geo_context__hpsa",
                          "geo_context__cdc_places")}
    ctx["data_layer"] = "geographic"
    ctx["source_name"] = ("Composite context: Census ACS 2020-2024 5-yr; CDC/ATSDR SVI 2022; USDA ERS RUCC 2023; "
                          "HRSA HPSA/MUA; CDC PLACES 2025")
    ctx["source_record_id"] = ctx["object_id"]
    ctx["source_version"] = " | ".join(f"{t}: {v['source_version']}" for t, v in versions.items())
    ctx["retrieved_at"] = max(v["retrieved_at"] for v in versions.values())
    ctx["source_geographic_resolution"] = ctx["geo_level"].map({"state": "state", "county": "county", "zcta": "zcta"})
    ctx["evidence_type"] = "composite_index"
    ctx["evidence_level"] = "population_context"
    ctx["provenance_notes"] = ("Wide join of context sources on FIPS/ZCTA (ecological). Not disease burden. State SVI "
                               "and state PLACES values are derived population-weighted means of county values; ZCTA "
                               "rows carry ACS + PLACES only. PR SVI is ranked within PR only.")
    return ctx


# ============================================================================ facility-layer counts
@dataclass
class Counts:
    county: pd.DataFrame   # index (condition_id, geo_id)
    state: pd.DataFrame


def _condition_core_groups() -> dict[str, list[str]]:
    rel = load_config("relevance")["condition_specialties"]
    return {c: list(v["core"]) for c, v in rel.items()}


def provider_counts(counties: pd.DataFrame, zpop: pd.DataFrame) -> tuple[Counts, dict]:
    groups = _condition_core_groups()
    psg = read_table("provider_specialty_groups", columns=["npi", "specialty_group", "entity_type"])
    psg = psg[psg["entity_type"].astype(str) == "1"]
    prov = read_table("providers", columns=["npi", "state_fips", "county_fips", "lat", "lon"]).set_index("npi")
    cty_rows, st_rows = [], []
    # radius structures
    pts = prov[prov["lat"].notna()][["lat", "lon"]].drop_duplicates().reset_index(drop=True)
    pts["pt"] = np.arange(len(pts))
    a_prov = radius_matrix(counties["lat"].to_numpy(), counties["lon"].to_numpy(), pts["lat"].to_numpy(),
                           pts["lon"].to_numpy(), PROVIDER_RADIUS_KM)
    zp = zpop.dropna(subset=["lat", "lon"])
    a_pop = radius_matrix(counties["lat"].to_numpy(), counties["lon"].to_numpy(), zp["lat"].to_numpy(),
                          zp["lon"].to_numpy(), PROVIDER_RADIUS_KM)
    pop50 = a_pop @ zp["total_population"].fillna(0).to_numpy()
    prov_pt = prov.reset_index().merge(pts, on=["lat", "lon"], how="left").set_index("npi")["pt"]
    meta = {}
    for cond, core in groups.items():
        for variant, gset in (("all", core), ("specialists", [g for g in core if g != "primary_care"])):
            npis = pd.Index(psg.loc[psg["specialty_group"].isin(gset), "npi"].unique())
            p = prov.reindex(npis)
            c_cnt = p.groupby("county_fips").size()
            s_cnt = p.groupby("state_fips").size()
            col = "relevant_providers_n" if variant == "all" else "relevant_specialists_n"
            cty_rows.append(pd.DataFrame({"condition_id": cond, "geo_id": c_cnt.index, "col": col, "v": c_cnt.values}))
            st_rows.append(pd.DataFrame({"condition_id": cond, "geo_id": s_cnt.index, "col": col, "v": s_cnt.values}))
            if variant == "all":
                ptc = prov_pt.reindex(npis).dropna().astype(int).value_counts()
                vec = np.zeros(len(pts))
                vec[ptc.index.to_numpy()] = ptc.to_numpy()
                within = a_prov @ vec
                cty_rows.append(pd.DataFrame({"condition_id": cond, "geo_id": counties["geo_id"].to_numpy(),
                                              "col": f"relevant_providers_within_{PROVIDER_RADIUS_KM}km_n", "v": within}))
                meta[cond] = {"groups": gset, "n_npis": int(len(npis)),
                              "n_geocoded": int(p["county_fips"].notna().sum())}
            else:
                meta[cond]["specialist_groups"] = gset
    cty = pd.concat(cty_rows).pivot_table(index=["condition_id", "geo_id"], columns="col", values="v", aggfunc="sum")
    st = pd.concat(st_rows).pivot_table(index=["condition_id", "geo_id"], columns="col", values="v", aggfunc="sum")
    popdf = pd.DataFrame({"geo_id": counties["geo_id"].to_numpy(), f"population_within_{PROVIDER_RADIUS_KM}km": pop50})
    return Counts(cty, st), {"meta": meta, "pop_within": popdf}


def hrsa_counts() -> tuple[pd.Series, pd.Series]:
    h = read_table("facilities__hrsa", columns=["site_bphc_number", "is_active", "is_service_delivery_site",
                                                "county_fips", "state_fips"])
    h = h[h["is_active"].astype(bool) & h["is_service_delivery_site"].astype(bool)]
    return h.groupby("county_fips")["site_bphc_number"].nunique(), h.groupby("state_fips")["site_bphc_number"].nunique()


def trial_counts(counties: pd.DataFrame) -> Counts:
    tc = read_table("trial_conditions", columns=["nct_id", "condition_id", "condition_literal_match"])
    tc = tc[tc["condition_literal_match"].astype(bool)].drop_duplicates(["nct_id", "condition_id"])
    ct = read_table("clinical_trials", columns=["nct_id", "overall_status"]).set_index("nct_id")["overall_status"]
    ts = read_table("trial_sites", columns=["nct_id", "country_group", "county_fips", "state_fips", "lat", "lon"])
    ts = ts[ts["country_group"] == "US"]
    ts = ts.assign(status=ts["nct_id"].map(ct))
    cty_rows, st_rows = [], []
    pts = ts[ts["lat"].notna()][["lat", "lon"]].drop_duplicates().reset_index(drop=True)
    pts["pt"] = np.arange(len(pts))
    tsp = ts.merge(pts, on=["lat", "lon"], how="left")
    amats = {r: radius_matrix(counties["lat"].to_numpy(), counties["lon"].to_numpy(), pts["lat"].to_numpy(),
                              pts["lon"].to_numpy(), r) for r in RADII_KM}
    cidx = {g: i for i, g in enumerate(counties["geo_id"])}
    nc = len(counties)
    for cond, g in tc.groupby("condition_id"):
        ids = set(g["nct_id"])
        s = tsp[tsp["nct_id"].isin(ids)]
        subsets = {"": s, "_recruiting": s[s["status"] == "RECRUITING"],
                   "_active": s[s["status"].isin(ACTIVE_STATUSES)], "_completed": s[s["status"] == "COMPLETED"]}
        for suf, sub in subsets.items():
            c = sub.groupby("county_fips")["nct_id"].nunique()
            cty_rows.append(pd.DataFrame({"condition_id": cond, "geo_id": c.index, "col": f"trials_in_geo{suf}_n",
                                          "v": c.values}))
            sc = sub.groupby("state_fips")["nct_id"].nunique()
            st_rows.append(pd.DataFrame({"condition_id": cond, "geo_id": sc.index, "col": f"trials_in_geo{suf}_n",
                                         "v": sc.values}))
        for suf, sub in (("", s), ("_recruiting", subsets["_recruiting"])):
            sub = sub.reset_index(drop=True)
            if sub.empty:
                continue
            # site-level incidence: county x site (radius from the county internal point; membership) and site x trial
            tidx = {t: i for i, t in enumerate(sorted(sub["nct_id"].unique()))}
            b = sparse.csr_matrix((np.ones(len(sub)), (np.arange(len(sub)), sub["nct_id"].map(tidx).to_numpy())),
                                  shape=(len(sub), len(tidx)))
            has_pt = sub["pt"].notna().to_numpy()
            site_pt = sub["pt"].fillna(0).astype(int).to_numpy()
            mem_r = sub["county_fips"].map(cidx)
            ok = mem_r.notna().to_numpy()
            member = sparse.csr_matrix((np.ones(ok.sum()), (mem_r[ok].astype(int).to_numpy(), np.flatnonzero(ok))),
                                       shape=(nc, len(sub)))
            for r, a in amats.items():
                rad = a[:, site_pt].multiply(has_pt[None, :]).tocsr()
                n = np.asarray(((rad @ b) > 0).sum(axis=1)).ravel()
                cty_rows.append(pd.DataFrame({"condition_id": cond, "geo_id": counties["geo_id"].to_numpy(),
                                              "col": f"trials_within_{r}km{suf}_n", "v": n}))
                uni = ((rad + member) > 0).astype(float)
                n = np.asarray(((uni @ b) > 0).sum(axis=1)).ravel()
                cty_rows.append(pd.DataFrame({"condition_id": cond, "geo_id": counties["geo_id"].to_numpy(),
                                              "col": f"trials_in_geo_or_within_{r}km{suf}_n", "v": n}))
    cty = pd.concat(cty_rows).pivot_table(index=["condition_id", "geo_id"], columns="col", values="v", aggfunc="sum")
    st = pd.concat(st_rows).pivot_table(index=["condition_id", "geo_id"], columns="col", values="v", aggfunc="sum")
    return Counts(cty, st)


def nih_counts() -> Counts:
    links = read_table("nih_project_conditions", columns=["appl_id", "condition_id", "match_tier",
                                                          "likely_false_positive"])
    links = links[(links["match_tier"] == "title_abstract") & ~links["likely_false_positive"].fillna(False).astype(bool)]
    links = links.drop_duplicates(["appl_id", "condition_id"])
    pj = read_table("nih_projects", columns=["appl_id", "core_project_num", "org_ipf_code", "org_name", "county_fips",
                                             "state_fips", "is_active"])
    pj["org_key"] = pj["org_ipf_code"].fillna("name:" + pj["org_name"].fillna(""))
    d = links.merge(pj, on="appl_id")
    d["active_core"] = d["core_project_num"].where(d["is_active"].fillna(False).astype(bool))
    out = {}
    for level, col in (("county", "county_fips"), ("state", "state_fips")):
        g = d.dropna(subset=[col]).groupby(["condition_id", col])
        f = pd.DataFrame({"nih_projects_n": g["appl_id"].nunique(), "nih_core_projects_n": g["core_project_num"].nunique(),
                          "nih_active_core_projects_n": g["active_core"].nunique(), "nih_orgs_n": g["org_key"].nunique()})
        f.index = f.index.set_names(["condition_id", "geo_id"])
        out[level] = f
    return Counts(out["county"], out["state"])


# ============================================================================ features
def _primary_burden(burden: pd.DataFrame) -> pd.DataFrame:
    p = burden[burden["measure_role"] == "primary"].copy()
    return p.drop_duplicates(["condition_id", "geo_level", "geo_id"])


def _unknown_reason_for_missing(row, defs) -> str:
    return UNKNOWN


def not_published_reasons(rows: pd.DataFrame, burden: pd.DataFrame) -> list[str]:
    """Why a geography has no value in a condition's defined primary measure, from what the source holds for it.

    rows: geo_condition_features rows (geo_id, geo_level, burden_measure_id). CMS MMD measures distinguish a county
    absent from every MMD file ingested (2012-2023, prevalence and context measures) from a county that CMS publishes
    in other years or for other cells; CDC Lyme distinguishes a county code absent from the county file from an
    omitted year. Nothing is inferred beyond what the tables hold (1:1 FIPS renames are bridged at ingestion, so a
    renamed county is not 'absent')."""
    if rows.empty:
        return []
    cms = burden[burden["source_table"] == "geo_condition_burden__cms_mmd"]
    cms_years = cms.groupby(["geo_level", "geo_id"])["year"].agg(lambda s: sorted({int(y) for y in s.dropna()}))
    try:
        ctx_geos = set(read_table("geo_context__cms_mmd", columns=["geo_id"])["geo_id"])
    except FileNotFoundError:
        ctx_geos = set()
    prim = cms[cms["year"] == PRIMARY_CMS_YEAR]
    prim_cells = set(zip(prim["geo_level"], prim["geo_id"], prim["measure_id"]))
    lyme = burden[burden["source_table"] == "geo_condition_burden__cdc_lyme"]
    lyme_years = lyme.groupby("geo_id")["year"].agg(lambda s: sorted({int(y) for y in s.dropna()}))
    out = []
    for r in rows.itertuples(index=False):
        mid, lvl, gid = str(r.burden_measure_id), str(r.geo_level), str(r.geo_id)
        if mid.startswith("cms_mmd_"):
            yrs = cms_years.get((lvl, gid))
            if yrs is None and gid not in ctx_geos:
                out.append(f"CMS MMD publishes no row for this {lvl} code in any claims year (2012-2023) or any "
                           "ingested measure: the entity is absent from the MMD files (CMS omits cells with fewer "
                           "than 11 beneficiaries and does not document other omissions).")
            elif (lvl, gid, mid) not in prim_cells and yrs and PRIMARY_CMS_YEAR not in yrs:
                out.append(f"CMS MMD publishes no {PRIMARY_CMS_YEAR} prevalence cell for this {lvl} (its prevalence "
                           f"rows cover {yrs[0]}-{yrs[-1]}, {len(yrs)} years).")
            else:
                out.append(f"The CMS MMD {PRIMARY_CMS_YEAR} cell for this measure is omitted for this {lvl} while "
                           "other MMD cells for it are published (CMS omits cells with fewer than 11 beneficiaries).")
        elif mid.startswith("lyme"):
            yrs = lyme_years.get(gid)
            out.append("The CDC Lyme county file has no row for this county code in any year." if yrs is None else
                       f"CDC Lyme publishes no {PRIMARY_LYME_YEAR} value for this county (rows for "
                       f"{yrs[0]}-{yrs[-1]}).")
        else:
            out.append(f"No row for this geography in the primary measure {mid} (not published by the source).")
    return out


def build_features(burden: pd.DataFrame, context: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    geos = _geo_names()
    base = geos[~geos["ct_legacy"]][["geo_id", "geo_level", "name", "state_fips", "state_abbr", "state_name", "lat", "lon"]]
    conds = target_conditions()
    df = base.merge(pd.DataFrame({"condition_id": conds}), how="cross")
    df["geo_name"] = np.where(df["geo_level"] == "county", df["name"] + ", " + df["state_name"], df["name"])
    df["county_fips"] = np.where(df["geo_level"] == "county", df["geo_id"], None)
    df["in_analysis_universe"] = df["state_fips"].isin(UNIVERSE_STATES)
    ctx = context[context["geo_level"].isin(["state", "county"])].set_index("geo_id")
    df["population_total"] = df["geo_id"].map(ctx["total_population"])
    df["population_adults_18plus"] = df["geo_id"].map(ctx["n_adults_18plus"])
    df["small_population_flag"] = df["population_total"] < load_config("scoring")["ranking"]["min_population"]
    # ---- burden
    defs = definitions_frame().set_index(["condition_id", "geo_level"])
    prim = _primary_burden(burden).set_index(["condition_id", "geo_level", "geo_id"])
    key = pd.MultiIndex.from_frame(df[["condition_id", "geo_level", "geo_id"]])
    pr = prim.reindex(key)
    dkey = pd.MultiIndex.from_frame(df[["condition_id", "geo_level"]])
    dd = defs.reindex(dkey)
    df["burden_defined_level"] = dd["burden_evidence_level"].to_numpy()
    df["burden_value"] = pr["value"].to_numpy()
    has = df["burden_value"].notna()
    df["burden_measure_id"] = np.where(has, pr["measure_id"].to_numpy(), dd["primary_measure_id"].to_numpy())
    df["burden_measure_label"] = pr["measure_label"].to_numpy()
    df["burden_metric_type"] = pr["metric_type"].to_numpy()
    df["burden_value_unit"] = pr["value_unit"].to_numpy()
    df["burden_evidence_level"] = np.where(has, pr["burden_evidence_level"].to_numpy(), "D")
    df["burden_source_resolution"] = pr["source_geographic_resolution"].to_numpy()
    df["burden_inherited"] = pr["inherited"].fillna(False).astype(bool).to_numpy()
    df["burden_ci_low"] = pr["ci_low"].to_numpy()
    df["burden_ci_high"] = pr["ci_high"].to_numpy()
    df["burden_year"] = pr["year"].to_numpy()
    df["burden_period"] = pr["period"].to_numpy()
    df["burden_source_name"] = pr["source_name"].to_numpy()
    df["burden_row_id"] = pr["burden_row_id"].to_numpy()
    df["burden_uncertainty_multiplier"] = df["burden_evidence_level"].map(LEVEL_MULT)
    reason = pd.Series(pd.NA, index=df.index, dtype="object")
    is_d = df["burden_defined_level"] == "D"
    reason[is_d] = dd["rationale"].to_numpy()[is_d.to_numpy()]
    miss = ~has & ~is_d
    reason[miss & (df["state_fips"] == "09") & df["condition_id"].isin(["me_cfs", "fibromyalgia", "migraine"])
           & (df["geo_level"] == "county")] = ("CMS MMD publishes Connecticut only on legacy counties (09001-09015); "
                                               "the 2024 planning regions have no CMS value (not remapped).")
    reason[miss & reason.isna() & ~df["in_analysis_universe"]] = "Source does not cover this territory."
    # not published: the measured, source-specific reason (entity absent from the source vs cell omitted)
    rest = miss & reason.isna()
    reason[rest] = not_published_reasons(df.loc[rest], burden)
    df["burden_unknown_reason"] = reason
    # ---- counts (same resolution only)
    df["burden_count"] = np.nan
    df["burden_count_kind"] = None
    lcs = (df["condition_id"] == "long_covid") & (df["geo_level"] == "state") & has
    df.loc[lcs, "burden_count"] = df.loc[lcs, "burden_value"] / 100 * df.loc[lcs, "population_adults_18plus"]
    df.loc[lcs, "burden_count_kind"] = "survey_prevalence_x_acs_adults (HPS period x ACS 2020-2024; approximate)"
    cnt = burden[(burden["condition_id"] == "lyme_disease") & (burden["measure_role"] == "count_companion")]
    cmap = cnt.set_index(["geo_level", "geo_id"])["value"]
    ly = (df["condition_id"] == "lyme_disease")
    lk = pd.MultiIndex.from_frame(df.loc[ly, ["geo_level", "geo_id"]])
    df.loc[ly, "burden_count"] = cmap.reindex(lk).to_numpy()
    df.loc[ly & df["burden_count"].notna(), "burden_count_kind"] = f"surveillance_cases_{PRIMARY_LYME_YEAR}"
    df["burden_count_note"] = np.where(df["burden_count"].notna(), "same-resolution count",
                                       np.where(df["burden_defined_level"] == "D", "no burden measure",
                                                "not defensible at this resolution (inherited value, proxy, or "
                                                "Medicare FFS denominators published only as size classes)"))
    # ---- county C-proxy column (+ state aggregate)
    px = burden[burden["measure_role"].eq("county_proxy") | (burden["derivation"].eq("aggregated_from_county")
                                                             & burden["measure_id"].str.endswith("_crude_state_agg"))]
    px = px.drop_duplicates(["condition_id", "geo_level", "geo_id"]).set_index(["condition_id", "geo_level", "geo_id"])
    pp = px.reindex(key)
    df["burden_proxy_measure_id"] = pp["measure_id"].to_numpy()
    df["burden_proxy_value"] = pp["value"].to_numpy()
    df["burden_proxy_level"] = np.where(pp["value"].notna().to_numpy(), "C", None)
    df["burden_proxy_ci_low"] = pp["ci_low"].to_numpy()
    df["burden_proxy_ci_high"] = pp["ci_high"].to_numpy()
    df["burden_proxy_derivation"] = pp["derivation"].to_numpy()
    df["burden_proxy_count_modeled"] = (pp["value"] / 100 * pp["denominator"]).to_numpy()
    df["burden_proxy_uncertainty_multiplier"] = np.where(pp["value"].notna().to_numpy(), LEVEL_MULT["C"], np.nan)
    # ---- vulnerability + rurality + HPSA
    for c in ["svi_overall", "svi_socioeconomic", "svi_household", "svi_minority", "svi_housing_transport",
              "svi_ranking_universe", "svi_derivation", "rucc_2023", "is_metro", "hpsa_pc_any_designation",
              "hpsa_pc_whole_county", "hpsa_pc_coverage", "hpsa_pc_max_score", "hpsa_mh_any_designation",
              "hpsa_mua_coverage"]:
        df[c] = df["geo_id"].map(ctx[c]) if c in ctx.columns else np.nan
    # ---- providers / HRSA / trials / NIH
    counties = geos[(geos["geo_level"] == "county") & ~geos["ct_legacy"]][["geo_id", "lat", "lon"]].reset_index(drop=True)
    zpop = read_table("zcta_centroids", columns=["zcta", "lat", "lon"]).merge(
        read_table("geo_context_zcta__acs", columns=["zcta", "total_population"]), on="zcta", how="left")
    pcounts, pmeta = provider_counts(counties, zpop)
    tcounts = trial_counts(counties)
    ncounts = nih_counts()
    for level, frames in (("county", (pcounts.county, tcounts.county, ncounts.county)),
                          ("state", (pcounts.state, tcounts.state, ncounts.state))):
        m = df["geo_level"] == level
        k = pd.MultiIndex.from_frame(df.loc[m, ["condition_id", "geo_id"]])
        for f in frames:
            for c in f.columns:
                if c not in df.columns:
                    df[c] = np.nan
                df.loc[m, c] = f[c].reindex(k).to_numpy()
    count_cols = [c for c in df.columns if (c.startswith("trials_") or c.startswith("nih_")
                                            or c.startswith("relevant_providers") or c.startswith("relevant_specialists"))
                  and c.endswith("_n")]
    inu = df["population_total"].notna()
    for c in count_cols:   # a missing count in a geography with population is a true zero
        df.loc[inu & df[c].isna(), c] = 0
    for c in [c for c in count_cols if "within_" in c]:
        df.loc[df["geo_level"] == "state", c] = np.nan   # radius counts are county-only
    pop_within = pmeta["pop_within"].set_index("geo_id")[f"population_within_{PROVIDER_RADIUS_KM}km"]
    df[f"population_within_{PROVIDER_RADIUS_KM}km"] = np.where(df["geo_level"] == "county", df["geo_id"].map(pop_within), np.nan)
    per100k = lambda n, p: np.where(p > 0, n / p * 1e5, np.nan)  # noqa: E731
    df["relevant_providers_per_100k"] = per100k(df["relevant_providers_n"], df["population_total"])
    df["relevant_specialists_per_100k"] = per100k(df["relevant_specialists_n"], df["population_total"])
    df[f"relevant_providers_per_100k_within_{PROVIDER_RADIUS_KM}km"] = per100k(
        df[f"relevant_providers_within_{PROVIDER_RADIUS_KM}km_n"], df[f"population_within_{PROVIDER_RADIUS_KM}km"])
    core = _condition_core_groups()
    df["relevant_specialty_groups"] = df["condition_id"].map(lambda c: ";".join(core.get(c, [])))
    df["provider_count_variant"] = ("distinct individual NPIs (entity type 1), any taxonomy code in a core group of "
                                    "configs/relevance.yaml condition_specialties, primary practice location; "
                                    "self-reported taxonomy, not evidence of treating the condition")
    hc, hs = hrsa_counts()
    df["hrsa_sites_n"] = np.where(df["geo_level"] == "county", df["geo_id"].map(hc), df["geo_id"].map(hs))
    df.loc[inu & df["hrsa_sites_n"].isna(), "hrsa_sites_n"] = 0
    df["hrsa_sites_per_100k"] = per100k(df["hrsa_sites_n"], df["population_total"])
    df["trial_precision_view"] = "literal (trial_conditions.condition_literal_match); any registration date"
    df["nih_precision_view"] = "title_abstract match, likely false positives excluded; FY2015-FY2026"
    # ---- diagnostic desert
    df = add_desert(df)
    # ---- ids + provenance
    df["object_id"] = "geo:" + df["geo_id"]
    df["feature_id"] = df["condition_id"] + "|" + df["geo_id"]
    df["data_layer"] = "derived"
    df["source_name"] = ("measure_it.geography.features: geo_condition_burden + geo_context + NPPES + HRSA + "
                         "ClinicalTrials.gov + NIH RePORTER (ecological composite)")
    df["source_record_id"] = df["feature_id"]
    df["source_version"] = "built " + utc_now_iso()
    df["retrieved_at"] = burden["retrieved_at"].max()
    df["source_geographic_resolution"] = df["geo_level"]
    df["evidence_type"] = "derived_score"
    df["evidence_level"] = df["burden_evidence_level"]
    df["provenance_notes"] = ("Ecological features: burden, vulnerability, providers, trials and NIH describe places, "
                              "not people, and are never person-linked. Inherited burden keeps source resolution "
                              "'state'. Research activity is not burden.")
    df = df.drop(columns=["name"])
    return df, pmeta


def _desert_components(g: pd.DataFrame, burden_col: str, prov_col: str, trial_col: str) -> pd.DataFrame:
    u = g[g["in_analysis_universe"]]
    out = pd.DataFrame(index=g.index)
    out["b"] = percentile_rank(u[burden_col])
    out["v"] = percentile_rank(u["svi_overall"])
    out["p"] = 1 - percentile_rank(u[prov_col])
    out["t"] = 1 - percentile_rank(u[trial_col])
    return out


def add_desert(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in ["desert_pct_burden", "desert_pct_vulnerability", "desert_pct_low_providers", "desert_pct_low_trials",
              "diagnostic_desert", "diagnostic_desert_access_only"]:
        df[c] = np.nan
    sens = {
        "diagnostic_desert_sens_specialists": ("burden_value", "relevant_specialists_per_100k", None),
        "diagnostic_desert_sens_providers_50km": ("burden_value", f"relevant_providers_per_100k_within_{PROVIDER_RADIUS_KM}km", None),
        "diagnostic_desert_sens_trials_100km": ("burden_value", "relevant_providers_per_100k", "trials_in_geo_or_within_100km_n"),
        "diagnostic_desert_sens_recruiting_100km": ("burden_value", "relevant_providers_per_100k",
                                                    "trials_in_geo_or_within_100km_recruiting_n"),
        "diagnostic_desert_sens_trials_50km_internal_point_only": ("burden_value", "relevant_providers_per_100k",
                                                                   "trials_within_50km_n"),
        "diagnostic_desert_sens_proxy_burden": ("burden_proxy_value", "relevant_providers_per_100k", None),
    }
    for c in sens:
        df[c] = np.nan
    for (cond, level), g in df.groupby(["condition_id", "geo_level"]):
        trial_col = "trials_in_geo_or_within_50km_n" if level == "county" else "trials_in_geo_n"
        comp = _desert_components(g, "burden_value", "relevant_providers_per_100k", trial_col)
        df.loc[g.index, "desert_pct_burden"] = comp["b"]
        df.loc[g.index, "desert_pct_vulnerability"] = comp["v"]
        df.loc[g.index, "desert_pct_low_providers"] = comp["p"]
        df.loc[g.index, "desert_pct_low_trials"] = comp["t"]
        df.loc[g.index, "diagnostic_desert"] = comp[["b", "v", "p", "t"]].mean(axis=1, skipna=False)
        df.loc[g.index, "diagnostic_desert_access_only"] = comp[["v", "p", "t"]].mean(axis=1, skipna=False)
        for c, (bcol, pcol, tcol) in sens.items():
            if level == "state" and ("50km" in c or "100km" in c or "internal_point" in c):
                continue
            if c.endswith("proxy_burden") and not (cond == "long_covid" and level == "county"):
                continue
            s = _desert_components(g, bcol, pcol, tcol or trial_col)
            df.loc[g.index, c] = s[["b", "v", "p", "t"]].mean(axis=1, skipna=False)
    reason = pd.Series(pd.NA, index=df.index, dtype="object")
    reason[~df["in_analysis_universe"]] = "outside analysis universe (50 states + DC): SVI/PLACES/HPS not poolable"
    reason[reason.isna() & (df["burden_defined_level"] == "D")] = "level D: no burden measure (see diagnostic_desert_access_only)"
    reason[reason.isna() & df["burden_value"].isna()] = "burden value missing for this geography"
    reason[reason.isna() & df["diagnostic_desert"].isna()] = "a non-burden component is missing"
    df["desert_unknown_reason"] = reason.where(df["diagnostic_desert"].isna())
    base_def = ("mean(pr(burden), pr(SVI overall), 1-pr(relevant providers per 100k), "
                "1-pr(relevant trials with a U.S. site in the county or within 50 km of its internal "
                "point [county] / in the state [state])); percentile ranks within 50 states + DC per "
                "condition and level; equal weights")
    inh = df["burden_inherited"].fillna(False).astype(bool)
    mod = (df["burden_source_name"].astype(str).str.contains("small-area", na=False)
           & ~inh) if "burden_source_name" in df.columns else pd.Series(False, index=df.index)
    df["desert_definition"] = np.where(
        inh, base_def + ". BURDEN COMPONENT IS THE INHERITED STATE VALUE: every county in a state gets the same "
                        "burden percentile, so this index ranks counties partly by their state's estimate; it is "
                        "not evidence of county-level burden",
        np.where(mod, base_def + ". BURDEN COMPONENT IS A MODELLED SMALL-AREA ESTIMATE (BRFSS 2023 MRP): its "
                                 "within-state differences come from county composition and covariates only, not "
                                 "from observed county prevalence", base_def))
    # reviewer addition: resolution of the burden component actually used in the index
    df["desert_burden_resolution"] = np.where(
        df["diagnostic_desert"].isna(), None,
        np.where(inh, "state (inherited; no within-state variation)",
                 np.where(mod, "county (modelled small-area estimate)",
                          df["burden_source_resolution"].astype(object))))
    return df


# ============================================================================ analysis: Test 4
def _series_frames(features: pd.DataFrame, burden: pd.DataFrame) -> list[dict]:
    """Pre-specified Test 4 series: dicts with id, level, frame (index geo_id; prev, count, pop, svi, desert)."""
    f = features[features["in_analysis_universe"]]
    out = []

    def mk(sid, label, cond, level, prev_col="burden_value", count_col="burden_count", desert=False):
        g = f[(f["condition_id"] == cond) & (f["geo_level"] == level)].set_index("geo_id")
        fr = pd.DataFrame({"prev": g[prev_col], "count": g[count_col] if count_col else np.nan,
                           "pop": g["population_total"], "svi": g["svi_overall"],
                           "desert": g["diagnostic_desert"] if desert else np.nan})
        out.append({"series_id": sid, "label": label, "condition_id": cond, "geo_level": level, "frame": fr})

    mk("lyme_county", "Lyme incidence 2023 (A) / reported cases", "lyme_disease", "county")
    mk("lyme_state", "Lyme incidence 2023, state aggregate (A) / cases", "lyme_disease", "state")
    mk("long_covid_state", "HPS current long COVID % (A) / adults x %", "long_covid", "state")
    mk("long_covid_county_proxy", "PLACES PHLTH crude (C proxy) / modeled count", "long_covid", "county",
       prev_col="burden_proxy_value", count_col="burden_proxy_count_modeled")
    mk("cms51_county", "CMS 51 fibromyalgia/pain/fatigue prevalence (B; C for ME/CFS)", "fibromyalgia", "county",
       count_col=None)
    mk("cms59_county", "CMS 59 migraine/headache prevalence (B)", "migraine", "county", count_col=None)
    for cond in ("long_covid", "me_cfs", "fibromyalgia", "migraine", "lyme_disease", "ptlds"):
        for level in ("county", "state"):
            g = f[(f["condition_id"] == cond) & (f["geo_level"] == level)].set_index("geo_id")
            if g["diagnostic_desert"].notna().sum() < 5:
                continue
            fr = pd.DataFrame({"prev": g["burden_value"], "count": np.nan, "pop": g["population_total"],
                               "svi": g["svi_overall"], "desert": g["diagnostic_desert"]})
            out.append({"series_id": f"desert_{cond}_{level}", "label": f"diagnostic desert, {cond}",
                        "condition_id": cond, "geo_level": level, "frame": fr, "desert_only": True})
    return out


def test4_series_legend(corr: pd.DataFrame) -> str:
    """One line per Test 4 series: its measure and burden evidence level, so a proxy series (e.g.
    long_covid_county_proxy = PLACES PHLTH, level C) is never read as the named condition's prevalence."""
    if corr is None or corr.empty or not {"series_id", "label"} <= set(corr.columns):
        return ""
    seen = corr[["series_id", "label"]].drop_duplicates("series_id")
    items = [f"`{s}` = {lab}" for s, lab in seen.itertuples(index=False) if not str(s).startswith("desert_")]
    return ("\nSeries (measure / count, burden evidence level in brackets; `variable` = prevalence is that measure's "
            "rate, a proxy where the level is C): " + "; ".join(items) + ". `desert_*` rows are the diagnostic desert "
            "index of the named condition.\n")


def run_test4(features: pd.DataFrame, burden: pd.DataFrame) -> dict[str, pd.DataFrame]:
    min_pop = load_config("scoring")["ranking"]["min_population"]
    top_n = load_config("scoring")["ranking"]["top_n"]
    corr_rows, ov_rows, part_rows, dec_rows, verdict_rows, top_rows = [], [], [], [], [], []
    small_rows = []
    for s in _series_frames(features, burden):
        fr, sid, lvl = s["frame"], s["series_id"], s["geo_level"]
        desert_only = s.get("desert_only", False)
        vars_ = {"desert": fr["desert"]} if desert_only else {"prevalence": fr["prev"], "count": fr["count"],
                                                              "desert": fr["desert"]}
        for vname, v in vars_.items():
            ok = v.notna() & fr["pop"].notna()
            if ok.sum() < 5:
                continue
            rho = spearman(v[ok], fr["pop"][ok])
            lo, hi = bootstrap_spearman_ci(v[ok], fr["pop"][ok])
            p = float(stats.spearmanr(v[ok], fr["pop"][ok]).pvalue)
            if lvl == "county":   # reviewer addition: resample whole states (spatial dependence)
                clo, chi = cluster_bootstrap_spearman_ci(v[ok], fr["pop"][ok], v[ok].index.str[:2].to_numpy())
            else:
                clo, chi = np.nan, np.nan
            corr_rows.append({"series_id": sid, "label": s["label"], "geo_level": lvl, "variable": vname,
                              "n": int(ok.sum()), "rho_with_population": rho, "ci95_low": lo, "ci95_high": hi,
                              "ci95_low_state_cluster": clo, "ci95_high_state_cluster": chi, "p_value": p})
            # reviewer addition: the pre-specified top-N test only looks for LARGE-population dominance. Check
            # the other tail: are the top-N geographies concentrated among the smallest populations?
            N_ = top_n if lvl == "county" else TOP_N_STATE
            fr3 = fr[ok & (fr["pop"] >= min_pop)] if lvl == "county" else fr[ok]
            if len(fr3) >= 2 * N_:
                vv, pp = v.loc[fr3.index], fr3["pop"]
                top = top_n_ids(vv, N_)
                smallest = top_n_ids(-pp, N_)
                q = pp.quantile(0.25)
                bottom = set(pp.index[pp <= q])
                M, K = len(fr3), len(bottom)
                k_small, k_q = len(top & smallest), len(top & bottom)
                small_rows.append({
                    "series_id": sid, "variable": vname, "geo_level": lvl, "M": M, "N": N_,
                    "population_floor": min_pop if lvl == "county" else 0,
                    "topN_overlap_with_smallest_N": k_small, "expected_overlap_smallest": N_ * N_ / M,
                    "p_hypergeom_smallest": float(stats.hypergeom.sf(k_small - 1, M, N_, N_)),
                    "bottom_quartile_population_max": float(q), "topN_in_bottom_pop_quartile": k_q,
                    "expected_in_bottom_quartile": N_ * K / M,
                    "p_hypergeom_bottom_quartile": float(stats.hypergeom.sf(k_q - 1, M, K, N_)),
                    "topN_in_bottom_pop_quartile_tie_averaged": float(
                        topn_inclusion_prob(vv, N_).reindex(sorted(bottom)).fillna(0).sum()),
                    "ties_at_topN_boundary": int(((topn_inclusion_prob(vv, N_) > 0)
                                                  & (topn_inclusion_prob(vv, N_) < 1)).sum()),
                    "topN_median_population": float(pp.loc[list(top)].median()),
                    "all_median_population": float(pp.median())})
        # top-N overlaps
        for floor, N in ((True, top_n if lvl == "county" else TOP_N_STATE), (True, 100 if lvl == "county" else None),
                         (False, top_n if lvl == "county" else None)):
            if N is None:
                continue
            fr2 = fr[fr["pop"] >= min_pop] if (floor and lvl == "county") else fr
            pairs = ([("desert", "population")] if desert_only else
                     [("prevalence", "population"), ("count", "population"), ("prevalence", "count"),
                      ("desert", "population")])
            cols = {"prevalence": "prev", "count": "count", "population": "pop", "desert": "desert"}
            for a, b in pairs:
                if fr2[cols[a]].notna().sum() < N or fr2[cols[b]].notna().sum() < N:
                    continue
                o = topn_overlap(fr2[cols[a]], fr2[cols[b]], N)
                ov_rows.append({"series_id": sid, "geo_level": lvl, "ranking_a": a, "ranking_b": b,
                                "population_floor": min_pop if (floor and lvl == "county") else 0, **o})
        # partial correlations
        lp = np.log(fr["pop"].where(fr["pop"] > 0))
        if not desert_only:
            if fr["count"].notna().sum() >= 5:
                part_rows.append({"series_id": sid, "x": "count", "y": "prevalence", "control": "log_population",
                                  "n": int((fr["count"].notna() & fr["prev"].notna() & lp.notna()).sum()),
                                  "raw_spearman": spearman(fr["count"], fr["prev"]),
                                  "partial_spearman": partial_spearman(fr["count"], fr["prev"], lp)})
            if fr["svi"].notna().sum() >= 5:
                part_rows.append({"series_id": sid, "x": "prevalence", "y": "svi_overall", "control": "log_population",
                                  "n": int((fr["svi"].notna() & fr["prev"].notna() & lp.notna()).sum()),
                                  "raw_spearman": spearman(fr["prev"], fr["svi"]),
                                  "partial_spearman": partial_spearman(fr["prev"], fr["svi"], lp)})
        if fr["desert"].notna().sum() >= 5:
            part_rows.append({"series_id": sid, "x": "desert", "y": "burden", "control": "log_population",
                              "n": int((fr["desert"].notna() & fr["prev"].notna() & lp.notna()).sum()),
                              "raw_spearman": spearman(fr["desert"], fr["prev"]),
                              "partial_spearman": partial_spearman(fr["desert"], fr["prev"], lp)})
        # log-count decomposition
        if not desert_only and fr["count"].notna().sum() >= 5:
            d = fr[(fr["count"] > 0) & (fr["prev"] > 0) & (fr["pop"] > 0)]
            if len(d) >= 5:
                lc_, lpop, lprev = np.log(d["count"]), np.log(d["pop"]), np.log(d["prev"])
                dec_rows.append({"series_id": sid, "n_count_positive": len(d),
                                 "n_count_zero": int((fr["count"] == 0).sum()),
                                 "r2_logcount_on_logpop": float(np.corrcoef(lc_, lpop)[0, 1] ** 2),
                                 "r2_logcount_on_logprev": float(np.corrcoef(lc_, lprev)[0, 1] ** 2)})
        # top lists (transparency)
        N = top_n if lvl == "county" else TOP_N_STATE
        fr2 = fr[fr["pop"] >= min_pop] if lvl == "county" else fr
        for rk in (["desert"] if desert_only else ["prevalence", "count", "population", "desert"]):
            col = {"prevalence": "prev", "count": "count", "population": "pop", "desert": "desert"}[rk]
            if fr2[col].notna().sum() < N:
                continue
            ids = fr2[col].dropna().sort_values(ascending=False).index[:N]
            for i, gid in enumerate(ids, 1):
                top_rows.append({"series_id": sid, "ranking": rk, "rank": i, "geo_id": gid, "value": fr2.loc[gid, col],
                                 "population": fr2.loc[gid, "pop"]})
    corr = pd.DataFrame(corr_rows)
    ov = pd.DataFrame(ov_rows)
    part = pd.DataFrame(part_rows)
    dec = pd.DataFrame(dec_rows)
    # verdicts per pre-specified rule
    for _, r in corr.iterrows():
        lvl = r["geo_level"]
        N = top_n if lvl == "county" else TOP_N_STATE
        o = ov[(ov["series_id"] == r["series_id"]) & (ov["ranking_a"] == r["variable"]) & (ov["ranking_b"] == "population")
               & (ov["N"] == N) & (ov["population_floor"] == (min_pop if lvl == "county" else 0))] if len(ov) else ov
        o = o.iloc[0] if len(o) else None
        above = (o is not None) and (o["p_hypergeom"] < 0.05)
        rho = r["rho_with_population"]
        if rho > 0.7 and above:
            v = "population-driven"
        elif abs(rho) < 0.3 and (o is None or not above):
            v = "burden-driven (not population-driven)"
        else:
            v = "mixed"
        verdict_rows.append({"series_id": r["series_id"], "variable": r["variable"], "geo_level": lvl,
                             "rho_with_population": rho, "topN": None if o is None else int(o["N"]),
                             "topN_overlap_with_population": None if o is None else o["overlap"],
                             "expected_overlap": None if o is None else o["expected_overlap"],
                             "p_hypergeom": None if o is None else o["p_hypergeom"], "verdict": v})
    ver = pd.DataFrame(verdict_rows)
    # reviewer additions: multiplicity (Holm across the verdict family) and the small-population tail
    ver["p_hypergeom_holm"] = holm(ver["p_hypergeom"].astype(float))
    small = pd.DataFrame(small_rows)
    if len(small):
        small["p_hypergeom_bottom_quartile_holm"] = holm(small["p_hypergeom_bottom_quartile"])
        ver = ver.merge(small[["series_id", "variable", "topN_in_bottom_pop_quartile", "expected_in_bottom_quartile",
                               "p_hypergeom_bottom_quartile", "p_hypergeom_bottom_quartile_holm"]],
                        on=["series_id", "variable"], how="left")
        sig = ver["p_hypergeom_bottom_quartile_holm"] < 0.05
        ver["reviewer_reading"] = np.select(
            [ver["verdict"].eq("population-driven"),
             sig & (ver["rho_with_population"] < 0),
             sig,
             ver["rho_with_population"].abs() >= 0.3],
            ["population-driven (pre-specified rule)",
             "not driven by LARGE population, but the top-N is concentrated in small-population geographies "
             "(inverse population association)",
             "top-N concentrated in small-population geographies although rho with population is not negative",
             "moderate association with population (|rho| >= 0.3) without top-N concentration in either tail"],
            default="no top-N concentration in either population tail after Holm correction (|rho| < 0.3)")
    geos = _geo_names().set_index("geo_id")
    top = pd.DataFrame(top_rows)
    if len(top):
        top["geo_name"] = top["geo_id"].map(geos["name"] + ", " + geos["state_abbr"].fillna(""))
    return {"test4_population_correlations": corr, "test4_topn_overlap": ov, "test4_partial_correlations": part,
            "test4_log_count_decomposition": dec, "test4_verdicts": ver,
            "test4_top_lists": top, "test4_small_population_check": small}


# ============================================================================ analysis: Test 6
def county_weights(county_ids: pd.Index) -> dict[str, sparse.csr_matrix]:
    import geopandas as gpd
    from ..config import PROCESSED
    poly = gpd.read_parquet(PROCESSED / "county_boundaries_2024.geoparquet").set_index("geo_id").reindex(county_ids)
    poly = poly[poly.geometry.notna()].to_crs(5070)
    idx = {g: i for i, g in enumerate(county_ids)}
    buf = poly.copy()
    buf["geometry"] = buf.geometry.buffer(CONTIGUITY_BUFFER_M)
    j = gpd.sjoin(buf[["geometry"]].reset_index(), buf[["geometry"]].reset_index(), predicate="intersects")
    j = j[j["geo_id_left"] != j["geo_id_right"]]
    r = j["geo_id_left"].map(idx).to_numpy()
    c = j["geo_id_right"].map(idx).to_numpy()
    n = len(county_ids)
    w_cont = sparse.csr_matrix((np.ones(len(r)), (r, c)), shape=(n, n))
    w_cont = ((w_cont + w_cont.T) > 0).astype(float)
    g = read_table("geographies").set_index("geo_id").reindex(county_ids)
    band = radius_matrix(g["lat"].to_numpy(), g["lon"].to_numpy(), g["lat"].to_numpy(), g["lon"].to_numpy(),
                         MORAN_BAND_KM).astype(float)
    band = band - sparse.diags(band.diagonal())
    band.eliminate_zeros()
    return {"contiguity": w_cont.tocsr(), f"distance_band_{MORAN_BAND_KM}km": band.tocsr().astype(float)}


def _moran_perm(x: np.ndarray, w: sparse.csr_matrix, groups: np.ndarray, rng, n_perm: int, scheme: str) -> np.ndarray:
    out = np.empty(n_perm)
    for b in range(n_perm):
        xp = rng.permutation(x) if scheme == "national" else permute_within_groups(x, groups, rng)
        out[b] = morans_i(xp, w)
    return out


def _corr_perm(x: np.ndarray, y: np.ndarray, groups: np.ndarray | None, rng, n_perm: int, scheme: str) -> np.ndarray:
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    out = np.empty(n_perm)
    for b in range(n_perm):
        xp = rng.permutation(rx) if scheme == "national" else permute_within_groups(rx, groups, rng)
        out[b] = np.corrcoef(xp, ry)[0, 1]
    return out


def run_test6(features: pd.DataFrame, burden: pd.DataFrame, context: pd.DataFrame) -> dict[str, pd.DataFrame]:
    rng = np.random.default_rng(config.SEED)
    cty = context[(context["geo_level"] == "county") & context["state_fips"].isin(UNIVERSE_STATES)].set_index("geo_id")
    fc = features[(features["geo_level"] == "county") & features["in_analysis_universe"]]
    bc = burden[burden["geo_level"] == "county"]

    def bvals(measure_id, year=None, cond=None, role=None):
        s = bc[bc["measure_id"] == measure_id]
        if year is not None:
            s = s[s["year"] == year]
        if cond is not None:
            s = s[s["condition_id"] == cond]
        if role is not None:
            s = s[s["measure_role"] == role]
        return s.drop_duplicates("geo_id").set_index("geo_id")["value"]

    cms51 = bvals("cms_mmd_51_prev_actual_all", PRIMARY_CMS_YEAR, "fibromyalgia", "primary")
    cms59 = bvals("cms_mmd_59_prev_actual_all", PRIMARY_CMS_YEAR, "migraine", "primary")
    lyme23 = bvals("lyme_incidence_per_100k", PRIMARY_LYME_YEAR, "lyme_disease")
    lyme19 = bvals("lyme_incidence_per_100k", ALT_LYME_PRE2022_YEAR, "lyme_disease")
    phlth = cty["places_PHLTH_crude"]
    arth = cty["places_ARTHRITIS_crude"]
    rels = [("R1", "CMS 51 prevalence (claims) vs PLACES ARTHRITIS crude (BRFSS model)", cms51, arth),
            ("R2", "CMS 51 prevalence vs PLACES PHLTH crude", cms51, phlth),
            ("R3", "CMS 59 migraine prevalence vs PLACES PHLTH crude", cms59, phlth),
            ("R4", "Positive control: Lyme incidence 2023 vs 2019 (same source)", lyme23, lyme19)]
    rel_rows, null_rows = [], []
    for rid, desc, x, y in rels:
        d = pd.DataFrame({"x": x, "y": y}).reindex(cty.index).dropna()
        groups = d.index.str[:2].to_numpy()
        obs = spearman(d["x"], d["y"])
        crh = crh_effective_n(d["x"], d["y"], cty.loc[d.index, "lat"], cty.loc[d.index, "lon"])
        for scheme in ("national", "within_state"):
            null = _corr_perm(d["x"].to_numpy(), d["y"].to_numpy(), groups, rng, N_PERM, scheme)
            rel_rows.append({"relationship_id": rid, "description": desc, "geo_level": "county", "n": len(d),
                             "observed_spearman": obs, "permutation_scheme": scheme, "n_perm": N_PERM,
                             **permutation_summary(obs, null, 0.0),
                             **(crh if scheme == "national" else {}),
                             "min_detectable_rho_80pct": min_detectable_rho(
                                 crh["n_eff_spatial"] if scheme == "national" else len(d))})
            null_rows.append(pd.DataFrame({"test": rid, "scheme": scheme, "null_value": null}))
    # state relationships
    sb = burden[(burden["geo_level"] == "state") & burden["geo_id"].isin(UNIVERSE_STATES)]
    # the HPS state estimate itself (latest non-suppressed period), whether or not it is the primary measure
    hsel = sb[(sb["condition_id"] == "long_covid") & (sb["source_table"] == "geo_condition_burden__cdc_long_covid")
              & (sb["measure_id"] == "lc_current_pct_all_adults") & (sb["derivation"] == "as_published")
              & sb["value"].notna()]
    hps = hsel.sort_values("period_end").groupby("geo_id").tail(1).set_index("geo_id")["value"]
    phl_st = sb[sb["measure_id"] == "places_PHLTH_crude_state_agg"].drop_duplicates("geo_id").set_index("geo_id")["value"]
    cms51_st = sb[(sb["measure_id"] == "cms_mmd_51_prev_actual_all") & (sb["year"] == PRIMARY_CMS_YEAR)
                  & (sb["condition_id"] == "fibromyalgia")].set_index("geo_id")["value"]
    sctx = context[context["geo_level"] == "state"].set_index("geo_id")
    for rid, desc, x, y in (("R5", "State: HPS current long COVID % vs PLACES PHLTH state aggregate", hps, phl_st),
                            ("R6", "State: HPS current long COVID % vs CMS 51 state prevalence", hps, cms51_st)):
        d = pd.DataFrame({"x": x, "y": y}).dropna()
        obs = spearman(d["x"], d["y"])
        null = _corr_perm(d["x"].to_numpy(), d["y"].to_numpy(), None, rng, N_PERM, "national")
        crh = crh_effective_n(d["x"], d["y"], sctx.reindex(d.index)["lat"], sctx.reindex(d.index)["lon"],
                              band_km=250.0, max_km=3000.0)
        rel_rows.append({"relationship_id": rid, "description": desc, "geo_level": "state", "n": len(d),
                         "observed_spearman": obs, "permutation_scheme": "national", "n_perm": N_PERM,
                         "spearman_p_value": float(stats.spearmanr(d["x"], d["y"]).pvalue),
                         **permutation_summary(obs, null, 0.0), **crh,
                         "min_detectable_rho_80pct": min_detectable_rho(crh["n_eff_spatial"])})
        null_rows.append(pd.DataFrame({"test": rid, "scheme": "national", "null_value": null}))
    # Moran's I
    ids = cty.index
    W = county_weights(ids)
    # the inherited HPS state value (an alternate row since 2026-10-07) and the primary county burden (the BRFSS
    # small-area estimate when long_covid_burden = brfss_sae) are tested separately
    lci = burden[(burden["condition_id"] == "long_covid") & (burden["geo_level"] == "county") & burden["inherited"]
                 & burden["source_measure_id"].eq("lc_current_pct_all_adults")]
    lc_inh = lci.drop_duplicates("geo_id").set_index("geo_id")["value"]
    lc_prim = fc[fc["condition_id"] == "long_covid"].set_index("geo_id")["burden_value"]
    prov = fc[fc["condition_id"] == "long_covid"].set_index("geo_id")["relevant_providers_per_100k"]
    des = fc[fc["condition_id"] == "me_cfs"].set_index("geo_id")["diagnostic_desert"]
    variables = {
        "log1p_lyme_incidence_2023": np.log1p(lyme23), "cms51_prevalence_2022": cms51, "cms59_prevalence_2022": cms59,
        "places_PHLTH_crude": phlth, "places_ARTHRITIS_crude": arth, "svi_overall": cty["svi_overall"],
        "log1p_relevant_provider_density_long_covid": np.log1p(prov),
        "diagnostic_desert_me_cfs": des, "long_covid_inherited_state_value": lc_inh,
    }
    if not lc_prim.reindex(lc_inh.index).equals(lc_inh):
        variables["long_covid_primary_county_burden"] = lc_prim
    mor_rows = []
    for vname, s in variables.items():
        x = s.reindex(ids)
        for wname, wfull in W.items():
            ok = x.notna().to_numpy()
            w = wfull[ok][:, ok]
            has_nb = np.asarray(w.sum(axis=1)).ravel() > 0
            xv = x.to_numpy()[ok][has_nb]
            w = row_standardize(w[has_nb][:, has_nb])
            keep = np.asarray(w.sum(axis=1)).ravel() > 0
            xv, w = xv[keep], row_standardize(w[keep][:, keep])
            groups = ids.to_numpy()[ok][has_nb][keep]
            groups = np.array([g[:2] for g in groups])
            obs = morans_i(xv, w)
            e_i = -1.0 / (len(xv) - 1)
            for scheme in ("national", "within_state"):
                null = _moran_perm(xv, w, groups, rng, N_PERM_MORAN, scheme)
                mor_rows.append({"variable": vname, "weights": wname, "n": len(xv), "observed_I": obs,
                                 "expected_I": e_i, "permutation_scheme": scheme, "n_perm": N_PERM_MORAN,
                                 **permutation_summary(obs, null, e_i)})
                null_rows.append(pd.DataFrame({"test": f"moran:{vname}:{wname}", "scheme": scheme, "null_value": null}))
    rel = pd.DataFrame(rel_rows)
    rel["collapsed_under_permutation"] = rel["observed_outside_null95"] & (rel["null_mean"].abs() < 0.05)
    mor = pd.DataFrame(mor_rows)
    mor["collapsed_under_permutation"] = mor["observed_outside_null95"] & ((mor["null_mean"] - mor["expected_I"]).abs() < 0.02)
    # reviewer additions: Holm within each family (relationships x scheme; Moran x scheme). The national family
    # also gets a Holm-adjusted spatial (CRH) p. Within-state rows use the null-centred p.
    for t in (rel, mor):
        t["p_holm"] = np.nan
        for scheme, idx in t.groupby("permutation_scheme").groups.items():
            col = "p_perm" if scheme == "national" else "p_perm_vs_null_centre"
            t.loc[idx, "p_holm"] = holm(t.loc[idx, col])
    if "p_spatial_crh" in rel.columns:
        nat = rel["permutation_scheme"] == "national"
        rel["p_spatial_crh_holm"] = np.nan
        rel.loc[nat, "p_spatial_crh_holm"] = holm(rel.loc[nat, "p_spatial_crh"])
        rel["spatially_adjusted_reading"] = np.where(
            ~nat, None,
            np.where(~rel["observed_outside_null95"], "no detectable structure (naive or adjusted)",
                     np.where(rel["p_spatial_crh_holm"] < 0.05,
                              "survives spatial (CRH effective-n) adjustment and Holm correction",
                              "significant only under the naive national shuffle; does NOT survive spatial (CRH) "
                              "adjustment + Holm")))
    rel["no_structure_value"] = 0.0
    mor["no_structure_value"] = mor["expected_I"]
    for t, obs_col in ((rel, "observed_spearman"), (mor, "observed_I")):
        nat = t["permutation_scheme"] == "national"
        out_ = t["observed_outside_null95"]
        t["within_state_component"] = np.where(nat, np.nan, t[obs_col] - t["null_mean"])
        t["interpretation"] = np.select(
            [nat & out_, nat & ~out_, ~nat & ~out_, ~nat & out_ & (t[obs_col] > t["null_p97_5"])],
            ["structure under the naive national shuffle: collapses to the no-structure value (see "
             "spatially_adjusted_reading)",
             "no detectable structure: observed value lies within the national permutation band",
             "fully carried by between-state differences: the within-state shuffle reproduces the observed value",
             "survives beyond between-state differences: observed exceeds the within-state-shuffle band "
             "(between-state share = retained_fraction)"],
            default="between-state component (null mean) exceeds the observed value: the within-state association "
                    "is weaker or opposite in sign")
    return {"test6_geography_relationships": rel, "test6_geography_morans_i": mor,
            "test6_geography_null_draws": pd.concat(null_rows, ignore_index=True)}


def crh_calibration(context: pd.DataFrame, n_points: int = 1500, n_sim: int = 200,
                    ranges_km: tuple = (150, 400)) -> pd.DataFrame:
    """False-positive rate of the naive (exchangeable) test vs the CRH-adjusted test for two INDEPENDENT spatially
    autocorrelated Gaussian fields (exponential covariance) on real county internal points (contiguous US sample).
    Reviewer addition: shows why the national label shuffle cannot, on its own, establish a cross-source relationship.
    The naive test is the Spearman t-test, which for n in the thousands matches the national permutation test.
    """
    rng = np.random.default_rng(config.SEED)
    c = context[(context["geo_level"] == "county") & context["state_fips"].isin(UNIVERSE_STATES)
                & ~context["state_fips"].isin(["02", "15"])].dropna(subset=["lat", "lon"])
    c = c.iloc[np.sort(rng.choice(len(c), size=min(n_points, len(c)), replace=False))]
    lat, lon = c["lat"].to_numpy(), c["lon"].to_numpy()
    d = haversine_km(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
    rows = []
    for r_km in ranges_km:
        chol = np.linalg.cholesky(np.exp(-d / r_km) + 1e-6 * np.eye(len(lat)))
        naive, adj, neff = 0, 0, []
        for _ in range(n_sim):
            x, y = chol @ rng.normal(size=len(lat)), chol @ rng.normal(size=len(lat))
            naive += stats.spearmanr(x, y).pvalue < 0.05
            o = crh_effective_n(x, y, lat, lon)
            adj += o["p_spatial_crh"] < 0.05
            neff.append(o["n_eff_spatial"])
        rows.append({"field_range_km": r_km, "n_points": len(lat), "n_sim": n_sim, "nominal_alpha": 0.05,
                     "false_positive_rate_naive": naive / n_sim, "false_positive_rate_crh": adj / n_sim,
                     "median_n_eff": float(np.median(neff))})
    return pd.DataFrame(rows)


# ============================================================================ coverage, sensitivity
def coverage_matrix(features: pd.DataFrame, burden: pd.DataFrame) -> pd.DataFrame:
    defs = definitions_frame()
    rows = []
    for _, d in defs.iterrows():
        cond, lvl = d["condition_id"], d["geo_level"]
        if lvl == "national":
            b = burden[(burden["condition_id"] == cond) & (burden["geo_level"] == "national")
                       & (burden["measure_role"] == "primary")]
            n_val, n_tot = int(b["value"].notna().sum()), 1
        else:
            f = features[(features["condition_id"] == cond) & (features["geo_level"] == lvl) & features["in_analysis_universe"]]
            n_val, n_tot = int(f["burden_value"].notna().sum()), len(f)
        best = d["burden_evidence_level"]
        label = best + (" (inherited from state)" if d["inherited"] else "")
        proxy = d["county_proxy_places_measure"]
        has_proxy = isinstance(proxy, str) and bool(proxy)
        rows.append({"condition_id": cond, "geo_level": lvl, "best_level": best, "best_level_label": label,
                     "primary_measure_id": d["primary_measure_id"] if isinstance(d["primary_measure_id"], str) else UNKNOWN,
                     "source_resolution": (d["source_geographic_resolution"]
                                           if isinstance(d["source_geographic_resolution"], str) else UNKNOWN),
                     "inherited": d["inherited"],
                     "county_proxy_column": f"places_{proxy}_crude (C)" if has_proxy else "none",
                     "n_geos_with_value": n_val, "n_geos_in_universe": n_tot,
                     "coverage_pct": 100 * n_val / n_tot if n_tot else np.nan})
    return pd.DataFrame(rows)


def desert_sensitivity(features: pd.DataFrame) -> pd.DataFrame:
    rows = []
    top_n = load_config("scoring")["ranking"]["top_n"]
    min_pop = load_config("scoring")["ranking"]["min_population"]
    f = features[features["in_analysis_universe"]]
    variants = [c for c in f.columns if c.startswith("diagnostic_desert_sens_")] + ["diagnostic_desert_access_only"]
    for (cond, lvl), g in f.groupby(["condition_id", "geo_level"]):
        if g["diagnostic_desert"].notna().sum() < 5:
            continue
        g = g.set_index("geo_id")
        gf = g[g["population_total"] >= min_pop] if lvl == "county" else g
        N = top_n if lvl == "county" else TOP_N_STATE
        for v in variants:
            if g[v].notna().sum() < 5:
                continue
            o = topn_overlap(gf["diagnostic_desert"], gf[v], N)
            rows.append({"condition_id": cond, "geo_level": lvl, "variant": v,
                         "n": int((g[v].notna() & g["diagnostic_desert"].notna()).sum()),
                         "spearman_with_primary": spearman(g["diagnostic_desert"], g[v]),
                         "topN": o["N"], "topN_overlap": o["overlap"], "topN_jaccard": o["jaccard"]})
    return pd.DataFrame(rows)


def primary_vs_alternates(burden: pd.DataFrame) -> pd.DataFrame:
    rows = []
    def pair(label, level, a, b):
        d = pd.DataFrame({"a": a, "b": b}).dropna()
        d = d[d.index.str[:2].isin(UNIVERSE_STATES)]
        rows.append({"comparison": label, "geo_level": level, "n": len(d), "spearman": spearman(d["a"], d["b"]),
                     "median_a": d["a"].median() if len(d) else np.nan, "median_b": d["b"].median() if len(d) else np.nan})
    b = burden[burden["in_canonical_geographies"] | (burden["geo_level"] == "national")]

    def sel(cond, level, mid, year=None, role=None):
        s = b[(b["condition_id"] == cond) & (b["geo_level"] == level) & (b["measure_id"] == mid)]
        if year is not None:
            s = s[s["year"] == year]
        if role is not None:
            s = s[s["measure_role"] == role]
        return s.drop_duplicates("geo_id").set_index("geo_id")["value"]
    pair("HPS long COVID: primary (latest non-suppressed) vs 2024 mean", "state",
         sel("long_covid", "state", "lc_current_pct_all_adults", role="primary"),
         sel("long_covid", "state", "lc_current_pct_all_adults_2024_mean"))
    for cond, code in (("fibromyalgia", 51), ("migraine", 59)):
        for lvl in ("county", "state"):
            prim = sel(cond, lvl, f"cms_mmd_{code}_prev_actual_all", PRIMARY_CMS_YEAR)
            pair(f"CMS {code}: 2022 actual vs 2022 age-standardized", lvl, prim,
                 sel(cond, lvl, f"cms_mmd_{code}_prev_agestd_all", PRIMARY_CMS_YEAR))
            pair(f"CMS {code}: 2022 actual vs 2023 preliminary actual", lvl, prim,
                 sel(cond, lvl, f"cms_mmd_{code}_prev_actual_all", ALT_CMS_PRELIM_YEAR))
            pair(f"CMS {code}: 2022 actual all ages vs age < 65", lvl, prim,
                 sel(cond, lvl, f"cms_mmd_{code}_prev_actual_lt65", PRIMARY_CMS_YEAR))
    lp = sel("lyme_disease", "county", "lyme_incidence_per_100k", PRIMARY_LYME_YEAR)
    pair("Lyme county: 2023 vs pooled 2022-2023", "county", lp, sel("lyme_disease", "county", "lyme_incidence_per_100k_2022_2023"))
    pair("Lyme county: 2023 vs 2019", "county", lp, sel("lyme_disease", "county", "lyme_incidence_per_100k", ALT_LYME_PRE2022_YEAR))
    pair("PLACES PHLTH crude vs age-adjusted (county)", "county", sel("long_covid", "county", "places_PHLTH_crude"),
         sel("long_covid", "county", "places_PHLTH_ageadj"))
    pair("PLACES PHLTH crude vs COGNITION crude (county)", "county", sel("long_covid", "county", "places_PHLTH_crude"),
         sel("long_covid", "county", "places_COGNITION_crude"))
    return pd.DataFrame(rows)


def trial_radius_qa(features: pd.DataFrame) -> pd.DataFrame:
    """How often the 50 km internal-point radius misses trials sited inside the county (deviation 1)."""
    f = features[(features["geo_level"] == "county") & features["in_analysis_universe"]]
    rows = []
    for cond, g in f.groupby("condition_id"):
        miss = g["trials_in_geo_n"] > g["trials_within_50km_n"]
        rows.append({"condition_id": cond, "n_counties": len(g),
                     "n_counties_with_in_county_trials": int((g["trials_in_geo_n"] > 0).sum()),
                     "n_counties_in_county_exceeds_within_50km": int(miss.sum()),
                     "n_counties_union_exceeds_within_50km": int((g["trials_in_geo_or_within_50km_n"]
                                                                  > g["trials_within_50km_n"]).sum()),
                     "max_in_county_minus_within_50km": float((g["trials_in_geo_n"] - g["trials_within_50km_n"]).max())})
    return pd.DataFrame(rows)


def missingness_table(features: pd.DataFrame) -> pd.DataFrame:
    f = features[features["in_analysis_universe"]]
    cols = ["burden_value", "burden_proxy_value", "svi_overall", "relevant_providers_per_100k",
            "trials_in_geo_or_within_50km_n", "hrsa_sites_per_100k", "diagnostic_desert", "diagnostic_desert_access_only"]
    rows = []
    for (cond, lvl), g in f.groupby(["condition_id", "geo_level"]):
        r = {"condition_id": cond, "geo_level": lvl, "n": len(g)}
        for c in cols:
            r[f"missing_{c}"] = int(g[c].isna().sum())
        rows.append(r)
    return pd.DataFrame(rows)


# ============================================================================ figures (drafts)
def draft_figures(cov: pd.DataFrame, t4: dict, t6: dict, features: pd.DataFrame) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    DRAFTS.mkdir(parents=True, exist_ok=True)
    paths = []
    # coverage matrix
    lv = {"A": 4, "B": 3, "C": 2, "D": 1}
    m = cov.pivot(index="condition_id", columns="geo_level", values="best_level")[["national", "state", "county"]]
    inh = cov.pivot(index="condition_id", columns="geo_level", values="inherited")[["national", "state", "county"]]
    fig, ax = plt.subplots(figsize=(5.5, 6))
    ax.imshow(m.replace(lv).astype(float).to_numpy(), cmap="Greens", vmin=0, vmax=4.5, aspect="auto")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, m.iat[i, j] + ("*" if inh.iat[i, j] else ""), ha="center", va="center")
    ax.set_xticks(range(3), m.columns)
    ax.set_yticks(range(len(m)), m.index)
    ax.set_title("Best burden evidence level (* = inherited state value)", fontsize=9)
    fig.tight_layout()
    p = DRAFTS / "geo_burden_coverage_matrix.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p.relative_to(config.PROJECT_ROOT)))
    # Test 4 rho with population
    c = t4["test4_population_correlations"]
    c = c.assign(lab=c["series_id"] + " : " + c["variable"]).sort_values("rho_with_population")
    fig, ax = plt.subplots(figsize=(7, max(3, 0.28 * len(c))))
    ax.errorbar(c["rho_with_population"], range(len(c)), xerr=[c["rho_with_population"] - c["ci95_low"],
                                                                c["ci95_high"] - c["rho_with_population"]], fmt="o", ms=3)
    ax.axvline(0, color="grey", lw=0.8)
    for x in (0.3, 0.7, -0.3):
        ax.axvline(x, color="grey", lw=0.5, ls=":")
    ax.set_yticks(range(len(c)), c["lab"], fontsize=7)
    ax.set_xlabel("Spearman rho with total population (95% bootstrap CI)")
    ax.set_title("Test 4: burden vs population size", fontsize=9)
    fig.tight_layout()
    p = DRAFTS / "geo_test4_population_rho.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p.relative_to(config.PROJECT_ROOT)))
    # Test 6 relationships: observed vs null
    rel = t6["test6_geography_relationships"]
    nulls = t6["test6_geography_null_draws"]
    rids = rel["relationship_id"].unique()
    fig, axes = plt.subplots(1, len(rids), figsize=(2.4 * len(rids), 2.6), sharey=False)
    for ax, rid in zip(np.atleast_1d(axes), rids):
        for scheme, col in (("national", "tab:blue"), ("within_state", "tab:orange")):
            nv = nulls[(nulls["test"] == rid) & (nulls["scheme"] == scheme)]["null_value"]
            if len(nv):
                ax.hist(nv, bins=30, alpha=0.6, color=col, label=scheme)
        ax.axvline(rel[rel["relationship_id"] == rid]["observed_spearman"].iloc[0], color="k", lw=1.5)
        ax.set_title(rid, fontsize=8)
        ax.tick_params(labelsize=6)
    np.atleast_1d(axes)[0].legend(fontsize=6)
    fig.suptitle("Test 6: observed Spearman (black) vs permuted-label nulls", fontsize=9)
    fig.tight_layout()
    p = DRAFTS / "geo_test6_relationship_nulls.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p.relative_to(config.PROJECT_ROOT)))
    # Moran's I summary
    mo = t6["test6_geography_morans_i"]
    mo = mo[mo["weights"] == "contiguity"]
    piv = mo.pivot(index="variable", columns="permutation_scheme", values="null_mean")
    obs = mo.drop_duplicates("variable").set_index("variable")["observed_I"]
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    y = np.arange(len(obs))
    ax.scatter(obs.to_numpy(), y, label="observed", color="k")
    ax.scatter(piv.loc[obs.index, "national"], y, label="national shuffle (mean)", marker="x")
    ax.scatter(piv.loc[obs.index, "within_state"], y, label="within-state shuffle (mean)", marker="+")
    ax.set_yticks(y, obs.index, fontsize=7)
    ax.set_xlabel("Moran's I (county contiguity)")
    ax.legend(fontsize=6)
    fig.tight_layout()
    p = DRAFTS / "geo_test6_morans_i.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(str(p.relative_to(config.PROJECT_ROOT)))
    # desert map (me_cfs, county)
    try:
        import geopandas as gpd
        from ..config import PROCESSED
        poly = gpd.read_parquet(PROCESSED / "county_boundaries_2024.geoparquet")
        f = features[(features["condition_id"] == "me_cfs") & (features["geo_level"] == "county")]
        poly = poly.merge(f[["geo_id", "diagnostic_desert"]], on="geo_id", how="inner")
        poly = poly[poly["geo_id"].str[:2].isin(set(UNIVERSE_STATES) - {"02", "15"})].to_crs(5070)
        fig, ax = plt.subplots(figsize=(9, 5.5))
        poly.plot(column="diagnostic_desert", ax=ax, cmap="magma_r", legend=True, linewidth=0,
                  missing_kwds={"color": "lightgrey"})
        ax.set_axis_off()
        ax.set_title("Draft: diagnostic desert index, ME/CFS (burden = CMS 51, level C proxy); grey = unknown", fontsize=9)
        fig.tight_layout()
        p = DRAFTS / "geo_desert_map_me_cfs.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        paths.append(str(p.relative_to(config.PROJECT_ROOT)))
    except Exception as e:  # map is a convenience; failure is reported, not hidden
        paths.append(f"map failed: {e}")
    return paths


# ============================================================================ results document
def _fmt(x, nd=2):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return f"{x:,}"
    return f"{x:.{nd}f}"


def _md(df: pd.DataFrame, cols: list[str], nd: int = 2) -> str:
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df[cols].iterrows():
        vals = []
        for v in r.to_numpy():
            if isinstance(v, (float, np.floating)):
                vals.append(_fmt(float(v), nd))
            elif isinstance(v, (bool, np.bool_)):
                vals.append(str(bool(v)))
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_results_md(stats_: dict) -> Path:
    """Generate results/GEOGRAPHY_RESULTS.md from the CSVs and tables this module wrote."""
    T = {p.stem: pd.read_csv(p, dtype={"geo_id": str}) for p in TABLES.glob("*.csv")
         if p.stem.startswith(("test4_", "test6_geography_", "geo_"))}
    cov, corr, ov = T["geo_burden_coverage_matrix"], T["test4_population_correlations"], T["test4_topn_overlap"]
    part, dec, ver = T["test4_partial_correlations"], T["test4_log_count_decomposition"], T["test4_verdicts"]
    rel, mor = T["test6_geography_relationships"], T["test6_geography_morans_i"]
    sens, alt = T["geo_desert_sensitivity"], T["geo_primary_vs_alternates"]
    top_n = load_config("scoring")["ranking"]["top_n"]
    s = stats_
    L = []
    L.append("# Geography results: burden, context, opportunity features, Test 4 and Test 6 (geography labels)\n")
    L.append(f"Generated by `{PRODUCER}` at {s['generated_at']} from files it wrote (tables under `results/tables/`, "
             "processed tables under `data/processed/`). Plan: `docs/ANALYSIS_PLAN_GEOGRAPHY.md`; definitions: "
             "`docs/BURDEN_DEFINITIONS.md`. All joins are ecological: nothing here describes an individual.\n")
    if s.get("reviewer_reading"):
        L.append("## 0. Independent review: corrections to the first version (2026-09-23)\n")
        L.append("The review re-ran the module end to end (identical tables) and recomputed headline numbers "
                 "independently from the input partitions and with `esda`/`libpysal` (Moran's I). The arithmetic "
                 "was correct. Several inferences were stronger than the evidence; they are corrected below and in "
                 "the sections that follow (all numbers pulled from the tables):\n")
        L.append(s["reviewer_reading"] + "\n")
    L.append("## 1. What was built\n")
    L.append(_md(pd.DataFrame(s["tables"]), ["table", "rows", "note"]))
    L.append("")
    L.append("## 2. Evidence-level coverage matrix (condition x geography level x best level)\n")
    L.append("`*` = the county value is the state estimate carried as inherited context (source resolution `state`), "
             "not county prevalence. Coverage = geographies in the 50 states + DC with a primary value.\n")
    cv = cov.copy()
    cv["coverage"] = cv["n_geos_with_value"].astype(int).astype(str) + " / " + cv["n_geos_in_universe"].astype(int).astype(str)
    cv["best"] = cv["best_level"] + np.where(cv["inherited"], "*", "")
    L.append(_md(cv, ["condition_id", "geo_level", "best", "primary_measure_id", "source_resolution",
                      "county_proxy_column", "coverage"]))
    nd = (cov["best_level"] == "D").groupby(cov["condition_id"]).all()
    L.append(f"\n**{int(nd.sum())} of {len(nd)} target conditions have no usable burden measure at any level "
             f"(level D): {', '.join(sorted(nd[nd].index))}.** Their burden is null and their diagnostic desert is "
             "null; only `diagnostic_desert_access_only` (not burden-informed) is defined for them. No level-A "
             "county measure exists for any condition except Lyme disease; long COVID at county level is inherited "
             "state context.\n")
    if s.get("coverage_reading"):
        L.append(s["coverage_reading"] + "\n")
    L.append("## 3. Primary vs pre-specified alternate burden measures\n")
    L.append(_md(alt, ["comparison", "geo_level", "n", "spearman", "median_a", "median_b"]))
    L.append("\n" + s.get("alternates_reading", "") + "\n")
    # Test 4
    L.append("## 4. Test 4 — are hotspots driven by burden or by population size?\n")
    L.append("Spearman rho with total population (ACS 2020-2024), 95% percentile-bootstrap CI (2,000 resamples). "
             "`ci95_*` resample geographies i.i.d. (pre-specified); `cl_*` resample whole states with their counties "
             "(reviewer addition; counties are spatially dependent, so the i.i.d. interval is too narrow):\n")
    cr = corr.rename(columns={"ci95_low_state_cluster": "cl_low", "ci95_high_state_cluster": "cl_high"})
    L.append(_md(cr, ["series_id", "variable", "n", "rho_with_population", "ci95_low", "ci95_high", "cl_low",
                      "cl_high"]))
    L.append(test4_series_legend(corr))
    L.append(f"\nTop-N overlap with the population ranking (counties: N = {top_n}, population >= "
             f"{load_config('scoring')['ranking']['min_population']:,}; states: N = {TOP_N_STATE}); expected = N^2/M under "
             "independence; P = hypergeometric P(overlap >= observed). `overlap` breaks ties by FIPS (pre-specified); "
             "`tie_avg` is the expected overlap under random tie-breaking and `ties_a` the number of values tied at "
             "the top-N boundary of ranking_a (CMS publishes integer percentages):\n")
    ovp = ov[(ov["N"].isin([top_n, TOP_N_STATE])) & ((ov["geo_level"] == "state") | (ov["population_floor"] > 0))]
    ovp = ovp.rename(columns={"overlap_tie_averaged": "tie_avg", "ties_at_boundary_a": "ties_a"})
    L.append(_md(ovp, ["series_id", "ranking_a", "ranking_b", "M", "N", "overlap", "tie_avg", "ties_a",
                       "expected_overlap", "p_hypergeom"], 3))
    L.append("\nPre-specified verdicts (population-driven: rho > 0.7 and top-N overlap above chance at P < 0.05; "
             "burden-driven: |rho| < 0.3 and overlap not above chance; otherwise mixed). The pre-specified rule only "
             "tests the LARGE-population tail; `in_bottom_q` (reviewer addition) is how many of the top-N fall in the "
             "bottom population quartile (expected `exp_bq`), with a Holm-adjusted hypergeometric P across the "
             "family (`p_bq_holm`), and `reviewer_reading` combines both tails:\n")
    vv = ver.rename(columns={"topN_in_bottom_pop_quartile": "in_bottom_q", "expected_in_bottom_quartile": "exp_bq",
                             "p_hypergeom_bottom_quartile_holm": "p_bq_holm"})
    L.append(_md(vv, ["series_id", "variable", "rho_with_population", "topN_overlap_with_population",
                      "expected_overlap", "verdict", "in_bottom_q", "exp_bq", "p_bq_holm", "reviewer_reading"]))
    L.append("\nPartial Spearman correlations controlling for log population:\n")
    L.append(_md(part, ["series_id", "x", "y", "n", "raw_spearman", "partial_spearman"]))
    if len(dec):
        L.append("\nShare of variance of log(count) explained by log(population) vs by log(prevalence) (counts > 0):\n")
        L.append(_md(dec, ["series_id", "n_count_positive", "n_count_zero", "r2_logcount_on_logpop",
                           "r2_logcount_on_logprev"]))
    L.append("\n**Reading (numbers pulled from the tables above).**\n\n" + s["test4_reading"] + "\n")
    # Test 6
    L.append("## 5. Test 6 — shuffled geography labels (negative controls)\n")
    L.append(f"Spearman relationships between independent sources, observed vs {N_PERM:,} permutations of the first "
             "variable (national shuffle; and within-state shuffle for counties). `p` is the permutation P used for "
             "the row: for the national shuffle the pre-specified P around 0; for the within-state shuffle the P "
             "around the null's own centre (`p_perm_vs_null_centre`, reviewer correction: the pre-specified P around "
             "0 does not test anything for a null centred on the between-state component, e.g. it gave 0.992 for R1 "
             "although the observed value lies below the within-state null band). The null mean of the "
             "within-state shuffle is the between-state component:\n")
    rl = rel.copy()
    rl["p"] = np.where(rl["permutation_scheme"] == "national", rl["p_perm"], rl["p_perm_vs_null_centre"])
    L.append(_md(rl, ["relationship_id", "description", "n", "observed_spearman", "permutation_scheme", "null_mean",
                      "null_p2_5", "null_p97_5", "p", "p_holm", "within_state_component", "interpretation"], 3))
    if "n_eff_spatial" in rel.columns:
        L.append("\nSpatially adjusted inference for the national-shuffle rows (reviewer addition). A label shuffle "
                 "treats counties as exchangeable, but both variables are spatially autocorrelated (Moran's I below), "
                 "so its null is too narrow. `n_eff` is the Clifford-Richardson-Hemon effective sample size (100 km "
                 "distance classes for counties, 250 km for states), `p_crh` the corresponding t-test P, `p_crh_holm` "
                 "Holm-adjusted across R1-R6, and `mdr80` the smallest |rho| detectable with 80% power at n_eff:\n")
        nat = rel[rel["permutation_scheme"] == "national"].rename(columns={
            "n_eff_spatial": "n_eff", "p_spatial_crh": "p_crh", "p_spatial_crh_holm": "p_crh_holm",
            "min_detectable_rho_80pct": "mdr80"})
        L.append(_md(nat, ["relationship_id", "n", "n_eff", "observed_spearman", "p_perm", "p_crh", "p_crh_holm",
                           "mdr80", "spatially_adjusted_reading"], 3))
    L.append(f"\nMoran's I (row-standardised; {N_PERM_MORAN} permutations; E[I] = -1/(n-1)). `p` as above (national: "
             "around E[I]; within-state: around the null centre); `p_holm` within each scheme (18 tests):\n")
    mo = mor.copy()
    mo["p"] = np.where(mo["permutation_scheme"] == "national", mo["p_perm"], mo["p_perm_vs_null_centre"])
    L.append(_md(mo, ["variable", "weights", "n", "observed_I", "permutation_scheme", "null_mean", "null_p2_5",
                      "null_p97_5", "p", "p_holm", "retained_fraction"], 3))
    L.append("\n**Reading (numbers pulled from the tables above).**\n\n" + s["test6_reading"] + "\n")
    L.append("## 6. Diagnostic desert: definition and sensitivity\n")
    L.append("`diagnostic_desert = mean(pr(burden), pr(SVI overall), 1 - pr(relevant providers per 100k), "
             "1 - pr(relevant trials sited in the county or within 50 km of its internal point))` (state: trials in "
             "the state), percentile ranks within 50 states + DC "
             "per condition and level. Every component is a column of `geo_condition_features`. Sensitivity of the "
             f"ranking to pre-specified variants (Spearman with the primary; top-{top_n} overlap among counties with "
             "population >= 10,000):\n")
    L.append(_md(sens, ["condition_id", "geo_level", "variant", "n", "spearman_with_primary", "topN", "topN_overlap"]))
    L.append("\n" + s["desert_reading"] + "\n")
    L.append(s.get("provider_reading", "") + "\n")
    L.append("## 7. Caveats that limit every number above\n")
    for c in s["caveats"]:
        L.append(f"* {c}")
    L.append("\n## 8. Files\n")
    for f in s["files"]:
        L.append(f"* `{f}`")
    path = config.RESULTS / "GEOGRAPHY_RESULTS.md"
    path.write_text("\n".join(L) + "\n")
    return path


def _readings(t4: dict, t6: dict, sens: pd.DataFrame, extra: dict) -> dict:
    """Plain-language readings assembled from the result tables (every number is pulled, never typed)."""
    corr, ver, ov = t4["test4_population_correlations"], t4["test4_verdicts"], t4["test4_topn_overlap"]
    part, dec = t4["test4_partial_correlations"], t4["test4_log_count_decomposition"]
    top_n = load_config("scoring")["ranking"]["top_n"]

    def c(sid, var):
        r = corr[(corr["series_id"] == sid) & (corr["variable"] == var)]
        return None if r.empty else r.iloc[0]

    def v(sid, var):
        r = ver[(ver["series_id"] == sid) & (ver["variable"] == var)]
        return None if r.empty else r.iloc[0]

    lines = []
    for sid, lab in (("long_covid_state", "Long COVID, state (HPS, level A)"),
                     ("lyme_county", "Lyme disease, county (level A)"),
                     ("lyme_state", "Lyme disease, state aggregate (level A)"),
                     ("long_covid_county_proxy", "PLACES PHLTH, county (level C proxy for long COVID / ME/CFS)"),
                     ("cms51_county", "CMS 51, county (fibromyalgia B / ME/CFS C)"),
                     ("cms59_county", "CMS 59, county (migraine B)")):
        a, va = c(sid, "prevalence"), v(sid, "prevalence")
        if a is None:
            continue
        t = (f"* **{lab}.** Prevalence/incidence vs population: rho = {a['rho_with_population']:.2f} "
             f"[{a['ci95_low']:.2f}, {a['ci95_high']:.2f}], top-N overlap with the population ranking "
             f"{_fmt(va['topN_overlap_with_population'], 0)} (expected {va['expected_overlap']:.2f}) -> {va['verdict']}.")
        b, vb = c(sid, "count"), v(sid, "count")
        if b is not None:
            t += (f" Count vs population: rho = {b['rho_with_population']:.2f} [{b['ci95_low']:.2f}, "
                  f"{b['ci95_high']:.2f}], top-N overlap {_fmt(vb['topN_overlap_with_population'], 0)} "
                  f"(expected {vb['expected_overlap']:.2f}, P = {vb['p_hypergeom']:.3g}) -> {vb['verdict']}.")
            d = dec[dec["series_id"] == sid]
            if len(d):
                d = d.iloc[0]
                t += (f" Variance of log(count) explained by log(population) {d['r2_logcount_on_logpop']:.2f} vs by "
                      f"log(prevalence) {d['r2_logcount_on_logprev']:.2f} (n = {int(d['n_count_positive']):,} with "
                      "count > 0).")
        else:
            t += " No count: not defensible at this resolution."
        lines.append(t)
    dv = ver[ver["variable"] == "desert"].drop_duplicates(["series_id"])
    dcorr = corr[corr["variable"] == "desert"].drop_duplicates("series_id")
    dov = ov[(ov["ranking_a"] == "desert") & (ov["ranking_b"] == "population")
             & (((ov["geo_level"] == "county") & (ov["population_floor"] > 0) & (ov["N"] == top_n))
                | ((ov["geo_level"] == "state") & (ov["N"] == TOP_N_STATE)))]
    lines.append(f"* **Diagnostic desert** ({len(dcorr)} condition x level series): rho with population ranges "
                 f"{dcorr['rho_with_population'].min():.2f} to {dcorr['rho_with_population'].max():.2f}; top-N overlap "
                 f"with the population ranking {int(dov['overlap'].min())} to {int(dov['overlap'].max())} "
                 f"({int((dov['p_hypergeom'] < 0.05).sum())} of {len(dov)} above chance at P < 0.05). Pre-specified "
                 "verdicts: " + ", ".join(f"{k} {n}" for k, n in dv["verdict"].value_counts().items()) + ".")
    sm = t4.get("test4_small_population_check")
    if sm is not None and len(sm):
        dsm = sm[sm["series_id"].str.startswith("desert_") & (sm["geo_level"] == "county")]
        sig = dsm[dsm["p_hypergeom_bottom_quartile_holm"] < 0.05]
        nonsig = dsm[dsm["p_hypergeom_bottom_quartile_holm"] >= 0.05]
        dc = corr[corr["series_id"].str.startswith("desert_") & (corr["geo_level"] == "county")]
        lines.append(
            "* **Reviewer correction: the county desert short lists are concentrated in small-population counties, "
            "which the pre-specified test cannot see** (it only compares with the LARGEST-population ranking). Of the "
            f"top {int(dsm['N'].iloc[0])} counties (population >= 10,000), "
            f"{int(dsm['topN_in_bottom_pop_quartile'].min())} to {int(dsm['topN_in_bottom_pop_quartile'].max())} fall "
            f"in the bottom population quartile (expected {dsm['expected_in_bottom_quartile'].iloc[0]:.1f}); "
            f"significant after Holm for {len(sig)} of {len(dsm)} county series ("
            + ", ".join(f"{r.series_id.replace('desert_', '').replace('_county', '')} {int(r.topN_in_bottom_pop_quartile)}"
                        f"/{int(r.N)}, P_Holm = {r.p_hypergeom_bottom_quartile_holm:.2g}" for r in sig.itertuples())
            + "); not significant for " + (", ".join(r.series_id.replace("desert_", "").replace("_county", "")
                                                     + f" {int(r.topN_in_bottom_pop_quartile)}/{int(r.N)}"
                                                     for r in nonsig.itertuples()) or "none")
            + f". Median population of the top-{int(dsm['N'].iloc[0])}: "
            f"{int(dsm['topN_median_population'].min()):,} to {int(dsm['topN_median_population'].max()):,} vs "
            f"{int(dsm['all_median_population'].iloc[0]):,} overall. County desert rho with population "
            f"{dc['rho_with_population'].min():.2f} to {dc['rho_with_population'].max():.2f} (state-cluster bootstrap "
            f"CIs spanning {dc['ci95_low_state_cluster'].min():.2f} to {dc['ci95_high_state_cluster'].max():.2f}). The access "
            "components drive this: the trial component is a raw count of accessible trials (not per capita) and "
            "per-capita relevant-provider density is lower in small counties. So 'not population-driven' holds only "
            "for the large-population tail; the county desert partly measures smallness/rurality. State-level desert "
            "series show no small-population concentration.")
        px_ = sm[(sm["series_id"] == "long_covid_county_proxy") & (sm["variable"] == "prevalence")]
        if len(px_):
            r = px_.iloc[0]
            lines.append(f"* The PLACES PHLTH county proxy has the same small-county concentration: "
                         f"{int(r['topN_in_bottom_pop_quartile'])}/{int(r['N'])} of its top counties in the bottom "
                         f"population quartile (P_Holm = {r['p_hypergeom_bottom_quartile_holm']:.2g}).")
        cm = ov[ov["series_id"].isin(["cms51_county", "cms59_county"]) & (ov["ranking_b"] == "population")
                & (ov["population_floor"] > 0) & (ov["N"] == top_n)]
        if len(cm) and "ties_at_boundary_a" in cm.columns:
            lines.append("* CMS top-N lists are partly decided by ties: CMS publishes integer percentages, so "
                         + "; ".join(f"{r.series_id} has {int(r.ties_at_boundary_a)} counties tied at the top-{int(r.N)} "
                                     f"boundary (tie-averaged overlap with population {r.overlap_tie_averaged:.2f} vs "
                                     f"{int(r.overlap)} with the FIPS tie-break)" for r in cm.itertuples())
                         + ". Their top-N memberships should not be read as a ranking.")
    pd_ = part[(part["x"] == "desert") & (part["y"] == "burden")].drop_duplicates("series_id")
    neg = pd_[pd_["partial_spearman"] < 0]
    lines.append(f"* Desert vs its own burden component, controlling for log population: partial rho "
                 f"{pd_['partial_spearman'].min():.2f} to {pd_['partial_spearman'].max():.2f}"
                 + (f"; negative for {', '.join(neg['series_id'])} (high-burden places there have more providers/"
                    "trials and lower SVI, so the access components pull the index the other way)." if len(neg) else "."))
    px = part[(part["series_id"] == "long_covid_county_proxy") & (part["y"] == "svi_overall")]
    if len(px):
        lines.append(f"* The long-COVID county C-proxy (PLACES PHLTH) tracks social vulnerability: Spearman with SVI "
                     f"{px['raw_spearman'].iloc[0]:.2f} (partial {px['partial_spearman'].iloc[0]:.2f}); used as burden it "
                     "would largely double-count the vulnerability component.")
    t4r = "\n".join(lines)

    rel, mor = t6["test6_geography_relationships"], t6["test6_geography_morans_i"]
    L = []
    nat = rel[rel["permutation_scheme"] == "national"]
    real = nat[nat["observed_outside_null95"]]
    none_ = nat[~nat["observed_outside_null95"]]
    L.append("* **Cross-source relationships, national shuffle (naive; see the spatial correction below).** Structure "
             "that collapses under the shuffle: "
             + (", ".join(f"{r.relationship_id} (rho {r.observed_spearman:.3f}, null mean {r.null_mean:.3f}, P = "
                          f"{r.p_perm:.3f})" for r in real.itertuples()) or "none")
             + ". No detectable relationship to begin with (observed inside the null band): "
             + (", ".join(f"{r.relationship_id} (rho {r.observed_spearman:.3f}, P = {r.p_perm:.3f})"
                          for r in none_.itertuples()) or "none") + ".")
    if "p_spatial_crh" in nat.columns:
        adj_ok = nat[nat["observed_outside_null95"] & (nat["p_spatial_crh_holm"] < 0.05)]
        adj_bad = nat[nat["observed_outside_null95"] & ~(nat["p_spatial_crh_holm"] < 0.05)]
        L.append("* **Reviewer correction: spatially adjusted inference.** The national shuffle ignores spatial "
                 "autocorrelation, so its P values are anti-conservative. With the Clifford-Richardson-Hemon effective "
                 "sample size and Holm correction across R1-R6: "
                 + ("; ".join(f"{r.relationship_id} survives (n = {int(r.n):,} -> n_eff = {r.n_eff_spatial:.0f}, "
                              f"P_CRH = {r.p_spatial_crh:.2g}, P_Holm = {r.p_spatial_crh_holm:.2g})"
                              for r in adj_ok.itertuples()) or "nothing survives")
                 + "; " + ("; ".join(f"{r.relationship_id} does NOT survive (rho {r.observed_spearman:.3f}; n = "
                                     f"{int(r.n):,} -> n_eff = {r.n_eff_spatial:.0f}, P_CRH = {r.p_spatial_crh:.2g}, "
                                     f"P_Holm = {r.p_spatial_crh_holm:.2g})" for r in adj_bad.itertuples())
                           or "no naive-significant relationship is lost")
                 + ". Only the same-source positive control is robust; no cross-source relationship (R1-R3, R5, R6) "
                   "survives. The null results are also weakly powered: the smallest |rho| detectable with 80% power "
                   "at n_eff is "
                 + ", ".join(f"{r.relationship_id} {r.min_detectable_rho_80pct:.2f}" for r in nat.itertuples()) + ".")
        cal = t6.get("test6_geography_crh_calibration")
        if cal is not None and len(cal):
            L.append("* **Calibration of the two tests** (`test6_geography_crh_calibration.csv`): for two INDEPENDENT "
                     f"spatially autocorrelated fields on {int(cal['n_points'].iloc[0]):,} real county points "
                     f"({int(cal['n_sim'].iloc[0])} simulations per field range), the naive exchangeable test rejects "
                     "at nominal 5% in "
                     + ", ".join(f"{100 * r.false_positive_rate_naive:.0f}%" for r in cal.itertuples())
                     + " of cases (field ranges " + ", ".join(f"{int(r.field_range_km)} km" for r in cal.itertuples())
                     + "), the CRH-adjusted test in "
                     + ", ".join(f"{100 * r.false_positive_rate_crh:.1f}%" for r in cal.itertuples())
                     + ". A national-shuffle P value alone therefore cannot establish a cross-source relationship "
                       "between spatially smooth layers.")
    ws = rel[rel["permutation_scheme"] == "within_state"]
    L.append("* **Within-state shuffle (counties).** The null mean is the between-state component; the P around the "
             "null centre tests for a within-state association (reviewer correction; the pre-specified P around 0 "
             "is not informative here). " + "; ".join(
                 f"{r.relationship_id}: between-state component {r.null_mean:.3f}, within-state component "
                 f"{r.within_state_component:.3f} (P = {r.p_perm_vs_null_centre:.3f}, Holm {r.p_holm:.3f})"
                 for r in ws.itertuples())
             + ". These within-state P values still treat counties within a state as exchangeable, so they share the "
               "spatial caveat above.")
    mn = mor[mor["permutation_scheme"] == "national"]
    mw = mor[mor["permutation_scheme"] == "within_state"]
    L.append(f"* **Moran's I.** Observed I ranges {mn['observed_I'].min():.3f} to {mn['observed_I'].max():.3f} across "
             f"{mn['variable'].nunique()} variables x {mn['weights'].nunique()} weight schemes; the national shuffle "
             f"collapses it to the no-structure value in {int(mn['collapsed_under_permutation'].sum())} of {len(mn)} cases "
             f"(null means {mn['null_mean'].min():.3f} to {mn['null_mean'].max():.3f}; "
             + (f"{int((mn['p_holm'] < 0.05).sum())} of {len(mn)} significant after Holm" if "p_holm" in mn else "")
             + f"). Under the within-state shuffle the retained fraction (between-state share) is "
             f"{mw['retained_fraction'].min():.2f} to {mw['retained_fraction'].max():.2f}"
             + (f"; the observed I exceeds the within-state null (Holm P < 0.05) for "
                f"{int(((mw['p_holm'] < 0.05) & (mw['observed_I'] > mw['null_mean'])).sum())} of {len(mw)} "
                "variable x weight cases, i.e. those variables have spatial structure beyond between-state "
                "differences; the exceptions are "
                + (", ".join(sorted(set(mw.loc[~((mw['p_holm'] < 0.05) & (mw['observed_I'] > mw['null_mean'])),
                                               'variable']))) or "none")
                if "p_holm" in mw else "") + ".")
    inh = mw[mw["variable"] == "long_covid_inherited_state_value"]
    if len(inh):
        L.append(f"* The inherited long-COVID county value keeps Moran's I = {inh['null_mean'].iloc[0]:.3f} "
                 f"(= observed {inh['observed_I'].iloc[0]:.3f}) after every within-state shuffle: it carries no "
                 "within-state information, which is why it is flagged inherited (`burden_inherited`). Reviewer "
                 "correction: it IS used, as pre-specified, as the burden component of the long-COVID county "
                 "diagnostic desert, where every county in a state receives the same burden percentile "
                 "(`desert_burden_resolution` = 'state (inherited; no within-state variation)'). That county desert "
                 "therefore ranks counties partly by their state's HPS estimate and must not be read as evidence of "
                 "county long-COVID burden; its within-state ordering comes only from SVI, provider density and trials.")
    pdn = mor[(mor["variable"] == "log1p_relevant_provider_density_long_covid")]
    if len(pdn):
        a = pdn[(pdn["permutation_scheme"] == "national") & (pdn["weights"] == "contiguity")].iloc[0]
        L.append(f"* Relevant-provider density has weak spatial autocorrelation (I = {a['observed_I']:.3f}, contiguity), "
                 "much weaker than the burden and vulnerability layers.")
    t6r = "\n".join(L)

    s = sens[sens["variant"] != "diagnostic_desert_access_only"]
    ao = sens[sens["variant"] == "diagnostic_desert_access_only"]
    dr = (f"**Reading.** Across {len(s)} condition x level x variant comparisons the Spearman correlation with the "
          f"primary desert is {s['spearman_with_primary'].median():.2f} (median; range "
          f"{s['spearman_with_primary'].min():.2f} to {s['spearman_with_primary'].max():.2f}), but top-N overlap ranges "
          f"{int(s['topN_overlap'].min())} to {int(s['topN_overlap'].max())} of N: the overall ordering is stable while "
          "the short list at the top is not. Dropping burden (access-only index) keeps "
          f"{int(ao['topN_overlap'].min())} to {int(ao['topN_overlap'].max())} of the top-N, so the burden component "
          "decides which places reach the top. ")
    qa = extra.get("trial_radius_qa")
    if qa is not None and len(qa):
        dr += (f"Deviation 1 check: the union trial measure exceeds the internal-point-only 50 km count in "
               f"{int(qa['n_counties_union_exceeds_within_50km'].min())} to "
               f"{int(qa['n_counties_union_exceeds_within_50km'].max())} counties per condition; the internal-point-only "
               "variant correlates with the primary at "
               f"{sens[sens['variant'].str.endswith('internal_point_only')]['spearman_with_primary'].min():.3f} or higher.")
    alt = extra.get("alternates")
    ar = ""
    if alt is not None and len(alt):
        h = alt[alt["comparison"].str.startswith("HPS")]
        u65 = alt[alt["comparison"].str.contains("age < 65") & (alt["geo_level"] == "county")]
        ar = ("**Reading.** " + (f"The state long-COVID ranking is period-sensitive: primary (latest non-suppressed "
                                 f"period) vs the 2024 mean, Spearman {h['spearman'].iloc[0]:.2f} (n = {int(h['n'].iloc[0])}); "
                                 "single HPS periods have wide CIs, so state ranks should be read with their CIs. "
                                 if len(h) else "")
              + (f"The CMS measures depend strongly on the age stratum (all ages vs < 65, county Spearman "
                 f"{u65['spearman'].min():.2f} to {u65['spearman'].max():.2f}), whereas actual vs age-standardized and "
                 "2022 vs 2023 agree closely (table)." if len(u65) else ""))
    f = extra.get("features")
    cov_r = ""
    if f is not None:
        u = f[f["in_analysis_universe"] & (f["geo_level"] == "county") & (f["condition_id"] == "migraine")]
        miss = u[u["burden_value"].isna()]
        by = miss["state_fips"].map(STATE_FIPS).fillna(miss["state_fips"]).value_counts()
        cov_r = (f"CMS county coverage gaps ({len(miss)} of {len(u):,} counties without a 2022 CMS value): "
                 + ", ".join(f"{k} {v}" for k, v in by.items())
                 + ". Connecticut's 9 planning regions are missing because CMS publishes legacy CT counties; the others "
                   "are not published by the source (see `burden_unknown_reason`).")
    # reviewer summary (every number pulled from the tables)
    rv = []
    sm = t4.get("test4_small_population_check")
    if sm is not None and len(sm):
        dsm = sm[sm["series_id"].str.startswith("desert_") & (sm["geo_level"] == "county")]
        rv.append(f"* **Test 4, diagnostic desert.** The pre-specified verdict 'not population-driven' tests only the "
                  f"large-population tail. The county desert top-{int(dsm['N'].iloc[0])} lists hold "
                  f"{int(dsm['topN_in_bottom_pop_quartile'].min())}-{int(dsm['topN_in_bottom_pop_quartile'].max())} "
                  f"bottom-population-quartile counties (expected {dsm['expected_in_bottom_quartile'].iloc[0]:.1f}; "
                  f"Holm P < 0.05 for {int((dsm['p_hypergeom_bottom_quartile_holm'] < 0.05).sum())} of {len(dsm)} "
                  "series). The county desert partly measures small population/rurality via its access components "
                  "(`test4_small_population_check.csv`, `reviewer_reading` in `test4_verdicts.csv`).")
    nat_ = rel[rel["permutation_scheme"] == "national"]
    if "p_spatial_crh_holm" in nat_.columns:
        lost = nat_[nat_["observed_outside_null95"] & ~(nat_["p_spatial_crh_holm"] < 0.05)]
        rv.append("* **Test 6, cross-source relationships.** The national shuffle ignores spatial autocorrelation. "
                  "After a Clifford-Richardson-Hemon effective-sample-size correction and Holm across R1-R6, "
                  + ("; ".join(f"{r.relationship_id} (rho {r.observed_spearman:.3f}, naive P {r.p_perm:.3f}) drops to "
                               f"P_CRH = {r.p_spatial_crh:.2g}, P_Holm = {r.p_spatial_crh_holm:.2g} (n_eff "
                               f"{r.n_eff_spatial:.0f} of {int(r.n):,})" for r in lost.itertuples()) or "no relationship changes")
                  + ". Only the same-source positive control R4 survives."
                  + ((" In simulation the naive test rejects "
                      + "/".join(f"{100 * r.false_positive_rate_naive:.0f}%" for r in
                                 t6["test6_geography_crh_calibration"].itertuples())
                      + " of independent smooth-field pairs at nominal 5% (CRH: "
                      + "/".join(f"{100 * r.false_positive_rate_crh:.1f}%" for r in
                                 t6["test6_geography_crh_calibration"].itertuples()) + ").")
                     if "test6_geography_crh_calibration" in t6 else ""))
    ws_ = rel[rel["permutation_scheme"] == "within_state"]
    if "p_perm_vs_null_centre" in ws_.columns and len(ws_):
        r1 = ws_[ws_["relationship_id"] == "R1"]
        rv.append("* **Test 6, within-state P values.** The pre-specified P (around 0) is uninformative for a null "
                  "centred on the between-state component"
                  + (f" (R1: {r1['p_perm'].iloc[0]:.3f} although the observed value lies below the null band; "
                     f"centred P = {r1['p_perm_vs_null_centre'].iloc[0]:.3f})" if len(r1) else "")
                  + ". Within-state rows now report the P around the null centre, and Moran's I rows likewise.")
    rv.append("* **Long-COVID county desert.** The earlier text said the inherited state value is 'never ranked as "
              "county burden'. That was false: it is the burden component of the long-COVID county desert, with "
              "no within-state variation. The index is now flagged (`desert_burden_resolution`, `desert_definition`) "
              "and the text corrected.")
    cc = corr[corr["ci95_low_state_cluster"].notna()] if "ci95_low_state_cluster" in corr else corr.iloc[0:0]
    if len(cc):
        ratio = ((cc["ci95_high_state_cluster"] - cc["ci95_low_state_cluster"]) / (cc["ci95_high"] - cc["ci95_low"]))
        iid_excl = (cc["ci95_low"] > 0) | (cc["ci95_high"] < 0)
        cl_excl = (cc["ci95_low_state_cluster"] > 0) | (cc["ci95_high_state_cluster"] < 0)
        rv.append(f"* **Uncertainty and multiplicity.** County Test 4 correlations now also carry state-cluster "
                  f"bootstrap CIs, {ratio.min():.1f}x to {ratio.max():.1f}x wider than the i.i.d. CIs (median "
                  f"{ratio.median():.1f}x). Of {int(iid_excl.sum())} county correlations whose i.i.d. CI excluded 0, "
                  f"{int((iid_excl & ~cl_excl).sum())} no longer do. Holm-adjusted P values are added for the "
                  "Test 4 verdict family, the Test 6 relationship families and the Moran's I families. Top-N overlaps "
                  "also report tie-averaged values, because CMS publishes integer percentages.")
    rv.append("* **Query layer.** Condition names resolve by unique prefix ('Lyme' -> lyme_disease), place names are "
              "accent-insensitive ('Dona Ana County, NM'), 'US' resolves to the national level, and a geography id "
              "at the wrong level or a missing value now returns the specific reason.")
    reviewer = "\n".join(rv)
    pr_ = ""
    if f is not None:
        st = f[(f["geo_level"] == "state") & f["in_analysis_universe"]]
        g = st.groupby("condition_id")[["relevant_specialists_n", "relevant_providers_n"]].sum()
        share = (g["relevant_specialists_n"] / g["relevant_providers_n"]).dropna()
        pr_ = (f"**Provider counts.** Primary care dominates the relevant-provider variant: non-primary-care specialists "
               f"are {100 * share.min():.1f}% to {100 * share.max():.1f}% of the relevant individual NPIs across the "
               f"{len(share)} conditions (state sums), so the primary density mostly measures primary-care supply; the "
               "specialists-only variant (`relevant_specialists_per_100k`) is kept for that reason.")
    return {"test4_reading": t4r, "test6_reading": t6r, "desert_reading": dr, "alternates_reading": ar,
            "provider_reading": pr_, "coverage_reading": cov_r, "reviewer_reading": reviewer}


# ============================================================================ orchestration
def _write_csv(df: pd.DataFrame, name: str) -> str:
    TABLES.mkdir(parents=True, exist_ok=True)
    p = TABLES / f"{name}.csv"
    df.to_csv(p, index=False)
    return str(p.relative_to(config.PROJECT_ROOT))


def build_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    burden = build_burden()
    write_table(burden, "geo_condition_burden", producer=PRODUCER,
                description="Harmonised long burden table (HPS, CDC Lyme, CMS MMD, PLACES level-C proxies, inherited "
                            "state long-COVID context) with evidence level, measure role and inherited flag")
    defs = definitions_frame().rename(columns={"source_geographic_resolution": "primary_source_resolution"})
    defs = defs.assign(data_layer="metadata", source_name="docs/BURDEN_DEFINITIONS.md (pre-specified 2026-09-23)",
                       source_record_id=defs["condition_id"] + "|" + defs["geo_level"], source_version="v1",
                       retrieved_at=utc_now_iso(), source_geographic_resolution="none",
                       evidence_type="curated_config", evidence_level=defs["burden_evidence_level"],
                       provenance_notes="Pre-specified primary burden measure per condition x geography level.")
    write_table(defs, "geo_burden_definitions", producer=PRODUCER,
                description="Pre-specified primary burden measure per condition x geography level (the proxy map)")
    context = build_context()
    write_table(context, "geo_context", producer=PRODUCER,
                description="One wide row per state / county / ZCTA: ACS, SVI, RUCC, HPSA/MUA, PLACES context")
    features, pmeta = build_features(burden, context)
    write_table(features, "geo_condition_features", producer=PRODUCER,
                description="One row per (state|county, condition): burden (+level, resolution, inherited), "
                            "vulnerability, relevant providers, HRSA, HPSA, trials, NIH, diagnostic_desert + components")
    return burden, context, features, pmeta


def run(skip_build: bool = False) -> dict:
    if skip_build:
        burden, context, features = (read_table(t) for t in ("geo_condition_burden", "geo_context",
                                                             "geo_condition_features"))
        pmeta = {"meta": {}}
    else:
        burden, context, features, pmeta = build_tables()
    files = []
    cov = coverage_matrix(features, burden)
    files.append(_write_csv(cov, "geo_burden_coverage_matrix"))
    alt = primary_vs_alternates(burden)
    files.append(_write_csv(alt, "geo_primary_vs_alternates"))
    sens = desert_sensitivity(features)
    files.append(_write_csv(sens, "geo_desert_sensitivity"))
    files.append(_write_csv(missingness_table(features), "geo_feature_missingness"))
    files.append(_write_csv(trial_radius_qa(features), "geo_trial_radius_qa"))
    top_n = load_config("scoring")["ranking"]["top_n"]
    min_pop = load_config("scoring")["ranking"]["min_population"]
    dtop = (features[features["in_analysis_universe"] & features["diagnostic_desert"].notna()
                     & ((features["geo_level"] == "state") | (features["population_total"] >= min_pop))]
            # stable sort with explicit tie-breakers: exact desert ties (common at state level) must not be ordered
            # by the sort algorithm, or a tie at the top_n boundary could change membership between builds
            .sort_values(["diagnostic_desert", "condition_id", "geo_level", "geo_id"],
                         ascending=[False, True, True, True], kind="mergesort")
            .groupby(["condition_id", "geo_level"]).head(top_n))
    files.append(_write_csv(dtop[["condition_id", "geo_level", "geo_id", "geo_name", "object_id", "diagnostic_desert",
                                  "desert_pct_burden", "desert_pct_vulnerability", "desert_pct_low_providers",
                                  "desert_pct_low_trials", "burden_value", "burden_evidence_level", "burden_inherited",
                                  "svi_overall", "relevant_providers_per_100k", "trials_in_geo_or_within_50km_n",
                                  "trials_in_geo_n", "population_total"]], "geo_desert_top25"))
    if pmeta.get("meta"):
        pm = pd.DataFrame([{"condition_id": k, **{kk: (";".join(vv) if isinstance(vv, list) else vv)
                                                  for kk, vv in v.items()}} for k, v in pmeta["meta"].items()])
        files.append(_write_csv(pm, "geo_provider_group_counts"))
    t4 = run_test4(features, burden)
    for k, v in t4.items():
        files.append(_write_csv(v, k))
    t6 = run_test6(features, burden, context)
    t6["test6_geography_crh_calibration"] = crh_calibration(context)
    for k, v in t6.items():
        files.append(_write_csv(v, k))
    figs = draft_figures(cov, t4, t6, features)
    readings = _readings(t4, t6, sens, {"trial_radius_qa": trial_radius_qa(features), "alternates": alt,
                                        "features": features})
    tables = []
    for t, note in (("geo_condition_burden", "long harmonised burden table"),
                    ("geo_burden_definitions", "pre-specified primary measure per condition x level"),
                    ("geo_context", "wide context per state/county/ZCTA"),
                    ("geo_condition_features", "one row per state|county x condition")):
        meta = json.loads((config.PROCESSED / f"{t}.meta.json").read_text())
        tables.append({"table": f"data/processed/{t}.parquet", "rows": meta["rows"], "note": note})
    b = burden
    role_counts = b["measure_role"].value_counts().to_dict()
    tables.append({"table": "geo_condition_burden rows by measure_role", "rows": len(b),
                   "note": "; ".join(f"{k}: {v:,}" for k, v in role_counts.items())})
    tables.append({"table": "geo_condition_burden inherited rows (long_covid county)",
                   "rows": int(b["inherited"].sum()), "note": "source_geographic_resolution = state"})
    legacy = int(((b["geo_match_status"] == "matched_ct_legacy")).sum())
    tables.append({"table": "geo_condition_burden rows on legacy CT counties (not in features)", "rows": legacy,
                   "note": "CMS 2012-2023 + Lyme 2001-2022; not remapped to planning regions"})
    caveats = [
        "Burden rows are ecological estimates for places; none describes any person in the person-level datasets.",
        "Long COVID has no county measure; the county value is the state HPS estimate (inherited, source resolution "
        "'state'). HPS is an experimental, low-response online survey that ended 2024-09-16; state CIs are wide.",
        "CMS MMD measures cover Medicare fee-for-service beneficiaries only (mostly 65+ and disability-entitled) and "
        "reflect coding practice; CMS 51 is a broad pain/fatigue composite (identical values for fibromyalgia and "
        "ME/CFS); CT planning regions have no CMS value.",
        "Lyme counts are reported cases by county of residence (undercount; 2022+ case definition; excludes cases with "
        "unknown county).",
        "PLACES values are modeled small-area estimates that borrow strength from state BRFSS and county covariates; "
        "Kentucky and Pennsylvania lack the 2023-BRFSS measures.",
        "Provider counts are self-reported NPPES taxonomies at ZIP-centroid locations; a specialty does not mean a "
        "provider evaluates or treats the condition. Trial sites are city centroids, so 50/100 km radii are approximate.",
        "Research activity (trials, NIH) is not burden; it follows academic and population density.",
    ]
    stats_ = {"generated_at": utc_now_iso(), "tables": tables, "caveats": caveats,
              "files": files + figs + ["docs/BURDEN_DEFINITIONS.md", "docs/ANALYSIS_PLAN_GEOGRAPHY.md"], **readings}
    md = write_results_md(stats_)
    return {"files": files, "figures": figs, "results_md": str(md)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--skip-build", action="store_true", help="reuse the processed tables; rerun analyses only")
    args = ap.parse_args()
    out = run(skip_build=args.skip_build)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
