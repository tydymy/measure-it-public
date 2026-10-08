"""Deployment-opportunity scoring (SPEC Phase 4): components, composite, Monte Carlo rank intervals.

    build(...)                        -> writes data/processed/deployment_opportunities.parquet
    score_combo(condition, measurement, level, n_draws=...) -> one frame (same engine, any combination)

Pre-specified in docs/ANALYSIS_PLAN_SCORING.md; method reference docs/SCORING.md.

For geography g, condition c (or a condition set) and measurement m (class or bundle):

    burden(g,c)              primary measure of geo_condition_features (evidence level, inherited flag, resolution)
    vulnerability(g)         CDC/ATSDR SVI overall
    diagnostic_desert(g,c)   the geography module's index (access-only variant for level D, labelled)
    clinic_capacity(g,m)     implementer-group individual providers + implementer facilities per capita,
                             in the county and within the default radius (50 km)
    research_readiness(g,c,m) distinct condition trials, NIH core projects and technology-experience trials at
                             facilities in the county or within 50 km (state: in the state)
    technology_saturation(g,m) facilities with technology experience per capita (reported; not in the default
                             formula; one sensitivity weight set)

Each component is rescaled to [0, 1] by percentile rank within the 50 states + DC (configs/scoring.yaml) and
every raw and normalised column is kept. The composite is a weighted mean of the available normalised components
(level D burden is excluded and the other weights renormalised). Ranks are computed among ranking-eligible
geographies (counties with population >= scoring.yaml ranking.min_population; all states).

Burden completeness (2026-09-24; docs/ANALYSIS_PLAN_SCORING.md deviation 7): a condition's or condition set's
burden is computed only where EVERY member with a defined burden measure at the level has a value (a single condition
is a set of one). Regions where such a member has no value are 'incomplete': they keep their non-burden components,
get no burden, composite or rank, and are listed separately. Members with no defined burden (level D everywhere) are
left out of the burden; if no member has one, the burden weight is removed and the other weights renormalised.

Uncertainty: burden is perturbed within its CI inflated by the evidence-level multiplier (CMS: binomial SE at
the size-class midpoint + integer rounding; inherited county values: a shared state error, further inflated by
the inherited-burden multiplier of configs/scoring.yaml because a state value carries no within-state
information, plus a per-county deviation tau), and weights are drawn from a Dirichlet around the default; the
5th-95th percentile of the rank over the draws is the rank interval.

Everything here is ecological: places, facilities and registries, never people. A ranked region is a candidate
deployment opportunity, not a validated diagnostic pathway.
"""
from __future__ import annotations

import json
import re
import sys
import time
import warnings
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

from ..config import PROCESSED, SEED, UNKNOWN, load_config, utc_now_iso
from ..facilities import matching as M
from ..facilities.registry import SEP, haversine_km
from ..provenance import add_provenance
from ..store import read_table, write_table
from . import metric_link as ML

PRODUCER = "measure_it.scoring.opportunity"
TABLE = "deployment_opportunities"
COMPONENTS = ["burden", "vulnerability", "diagnostic_desert", "clinic_capacity", "research_readiness"]
PCT = {c: f"{c}_pct" for c in COMPONENTS}
# the metric -> translation link (docs/ANALYSIS_PLAN_METRIC_LINK.md, 2026-09-25): two absolute [0, 1] components used
# only by weight sets that name them (evidence_weighted); the five components above and every older set are unchanged
EVIDENCE_COMPONENTS = list(ML.EVIDENCE_COMPONENTS)
EVIDENCE_WEIGHT_SET = "evidence_weighted"
YIELD_BASES = ("count", "percentile", "level_D")
SINGLE_CONDITIONS = ["long_covid", "me_cfs", "pots", "dysautonomia"]
QUERY_SET_ID = "long_covid_or_me_cfs"
MEASUREMENTS = ["wearable_autonomic_activity_monitoring", "nailfold_capillaroscopy", "autonomic_function_testing",
                "exercise_capacity_testing"]
PRIMARY = {"condition": QUERY_SET_ID, "measurement": "wearable_autonomic_activity_monitoring", "level": "county"}
LEVELS = ("county", "state")
MC_DRAWS = 1000
MC_CHUNK = 200
DIRICHLET_CONCENTRATION = 20.0
LEVEL_ORDER = {"A": 0, "B": 1, "C": 2}
# complete = a burden only where every member with a defined burden has a value (since 2026-09-24);
# available = mean over the members with a value in the row (the rule before 2026-09-24; kept for the before/after)
BURDEN_RULES = ("complete", "available")
# all = relevant-provider density as built by the geography module / configs; specialists_only = sensitivity S6
PROVIDER_VARIANTS = ("all", "specialists_only")
SPECIALIST_EXCLUDED_GROUP = "primary_care"
# Medicare FFS size classes published by CMS MMD -> geometric midpoint (10,000+: the lower bound, conservative)
SIZE_CLASS_N = {"11-499": 74.0, "500-999": 707.0, "1,000-4,999": 2236.0, "5,000-9,999": 7071.0, "10,000+": 10000.0}
GUARDRAIL = ("A ranked region is a candidate deployment opportunity for pilot evaluation, not a validated diagnostic "
             "pathway. All joins are ecological (places, facilities, registries); nothing describes a person.")


# ------------------------------------------------------------------------------------------------------------------
# configuration helpers (pure)
# ------------------------------------------------------------------------------------------------------------------

def scoring_cfg() -> dict:
    return load_config("scoring")


def radius_km() -> float:
    return float(scoring_cfg().get("clinic_matching", {}).get("default_radius_km", 50))


def min_population() -> int:
    return int(scoring_cfg().get("ranking", {}).get("min_population", 10000))


def default_weights() -> dict:
    return dict(scoring_cfg()["default_weights"])


def weight_sets(include_extra: bool = True) -> dict[str, dict]:
    """Named weight sets from configs/scoring.yaml. include_extra=False drops the sets that use a component
    outside COMPONENTS (saturation_adjusted, which enters low_saturation = 1 - technology_saturation_pct)."""
    ws = {k: dict(v) for k, v in scoring_cfg()["weight_sets"].items()}
    if not include_extra:
        ws = {k: v for k, v in ws.items() if set(v) <= set(COMPONENTS)}
    return ws


def multiplier(level) -> float:
    m = scoring_cfg()["evidence_level_uncertainty_multiplier"].get(str(level))
    return float(m) if m is not None else np.nan


def inherited_multiplier() -> float:
    """Uncertainty multiplier for a state value inherited by a county (configs/scoring.yaml; 1.0 if absent)."""
    return float(scoring_cfg().get("inherited_burden_uncertainty_multiplier", 1.0))


def condition_sets() -> dict[str, dict]:
    dc = load_config("conditions")["demo_cluster"]
    return {
        dc["id"]: {"label": dc["label"], "members": list(dc["members"]), "set_kind": "demo_cluster",
                   "aliases": ["demo cluster", "autonomic activity invisible illness", dc["label"]]},
        QUERY_SET_ID: {"label": "Long COVID or ME/CFS", "members": ["long_covid", "me_cfs"],
                       "set_kind": "query_condition_set",
                       "aliases": ["Long COVID or ME/CFS", "Long COVID / ME/CFS", "long covid or mecfs",
                                   "Long COVID and ME/CFS", "long_covid+me_cfs"]},
    }


def condition_spec(condition_id: str) -> dict:
    """Canonical condition id or condition-set id -> {condition_id, kind, members, label}."""
    sets = condition_sets()
    if condition_id in sets:
        s = sets[condition_id]
        return {"condition_id": condition_id, "kind": "condition_set", "members": list(s["members"]),
                "label": s["label"], "set_kind": s["set_kind"]}
    conds = {c["id"]: c for c in load_config("conditions")["conditions"]}
    if condition_id in conds:
        return {"condition_id": condition_id, "kind": "condition", "members": [condition_id],
                "label": conds[condition_id]["preferred_name"], "set_kind": None}
    if "+" in condition_id:  # ad hoc set of canonical ids, e.g. "fibromyalgia+me_cfs"
        mem = sorted(set(condition_id.split("+")))
        if all(m in conds for m in mem):
            return {"condition_id": "+".join(mem), "kind": "condition_set", "members": mem,
                    "label": " or ".join(conds[m]["preferred_name"] for m in mem), "set_kind": "ad_hoc"}
    raise KeyError(f"unknown condition or condition set {condition_id!r}")


def measurement_spec(measurement_id: str) -> dict:
    """Measurement class/bundle -> members, implementer groups, adapter (the adapter swap goes through here)."""
    m = M.resolve_measurement(measurement_id)
    if m["status"] != "matched":
        raise KeyError(f"measurement {measurement_id!r}: {m.get('reason')}")
    out = dict(m)
    out["adapter_name"], out["adapter_status"], out["adapter_id"] = UNKNOWN, UNKNOWN, UNKNOWN
    try:
        from ..measurements.adapters import adapter_for_bundle, resolve_bundle
        if m["kind"] == "bundle":
            rb = resolve_bundle(m["measurement_id"])
            out["adapter_name"] = rb.get("adapter_name", UNKNOWN)
            a = adapter_for_bundle(m["measurement_id"])
            if a is not None:
                d = a.describe_measurement()
                out["adapter_status"], out["adapter_id"] = d.status, d.adapter_id
            else:
                out["adapter_status"] = "no adapter (adapter: null in configs/relevance.yaml)"
    except Exception as e:  # adapter registry is optional for scoring
        out["adapter_status"] = f"adapter lookup failed: {type(e).__name__}"
    return out


def percentile_rank(x) -> np.ndarray:
    """pr(x) = rank(average ties) / n over non-null values; NaN stays NaN (same as geography.features)."""
    s = pd.Series(np.asarray(x, dtype=float))
    return s.rank(pct=True, method="average").to_numpy()


def pct_rank_rows(X: np.ndarray) -> np.ndarray:
    """Row-wise percentile rank for a (draws x n) matrix whose NaN pattern is the same in every row."""
    X = np.atleast_2d(np.asarray(X, dtype=float))
    out = np.full(X.shape, np.nan)
    ok = ~np.isnan(X).any(axis=0)
    if ok.sum() == 0:
        return out
    out[:, ok] = rankdata(X[:, ok], axis=1, method="average") / ok.sum()
    return out


def weighted_composite(comps: dict[str, np.ndarray], weights: dict[str, float] | np.ndarray,
                       order: list[str] | None = None) -> np.ndarray:
    """Weighted mean over AVAILABLE components (NaN component -> its weight removed, the rest renormalised).

    comps: name -> array (n,) or (draws, n). weights: dict name -> w, or array (draws, K) in `order`."""
    order = order or list(comps)
    arrs = [np.atleast_2d(np.asarray(comps[k], dtype=float)) for k in order]
    D = max(a.shape[0] for a in arrs)
    X = np.stack([np.broadcast_to(a, (D, a.shape[1])) for a in arrs], axis=2)  # D x n x K
    if isinstance(weights, dict):
        W = np.array([[float(weights.get(k, 0.0)) for k in order]])
    else:
        W = np.asarray(weights, dtype=float)
    W = W[:, None, :]  # D x 1 x K
    avail = ~np.isnan(X)
    num = np.where(avail, X, 0.0) * W
    den = (avail * W).sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = num.sum(axis=2) / den
    c[den <= 0] = np.nan
    return c


def set_burden_member_mean(B: np.ndarray, defined, rule: str = "complete") -> np.ndarray:
    """Raw (set) burden from member burden percentiles B (members x draws x n).

    complete  (since 2026-09-24): mean over the members with a defined burden (`defined[k]`), NaN wherever one of
              them has no value (an incomplete region), all NaN when no member has a defined burden;
    available (the rule before 2026-09-24): mean over the members with a value in the row, so a row missing one
              member rests on the others alone. Kept only for the before/after comparison."""
    B = np.asarray(B, dtype=float)
    defined = np.asarray(defined, dtype=bool)
    if rule == "available":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmean(B, axis=0)
    if rule != "complete":
        raise ValueError(f"unknown burden rule {rule!r}; use one of {BURDEN_RULES}")
    if not defined.any():
        return np.full(B.shape[1:], np.nan)
    return B[defined].mean(axis=0)


def rank_rows(C: np.ndarray, eligible: np.ndarray) -> np.ndarray:
    """Rank 1 = highest composite among eligible columns; ties by column order (rows are sorted by geo_id)."""
    C = np.atleast_2d(np.asarray(C, dtype=float))
    R = np.full(C.shape, np.nan)
    idx = np.flatnonzero(eligible)
    if len(idx) == 0:
        return R
    sub = C[:, idx]
    key = np.where(np.isnan(sub), np.inf, -sub)
    order = np.argsort(key, axis=1, kind="stable")
    ranks = np.empty(order.shape, dtype=float)
    np.put_along_axis(ranks, order, np.broadcast_to(np.arange(1, len(idx) + 1, dtype=float), order.shape), axis=1)
    ranks[np.isnan(sub)] = np.nan
    R[:, idx] = ranks
    return R


def size_class_n(desc) -> float:
    """'Medicare FFS beneficiaries; size class 1,000-4,999' -> 2236 (geometric midpoint)."""
    if desc is None or (isinstance(desc, float) and np.isnan(desc)):
        return np.nan
    m = re.search(r"size class\s+([0-9,]+\s*-\s*[0-9,]+|[0-9,]+\+)", str(desc))
    if not m:
        return np.nan
    return SIZE_CLASS_N.get(m.group(1).replace(" ", ""), np.nan)


def cms_sd_pp(value_pct, n) -> np.ndarray:
    """Binomial SE (percentage points) at size-class n plus the integer-rounding variance 1/12."""
    p = np.clip(np.asarray(value_pct, dtype=float) / 100.0, 0, 1)
    n = np.asarray(n, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.sqrt(p * (1 - p) / n * 1e4 + 1.0 / 12.0)


def between_state_tau(values, ses) -> float:
    """Method-of-moments between-state SD: sqrt(max(0, var(values) - mean(se^2)))."""
    v = np.asarray(values, dtype=float)
    s = np.asarray(ses, dtype=float)
    ok = ~np.isnan(v) & ~np.isnan(s)
    if ok.sum() < 3:
        return np.nan
    return float(np.sqrt(max(0.0, np.var(v[ok], ddof=1) - np.mean(s[ok] ** 2))))


# ------------------------------------------------------------------------------------------------------------------
# spatial base: locations, facility and provider arrays, geography x location incidence
# ------------------------------------------------------------------------------------------------------------------

@dataclass
class LevelGeo:
    level: str
    geo: pd.DataFrame                 # universe rows sorted by geo_id
    loc_in: sparse.csr_matrix         # n_geo x n_loc (1 = location inside the geography)
    loc_50: sparse.csr_matrix | None  # n_geo x n_loc (1 = within radius of the internal point); None for states
    pool: sparse.csr_matrix           # in OR within radius (state: in)


@dataclass
class Spatial:
    loc: pd.DataFrame                 # lat, lon, county_fips, state
    fac_loc: np.ndarray               # facility index -> location id
    npi_loc: np.ndarray               # individual NPI index -> location id
    npi_bits: np.ndarray              # int64 bitmask of specialty groups per NPI
    group_bit: dict                   # group -> bit
    base: M.Base                      # facilities base shared with find_candidate_clinics
    levels: dict = field(default_factory=dict)
    input_versions: dict = field(default_factory=dict)


def _meta(name: str) -> dict:
    p = PROCESSED / f"{name}.meta.json"
    return json.loads(p.read_text()) if p.exists() else {}


def _universe(level: str) -> pd.DataFrame:
    cols = ["geo_id", "geo_level", "geo_name", "state_fips", "state_abbr", "lat", "lon", "population_total",
            "population_adults_18plus", "population_within_50km", "small_population_flag", "in_analysis_universe"]
    f = read_table("geo_condition_features", columns=cols + ["condition_id"])
    f = f[(f["condition_id"] == "long_covid") & (f["geo_level"] == level) & f["in_analysis_universe"]]
    g = f.drop(columns=["condition_id"]).sort_values("geo_id").reset_index(drop=True)
    if level == "county":
        g["eligible"] = g["population_total"].fillna(0) >= min_population()
    else:
        g["eligible"] = True
    return g


@lru_cache(maxsize=1)
def spatial() -> Spatial:
    t0 = time.time()
    b = M.load_base()
    fac = b.fac
    p = read_table("providers", columns=["npi", "entity_type", "state", "county_fips", "lat", "lon",
                                         "specialty_groups"])
    p = p[(p["entity_type"].astype(str) == "1") & p["lat"].notna() & p["specialty_groups"].fillna("").ne("")]
    p = p.reset_index(drop=True)
    groups = sorted({g for s in p["specialty_groups"].unique() for g in str(s).split(SEP) if g} | set(b.groups))
    group_bit = {g: i for i, g in enumerate(groups)}
    assert len(groups) < 63
    bits = np.zeros(len(p), dtype=np.int64)
    for s, idx in p.groupby("specialty_groups").indices.items():
        v = 0
        for g in str(s).split(SEP):
            if g:
                v |= 1 << group_bit[g]
        bits[idx] = v
    # unified location table
    fl = pd.DataFrame({"lat": fac["lat"].round(6).to_numpy(), "lon": fac["lon"].round(6).to_numpy(),
                       "county_fips": fac["county_fips"].astype(object).where(fac["county_fips"].notna(), "").to_numpy(),
                       "state": fac["state"].astype(object).where(fac["state"].notna(), "").to_numpy()})
    pl = pd.DataFrame({"lat": p["lat"].round(6).to_numpy(), "lon": p["lon"].round(6).to_numpy(),
                       "county_fips": p["county_fips"].astype(object).where(p["county_fips"].notna(), "").to_numpy(),
                       "state": p["state"].astype(object).where(p["state"].notna(), "").to_numpy()})
    allloc = pd.concat([fl, pl], ignore_index=True)
    codes, uniq = pd.factorize(pd.MultiIndex.from_frame(allloc))
    loc = uniq.to_frame(index=False)
    loc.columns = ["lat", "lon", "county_fips", "state"]
    fac_loc = codes[: len(fl)].astype(np.int64)
    npi_loc = codes[len(fl):].astype(np.int64)
    S = Spatial(loc=loc, fac_loc=fac_loc, npi_loc=npi_loc, npi_bits=bits, group_bit=group_bit, base=b)
    for level in LEVELS:
        S.levels[level] = _level_geo(level, loc)
    S.input_versions = {n: _meta(n).get("created_at", UNKNOWN) for n in
                        ("geo_condition_features", "geo_condition_burden", "facilities", "facility_trials",
                         "facility_nih_projects", "providers", "geographies")}
    # progress goes to stderr: a tool's stdout (ask --json, python -m measure_it.tools) must stay parseable
    print(f"[scoring] spatial base: {len(fac):,} facilities, {len(p):,} individual NPIs, {len(loc):,} locations "
          f"({time.time() - t0:.1f}s)", file=sys.stderr, flush=True)
    return S


def _incidence(rows: np.ndarray, cols: np.ndarray, shape) -> sparse.csr_matrix:
    m = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=shape)
    m.sum_duplicates()
    m.data[:] = 1.0
    return m


def _level_geo(level: str, loc: pd.DataFrame) -> LevelGeo:
    geo = _universe(level)
    n_geo, n_loc = len(geo), len(loc)
    if level == "county":
        gi = pd.Series(np.arange(n_geo), index=geo["geo_id"])
        r = loc["county_fips"].map(gi)
        ok = r.notna().to_numpy()
        loc_in = _incidence(r[ok].to_numpy(int), np.flatnonzero(ok), (n_geo, n_loc))
        from sklearn.neighbors import BallTree
        R = 6371.0088
        rad = radius_km()
        tree = BallTree(np.radians(loc[["lat", "lon"]].to_numpy(float)), metric="haversine")
        pts = np.radians(geo[["lat", "lon"]].to_numpy(float))
        hits = tree.query_radius(pts, r=(rad + 0.5) / R)
        rows, cols = [], []
        for i, h in enumerate(hits):
            if len(h):
                d = haversine_km(geo["lat"].iloc[i], geo["lon"].iloc[i], loc["lat"].to_numpy()[h],
                                 loc["lon"].to_numpy()[h])
                h = h[d <= rad]
                rows.append(np.full(len(h), i))
                cols.append(h)
        loc_50 = _incidence(np.concatenate(rows), np.concatenate(cols), (n_geo, n_loc))
        pool = ((loc_in + loc_50) > 0).astype(np.float32).tocsr()
    else:
        gi = pd.Series(np.arange(n_geo), index=geo["state_abbr"])
        r = loc["state"].map(gi)
        ok = r.notna().to_numpy()
        loc_in = _incidence(r[ok].to_numpy(int), np.flatnonzero(ok), (n_geo, n_loc))
        loc_50 = None
        pool = loc_in
    return LevelGeo(level=level, geo=geo, loc_in=loc_in, loc_50=loc_50, pool=pool)


# ------------------------------------------------------------------------------------------------------------------
# facility / provider measures (clinic_capacity, research_readiness, technology_saturation)
# ------------------------------------------------------------------------------------------------------------------

def _group_mask(S: Spatial, groups) -> tuple[np.ndarray, int]:
    gi = {g: i for i, g in enumerate(S.base.groups)}
    fcols = [gi[g] for g in groups if g in gi]
    fac_impl = S.base.G[:, fcols].any(axis=1) if fcols else np.zeros(len(S.base.fac), bool)
    bits = 0
    for g in groups:
        if g in S.group_bit:
            bits |= 1 << S.group_bit[g]
    return fac_impl, bits


@lru_cache(maxsize=64)
def _item_incidence(kind: str, key: tuple) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(facility index, item code, item labels) for condition trials / NIH core projects / technology trials."""
    b = spatial().base
    if kind == "condition_trials":
        d = b.ft[M._has(b.ft["condition_ids_literal"], list(key))][["fidx", "nct_id"]].drop_duplicates()
        item = d["nct_id"].astype(str)
    elif kind == "condition_nih":
        d = b.fn[M._has(b.fn["condition_ids_precision"], list(key))][["fidx", "core_project_num"]].drop_duplicates()
        item = d["core_project_num"].astype(str)
    elif kind == "tech_trials":
        d = b.ft[M._has(b.ft["measurement_classes"], list(key))][["fidx", "nct_id"]].drop_duplicates()
        item = d["nct_id"].astype(str)
    elif kind == "tech_trials_condition":
        cond, members = key
        ft = b.ft[M._has(b.ft["measurement_classes"], list(members))]
        d = ft[M._has(ft["condition_ids_literal"], list(cond))][["fidx", "nct_id"]].drop_duplicates()
        item = d["nct_id"].astype(str)
    else:
        raise ValueError(kind)
    codes, labels = pd.factorize(item)
    return d["fidx"].to_numpy(np.int64), codes.astype(np.int64), np.asarray(labels)


def _nih_appl_ids(cond_ids: tuple) -> pd.DataFrame:
    b = spatial().base
    return b.fn[M._has(b.fn["condition_ids_precision"], list(cond_ids))][["fidx", "appl_id", "core_project_num"]]


def _distinct_per_geo(pool: sparse.csr_matrix, fac_loc_eff: np.ndarray, fidx: np.ndarray, codes: np.ndarray,
                      n_items: int) -> tuple[np.ndarray, sparse.csr_matrix]:
    n_loc = pool.shape[1]
    if len(fidx) == 0:
        return np.zeros(pool.shape[0]), sparse.csr_matrix((pool.shape[0], max(n_items, 1)))
    LI = _incidence(fac_loc_eff[fidx], codes, (n_loc, max(n_items, 1)))
    R = (pool @ LI).tocsr()
    return np.diff(R.indptr).astype(float), R


def facility_measures(S: Spatial, level: str, cond_ids: tuple, meas: dict, fac_perm: np.ndarray | None = None,
                      npi_perm: np.ndarray | None = None, exclude_groups: tuple = (),
                      keep_incidence: bool = False) -> dict:
    """Raw counts and rates for clinic_capacity / research_readiness / technology_saturation.

    fac_perm / npi_perm: permutation (facility j occupies slot perm[j]) for the Test 6 location shuffle; the
    facility's attributes stay with it, its location tuple is the slot's. None = observed locations."""
    L = S.levels[level]
    geo = L.geo
    n_loc = len(S.loc)
    fac_loc = S.fac_loc if fac_perm is None else S.fac_loc[fac_perm]
    npi_loc = S.npi_loc if npi_perm is None else S.npi_loc[npi_perm]
    groups = [g for g in meas["implementer_groups"] if g not in set(exclude_groups)]
    fac_impl, bits = _group_mask(S, groups)
    npi_impl = (S.npi_bits & bits) != 0
    pop = geo["population_total"].to_numpy(float)
    pop = np.where(pop > 0, pop, np.nan)  # a rate needs a positive denominator (NaN otherwise, not inf)
    out = {}
    fl_counts = np.bincount(fac_loc[fac_impl], minlength=n_loc).astype(float)
    nl_counts = np.bincount(npi_loc[npi_impl], minlength=n_loc).astype(float)
    all_fac = np.bincount(fac_loc, minlength=n_loc).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["cc_impl_providers_in_geo_n"] = L.loc_in @ nl_counts
        out["cc_impl_providers_in_geo_per_100k"] = out["cc_impl_providers_in_geo_n"] / pop * 1e5
        out["cc_impl_facilities_in_geo_n"] = L.loc_in @ fl_counts
        out["cc_impl_facilities_in_geo_per_100k"] = out["cc_impl_facilities_in_geo_n"] / pop * 1e5
        if L.loc_50 is not None:
            pop50 = geo["population_within_50km"].to_numpy(float)
            pop50 = np.where(pop50 > 0, pop50, np.nan)
            out["cc_impl_providers_within_50km_n"] = L.loc_50 @ nl_counts
            out["cc_impl_providers_within_50km_per_100k"] = out["cc_impl_providers_within_50km_n"] / pop50 * 1e5
            out["cc_impl_facilities_within_50km_n"] = L.loc_50 @ fl_counts
            out["cc_impl_facilities_within_50km_per_100k"] = out["cc_impl_facilities_within_50km_n"] / pop50 * 1e5
        else:
            pop50 = pop
        out["pool_facilities_n"] = L.pool @ all_fac
        # research readiness (distinct items over the pool)
        inc = {}
        for name, kind, key in (("rr_condition_trials_pool_n", "condition_trials", cond_ids),
                                ("rr_condition_nih_core_pool_n", "condition_nih", cond_ids),
                                ("rr_tech_trials_pool_n", "tech_trials", tuple(meas["members"])),
                                ("rr_tech_trials_condition_pool_n", "tech_trials_condition",
                                 (cond_ids, tuple(meas["members"])))):
            fidx, codes, labels = _item_incidence(kind, key)
            out[name], R = _distinct_per_geo(L.pool, fac_loc, fidx, codes, len(labels))
            if keep_incidence:
                inc[name] = (R, labels)
        # sites: facilities with >= 1 condition trial; technology-experience facilities (saturation)
        fidx, _, _ = _item_incidence("condition_trials", cond_ids)
        site = np.zeros(len(S.fac_loc), bool)
        site[np.unique(fidx)] = True
        out["rr_condition_trial_sites_pool_n"] = L.pool @ np.bincount(fac_loc[site], minlength=n_loc).astype(float)
        fidx, _, _ = _item_incidence("tech_trials", tuple(meas["members"]))
        tech = np.zeros(len(S.fac_loc), bool)
        tech[np.unique(fidx)] = True
        out["ts_tech_facilities_pool_n"] = L.pool @ np.bincount(fac_loc[tech], minlength=n_loc).astype(float)
        out["technology_saturation"] = out["ts_tech_facilities_pool_n"] / pop50 * 1e5
    for k, v in list(out.items()):
        out[k] = np.asarray(v, dtype=float).ravel()
    if keep_incidence:
        out["_incidence"] = inc
    return out


def condition_core_groups(cid: str) -> list[str]:
    """Core specialty groups of a condition (configs/relevance.yaml condition_specialties; curated)."""
    rel = load_config("relevance")["condition_specialties"]
    if cid not in rel:
        raise KeyError(f"no condition_specialties entry for {cid!r}")
    return list(rel[cid]["core"])


def relevant_provider_counts(cid: str, level: str, exclude_groups: tuple = ()) -> tuple[np.ndarray, np.ndarray]:
    """Relevant-provider density of the diagnostic desert, counted here from `providers` (spatial base).

    Distinct individual NPIs (entity type 1) with >= 1 taxonomy in the condition's core groups minus exclude_groups,
    primary practice location in the geography (county: its county FIPS; state: its state); returns (count, per 100k
    ACS population; NaN for a non-positive population). With exclude_groups=('primary_care',) this is the geography
    module's `relevant_specialists_n` (sensitivity S6); without, its `relevant_providers_n`."""
    return provider_counts_for_groups([g for g in condition_core_groups(cid) if g not in set(exclude_groups)], level)


@lru_cache(maxsize=1)
def _npi_state_groups() -> pd.DataFrame:
    """Individual NPIs with a specialty group and a state (including the ~2% without a county geocode), exploded to
    one row per NPI x group; the geography module counts state providers on this basis (by state_fips)."""
    p = read_table("providers", columns=["npi", "entity_type", "state_fips", "specialty_groups"])
    p = p[(p["entity_type"].astype(str) == "1") & p["state_fips"].notna() & p["specialty_groups"].fillna("").ne("")]
    x = p[["npi", "state_fips"]].assign(group=p["specialty_groups"].astype(str).str.split(SEP)).explode("group")
    return x[x["group"].fillna("").ne("")].reset_index(drop=True)


def provider_counts_for_groups(groups: list[str], level: str) -> tuple[np.ndarray, np.ndarray]:
    """Distinct individual NPIs with >= 1 taxonomy in `groups`, located in each geography, and per 100k.

    County: primary practice location in the county (the spatial base, as clinic_capacity). State: the NPI's state
    FIPS, which also counts NPIs without a county geocode (the geography module's state basis)."""
    S = spatial()
    L = S.levels[level]
    if level == "state":
        x = _npi_state_groups()
        cnt = x[x["group"].isin(list(groups))].drop_duplicates("npi").groupby("state_fips").size()
        n = cnt.reindex(L.geo["state_fips"].astype(str)).fillna(0).to_numpy(float)
    else:
        bits = 0
        for g in groups:
            if g in S.group_bit:
                bits |= 1 << S.group_bit[g]
        sel = (S.npi_bits & bits) != 0
        n = np.asarray(L.loc_in @ np.bincount(S.npi_loc[sel], minlength=len(S.loc)).astype(float)).ravel()
    pop = L.geo["population_total"].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(pop > 0, n / pop * 1e5, np.nan)
    return n, rate


CC_RATE_COLS = {"county": ["cc_impl_providers_in_geo_per_100k", "cc_impl_providers_within_50km_per_100k",
                           "cc_impl_facilities_in_geo_per_100k", "cc_impl_facilities_within_50km_per_100k"],
                "state": ["cc_impl_providers_in_geo_per_100k", "cc_impl_facilities_in_geo_per_100k"]}
RR_COLS = ["rr_condition_trials_pool_n", "rr_condition_nih_core_pool_n", "rr_tech_trials_pool_n"]


def capacity_readiness(fm: dict, level: str) -> dict:
    """Sub-percentiles and components from facility_measures output."""
    out = {}
    cc = []
    for c in CC_RATE_COLS[level]:
        out[f"{c}_pct"] = percentile_rank(fm[c])
        cc.append(out[f"{c}_pct"])
    out["clinic_capacity"] = np.nanmean(np.vstack(cc), axis=0) if cc else np.full(len(fm["pool_facilities_n"]), np.nan)
    rr = []
    for c in RR_COLS:
        out[f"{c}_pct"] = percentile_rank(fm[c])
        rr.append(out[f"{c}_pct"])
    out["research_readiness"] = np.mean(np.vstack(rr), axis=0)
    out["clinic_capacity_pct"] = percentile_rank(out["clinic_capacity"])
    out["research_readiness_pct"] = percentile_rank(out["research_readiness"])
    out["technology_saturation_pct"] = percentile_rank(fm["technology_saturation"])
    return out


def capacity_basis_for(measurement_id: str) -> str:
    """'measurement_activity' when configs/scoring.yaml clinic_capacity_basis is 'activity_where_dedicated' and the
    measurement has a dedicated billing code (facilities.activity.capacity_basis), else 'implementer_density'."""
    if scoring_cfg().get("clinic_capacity_basis", "implementer_density") != "activity_where_dedicated":
        return "implementer_density"
    from ..facilities import activity as A
    try:
        return "measurement_activity" if A.capacity_basis(measurement_id) == "dedicated" else "implementer_density"
    except Exception as e:  # activity tables not built: keep density and say so
        print(f"[scoring] activity capacity unavailable for {measurement_id} ({e}); implementer density used")
        return "implementer_density"


def apply_capacity_basis(cr: dict, measurement_id: str, level: str, geo: pd.DataFrame) -> dict:
    """Replace clinic_capacity by measurement-specific clinician activity where the basis says so; the density
    version stays in clinic_capacity_density(_pct)."""
    cr["clinic_capacity_density"] = cr["clinic_capacity"]
    cr["clinic_capacity_density_pct"] = cr["clinic_capacity_pct"]
    basis = capacity_basis_for(measurement_id)
    cr["clinic_capacity_basis"] = basis
    if basis == "measurement_activity":
        from ..facilities import activity as A
        pct = A.capacity_pct_for(measurement_id, level, geo["geo_id"].to_numpy())
        cr["clinic_capacity"] = pct
        cr["clinic_capacity_pct"] = percentile_rank(pct)
    return cr


# ------------------------------------------------------------------------------------------------------------------
# burden / desert per member condition
# ------------------------------------------------------------------------------------------------------------------

FEATURE_COLS = ["geo_id", "geo_level", "condition_id", "burden_value", "burden_measure_id", "burden_measure_label",
                "burden_value_unit", "burden_evidence_level", "burden_defined_level", "burden_source_resolution",
                "burden_inherited", "burden_ci_low", "burden_ci_high", "burden_period", "burden_year",
                "burden_source_name", "burden_row_id", "burden_uncertainty_multiplier", "burden_unknown_reason",
                "burden_proxy_measure_id", "burden_proxy_value", "burden_proxy_ci_low", "burden_proxy_ci_high",
                "svi_overall", "desert_pct_burden", "desert_pct_vulnerability", "desert_pct_low_providers",
                "desert_pct_low_trials", "diagnostic_desert", "diagnostic_desert_access_only",
                "desert_burden_resolution", "desert_unknown_reason", "relevant_providers_per_100k",
                "trials_in_geo_or_within_50km_n", "trials_in_geo_n", "nih_core_projects_n", "object_id",
                "source_version", "retrieved_at", "in_analysis_universe"]


@lru_cache(maxsize=1)
def _features() -> pd.DataFrame:
    f = read_table("geo_condition_features", columns=FEATURE_COLS)
    return f[f["in_analysis_universe"]].copy()


@lru_cache(maxsize=1)
def _denominators() -> pd.Series:
    b = read_table("geo_condition_burden", columns=["burden_row_id", "denominator_description"])
    return b.drop_duplicates("burden_row_id").set_index("burden_row_id")["denominator_description"]


@dataclass
class Member:
    condition_id: str
    frame: pd.DataFrame     # aligned to LevelGeo.geo
    value: np.ndarray
    sd: np.ndarray
    sd_basis: np.ndarray
    inherited: np.ndarray
    state_idx: np.ndarray
    tau: float
    inh_mult: float         # inherited-burden multiplier already applied to sd[inherited]
    level_d: np.ndarray     # no usable burden for the row
    v: np.ndarray
    p: np.ndarray
    t: np.ndarray
    desert: np.ndarray       # burden-informed desert (NaN where not defined)
    desert_access: np.ndarray
    defined: bool = True     # the condition has a burden measure at this level (burden_defined_level != D)
    defined_level: str = "D"  # its defined evidence level (A/B/C; D = none)
    missing_reason: np.ndarray | None = None  # why a defined burden has no value in a row (None elsewhere)
    provider_variant: str = "all"


def member_data(cid: str, level: str, geo: pd.DataFrame, burden_override: str | None = None,
                provider_variant: str = "all") -> Member:
    """burden_override='proxy': use the PLACES C-proxy column (sensitivity S2 only).
    provider_variant='specialists_only': the desert's provider term counts the condition's core groups minus
    primary_care (sensitivity S6); 'all' = the geography module's term."""
    if provider_variant not in PROVIDER_VARIANTS:
        raise ValueError(f"unknown provider_variant {provider_variant!r}; use one of {PROVIDER_VARIANTS}")
    f = _features()
    f = f[(f["condition_id"] == cid) & (f["geo_level"] == level)].set_index("geo_id").reindex(geo["geo_id"])
    if f["condition_id"].isna().all():
        raise KeyError(f"no geo_condition_features rows for {cid} at {level}")
    f = f.reset_index()
    value = f["burden_value"].to_numpy(float)
    lev = f["burden_evidence_level"].astype(object).to_numpy()
    lo, hi = f["burden_ci_low"].to_numpy(float), f["burden_ci_high"].to_numpy(float)
    inherited = f["burden_inherited"].fillna(False).astype(bool).to_numpy()
    if burden_override == "proxy":
        value = f["burden_proxy_value"].to_numpy(float)
        lo, hi = f["burden_proxy_ci_low"].to_numpy(float), f["burden_proxy_ci_high"].to_numpy(float)
        lev = np.where(np.isnan(value), "D", "C").astype(object)
        inherited = np.zeros(len(value), bool)
    mult = np.array([multiplier(x) for x in lev])
    has_ci = ~np.isnan(lo) & ~np.isnan(hi)
    sd = np.where(has_ci, (hi - lo) / (2 * 1.959964), np.nan)
    basis = np.where(has_ci, "published 95% CI", "").astype(object)
    if burden_override != "proxy":
        n = f["burden_row_id"].map(_denominators()).map(size_class_n).to_numpy(float)
        cms = ~has_ci & ~np.isnan(n) & ~np.isnan(value)
        sd = np.where(cms, cms_sd_pp(value, n), sd)
        basis = np.where(cms, "binomial SE at CMS FFS size-class midpoint n=" + pd.Series(n).map(
            lambda x: f"{x:,.0f}" if x == x else "").to_numpy() + " + integer rounding", basis)
    sd = sd * mult
    level_d = np.isnan(value) | (lev == "D")
    inh_mult = inherited_multiplier()
    sd = np.where(inherited & ~level_d, sd * inh_mult, sd)
    basis = np.where(level_d, "level D: no burden", np.where(basis == "", "no uncertainty information", basis))
    basis = np.where(inherited & ~level_d, basis + f" of the state estimate (shared by the state's counties) x "
                     f"inherited-burden multiplier {inh_mult:g} + per-county deviation tau", basis)
    st = geo["state_fips"].astype(str).to_numpy()
    state_idx = pd.factorize(st)[0]
    tau = np.nan
    if inherited.any():
        sf = _features()
        s = sf[(sf["condition_id"] == cid) & (sf["geo_level"] == "state")]
        se = (s["burden_ci_high"] - s["burden_ci_low"]) / (2 * 1.959964)
        tau = between_state_tau(s["burden_value"], se)
    v = f["desert_pct_vulnerability"].to_numpy(float)
    p = f["desert_pct_low_providers"].to_numpy(float)
    t = f["desert_pct_low_trials"].to_numpy(float)
    desert = f["diagnostic_desert"].to_numpy(float)
    desert_access = f["diagnostic_desert_access_only"].to_numpy(float)
    if provider_variant == "specialists_only":
        _, rate = relevant_provider_counts(cid, level, exclude_groups=(SPECIALIST_EXCLUDED_GROUP,))
        p = 1.0 - percentile_rank(rate)
        desert_access = (v + p + t) / 3.0
    if burden_override == "proxy" or provider_variant == "specialists_only":
        desert = np.where(level_d, np.nan, (percentile_rank(value) + v + p + t) / 4)
    # a defined burden: the condition has a burden measure at this level (for the S2 proxy: any proxy value)
    if burden_override == "proxy":
        defined = bool((~np.isnan(value)).any())
        defined_level = "C" if defined else "D"
        reason = np.full(len(value), "no PLACES PHLTH proxy value for this geography", dtype=object)
    else:
        dl = f["burden_defined_level"].astype(object).where(f["burden_defined_level"].notna(), "D").astype(str)
        dls = [x for x in dl.unique() if x in LEVEL_ORDER]
        defined = bool(dls)
        defined_level = _least_direct(dls) if dls else "D"
        reason = f["burden_unknown_reason"].astype(object).where(
            f["burden_unknown_reason"].notna(),
            "no value for this geography in the primary measure (no source reason recorded)").to_numpy()
    missing_reason = np.where(level_d & defined, reason, None)
    return Member(condition_id=cid, frame=f, value=np.where(level_d, np.nan, value), sd=sd, sd_basis=basis,
                  inherited=inherited & ~level_d, state_idx=state_idx, tau=tau, inh_mult=inh_mult, level_d=level_d,
                  v=v, p=p, t=t, desert=np.where(level_d, np.nan, desert), desert_access=desert_access,
                  defined=defined, defined_level=defined_level, missing_reason=missing_reason,
                  provider_variant=provider_variant)


def burden_complete_mask(members: list, rule: str = "complete") -> np.ndarray:
    """True where every member with a defined burden has a value (rule 'complete'); all True for 'available'."""
    n = len(members[0].value)
    ok = np.ones(n, bool)
    if rule == "available":
        return ok
    if rule != "complete":
        raise ValueError(f"unknown burden rule {rule!r}; use one of {BURDEN_RULES}")
    for m in members:
        if m.defined:
            ok &= ~m.level_d
    return ok


# ------------------------------------------------------------------------------------------------------------------
# a scored combination and its evaluation (observed, Monte Carlo, permutations)
# ------------------------------------------------------------------------------------------------------------------

@dataclass
class Combo:
    cond: dict
    meas: dict
    level: str
    geo: pd.DataFrame
    members: list
    fm: dict                  # facility measures (raw)
    cr: dict                  # capacity/readiness percentiles and components
    vulnerability_pct: np.ndarray
    eligible: np.ndarray      # ranked rows: population-eligible AND complete burden
    pop_eligible: np.ndarray | None = None      # population >= min_population (county); all states
    burden_complete: np.ndarray | None = None   # every member with a defined burden has a value
    burden_rule: str = "complete"
    provider_variant: str = "all"
    exclude_groups: tuple = ()
    perf: dict | None = None                    # scored measurement_performance row (metric link)
    reach: pd.DataFrame | None = None           # reach, reach_population, ... aligned to geo
    yield_basis: str = "percentile"             # count / percentile / level_D (plan section 4.1)
    adults: np.ndarray | None = None            # ACS adults 18+


def performance_known(cb: "Combo") -> bool:
    return bool(cb.perf) and cb.perf.get("performance_status") == "known"


def yield_basis(members: list) -> str:
    """count: single condition, level A/B burden, not inherited, percent of all adults (a defensible count at this
    resolution); level_D: no member has a burden measure; percentile otherwise (docs/ANALYSIS_PLAN_METRIC_LINK.md 4.1)."""
    if not any(m.defined for m in members):
        return "level_D"
    if len(members) == 1:
        m = members[0]
        f = m.frame
        ok = ~m.level_d
        unit = f["burden_value_unit"].astype(str)[ok]
        label = f["burden_measure_label"].astype(str)[ok].str.lower()
        # a modelled small-area county value (BRFSS MRP, since 2026-10-07) is not counted as people at county level:
        # its within-state differences are not observed (results/BRFSS_SAE_RESULTS.md); state aggregates are
        modelled_county = (f["burden_source_name"].astype(str).str.contains("small-area")[ok].any()
                           and "county" in set(f["burden_source_resolution"].astype(str)[ok]))
        if (m.defined_level in ("A", "B") and not m.inherited.any() and ok.any() and (unit == "percent").all()
                and label.str.contains("of all adults").all() and not modelled_county):
            return "count"
    return "percentile"


def make_combo(condition_id: str, measurement_id: str, level: str, *, exclude_groups: tuple = (),
               burden_override: str | None = None, keep_incidence: bool = False, burden_rule: str = "complete",
               provider_variant: str = "all") -> Combo:
    """One scored combination. burden_rule 'complete' (default) ranks only regions where every member with a defined
    burden has a value; 'available' is the rule before 2026-09-24 (before/after comparison only).
    provider_variant 'specialists_only' + exclude_groups=('primary_care',) is sensitivity S6."""
    if burden_rule not in BURDEN_RULES:
        raise ValueError(f"unknown burden rule {burden_rule!r}; use one of {BURDEN_RULES}")
    S = spatial()
    cond = condition_spec(condition_id)
    meas = measurement_spec(measurement_id)
    L = S.levels[level]
    geo = L.geo
    # burden_override ('proxy', sensitivity S2) applies to long_covid only
    members = [member_data(c, level, geo, burden_override=burden_override if c == "long_covid" else None,
                           provider_variant=provider_variant)
               for c in cond["members"]]
    fm = facility_measures(S, level, tuple(cond["members"]), meas, exclude_groups=exclude_groups,
                           keep_incidence=keep_incidence)
    cr = capacity_readiness(fm, level)
    apply_capacity_basis(cr, meas["measurement_id"], level, geo)
    vul = percentile_rank(geo_svi(level, geo))
    pop_elig = geo["eligible"].to_numpy(bool)
    complete = burden_complete_mask(members, burden_rule)
    perf = ML.performance_for(cond["condition_id"], meas["measurement_id"])
    reach = ML.reach_for(level, [g for g in meas["implementer_groups"] if g not in set(exclude_groups)], geo)
    return Combo(cond=cond, meas=meas, level=level, geo=geo, members=members, fm=fm, cr=cr, vulnerability_pct=vul,
                 eligible=pop_elig & complete, pop_eligible=pop_elig, burden_complete=complete,
                 burden_rule=burden_rule, provider_variant=provider_variant, exclude_groups=tuple(exclude_groups),
                 perf=perf, reach=reach, yield_basis=yield_basis(members),
                 adults=geo["population_adults_18plus"].to_numpy(float))


def geo_svi(level: str, geo: pd.DataFrame) -> np.ndarray:
    f = _features()
    s = f[(f["condition_id"] == "long_covid") & (f["geo_level"] == level)].set_index("geo_id")["svi_overall"]
    return s.reindex(geo["geo_id"]).to_numpy(float)


def evaluate(cb: Combo, *, member_values: list | None = None, vulnerability_pct: np.ndarray | None = None,
             capacity_pct: np.ndarray | None = None, readiness_pct: np.ndarray | None = None,
             saturation_pct: np.ndarray | None = None, weights=None, order: list[str] | None = None,
             perf_draws: dict | None = None) -> dict:
    """Burden -> burden_pct (members, set mean) -> desert -> desert_pct -> composite, for 1 or many draws.

    member_values: per member (n,) or (draws, n) burden values (default observed). Other overrides replace a
    normalised component (Test 6). weights: dict or (draws, K) array over `order` (default equal weights)."""
    n = len(cb.geo)
    vals = member_values if member_values is not None else [m.value for m in cb.members]
    v = cb.vulnerability_pct if vulnerability_pct is None else vulnerability_pct
    bps, deserts, access = [], [], []
    for m, x in zip(cb.members, vals):
        x2 = np.atleast_2d(np.asarray(x, dtype=float))
        bp = pct_rank_rows(x2)
        bps.append(bp)
        vv = m.v if vulnerability_pct is None else v
        d = (bp + vv + m.p + m.t) / 4.0
        deserts.append(d)
        access.append(np.atleast_2d((vv + m.p + m.t) / 3.0))
    D = max(b.shape[0] for b in bps)
    B = np.stack([np.broadcast_to(b, (D, n)) for b in bps])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        # Burden completeness (deviation 7, 2026-09-24): mean over the members with a defined burden, NaN where one
        # of them has no value (incomplete region: excluded from the ranking); 'available' = the old rule.
        burden_member_mean = set_burden_member_mean(B, [m.defined for m in cb.members], cb.burden_rule)
        # A condition set's burden is an average of member percentiles; like every other averaged component
        # (desert, clinic_capacity, research_readiness) it is re-normalised with pr() before weighting
        # (plan section 3; configs/scoring.yaml). Without this the set burden has a compressed spread (SD ~0.21
        # vs ~0.28 for the other components) and silently carries less than its nominal weight. For a single
        # condition pr(pr(x)) == pr(x), so nothing changes there. (Reviewer correction 2026-09-23.)
        burden_pct = burden_member_mean if len(cb.members) == 1 else pct_rank_rows(burden_member_mean)
        DS = np.stack([np.broadcast_to(d, (D, n)) for d in deserts])
        desert = np.nanmean(DS, axis=0)
        AC = np.stack([np.broadcast_to(a, (D, n)) for a in access])
        desert_access = np.nanmean(AC, axis=0)
    no_informed = np.isnan(desert)
    desert = np.where(no_informed, desert_access, desert)
    desert_pct = pct_rank_rows(desert)
    comps = {"burden": burden_pct, "vulnerability": v,
             "diagnostic_desert": desert_pct,
             "clinic_capacity": cb.cr["clinic_capacity_pct"] if capacity_pct is None else capacity_pct,
             "research_readiness": cb.cr["research_readiness_pct"] if readiness_pct is None else readiness_pct}
    sat = cb.cr["technology_saturation_pct"] if saturation_pct is None else saturation_pct
    comps["low_saturation"] = 1.0 - np.asarray(sat, dtype=float)
    ev_comps, pB = evidence_components(cb, burden_pct, vals, perf_draws)
    comps.update(ev_comps)
    order = order or COMPONENTS
    w = weights if weights is not None else default_weights()
    comp = mask_incomplete(cb, weighted_composite({k: comps[k] for k in order}, w, order))
    return {"burden_pct": burden_pct, "burden_member_mean": burden_member_mean,
            "diagnostic_desert": desert, "desert_access_only_used": no_informed,
            "diagnostic_desert_pct": desert_pct, "components": comps, "composite": comp, "yield_burden_pct": pB}


def burden_count(cb: Combo, values=None) -> np.ndarray | None:
    """Count path only: burden value (% of all adults) / 100 x ACS adults 18+ (NaN elsewhere); None otherwise."""
    if cb.yield_basis != "count":
        return None
    v = cb.members[0].value if values is None else values
    return np.atleast_2d(np.asarray(v, dtype=float)) / 100.0 * cb.adults[None, :]


def evidence_components(cb: Combo, burden_pct: np.ndarray, member_values=None,
                        perf_draws: dict | None = None) -> tuple[dict, np.ndarray | None]:
    """measurement_evidence (constant, absolute) and expected_yield = pr(B) x sensitivity x reach (absolute, not
    re-normalised: docs/ANALYSIS_PLAN_METRIC_LINK.md sections 2.6 and 4.3). NaN when performance is not known (the
    evidence-weighted composite is then masked entirely) or when there is no burden (level D)."""
    n = len(cb.geo)
    perf = cb.perf or {}
    me = perf.get("measurement_evidence")
    me = float(me) if me is not None and pd.notna(me) else np.nan
    out = {"measurement_evidence": np.full(n, me if performance_known(cb) else np.nan)}
    if not performance_known(cb) or cb.yield_basis == "level_D" or cb.reach is None:
        out["expected_yield"] = np.full(n, np.nan)
        return out, None
    sens = float(perf["op_sensitivity"]) if perf_draws is None else np.asarray(perf_draws["sens"], float)[:, None]
    if cb.yield_basis == "count":
        vals = cb.members[0].value if member_values is None else member_values[0]
        pB = pct_rank_rows(burden_count(cb, vals))
    else:
        pB = np.atleast_2d(burden_pct)
    out["expected_yield"] = pB * sens * cb.reach["reach"].to_numpy(float)[None, :]
    return out, pB


def mask_evidence(cb: Combo, C: np.ndarray, comps: dict) -> np.ndarray:
    """UNKNOWN rule (plan section 5): no evidence-weighted composite when performance is not known; with known
    performance, a region without an expected yield (reach or burden missing) gets none either, except that a level-D
    condition (no burden at all) drops the yield weight like the burden weight."""
    C = np.array(np.atleast_2d(C), dtype=float, copy=True)
    if not performance_known(cb):
        return np.full(C.shape, np.nan)
    if cb.yield_basis != "level_D":
        y = np.broadcast_to(np.atleast_2d(np.asarray(comps["expected_yield"], dtype=float)), C.shape)
        C[np.isnan(y)] = np.nan
    return C


def uses_evidence(weights: dict) -> bool:
    return any(float(weights.get(k, 0) or 0) > 0 for k in EVIDENCE_COMPONENTS)


def mask_incomplete(cb: Combo, C: np.ndarray) -> np.ndarray:
    """Composite -> NaN in incomplete-burden regions (no composite and no rank there; components are kept)."""
    if cb.burden_complete is None or cb.burden_complete.all():
        return C
    C = np.array(C, dtype=float, copy=True)
    C[..., ~cb.burden_complete] = np.nan
    return C


def named_composites(cb: Combo, ev: dict | None = None) -> dict[str, np.ndarray]:
    ev = ev or evaluate(cb)
    out = {}
    for name, w in weight_sets().items():
        order = [k for k in w]
        c = mask_incomplete(cb, weighted_composite({k: ev["components"][k] for k in order}, w, order))
        if uses_evidence(w):
            c = mask_evidence(cb, c, ev["components"])
        out[name] = c[0]
    return out


def draw_burden(cb: Combo, rng: np.random.Generator, D: int, tau_scale: float = 1.0,
                inherited_multiplier: float | None = None) -> list[np.ndarray]:
    """Burden draws per member. inherited_multiplier overrides the configured inherited-burden multiplier
    (sensitivity S1b); None uses the configured value already in Member.sd."""
    out = []
    n = len(cb.geo)
    for m in cb.members:
        sd = np.nan_to_num(m.sd, nan=0.0)
        e = rng.standard_normal((D, n)) * sd
        if m.inherited.any():
            n_states = int(m.state_idx.max()) + 1
            es = rng.standard_normal((D, n_states))
            inh = m.inherited
            sd_inh = sd[inh] if inherited_multiplier is None else sd[inh] / m.inh_mult * inherited_multiplier
            e[:, inh] = es[:, m.state_idx[inh]] * sd_inh
            tau = 0.0 if np.isnan(m.tau) else m.tau * tau_scale
            e[:, inh] += rng.standard_normal((D, int(inh.sum()))) * tau
        x = np.maximum(m.value[None, :] + e, 0.0)
        x[:, np.isnan(m.value)] = np.nan
        out.append(x)
    return out


def monte_carlo(cb: Combo, n_draws: int = MC_DRAWS, seed: int = SEED, burden_noise: bool = True,
                weight_noise: bool = True, tau_scale: float = 1.0,
                inherited_multiplier: float | None = None) -> np.ndarray:
    """Rank matrix (draws x n) under burden and/or weight uncertainty; NaN for non-eligible rows."""
    rng = np.random.default_rng(seed)
    w0 = default_weights()
    alpha = np.array([w0[k] for k in COMPONENTS]) * DIRICHLET_CONCENTRATION
    ranks = np.full((n_draws, len(cb.geo)), np.nan, dtype=np.float32)
    done = 0
    while done < n_draws:
        D = min(MC_CHUNK, n_draws - done)
        vals = draw_burden(cb, rng, D, tau_scale, inherited_multiplier) if burden_noise else None
        W = rng.dirichlet(alpha, size=D) if weight_noise else np.array([[w0[k] for k in COMPONENTS]])
        ev = evaluate(cb, member_values=vals, weights=W, order=COMPONENTS)
        C = ev["composite"]
        if C.shape[0] == 1 and D > 1:
            C = np.repeat(C, D, axis=0)
        ranks[done:done + D] = rank_rows(C, cb.eligible)
        done += D
    return ranks


def monte_carlo_evidence(cb: Combo, n_draws: int = MC_DRAWS, seed: int = SEED + 11) -> dict | None:
    """Evidence-weighted Monte Carlo (docs/ANALYSIS_PLAN_METRIC_LINK.md 4.4): burden draws as in monte_carlo(),
    operating-point sensitivity / specificity on the logit scale from their 95% CIs, weights ~ Dirichlet(20 x
    evidence_weighted). Returns the rank matrix and the drawn yields / false positives; None when performance is not
    known (the combination is not ranked under evidence_weighted)."""
    ws = weight_sets().get(EVIDENCE_WEIGHT_SET)
    if ws is None or not performance_known(cb):
        return None
    perf = cb.perf
    rng = np.random.default_rng(seed)
    order = [k for k in ws]
    alpha = np.array([float(ws[k]) for k in order]) * DIRICHLET_CONCENTRATION
    n = len(cb.geo)
    reach = cb.reach["reach"].to_numpy(float)[None, :]
    adults = cb.adults[None, :]
    ranks = np.full((n_draws, n), np.nan, dtype=np.float32)
    tp = np.full((n_draws, n), np.nan, dtype=np.float32)
    fp = np.full((n_draws, n), np.nan, dtype=np.float32)
    done = 0
    while done < n_draws:
        D = min(MC_CHUNK, n_draws - done)
        vals = draw_burden(cb, rng, D)
        sens = ML.draw_logit(float(perf["op_sensitivity"]), float(perf["op_sensitivity_ci_low"]),
                             float(perf["op_sensitivity_ci_high"]), rng, D)
        spec = ML.draw_logit(float(perf["op_specificity"]), _f(perf.get("op_specificity_ci_low")),
                             _f(perf.get("op_specificity_ci_high")), rng, D)
        W = rng.dirichlet(alpha, size=D)
        ev = evaluate(cb, member_values=vals, weights=W, order=order, perf_draws={"sens": sens, "spec": spec})
        C = mask_evidence(cb, ev["composite"], ev["components"])
        ranks[done:done + D] = rank_rows(C, cb.eligible)
        if cb.yield_basis == "count":
            cnt = burden_count(cb, vals[0])
            tp[done:done + D] = cnt * sens[:, None] * reach
            fp[done:done + D] = (1 - spec)[:, None] * (adults - cnt) * reach
        elif cb.yield_basis == "percentile":
            tp[done:done + D] = np.atleast_2d(ev["burden_pct"]) * sens[:, None] * reach
            fp[done:done + D] = (1 - spec)[:, None] * adults * reach
        else:
            fp[done:done + D] = (1 - spec)[:, None] * adults * reach
        done += D
    return {"ranks": ranks, "yield": tp, "false_positives": fp}


def _f(v) -> float:
    return float(v) if v is not None and pd.notna(v) else np.nan


def summarize_ranks(R: np.ndarray, prefix: str) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q = np.nanpercentile(R, [5, 50, 95], axis=0)
        top10 = np.nanmean(np.where(np.isnan(R), np.nan, R <= 10), axis=0)
        top25 = np.nanmean(np.where(np.isnan(R), np.nan, R <= 25), axis=0)
    return {f"{prefix}_p05": q[0], f"{prefix}_p50": q[1], f"{prefix}_p95": q[2],
            f"p_top10_{prefix.replace('rank_', '')}": top10, f"p_top25_{prefix.replace('rank_', '')}": top25}


# ------------------------------------------------------------------------------------------------------------------
# frame assembly
# ------------------------------------------------------------------------------------------------------------------

def _least_direct(levels) -> str:
    ls = [l for l in levels if l in LEVEL_ORDER]
    return max(ls, key=LEVEL_ORDER.get) if ls else "D"


def combo_frame(cb: Combo, n_draws: int = MC_DRAWS, seed: int = SEED, mc: bool = True) -> pd.DataFrame:
    geo = cb.geo
    n = len(geo)
    ev = evaluate(cb)
    cond, meas = cb.cond, cb.meas
    df = pd.DataFrame({
        "object_id": [f"opportunity:{cond['condition_id']}|{meas['measurement_id']}|{g}" for g in geo["geo_id"]],
        "condition_id": cond["condition_id"], "condition_kind": cond["kind"],
        "condition_label": cond["label"], "condition_members": SEP.join(cond["members"]),
        "measurement_id": meas["measurement_id"], "measurement_label": meas["label"],
        "measurement_kind": meas["kind"], "measurement_members": SEP.join(meas["members"]),
        "implementer_groups": SEP.join(meas["implementer_groups"]),
        "measurement_adapter": meas.get("adapter_name", UNKNOWN), "measurement_adapter_id": meas.get("adapter_id", UNKNOWN),
        "measurement_adapter_status": meas.get("adapter_status", UNKNOWN),
        "geo_level": cb.level, "geo_id": geo["geo_id"].to_numpy(), "geo_object_id": "geo:" + geo["geo_id"],
        "geo_name": geo["geo_name"].to_numpy(), "state_fips": geo["state_fips"].to_numpy(),
        "state_abbr": geo["state_abbr"].to_numpy(), "population_total": geo["population_total"].to_numpy(),
        "population_within_50km": geo["population_within_50km"].to_numpy() if cb.level == "county" else np.nan,
        "small_population_flag": ~geo["eligible"].to_numpy(bool) if cb.level == "county" else False,
        "rank_eligible": cb.eligible,
    })
    # ---- burden
    single = len(cb.members) == 1
    complete = cb.burden_complete if cb.burden_complete is not None else np.ones(n, bool)
    has_val = np.array([[not m.level_d[i] for m in cb.members] for i in range(n)]).reshape(n, len(cb.members))
    defined = np.array([m.defined for m in cb.members])
    # contributing = members whose percentile enters the (set) burden of the row; none in an incomplete row
    contrib = has_val & defined[None, :] & complete[:, None]
    missing = ~has_val & defined[None, :]          # a defined member without a value in the row
    if single:
        m = cb.members[0]
        f = m.frame
        df["burden_value"] = m.value
        for c in ("burden_measure_id", "burden_measure_label", "burden_value_unit", "burden_period", "burden_year",
                  "burden_source_name", "burden_row_id", "burden_ci_low", "burden_ci_high"):
            df[c] = f[c].to_numpy()
        df["burden_evidence_level"] = np.where(m.level_d, "D", f["burden_evidence_level"].astype(object).to_numpy())
        df["burden_inherited"] = m.inherited
        df["burden_source_resolution"] = np.where(m.level_d, UNKNOWN, f["burden_source_resolution"].astype(object))
        df["burden_uncertainty_multiplier"] = [multiplier(x) for x in df["burden_evidence_level"]]
        df["burden_inherited_uncertainty_multiplier"] = np.where(m.inherited, m.inh_mult, 1.0)
        df["burden_mc_sd"] = np.where(m.level_d, np.nan, m.sd)
        df["burden_mc_sd_basis"] = m.sd_basis
        df["burden_unknown_reason"] = np.where(m.level_d, f["burden_unknown_reason"].fillna(
            "level D: no usable burden estimate").astype(object), None)
    else:
        df["burden_value"] = np.nan
        for c in ("burden_measure_id", "burden_measure_label", "burden_value_unit", "burden_period", "burden_year",
                  "burden_source_name", "burden_row_id", "burden_ci_low", "burden_ci_high"):
            df[c] = None if c not in ("burden_ci_low", "burden_ci_high", "burden_year") else np.nan
        levs = [m.frame["burden_evidence_level"].astype(object).to_numpy() for m in cb.members]
        df["burden_evidence_level"] = [_least_direct([levs[k][i] for k in range(len(cb.members)) if contrib[i, k]])
                                       for i in range(n)]
        df["burden_inherited"] = np.array([any(m.inherited[i] for m in cb.members) for i in range(n)])
        df["burden_source_resolution"] = [SEP.join(sorted({str(m.frame["burden_source_resolution"].iloc[i])
                                                            for k, m in enumerate(cb.members) if contrib[i, k]}))
                                          or UNKNOWN for i in range(n)]
        df["burden_uncertainty_multiplier"] = [multiplier(x) for x in df["burden_evidence_level"]]
        df["burden_inherited_uncertainty_multiplier"] = np.where(df["burden_inherited"], inherited_multiplier(), 1.0)
        df["burden_mc_sd"] = np.nan
        df["burden_mc_sd_basis"] = "per member (see member_components)"
        df["burden_unknown_reason"] = np.where(
            ~complete, "incomplete set burden (see burden_incomplete_reason)",
            np.where(~contrib.any(axis=1), "no member has a usable burden estimate", None))
        df["burden_value_note"] = ("condition set: raw set burden = burden_set_member_mean_pct, the mean of the "
                                   "burden percentiles of the members with a defined (non-D) burden measure at this "
                                   "level, computed only where every such member has a value; burden_pct = "
                                   "pr(burden_set_member_mean_pct) over those complete regions (re-normalised like "
                                   "every averaged component); members' raw values in member_components. Regions "
                                   "where a defined member has no value are incomplete (burden_incomplete): no set "
                                   "burden, composite or rank")
    df["burden_members_contributing"] = [SEP.join(m.condition_id for k, m in enumerate(cb.members) if contrib[i, k])
                                         for i in range(n)]
    # members with no burden measure at this level (level D everywhere): never part of the burden
    df["burden_members_excluded_level_D"] = SEP.join(m.condition_id for m in cb.members if not m.defined)
    # members with a defined burden but no value in this region (-> incomplete burden)
    df["burden_members_missing"] = [SEP.join(m.condition_id for k, m in enumerate(cb.members) if missing[i, k])
                                    for i in range(n)]
    df["burden_incomplete"] = ~complete
    df["burden_incomplete_reason"] = [incomplete_reason(cb, i) if not complete[i] else None for i in range(n)]
    df["burden_missing_source_reason"] = [
        "; ".join(f"{m.condition_id}: {m.missing_reason[i]}" for k, m in enumerate(cb.members) if missing[i, k]) or None
        for i in range(n)]
    df["burden_completeness_rule"] = (
        "complete: burden only where every member with a defined burden has a value (docs/ANALYSIS_PLAN_SCORING.md "
        "deviation 7)" if cb.burden_rule == "complete" else "available: mean over members with a value (before "
                                                             "2026-09-24)")
    # burden weight removed and the other weights renormalised: no member has a defined burden (level D)
    df["burden_excluded_level_D"] = ~contrib.any(axis=1) & complete
    df["burden_weight_renormalised"] = df["burden_excluded_level_D"]
    df["burden_set_member_mean_pct"] = np.nan if single else ev["burden_member_mean"][0]
    df["burden_pct"] = ev["burden_pct"][0]
    # ---- vulnerability
    df["svi_overall"] = geo_svi(cb.level, geo)
    df["vulnerability_pct"] = cb.vulnerability_pct
    # ---- desert
    if single:
        f = cb.members[0].frame
        for c in ("desert_pct_burden", "desert_pct_vulnerability", "desert_pct_low_providers", "desert_pct_low_trials",
                  "diagnostic_desert_access_only", "desert_burden_resolution", "relevant_providers_per_100k",
                  "trials_in_geo_or_within_50km_n" if cb.level == "county" else "trials_in_geo_n"):
            df[c] = f[c].to_numpy()
    df["diagnostic_desert"] = ev["diagnostic_desert"][0]
    df["diagnostic_desert_variant"] = np.where(
        ev["desert_access_only_used"][0], "access_only (level D; not burden-informed)",
        "burden_informed" if single else "mean of members' burden-informed deserts")
    df["diagnostic_desert_pct"] = ev["diagnostic_desert_pct"][0]
    # ---- clinic capacity, research readiness, saturation (raw + sub-percentiles + components)
    for k, v in cb.fm.items():
        if not k.startswith("_"):
            df[k] = v
    for k, v in cb.cr.items():
        df[k] = v
    # ---- composites and ranks
    comps = named_composites(cb, ev)
    for name, c in comps.items():
        df[f"composite_{name}"] = c
        df[f"rank_{name}"] = rank_rows(c, cb.eligible)[0]
    df["composite"] = df["composite_equal"]
    df["rank"] = df["rank_equal"]
    df["rank_equal_incl_small"] = rank_rows(comps["equal"], np.ones(n, bool))[0]
    df["weights_default"] = json.dumps(default_weights())
    # ---- member components (sets)
    if not single:
        mc = []
        for i in range(n):
            d = {}
            for m in cb.members:
                f = m.frame
                d[m.condition_id] = {
                    "burden_value": None if np.isnan(m.value[i]) else float(m.value[i]),
                    "burden_evidence_level": "D" if m.level_d[i] else str(f["burden_evidence_level"].iloc[i]),
                    "burden_defined_at_level": bool(m.defined),
                    "contributes_to_set_burden": bool(contrib[i, cb.members.index(m)]),
                    "burden_inherited": bool(m.inherited[i]),
                    # a level-D member has no measure: UNKNOWN / None, never the string 'nan'
                    "burden_source_resolution": (UNKNOWN if pd.isna(f["burden_source_resolution"].iloc[i])
                                                 else str(f["burden_source_resolution"].iloc[i])),
                    "burden_ci": [None if pd.isna(f["burden_ci_low"].iloc[i]) else float(f["burden_ci_low"].iloc[i]),
                                  None if pd.isna(f["burden_ci_high"].iloc[i]) else float(f["burden_ci_high"].iloc[i])],
                    "burden_measure_id": (None if pd.isna(f["burden_measure_id"].iloc[i])
                                          else str(f["burden_measure_id"].iloc[i])),
                    "burden_source_name": (None if pd.isna(f["burden_source_name"].iloc[i])
                                           else str(f["burden_source_name"].iloc[i])),
                    "burden_mc_sd": None if np.isnan(m.sd[i]) else round(float(m.sd[i]), 4),
                    "diagnostic_desert": None if np.isnan(m.desert[i]) else round(float(m.desert[i]), 5),
                    "diagnostic_desert_access_only": None if np.isnan(m.desert_access[i]) else round(float(m.desert_access[i]), 5),
                    "desert_burden_resolution": None if pd.isna(f["desert_burden_resolution"].iloc[i]) else str(f["desert_burden_resolution"].iloc[i]),
                }
            mc.append(json.dumps(d))
        df["member_components"] = mc
    else:
        df["member_components"] = None
    # ---- Monte Carlo
    if mc:
        R = monte_carlo(cb, n_draws, seed)
        for k, v in summarize_ranks(R, "rank_mc").items():
            df[k] = v
        Rb = monte_carlo(cb, n_draws, seed + 1, burden_noise=True, weight_noise=False)
        s = summarize_ranks(Rb, "rank_mc_burden")
        df["rank_mc_burden_p05"], df["rank_mc_burden_p95"] = s["rank_mc_burden_p05"], s["rank_mc_burden_p95"]
        Rw = monte_carlo(cb, n_draws, seed + 2, burden_noise=False, weight_noise=True)
        s = summarize_ranks(Rw, "rank_mc_weights")
        df["rank_mc_weights_p05"], df["rank_mc_weights_p95"] = s["rank_mc_weights_p05"], s["rank_mc_weights_p95"]
        df["mc_n_draws"] = n_draws
        df["mc_seed"] = seed
        df["mc_dirichlet_concentration"] = DIRICHLET_CONCENTRATION
        df["mc_tau_inherited"] = next((m.tau for m in cb.members if m.inherited.any()), np.nan)
        df["mc_inherited_multiplier"] = inherited_multiplier()
    df = pd.concat([df, metric_link_columns(cb, df, ev, n_draws=n_draws if mc else 0)], axis=1).copy()
    df["uncertainties"] = [json.dumps(u) for u in row_uncertainties(cb, df)]
    return df


PERF_COLS = {  # deployment_opportunities column -> measurement_performance field
    "measurement_performance_status": "performance_status", "measurement_performance_tier": "quality_tier",
    "measurement_performance_tier_label": "quality_tier_label",
    "measurement_performance_record_id": "selected_record_id", "measurement_performance_source_kind": "source_kind",
    "measurement_performance_dataset": "dataset_id", "measurement_performance_label_basis": "label_basis",
    "measurement_performance_comparator": "comparator", "measurement_performance_comparator_kind": "comparator_kind",
    "measurement_performance_auroc": "auroc", "measurement_performance_auroc_ci_low": "auroc_ci_low",
    "measurement_performance_auroc_ci_high": "auroc_ci_high", "measurement_performance_n_cases": "n_cases",
    "measurement_performance_n_controls": "n_controls", "op_sensitivity": "op_sensitivity",
    "op_sensitivity_ci_low": "op_sensitivity_ci_low", "op_sensitivity_ci_high": "op_sensitivity_ci_high",
    "op_specificity": "op_specificity", "op_specificity_ci_low": "op_specificity_ci_low",
    "op_specificity_ci_high": "op_specificity_ci_high", "measurement_evidence_undiscounted":
    "measurement_evidence_undiscounted", "measurement_evidence_conflict": "evidence_conflict",
    "measurement_evidence_conflict_note": "evidence_conflict_note", "measurement_performance_note": "performance_note",
}
YIELD_BASIS_NOTE = {
    "count": ("count path: burden value (% of all adults) x ACS adults 18+ is a defensible count at this resolution; "
              "expected_detectable_cases = count x sensitivity x reach"),
    "percentile": ("prevalence-percentile path: no defensible case count at this resolution (inherited state value, "
                   "proxy burden or condition set), so expected_yield_index = burden percentile x sensitivity x reach "
                   "is an index, not a count; expected_false_positives_upper treats every reached adult as a "
                   "non-case"),
    "level_D": "no burden measure (level D): no expected yield; the yield weight is removed like the burden weight",
}


def metric_link_columns(cb: Combo, base: pd.DataFrame, ev: dict, n_draws: int = MC_DRAWS) -> pd.DataFrame:
    """Performance record, reach, expected yield / false positives (every factor a column), the evidence-weighted
    status and its Monte Carlo (docs/ANALYSIS_PLAN_METRIC_LINK.md). Returns the new columns (index of `base`)."""
    perf = cb.perf or {}
    n = len(base)
    df = {}
    df["measurement_performance_object_id"] = perf.get("object_id")
    for col, key in PERF_COLS.items():
        v = perf.get(key)
        df[col] = v if v is not None and not (isinstance(v, float) and np.isnan(v)) else (
            np.nan if col not in ("measurement_performance_status",) else UNKNOWN)
    df["measurement_performance_status"] = perf.get("performance_status", UNKNOWN)
    df["measurement_performance_source_object_ids"] = perf.get("source_object_ids", "[]")
    known = performance_known(cb)
    df["measurement_evidence"] = np.asarray(ev["components"]["measurement_evidence"], dtype=float).ravel()
    r = cb.reach
    df["reach"] = r["reach"].to_numpy(float)
    df["reach_zcta_population_total"] = r["reach_zcta_population_total"].to_numpy(float)
    df["reach_population"] = r["reach_population"].to_numpy(float)
    df["reach_zcta_n"] = r["reach_zcta_n"].to_numpy(float)
    df["reach_radius_km"] = radius_km()
    df["adults_18plus"] = cb.adults
    df["reached_adults"] = cb.adults * df["reach"]
    df["expected_yield_basis"] = cb.yield_basis
    df["expected_yield_basis_note"] = YIELD_BASIS_NOTE[cb.yield_basis]
    se = _f(perf.get("op_sensitivity")) if known else np.nan
    sp = _f(perf.get("op_specificity")) if known else np.nan
    reach = df["reach"]
    nan = np.full(n, np.nan)
    df["burden_count"] = burden_count(cb)[0] if cb.yield_basis == "count" else nan
    if cb.yield_basis == "count":
        cnt = df["burden_count"]
        tp = cnt * se * reach
        fp = (1 - sp) * (cb.adults - cnt) * reach
        df["expected_detectable_cases"] = tp
        df["expected_false_positives"] = fp
        with np.errstate(invalid="ignore", divide="ignore"):
            df["expected_ppv"] = tp / (tp + fp)
            df["false_positives_per_detected_case"] = fp / tp
    else:
        for c in ("expected_detectable_cases", "expected_false_positives", "expected_ppv",
                  "false_positives_per_detected_case"):
            df[c] = nan
    df["expected_yield_index"] = (ev["burden_pct"][0] * se * reach) if cb.yield_basis == "percentile" else nan
    df["expected_false_positives_upper"] = ((1 - sp) * df["reached_adults"] if cb.yield_basis != "count" else nan)
    df["expected_yield"] = np.asarray(ev["components"]["expected_yield"], dtype=float).reshape(-1, n)[0]
    df["expected_yield_excluded_level_D"] = known and cb.yield_basis == "level_D"
    st = perf.get("performance_status", UNKNOWN)
    if st == "known":
        status = np.where(base["composite_evidence_weighted"].notna() if "composite_evidence_weighted" in base else False,
                          "ranked" if cb.yield_basis != "level_D" else
                          "ranked without expected_yield (level D: no burden measure; weight renormalised)",
                          "not ranked in this region (incomplete burden, no reach value or outside the universe)")
    elif st == "partial":
        status = np.full(n, "not ranked: measurement performance partial (AUROC only; no operating point at >= 90% "
                            "specificity, so expected yield is UNKNOWN)")
    else:
        status = np.full(n, "not ranked: measurement performance UNKNOWN (no record); never scored 0 or 1")
    df["evidence_weighted_status"] = status
    df["expected_yield_language"] = ML.YIELD_LANGUAGE + " " + ML.HEALTHY_CAVEAT
    mcs = monte_carlo_evidence(cb, n_draws) if n_draws else None
    for c in ("rank_mc_ew_p05", "rank_mc_ew_p50", "rank_mc_ew_p95", "p_top10_mc_ew", "p_top25_mc_ew",
              "expected_yield_mc_p05", "expected_yield_mc_p50", "expected_yield_mc_p95",
              "expected_false_positives_mc_p05", "expected_false_positives_mc_p50",
              "expected_false_positives_mc_p95"):
        df[c] = nan
    if mcs is not None:
        s = summarize_ranks(mcs["ranks"], "rank_mc_ew")
        for k, v in s.items():
            df[k] = v
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for key, pre in (("yield", "expected_yield_mc"), ("false_positives", "expected_false_positives_mc")):
                q = np.nanpercentile(mcs[key], [5, 50, 95], axis=0)
                df[f"{pre}_p05"], df[f"{pre}_p50"], df[f"{pre}_p95"] = q[0], q[1], q[2]
        df["mc_ew_n_draws"] = n_draws
    df["expected_yield_mc_quantity"] = {"count": "expected_detectable_cases", "percentile": "expected_yield_index",
                                        "level_D": UNKNOWN}[cb.yield_basis]
    return pd.DataFrame({k: (np.full(n, v, dtype=object) if isinstance(v, str) or v is None else v)
                         for k, v in df.items()}, index=base.index)


def incomplete_reason(cb: Combo, i: int) -> str:
    """Why region i has no (set) burden: a member with a defined burden measure has no value there."""
    miss = [m for m in cb.members if m.defined and m.level_d[i]]
    parts = [f"{m.condition_id} has a defined burden measure at {cb.level} level (level {m.defined_level}) but no "
             f"value for this region ({str(m.missing_reason[i] if m.missing_reason is not None else UNKNOWN).rstrip('.')})"
             for m in miss]
    if len(cb.members) > 1:
        return ("incomplete set burden: " + "; ".join(parts) + ". The set burden is computed only where every member "
                "with a defined burden has a value, so this region is excluded from the set's ranking and listed "
                "separately; resting its burden on the other member(s) alone would rank a different quantity.")
    return ("incomplete burden: " + "; ".join(parts) + ". The region is excluded from the ranking and listed "
            "separately rather than ranked on the other components with the burden weight renormalised.")


def row_uncertainties(cb: Combo, df: pd.DataFrame) -> list[list[str]]:
    out = []
    perf_lines = performance_uncertainties(cb)
    tau = next((m.tau for m in cb.members if m.inherited.any()), np.nan)
    mult_c = multiplier("C")
    any_defined = any(m.defined for m in cb.members)
    is_set = len(cb.members) > 1
    for i in range(len(df)):
        u = []
        complete = cb.burden_complete is None or bool(cb.burden_complete[i])
        if not complete:
            u.append(incomplete_reason(cb, i) + " Its non-burden components are kept.")
        for m in cb.members:
            f = m.frame
            if m.level_d[i]:
                if m.defined and cb.burden_rule == "complete":
                    pass  # incomplete region: stated once above
                elif not is_set or not any_defined:
                    tail = ("; no member of the set has a usable burden estimate" if is_set else "") + \
                           "; burden weight removed and the other weights renormalised (the ranking carries no burden " \
                           "information)"
                    u.append(f"No usable burden estimate (level D) for {m.condition_id}: excluded from burden-weighted "
                             "ranking" + tail)
                else:
                    u.append(f"No usable burden estimate (level D) for {m.condition_id}: no burden measure is defined "
                             f"for it at {cb.level} level, so it is not part of the set burden (the other members' "
                             "burden is used in every ranked region alike)")
            elif m.inherited[i]:
                u.append(f"{m.condition_id}: burden is the state HPS estimate for {df['state_abbr'].iloc[i]} inherited "
                         "by every county of the state (source resolution: state). It carries no within-state "
                         "variation and is not county prevalence; the Monte Carlo inflates the state estimate's "
                         f"uncertainty x{m.inh_mult:g} (inherited-burden multiplier) and adds a per-county deviation "
                         f"(tau = {tau:.2f} percentage points, assumed equal to the between-state SD).")
            elif "small-area" in str(f["burden_source_name"].iloc[i]):
                u.append(f"{m.condition_id}: burden is a modelled small-area estimate (BRFSS 2023 multilevel "
                         "regression and post-stratification), not observed county prevalence; its within-state "
                         "differences come from county age/sex/race/education composition, urbanicity and covariates "
                         "only (results/BRFSS_SAE_RESULTS.md); Kentucky and Pennsylvania values are model predictions "
                         "(absent from the 2023 public file).")
            elif str(f["burden_evidence_level"].iloc[i]) == "C":
                u.append(f"{m.condition_id}: burden is a level-C proxy ({f['burden_measure_label'].iloc[i]}); "
                         f"not the condition's prevalence; uncertainty band inflated x{mult_c:.1f}.")
            if cb.level == "county" and f["desert_burden_resolution"].iloc[i] is not None and \
                    "inherited" in str(f["desert_burden_resolution"].iloc[i]):
                u.append(f"{m.condition_id}: the diagnostic desert's burden term is the inherited state value.")
        if cb.level == "county" and bool(df["small_population_flag"].iloc[i]):
            u.append(f"Small population (< {min_population():,}): scored but not ranked (small denominators).")
        if "rank_mc_p05" in df:
            lo, hi = df["rank_mc_p05"].iloc[i], df["rank_mc_p95"].iloc[i]
            if lo == lo:
                u.append(f"Rank interval (5th-95th percentile over {int(df['mc_n_draws'].iloc[i])} draws of burden "
                         f"and weights): {lo:.0f}-{hi:.0f}.")
        u += perf_lines
        out.append(u)
    return out


def performance_uncertainties(cb: Combo) -> list[str]:
    """Row-independent statements about the measurement's performance evidence (plan section 5, rule 5)."""
    perf = cb.perf or {}
    st = perf.get("performance_status", UNKNOWN)
    out = []
    if st == "known":
        out.append(f"Measurement performance ({ML.tier_text(perf['quality_tier'])}): {perf.get('selected_record_id')} "
                   f"AUROC {perf['auroc']:.3f} ({perf['auroc_ci_low']:.3f}-{perf['auroc_ci_high']:.3f}), "
                   f"{perf.get('label_basis')} vs {perf.get('comparator')}; sensitivity {perf['op_sensitivity']:.2f} "
                   f"({perf['op_sensitivity_ci_low']:.2f}-{perf['op_sensitivity_ci_high']:.2f}) at specificity "
                   f"{perf['op_specificity']:.2f}. {ML.HEALTHY_CAVEAT}")
        out.append(ML.YIELD_LANGUAGE + " " + YIELD_BASIS_NOTE[cb.yield_basis] + ".")
    elif st == "partial":
        out.append(f"Measurement performance partial ({ML.tier_text(perf['quality_tier'])}): "
                   f"{perf.get('selected_record_id')} AUROC {perf['auroc']:.3f} but no operating point at >= 90% "
                   "specificity, so the expected yield is UNKNOWN; not ranked under evidence_weighted. The "
                   "equal-weight rank carries no measurement-performance information.")
    else:
        out.append("Measurement performance UNKNOWN / NOT AVAILABLE for this condition and measurement (no record in "
                   "this project's results or the curated published evidence): not ranked under evidence_weighted, "
                   "never scored 0 or 1. The equal-weight rank carries no measurement-performance information.")
    if perf.get("evidence_conflict"):
        out.append(f"Evidence conflict: {perf.get('evidence_conflict_note')}.")
    return out


def _provenance(df: pd.DataFrame, S: Spatial) -> pd.DataFrame:
    built = utc_now_iso()
    versions = "; ".join(f"{k}@{v}" for k, v in S.input_versions.items())
    df = add_provenance(
        df, data_layer="derived", source_name=f"{PRODUCER} (geo_condition_features, facilities, facility_trials, "
                                              "facility_nih_projects, providers, measurement_performance, "
                                              "geo_context_zcta__acs, configs/scoring.yaml, configs/relevance.yaml)",
        source_version=f"scoring build {built}; inputs: {versions}", retrieved_at=built,
        evidence_type="derived_score", source_record_id="object_id", source_geographic_resolution="county",
        evidence_level=df["burden_evidence_level"].astype(str).to_numpy(),
        provenance_notes=("Ecological derived score (places, facilities, registries; never people). Components are "
                          "kept; composite is for demonstration and must be read with its components, rank interval "
                          "and weight sensitivity. Curated relevance maps (implementer groups, condition specialties) "
                          "are assumptions. Registered trials/grants are research-activity signals, not evidence that "
                          "a measurement works. A ranked region is a candidate deployment opportunity."))
    res = np.where((df["geo_level"] == "state") | df["burden_inherited"].astype(bool), "state", "county")
    df["source_geographic_resolution"] = res
    return df


def score_combo(condition_id: str, measurement_id: str, level: str = "county", n_draws: int = MC_DRAWS,
                mc: bool = True, **kw) -> pd.DataFrame:
    cb = make_combo(condition_id, measurement_id, level, **kw)
    return _provenance(combo_frame(cb, n_draws=n_draws, mc=mc), spatial())


def add_joint_ranks(df: pd.DataFrame, weight_sets_=("equal", EVIDENCE_WEIGHT_SET)) -> None:
    """rank_<ws>_joint: per condition x level, every ranking-eligible region x measurement pair with a composite,
    ranked by the composite (ties by geo_id, then measurement_id). Under evidence_weighted, measurements without known
    performance have no pairs (docs/ANALYSIS_PLAN_METRIC_LINK.md section 5)."""
    for ws in weight_sets_:
        col = f"composite_{ws}"
        out = pd.Series(np.nan, index=df.index)
        if col in df:
            for _, g in df.groupby(["condition_id", "geo_level"]):
                ok = g[g["rank_eligible"].astype(bool) & g[col].notna()]
                ok = ok.sort_values([col, "geo_id", "measurement_id"], ascending=[False, True, True])
                out.loc[ok.index] = np.arange(1, len(ok) + 1, dtype=float)
        df[f"rank_{ws}_joint"] = out


def all_condition_ids() -> list[str]:
    return SINGLE_CONDITIONS + list(condition_sets())


def build(conditions: list[str] | None = None, measurements: list[str] | None = None, levels=LEVELS,
          n_draws: int = MC_DRAWS, write: bool = True) -> pd.DataFrame:
    t0 = time.time()
    S = spatial()
    conditions = conditions or all_condition_ids()
    measurements = measurements or MEASUREMENTS
    frames = []
    for level in levels:
        for mid in measurements:
            for cid in conditions:
                cb = make_combo(cid, mid, level)
                frames.append(combo_frame(cb, n_draws=n_draws))
                print(f"[scoring] {cid} x {mid} x {level}: {len(frames[-1])} rows ({time.time() - t0:.0f}s)",
                      file=sys.stderr, flush=True)
    df = pd.concat(frames, ignore_index=True)
    add_joint_ranks(df)
    df = _provenance(df, S)
    assert df["object_id"].is_unique
    if write:
        write_table(df, TABLE, producer=PRODUCER,
                    description="Deployment-opportunity components, composites under every weight set, and Monte "
                                "Carlo rank intervals per condition (or condition set) x measurement x geography. "
                                "docs/SCORING.md")
    return df


if __name__ == "__main__":
    build()
