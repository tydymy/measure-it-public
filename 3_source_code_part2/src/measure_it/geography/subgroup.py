"""County and state burden of a device-defined subgroup (HARMONIZATION_CONTRACT section 2 and 4).

For a computable phenotype with `subgroup_strata_fractions` (fraction of condition cases that are device-positive, by
age band x sex, with a pooled `all`/`all` row) this module estimates how many people in each county and state belong
to the subgroup:

    long_covid    adults in subgroup = sum over ACS adult age x sex cells of
                  ACS adults x SAE Long COVID prevalence for that county, age7 band and sex (geography.sae draws)
                  x subgroup fraction for the matching age band x sex (pooled fraction where the cell is suppressed);
                  burden_basis = "brfss_2023_mrp_sae".
    lyme_disease  annual reported Lyme cases in the subgroup = 2023 county reported cases (CDC surveillance, county of
                  residence, all ages; 2022 case definition) x pooled subgroup fraction (the county series has no age
                  or sex); burden_basis = "cdc_lyme_reported_cases_2023".
    ptlds         derived: 2023 reported Lyme cases x published PTLD fraction among treated early Lyme
                  (Aucott et al. 2022, Int J Infect Dis 116:230-237, PMID 35066160: 13.7% of 234 = 32/234 met PTLD
                  criteria vs 4.1% of 49 without Lyme) x pooled subgroup fraction; burden_basis =
                  "derived_ptlds_from_lyme_cases_aucott2022" (evidence level C).
    any other     the condition's primary burden in `geo_condition_features`: counts only where that table defines a
                  count at the row's own resolution; inherited state values give state rows only (guardrail 5).

Uncertainty: Monte Carlo over the SAE posterior draws (or Gamma(k, 1) draws for reported counts k, a Jeffreys beta
for the PTLD fraction) and over the subgroup fractions (Jeffreys beta from n_positive / n_cases per stratum; one draw
per stratum shared by all geographies, because the fraction is a cohort quantity). 95% intervals are the 2.5/97.5
percentiles.

Output: geo_subgroup_burden (non-partition; data_layer "derived"; one row per phenotype x geography). If
`subgroup_strata_fractions` does not exist yet, nothing is written and the reason is printed.

Run: uv run python -m measure_it.geography.subgroup
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from .. import config
from ..config import UNKNOWN, utc_now_iso
from ..geography.crosswalk import STATE_FIPS
from ..provenance import add_provenance
from ..config import PROCESSED
from ..store import read_table, table_exists, write_table
from . import sae

PRODUCER = "measure_it.geography.subgroup"
TABLE = "geo_subgroup_burden"
FRACTIONS_TABLE = "subgroup_strata_fractions"
N_DRAWS = 500
LYME_YEAR = 2023
PTLDS_SOURCE = {"citation": "Aucott JN, Yang T, Yoon I, Powell D, Geller SA, Rebman AW. Risk of post-treatment Lyme "
                            "disease in patients with ideally-treated early Lyme disease: a prospective cohort study. "
                            "Int J Infect Dis. 2022;116:230-237. doi:10.1016/j.ijid.2022.01.033. PMID 35066160",
                "n_positive": 32, "n_total": 234, "point": 0.137,
                "comparison_without_lyme": "4.1% of 49 participants without prior Lyme met the same criteria",
                "attributable_excess": 0.137 - 0.041}
LYME_CAVEATS = ("reported surveillance cases (under-ascertained; CDC estimates true annual incidence several-fold "
                "higher); 2022+ case definition: high-incidence states report on laboratory evidence alone, so 2023 "
                "is not comparable with pre-2022 years and high- vs low-incidence states are ascertained differently; "
                "all ages (not adults); county of residence")


# =================================================================================================== pure helpers
def parse_age_band(band: str) -> tuple[int, int] | None:
    """'40-49' -> (40, 49); '90+' -> (90, 200); 'all' -> None."""
    b = str(band).strip().lower()
    if b in ("all", "", "nan", "none"):
        return None
    m = re.fullmatch(r"(\d+)\s*\+", b)
    if m:
        return int(m.group(1)), 200
    m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", b)
    if m:
        return int(m.group(1)), int(m.group(2))
    raise ValueError(f"unparsed age band {band!r}")


def _beta_params(row) -> tuple[float, float] | None:
    n, k = row.get("n_cases"), row.get("n_positive")
    if n is not None and k is not None and pd.notna(n) and pd.notna(k) and n > 0:
        return float(k) + 0.5, float(n) - float(k) + 0.5          # Jeffreys
    f, lo, hi = row.get("fraction"), row.get("ci_low"), row.get("ci_high")
    if f is None or pd.isna(f):
        return None
    if lo is not None and hi is not None and pd.notna(lo) and pd.notna(hi) and hi > lo and 0 < f < 1:
        sd = (hi - lo) / 3.92
        common = f * (1 - f) / sd ** 2 - 1
        if common > 0:
            return f * common, (1 - f) * common
    return None


def fraction_draws(fr: pd.DataFrame, n: int, rng: np.random.Generator) -> dict[tuple, np.ndarray]:
    """(age_band, sex) -> fraction draws for every usable (not suppressed, fraction present) stratum."""
    out = {}
    for _, r in fr.iterrows():
        sup = r.get("suppressed", False)
        if (pd.notna(sup) and bool(sup)) or pd.isna(r.get("fraction")):
            continue
        ab = _beta_params(r)
        key = (str(r["age_band"]).lower(), str(r["sex"]).lower())
        out[key] = rng.beta(*ab, size=n) if ab else np.full(n, float(r["fraction"]))
    if ("all", "all") not in out:
        raise ValueError("subgroup_strata_fractions needs a usable pooled row (age_band='all', sex='all')")
    return out


def stratum_key(age_low: int, sex: str, available: dict) -> tuple[tuple, str]:
    """Fraction stratum for an ACS cell (lower age bound, sex): exact band x sex, then band x all, all x sex, pooled.

    The ACS cell is assigned to the band containing its lower bound (ACS 85+ goes to the band containing 85)."""
    bands = {k[0]: parse_age_band(k[0]) for k in available if k[0] != "all"}
    band = next((b for b, rng_ in bands.items() if rng_ and rng_[0] <= age_low <= rng_[1]), None)
    for key, how in (((band, sex), "age_band_x_sex"), ((band, "all"), "age_band"), (("all", sex), "sex"),
                     (("all", "all"), "pooled")):
        if key in available:
            return key, how
    raise KeyError("no pooled fraction")


def subgroup_counts(N: np.ndarray, prev: np.ndarray, frac: np.ndarray) -> np.ndarray:
    """N (C, F) population, prev (C, F, D) prevalence draws, frac (F, D) fraction draws -> (C, D) subgroup counts."""
    return np.einsum("cf,cfd,fd->cd", N, prev, frac)


def summarise(draws: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return draws.mean(1), np.quantile(draws, 0.025, axis=1), np.quantile(draws, 0.975, axis=1)


def _rows(ids, level, base, sub, pop, extra: dict) -> pd.DataFrame:
    bm, bl, bh = summarise(base)
    sm, sl, sh = summarise(sub)
    rate = sub / np.where(pop > 0, pop, np.nan)[:, None] * 1e5
    rm, rl, rh = summarise(np.nan_to_num(rate))
    return pd.DataFrame({"geo_id": ids, "geo_level": level, "base_population": pop,
                         "base_burden": bm, "base_burden_ci_low": bl, "base_burden_ci_high": bh,
                         "subgroup_estimate": sm, "subgroup_ci_low": sl, "subgroup_ci_high": sh,
                         "subgroup_rate_per_100k": rm, "subgroup_rate_ci_low": rl, "subgroup_rate_ci_high": rh,
                         **extra})


def _aggregate(ids: np.ndarray, base: np.ndarray, sub: np.ndarray, pop: np.ndarray, key: np.ndarray):
    keys = pd.unique(key)
    B = np.vstack([base[key == k].sum(0) for k in keys])
    S = np.vstack([sub[key == k].sum(0) for k in keys])
    P = np.array([pop[key == k].sum() for k in keys])
    return keys, B, S, P


# =================================================================================================== bases
def long_covid_rows(fr: pd.DataFrame, fdraw: dict, measure: str = sae.PRIMARY) -> tuple[pd.DataFrame, dict]:
    draws = sae.load_draws(measure)
    if draws is None:
        return pd.DataFrame(), {"long_covid": "SAE draws missing (run measure_it.geography.sae)"}
    cells = read_table(sae.CELLS_TABLE)
    geo_ids = draws["geo_id"].astype(str)
    D = min(N_DRAWS, draws["cell_draws"].shape[2])
    cell_draws = draws["cell_draws"][:, :, :D]
    labels = list(draws["cell_labels"].astype(str))
    fine = cells.pivot_table(index="geo_id", columns=["sex", "age_low"], values="population", aggfunc="sum")
    fine = fine.reindex(geo_ids).fillna(0.0)
    cols = list(fine.columns)
    N = fine.to_numpy(float)
    age7 = cells.drop_duplicates(["sex", "age_low"]).set_index(["sex", "age_low"])["age7"]
    cell_idx = np.array([labels.index(f"{age7[(s, a)]}|{s}") for s, a in cols])
    prev = cell_draws[:, cell_idx, :]                                    # (C, F, D)
    keys, hows = zip(*(stratum_key(int(a), s, fdraw) for s, a in cols))
    frac = np.vstack([fdraw[k][:D] for k in keys])                        # (F, D)
    sub = subgroup_counts(N, prev, frac)
    base = np.einsum("cf,cfd->cd", N, prev)
    pop = N.sum(1)
    w = N.sum(0)
    fallback_share = float(w[np.array(hows) != "age_band_x_sex"].sum() / w.sum())
    extra = {"burden_basis": "brfss_2023_mrp_sae",
             "base_metric": "adults 18+ with current Long COVID (modelled prevalent count)",
             "base_population_kind": "adults 18+ (ACS 2020-2024)",
             "base_source_table": sae.TABLE, "base_measure_id": measure, "base_evidence_level": "A (modelled)",
             "fraction_basis": "age band x sex strata (" + ", ".join(sorted(set(hows))) + ")",
             "share_population_on_fallback_fraction": fallback_share,
             "assumptions": "SAE age7 x sex prevalence applies to every ACS age cell inside the band; subgroup "
                            "fraction of the cohort applies to every geography (no geographic variation in the "
                            "device-positive share is modelled)"}
    county = _rows(geo_ids, "county", base, sub, pop, extra)
    st = np.array([g[:2] for g in geo_ids])
    k, B_, S_, P_ = _aggregate(geo_ids, base, sub, pop, st)
    state = _rows(k, "state", B_, S_, P_, extra)
    k, B_, S_, P_ = _aggregate(geo_ids, base, sub, pop, np.array(["US"] * len(geo_ids)))
    nat = _rows(k, "national", B_, S_, P_, extra)
    return pd.concat([county, state, nat], ignore_index=True), {}


def lyme_base(rng: np.random.Generator, D: int) -> tuple[pd.DataFrame, np.ndarray]:
    ly = read_table("geo_condition_burden__cdc_lyme")
    c = ly[(ly["year"] == LYME_YEAR) & (ly["measure_id"] == "lyme_reported_cases") & (ly["geo_level"] == "county")
           & ly["in_canonical_geographies"].astype(bool)].copy()
    c = c.dropna(subset=["value"]).sort_values("geo_id")
    cases = c["value"].to_numpy(float)
    # Poisson noise around the reported count: Gamma(k, 1) draws (posterior of the Poisson mean under the 1/lambda
    # prior; mean = the reported count k). A county reporting 0 cases keeps 0 (no imputed cases).
    draws = rng.gamma(np.where(cases > 0, cases, 1.0)[:, None], 1.0, size=(len(cases), D))
    draws[cases == 0] = 0.0
    return c, draws


def lyme_rows(fdraw: dict, rng, ptlds: bool = False) -> tuple[pd.DataFrame, dict]:
    D = N_DRAWS
    c, base = lyme_base(rng, D)
    frac = fdraw[("all", "all")][:D]
    if ptlds:
        pt = rng.beta(PTLDS_SOURCE["n_positive"] + 0.5, PTLDS_SOURCE["n_total"] - PTLDS_SOURCE["n_positive"] + 0.5, D)
        base = base * pt[None, :]
    sub = base * frac[None, :]
    pop = c["denominator"].to_numpy(float)
    common = {"base_population_kind": "residents, all ages (Census PEP July 2023, cdc_lyme denominator)",
              "base_source_table": "geo_condition_burden__cdc_lyme", "base_measure_id": "lyme_reported_cases",
              "fraction_basis": "pooled (age_band='all', sex='all'): county case counts have no age or sex",
              "share_population_on_fallback_fraction": 1.0}
    if ptlds:
        extra = {**common, "burden_basis": "derived_ptlds_from_lyme_cases_aucott2022",
                 "base_metric": "derived annual PTLD cases among reported Lyme cases (reported cases x 13.7%)",
                 "base_evidence_level": "C (derived: surveillance count x published cohort fraction)",
                 "assumptions": (f"PTLD fraction among treated early Lyme 13.7% (32/234; Jeffreys beta draws), "
                                 f"{PTLDS_SOURCE['citation']}; {PTLDS_SOURCE['comparison_without_lyme']} (attributable "
                                 f"excess ~{100 * PTLDS_SOURCE['attributable_excess']:.1f}%, a lower sensitivity bound); "
                                 "a single ideally-treated Maryland cohort applied to all reported cases; "
                                 + LYME_CAVEATS)}
    else:
        extra = {**common, "burden_basis": f"cdc_lyme_reported_cases_{LYME_YEAR}",
                 "base_metric": f"reported Lyme disease cases, {LYME_YEAR} (annual incident, all ages)",
                 "base_evidence_level": "A (surveillance count)",
                 "assumptions": "Gamma(k, 1) draws for the Poisson mean of the reported count k; " + LYME_CAVEATS}
    ids = c["geo_id"].to_numpy()
    county = _rows(ids, "county", base, sub, pop, extra)
    county["lyme_state_incidence_category"] = c["state_incidence_category"].to_numpy()
    county["lyme_state_cases_missing_from_county_file_pct"] = c["state_cases_missing_from_county_file_pct"].to_numpy()
    k, B_, S_, P_ = _aggregate(ids, base, sub, pop, c["state_fips"].to_numpy())
    state = _rows(k, "state", B_, S_, P_, {**extra, "assumptions": extra["assumptions"] + "; state = sum of its "
                                           "counties (cases the state could not assign to a county are not included)"})
    k, B_, S_, P_ = _aggregate(ids, base, sub, pop, np.array(["US"] * len(ids)))
    nat = _rows(k, "national", B_, S_, P_, extra)
    return pd.concat([county, state, nat], ignore_index=True), {}


def generic_rows(condition_id: str, fdraw: dict) -> tuple[pd.DataFrame, dict]:
    """Fallback for a condition without a modelled base: geo_condition_features primary burden, pooled fraction."""
    if not table_exists("geo_condition_features"):
        return pd.DataFrame(), {condition_id: "geo_condition_features missing"}
    f = read_table("geo_condition_features")
    f = f[(f["condition_id"] == condition_id) & f["geo_level"].isin(["county", "state"])].copy()
    if f.empty or f["burden_value"].isna().all():
        return pd.DataFrame(), {condition_id: "no burden measure (level D)"}
    frac = fdraw[("all", "all")]
    fm, fl, fh = frac.mean(), np.quantile(frac, 0.025), np.quantile(frac, 0.975)
    inherited = f["burden_inherited"].fillna(False).astype(bool)
    count = f["burden_count"].where(~inherited)
    out = pd.DataFrame({
        "geo_id": f["geo_id"], "geo_level": f["geo_level"],
        "base_population": f.get("population_adults_18plus"),
        "base_burden": count, "base_burden_ci_low": np.nan, "base_burden_ci_high": np.nan,
        "subgroup_estimate": count * fm, "subgroup_ci_low": count * fl, "subgroup_ci_high": count * fh,
        "subgroup_rate_per_100k": np.nan, "subgroup_rate_ci_low": np.nan, "subgroup_rate_ci_high": np.nan,
        "burden_basis": np.where(inherited, "inherited_state_value_no_county_count",
                                 "geo_condition_features_primary_burden:" + f["burden_measure_id"].astype(str)),
        "base_metric": f["burden_count_kind"].fillna(UNKNOWN), "base_population_kind": "adults 18+ (ACS)",
        "base_source_table": "geo_condition_features", "base_measure_id": f["burden_measure_id"],
        "base_evidence_level": f["burden_evidence_level"],
        "fraction_basis": "pooled (age_band='all', sex='all')", "share_population_on_fallback_fraction": 1.0,
        "assumptions": np.where(count.isna(), "no count defensible at this resolution (" +
                                f["burden_count_note"].fillna("no count").astype(str) + ")",
                                "count as defined in geo_condition_features x pooled fraction"),
    })
    return out, {}


# =================================================================================================== build
def build(fractions: pd.DataFrame, seed: int = config.SEED) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    frames, notes = [], {}
    for pid, fr in fractions.groupby("phenotype_id"):
        cond = str(fr["condition_id"].iloc[0])
        fdraw = fraction_draws(fr, N_DRAWS, rng)
        pooled = fr[(fr["age_band"].astype(str).str.lower() == "all") & (fr["sex"].astype(str).str.lower() == "all")]
        if cond == "long_covid":
            t, n = long_covid_rows(fr, fdraw)
            if t.empty:
                t, n2 = generic_rows(cond, fdraw)
                n.update(n2)
        elif cond in ("lyme_disease", "ptlds"):
            t, n = lyme_rows(fdraw, rng, ptlds=cond == "ptlds")
        else:
            t, n = generic_rows(cond, fdraw)
        notes.update({f"{pid}:{k}": v for k, v in n.items()})
        if t.empty:
            continue
        t.insert(0, "condition_id", cond)
        t.insert(0, "phenotype_id", pid)
        t["pooled_fraction"] = float(pooled["fraction"].iloc[0])
        t["pooled_fraction_n_cases"] = pooled["n_cases"].iloc[0] if "n_cases" in pooled else np.nan
        frames.append(t)
    if not frames:
        return pd.DataFrame(), notes
    out = pd.concat(frames, ignore_index=True)
    out["state_fips"] = np.where(out["geo_level"] == "national", None, out["geo_id"].astype(str).str[:2])
    out["state_abbr"] = out["state_fips"].map(STATE_FIPS)
    out["n_draws"] = N_DRAWS
    out["object_id"] = "geo:" + out["geo_id"].astype(str)
    return out, notes


def _drop_stale() -> None:
    """No usable phenotype: remove an earlier geo_subgroup_burden so the table never outlives its inputs."""
    for suffix in (".parquet", ".meta.json"):
        (PROCESSED / f"{TABLE}{suffix}").unlink(missing_ok=True)


def run() -> dict:
    if not table_exists(FRACTIONS_TABLE):
        msg = f"{FRACTIONS_TABLE} not found (produced by `measure-it byod subgroup`); geo_subgroup_burden not written"
        print(msg)
        _drop_stale()
        return {"written": False, "reason": msg}
    fr = read_table(FRACTIONS_TABLE)
    # demo-tagged phenotypes (the synthetic tutorial) are left out unless MEASURE_IT_BYOD_DEMO=1, as in the metric
    # link (byod/records.py); the flag travels on the aggregate fractions table itself
    import os
    if "demo" in fr.columns and os.environ.get("MEASURE_IT_BYOD_DEMO") != "1":
        fr = fr[~fr["demo"].fillna(False).astype(bool)]
    out, notes = build(fr) if not fr.empty else (pd.DataFrame(), {"reason": "only demo phenotypes"})
    if out.empty:
        print(json.dumps(notes, indent=2, default=str))
        _drop_stale()
        return {"written": False, "notes": notes}
    g = read_table("geographies").set_index("geo_id")["name"]
    out["geo_name"] = out["geo_id"].map(g).fillna(out["geo_id"].map({"US": "United States (50 states + DC)"}))
    ev = {"brfss_2023_mrp_sae": "modeled_small_area_estimate"}
    out = add_provenance(
        out, data_layer="derived", source_name="Measure-It subgroup burden (SAE / surveillance x cohort fraction)",
        source_version=f"{PRODUCER} {utc_now_iso()}", retrieved_at=utc_now_iso(),
        evidence_type="modeled_small_area_estimate",
        source_record_id=lambda d: d["phenotype_id"] + "|" + d["geo_level"] + "|" + d["geo_id"].astype(str),
        source_geographic_resolution="county",
        provenance_notes="Estimated people in a device-defined subgroup = base burden x cohort subgroup fraction; "
                         "see burden_basis and assumptions. Ecological, aggregate only.")
    basis = out["burden_basis"].astype(str)
    out["evidence_type"] = np.select([basis.isin(list(ev)), basis.str.startswith("cdc_lyme")],
                                     ["modeled_small_area_estimate", "surveillance_case_count"],
                                     "derived_evidence_summary")
    out["evidence_level"] = out["base_evidence_level"].astype(str).str[0]
    out["source_geographic_resolution"] = np.where(
        out["burden_basis"].eq("inherited_state_value_no_county_count") | out["geo_level"].eq("state"), "state",
        np.where(out["geo_level"].eq("national"), "national", "county"))
    write_table(out, TABLE, producer=PRODUCER, description="device-defined subgroup burden by county/state")
    summary = {"written": True, "rows": len(out), "phenotypes": sorted(out["phenotype_id"].unique()), "notes": notes,
               "national": out[out["geo_level"] == "national"][["phenotype_id", "subgroup_estimate", "subgroup_ci_low",
                                                                 "subgroup_ci_high"]].to_dict("records")}
    print(json.dumps(summary, indent=2, default=str))
    return summary


if __name__ == "__main__":
    run()
