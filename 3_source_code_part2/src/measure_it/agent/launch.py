"""Launch agent part C: the launch ranking ("greatest impact", docs/AGENT_CONTRACT.md section 4).

For a condition x measurement:

(a) engine ranking. The scored grid (`deployment_opportunities`, county level) under `evidence_weighted` when the
    measurement has a KNOWN performance record for the condition (after `byod deploy`, or a public record), else under
    `equal`. Every component is reported (burden_pct, vulnerability_pct, diagnostic_desert_pct, clinic_capacity_pct
    with its measurement-specific basis, research_readiness_pct, and measurement_evidence / expected_yield when
    known), with the engine's Monte Carlo rank interval (equal weights: burden + Dirichlet weight draws, 1000 draws;
    evidence_weighted: its own MC).
(b) device-subgroup re-score, when a device subgroup exists for the condition: the same engine combination
    (opportunity.make_combo + opportunity.evaluate, as scoring/stage_layers.py does) with the burden member values of
    the subgroup's condition replaced by the estimated device-subgroup rate per 100k adults per county
    (geography.subgroup.build, computed in memory for this phenotype). Components, composite and rank are reported
    under the same weight set as (a); an equal-weight Monte Carlo interval draws Dirichlet weights and subgroup rates
    (normal from the subgroup's 95% interval). The estimated number of adults in the subgroup is reported next to it.

The launch ranking ("greatest impact") is (b) where a subgroup exists, else (a); the result names which.

Conditions outside the scored grid (e.g. lyme_disease, ptlds) use the BURDEN-ONLY fallback the outreach module
implements (primary county burden measure), and with a subgroup a SUBGROUP-BURDEN-ONLY ranking (estimated subgroup
cases per county); both are labelled as burden-only, not deployment composites.

For the top N counties of each ranking: billing clinicians (geo_measurement_capacity), measurement experience sites
(measurement_experience_sites) and candidate facilities (clinic_registry: clinic candidates whose specialty groups
include one of the measurement's implementer groups); names only in the outreach workbook.
"""
from __future__ import annotations

import json
import re
import warnings
from unittest import mock

import numpy as np
import pandas as pd

from ..config import SEED, UNKNOWN
from ..store import read_table, table_exists

ENGINE_COLS = ["geo_id", "geo_name", "state_abbr", "population_total", "adults_18plus", "rank_equal",
               "rank_evidence_weighted", "composite_equal", "composite_evidence_weighted", "burden_value",
               "burden_measure_label", "burden_evidence_level", "burden_source_name", "burden_pct",
               "vulnerability_pct", "diagnostic_desert_pct", "clinic_capacity_pct", "clinic_capacity_basis",
               "research_readiness_pct", "measurement_evidence", "expected_yield", "expected_detectable_cases",
               "rank_mc_p05", "rank_mc_p50", "rank_mc_p95", "p_top10_mc", "rank_mc_ew_p05", "rank_mc_ew_p50",
               "rank_mc_ew_p95", "p_top10_mc_ew", "measurement_performance_status",
               "measurement_performance_record_id", "uncertainties"]
CAPACITY_NOTE = {
    "measurement_activity": "clinic_capacity = measurement-specific clinician activity (Medicare FFS billing of a "
                            "dedicated code; facilities.activity)",
    "implementer_density": "clinic_capacity = implementer-group provider and facility density (the measurement has no "
                           "dedicated billing code)",
}
MC_DRAWS = 1000


# ------------------------------------------------------------------------------------------------------------------
# grid membership, subgroup burden
# ------------------------------------------------------------------------------------------------------------------
def in_scored_grid(condition: str, measurement: str) -> bool:
    from ..scoring import opportunity as O
    return condition in O.all_condition_ids() and measurement in O.MEASUREMENTS


def subgroup_burden(ds: str, phenotype_id: str) -> tuple[pd.DataFrame | None, dict]:
    """geo_subgroup_burden rows for ONE phenotype, computed in memory (the shared table is not written)."""
    from ..byod import common as K
    from ..geography import subgroup as G
    if not table_exists(K.STRATA_TABLE):
        return None, {"reason": f"{K.STRATA_TABLE} missing"}
    fr = read_table(K.STRATA_TABLE)
    fr = fr[(fr["dataset_id"] == ds) & (fr["phenotype_id"] == phenotype_id)]
    if fr.empty:
        return None, {"reason": f"no strata fractions for {phenotype_id}"}
    out, notes = G.build(fr.reset_index(drop=True))
    if out.empty:
        return None, notes
    names = read_table("geographies", columns=["geo_id", "name"]).drop_duplicates("geo_id").set_index("geo_id")["name"]
    out["geo_name"] = out["geo_id"].map(names).fillna(out["geo_id"].map({"US": "United States (50 states + DC)"}))
    return out, notes


# ------------------------------------------------------------------------------------------------------------------
# (a) engine ranking
# ------------------------------------------------------------------------------------------------------------------
def _clean(v):
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def _records(df: pd.DataFrame) -> list[dict]:
    return [{k: _clean(v) for k, v in r.items()} for r in df.to_dict("records")]


def engine_ranking(condition: str, measurement: str, top: int = 10) -> dict:
    """(a) the engine's own county ranking (scored grid) or the outreach module's burden-only fallback."""
    if not in_scored_grid(condition, measurement):
        return burden_only_ranking(condition, measurement, top)
    d = read_table("deployment_opportunities")
    d = d[(d["condition_id"] == condition) & (d["measurement_id"] == measurement) & (d["geo_level"] == "county")]
    known = bool((d["measurement_performance_status"] == "known").any() and d["rank_evidence_weighted"].notna().any())
    ws = "evidence_weighted" if known else "equal"
    rank_col = f"rank_{ws}"
    top_rows = d[d[rank_col] <= top].sort_values(rank_col)
    cols = [c for c in ENGINE_COLS if c in d.columns]
    rows = top_rows[cols].copy()
    rows.insert(0, "rank", top_rows[rank_col].astype(int).to_numpy())
    rows["composite"] = top_rows[f"composite_{ws}"].to_numpy()
    if ws == "evidence_weighted":
        rows["mc_rank_p05"], rows["mc_rank_p95"] = top_rows["rank_mc_ew_p05"], top_rows["rank_mc_ew_p95"]
        rows["mc_p_top10"] = top_rows["p_top10_mc_ew"]
    else:
        rows["mc_rank_p05"], rows["mc_rank_p95"] = top_rows["rank_mc_p05"], top_rows["rank_mc_p95"]
        rows["mc_p_top10"] = top_rows["p_top10_mc"]
    rows["uncertainties"] = rows["uncertainties"].map(lambda s: json.loads(s) if isinstance(s, str) else [])
    basis = str(d["clinic_capacity_basis"].dropna().iloc[0]) if d["clinic_capacity_basis"].notna().any() else UNKNOWN
    perf = d.iloc[0] if len(d) else None
    return {
        "basis": "deployment_opportunities", "kind": "deployment_composite", "weight_set": ws,
        "why": ("the measurement has a KNOWN performance record for this condition: evidence_weighted ranking"
                if known else "no known performance record for this condition x measurement (UNKNOWN is never "
                              "imputed): equal-weight ranking"),
        "n_ranked": int(d[rank_col].notna().sum()), "clinic_capacity_basis": basis,
        "clinic_capacity_note": CAPACITY_NOTE.get(basis, basis),
        "monte_carlo": ("evidence_weighted Monte Carlo (burden, sensitivity/specificity and Dirichlet weight draws)"
                        if known else "equal-weight Monte Carlo (burden draws + Dirichlet(20 x equal) weights, "
                                      "1000 draws)") + "; 5th-95th percentile of the rank and P(top 10)",
        "performance": None if perf is None else {
            "status": perf["measurement_performance_status"], "record_id": _clean(perf["measurement_performance_record_id"]),
            "tier_label": _clean(perf.get("measurement_performance_tier_label")),
            "auroc": _clean(perf.get("measurement_performance_auroc")),
            "auroc_ci_low": _clean(perf.get("measurement_performance_auroc_ci_low")),
            "auroc_ci_high": _clean(perf.get("measurement_performance_auroc_ci_high")),
            "op_sensitivity": _clean(perf.get("op_sensitivity")), "op_specificity": _clean(perf.get("op_specificity")),
            "measurement_evidence": _clean(perf.get("measurement_evidence"))},
        "counties": _records(rows),
    }


def burden_only_ranking(condition: str, measurement: str, top: int) -> dict:
    """Outside the scored grid: outreach.select_counties' BURDEN-ONLY fallback on the primary county burden measure
    (phenotype_id that matches nothing, so a shared geo_subgroup_burden table from another dataset is never used)."""
    from .. import outreach as OR
    cond = OR.resolve_condition(condition)
    meas = OR.resolve_measurement(measurement)
    c, basis = OR.select_counties(cond, meas, top, phenotype_id="__launch_primary_burden_only__")
    c = c.rename(columns={"county_fips": "geo_id"})
    return {"basis": basis["basis"], "kind": "burden_only", "weight_set": None,
            "why": f"{condition} x {measurement} is outside the scored grid: {basis['description']}",
            "universe": basis.get("universe"), "n_ranked": None, "clinic_capacity_basis": None,
            "monte_carlo": "not available (burden-only ranking)", "performance": None, "counties": _records(c)}


# ------------------------------------------------------------------------------------------------------------------
# (b) device-subgroup re-score
# ------------------------------------------------------------------------------------------------------------------
def subgroup_rescore(condition: str, measurement: str, sub: pd.DataFrame, top: int = 10, n_draws: int = MC_DRAWS,
                     seed: int = SEED) -> dict:
    """Re-score the engine combination with the subgroup's condition burden replaced by the subgroup rate."""
    from ..scoring import opportunity as O
    county = sub[sub["geo_level"] == "county"].drop_duplicates("geo_id").set_index("geo_id")
    sub_cond = str(sub["condition_id"].iloc[0])
    cb = O.make_combo(condition, measurement, "county")
    gid = cb.geo["geo_id"].astype(str).to_numpy()
    rate = county["subgroup_rate_per_100k"].reindex(gid).to_numpy(float)
    lo = county["subgroup_rate_ci_low"].reindex(gid).to_numpy(float)
    hi = county["subgroup_rate_ci_high"].reindex(gid).to_numpy(float)
    est = county["subgroup_estimate"].reindex(gid).to_numpy(float)
    est_lo = county["subgroup_ci_low"].reindex(gid).to_numpy(float)
    est_hi = county["subgroup_ci_high"].reindex(gid).to_numpy(float)
    replaced = [i for i, m in enumerate(cb.members) if m.condition_id == sub_cond]
    if not replaced:
        return {"status": "skipped", "reason": f"the device subgroup is defined within {sub_cond}, which is not a "
                                               f"member of {condition}"}
    vals = [rate if i in replaced else m.value for i, m in enumerate(cb.members)]
    ev = O.evaluate(cb, member_values=vals)
    ws = "evidence_weighted" if O.performance_known(cb) else "equal"
    comp = O.named_composites(cb, ev)[ws]
    base = O.named_composites(cb)[ws]
    elig = cb.eligible & ~np.isnan(comp)
    elig_base = cb.eligible & ~np.isnan(base)
    rank = O.rank_rows(comp, elig)[0]
    rank_base = O.rank_rows(base, elig_base)[0]
    # equal-weight Monte Carlo: Dirichlet weights x subgroup-rate draws (normal from the 95% interval, >= 0)
    rng = np.random.default_rng(seed)
    w0 = O.default_weights()
    alpha = np.array([w0[k] for k in O.COMPONENTS]) * O.DIRICHLET_CONCENTRATION
    sd = np.where(np.isfinite(hi - lo), (hi - lo) / 3.92, 0.0)
    R = np.full((n_draws, len(gid)), np.nan, dtype=np.float32)
    done = 0
    while done < n_draws:
        D = min(O.MC_CHUNK, n_draws - done)
        draw = np.maximum(rate[None, :] + rng.standard_normal((D, len(gid))) * sd[None, :], 0.0)
        draw[:, np.isnan(rate)] = np.nan
        v2 = [draw if i in replaced else m.value for i, m in enumerate(cb.members)]
        W = rng.dirichlet(alpha, size=D)
        C = O.evaluate(cb, member_values=v2, weights=W, order=O.COMPONENTS)["composite"]
        R[done:done + D] = O.rank_rows(C, elig)
        done += D
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q05, q50, q95 = np.nanpercentile(R, [5, 50, 95], axis=0)
        p10 = np.nanmean(np.where(np.isnan(R), np.nan, R <= top), axis=0)
    comps = ev["components"]
    order = np.argsort(np.where(np.isnan(rank), np.inf, rank), kind="stable")
    top_idx = [i for i in order if rank[i] <= top]
    names = cb.geo["geo_name"].to_numpy()
    st = cb.geo["state_abbr"].to_numpy() if "state_abbr" in cb.geo else np.array([""] * len(gid))
    pop = cb.geo["population_total"].to_numpy(float) if "population_total" in cb.geo else np.full(len(gid), np.nan)
    rows = []
    for i in top_idx:
        rows.append({
            "rank": int(rank[i]), "geo_id": gid[i], "geo_name": names[i], "state_abbr": st[i],
            "population_total": _clean(pop[i]), "adults_18plus": _clean(cb.adults[i]),
            "rank_engine_same_weights": int(rank_base[i]) if np.isfinite(rank_base[i]) else None,
            "subgroup_estimate": _clean(est[i]), "subgroup_ci_low": _clean(est_lo[i]),
            "subgroup_ci_high": _clean(est_hi[i]), "subgroup_rate_per_100k": _clean(rate[i]),
            "subgroup_rate_ci_low": _clean(lo[i]), "subgroup_rate_ci_high": _clean(hi[i]),
            "burden_pct": _clean(np.atleast_2d(comps["burden"])[0][i]),
            "vulnerability_pct": _clean(np.atleast_2d(comps["vulnerability"])[0][i]),
            "diagnostic_desert_pct": _clean(np.atleast_2d(comps["diagnostic_desert"])[0][i]),
            "clinic_capacity_pct": _clean(np.atleast_2d(comps["clinic_capacity"])[0][i]),
            "research_readiness_pct": _clean(np.atleast_2d(comps["research_readiness"])[0][i]),
            "measurement_evidence": _clean(np.atleast_2d(comps["measurement_evidence"])[0][i]),
            "expected_yield": _clean(np.atleast_2d(comps["expected_yield"])[0][i]),
            "composite": _clean(comp[i]), "mc_rank_p05": _clean(q05[i]), "mc_rank_p50": _clean(q50[i]),
            "mc_rank_p95": _clean(q95[i]), "mc_p_top10": _clean(p10[i])})
    both = elig & elig_base
    rho = float(pd.Series(comp[both]).corr(pd.Series(base[both]), method="spearman")) if both.sum() > 2 else None
    t_new = set(np.flatnonzero(rank <= top))
    t_base = set(np.flatnonzero(rank_base <= top))
    return {
        "status": "done", "basis": "deployment_composite_with_subgroup_burden", "kind": "deployment_composite",
        "weight_set": ws, "subgroup_condition": sub_cond,
        "replaced_members": [cb.members[i].condition_id for i in replaced],
        "kept_members": [m.condition_id for i, m in enumerate(cb.members) if i not in replaced],
        "burden_replacement": "estimated device-subgroup rate per 100k adults (geography.subgroup: base burden x "
                              "cohort subgroup fraction by age band x sex), percentile-ranked like the burden it "
                              "replaces",
        "burden_basis": sorted(sub["burden_basis"].astype(str).unique()),
        "n_ranked": int(np.isfinite(rank).sum()), "clinic_capacity_basis": cb.cr.get("clinic_capacity_basis"),
        "clinic_capacity_note": CAPACITY_NOTE.get(cb.cr.get("clinic_capacity_basis"), ""),
        "monte_carlo": f"equal-weight Monte Carlo of the re-score ({n_draws} draws: Dirichlet(20 x equal) weights x "
                       "subgroup-rate draws, normal from the 95% interval, independent per county: conservative, since the "
                       "shared cohort fraction moves every county together); 5th-95th percentile of the rank and "
                       f"P(top {top})",
        "effect": {"spearman_vs_engine": rho, "top_overlap": len(t_new & t_base), "top_n": top,
                   "entering": [names[i] for i in sorted(t_new - t_base, key=lambda i: rank[i])],
                   "leaving": [names[i] for i in sorted(t_base - t_new, key=lambda i: rank_base[i])]},
        "counties": rows,
    }


def subgroup_burden_only_ranking(sub: pd.DataFrame, top: int) -> dict:
    """Outside the scored grid with a subgroup: rank counties by the estimated number of people in the subgroup
    (population >= scoring.yaml ranking.min_population), as outreach's geo_subgroup_burden fallback does."""
    from .. import outreach as OR
    geo = OR._county_frame()
    geo = geo[geo["population_total"].fillna(0) >= OR.min_population()]
    s = sub[sub["geo_level"] == "county"][["geo_id", "subgroup_estimate", "subgroup_ci_low", "subgroup_ci_high",
                                           "subgroup_rate_per_100k", "base_burden"]]
    m = geo.merge(s.rename(columns={"geo_id": "county_fips"}), on="county_fips", how="inner")
    m = m.dropna(subset=["subgroup_estimate"]).sort_values(["subgroup_estimate", "county_fips"],
                                                           ascending=[False, True]).reset_index(drop=True)
    m["rank"] = np.arange(1, len(m) + 1)
    m = m[m["rank"] <= top].rename(columns={"county_fips": "geo_id"})
    cols = ["rank", "geo_id", "geo_name", "state_abbr", "population_total", "subgroup_estimate", "subgroup_ci_low",
            "subgroup_ci_high", "subgroup_rate_per_100k", "base_burden"]
    return {"status": "done", "basis": "geo_subgroup_burden", "kind": "burden_only", "weight_set": None,
            "why": "SUBGROUP-BURDEN-ONLY ranking (estimated people in the device subgroup per county); the condition "
                   "is outside the scored grid, so there is no deployment composite",
            "burden_basis": sorted(sub["burden_basis"].astype(str).unique()), "n_ranked": None,
            "monte_carlo": "subgroup estimate 95% interval per county (no rank Monte Carlo)", "counties": _records(m[cols])}


# ------------------------------------------------------------------------------------------------------------------
# contacts per county
# ------------------------------------------------------------------------------------------------------------------
def contact_counts(fips: list[str], measurement: str) -> list[dict]:
    """Billing clinicians, experience sites and candidate facilities per county (counts only)."""
    from ..scoring.opportunity import measurement_spec
    try:
        spec = measurement_spec(measurement)
        groups = set(spec.get("implementer_groups", []))
        mid = spec["measurement_id"]
    except KeyError:
        groups, mid = set(), measurement
    out = {f: {"geo_id": f} for f in fips}
    if table_exists("geo_measurement_capacity"):
        cap = read_table("geo_measurement_capacity")
        cap = cap[(cap["measurement_id"] == mid) & (cap["geo_level"] == "county") & cap["geo_id"].isin(fips)]
        for r in cap.itertuples():
            out[r.geo_id].update({"capacity_basis": r.capacity_basis,
                                  "billing_clinicians_n": _clean(r.active_clinicians_n),
                                  "billing_clinicians_within_50km_n": _clean(r.active_clinicians_within_50km_n),
                                  "billing_clinicians_per_100k_adults": _clean(r.active_clinicians_per_100k_adults),
                                  "billing_data_year": _clean(r.data_year)})
    if table_exists("measurement_experience_sites"):
        e = read_table("measurement_experience_sites", columns=["measurement_id", "county_fips", "facility_id"])
        e = e[(e["measurement_id"] == mid) & e["county_fips"].isin(fips)]
        n = e.groupby("county_fips")["facility_id"].nunique()
        for f in fips:
            out[f]["experience_sites_n"] = int(n.get(f, 0))
    if table_exists("clinic_registry"):
        c = read_table("clinic_registry", columns=["county_fips", "facility_id", "is_clinic_candidate",
                                                   "specialty_groups"])
        c = c[c["county_fips"].isin(fips) & c["is_clinic_candidate"].fillna(False).astype(bool)]
        sg = c["specialty_groups"].fillna("").astype(str).map(lambda s: set(x for x in re.split(r"[|;,]", s) if x))
        hit = sg.map(lambda s: bool(s & groups))
        n = c[hit].groupby("county_fips")["facility_id"].nunique()
        for f in fips:
            out[f]["candidate_facilities_n"] = int(n.get(f, 0))
    return [out[f] for f in fips]


# ------------------------------------------------------------------------------------------------------------------
# the launch ranking
# ------------------------------------------------------------------------------------------------------------------
def launch_ranking(condition: str, measurement: str, top: int = 10, subgroup: pd.DataFrame | None = None,
                   phenotype: dict | None = None, log=print) -> dict:
    """(a) engine ranking, (b) subgroup re-score where a subgroup exists; contacts for the top counties."""
    scored = in_scored_grid(condition, measurement)
    eng = engine_ranking(condition, measurement, top)
    sub_res = None
    if subgroup is not None and len(subgroup):
        if scored:
            sub_res = subgroup_rescore(condition, measurement, subgroup, top)
        else:
            from .. import outreach as OR
            members = OR.resolve_condition(condition)["members"]
            sc = str(subgroup["condition_id"].iloc[0])
            sub_res = (subgroup_burden_only_ranking(subgroup, top) if sc in members else
                       {"status": "skipped", "reason": f"the device subgroup is defined within {sc}, not within "
                                                       f"{condition}"})
        if sub_res.get("status") == "done":
            sub_res["phenotype_id"] = (phenotype or {}).get("phenotype_id")
            sub_res["synthetic"] = bool((phenotype or {}).get("synthetic"))
    primary = "subgroup" if sub_res and sub_res.get("status") == "done" else "engine"
    fips = list(dict.fromkeys([c["geo_id"] for c in eng["counties"]]
                              + [c["geo_id"] for c in ((sub_res or {}).get("counties") or [])]))
    contacts = {r["geo_id"]: r for r in contact_counts(fips, measurement)}
    for block in (eng, sub_res or {}):
        for c in block.get("counties") or []:
            c["contacts"] = {k: v for k, v in contacts.get(c["geo_id"], {}).items() if k != "geo_id"}
    definition = (
        "Greatest impact = the device-subgroup re-score: the engine's composite with the condition burden replaced by "
        "the estimated device-subgroup rate per 100k adults per county"
        + (f" under the {sub_res.get('weight_set')} weight set" if primary == "subgroup" and sub_res.get("weight_set")
           else " (subgroup-burden-only: the condition is outside the scored grid)")
        if primary == "subgroup" else
        ("Greatest impact = the engine's " + (f"{eng['weight_set']} ranking ({eng['why']})" if eng.get("weight_set")
                                              else f"burden-only fallback ({eng['why']})")
         + "; no device subgroup was available to replace the burden"))
    caveats = launch_caveats(condition, measurement, eng, sub_res, phenotype)
    return {"condition": condition, "measurement": measurement, "scored_grid": scored, "top": top,
            "primary": primary, "definition": definition, "engine": eng, "subgroup": sub_res, "caveats": caveats}


def launch_caveats(condition, measurement, eng, sub_res, phenotype) -> list[str]:
    out = []
    if eng.get("kind") == "burden_only" or (sub_res or {}).get("kind") == "burden_only":
        out.append(f"{condition} x {measurement} is outside the scored grid: the rankings are BURDEN-ONLY (no "
                   "vulnerability, desert, capacity or readiness components).")
    srcs = " ".join(str(c.get("burden_source_name") or "") for c in eng.get("counties", []))
    if "small-area" in srcs or condition in ("long_covid", "long_covid_or_me_cfs"):
        out.append("Long COVID county burden is the BRFSS 2023 small-area model: within-state differences come from "
                   "county composition and covariates (modelled, not observed county prevalence).")
    basis = eng.get("clinic_capacity_basis") or (sub_res or {}).get("clinic_capacity_basis")
    if basis == "implementer_density":
        out.append(f"{measurement} has no dedicated billing code (nailfold capillaroscopy has no CPT/HCPCS code): "
                   "clinic_capacity is implementer-group density, and billing-clinician counts use related or proxy "
                   "codes only.")
    out.append("Clinician activity and billing counts cover Medicare fee-for-service only.")
    if (eng.get("performance") or {}).get("status") != "known":
        out.append("No known performance record for this condition x measurement in the rankings: the engine ranking "
                   "is equal-weight and carries no device-accuracy information.")
    if sub_res and sub_res.get("status") == "done":
        out.append("The subgroup re-score applies one cohort's device-positive fraction (by age band x sex) to every "
                   "county: no geographic variation in the device-positive share is modelled.")
        if (phenotype or {}).get("synthetic") or sub_res.get("synthetic"):
            out.append("The device subgroup comes from SYNTHETIC data: its re-score demonstrates the mechanics only.")
        if any("lyme" in b for b in sub_res.get("burden_basis", [])):
            out.append("Lyme burden is reported surveillance cases (under-ascertained; 2022+ case definition); PTLDS "
                       "is derived as 13.7% of reported cases (Aucott 2022).")
    elif sub_res:
        out.append(f"Device-subgroup re-score not computed: {sub_res.get('reason')}")
    return out


# ------------------------------------------------------------------------------------------------------------------
# outreach workbook (written under the run folder)
# ------------------------------------------------------------------------------------------------------------------
def outreach_workbook(condition: str, measurement: str, launch: dict, out_dir, top: int = 10,
                      per_county: int = 20) -> dict:
    """measure_it.outreach.build_outreach for the launch ranking's counties, written under out_dir (never committed).

    Engine ranking on the scored grid: build_outreach's own county selection (same weight set). Otherwise the launch
    ranking's counties are handed to build_outreach (its select_counties is replaced for this call only)."""
    from .. import outreach as OR
    from .orchestrate import local_only
    from pathlib import Path
    out_dir = local_only(Path(out_dir))
    prim = launch["subgroup"] if launch["primary"] == "subgroup" else launch["engine"]
    if prim.get("basis") == "deployment_opportunities":
        return OR.build_outreach(condition, measurement, top, per_county, weight_set=prim["weight_set"],
                                 out_dir=out_dir)
    cty = pd.DataFrame(prim["counties"]).rename(columns={"geo_id": "county_fips"})
    keep = ["rank", "county_fips", "geo_name", "state_abbr", "population_total"] + [
        c for c in ("composite", "subgroup_estimate", "subgroup_rate_per_100k", "burden_pct", "vulnerability_pct",
                    "diagnostic_desert_pct", "clinic_capacity_pct", "research_readiness_pct", "burden_score",
                    "mc_rank_p05", "mc_rank_p95") if c in cty.columns]
    cty = cty[keep]
    label = ("LAUNCH ranking (measure-it launch): " + launch["definition"]
             + (f"; phenotype {prim.get('phenotype_id')}" if prim.get("phenotype_id") else "")
             + (" (SYNTHETIC phenotype: test data, not a real estimate)" if prim.get("synthetic") else ""))
    basis = {"basis": "launch_" + str(prim.get("basis")), "description": label}
    with mock.patch.object(OR, "select_counties", lambda *a, **k: (cty.copy(), dict(basis))):
        return OR.build_outreach(condition, measurement, top, per_county, out_dir=out_dir)
