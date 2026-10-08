"""Deployment recommendations (SPEC Phases 5-6 and the deployment-opportunity output schema).

    rank_deployment_opportunities(condition, measurement, geography_level='county', top_n=10,
                                  weight_set='equal') -> dict

Returns, per ranked region, the SPEC output schema:
    condition, phenotype, measurement, technology {name, regulatory_context, technology_evidence},
    geography {name, fips, burden, burden_evidence_level, vulnerability, diagnostic_desert, ...components},
    candidate_sites, research_evidence, molecular_context, uncertainties, provenance, recommended_next_step
Every region is a CANDIDATE deployment opportunity for pilot evaluation, never a validated diagnostic pathway.
Missing data come back as config.UNKNOWN with a reason; nothing is filled in.

Aliases: 'Long COVID or ME/CFS' -> condition set long_covid_or_me_cfs; any 'A or B' of resolvable conditions ->
an ad hoc set scored with the same engine; 'wearable autonomic monitoring' -> the wearable bundle.

`run()` writes the SPEC demo query: processed table `deployment_candidates` (one row per recommendation incl.
its JSON), results/example_deployment_recommendation.json, results/tables/demo_top10_*.csv and draft maps.
"""
from __future__ import annotations

import json
import re
import time
from functools import lru_cache

import numpy as np
import pandas as pd

from ..config import MAPS, RESULTS, TABLES, UNKNOWN, utc_now_iso
from ..facilities import matching as M
from ..facilities.registry import SEP
from ..ontology.normalize import norm_text
from ..provenance import add_provenance
from ..store import processed_path, read_table, table_exists, write_table
from . import metric_link as ML
from . import opportunity as O

PRODUCER = "measure_it.scoring.recommend"
DEMO_QUERY = ("Find 10 U.S. regions where wearable monitoring of autonomic/activity abnormalities in Long COVID or "
              "ME/CFS would be useful to evaluate, and identify clinics/research sites that could plausibly "
              "participate.")
DEMO_ARGS = {"condition": "Long COVID or ME/CFS", "measurement": "wearable autonomic monitoring",
             "geography_level": "county", "top_n": 10, "weight_set": "equal"}
# the primary alternative (docs/ANALYSIS_PLAN_METRIC_LINK.md): the same query ranked with the measurement's own evidence
DEMO_ARGS_EVIDENCE = {**DEMO_ARGS, "weight_set": "evidence_weighted"}
MAX_IDS = 25
MOLECULAR_LABEL = ("condition-level molecular enrichment (public databases; not participant-linked; not patient "
                   "multi-omics; not mechanism)")
GENERAL_UNCERTAINTIES = [
    "Ecological ranking: burden, vulnerability, providers, trials and grants describe places and registries, not "
    "people; no individual's location or condition is inferred.",
    "Curated assumptions: condition specialties and measurement implementer groups (configs/relevance.yaml); an NPPES "
    "taxonomy does not mean a clinician evaluates or treats the condition.",
    "The composite re-uses information: the diagnostic desert contains burden and vulnerability, and its access "
    "terms (few providers/trials) oppose clinic_capacity and research_readiness.",
    "Registered trials and NIH grants are research-activity signals, not evidence that the measurement works; an FDA "
    "record is a deployment-readiness signal, not evidence that a device detects the condition.",
    "Facility and provider locations are ZIP/ZCTA or city centroids; the 50 km radius is measured from the county's "
    "Census internal point, a weak proxy for where people live in large counties.",
]


# ------------------------------------------------------------------------------------------------------------------
# JSON helpers
# ------------------------------------------------------------------------------------------------------------------

def _py(v):
    if isinstance(v, dict):
        return {str(k): _py(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, set)):
        return [_py(x) for x in v]
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.ndarray):
        return [_py(x) for x in v.tolist()]
    if v is pd.NA or v is pd.NaT:
        return None
    return v


def _num(v, nd=4):
    try:
        if v is None or pd.isna(v):
            return None
    except (TypeError, ValueError):
        return v
    return round(float(v), nd)


def _or_unknown(v):
    return UNKNOWN if v is None or (isinstance(v, float) and not np.isfinite(v)) else v


# ------------------------------------------------------------------------------------------------------------------
# resolution
# ------------------------------------------------------------------------------------------------------------------

def resolve_condition_query(text: str) -> dict:
    """Condition id / alias / 'A or B' -> condition or condition-set id scored by the engine."""
    q = str(text or "").strip()
    if not q:
        return {"status": UNKNOWN, "query": q, "reason": "empty condition"}
    sets = O.condition_sets()
    if q in sets:
        return {"status": "matched", "query": q, "condition_id": q, "match_reason": "condition-set id"}
    nq = norm_text(q)
    for sid, s in sets.items():
        if nq in {norm_text(a) for a in s["aliases"]} | {norm_text(s["label"])}:
            return {"status": "matched", "query": q, "condition_id": sid, "match_reason": "condition-set alias"}
    parts = [p for p in re.split(r"\s+or\s+|\s+and\s+|\s*\+\s*|\s*;\s*|\s*,\s*", q, flags=re.I) if p.strip()]
    if len(parts) > 1:
        ids, bad = [], []
        for p in parts:
            r = M.resolve_condition(p)
            (ids.append(r["condition_id"]) if r["status"] == "matched" else bad.append({"part": p, **r}))
        if bad:
            return {"status": UNKNOWN, "query": q, "reason": "a part of the condition query did not resolve",
                    "unresolved": bad}
        ids = sorted(set(ids))
        for sid, s in sets.items():
            if set(s["members"]) == set(ids):
                return {"status": "matched", "query": q, "condition_id": sid,
                        "match_reason": f"parts resolved to {ids} = members of condition set {sid}"}
        if len(ids) == 1:
            return {"status": "matched", "query": q, "condition_id": ids[0], "match_reason": "parts resolved to one id"}
        return {"status": "matched", "query": q, "condition_id": "+".join(ids),
                "match_reason": f"ad hoc condition set of {ids} (same aggregation as the demo cluster)"}
    r = M.resolve_condition(q)
    if r["status"] != "matched":
        return {"status": UNKNOWN, "query": q, "reason": r.get("reason"), "candidates": r.get("candidates", [])}
    return {"status": "matched", "query": q, "condition_id": r["condition_id"], "match_reason": r.get("match_reason")}


# ------------------------------------------------------------------------------------------------------------------
# context (condition-level; computed once per query)
# ------------------------------------------------------------------------------------------------------------------

@lru_cache(maxsize=32)
def phenotype_context(cid: str) -> dict:
    from ..wearables.signatures import get_patient_phenotype_signature
    try:
        r = get_patient_phenotype_signature(cid)
    except Exception as e:  # pragma: no cover
        return {"condition_id": cid, "status": UNKNOWN, "reason": f"{type(e).__name__}: {e}"}
    if r.get("status") != "found":
        return {"condition_id": cid, "status": UNKNOWN, "reason": r.get("reason", "no signature")}
    ds = []
    for d in r.get("datasets", []):
        ds.append({
            "dataset_id": d.get("dataset_id"), "phenotype_id": d.get("phenotype_id"),
            "phenotype_label": d.get("phenotype_label"), "evidence_level": d.get("evidence_level"),
            "is_proxy": d.get("is_proxy"), "n_cases": d.get("n_cases"), "n_controls": d.get("n_controls"),
            "n_features_tested": d.get("n_features_tested"),
            "n_features_fdr_significant": d.get("n_features_fdr_significant"),
            "top_features": [{k: f.get(k) for k in ("feature", "feature_label", "effect_measure", "effect_size",
                                                     "ci_low", "ci_high", "q_value_bh", "object_id") if k in f}
                             for f in (d.get("top_features") or [])[:3]],
            "null_results": (d.get("null_results") or {}).get("summary"),
            "model_level_results": [{k: m.get(k) for k in ("feature", "effect_measure", "effect_size", "ci_low",
                                                          "ci_high", "object_id")}
                                    for m in ((d.get("model_level_results") or {}).get("rows") or [])[:3]],
        })
    return {"condition_id": cid, "status": "found", "datasets": ds, "interpretation": r.get("interpretation"),
            "note": "Public person-level cohorts; group differences, not diagnostic accuracy, and never an individual "
                    "finding. These participants are not the people of any ranked region."}


@lru_cache(maxsize=32)
def molecular_context(cid: str) -> dict:
    from ..omics.query import get_molecular_context
    try:
        r = get_molecular_context(cid, top_n=5)
    except Exception as e:  # pragma: no cover
        return {"condition_id": cid, "label": MOLECULAR_LABEL, "status": UNKNOWN, "reason": f"{type(e).__name__}: {e}"}
    if r.get("status") != "available":
        return {"condition_id": cid, "label": MOLECULAR_LABEL, "status": UNKNOWN,
                "reason": r.get("reason", r.get("status"))}
    ev = {}
    for etype, v in (r.get("evidence") or {}).items():
        tops = []
        for ent_type, lst in (v.get("top_entities") or {}).items():
            for e in lst[:3]:
                tops.append({"entity_type": ent_type, "entity_id": e.get("entity_id"),
                             "entity_label": e.get("entity_label"), "sources": e.get("sources")})
        ev[etype] = {"n_records": v.get("n_records"), "n_distinct_entities": v.get("n_distinct_entities"),
                     "counts_by_source": v.get("counts_by_source"), "top_entities": tops[:5]}
    return {"condition_id": cid, "label": MOLECULAR_LABEL, "layer": r.get("layer"), "status": "available",
            "evidence_summary": ev, "cross_source_genetic_agreement": r.get("cross_source_genetic_agreement"),
            "caveats": r.get("caveats"),
            "provenance": [{k: p.get(k) for k in ("source_name", "source_version", "retrieved_at")}
                           for p in (r.get("provenance") or [])]}


def _signal_verdict(ps):
    """Verdict of a phenotype_signal_strength dimension (dict with 'value' / 'status', or a plain string)."""
    if isinstance(ps, dict):
        v = ps.get("value", ps.get("tier"))
        return v if v is not None else ps.get("status", UNKNOWN)
    return ps if ps is not None else UNKNOWN


@lru_cache(maxsize=32)
def technology_context(mid: str, members: tuple, cond_ids: tuple) -> dict:
    from ..measurements.query import get_measurement_evidence
    from ..measurements.regulatory import get_regulatory_context
    reg = []
    for m in members:
        r = get_regulatory_context(m)
        if r.get("status") != "matched":
            reg.append({"measurement_id": m, "status": UNKNOWN, "reason": r.get("reason")})
            continue
        s = next((x for x in r.get("summary", []) if x["measurement_id"] == m), None)
        reg.append({"measurement_id": m, "status": "matched", "summary": s,
                    "example_fda_records": [{k: d.get(k) for k in ("clearance_id", "device_name", "applicant",
                                                                  "decision_date", "product_code")}
                                            for d in (r.get("matching_fda_records") or [])[:3]],
                    "regulatory_note": r.get("regulatory_note")})
    tev = []
    for cid in cond_ids:
        r = get_measurement_evidence(mid, cid)
        if r.get("status") != "ok":
            tev.append({"condition_id": cid, "status": UNKNOWN, "reason": r.get("reason", r.get("status"))})
            continue
        for e in r.get("evidence", []):
            d = e.get("dimensions", {})
            tev.append({"condition_id": cid, "measurement_id": e.get("measurement_id"),
                        "evidence_object_id": e.get("evidence_object_id"),
                        "measurement_evidence_strength": {k: (d.get("measurement_evidence_strength") or {}).get(k)
                                                          for k in ("tier", "n_trials_objective",
                                                                    "n_trials_outcome_measure", "n_nih_core_projects")},
                        "technology_maturity": (d.get("technology_maturity") or {}).get("tier"),
                        "regulatory_visibility": (d.get("regulatory_visibility") or {}).get("tier"),
                        "deployment_complexity": (d.get("deployment_complexity") or {}).get("value"),
                        # the dimension is a dict whose verdict is 'value' (null_result, supported, ... or UNKNOWN);
                        # reading a 'tier' key here returned None for every row
                        "phenotype_signal_strength": _signal_verdict(d.get("phenotype_signal_strength")),
                        "example_trial_object_ids": (e.get("trial_object_ids") or [])[:5],
                        "example_grant_object_ids": (e.get("grant_object_ids") or [])[:5]})
    return {"regulatory_context": reg, "technology_evidence": tev,
            "guardrail": "Trials/grants show what researchers deploy, not that a measurement works; FDA records are "
                         "deployment-readiness signals, not evidence that a device diagnoses the target illness."}


PERF_KEYS = ("object_id", "performance_status", "quality_tier", "quality_tier_label", "selected_record_id",
             "source_kind", "dataset_id", "analysis", "label_basis", "comparator", "comparator_kind", "n_cases",
             "n_controls", "auroc", "auroc_ci_low", "auroc_ci_high", "auroc_basis", "op_sensitivity",
             "op_sensitivity_ci_low", "op_sensitivity_ci_high", "op_specificity", "op_specificity_ci_low",
             "op_specificity_ci_high", "op_basis", "measurement_evidence", "measurement_evidence_undiscounted",
             "evidence_conflict", "evidence_conflict_note", "performance_note")


def performance_block(cid: str, mid: str) -> dict:
    """The scored measurement_performance record of (condition or set, measurement): the metric -> translation link
    (docs/ANALYSIS_PLAN_METRIC_LINK.md). UNKNOWN stays UNKNOWN; tier 3 rows are published claims, not reproduced."""
    try:
        p = ML.performance_for(cid, mid)
    except KeyError as e:
        return {"performance_status": UNKNOWN, "reason": str(e)}
    out = {}
    for k in PERF_KEYS:
        v = p.get(k)
        out[k] = _num(v) if isinstance(v, (float, np.floating)) else (_py(v) if v is not None else None)
    for k in ("source_object_ids", "caveats"):
        try:
            out[k] = json.loads(p.get(k) or "[]")
        except (TypeError, ValueError):
            out[k] = []
    if out.get("performance_status") != "known":
        out["performance_status"] = p.get("performance_status") or UNKNOWN
    out["target_specificity"] = ML.TARGET_SPECIFICITY
    out["tier_text"] = ML.tier_text(p.get("quality_tier") or 4)
    out["healthy_control_caveat"] = ML.HEALTHY_CAVEAT
    out["expected_yield_language"] = ML.YIELD_LANGUAGE
    return out


# ------------------------------------------------------------------------------------------------------------------
# region-level evidence
# ------------------------------------------------------------------------------------------------------------------

def research_evidence(cb: O.Combo, i: int) -> dict:
    """trial:/nih: object ids at facilities in the region's pool (county or within 50 km; state: in state)."""
    S = O.spatial()
    b = S.base
    L = S.levels[cb.level]
    locs = L.pool[i].indices
    fac = np.flatnonzero(np.isin(S.fac_loc, locs))
    cond = list(cb.cond["members"])
    ft = b.ft[b.ft["fidx"].isin(fac)]
    tc = ft[M._has(ft["condition_ids_literal"], cond)]
    tt = ft[M._has(ft["measurement_classes"], list(cb.meas["members"]))]
    fn = b.fn[b.fn["fidx"].isin(fac)]
    fn = fn[M._has(fn["condition_ids_precision"], cond)]
    trials = sorted(tc["nct_id"].unique())
    tech = sorted(tt["nct_id"].unique())
    nih = fn.sort_values(["fiscal_year", "appl_id"], ascending=[False, True])
    return {
        "pool": "facilities in the county or within 50 km of its internal point" if cb.level == "county"
                else "facilities in the state",
        "condition_trials": {"n": len(trials), "object_ids": ["trial:" + x for x in trials[:MAX_IDS]]},
        "technology_experience_trials": {"n": len(tech), "object_ids": ["trial:" + x for x in tech[:MAX_IDS]],
                                         "note": "trials (any target condition) whose registered text mentions a "
                                                 "member measurement class"},
        "nih_projects": {"n_core_projects": int(nih["core_project_num"].nunique()),
                         "n_appl_ids": int(nih["appl_id"].nunique()),
                         "object_ids": ["nih:" + str(x) for x in nih["appl_id"].drop_duplicates().head(MAX_IDS)],
                         "examples": [{"object_id": f"nih:{r.appl_id}", "core_project_num": r.core_project_num,
                                       "fiscal_year": _py(r.fiscal_year), "title": r.project_title}
                                      for r in nih.drop_duplicates("core_project_num").head(5).itertuples()]},
        "note": "Registration and funding are research-activity signals, not evidence that a measurement works.",
    }


def candidate_sites(geo_id: str, cond_ids: list[str], mid: str) -> tuple[list[dict], list[dict]]:
    """find_candidate_clinics per member condition, merged by facility (member conditions recorded)."""
    merged: dict[str, dict] = {}
    status = []
    for cid in cond_ids:
        r = M.find_candidate_clinics(geo_id, cid, mid)
        status.append({"condition_id": cid, "status": r.get("status"), "reason": r.get("reason"),
                       "n_facilities_in_pool": r.get("n_facilities_in_pool"), "n_eligible": r.get("n_eligible"),
                       "characteristics_absent_in_pool": r.get("characteristics_absent_in_pool"),
                       "absence_note": r.get("absence_note")})
        for c in r.get("candidates", []):
            fid = c["facility_id"]
            if fid not in merged:
                merged[fid] = {k: c.get(k) for k in ("facility_id", "object_id", "facility_name", "primary_kind",
                                                     "city", "state", "zip5", "county_fips", "lat", "lon",
                                                     "geocode_precision", "distance_km", "inside_geography")}
                merged[fid].update({"selected_for_conditions": [], "selected_via": {}, "characteristics": {},
                                    "reasons": [], "source_object_ids": {}, "framing": c["framing"]})
            m = merged[fid]
            m["selected_for_conditions"].append(cid)
            m["selected_via"][cid] = c["selected_via"]
            m["characteristics"][cid] = {k: {"value": v["value"], "rank": v["rank"], "n_ranked": v["n_ranked"]}
                                         for k, v in c["characteristics"].items()}
            for rs in c["reasons"]:
                if rs not in m["reasons"]:
                    m["reasons"].append(rs)
            for k, v in c["source_object_ids"].items():
                if isinstance(v, list):
                    m["source_object_ids"].setdefault(k, [])
                    m["source_object_ids"][k] = sorted(set(m["source_object_ids"][k]) | set(v))
    return list(merged.values()), status


# ------------------------------------------------------------------------------------------------------------------
# main query
# ------------------------------------------------------------------------------------------------------------------

def _unknown(reason: str, **extra) -> dict:
    return {"status": UNKNOWN, "reason": reason, "recommendations": [], **extra}


def _rows(cid: str, mid: str, level: str) -> tuple[pd.DataFrame, str]:
    if table_exists(O.TABLE):
        import pyarrow.parquet as pq
        t = pq.read_table(processed_path(O.TABLE), filters=[("condition_id", "==", cid), ("measurement_id", "==", mid),
                                                            ("geo_level", "==", level)])
        df = t.to_pandas()
        if len(df):
            return df, "deployment_opportunities (precomputed)"
    return O.score_combo(cid, mid, level), "computed on the fly with measure_it.scoring.opportunity (same engine)"


def rank_deployment_opportunities(condition: str, measurement: str, geography_level: str = "county",
                                  top_n: int = 10, weight_set: str = "equal", include_small_population: bool = False,
                                  with_context: bool = True) -> dict:
    """Top-N candidate deployment opportunities with the SPEC output schema per region."""
    t0 = time.time()
    query = {"condition": condition, "measurement": measurement, "geography_level": geography_level,
             "top_n": top_n, "weight_set": weight_set}
    if geography_level not in O.LEVELS:
        return _unknown(f"geography_level must be one of {O.LEVELS}", query=query)
    try:
        top_ok = int(top_n) == top_n and int(top_n) >= 1
    except (TypeError, ValueError):
        top_ok = False
    if not top_ok:
        return _unknown("top_n must be a positive integer", query=query)
    top_n = int(top_n)
    ws = O.weight_sets()
    if weight_set not in ws:
        return _unknown(f"unknown weight_set; available: {sorted(ws)}", query=query)
    rc = resolve_condition_query(condition)
    if rc["status"] != "matched":
        return _unknown(f"condition: {rc.get('reason')}", query=query, condition_resolution=rc)
    rm = M.resolve_measurement(measurement)
    if rm["status"] != "matched":
        return _unknown(f"measurement: {rm.get('reason')}", query=query, measurement_resolution=rm)
    cid, mid = rc["condition_id"], rm["measurement_id"]
    try:
        spec = O.condition_spec(cid)
    except KeyError as e:
        return _unknown(str(e), query=query)
    perf = performance_block(cid, mid)
    if O.uses_evidence(ws[weight_set]) and perf.get("performance_status") != "known":
        return _unknown(f"not ranked under '{weight_set}': measurement performance "
                        f"{perf.get('performance_status')} for {cid} x {mid} (docs/ANALYSIS_PLAN_METRIC_LINK.md "
                        "section 5: an UNKNOWN or partial performance record is never scored 0 or 1; use weight_set "
                        "'equal', whose rank carries no measurement-performance information)", query=query,
                        condition_resolution=rc, measurement_performance=perf)
    df, basis = _rows(cid, mid, geography_level)
    rank_col = f"rank_{weight_set}"
    if include_small_population and geography_level == "county":
        # every scored county, same composite, ties by geo_id (for 'equal' this equals rank_equal_incl_small)
        rank_col = f"rank_{weight_set}_incl_small"
        df = df.sort_values("geo_id").reset_index(drop=True)
        df[rank_col] = O.rank_rows(df[f"composite_{weight_set}"].to_numpy(float), np.ones(len(df), bool))[0]
    ranked = df[df[rank_col].notna()].sort_values(rank_col)
    n_ranked = int(len(ranked))
    if n_ranked == 0:
        return _unknown(f"no region has a defined {rank_col} (e.g. burden_only for a level-D condition)",
                        query=query, condition_resolution=rc)
    top = ranked.head(int(top_n))
    cb = O.make_combo(cid, mid, geography_level) if with_context else None
    gi = {g: i for i, g in enumerate(cb.geo["geo_id"])} if cb is not None else {}
    meas = O.measurement_spec(mid)
    ctx = {}
    if with_context:
        ctx["phenotype"] = [phenotype_context(c) for c in spec["members"]]
        ctx["molecular"] = [molecular_context(c) for c in spec["members"]]
        ctx["technology"] = dict(technology_context(mid, tuple(meas["members"]), tuple(spec["members"])))
        ctx["technology"]["measurement_performance"] = perf
        ctx["technology"]["member_condition_performance"] = (
            [performance_block(c, mid) for c in spec["members"]] if len(spec["members"]) > 1 else [])
    stored = basis.startswith("deployment_opportunities")
    member_rows = stored_member_rows(spec["members"], mid, geography_level) if not stored else set()
    recs = [recommendation(r, spec, meas, weight_set, rank_col, n_ranked, ctx, cb, gi.get(r["geo_id"]),
                           stored=stored, member_rows=member_rows)
            for r in top.to_dict("records")]
    excluded = excluded_incomplete(df, geography_level, include_small_population, stored=stored)
    meta = json.loads((processed_path(O.TABLE).with_suffix(".meta.json")).read_text()) if table_exists(O.TABLE) else {}
    return _py({
        "status": "ok", "query": query, "condition_resolution": rc,
        "measurement_resolution": {k: rm[k] for k in ("measurement_id", "label", "kind", "members",
                                                      "implementer_groups", "match_reason")},
        "weight_set": weight_set, "weights": ws[weight_set], "rank_column": rank_col, "n_regions_ranked": n_ranked,
        "ranking_basis": basis,
        "opportunity_rows_stored": stored,
        "ranking_universe": ("counties of the 50 states + DC with population >= "
                             f"{O.min_population():,}" if geography_level == "county" and not include_small_population
                             else f"all {geography_level} geographies of the 50 states + DC"),
        "method": "docs/SCORING.md; composite = weighted mean of percentile-normalised components; rank interval = "
                  "5th-95th percentile over Monte Carlo draws of burden (CI x evidence-level multiplier) and weights "
                  "(Dirichlet around the default)",
        "framing": O.GUARDRAIL, "general_uncertainties": GENERAL_UNCERTAINTIES,
        "deployment_opportunities_version": meta.get("created_at", basis),
        "recommendations": recs,
        # after the recommendations, so the envelope's first provenance id stays the stored opportunity row
        "measurement_performance": perf,
        "excluded_incomplete_burden": excluded,
        "excluded_incomplete_burden_note": EXCLUDED_NOTE,
        "runtime_s": round(time.time() - t0, 1),
    })


ON_THE_FLY_NOTE = ("computed on the fly with measure_it.scoring.opportunity (same engine); not a stored row, so "
                   "trace_evidence cannot resolve this opportunity id. Trace the component ids in provenance.object_ids "
                   "instead: geo:<fips> (burden, vulnerability and desert inputs per member condition), "
                   "condition:<member>, measurement:<id> and, where stored, each member condition's own opportunity row.")


def stored_member_rows(members: list[str], mid: str, level: str) -> set[str]:
    """Object ids of the member conditions' stored opportunity rows (for an on-the-fly set or measurement)."""
    if not table_exists(O.TABLE) or len(members) < 2:
        return set()
    import pyarrow.parquet as pq
    t = pq.read_table(processed_path(O.TABLE), columns=["object_id"],
                      filters=[("condition_id", "in", list(members)), ("measurement_id", "==", mid),
                               ("geo_level", "==", level)])
    return set(t.column("object_id").to_pylist())


EXCLUDED_NOTE = ("Regions not ranked because their burden is incomplete: a member condition with a defined burden "
                 "measure at this level has no value there (e.g. CMS does not publish the 2024 Connecticut planning "
                 "regions). Their burden is not computed from the other member(s) alone, which would rank a different "
                 "quantity (docs/ANALYSIS_PLAN_SCORING.md deviation 7). They may still be deployment opportunities; "
                 "this ranking cannot place them. Their non-burden components are shown.")


def excluded_incomplete(df: pd.DataFrame, level: str, include_small_population: bool = False,
                        stored: bool = True) -> list[dict]:
    """Incomplete-burden regions of the ranked combination (listed separately, never ranked). object_id is the stored
    opportunity row; None for an on-the-fly combination (no stored row; trace geo_object_id instead)."""
    if "burden_incomplete" not in df.columns:
        return []
    inc = df[df["burden_incomplete"].fillna(False).astype(bool)]
    if level == "county" and not include_small_population:
        inc = inc[~inc["small_population_flag"].fillna(False).astype(bool)]
    out = []
    for r in inc.sort_values("geo_id").to_dict("records"):
        out.append({
            "name": r["geo_name"], "fips": r["geo_id"], "state": r["state_abbr"],
            "object_id": r["object_id"] if stored else None,
            **({} if stored else {"object_id_note": "computed on the fly; not a stored row"}),
            "geo_object_id": r["geo_object_id"], "population": _num(r.get("population_total"), 0),
            "small_population": bool(r.get("small_population_flag")),
            "members_without_value": [x for x in str(r.get("burden_members_missing") or "").split(SEP) if x],
            "reason": r.get("burden_incomplete_reason"),
            "components_kept": {k: _num(r.get(k)) for k in ("vulnerability_pct", "diagnostic_desert_pct",
                                                            "clinic_capacity_pct", "research_readiness_pct")},
            "burden": UNKNOWN, "composite": UNKNOWN, "rank": UNKNOWN})
    return out


def next_step_text(geo_block: dict, weight_set: str, rank_col: str, n_ranked: int, meas: dict, spec: dict,
                   geo_name: str, site_names: list[str]) -> str:
    """recommended_next_step: the rank names its weight set and universe, and the Monte Carlo interval says that it is
    drawn around the EQUAL weights over the population-eligible regions (docs/SCORING.md section 4), so a rank under
    another weight set, or among counties incl. small ones, is never presented as if it were the interval's rank."""
    rank = geo_block["rank"]
    lo, hi = geo_block["rank_interval_5_95"]
    incl_small = rank_col.endswith("_incl_small")
    universe = " among all counties incl. small ones" if incl_small else ""
    rank_txt = (f"rank {rank:.0f} of {n_ranked}{universe} under the '{weight_set}' weights" if rank is not None
                else "unranked")
    if lo is None:
        mc = "no Monte Carlo rank interval (the region is outside the Monte Carlo ranking universe)"
    elif weight_set == "equal" and not incl_small:
        mc = f"5th-95th percentile rank interval {lo:.0f}-{hi:.0f} over Monte Carlo draws of burden and weights"
    elif weight_set == O.EVIDENCE_WEIGHT_SET and not incl_small and \
            str(geo_block.get("monte_carlo_basis", "")).startswith("rank_interval_5_95 and p_top10_monte_carlo are "
                                                                     "drawn around the EVIDENCE_WEIGHTED"):
        mc = (f"5th-95th percentile rank interval {lo:.0f}-{hi:.0f} over Monte Carlo draws of burden, the performance "
              "record's sensitivity / specificity and the evidence-weighted weights")
    else:
        eq = geo_block["ranks_under_weight_sets"].get("equal")
        target = f"the '{weight_set}'-weight rank" + (" among all counties incl. small ones" if incl_small else "")
        mc = (f"the Monte Carlo 5th-95th percentile rank interval {lo:.0f}-{hi:.0f} belongs to the equal-weight rank "
              f"among population-eligible regions ({'UNKNOWN' if eq is None else f'{eq:.0f}'}), drawn around the "
              f"equal weights; it is not an interval for {target}")
    partners = '; '.join(site_names) if site_names else 'UNKNOWN / NOT AVAILABLE (no candidate facility returned)'
    if rank is None and lo is None:
        return f"Candidate deployment opportunity ({rank_txt}); pilot evaluation only; not a validated diagnostic pathway."
    ey = geo_block.get("expected_yield") or {}
    perf_txt = ""
    if ey.get("performance_status") == "known":
        if ey.get("basis") == "count":
            perf_txt = (f" As a planning estimate (not a prediction of diagnoses), screening the reached adults at the "
                        f"scored operating point would flag about {_fmtn(ey.get('expected_detectable_cases'))} "
                        f"detectable cases and {_fmtn(ey.get('expected_false_positives'))} false positives (PPV "
                        f"{ey.get('expected_ppv')}), with performance measured against healthy controls, which "
                        "overstates real-world performance.")
        elif ey.get("basis") == "percentile":
            perf_txt = (" No defensible case count exists at this resolution, so the expected yield is a relative "
                        f"index ({ey.get('expected_yield_index')}); false positives could reach "
                        f"{_fmtn(ey.get('expected_false_positives_upper'))} among the reached adults (upper bound), a "
                        "reminder that where prevalence is low most flags are false positives.")
    else:
        perf_txt = (" The measurement's performance for this condition is " + str(ey.get("performance_status") or
                                                                              UNKNOWN)
                    + " in this project's evidence, so its expected yield is UNKNOWN; this rank carries no "
                      "measurement-performance information.")
    return (f"Candidate deployment opportunity ({rank_txt}; {mc}). A pilot evaluation of {meas['label'].lower()} for "
            f"{spec['label']} in {geo_name} could test feasibility and whether the measurement captures a measurable "
            f"phenotype locally, working with candidate partners such as {partners}, which have characteristics "
            "suggesting they may be viable implementation or study partners." + perf_txt + " This is a candidate "
            "deployment opportunity, not a validated diagnostic pathway; check the listed uncertainties first.")


def _fmtn(v) -> str:
    return f"{v:,.0f}" if isinstance(v, (int, float)) and np.isfinite(v) else str(v)


def expected_yield_block(r: dict) -> dict:
    """Expected yield and false positives of the region at the scored operating point (every factor shown)."""
    basis = r.get("expected_yield_basis")
    st = r.get("measurement_performance_status")
    out = {"basis": basis, "basis_note": r.get("expected_yield_basis_note"), "performance_status": st,
           "operating_point": {"sensitivity": _num(r.get("op_sensitivity"), 3),
                               "sensitivity_ci_95": [_num(r.get("op_sensitivity_ci_low"), 3),
                                                     _num(r.get("op_sensitivity_ci_high"), 3)],
                               "specificity": _num(r.get("op_specificity"), 3),
                               "target_specificity": ML.TARGET_SPECIFICITY},
           "reach": _num(r.get("reach"), 3), "reach_radius_km": _num(r.get("reach_radius_km"), 0),
           "adults_18plus": _num(r.get("adults_18plus"), 0), "reached_adults": _num(r.get("reached_adults"), 0),
           "expected_yield_component": _num(r.get("expected_yield"))}
    if st != "known":
        out["expected_detectable_cases"] = UNKNOWN
        out["expected_false_positives"] = UNKNOWN
        out["note"] = ("measurement performance " + str(st) + ": no expected yield or false positives (nothing is "
                       "imputed from another condition or measurement)")
    elif basis == "count":
        out.update({"burden_count": _num(r.get("burden_count"), 0),
                    "expected_detectable_cases": _num(r.get("expected_detectable_cases"), 0),
                    "expected_detectable_cases_interval_5_95": [_num(r.get("expected_yield_mc_p05"), 0),
                                                                _num(r.get("expected_yield_mc_p95"), 0)],
                    "expected_false_positives": _num(r.get("expected_false_positives"), 0),
                    "expected_false_positives_interval_5_95": [_num(r.get("expected_false_positives_mc_p05"), 0),
                                                               _num(r.get("expected_false_positives_mc_p95"), 0)],
                    "expected_ppv": _num(r.get("expected_ppv"), 3),
                    "false_positives_per_detected_case": _num(r.get("false_positives_per_detected_case"), 2)})
    elif basis == "percentile":
        out.update({"expected_detectable_cases": "UNKNOWN / NOT AVAILABLE (no defensible case count at this "
                                                 "resolution; prevalence-percentile path)",
                    "expected_yield_index": _num(r.get("expected_yield_index")),
                    "expected_yield_index_interval_5_95": [_num(r.get("expected_yield_mc_p05")),
                                                           _num(r.get("expected_yield_mc_p95"))],
                    "expected_false_positives_upper": _num(r.get("expected_false_positives_upper"), 0),
                    "expected_false_positives_upper_interval_5_95": [_num(r.get("expected_false_positives_mc_p05"), 0),
                                                                     _num(r.get("expected_false_positives_mc_p95"), 0)]})
    else:
        out.update({"expected_detectable_cases": UNKNOWN, "expected_false_positives_upper":
                    _num(r.get("expected_false_positives_upper"), 0)})
    out["language"] = ML.YIELD_LANGUAGE + " " + ML.HEALTHY_CAVEAT
    return out


def recommendation(r: dict, spec: dict, meas: dict, weight_set: str, rank_col: str, n_ranked: int, ctx: dict,
                   cb: O.Combo | None, i: int | None, stored: bool = True, member_rows: set | None = None) -> dict:
    members = spec["members"]
    unc = json.loads(r["uncertainties"]) if isinstance(r.get("uncertainties"), str) else []
    mc = json.loads(r["member_components"]) if isinstance(r.get("member_components"), str) else None
    for comp in (mc or {}).values():   # rows built before opportunity.py stopped writing str(NaN): no 'nan' strings
        if isinstance(comp, dict):
            if comp.get("burden_source_resolution") == "nan":
                comp["burden_source_resolution"] = UNKNOWN
            if comp.get("burden_measure_id") == "nan":
                comp["burden_measure_id"] = None
    if len(members) == 1:
        burden = {"value": _or_unknown(_num(r.get("burden_value"))), "unit": r.get("burden_value_unit"),
                  "measure_id": r.get("burden_measure_id"), "measure_label": r.get("burden_measure_label"),
                  "ci_95": [_num(r.get("burden_ci_low")), _num(r.get("burden_ci_high"))],
                  "period": r.get("burden_period"), "percentile": _num(r.get("burden_pct")),
                  "inherited": bool(r.get("burden_inherited")),
                  "source_resolution": r.get("burden_source_resolution")}
    else:
        burden = {"percentile": _num(r.get("burden_pct")),
                  "member_mean_percentile": _num(r.get("burden_set_member_mean_pct")), "members": mc,
                  "aggregation": "percentile rank (re-normalised over the complete regions) of the mean of the burden "
                                 "percentiles of the members with a defined (non-D) burden measure; computed only "
                                 "where every such member has a value (regions where one is missing are not ranked)",
                  "members_without_usable_burden": [x for x in str(r.get("burden_members_excluded_level_D") or "")
                                                    .split(SEP) if x],
                  "members_without_usable_burden_note": "no burden measure is defined for these members at this "
                                                        "level (level D); they are not part of the set burden",
                  "inherited": bool(r.get("burden_inherited"))}
    srcs = [str(r.get("burden_source_name") or "")] + [str(c.get("burden_source_name") or "")
                                                       for c in (mc or {}).values() if isinstance(c, dict)]
    burden["modelled_small_area"] = any("small-area" in x for x in srcs)
    level_note = None
    if burden["modelled_small_area"] and r.get("geo_level") == "county":
        level_note = ("the Long COVID burden is a modelled small-area estimate (BRFSS 2023 MRP), not observed county "
                      "prevalence; within-state differences come from county composition and covariates only")
    if bool(r.get("burden_inherited")) and r.get("geo_level") == "county":
        if len(members) == 1:
            level_note = ("the burden is the state estimate inherited by the county (source resolution: state); its "
                          "evidence level describes the state measure, not county prevalence")
        else:
            level_note = ("the set's evidence level is its least direct contributing member; the long-COVID member's "
                          "burden is the state estimate inherited by the county (source resolution: state; its level "
                          "describes the state measure, not county prevalence)")
    if not bool(r.get("rank_eligible", True)):
        unc = unc + ["Small-population county: shown only because small counties were included."]
    incl_small = rank_col.endswith("_incl_small")
    ew_mc = weight_set == O.EVIDENCE_WEIGHT_SET and not incl_small and r.get("rank_mc_ew_p05") is not None \
        and pd.notna(r.get("rank_mc_ew_p05"))
    mc_lo, mc_hi, mc_p10 = (("rank_mc_ew_p05", "rank_mc_ew_p95", "p_top10_mc_ew") if ew_mc else
                            ("rank_mc_p05", "rank_mc_p95", "p_top10_mc"))
    geo_block = {
        "name": r["geo_name"], "fips": r["geo_id"], "level": r["geo_level"], "object_id": r["geo_object_id"],
        "state": r["state_abbr"], "population": _num(r.get("population_total"), 0),
        "burden": burden, "burden_evidence_level": r.get("burden_evidence_level"),
        "burden_evidence_level_note": level_note,
        "burden_excluded_level_D": bool(r.get("burden_excluded_level_D")),
        "vulnerability": {"svi_overall": _num(r.get("svi_overall")), "percentile": _num(r.get("vulnerability_pct"))},
        "diagnostic_desert": {"index": _num(r.get("diagnostic_desert")), "percentile": _num(r.get("diagnostic_desert_pct")),
                              "variant": r.get("diagnostic_desert_variant")},
        "clinic_capacity": {"index": _num(r.get("clinic_capacity")), "percentile": _num(r.get("clinic_capacity_pct")),
                            "implementer_providers_in_county_per_100k": _num(r.get("cc_impl_providers_in_geo_per_100k"), 1),
                            "implementer_providers_within_50km_per_100k": _num(r.get("cc_impl_providers_within_50km_per_100k"), 1),
                            "implementer_facilities_in_county_per_100k": _num(r.get("cc_impl_facilities_in_geo_per_100k"), 1),
                            "implementer_facilities_within_50km_per_100k": _num(r.get("cc_impl_facilities_within_50km_per_100k"), 1)},
        "research_readiness": {"index": _num(r.get("research_readiness")),
                               "percentile": _num(r.get("research_readiness_pct")),
                               "condition_trials_in_pool": _num(r.get("rr_condition_trials_pool_n"), 0),
                               "condition_nih_core_projects_in_pool": _num(r.get("rr_condition_nih_core_pool_n"), 0),
                               "technology_experience_trials_in_pool": _num(r.get("rr_tech_trials_pool_n"), 0)},
        "technology_saturation": {"facilities_with_technology_experience_per_100k": _num(r.get("technology_saturation"), 2),
                                  "percentile": _num(r.get("technology_saturation_pct")),
                                  "note": "reported only; not in the default formula"},
        "composite": _num(r.get(f"composite_{weight_set}")), "rank": _num(r.get(rank_col), 0),
        "n_regions_ranked": n_ranked,
        "rank_basis": {"weight_set": weight_set, "rank_column": rank_col,
                       "universe": ("all counties incl. small ones" if rank_col.endswith("_incl_small") else
                                    "population-eligible regions with a complete burden")},
        "rank_interval_5_95": [_num(r.get(mc_lo), 0), _num(r.get(mc_hi), 0)],
        "p_top10_monte_carlo": _num(r.get(mc_p10), 3),
        "monte_carlo_basis": (("rank_interval_5_95 and p_top10_monte_carlo are drawn around the EVIDENCE_WEIGHTED "
                               "weights (Dirichlet) with burden uncertainty and the sensitivity / specificity CIs of the "
                               "performance record, among population-eligible regions with a complete burden; they "
                               "describe the evidence-weighted rank") if ew_mc else
                              ("rank_interval_5_95 and p_top10_monte_carlo are drawn around the EQUAL weights (Dirichlet) "
                               "with burden uncertainty, among population-eligible regions with a complete burden; they "
                               "describe the equal-weight rank"
                               + ("" if weight_set == "equal" and not incl_small else
                                  f", not the rank under '{weight_set}'"
                                  + (" among all counties incl. small ones" if incl_small else "")))),
        "ranks_under_weight_sets": {k: _num(r.get(f"rank_{k}"), 0) for k in O.weight_sets()},
        "measurement_evidence": {"value": _num(r.get("measurement_evidence")),
                                 "performance_status": r.get("measurement_performance_status"),
                                 "performance_tier": _num(r.get("measurement_performance_tier"), 0),
                                 "object_id": r.get("measurement_performance_object_id"),
                                 "note": "constant for every region of this condition x measurement: it changes "
                                         "which measurement is preferred, not the order of regions"},
        "expected_yield": expected_yield_block(r),
        "evidence_weighted": {"status": r.get("evidence_weighted_status"),
                              "composite": _num(r.get(f"composite_{O.EVIDENCE_WEIGHT_SET}")),
                              "rank": _num(r.get(f"rank_{O.EVIDENCE_WEIGHT_SET}"), 0),
                              "rank_interval_5_95": [_num(r.get("rank_mc_ew_p05"), 0), _num(r.get("rank_mc_ew_p95"), 0)],
                              "p_top10_monte_carlo": _num(r.get("p_top10_mc_ew"), 3),
                              "joint_rank_over_region_x_measurement": _num(r.get(f"rank_{O.EVIDENCE_WEIGHT_SET}_joint"), 0)},
    }
    sites, site_status, rev = [], [], {}
    if cb is not None and i is not None:
        sites, site_status = candidate_sites(r["geo_id"], members, meas["measurement_id"])
        rev = research_evidence(cb, i)
        for s in site_status:
            if s.get("absence_note"):
                unc.append(f"{s['condition_id']}: {s['absence_note']}")
    phen = ctx.get("phenotype", [])
    mol = ctx.get("molecular", [])
    tech = ctx.get("technology", {})
    member_opps = [f"opportunity:{c}|{meas['measurement_id']}|{r['geo_id']}" for c in members]
    member_opps = [o for o in member_opps if o in (member_rows or set())]
    perf_ids = [r.get("measurement_performance_object_id")] if r.get("measurement_performance_object_id") else []
    try:
        perf_ids += [i for i in json.loads(r.get("measurement_performance_source_object_ids") or "[]")
                     if not str(i).startswith("published_evidence:")]
    except (TypeError, ValueError):
        pass
    prov_ids = ([r["object_id"]] if stored else []) + [r["geo_object_id"], f"measurement:{meas['measurement_id']}"] + \
        perf_ids + \
        [f"condition:{c}" for c in members] + member_opps + [s["object_id"] for s in sites]
    prov_ids += (rev.get("condition_trials", {}).get("object_ids", []) +
                 rev.get("technology_experience_trials", {}).get("object_ids", []) +
                 rev.get("nih_projects", {}).get("object_ids", []))
    for p in phen:
        for d in p.get("datasets", []):
            prov_ids += [f["object_id"] for f in d.get("top_features", []) if f.get("object_id")]
    site_names = [s["facility_name"] for s in sites[:3]]
    next_step = next_step_text(geo_block, weight_set, rank_col, n_ranked, meas, spec, r["geo_name"], site_names)
    return {
        "condition": {"condition_id": spec["condition_id"], "label": spec["label"], "kind": spec["kind"],
                      "members": members, "object_ids": [f"condition:{c}" for c in members]},
        "phenotype": phen if phen else UNKNOWN,
        "measurement": {"measurement_id": meas["measurement_id"], "label": meas["label"], "kind": meas["kind"],
                        "members": meas["members"], "implementer_groups": meas["implementer_groups"],
                        "deployment_complexity": meas["deployment_complexity"],
                        "adapter": meas.get("adapter_name"), "adapter_status": meas.get("adapter_status"),
                        "object_id": f"measurement:{meas['measurement_id']}"},
        "technology": {"name": meas["label"], "regulatory_context": tech.get("regulatory_context", UNKNOWN),
                       "technology_evidence": tech.get("technology_evidence", []),
                       "measurement_performance": tech.get("measurement_performance", UNKNOWN),
                       "member_condition_performance": tech.get("member_condition_performance", []),
                       "expected_yield_here": geo_block["expected_yield"],
                       "guardrail": tech.get("guardrail")},
        "geography": geo_block,
        "candidate_sites": sites if sites else UNKNOWN,
        "candidate_site_query_status": site_status,
        "research_evidence": rev if rev else UNKNOWN,
        "molecular_context": mol if mol else UNKNOWN,
        "uncertainties": unc,
        "provenance": {"object_ids": list(dict.fromkeys(prov_ids)),
                       "opportunity_row": {"object_id": r["object_id"], "stored": bool(stored),
                                           "table": O.TABLE if stored else None,
                                           "note": ("stored row; trace_evidence resolves it" if stored
                                                    else ON_THE_FLY_NOTE)},
                       "sources": [r.get("source_name")], "source_version": r.get("source_version"),
                       "retrieved_at": r.get("retrieved_at"),
                       "tools": ["measure_it.scoring.opportunity", "facilities.matching.find_candidate_clinics",
                                 "wearables.signatures.get_patient_phenotype_signature",
                                 "omics.query.get_molecular_context", "measurements.regulatory.get_regulatory_context",
                                 "measurements.query.get_measurement_evidence"]},
        "recommended_next_step": next_step,
    }


# ------------------------------------------------------------------------------------------------------------------
# demo outputs
# ------------------------------------------------------------------------------------------------------------------

def _demo_frames(res: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    regions, sites, ev = [], [], []
    for rec in res["recommendations"]:
        g = rec["geography"]
        regions.append({"rank": g["rank"], "geo_id": g["fips"], "geo_name": g["name"], "state": g["state"],
                        "population": g["population"], "composite": g["composite"],
                        "rank_interval_p05": g["rank_interval_5_95"][0], "rank_interval_p95": g["rank_interval_5_95"][1],
                        "p_top10_monte_carlo": g["p_top10_monte_carlo"],
                        "burden_evidence_level": g["burden_evidence_level"],
                        "burden_pct": g["burden"].get("percentile"),
                        "burden_inherited": g["burden"]["inherited"],
                        "vulnerability_pct": g["vulnerability"]["percentile"],
                        "diagnostic_desert_pct": g["diagnostic_desert"]["percentile"],
                        "clinic_capacity_pct": g["clinic_capacity"]["percentile"],
                        "research_readiness_pct": g["research_readiness"]["percentile"],
                        "technology_saturation_pct": g["technology_saturation"]["percentile"],
                        "condition_trials_in_pool": g["research_readiness"]["condition_trials_in_pool"],
                        "condition_nih_core_projects_in_pool": g["research_readiness"]["condition_nih_core_projects_in_pool"],
                        "technology_experience_trials_in_pool": g["research_readiness"]["technology_experience_trials_in_pool"],
                        **{f"rank_{k}": v for k, v in g["ranks_under_weight_sets"].items()},
                        "measurement_performance_status": g["measurement_evidence"]["performance_status"],
                        "measurement_performance_tier": g["measurement_evidence"]["performance_tier"],
                        "measurement_evidence": g["measurement_evidence"]["value"],
                        "reach": g["expected_yield"].get("reach"),
                        "expected_yield_basis": g["expected_yield"].get("basis"),
                        "expected_yield_component": g["expected_yield"].get("expected_yield_component"),
                        "expected_yield_index": g["expected_yield"].get("expected_yield_index"),
                        "expected_detectable_cases": g["expected_yield"].get("expected_detectable_cases"),
                        "expected_false_positives_upper": g["expected_yield"].get("expected_false_positives_upper"),
                        "expected_false_positives": g["expected_yield"].get("expected_false_positives"),
                        "rank_evidence_weighted_interval": g["evidence_weighted"]["rank_interval_5_95"],
                        "n_candidate_sites": len(rec["candidate_sites"]) if isinstance(rec["candidate_sites"], list) else 0,
                        "object_id": rec["provenance"]["opportunity_row"]["object_id"]})
        for s in rec["candidate_sites"] if isinstance(rec["candidate_sites"], list) else []:
            sites.append({"region_rank": g["rank"], "geo_id": g["fips"], "geo_name": g["name"],
                          "facility_object_id": s["object_id"], "facility_name": s["facility_name"],
                          "primary_kind": s["primary_kind"], "city": s["city"], "state": s["state"],
                          "distance_km": s["distance_km"], "inside_geography": s["inside_geography"],
                          "selected_for_conditions": SEP.join(s["selected_for_conditions"]),
                          "selected_via": json.dumps(s["selected_via"]), "reasons": " | ".join(s["reasons"]),
                          "framing": s["framing"]})
        rv = rec["research_evidence"] if isinstance(rec["research_evidence"], dict) else {}
        for kind, key in (("condition_trial", "condition_trials"), ("technology_experience_trial",
                                                                   "technology_experience_trials"),
                          ("nih_project", "nih_projects")):
            for oid in (rv.get(key) or {}).get("object_ids", []):
                ev.append({"region_rank": g["rank"], "geo_id": g["fips"], "geo_name": g["name"], "kind": kind,
                           "object_id": oid})
    return pd.DataFrame(regions), pd.DataFrame(sites), pd.DataFrame(ev)


def deployment_candidates_table(res: dict, query_text: str) -> pd.DataFrame:
    rows = []
    for rec in res["recommendations"]:
        g = rec["geography"]
        rows.append({"object_id": rec["provenance"]["opportunity_row"]["object_id"], "query_text": query_text,
                     "query_args": json.dumps(res["query"]), "condition_id": rec["condition"]["condition_id"],
                     "measurement_id": rec["measurement"]["measurement_id"], "geo_level": g["level"],
                     "geo_id": g["fips"], "geo_object_id": g["object_id"], "geo_name": g["name"],
                     "weight_set": res["weight_set"], "rank": g["rank"], "composite": g["composite"],
                     "rank_interval_p05": g["rank_interval_5_95"][0], "rank_interval_p95": g["rank_interval_5_95"][1],
                     "burden_evidence_level": g["burden_evidence_level"],
                     "burden_inherited": bool(g["burden"].get("inherited")),
                     "rank_evidence_weighted": g["evidence_weighted"]["rank"],
                     "evidence_weighted_status": g["evidence_weighted"]["status"],
                     "measurement_performance_object_id": g["measurement_evidence"]["object_id"],
                     "measurement_performance_status": g["measurement_evidence"]["performance_status"],
                     "measurement_performance_tier": g["measurement_evidence"]["performance_tier"],
                     "measurement_evidence": g["measurement_evidence"]["value"],
                     "expected_yield_basis": g["expected_yield"].get("basis"),
                     "expected_yield_index": g["expected_yield"].get("expected_yield_index"),
                     "expected_false_positives_upper": g["expected_yield"].get("expected_false_positives_upper"),
                     "n_candidate_sites": len(rec["candidate_sites"]) if isinstance(rec["candidate_sites"], list) else 0,
                     "candidate_site_object_ids": SEP.join(s["object_id"] for s in rec["candidate_sites"])
                     if isinstance(rec["candidate_sites"], list) else "",
                     "n_uncertainties": len(rec["uncertainties"]),
                     "recommended_next_step": rec["recommended_next_step"],
                     "recommendation_json": json.dumps(rec)})
    df = pd.DataFrame(rows)
    df = add_provenance(
        df, data_layer="derived", source_name=f"{PRODUCER} (deployment_opportunities + facilities.matching + "
                                              "wearables.signatures + omics.query + measurements.regulatory/query)",
        source_version=f"recommend {utc_now_iso()}; {res.get('deployment_opportunities_version')}",
        retrieved_at=utc_now_iso(), evidence_type="derived_score", source_record_id="object_id",
        source_geographic_resolution="county" if res["query"]["geography_level"] == "county" else "state",
        evidence_level=df["burden_evidence_level"].astype(str).to_numpy(),
        provenance_notes="Candidate deployment opportunities for pilot evaluation, not validated diagnostic pathways. "
                         "Ecological joins only; candidate sites are facilities with characteristics suggesting they "
                         "may be viable implementation or study partners.")
    # same rule as deployment_opportunities (CONVENTIONS 3.5): a county row whose burden is a state value inherited
    # by the county keeps source_geographic_resolution = state
    df["source_geographic_resolution"] = np.where((df["geo_level"] == "state") | df["burden_inherited"], "state",
                                                  "county")
    return df


def draft_maps(df_primary: pd.DataFrame, top_ids: list[str]) -> list[str]:
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = []
    poly = gpd.read_parquet(processed_path("county_boundaries_2024").with_suffix(".geoparquet"))
    gid = next(c for c in ("geo_id", "GEOID", "geoid") if c in poly.columns)
    poly = poly.rename(columns={gid: "geo_id"})
    poly = poly[~poly["geo_id"].str[:2].isin(["02", "15", "60", "66", "69", "72", "78"])].to_crs(5070)
    d = poly.merge(df_primary[["geo_id", "composite_equal", "rank_equal", "rank_burden_only", "rank_mc_p05",
                               "rank_mc_p95", "rank_eligible"]], on="geo_id", how="left")
    MAPS_D = MAPS / "drafts"
    MAPS_D.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(13, 8))
    d.plot(column="composite_equal", ax=ax, cmap="viridis", legend=True, linewidth=0,
           missing_kwds={"color": "#dddddd"}, legend_kwds={"shrink": 0.5, "label": "composite (equal weights)"})
    top = d[d["geo_id"].isin(top_ids)]
    top.boundary.plot(ax=ax, color="red", linewidth=1.6)
    for r in top.itertuples():
        p = r.geometry.representative_point()
        ax.annotate(f"{int(r.rank_equal)}", (p.x, p.y), color="red", fontsize=8, ha="center", weight="bold")
    ax.set_axis_off()
    ax.set_title("DRAFT - Candidate deployment opportunities: Long COVID or ME/CFS x wearable autonomic/activity "
                 "monitoring (county, equal weights)\nRed = top 10 (counties with population >= 10,000); grey = "
                 "no composite (incl. incomplete burden: no CMS ME/CFS value, e.g. the Connecticut planning regions). "
                 "Alaska/Hawaii omitted in this draft.\nEcological ranking; long-COVID burden is the state estimate "
                 "inherited by counties.", fontsize=9)
    p = MAPS_D / "demo_opportunity_county.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    out.append(str(p))
    # Test 5: burden-only vs full top-25
    fig, ax = plt.subplots(figsize=(13, 8))
    d.plot(ax=ax, color="#eeeeee", linewidth=0)
    full = d["rank_equal"] <= 25
    bo = d["rank_burden_only"] <= 25
    d[full & bo].plot(ax=ax, color="#6a3d9a", linewidth=0)
    d[full & ~bo].plot(ax=ax, color="#1f78b4", linewidth=0)
    d[~full & bo].plot(ax=ax, color="#e31a1c", linewidth=0)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#6a3d9a", label="top-25 in both"),
                       Patch(color="#1f78b4", label="top-25 full ranking only (enters)"),
                       Patch(color="#e31a1c", label="top-25 burden-only only (leaves)")], loc="lower left", fontsize=8)
    ax.set_axis_off()
    ax.set_title("DRAFT - Test 5: burden-only vs full (equal-weight) top-25 counties, Long COVID or ME/CFS x wearable "
                 "monitoring. Burden-only ties are broken by FIPS (see test5 tables for tie-averaged overlap).",
                 fontsize=9)
    p = MAPS_D / "test5_burden_only_vs_full_top25.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    out.append(str(p))
    # rank-interval width
    fig, ax = plt.subplots(figsize=(13, 8))
    d["interval_width"] = d["rank_mc_p95"] - d["rank_mc_p05"]
    d.plot(column="interval_width", ax=ax, cmap="magma_r", legend=True, linewidth=0, missing_kwds={"color": "#dddddd"},
           legend_kwds={"shrink": 0.5, "label": "5th-95th percentile rank interval width"})
    ax.set_axis_off()
    ax.set_title("DRAFT - Monte Carlo rank-interval width (burden CI x evidence multiplier + Dirichlet weights), "
                 "Long COVID or ME/CFS x wearable monitoring", fontsize=9)
    p = MAPS_D / "demo_rank_interval_width.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    out.append(str(p))
    return out


def run(write: bool = True) -> dict:
    res = rank_deployment_opportunities(**DEMO_ARGS)
    res["query_text"] = DEMO_QUERY
    if not write or res.get("status") != "ok":
        return res
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "example_deployment_recommendation.json").write_text(json.dumps(res, indent=2, default=str))
    regions, sites, ev = _demo_frames(res)
    TABLES.mkdir(parents=True, exist_ok=True)
    regions.to_csv(TABLES / "demo_top10_regions.csv", index=False)
    # the same query under the pre-declared primary alternative (docs/ANALYSIS_PLAN_METRIC_LINK.md)
    res_ew = rank_deployment_opportunities(**DEMO_ARGS_EVIDENCE)
    res_ew["query_text"] = DEMO_QUERY + " [weight set: evidence_weighted]"
    (RESULTS / "example_deployment_recommendation_evidence_weighted.json").write_text(
        json.dumps(res_ew, indent=2, default=str))
    if res_ew.get("status") == "ok":
        _demo_frames(res_ew)[0].to_csv(TABLES / "demo_top10_regions_evidence_weighted.csv", index=False)
    res["evidence_weighted_demo"] = {"status": res_ew.get("status"), "reason": res_ew.get("reason"),
                                     "file": "results/example_deployment_recommendation_evidence_weighted.json"}
    sites.to_csv(TABLES / "demo_top10_candidate_sites.csv", index=False)
    ev.to_csv(TABLES / "demo_top10_research_evidence.csv", index=False)
    dc = deployment_candidates_table(res, DEMO_QUERY)
    write_table(dc, "deployment_candidates", producer=PRODUCER,
                description="One row per recommendation of the SPEC demo query, with the full SPEC-schema JSON.")
    cid = res["condition_resolution"]["condition_id"]
    mid = res["measurement_resolution"]["measurement_id"]
    df, _ = _rows(cid, mid, "county")
    maps = draft_maps(df, [r["geography"]["fips"] for r in res["recommendations"]])
    print(f"[recommend] {len(res['recommendations'])} recommendations; maps: {maps}", flush=True)
    return res


if __name__ == "__main__":
    run()
