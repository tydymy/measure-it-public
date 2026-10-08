"""Page 5: Deployment Recommendation - the SPEC demo query answered by rank_deployment_opportunities: a transparent
ranked table (every component, rank interval, ranks under other weight sets), a map with candidate sites, the evidence
with trace_evidence, the limitations and the recommendation JSON."""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from components import charts as C  # noqa: E402
from components import data as D  # noqa: E402
from components import maps as M  # noqa: E402
from components.common import (FRAMING, ToolLog, caveat_list, fmt, json_bytes, level_badge,  # noqa: E402
                               limitations_panel, ok, records_df, setup_page, show_df, trace_explorer, unknown_box)
from measure_it.config import UNKNOWN, load_config  # noqa: E402
from measure_it.scoring.recommend import DEMO_ARGS, DEMO_QUERY  # noqa: E402

COMPONENTS = ["burden", "vulnerability", "diagnostic_desert", "clinic_capacity", "research_readiness"]

setup_page("Deployment Recommendation",
           "Where should a candidate measurement be deployed to generate the most useful clinical evidence, and which "
           "facilities could plausibly participate? Answered by rank_deployment_opportunities, one row per region.")
log = ToolLog()
stamp = D.stamp()
cat = D.catalog(stamp)

st.markdown(f"**Example query (SPEC):** _{DEMO_QUERY}_")
with st.expander("Query parameters", expanded=True):
    q1, q2, q3, q4 = st.columns([1.6, 1.6, 0.9, 1.1])
    cond = q1.text_input("Condition or condition set", value=DEMO_ARGS["condition"], key="rec_condition",
                         help="A condition, 'A or B', or a condition-set alias ('Long COVID or ME/CFS', 'demo cluster')")
    meas = q2.text_input("Measurement class or bundle", value=DEMO_ARGS["measurement"], key="rec_measurement")
    level = q3.radio("Level", ["county", "state"], horizontal=True, key="rec_level")
    ws_names = list(cat["weight_sets"])
    ws = q4.selectbox("Weight set", ws_names, index=ws_names.index(DEMO_ARGS["weight_set"]), key="rec_weights")
    q5, q6, q7 = st.columns([1, 1, 1.4])
    top_n = q5.number_input("Regions", 1, 25, int(DEMO_ARGS["top_n"]), key="rec_topn")
    n_sites = q6.number_input("Candidate sites per region (map)", 1, 10, 3, key="rec_sites")
    small = q7.checkbox(f"Include counties below {cat['min_population']:,} residents", value=False, key="rec_small",
                        disabled=level != "county")

with st.spinner("Ranking candidate deployment opportunities ..."):
    env = log.call("rank_deployment_opportunities", condition=cond, measurement=meas, geography_level=level,
                   top_n=int(top_n), weight_set=ws, include_small_population=bool(small), detail="full")
if not ok(env):
    unknown_box(env, "Deployment recommendation")
    _perf = (env.get("data") or {}).get("measurement_performance")
    if isinstance(_perf, dict):
        st.warning("Measurement performance for this condition and measurement: "
                   f"{_perf.get('performance_status')}. {_perf.get('performance_note') or ''} Under the "
                   "'evidence_weighted' weights such a measurement is not ranked (never scored 0 or 1); choose the "
                   "'equal' weights for a ranking that carries no measurement-performance information.",
                   icon=":material/help:")
    limitations_panel(log, ["No ranking was produced for this query; nothing is filled in."])
    st.stop()

d = env["data"]
recs = d["recommendations"]
w = d.get("weights") or {}
st.info(FRAMING, icon=":material/info:")
mc = st.columns(4)
mc[0].metric("Regions ranked", fmt(d.get("n_regions_ranked"), 0))
mc[1].metric("Condition", (d.get("condition_resolution") or {}).get("condition_id", UNKNOWN))
mc[2].metric("Measurement", (d.get("measurement_resolution") or {}).get("measurement_id", UNKNOWN))
mc[3].metric("Weight set", ws)

# ------------------------------------------------------------------------------------------------ metric link
perf = d.get("measurement_performance") or {}
st.subheader("The measurement's own evidence (metric link)")
pm = st.columns(5)
pm[0].metric("Performance", str(perf.get("performance_status") or UNKNOWN).replace(" / NOT AVAILABLE", ""))
pm[1].metric("Tier", fmt(perf.get("quality_tier"), 1) if (perf.get("quality_tier") or 4) < 4 else "none")  # 2.5 = BYOD
pm[2].metric("AUROC (lower 95%)", f"{fmt(perf.get('auroc'), 3)} ({fmt(perf.get('auroc_ci_low'), 3)})"
             if perf.get("auroc") is not None else UNKNOWN)
pm[3].metric("Sensitivity @ specificity", f"{fmt(perf.get('op_sensitivity'))} @ {fmt(perf.get('op_specificity'))}"
             if perf.get("op_sensitivity") is not None else UNKNOWN)
pm[4].metric("measurement_evidence", fmt(perf.get("measurement_evidence"), 3)
             if perf.get("measurement_evidence") is not None else UNKNOWN)
if perf.get("performance_status") == "known":
    st.caption(f"Record {perf.get('selected_record_id')} ({perf.get('tier_text')}): {perf.get('label_basis')} vs "
               f"{perf.get('comparator')}; n = {fmt(perf.get('n_cases'), 0)} / {fmt(perf.get('n_controls'), 0)}. "
               f"{perf.get('healthy_control_caveat')} {perf.get('expected_yield_language')}")
else:
    st.warning(f"{perf.get('performance_note') or UNKNOWN} This ranking carries no measurement-performance "
               "information; under the 'evidence_weighted' weights this measurement is not ranked.",
               icon=":material/help:")
if perf.get("evidence_conflict"):
    st.warning(f"Evidence conflict: {perf.get('evidence_conflict_note')}", icon=":material/warning:")
st.caption(f"Universe: {d.get('ranking_universe')}. Basis: {d.get('ranking_basis')} "
           f"(version {d.get('deployment_opportunities_version')}). Weights: "
           + ", ".join(f"{k} {v:.2f}" for k, v in w.items()) + f". Method: {d.get('method')}")

# ------------------------------------------------------------------------------------------------ ranked table
# The Monte Carlo rank interval and P(top 10) are always drawn around the EQUAL weights (docs/SCORING.md section 4);
# under another weight set the rank can fall outside the interval, so the columns name their weights.
MC_TAG = "" if ws in ("equal", "evidence_weighted") else " (equal-weight MC)"
RI, PT = "rank interval 5-95" + MC_TAG, "P(top 10)" + MC_TAG
rows = []
for r in recs:
    g = r["geography"]
    b = g.get("burden") or {}
    members = b.get("members") or {}
    mem_txt = "; ".join(f"{c}: {fmt(m.get('burden_value'))} (level {m.get('burden_evidence_level') or 'D'}"
                        f"{', inherited state value' if m.get('burden_inherited') else ''})" for c, m in members.items())
    if not members:
        mem_txt = (f"{fmt(b.get('value'))} {b.get('unit') or ''} (level {g.get('burden_evidence_level') or 'D'}"
                   f"{', inherited state value' if b.get('inherited') else ''})")
    ranks = [v for v in (g.get("ranks_under_weight_sets") or {}).values() if isinstance(v, (int, float))]
    lo, hi = (g.get("rank_interval_5_95") or [None, None])[:2]
    rows.append({   # component percentiles right after the composite, so they are visible without scrolling
        "rank": g.get("rank"), "region": g.get("name"), "composite": g.get("composite"),
        "burden pct": b.get("percentile"),
        "vulnerability pct": (g.get("vulnerability") or {}).get("percentile"),
        "desert pct": (g.get("diagnostic_desert") or {}).get("percentile"),
        "capacity pct": (g.get("clinic_capacity") or {}).get("percentile"),
        "readiness pct": (g.get("research_readiness") or {}).get("percentile"),
        "measurement evidence": (g.get("measurement_evidence") or {}).get("value"),
        "expected yield": (g.get("expected_yield") or {}).get("expected_yield_component"),
        RI: f"{fmt(lo, 0)}-{fmt(hi, 0)}", PT: g.get("p_top10_monte_carlo"),
        "expected detectable cases": (g.get("expected_yield") or {}).get("expected_detectable_cases")
        if isinstance((g.get("expected_yield") or {}).get("expected_detectable_cases"), (int, float)) else
        ("index " + fmt((g.get("expected_yield") or {}).get("expected_yield_index"), 3)
         if (g.get("expected_yield") or {}).get("expected_yield_index") is not None else UNKNOWN),
        "expected false positives": (g.get("expected_yield") or {}).get("expected_false_positives")
        if isinstance((g.get("expected_yield") or {}).get("expected_false_positives"), (int, float)) else
        (g.get("expected_yield") or {}).get("expected_false_positives_upper"),
        "PPV": (g.get("expected_yield") or {}).get("expected_ppv"),
        "reach": (g.get("expected_yield") or {}).get("reach"),
        "evidence-weighted rank": (g.get("evidence_weighted") or {}).get("rank"),
        "burden level": (g.get("burden_evidence_level") or "D") + (" (inherited state value)" if b.get("inherited")
                                                                   else ""),
        "burden (members)": mem_txt,
        "state": g.get("state"), "fips": g.get("fips"), "population": g.get("population"),
        "SVI": (g.get("vulnerability") or {}).get("svi_overall"),
        "desert index": (g.get("diagnostic_desert") or {}).get("index"),
        "capacity index": (g.get("clinic_capacity") or {}).get("index"),
        "readiness index": (g.get("research_readiness") or {}).get("index"),
        "condition trials in reach": (g.get("research_readiness") or {}).get("condition_trials_in_pool"),
        "NIH core projects in reach": (g.get("research_readiness") or {}).get("condition_nih_core_projects_in_pool"),
        "technology saturation pct (reported only)": (g.get("technology_saturation") or {}).get("percentile"),
        "rank range across weight sets": f"{min(ranks):.0f}-{max(ranks):.0f}" if ranks else UNKNOWN,
        # stored rankings cite their deployment_opportunities row; an on-the-fly ranking has no stored row
        "object_id": (((r.get("provenance") or {}).get("opportunity_row") or {}).get("object_id")
                      if ((r.get("provenance") or {}).get("opportunity_row") or {}).get("stored")
                      else "computed on the fly; not a stored row (trace geo / condition / measurement ids)"),
        "_lo": lo, "_hi": hi})
tdf = pd.DataFrame(rows)
st.subheader("Ranked candidate deployment opportunities")
show_df(tdf.drop(columns=["_lo", "_hi"]), key="rec_table", column_config={
    "composite": st.column_config.NumberColumn(format="%.3f"),
    PT: st.column_config.NumberColumn(format="%.2f"),
    **{c: st.column_config.NumberColumn(format="%.2f") for c in tdf.columns if c.endswith("pct") or c.endswith(
        "index") or c in ("SVI", "measurement evidence", "expected yield", "reach", "PPV")},
    "expected false positives": st.column_config.NumberColumn(format="%d"),
    "population": st.column_config.NumberColumn(format="%d")})
# method numbers come from the ranking output and the scoring config, not typed into the page
_draws = next((m.group(1) for rec in recs for u in rec.get("uncertainties") or []
               for m in [re.search(r"over (\d+) draws", str(u))] if m), None)
_inh = load_config("scoring").get("inherited_burden_uncertainty_multiplier")
st.caption("pct = percentile rank among the ranked regions (0-1). The composite is the weighted mean of the component "
           "percentiles; the rank interval is the 5th-95th percentile rank over "
           + (f"{int(_draws):,} " if _draws else "") + "Monte Carlo draws of burden "
           "(CI x evidence-level multiplier" + (f"; x{_inh:g} for a state value inherited by a county" if _inh else "")
           + ") and weights. The 'rank range across weight sets' is the best and worst rank under the "
           f"{len(ws_names)} named weight sets."
           + ("" if ws in ("equal", "evidence_weighted") else
              f" The Monte Carlo is drawn around the equal weights, so under the '{ws}' "
              "weights a region's rank can fall outside its interval; the interval and P(top 10) are labelled "
              "'equal-weight MC'.")
           + (" Under 'evidence_weighted' the interval is drawn around those weights and also varies the performance "
              "record's sensitivity / specificity." if ws == "evidence_weighted" else "")
           + " Metric link: 'measurement evidence' is the same for every region (it separates measurements, not "
             "regions); 'expected yield' = burden percentile x sensitivity x reach; 'expected detectable cases' is a "
             "count only where a defensible count exists (otherwise an index), and 'expected false positives' is then "
             "an upper bound. These are planning estimates for a candidate pilot, not predictions of diagnoses.")

# ------------------------------------------------------------------------------------------------ charts
lvl_note = sorted({r["burden level"] for r in rows})
cc = st.columns(2)
with cc[0]:
    parts, comp_rows, gap = [], [], []
    wsum_all = {k: float(v) for k, v in w.items()}
    for r in recs:
        g = r["geography"]
        pct = {"burden": (g.get("burden") or {}).get("percentile"),
               "vulnerability": (g.get("vulnerability") or {}).get("percentile"),
               "diagnostic_desert": (g.get("diagnostic_desert") or {}).get("percentile"),
               "clinic_capacity": (g.get("clinic_capacity") or {}).get("percentile"),
               "research_readiness": (g.get("research_readiness") or {}).get("percentile")}
        pct["measurement_evidence"] = (g.get("measurement_evidence") or {}).get("value")
        pct["expected_yield"] = (g.get("expected_yield") or {}).get("expected_yield_component")
        if "low_saturation" in wsum_all:
            ts = (g.get("technology_saturation") or {}).get("percentile")
            pct["low_saturation"] = None if ts is None else 1 - ts
        avail = {k: v for k, v in pct.items() if isinstance(v, (int, float)) and wsum_all.get(k, 0) > 0}
        den = sum(wsum_all[k] for k in avail) or np.nan
        row = {"region": f"{g.get('rank'):.0f}. {g.get('name')}"}
        for k in wsum_all:
            row[k] = wsum_all[k] * avail[k] / den if k in avail else 0.0
        comp_rows.append(row)
        gap.append(abs(sum(row[k] for k in wsum_all) - (g.get("composite") or np.nan)))
    cdf = pd.DataFrame(comp_rows)
    short = {"burden": "burden", "vulnerability": "vulnerability", "diagnostic_desert": "desert",
             "clinic_capacity": "capacity", "research_readiness": "readiness", "low_saturation": "low saturation",
             "measurement_evidence": "measurement evidence", "expected_yield": "expected yield"}
    parts = [(k, short.get(k, k)) for k in wsum_all if wsum_all[k] > 0]
    st.markdown("**What each component contributes** (weight x percentile / sum of weights of the available "
                "components; the bar length is the composite)")
    st.plotly_chart(C.stacked_contributions(cdf, "region", parts, xtitle=f"composite ({ws} weights)"),
                    key="rec_contrib")
    if gap and np.nanmax(gap) > 0.01:
        st.caption(f"Note: the stacked parts differ from the tool's composite by up to {np.nanmax(gap):.3f}; the table "
                   "shows the tool's value.")
with cc[1]:
    st.markdown("**Rank and its 5th-95th percentile interval** (Monte Carlo over burden and weights"
                + (")" if ws == "equal" else f", drawn around the equal weights; the dot is the rank under '{ws}')"))
    ri = tdf.assign(label=[f"{r:.0f}. {n}" for r, n in zip(tdf["rank"], tdf["region"])]).dropna(subset=["_lo", "_hi"])
    if len(ri):
        st.plotly_chart(C.error_bars(ri, "rank", "label", "_lo", "_hi", xtitle="rank (1 = top); interval 5th-95th pct",
                                     hover_cols=[RI, PT, "burden level"]),
                        key="rec_intervals")
    st.caption(f"Burden evidence level(s) in this list: {', '.join(lvl_note)}. A wide interval means the region's place "
               "in the short list is not stable.")

# ------------------------------------------------------------------------------------------------ map
st.subheader("Map")
fig, map_note = M.ranking_map(d, ws, level, cat, stamp, log=log, height=560)
site_rows = []
for r in recs:
    for s in (r.get("candidate_sites") if isinstance(r.get("candidate_sites"), list) else [])[:int(n_sites)]:
        # candidate_sites carry the facility's lat / lon / geocode_precision (rank_deployment_opportunities)
        site_rows.append({"facility_id": s.get("facility_id") or str(s.get("object_id", "")).split(":", 1)[-1],
                          "facility": s.get("facility_name"), "region": r["geography"]["name"],
                          "selected via": "; ".join(f"{k}: {v}" for k, v in (s.get("selected_via") or {}).items()),
                          "distance_km": s.get("distance_km"), "object_id": s.get("object_id"),
                          "lat": s.get("lat"), "lon": s.get("lon"), "geocode_precision": s.get("geocode_precision")})
sdf = pd.DataFrame(site_rows)
if len(sdf):
    M.add_points(fig, sdf.dropna(subset=["lat", "lon"]), "sites",
                 ["facility", "region", "selected via", "distance_km", "geocode_precision"], base_size=9,
                 name=f"Candidate sites (first {int(n_sites)} per region)")
st.plotly_chart(fig, key="rec_map")
st.caption(map_note + "; orange diamonds = candidate sites (facilities with characteristics suggesting they may be "
           "viable implementation or study partners), at ZIP/ZCTA or city centroids.")

# ------------------------------------------------------------------------------------------------ region detail
st.subheader("Evidence for one region")
opts = list(range(len(recs)))
i = st.selectbox("Region", opts, key="rec_region",
                 format_func=lambda k: f"{recs[k]['geography']['rank']:.0f}. {recs[k]['geography']['name']}")
rec = recs[i]
g = rec["geography"]
st.markdown(level_badge(g.get("burden_evidence_level")) + (f" <span style='font-size:0.85rem'>"
                                                            f"{g.get('burden_evidence_level_note')}</span>"
                                                            if g.get("burden_evidence_level_note") else ""),
            unsafe_allow_html=True)
st.success(rec.get("recommended_next_step") or UNKNOWN, icon=":material/arrow_forward:")
t1, t2, t3, t4 = st.tabs(["Candidate sites", "Research evidence", "Uncertainties", "Provenance"])
with t1:
    cs = rec.get("candidate_sites")
    if isinstance(cs, list) and cs:
        srows = []
        for s in cs:
            row = {"facility": s.get("facility_name"), "kind": s.get("primary_kind"), "city": s.get("city"),
                   "km": s.get("distance_km"), "inside": s.get("inside_geography"),
                   "selected via": "; ".join(f"{k}: {v}" for k, v in (s.get("selected_via") or {}).items())}
            for cid_, chars in (s.get("characteristics") or {}).items():
                for k, v in (chars or {}).items():
                    if isinstance(v, dict) and v.get("value") not in (None, 0):
                        row[f"{cid_}: {k}"] = f"{v.get('value')} (rank {v.get('rank') if v.get('rank') else '-'} of "
                        row[f"{cid_}: {k}"] += f"{v.get('n_ranked')})"
            row["object_id"] = s.get("object_id")
            srows.append(row)
        show_df(pd.DataFrame(srows), key="rec_sites_table")
        with st.expander("Reasons per facility"):
            for s in cs:
                st.markdown(f"**{s.get('facility_name')}**: " + " ".join(map(str, s.get("reasons") or [])))
    else:
        st.warning(f"Candidate sites: {UNKNOWN}")
    for q in rec.get("candidate_site_query_status") or []:
        if q.get("absence_note"):
            st.info(f"{q.get('condition_id')}: {q['absence_note']}")
with t2:
    rev = rec.get("research_evidence")
    if isinstance(rev, dict):
        st.caption(f"Pool: {rev.get('pool')}. {rev.get('note') or ''}")
        for k in ("condition_trials", "technology_experience_trials", "nih_projects"):
            v = rev.get(k) or {}
            n = v.get("n", v.get("n_core_projects"))
            st.markdown(f"**{k.replace('_', ' ')}**: {fmt(n, 0)} " + (f"({v.get('note')})" if v.get("note") else ""))
            if v.get("object_ids"):
                st.caption(", ".join(v["object_ids"][:25]))
            if v.get("examples"):
                show_df(records_df(v["examples"]), key=f"rec_ev_{k}")
    else:
        st.warning(f"Research evidence: {UNKNOWN}")
with t3:
    caveat_list(rec.get("uncertainties") or [])
with t4:
    prov = rec.get("provenance") or {}
    st.caption(f"Sources: {'; '.join(map(str, prov.get('sources') or []))}. Version: {prov.get('source_version')}. "
               f"Tools: {', '.join(prov.get('tools') or [])}.")
    trace_explorer(prov.get("object_ids") or [], key="rec_trace", log=log,
                   label=f"Trace evidence: {len(prov.get('object_ids') or [])} object ids cited by this region")

# ------------------------------------------------------------------------------------------------ shared context
st.subheader("Condition- and measurement-level context (identical for every region)")
x1, x2, x3 = st.tabs(["Phenotype (public cohorts)", "Condition-level molecular enrichment", "Technology"])
with x1:
    st.caption("Group statistics from public person-level cohorts; these participants are not the people of any "
               "ranked region.")
    prow = []
    for p in rec.get("phenotype") if isinstance(rec.get("phenotype"), list) else []:
        for ds in p.get("datasets") or []:
            tf = (ds.get("top_features") or [{}])
            nr = ds.get("null_results")
            prow.append({"condition": p.get("condition_id"), "dataset": ds.get("dataset_id"),
                         "proxy": ds.get("is_proxy"), "cases": ds.get("n_cases"), "controls": ds.get("n_controls"),
                         "tested": ds.get("n_features_tested"), "FDR-significant": ds.get("n_features_fdr_significant"),
                         "null results": nr if isinstance(nr, str) else (nr or {}).get("summary", "")
                         if isinstance(nr, dict) else "",
                         "top feature": (tf[0] or {}).get("feature") if tf else None,
                         "effect": (tf[0] or {}).get("effect_size") if tf else None,
                         "phenotype": ds.get("phenotype_label")})
    show_df(pd.DataFrame(prow), key="rec_phen")
with x2:
    st.info("Condition-level molecular enrichment: public databases, not participant-linked, not patient multi-omics, "
            "not mechanism.", icon=":material/biotech:")
    mrow = []
    for m in rec.get("molecular_context") if isinstance(rec.get("molecular_context"), list) else []:
        for et, v in (m.get("evidence_summary") or {}).items():
            mrow.append({"condition": m.get("condition_id"), "evidence type": et,
                         **{k: (", ".join(f"{a} ({b})" for a, b in v[k].items()) if isinstance(v.get(k), dict)
                                else v.get(k)) for k in v}})
    show_df(records_df(mrow), key="rec_mol")
with x3:
    tech = rec.get("technology") or {}
    st.caption(tech.get("guardrail") or "")
    rrow = [{"class": r.get("measurement_id"), "status": r.get("status"),
             **{k: v for k, v in (r.get("summary") or {}).items() if not isinstance(v, (list, dict))}}
            for r in tech.get("regulatory_context") or [] if isinstance(r, dict)]
    show_df(pd.DataFrame(rrow), key="rec_reg")
    trow = [{"condition": e.get("condition_id"), "class": e.get("measurement_id"),
             "evidence tier": (e.get("measurement_evidence_strength") or {}).get("tier"),
             "technology maturity": e.get("technology_maturity"), "regulatory visibility": e.get("regulatory_visibility"),
             "deployment complexity": e.get("deployment_complexity"), "object_id": e.get("evidence_object_id")}
            for e in tech.get("technology_evidence") or [] if isinstance(e, dict)]
    show_df(pd.DataFrame(trow), key="rec_techev")
    st.markdown("**Measurement performance** (the record the evidence-weighted ranking uses; members of a set shown "
                "for context only, never substituted):")
    perows = [tech.get("measurement_performance")] + list(tech.get("member_condition_performance") or [])
    show_df(pd.DataFrame([{"object_id": x.get("object_id"), "status": x.get("performance_status"),
                           "tier": x.get("quality_tier"), "record": x.get("selected_record_id"),
                           "AUROC": x.get("auroc"), "AUROC lower 95%": x.get("auroc_ci_low"),
                           "sensitivity": x.get("op_sensitivity"), "specificity": x.get("op_specificity"),
                           "comparator": x.get("comparator"), "measurement_evidence": x.get("measurement_evidence")}
                          for x in perows if isinstance(x, dict)]), key="rec_perf")

exc = d.get("excluded_incomplete_burden") or []
if exc:
    with st.expander(f"Regions not ranked because their burden is incomplete ({len(exc)})"):
        st.caption(d.get("excluded_incomplete_burden_note") or "")
        show_df(records_df(exc, ["name", "fips", "state", "population", "members_without_value", "reason",
                                 "components_kept"]), key="rec_excluded")

# ------------------------------------------------------------------------------------------------ downloads
st.subheader("Recommendation JSON")
compact = log.call("rank_deployment_opportunities", condition=cond, measurement=meas, geography_level=level,
                   top_n=int(top_n), weight_set=ws, include_small_population=bool(small), detail="compact")
dc = st.columns(2)
dc[0].download_button("Download recommendation JSON (compact, SPEC schema per region)", json_bytes(compact),
                      file_name=f"deployment_recommendation_{level}_{ws}.json", mime="application/json",
                      key="rec_dl_compact", disabled=not ok(compact))
dc[1].download_button("Download full-detail JSON (all sites and context)", json_bytes(env),
                      file_name=f"deployment_recommendation_{level}_{ws}_full.json", mime="application/json",
                      key="rec_dl_full")
st.caption("Same envelopes as the MCP tool and the API route /tools/rank_deployment_opportunities: status, data (SPEC "
           "deployment-opportunity schema per region), caveats, provenance (object ids + sources) and truncation.")

limitations_panel(log, list(d.get("general_uncertainties") or []) + [
    d.get("framing") or "",
    "Burden: the set's evidence level is its least direct member; a long-COVID county value is the state estimate "
    "inherited by the county; ME/CFS burden is a level-C claims proxy; POTS and dysautonomia have no usable burden "
    "(level D).",
    "The short list is weight-sensitive and its rank intervals are wide (results/SCORING_RESULTS.md, Tests 5-7).",
    "Measurement performance comes from case-control comparisons against healthy or recovered controls, which "
    "overstate real-world performance; expected yields are planning estimates for a candidate pilot, not predictions "
    "of diagnoses (results/SCORING_RESULTS.md section 12).",
])
