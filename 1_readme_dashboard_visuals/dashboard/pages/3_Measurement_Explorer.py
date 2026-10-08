"""Page 3: Measurement Explorer - one measurement class or bundle: the conditions where research deploys it, the
physiological signals, public supporting datasets, relevant trials, FDA device classes and regions with unmet need."""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from components import charts as C  # noqa: E402
from components import data as D  # noqa: E402
from components import maps as M  # noqa: E402
from components.common import (FRAMING, ToolLog, badge, fmt, limitations_panel, ok, records_df, sections,  # noqa: E402
                               setup_page, show_df, trace_explorer, unknown_box)
from measure_it.config import UNKNOWN  # noqa: E402

setup_page("Measurement Explorer",
           "Type a measurement class or bundle (e.g. 'wearable autonomic monitoring', 'CPET', 'nailfold "
           "capillaroscopy', 'hrv'). Research activity, datasets, trials and FDA records are signals of where the "
           "measurement is being deployed, not evidence that it works.")
log = ToolLog()
stamp = D.stamp()
cat = D.catalog(stamp)
labels = cat["labels"]


def parse(v):
    if isinstance(v, (list, dict)):
        return v
    if not isinstance(v, str) or not v.strip():
        return None
    for f in (json.loads, ast.literal_eval):
        try:
            return f(v)
        except (ValueError, SyntaxError):
            continue
    return None


c1, c2, c3, c4 = st.columns([2, 1.6, 1.1, 1.2])
text = c1.text_input("Measurement", value="wearable autonomic monitoring", key="meas_text")
cond_opts = [c["id"] for c in cat["conditions"]] + [s["id"] for s in cat["sets"]]
ctx_cond = c2.selectbox("Condition for trials and unmet need", cond_opts, key="meas_condition",
                        index=cond_opts.index("long_covid"), format_func=lambda c: labels.get(c, c))
level = c3.radio("Unmet-need level", ["county", "state"], horizontal=True, key="meas_level")
ws_names = list(cat["weight_sets"])
ws = c4.selectbox("Unmet-need weight set", ws_names, key="meas_weights",
                  index=ws_names.index("access_gap") if "access_gap" in ws_names else 0)

members_ctx = cat["members"].get(ctx_cond, [ctx_cond])
res_env = log.call("get_measurement_evidence", measurement=text, condition=members_ctx[0], max_object_ids=10)
if not ok(res_env) or not isinstance((res_env.get("data") or {}).get("measurement"), dict):
    unknown_box(res_env, f"Measurement '{text}'")
    limitations_panel(log, ["The measurement did not resolve to a registry class or bundle; nothing below is "
                            "guessed."])
    st.stop()

mres = res_env["data"]["measurement"]
mid = (mres.get("measurement_ids") or [None])[0]
member_ids = mres.get("member_ids") or [mid]
tr = log.call("trace_evidence", object_id=f"measurement:{mid}", max_rows=1)
row = ((tr.get("data") or {}).get("rows") or [{}])[0] if ok(tr) else {}
st.markdown(f"### {row.get('name') or mid}  \n"
            + " ".join(badge(b) for b in [f"measurement:{mid}", mres.get("kind"), row.get("modality_family"),
                                          f"adapter: {row.get('adapter') or 'none'}"] if b)
            + f"  \n<span style='opacity:0.7;font-size:0.85rem'>{mres.get('match_reason')}. Members: "
              f"{', '.join(member_ids)}. Adapter status: {row.get('adapter_status') or UNKNOWN}.</span>",
            unsafe_allow_html=True)

tabs = sections(["Performance evidence (metric link)", "Conditions where it is studied", "Signals & public datasets",
                "Relevant trials", "FDA device classes", "Regions with unmet need"])

# ------------------------------------------------------------------------------------------------ performance
with tabs[0]:
    st.caption("measurement_performance (docs/ANALYSIS_PLAN_METRIC_LINK.md): the performance record the ranking uses for "
               "each condition x this measurement. Tier 1 = this project's own computation on a clinical / study case "
               "definition; tier 2 = own computation on a proxy or self-reported label; tier 3 = published claim, not "
               "reproduced; UNKNOWN = no record (never imputed from another condition or measurement). Sensitivity is at "
               ">= 0.90 specificity. Every record compares cases with healthy or recovered controls, which overstates "
               "real-world performance.")
    pf = D.measurement_performance(stamp)
    pf = pf[pf["measurement_id"] == mid] if len(pf) else pf
    if not len(pf):
        st.warning(f"No measurement_performance rows for {mid}: {UNKNOWN} (the grid covers the scored conditions x the "
                   "four bundles, their member classes and blood biomarkers / metabolomics).")
    else:
        known = pf[pf["performance_status"] == "known"]
        mcol = st.columns(4)
        mcol[0].metric("Conditions with a known record", f"{len(known)} of {len(pf)}")
        mcol[1].metric("Best tier present", fmt(pf.loc[pf["quality_tier"] < 4, "quality_tier"].min(), 1)
                       if (pf["quality_tier"] < 4).any() else UNKNOWN)
        mcol[2].metric("Highest measurement_evidence", fmt(pf["measurement_evidence"].max(), 3)
                       if pf["measurement_evidence"].notna().any() else UNKNOWN)
        mcol[3].metric("UNKNOWN cells", int((pf["performance_status"] == UNKNOWN).sum()))
        prow = [{"condition": labels.get(r.condition_id, r.condition_id), "status": r.performance_status,
                 "tier": fmt(r.quality_tier, 1) if r.quality_tier < 4 else "none",
                 "record": r.selected_record_id or UNKNOWN,
                 "AUROC (95% CI)": (f"{r.auroc:.3f} ({r.auroc_ci_low:.3f}-{r.auroc_ci_high:.3f})"
                                    if pd.notna(r.auroc) else UNKNOWN),
                 "sensitivity @ spec": (f"{r.op_sensitivity:.2f} ({r.op_sensitivity_ci_low:.2f}-"
                                        f"{r.op_sensitivity_ci_high:.2f}) @ {r.op_specificity:.2f}"
                                        if pd.notna(r.op_sensitivity) else UNKNOWN),
                 "measurement_evidence": r.measurement_evidence, "comparator": r.comparator or UNKNOWN,
                 "label basis": r.label_basis or UNKNOWN, "conflict": r.evidence_conflict_note or "",
                 "object_id": r.object_id} for r in pf.itertuples()]
        pdf = pd.DataFrame(prow)
        plot = pdf.dropna(subset=["measurement_evidence"])
        if len(plot):
            st.plotly_chart(C.hbar(plot, "measurement_evidence", "condition",
                                   xtitle="measurement_evidence = tier factor x (AUROC lower 95% bound - 0.5) / 0.5",
                                   hover_cols=["record", "AUROC (95% CI)", "sensitivity @ spec", "tier"],
                                   text="measurement_evidence"), key="meas_perf_bar")
        show_df(pdf, key="meas_perf", column_config={"measurement_evidence": st.column_config.NumberColumn(format="%.3f")})
        for r in pf[pf["evidence_conflict"].fillna(False).astype(bool)].itertuples():
            st.warning(f"{labels.get(r.condition_id, r.condition_id)}: {r.evidence_conflict_note}")
        st.info("An UNKNOWN or partial record means this measurement is not ranked under the 'evidence_weighted' "
                "weights for that condition; its 'equal'-weight rank carries no measurement-performance information.",
                icon=":material/help:")

# ------------------------------------------------------------------------------------------------ conditions
with tabs[1]:
    st.caption("get_measurement_evidence(measurement, condition) for every registry condition: registered trials "
               "that describe objective use of the measurement (literal condition match) and NIH core projects.")
    rows, unknown_rows = [], []
    for c in cat["conditions"]:
        env = res_env if c["id"] == members_ctx[0] else log.call("get_measurement_evidence", measurement=text,
                                                                  condition=c["id"], max_object_ids=10)
        if not ok(env):
            unknown_rows.append(f"**{c['label']}**: {env.get('status')}. {env.get('reason')}")
            continue
        evs = env["data"].get("evidence") or []
        e = next((x for x in evs if x.get("measurement_id") == mid), evs[0] if evs else {})
        cnt, dm = e.get("counts") or {}, e.get("dimensions") or {}
        pss = dm.get("phenotype_signal_strength") or {}
        rows.append({"condition": c["label"], "trials (objective use)": cnt.get("n_trials_objective"),
                     "as outcome measure": cnt.get("n_trials_outcome_measure"),
                     "recruiting": cnt.get("n_trials_recruiting"), "NIH core projects": cnt.get("n_nih_core_projects"),
                     "evidence tier": (dm.get("measurement_evidence_strength") or {}).get("tier"),
                     "phenotype signal": pss.get("status") or pss.get("value"),
                     "object_id": e.get("evidence_object_id")})
    cdf = pd.DataFrame(rows)
    if "trials (objective use)" in cdf:
        cdf = cdf.sort_values("trials (objective use)", ascending=False, na_position="last")
        plot = cdf.dropna(subset=["trials (objective use)"])
        if len(plot):
            st.plotly_chart(C.hbar(plot, "trials (objective use)", "condition",
                                   xtitle="registered trials describing objective use (literal condition match)",
                                   hover_cols=["recruiting", "NIH core projects", "evidence tier"],
                                   text="trials (objective use)"), key="meas_cond_bar")
    show_df(cdf, key="meas_conditions")
    if unknown_rows:
        st.markdown("Not in the table (nothing filled in):  \n" + "  \n".join(unknown_rows))

# ------------------------------------------------------------------------------------------------ signals / datasets
with tabs[2]:
    sig = parse(row.get("signals")) or []
    st.markdown("**Physiological signals** (measurement registry): "
                + (" ".join(badge(s) for s in sig) if sig else UNKNOWN), unsafe_allow_html=True)
    ds = parse(row.get("public_person_level_datasets")) or []
    st.markdown("**Public person-level datasets that contain this measurement** (group statistics only; never linked "
                "to places):")
    show_df(records_df(ds), key="meas_datasets")
    res = []
    for e in res_env["data"].get("evidence") or []:
        for r in e.get("phenotype_signal_results") or []:
            res.append({"class": e.get("measurement_id"), "dataset": r.get("dataset_id"), "label": r.get("label_id"),
                        "label class": r.get("label_class"), "role": r.get("role"), "cases": r.get("n_cases"),
                        "controls": r.get("n_controls"), "features tested": r.get("n_features_tested"),
                        "FDR-significant": r.get("n_features_fdr"), "best feature": r.get("best_feature"),
                        "best effect": fmt(r.get("best_effect")), "model": r.get("model_metric"),
                        "model estimate": fmt(r.get("model_estimate")), "verdict": r.get("verdict"),
                        "object_id": r.get("object_id")})
    st.markdown(f"**Phenotype-signal results in public data** (context condition: {labels.get(ctx_cond, ctx_cond)}; "
                "`role` says when a dataset is adjacent evidence, e.g. an acute-infection cohort, not the condition):")
    show_df(pd.DataFrame(res).drop_duplicates(subset=["object_id"]) if res else pd.DataFrame(), key="meas_signal")

# ------------------------------------------------------------------------------------------------ trials
with tabs[3]:
    tri = log.call("find_relevant_trials", condition=ctx_cond, measurement=text, max_trials=25)
    st.caption("find_relevant_trials: ClinicalTrials.gov trials registered for the condition (literal condition "
               "match) whose registered text mentions a member measurement class. A registration is a "
               "research-activity signal.")
    if not ok(tri):
        unknown_box(tri, "Relevant trials")
    else:
        s = tri["data"]["summary"]
        mc = st.columns(4)
        mc[0].metric("Trials", s.get("n_trials"))
        mc[1].metric("Recruiting", s.get("n_recruiting"))
        mc[2].metric("Interventional", (s.get("by_study_type") or {}).get("INTERVENTIONAL", 0))
        mc[3].metric("Observational", (s.get("by_study_type") or {}).get("OBSERVATIONAL", 0))
        show_df(records_df(tri["data"]["trials"], ["object_id", "brief_title", "overall_status", "study_type",
                                                   "phases", "start_date", "lead_sponsor", "n_us_locations",
                                                   "measurement_classes_matched", "measurement_match_example"]),
                key="meas_trials")

# ------------------------------------------------------------------------------------------------ FDA
with tabs[4]:
    st.caption("get_regulatory_context per member class: mapped FDA product codes (curated mapping with a confidence), "
               "device class, regulation number and decision counts. A regulatory record is a deployment-readiness "
               "signal, not evidence that the device detects the condition.")
    reg_rows = []
    for m in member_ids:
        reg = log.call("get_regulatory_context", technology=m, max_records=3)
        if not ok(reg):
            reg_rows.append({"member": m, "product code": UNKNOWN, "device": reg.get("reason")})
            continue
        for r in reg["data"].get("regulatory_status") or []:
            reg_rows.append({"member": m, "product code": r.get("product_code"), "device": r.get("device_name"),
                             "class": r.get("device_class"), "regulation": r.get("regulation_number"),
                             "specialty": r.get("medical_specialty"), "510(k)": r.get("n_510k"),
                             "De Novo": r.get("n_denovo"), "PMA": r.get("n_pma"),
                             "U.S. establishments": r.get("n_registered_establishments_us"),
                             "mapping confidence": r.get("mapping_confidence"),
                             "latest decision": r.get("most_recent_decision_date"),
                             "object_id": r.get("product_code_object_id")})
    rdf = pd.DataFrame(reg_rows)
    if len(rdf) and "class" in rdf:
        cls = rdf.dropna(subset=["class"]).groupby("class").size()
        st.markdown("Device classes among mapped product codes: "
                    + ", ".join(f"class {k}: {v} code(s)" for k, v in cls.items()))
    show_df(rdf, key="meas_fda")

# ------------------------------------------------------------------------------------------------ unmet need
with tabs[5]:
    st.caption(f"rank_deployment_opportunities({labels.get(ctx_cond, ctx_cond)}, {mid}, {level}) under the "
               f"'{ws}' weights (access_gap puts 0.4 on the diagnostic desert and 0 on research readiness). "
               "Combinations outside the precomputed table are scored on the fly with the same engine.")
    st.markdown(FRAMING)
    with st.spinner("Ranking candidate deployment opportunities ..."):
        rk = log.call("rank_deployment_opportunities", condition=ctx_cond, measurement=mid, geography_level=level,
                      top_n=10, weight_set=ws)
    if not ok(rk):
        unknown_box(rk, "Regions with unmet need")
    else:
        # the Monte Carlo rank interval is always drawn around the EQUAL weights (docs/SCORING.md); under another
        # weight set the rank can lie outside it, so the column says which weights it belongs to
        ri_col = "rank interval 5-95" + ("" if ws in ("equal", "evidence_weighted") else " (equal-weight MC)")
        rows = []
        for r in rk["data"]["recommendations"]:
            g = r["geography"]
            rows.append({"rank": g.get("rank"), "region": g.get("name"), "composite": g.get("composite"),
                         "burden pct": (g.get("burden") or {}).get("percentile"),
                         "SVI pct": (g.get("vulnerability") or {}).get("percentile"),
                         "desert pct": (g.get("diagnostic_desert") or {}).get("percentile"),
                         "capacity pct": (g.get("clinic_capacity") or {}).get("percentile"),
                         "readiness pct": (g.get("research_readiness") or {}).get("percentile"),
                         "burden level": (g.get("burden_evidence_level") or "D")
                         + (" (inherited state value)" if (g.get("burden") or {}).get("inherited") else ""),
                         ri_col: "-".join(fmt(x, 0) for x in g.get("rank_interval_5_95") or []) or UNKNOWN,
                         "performance": (g.get("measurement_evidence") or {}).get("performance_status"),
                         "measurement evidence": (g.get("measurement_evidence") or {}).get("value"),
                         "expected yield": (g.get("expected_yield") or {}).get("expected_yield_component"),
                         "expected false positives (upper)": (g.get("expected_yield") or {}).get(
                             "expected_false_positives_upper", (g.get("expected_yield") or {}).get(
                                 "expected_false_positives")),
                         "geo_id": g.get("fips")})
        udf = pd.DataFrame(rows)
        st.caption(f"{rk['data'].get('n_regions_ranked')} regions ranked ({rk['data'].get('ranking_universe')}); "
                   f"basis: {rk['data'].get('ranking_basis')}.")
        fig, note = M.ranking_map(rk["data"], ws, level, cat, stamp, log=log, height=440)
        st.plotly_chart(fig, key="meas_unmet_map")
        st.caption(note + ".")
        show_df(udf, key="meas_unmet", column_config={
            c: st.column_config.NumberColumn(format="%.3f" if c == "composite" else "%.2f")
            for c in udf.columns if c == "composite" or c.endswith("pct")})
        if ws not in ("equal", "evidence_weighted"):
            st.caption(f"The rank is under the '{ws}' weights; the 5th-95th percentile rank interval comes from Monte "
                       "Carlo draws around the equal weights, so a rank can fall outside it.")
        st.caption("Expected yield and false positives are planning estimates for a candidate pilot, not predictions "
                   "of diagnoses; the false-positive column is an upper bound where no defensible case count exists.")

ids = [f"measurement:{mid}"]
for e in log.calls:
    ids += (e.get("provenance") or {}).get("object_ids") or []
trace_explorer(ids, key="meas_trace", log=log)

limitations_panel(log, [
    "Counts of registered trials and NIH projects come from text mining of registered text (audited precision per "
    "class in the measurement registry); absence of a mention is not evidence of non-use.",
    "A trial using a measurement is not evidence the measurement works; an FDA product code is a deployment-readiness "
    "signal, not evidence that the device detects or diagnoses the condition.",
    "Public datasets are person-level cohorts summarised as group statistics; some are adjacent evidence (e.g. acute "
    "COVID-19 cohorts) rather than the named condition.",
    "Regions with unmet need are candidate deployment opportunities under the chosen weights; burden carries its "
    "evidence level and state values on counties are inherited context.",
])
