"""Page 2: Condition Explorer - one condition through the engine: ontology, phenotype signatures (incl. null results),
wearable abnormalities from public cohorts, condition-level molecular enrichment, candidate measurements with the five
separate dimensions, FDA technologies and burden sources with evidence levels."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from components import charts as C  # noqa: E402
from components import data as D  # noqa: E402
from components.common import (EVIDENCE_LEVELS, ToolLog, badge, fmt, is_unknown, level_badge,  # noqa: E402
                               limitations_panel, ok, records_df, sections, setup_page, show_df, trace_explorer,
                               unknown_box)
from measure_it.config import UNKNOWN  # noqa: E402

setup_page("Condition Explorer",
           "Type a condition (name, alias, acronym, ICD-10-CM code or MONDO id). Everything below comes from the tool "
           "facade; each section says which tool answered.")
log = ToolLog()
stamp = D.stamp()

c1, c2, c3 = st.columns([2, 2, 1])
text = c1.text_input("Condition", value="Long COVID", key="cond_text")
axes = D.phenotype_axes(stamp)
axis_ids = [""] + [a for a, _ in axes]
axis_lbl = dict(axes)
phen = c2.selectbox("Phenotype focus for candidate measurements (optional)", axis_ids, key="cond_phenotype",
                    format_func=lambda a: "(none: all phenotypes)" if not a else f"{axis_lbl.get(a, a)} ({a})")
top_n = c3.number_input("Top features / candidates", 3, 20, 8, key="cond_topn")

norm = log.call("normalize_condition", condition_or_code=text)
nd = norm.get("data") or {}
matches = nd.get("matches") or []
status = nd.get("match_status")
if not ok(norm) or not matches:
    unknown_box(norm, f"Condition '{text}'")
    srch = log.call("search_condition", condition=text)
    if ok(srch):
        st.markdown("**Closest registry conditions (search_condition):**")
        show_df(records_df(srch["data"]["candidates"], ["canonical_condition_id", "preferred_name", "score",
                                                        "match_reason"]), key="cond_search")
    limitations_panel(log, ["The condition did not resolve to one registry condition; nothing below is guessed."])
    st.stop()

if status in ("ambiguous", "partial"):
    # the tool facade never uses a partial or ambiguous match (tools.resolve_condition_arg); neither does this page
    st.warning(f"**Condition '{text}': {UNKNOWN}**. match_status is **{status}**: the text is not one unambiguous "
               "registry condition, so no condition is chosen for you. Type one of the candidates below (name, id or "
               "code).")
    show_df(records_df(matches + list(nd.get("other_candidates") or []),
                       ["canonical_condition_id", "preferred_name", "match_reason", "confidence"]),
            key="cond_candidates_ambiguous")
    limitations_panel(log, ["The condition did not resolve to one registry condition; nothing below is guessed."])
    st.stop()
m = matches[0]
cid = m["canonical_condition_id"]
st.markdown(f"### {m.get('preferred_name')}  \n"
            + " ".join(badge(b) for b in [f"condition:{cid}", m.get("primary_mondo_id") or "",
                                          *(m.get("icd10cm_codes") or [])[:3], *(m.get("mesh_ids") or [])[:2]] if b)
            + f"  \n<span style='opacity:0.7;font-size:0.85rem'>{m.get('match_reason')}</span>",
            unsafe_allow_html=True)

tabs = sections(["Phenotype features", "Wearable abnormalities", "Molecular enrichment", "Candidate measurements",
                "Technologies (FDA)", "Burden sources"])

# ------------------------------------------------------------------------------------------------ phenotype
sig = log.call("get_patient_phenotype_signature", condition=cid, top_n=int(top_n))
with tabs[0]:
    st.caption("get_patient_phenotype_signature: group differences in public person-level datasets (NHANES 2011-2014, "
               "Stanford wearable cohorts). Group statistics only; these participants are never located or linked to "
               "any place, facility or trial.")
    if not ok(sig):
        unknown_box(sig, "Phenotype signatures")
    else:
        ds = sig["data"].get("datasets") or []
        rows = [{"dataset": d.get("dataset_id"), "phenotype": d.get("phenotype_label"),
                 "label basis": d.get("label_basis"), "proxy label": d.get("is_proxy"),
                 "cases": d.get("n_cases"), "controls": d.get("n_controls"), "effect measure": d.get("effect_measure"),
                 "features tested": d.get("n_features_tested"),
                 "FDR-significant": d.get("n_features_fdr_significant"),
                 "null results": (d.get("null_results") or {}).get("summary"),
                 "phenotype-level null": (d.get("null_results") or {}).get("phenotype_level_null")} for d in ds]
        show_df(pd.DataFrame(rows), key="sig_datasets")
        for d in ds:
            with st.expander(f"{d.get('dataset_id')} | {d.get('phenotype_label')}"):
                st.markdown(f"**Definition:** {d.get('phenotype_definition')}")
                if d.get("caveats"):
                    st.caption(str(d.get("caveats")))
                feats = records_df(d.get("top_features"), ["feature_label", "feature", "effect_size", "ci_low",
                                                           "ci_high", "q_value_bh", "object_id"])
                if len(feats):
                    st.markdown("FDR-significant features (BH q < 0.05):")
                    show_df(feats, key=f"sig_feats_{d.get('dataset_id')}_{d.get('phenotype_id')}")
                else:
                    st.markdown("**No feature passed BH-FDR q < 0.05** (a null result, reported as such).")
                nr = d.get("null_results") or {}
                if nr.get("features"):
                    st.caption("Null (not significant) features: " + ", ".join(map(str, nr["features"])))
                mlr = (d.get("model_level_results") or {}).get("rows") or []
                if mlr:
                    st.markdown("Model-level results stored by the producer (not per-feature signatures):")
                    show_df(records_df(mlr), key=f"sig_models_{d.get('dataset_id')}_{d.get('phenotype_id')}")

with tabs[1]:
    st.caption("Wearable features that differ between cases and controls in public cohorts (effect size with 95% CI; "
               "only BH-FDR-significant features are drawn, null features are listed below).")
    if not ok(sig):
        unknown_box(sig, "Wearable abnormalities")
    else:
        pts, nulls = [], []
        for d in sig["data"].get("datasets") or []:
            plabel = str(d.get("phenotype_label") or d.get("phenotype_id"))
            grp = f"{d.get('dataset_id')}: {plabel[:90]}" + (" [proxy label]" if d.get("is_proxy") else "")
            for f in d.get("top_features") or []:
                if C.finite(f.get("effect_size")) == C.finite(f.get("effect_size")):
                    pts.append({"group": grp, "feature": f.get("feature_label") or f.get("feature"),
                                "effect": f.get("effect_size"), "ci_low": f.get("ci_low"), "ci_high": f.get("ci_high"),
                                "q": fmt(f.get("q_value_bh"), 4), "measure": d.get("effect_measure"),
                                "id": f.get("object_id")})
            nulls += [{"group": grp, "null features": ", ".join(map(str, (d.get("null_results") or {})
                                                                     .get("features") or []))}]
        if pts:
            fdf = pd.DataFrame(pts)
            fdf["feature"] = fdf["feature"].astype(str).str.slice(0, 60)
            st.plotly_chart(C.forest(fdf, "group", "feature", hover_cols=["group", "measure", "q", "id"],
                                     xtitle="effect size (see effect measure; 0 = no difference)"),
                            key="cond_forest")
        else:
            st.info(f"No FDR-significant wearable feature for this condition in the public cohorts processed here "
                    f"({UNKNOWN} as a positive signal; null results below).")
        for d in sig["data"].get("datasets") or []:
            if not d.get("n_features_fdr_significant"):
                st.warning(f"**Null result:** {d.get('dataset_id')}: {d.get('phenotype_label')}: 0 of "
                           f"{d.get('n_features_tested')} tested wearable features passed BH-FDR q < 0.05 "
                           f"({fmt(d.get('n_cases'), 0)} cases, {fmt(d.get('n_controls'), 0)} controls).")
        n_sig = len(pts)
        st.caption(f"{n_sig} FDR-significant feature(s) drawn. Datasets and their null (not significant) features:")
        show_df(pd.DataFrame(nulls), key="cond_nulls")

# ------------------------------------------------------------------------------------------------ molecular
mol = log.call("get_molecular_context", condition=cid, top_n=5)
with tabs[2]:
    st.info("**Condition-level molecular enrichment** from public databases (Open Targets, GWAS Catalog, GEO/SRA, "
            "mapMECFS supplements). Not measured on any participant, not participant-linked, **not patient "
            "multi-omics** and not a mechanism.", icon=":material/biotech:")
    if not ok(mol):
        unknown_box(mol, "Molecular context")
    else:
        md = mol["data"]
        ev = md.get("evidence") or {}
        erows = [{"evidence type": k, "records": v.get("n_records"), "distinct entities": v.get("n_distinct_entities"),
                  "sources": "; ".join(f"{s} ({n})" for s, n in (v.get("counts_by_source") or {}).items())}
                 for k, v in ev.items() if isinstance(v, dict)]
        show_df(pd.DataFrame(erows), key="mol_evidence")
        top = []
        for k, v in ev.items():
            for et, ents in ((v or {}).get("top_entities") or {}).items():
                for e in ents[:5]:
                    top.append({"evidence type": k, "entity type": et, "id": e.get("entity_id"),
                                "label": e.get("entity_label"), "records": e.get("n_records"),
                                "source score": e.get("max_source_score"), "rank basis": e.get("rank_basis")})
        with st.expander(f"Top entities per evidence type ({len(top)})"):
            show_df(pd.DataFrame(top), key="mol_top")
        coh = md.get("coherence") or {}
        st.markdown(f"**Test 3 coherence** (status: {coh.get('status')}): physiological systems supported by "
                    "pathway enrichment of the condition's public gene lists.")
        sysd = records_df(coh.get("supported_systems"), ["physiological_system", "source", "fdr", "n_enriched_pathways",
                                                         "condition_specific", "post_hoc_flagged",
                                                         "top_enriched_pathways"])
        if "fdr" in sysd:
            sysd["fdr"] = [f"{v:.1e}" if isinstance(v, (int, float)) else v for v in sysd["fdr"]]
        show_df(sysd, key="mol_systems")
        st.caption(coh.get("post_hoc_flag_note") or "")
        agr = md.get("cross_source_genetic_agreement") or {}
        if agr:
            st.caption(f"Cross-source genetic agreement: Open Targets genetic genes {agr.get('open_targets_genetic_genes')}"
                       f", GWAS Catalog mapped genes {agr.get('gwas_catalog_mapped_genes')}, in both "
                       f"{agr.get('in_both')}.")
        with st.expander("Gene-list agreement between sources"):
            show_df(records_df(coh.get("gene_agreement")), key="mol_agree")

# ------------------------------------------------------------------------------------------------ measurements
disc = log.call("discover_candidate_measurements", condition=cid, phenotype=phen or None, top_n=int(top_n) + 4)
cand_ids = []
with tabs[3]:
    st.caption("discover_candidate_measurements: candidates in the documented research-activity order. The five "
               "dimensions are kept separate (never one score); the order is research activity, not validity.")
    if not ok(disc):
        unknown_box(disc, "Candidate measurements")
    else:
        dd = disc["data"]
        if phen and isinstance(dd.get("phenotype"), dict):
            ph = dd["phenotype"]
            st.markdown(f"Phenotype focus: **{ph.get('axis_label') or phen}** ({ph.get('status')})")
        rows = []
        for c in dd.get("candidates") or []:
            dm = c.get("dimensions") or {}
            pss = dm.get("phenotype_signal_strength") or {}
            mes = dm.get("measurement_evidence_strength") or {}
            tm = dm.get("technology_maturity") or {}
            rv = dm.get("regulatory_visibility") or {}
            dc = dm.get("deployment_complexity") or {}
            cand_ids.append(c.get("measurement_id"))
            name = str(c.get("name") or c.get("measurement_id"))
            # the five dimensions first (tier / status words), the counts behind them after, so all five are visible
            rows.append({
                "rank": c.get("rank"), "measurement": name.split(" (")[0],
                "phenotype signal strength": str(pss.get("status") or pss.get("value")).replace("_", " "),
                "measurement evidence strength": mes.get("tier"),
                "technology maturity": tm.get("tier"),
                "regulatory visibility": rv.get("tier"),
                "deployment complexity": (f"{dc.get('value')} of 5" if not is_unknown(dc.get("value"))
                                          else UNKNOWN),
                "trials (objective use)": mes.get("n_trials_objective"),
                "NIH core projects": mes.get("n_nih_core_projects"),
                "FDA decisions": tm.get("n_decisions_high_medium_confidence"),
                "family": c.get("modality_family"),
                "text-mining precision": fmt((c.get("text_mining_precision") or {}).get("strict")),
                "full name": name, "object_id": c.get("evidence_object_id")})
        show_df(pd.DataFrame(rows), key="cond_candidates")
        st.caption("Deployment complexity is a curated ordinal (1 = home/consumer device, 5 = specialist lab "
                   "procedure; configs/relevance.yaml). Phenotype signal strength says whether a public person-level "
                   "dataset with this condition's label measured the class (a null result is shown as such).")
        st.caption(dd.get("ordering_note") or "")
        with st.expander("Reasons per candidate"):
            for c in dd.get("candidates") or []:
                st.markdown(f"**{c.get('name')}**: " + " ".join(map(str, c.get("reasons") or [])))

# ------------------------------------------------------------------------------------------------ FDA
with tabs[4]:
    st.caption("get_regulatory_context per candidate measurement class. An FDA record is a deployment-readiness "
               "signal, not evidence that a device detects or diagnoses the condition.")
    reg_rows, reg_sum = [], []
    for mid in cand_ids[:6]:
        reg = log.call("get_regulatory_context", technology=mid, max_records=3)
        if not ok(reg):
            reg_rows.append({"measurement": mid, "product code": UNKNOWN, "device": reg.get("reason")})
            reg_sum.append({"measurement": mid, "product codes": UNKNOWN, "reason": reg.get("reason")})
            continue
        rs = reg["data"].get("regulatory_status") or []
        sm = next((x for x in reg["data"].get("summary") or [] if isinstance(x, dict)), {})
        classes = sorted({str(r.get("device_class")) for r in rs if not is_unknown(r.get("device_class"))})
        reg_sum.append({"measurement": mid, "product codes": len(rs),
                        "device classes": ", ".join(classes) or UNKNOWN,
                        "best mapping confidence": sm.get("confidence_of_best_mapping"),
                        "510(k)": sm.get("n_510k_total"), "De Novo": sm.get("n_denovo_total"),
                        "PMA": sm.get("n_pma_total"),
                        "example codes": ", ".join(f"{r.get('product_code')} ({str(r.get('device_name'))[:40]})"
                                                   for r in rs[:3])})
        for r in rs:
            reg_rows.append({"measurement": mid, "product code": r.get("product_code"),
                             "device": r.get("device_name"), "class": r.get("device_class"),
                             "regulation": r.get("regulation_number"), "510(k)": r.get("n_510k"),
                             "De Novo": r.get("n_denovo"), "PMA": r.get("n_pma"),
                             "mapping confidence": r.get("mapping_confidence"),
                             "latest decision": r.get("most_recent_decision_date"),
                             "object_id": r.get("product_code_object_id")})
    show_df(pd.DataFrame(reg_sum), key="cond_fda_summary")
    with st.expander(f"Every mapped product code ({len(reg_rows)})"):
        show_df(pd.DataFrame(reg_rows), key="cond_fda")

# ------------------------------------------------------------------------------------------------ burden
with tabs[5]:
    st.caption("get_condition_burden at each level. Evidence levels: "
               + "; ".join(f"{k} = {v}" for k, v in EVIDENCE_LEVELS.items()) + ".")
    brow = []
    for lvl in ("national", "state", "county"):
        b = log.call("get_condition_burden", condition=cid, geography_level=lvl, max_rows=60, order="value_desc")
        bd = b.get("data") or {}
        s = bd.get("summary") or {}
        brow.append({"level": lvl, "status": b.get("status"), "evidence level": bd.get("burden_evidence_level"),
                     "primary measure": bd.get("primary_measure_id"),
                     "source resolution": bd.get("primary_source_resolution"), "inherited": bd.get("inherited"),
                     "rows": s.get("n_rows"), "min": s.get("min"), "median": s.get("median"), "max": s.get("max"),
                     "units": ", ".join(s.get("value_units") or []), "rationale / reason": bd.get("rationale")
                     or b.get("reason")})
        if lvl == "state":
            state_env = b
    lvl0 = brow[1]["evidence level"] or "D"
    st.markdown(level_badge(lvl0), unsafe_allow_html=True)
    show_df(pd.DataFrame(brow), key="cond_burden")
    if ok(state_env) and state_env["data"].get("rows"):
        sr = pd.DataFrame(state_env["data"]["rows"])
        sr = sr[pd.to_numeric(sr["value"], errors="coerce").notna()].head(15)
        if len(sr):
            sr["ci_low"] = pd.to_numeric(sr.get("ci_low"), errors="coerce").fillna(sr["value"])
            sr["ci_high"] = pd.to_numeric(sr.get("ci_high"), errors="coerce").fillna(sr["value"])
            st.markdown(f"Highest state values (level {lvl0}; {sr['value_unit'].iloc[0]}; 95% CI where published):")
            st.plotly_chart(C.error_bars(sr, "value", "geo_name", "ci_low", "ci_high",
                                         xtitle=f"{sr['measure_id'].iloc[0]} ({sr['value_unit'].iloc[0]})",
                                         hover_cols=["period", "source_name", "burden_evidence_level"]),
                            key="cond_state_burden")
    for c in (state_env.get("caveats") or [])[:4]:
        st.caption(c)

ids = []
for e in log.calls:
    ids += (e.get("provenance") or {}).get("object_ids") or []
trace_explorer([f"condition:{cid}"] + ids, key="cond_trace", log=log)

limitations_panel(log, [
    "Phenotype signatures are group statistics from public cohorts; several labels are proxies (e.g. an ME/CFS-like "
    "symptom proxy in NHANES) and some condition labels give null results (reported, not hidden).",
    "Molecular context is condition-level molecular enrichment from public databases; it is never linked to the "
    "wearable participants and is never patient multi-omics.",
    "Candidate measurements are ordered by research activity (registered trials, NIH projects), not by validity; "
    "trial registrations and FDA records are research-activity and deployment-readiness signals, not evidence of "
    "efficacy or diagnostic accuracy.",
    "Burden: the evidence level is shown with every burden value; a state value attached to counties is inherited "
    "context, never county prevalence.",
])
