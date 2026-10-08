"""Page 6: Your data - user-supplied (bring-your-own-data) datasets on this machine: manifest, locked plan, evaluation
metrics, the performance record that crosses into the ranking, and how the ranking moved (docs/BRING_YOUR_OWN_DATA.md).

Only aggregate results are shown; no person row or participant id is displayed. The embedding panel shows the dataset's
own Digital Phenotype Vector PCA as a picture, never pooled with another dataset."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from components.common import ToolLog, badge, caveat_list, limitations_panel, show_df, unknown_box  # noqa: E402
from components.common import setup_page  # noqa: E402
from measure_it.byod import common as K  # noqa: E402
from measure_it.byod import tools as BT  # noqa: E402
from measure_it.config import UNKNOWN  # noqa: E402

setup_page("Your data",
           "User-supplied datasets added with `measure-it byod` (validate -> ingest -> evaluate -> deploy -> remove). "
           "Local-only; only an aggregate performance record reaches the ranking.", icon=":material/upload_file:")
log = ToolLog()
env = BT.list_user_datasets()
log.calls.append(env)
rows = (env.get("data") or {}).get("datasets") or []
if env.get("status") != "ok" or not rows:
    unknown_box(env, "User datasets")
    st.markdown("Add one with the synthetic walkthrough (docs/WALKTHROUGH.md):\n\n"
                "```bash\nuv run python examples/byod_synthetic/make_example.py\n"
                "uv run measure-it byod validate examples/byod_synthetic/synthetic_demo\n"
                "uv run measure-it byod ingest examples/byod_synthetic/synthetic_demo\n"
                "uv run measure-it byod evaluate examples/byod_synthetic/synthetic_demo\n"
                "uv run measure-it byod deploy examples/byod_synthetic/synthetic_demo --demo\n```")
    limitations_panel(log, ["No user dataset is ingested on this machine; nothing below is filled in."])
    st.stop()

st.subheader("User datasets on this machine")
show_df(pd.DataFrame(rows)[["dataset_id", "title", "owner", "status", "synthetic", "demo", "comparator_type",
                            "primary_condition_id", "measurement_class", "n_cases_primary", "n_controls_primary",
                            "plan_sha256"]], key="byod_list")
ids = [r["dataset_id"] for r in rows]
ds = st.selectbox("Dataset", ids, key="byod_ds")
m_env = BT.get_user_dataset_metric(ds)
log.calls.append(m_env)
d = m_env.get("data") or {}
meta = d.get("dataset") or {}
flags = [badge("SYNTHETIC") if meta.get("synthetic") else "", badge("demo: excluded unless --demo")
         if meta.get("demo") else "", badge(f"comparator: {meta.get('comparator_type')}"),
         badge("tier 2.5 user_supplied_own_computation")]
st.markdown(" ".join(f for f in flags if f), unsafe_allow_html=True)
if meta.get("comparator_type") == "healthy":
    st.warning("Comparator: healthy controls. Performance against healthy controls is optimistic: it overstates "
               "real-world performance, where the question is this illness versus other causes of the same symptoms.")

tab_plan, tab_metrics, tab_rank, tab_fig = st.tabs(["1. Locked plan", "2. Metrics", "3. Ranking change",
                                                    "4. Heatmap and embedding"])
with tab_plan:
    lock = d.get("locked_plan")
    if isinstance(lock, dict):
        st.markdown(f"Locked {lock.get('locked_at')} before any outcome was computed; sha256 `{lock.get('sha256')}`"
                    + (f"; amends `{lock['amends'][:12]}`" if lock.get("amends") else ""))
        st.json(lock.get("plan", {}), expanded=False)
    else:
        st.caption(f"{UNKNOWN}: not evaluated yet (measure-it byod evaluate).")
with tab_metrics:
    models = d.get("models")
    if isinstance(models, list) and models:
        cols = ["model", "role", "n_cases", "n_controls", "auroc", "auroc_ci_low", "auroc_ci_high", "auprc",
                "op_sensitivity", "op_specificity", "perm_p", "n_perm", "decision"]
        show_df(pd.DataFrame(models)[[c for c in cols if c in pd.DataFrame(models).columns]], key="byod_models")
    else:
        st.caption(f"{UNKNOWN}: no evaluation.")
    rec = d.get("performance_record")
    if isinstance(rec, dict):
        st.markdown(f"**Performance record** `{rec.get('record_id')}`: {rec.get('target_condition')} x "
                    f"{rec.get('primary_class')}; AUROC {rec.get('auroc'):.3f} ({rec.get('auroc_ci_low'):.3f}-"
                    f"{rec.get('auroc_ci_high'):.3f}); sensitivity {rec.get('op_sensitivity'):.2f} at specificity "
                    f"{rec.get('op_specificity'):.2f}; measurement-only {rec.get('measurement_only')}.")
        caveat_list(rec.get("caveats") or [])
with tab_rank:
    rep = d.get("deploy_report")
    if isinstance(rep, dict):
        st.markdown(f"Deployed {rep.get('at')} (demo mode {rep.get('demo_mode')}). {rep.get('framing')}")
        cells = pd.DataFrame([{"condition": c["condition_id"], "measurement": c["measurement_id"],
                               "before": f"{c['before'].get('performance_status')} / "
                                         f"{c['before'].get('selected_record_id') or '-'}",
                               "after": f"{c['after'].get('performance_status')} / "
                                        f"{c['after'].get('selected_record_id') or '-'}", "why": c["why"]}
                              for c in rep.get("relevant_cells", [])])
        show_df(cells, key="byod_cells")
        for ch in rep.get("ranking_changes", []):
            st.markdown(f"**{ch['condition_id']} x {ch['measurement_id']} ({ch['geo_level']})**: regions ranked "
                        f"under evidence_weighted {ch['n_ranked_evidence_weighted_before']} -> "
                        f"{ch['n_ranked_evidence_weighted_after']}; best joint rank {ch['best_joint_rank_before']} -> "
                        f"{ch['best_joint_rank_after']}; equal-weight ranks unchanged: "
                        f"{ch['equal_weight_ranks_unchanged']}")
            show_df(pd.DataFrame(ch.get("top10_evidence_weighted_after") or []), key=f"top_{ch['geo_level']}_"
                                                                                     f"{ch['measurement_id']}")
    else:
        st.caption(f"{UNKNOWN}: not deployed yet (measure-it byod deploy).")
with tab_fig:
    out = K.results_path(ds)
    for name, cap in (("heatmap_measurement_evidence.png", "measurement_evidence by condition x measurement after "
                                                             "the deploy (red outline: the user record is selected)"),
                      ("embedding_dpv_pca.png", "the dataset's own Digital Phenotype Vector PCA (within-dataset; "
                                                "labels colour points only)")):
        p = out / name
        if p.exists():
            st.image(str(p), caption=cap)
        else:
            st.caption(f"{UNKNOWN}: {name} not generated (run measure-it byod deploy).")
    hp = out / "heatmap_measurement_evidence.csv"
    if hp.exists():
        with st.expander("Heatmap data (CSV)"):
            show_df(pd.read_csv(hp), key="byod_heat")

limitations_panel(log, [
    "User-supplied, local-only data: not public and not re-analysable by others; de-identification, labels and "
    "permission are attested by the data owner (docs/BRING_YOUR_OWN_DATA.md).",
    "Only the aggregate performance record crosses into the ranking; person rows never carry geography and are never "
    "joined to a place, facility or another dataset.",
    "A demo (SYNTHETIC) record enters rankings only after `measure-it byod deploy --demo`.",
    "Every ranked region is a candidate deployment opportunity for pilot evaluation, not a validated diagnostic "
    "pathway."])
