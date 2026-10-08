"""Measure It to Cure It (public-data alpha): dashboard home.

    uv run measure-it dashboard                 # or: uv run streamlit run dashboard/app.py

Five pages (dashboard/pages/): National Opportunity Map, Condition Explorer, Measurement Explorer, Geography Explorer,
Deployment Recommendation; plus page 6, Your data (user-supplied datasets, docs/BRING_YOUR_OWN_DATA.md). Every page calls the shared tool facade `measure_it.tools` in-process (the same functions
the MCP server and the FastAPI app expose); see docs/DASHBOARD.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from components.common import FRAMING, ToolLog, limitations_panel, ok, setup_page, show_df, unknown_box  # noqa: E402
from measure_it.mcp.server import DATA_LAYERS_MD  # noqa: E402

PAGES = [
    ("pages/1_National_Opportunity_Map.py", "National Opportunity Map", ":material/map:",
     "Burden (with evidence level), diagnostic deserts and candidate deployment opportunities by county or state; "
     "overlays for HRSA health centers, trial sites, NIH-funded organisations and specialist density."),
    ("pages/2_Condition_Explorer.py", "Condition Explorer", ":material/stethoscope:",
     "e.g. Long COVID: phenotype signatures incl. null results, wearable abnormalities from public cohorts, "
     "condition-level molecular enrichment, candidate measurements (five separate dimensions), FDA technologies, "
     "burden sources."),
    ("pages/3_Measurement_Explorer.py", "Measurement Explorer", ":material/sensors:",
     "e.g. wearable autonomic monitoring: conditions where it is studied, signals, public datasets, trials, FDA device "
     "classes, regions with unmet need."),
    ("pages/4_Geography_Explorer.py", "Geography Explorer", ":material/location_on:",
     "e.g. San Diego County, California: burden by condition with levels, vulnerability, providers, specialists, "
     "trials, NIH research, candidate technologies, nearby health centers."),
    ("pages/5_Deployment_Recommendation.py", "Deployment Recommendation", ":material/checklist:",
     "The demo query: a transparent ranked table with every component and rank interval, a map, candidate sites, "
     "evidence with trace_evidence, limitations and the recommendation JSON."),
    ("pages/6_Your_Data.py", "Your data", ":material/upload_file:",
     "User-supplied datasets (measure-it byod): manifest, locked plan, metrics, the performance record that reaches "
     "the ranking and how the ranking moved; local-only, aggregate results only."),
]


def home() -> None:
    setup_page("Measure It to Cure It: public-data alpha",
               "A modular measurement-validation and deployment engine for invisible illnesses, built only on public data.")
    log = ToolLog()

    st.markdown(
        "> Invisible illnesses are heterogeneous and poorly represented by diagnosis codes alone. Objective physiological "
        "measurements can reveal latent phenotypes. Public biomedical data can identify which measurements may be "
        "informative, while U.S. federal open data can identify the regions, facilities, and clinicians where those "
        "measurements could have the greatest deployment value.")

    steps = ["Person-level public data", "Objective wearable phenotype", "Condition ontology",
             "Condition-level molecular enrichment", "Candidate measurement / technology", "Geographic burden & unmet need",
             "Clinic / research-site readiness", "Candidate deployment opportunity"]
    cols = st.columns(len(steps))
    for i, (c, s) in enumerate(zip(cols, steps)):
        c.markdown(f"<div style='border:1px solid rgba(128,128,128,0.35);border-radius:8px;padding:8px 6px;"
                   f"min-height:118px;margin-bottom:8px;font-size:0.86rem'><span style='opacity:0.6'>{i + 1}</span><br>{s}</div>",
                   unsafe_allow_html=True)
    st.caption("The measurement is swappable: wearables today (WearableAdapter); a CAPRIO nailfold-capillaroscopy adapter "
               "later (CapillaroscopyAdapter stub, no private data) feeds the same ontology, geography, facility and "
               "agent layers.")

    st.subheader("Pages")
    for i in range(0, len(PAGES), 3):
        cs = st.columns(3)
        for c, (path, label, icon, desc) in zip(cs, PAGES[i:i + 3]):
            with c.container(border=True):
                try:
                    st.page_link(path, label=label, icon=icon)
                except Exception:  # noqa: BLE001 - page registry unavailable (e.g. a bare script run)
                    st.markdown(f"**{label}**")
                st.caption(desc)

    st.info(FRAMING, icon=":material/info:")

    st.subheader("Public sources in this build")
    env = log.call("list_sources")
    if ok(env):
        d = env["data"]
        m = st.columns(len(d["by_data_layer"]) + 1)
        m[0].metric("Sources", d["n_sources"])
        for c, (k, v) in zip(m[1:], sorted(d["by_data_layer"].items())):
            c.metric(k.replace("_", " "), v)
        src = pd.DataFrame(d["sources"])
        cols = [c for c in ["source_id", "name", "publisher", "data_layer", "status", "source_version", "retrieved_at",
                            "audit"] if c in src]
        with st.expander("SOURCE_REGISTRY.yaml (each source has data/raw/<source>/DATA_AUDIT.md)"):
            show_df(src[cols], key="home_sources")
    else:
        unknown_box(env, "Source registry")

    st.subheader("Four data layers, never joined as the same people")
    st.markdown(DATA_LAYERS_MD.split("\n", 2)[2])

    limitations_panel(log, [
        "Public-data alpha: every number on these pages comes from a stored table through the tool facade; nothing is "
        "generated by a language model. Missing data are shown as UNKNOWN / NOT AVAILABLE with the tool's reason.",
        "Candidate deployment opportunities rest on ecological joins, curated assumptions (configs/relevance.yaml) and "
        "proxy or inherited burden; the short list is weight-sensitive (results/SCORING_RESULTS.md).",
        "Person-level results (NHANES, Stanford wearable cohorts, MapMECFS) are group statistics from public cohorts; "
        "several labels are proxies (e.g. an ME/CFS-like symptom proxy), and some results are null.",
    ])


nav = st.navigation(
    [st.Page(home, title="Home", icon=":material/home:", default=True)]
    + [st.Page(path, title=label, icon=icon, url_path=Path(path).stem.split("_", 1)[1]) for path, label, icon, _ in PAGES])
nav.run()
