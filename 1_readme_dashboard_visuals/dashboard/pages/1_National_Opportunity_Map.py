"""Page 1: National Opportunity Map - burden (with evidence level), diagnostic desert or candidate deployment
opportunity as a county/state choropleth, with independent overlays for HRSA health-center sites, trial sites, NIH-funded
organisations and specialist density."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from components import data as D  # noqa: E402
from components import maps as M  # noqa: E402
from components.common import (FRAMING, ToolLog, burden_caveat_box, limitations_panel, ok, setup_page,  # noqa: E402
                               show_df, unknown_box)

LAYERS = ["Deployment opportunity", "Diagnostic desert", "Burden (with evidence level)"]

setup_page("National Opportunity Map",
           "Where burden, measurement deserts and candidate deployment opportunities sit, with the facility and research "
           "infrastructure that could act on them. Choose a layer, then toggle the overlays independently.")
log = ToolLog()
stamp = D.stamp()
cat = D.catalog(stamp)
labels = cat["labels"]

# ------------------------------------------------------------------------------------------------ controls
scored = cat["scored_conditions"]
cond_opts = scored + [c["id"] for c in cat["conditions"] if c["id"] not in scored]
c1, c2, c3, c4, c5 = st.columns([1.1, 1.6, 1.5, 1.5, 1.2])
level = c1.radio("Geography level", ["county", "state"], horizontal=True, key="map_level")
layer = c2.radio("Choropleth", LAYERS, key="map_layer")
cond = c3.selectbox("Condition or condition set", cond_opts, key="map_condition",
                    index=cond_opts.index("long_covid_or_me_cfs") if "long_covid_or_me_cfs" in cond_opts else 0,
                    format_func=lambda c: labels.get(c, c) + ("" if c in scored else "  (burden / desert only)"))
meas_ids = list(cat["scored_measurements"])
meas = c4.selectbox("Measurement bundle", meas_ids, key="map_measurement",
                    format_func=lambda m: cat["scored_measurements"].get(m, m))
ws_names = list(cat["weight_sets"])
ws = c5.selectbox("Weight set", ws_names, key="map_weights", index=ws_names.index("equal") if "equal" in ws_names else 0)

o1, o2, o3, o4, o5 = st.columns(5)
show_hrsa = o1.checkbox("HRSA health-center sites", value=False, key="ov_hrsa")
show_trials = o2.checkbox("Trial sites", value=True, key="ov_trials")
active_only = o3.checkbox("Trial sites: active/recruiting only", value=False, key="ov_trials_active",
                          disabled=not show_trials)
show_nih = o4.checkbox("NIH-funded organisations", value=False, key="ov_nih")
show_spec = o5.checkbox("Specialist density", value=False, key="ov_spec")

members = cat["members"].get(cond, [cond])
w = cat["weight_sets"][ws]
st.caption(f"Weights ({ws}): " + ", ".join(f"{k} {v:.2f}" for k, v in w.items())
           + f". Measurement adapter: {cat['measurement_adapter_status'].get(meas, 'n/a')}.")

# ------------------------------------------------------------------------------------------------ caveats (always)
burden_caveat_box(log, members, level, labels)

# ------------------------------------------------------------------------------------------------ choropleth data
fc = D.geojson(level, 0.02 if level == "county" else None)
state_fc = D.geojson("state") if level == "county" else None
highlight, table, value_col, title = [], pd.DataFrame(), None, ""
if layer == "Deployment opportunity":
    if cond not in scored:
        st.warning(f"**Deployment opportunity: UNKNOWN / NOT AVAILABLE** for {labels.get(cond, cond)}: "
                   "deployment_opportunities is precomputed for " + ", ".join(labels.get(c, c) for c in scored)
                   + ". The Deployment Recommendation page computes other combinations on the fly with the same "
                     "engine.")
        df = pd.DataFrame(columns=["geo_id"])
    else:
        d = D.opportunity_rows(cond, meas, level, stamp)
        log.table("deployment_opportunities", f"composite_{ws} / rank_{ws} per geography for the choropleth")
        rank = d[f"rank_{ws}"]
        n_ranked = int(rank.notna().sum())
        comp = d[f"composite_{ws}"].where(rank.notna())
        why = np.where(d["burden_incomplete"].fillna(False).astype(bool), "not ranked: incomplete member burden",
                       np.where(d["small_population_flag"].fillna(False).astype(bool),
                                f"not ranked: population < {cat['min_population']:,}",
                                "not ranked under this weight set"))
        df = pd.DataFrame({
            "geo_id": d["geo_id"], "Region": d["geo_name"], "value": comp,
            "Rank": [f"{r:.0f} of {n_ranked}" if pd.notna(r) else "not ranked" for r in rank],
            "Composite": d[f"composite_{ws}"],
            "Rank interval 5-95 (equal-weight MC)": [f"{a:.0f}-{b:.0f}" if pd.notna(a) else "n/a"
                                                     for a, b in zip(d["rank_mc_p05"], d["rank_mc_p95"])],
            "P(top 10, equal-weight MC)": d["p_top10_mc"],
            "Burden level": d["burden_evidence_level"].fillna("D") + np.where(
                d["burden_inherited"].fillna(False).astype(bool), " (inherited state value)", ""),
            "Burden pct": d["burden_pct"], "SVI pct": d["vulnerability_pct"], "Desert pct": d["diagnostic_desert_pct"],
            "Capacity pct": d["clinic_capacity_pct"], "Readiness pct": d["research_readiness_pct"],
            "Why no value": why, "rank_num": rank, "object_id": d["object_id"]})
        value_col, title = "value", f"Composite ({ws})"
        env = log.call("rank_deployment_opportunities", condition=cond, measurement=meas, geography_level=level,
                       top_n=10, weight_set=ws)
        if ok(env):
            highlight = [r["geography"]["fips"] for r in env["data"]["recommendations"]]
        else:
            unknown_box(env, "Top-10 candidate deployment opportunities")
        # component percentiles right after the composite, so they are visible without scrolling
        table = df[df["rank_num"].notna()].sort_values("rank_num")[
            ["Rank", "Region", "Composite", "Burden pct", "SVI pct", "Desert pct", "Capacity pct", "Readiness pct",
             "Rank interval 5-95 (equal-weight MC)", "P(top 10, equal-weight MC)", "Burden level", "geo_id",
             "object_id"]]
elif layer == "Diagnostic desert":
    d = D.condition_layer(cond, level, stamp)
    log.table("geo_condition_features" if cond not in {s["id"] for s in cat["sets"]} else "deployment_opportunities",
              "diagnostic_desert per geography")
    variant = d["desert_variant"].iloc[0] if len(d) else ""
    df = pd.DataFrame({"geo_id": d.get("geo_id", []), "Region": d.get("geo_name", []), "value": d.get("desert", []),
                       "Desert index": d.get("desert", []), "Variant": d.get("desert_variant", []),
                       "Burden level": d.get("burden_level", []),
                       "Inherited state value": d.get("inherited", []), "Burden": d.get("burden_detail", [])})
    value_col, title = "value", "Diagnostic desert"
    st.caption(f"Diagnostic desert = high burden + high vulnerability + few relevant providers + few relevant trials "
               f"(geography module; raw components in geo_condition_features). Variant shown: {variant}.")
    table = df.dropna(subset=["value"]).sort_values("value", ascending=False).drop(columns=["value"])
else:
    d = D.condition_layer(cond, level, stamp)
    log.table("geo_condition_features" if cond not in {s["id"] for s in cat["sets"]} else "deployment_opportunities",
              "burden value, evidence level, inherited flag and source resolution per geography")
    kind = d["burden_kind"].iloc[0] if len(d) else "burden"
    inherited_share = float(d["inherited"].mean()) if len(d) else 0.0
    if level == "county" and inherited_share > 0:
        # a state value carried by every county of the state: label it on the colour bar itself, never as prevalence
        kind = ("inherited STATE value on counties" + kind.removeprefix("burden value")
                if inherited_share > 0.5 and kind.startswith("burden value")
                else f"{kind}; incl. inherited state value")
    df = pd.DataFrame({"geo_id": d.get("geo_id", []), "Region": d.get("geo_name", []), "value": d.get("burden", []),
                       "Burden": d.get("burden", []), "Evidence level": d.get("burden_level", []),
                       "Inherited state value": d.get("inherited", []), "Source resolution": d.get("resolution", []),
                       "Detail": d.get("burden_detail", []), "Measure": d.get("measure_label", []),
                       "Why no value": d.get("unknown_reason", [])})
    modelled_share = (float(d["measure_label"].astype(str).str.contains("modelled").mean())
                      if len(d) and "measure_label" in d else 0.0)
    if level == "county" and modelled_share > 0.5 and kind.startswith("burden value"):
        kind = "MODELLED small-area estimate on counties" + kind.removeprefix("burden value")
    value_col, title = "value", kind
    table = df.dropna(subset=["value"]).sort_values("value", ascending=False).drop(columns=["value", "Why no value"])
    if level == "county" and modelled_share > 0:
        st.warning(f"**{labels.get(cond, cond)} at county level:** {modelled_share:.0%} of counties carry a "
                   "**modelled** small-area estimate (BRFSS 2023 multilevel regression and post-stratification): "
                   "within-state differences come from county composition and covariates only. The colours are "
                   "**not observed county prevalence**.", icon=":material/warning:")
    if level == "county" and inherited_share > 0:
        st.warning(f"**{labels.get(cond, cond)} at county level:** {inherited_share:.0%} of counties carry their "
                   "**state** estimate (inherited context: identical for every county of a state, no within-state "
                   "variation). The colours and the table below are **not county prevalence**.",
                   icon=":material/warning:")

# ------------------------------------------------------------------------------------------------ figure
hover = [c for c in df.columns if c not in ("geo_id", "value", "rank_num", "object_id", "Why no value")]
if value_col is not None and len(df):
    hv = df.copy()
    if "Detail" in hv:
        hv["Detail"] = "<br>" + hv["Detail"].astype(str)
    if "Burden" in hv and hv["Burden"].dtype == object:
        hv["Burden"] = "<br>" + hv["Burden"].astype(str)
    fig = M.choropleth_map(fc, hv, value_col, hover, title,
                           nodata_hover=["Region", "Why no value"] if "Why no value" in hv else ["Region"],
                           state_fc=state_fc, highlight=highlight)
else:
    fig = M.choropleth_map(fc, pd.DataFrame({"geo_id": [], "value": []}), "value", ["geo_id"], title or "",
                           state_fc=state_fc)

if show_hrsa:
    h = D.hrsa_points(stamp)
    log.table("facilities__hrsa", "HRSA health-center service delivery sites (source lat/lon)")
    M.add_points(fig, h.rename(columns={"facility_name": "Site", "health_center_type": "type"}), "hrsa",
                 ["Site", "city", "state_abbr", "type", "object_id"], base_size=3.5)
if show_trials:
    tp = D.trial_site_points(tuple(members), active_only, stamp)
    log.table("trial_sites + trial_conditions + clinical_trials", "U.S. sites of trials registered for the condition")
    M.add_points(fig, tp.rename(columns={"city": "City", "n_trials": "trials", "n_recruiting": "recruiting sites",
                                         "example_ids": "e.g."}), "trials",
                 ["City", "state", "trials", "recruiting sites", "e.g."], size_col="trials", base_size=4,
                 max_size=14)
if show_nih:
    npnt = D.nih_org_points(tuple(members), stamp)
    log.table("nih_projects + nih_project_conditions", "NIH awardee organisations (title/abstract condition match)")
    M.add_points(fig, npnt.rename(columns={"org_name": "Organisation", "n_core_projects": "core projects",
                                           "n_active_core_projects": "active core projects"}), "nih",
                 ["Organisation", "city", "state", "core projects", "active core projects"], size_col="core projects",
                 base_size=6, max_size=22)
spec_groups = []
if show_spec:
    sp, spec_groups = D.specialist_density(tuple(members), level, stamp)
    log.table("providers (via scoring.opportunity.provider_counts_for_groups)",
              "individual NPIs in the condition's core specialty groups other than primary care, per 100k")
    sp = sp[sp["per_100k"].notna() & (sp["n"] > 0)]
    M.add_points(fig, sp.rename(columns={"geo_name": "Region", "n": "specialists", "per_100k": "per 100k"}),
                 "specialists", ["Region", "specialists", "per 100k"], size_col="per 100k",
                 base_size=2.5 if level == "county" else 6, max_size=14 if level == "county" else 30,
                 name="Specialist density (per 100k)")

st.plotly_chart(fig, key="national_map", config={"displaylogo": False, "scrollZoom": False})
legend_bits = ["grey = no value (UNKNOWN / NOT AVAILABLE; hover gives the reason)"]
if highlight:
    legend_bits.append("red outline = top 10 under the selected weights (rank_deployment_opportunities)")
if show_spec and spec_groups:
    legend_bits.append("specialist groups: " + ", ".join(spec_groups))
st.caption("; ".join(legend_bits) + ". Facility and trial locations are ZIP/ZCTA or city centroids.")

# ------------------------------------------------------------------------------------------------ table
st.subheader("Ranked table" if layer == "Deployment opportunity" else "Highest values")
if layer == "Deployment opportunity":
    st.markdown(FRAMING)
show_df(table.head(50), key="map_table", column_config={
    c: st.column_config.NumberColumn(format="%.3f" if c == "Composite" else "%.2f")
    for c in table.columns if c == "Composite" or c.endswith("pct") or c.startswith("P(top")})
if len(table):
    st.download_button("Download this layer (CSV)", table.to_csv(index=False).encode(),
                       file_name=f"{layer.split()[0].lower()}_{cond}_{meas if layer.startswith('Deploy') else ''}"
                                 f"_{level}_{ws}.csv".replace("__", "_"), mime="text/csv", key="map_csv")

limitations_panel(log, [
    "Ecological map: every value describes a place or a registry, never a person; no participant of any public cohort "
    "is located.",
    "Burden layers show the primary measure with its evidence level (A direct, B coded condition, C proxy, D none). A "
    "state value shown on counties is inherited context, not county prevalence; a condition set takes its least "
    "direct member's level.",
    "The composite is a weighted mean of percentile-normalised components; the rank interval comes from Monte Carlo "
    "draws of burden and weights around the equal weights (docs/SCORING.md). Which regions make a short list depends "
    "on the weights (results/SCORING_RESULTS.md, Test 7).",
    "Trial sites are ClinicalTrials.gov city-level geopoints of trials registered for the condition (literal condition "
    "match); NIH organisations are RePORTER awardees with title/abstract matches. Both are research-activity signals, "
    "not evidence that a measurement works.",
    "Specialist density counts self-reported NPPES taxonomies in curated specialty groups (configs/relevance.yaml); a "
    "taxonomy does not mean a clinician evaluates or treats the condition.",
    "Boundaries are simplified for display (measure_it.api.geo).",
])
