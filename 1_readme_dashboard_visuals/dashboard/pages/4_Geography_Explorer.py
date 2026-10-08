"""Page 4: Geography Explorer - one county or state: burden by condition with evidence levels, vulnerability,
provider density and relevant specialists, active/completed trials, NIH-funded research, candidate technologies and
candidate facilities, and nearby health centers on a map."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from components import charts as C  # noqa: E402
from components import data as D  # noqa: E402
from components import maps as M  # noqa: E402
from components.common import (EVIDENCE_LEVELS, ToolLog, fmt, is_unknown, limitations_panel, ok,  # noqa: E402
                               records_df, sections, setup_page, show_df, trace_explorer, unknown_box)
from measure_it.config import UNKNOWN  # noqa: E402

setup_page("Geography Explorer",
           "Type a county or state ('San Diego County, California', 'Cook County, IL', 'Mississippi', FIPS '06073'). "
           "Everything here describes the place and its registries (ecological), never a person.")
log = ToolLog()
stamp = D.stamp()
cat = D.catalog(stamp)
labels = cat["labels"]

c1, c2, c3, c4 = st.columns([2, 1.5, 1.6, 1])
text = c1.text_input("Geography", value="San Diego County, California", key="geo_text")
cond_ids = [c["id"] for c in cat["conditions"]]
cond = c2.selectbox("Condition (facilities, technologies)", cond_ids, index=cond_ids.index("long_covid"),
                    key="geo_condition", format_func=lambda c: labels.get(c, c))
meas_ids = list(cat["scored_measurements"])
meas = c3.selectbox("Measurement bundle", meas_ids, key="geo_measurement",
                    format_func=lambda m: cat["scored_measurements"].get(m, m))
radius = c4.slider("Radius (km)", 10, 150, 50, step=10, key="geo_radius")

ctx = log.call("get_geographic_context", geography_id=text)
if not ok(ctx):
    unknown_box(ctx, f"Geography '{text}'")
    limitations_panel(log, ["The geography did not resolve to one state, county or ZCTA; nothing below is guessed."])
    st.stop()

d = ctx["data"]
g = d["geography"]
gid, glevel = g["geo_id"], g["geo_level"]
pop = (d.get("population") or {}).get("total_population")
vul = d.get("vulnerability") or {}
rur = d.get("rurality") or {}
st.markdown(f"### {g.get('name')}  \n`geo:{gid}` | {glevel} | {g.get('state_name')}")
mc = st.columns(5)
mc[0].metric("Population (ACS 2020-2024)", fmt(pop, 0))
mc[1].metric("SVI overall (percentile)", fmt(vul.get("svi_overall")))
mc[2].metric("HRSA health-center sites", fmt((d.get("hrsa_health_centers") or {}).get("n_service_delivery_sites"), 0))
mc[3].metric("Trials with a site here", fmt((d.get("trials") or {}).get("distinct_relevant_trials_with_site_here"), 0))
mc[4].metric("NIH core projects", fmt((d.get("nih_projects") or {}).get("core_projects"), 0))
if rur.get("rucc_description"):
    st.caption(f"Rurality: RUCC {rur.get('rucc_2023')} ({rur.get('rucc_description')}). "
               f"Shortage designations: primary care {(d.get('shortage_designations') or {}).get('pc_coverage')}, "
               f"MUA {(d.get('shortage_designations') or {}).get('mua_coverage')}.")

tabs = sections(["Burden by condition", "Vulnerability & context", "Providers & specialists", "Trials & NIH research",
                "Candidate technologies & facilities", "Nearby health centers (map)"])

# ------------------------------------------------------------------------------------------------ burden
with tabs[0]:
    st.caption("Evidence levels: " + "; ".join(f"{k} = {v}" for k, v in EVIDENCE_LEVELS.items())
               + ". A state value attached to a county is inherited context, never county prevalence.")
    rows = []
    for b in d.get("burden_by_condition") or []:
        v = b.get("burden_value")
        unit, label = "", b.get("burden_measure_id")
        if not is_unknown(v) and glevel in ("state", "county"):   # unit + label of this geography's burden row
            be = log.call("get_condition_burden", condition=b["condition_id"], geography_level=glevel, geo_id=gid,
                          max_rows=1)
            br = ((be.get("data") or {}).get("rows") or [{}])[0] if ok(be) else {}
            unit, label = br.get("value_unit") or "", br.get("measure_label") or label
        val = fmt(v, 3) if isinstance(v, float) and abs(v) < 1 else fmt(v, 1)
        if b.get("burden_inherited") is True and not is_unknown(v):   # the value itself says it is a state value
            val = f"{val} ({g.get('state_abbr') or 'state'} STATE value, inherited)"
        modelled = "small-area" in str(b.get("burden_source_name") or "") and glevel == "county"
        if modelled and not is_unknown(v):                            # ... or a modelled county estimate
            val = f"{val} (modelled small-area estimate)"
        rows.append({"condition": labels.get(b["condition_id"], b["condition_id"]),
                     "burden": val, "unit": unit,
                     "evidence level": b.get("burden_evidence_level"),
                     "inherited state value": b.get("burden_inherited"),
                     "modelled estimate": modelled,
                     "source resolution": b.get("burden_source_resolution"),
                     "95% CI": (f"{b.get('burden_ci_low')}-{b.get('burden_ci_high')}"
                                if not is_unknown(b.get("burden_ci_low")) else ""),
                     "measure": label, "period": b.get("burden_period"),
                     "diagnostic desert": b.get("diagnostic_desert"),
                     "desert (access only)": b.get("diagnostic_desert_access_only"),
                     "why unknown": b.get("burden_unknown_reason") if is_unknown(v) else ""})
    bdf = pd.DataFrame(rows)
    if len(bdf):   # usable burden first (A, B, C), then level D
        bdf = bdf.sort_values(["evidence level", "condition"], key=lambda s: s.fillna("D") if s.name != "condition"
                              else s).reset_index(drop=True)
    show_df(bdf, key="geo_burden")
    n_mod = int(bdf["modelled estimate"].fillna(False).astype(bool).sum()) if len(bdf) else 0
    if n_mod and glevel == "county":
        st.warning(f"{n_mod} condition(s) here carry a **modelled** small-area estimate (BRFSS 2023 MRP), not "
                   "observed county prevalence: within-state differences come from county composition and "
                   "covariates only.", icon=":material/warning:")
    n_inh = int(bdf["inherited state value"].fillna(False).astype(bool).sum()) if len(bdf) else 0
    if n_inh and glevel == "county":
        st.warning(f"{n_inh} condition(s) here carry the **state** estimate inherited by the county (no within-state "
                   "variation; not county prevalence).")

# ------------------------------------------------------------------------------------------------ vulnerability
with tabs[1]:
    th = [("svi_overall", "overall"), ("svi_socioeconomic", "socioeconomic"), ("svi_household", "household"),
          ("svi_minority", "racial/ethnic minority status"), ("svi_housing_transport", "housing / transportation")]
    vdf = pd.DataFrame([{"theme": lbl, "percentile": vul.get(k)} for k, lbl in th
                        if isinstance(vul.get(k), (int, float))])
    if len(vdf):
        st.plotly_chart(C.hbar(vdf, "percentile", "theme", xtitle=f"CDC/ATSDR SVI percentile "
                                                                 f"({vul.get('svi_ranking_universe')}; release "
                                                                 f"{fmt(vul.get('svi_release'), 0)})",
                               text="percentile", height=240), key="geo_svi")
        st.caption("Vulnerability is social context, not disease prevalence.")
    else:
        st.caption(f"SVI: {UNKNOWN}")
    a, p = d.get("acs_context") or {}, d.get("places_context") or {}
    cc = st.columns(2)
    with cc[0]:
        st.markdown("**ACS context**")
        show_df(pd.DataFrame({"variable": list(a), "value": [fmt(v) for v in a.values()]}), key="geo_acs")
    with cc[1]:
        st.markdown(f"**CDC PLACES** ({d.get('places_derivation') or ''}; crude %)")
        show_df(pd.DataFrame({"measure": list(p), "value": [fmt(v) for v in p.values()]}), key="geo_places")

# ------------------------------------------------------------------------------------------------ providers
with tabs[2]:
    pr = records_df(d.get("providers_by_specialty_group"), ["label", "specialty_group", "group_type",
                                                             "n_individual_providers", "n_individual_physicians",
                                                             "n_organizations", "individual_providers_per_100k"])
    if len(pr):
        ind = pr[pr["group_type"].ne("facility")].sort_values("individual_providers_per_100k", ascending=False)
        ind = ind[pd.to_numeric(ind["individual_providers_per_100k"], errors="coerce") > 0]
        if len(ind):
            st.plotly_chart(C.hbar(ind.head(20), "individual_providers_per_100k", "label",
                                   xtitle="individual NPPES providers per 100k residents (self-reported taxonomy)",
                                   hover_cols=["n_individual_providers", "n_individual_physicians"]),
                            key="geo_providers")
    show_df(pr, key="geo_provider_table")
    st.markdown("**Relevant specialists per condition** (the condition's core specialty groups other than primary "
                "care; configs/relevance.yaml):")
    spec = pd.DataFrame([{"condition": labels.get(b["condition_id"], b["condition_id"]),
                          "specialists (excl. primary care)": b.get("relevant_specialists_n"),
                          "specialists per 100k": b.get("relevant_specialists_per_100k"),
                          "incl. primary care": b.get("relevant_providers_n"),
                          "incl. primary care per 100k": b.get("relevant_providers_per_100k"),
                          "specialty groups": str(b.get("relevant_specialty_groups") or "").replace(";", ", ")}
                         for b in d.get("burden_by_condition") or []])
    show_df(spec, key="geo_specialists", column_config={
        c: st.column_config.NumberColumn(format="%.1f") for c in ("specialists per 100k", "incl. primary care per 100k")})

# ------------------------------------------------------------------------------------------------ trials / NIH
with tabs[3]:
    t, n = d.get("trials") or {}, d.get("nih_projects") or {}
    tc = st.columns(4)
    tc[0].metric("Relevant trials with a site here", fmt(t.get("distinct_relevant_trials_with_site_here"), 0))
    tc[1].metric("Active", fmt(t.get("active"), 0))
    tc[2].metric("Recruiting", fmt(t.get("recruiting"), 0))
    tc[3].metric("Completed", fmt(t.get("completed"), 0))
    bc = t.get("by_condition") or {}
    if bc:
        tdf = pd.DataFrame([{"condition": labels.get(k, k), "trials": v} for k, v in bc.items()])
        st.plotly_chart(C.hbar(tdf.sort_values("trials", ascending=False), "trials", "condition",
                               xtitle=f"trials with a site here ({t.get('precision_view')})", text="trials"),
                        key="geo_trials")
    if t.get("example_nct_ids"):
        st.caption("Examples: " + ", ".join(t["example_nct_ids"][:10]))
    st.markdown(f"**NIH-funded research** ({n.get('precision_view')}): {fmt(n.get('distinct_projects_appl_id'), 0)} "
                f"award records, {fmt(n.get('core_projects'), 0)} core projects ({fmt(n.get('active_core_projects'), 0)} "
                f"active), {fmt(n.get('organizations'), 0)} organisations.")
    show_df(records_df(n.get("top_organizations")), key="geo_nih")

# ------------------------------------------------------------------------------------------------ technologies / facilities
clin = log.call("find_candidate_clinics", geography_id=gid, condition=cond, measurement=meas,
                radius_km=float(radius), max_sites=15)
with tabs[4]:
    disc = log.call("discover_candidate_measurements", condition=cond, top_n=8)
    st.markdown(f"**Candidate technologies for {labels.get(cond, cond)}** (condition-level, from "
                "discover_candidate_measurements; the same list applies to every geography):")
    if ok(disc):
        show_df(pd.DataFrame([{
            "rank": c.get("rank"), "measurement": c.get("name"),
            "evidence tier": ((c.get("dimensions") or {}).get("measurement_evidence_strength") or {}).get("tier"),
            "technology maturity": ((c.get("dimensions") or {}).get("technology_maturity") or {}).get("tier"),
            "regulatory visibility": ((c.get("dimensions") or {}).get("regulatory_visibility") or {}).get("tier"),
            "deployment complexity": ((c.get("dimensions") or {}).get("deployment_complexity") or {}).get("value"),
            "phenotype signal": ((c.get("dimensions") or {}).get("phenotype_signal_strength") or {}).get("status")}
            for c in disc["data"].get("candidates") or []]), key="geo_tech")
    else:
        unknown_box(disc, "Candidate technologies")
    st.markdown(f"**Candidate facilities within {radius} km** (find_candidate_clinics): facilities with "
                "characteristics suggesting they may be viable implementation or study partners; six characteristics "
                "ranked separately, never one score.")
    if ok(clin):
        cd = clin["data"]
        st.caption(f"Pool: {cd.get('n_facilities_in_pool')} facilities, {cd.get('n_eligible')} eligible. With each "
                   "characteristic: " + ", ".join(f"{k} {v}" for k, v in (cd.get("n_with_characteristic") or {}).items()))
        crow = []
        for c in cd.get("candidates") or []:
            ch = c.get("characteristics") or {}
            crow.append({"facility": c.get("facility_name"), "kind": c.get("primary_kind"), "city": c.get("city"),
                         "km": c.get("distance_km"), "selected via": c.get("selected_via"),
                         **{k: f"{v.get('value')} (rank {v.get('rank') if v.get('rank') is not None else '-'})"
                            for k, v in ch.items() if isinstance(v, dict)},
                         "object_id": c.get("object_id")})
        show_df(pd.DataFrame(crow), key="geo_clinics")
        with st.expander("Reasons per facility"):
            for c in cd.get("candidates") or []:
                st.markdown(f"**{c.get('facility_name')}**: " + " ".join(map(str, c.get("reasons") or [])))
        if cd.get("absence_note"):
            st.info(cd["absence_note"])
    else:
        unknown_box(clin, "Candidate facilities")

# ------------------------------------------------------------------------------------------------ map
with tabs[5]:
    lat, lon = g.get("lat"), g.get("lon")
    if glevel == "county":
        area = D.subset_geojson(D.geojson("county", 0.002, gid[:2]), [gid])
        context = D.geojson("county", 0.005)       # neighbours drawn from the same boundaries (no basemap coastline)
        near = D.hrsa_near(lat, lon, float(radius), stamp, county_fips=gid)
        near_note = f"HRSA sites inside the county or within {radius} km of its internal point"
    else:
        area = D.subset_geojson(D.geojson("state", 0.005), [gid])
        context = D.geojson("state", 0.005)
        near = D.hrsa_points(stamp)
        near = near[near["state_abbr"].eq(g.get("state_abbr"))].reset_index(drop=True)
        near_note = "HRSA sites in the state"
    log.table("facilities__hrsa", "HRSA health-center sites near the geography (source lat/lon)")
    pts = [(near.rename(columns={"facility_name": "site", "health_center_type": "type"}), "hrsa",
            ["site", "city", "type", "location_setting", "distance_km", "object_id"] if "distance_km" in near
            else ["site", "city", "type", "location_setting", "object_id"])]
    if ok(clin):
        cl = pd.DataFrame(clin["data"].get("candidates") or [])
        if len(cl) and {"lat", "lon"} <= set(cl.columns):
            pts.append((cl.rename(columns={"facility_name": "facility"}), "clinics",
                        ["facility", "city", "selected_via", "distance_km", "object_id"]))
    st.plotly_chart(M.local_map(area, pts, height=560, area_hover=g.get("name", ""), context_fc=context),
                    key="geo_map")
    st.caption(f"{len(near)} {near_note} (orange circles); candidate facilities from find_candidate_clinics (violet "
               "squares). Facility locations are the source's coordinates or ZIP/ZCTA centroids.")
    show_df(near.head(40)[[c for c in ["facility_name", "city", "health_center_type", "location_setting",
                                       "distance_km", "object_id"] if c in near]], key="geo_hrsa_table")

ids = [f"geo:{gid}"]
for e in log.calls:
    ids += (e.get("provenance") or {}).get("object_ids") or []
trace_explorer(ids, key="geo_trace", log=log)

limitations_panel(log, [
    "Ecological context: burden, vulnerability, providers, trials and grants describe the place and its registries; "
    "no individual's location or condition is inferred.",
    "Burden rows carry their evidence level; long COVID county values are the state Household Pulse Survey estimate "
    "(inherited, not county prevalence); level-D conditions have no usable burden.",
    "Provider counts are self-reported NPPES taxonomies at ZIP-centroid locations; a specialty does not mean a "
    "provider evaluates or treats a condition.",
    "Candidate technologies are condition-level (not specific to this geography); candidate facilities are "
    "implementation or study-partner candidates with separate characteristic ranks, not a quality ranking.",
])
