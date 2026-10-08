"""Cached, read-only data for the dashboard's map layers and selectors.

The per-entity answers on every page come from `measure_it.tools`. A choropleth or a national point layer needs every
county / site at once, which the tools deliberately cap (<= 100 list items), so these loaders read the same processed
tables the tools read (deployment_opportunities, geo_condition_features, facilities__hrsa, trial_sites,
nih_projects, facilities) and the simplified boundaries of `measure_it.api.geo`. They select and aggregate rows; they
never compute a score. Every cached function takes the tool facade's data stamp, so a rebuilt table invalidates it.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import streamlit as st

from measure_it import tools as T
from measure_it.config import PROCESSED, UNKNOWN, load_config

SEP = "|"


def stamp() -> float:
    return T.data_stamp()


def _read(name: str, columns: list[str] | None = None, filters=None) -> pd.DataFrame:
    path = PROCESSED / f"{name}.parquet"
    if not path.exists():
        return pd.DataFrame(columns=columns or [])
    if columns:
        have = set(pq.ParquetFile(path).schema.names)
        columns = [c for c in columns if c in have]
    return pq.read_table(path, columns=columns, filters=filters).to_pandas()


# ------------------------------------------------------------------------------------------------ catalog

@st.cache_data(show_spinner=False)
def catalog(stamp: float) -> dict:
    """Selector options: registry conditions, condition sets, scored combinations, measurements, weight sets."""
    from measure_it.scoring import opportunity as O
    conds = [{"id": c["id"], "label": c["preferred_name"], "kind": "condition"}
             for c in load_config("conditions")["conditions"]]
    sets = [{"id": k, "label": v["label"], "kind": "condition_set", "members": v["members"]}
            for k, v in O.condition_sets().items()]
    scored = _read("deployment_opportunities", ["condition_id", "condition_label", "measurement_id",
                                                "measurement_label", "measurement_adapter_status"])
    scored = scored.drop_duplicates(["condition_id", "measurement_id"])
    meas = _read("measurement_registry", ["measurement_id", "name", "kind", "adapter", "adapter_status"])
    return {
        "conditions": conds, "sets": sets,
        "labels": {**{c["id"]: c["label"] for c in conds}, **{s["id"]: s["label"] for s in sets}},
        "members": {**{c["id"]: [c["id"]] for c in conds}, **{s["id"]: list(s["members"]) for s in sets}},
        "scored_conditions": list(dict.fromkeys(scored["condition_id"])),
        "scored_measurements": {r.measurement_id: r.measurement_label
                                for r in scored.drop_duplicates("measurement_id").itertuples()},
        "measurement_adapter_status": {r.measurement_id: r.measurement_adapter_status
                                       for r in scored.drop_duplicates("measurement_id").itertuples()},
        "measurements": meas.to_dict("records"),
        "weight_sets": O.weight_sets(),
        "min_population": O.min_population(),
    }


# ------------------------------------------------------------------------------------------------ choropleth layers

OPP_BASE = ["object_id", "geo_id", "geo_name", "state_abbr", "population_total", "small_population_flag",
            "rank_eligible", "burden_value", "burden_value_unit", "burden_measure_label", "burden_pct",
            "burden_evidence_level", "burden_inherited", "burden_source_resolution", "burden_incomplete",
            "burden_incomplete_reason", "burden_excluded_level_D", "member_components", "svi_overall",
            "vulnerability_pct", "diagnostic_desert", "diagnostic_desert_pct", "diagnostic_desert_variant",
            "clinic_capacity_pct", "research_readiness_pct", "technology_saturation_pct", "rank_mc_p05",
            "rank_mc_p95", "p_top10_mc", "measurement_adapter_status",
            # metric -> translation link (docs/ANALYSIS_PLAN_METRIC_LINK.md)
            "measurement_performance_status", "measurement_performance_tier", "measurement_evidence", "expected_yield",
            "expected_yield_basis", "expected_yield_index", "expected_detectable_cases", "expected_false_positives",
            "expected_false_positives_upper", "reach", "evidence_weighted_status", "rank_mc_ew_p05", "rank_mc_ew_p95"]


@st.cache_data(show_spinner=False)
def opportunity_rows(condition_id: str, measurement_id: str, level: str, stamp: float) -> pd.DataFrame:
    """deployment_opportunities rows of one precomputed (condition or set) x measurement x level combination."""
    ws = catalog(stamp)["weight_sets"]
    cols = OPP_BASE + [f"composite_{w}" for w in ws] + [f"rank_{w}" for w in ws]
    return _read("deployment_opportunities", cols, filters=[("condition_id", "==", condition_id),
                                                           ("measurement_id", "==", measurement_id),
                                                           ("geo_level", "==", level)])


def _member_text(mc: str | None) -> str:
    try:
        d = json.loads(mc) if isinstance(mc, str) else {}
    except ValueError:
        return ""
    out = []
    for cid, m in d.items():
        v = m.get("burden_value")
        lvl = m.get("burden_evidence_level") or "D"
        inh = " inherited STATE value" if m.get("burden_inherited") else ""
        out.append(f"{cid}: {UNKNOWN if v is None else v} (level {lvl}{inh})")
    return "<br>".join(out)


@st.cache_data(show_spinner=False)
def condition_layer(condition_id: str, level: str, stamp: float) -> pd.DataFrame:
    """Burden (with evidence level, inherited flag, source resolution) and diagnostic desert per geography.

    Single conditions: geo_condition_features. Condition sets: deployment_opportunities (the set's burden is the
    percentile of its members' mean burden percentile; the desert is the mean of the members' deserts)."""
    cat = catalog(stamp)
    if condition_id in {s["id"] for s in cat["sets"]}:
        mid = next(iter(cat["scored_measurements"]), None)
        d = opportunity_rows(condition_id, mid, level, stamp) if mid else pd.DataFrame()
        if not len(d):
            return pd.DataFrame()
        out = pd.DataFrame({
            "geo_id": d["geo_id"], "geo_name": d["geo_name"], "burden": d["burden_pct"],
            "burden_kind": "set burden percentile (0-1)",
            "burden_level": d["burden_evidence_level"].fillna("D"),
            "inherited": d["burden_inherited"].fillna(False).astype(bool),
            "resolution": d["burden_source_resolution"].fillna(UNKNOWN),
            "burden_detail": d["member_components"].map(_member_text),
            "unknown_reason": np.where(d["burden_incomplete"].fillna(False).astype(bool),
                                       d["burden_incomplete_reason"].fillna("incomplete member burden"), ""),
            "desert": d["diagnostic_desert"], "desert_variant": d["diagnostic_desert_variant"].fillna(""),
            "measure_label": "mean of member burden percentiles"})
        return out.reset_index(drop=True)
    g = _read("geo_condition_features",
              ["geo_id", "geo_name", "burden_value", "burden_value_unit", "burden_measure_label",
               "burden_evidence_level", "burden_inherited", "burden_source_resolution", "burden_ci_low",
               "burden_ci_high", "burden_period", "burden_unknown_reason", "diagnostic_desert",
               "diagnostic_desert_access_only", "desert_unknown_reason", "in_analysis_universe"],
              filters=[("condition_id", "==", condition_id), ("geo_level", "==", level)])
    if not len(g):
        return pd.DataFrame()
    desert = pd.to_numeric(g["diagnostic_desert"], errors="coerce")
    access = pd.to_numeric(g["diagnostic_desert_access_only"], errors="coerce")
    use_access = desert.isna().all() and access.notna().any()
    unit = g["burden_value_unit"].dropna().astype(str)
    ci = [f" (95% CI {lo}-{hi})" if pd.notna(lo) and pd.notna(hi) else ""
          for lo, hi in zip(g["burden_ci_low"], g["burden_ci_high"])]
    out = pd.DataFrame({
        "geo_id": g["geo_id"], "geo_name": g["geo_name"],
        "burden": pd.to_numeric(g["burden_value"], errors="coerce"),
        "burden_kind": f"burden value ({unit.iloc[0] if len(unit) else 'unit unknown'})",
        "burden_level": g["burden_evidence_level"].fillna("D"),
        "inherited": g["burden_inherited"].fillna(False).astype(bool),
        "resolution": g["burden_source_resolution"].fillna(UNKNOWN),
        "burden_detail": [f"{v}{c} {p or ''}".strip() if pd.notna(v) else ""
                          for v, c, p in zip(g["burden_value"], ci, g["burden_period"])],
        "unknown_reason": g["burden_unknown_reason"].fillna(""),
        "desert": access if use_access else desert,
        "desert_variant": ("access-only (no burden term: burden level D)" if use_access
                           else "burden-informed diagnostic desert"),
        "measure_label": g["burden_measure_label"].fillna(UNKNOWN)})
    return out.reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ point layers

@st.cache_data(show_spinner=False)
def hrsa_points(stamp: float) -> pd.DataFrame:
    """HRSA health-center service delivery sites with coordinates (facility layer)."""
    h = _read("facilities__hrsa", ["site_bphc_number", "facility_name", "city", "state_abbr", "county_fips",
                                   "health_center_type", "location_setting", "is_service_delivery_site", "is_active",
                                   "lat", "lon", "geocode_method"])
    h = h[h["is_service_delivery_site"].fillna(False).astype(bool) & h["lat"].notna() & h["lon"].notna()]
    h = h.assign(object_id="hrsa_site:" + h["site_bphc_number"].astype(str))
    return h.drop(columns=["is_service_delivery_site"]).reset_index(drop=True)


@st.cache_data(show_spinner=False)
def trial_site_points(members: tuple[str, ...], active_only: bool, stamp: float) -> pd.DataFrame:
    """U.S. sites of trials registered for the condition(s) (literal condition match), aggregated to the site's
    geocode (ClinicalTrials.gov city-level geopoints)."""
    tc = _read("trial_conditions", ["nct_id", "condition_id", "condition_literal_match"],
               filters=[("condition_id", "in", list(members))])
    ncts = sorted(set(tc.loc[tc["condition_literal_match"].fillna(False).astype(bool), "nct_id"]))
    if not ncts:
        return pd.DataFrame(columns=["lat", "lon", "n_trials"])
    ct = _read("clinical_trials", ["nct_id", "overall_status", "brief_title"], filters=[("nct_id", "in", ncts)])
    if active_only:
        ct = ct[ct["overall_status"].isin(["RECRUITING", "NOT_YET_RECRUITING", "ENROLLING_BY_INVITATION",
                                           "ACTIVE_NOT_RECRUITING"])]
    s = _read("trial_sites", ["nct_id", "facility", "city", "state_abbr", "country_group", "lat", "lon"],
              filters=[("nct_id", "in", sorted(ct["nct_id"]))])
    s = s[(s["country_group"] == "US") & s["lat"].notna() & s["lon"].notna()].merge(ct, on="nct_id", how="left")
    if not len(s):
        return pd.DataFrame(columns=["lat", "lon", "n_trials"])
    s["recruiting"] = s["overall_status"].eq("RECRUITING")
    g = s.groupby(["lat", "lon"], as_index=False).agg(
        city=("city", "first"), state=("state_abbr", "first"), n_trials=("nct_id", "nunique"),
        n_sites=("nct_id", "size"), n_recruiting=("recruiting", "sum"),
        example_ids=("nct_id", lambda x: ", ".join(sorted(set(x))[:4])))
    return g.sort_values("n_trials", ascending=False).reset_index(drop=True)


@st.cache_data(show_spinner=False)
def nih_org_points(members: tuple[str, ...], stamp: float) -> pd.DataFrame:
    """NIH RePORTER awardee organisations with projects whose title/abstract names the condition(s) (likely false
    positives excluded, as the geography module's precision view), aggregated per organisation location."""
    pc = _read("nih_project_conditions", ["appl_id", "condition_id", "matched_in_title", "matched_in_abstract",
                                          "likely_false_positive"], filters=[("condition_id", "in", list(members))])
    keep = (pc["matched_in_title"].fillna(False).astype(bool) | pc["matched_in_abstract"].fillna(False).astype(bool)) \
        & ~pc["likely_false_positive"].fillna(False).astype(bool)
    appl = sorted(set(pc.loc[keep, "appl_id"].astype(str)))
    if not appl:
        return pd.DataFrame(columns=["lat", "lon", "n_core_projects"])
    p = _read("nih_projects", ["appl_id", "core_project_num", "org_name", "org_city", "org_state", "lat", "lon",
                               "is_active", "fiscal_year"])
    p = p[p["appl_id"].astype(str).isin(appl) & p["lat"].notna() & p["lon"].notna()]
    if not len(p):
        return pd.DataFrame(columns=["lat", "lon", "n_core_projects"])
    p = p.assign(active_core=np.where(p["is_active"].fillna(False).astype(bool), p["core_project_num"], None))
    g = p.groupby(["org_name", "lat", "lon"], as_index=False).agg(
        city=("org_city", "first"), state=("org_state", "first"), n_core_projects=("core_project_num", "nunique"),
        n_active_core_projects=("active_core", "nunique"), n_award_records=("appl_id", "nunique"),
        fy_min=("fiscal_year", "min"), fy_max=("fiscal_year", "max"))
    return g.sort_values("n_core_projects", ascending=False).reset_index(drop=True)


@st.cache_data(show_spinner=False)
def specialist_density(members: tuple[str, ...], level: str, stamp: float) -> tuple[pd.DataFrame, list[str]]:
    """Distinct individual NPIs with a taxonomy in the condition(s)' core specialty groups other than primary care, per
    100k residents (the geography module's relevant_specialists definition, via the scoring engine's
    provider_counts_for_groups). Returns (rows with lat/lon, groups used)."""
    from measure_it.scoring import opportunity as O
    groups = []
    for cid in members:
        try:
            groups += [g for g in O.condition_core_groups(cid) if g != "primary_care"]
        except KeyError:
            continue
    groups = sorted(set(groups))
    if not groups:
        return pd.DataFrame(columns=["geo_id", "lat", "lon", "n", "per_100k"]), []
    n, rate = O.provider_counts_for_groups(groups, level)
    geo = O.spatial().levels[level].geo
    out = pd.DataFrame({"geo_id": geo["geo_id"].astype(str).to_numpy(), "geo_name": geo["geo_name"].to_numpy(),
                        "lat": geo["lat"].to_numpy(float), "lon": geo["lon"].to_numpy(float),
                        "population": geo["population_total"].to_numpy(float), "n": n, "per_100k": rate})
    return out, groups


@st.cache_data(show_spinner=False)
def phenotype_axes(stamp: float) -> list[tuple[str, str]]:
    """(axis_id, axis_label) of the phenotype axes (selector for discover_candidate_measurements)."""
    a = _read("phenotype_axes", ["axis_id", "axis_label"])
    return list(zip(a["axis_id"], a["axis_label"]))


def hrsa_near(lat: float, lon: float, radius_km: float, stamp: float, county_fips: str | None = None) -> pd.DataFrame:
    """HRSA sites within radius_km of a point (great-circle) or inside a county."""
    h = hrsa_points(stamp)
    if not len(h):
        return h
    r = 6371.0088
    p1, p2 = np.radians(lat), np.radians(h["lat"].to_numpy(float))
    dphi, dlmb = p2 - p1, np.radians(h["lon"].to_numpy(float) - lon)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    dist = 2 * r * np.arcsin(np.sqrt(a))
    out = h.assign(distance_km=np.round(dist, 1))
    keep = out["distance_km"] <= radius_km
    if county_fips:
        keep |= out["county_fips"].astype(str).eq(county_fips)
    return out[keep].sort_values("distance_km").reset_index(drop=True)


@st.cache_data(show_spinner=False)
def geographies(stamp: float) -> pd.DataFrame:
    return _read("geographies", ["geo_id", "geo_level", "name", "state_abbr", "state_name", "lat", "lon"])


@st.cache_resource(show_spinner=False)
def geojson(level: str, tolerance: float | None = None, state: str | None = None) -> dict:
    """Simplified display boundaries (measure_it.api.geo; cached once per process, shared, never mutated)."""
    from measure_it.api.geo import boundaries_geojson
    return boundaries_geojson(level, tolerance, state)


def subset_geojson(fc: dict, ids) -> dict:
    """FeatureCollection restricted to ids (keeps each Plotly trace's embedded geometry small)."""
    want = set(map(str, ids))
    return {"type": "FeatureCollection", "features": [f for f in fc["features"] if f["id"] in want]}


@st.cache_data(show_spinner=False)
def measurement_performance(stamp: float) -> pd.DataFrame:
    """measurement_performance: the scored performance record per condition x measurement (metric link)."""
    return _read("measurement_performance", ["object_id", "condition_id", "measurement_id", "measurement_kind",
                                             "performance_status", "quality_tier", "quality_tier_label",
                                             "selected_record_id", "source_kind", "dataset_id", "label_basis",
                                             "comparator", "comparator_kind", "n_cases", "n_controls", "auroc",
                                             "auroc_ci_low", "auroc_ci_high", "op_sensitivity",
                                             "op_sensitivity_ci_low", "op_sensitivity_ci_high", "op_specificity",
                                             "measurement_evidence", "evidence_conflict", "evidence_conflict_note",
                                             "performance_note", "n_records"])
