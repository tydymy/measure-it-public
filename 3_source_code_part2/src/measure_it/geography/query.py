"""Geographic query functions for the MCP / API layer (SPEC agentic layer).

get_condition_burden(condition, geography_level)  -> burden rows with value, evidence level, source resolution,
                                                     inherited flag, source and CI; UNKNOWN + reason for level D
get_geographic_context(geography_id)              -> one geography (FIPS, 'geo:<fips>' or a name such as
                                                     'San Diego County, California'): burden across target
                                                     conditions, vulnerability, ACS context, provider densities by
                                                     specialty group, HRSA sites, trials, NIH projects
resolve_geography(text), resolve_condition(text)  -> helpers

Everything returned is JSON-serialisable (None / python scalars). Missing data are returned as
config.UNKNOWN with a reason, never guessed. All values describe places (ecological), not people.
Tables are produced by measure_it.geography.features.
"""
from __future__ import annotations

import math
import re
import unicodedata
from functools import lru_cache

import numpy as np
import pandas as pd

from ..config import PROCESSED, UNKNOWN, load_config
from ..store import read_table
from .crosswalk import STATE_ABBR_TO_FIPS

LEVELS = ("national", "state", "county")
COUNTY_SUFFIXES = (" county", " parish", " borough", " census area", " city and borough", " municipality",
                   " planning region", " municipio")
ACTIVE_STATUSES = ("RECRUITING", "NOT_YET_RECRUITING", "ENROLLING_BY_INVITATION", "ACTIVE_NOT_RECRUITING")

CAVEAT_ECOLOGICAL = ("Ecological, place-level data: nothing here describes an individual or where any "
                     "participant lives.")
CAVEAT_INHERITED = ("long_covid county values are the STATE Household Pulse Survey estimate carried as inherited "
                    "context (source resolution 'state'); they are not county prevalence.")
CAVEAT_MODELLED = ("long_covid county values are a MODELLED small-area estimate (BRFSS 2023 multilevel regression and "
                   "post-stratification, results/BRFSS_SAE_RESULTS.md), not observed county prevalence: within-state "
                   "differences come from county composition and covariates only.")
CAVEAT_PROVIDERS = ("Provider counts are self-reported NPPES taxonomies at ZIP-centroid locations; a specialty does "
                    "not mean a provider evaluates or treats the condition.")
CAVEAT_TRIALS = ("A trial or grant is a research-activity signal, not evidence that a measurement or treatment works; "
                 "trial sites are geocoded to city centroids.")


# ============================================================================ table cache
@lru_cache(maxsize=None)
def _cached(name: str, mtime: float, columns: tuple | None) -> pd.DataFrame:
    return read_table(name, columns=list(columns) if columns else None)


def _t(name: str, columns: list[str] | None = None) -> pd.DataFrame:
    p = PROCESSED / f"{name}.parquet"
    return _cached(name, p.stat().st_mtime if p.exists() else 0.0, tuple(columns) if columns else None)


def _py(v):
    """Python scalar for JSON; NaN/NA -> None."""
    if v is None:
        return None
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        f = float(v)
        return None if math.isnan(f) else f
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    return v


def _or_unknown(v):
    v = _py(v)
    return UNKNOWN if v is None else v


def _record(row: pd.Series, cols: list[str]) -> dict:
    return {c: _py(row.get(c)) for c in cols}


# ============================================================================ resolution
def _norm(s: str) -> str:
    """Lower-case, accent-folded ('Doña' == 'dona'), dots removed, whitespace collapsed."""
    s = unicodedata.normalize("NFKD", str(s))
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    return re.sub(r"\s+", " ", s.strip().lower().replace(".", ""))


def _condition_names() -> list[tuple[str, set[str]]]:
    out = []
    for c in load_config("conditions")["conditions"]:
        names = [c["id"], c.get("preferred_name", "")] + list(c.get("aliases") or []) + list(c.get("search_terms") or [])
        names += [c["id"].replace("_", " ")]
        out.append((c["id"], {_norm(n) for n in names if n}))
    return out


def resolve_condition(text: str) -> str | None:
    """Map an id, preferred name, alias or search phrase to a canonical condition id.

    Exact (case/accent-insensitive) match first; otherwise a whole-word prefix match ('Lyme' -> 'lyme disease')
    accepted only when exactly one condition matches. Returns None when nothing (or more than one) matches.
    """
    if text is None:
        return None
    t = _norm(str(text).replace("condition:", ""))
    if not t:
        return None
    names = _condition_names()
    for cid, ns in names:
        if t in ns:
            return cid
    if len(t) >= 3:
        hits = {cid for cid, ns in names if any(n.startswith(t + " ") for n in ns)}
        if len(hits) == 1:
            return hits.pop()
    return None


def _state_lookup() -> dict[str, str]:
    g = _t("geographies")
    st = g[g["geo_level"] == "state"]
    out = {_norm(n): f for n, f in zip(st["name"], st["geo_id"])}
    out.update({_norm(a): STATE_ABBR_TO_FIPS[a] for a in STATE_ABBR_TO_FIPS})
    return out


def resolve_geography(text: str) -> dict:
    """Resolve a FIPS code, 'geo:<fips>', 'zcta:<zcta>' or a place name to one geography.

    Returns {'status': 'ok', 'geo_id', 'geo_level', ...} or {'status': 'unknown'|'ambiguous', 'reason', 'candidates'}.
    """
    if text is None or not str(text).strip():
        return {"status": "unknown", "reason": "empty geography id"}
    raw = str(text).strip()
    g = _t("geographies")
    s = re.sub(r"^geo:", "", raw, flags=re.I).strip()
    if _norm(s) in {"us", "usa", "united states", "united states of america", "national", "nation"}:
        return {"status": "ok", "geo_id": "US", "object_id": "geo:US", "geo_level": "national",
                "name": "United States", "state_fips": None, "state_abbr": None, "state_name": None,
                "lat": None, "lon": None, "vintage": "national"}
    zm = re.match(r"^zcta:?\s*(\d{5})$", s, flags=re.I)
    if zm:
        return _zcta_hit(zm.group(1))
    if re.fullmatch(r"\d{1,2}", s):
        s = s.zfill(2)
    if re.fullmatch(r"\d{2}|\d{5}", s):
        hit = g[g["geo_id"] == s]
        if len(hit):
            return _hit(hit.iloc[0])
        if len(s) == 5:
            z = _zcta_hit(s)
            if z["status"] == "ok":
                z["note"] = f"{s} is not a county FIPS code; interpreted as a ZCTA."
                return z
        return {"status": "unknown", "reason": f"no state/county with FIPS {s!r} in geographies (2024 vintage)"}
    states = _state_lookup()
    parts = [p.strip() for p in s.split(",")]
    if len(parts) == 1 and _norm(parts[0]) in states:
        return _hit(g[g["geo_id"] == states[_norm(parts[0])]].iloc[0])
    place = _norm(parts[0])
    state_fips = states.get(_norm(parts[-1])) if len(parts) > 1 else None
    if len(parts) > 1 and state_fips is None:
        return {"status": "unknown", "reason": f"state {parts[-1]!r} not recognised"}
    cty = g[(g["geo_level"] == "county")]
    if state_fips:
        cty = cty[cty["state_fips"] == state_fips]
    names = cty["name"].map(_norm)
    exact = cty[names == place]
    if exact.empty:
        stripped = names.map(lambda n: next((n[: -len(x)] for x in COUNTY_SUFFIXES if n.endswith(x)), n))
        base = place
        for x in COUNTY_SUFFIXES:
            if base.endswith(x):
                base = base[: -len(x)]
        exact = cty[(stripped == base) | (names == base + " city")]
    exact = exact.sort_values("ct_legacy")   # prefer the 2024 vintage
    if len(exact) > 1 and exact["ct_legacy"].any():
        exact = exact[~exact["ct_legacy"]] if (~exact["ct_legacy"]).any() else exact
    if len(exact) == 1:
        return _hit(exact.iloc[0])
    if len(exact) > 1:
        return {"status": "ambiguous", "reason": f"{raw!r} matches {len(exact)} geographies; pass a FIPS code",
                "candidates": [{"geo_id": r.geo_id, "name": f"{r.name}, {r.state_name}"} for r in exact.itertuples()]}
    return {"status": "unknown", "reason": f"no county or state named {raw!r} in geographies (2024 vintage)"}


def _hit(r: pd.Series) -> dict:
    name = r["name"] if r["geo_level"] == "state" else f"{r['name']}, {r['state_name']}"
    out = {"status": "ok", "geo_id": r["geo_id"], "object_id": f"geo:{r['geo_id']}", "geo_level": r["geo_level"],
           "name": name, "state_fips": _py(r["state_fips"]), "state_abbr": _py(r["state_abbr"]),
           "state_name": _py(r["state_name"]), "lat": _py(r["lat"]), "lon": _py(r["lon"]),
           "vintage": _py(r["vintage"])}
    if bool(r.get("ct_legacy", False)):
        out["note"] = ("Legacy Connecticut county (pre-2022). The canonical 2024 geography is the planning region; "
                       "most sources publish on planning regions and values are not remapped.")
    return out


def _zcta_hit(z: str) -> dict:
    ctx = _t("geo_context")
    r = ctx[(ctx["geo_level"] == "zcta") & (ctx["geo_id"] == z)]
    if r.empty:
        return {"status": "unknown", "reason": f"ZCTA {z} not in the 2024 Gazetteer"}
    r = r.iloc[0]
    return {"status": "ok", "geo_id": z, "object_id": f"geo:zcta:{z}", "geo_level": "zcta", "name": f"ZCTA {z}",
            "state_fips": _py(r["state_fips"]), "state_abbr": _py(r["state_abbr"]), "county_fips": _py(r["county_fips"]),
            "lat": _py(r["lat"]), "lon": _py(r["lon"]), "vintage": _py(r["vintage"])}


# ============================================================================ get_condition_burden
BURDEN_ROW_COLS = ["geo_id", "geo_name", "value", "value_unit", "ci_low", "ci_high", "ci_level", "year", "period",
                   "burden_evidence_level", "source_geographic_resolution", "inherited", "measure_id", "measure_label",
                   "metric_type", "measure_role", "derivation", "proxy_for_condition", "suppressed", "source_name",
                   "source_version", "retrieved_at", "burden_row_id"]


def get_condition_burden(condition: str, geography_level: str = "state", *, geo_id: str | None = None,
                         include_alternates: bool = False) -> dict:
    """Burden estimates for a condition at a geography level ('national', 'state' or 'county').

    Returns the pre-specified primary measure (docs/BURDEN_DEFINITIONS.md) for every geography at that level,
    each with value, evidence level (A/B/C/D), source resolution, inherited flag, source and CI.
    Level-D conditions return value UNKNOWN and the reason. `include_alternates` adds the county C-proxy column
    and the pre-specified sensitivity alternates, each labelled with its measure_role.
    """
    cid = resolve_condition(condition)
    if cid is None:
        return {"status": "unknown", "condition_input": condition, "value": UNKNOWN,
                "reason": f"condition {condition!r} not in configs/conditions.yaml (ids, names, aliases)"}
    lvl = str(geography_level).strip().lower()
    if lvl not in LEVELS:
        return {"status": "unknown", "condition_id": cid, "value": UNKNOWN,
                "reason": f"geography_level must be one of {LEVELS}; got {geography_level!r}"}
    defs = _t("geo_burden_definitions")
    d = defs[(defs["condition_id"] == cid) & (defs["geo_level"] == lvl)].iloc[0]
    out = {"status": "ok", "condition_id": cid, "condition_input": condition, "geography_level": lvl,
           "burden_evidence_level": d["burden_evidence_level"],
           "primary_measure_id": _or_unknown(d["primary_measure_id"]),
           "primary_source_resolution": _or_unknown(d["primary_source_resolution"]),
           "inherited": bool(d["inherited"]), "rationale": d["rationale"],
           "uncertainty_multiplier": _py(d["uncertainty_multiplier"]),
           "alternates": _py(d["alternates"]) or "", "caveats": [CAVEAT_ECOLOGICAL]}
    if d["burden_evidence_level"] == "D":
        out.update({"status": "unknown", "value": UNKNOWN, "reason": d["rationale"], "rows": [], "n_rows": 0})
        return out
    b = _t("geo_condition_burden")
    roles = ["primary", "county_proxy", "alternate"] if include_alternates else ["primary"]
    rows = b[(b["condition_id"] == cid) & (b["geo_level"] == lvl) & b["measure_role"].isin(roles)]
    if geo_id is not None:
        gr = resolve_geography(geo_id)
        if gr.get("status") != "ok":
            out.update({"status": "unknown", "value": UNKNOWN, "reason": gr.get("reason"), "rows": [], "n_rows": 0,
                        **({"candidates": gr["candidates"]} if gr.get("candidates") else {})})
            return out
        if gr["geo_level"] != lvl:
            out.update({"status": "unknown", "value": UNKNOWN, "rows": [], "n_rows": 0,
                        "reason": f"{geo_id!r} resolves to a {gr['geo_level']} ({gr['geo_id']}), but "
                                  f"geography_level={lvl!r}; pass geography_level={gr['geo_level']!r}"})
            return out
        rows = rows[rows["geo_id"] == gr["geo_id"]]
        if rows.empty:
            reason = "no value for this condition/geography in the primary measure (see burden_unknown_reason in geo_condition_features)"
            if lvl in ("state", "county"):
                fz = _t("geo_condition_features", ["geo_id", "geo_level", "condition_id", "burden_unknown_reason"])
                hit = fz[(fz["geo_id"] == gr["geo_id"]) & (fz["geo_level"] == lvl) & (fz["condition_id"] == cid)]
                if len(hit) and isinstance(hit["burden_unknown_reason"].iloc[0], str):
                    reason = hit["burden_unknown_reason"].iloc[0]
            out.update({"status": "unknown", "value": UNKNOWN, "rows": [], "n_rows": 0, "geo_id": gr["geo_id"],
                        "reason": reason})
            return out
    rows = rows.sort_values(["measure_role", "geo_id"])
    recs = []
    for _, r in rows.iterrows():
        rec = _record(r, BURDEN_ROW_COLS)
        rec["object_id"] = f"geo:{r['geo_id']}" if r["geo_level"] != "national" else "geo:US"
        if rec["value"] is None:
            rec["value"] = UNKNOWN
            rec["value_missing_reason"] = ("suppressed by the source" if rec.get("suppressed") else
                                           "not published for this geography")
        for c in ("ci_low", "ci_high"):
            if rec[c] is None:
                rec[c] = UNKNOWN
        recs.append(rec)
    if not recs:
        out.update({"status": "unknown", "value": UNKNOWN, "rows": [], "n_rows": 0,
                    "reason": "no value for this condition/geography in the primary measure (see burden_unknown_reason in geo_condition_features)"})
        return out
    if cid == "long_covid" and lvl == "county":
        derivs = {str(r.get("derivation")) for r in recs}
        if "inherited_from_state" in derivs:
            out["caveats"].append(CAVEAT_INHERITED)
        if "mrp_small_area_estimate" in derivs:
            out["caveats"].append(CAVEAT_MODELLED)
    if d["burden_evidence_level"] == "C":
        out["caveats"].append("Level C: symptom/comorbidity/antecedent proxy, not a measure of the condition itself.")
    if lvl != "national" and cid in ("me_cfs", "fibromyalgia", "migraine"):
        out["caveats"].append("Medicare fee-for-service beneficiaries only (CMS MMD coded-condition prevalence).")
    out["rows"] = recs
    out["n_rows"] = len(recs)
    return out


# ============================================================================ get_geographic_context
ACS_KEYS = ["total_population", "n_adults_18plus", "median_age", "pct_age_65_plus", "pct_with_disability",
            "pct_uninsured", "pct_below_poverty", "median_household_income", "pct_households_no_vehicle",
            "pct_households_broadband", "pop_density_per_km2", "acs_vintage"]
SVI_KEYS = ["svi_overall", "svi_socioeconomic", "svi_household", "svi_minority", "svi_housing_transport",
            "svi_ranking_universe", "svi_derivation", "svi_release"]
PLACES_KEYS = ["places_PHLTH_crude", "places_GHLTH_crude", "places_DISABILITY_crude", "places_COGNITION_crude",
               "places_MOBILITY_crude", "places_DEPRESSION_crude", "places_SLEEP_crude", "places_LPA_crude",
               "places_ACCESS2_crude", "places_CHECKUP_crude", "places_LACKTRPT_crude"]
FEATURE_KEYS = ["burden_value", "burden_measure_id", "burden_evidence_level", "burden_source_resolution",
                "burden_inherited", "burden_ci_low", "burden_ci_high", "burden_year", "burden_period",
                "burden_source_name", "burden_unknown_reason", "burden_count", "burden_count_kind",
                "burden_proxy_measure_id", "burden_proxy_value", "burden_proxy_level", "burden_uncertainty_multiplier",
                "relevant_specialty_groups", "relevant_providers_n", "relevant_providers_per_100k",
                "relevant_specialists_n", "relevant_specialists_per_100k", "trials_in_geo_n", "trials_in_geo_active_n",
                "trials_in_geo_recruiting_n", "trials_in_geo_completed_n", "trials_in_geo_or_within_50km_n",
                "trials_within_50km_n", "trials_within_100km_n", "nih_projects_n", "nih_core_projects_n",
                "nih_active_core_projects_n", "nih_orgs_n", "diagnostic_desert", "desert_pct_burden",
                "desert_pct_vulnerability", "desert_pct_low_providers", "desert_pct_low_trials",
                "diagnostic_desert_access_only", "desert_unknown_reason", "desert_burden_resolution"]


def _providers_by_group(geo: dict, population) -> list[dict]:
    pdc = _t("provider_density_county", ["county_fips", "state_abbr", "specialty_group", "group_label", "group_type",
                                         "n_individual_providers", "n_individual_physicians", "n_organizations"])
    if geo["geo_level"] == "county":
        d = pdc[pdc["county_fips"] == geo["geo_id"]]
        note = "county of the provider's primary practice ZIP (ZCTA internal point)"
    else:
        d = pdc[pdc["county_fips"].str[:2] == geo["geo_id"]]
        d = d.groupby(["specialty_group", "group_label", "group_type"], as_index=False)[
            ["n_individual_providers", "n_individual_physicians", "n_organizations"]].sum()
        note = "sum over the state's counties (geocoded NPIs only)"
    out = []
    for _, r in d.sort_values(["group_type", "specialty_group"]).iterrows():
        per = (r["n_individual_providers"] / population * 1e5) if population else None
        out.append({"specialty_group": r["specialty_group"], "label": r["group_label"], "group_type": r["group_type"],
                    "n_individual_providers": int(r["n_individual_providers"]),
                    "n_individual_physicians": int(r["n_individual_physicians"]),
                    "n_organizations": int(r["n_organizations"]),
                    "individual_providers_per_100k": _or_unknown(per), "placement": note})
    return out


def _trials_overall(geo: dict) -> dict:
    tc = _t("trial_conditions", ["nct_id", "condition_id", "condition_literal_match"])
    tc = tc[tc["condition_literal_match"].astype(bool)]
    ts = _t("trial_sites", ["nct_id", "country_group", "county_fips", "state_fips"])
    ct = _t("clinical_trials", ["nct_id", "overall_status"]).set_index("nct_id")["overall_status"]
    col = "county_fips" if geo["geo_level"] == "county" else "state_fips"
    ids = set(ts[(ts["country_group"] == "US") & (ts[col] == geo["geo_id"])]["nct_id"]) & set(tc["nct_id"])
    st = ct.reindex(sorted(ids))
    by = tc[tc["nct_id"].isin(ids)].groupby("condition_id")["nct_id"].nunique().to_dict()
    return {"distinct_relevant_trials_with_site_here": len(ids),
            "active": int(st.isin(ACTIVE_STATUSES).sum()), "recruiting": int((st == "RECRUITING").sum()),
            "completed": int((st == "COMPLETED").sum()), "by_condition": {k: int(v) for k, v in by.items()},
            "precision_view": "literal condition match (trial_conditions.condition_literal_match)",
            "example_nct_ids": [f"trial:{n}" for n in sorted(ids)[:10]]}


def _nih_overall(geo: dict) -> dict:
    links = _t("nih_project_conditions", ["appl_id", "condition_id", "match_tier", "likely_false_positive"])
    links = links[(links["match_tier"] == "title_abstract") & ~links["likely_false_positive"].fillna(False).astype(bool)]
    pj = _t("nih_projects", ["appl_id", "core_project_num", "org_name", "org_ipf_code", "county_fips", "state_fips",
                             "is_active"])
    col = "county_fips" if geo["geo_level"] == "county" else "state_fips"
    p = pj[(pj[col] == geo["geo_id"]) & pj["appl_id"].isin(links["appl_id"])]
    orgs = p["org_ipf_code"].fillna("name:" + p["org_name"].fillna("")).nunique()
    active = p[p["is_active"].fillna(False).astype(bool)]["core_project_num"].nunique()
    top_orgs = p["org_name"].value_counts().head(5)
    return {"distinct_projects_appl_id": int(p["appl_id"].nunique()), "core_projects": int(p["core_project_num"].nunique()),
            "active_core_projects": int(active), "organizations": int(orgs),
            "top_organizations": [{"org_name": k, "n_award_records": int(v)} for k, v in top_orgs.items()],
            "precision_view": "title/abstract match; likely false positives excluded; FY2015-FY2026"}


def _hrsa(geo: dict, population) -> dict:
    h = _t("facilities__hrsa", ["site_bphc_number", "facility_name", "city", "is_active", "is_service_delivery_site",
                                "county_fips", "state_fips", "site_type", "location_setting"])
    col = "county_fips" if geo["geo_level"] == "county" else "state_fips"
    h = h[(h[col] == geo["geo_id"]) & h["is_active"].astype(bool) & h["is_service_delivery_site"].astype(bool)]
    n = int(h["site_bphc_number"].nunique())
    out = {"n_service_delivery_sites": n,
           "sites_per_100k": _or_unknown(n / population * 1e5 if population else None),
           "note": "HRSA health-center program sites (FQHC + look-alike); no service-line data in the source"}
    if geo["geo_level"] == "county":
        out["sites"] = [{"object_id": f"hrsa_site:{r.site_bphc_number}", "facility_name": r.facility_name,
                         "city": r.city, "location_setting": r.location_setting}
                        for r in h.sort_values("facility_name").head(15).itertuples()]
    return out


def get_geographic_context(geography_id: str) -> dict:
    """Context for one state or county (FIPS, 'geo:<fips>' or a name like 'San Diego County, California')."""
    geo = resolve_geography(geography_id)
    if geo.get("status") != "ok":
        return {"status": geo.get("status", "unknown"), "geography_input": geography_id, "value": UNKNOWN, **geo}
    ctx = _t("geo_context")
    lvl = geo["geo_level"]
    if lvl == "national":
        return {"status": "unknown", "geography": geo, "value": UNKNOWN,
                "reason": "geographic context is built for states, counties and ZCTAs only; for national burden use "
                          "get_condition_burden(condition, 'national')"}
    row = ctx[(ctx["geo_level"] == lvl) & (ctx["geo_id"] == geo["geo_id"])]
    if lvl == "zcta":
        r = row.iloc[0] if len(row) else pd.Series(dtype=object)
        return {"status": "ok", "geography": geo,
                "acs_context": {k: _or_unknown(r.get(k)) for k in ACS_KEYS if k in r.index},
                "places_context": {k: _or_unknown(r.get(k)) for k in PLACES_KEYS if k in r.index},
                "note": "ZCTA rows carry ACS and PLACES context only; burden, vulnerability and facility features "
                        "are county-level (see the containing county).",
                "caveats": [CAVEAT_ECOLOGICAL]}
    if row.empty:
        return {"status": "unknown", "geography": geo, "value": UNKNOWN,
                "reason": "geography has no geo_context row (legacy Connecticut county or not built)"}
    r = row.iloc[0]
    pop = _py(r.get("total_population"))
    feats = _t("geo_condition_features")
    f = feats[(feats["geo_level"] == lvl) & (feats["geo_id"] == geo["geo_id"])]
    burden = []
    for _, fr in f.sort_values("condition_id").iterrows():
        rec = {"condition_id": fr["condition_id"]}
        for k in FEATURE_KEYS:
            v = _py(fr.get(k))
            rec[k] = UNKNOWN if v is None and not k.endswith("reason") else v
        if rec["burden_evidence_level"] == "D":
            rec["burden_value"] = UNKNOWN
        if "small-area" in str(fr.get("burden_source_name") or "") and lvl == "county":
            rec["burden_value_note"] = ("MODELLED small-area estimate (BRFSS 2023 MRP), not county prevalence: "
                                        "within-state differences come from county composition and covariates only.")
        if rec.get("burden_inherited") is True:
            rec["burden_value_note"] = ("STATE estimate carried as inherited context (source resolution 'state'); "
                                        "not county prevalence. The county diagnostic_desert uses it as its burden "
                                        "component, so that index ranks counties partly by their state's value.")
        burden.append(rec)
    min_pop = load_config("scoring")["ranking"]["min_population"]
    rucc = {k: _or_unknown(r.get(k)) for k in ("rucc_2023", "rucc_description", "metro_status", "is_metro")} \
        if lvl == "county" else {"note": UNKNOWN + " (RUCC is a county classification)"}
    hpsa = {k.replace("hpsa_", ""): _or_unknown(r.get(k)) for k in r.index if k.startswith("hpsa_")} \
        if lvl == "county" else {"note": UNKNOWN + " (HPSA/MUA summarised by county only)"}
    caveats = [CAVEAT_ECOLOGICAL, CAVEAT_MODELLED, CAVEAT_INHERITED, CAVEAT_PROVIDERS, CAVEAT_TRIALS]
    if pop is not None and pop < min_pop:
        caveats.append(f"Small population ({int(pop):,} < {min_pop:,}): rates are unstable.")
    if geo.get("state_fips") in ("60", "66", "69", "72", "78"):
        caveats.append("Territory: outside the 50 states + DC analysis universe; SVI/PLACES/HPS are missing or "
                       "not poolable, so diagnostic_desert is UNKNOWN.")
    return {
        "status": "ok", "geography": geo,
        "population": {"total_population": _or_unknown(pop), "adults_18plus": _or_unknown(r.get("n_adults_18plus")),
                       "source": "ACS 2020-2024 5-year", "small_population_flag": (pop is not None and pop < min_pop)},
        "burden_by_condition": burden,
        "vulnerability": {k: _or_unknown(r.get(k)) for k in SVI_KEYS},
        "acs_context": {k: _or_unknown(r.get(k)) for k in ACS_KEYS},
        "places_context": {k: _or_unknown(r.get(k)) for k in PLACES_KEYS},
        "places_derivation": _or_unknown(r.get("places_derivation")),
        "rurality": rucc, "shortage_designations": hpsa,
        "providers_by_specialty_group": _providers_by_group(geo, pop),
        "hrsa_health_centers": _hrsa(geo, pop),
        "trials": _trials_overall(geo),
        "nih_projects": _nih_overall(geo),
        "caveats": caveats,
        "provenance": ["geo_condition_features", "geo_context", "geo_condition_burden", "provider_density_county",
                       "facilities__hrsa", "trial_sites", "trial_conditions", "clinical_trials", "nih_projects",
                       "nih_project_conditions"],
    }
