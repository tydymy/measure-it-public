"""Clinic matching (SPEC Phase 5): candidate implementation / study partners around a geography.

    find_candidate_clinics(geography_id, condition, measurement, radius_km=None, max_sites=None) -> dict

Candidates are resolved facilities (`facilities`, built by measure_it.facilities.registry) within `radius_km`
of the geography's Census internal point, plus any inside a county geography. Six characteristics are ranked
SEPARATELY and never combined into one score:

    clinical_specialty_match       condition specialties + measurement implementer groups (configs/relevance.yaml,
                                   curated assumptions) present in the facility's self-reported NPPES taxonomies
    relevant_trial_history         U.S. ClinicalTrials.gov trials of the condition at the facility (literal view)
    NIH_research_activity          NIH RePORTER core projects of the facility's organisation for the condition
    technology_experience          trials at the facility whose registered text mentions the measurement
    community_access               FQHC/look-alike, community-type facility, county primary-care HPSA, nonmetro county
    distance_to_target_population  km to the geography's internal point

The returned list is a round-robin over the six rankings (1st of each, then 2nd of each, ...). Every candidate
carries its raw values, its rank on each characteristic, plain-language reasons and the framing sentence below.
It never says "best". A taxonomy code does not mean a facility treats a condition; a registered trial is not
evidence that a measurement works; nothing here locates any patient.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from ..config import UNKNOWN, load_config
from ..ontology.normalize import norm_text, normalize_condition
from ..store import processed_path, read_table
from .registry import SEP, haversine_km

FRAMING = "This facility has characteristics suggesting it may be a viable implementation or study partner."
CHARACTERISTICS = ["clinical_specialty_match", "relevant_trial_history", "NIH_research_activity",
                   "technology_experience", "community_access", "distance_to_target_population"]
GUARDRAILS = [
    "Candidates are facilities with characteristics that may make them viable implementation or study partners; "
    "the list is not a ranking of quality and never identifies a 'best' clinic.",
    "Specialty groups are self-reported NPPES taxonomies; they do not mean a facility evaluates or treats the condition. "
    "Organisation NPIs under physician taxonomies are group practices.",
    "Condition specialties and measurement implementer groups are curated assumptions (configs/relevance.yaml).",
    "A trial that registered a measurement is a research-activity signal, not evidence that the measurement works.",
    "NIH awards are a research-capability signal, not patient burden; award obligations are not expenditure.",
    "Facility locations are ZIP centroids (NPPES, RePORTER), HRSA address points, or ClinicalTrials.gov city/ZIP "
    "centroids; distances are approximate. No patient location is inferred.",
    "County HPSA and RUCC flags describe the facility's county, not the facility.",
]
CHARACTERISTIC_DEFINITIONS = {
    "clinical_specialty_match": "Condition specialty groups and measurement implementer groups (relevance.yaml) present "
                                "in the facility's NPPES taxonomy groups (HRSA roster adds fqhc). Ranked: matches both "
                                "kinds first, then number of distinct matched groups.",
    "relevant_trial_history": "Distinct U.S. ClinicalTrials.gov trials of the condition registered at the facility "
                              "(condition named literally in the record). Ranked by count, then open-status count, "
                              "then latest start year.",
    "NIH_research_activity": "Distinct NIH RePORTER core projects of the facility's awardee organisation whose title/"
                             "abstract mention the condition (likely false positives excluded). Ranked by core "
                             "projects, then active core projects, then award obligations.",
    "technology_experience": "Distinct trials at the facility (any target condition) whose registered titles/summary/"
                             "keywords/interventions/outcomes match the measurement's configs/measurements.yaml "
                             "patterns. Ranked by count, then count for this condition.",
    "community_access": "Number of flags: FQHC/look-alike (HRSA roster or NPPES taxonomy); community/rural-health/"
                        "critical-access/public-health taxonomy; county primary-care HPSA area/population "
                        "designation; nonmetro county (RUCC 4-9).",
    "distance_to_target_population": "Great-circle km from the facility location to the geography's Census internal "
                                     "point (ascending). Not ranked for state geographies.",
}


def _cfg_matching() -> dict:
    return load_config("scoring").get("clinic_matching", {})


# --------------------------------------------------------------------------------------------------------------
# resolvers
# --------------------------------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _geographies() -> pd.DataFrame:
    return read_table("geographies", columns=["geo_id", "geo_level", "name", "state_fips", "state_abbr", "state_name",
                                              "county_fips", "lat", "lon", "ct_legacy"])


def resolve_geography(geography_id) -> dict:
    """County FIPS ('06073', 'geo:06073'), state FIPS/abbreviation/name, or 'San Diego County, California'."""
    g = _geographies()
    q = str(geography_id or "").strip()
    if not q:
        return {"status": UNKNOWN, "query": q, "reason": "empty geography"}
    s = q[4:] if q.lower().startswith("geo:") else q
    note = ""
    if re.fullmatch(r"\d{4}", s):
        s, note = "0" + s, "4-digit FIPS zero-padded to 5 digits"
    if re.fullmatch(r"\d", s):
        s, note = "0" + s, "1-digit state FIPS zero-padded"
    row = g[g["geo_id"] == s]
    if row.empty and re.fullmatch(r"[A-Za-z]{2}", s):
        row = g[(g["geo_level"] == "state") & (g["state_abbr"] == s.upper())]
    if row.empty:
        key = norm_text(s)
        st = g[g["geo_level"] == "state"]
        row = st[st["state_name"].map(norm_text) == key]
        if row.empty and "," in s:
            county, state = (x.strip() for x in s.rsplit(",", 1))
            srow = st[(st["state_name"].map(norm_text) == norm_text(state)) | (st["state_abbr"] == state.upper())]
            if not srow.empty:
                cand = g[(g["geo_level"] == "county") & (g["state_fips"] == srow["state_fips"].iloc[0]) & ~g["ct_legacy"]]
                nm = cand["name"].map(norm_text)
                ck = norm_text(county)
                row = cand[nm == ck]
                if row.empty:
                    row = cand[nm.str.replace(r" (county|parish|borough|census area|municipality|city and borough)$",
                                              "", regex=True) == re.sub(
                        r" (county|parish|borough|census area|municipality|city and borough)$", "", ck)]
    if row.empty:
        return {"status": UNKNOWN, "query": q, "reason": "not a state or 2024 county geo_id, abbreviation or name "
                                                          "in the geographies table"}
    if len(row) > 1:
        return {"status": UNKNOWN, "query": q, "reason": "ambiguous geography",
                "candidates": row["geo_id"].tolist()}
    r = row.iloc[0]
    out = {"status": "matched", "query": q, "geo_id": r["geo_id"], "object_id": f"geo:{r['geo_id']}",
           "geo_level": r["geo_level"], "name": r["name"], "state_abbr": r["state_abbr"],
           "state_name": r["state_name"], "lat": float(r["lat"]), "lon": float(r["lon"]),
           "ct_legacy": bool(r["ct_legacy"])}
    if note:
        out["note"] = note
    if out["ct_legacy"]:
        out["note"] = ("legacy Connecticut county (pre-2022); facilities carry 2024 planning-region counties, so only "
                       "the radius search applies")
    return out


def resolve_condition(condition) -> dict:
    ids = {c["id"]: c for c in load_config("conditions")["conditions"]}
    q = str(condition or "").strip()
    if q in ids:
        return {"status": "matched", "query": q, "condition_id": q, "preferred_name": ids[q]["preferred_name"],
                "object_id": f"condition:{q}", "match_reason": "canonical condition id"}
    try:
        r = normalize_condition(q)
    except FileNotFoundError as e:  # ontology tables not built
        return {"status": UNKNOWN, "query": q, "reason": f"ontology tables unavailable: {e}"}
    if r.get("status") == "matched" and r["matches"]:
        m = r["matches"][0]
        cid = m["canonical_condition_id"]
        return {"status": "matched", "query": q, "condition_id": cid,
                "preferred_name": ids.get(cid, {}).get("preferred_name", m.get("preferred_name")),
                "object_id": f"condition:{cid}", "match_reason": m.get("match_reason")}
    reason = r.get("reason") or f"normalize_condition status {r.get('status')}: not an unambiguous match"
    cands = [m["canonical_condition_id"] for m in r.get("matches", []) + r.get("other_candidates", [])]
    return {"status": UNKNOWN, "query": q, "reason": reason, "candidates": cands}


def resolve_measurement(measurement) -> dict:
    """A measurement class id/name/alias or bundle id/label/alias -> member classes and implementer groups.

    The shared resolver (measure_it.measurements.resolve, single=True: one class or bundle is required here)."""
    from ..measurements.resolve import resolve_measurement as _resolve
    return _resolve(measurement, single=True)


def condition_groups(condition_id: str) -> list[str]:
    return list(load_config("relevance").get("condition_specialties", {}).get(condition_id, {}).get("core", []))


# --------------------------------------------------------------------------------------------------------------
# base data
# --------------------------------------------------------------------------------------------------------------

@dataclass
class Base:
    fac: pd.DataFrame            # clinic candidates (geocoded), positional index 0..n-1
    groups: list                 # group ids (columns of G)
    G: np.ndarray                # n x len(groups) bool membership
    ft: pd.DataFrame             # facility_trials with fidx
    fn: pd.DataFrame             # facility_nih_projects with fidx
    lat: np.ndarray
    lon: np.ndarray
    county: np.ndarray
    state: np.ndarray
    hpsa: np.ndarray             # county primary-care HPSA area designation (location attribute)
    nonmetro: np.ndarray         # county nonmetro (location attribute)
    fqhc: np.ndarray             # facility attribute
    community_type: np.ndarray   # facility attribute


def _stamp() -> tuple:
    return tuple(processed_path(n).stat().st_mtime if processed_path(n).exists() else 0
                 for n in ("facilities", "facility_trials", "facility_nih_projects"))


@lru_cache(maxsize=2)
def _base(stamp: tuple) -> Base:
    cols = ["facility_id", "object_id", "facility_name", "primary_kind", "facility_kinds", "sources", "address", "city",
            "state", "zip5", "county_fips", "state_fips", "lat", "lon", "geocode_precision", "specialty_groups",
            "is_fqhc_or_lookalike", "fqhc_evidence", "is_community_type_facility", "is_hrsa_site",
            "hrsa_is_look_alike", "hrsa_location_type", "county_pc_hpsa_area", "county_pc_hpsa_coverage", "county_nonmetro", "rucc_2023",
            "npi_object_ids", "hrsa_site_object_ids", "resolution_tier", "match_confidence_min", "is_clinic_candidate",
            "source_version", "retrieved_at"]
    f = read_table("facilities", columns=cols)
    f = f[f["is_clinic_candidate"].astype(bool) & f["lat"].notna()].reset_index(drop=True)
    sg = f["specialty_groups"].fillna("").astype(str)
    groups = sorted({g for s in sg.unique() for g in s.split(SEP) if g} | {"fqhc"})
    G = np.zeros((len(f), len(groups)), dtype=bool)
    gi = {g: i for i, g in enumerate(groups)}
    for r, s in enumerate(sg):
        for g in s.split(SEP):
            if g:
                G[r, gi[g]] = True
    fq = f["is_fqhc_or_lookalike"].astype(bool).to_numpy()
    G[:, gi["fqhc"]] |= fq
    idx = pd.Series(np.arange(len(f)), index=f["facility_id"])
    ft = read_table("facility_trials", columns=["facility_id", "nct_id", "condition_ids_literal", "condition_ids_all",
                                                "measurement_classes", "site_recruiting", "overall_recruiting",
                                                "is_open_status", "start_year", "is_device_or_diagnostic",
                                                "overall_status", "brief_title"])
    ft["fidx"] = ft["facility_id"].map(idx)
    ft = ft[ft["fidx"].notna()].copy()
    ft["fidx"] = ft["fidx"].astype(int)
    fn = read_table("facility_nih_projects", columns=["facility_id", "appl_id", "core_project_num", "project_num",
                                                      "fiscal_year", "is_subproject", "is_active", "award_amount",
                                                      "condition_ids_precision", "condition_ids_all",
                                                      "project_title"])
    fn["fidx"] = fn["facility_id"].map(idx)
    fn = fn[fn["fidx"].notna()].copy()
    fn["fidx"] = fn["fidx"].astype(int)
    return Base(fac=f, groups=groups, G=G, ft=ft, fn=fn,
                lat=f["lat"].to_numpy(float), lon=f["lon"].to_numpy(float),
                county=f["county_fips"].astype(object).where(f["county_fips"].notna(), "").to_numpy(),
                state=f["state"].astype(object).where(f["state"].notna(), "").to_numpy(),
                hpsa=f["county_pc_hpsa_area"].astype(bool).to_numpy(),
                nonmetro=f["county_nonmetro"].astype(bool).to_numpy(), fqhc=fq,
                community_type=f["is_community_type_facility"].astype(bool).to_numpy())


def load_base() -> Base:
    return _base(_stamp())


def _has(s: pd.Series, ids) -> pd.Series:
    ids = [i for i in ids if i]
    if not ids:
        return pd.Series(False, index=s.index)
    return s.fillna("").str.contains(r"(?:^|\|)(?:" + "|".join(map(re.escape, ids)) + r")(?:\||$)", regex=True)


@dataclass
class Features:
    n_cond_groups: np.ndarray
    n_meas_groups: np.ndarray
    n_union_groups: np.ndarray
    trials_cond: np.ndarray
    trials_cond_open: np.ndarray
    trials_cond_latest: np.ndarray
    nih_core: np.ndarray
    nih_active_core: np.ndarray
    nih_obligations: np.ndarray
    tech: np.ndarray
    tech_cond: np.ndarray
    intrinsic_access: np.ndarray   # fqhc + community-type facility (facility attributes)
    eligible: np.ndarray
    cond_groups: list
    meas_groups: list


def _count_by(fidx: pd.Series, n: int) -> np.ndarray:
    out = np.zeros(n, dtype=np.int64)
    if len(fidx):
        vc = fidx.value_counts()
        out[vc.index.to_numpy(int)] = vc.to_numpy()
    return out


@lru_cache(maxsize=16)
def _features_cached(stamp: tuple, condition_id: str, members: tuple, meas_groups: tuple) -> Features:
    from ..ingestion.nih_reporter import obligation_total

    b = _base(stamp)
    n = len(b.fac)
    gi = {g: i for i, g in enumerate(b.groups)}
    cgroups = [g for g in condition_groups(condition_id) if g in gi]
    mgroups = [g for g in meas_groups if g in gi]
    Gc = b.G[:, [gi[g] for g in cgroups]] if cgroups else np.zeros((n, 0), bool)
    Gm = b.G[:, [gi[g] for g in mgroups]] if mgroups else np.zeros((n, 0), bool)
    union = sorted(set(cgroups) | set(mgroups))
    Gu = b.G[:, [gi[g] for g in union]] if union else np.zeros((n, 0), bool)

    ft = b.ft
    tc = ft[_has(ft["condition_ids_literal"], [condition_id])].drop_duplicates(["fidx", "nct_id"])
    trials_cond = _count_by(tc["fidx"], n)
    trials_open = _count_by(tc.loc[tc["is_open_status"].astype(bool), "fidx"], n)
    latest = np.zeros(n, dtype=np.int64)
    if len(tc):
        ly = tc.groupby("fidx")["start_year"].max().dropna()
        latest[ly.index.to_numpy(int)] = ly.to_numpy(int)
    tm = ft[_has(ft["measurement_classes"], list(members))].drop_duplicates(["fidx", "nct_id"])
    tech = _count_by(tm["fidx"], n)
    tech_cond = _count_by(tm.loc[_has(tm["condition_ids_literal"], [condition_id]), "fidx"], n)

    fn = b.fn[_has(b.fn["condition_ids_precision"], [condition_id])]
    core = fn.drop_duplicates(["fidx", "core_project_num"])
    nih_core = _count_by(core["fidx"], n)
    act = fn[fn["is_active"].fillna(False).astype(bool)].drop_duplicates(["fidx", "core_project_num"])
    nih_active = _count_by(act["fidx"], n)
    obl = np.full(n, np.nan)
    for i, g in fn.groupby("fidx"):
        if g["award_amount"].notna().any():
            obl[int(i)] = obligation_total(g)
    intrinsic = b.fqhc.astype(int) + b.community_type.astype(int)
    n_union = Gu.sum(1)
    eligible = (n_union > 0) | (trials_cond > 0) | (nih_core > 0) | (tech > 0) | (intrinsic > 0)
    return Features(n_cond_groups=Gc.sum(1), n_meas_groups=Gm.sum(1), n_union_groups=n_union,
                    trials_cond=trials_cond, trials_cond_open=trials_open, trials_cond_latest=latest,
                    nih_core=nih_core, nih_active_core=nih_active, nih_obligations=obl, tech=tech,
                    tech_cond=tech_cond, intrinsic_access=intrinsic, eligible=eligible,
                    cond_groups=cgroups, meas_groups=mgroups)


def features(condition_id: str, measurement: dict | None) -> Features:
    members = tuple(measurement["members"]) if measurement else ()
    mg = tuple(measurement["implementer_groups"]) if measurement else ()
    return _features_cached(_stamp(), condition_id, members, mg)


# --------------------------------------------------------------------------------------------------------------
# core selection (shared by find_candidate_clinics and the Test 6 shuffle)
# --------------------------------------------------------------------------------------------------------------

def characteristic_orders(F: Features, pool: np.ndarray, dist: np.ndarray, access: np.ndarray,
                          rank_distance: bool = True) -> dict[str, np.ndarray]:
    """For each characteristic, eligible pool facilities that have it, in rank order (ties: distance, index)."""
    e = pool[F.eligible[pool]]
    d = dist[e]
    orders = {}

    def order(mask, *keys_primary_first):
        sel = e[mask]
        if len(sel) == 0:
            return sel
        keys = [sel, dist[sel]] + [k[mask] for k in reversed(keys_primary_first)]
        return sel[np.lexsort(keys)]

    both = ((F.n_cond_groups[e] > 0) & (F.n_meas_groups[e] > 0)).astype(int)
    orders["clinical_specialty_match"] = order(F.n_union_groups[e] > 0, -both, -F.n_union_groups[e])
    orders["relevant_trial_history"] = order(F.trials_cond[e] > 0, -F.trials_cond[e], -F.trials_cond_open[e],
                                             -F.trials_cond_latest[e])
    orders["NIH_research_activity"] = order(F.nih_core[e] > 0, -F.nih_core[e], -F.nih_active_core[e],
                                            -np.nan_to_num(F.nih_obligations[e]))
    orders["technology_experience"] = order(F.tech[e] > 0, -F.tech[e], -F.tech_cond[e])
    orders["community_access"] = order(access[e] > 0, -access[e])
    orders["distance_to_target_population"] = (e[np.lexsort((e, d))] if rank_distance else e[:0])
    return orders


def round_robin(orders: dict[str, np.ndarray], max_sites: int) -> tuple[list, dict]:
    selected, via = [], {}
    pos = {c: 0 for c in CHARACTERISTICS}
    while len(selected) < max_sites:
        progressed = False
        for c in CHARACTERISTICS:
            o = orders.get(c, np.array([], dtype=int))
            while pos[c] < len(o) and int(o[pos[c]]) in via:
                pos[c] += 1
            if pos[c] < len(o):
                f = int(o[pos[c]])
                via[f] = c
                selected.append(f)
                pos[c] += 1
                progressed = True
                if len(selected) >= max_sites:
                    break
        if not progressed:
            break
    return selected, via


def select(F: Features, lat: np.ndarray, lon: np.ndarray, county: np.ndarray, state: np.ndarray,
           hpsa: np.ndarray, nonmetro: np.ndarray, geo: dict, radius_km: float, max_sites: int) -> dict:
    """Pool -> orders -> round-robin selection. Location arrays may be permuted (Test 6)."""
    dist = haversine_km(geo["lat"], geo["lon"], lat, lon)
    if geo["geo_level"] == "state":
        pool_mask = state == geo["state_abbr"]
        rank_distance = False
    else:
        pool_mask = dist <= radius_km
        if not geo.get("ct_legacy"):
            pool_mask |= county == geo["geo_id"]
        rank_distance = True
    pool = np.flatnonzero(pool_mask)
    access = F.intrinsic_access + hpsa.astype(int) + nonmetro.astype(int)
    orders = characteristic_orders(F, pool, dist, access, rank_distance)
    sel, via = round_robin(orders, max_sites)
    return {"pool": pool, "dist": dist, "access": access, "orders": orders, "selected": sel, "via": via,
            "n_eligible": int(F.eligible[pool].sum())}


# --------------------------------------------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------------------------------------------

def _fmt_usd(x) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "award amount not reported"
    return f"${x / 1e6:,.1f}M in award obligations"


def _ids(s: str, k: int) -> list[str]:
    v = [x for x in str(s or "").split(SEP) if x]
    return v[:k]


def _unknown(reason: str, **extra) -> dict:
    return {"status": UNKNOWN, "reason": reason, "candidates": [], "framing": FRAMING, "guardrails": GUARDRAILS,
            **extra}


def find_candidate_clinics(geography_id, condition, measurement=None, radius_km: float | None = None,
                           max_sites: int | None = None) -> dict:
    """Candidate implementation / study partner facilities for (geography, condition, measurement).

    Returns a dict with status ('ok', 'no_eligible_candidates' or UNKNOWN with a reason), the resolved inputs,
    pool sizes, and `candidates`: each with raw values and a separate rank for each of the six characteristics,
    plain-language reasons, object ids and the framing sentence. Never a single combined score.
    """
    cfg = _cfg_matching()
    radius_km = float(radius_km if radius_km is not None else cfg.get("default_radius_km", 50))
    max_sites = int(max_sites if max_sites is not None else cfg.get("max_sites_per_geography", 15))
    geo = resolve_geography(geography_id)
    if geo["status"] != "matched":
        return _unknown(f"geography: {geo['reason']}", geography=geo)
    cond = resolve_condition(condition)
    if cond["status"] != "matched":
        return _unknown(f"condition: {cond['reason']}", geography=geo, condition=cond)
    meas = None
    if measurement is not None and str(measurement).strip():
        meas = resolve_measurement(measurement)
        if meas["status"] != "matched":
            return _unknown(f"measurement: {meas['reason']}", geography=geo, condition=cond, measurement=meas)
    try:
        b = load_base()
    except FileNotFoundError as e:
        return _unknown(f"facility tables not built: {e}", geography=geo, condition=cond, measurement=meas)
    cid = cond["condition_id"]
    F = features(cid, meas)
    res = select(F, b.lat, b.lon, b.county, b.state, b.hpsa, b.nonmetro, geo, radius_km, max_sites)
    base_info = {"geography": geo, "condition": {**cond, "specialty_groups": condition_groups(cid)},
                 "measurement": meas if meas else {"status": "not_requested",
                                                   "note": "technology_experience and measurement implementer "
                                                           "groups are not evaluated"},
                 "radius_km": radius_km, "max_sites": max_sites, "framing": FRAMING, "guardrails": GUARDRAILS,
                 "characteristic_definitions": CHARACTERISTIC_DEFINITIONS,
                 "ranking_method": "Each characteristic is ranked separately among eligible facilities in the pool; "
                                   "the list is a round-robin over the six rankings (first of each, then second of "
                                   "each, ...), ties broken by distance. No combined score is computed.",
                 "provenance": {"tables": ["facilities", "facility_trials", "facility_nih_projects",
                                           "geographies", "configs/relevance.yaml", "configs/measurements.yaml"],
                                "facility_source_version": str(b.fac["source_version"].iloc[0]) if len(b.fac) else UNKNOWN}}
    pool = res["pool"]
    if len(pool) == 0:
        where = (f"inside {geo['name']}" if geo["geo_level"] == "state"
                 else f"within {radius_km:g} km of the internal point of {geo['name']} or inside it")
        return _unknown(f"no geocoded facility {where}", **base_info, n_facilities_in_pool=0)
    n_with = {c: int(len(res["orders"][c])) for c in CHARACTERISTICS}
    absent = [c for c in ("relevant_trial_history", "NIH_research_activity", "technology_experience")
              if n_with[c] == 0 and (c != "technology_experience" or meas)]
    out = {"status": "ok", **base_info, "n_facilities_in_pool": int(len(pool)), "n_eligible": res["n_eligible"],
           "n_with_characteristic": n_with, "characteristics_absent_in_pool": absent}
    if absent:
        what = cond["preferred_name"] + (f" x {meas['label']}" if meas else "")
        out["absence_note"] = (f"For {what} [condition:{cid}], no facility in the pool has " + ", ".join(absent)
                               + ": a possible measurement/research desert around this geography (absence in these "
                               "public registries, not proof that no such activity exists).")
    if not res["selected"]:
        out["status"] = "no_eligible_candidates"
        out["reason"] = ("facilities exist in the pool but none has a facility-level relevance signal (specialty "
                         "match, trial history, NIH activity, technology experience or community-health type) for "
                         "this condition/measurement")
        out["candidates"] = []
        return out
    ranks = {c: {int(f): i + 1 for i, f in enumerate(res["orders"][c])} for c in CHARACTERISTICS}
    cands = [_describe(i, k, b, F, res, ranks, n_with, geo, cond, meas) for k, i in enumerate(res["selected"])]
    out["candidates"] = cands
    out["n_returned"] = len(cands)
    return out


def _describe(i: int, order: int, b: Base, F: Features, res: dict, ranks: dict, n_with: dict, geo: dict,
              cond: dict, meas: dict | None) -> dict:
    r = b.fac.iloc[i]
    cid = cond["condition_id"]
    cname = cond["preferred_name"]
    gi = {g: j for j, g in enumerate(b.groups)}
    cg = [g for g in F.cond_groups if b.G[i, gi[g]]]
    mg = [g for g in F.meas_groups if b.G[i, gi[g]]]
    d = float(res["dist"][i])
    inside = (geo["geo_level"] == "state" and b.state[i] == geo["state_abbr"]) or (b.county[i] == geo["geo_id"])
    ft = b.ft[b.ft["fidx"] == i]
    fn = b.fn[b.fn["fidx"] == i]
    tc = ft[_has(ft["condition_ids_literal"], [cid])]
    tm = ft[_has(ft["measurement_classes"], meas["members"])] if meas else ft.iloc[0:0]
    fnc = fn[_has(fn["condition_ids_precision"], [cid])]
    reasons = []
    kind = r["primary_kind"]
    if cg or mg:
        txt = []
        if cg:
            txt.append(f"taxonomy groups that commonly see {cname}: {', '.join(cg)}")
        if mg and meas:
            txt.append(f"groups that could implement {meas['label']}: {', '.join(mg)}")
        gp = " (organisation NPI of a group practice)" if kind == "nppes_group_practice" else ""
        reasons.append("Self-reported NPPES " + "; ".join(txt) + gp +
                       ". A taxonomy code does not mean the facility evaluates or treats " + cname + ".")
    if F.trials_cond[i] > 0:
        ex = ", ".join("trial:" + x for x in tc["nct_id"].drop_duplicates().head(3))
        reasons.append(f"Registered site of {int(F.trials_cond[i])} U.S. {cname} trial(s) in ClinicalTrials.gov "
                       f"({int(F.trials_cond_open[i])} with open status; latest start year "
                       f"{int(F.trials_cond_latest[i]) if F.trials_cond_latest[i] else 'not reported'}; e.g. {ex}). "
                       "Registration is a research-activity signal, not evidence of effectiveness.")
    if F.nih_core[i] > 0:
        reasons.append(f"Awardee organisation of {int(F.nih_core[i])} NIH RePORTER core project(s) whose title/"
                       f"abstract mention {cname} ({int(F.nih_active_core[i])} active; "
                       f"{_fmt_usd(F.nih_obligations[i])}, FY2015-FY2026, obligations not expenditure).")
    if meas and F.tech[i] > 0:
        cls = sorted({c for s in tm["measurement_classes"] for c in str(s).split(SEP) if c in meas["members"]})
        reasons.append(f"Registered site of {int(F.tech[i])} trial(s) of any of the registry's target conditions (not "
                       f"only {cname}) whose registered text mentions {', '.join(cls)}; {int(F.tech_cond[i])} of "
                       f"them are {cname} trials. Experience with the measurement in a study, not evidence that it "
                       "works.")
    acc = []
    if r["is_hrsa_site"]:
        acc.append("HRSA " + ("FQHC look-alike" if r["hrsa_is_look_alike"] else "health center program") + " site"
                   + (" (mobile van; location is its registered address)" if r["hrsa_location_type"] == "Mobile Van"
                      else ""))
    elif r["is_fqhc_or_lookalike"]:
        acc.append("self-reported FQHC taxonomy (NPPES)")
    if r["is_community_type_facility"]:
        acc.append("community / rural-health / critical-access / public-health clinic taxonomy")
    if res["access"][i] - F.intrinsic_access[i] > 0:
        if b.hpsa[i]:
            acc.append(f"its county has a primary-care HPSA designation ({r['county_pc_hpsa_coverage']})")
        if b.nonmetro[i]:
            acc.append(f"its county is nonmetro (RUCC {r['rucc_2023']})")
    if acc:
        reasons.append("Community access: " + "; ".join(acc) + ".")
    loc = (f"{d:.1f} km from the internal point of {geo['name']}" + (" (inside it)" if inside else "") +
           f"; location precision: {r['geocode_precision']}.")
    reasons.append(loc[0].upper() + loc[1:])

    def ch(c, value, detail):
        return {"value": value, "rank": ranks[c].get(i), "n_ranked": n_with[c], "detail": detail}

    chars = {
        "clinical_specialty_match": ch("clinical_specialty_match", int(F.n_union_groups[i]),
                                       {"condition_groups_matched": cg, "measurement_groups_matched": mg}),
        "relevant_trial_history": ch("relevant_trial_history", int(F.trials_cond[i]),
                                     {"open_status": int(F.trials_cond_open[i]),
                                      "latest_start_year": int(F.trials_cond_latest[i]) or None}),
        "NIH_research_activity": ch("NIH_research_activity", int(F.nih_core[i]),
                                    {"active_core_projects": int(F.nih_active_core[i]),
                                     "award_obligations_usd": (None if np.isnan(F.nih_obligations[i])
                                                               else float(F.nih_obligations[i]))}),
        "technology_experience": ch("technology_experience", int(F.tech[i]) if meas else UNKNOWN,
                                    {"for_this_condition": int(F.tech_cond[i]) if meas else UNKNOWN}),
        "community_access": ch("community_access", int(res["access"][i]),
                               {"fqhc_or_lookalike": bool(r["is_fqhc_or_lookalike"]),
                                "fqhc_evidence": r["fqhc_evidence"] or None,
                                "community_type_facility": bool(r["is_community_type_facility"]),
                                "county_pc_hpsa_area": bool(b.hpsa[i]), "county_nonmetro": bool(b.nonmetro[i])}),
        "distance_to_target_population": ch("distance_to_target_population", round(d, 2),
                                            {"inside_geography": bool(inside),
                                             "geocode_precision": r["geocode_precision"]}),
    }
    return {
        "selection_order": order + 1, "selected_via": res["via"][i],
        "facility_id": r["facility_id"], "object_id": r["object_id"], "facility_name": r["facility_name"],
        "primary_kind": kind, "facility_kinds": r["facility_kinds"], "sources": r["sources"],
        "address": r["address"], "city": r["city"], "state": r["state"], "zip5": r["zip5"],
        "county_fips": r["county_fips"], "lat": float(r["lat"]), "lon": float(r["lon"]),
        "geocode_precision": r["geocode_precision"], "distance_km": round(d, 2), "inside_geography": bool(inside),
        "characteristics": chars, "reasons": reasons, "framing": FRAMING,
        "resolution_tier": r["resolution_tier"],
        "source_object_ids": {
            "npi": _ids(r["npi_object_ids"], 5), "n_npi": len(_ids(r["npi_object_ids"], 10 ** 6)),
            "hrsa_site": _ids(r["hrsa_site_object_ids"], 5),
            "trial": sorted("trial:" + tc["nct_id"].unique())[:10] if len(tc) else [],
            "n_trial_condition": int(tc["nct_id"].nunique()),
            "nih": sorted("nih:" + fnc["appl_id"].astype(str).unique())[:10] if len(fnc) else [],
            "n_nih_condition": int(fnc["appl_id"].nunique()),
        },
    }


def candidates_frame(result: dict) -> pd.DataFrame:
    """Flatten a find_candidate_clinics result to one row per candidate (values and ranks side by side)."""
    rows = []
    for c in result.get("candidates", []):
        row = {k: c[k] for k in ("selection_order", "selected_via", "facility_id", "object_id", "facility_name",
                                 "primary_kind", "city", "state", "zip5", "county_fips", "geocode_precision",
                                 "distance_km", "inside_geography", "resolution_tier")}
        for ch, v in c["characteristics"].items():
            row[f"{ch}__value"] = v["value"]
            row[f"{ch}__rank"] = v["rank"]
            row[f"{ch}__n_ranked"] = v["n_ranked"]
        row["reasons"] = " | ".join(c["reasons"])
        row["framing"] = c["framing"]
        rows.append(row)
    return pd.DataFrame(rows)
