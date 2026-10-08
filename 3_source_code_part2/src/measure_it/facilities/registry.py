"""Unified facility table, clinic_registry and research_site_registry (SPEC §F, §G).

Sources (all already processed upstream; read their DATA_AUDIT.md):
  * NPPES organisation NPIs in a specialty/facility/research group (`providers`). Organisation NPIs under
    physician taxonomies are GROUP PRACTICES, not facilities; taxonomies are self-reported.
  * HRSA health-center service-delivery and look-alike sites (`facilities__hrsa`), point geocodes.
  * U.S. ClinicalTrials.gov facilities (`trial_sites`), excluding sponsor placeholders ("Research Site"),
    unnamed sites and non-facility locations (participants' homes, virtual/online studies).
  * NIH RePORTER awardee organisations (`nih_projects`).

Entity resolution is conservative: a link needs a normalised-name agreement AND the same ZIP5, or the same
city+state within 25 km (with >= 2 informative name tokens); never name-only, never across states. Every link
keeps `match_method` and `match_confidence`; union-find clustering allows at most one HRSA site and one NIH
organisation per facility and no facility wider than 25 km. `facility_id` is the anchor record id; its object id
is `facility:<facility_id>`.

Outputs (data/processed): facilities, facility_source_links, facility_trials, facility_nih_projects,
clinic_registry, research_site_registry. Validation tables go to results/tables/facility_*.csv.

Reproduce: `uv run python -m measure_it.facilities.registry`
"""
from __future__ import annotations

import argparse
import json
import re
import unicodedata
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from rapidfuzz import fuzz
from rapidfuzz.process import cdist

from ..config import SEED, TABLES, UNKNOWN, load_config, utc_now_iso
from ..geography.crosswalk import STATE_ABBR_TO_FIPS, zip_to_geo
from ..provenance import PROVENANCE_COLUMNS
from ..store import read_table, write_table

PRODUCER = "facilities.registry"
SEP = "|"
MAX_SPAN_KM = 25.0
WORKERS = 6

# --------------------------------------------------------------------------------------------------------------
# name normalisation
# --------------------------------------------------------------------------------------------------------------

ABBREV = {
    "univ": "university", "ctr": "center", "cntr": "center", "centre": "center", "centers": "center",
    "hosp": "hospital", "med": "medical", "hlth": "health", "inst": "institute", "natl": "national",
    "assn": "association", "assoc": "associates", "svcs": "services", "svc": "service", "sys": "system",
    "dept": "department", "sch": "school", "coll": "college", "mem": "memorial", "reg": "regional",
    "cmty": "community", "comm": "community", "fam": "family", "mt": "mount", "ft": "fort", "st": "saint",
    "ste": "sainte", "children": "childrens", "chldrns": "childrens", "hlthcare": "healthcare",
    "clin": "clinic", "clinics": "clinic", "hospitals": "hospital", "labs": "laboratory", "lab": "laboratory",
    "rsch": "research", "res": "research", "sciences": "science", "scis": "science", "sci": "science",
}
LEGAL_SUFFIXES = {"inc", "llc", "pc", "pa", "pllc", "ltd", "corp", "corporation", "co", "company", "lp", "llp",
                  "plc", "incorporated", "limited", "sc", "psc", "ps", "md", "apc", "apmc", "lc", "pllp", "dba"}
STOP = {"of", "the", "and", "at", "in", "for", "a", "an", "on", "to", "by", "dba", "aka"}
# facility-type words that do not identify a particular institution
GENERIC = {
    "clinic", "clinica", "center", "medical", "medicine", "health", "healthcare", "hospital", "care", "services",
    "service", "group", "groups", "associates", "associate", "association", "practice", "practices", "family",
    "physicians", "physician", "research", "clinical", "institute", "institutes", "system", "systems",
    "network", "networks", "community", "department", "office", "offices", "site", "sites", "university",
    "college", "school", "partners", "specialists", "specialty", "specialties", "primary", "general", "regional",
    "foundation", "trial", "trials", "study", "studies", "professional", "professionals", "management",
    "consultants", "program", "programs", "unit", "division", "campus", "main", "north", "south", "east", "west",
    "central", "national", "county", "memorial", "saint", "sainte", "childrens", "science", "laboratory",
    "healthsystem", "affiliates", "affiliated", "corporate", "llc", "inc", "pc", "center", "centers", "new",
}
NONFACILITY_RE = re.compile(
    r"^\s*site\s+reference\s+id\b|\bcontact\b.*\bfor\b.*\blocations?\b|"
    r"^\s*(?:site|location|center|centre)\s*(?:#|no|number)?\s*\d+\s*$|"
    r"\b(?:subjects?|participants?|patients?)'?s?\s+(?:own\s+)?(?:homes?|residences?)\b|\bvirtual\b|"
    r"\bremote\s+(?:trial|study|site)\b|\bonline\b|\bdecentrali[sz]ed\b|\bhome[- ]based\s+study\b|"
    r"\btelemedicine\s+study\b|\bno\s+physical\s+site\b|digital\s+research\s+platform|\bproofpilot\b",
    re.IGNORECASE)


def _ascii_lower(text) -> str:
    if text is None or (isinstance(text, float) and np.isnan(text)) or text is pd.NA:
        return ""
    return unicodedata.normalize("NFKD", str(text)).encode("ascii", "ignore").decode().lower()


def name_tokens(text) -> list[str]:
    """Normalised tokens of a facility/organisation name (abbreviations expanded, legal suffixes dropped)."""
    s = _ascii_lower(text)
    s = re.sub(r"['’`.]", "", s)         # children's -> childrens; P.A. -> pa; St. -> st
    s = s.replace("&", " and ").replace("+", " and ")
    s = re.sub(r"\bd/?b/?a\b.*$", "", s)      # keep the legal name before "d/b/a ..."
    s = re.sub(r"[^a-z0-9]+", " ", s)
    toks = [ABBREV.get(t, t) for t in s.split()]
    while toks and toks[-1] in LEGAL_SUFFIXES:   # legal suffixes only at the end ("... PA", "... INC")
        toks.pop()
    while toks and toks[0] == "the":
        toks.pop(0)
    return toks


def norm_name(text) -> str:
    return " ".join(name_tokens(text))


def informative_tokens(tokens) -> frozenset:
    """Tokens that can identify an institution: not generic facility words, not stop words, >= 2 chars."""
    return frozenset(t for t in tokens if t not in GENERIC and t not in STOP and len(t) >= 2)


def norm_city(text) -> str:
    s = re.sub(r"['’`.]", "", _ascii_lower(text))
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return " ".join(ABBREV.get(t, t) if t in ("st", "ft", "mt", "ste") else t for t in s.split())


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = (np.radians(np.asarray(x, dtype=float)) for x in (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


# --------------------------------------------------------------------------------------------------------------
# link rules
# --------------------------------------------------------------------------------------------------------------

CONFIDENCE = {
    "shared_npi": 0.99,
    "nppes_same_name_same_zip": 0.95,
    "exact_name_same_zip": 0.95,
    "hrsa_org_name_same_zip": 0.85,
    "fuzzy_name_same_zip": 0.85,
    "contained_name_same_zip": 0.80,
    "exact_name_same_city": 0.80,
    "fuzzy_name_same_city": 0.75,
    "contained_name_same_city": 0.70,
}
FUZZY_ZIP, FUZZY_CITY = 92.0, 95.0


def confidence_tier(c: float) -> str:
    if pd.isna(c):
        return "single_record"
    return "high" if c >= 0.9 else ("medium" if c >= 0.75 else "low")


def _digit_tokens(norm: str) -> frozenset:
    return frozenset(t for t in norm.split() if any(ch.isdigit() for ch in t))


def _in_order(short_seq: list, long_seq: list) -> bool:
    """short_seq is a subsequence of long_seq restricted to short's tokens (same order, gaps allowed)."""
    it = iter([t for t in long_seq if t in set(short_seq)])
    return all(t in it for t in short_seq)


def classify_pair(na: str, nb: str, ia: frozenset, ib: frozenset, level: str, place_tokens: frozenset,
                  sort_ratio: float | None = None, city_tokens: frozenset = frozenset()) -> str | None:
    """Apply the name rules to one pair already known to share a ZIP5 (level='zip') or a city+state within
    MAX_SPAN_KM (level='city'). Returns the method name or None. Pure function (unit-tested).

    Fuzzy and containment links are refused when the names carry different numbers ('Site 146' vs 'Site 147' are
    different sites). Containment needs the shorter name's informative tokens in the same order in the longer name,
    extra informative tokens that are only place names, and at least one informative token that is not the city."""
    if not na or not nb:
        return None
    min_info = 1 if level == "zip" else 2
    if na == nb:
        return f"exact_name_same_{level}" if len(ia) >= min_info else None
    if _digit_tokens(na) != _digit_tokens(nb):
        return None
    shared = ia & ib
    if sort_ratio is None:
        sort_ratio = fuzz.token_sort_ratio(na, nb)
    thr = FUZZY_ZIP if level == "zip" else FUZZY_CITY
    if sort_ratio >= thr and len(shared) >= min_info and min(len(ia), len(ib)) >= min_info:
        return f"fuzzy_name_same_{level}"
    (short, ns), (long_, nl) = sorted([(ia, na), (ib, nb)], key=lambda x: (len(x[0]), x[1]))
    if (len(short) >= min_info and (short - city_tokens) and short <= long_ and (long_ - short) <= place_tokens
            and _in_order([t for t in ns.split() if t in short], [t for t in nl.split() if t in short])):
        return f"contained_name_same_{level}"
    return None


# --------------------------------------------------------------------------------------------------------------
# source records
# --------------------------------------------------------------------------------------------------------------

STATE_NAME_TO_ABBR: dict[str, str] = {}


def _state_names() -> dict[str, str]:
    if not STATE_NAME_TO_ABBR:
        g = read_table("geographies", columns=["geo_level", "state_name", "state_abbr"])
        g = g[g["geo_level"] == "state"]
        STATE_NAME_TO_ABBR.update({str(n).lower(): a for n, a in zip(g["state_name"], g["state_abbr"])})
    return STATE_NAME_TO_ABBR


def _s(x) -> str | None:
    if x is None or x is pd.NA or (isinstance(x, float) and np.isnan(x)):
        return None
    x = str(x).strip()
    return x or None


def _join(values, sep: str = SEP) -> str:
    vals = [str(v) for v in values if _s(v) is not None]
    return sep.join(dict.fromkeys(vals))


def load_nppes_records() -> pd.DataFrame:
    cols = ["npi", "entity_type", "name", "organization_name", "is_organization_subpart", "parent_organization_name",
            "primary_taxonomy_code", "primary_taxonomy_display_name", "specialty_groups", "primary_specialty_groups",
            "group_types", "address_line1", "city", "state", "state_fips", "zip5", "lat", "lon", "county_fips",
            "geocode_method", "source_geographic_resolution", "last_update_date", "in_50_states_dc",
            "source_name", "source_version", "retrieved_at"]
    p = read_table("providers", columns=cols)
    p = p[p["entity_type"] == 2].copy()
    gt = p["group_types"].fillna("")
    p["kind"] = np.where(gt.str.contains("facility"), "nppes_facility",
                         np.where(gt.str.contains("specialty"), "nppes_group_practice", "nppes_research"))
    p["record_id"] = "npi:" + p["npi"].astype(str)
    p["source"] = "nppes"
    p["native_id"] = p["npi"].astype(str)
    p["alt_name"] = None
    p["precision"] = np.where(p["lat"].notna(), "zcta_centroid", "none")
    return p


def load_hrsa_records() -> pd.DataFrame:
    h = read_table("facilities__hrsa")
    h = h.copy()
    h["record_id"] = "hrsa_site:" + h["site_bphc_number"].astype(str)
    h["source"] = "hrsa"
    h["native_id"] = h["site_bphc_number"].astype(str)
    h["name"] = h["facility_name"]
    h["alt_name"] = h["health_center_name"]
    h["state"] = h["state_abbr"]
    h["kind"] = np.where(h["is_look_alike"].astype(bool), "hrsa_look_alike_site", "hrsa_health_center_site")
    res = h["source_geographic_resolution"].astype(str)
    h["precision"] = np.where(h["lat"].isna(), "none", np.where(res == "point", "point", "zcta_centroid"))
    return h


def load_ctgov_site_rows() -> tuple[pd.DataFrame, dict]:
    """U.S. identifiable ClinicalTrials.gov site rows (one per NCT x site) with the upstream site_key."""
    from ..ingestion.clinicaltrials import site_key

    s = read_table("trial_sites")
    stats = {"trial_site_rows_us": int((s["country_group"] == "US").sum())}
    s = s[s["country_group"] == "US"].copy()
    stats["us_rows_sponsor_placeholder"] = int(s["facility_generic"].astype(bool).sum())
    stats["us_rows_facility_missing"] = int(s["facility_missing"].astype(bool).sum())
    s = s[s["facility_identifiable"].astype(bool)]
    nonfac = s["facility"].fillna("").str.contains(NONFACILITY_RE)
    stats["us_rows_nonfacility_location"] = int(nonfac.sum())
    stats["us_rows_nonfacility_examples"] = sorted(s.loc[nonfac, "facility"].unique().tolist())[:15]
    s = s[~nonfac]
    # a "facility" name that is only a place (e.g. "Virginia Beach" in Virginia Beach) identifies no facility
    names = _state_names()
    extra = {"usa", "us", "united", "states", "america", "metro", "area", "greater"}
    place_only = [bool(set(name_tokens(f))) and set(name_tokens(f)) <= (set(norm_city(c).split()) | set(
        str(st).lower().split()) | set(names.get(str(st).lower(), "").lower().split()) | extra)
        for f, c, st in zip(s["facility"], s["city"], s["state"])]
    place_only = np.array(place_only, dtype=bool)
    stats["us_rows_place_name_only"] = int(place_only.sum())
    stats["us_rows_place_name_only_examples"] = sorted(s.loc[place_only, "facility"].unique().tolist())[:15]
    s = s[~place_only]
    stats["us_rows_identifiable_kept"] = int(len(s))
    s["site_key"] = [site_key(*r) for r in s[["facility", "city", "state", "country"]].itertuples(index=False)]
    return s, stats


def load_ctgov_records(site_rows: pd.DataFrame) -> pd.DataFrame:
    """One record per site_key; ZIP-centroid location when the ZIP resolves in the registered state,
    otherwise the ClinicalTrials.gov geoPoint (a city centroid)."""
    def mode(x):
        x = x.dropna()
        return x.mode().iloc[0] if len(x) else None

    g = site_rows.groupby("site_key", sort=True)
    rec = pd.DataFrame({
        "facility": g["facility"].agg(mode), "city": g["city"].agg(mode), "state_name": g["state"].agg(mode),
        "zip_raw": g["zip"].agg(mode), "geo_lat": g["lat"].agg(mode), "geo_lon": g["lon"].agg(mode),
        "geo_res": g["source_geographic_resolution"].agg(mode), "county_fips_ctgov": g["county_fips"].agg(mode),
        "state_fips_ctgov": g["state_fips"].agg(mode), "retrieved_at": g["retrieved_at"].max(),
        "source_version": g["source_version"].agg(mode), "source_name": g["source_name"].agg(mode),
        "n_site_rows": g.size(),
    }).reset_index()
    names = _state_names()
    rec["state"] = rec["state_name"].map(lambda x: names.get(str(x).lower()) if _s(x) else None)
    z = zip_to_geo(rec["zip_raw"].fillna(""))
    rec["zip5"] = z["zip5"].values
    ok = z["lat"].notna().values & ((z["state_fips"].values == rec["state_fips_ctgov"].values) |
                                    pd.isna(z["state_fips"].values))
    rec["lat"] = np.where(ok, z["lat"].values, rec["geo_lat"].values).astype(float)
    rec["lon"] = np.where(ok, z["lon"].values, rec["geo_lon"].values).astype(float)
    rec["precision"] = np.where(ok, "zcta_centroid",
                                np.where(rec["geo_lat"].notna(),
                                         np.where(rec["geo_res"] == "zcta_centroid", "zcta_centroid", "city_centroid"),
                                         "none"))
    rec["geocode_method"] = np.where(ok, "ctgov_zip5_as_zcta_internal_point",
                                     np.where(rec["geo_lat"].notna(), "ctgov_geopoint_city_level", "not_geocoded"))
    zc = pd.Series(z["county_fips"].values)
    rec["county_fips"] = np.where(ok & zc.notna().values, zc.values, rec["county_fips_ctgov"].values)
    rec["state_fips"] = rec["state_fips_ctgov"]
    rec["record_id"] = "ctgov_site:" + rec["site_key"]
    rec["native_id"] = rec["site_key"]
    rec["source"] = "ctgov"
    rec["name"] = rec["facility"]
    rec["alt_name"] = None
    rec["kind"] = "ctgov_trial_site"
    return rec


def load_nih_records() -> tuple[pd.DataFrame, pd.DataFrame]:
    n = read_table("nih_projects")
    n = n.copy()
    n["org_key"] = n["org_ipf_code"].astype("string").fillna("name:" + n["org_name"].astype("string").fillna(UNKNOWN))
    us = n["org_country"].fillna("").str.upper() == "UNITED STATES"
    n = n[us]
    n = n.sort_values(["fiscal_year", "appl_id"])
    g = n.groupby("org_key", sort=True)
    last = g.tail(1).set_index("org_key")
    rec = pd.DataFrame({
        "org_key": last.index, "name": last["org_name"].values, "city": last["org_city"].values,
        "state": last["org_state"].values, "zip5": last["zip5"].values, "lat": last["lat"].values,
        "lon": last["lon"].values, "county_fips": last["county_fips"].values, "state_fips": last["state_fips"].values,
        "geocode_method": last["geocode_method"].values, "res": last["source_geographic_resolution"].values,
        "org_type_name": last["org_type_name"].values, "retrieved_at": last["retrieved_at"].values,
        "source_version": last["source_version"].values, "source_name": last["source_name"].values,
    })
    rec["n_appl_ids"] = g.size().reindex(rec["org_key"]).values
    rec["record_id"] = "nih_org:" + rec["org_key"].astype(str)
    rec["native_id"] = rec["org_key"].astype(str)
    rec["source"] = "nih"
    rec["alt_name"] = None
    rec["kind"] = "nih_awardee_org"
    rec["precision"] = np.where(rec["lat"].isna(), "none", np.where(rec["res"] == "point", "point", "zcta_centroid"))
    return rec, n


NODE_COLS = ["record_id", "source", "native_id", "kind", "name", "alt_name", "city", "state", "zip5", "lat", "lon",
             "precision", "county_fips", "state_fips"]


def _prepare_nodes(frames: list[pd.DataFrame]) -> pd.DataFrame:
    nodes = pd.concat([f[NODE_COLS] for f in frames], ignore_index=True)
    for c in ("name", "alt_name", "city", "state", "zip5", "county_fips", "state_fips"):
        nodes[c] = nodes[c].astype(object).where(nodes[c].notna(), None)
    nodes["state"] = nodes["state"].map(lambda x: str(x).upper() if _s(x) else None)
    nodes["zip5"] = nodes["zip5"].map(lambda x: str(x)[:5] if _s(x) and re.match(r"^\d{5}", str(x)) else None)
    toks = nodes["name"].map(name_tokens)
    nodes["norm"] = toks.map(" ".join)
    nodes["info"] = toks.map(informative_tokens)
    nodes["alt_norm"] = nodes["alt_name"].map(lambda x: norm_name(x) if _s(x) else "")
    nodes["city_n"] = nodes["city"].map(norm_city)
    nodes["lat"] = pd.to_numeric(nodes["lat"], errors="coerce").astype(float)
    nodes["lon"] = pd.to_numeric(nodes["lon"], errors="coerce").astype(float)
    return nodes


# --------------------------------------------------------------------------------------------------------------
# candidate edges
# --------------------------------------------------------------------------------------------------------------

def _place_tokens(state_abbr: str | None, *cities: str) -> frozenset:
    toks = set()
    for c in cities:
        toks.update((c or "").split())
    if state_abbr:
        toks.add(state_abbr.lower())
        fips = STATE_ABBR_TO_FIPS.get(state_abbr)
        for n, a in _state_names().items():
            if a == state_abbr:
                toks.update(n.split())
        if fips:
            toks.add(fips)
    return frozenset(toks)


def nppes_group_edges(nodes: pd.DataFrame) -> pd.DataFrame:
    """Within-NPPES: identical normalised name and ZIP5 -> same facility (subparts/duplicate NPIs)."""
    p = nodes[(nodes["source"] == "nppes") & (nodes["norm"] != "") & nodes["zip5"].notna()]
    p = p.assign(npi=p["native_id"].astype(np.int64))
    anchor = p.groupby(["norm", "zip5"])["npi"].transform("min")
    m = p[p["npi"] != anchor]
    return pd.DataFrame({"a": m["record_id"].values, "b": "npi:" + anchor[m.index].astype(str).values,
                         "method": "nppes_same_name_same_zip", "dist_km": 0.0})


def shared_npi_edges(nodes: pd.DataFrame, hrsa: pd.DataFrame) -> pd.DataFrame:
    npis = set(nodes.loc[nodes["source"] == "nppes", "native_id"])
    h = hrsa[hrsa["site_npi"].notna()]
    h = h[h["site_npi"].astype(str).isin(npis)]
    return pd.DataFrame({"a": h["record_id"].values, "b": "npi:" + h["site_npi"].astype(str).values,
                         "method": "shared_npi", "dist_km": np.nan})


def _nppes_group_representatives(nodes: pd.DataFrame, group_edges: pd.DataFrame) -> pd.DataFrame:
    """Cross-source matching compares against one node per NPPES (name, ZIP) group."""
    members = set(group_edges["a"])
    return nodes[~((nodes["source"] == "nppes") & nodes["record_id"].isin(members))]


def _pair_allowed(sa: str, sb: str) -> bool:
    if sa == sb:
        return sa == "ctgov"   # ctgov spelling variants; NPPES handled by grouping; never HRSA-HRSA / NIH-NIH
    return True


def exact_zip_edges(nodes: pd.DataFrame) -> pd.DataFrame:
    base = nodes[(nodes["norm"] != "") & nodes["zip5"].notna()]
    names = base[["record_id", "source", "norm", "zip5", "info"]].assign(method="exact_name_same_zip")
    alt = base[base["alt_norm"] != ""][["record_id", "source", "alt_norm", "zip5"]].rename(columns={"alt_norm": "norm"})
    alt = alt.assign(info=alt["norm"].map(lambda s: informative_tokens(s.split())), method="hrsa_org_name_same_zip")
    left = pd.concat([names, alt], ignore_index=True)
    left = left[left["info"].map(len) >= 1]
    right = names[["record_id", "source", "norm", "zip5"]]
    m = left.merge(right, on=["norm", "zip5"], suffixes=("_a", "_b"))
    m = m[m["record_id_a"] < m["record_id_b"]]
    m = m[[_pair_allowed(a, b) for a, b in zip(m["source_a"], m["source_b"])]]
    m = m[~((m["method"] == "hrsa_org_name_same_zip") & (m["source_b"] == "hrsa"))]
    # the HRSA organisation name can sit on either side of the merge
    alt_r = alt[["record_id", "source", "norm", "zip5"]]
    m2 = names.merge(alt_r, on=["norm", "zip5"], suffixes=("_a", "_b"))
    m2 = m2[(m2["record_id_a"] < m2["record_id_b"]) & (m2["source_a"] != "hrsa") & (m2["info"].map(len) >= 1)]
    m2 = m2.assign(method="hrsa_org_name_same_zip")
    out = pd.concat([m[["record_id_a", "record_id_b", "method"]], m2[["record_id_a", "record_id_b", "method"]]])
    out = out.rename(columns={"record_id_a": "a", "record_id_b": "b"}).assign(dist_km=0.0)
    return out.drop_duplicates(["a", "b", "method"])


def _block_pairs(block: pd.DataFrame, level: str, query_sources: frozenset | None = None) -> list[tuple]:
    """Fuzzy/containment candidates inside one ZIP or city block (only pairs involving a non-NPPES record)."""
    out = []
    src = block["source"].to_numpy()
    norms = block["norm"].tolist()
    info = block["info"].tolist()
    rid = block["record_id"].to_numpy()
    zips = block["zip5"].to_numpy()
    lat = block["lat"].to_numpy()
    lon = block["lon"].to_numpy()
    q = np.flatnonzero(src != "nppes") if query_sources is None else np.flatnonzero(np.isin(src, list(query_sources)))
    if len(q) == 0 or len(block) < 2:
        return out
    ratios = cdist([norms[i] for i in q], norms, scorer=fuzz.token_sort_ratio, workers=1)
    setr = cdist([norms[i] for i in q], norms, scorer=fuzz.token_set_ratio, workers=1)
    cities = frozenset(block["city_n"].iloc[0].split()) if level == "city" else frozenset()
    for qi, i in enumerate(q):
        cand = np.flatnonzero((ratios[qi] >= min(FUZZY_ZIP, FUZZY_CITY)) | (setr[qi] >= 100))
        for j in cand:
            if j == i or not norms[i] or not norms[j] or norms[i] == norms[j]:
                continue
            if src[j] != "nppes" and rid[j] < rid[i] and (query_sources is None or src[j] in query_sources):
                continue  # each unordered pair once
            if not _pair_allowed(src[i], src[j]):
                continue
            if level == "city":
                if zips[i] is not None and zips[i] == zips[j]:
                    continue  # handled at ZIP level
                d = float(haversine_km(lat[i], lon[i], lat[j], lon[j])) if not (
                    np.isnan(lat[i]) or np.isnan(lat[j])) else np.nan
                if not np.isnan(d) and d > MAX_SPAN_KM:
                    continue
            else:
                d = 0.0
            ctoks = frozenset(block["city_n"].iloc[i].split()) | frozenset(block["city_n"].iloc[j].split())
            place = _place_tokens(block["state"].iloc[0], block["city_n"].iloc[i], block["city_n"].iloc[j]) | cities
            meth = classify_pair(norms[i], norms[j], info[i], info[j], level, place, float(ratios[qi, j]), ctoks)
            if meth:
                out.append((rid[i], rid[j], meth, d))
    return out


def fuzzy_edges(nodes: pd.DataFrame, level: str, query_sources: frozenset | None = None) -> pd.DataFrame:
    if level == "zip":
        base = nodes[(nodes["norm"] != "") & nodes["zip5"].notna()]
        keys = ["zip5"]
    else:
        base = nodes[(nodes["norm"] != "") & (nodes["city_n"] != "") & nodes["state"].notna()]
        keys = ["state", "city_n"]
    isq = (base["source"] != "nppes") if query_sources is None else base["source"].isin(list(query_sources))
    has_q = isq.groupby([base[k] for k in keys]).transform("any")
    base = base[has_q]
    rows = []
    for _, blk in base.groupby(keys, sort=True):
        if len(blk) > 1:
            rows.extend(_block_pairs(blk.reset_index(drop=True), level, query_sources))
    return pd.DataFrame(rows, columns=["a", "b", "method", "dist_km"])


def exact_city_edges(nodes: pd.DataFrame) -> pd.DataFrame:
    base = nodes[(nodes["norm"] != "") & (nodes["city_n"] != "") & nodes["state"].notna()]
    base = base[base["info"].map(len) >= 2]
    has_q = (base["source"] != "nppes").groupby([base["norm"], base["state"], base["city_n"]]).transform("any")
    base = base[has_q]
    m = base.merge(base, on=["norm", "state", "city_n"], suffixes=("_a", "_b"))
    m = m[m["record_id_a"] != m["record_id_b"]]
    m = m[(m["source_a"] != "nppes") | (m["source_b"] != "nppes")]
    m = m[~((m["source_a"] != "nppes") & (m["source_b"] != "nppes") & (m["record_id_a"] > m["record_id_b"]))]
    m = m[~((m["source_a"] == "nppes") & (m["source_b"] != "nppes"))]  # keep orientation non-NPPES -> other
    m = m[[_pair_allowed(a, b) for a, b in zip(m["source_a"], m["source_b"])]]
    m = m[~(m["zip5_a"].notna() & (m["zip5_a"] == m["zip5_b"]))]
    d = haversine_km(m["lat_a"], m["lon_a"], m["lat_b"], m["lon_b"])
    m = m.assign(dist_km=d)
    m = m[~(m["dist_km"] > MAX_SPAN_KM)]
    return pd.DataFrame({"a": m["record_id_a"].values, "b": m["record_id_b"].values,
                         "method": "exact_name_same_city", "dist_km": m["dist_km"].values})


def mutual_best(edges: pd.DataFrame, source_of: dict) -> pd.DataFrame:
    """City-level edges: keep an edge only if it is each endpoint's best link to the other endpoint's source."""
    if edges.empty:
        return edges
    e = edges.copy()
    e["conf"] = e["method"].map(CONFIDENCE)
    e["d"] = e["dist_km"].fillna(0.0)
    e["sa"] = e["a"].map(source_of)
    e["sb"] = e["b"].map(source_of)
    both = pd.concat([
        e.assign(node=e["a"], other_src=e["sb"]), e.assign(node=e["b"], other_src=e["sa"])], ignore_index=False)
    both["eid"] = both.index
    both = both.sort_values(["node", "other_src", "conf", "d", "a", "b"], ascending=[True, True, False, True, True, True])
    best = both.groupby(["node", "other_src"]).head(1)["eid"]
    counts = best.value_counts()
    keep = counts[counts == 2].index
    return edges.loc[keep]


def candidate_edges(nodes: pd.DataFrame, hrsa: pd.DataFrame, *, use_shared_npi: bool = True) -> pd.DataFrame:
    grp = nppes_group_edges(nodes)
    reps = _nppes_group_representatives(nodes, grp)
    parts = [grp]
    if use_shared_npi:
        sn = shared_npi_edges(nodes, hrsa)
        # point the HRSA site at the NPPES group anchor of that NPI
        to_anchor = dict(zip(grp["a"], grp["b"]))
        sn["b"] = sn["b"].map(lambda r: to_anchor.get(r, r))
        parts.append(sn)
    parts.append(exact_zip_edges(reps))
    parts.append(fuzzy_edges(reps, "zip"))
    city = pd.concat([exact_city_edges(reps), fuzzy_edges(reps, "city")], ignore_index=True)
    source_of = dict(zip(nodes["record_id"], nodes["source"]))
    parts.append(mutual_best(city, source_of))
    e = pd.concat(parts, ignore_index=True)
    e["conf"] = e["method"].map(CONFIDENCE)
    e = e.sort_values(["conf"], ascending=False).drop_duplicates(["a", "b"])
    return e.reset_index(drop=True)


# --------------------------------------------------------------------------------------------------------------
# constrained clustering
# --------------------------------------------------------------------------------------------------------------

@dataclass
class _Cluster:
    n_hrsa: int = 0
    n_nih: int = 0
    coords: list = field(default_factory=list)
    states: set = field(default_factory=set)


def cluster(nodes: pd.DataFrame, edges: pd.DataFrame) -> tuple[dict, dict, dict]:
    """Union-find over edges in descending confidence, with the HRSA/NIH/state/span constraints.

    Returns (root_of, accepted_edge_of_record, rejected_counts)."""
    parent: dict[str, str] = {}
    info: dict[str, _Cluster] = {}
    src = dict(zip(nodes["record_id"], nodes["source"]))
    states = nodes["state"] if "state" in nodes.columns else pd.Series([None] * len(nodes), index=nodes.index)
    for rid, s, la, lo, st in zip(nodes["record_id"], nodes["source"], nodes["lat"], nodes["lon"], states):
        parent[rid] = rid
        info[rid] = _Cluster(n_hrsa=int(s == "hrsa"), n_nih=int(s == "nih"),
                             coords=[] if np.isnan(la) else [(round(la, 4), round(lo, 4))],
                             states={st} if _s(st) else set())

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    e = edges.assign(d=edges["dist_km"].fillna(0.0)).sort_values(["conf", "d", "a", "b"],
                                                                     ascending=[False, True, True, True])
    accepted: dict[str, tuple] = {}
    rejected = defaultdict(int)
    for a, b, meth, conf, d in zip(e["a"], e["b"], e["method"], e["conf"], e["dist_km"]):
        if a not in parent or b not in parent:
            continue
        ra, rb = find(a), find(b)
        if ra == rb:
            continue
        ca, cb = info[ra], info[rb]
        if ca.n_hrsa + cb.n_hrsa > 1:
            rejected["two_hrsa_sites"] += 1
            continue
        if ca.n_nih + cb.n_nih > 1:
            rejected["two_nih_orgs"] += 1
            continue
        if len(ca.states | cb.states) > 1:
            rejected["cross_state"] += 1
            continue
        if ca.coords and cb.coords:
            A = np.array(ca.coords)
            B = np.array(cb.coords)
            dd = haversine_km(A[:, None, 0], A[:, None, 1], B[None, :, 0], B[None, :, 1])
            if float(np.max(dd)) > MAX_SPAN_KM:
                rejected["span_over_25km"] += 1
                continue
        # union (attach smaller to larger)
        if len(ca.coords) < len(cb.coords):
            ra, rb, ca, cb = rb, ra, cb, ca
        parent[rb] = ra
        ca.n_hrsa += cb.n_hrsa
        ca.n_nih += cb.n_nih
        ca.coords = list(dict.fromkeys(ca.coords + cb.coords))
        ca.states |= cb.states
        # record the linking edge for both endpoints (first accepted edge wins)
        for x, y in ((a, b), (b, a)):
            accepted.setdefault(x, (meth, conf, y, d))
    root = {r: find(r) for r in parent}
    _ = src
    return root, accepted, dict(rejected)


SOURCE_PRIORITY = {"hrsa": 0, "nppes": 1, "nih": 2, "ctgov": 3}
PRECISION_RANK = {"point": 0, "zcta_centroid": 1, "city_centroid": 2, "none": 3}


def _facility_id(record_id: str) -> str:
    src, native = record_id.split(":", 1)
    if src == "ctgov_site":
        import hashlib
        return "ctgov-" + hashlib.sha1(native.encode()).hexdigest()[:12]
    if src == "nih_org":
        return "nihorg-" + re.sub(r"[^A-Za-z0-9]+", "-", native).strip("-")[:40]
    if src == "hrsa_site":
        return "hrsa-" + native
    return f"{src}-{native}"


def assign_facilities(nodes: pd.DataFrame, root: dict) -> pd.DataFrame:
    n = nodes[["record_id", "source", "native_id", "precision"]].copy()
    n["root"] = n["record_id"].map(root)
    n["prio"] = n["source"].map(SOURCE_PRIORITY)
    n["npi_num"] = np.where(n["source"] == "nppes", pd.to_numeric(n["native_id"], errors="coerce"), np.nan)
    n = n.sort_values(["root", "prio", "npi_num", "record_id"])
    anchor = n.groupby("root")["record_id"].first()
    n["anchor_record_id"] = n["root"].map(anchor)
    n["facility_id"] = n["anchor_record_id"].map(_facility_id)
    return n.set_index("record_id")[["facility_id", "anchor_record_id"]]


# --------------------------------------------------------------------------------------------------------------
# trial and NIH link tables
# --------------------------------------------------------------------------------------------------------------

OPEN_STATUSES = {"RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"}


def trial_measurement_classes() -> pd.DataFrame:
    """nct_id x measurement class ids whose configs/measurements.yaml patterns match the registered text
    (titles, brief summary, keywords, intervention name/description/other names, outcome measure/description;
    eligibility criteria are not searched) - the same fields and regexes as clinicaltrials.trials_for_condition.
    A match means the trial REGISTERED the measurement, not that it works."""
    from ..ingestion.clinicaltrials import _text_hits, measurement_regex

    ct = read_table("clinical_trials", columns=["nct_id", "brief_title", "official_title", "brief_summary", "keywords"])
    iv = read_table("trial_interventions", columns=["nct_id", "intervention_name", "intervention_description",
                                                    "other_names"])
    oc = read_table("trial_outcomes", columns=["nct_id", "measure", "description"])
    fields = [(ct, c) for c in ("brief_title", "official_title", "brief_summary", "keywords")]
    fields += [(iv, c) for c in ("intervention_name", "intervention_description", "other_names")]
    fields += [(oc, c) for c in ("measure", "description")]
    rows = []
    for m in load_config("measurements")["measurement_classes"]:
        inc, exc, _ = measurement_regex(m["id"])
        hit_ids: set[str] = set()
        for frame, col in fields:
            h = _text_hits(frame[col].astype(object), inc, exc)
            hit_ids.update(frame.loc[h.notna().values, "nct_id"])
        rows.extend((n, m["id"]) for n in hit_ids)
    return pd.DataFrame(rows, columns=["nct_id", "measurement_class"])


def build_facility_trials(site_rows: pd.DataFrame, fac_of: pd.Series) -> pd.DataFrame:
    ct = read_table("clinical_trials", columns=["nct_id", "overall_status", "study_type", "start_year",
                                                "last_update_post_date", "primary_completion_date",
                                                "intervention_types", "brief_title", "lead_sponsor_class",
                                                "source_version", "retrieved_at"])
    tc = read_table("trial_conditions", columns=["nct_id", "condition_id", "condition_literal_match"])
    lit = tc[tc["condition_literal_match"].astype(bool)].groupby("nct_id")["condition_id"].agg(
        lambda x: _join(sorted(set(x))))
    allc = tc.groupby("nct_id")["condition_id"].agg(lambda x: _join(sorted(set(x))))
    meas = trial_measurement_classes()
    mc = meas.groupby("nct_id")["measurement_class"].agg(lambda x: _join(sorted(set(x))))
    s = site_rows[["nct_id", "site_key", "facility", "site_status"]].copy()
    s["facility_id"] = ("ctgov_site:" + s["site_key"]).map(fac_of)
    s["site_recruiting"] = s["site_status"].fillna("") == "RECRUITING"
    g = s.groupby(["facility_id", "nct_id"], sort=True).agg(
        site_keys=("site_key", lambda x: _join(sorted(set(x)))),
        facility_names_registered=("facility", lambda x: _join(sorted(set(x)))),
        site_status=("site_status", lambda x: _join(sorted(set(x.dropna())))),
        site_recruiting=("site_recruiting", "any")).reset_index()
    out = g.merge(ct, on="nct_id", how="left")
    out["condition_ids_literal"] = out["nct_id"].map(lit).fillna("")
    out["condition_ids_all"] = out["nct_id"].map(allc).fillna("")
    out["measurement_classes"] = out["nct_id"].map(mc).fillna("")
    it = out["intervention_types"].fillna("")
    out["is_device_or_diagnostic"] = it.str.contains(r"\bDEVICE\b|\bDIAGNOSTIC_TEST\b")
    out["overall_recruiting"] = out["overall_status"] == "RECRUITING"
    out["is_open_status"] = out["overall_status"].isin(OPEN_STATUSES)
    out["object_id"] = "trial:" + out["nct_id"]
    out["facility_object_id"] = "facility:" + out["facility_id"]
    out = out.rename(columns={"source_version": "_sv", "retrieved_at": "_ra"})
    out["data_layer"] = "facility"
    out["source_name"] = "ClinicalTrials.gov API v2 (U.S. NLM) via facility entity resolution"
    out["source_record_id"] = out["nct_id"] + "@" + out["facility_id"]
    out["source_version"] = out["_sv"]
    out["retrieved_at"] = out["_ra"]
    out["source_geographic_resolution"] = "none"
    out["evidence_type"] = "trial_registry_record"
    out["evidence_level"] = np.where(out["condition_ids_literal"] != "", "condition_literal_match",
                                     "mesh_or_synonym_expansion_only")
    out["provenance_notes"] = ("U.S. identifiable site of a target-condition trial, attached to a resolved facility. "
                               "Registration is a research-activity signal, not evidence a measurement or "
                               "intervention works. measurement_classes = configs/measurements.yaml patterns found in "
                               "registered titles/summary/keywords/interventions/outcomes.")
    return out.drop(columns=["_sv", "_ra"])


def build_facility_nih(nih_rows: pd.DataFrame, fac_of: pd.Series) -> pd.DataFrame:
    pc = read_table("nih_project_conditions", columns=["appl_id", "condition_id", "match_tier", "likely_false_positive"])
    prec = pc[(pc["match_tier"] == "title_abstract") & ~pc["likely_false_positive"].fillna(False).astype(bool)]
    cp = prec.groupby("appl_id")["condition_id"].agg(lambda x: _join(sorted(set(x))))
    ca = pc.groupby("appl_id")["condition_id"].agg(lambda x: _join(sorted(set(x))))
    n = nih_rows.copy()
    n["facility_id"] = ("nih_org:" + n["org_key"].astype(str)).map(fac_of)
    keep = ["facility_id", "appl_id", "org_key", "org_name", "project_num", "core_project_num", "fiscal_year",
            "is_subproject", "is_active", "award_amount", "activity_code", "agency_code", "pi_profile_ids",
            "contact_pi_name", "n_pis", "project_title", "augmented_terms_hit", "project_start_date",
            "project_end_date", "source_version", "retrieved_at"]
    out = n[keep].copy()
    out["condition_ids_precision"] = out["appl_id"].map(cp).fillna("")
    out["condition_ids_all"] = out["appl_id"].map(ca).fillna("")
    out["object_id"] = "nih:" + out["appl_id"].astype(str)
    out["facility_object_id"] = "facility:" + out["facility_id"]
    out = out.rename(columns={"source_version": "_sv", "retrieved_at": "_ra"})
    out["data_layer"] = "facility"
    out["source_name"] = "NIH RePORTER Project API v2 via facility entity resolution"
    out["source_record_id"] = out["appl_id"].astype(str)
    out["source_version"] = out["_sv"]
    out["retrieved_at"] = out["_ra"]
    out["source_geographic_resolution"] = "none"
    out["evidence_type"] = "grant_record"
    out["evidence_level"] = np.where(out["condition_ids_precision"] != "", "title_abstract_match",
                                     "terms_only_or_false_positive")
    out["provenance_notes"] = ("RePORTER fiscal-year award record of the facility's awardee organisation; "
                               "research-capability signal, not burden. condition_ids_precision = title/abstract "
                               "matches without likely false positives.")
    return out.drop(columns=["_sv", "_ra"])


# --------------------------------------------------------------------------------------------------------------
# facility table
# --------------------------------------------------------------------------------------------------------------

def _county_context() -> pd.DataFrame:
    h = read_table("geo_context__hpsa", columns=["county_fips", "pc_coverage", "mua_coverage"])
    r = read_table("geo_context__rucc", columns=["geo_id", "rucc_2023", "is_metro"]).rename(columns={"geo_id": "county_fips"})
    c = h.merge(r, on="county_fips", how="outer")
    c["county_pc_hpsa_area"] = c["pc_coverage"].isin(["whole_county_geographic", "whole_county_population",
                                                       "partial_county"])
    c["county_nonmetro"] = c["rucc_2023"].notna() & (c["rucc_2023"] >= 4)
    return c.rename(columns={"pc_coverage": "county_pc_hpsa_coverage", "mua_coverage": "county_mua_coverage"})


def _colocated_specialists() -> pd.Series:
    """Individual NPIs (entity type 1) per ZIP5 and SPEC specialty group, as 'group:n|...' strings.
    Co-location in a ZIP is NOT affiliation with any facility."""
    psg = read_table("provider_specialty_groups", columns=["npi", "specialty_group", "entity_type", "spec_listed"])
    psg = psg[(psg["entity_type"] == 1) & psg["spec_listed"].astype(bool)]
    z = read_table("providers", columns=["npi", "zip5"])
    psg = psg.merge(z, on="npi", how="inner")
    cnt = psg.groupby(["zip5", "specialty_group"]).size().reset_index(name="n")
    cnt = cnt.sort_values(["zip5", "n", "specialty_group"], ascending=[True, False, True])
    cnt["s"] = cnt["specialty_group"].astype(str) + ":" + cnt["n"].astype(str)
    return cnt.groupby("zip5")["s"].agg(SEP.join)


def build_facilities(nodes: pd.DataFrame, fac: pd.DataFrame, accepted: dict, extras: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (facility_source_links, facilities)."""
    n = nodes.join(fac, on="record_id")
    acc = pd.DataFrame.from_dict(accepted, orient="index", columns=["match_method", "match_confidence",
                                                                   "linked_to_record_id", "link_distance_km"])
    n = n.join(acc, on="record_id")
    is_anchor = n["record_id"] == n["anchor_record_id"]
    n.loc[n["match_method"].isna(), "match_method"] = np.where(is_anchor[n["match_method"].isna()], "anchor_single",
                                                               "unlinked")
    sizes = n.groupby("facility_id")["record_id"].transform("size")
    n.loc[sizes == 1, "match_method"] = "single_record"
    n.loc[(sizes > 1) & is_anchor & n["linked_to_record_id"].isna(), "match_method"] = "anchor"
    n["match_confidence"] = np.where(n["match_method"].isin(["single_record", "anchor"]), 1.0, n["match_confidence"])

    # per-source extras
    hrsa, nppes, ctgov, nih = extras["hrsa"], extras["nppes"], extras["ctgov"], extras["nih"]
    n = n.merge(nppes[["record_id", "specialty_groups", "primary_specialty_groups", "group_types",
                       "primary_taxonomy_display_name", "is_organization_subpart", "address_line1",
                       "last_update_date"]], on="record_id", how="left")
    n = n.merge(hrsa[["record_id", "health_center_name", "is_look_alike", "is_service_delivery_site", "site_type",
                      "location_setting", "location_type", "address", "site_npi"]], on="record_id", how="left")
    n = n.merge(nih[["record_id", "org_type_name", "n_appl_ids"]], on="record_id", how="left")
    n = n.merge(ctgov[["record_id", "n_site_rows", "geocode_method"]].rename(
        columns={"geocode_method": "ctgov_geocode_method"}), on="record_id", how="left")
    n["address"] = n["address"].where(n["address"].notna(), n["address_line1"])

    links = n.copy()
    links["object_id"] = np.select(
        [links["source"] == "nppes", links["source"] == "hrsa"],
        ["npi:" + links["native_id"].astype(str), "hrsa_site:" + links["native_id"].astype(str)], default="")
    links["facility_object_id"] = "facility:" + links["facility_id"]

    # facility-level aggregation
    n["prec_rank"] = n["precision"].map(PRECISION_RANK)
    n["prio"] = n["source"].map(SOURCE_PRIORITY)
    loc = n.sort_values(["facility_id", "prec_rank", "prio", "record_id"]).groupby("facility_id").first()
    disp = n.assign(dprio=n["source"].map({"hrsa": 0, "ctgov": 1, "nih": 2, "nppes": 3}),
                    rows=-n["n_site_rows"].fillna(0)).sort_values(
        ["facility_id", "dprio", "rows", "record_id"]).groupby("facility_id")["name"].first()

    def agg_join(col, sub=None, prefix=""):
        d = n if sub is None else n[n["source"] == sub]
        return d.groupby("facility_id")[col].agg(lambda x: _join(sorted(prefix + str(v) for v in x.dropna())))

    groups = n.groupby("facility_id")["specialty_groups"].agg(
        lambda x: _join(sorted({g for v in x.dropna() for g in str(v).split("|") if g})))
    f = pd.DataFrame(index=loc.index)
    f["facility_name"] = disp
    f["member_names"] = n.groupby("facility_id")["name"].agg(lambda x: _join(list(dict.fromkeys(x.dropna()))[:8]))
    f["facility_kinds"] = n.groupby("facility_id")["kind"].agg(lambda x: _join(sorted(set(x))))
    f["sources"] = n.groupby("facility_id")["source"].agg(lambda x: _join(sorted(set(x), key=SOURCE_PRIORITY.get)))
    f["n_source_records"] = n.groupby("facility_id").size()
    f["n_npis"] = n[n["source"] == "nppes"].groupby("facility_id").size().reindex(f.index).fillna(0).astype(int)
    f["npi_object_ids"] = agg_join("native_id", "nppes", "npi:").reindex(f.index).fillna("")
    f["hrsa_site_object_ids"] = agg_join("native_id", "hrsa", "hrsa_site:").reindex(f.index).fillna("")
    f["ctgov_site_keys"] = agg_join("native_id", "ctgov").reindex(f.index).fillna("")
    f["nih_org_keys"] = agg_join("native_id", "nih").reindex(f.index).fillna("")
    f["specialty_groups"] = groups
    hr = n[n["source"] == "hrsa"].groupby("facility_id")
    f["is_hrsa_site"] = f.index.isin(hr.groups.keys())
    f["hrsa_is_look_alike"] = hr["is_look_alike"].any().reindex(f.index).fillna(False).astype(bool)
    f["hrsa_service_delivery"] = hr["is_service_delivery_site"].any().reindex(f.index).fillna(False).astype(bool)
    f["hrsa_site_type"] = hr["site_type"].first().reindex(f.index)
    f["hrsa_location_setting"] = hr["location_setting"].first().reindex(f.index)
    f["hrsa_location_type"] = hr["location_type"].first().reindex(f.index)
    f["hrsa_health_center_name"] = hr["health_center_name"].first().reindex(f.index)
    for c in ("address", "city", "state", "zip5", "county_fips", "state_fips", "lat", "lon", "precision"):
        f[c] = loc[c]
    f["location_source_record_id"] = loc.index.map(
        n.sort_values(["facility_id", "prec_rank", "prio", "record_id"]).groupby("facility_id")["record_id"].first())
    f = f.rename(columns={"precision": "geocode_precision"})
    methods = n[~n["match_method"].isin(["single_record", "anchor"])].groupby("facility_id")["match_method"].agg(
        lambda x: _join(sorted(set(x))))
    f["match_methods"] = methods.reindex(f.index).fillna("")
    minconf = n[~n["match_method"].isin(["single_record", "anchor"])].groupby("facility_id")["match_confidence"].min()
    f["match_confidence_min"] = minconf.reindex(f.index)
    f["resolution_tier"] = [("single_record" if k == 1 else confidence_tier(c))
                            for k, c in zip(f["n_source_records"], f["match_confidence_min"])]
    gl = groups.reindex(f.index).fillna("")
    f["is_fqhc_or_lookalike"] = f["is_hrsa_site"] | gl.str.contains(r"(?:^|\|)fqhc(?:\||$)")
    f["fqhc_evidence"] = np.where(f["is_hrsa_site"], "hrsa_site_roster",
                                  np.where(f["is_fqhc_or_lookalike"], "nppes_self_reported_taxonomy", ""))
    f["is_community_type_facility"] = gl.str.contains(
        r"(?:^|\|)(?:community_health_clinic|rural_health_clinic|critical_access_hospital|public_health_clinic)(?:\||$)")
    kinds = f["facility_kinds"]
    f["is_group_practice_only"] = kinds == "nppes_group_practice"
    f["primary_kind"] = np.select(
        [f["is_hrsa_site"], kinds.str.contains("nppes_facility"), kinds.str.contains("nppes_group_practice"),
         kinds.str.contains("nih_awardee_org"), kinds.str.contains("ctgov_trial_site"),
         kinds.str.contains("nppes_research")],
        ["hrsa_health_center_site", "nppes_facility", "nppes_group_practice", "nih_awardee_org", "ctgov_trial_site",
         "nppes_research"], default="other")
    # a facility whose only clinical record is an HRSA administrative site is not a clinic candidate
    f["is_clinic_candidate"] = (f["geocode_precision"] != "none") & ~(
        f["is_hrsa_site"] & ~f["hrsa_service_delivery"] & (f["sources"] == "hrsa"))
    ctx = _county_context().set_index("county_fips")
    f = f.join(ctx[["county_pc_hpsa_coverage", "county_pc_hpsa_area", "county_mua_coverage", "rucc_2023",
                    "county_nonmetro"]], on="county_fips")
    f["county_pc_hpsa_area"] = f["county_pc_hpsa_area"].fillna(False).astype(bool)
    f["county_nonmetro"] = f["county_nonmetro"].fillna(False).astype(bool)
    f = f.reset_index()
    f["object_id"] = "facility:" + f["facility_id"]
    return links, f


def _provenance_facility(df: pd.DataFrame, source_names: dict, source_versions: dict, retrieved: dict,
                         evidence_type: str, notes: str) -> pd.DataFrame:
    df = df.copy()
    srcs = df["sources"].str.split(SEP)
    df["data_layer"] = "facility"
    df["source_name"] = srcs.map(lambda s: "; ".join(source_names[x] for x in s))
    df["source_record_id"] = df["facility_id"]
    df["source_version"] = srcs.map(lambda s: "; ".join(source_versions[x] for x in s))
    df["retrieved_at"] = srcs.map(lambda s: max(retrieved[x] for x in s))
    df["source_geographic_resolution"] = df["geocode_precision"]
    df["evidence_type"] = evidence_type
    df["evidence_level"] = "entity_resolution:" + df["resolution_tier"]
    df["provenance_notes"] = notes
    return df


# --------------------------------------------------------------------------------------------------------------
# registries
# --------------------------------------------------------------------------------------------------------------

def build_research_registry(f: pd.DataFrame, ft: pd.DataFrame, fn: pd.DataFrame) -> pd.DataFrame:
    from ..ingestion.nih_reporter import obligation_total

    rows = []
    ft_g = dict(tuple(ft.groupby("facility_id")))
    fn_g = dict(tuple(fn.groupby("facility_id")))
    ids = sorted(set(ft_g) | set(fn_g))
    for fid in ids:
        t = ft_g.get(fid)
        p = fn_g.get(fid)
        r = {"facility_id": fid}
        if t is not None:
            lit = t["condition_ids_literal"].str.split(SEP).explode()
            lit = lit[lit.notna() & (lit != "")]
            allc = t["condition_ids_all"].str.split(SEP).explode()
            allc = allc[allc.notna() & (allc != "")]
            meas = t["measurement_classes"].str.split(SEP).explode()
            meas = meas[meas.notna() & (meas != "")]
            dev_meas = t[t["is_device_or_diagnostic"] & (t["measurement_classes"] != "")]
            r.update({
                "n_trials": int(t["nct_id"].nunique()),
                "n_trials_literal": int(t.loc[t["condition_ids_literal"] != "", "nct_id"].nunique()),
                "trials_by_condition_literal": _join(f"{k}:{v}" for k, v in lit.value_counts().sort_index().items()),
                "trials_by_condition_all": _join(f"{k}:{v}" for k, v in allc.value_counts().sort_index().items()),
                "condition_ids_literal": _join(sorted(set(lit))),
                "n_trials_site_recruiting": int(t.loc[t["site_recruiting"], "nct_id"].nunique()),
                "n_trials_overall_recruiting": int(t.loc[t["overall_recruiting"], "nct_id"].nunique()),
                "n_trials_open_status": int(t.loc[t["is_open_status"], "nct_id"].nunique()),
                "n_interventional": int((t["study_type"] == "INTERVENTIONAL").sum()),
                "n_observational": int((t["study_type"] == "OBSERVATIONAL").sum()),
                "first_start_year": t["start_year"].min(),
                "latest_start_year": t["start_year"].max(),
                "latest_update_post_date": t["last_update_post_date"].dropna().max() if t[
                    "last_update_post_date"].notna().any() else None,
                "n_device_or_diagnostic_trials": int(t.loc[t["is_device_or_diagnostic"], "nct_id"].nunique()),
                "n_trials_with_measurement": int(t.loc[t["measurement_classes"] != "", "nct_id"].nunique()),
                "measurement_classes_studied": _join(f"{k}:{v}" for k, v in meas.value_counts().sort_index().items()),
                "n_device_diag_trials_with_measurement": int(dev_meas["nct_id"].nunique()),
                "trial_object_ids": _join(sorted("trial:" + t["nct_id"].unique())),
            })
        if p is not None:
            pr = p["condition_ids_precision"].str.split(SEP).explode()
            pr = pr[pr.notna() & (pr != "")]
            prec = p[p["condition_ids_precision"] != ""]
            pis = {i.strip() for s in p["pi_profile_ids"].dropna() for i in str(s).split(";") if i.strip()}
            act = p[p["is_active"].fillna(False).astype(bool)]
            r.update({
                "n_nih_projects": int(p["appl_id"].nunique()),
                "n_nih_projects_precision": int(prec["appl_id"].nunique()),
                "nih_projects_by_condition_precision": _join(
                    f"{k}:{v}" for k, v in pr.value_counts().sort_index().items()),
                "nih_condition_ids_precision": _join(sorted(set(pr))),
                "n_nih_core_projects": int(p["core_project_num"].nunique()),
                "n_nih_core_projects_precision": int(prec["core_project_num"].nunique()),
                "n_nih_active_core_projects": int(act["core_project_num"].nunique()),
                "nih_award_obligations_usd": (obligation_total(p) if p["award_amount"].notna().any() else np.nan),
                "nih_award_obligations_usd_precision": (obligation_total(prec) if prec["award_amount"].notna().any()
                                                        else np.nan),
                "n_nih_award_amount_missing": int(p["award_amount"].isna().sum()),
                "n_nih_pis": len(pis),
                "nih_first_fiscal_year": int(p["fiscal_year"].min()),
                "nih_last_fiscal_year": int(p["fiscal_year"].max()),
                "nih_activity_codes": _join(sorted(p["activity_code"].dropna().unique())),
                "nih_measurement_terms": _join(sorted({x.strip() for s in p["augmented_terms_hit"].dropna()
                                                       for x in re.split(r"[;|]", str(s)) if x.strip()})),
                "nih_object_ids": _join(sorted("nih:" + p["appl_id"].astype(str).unique())),
            })
        rows.append(r)
    rr = pd.DataFrame(rows)
    int_cols = ["n_trials", "n_trials_literal", "n_trials_site_recruiting", "n_trials_overall_recruiting",
                "n_trials_open_status", "n_interventional", "n_observational", "n_device_or_diagnostic_trials",
                "n_trials_with_measurement", "n_device_diag_trials_with_measurement", "n_nih_projects",
                "n_nih_projects_precision", "n_nih_core_projects", "n_nih_core_projects_precision",
                "n_nih_active_core_projects", "n_nih_award_amount_missing", "n_nih_pis"]
    for c in int_cols:
        rr[c] = rr[c].fillna(0).astype(int)
    for c in ("first_start_year", "latest_start_year", "nih_first_fiscal_year", "nih_last_fiscal_year"):
        rr[c] = pd.to_numeric(rr[c], errors="coerce").astype("Int64")
    for c in ("trials_by_condition_literal", "trials_by_condition_all", "condition_ids_literal",
              "measurement_classes_studied", "trial_object_ids", "nih_projects_by_condition_precision",
              "nih_condition_ids_precision", "nih_activity_codes", "nih_measurement_terms", "nih_object_ids"):
        rr[c] = rr[c].fillna("")
    keep = ["facility_id", "object_id", "facility_name", "primary_kind", "facility_kinds", "sources", "address",
            "city", "state", "zip5", "county_fips", "state_fips", "lat", "lon", "geocode_precision",
            "resolution_tier", "match_methods", "match_confidence_min", "npi_object_ids", "hrsa_site_object_ids",
            "ctgov_site_keys", "nih_org_keys"]
    rr = f[keep].merge(rr, on="facility_id", how="inner")
    # Added at review: a facility enters this table with ANY ingested trial or NIH record, including trials that match
    # a target condition only through MeSH/synonym expansion and NIH links that are terms-only or likely false
    # positives. This flag marks rows with precision-view evidence (a literal-match trial or a title/abstract NIH match).
    # The aggregate recruiting / recency / device / measurement counts span all of the facility's ingested trials.
    rr["has_precision_view_record"] = (rr["n_trials_literal"] > 0) | (rr["n_nih_projects_precision"] > 0)
    return rr


def build_clinic_registry(f: pd.DataFrame, rr: pd.DataFrame, coloc: pd.Series) -> pd.DataFrame:
    c = f[f["sources"].str.contains("nppes|hrsa")].copy()
    c["colocated_individual_npis_same_zip"] = c["zip5"].map(coloc).fillna("")
    r = rr.set_index("facility_id")
    c["has_research_record"] = c["facility_id"].isin(r.index)
    c["n_trials"] = c["facility_id"].map(r["n_trials"]).fillna(0).astype(int)
    c["n_nih_projects"] = c["facility_id"].map(r["n_nih_projects"]).fillna(0).astype(int)
    return c


# --------------------------------------------------------------------------------------------------------------
# validation (V1 ground truth, V2 location shuffle)
# --------------------------------------------------------------------------------------------------------------

def validate_hrsa_npi(nodes: pd.DataFrame, hrsa: pd.DataFrame, edges_no_npi: pd.DataFrame) -> dict:
    """V1: HRSA site_npi as ground truth for name/location linking (shared_npi edges disabled)."""
    grp = edges_no_npi[edges_no_npi["method"] == "nppes_same_name_same_zip"]
    to_anchor = dict(zip(grp["a"], grp["b"]))
    root, _, _ = cluster(nodes, edges_no_npi)
    npis = set(nodes.loc[nodes["source"] == "nppes", "native_id"])
    h = hrsa[hrsa["site_npi"].notna() & hrsa["site_npi"].astype(str).isin(npis)]
    norm_of = dict(zip(nodes["record_id"], nodes["norm"]))
    members = defaultdict(list)
    for rid, rt in root.items():
        if rid.startswith("npi:"):
            members[rt].append(rid)
    strict = lenient = linked = 0
    for rid, npi in zip(h["record_id"], h["site_npi"].astype(str)):
        true_rec = "npi:" + npi
        true_anchor = to_anchor.get(true_rec, true_rec)
        mem = members.get(root[rid], [])
        if not mem:
            continue
        linked += 1
        anchors = {to_anchor.get(m, m) for m in mem}
        if true_anchor in anchors:
            strict += 1
            lenient += 1
        elif norm_of.get(true_rec) and any(norm_of.get(m) == norm_of.get(true_rec) for m in mem):
            lenient += 1
    n = len(h)
    return {"validation": "V1_hrsa_site_npi_ground_truth", "n_hrsa_sites_with_npi_in_nppes_orgs": n,
            "n_linked_to_any_nppes_by_name_location": linked,
            "n_linked_to_true_npi_group": strict, "n_linked_to_same_name_as_true_npi": lenient,
            "precision_strict": strict / linked if linked else np.nan,
            "precision_lenient": lenient / linked if linked else np.nan,
            "recall_strict": strict / n if n else np.nan}


def location_shuffle_false_links(nodes: pd.DataFrame, n_draws: int = 3) -> pd.DataFrame:
    """V2: permute location fields of ClinicalTrials.gov and NIH records, count surviving cross-source links."""
    grp = nppes_group_edges(nodes)
    reps = _nppes_group_representatives(nodes, grp)
    loc_cols = ["city", "state", "zip5", "lat", "lon", "city_n", "county_fips", "state_fips"]
    source_of = dict(zip(nodes["record_id"], nodes["source"]))

    qs = frozenset({"ctgov", "nih"})

    def links(frame):
        e = pd.concat([exact_zip_edges(frame), fuzzy_edges(frame, "zip", qs),
                       mutual_best(pd.concat([exact_city_edges(frame), fuzzy_edges(frame, "city", qs)],
                                             ignore_index=True), source_of)], ignore_index=True)
        sa, sb = e["a"].map(source_of), e["b"].map(source_of)
        e = e[(sa.isin(["ctgov", "nih"]) | sb.isin(["ctgov", "nih"])) & (sa != sb)]
        return e.groupby("method").size()

    rows = [links(reps).rename("observed")]
    rng = np.random.default_rng(SEED)
    q = reps["source"].isin(["ctgov", "nih"])
    for d in range(n_draws):
        sh = reps.copy()
        idx = np.flatnonzero(q.values)
        perm = rng.permutation(idx)
        for c in loc_cols:
            col = sh[c].to_numpy(copy=True)
            col[idx] = reps[c].to_numpy()[perm]
            sh[c] = col
        rows.append(links(sh).rename(f"shuffle_{d + 1}"))
    out = pd.concat(rows, axis=1).fillna(0).astype(int)
    sh_cols = [c for c in out.columns if c.startswith("shuffle_")]
    out["shuffle_mean"] = out[sh_cols].mean(axis=1)
    out["estimated_false_link_share"] = out["shuffle_mean"] / out["observed"].replace(0, np.nan)
    out.loc["ALL"] = out.sum(numeric_only=True)
    out.loc["ALL", "estimated_false_link_share"] = out.loc["ALL", "shuffle_mean"] / out.loc["ALL", "observed"]
    return out.reset_index().rename(columns={"index": "method"})


# --------------------------------------------------------------------------------------------------------------
# run
# --------------------------------------------------------------------------------------------------------------

def run(skip_validation: bool = False) -> dict:
    t0 = utc_now_iso()
    TABLES.mkdir(parents=True, exist_ok=True)
    nppes = load_nppes_records()
    hrsa = load_hrsa_records()
    site_rows, ct_stats = load_ctgov_site_rows()
    ctgov = load_ctgov_records(site_rows)
    nih, nih_rows = load_nih_records()
    nodes = _prepare_nodes([nppes, hrsa, ctgov, nih])
    print(f"records: {nodes['source'].value_counts().to_dict()}", flush=True)

    edges = candidate_edges(nodes, hrsa)
    print(f"candidate edges: {edges['method'].value_counts().to_dict()}", flush=True)
    root, accepted, rejected = cluster(nodes, edges)
    fac = assign_facilities(nodes, root)
    extras = {"nppes": nppes, "hrsa": hrsa, "ctgov": ctgov, "nih": nih}
    links, f = build_facilities(nodes, fac, accepted, extras)
    fac_of = fac["facility_id"]

    ft = build_facility_trials(site_rows, fac_of)
    fn = build_facility_nih(nih_rows, fac_of)
    rr = build_research_registry(f, ft, fn)
    coloc = _colocated_specialists()
    cr = build_clinic_registry(f, rr, coloc)

    names = {"nppes": str(nppes["source_name"].iloc[0]), "hrsa": str(hrsa["source_name"].iloc[0]),
             "ctgov": str(ctgov["source_name"].iloc[0]), "nih": str(nih["source_name"].iloc[0])}
    versions = {"nppes": str(nppes["source_version"].iloc[0]), "hrsa": str(hrsa["source_version"].iloc[0]),
                "ctgov": str(ctgov["source_version"].iloc[0]), "nih": str(nih["source_version"].iloc[0])}
    retrieved = {"nppes": str(nppes["retrieved_at"].max()), "hrsa": str(hrsa["retrieved_at"].max()),
                 "ctgov": str(ctgov["retrieved_at"].max()), "nih": str(nih["retrieved_at"].max())}

    f_out = _provenance_facility(
        f, names, versions, retrieved, "facility_registry",
        "Resolved facility (conservative name + ZIP/city+distance entity resolution across NPPES organisation NPIs, "
        "HRSA sites, U.S. ClinicalTrials.gov facilities and NIH RePORTER organisations). NPPES organisations under "
        "physician taxonomies are group practices; taxonomy is self-reported and does not mean the facility treats "
        "any condition. geocode_precision gives the location precision (ClinicalTrials.gov geoPoints are city "
        "centroids).")
    write_table(f_out, "facilities", producer=PRODUCER,
                description="Unified resolved facilities (NPPES orgs, HRSA sites, ClinicalTrials.gov US facilities, "
                            "NIH RePORTER orgs) with match method/confidence and every source id")

    links_out = links.copy()
    links_out["sources"] = links_out["source"]
    links_out["resolution_tier"] = links_out["match_confidence"].map(
        lambda c: "single_record" if pd.isna(c) else confidence_tier(c))
    links_out["geocode_precision"] = links_out["precision"]
    links_out = _provenance_facility(links_out, names, versions, retrieved, "facility_registry",
                                     "One source record and the resolved facility it was assigned to; match_method/"
                                     "match_confidence describe the link that attached it (anchor/single_record = "
                                     "the facility's own record).")
    links_out["source_record_id"] = links_out["record_id"]
    links_out["evidence_level"] = "match:" + links_out["match_method"].astype(str)
    keep_links = ["record_id", "object_id", "facility_id", "facility_object_id", "source", "native_id", "kind", "name",
                  "alt_name", "norm", "address", "city", "state", "zip5", "county_fips", "lat", "lon",
                  "geocode_precision", "match_method", "match_confidence", "linked_to_record_id", "link_distance_km",
                  "specialty_groups", "primary_taxonomy_display_name", "is_organization_subpart", "site_npi",
                  "is_service_delivery_site", "org_type_name", "n_appl_ids", "n_site_rows"] + PROVENANCE_COLUMNS
    links_out = links_out[keep_links]
    for c in ("alt_name", "site_npi", "is_organization_subpart", "org_type_name"):
        links_out[c] = links_out[c].astype("string")
    links_out["is_service_delivery_site"] = links_out["is_service_delivery_site"].astype("boolean")
    write_table(links_out, "facility_source_links", producer=PRODUCER,
                description="Source record -> resolved facility with match_method, match_confidence and linked record")
    write_table(ft, "facility_trials", producer=PRODUCER,
                description="Resolved facility x U.S. target-condition trial (literal and MeSH-expanded condition ids, "
                            "status, measurement classes registered)")
    write_table(fn, "facility_nih_projects", producer=PRODUCER,
                description="Resolved facility x NIH RePORTER appl_id (precision-view and all-tier condition ids)")

    rr_out = rr.merge(f_out[["facility_id", "sources"]], on="facility_id", how="left", suffixes=("", "_f"))
    rr_out = _provenance_facility(
        rr_out.drop(columns=["sources_f"]), names, versions, retrieved, "facility_registry",
        "Research readiness per resolved facility: ClinicalTrials.gov site history for target-condition trials "
        "(literal-match and MeSH-expanded views separate) and NIH RePORTER awards of the facility's awardee "
        "organisation (obligations counted once per appl_id without parent/subproject double counting). "
        "Research readiness is NOT patient burden and no burden variable is included. A trial registering a "
        "measurement is not evidence the measurement works. Rows with has_precision_view_record=False have only "
        "MeSH/synonym-expanded trials or terms-only/likely-false-positive NIH links; aggregate recruiting/recency/"
        "device/measurement counts span all ingested trials (per-condition literal counts are in "
        "trials_by_condition_literal).")
    write_table(rr_out, "research_site_registry", producer=PRODUCER,
                description="Research-oriented registry: per resolved facility trial history by condition, "
                            "recruiting, recency, device/diagnostic and measurement trials, NIH projects, "
                            "active projects, obligations, PIs, topics")
    cr_out = _provenance_facility(
        cr, names, versions, retrieved, "facility_registry",
        "Implementation-oriented registry: resolved facilities with an NPPES organisation NPI or an HRSA site. "
        "specialty_groups are self-reported NPPES taxonomies (HRSA roster adds fqhc) and do NOT mean the facility "
        "evaluates or treats any condition. County HPSA/RUCC flags are county context of the location, not "
        "facility attributes. colocated_individual_npis_same_zip counts individual NPIs sharing the ZIP "
        "(co-location, not affiliation).")
    write_table(cr_out, "clinic_registry", producer=PRODUCER,
                description="Implementation-oriented clinic registry (NPPES organisations and HRSA sites, resolved)")

    # ---------------- summaries / validation tables
    summ = []
    for s, g in nodes.groupby("source"):
        summ.append({"table": "source_records", "key": s, "n": len(g),
                     "n_geocoded": int((g["precision"] != "none").sum())})
    for k, v in ct_stats.items():
        if not isinstance(v, list):
            summ.append({"table": "ctgov_site_filter", "key": k, "n": v, "n_geocoded": np.nan})
    for m, v in edges["method"].value_counts().items():
        summ.append({"table": "candidate_edges", "key": m, "n": int(v), "n_geocoded": np.nan})
    for m, v in links["match_method"].value_counts().items():
        summ.append({"table": "record_match_method", "key": m, "n": int(v), "n_geocoded": np.nan})
    for k, v in rejected.items():
        summ.append({"table": "edges_rejected_by_constraint", "key": k, "n": int(v), "n_geocoded": np.nan})
    for k, v in f["sources"].value_counts().items():
        summ.append({"table": "facilities_by_source_combination", "key": k, "n": int(v), "n_geocoded": np.nan})
    for k, v in f["geocode_precision"].value_counts().items():
        summ.append({"table": "facilities_by_geocode_precision", "key": k, "n": int(v), "n_geocoded": np.nan})
    for k, v in f["resolution_tier"].value_counts().items():
        summ.append({"table": "facilities_by_resolution_tier", "key": k, "n": int(v), "n_geocoded": np.nan})
    for k, v in f["primary_kind"].value_counts().items():
        summ.append({"table": "facilities_by_primary_kind", "key": k, "n": int(v), "n_geocoded": np.nan})
    summ += [{"table": "outputs", "key": "facilities", "n": len(f), "n_geocoded": int((f["geocode_precision"] != "none").sum())},
             {"table": "outputs", "key": "clinic_registry", "n": len(cr), "n_geocoded": int((cr["geocode_precision"] != "none").sum())},
             {"table": "outputs", "key": "research_site_registry", "n": len(rr), "n_geocoded": int((rr["geocode_precision"] != "none").sum())},
             {"table": "outputs", "key": "facility_trials", "n": len(ft), "n_geocoded": np.nan},
             {"table": "outputs", "key": "facility_nih_projects", "n": len(fn), "n_geocoded": np.nan},
             {"table": "outputs", "key": "clinic_candidates_geocoded", "n": int(f["is_clinic_candidate"].sum()), "n_geocoded": np.nan}]
    # multi-source facilities: how many research facilities also carry a clinical (NPPES/HRSA) record
    rsrc = rr_out["sources"]
    summ += [{"table": "research_sites", "key": "with_nppes_or_hrsa_member", "n": int(rsrc.str.contains("nppes|hrsa").sum()), "n_geocoded": np.nan},
             {"table": "research_sites", "key": "with_trials", "n": int((rr["n_trials"] > 0).sum()), "n_geocoded": np.nan},
             {"table": "research_sites", "key": "with_nih", "n": int((rr["n_nih_projects"] > 0).sum()), "n_geocoded": np.nan},
             {"table": "research_sites", "key": "with_trials_and_nih", "n": int(((rr["n_trials"] > 0) & (rr["n_nih_projects"] > 0)).sum()), "n_geocoded": np.nan},
             {"table": "research_sites", "key": "with_precision_view_record", "n": int(rr["has_precision_view_record"].sum()), "n_geocoded": np.nan},
             {"table": "research_sites", "key": "only_mesh_expanded_trials_or_terms_only_nih", "n": int((~rr["has_precision_view_record"]).sum()), "n_geocoded": np.nan}]
    # consistency checks against the upstream tables
    us_ident_nct = site_rows["nct_id"].nunique()
    summ += [{"table": "consistency", "key": "us_identifiable_trials_upstream", "n": int(us_ident_nct), "n_geocoded": np.nan},
             {"table": "consistency", "key": "us_identifiable_trials_in_facility_trials", "n": int(ft["nct_id"].nunique()), "n_geocoded": np.nan},
             {"table": "consistency", "key": "nih_us_appl_ids_upstream", "n": int(nih_rows["appl_id"].nunique()), "n_geocoded": np.nan},
             {"table": "consistency", "key": "nih_appl_ids_in_facility_nih_projects", "n": int(fn["appl_id"].nunique()), "n_geocoded": np.nan}]
    pd.DataFrame(summ).to_csv(TABLES / "facility_resolution_summary.csv", index=False)

    # sample of links per method for human review
    rng = np.random.default_rng(SEED)
    lk = links[~links["match_method"].isin(["single_record", "anchor"])]
    rec = nodes.set_index("record_id")
    samp = []
    for m, g in lk.groupby("match_method"):
        take = g.iloc[rng.permutation(len(g))[:25]]
        for _, r in take.iterrows():
            o = rec.loc[r["linked_to_record_id"]] if r["linked_to_record_id"] in rec.index else None
            samp.append({"match_method": m, "match_confidence": r["match_confidence"],
                         "record_id": r["record_id"], "name": r["name"], "city": r["city"], "state": r["state"],
                         "zip5": r["zip5"], "linked_to_record_id": r["linked_to_record_id"],
                         "linked_name": None if o is None else o["name"],
                         "linked_city": None if o is None else o["city"],
                         "linked_zip5": None if o is None else o["zip5"],
                         "link_distance_km": r["link_distance_km"], "facility_id": r["facility_id"]})
    pd.DataFrame(samp).to_csv(TABLES / "facility_match_sample.csv", index=False)

    val_rows = []
    if not skip_validation:
        e_no = candidate_edges(nodes, hrsa, use_shared_npi=False)
        v1 = validate_hrsa_npi(nodes, hrsa, e_no)
        val_rows.append(v1)
        v2 = location_shuffle_false_links(nodes)
        v2.to_csv(TABLES / "facility_match_shuffle_control.csv", index=False)
        allrow = v2[v2["method"] == "ALL"].iloc[0]
        val_rows.append({"validation": "V2_location_shuffle_false_links_ctgov_nih",
                         "observed_cross_source_links": int(allrow["observed"]),
                         "shuffle_mean_links": float(allrow["shuffle_mean"]),
                         "estimated_false_link_share": float(allrow["estimated_false_link_share"])})
        pd.DataFrame(val_rows).to_csv(TABLES / "facility_match_validation.csv", index=False)
    stats = {"started": t0, "finished": utc_now_iso(), "n_facilities": len(f), "n_clinic_registry": len(cr),
             "n_research_site_registry": len(rr), "n_facility_trials": len(ft), "n_facility_nih": len(fn),
             "rejected": rejected, "validation": val_rows, "ctgov_filter": ct_stats}
    (TABLES / "facility_build_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    print(json.dumps(stats, indent=2, default=str))
    return stats


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Build facilities, clinic_registry and research_site_registry")
    ap.add_argument("--skip-validation", action="store_true")
    a = ap.parse_args(argv)
    run(skip_validation=a.skip_validation)


if __name__ == "__main__":
    main()
