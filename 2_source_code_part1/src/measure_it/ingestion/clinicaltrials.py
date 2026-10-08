"""ClinicalTrials.gov (API v2) ingestion for the target invisible-illness conditions.

Pipeline
--------
1. Record the API and data version: GET /api/v2/version -> {apiVersion, dataTimestamp}.
   The dataTimestamp (UTC) is this source's `source_version`.
2. For every condition in configs/conditions.yaml run ONE query per search term:
   ``query.cond="<term>"`` (phrase-quoted), all overall statuses, all years,
   pageSize=1000, following nextPageToken. Bare acronyms (e.g. "PASC", "ME/CFS")
   are never searched; they are logged as skipped in query_log.json.
3. Union and dedupe the matched NCT IDs, then fetch the full records in batches
   with ``filter.ids`` (fields=protocolSection,derivedSection,hasResults).
   Every API page is saved under data/raw/clinicaltrials_gov/ via download_file
   (so MANIFEST.json holds url, bytes, sha256 and retrieval time per page).
4. Write long tables: clinical_trials, trial_conditions, trial_interventions,
   trial_outcomes, trial_sites, plus a per-site-key history table
   ctgov_facility_summary (SPEC section G input; it was named
   research_site_registry__clinicaltrials_gov until 2026-09-23, which collided with
   store.union_partitions("research_site_registry"); the canonical
   research_site_registry is built by measure_it.facilities.registry).

Query semantics (measured against the live API on 2026-09-23)
-------------------------------------------------------------
* query.cond searches the "ConditionSearch" area: Condition, BriefTitle,
  OfficialTitle, ConditionMeshTerm, ConditionAncestorTerm and Keyword, all with
  synonym expansion (api_search_areas.json is saved as a raw file).
* Unquoted multi-word terms are AND-ed word by word: "post-COVID condition"
  returned 864 studies unquoted vs 35 quoted; the unquoted extras included
  "long-term care" COVID vaccine trials. Every quoted hit is also an unquoted hit,
  so phrases are always quoted. Quoting does cost some recall (e.g. trials registered
  as "Long COVID-19" are missed by "long COVID"); quoting_check() measures how many
  unquoted-only hits name a condition term, and the audit reports it.
* Even quoted phrases are expanded through MeSH ancestors and synonyms
  ("autonomic dysfunction" returns multiple system atrophy and Parkinson trials).
  trial_conditions therefore records `match_fields` / `literal_match`: whether
  the search term literally appears in the record's conditions, keywords or
  titles, or the study was matched only through the API's expansion.

Guardrail
---------
Registration of a trial that uses, measures or tests something is NOT evidence
that the measurement works, that a device diagnoses the condition, or that an
intervention is effective. Every row carries this note in provenance_notes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests

from ..config import PROCESSED, UNKNOWN, load_config, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import processed_path, read_table, write_table

SOURCE_ID = "clinicaltrials_gov"
SOURCE_NAME = "ClinicalTrials.gov API v2 (U.S. National Library of Medicine)"
API_BASE = "https://clinicaltrials.gov/api/v2"
STUDIES_URL = f"{API_BASE}/studies"
VERSION_URL = f"{API_BASE}/version"
SEARCH_AREAS_URL = f"{API_BASE}/studies/search-areas"
OAS_URL = "https://clinicaltrials.gov/api/oas/v2"
LANDING_URL = "https://clinicaltrials.gov/data-api/api"
RECORD_FIELDS = "protocolSection,derivedSection,hasResults"
PAGE_SIZE = 1000          # API maximum; larger values are coerced down to 1000
ID_BATCH_SIZE = 200       # NCT IDs per filter.ids request (~2.9 kB URL)
MIN_INTERVAL_S = 1.3      # polite spacing between live requests to clinicaltrials.gov
PRODUCER = "ingestion.clinicaltrials.run"
MODULE = "measure_it.ingestion.clinicaltrials"
LIST_SEP = " | "
GEOPOINT_METHOD = "ctgov_geopoint_city_level"

GUARDRAIL_NOTE = ("Trial registry record: registration of a trial that uses or measures something is not "
                  "evidence that the measurement works or that the intervention is effective.")

US_COUNTRY = "United States"
# ClinicalTrials.gov lists U.S. territories as their own countries.
US_TERRITORY_STATE_FIPS = {
    "Puerto Rico": "72", "Guam": "66", "Virgin Islands (U.S.)": "78", "U.S. Virgin Islands": "78",
    "American Samoa": "60", "Northern Mariana Islands": "69",
}

# Facility names that identify no facility (sponsor placeholders, non-facility locations). Flagged, not dropped.
GENERIC_FACILITY_RE = re.compile(
    r"^\s*[\"']?(?:"
    # "<Sponsor> Investigational Site", "Amgen Research Site" (sponsor prefix of <=2 words, never "clinical")
    r"(?:(?!clinical\b)[\w&.\-]+\s+){0,2}"
    r"(?:investigational|investigative|investigator|research|study|trial|clinical(?:\s+trial)?)\s+sites?"
    # "Research Site", "Clinical Research Site", "Clinical Trial Site", "Study Site"
    r"|(?:clinical\s+|study\s+|trial\s+)?(?:research|study|trial|clinical)\s+sites?"
    r"|local\s+institution"
    r"|sites?"
    r"|for\s+additional\s+information.*"
    r"|all\s+sites\s+listed\s+under.*"
    # facilities-review additions (2026-09-23): "Site Reference ID/Investigator# 56266", "SIte reference ID 304"
    r"|site\s+reference\s+id(?:\s*/\s*investigator\s*#?)?"
    # "Contact Medtronic for Exact Location(s)", "Contact <sponsor> for site locations"
    r"|contact\b.{0,80}?\bfor\b.{0,40}?\blocations?\b.*"
    # virtual / remote / online / decentralised studies and participants' homes are not facilities
    r"|.*\b(?:virtual|online|decentrali[sz]ed|proofpilot)\b.*"
    r"|.*\bremote\s+(?:trial|study|site|participation)\b.*"
    r"|.*\b(?:home[- ]based|telemedicine|telehealth)\s+(?:study|trial)\b.*"
    r"|.*\bno\s+physical\s+site\b.*"
    r"|.*\bdigital\s+research\s+platform\b.*"
    r"|.*\b(?:subjects?|participants?|patients?)'?s?\s+(?:own\s+)?(?:homes?|residences?)\b.*"
    r")"
    # optional trailing site numbers: "Site 101", "Investigational Site Number : 8400001", "Local Institution - 0012"
    r"(?:(?:\s*(?:number|no\.?|nr\.?|#|:|-|–))*\s*[\w\-/.]*\d[\w\-/.]*)*\s*[.,]?\s*$",
    re.IGNORECASE,
)

_last_live_request = [0.0]


# --------------------------------------------------------------------------------------------
# small pure helpers (unit-tested)
# --------------------------------------------------------------------------------------------

def is_bare_acronym(term: str) -> bool:
    """True for single-token, mostly-uppercase terms such as 'PASC', 'ME/CFS', 'hEDS'.

    Such terms are never sent as searches: they collide with unrelated meanings.
    """
    t = term.strip()
    if not t or re.search(r"\s", t):
        return False
    upper = sum(c.isupper() for c in t)
    lower = sum(c.islower() for c in t)
    return upper >= 2 and upper >= lower


def query_string(term: str) -> str:
    """The exact query.cond value sent for a search term (phrase-quoted Essie expression)."""
    t = term.strip().replace('"', "")
    return f'"{t}"'


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def _is_missing(value) -> bool:
    """None, NaN, pd.NA or NaT (values read back from a DataFrame are NaN, not None)."""
    if value is None:
        return True
    if isinstance(value, str):
        return False
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def normalize_text(text: str | None) -> str:
    """Lower-case; hyphens, slashes, underscores and punctuation become spaces; collapse whitespace.

    Missing values (None/NaN/NA) normalize to "" -- never to the string "nan".
    """
    if _is_missing(text):
        return ""
    s = str(text).lower()
    s = re.sub(r"[‐-―\-_/,;:()\[\]'\"’.]", " ", s)
    return re.sub(r"\s+", " ", s).strip()


# Spelling variants of one disease/virus token that a literal check must not treat as different words:
# "COVID-19"/"COVID 19"/"COVID19" -> "covid"; "SARS-CoV2"/"SARS CoV 2" -> "sars cov 2".
_TOKEN_CANON = [(re.compile(r"(?<![a-z0-9])covid ?19(?![a-z0-9])"), "covid"),
                (re.compile(r"(?<![a-z0-9])sars ?cov ?2(?![a-z0-9])"), "sars cov 2")]


def _match_text(text) -> str:
    s = normalize_text(text)
    for rx, rep in _TOKEN_CANON:
        s = rx.sub(rep, s)
    return s


def literal_term_in(term: str, text: str | None) -> bool:
    """Does the (normalized) term occur in text, starting at a word boundary?

    Normalization folds case, punctuation/hyphens and the COVID-19 / SARS-CoV-2 spelling variants, so
    "long COVID" is found in "Long COVID-19" and "post-COVID-19 condition" in "Post-COVID Condition".
    """
    nt = _match_text(term)
    if not nt:
        return False
    return re.search(r"(?<![a-z0-9])" + re.escape(nt), _match_text(text)) is not None


def join_list(values) -> str | None:
    vals = [str(v).strip() for v in (values or []) if v is not None and str(v).strip()]
    return LIST_SEP.join(vals) if vals else None


def split_list(value) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return [v for v in str(value).split(LIST_SEP) if v]


def date_year(date_str: str | None) -> int | None:
    m = re.match(r"^(\d{4})", str(date_str or ""))
    return int(m.group(1)) if m else None


def is_generic_facility(name: str | None, city: str | None = None, state: str | None = None) -> bool:
    """True for a name that identifies no facility: a sponsor placeholder or non-facility location
    (GENERIC_FACILITY_RE), or, when the site's city is given, a name that is only the city
    (optionally with its state), e.g. "Denver" in Denver, Colorado."""
    if _is_missing(name) or not str(name).strip():
        return False
    if GENERIC_FACILITY_RE.match(str(name)) is not None:
        return True
    if city is not None and not _is_missing(city):
        n, c = normalize_text(name), normalize_text(city)
        st = "" if _is_missing(state) else normalize_text(state)
        if n and c and n in {c, f"{c} {st}".strip(), f"{c} usa", f"{c} {st} usa".strip()}:
            return True
    return False


def site_key(facility, city, state, country) -> str:
    return "|".join(normalize_text(x) for x in (facility, city, state, country))


def _get(d, *path, default=None):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


# --------------------------------------------------------------------------------------------
# retrieval
# --------------------------------------------------------------------------------------------

def _studies_url(params: dict) -> str:
    return f"{STUDIES_URL}?{urlencode(params)}"


def _fetch(url: str, filename: str, *, refresh: bool = False) -> Path:
    """download_file with a polite throttle and retry on 429/5xx; skips files already in the manifest."""
    dest = raw_dir(SOURCE_ID) / filename
    if dest.exists() and not refresh and filename in load_manifest(SOURCE_ID)["files"]:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    last_exc: Exception | None = None
    for attempt in range(6):
        wait = _last_live_request[0] + MIN_INTERVAL_S - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        _last_live_request[0] = time.monotonic()
        try:
            return download_file(url, SOURCE_ID, filename, refresh=refresh, timeout=300)
        except requests.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else None
            if status in (429, 500, 502, 503, 504):
                last_exc = exc
                time.sleep(min(120, 5 * 2 ** attempt))
                continue
            raise
        except (requests.ConnectionError, requests.Timeout) as exc:
            last_exc = exc
            time.sleep(min(120, 5 * 2 ** attempt))
    raise RuntimeError(f"giving up on {url}: {last_exc}")


def fetch_metadata(refresh: bool = False) -> dict:
    """Version (data timestamp), search-area definitions and the OpenAPI spec, saved as raw files."""
    vpath = _fetch(VERSION_URL, "api_version.json", refresh=refresh)
    _fetch(SEARCH_AREAS_URL, "api_search_areas.json", refresh=refresh)
    _fetch(OAS_URL, "api_oas_v2.yaml", refresh=refresh)
    v = json.loads(vpath.read_text())
    m = load_manifest(SOURCE_ID)["files"]["api_version.json"]
    return {"api_version": v.get("apiVersion"), "data_timestamp": v.get("dataTimestamp"),
            "version_retrieved_at": m["retrieved_at"]}


def search_plan(conditions: list[dict] | None = None) -> pd.DataFrame:
    """One row per (condition, search term): the exact query string, or why it was skipped."""
    conds = conditions if conditions is not None else load_config("conditions")["conditions"]
    rows = []
    for c in conds:
        for term in c.get("search_terms", []) or []:
            skip = is_bare_acronym(term)
            rows.append({
                "condition_id": c["id"], "preferred_name": c.get("preferred_name"), "term": term,
                "query_param": "query.cond", "query_string": None if skip else query_string(term),
                "skipped": skip, "skip_reason": "bare acronym (never searched)" if skip else "",
            })
    return pd.DataFrame(rows)


def fetch_term(condition_id: str, term: str, *, refresh: bool = False) -> dict:
    """Run query.cond="<term>" across all statuses/years, following nextPageToken; returns matched IDs."""
    base = {"query.cond": query_string(term), "fields": "NCTId", "pageSize": PAGE_SIZE, "countTotal": "true"}
    for attempt in range(2):
        ids: list[str] = []
        files: list[str] = []
        urls: list[str] = []
        total = None
        token = None
        page = 1
        try:
            while True:
                params = dict(base)
                if token:
                    params["pageToken"] = token
                url = _studies_url(params)
                fname = f"queries/{condition_id}__{slugify(term)}__p{page:03d}.json"
                path = _fetch(url, fname, refresh=refresh or attempt > 0)
                d = json.loads(path.read_text())
                if page == 1:
                    total = d.get("totalCount")
                ids += [s["protocolSection"]["identificationModule"]["nctId"] for s in d.get("studies", [])]
                files.append(fname)
                urls.append(url)
                token = d.get("nextPageToken")
                if not token:
                    break
                page += 1
            break
        except requests.HTTPError:
            # A stale page token from an interrupted earlier run: redo the whole chain once.
            if attempt == 0 and page > 1:
                continue
            raise
    retrieved = load_manifest(SOURCE_ID)["files"][files[0]]["retrieved_at"]
    return {"condition_id": condition_id, "term": term, "query_param": "query.cond",
            "query_string": query_string(term), "request_url_page1": urls[0], "n_pages": len(files),
            "total_count": total, "n_ids_retrieved": len(ids), "n_unique_ids": len(set(ids)),
            "raw_files": files, "retrieved_at": retrieved, "nct_ids": sorted(set(ids))}


QUOTING_CHECK_TERMS = [("long_covid", "post-COVID condition"), ("long_covid", "long COVID"),
                       ("dysautonomia", "autonomic dysfunction")]


EXTRA_FIELDS = "NCTId,BriefTitle,OfficialTitle,Acronym,Condition,Keyword"


def _light_trial(study: dict) -> dict:
    """The fields match_fields() reads, from a light (titles/conditions/keywords) API record."""
    ps = study.get("protocolSection", {})
    idm, cnd = ps.get("identificationModule", {}), ps.get("conditionsModule", {})
    return {"nct_id": idm.get("nctId"), "brief_title": idm.get("briefTitle"), "official_title": idm.get("officialTitle"),
            "acronym": idm.get("acronym"), "conditions": join_list(cnd.get("conditions")),
            "keywords": join_list(cnd.get("keywords"))}


def quoting_check(term_results: dict[tuple[str, str], set], *, cond_unions: dict[str, set] | None = None,
                  cond_terms: dict[str, list[str]] | None = None, refresh: bool = False) -> list[dict]:
    """Why phrases are quoted, and what quoting costs: unquoted vs quoted query.cond result sets for a few terms.

    Precision side: n_unquoted vs n_quoted. Recall side: the unquoted-only hits that NO quoted term of the
    same condition retrieved are fetched (titles, conditions, keywords) and counted when one of the
    condition's searched terms appears literally in them -- trials that quoting lost although they name
    the condition. All pages are saved under queries_unquoted_check/. Not used to select trials.
    """
    out = []
    for cid, term in QUOTING_CHECK_TERMS:
        quoted = term_results.get((cid, term))
        if quoted is None:
            continue
        ids: set[str] = set()
        token, page, total = None, 1, None
        while True:
            params = {"query.cond": term, "fields": "NCTId", "pageSize": PAGE_SIZE, "countTotal": "true"}
            if token:
                params["pageToken"] = token
            path = _fetch(_studies_url(params), f"queries_unquoted_check/{slugify(term)}__p{page:03d}.json",
                          refresh=refresh)
            d = json.loads(path.read_text())
            total = d.get("totalCount") if page == 1 else total
            ids |= {s["protocolSection"]["identificationModule"]["nctId"] for s in d.get("studies", [])}
            token = d.get("nextPageToken")
            if not token:
                break
            page += 1
        row = {"term": term, "unquoted_query": term, "quoted_query": query_string(term),
               "n_unquoted": len(ids), "n_quoted": len(quoted),
               "n_unquoted_not_in_quoted": len(ids - quoted), "n_quoted_not_in_unquoted": len(quoted - ids)}
        if cond_unions is not None and cond_terms is not None and cid in cond_unions:
            extra = sorted(ids - cond_unions[cid])
            n_lit, examples = 0, []
            for i in range(0, len(extra), ID_BATCH_SIZE):
                batch = extra[i:i + ID_BATCH_SIZE]
                h = hashlib.sha1(",".join(batch).encode()).hexdigest()[:16]
                params = {"filter.ids": ",".join(batch), "fields": EXTRA_FIELDS, "pageSize": PAGE_SIZE}
                path = _fetch(_studies_url(params),
                              f"queries_unquoted_check/{slugify(term)}__outside_{cid}_{h}.json", refresh=refresh)
                for st in json.loads(path.read_text()).get("studies", []):
                    lt = _light_trial(st)
                    if condition_literal_match(cond_terms[cid], lt):
                        n_lit += 1
                        if len(examples) < 5:
                            examples.append(f"{lt['nct_id']} \"{(lt['brief_title'] or '')[:70]}\"")
            row.update({"n_unquoted_only_outside_condition": len(extra),
                        "n_of_those_naming_condition_term": n_lit,
                        "examples_lost_by_quoting": "; ".join(examples)})
        out.append(row)
    return out


def fetch_records(nct_ids, *, refresh: bool = False) -> tuple[list[str], dict]:
    """Fetch full records for the NCT IDs in batches; returns (raw file names, coverage stats)."""
    ids = sorted(set(nct_ids))
    files: list[str] = []
    returned: set[str] = set()
    for i in range(0, len(ids), ID_BATCH_SIZE):
        batch = ids[i:i + ID_BATCH_SIZE]
        h = hashlib.sha1(",".join(batch).encode()).hexdigest()[:16]
        params = {"filter.ids": ",".join(batch), "fields": RECORD_FIELDS, "pageSize": PAGE_SIZE,
                  "countTotal": "true"}
        fname = f"records/batch_{h}.json"
        path = _fetch(_studies_url(params), fname, refresh=refresh)
        d = json.loads(path.read_text())
        if d.get("nextPageToken"):
            raise RuntimeError(f"{fname}: unexpected second page for a {len(batch)}-id batch")
        returned |= {s["protocolSection"]["identificationModule"]["nctId"] for s in d.get("studies", [])}
        files.append(fname)
    stats = {"requested": len(ids), "returned": len(returned & set(ids)),
             "missing": sorted(set(ids) - returned), "unexpected_extra": sorted(returned - set(ids))}
    return files, stats


# --------------------------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------------------------

def parse_study(study: dict) -> dict:
    """Flatten one API v2 study into {'trial': row, 'interventions': [...], 'outcomes': [...], 'sites': [...]}."""
    ps = study.get("protocolSection", {})
    ds = study.get("derivedSection", {})
    idm = ps.get("identificationModule", {})
    stm = ps.get("statusModule", {})
    spm = ps.get("sponsorCollaboratorsModule", {})
    dsc = ps.get("descriptionModule", {})
    cnd = ps.get("conditionsModule", {})
    dgn = ps.get("designModule", {})
    aim = ps.get("armsInterventionsModule", {})
    out = ps.get("outcomesModule", {})
    elg = ps.get("eligibilityModule", {})
    clm = ps.get("contactsLocationsModule", {})
    nct = idm.get("nctId")

    collaborators = spm.get("collaborators", []) or []
    interventions = aim.get("interventions", []) or []
    arm_groups = aim.get("armGroups", []) or []
    locations = clm.get("locations", []) or []
    countries = sorted({loc.get("country") for loc in locations if loc.get("country")})

    trial = {
        "nct_id": nct,
        "brief_title": idm.get("briefTitle"),
        "official_title": idm.get("officialTitle"),
        "acronym": idm.get("acronym"),
        "brief_summary": dsc.get("briefSummary"),
        "overall_status": stm.get("overallStatus"),
        "last_known_status": stm.get("lastKnownStatus"),
        "why_stopped": stm.get("whyStopped"),
        "start_date": _get(stm, "startDateStruct", "date"),
        "start_date_type": _get(stm, "startDateStruct", "type"),
        "primary_completion_date": _get(stm, "primaryCompletionDateStruct", "date"),
        "primary_completion_date_type": _get(stm, "primaryCompletionDateStruct", "type"),
        "completion_date": _get(stm, "completionDateStruct", "date"),
        "completion_date_type": _get(stm, "completionDateStruct", "type"),
        "study_first_post_date": _get(stm, "studyFirstPostDateStruct", "date"),
        "last_update_post_date": _get(stm, "lastUpdatePostDateStruct", "date"),
        "status_verified_date": stm.get("statusVerifiedDate"),
        "start_year": date_year(_get(stm, "startDateStruct", "date")),
        "study_type": dgn.get("studyType"),
        "phases": join_list(dgn.get("phases")),
        "allocation": _get(dgn, "designInfo", "allocation"),
        "intervention_model": _get(dgn, "designInfo", "interventionModel"),
        "primary_purpose": _get(dgn, "designInfo", "primaryPurpose"),
        "masking": _get(dgn, "designInfo", "maskingInfo", "masking"),
        "observational_model": _get(dgn, "designInfo", "observationalModel"),
        "time_perspective": _get(dgn, "designInfo", "timePerspective"),
        "enrollment_count": _get(dgn, "enrollmentInfo", "count"),
        "enrollment_type": _get(dgn, "enrollmentInfo", "type"),
        "lead_sponsor": _get(spm, "leadSponsor", "name"),
        "lead_sponsor_class": _get(spm, "leadSponsor", "class"),
        "collaborators": join_list(c.get("name") for c in collaborators),
        "collaborator_classes": join_list(c.get("class") for c in collaborators),
        "n_collaborators": len(collaborators),
        "conditions": join_list(cnd.get("conditions")),
        "keywords": join_list(cnd.get("keywords")),
        "condition_mesh_ids": join_list(m.get("id") for m in _get(ds, "conditionBrowseModule", "meshes", default=[]) or []),
        "condition_mesh_terms": join_list(m.get("term") for m in _get(ds, "conditionBrowseModule", "meshes", default=[]) or []),
        "condition_mesh_ancestors": join_list(m.get("term") for m in _get(ds, "conditionBrowseModule", "ancestors", default=[]) or []),
        "intervention_types": join_list(sorted({i.get("type") for i in interventions if i.get("type")})),
        "intervention_names": join_list(i.get("name") for i in interventions),
        "n_interventions": len(interventions),
        "arm_groups_json": json.dumps(
            [{"label": a.get("label"), "type": a.get("type"), "description": a.get("description"),
              "intervention_names": a.get("interventionNames", [])} for a in arm_groups],
            ensure_ascii=False) if arm_groups else None,
        "n_arm_groups": len(arm_groups),
        "n_primary_outcomes": len(out.get("primaryOutcomes", []) or []),
        "n_secondary_outcomes": len(out.get("secondaryOutcomes", []) or []),
        "n_other_outcomes": len(out.get("otherOutcomes", []) or []),
        "eligibility_criteria": elg.get("eligibilityCriteria"),
        "healthy_volunteers": elg.get("healthyVolunteers"),
        "sex": elg.get("sex"),
        "minimum_age": elg.get("minimumAge"),
        "maximum_age": elg.get("maximumAge"),
        "std_ages": join_list(elg.get("stdAges")),
        "n_locations": len(locations),
        "n_us_locations": sum(1 for loc in locations if loc.get("country") == US_COUNTRY),
        "location_countries": join_list(countries),
        "has_results": bool(study.get("hasResults", False)),
        "ctgov_version_holder": _get(ds, "miscInfoModule", "versionHolder"),
        "study_url": f"https://clinicaltrials.gov/study/{nct}",
    }

    iv_rows = [{
        "nct_id": nct, "intervention_idx": k, "intervention_type": i.get("type"),
        "intervention_name": i.get("name"), "intervention_description": i.get("description"),
        "other_names": join_list(i.get("otherNames")), "arm_group_labels": join_list(i.get("armGroupLabels")),
    } for k, i in enumerate(interventions)]

    oc_rows = []
    for otype, key in (("primary", "primaryOutcomes"), ("secondary", "secondaryOutcomes"), ("other", "otherOutcomes")):
        for k, o in enumerate(out.get(key, []) or []):
            oc_rows.append({"nct_id": nct, "outcome_type": otype, "outcome_idx": k, "measure": o.get("measure"),
                            "description": o.get("description"), "time_frame": o.get("timeFrame")})

    site_rows = [{
        "nct_id": nct, "site_idx": k, "facility": loc.get("facility"), "city": loc.get("city"),
        "state": loc.get("state"), "zip": loc.get("zip"), "country": loc.get("country"),
        "site_status": loc.get("status"),
        "geopoint_lat": _get(loc, "geoPoint", "lat"), "geopoint_lon": _get(loc, "geoPoint", "lon"),
    } for k, loc in enumerate(locations)]   # contacts (named people, phones, emails) are deliberately not kept

    return {"trial": trial, "interventions": iv_rows, "outcomes": oc_rows, "sites": site_rows}


def match_fields(term: str, trial: dict) -> list[str]:
    """Record fields in which the search term appears literally (after normalization)."""
    hits = []
    for field, values in (("condition", split_list(trial.get("conditions"))),
                          ("keyword", split_list(trial.get("keywords"))),
                          ("brief_title", [trial.get("brief_title")]),
                          ("official_title", [trial.get("official_title")]),
                          ("acronym", [trial.get("acronym")]),
                          ("condition_mesh_term", split_list(trial.get("condition_mesh_terms"))),
                          ("condition_mesh_ancestor", split_list(trial.get("condition_mesh_ancestors")))):
        if any(literal_term_in(term, v) for v in values):
            hits.append(field)
    return hits


LITERAL_FIELDS = {"condition", "keyword", "brief_title", "official_title", "acronym"}


def condition_literal_match(terms, trial: dict) -> bool:
    """True when ANY of a condition's searched terms appears literally in the record's
    conditions, keywords, titles or acronym -- not only the term whose query retrieved it.

    (A trial retrieved by "long COVID" through synonym expansion whose conditions list
    "Post-COVID Conditions" is a literal match for the long_covid condition.)
    """
    return any(LITERAL_FIELDS & set(match_fields(t, trial)) for t in terms)


# --------------------------------------------------------------------------------------------
# geography
# --------------------------------------------------------------------------------------------

def _state_name_to_fips() -> dict[str, str]:
    g = read_table("geographies", columns=["geo_level", "state_name", "state_fips", "ct_legacy"])
    st = g[(g["geo_level"] == "state")]
    return {normalize_text(n): f for n, f in zip(st["state_name"], st["state_fips"])}


def geocode_sites(sites: pd.DataFrame) -> pd.DataFrame:
    """Add country flags, state/county FIPS and geocode method columns to trial sites.

    lat/lon: the API geoPoint when present (resolution 'city_centroid': measured to be a city-level
    geocode, hence geocode_method='ctgov_geopoint_city_level'); else, for U.S. and territory
    sites with a ZIP, the ZCTA internal point (resolution 'zcta_centroid').

    county_fips, in order:
      1. the county of the site's own ZIP (ZIP5 -> ZCTA internal point, crosswalk.zip_to_geo) when it
         resolves and lies in the registered state -- except
      2. when the geoPoint's county covers >= 50% of that ZCTA's land (2020 ZCTA-county relationship
         file): then the city centroid and most of the ZIP area agree, and the geoPoint county is used
         (a ZCTA's internal point can fall in a minority county, e.g. ZCTA 32174 Ormond Beach);
      3. else point-in-polygon of the geoPoint against the 2024 cartographic county boundaries.
    ZIP comes first because geoPoint is the city centroid: sites whose postal city is "Atlanta" or
    "St Louis" but which lie in DeKalb or St. Louis County would otherwise all land in the core-city
    county. Both candidates are kept (county_fips_zip, county_fips_pip) with county_methods_agree.
    """
    import geopandas as gpd

    from ..geography.crosswalk import STATE_FIPS, zip_to_geo

    s = sites.copy()
    s["country_group"] = np.select(
        [s["country"] == US_COUNTRY, s["country"].isin(list(US_TERRITORY_STATE_FIPS))],
        ["US", "US_territory"], default="non_US")
    s.loc[s["country"].isna(), "country_group"] = "unknown_country"
    us_like = s["country_group"].isin(["US", "US_territory"])

    name_to_fips = _state_name_to_fips()
    rep = pd.Series(pd.NA, index=s.index, dtype="object")
    m_us = s["country_group"] == "US"
    rep[m_us] = s.loc[m_us, "state"].map(lambda x: name_to_fips.get(normalize_text(x)) if pd.notna(x) else None)
    m_ter = s["country_group"] == "US_territory"
    rep[m_ter] = s.loc[m_ter, "country"].map(US_TERRITORY_STATE_FIPS)
    s["reported_state_fips"] = rep.where(rep.notna(), None)

    # point-in-polygon for U.S./territory geoPoints
    counties = gpd.read_parquet(PROCESSED / "county_boundaries_2024.geoparquet")[["geo_id", "geometry"]]
    has_pt = s["geopoint_lat"].notna() & s["geopoint_lon"].notna()
    pm = us_like & has_pt
    s["county_fips_pip"] = None
    if pm.any():
        pts = gpd.GeoDataFrame(
            {"_i": s.index[pm]},
            geometry=gpd.points_from_xy(s.loc[pm, "geopoint_lon"], s.loc[pm, "geopoint_lat"]), crs=4326,
        ).to_crs(counties.crs)
        j = gpd.sjoin(pts, counties, how="left", predicate="within").drop_duplicates("_i")
        s.loc[j["_i"].values, "county_fips_pip"] = j["geo_id"].values
    s["county_fips_pip"] = s["county_fips_pip"].where(s["county_fips_pip"].notna(), None)

    # ZIP -> ZCTA internal point -> county
    s["county_fips_zip"] = None
    s["zcta_lat"] = np.nan
    s["zcta_lon"] = np.nan
    zm = us_like & s["zip"].notna()
    if zm.any():
        zg = zip_to_geo(s.loc[zm, "zip"])
        s.loc[zm, "county_fips_zip"] = zg["county_fips"].values
        s.loc[zm, "zcta_lat"] = zg["lat"].values
        s.loc[zm, "zcta_lon"] = zg["lon"].values
    s["county_fips_zip"] = s["county_fips_zip"].where(s["county_fips_zip"].notna(), None)

    # A ZIP whose county lies in another state than the registered one is treated as a typo.
    zip_state_ok = s["county_fips_zip"].notna() & (
        s["reported_state_fips"].isna() | (s["county_fips_zip"].str[:2] == s["reported_state_fips"]))
    s["zip_state_conflict"] = s["county_fips_zip"].notna() & ~zip_state_ok
    # geoPoint county's share of the ZIP's ZCTA land area (2020 relationship file; CT 2024 regions never match)
    xw = read_table("zcta_county_crosswalk", columns=["zcta", "county_fips_2020", "land_share_of_zcta"])
    share = dict(zip(zip(xw["zcta"], xw["county_fips_2020"]), xw["land_share_of_zcta"]))
    z5 = s["zip"].astype(str).str.extract(r"^(\d{5})")[0]
    s["pip_county_share_of_zcta"] = [
        share.get((z, c)) if isinstance(z, str) and isinstance(c, str) else None
        for z, c in zip(z5, s["county_fips_pip"])]
    s["pip_county_share_of_zcta"] = pd.to_numeric(s["pip_county_share_of_zcta"])
    pip_majority = (zip_state_ok & s["county_fips_pip"].notna() & (s["county_fips_pip"] != s["county_fips_zip"])
                    & (s["pip_county_share_of_zcta"] >= 0.5))
    s["county_fips"] = s["county_fips_zip"].where(zip_state_ok & ~pip_majority, s["county_fips_pip"])
    s["county_fips"] = s["county_fips"].where(s["county_fips"].notna(), None)
    s["county_assignment_method"] = np.select(
        [pip_majority, zip_state_ok, s["county_fips_pip"].notna()],
        ["geopoint_county_majority_of_zcta", "zip5_as_zcta_internal_point", "geopoint_point_in_polygon_cb_2024_5m"],
        default="none")
    both = s["county_fips_pip"].notna() & s["county_fips_zip"].notna()
    s["county_methods_agree"] = pd.Series(pd.NA, index=s.index, dtype="boolean")
    s.loc[both, "county_methods_agree"] = (s.loc[both, "county_fips_pip"] == s.loc[both, "county_fips_zip"]).values

    use_zcta = ~has_pt & us_like & s["zcta_lat"].notna()
    s["lat"] = s["geopoint_lat"].where(has_pt, s["zcta_lat"].where(use_zcta))
    s["lon"] = s["geopoint_lon"].where(has_pt, s["zcta_lon"].where(use_zcta))
    s["geocode_method"] = np.select([has_pt, use_zcta], [GEOPOINT_METHOD, "zip5_as_zcta_internal_point"],
                                    default="not_geocoded")
    s["state_fips"] = s["county_fips"].str[:2].where(s["county_fips"].notna(), s["reported_state_fips"])
    s["state_fips"] = s["state_fips"].where(s["state_fips"].notna(), None)
    s["state_abbr"] = s["state_fips"].map(STATE_FIPS)
    s["pip_state_matches_reported"] = pd.Series(pd.NA, index=s.index, dtype="boolean")
    chk = s["county_fips_pip"].notna() & s["reported_state_fips"].notna()
    s.loc[chk, "pip_state_matches_reported"] = (
        s.loc[chk, "county_fips_pip"].str[:2] == s.loc[chk, "reported_state_fips"]).values
    s["source_geographic_resolution"] = np.select(
        [has_pt, use_zcta, us_like & s["state_fips"].notna()], ["city_centroid", "zcta_centroid", "state"],
        default="none")
    s = s.drop(columns=["zcta_lat", "zcta_lon"])
    return s


# --------------------------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------------------------

def _load_records(files: list[str]) -> tuple[list[dict], dict[str, str]]:
    manifest = load_manifest(SOURCE_ID)["files"]
    studies, retrieved = [], {}
    seen = set()
    for f in files:
        d = json.loads((raw_dir(SOURCE_ID) / f).read_text())
        for s in d.get("studies", []):
            nct = s["protocolSection"]["identificationModule"]["nctId"]
            if nct in seen:
                continue
            seen.add(nct)
            studies.append(s)
            retrieved[nct] = manifest[f]["retrieved_at"]
    return studies, retrieved


def _provenance(df: pd.DataFrame, *, version: str, retrieved: pd.Series | str, record_id,
                resolution="none", notes: str = GUARDRAIL_NOTE) -> pd.DataFrame:
    out = add_provenance(df, data_layer="facility", source_name=SOURCE_NAME, source_version=version,
                         retrieved_at=retrieved if isinstance(retrieved, str) else "pending",
                         evidence_type="trial_registry_record", source_record_id=record_id,
                         source_geographic_resolution="none", provenance_notes=notes)
    if not isinstance(retrieved, str):
        out["retrieved_at"] = retrieved.values
    if isinstance(resolution, pd.Series):
        out["source_geographic_resolution"] = resolution.values
    else:
        out["source_geographic_resolution"] = resolution
    return out


def build_site_history(trials: pd.DataFrame, sites: pd.DataFrame, conds: pd.DataFrame) -> pd.DataFrame:
    """Per-facility trial history (SPEC section G): counts, recruiting, conditions, recency, device experience.

    A facility is keyed by normalized (facility, city, state, country); spelling variants of one
    institution stay separate rows, and sponsor placeholder names are flagged `facility_generic`.
    Sites registered without a facility name share one key per city ("|city|state|country"); such a row
    is a location bucket, not a facility, and is flagged `facility_missing`. Use `facility_identifiable`
    (neither generic nor missing) for facility-level matching.
    """
    t = trials[["nct_id", "overall_status", "study_type", "start_year", "last_update_post_date",
                "intervention_types"]].copy()
    t["device_or_diagnostic"] = t["intervention_types"].fillna("").str.contains(r"\bDEVICE\b|\bDIAGNOSTIC_TEST\b")
    cmap = conds.groupby("nct_id")["condition_id"].agg(lambda x: sorted(set(x)))
    s = sites.merge(t, on="nct_id", how="left")
    s["site_key"] = [site_key(*r) for r in s[["facility", "city", "state", "country"]].itertuples(index=False)]
    s["condition_ids"] = s["nct_id"].map(cmap)

    def mode(x):
        x = x.dropna()
        return x.mode().iloc[0] if len(x) else None

    rows = []
    for key, g in s.groupby("site_key", sort=False):
        ncts = sorted(set(g["nct_id"]))
        cids = sorted({c for lst in g["condition_ids"] if isinstance(lst, list) for c in lst})
        rows.append({
            "site_key": key, "facility": mode(g["facility"]), "city": mode(g["city"]), "state": mode(g["state"]),
            "zip": mode(g["zip"]), "country": mode(g["country"]), "country_group": mode(g["country_group"]),
            "lat": mode(g["lat"]), "lon": mode(g["lon"]), "geocode_method": mode(g["geocode_method"]),
            "county_fips": mode(g["county_fips"]), "state_fips": mode(g["state_fips"]),
            "facility_generic": bool(g["facility_generic"].any()),
            "facility_missing": bool(g["facility_missing"].all()),
            "facility_identifiable": bool(g["facility_identifiable"].all()),
            "n_trials": len(ncts),
            "n_trials_site_recruiting": int(g.loc[g["site_status"] == "RECRUITING", "nct_id"].nunique()),
            "n_trials_overall_recruiting": int(g.loc[g["overall_status"] == "RECRUITING", "nct_id"].nunique()),
            "n_interventional": int(g.loc[g["study_type"] == "INTERVENTIONAL", "nct_id"].nunique()),
            "n_observational": int(g.loc[g["study_type"] == "OBSERVATIONAL", "nct_id"].nunique()),
            "n_device_or_diagnostic_trials": int(g.loc[g["device_or_diagnostic"].fillna(False), "nct_id"].nunique()),
            "condition_ids": join_list(cids), "n_conditions": len(cids),
            "first_start_year": g["start_year"].min() if g["start_year"].notna().any() else None,
            "latest_start_year": g["start_year"].max() if g["start_year"].notna().any() else None,
            "latest_update_post_date": g["last_update_post_date"].dropna().max() if g["last_update_post_date"].notna().any() else None,
            "nct_ids": join_list(ncts),
            "source_geographic_resolution": mode(g["source_geographic_resolution"]) or "none",
            "retrieved_at": g["retrieved_at"].max(),
        })
    out = pd.DataFrame(rows)
    for c in ("first_start_year", "latest_start_year"):
        out[c] = pd.to_numeric(out[c]).astype("Int64")
    return out


def run(refresh: bool = False) -> dict:
    """Fetch, parse, geocode and write all ClinicalTrials.gov tables; returns measured stats."""
    meta = fetch_metadata(refresh=refresh)
    version = f"ClinicalTrials.gov API {meta['api_version']}; dataTimestamp {meta['data_timestamp']} UTC"

    plan = search_plan()
    log = []
    for r in plan.itertuples(index=False):
        if r.skipped:
            log.append({"condition_id": r.condition_id, "term": r.term, "query_param": r.query_param,
                        "query_string": None, "skipped": True, "skip_reason": r.skip_reason})
            continue
        res = fetch_term(r.condition_id, r.term, refresh=refresh)
        res["skipped"] = False
        res["skip_reason"] = ""
        if res["total_count"] is not None and res["n_unique_ids"] != res["total_count"]:
            res["count_warning"] = f"retrieved {res['n_unique_ids']} unique of totalCount {res['total_count']}"
        log.append(res)
        print(f"[ctgov] {r.condition_id:<24} {r.query_string:<45} total={res['total_count']}", flush=True)

    # the searched (non-skipped) terms of each condition, and each condition's union of matched NCT IDs
    cond_terms: dict[str, list[str]] = {}
    cond_unions: dict[str, set] = {}
    for e in log:
        if not e["skipped"]:
            cond_terms.setdefault(e["condition_id"], []).append(e["term"])
            cond_unions.setdefault(e["condition_id"], set()).update(e["nct_ids"])
    qcheck = quoting_check({(e["condition_id"], e["term"]): set(e["nct_ids"]) for e in log if not e["skipped"]},
                           cond_unions=cond_unions, cond_terms=cond_terms, refresh=refresh)
    all_ids = sorted({i for e in log if not e["skipped"] for i in e["nct_ids"]})
    print(f"[ctgov] union of matched NCT IDs: {len(all_ids)}", flush=True)
    rec_files, coverage = fetch_records(all_ids, refresh=refresh)
    studies, retrieved = _load_records(rec_files)

    parsed = [parse_study(s) for s in studies]
    trials = pd.DataFrame([p["trial"] for p in parsed])
    trials = trials[trials["nct_id"].isin(all_ids)].reset_index(drop=True)
    interventions = pd.DataFrame([r for p in parsed for r in p["interventions"]])
    outcomes = pd.DataFrame([r for p in parsed for r in p["outcomes"]])
    sites = pd.DataFrame([r for p in parsed for r in p["sites"]])
    for df in (interventions, outcomes, sites):
        df.drop(df.index[~df["nct_id"].isin(all_ids)], inplace=True)
    trials["enrollment_count"] = pd.to_numeric(trials["enrollment_count"]).astype("Int64")
    trials["start_year"] = pd.to_numeric(trials["start_year"]).astype("Int64")

    # --- trial_conditions: one row per NCT x condition x matched term
    tr_by_id = trials.set_index("nct_id").to_dict("index")
    names = {c["id"]: c.get("preferred_name") for c in load_config("conditions")["conditions"]}
    tc_rows = []
    cond_lit_cache: dict[tuple[str, str], bool] = {}
    for e in log:
        if e["skipped"]:
            continue
        for nct in e["nct_ids"]:
            if nct not in tr_by_id:
                continue
            hits = match_fields(e["term"], tr_by_id[nct])
            ck = (nct, e["condition_id"])
            if ck not in cond_lit_cache:
                cond_lit_cache[ck] = condition_literal_match(cond_terms[e["condition_id"]], tr_by_id[nct])
            tc_rows.append({"nct_id": nct, "condition_id": e["condition_id"],
                            "condition_preferred_name": names.get(e["condition_id"]),
                            "matched_term": e["term"], "query_param": "query.cond",
                            "query_string": e["query_string"], "match_fields": join_list(hits),
                            "literal_match": bool(LITERAL_FIELDS & set(hits)),
                            "condition_literal_match": cond_lit_cache[ck],
                            "query_retrieved_at": e["retrieved_at"]})
    tconds = pd.DataFrame(tc_rows)

    agg = tconds.groupby("nct_id").agg(
        matched_condition_ids=("condition_id", lambda x: join_list(sorted(set(x)))),
        matched_terms=("matched_term", lambda x: join_list(sorted(set(x)))),
        n_matched_conditions=("condition_id", "nunique"),
        any_literal_match=("condition_literal_match", "any"))
    trials = trials.merge(agg, left_on="nct_id", right_index=True, how="left")

    # --- sites
    sites = geocode_sites(sites)
    sites["facility_generic"] = pd.Series(
        [is_generic_facility(f, c, st) for f, c, st in zip(sites["facility"], sites["city"], sites["state"])],
        index=sites.index, dtype=bool)
    sites["facility_missing"] = sites["facility"].map(lambda x: _is_missing(x) or not str(x).strip()).astype(bool)
    # a site names an identifiable facility only when the name is present and is not a sponsor placeholder
    sites["facility_identifiable"] = ~(sites["facility_generic"] | sites["facility_missing"])
    sites["retrieved_at"] = sites["nct_id"].map(retrieved)
    rec_ret = trials["nct_id"].map(retrieved)

    # --- provenance + write
    geo_note = ("lat/lon = ClinicalTrials.gov geoPoint (ClinicalTrials.gov's own geocode; measured to be the city "
                "centroid, not the facility address) when geocode_method='ctgov_geopoint_city_level'; else ZIP->ZCTA "
                "internal point. county_fips = county of the site's ZIP (ZCTA internal point) when it resolves in the "
                "registered state, unless the geoPoint county holds >=50% of that ZCTA's land; else point-in-polygon of "
                "geoPoint against cb_2024_5m counties (see county_assignment_method). " + GUARDRAIL_NOTE)
    trials_p = _provenance(trials, version=version, retrieved=rec_ret, record_id="nct_id")
    tconds_p = _provenance(tconds, version=version, retrieved=tconds["query_retrieved_at"],
                           record_id=lambda d: d["nct_id"] + ":" + d["condition_id"] + ":" + d["matched_term"],
                           notes="Matched by ClinicalTrials.gov query.cond (Essie: synonyms + MeSH ancestors); "
                                 "literal_match=False means this row's term is not in the record's conditions/"
                                 "keywords/titles/acronym; condition_literal_match=False means none of the "
                                 "condition's searched terms is. " + GUARDRAIL_NOTE)
    iv_p = _provenance(interventions.reset_index(drop=True), version=version,
                       retrieved=interventions["nct_id"].map(retrieved).reset_index(drop=True),
                       record_id=lambda d: d["nct_id"] + ":intervention:" + d["intervention_idx"].astype(str))
    oc_p = _provenance(outcomes.reset_index(drop=True), version=version,
                       retrieved=outcomes["nct_id"].map(retrieved).reset_index(drop=True),
                       record_id=lambda d: d["nct_id"] + ":" + d["outcome_type"] + "_outcome:" + d["outcome_idx"].astype(str))
    sites = sites.reset_index(drop=True)
    res_col = sites.pop("source_geographic_resolution")
    ret_col = sites.pop("retrieved_at")
    sites_p = _provenance(sites, version=version, retrieved=ret_col,
                          record_id=lambda d: d["nct_id"] + ":site:" + d["site_idx"].astype(str),
                          resolution=res_col, notes=geo_note)

    hist = build_site_history(trials, sites_p.assign(source_geographic_resolution=res_col.values), tconds)
    h_res = hist.pop("source_geographic_resolution")
    h_ret = hist.pop("retrieved_at")
    hist_p = _provenance(hist, version=version, retrieved=h_ret, record_id="site_key", resolution=h_res,
                         notes="Aggregated from trial_sites for trials matching the target conditions; facility "
                               "identity = normalized (facility, city, state, country), so spelling variants are "
                               "separate rows. Research readiness signal, not patient burden. " + GUARDRAIL_NOTE)

    write_table(trials_p, "clinical_trials", producer=PRODUCER,
                description="One row per NCT study matching any target-condition search term (ClinicalTrials.gov API v2)")
    write_table(tconds_p, "trial_conditions", producer=PRODUCER,
                description="NCT x condition_id x matched search term, with the exact query string and literal-match flag")
    write_table(iv_p, "trial_interventions", producer=PRODUCER, description="One row per trial intervention")
    write_table(oc_p, "trial_outcomes", producer=PRODUCER,
                description="One row per registered primary/secondary/other outcome measure")
    write_table(sites_p, "trial_sites", producer=PRODUCER,
                description="One row per NCT x location, with geoPoint/ZIP geocode and 2024 county FIPS")
    write_table(hist_p, "ctgov_facility_summary", producer=PRODUCER,
                description="Per-site-key ClinicalTrials.gov site history for target-condition trials (SPEC G input; "
                            "formerly research_site_registry__clinicaltrials_gov)")
    stale = processed_path("research_site_registry__clinicaltrials_gov")  # pre-2026-09-23 name of this table
    if stale.exists():
        stale.unlink()
        stale.with_suffix(".meta.json").unlink(missing_ok=True)

    # --- query log (raw-side record of every query string sent)
    qlog = {"source_id": SOURCE_ID, "generated_at": utc_now_iso(), "api": meta, "record_fields": RECORD_FIELDS,
            "id_batch_size": ID_BATCH_SIZE, "record_coverage": coverage, "quoting_check": qcheck,
            "queries": [{k: v for k, v in e.items() if k != "nct_ids"} for e in log]}
    (raw_dir(SOURCE_ID) / "query_log.json").write_text(json.dumps(qlog, indent=2))

    stats = compute_stats(meta, log, coverage, trials_p, tconds_p, iv_p, oc_p, sites_p, hist_p)
    stats["quoting_check"] = qcheck
    write_audit(stats)
    write_registry(stats)
    print(json.dumps({k: v for k, v in stats.items() if isinstance(v, (int, str, float))}, indent=1))
    return stats


# --------------------------------------------------------------------------------------------
# audit + registry
# --------------------------------------------------------------------------------------------

def _miss(df: pd.DataFrame, cols: list[str], mask: pd.Series | None = None) -> list[dict]:
    d = df if mask is None else df[mask]
    out = []
    for c in cols:
        v = d[c]
        n_missing = int((v.isna() | (v.astype(str).str.strip() == "")).sum())
        out.append({"variable": c, "n": len(d), "missing": n_missing,
                    "pct": round(100 * n_missing / len(d), 1) if len(d) else None})
    return out


# MeSH descriptors whose share of a condition's non-literal matches the audit reports by name
NAMED_EXPANSION_MESH = {"dysautonomia": ["Autonomic Nervous System Diseases"],
                        "mcas": ["Mast Cell Activation Disorders", "Mastocytosis"]}


def compute_stats(meta, log, coverage, trials, tconds, iv, oc, sites, hist) -> dict:
    us = sites[sites["country_group"] == "US"]
    us_ct = us["reported_state_fips"] == "09"
    snapshot_date = str(meta["data_timestamp"] or "")[:10]
    # Connecticut county vintages from the shared geographies table (2024 planning regions vs legacy counties)
    geo = read_table("geographies", columns=["geo_id", "geo_level", "state_fips", "ct_legacy"])
    ct_geo = geo[(geo["geo_level"] == "county") & (geo["state_fips"] == "09")]
    ct_2024 = set(ct_geo.loc[~ct_geo["ct_legacy"].astype(bool), "geo_id"])
    ct_legacy = set(ct_geo.loc[ct_geo["ct_legacy"].astype(bool), "geo_id"])
    us_pt = us[us["geocode_method"] == GEOPOINT_METHOD]
    # How precise is geoPoint? Distinct coordinates vs distinct (city, state) among U.S. geoPoint sites.
    cityset = us_pt.assign(_c=us_pt["city"].map(normalize_text) + "|" + us_pt["state"].map(normalize_text))
    coords = cityset.groupby(["geopoint_lat", "geopoint_lon"])
    fac_per_coord = coords["facility"].agg(lambda x: x.map(normalize_text).nunique())
    per_cond = tconds.groupby("condition_id").agg(n_trials=("nct_id", "nunique"))
    per_cond["n_trials_literal"] = (tconds[tconds["condition_literal_match"]]
                                    .groupby("condition_id")["nct_id"].nunique())
    per_cond["n_trials_literal"] = per_cond["n_trials_literal"].fillna(0).astype(int)
    per_cond = per_cond.reset_index()
    status = trials["overall_status"].value_counts(dropna=False).rename_axis("overall_status").reset_index(name="n_trials")
    cond_status = (tconds[["nct_id", "condition_id"]].drop_duplicates()
                   .merge(trials[["nct_id", "overall_status"]], on="nct_id")
                   .pivot_table(index="condition_id", columns="overall_status", values="nct_id",
                                aggfunc="nunique", fill_value=0))
    stype = trials["study_type"].value_counts(dropna=False).rename_axis("study_type").reset_index(name="n_trials")
    both = us["county_methods_agree"].notna()
    # what the API's expansion adds: per condition, trials with no literal term match
    from collections import Counter
    t_by = trials.set_index("nct_id")
    mesh_df = Counter(m for a, b in zip(trials["condition_mesh_terms"], trials["condition_mesh_ancestors"])
                      for m in set(split_list(a)) | set(split_list(b)))
    expansion = []
    named_counts: dict[str, str] = {}
    for cid, g in tconds.groupby("condition_id"):
        lit = g.groupby("nct_id")["condition_literal_match"].any()
        nl = t_by.loc[lit.index[~lit.values]]
        if nl.empty:
            continue
        conds_c = Counter(c.lower() for v in nl["conditions"] for c in split_list(v))
        mesh_c = Counter(m for a, b in zip(nl["condition_mesh_terms"], nl["condition_mesh_ancestors"])
                         for m in set(split_list(a)) | set(split_list(b)))
        # most specific MeSH (lowest frequency across all trials) carried by >= half of the non-literal trials
        shared = [m for m, c in mesh_c.items() if c >= 0.5 * len(nl)]
        spec = min(shared, key=lambda m: (mesh_df[m], m)) if shared else None
        for named in NAMED_EXPANSION_MESH.get(cid, []):
            named_counts[f"{cid}:{named}"] = f"{mesh_c.get(named, 0)} of {len(nl)}"
        expansion.append({"condition_id": cid, "n_trials": int(len(lit)), "n_non_literal": int(len(nl)),
                          "top_registered_conditions_of_non_literal": "; ".join(f"{k} ({v})" for k, v in conds_c.most_common(4)),
                          "most_specific_mesh_shared_by_half": (f"{spec} ({mesh_c[spec]} of {len(nl)})" if spec else None)})
    # per-term diagnostics: literal hits, and hits no other term of the same condition found
    term_rows = []
    for e in log:
        row = {k: e.get(k) for k in ("condition_id", "term", "query_string", "total_count", "n_unique_ids",
                                     "n_pages", "skipped", "skip_reason")}
        if not e["skipped"]:
            mine = tconds[(tconds["condition_id"] == e["condition_id"]) & (tconds["matched_term"] == e["term"])]
            others = set(tconds.loc[(tconds["condition_id"] == e["condition_id"])
                                    & (tconds["matched_term"] != e["term"]), "nct_id"])
            row["n_literal"] = int(mine["literal_match"].sum())
            row["n_only_this_term"] = int(len(set(mine["nct_id"]) - others))
        term_rows.append(row)
    return {
        "api_version": meta["api_version"], "data_timestamp": meta["data_timestamp"],
        "version_retrieved_at": meta["version_retrieved_at"],
        "records_retrieved_min": str(trials["retrieved_at"].min()), "records_retrieved_max": str(trials["retrieved_at"].max()),
        "n_queries_run": sum(1 for e in log if not e["skipped"]), "n_terms_skipped": sum(1 for e in log if e["skipped"]),
        "query_log": term_rows,
        "expansion": expansion,
        "expansion_named": named_counts,
        "record_coverage": coverage,
        "n_trials": int(trials["nct_id"].nunique()), "n_trial_condition_rows": len(tconds),
        "n_trial_condition_pairs": int(tconds[["nct_id", "condition_id"]].drop_duplicates().shape[0]),
        "n_trials_multi_condition": int((trials["n_matched_conditions"] > 1).sum()),
        "n_trials_no_literal_match": int((~trials["any_literal_match"].astype(bool)).sum()),
        "per_condition": per_cond.to_dict("records"), "per_status": status.to_dict("records"),
        "condition_by_status": cond_status, "per_study_type": stype.to_dict("records"),
        "n_interventions": len(iv), "intervention_types": iv["intervention_type"].value_counts(dropna=False).to_dict(),
        "n_outcomes": len(oc), "outcome_types": oc["outcome_type"].value_counts().to_dict(),
        "n_sites": len(sites), "n_trials_with_sites": int(sites["nct_id"].nunique()),
        "n_trials_without_sites": int((trials["n_locations"] == 0).sum()),
        "country_groups": sites["country_group"].value_counts().to_dict(),
        "n_countries": int(sites["country"].nunique()),
        "top_countries": sites["country"].value_counts().head(10).to_dict(),
        "n_us_sites": len(us), "n_us_sites_geopoint": len(us_pt),
        "n_us_sites_zip_only_geocode": int((us["geocode_method"] == "zip5_as_zcta_internal_point").sum()),
        "n_us_sites_not_geocoded": int((us["geocode_method"] == "not_geocoded").sum()),
        "n_sites_no_geopoint": int(sites["geopoint_lat"].isna().sum()),
        "n_sites_no_coordinates_any": int(sites["lat"].isna().sum()),
        "n_us_sites_with_county": int(us["county_fips"].notna().sum()),
        "n_us_sites_county_pip": int((us["county_assignment_method"] == "geopoint_point_in_polygon_cb_2024_5m").sum()),
        "n_us_zip_state_conflict": int(us["zip_state_conflict"].sum()),
        "n_us_county_pip_majority": int((us["county_assignment_method"] == "geopoint_county_majority_of_zcta").sum()),
        # CT excluded: the 2020 relationship file has legacy CT counties, so a 2024 CT planning-region
        # geoPoint county never overlaps a ZCTA there (vintage mismatch, not evidence of non-overlap)
        "n_us_disagree_pip_outside_zcta": int(((us["county_methods_agree"] == False)  # noqa: E712
                                               & us["pip_county_share_of_zcta"].isna() & ~us_ct).sum()),
        "n_us_disagree_pip_overlaps_zcta": int(((us["county_methods_agree"] == False)  # noqa: E712
                                                & us["pip_county_share_of_zcta"].notna()).sum()),
        "county_disagreement_top": (us[us["county_methods_agree"] == False]  # noqa: E712
                                    .groupby(["city", "state", "county_fips_pip", "county_fips_zip"]).size()
                                    .sort_values(ascending=False).head(8).rename("n_sites").reset_index()
                                    .to_dict("records")),
        "n_us_sites_county_zip": int((us["county_assignment_method"] == "zip5_as_zcta_internal_point").sum()),
        "n_ct_sites": int(us_ct.sum()),
        "n_ct_sites_with_county": int(us.loc[us_ct, "county_fips"].notna().sum()),
        "n_ct_sites_2024_region": int(us.loc[us_ct, "county_fips"].isin(ct_2024).sum()),
        "n_ct_sites_legacy_county": int(us.loc[us_ct, "county_fips"].isin(ct_legacy).sum()),
        "n_ct_disagree": int((us_ct & (us["county_methods_agree"] == False)).sum()),  # noqa: E712
        "n_ct_share_available": int(us.loc[us_ct, "pip_county_share_of_zcta"].notna().sum()),
        "n_us_geopoint_outside_polygons": int((us_pt["county_fips_pip"].isna()).sum()),
        "n_us_both_county_methods": int(both.sum()),
        "n_us_county_methods_agree": int(us.loc[both, "county_methods_agree"].astype(bool).sum()),
        "n_us_pip_state_mismatch": int((us["pip_state_matches_reported"] == False).sum()),  # noqa: E712
        "n_us_state_unmapped": int(us["reported_state_fips"].isna().sum()),
        "n_us_counties_with_sites": int(us["county_fips"].nunique()),
        "n_territory_sites": int((sites["country_group"] == "US_territory").sum()),
        "n_us_distinct_geopoints": int(len(fac_per_coord)),
        "n_us_distinct_city_state_geopoint_sites": int(cityset["_c"].nunique()),
        "n_us_geopoints_shared_by_multiple_facilities": int((fac_per_coord > 1).sum()),
        "n_us_sites_at_shared_geopoint": int(sum(
            fac_per_coord.to_dict()[k] > 1 for k in zip(us_pt["geopoint_lat"], us_pt["geopoint_lon"]))),
        "n_sites_generic_facility": int(sites["facility_generic"].sum()),
        "n_sites_missing_facility": int(sites["facility_missing"].sum()),
        "n_sites_identifiable_facility": int(sites["facility_identifiable"].sum()),
        "n_research_sites": len(hist), "n_research_sites_us": int((hist["country_group"] == "US").sum()),
        "n_research_sites_unnamed_bucket": int(hist["facility_missing"].sum()),
        "n_research_sites_generic": int(hist["facility_generic"].sum()),
        "n_research_sites_identifiable": int(hist["facility_identifiable"].sum()),
        "n_research_sites_identifiable_us": int((hist["facility_identifiable"] & (hist["country_group"] == "US")).sum()),
        "miss_trials": _miss(trials, ["brief_summary", "official_title", "start_date", "primary_completion_date",
                                      "completion_date", "phases", "allocation", "enrollment_count",
                                      "lead_sponsor_class", "conditions", "keywords", "condition_mesh_terms",
                                      "eligibility_criteria", "sex", "minimum_age", "maximum_age",
                                      "intervention_names", "last_update_post_date"]),
        "miss_trials_interventional": _miss(trials, ["phases", "allocation", "intervention_names"],
                                            trials["study_type"] == "INTERVENTIONAL"),
        "miss_sites": _miss(sites, ["facility", "city", "state", "zip", "country", "site_status", "geopoint_lat"]),
        "miss_us_sites": _miss(sites, ["facility", "city", "state", "zip", "geopoint_lat", "county_fips"],
                               sites["country_group"] == "US"),
        "miss_outcomes": _miss(oc, ["measure", "description", "time_frame"]),
        "miss_interventions": _miss(iv, ["intervention_type", "intervention_name", "intervention_description",
                                         "other_names"]),
        "row_counts": {"clinical_trials": len(trials), "trial_conditions": len(tconds),
                       "trial_interventions": len(iv), "trial_outcomes": len(oc), "trial_sites": len(sites),
                       "ctgov_facility_summary": len(hist)},
        "far_future_start": trials.loc[trials["start_year"] > 2035, ["nct_id", "start_date", "overall_status"]]
                                  .to_dict("records"),
        "snapshot_date": snapshot_date,
        "n_start_after_today": int((trials["start_date"].fillna("") > snapshot_date).sum()),
        "n_start_missing": int(trials["start_date"].isna().sum()),
        "start_year_min": int(trials["start_year"].min()) if trials["start_year"].notna().any() else None,
        "start_year_max": int(trials["start_year"].max()) if trials["start_year"].notna().any() else None,
    }


def _md_table(records: list[dict], cols: list[str] | None = None) -> str:
    if not records:
        return "_none_\n"
    cols = cols or list(records[0].keys())
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in records:
        lines.append("| " + " | ".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + " |")
    return "\n".join(lines) + "\n"


def _pct(a: int, b: int) -> str:
    return f"{a:,} ({100 * a / b:.1f}%)" if b else f"{a:,}"


def write_audit(st: dict) -> Path:
    manifest = load_manifest(SOURCE_ID)["files"]
    n_q = sum(1 for k in manifest if k.startswith("queries/"))
    n_r = sum(1 for k in manifest if k.startswith("records/"))
    bytes_r = sum(v["bytes"] for k, v in manifest.items() if k.startswith("records/"))
    cbs = st["condition_by_status"]
    cbs_records = cbs.reset_index().to_dict("records") if len(cbs) else []
    us = st["n_us_sites"]
    qlog = [dict(r, query_string=f"`{r['query_string']}`" if r["query_string"] else "—") for r in st["query_log"]]
    text = f"""# DATA AUDIT — ClinicalTrials.gov (API v2)

| Field | Value |
|---|---|
| source_id | {SOURCE_ID} |
| Source (dataset/API name, exact files/endpoints) | ClinicalTrials.gov API v2: `GET {STUDIES_URL}` (query.cond per search term; filter.ids for full records), `GET {VERSION_URL}`, `GET {SEARCH_AREAS_URL}`, `GET {OAS_URL}` |
| Publishing organization | U.S. National Library of Medicine (NLM), National Institutes of Health |
| Retrieval date (UTC) | version endpoint {st['version_retrieved_at']}; study records {st['records_retrieved_min']} to {st['records_retrieved_max']} |
| Source version / release | API {st['api_version']}; dataTimestamp **{st['data_timestamp']}** (UTC, from /api/v2/version) |
| Source update date / cadence | dataTimestamp {st['data_timestamp']}; the registry is continuously updated by sponsors, so counts drift between runs (re-run with `--refresh` to update) |
| License / access conditions | Open API, no key or registration. The ClinicalTrials.gov terms page (https://clinicaltrials.gov/about-site/terms-conditions) is JavaScript-rendered and its text could not be retrieved programmatically on 2026-09-23. NLM's data terms (https://www.nlm.nih.gov/databases/download/terms_and_conditions.html) ask users to credit "Courtesy of the U.S. National Library of Medicine", not imply NLM endorsement, and disclose when redistributed data are not current; that page does not name ClinicalTrials.gov explicitly. |
| Unit of observation | Registered study (NCT ID); long tables per condition match, intervention, outcome measure and site |
| Sample size (actual, as ingested) | {st['n_trials']:,} unique studies; {st['n_trial_condition_pairs']:,} study x condition pairs; {st['n_interventions']:,} interventions; {st['n_outcomes']:,} outcome measures; {st['n_sites']:,} study x site rows ({us:,} U.S.) |
| Geography (resolution, vintage) | Site-level: ClinicalTrials.gov geoPoint (measured city-level, see below) or ZIP->ZCTA internal point; county = 2024 cartographic boundaries (cb_2024_us_county_5m) |
| Person-level? | no (registry of studies; no participant data) |
| Geographic? | yes (site locations only; no patient locations) |
| Omics? | no (omics appear only as registered outcome text) |
| Wearable? | no (wearables appear only as registered intervention/outcome text) |
| True participant linkage across modalities? | no — not applicable; trial records hold no participant-level data |

> **Guardrail.** Presence of a measurement in a registered trial (as an intervention, device or outcome) is
> **not** evidence that the measurement works, that a device diagnoses the target illness, or that an
> intervention is effective. Registration is a research-activity and deployment-interest signal only.

## Files / endpoints retrieved
All files are under `data/raw/{SOURCE_ID}/`; `MANIFEST.json` holds URL, bytes, sha256 and retrieval time for each.

* `api_version.json` — `{VERSION_URL}`
* `api_search_areas.json` — `{SEARCH_AREAS_URL}` (defines which fields query.cond searches)
* `api_oas_v2.yaml` — `{OAS_URL}` (OpenAPI spec: pageSize max 1000, nextPageToken paging)
* `queries/<condition>__<term>__pNNN.json` — {n_q} pages. Parameters: `query.cond="<term>"`, `fields=NCTId`,
  `pageSize=1000`, `countTotal=true`, then `pageToken=<nextPageToken>`. No status, date or type filter.
* `records/batch_<sha1>.json` — {n_r} pages, {bytes_r / 1e6:.1f} MB. Parameters: `filter.ids=<up to {ID_BATCH_SIZE} NCT IDs>`,
  `fields={RECORD_FIELDS}`, `pageSize=1000`.
* `query_log.json` — every query string sent, its totalCount, retrieved count and raw file names.

Record coverage: requested {st['record_coverage']['requested']:,} NCT IDs, returned {st['record_coverage']['returned']:,};
missing {len(st['record_coverage']['missing'])}; unexpected extras {len(st['record_coverage']['unexpected_extra'])}.

### Query strings (one query per configured search term; `query.cond`)
Terms are phrase-quoted. An unquoted query ANDs the words anywhere in conditions/titles/keywords (e.g. pulling in
"long-term care" COVID-vaccine trials for `long COVID`). Quoted vs unquoted, measured in this run (unquoted pages saved
under `queries_unquoted_check/`, used only for this comparison):

{_md_table(st.get('quoting_check', []), ['term', 'quoted_query', 'n_quoted', 'n_unquoted', 'n_unquoted_not_in_quoted', 'n_quoted_not_in_unquoted', 'n_unquoted_only_outside_condition', 'n_of_those_naming_condition_term'])}
**Quoting costs some recall.** `n_unquoted_only_outside_condition` = unquoted hits that no quoted term of the same
condition retrieved; `n_of_those_naming_condition_term` = how many of them nevertheless carry one of the condition's
searched terms literally in conditions/keywords/titles/acronym (records fetched and saved under
`queries_unquoted_check/*__outside_*`). These trials are **not** in the tables. The count is a lower bound: registrations
such as "COVID-19 Long-Hauler Syndrome" name the condition without containing a configured term. Examples:
{'; '.join(dict.fromkeys(e for r in st.get('quoting_check', []) for e in (r.get('examples_lost_by_quoting') or '').split('; ') if e)) or 'none'}.
Quoting was kept because the unquoted precision loss is much larger (see `n_unquoted` vs `n_quoted`); treat the
trial sets as a high-precision sample, not a complete census.

Bare acronyms in configs/conditions.yaml are never searched.

{_md_table(qlog, ["condition_id", "term", "query_string", "total_count", "n_unique_ids", "n_pages", "n_literal", "n_only_this_term", "skipped", "skip_reason"])}
`n_literal` = hits whose record contains the term itself in conditions, keywords, titles or acronym;
`n_only_this_term` = hits that no other term of the same condition returned. Identical `total_count` across a
condition's terms (e.g. all three dysautonomia phrases) shows the API mapped them to the same concept.
## Key variables
| Table.column | Meaning | Downstream use |
|---|---|---|
| clinical_trials.nct_id | ClinicalTrials.gov identifier (primary key) | joins every trial_* table |
| clinical_trials.overall_status / start_date / primary_completion_date / completion_date / last_update_post_date | registry status and dates (strings as registered: YYYY-MM or YYYY-MM-DD; *_type = ACTUAL/ESTIMATED) | recency, recruiting activity |
| clinical_trials.study_type / phases / allocation / primary_purpose | design | trial maturity |
| clinical_trials.enrollment_count (+ enrollment_type) | planned (ESTIMATED) or actual enrollment | scale |
| clinical_trials.lead_sponsor / lead_sponsor_class / collaborators | sponsor (class: INDUSTRY, NIH, FED, OTHER, ...) | research readiness |
| clinical_trials.conditions / keywords / condition_mesh_terms / condition_mesh_ancestors | registered conditions and NLM-derived MeSH | ontology linking |
| clinical_trials.eligibility_criteria / sex / minimum_age / maximum_age | eligibility (free text + structured) | population description only |
| clinical_trials.arm_groups_json | arm label, type, description, intervention names | design detail |
| trial_conditions.condition_id / matched_term / query_string | which configured condition and term matched | condition filter |
| trial_conditions.match_fields / literal_match | where this row's term literally occurs (condition, keyword, title, acronym, MeSH term, MeSH ancestor); literal_match=False = this term matched only through API synonym/MeSH expansion | per-term diagnostics |
| trial_conditions.condition_literal_match | any of the condition's searched terms literally in conditions/keywords/titles/acronym | precision filter (literal_only) |
| trial_interventions.intervention_type / name / description / other_names | registered interventions (DEVICE, DIAGNOSTIC_TEST, DRUG, BEHAVIORAL, ...) | measurement/technology discovery (text mining) |
| trial_outcomes.outcome_type / measure / description / time_frame | registered primary/secondary/other outcomes | objective-endpoint discovery (text mining) |
| trial_sites.facility / city / state / zip / country / site_status | registered location (contacts dropped deliberately) | research-site layer |
| trial_sites.lat / lon / geocode_method / source_geographic_resolution | geoPoint (`city_centroid`, geocode_method `ctgov_geopoint_city_level`) or ZCTA internal point (`zcta_centroid`) | maps (city-level precision) |
| trial_sites.county_fips / county_assignment_method | 2024 county: the site ZIP's county when it resolves in the registered state (unless the geoPoint county holds >=50% of the ZCTA's land), else point-in-polygon of geoPoint | ecological joins (site counts per county) |
| trial_sites.county_fips_zip / county_fips_pip / county_methods_agree / zip_state_conflict / pip_county_share_of_zcta | both county candidates and the evidence used to choose | QA |
| trial_sites.facility_generic / facility_missing / facility_identifiable | sponsor placeholder or non-facility names ("Research Site", "Local Institution", "Site 101", "Site Reference ID/Investigator# 5870", "Contact <sponsor> for exact locations", virtual/remote/online sites, a name that is only the city); no facility name; neither | facility matching uses facility_identifiable |

## Sample sizes

### Trials per condition (unique NCT IDs; a trial can match several conditions)
{_md_table(st['per_condition'], ['condition_id', 'n_trials', 'n_trials_literal'])}
`n_trials_literal` = trials where at least one of the condition's searched terms appears literally in the record's
conditions, keywords, titles or acronym (`trial_conditions.condition_literal_match`; case, punctuation, hyphens and the
COVID-19/SARS-CoV-2 spelling variants are folded). The rest are matched only through ClinicalTrials.gov's synonym and
MeSH expansion; for `long_covid` most of them are registered under the MeSH term "Post-Acute COVID-19 Syndrome" and
are genuine long-COVID trials, whereas for `dysautonomia` and `mcas` most are other diseases (below). What that expansion adds, per condition (counts of trials; last column = the MeSH
descriptor or ancestor carried by >= half of the non-literal trials that is rarest across all ingested trials):

{_md_table(st['expansion'])}
Read this before using a condition's trial set. For `dysautonomia`, {st['expansion_named'].get('dysautonomia:Autonomic Nervous System Diseases')}
non-literal trials carry the MeSH ancestor "Autonomic Nervous System Diseases" (multiple system atrophy, orthostatic
hypotension, CRPS, Parkinson disease). For `mcas`, {st['expansion_named'].get('mcas:Mast Cell Activation Disorders')} non-literal trials carry
"Mast Cell Activation Disorders" and {st['expansion_named'].get('mcas:Mastocytosis')} carry "Mastocytosis": they are mastocytosis and
hematologic-malignancy trials, which configs/conditions.yaml says to treat as a related, distinct entity. Filter with
`trial_conditions.condition_literal_match` (or `trials_for_condition(..., literal_only=True)`) where precision matters.

* Unique trials: **{st['n_trials']:,}**; trials matching >1 condition: {st['n_trials_multi_condition']:,};
  trials with no literal match for any term: {st['n_trials_no_literal_match']:,}.
* Registered start years span {st['start_year_min']}–{st['start_year_max']}; {st['n_start_after_today']:,} trials have a
  start date after the snapshot date {st['snapshot_date']} (planned/estimated) and {st['n_start_missing']:,} have none. Start years after 2035 are
  registry placeholders: {st['far_future_start']}.

### Trials per overall status
{_md_table(st['per_status'], ['overall_status', 'n_trials'])}
### Trials per study type
{_md_table(st['per_study_type'], ['study_type', 'n_trials'])}
### Condition x status (unique trials)
{_md_table(cbs_records)}
### Interventions and outcomes
* Intervention rows: {st['n_interventions']:,}; by type: {st['intervention_types']}
* Outcome rows: {st['n_outcomes']:,}; by type: {st['outcome_types']}

### Sites
* Study x site rows: **{st['n_sites']:,}** across {st['n_countries']} countries; trials with at least one site:
  {st['n_trials_with_sites']:,}; trials with **no** registered site: {st['n_trials_without_sites']:,}.
* Country groups: {st['country_groups']}; top countries: {st['top_countries']}
* U.S. sites: **{us:,}**. With API geoPoint (resolution `city_centroid`): {_pct(st['n_us_sites_geopoint'], us)}; no geoPoint but
  ZIP-geocoded to ZCTA centroid: {_pct(st['n_us_sites_zip_only_geocode'], us)}; not geocoded: {_pct(st['n_us_sites_not_geocoded'], us)}.
* Sites lacking an API geoPoint (all countries): {_pct(st['n_sites_no_geopoint'], st['n_sites'])}; lacking any coordinates after
  ZIP fallback: {_pct(st['n_sites_no_coordinates_any'], st['n_sites'])} (non-U.S. sites get no ZIP fallback).
* U.S. sites with a 2024 county FIPS: {_pct(st['n_us_sites_with_county'], us)} (from the site ZIP {st['n_us_sites_county_zip']:,};
  geoPoint county holding >=50% of the ZIP's land {st['n_us_county_pip_majority']:,}; geoPoint point-in-polygon because the ZIP
  was unusable {st['n_us_sites_county_pip']:,}); ZIPs resolving to a county outside the registered state
  (treated as typos, geoPoint used): {st['n_us_zip_state_conflict']:,}; U.S. geoPoints outside every county polygon:
  {st['n_us_geopoint_outside_polygons']:,}; distinct counties with >=1 U.S. site: {st['n_us_counties_with_sites']:,}.
* Point-in-polygon state vs registered state disagree: {st['n_us_pip_state_mismatch']:,} U.S. sites; U.S. sites whose
  state name did not map to a FIPS code: {st['n_us_state_unmapped']:,}. U.S.-territory sites (listed by
  ClinicalTrials.gov as separate countries): {st['n_territory_sites']:,}.
* **geoPoint precision (measured).** Among {st['n_us_sites_geopoint']:,} U.S. geoPoint sites there are
  {st['n_us_distinct_geopoints']:,} distinct coordinates for {st['n_us_distinct_city_state_geopoint_sites']:,} distinct
  city/state pairs; {st['n_us_geopoints_shared_by_multiple_facilities']:,} coordinates are shared by more than one
  distinct facility name, covering {st['n_us_sites_at_shared_geopoint']:,} site rows. geoPoint is ClinicalTrials.gov's
  own geocode of the city, not the facility address (e.g. "860 Fifth Avenue, New York" is placed at the New York City
  centroid 40.71427, -74.00597). Where both methods apply, the geoPoint point-in-polygon county and the ZIP county
  agree for {_pct(st['n_us_county_methods_agree'], st['n_us_both_county_methods'])} U.S. sites. The disagreements
  concentrate where a postal city name covers several counties or an independent city; largest groups:

{_md_table(st['county_disagreement_top'])}
  In these groups the ZIP identifies the facility's own postal area (e.g. Emory University, ZIP 30322 -> DeKalb 13089;
  University of Colorado Anschutz, ZIP 80045 -> Adams 08001), whereas geoPoint places every "Atlanta"/"Aurora" site
  at the core-city centroid (Fulton 13121 / Arapahoe 08005). Of the disagreeing sites, {st['n_us_disagree_pip_outside_zcta']:,}
  have a geoPoint county that does not overlap the site's ZCTA at all (2020 ZCTA-county relationship file), so the
  city centroid is certainly not the facility's county; in the other {st['n_us_disagree_pip_overlaps_zcta']:,} the ZCTA straddles both
  counties. **county_fips therefore uses the ZIP county first**, except when the geoPoint county holds >=50% of the
  ZCTA's land (ZCTA internal points can fall in a minority county, e.g. ZCTA 32174 Ormond Beach is 70.6% Volusia but its
  internal point lies in Flagler), and falls back to point-in-polygon when the ZIP does not resolve (PO-box/unique ZIPs,
  cross-state typos, no ZIP). Land share is area, not population, and uses 2020 counties (CT regions never match).
* **Connecticut (county vintage).** {st['n_ct_sites']:,} U.S. sites are in Connecticut; {st['n_ct_sites_with_county']:,} have a
  county_fips, of which {st['n_ct_sites_2024_region']:,} are 2024 planning regions (09110-09190) and
  {st['n_ct_sites_legacy_county']:,} are legacy counties (09001-09015): both the ZIP route (zcta_centroids, 2024 polygons)
  and point-in-polygon (cb_2024) yield planning regions, and nothing is mapped between vintages. The 2020 ZCTA-county
  relationship file uses legacy CT counties, so `pip_county_share_of_zcta` is null for CT (non-null for
  {st['n_ct_share_available']:,} CT sites) and the >=50% majority rule never applies there; the {st['n_ct_disagree']:,} CT
  ZIP-vs-geoPoint disagreements keep the ZIP county and are excluded from the "does not overlap" count above.
* Facility names that are sponsor placeholders (`facility_generic`): {_pct(st['n_sites_generic_facility'], st['n_sites'])};
  facility missing (`facility_missing`): {st['n_sites_missing_facility']:,}; named, non-placeholder facility
  (`facility_identifiable`): {_pct(st['n_sites_identifiable_facility'], st['n_sites'])}.
* Per-facility site history (`ctgov_facility_summary`): {st['n_research_sites']:,} keys
  ({st['n_research_sites_us']:,} U.S.); {st['n_research_sites_identifiable']:,} are identifiable facilities
  ({st['n_research_sites_identifiable_us']:,} U.S.), {st['n_research_sites_generic']:,} are sponsor-placeholder names and
  {st['n_research_sites_unnamed_bucket']:,} are **unnamed-site buckets** (all sites registered without a facility name in
  one city share a key, so their n_trials is a count per city, not per facility). Use `facility_identifiable` for
  facility-level matching.

## Missingness
Measured on the processed tables (missing = null or empty string).

### clinical_trials (all {st['n_trials']:,} trials)
{_md_table(st['miss_trials'])}
Interventional trials only (phase/allocation are not applicable to observational studies):
{_md_table(st['miss_trials_interventional'])}
### trial_sites (all rows)
{_md_table(st['miss_sites'])}
`site_status` is only populated by ClinicalTrials.gov for studies whose recruitment is ongoing.

### trial_sites (U.S. rows)
{_md_table(st['miss_us_sites'])}
### trial_outcomes
{_md_table(st['miss_outcomes'])}
### trial_interventions
{_md_table(st['miss_interventions'])}
## Linkage strategy
* `nct_id` joins clinical_trials, trial_conditions, trial_interventions, trial_outcomes and trial_sites.
* `trial_conditions.condition_id` = `configs/conditions.yaml` ids (the ontology layer's canonical_condition_id).
  `condition_mesh_ids` give MeSH descriptor IDs for ontology cross-checks.
* `trial_sites.county_fips` / `state_fips` join to `geographies` (2024 vintage) **ecologically only**: a site in a
  county says where research infrastructure is, not where participants live. Nothing here locates any patient.
* Trial sites and NIH RePORTER / NPPES / HRSA facilities have no shared identifier; matching facilities across
  sources needs name + location fuzzy matching and must be reported as approximate.
* Measurement/technology discovery is done by text-mining interventions and outcomes (helper
  `trials_for_condition(condition_id, measurement_pattern)` accepts a `configs/measurements.yaml` class id or a regex).
  A match means the trial registered that measurement, not that it works.

## Limitations and caveats
* **Registration != evidence.** Trials describe planned work; status, enrollment and dates are sponsor-reported and
  many are ESTIMATED. Completed trials may have unpublished or null results.
* **Search recall/precision.** query.cond matches synonyms and MeSH ancestors; see `literal_match`. Conversely, trials
  that describe a condition only in the summary (not conditions/titles/keywords) are not found. Acronym-only
  registrations (e.g. conditions listed only as "PASC" or "ME/CFS") are found only if the API's synonym expansion
  maps the phrase, since bare acronyms are never searched.
* **Condition grouping.** `post_infectious_syndrome` is a grouping; `dysautonomia` and `lyme_disease` are broad and
  include parent/child MeSH matches. A trial can match several conditions; counts per condition are not additive.
* **Geography.** geoPoint is city-level; lat/lon are therefore city centroids, labelled
  `source_geographic_resolution = city_centroid` (labelled `point` before 2026-09-23, when the shared vocabulary had
  no city-centroid level). County comes from the ZIP where possible; a ZIP->ZCTA
  mapping is itself approximate (ZIP != ZCTA; ZCTAs can straddle counties; the ZCTA internal point decides). Non-U.S. sites are kept with their geoPoint and no county. U.S. territories are listed as
  their own countries by ClinicalTrials.gov and are flagged `US_territory`.
* **Facility identity.** One institution appears under many spellings; sponsor placeholders ("Research Site",
  "Site Reference ID/Investigator# NNNNN", "Contact <sponsor> for exact locations"), virtual/remote/online sites and
  names that are only the city hide the facility (all flagged `facility_generic`). Per-facility counts are therefore approximate lower bounds per spelling.
* **Snapshot.** Counts reflect dataTimestamp {st['data_timestamp']}; the registry changes daily.

## Processed outputs
{_md_table([{'table': k, 'rows': f'{v:,}'} for k, v in st['row_counts'].items()])}
All tables carry the provenance columns (`data_layer=facility`, `evidence_type=trial_registry_record`,
`source_version` = API version + dataTimestamp). Site rows carry per-row `source_geographic_resolution`
(`city_centroid`, `zcta_centroid`, `state`, `none`).

## Reproduce
`uv run python -m {MODULE}` (add `--refresh` to re-query the live API; otherwise saved raw pages are reused).
Tests: `uv run pytest tests/test_clinicaltrials.py`.
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


def write_registry(st: dict) -> Path:
    return write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "ClinicalTrials.gov API v2",
        "publisher": "U.S. National Library of Medicine (NLM), NIH",
        "landing_url": LANDING_URL,
        "access_urls": [STUDIES_URL, VERSION_URL, SEARCH_AREAS_URL, OAS_URL],
        "license": ("Open U.S. government registry; NLM data terms ask for the credit line 'Courtesy of the U.S. "
                    "National Library of Medicine', no implied endorsement, and disclosure when redistributed data "
                    "are not current (https://www.nlm.nih.gov/databases/download/terms_and_conditions.html). The "
                    "ClinicalTrials.gov terms page is JavaScript-rendered and was not machine-read."),
        "access_conditions": "open API, no key or registration",
        "retrieved_at": st["records_retrieved_max"],
        "source_version": f"API {st['api_version']}; dataTimestamp {st['data_timestamp']} UTC",
        "update_date": st["data_timestamp"],
        "data_layer": "facility",
        "unit_of_observation": "registered study (NCT ID); long tables per condition match, intervention, outcome, site",
        "sample_size": {
            "trials": st["n_trials"], "trial_condition_pairs": st["n_trial_condition_pairs"],
            "interventions": st["n_interventions"], "outcomes": st["n_outcomes"], "sites": st["n_sites"],
            "us_sites": st["n_us_sites"], "us_sites_with_geopoint": st["n_us_sites_geopoint"],
            "us_sites_with_county": st["n_us_sites_with_county"], "research_site_keys": st["n_research_sites"],
            "queries_run": st["n_queries_run"], "terms_skipped_bare_acronym": st["n_terms_skipped"],
            "trials_per_condition": {r["condition_id"]: int(r["n_trials"]) for r in st["per_condition"]},
        },
        "geographic_resolution": ("site point (ClinicalTrials.gov geoPoint, city-level) or ZCTA centroid; county_fips "
                                  "(2024) from the site ZIP's ZCTA first, else geoPoint point-in-polygon"),
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (no participant data)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": list(st["row_counts"].keys()),
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": MODULE,
        "limitations": [
            "Registration is not evidence that a measurement works or an intervention is effective",
            "query.cond expands synonyms and MeSH ancestors; use trial_conditions.condition_literal_match for precision",
            "Phrase-quoted search misses some registrations that name the condition (e.g. 'Long COVID-19'); "
            "the number is measured in DATA_AUDIT.md (quoting check) -- the trial sets are not a complete census",
            "ctgov_facility_summary rows with facility_missing are per-city buckets of unnamed sites, not facilities",
            "Bare acronyms (PASC, ME/CFS) are never searched",
            "geoPoint is a city-level geocode, not the facility address",
            "Facility names vary in spelling; sponsor placeholder names are flagged facility_generic",
            f"Snapshot at dataTimestamp {st['data_timestamp']}",
        ],
    })


# --------------------------------------------------------------------------------------------
# query helpers (for the MCP layer)
# --------------------------------------------------------------------------------------------

@lru_cache(maxsize=4)
def _tables(stamp: tuple) -> dict[str, pd.DataFrame]:
    return {n: read_table(n) for n in ("clinical_trials", "trial_conditions", "trial_interventions",
                                       "trial_outcomes", "trial_sites")}


def _load_tables() -> dict[str, pd.DataFrame]:
    names = ("clinical_trials", "trial_conditions", "trial_interventions", "trial_outcomes", "trial_sites")
    stamp = tuple(processed_path(n).stat().st_mtime if processed_path(n).exists() else 0 for n in names)
    return _tables(stamp)


def measurement_regex(measurement_pattern: str) -> tuple[re.Pattern, re.Pattern | None, str]:
    """Resolve a measurements.yaml class id to its include/exclude regexes, else treat the input as a regex."""
    for m in load_config("measurements").get("measurement_classes", []):
        if m["id"] == measurement_pattern:
            inc = re.compile("|".join(f"(?:{p})" for p in m["patterns"]), re.IGNORECASE)
            exc = (re.compile("|".join(f"(?:{p})" for p in m["exclude_patterns"]), re.IGNORECASE)
                   if m.get("exclude_patterns") else None)
            return inc, exc, f"measurement class {m['id']}"
    return re.compile(measurement_pattern, re.IGNORECASE), None, "regex"


def _text_hits(texts: pd.Series, inc: re.Pattern, exc: re.Pattern | None) -> pd.Series:
    def hit(t):
        if t is None or (isinstance(t, float) and np.isnan(t)):
            return None
        t = str(t)
        if exc is not None:
            t = exc.sub(" ", t)
        m = inc.search(t)
        return m.group(0) if m else None
    return texts.map(hit)


def trials_for_condition(condition_id: str, measurement_pattern: str | None = None, *,
                         literal_only: bool = False, statuses: list[str] | None = None,
                         us_sites_only: bool = False) -> pd.DataFrame:
    """Trials registered for a target condition, optionally those that mention a measurement.

    condition_id: an id from configs/conditions.yaml (e.g. 'long_covid').
    measurement_pattern: a measurement class id from configs/measurements.yaml (its patterns and
        exclude_patterns are applied) or a case-insensitive regular expression. Searched fields: titles,
        brief summary, keywords, intervention name/description/other names, outcome measure/description.
        Eligibility criteria are not searched (they mention tests as exclusions).
    literal_only: keep only trials in which one of the condition's searched terms appears literally in
        conditions/keywords/titles/acronym (trial_conditions.condition_literal_match).
    statuses: optional list of overall_status values (e.g. ['RECRUITING']).

    Returns one row per trial with `measurement_match_fields` / `measurement_match_example` when a
    pattern is given. Presence in a trial is NOT evidence that the measurement works.
    """
    known = {c["id"] for c in load_config("conditions")["conditions"]}
    if condition_id not in known:
        raise ValueError(f"unknown condition_id {condition_id!r}; expected one of {sorted(known)}")
    t = _load_tables()
    tc = t["trial_conditions"]
    tc = tc[tc["condition_id"] == condition_id]
    lit_col = "condition_literal_match" if "condition_literal_match" in tc.columns else "literal_match"
    if literal_only:
        tc = tc[tc[lit_col].astype(bool)]
    per = tc.groupby("nct_id").agg(condition_terms=("matched_term", lambda x: join_list(sorted(set(x)))),
                                   condition_literal_match=(lit_col, "any"))
    cols = ["nct_id", "brief_title", "official_title", "overall_status", "study_type", "phases", "start_date",
            "primary_completion_date", "completion_date", "last_update_post_date", "enrollment_count",
            "lead_sponsor", "lead_sponsor_class", "conditions", "keywords", "intervention_types",
            "intervention_names", "n_locations", "n_us_locations", "has_results", "study_url",
            "source_version", "retrieved_at"]
    df = t["clinical_trials"][cols + ["brief_summary"]].merge(per, left_on="nct_id", right_index=True)
    if statuses:
        df = df[df["overall_status"].isin(statuses)]
    if us_sites_only:
        df = df[df["n_us_locations"] > 0]
    if measurement_pattern:
        inc, exc, kind = measurement_regex(measurement_pattern)
        ids = set(df["nct_id"])
        fields: dict[str, pd.Series] = {
            "brief_title": df.set_index("nct_id")["brief_title"],
            "official_title": df.set_index("nct_id")["official_title"],
            "brief_summary": df.set_index("nct_id")["brief_summary"],
            "keywords": df.set_index("nct_id")["keywords"],
        }
        iv = t["trial_interventions"]
        iv = iv[iv["nct_id"].isin(ids)]
        oc = t["trial_outcomes"]
        oc = oc[oc["nct_id"].isin(ids)]
        long = [pd.DataFrame({"nct_id": s.index, "field": f, "hit": _text_hits(s, inc, exc).values})
                for f, s in fields.items()]
        for f, frame, col in (("intervention_name", iv, "intervention_name"),
                              ("intervention_description", iv, "intervention_description"),
                              ("intervention_other_names", iv, "other_names"),
                              ("outcome_measure", oc, "measure"), ("outcome_description", oc, "description")):
            long.append(pd.DataFrame({"nct_id": frame["nct_id"].values, "field": f,
                                      "hit": _text_hits(frame[col], inc, exc).values}))
        hits = pd.concat(long, ignore_index=True).dropna(subset=["hit"])
        agg = hits.groupby("nct_id").agg(measurement_match_fields=("field", lambda x: join_list(sorted(set(x)))),
                                         measurement_match_example=("hit", "first"))
        df = df.merge(agg, left_on="nct_id", right_index=True)
        df["measurement_query"] = f"{kind}: {measurement_pattern}"
    df = df.drop(columns=["brief_summary"]).sort_values("last_update_post_date", ascending=False)
    df = df.reset_index(drop=True)
    df.attrs["guardrail"] = GUARDRAIL_NOTE
    df.attrs["condition_id"] = condition_id
    return df


def trial_sites_for_condition(condition_id: str, measurement_pattern: str | None = None, *,
                              us_only: bool = True, literal_only: bool = False,
                              statuses: list[str] | None = None, exclude_generic: bool = False) -> pd.DataFrame:
    """Sites of the trials returned by trials_for_condition (same arguments), with county/state FIPS.

    exclude_generic: drop sites that do not identify a facility (sponsor placeholder names such as
    "Research Site" and sites registered without a facility name), i.e. keep facility_identifiable.
    """
    trials = trials_for_condition(condition_id, measurement_pattern, literal_only=literal_only, statuses=statuses)
    s = _load_tables()["trial_sites"]
    s = s[s["nct_id"].isin(set(trials["nct_id"]))]
    if us_only:
        s = s[s["country_group"] == "US"]
    if exclude_generic:
        s = s[s["facility_identifiable"].astype(bool)] if "facility_identifiable" in s.columns else s[~s["facility_generic"]]
    keep = ["nct_id", "facility", "city", "state", "zip", "country", "country_group", "site_status", "lat", "lon",
            "geocode_method", "source_geographic_resolution", "county_fips", "county_assignment_method",
            "state_fips", "facility_generic", "facility_missing", "facility_identifiable"]
    out = s[keep].merge(trials[["nct_id", "brief_title", "overall_status", "study_type", "start_date"]],
                        on="nct_id", how="left").reset_index(drop=True)
    out.attrs["guardrail"] = GUARDRAIL_NOTE
    return out


def get_trial(nct_id: str) -> dict:
    """One trial with its conditions, interventions, outcomes and sites, or UNKNOWN when absent."""
    t = _load_tables()
    row = t["clinical_trials"][t["clinical_trials"]["nct_id"] == nct_id]
    if row.empty:
        return {"nct_id": nct_id, "status": UNKNOWN}
    rec = row.iloc[0].to_dict()
    for name in ("trial_conditions", "trial_interventions", "trial_outcomes", "trial_sites"):
        rec[name] = t[name][t[name]["nct_id"] == nct_id].to_dict("records")
    rec["guardrail"] = GUARDRAIL_NOTE
    return rec


def main() -> None:
    ap = argparse.ArgumentParser(description="Ingest ClinicalTrials.gov API v2 for the target conditions")
    ap.add_argument("--refresh", action="store_true", help="re-query the live API instead of reusing raw pages")
    args = ap.parse_args()
    run(refresh=args.refresh)


if __name__ == "__main__":
    main()
