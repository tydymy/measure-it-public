"""NIH RePORTER (API v2) grant search for the target invisible-illness conditions.

What this module does
---------------------
1. For every condition in configs/conditions.yaml and every one of its
   ``search_terms`` it runs one RePORTER ``advanced_text_search`` over project
   title, abstract and terms for fiscal years 2015-2026 (the exact criteria
   JSON of every request is saved in data/raw/nih_reporter/queries.json).
2. For each demo-cluster condition (long_covid, me_cfs, pots, dysautonomia) it
   runs the SPEC section D "measurement-augmented" queries:
   ``(<any condition search term>) AND <measurement term>`` for wearable,
   biomarker, diagnostic, sensor, imaging, autonomic, microvascular, capillary,
   digital health and point of care.
3. It writes three processed tables:

   nih_projects            one row per appl_id (a fiscal-year award record),
                           organization ZIP geocoded with crosswalk.zip_to_geo
   nih_project_conditions  appl_id x condition_id x matched condition query
   nih_project_query_hits  appl_id x augmented (condition AND measurement) query

Interpretation guardrails
-------------------------
* A grant is a *research-capability signal* for an organization. It is not
  patient burden, not clinical capacity, and not evidence that a measurement
  works (CONVENTIONS section 3.7).
* ``award_amount`` is the RePORTER total-cost award obligation of one
  fiscal-year application. It is not spending, and multi-year funded awards
  are obligated in one year.
* Multi-component awards (P01/U54/U19...) appear as a parent record plus
  subproject records with their own appl_ids. RePORTER documents the parent
  ``award_amount`` as including all its subprojects, and the parent amount
  always equals the sum of its IC fundings. The subproject amounts do not
  always sum to it, though: measured over the retrieved records, a few parents
  are *smaller* than their subprojects (see DATA_AUDIT.md). ``obligation_total``
  therefore counts the parent's amount whenever the parent is present, and the
  subprojects only when it is not. Totals always count each appl_id once,
  however many conditions or queries matched it.
* RePORTER's ``terms`` field is a thesaurus-expanded concept list: a phrase can
  match because a synonym concept was assigned even when the phrase is not in
  the title or abstract. Each condition link therefore records where the
  phrase was found locally (``match_tier``).

Reproduce: ``uv run python -m measure_it.ingestion.nih_reporter``
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import re
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import UNKNOWN, load_config, raw_dir, utc_now_iso
from ..download import _record as _manifest_record  # shared, locked MANIFEST writer
from ..download import download_file, sha256_file
from ..http import post
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, table_exists, write_table

SOURCE_ID = "nih_reporter"
SOURCE_NAME = "NIH RePORTER Project API v2 (POST /v2/projects/search)"
API_URL = "https://api.reporter.nih.gov/v2/projects/search"
LANDING_URL = "https://api.reporter.nih.gov/"
DOC_URLS = {
    "RePORTER_Project_API_V2_data_elements.pdf":
        "https://api.reporter.nih.gov/documents/Data%20Elements%20for%20RePORTER%20Project%20API_V2.pdf",
    "reporter_swagger_v2.json": "https://api.reporter.nih.gov/swagger/v2/swagger.json",
}
MODULE = "measure_it.ingestion.nih_reporter"

FISCAL_YEARS = list(range(2015, 2027))
SEARCH_FIELD = "projecttitle,abstracttext,terms"
OPERATOR = "advanced"
PAGE_LIMIT = 500          # API maximum
MAX_OFFSET = 14999        # API maximum offset
MAX_RETRIEVABLE = 15000   # offsets 0, 500, ..., 14500 -> 15,000 records per criteria set
SORT = {"sort_field": "appl_id", "sort_order": "asc"}  # deterministic paging

# SPEC section D measurement terms. The RePORTER text search does no stemming
# (measured 2026-09-23, FY2015-2026: "wearable" 5,902 hits, "wearables" 2,728,
# either 5,964), so each term is searched as an OR of its listed variants.
AUGMENT_TERMS: dict[str, list[str]] = {
    "wearable": ["wearable", "wearables"],
    "biomarker": ["biomarker", "biomarkers"],
    "diagnostic": ["diagnostic", "diagnostics"],
    "sensor": ["sensor", "sensors"],
    "imaging": ["imaging"],
    "autonomic": ["autonomic"],
    "microvascular": ["microvascular"],
    "capillary": ["capillary", "capillaries"],
    "digital health": ["digital health"],
    "point of care": ["point of care"],
}

# NIH RCDC spending categories whose names match a target condition. Names are
# the exact strings observed in spending_categories_desc of the retrieved
# records (see DATA_AUDIT.md); conditions without an RCDC category are absent.
# RCDC categorisation is NIH's own curated classification, used here only as
# an agreement check on the text search, never as a filter.
# Observed 2026-09-23 among the 294 distinct category names on retrieved records.
# No condition-specific category exists for dysautonomia, eds_hsd, mcas, ibs,
# gastroparesis or post_infectious_syndrome (only broad ones such as
# "Digestive Diseases"), so those are left unmapped rather than proxied.
# ptlds maps to its parent disease category "Lyme Disease" (broader).
RCDC_CATEGORY_FOR_CONDITION: dict[str, str] = {
    "long_covid": "Post-Acute Sequelae of SARS-CoV-2 infection (PASC) including Long COVID",
    "me_cfs": "Chronic Fatigue Syndrome (ME/CFS)",
    "pots": "Postural Orthostatic Tachycardia Syndrome",
    "fibromyalgia": "Fibromyalgia",
    "lyme_disease": "Lyme Disease",
    "ptlds": "Lyme Disease",
    "migraine": "Migraine",
    "endometriosis": "Endometriosis",
}

PROVENANCE_NOTE_PROJECTS = (
    "RePORTER fiscal-year award record (appl_id). award_amount = total-cost award obligation for that "
    "fiscal year, not expenditure; multi-project awards have a parent record plus subproject records, and "
    "totals count the parent's amount instead of its subprojects' (see obligation_total). "
    "lat/lon = Census ZCTA internal point of the organization ZIP (zcta_centroid), not an address; when the "
    "ZIP has no ZCTA, RePORTER's own organization point verified to lie in the org's state (point). "
    "reporter_lat/reporter_lon always keep RePORTER's coordinate. Research-capability signal "
    "only: not patient burden and not evidence that any measurement works."
)
PROVENANCE_NOTE_LINKS = (
    "Condition assignment = RePORTER advanced_text_search hit over title/abstract/terms (thesaurus-expanded "
    "terms field), not a curated disease classification. match_tier records where the phrase was found locally."
)


# --------------------------------------------------------------------------- query building
@dataclass
class QuerySpec:
    query_id: str
    query_type: str                 # "condition_term" | "augmented"
    condition_id: str
    search_term: str                # condition search term, or the SPEC augment term
    search_text: str                # exact advanced_text_search text sent to the API
    fiscal_years: list[int] = field(default_factory=lambda: list(FISCAL_YEARS))

    def criteria(self, fiscal_years: list[int] | None = None) -> dict:
        return make_criteria(self.search_text, fiscal_years if fiscal_years is not None else self.fiscal_years)


def phrase(term: str) -> str:
    """Quote a term as an exact phrase for the RePORTER text search."""
    t = term.replace('"', " ").strip()
    if not t:
        raise ValueError("empty search term")
    return f'"{t}"'


def or_group(terms: list[str]) -> str:
    """("a" OR "b") for several terms; a bare quoted phrase for one."""
    if not terms:
        raise ValueError("no terms")
    quoted = [phrase(t) for t in terms]
    return quoted[0] if len(quoted) == 1 else "(" + " OR ".join(quoted) + ")"


def make_criteria(search_text: str, fiscal_years: list[int]) -> dict:
    return {
        "fiscal_years": list(fiscal_years),
        "advanced_text_search": {"operator": OPERATOR, "search_field": SEARCH_FIELD, "search_text": search_text},
    }


def _slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def build_condition_queries(conditions: list[dict]) -> list[QuerySpec]:
    out = []
    for c in conditions:
        for term in c.get("search_terms", []):
            out.append(QuerySpec(
                query_id=f"cond:{c['id']}:{_slug(term)}", query_type="condition_term",
                condition_id=c["id"], search_term=term, search_text=phrase(term)))
    return out


def build_augmented_queries(conditions: list[dict], members: list[str]) -> list[QuerySpec]:
    by_id = {c["id"]: c for c in conditions}
    out = []
    for cid in members:
        cond_text = or_group(by_id[cid]["search_terms"])
        for term, variants in AUGMENT_TERMS.items():
            out.append(QuerySpec(
                query_id=f"aug:{cid}:{_slug(term)}", query_type="augmented", condition_id=cid,
                search_term=term, search_text=f"{cond_text} AND {or_group(variants)}"))
    return out


def all_queries() -> list[QuerySpec]:
    cfg = load_config("conditions")
    return build_condition_queries(cfg["conditions"]) + build_augmented_queries(
        cfg["conditions"], cfg["demo_cluster"]["members"])


# --------------------------------------------------------------------------- fetching
def page_offsets(total: int, limit: int = PAGE_LIMIT, max_offset: int = MAX_OFFSET) -> list[int]:
    """Offsets needed to page through `total` records; raises if the API offset cap is exceeded."""
    if total > max_offset + 1:
        raise ValueError(f"{total} records exceed the API offset cap ({max_offset}); split the query")
    return list(range(0, total, limit)) or [0]


def needs_fiscal_year_split(total: int) -> bool:
    return total > MAX_RETRIEVABLE


def _post_page(criteria: dict, offset: int, refresh: bool, poster=None) -> tuple[dict, str, bool]:
    body = {"criteria": criteria, "offset": offset, "limit": PAGE_LIMIT, **SORT}
    if poster is not None:
        return poster(body), utc_now_iso(), False
    r = post(API_URL, json_body=body, expect_json=True, refresh=refresh)
    r.raise_for_status()
    return r.json(), r.fetched_at, r.from_cache


def fetch_criteria(criteria: dict, *, refresh: bool = False, poster=None) -> tuple[list[dict], list[dict]]:
    """All records for one criteria set. Splits by fiscal year if the total exceeds the offset cap.

    Returns (records, request_log). `poster` (body -> response JSON) is for tests.
    """
    first, fetched_at, cached = _post_page(criteria, 0, refresh, poster)
    total = int(first["meta"]["total"])
    fys = criteria.get("fiscal_years") or []
    if needs_fiscal_year_split(total):
        if len(fys) <= 1:
            raise RuntimeError(f"a single fiscal year returns {total} records (> {MAX_RETRIEVABLE}); narrow the query")
        recs, log = [], [{"criteria": criteria, "offset": 0, "total": total, "split_by_fiscal_year": True,
                          "fetched_at": fetched_at, "from_cache": cached}]
        for fy in fys:
            r, lg = fetch_criteria({**criteria, "fiscal_years": [fy]}, refresh=refresh, poster=poster)
            recs += r
            log += lg
        return recs, log
    recs = list(first.get("results") or [])
    log = [{"criteria": criteria, "offset": 0, "total": total, "n_results": len(recs),
            "search_id": first["meta"].get("search_id"),
            "reporter_url": (first["meta"].get("properties") or {}).get("URL"),
            "fetched_at": fetched_at, "from_cache": cached}]
    for off in page_offsets(total)[1:]:
        page, fetched_at, cached = _post_page(criteria, off, refresh, poster)
        res = list(page.get("results") or [])
        recs += res
        log.append({"criteria": criteria, "offset": off, "total": int(page["meta"]["total"]),
                    "n_results": len(res), "fetched_at": fetched_at, "from_cache": cached})
        if len(res) < PAGE_LIMIT:
            break
    return recs, log


def fetch_all(queries: list[QuerySpec], *, refresh: bool = False, poster=None, verbose: bool = True):
    """Run every query. Returns (records by appl_id, hits [(query_id, appl_id)], query log)."""
    records: dict[int, dict] = {}
    hits: list[tuple[str, int]] = []
    qlog = []
    for i, q in enumerate(queries):
        recs, log = fetch_criteria(q.criteria(), refresh=refresh, poster=poster)
        ids = []
        for r in recs:
            aid = int(r["appl_id"])
            records.setdefault(aid, r)
            ids.append(aid)
        uniq = sorted(set(ids))
        hits += [(q.query_id, a) for a in uniq]
        top = log[0]
        qlog.append({**asdict(q), "request_body_template": {"criteria": q.criteria(), "limit": PAGE_LIMIT, **SORT},
                     "api_total": top["total"], "n_records_retrieved": len(recs), "n_unique_appl_ids": len(uniq),
                     "split_by_fiscal_year": bool(top.get("split_by_fiscal_year")),
                     "search_id": top.get("search_id"), "reporter_url": top.get("reporter_url"),
                     "requests": log})
        if top["total"] != len(uniq) and not top.get("split_by_fiscal_year"):
            qlog[-1]["warning"] = f"api total {top['total']} != unique appl_ids retrieved {len(uniq)}"
        if verbose:
            print(f"[{i + 1}/{len(queries)}] {q.query_id}: total={top['total']} unique={len(uniq)}", flush=True)
    return records, hits, qlog


# --------------------------------------------------------------------------- raw persistence
def raw_paths() -> dict[str, Path]:
    d = raw_dir(SOURCE_ID)
    return {"records": d / "records.jsonl.gz", "hits": d / "query_hits.csv", "queries": d / "queries.json"}


def save_raw(records: dict[int, dict], hits: list[tuple[str, int]], qlog: list[dict]) -> str:
    p = raw_paths()
    # mtime=0 keeps the gzip header (and so the MANIFEST sha256) identical across re-runs of the same data
    with open(p["records"], "wb") as raw, gzip.GzipFile(filename="records.jsonl", mode="wb", fileobj=raw,
                                                         mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8") as fh:
        for aid in sorted(records):
            fh.write(json.dumps(records[aid], sort_keys=True) + "\n")
    pd.DataFrame(hits, columns=["query_id", "appl_id"]).sort_values(["query_id", "appl_id"]).to_csv(p["hits"], index=False)
    fetched = [r["fetched_at"] for q in qlog for r in q["requests"] if r.get("fetched_at")]
    retrieved_at = max(fetched) if fetched else utc_now_iso()
    p["queries"].write_text(json.dumps({
        "source_id": SOURCE_ID, "api_url": API_URL, "method": "POST", "fiscal_years": FISCAL_YEARS,
        "search_field": SEARCH_FIELD, "operator": OPERATOR, "page_limit": PAGE_LIMIT, "sort": SORT,
        "retrieved_at_first": min(fetched) if fetched else None, "retrieved_at_last": retrieved_at,
        "n_queries": len(qlog), "queries": qlog}, indent=1, default=str))
    for key, path in p.items():
        _manifest_record(SOURCE_ID, path.name, {
            "url": API_URL, "final_url": API_URL, "bytes": path.stat().st_size, "sha256": sha256_file(path),
            "retrieved_at": retrieved_at, "last_modified": None, "etag": None,
            "content_type": {"records": "application/gzip", "hits": "text/csv", "queries": "application/json"}[key],
            "note": f"assembled from POST {API_URL} responses (exact request bodies in queries.json)",
        })
    return retrieved_at


def load_raw() -> tuple[dict[int, dict], list[tuple[str, int]], dict]:
    p = raw_paths()
    records = {}
    with gzip.open(p["records"], "rt", encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            records[int(r["appl_id"])] = r
    h = pd.read_csv(p["hits"])
    hits = list(zip(h["query_id"], h["appl_id"].astype(int)))
    return records, hits, json.loads(p["queries"].read_text())


def fetch_docs() -> None:
    (raw_dir(SOURCE_ID) / "docs").mkdir(parents=True, exist_ok=True)
    for fname, url in DOC_URLS.items():
        try:
            download_file(url, SOURCE_ID, f"docs/{fname}", allow_html=False)
        except Exception as exc:  # documentation only; never blocks ingestion
            print(f"doc download failed {url}: {exc}")


# --------------------------------------------------------------------------- flattening
def parse_terms(terms: str | None) -> str | None:
    """'<A><B C>' -> 'A; B C'."""
    if not terms:
        return None
    parts = [t.strip() for t in re.findall(r"<([^<>]*)>", terms) if t.strip()]
    return "; ".join(parts) if parts else terms.strip()


def _date(s):
    return pd.to_datetime(s, errors="coerce").date() if s else None


def _num(x):
    return float(x) if x is not None else np.nan


def flatten_record(r: dict) -> dict:
    org = r.get("organization") or {}
    otype = r.get("organization_type") or {}
    pis = r.get("principal_investigators") or []
    contact = [p for p in pis if p.get("is_contact_pi")]
    ic_admin = r.get("agency_ic_admin") or {}
    fund = r.get("agency_ic_fundings") or []
    geo = r.get("geo_lat_lon") or {}
    ss = r.get("full_study_section") or {}
    sub = r.get("subproject_id")
    zipc = org.get("org_zipcode")
    return {
        "appl_id": int(r["appl_id"]),
        "subproject_id": str(sub) if sub not in (None, "") else None,
        "is_subproject": sub not in (None, ""),
        "project_num": r.get("project_num"),
        "core_project_num": r.get("core_project_num"),
        "fiscal_year": int(r["fiscal_year"]) if r.get("fiscal_year") is not None else None,
        "project_title": r.get("project_title"),
        "abstract_text": r.get("abstract_text"),
        "public_health_relevance": r.get("phr_text"),
        "terms": parse_terms(r.get("terms")),
        "pref_terms": r.get("pref_terms"),
        "pi_names": "; ".join(p.get("full_name") or "" for p in pis) or None,
        "pi_profile_ids": "; ".join(str(p.get("profile_id")) for p in pis if p.get("profile_id") is not None) or None,
        "n_pis": len(pis),
        "contact_pi_name": r.get("contact_pi_name") or (contact[0].get("full_name") if contact else None),
        "contact_pi_profile_id": str(contact[0]["profile_id"]) if contact and contact[0].get("profile_id") is not None else None,
        "principal_investigators": json.dumps(pis, sort_keys=True) if pis else None,
        "program_officer_names": "; ".join(p.get("full_name") or "" for p in (r.get("program_officers") or [])) or None,
        "org_name": org.get("org_name"),
        "org_ipf_code": str(org["org_ipf_code"]) if org.get("org_ipf_code") not in (None, "") else None,
        "org_uei": org.get("primary_uei"),
        "org_city": org.get("org_city") or org.get("city"),
        "org_state": org.get("org_state"),
        "org_zipcode": str(zipc) if zipc not in (None, "") else None,
        "org_country": org.get("org_country") or org.get("country"),
        "org_dept_type": org.get("dept_type"),
        "org_type_name": otype.get("name"),
        "org_type_code": otype.get("code"),
        "cong_dist": r.get("cong_dist"),
        "reporter_lat": _num(geo.get("lat")),
        "reporter_lon": _num(geo.get("lon")),
        "award_amount": _num(r.get("award_amount")),
        "direct_cost_amt": _num(r.get("direct_cost_amt")),
        "indirect_cost_amt": _num(r.get("indirect_cost_amt")),
        "agency_code": r.get("agency_code"),
        "admin_ic_code": ic_admin.get("code"),
        "admin_ic_abbreviation": ic_admin.get("abbreviation"),
        "admin_ic_name": ic_admin.get("name"),
        "funding_ics": "; ".join(sorted({f.get("abbreviation") or f.get("code") or "" for f in fund})) or None,
        "ic_fundings": json.dumps(fund, sort_keys=True) if fund else None,
        "activity_code": r.get("activity_code"),
        "award_type": r.get("award_type"),
        "funding_mechanism": r.get("funding_mechanism"),
        "mechanism_code_dc": r.get("mechanism_code_dc"),
        "cfda_code": r.get("cfda_code"),
        "opportunity_number": r.get("opportunity_number"),
        "study_section": ss.get("name"),
        "project_start_date": _date(r.get("project_start_date")),
        "project_end_date": _date(r.get("project_end_date")),
        "budget_start": _date(r.get("budget_start")),
        "budget_end": _date(r.get("budget_end")),
        "award_notice_date": _date(r.get("award_notice_date")),
        "is_active": bool(r["is_active"]) if r.get("is_active") is not None else None,
        "is_new": bool(r["is_new"]) if r.get("is_new") is not None else None,
        "covid_response": "; ".join(r["covid_response"]) if r.get("covid_response") else None,
        "arra_funded": r.get("arra_funded"),
        "spending_categories": "; ".join(str(x) for x in (r.get("spending_categories") or [])) or None,
        "spending_categories_desc": r.get("spending_categories_desc"),
        "project_detail_url": r.get("project_detail_url"),
        "date_added": r.get("date_added"),
    }


# --------------------------------------------------------------------------- local phrase matching
_SEP = r"[\s\-/‐-―−]+"   # space, hyphen, slash, unicode dashes (e.g. "Ehlers–Danlos")

# Ambiguous acronyms among the configured search terms, with the context that
# must co-occur for the hit to count as the condition. Measured 2026-09-23:
# 20 long_covid appl_ids matched only "PASC": 2 mention COVID, 2 are the
# RECOVER-funded iDRAW project (opportunity OTA-21-015A, text never says
# COVID), and 16 are about pancreatic stellate cells (PaSC) or a phosphorylated
# inflammasome adaptor (pASC). The RECOVER opportunity number is therefore
# accepted as context (opportunity_number is searched too). "ME/CFS" needs no
# context (its 4 acronym-only hits were ME/CFS projects on inspection).
ACRONYM_CONTEXT: dict[str, str] = {
    "PASC": r"covid|sars[\s\-]?cov[\s\-]?2|coronavirus|post[\s\-]?acute sequelae|OTA-21-015",
}


def phrase_regex(term: str) -> re.Pattern:
    """Case-insensitive whole-phrase regex; space, hyphen, slash and dashes are interchangeable separators."""
    words = [w for w in re.split(_SEP, term.strip()) if w]
    return re.compile(r"(?<![A-Za-z0-9])" + _SEP.join(map(re.escape, words)) + r"(?![A-Za-z0-9])", re.I)


def acronym_context(term: str, *texts) -> bool | None:
    """True/False if `term` is an ambiguous acronym and its context is/is not present; None otherwise."""
    pat = ACRONYM_CONTEXT.get(term)
    if pat is None:
        return None
    rx = re.compile(pat, re.I)
    return any(isinstance(t, str) and rx.search(t) for t in texts)


def match_fields(terms: list[str], title, abstract, phr, terms_field) -> dict:
    pats = [phrase_regex(t) for t in terms]

    def hit(text):
        return isinstance(text, str) and bool(text) and any(p.search(text) for p in pats)
    in_title = hit(title)
    in_abs = hit(abstract) or hit(phr)
    in_terms = hit(terms_field)
    tier = "title_abstract" if (in_title or in_abs) else ("terms_only" if in_terms else "not_found_locally")
    return {"matched_in_title": in_title, "matched_in_abstract": in_abs, "matched_in_terms": in_terms, "match_tier": tier}


# --------------------------------------------------------------------------- funding totals
def obligation_total(df: pd.DataFrame, amount_col: str = "award_amount") -> float:
    """Sum award obligations once per appl_id without parent/subproject double counting.

    Rows are de-duplicated on appl_id. Within one (project_num, fiscal_year),
    when the parent record is present with an amount, only the parent is
    counted: RePORTER defines it as the total for all subprojects, and it equals
    the sum of the award's IC fundings. This holds even when the subproject
    amounts sum to more than the parent, a RePORTER inconsistency measured in
    DATA_AUDIT.md. Otherwise the matched subprojects are summed.
    """
    if df.empty:
        return 0.0
    d = df.drop_duplicates("appl_id")
    total = 0.0
    for _, g in d.groupby([d["project_num"].fillna("?"), d["fiscal_year"].fillna(-1)], sort=False):
        parent = g[~g["is_subproject"].astype(bool) & g[amount_col].notna()]
        total += float(parent[amount_col].sum()) if len(parent) else float(g[amount_col].fillna(0).sum())
    return total


def parent_subproject_consistency(df: pd.DataFrame, amount_col: str = "award_amount", tol: float = 1.0) -> dict:
    """Compare parent amounts with the sum of their subprojects present in df.

    Over (project_num, fiscal_year) groups holding a parent with an amount and
    at least one subproject, counts groups where the parent equals the subproject
    sum, exceeds it (other subprojects were not retrieved), or is smaller (the
    subprojects sum to more than the award, a RePORTER inconsistency), and the
    total excess in the smaller-parent groups.
    """
    d = df.drop_duplicates("appl_id")
    sub = d["is_subproject"].astype(bool)
    par = d[~sub & d[amount_col].notna()].groupby(["project_num", "fiscal_year"])[amount_col].sum()
    subs = d[sub].groupby(["project_num", "fiscal_year"])[amount_col].sum()
    both = par.to_frame("parent").join(subs.to_frame("subs"), how="inner")
    diff = both["parent"] - both["subs"]
    smaller = diff < -tol
    return {"groups": int(len(both)), "equal": int((diff.abs() <= tol).sum()), "parent_larger": int((diff > tol).sum()),
            "parent_smaller": int(smaller.sum()), "parent_smaller_excess": float(-diff[smaller].sum()),
            "parent_smaller_groups": [f"{pn} FY{fy}" for pn, fy in both.index[smaller]]}


def subproject_overlap_rows(df: pd.DataFrame) -> int:
    """Number of subproject rows whose parent (same project_num + FY) is also in df."""
    d = df.drop_duplicates("appl_id")
    parents = set(zip(d.loc[~d["is_subproject"].astype(bool), "project_num"],
                      d.loc[~d["is_subproject"].astype(bool), "fiscal_year"]))
    subs = d[d["is_subproject"].astype(bool)]
    return int(sum((pn, fy) in parents for pn, fy in zip(subs["project_num"], subs["fiscal_year"])))


# --------------------------------------------------------------------------- geocoding
GEO_ZCTA = "zip5_as_zcta_internal_point"
GEO_REPORTER = "reporter_org_point_state_verified"


def geocode_orgs(proj: pd.DataFrame) -> pd.DataFrame:
    """Organization location for each row.

    1. U.S. organizations only (a foreign postal code such as Nigeria's 900246
       must not be read as a U.S. ZIP): ZIP5 -> ZCTA internal point via
       crosswalk.zip_to_geo (resolution zcta_centroid).
    2. When the ZIP has no ZCTA (unique/PO-box ZIPs such as UCSF 94143 or Yale
       06520), fall back to RePORTER's own organization coordinate, accepted
       only if it falls inside a 2024 county of the organization's own state
       (resolution point, CONVENTIONS section 4).
    Everything else stays ungeocoded with a reason.
    """
    import geopandas as gpd

    from ..config import PROCESSED
    from ..geography.crosswalk import STATE_FIPS, zip_to_geo

    us = proj["org_country"].fillna("").str.upper().eq("UNITED STATES")
    geo = zip_to_geo(proj["org_zipcode"].where(us, "").fillna(""))
    out = pd.DataFrame(index=proj.index)
    out["zip5"] = geo["zip5"].where(us)
    for c in ("lat", "lon", "county_fips", "state_fips"):
        out[c] = geo[c].where(us)
    zok = us & (geo["geocode_method"] == GEO_ZCTA)
    out["geocode_method"] = np.where(zok, GEO_ZCTA, "not_geocoded")
    out["geocode_failure_reason"] = np.select(
        [zok, ~us, us & out["zip5"].isna()], [None, "non_us_organization", "no_zip_in_record"],
        default="zip_has_no_zcta")
    out["zip_geocode_failed"] = us & ~zok & out["zip5"].notna()
    out["source_geographic_resolution"] = np.where(zok, "zcta_centroid", "none")
    out["county_assignment"] = np.where(zok & out["county_fips"].notna(), "zcta_point_in_2024_county", None)
    # A few lakefront ZCTA internal points (e.g. 60611, Chicago) fall outside the
    # clipped 5m county polygons, so zcta_centroids has no county for them. Use
    # the Census 2020 ZCTA-county relationship (largest land share) instead;
    # skip Connecticut, whose 2020 counties are not 2024 counties.
    miss = zok & out["county_fips"].isna()
    if miss.any() and table_exists("zcta_county_crosswalk"):
        xw = read_table("zcta_county_crosswalk", columns=["zcta", "county_fips_2020", "land_share_of_zcta"])
        best = (xw.sort_values("land_share_of_zcta", ascending=False).drop_duplicates("zcta")
                .set_index("zcta")["county_fips_2020"])
        fill = out.loc[miss, "zip5"].map(best)
        fill = fill[fill.notna() & ~fill.str.startswith("09")]
        out.loc[fill.index, "county_fips"] = fill
        out.loc[fill.index, "state_fips"] = fill.str[:2]
        out.loc[fill.index, "county_assignment"] = "zcta_county_rel2020_max_land_share"

    need = us & ~zok & proj["reporter_lat"].notna() & proj["reporter_lon"].notna()
    if need.any():
        counties = gpd.read_parquet(PROCESSED / "county_boundaries_2024.geoparquet")[["geo_id", "state_abbr", "geometry"]]
        pts = gpd.GeoDataFrame(index=proj.index[need], crs=4326, geometry=gpd.points_from_xy(
            proj.loc[need, "reporter_lon"], proj.loc[need, "reporter_lat"])).to_crs(counties.crs)
        j = gpd.sjoin(pts, counties, how="left", predicate="within")
        j = j[~j.index.duplicated(keep="first")]
        ok = j["geo_id"].notna() & (j["state_abbr"] == proj.loc[j.index, "org_state"])
        idx = j.index[ok.values]
        out.loc[idx, "lat"] = proj.loc[idx, "reporter_lat"]
        out.loc[idx, "lon"] = proj.loc[idx, "reporter_lon"]
        out.loc[idx, "county_fips"] = j.loc[idx, "geo_id"]
        out.loc[idx, "state_fips"] = j.loc[idx, "geo_id"].str[:2]
        out.loc[idx, "geocode_method"] = GEO_REPORTER
        out.loc[idx, "county_assignment"] = "reporter_point_in_2024_county"
        out.loc[idx, "geocode_failure_reason"] = None
        out.loc[idx, "source_geographic_resolution"] = "point"
        bad = j.index[~ok.values]
        out.loc[bad, "geocode_failure_reason"] = out.loc[bad, "geocode_failure_reason"] + "; reporter_point_not_in_org_state"
    fips_abbr = out["state_fips"].map(STATE_FIPS)
    out["geo_state_matches_org_state"] = np.where(out["state_fips"].notna(), fips_abbr == proj["org_state"], None)
    return out


# --------------------------------------------------------------------------- table building
def build_tables(records: dict[int, dict], hits: list[tuple[str, int]], qmeta: dict, retrieved_at: str):
    conds = load_config("conditions")["conditions"]
    terms_by_cond = {c["id"]: c["search_terms"] for c in conds}
    qinfo = {q["query_id"]: q for q in qmeta["queries"]}
    version = f"RePORTER Project API v2.0; FY{FISCAL_YEARS[0]}-FY{FISCAL_YEARS[-1]}; retrieved {retrieved_at[:10]}"

    proj = pd.DataFrame([flatten_record(r) for r in records.values()]).sort_values("appl_id").reset_index(drop=True)
    geo = geocode_orgs(proj)
    for c in geo.columns:
        if c != "source_geographic_resolution":
            proj[c] = geo[c].values
    hit_df = pd.DataFrame(hits, columns=["query_id", "appl_id"])
    hit_df["query_type"] = hit_df["query_id"].map(lambda q: qinfo[q]["query_type"])
    hit_df["condition_id"] = hit_df["query_id"].map(lambda q: qinfo[q]["condition_id"])
    cond_hits = hit_df[hit_df["query_type"] == "condition_term"]
    matched = cond_hits.groupby("appl_id")["condition_id"].agg(lambda s: "; ".join(sorted(set(s))))
    proj["matched_condition_ids"] = proj["appl_id"].map(matched)
    aug_hits = hit_df[hit_df["query_type"] == "augmented"]
    aug_by = aug_hits.assign(t=aug_hits["query_id"].map(lambda q: qinfo[q]["search_term"])).groupby("appl_id")["t"]
    proj["augmented_terms_hit"] = proj["appl_id"].map(aug_by.agg(lambda s: "; ".join(sorted(set(s)))))

    proj = add_provenance(
        proj, data_layer="facility", source_name=SOURCE_NAME, source_version=version, retrieved_at=retrieved_at,
        evidence_type="grant_record", source_record_id="appl_id",
        source_geographic_resolution="zcta_centroid", evidence_level="award_record",
        provenance_notes=PROVENANCE_NOTE_PROJECTS)
    proj["source_geographic_resolution"] = geo["source_geographic_resolution"].values

    text_cols = proj.set_index("appl_id")[["project_title", "abstract_text", "public_health_relevance", "terms",
                                           "opportunity_number"]]
    rcdc = proj.set_index("appl_id")["spending_categories_desc"].fillna("")

    def _link_rows(df, term_lookup):
        rows = []
        for q, aid, cid in zip(df["query_id"], df["appl_id"], df["condition_id"]):
            t = text_cols.loc[aid]
            rows.append(match_fields(term_lookup(q), t["project_title"], t["abstract_text"],
                                     t["public_health_relevance"], t["terms"]))
        return pd.DataFrame(rows, index=df.index)

    pc = cond_hits.copy()
    pc["search_term"] = pc["query_id"].map(lambda q: qinfo[q]["search_term"])
    pc["search_text"] = pc["query_id"].map(lambda q: qinfo[q]["search_text"])
    pc = pc.join(_link_rows(pc, lambda q: [qinfo[q]["search_term"]]))
    pc["rcdc_category"] = pc["condition_id"].map(RCDC_CATEGORY_FOR_CONDITION)
    pc["rcdc_category_assigned"] = [
        (cat in [s.strip() for s in rcdc.loc[a].split(";")]) if isinstance(cat, str) else None
        for a, cat in zip(pc["appl_id"], pc["rcdc_category"])]
    pc["acronym_context_found"] = [
        acronym_context(term, *text_cols.loc[a].tolist()) for term, a in zip(pc["search_term"], pc["appl_id"])]
    pc["likely_false_positive"] = pc["acronym_context_found"].eq(False)
    pc["fiscal_year"] = pc["appl_id"].map(proj.set_index("appl_id")["fiscal_year"])
    pc = pc.drop(columns=["query_type"]).sort_values(["condition_id", "query_id", "appl_id"]).reset_index(drop=True)
    pc = add_provenance(
        pc, data_layer="facility", source_name=SOURCE_NAME, source_version=version, retrieved_at=retrieved_at,
        evidence_type="text_mined", source_record_id=lambda d: d["appl_id"].astype(str) + ":" + d["query_id"],
        source_geographic_resolution="none", evidence_level=pc["match_tier"], provenance_notes=PROVENANCE_NOTE_LINKS)

    qh = aug_hits.copy()
    qh["augmented_term"] = qh["query_id"].map(lambda q: qinfo[q]["search_term"])
    qh["search_text"] = qh["query_id"].map(lambda q: qinfo[q]["search_text"])
    m_aug = _link_rows(qh, lambda q: AUGMENT_TERMS[qinfo[q]["search_term"]])
    m_cond = _link_rows(qh, lambda q: terms_by_cond[qinfo[q]["condition_id"]])
    qh["augmented_term_match_tier"] = m_aug["match_tier"]
    qh["condition_term_match_tier"] = m_cond["match_tier"]
    qh["both_in_title_abstract"] = (m_aug["match_tier"] == "title_abstract") & (m_cond["match_tier"] == "title_abstract")
    fp = pc.groupby(["condition_id", "appl_id"])["likely_false_positive"].all()
    qh["condition_likely_false_positive"] = [bool(fp.get((c, a), False)) for c, a in zip(qh["condition_id"], qh["appl_id"])]
    qh["fiscal_year"] = qh["appl_id"].map(proj.set_index("appl_id")["fiscal_year"])
    qh = qh.drop(columns=["query_type"]).sort_values(["condition_id", "augmented_term", "appl_id"]).reset_index(drop=True)
    qh = add_provenance(
        qh, data_layer="facility", source_name=SOURCE_NAME, source_version=version, retrieved_at=retrieved_at,
        evidence_type="text_mined", source_record_id=lambda d: d["appl_id"].astype(str) + ":" + d["query_id"],
        source_geographic_resolution="none",
        evidence_level=np.where(qh["both_in_title_abstract"], "both_in_title_abstract", "includes_terms_field_match"),
        provenance_notes=PROVENANCE_NOTE_LINKS + " A hit is a grant mentioning the measurement term, not evidence "
                         "that the measurement works for the condition.")
    return proj, pc, qh


# --------------------------------------------------------------------------- helpers for MCP / scoring
def _conditions() -> list[str]:
    return [c["id"] for c in load_config("conditions")["conditions"]]


def projects_for_condition(condition_id: str, *, fiscal_years: list[int] | None = None,
                           active_only: bool = False, include_terms_only: bool = False,
                           include_likely_false_positives: bool = False) -> pd.DataFrame:
    """RePORTER projects (one row per appl_id) whose text search matched `condition_id`.

    Default = precision view: the condition phrase occurs in the title, abstract
    or narrative, and ambiguous-acronym hits without context are dropped.
    `include_terms_only=True` adds projects matched only through RePORTER's
    thesaurus-expanded terms field. That tier agrees with NIH's own RCDC
    category far less often (measured per condition in DATA_AUDIT.md), so it is
    off by default. Adds matched_search_terms and best_match_tier.
    Empty frame when nothing matched; ValueError for an unknown id.
    """
    if condition_id not in _conditions():
        raise ValueError(f"unknown condition_id {condition_id!r}; known: {_conditions()}")
    links = read_table("nih_project_conditions")
    links = links[links["condition_id"] == condition_id]
    if not include_likely_false_positives:
        links = links[~links["likely_false_positive"].fillna(False).astype(bool)]
    rank = {"title_abstract": 0, "terms_only": 1, "not_found_locally": 2}
    agg = links.assign(r=links["match_tier"].map(rank)).groupby("appl_id").agg(
        matched_search_terms=("search_term", lambda s: "; ".join(sorted(set(s)))),
        best_rank=("r", "min"))
    agg["best_match_tier"] = agg["best_rank"].map({v: k for k, v in rank.items()})
    proj = read_table("nih_projects")
    out = proj.merge(agg.drop(columns="best_rank"), left_on="appl_id", right_index=True, how="inner")
    if fiscal_years:
        out = out[out["fiscal_year"].isin(fiscal_years)]
    if active_only:
        out = out[out["is_active"].fillna(False).astype(bool)]
    if not include_terms_only:
        out = out[out["best_match_tier"] == "title_abstract"]
    out.insert(0, "condition_id", condition_id)
    return out.sort_values(["fiscal_year", "appl_id"], ascending=[False, True]).reset_index(drop=True)


def research_centers_for_condition(condition_id: str, *, fiscal_years: list[int] | None = None,
                                   include_terms_only: bool = False, include_likely_false_positives: bool = False,
                                   top_n: int | None = None) -> pd.DataFrame:
    """Organizations holding RePORTER projects matched to `condition_id`, aggregated per organization.

    A research-capability signal (grants, active projects, investigators,
    award obligations), not patient burden or clinical capacity. Award totals
    count each appl_id once and do not double count subprojects. Filters as in
    projects_for_condition (precision view by default).
    """
    p = projects_for_condition(condition_id, fiscal_years=fiscal_years, include_terms_only=include_terms_only,
                               include_likely_false_positives=include_likely_false_positives)
    cols = ["condition_id", "org_key", "org_name", "org_ipf_code", "org_city", "org_state", "org_zip5", "org_country",
            "lat", "lon", "county_fips", "state_fips", "geocode_method", "n_projects", "n_core_projects",
            "n_active_projects", "first_fiscal_year", "last_fiscal_year", "award_obligations_total_usd",
            "n_award_amount_missing", "n_investigators", "contact_pis", "activity_codes", "funding_ics",
            "n_title_abstract_matches", "recent_titles", "data_layer", "source_name", "source_record_id",
            "source_version", "retrieved_at", "source_geographic_resolution", "evidence_type", "evidence_level",
            "provenance_notes"]
    if p.empty:
        return pd.DataFrame(columns=cols)
    p = p.copy()
    p["org_key"] = p["org_ipf_code"].fillna("name:" + p["org_name"].fillna(UNKNOWN))
    rows = []
    for key, g in p.groupby("org_key", sort=False):
        latest = g.sort_values("fiscal_year").iloc[-1]
        pis = Counter(x for x in g["contact_pi_name"].dropna())
        inv = {i.strip() for s in g["pi_profile_ids"].dropna() for i in s.split(";") if i.strip()}
        titles = list(dict.fromkeys(g.sort_values("fiscal_year", ascending=False)["project_title"].dropna()))[:3]
        rows.append({
            "condition_id": condition_id, "org_key": key, "org_name": latest["org_name"],
            "org_ipf_code": latest["org_ipf_code"], "org_city": latest["org_city"], "org_state": latest["org_state"],
            "org_zip5": latest["zip5"], "org_country": latest["org_country"], "lat": latest["lat"],
            "lon": latest["lon"], "county_fips": latest["county_fips"], "state_fips": latest["state_fips"],
            "geocode_method": latest["geocode_method"], "n_projects": int(g["appl_id"].nunique()),
            "n_core_projects": int(g["core_project_num"].nunique()),
            "n_active_projects": int(g.loc[g["is_active"].fillna(False).astype(bool), "core_project_num"].nunique()),
            "first_fiscal_year": int(g["fiscal_year"].min()), "last_fiscal_year": int(g["fiscal_year"].max()),
            # NaN (not 0) when no record reports an amount, e.g. VA awards
            "award_obligations_total_usd": obligation_total(g) if g["award_amount"].notna().any() else np.nan,
            "n_award_amount_missing": int(g["award_amount"].isna().sum()),
            "n_investigators": len(inv),
            "contact_pis": "; ".join(n for n, _ in pis.most_common(5)) or UNKNOWN,
            "activity_codes": "; ".join(sorted(g["activity_code"].dropna().unique())),
            "funding_ics": "; ".join(sorted({i.strip() for s in g["funding_ics"].dropna() for i in s.split(";")})),
            "n_title_abstract_matches": int((g["best_match_tier"] == "title_abstract").sum()),
            "recent_titles": " | ".join(titles),
        })
    out = pd.DataFrame(rows)
    out["data_layer"] = "facility"
    out["source_name"] = SOURCE_NAME
    out["source_record_id"] = condition_id + ":" + out["org_key"].astype(str)
    out["source_version"] = p["source_version"].iloc[0]
    out["retrieved_at"] = p["retrieved_at"].iloc[0]
    out["source_geographic_resolution"] = out["geocode_method"].map(
        {GEO_ZCTA: "zcta_centroid", GEO_REPORTER: "point"}).fillna("none")
    out["evidence_type"] = "grant_record"
    out["evidence_level"] = ("aggregated_award_records:"
                             + ("all_match_tiers" if include_terms_only else "title_abstract_matches"))
    out["provenance_notes"] = ("Aggregated RePORTER award records matched by text search; research-capability "
                               "signal only. n_active_projects counts distinct core projects flagged is_active. "
                               "Award obligations are total-cost obligations, FY-summed, not expenditure.")
    out = out.sort_values(["n_core_projects", "award_obligations_total_usd"], ascending=False).reset_index(drop=True)
    return (out.head(top_n) if top_n else out)[cols]


# --------------------------------------------------------------------------- audit
def _pct(n, d):
    return f"{n:,} ({100 * n / d:.1f}%)" if d else f"{n:,}"


def _md_table(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r.values) + " |")
    return "\n".join(lines)


def compute_stats(proj: pd.DataFrame, pc: pd.DataFrame, qh: pd.DataFrame, qmeta: dict) -> dict:
    conds = load_config("conditions")["conditions"]
    s: dict = {}
    s["n_projects"] = int(len(proj))
    s["n_core_projects"] = int(proj["core_project_num"].nunique())
    s["n_orgs"] = int(proj["org_ipf_code"].fillna(proj["org_name"]).nunique())
    s["n_queries"] = int(qmeta["n_queries"])
    s["n_condition_queries"] = sum(q["query_type"] == "condition_term" for q in qmeta["queries"])
    s["n_augmented_queries"] = sum(q["query_type"] == "augmented" for q in qmeta["queries"])
    s["n_requests"] = sum(len(q["requests"]) for q in qmeta["queries"])
    s["any_fy_split"] = any(q["split_by_fiscal_year"] for q in qmeta["queries"])
    s["zero_hit_queries"] = [q["query_id"] for q in qmeta["queries"] if q["api_total"] == 0]
    s["query_warnings"] = [(q["query_id"], q["warning"]) for q in qmeta["queries"] if q.get("warning")]
    s["n_subprojects"] = int(proj["is_subproject"].sum())
    s["n_active"] = int(proj["is_active"].fillna(False).astype(bool).sum())
    s["union_obligations"] = obligation_total(proj)
    s["union_naive_sum"] = float(proj["award_amount"].sum())
    s["union_sub_overlap"] = subproject_overlap_rows(proj)
    s["parent_sub"] = parent_subproject_consistency(proj)
    # award_amount vs the per-IC total_cost list, and vs direct + indirect
    ic_sum = proj["ic_fundings"].map(
        lambda j: sum(float(f.get("total_cost") or 0) for f in json.loads(j)) if isinstance(j, str) else np.nan)
    m = proj["award_amount"].notna() & ic_sum.notna()
    s["award_eq_ic_sum"] = (int(((proj["award_amount"] - ic_sum).abs() <= 1)[m].sum()), int(m.sum()))
    m = proj[["award_amount", "direct_cost_amt", "indirect_cost_amt"]].notna().all(axis=1)
    di = (proj["award_amount"] - proj["direct_cost_amt"] - proj["indirect_cost_amt"]).abs() > 1
    s["award_ne_direct_indirect"] = (int(di[m].sum()), int(m.sum()))
    s["max_date_added"] = str(proj["date_added"].dropna().max())
    s["agency"] = proj["agency_code"].fillna("(missing)").value_counts().to_dict()
    s["country"] = proj["org_country"].fillna("(missing)").value_counts().to_dict()

    text_idx = proj.set_index("appl_id")
    rank = {"title_abstract": 0, "terms_only": 1, "not_found_locally": 2}
    inv_rank = {v: k for k, v in rank.items()}
    rows = []
    for c in conds:
        cid = c["id"]
        links = pc[pc["condition_id"] == cid]
        ids = links["appl_id"].unique()
        sub = proj[proj["appl_id"].isin(ids)]
        tier = links.assign(r=links["match_tier"].map(rank)).groupby("appl_id")["r"].min().map(inv_rank)
        keep = links[~links["likely_false_positive"].astype(bool)]
        ktier = keep.assign(r=keep["match_tier"].map(rank)).groupby("appl_id")["r"].min().map(inv_rank)
        prec_ids = ktier[ktier == "title_abstract"].index
        prec = proj[proj["appl_id"].isin(prec_ids)]
        n_fp = int(len(set(ids) - set(keep["appl_id"])))
        rc = links.drop_duplicates("appl_id")["rcdc_category_assigned"]
        rows.append({
            "condition_id": cid,
            "search_terms": len(c["search_terms"]),
            "zero_hit_terms": int(sum(q["api_total"] == 0 for q in qmeta["queries"]
                                      if q["query_type"] == "condition_term" and q["condition_id"] == cid)),
            "projects (appl_id)": int(len(ids)),
            "core projects": int(sub["core_project_num"].nunique()),
            "active (appl_id)": int(sub["is_active"].fillna(False).astype(bool).sum()),
            "active core projects": int(sub.loc[sub["is_active"].fillna(False).astype(bool), "core_project_num"].nunique()),
            "title/abstract match": int((tier == "title_abstract").sum()),
            "terms-field only": int((tier == "terms_only").sum()),
            "not found locally": int((tier == "not_found_locally").sum()),
            "likely false positive": n_fp,
            "RCDC category": RCDC_CATEGORY_FOR_CONDITION.get(cid, "none"),
            "RCDC assigned": (int(rc.fillna(False).astype(bool).sum()) if cid in RCDC_CATEGORY_FOR_CONDITION else "n/a"),
            "orgs": int(sub["org_ipf_code"].fillna(sub["org_name"]).nunique()),
            "award obligations USD (dedup)": round(obligation_total(sub)),
            "naive sum USD": round(float(sub["award_amount"].sum())),
            "award missing": int(sub["award_amount"].isna().sum()),
            "precision view: projects": int(len(prec)),
            "precision view: core projects": int(prec["core_project_num"].nunique()),
            "precision view: active (appl_id)": int(prec["is_active"].fillna(False).astype(bool).sum()),
            "precision view: orgs": int(prec["org_ipf_code"].fillna(prec["org_name"]).nunique()),
            "precision view: award obligations USD (dedup)": round(obligation_total(prec)),
        })
    s["per_condition"] = pd.DataFrame(rows)
    main_cols = [c for c in s["per_condition"].columns if not c.startswith("precision view")]
    s["per_condition_main"] = s["per_condition"][main_cols]
    s["per_condition_precision"] = s["per_condition"][["condition_id"] + [c for c in s["per_condition"].columns
                                                                          if c.startswith("precision view")]]

    # RCDC agreement by match tier (NIH-administered records, conditions with a category)
    x = pc.assign(r=pc["match_tier"].map(rank)).groupby(["condition_id", "appl_id"]).agg(
        r=("r", "min"), rc=("rcdc_category_assigned", "first")).reset_index()
    x = x[x["rc"].notna()].merge(proj[["appl_id", "agency_code", "fiscal_year"]], on="appl_id")
    x = x[x["agency_code"] == "NIH"]
    # only fiscal years in which NIH assigned that category to at least one
    # retrieved record (RCDC is assigned after a fiscal year closes: no FY2026
    # record has any category; the long COVID category starts in FY2022)
    cat_years = {}
    for cid, cat in RCDC_CATEGORY_FOR_CONDITION.items():
        has = proj["spending_categories_desc"].fillna("").map(lambda d: cat in [t.strip() for t in d.split(";")])
        cat_years[cid] = set(proj.loc[has, "fiscal_year"])
    x = x[[fy in cat_years.get(c, set()) for c, fy in zip(x["condition_id"], x["fiscal_year"])]]
    s["rcdc_years"] = {c: f"FY{min(v)}-FY{max(v)}" for c, v in cat_years.items() if v}
    nih = proj["agency_code"] == "NIH"
    s["rcdc_fy2026"] = (int((nih & (proj["fiscal_year"] == 2026) & proj["spending_categories_desc"].notna()).sum()),
                        int((nih & (proj["fiscal_year"] == 2026)).sum()))
    s["rcdc_non_nih"] = (int((~nih & proj["spending_categories_desc"].notna()).sum()), int((~nih).sum()))
    x["tier"] = x["r"].map(inv_rank)
    rt = x.groupby(["condition_id", "tier"])["rc"].agg(n="size", rcdc_assigned="sum").reset_index()
    rt["rcdc_assigned"] = rt["rcdc_assigned"].astype(int)
    rt["pct_rcdc"] = (100 * rt["rcdc_assigned"] / rt["n"]).round(1)
    s["rcdc_by_tier"] = rt
    wide = rt.pivot(index="condition_id", columns="tier", values="pct_rcdc")
    if {"terms_only", "title_abstract"} <= set(wide.columns):
        cmp = wide.dropna(subset=["terms_only", "title_abstract"])
        s["rcdc_lower"] = sorted(cmp.index[cmp["terms_only"] < cmp["title_abstract"]])
        s["rcdc_not_lower"] = sorted(cmp.index[cmp["terms_only"] >= cmp["title_abstract"]])
    else:
        s["rcdc_lower"], s["rcdc_not_lower"] = [], []

    fy = pc.drop_duplicates(["condition_id", "appl_id"]).pivot_table(
        index="condition_id", columns="fiscal_year", values="appl_id", aggfunc="count", fill_value=0)
    s["per_condition_fy"] = fy.reset_index()
    s["overall_fy"] = proj["fiscal_year"].value_counts().sort_index().to_dict()

    per_term = [{"query_id": q["query_id"], "search_text": q["search_text"], "api_total": q["api_total"],
                 "retrieved unique": q["n_unique_appl_ids"]}
                for q in qmeta["queries"] if q["query_type"] == "condition_term"]
    s["per_term"] = pd.DataFrame(per_term)

    # augmented queries + local re-check of the API's boolean logic: among the
    # condition's hits, which records contain a measurement variant anywhere
    # in title/abstract/narrative/terms?
    aug_rows = []
    for q in (q for q in qmeta["queries"] if q["query_type"] == "augmented"):
        cid, term = q["condition_id"], q["search_term"]
        api = set(qh.loc[qh["query_id"] == q["query_id"], "appl_id"])
        cset = pc.loc[pc["condition_id"] == cid, "appl_id"].unique()
        pats = [phrase_regex(v) for v in AUGMENT_TERMS[term]]
        local = {a for a in cset if any(isinstance(t, str) and p.search(t) for p in pats for t in (
            text_idx.at[a, "project_title"], text_idx.at[a, "abstract_text"],
            text_idx.at[a, "public_health_relevance"], text_idx.at[a, "terms"]))}
        sub = qh[qh["query_id"] == q["query_id"]]
        aug_rows.append({"condition_id": cid, "augmented_term": term, "api_hits": len(api),
                         "both_in_title_abstract": int(sub["both_in_title_abstract"].sum()),
                         "condition_likely_false_positive": int(sub["condition_likely_false_positive"].sum()),
                         "local_regex_hits": len(local), "api_and_local": len(api & local),
                         "api_only": len(api - local), "local_only": len(local - api)})
    s["augmented"] = pd.DataFrame(aug_rows)
    cond_ids = {(cid, a) for cid, a in zip(pc["condition_id"], pc["appl_id"])}
    s["aug_not_in_condition_set"] = int(sum((cid, a) not in cond_ids for cid, a in zip(qh["condition_id"], qh["appl_id"])))

    # multi-condition overlap and double-count check
    per_app = pc.drop_duplicates(["appl_id", "condition_id"]).groupby("appl_id")["condition_id"].nunique()
    s["multi_condition_appl_ids"] = int((per_app > 1).sum())
    s["sum_of_condition_obligations"] = float(s["per_condition"]["award obligations USD (dedup)"].sum())

    # ambiguous-acronym check (ACRONYM_CONTEXT)
    s["acronym_checks"] = []
    for term in ACRONYM_CONTEXT:
        lk = pc[pc["search_term"] == term]
        for cid, g in lk.groupby("condition_id"):
            other = set(pc.loc[(pc["condition_id"] == cid) & (pc["search_term"] != term), "appl_id"])
            only = set(g["appl_id"]) - other
            fp = set(g.loc[g["likely_false_positive"].astype(bool), "appl_id"])
            s["acronym_checks"].append({"term": term, "condition_id": cid, "hits": int(g["appl_id"].nunique()),
                                        "acronym_only_hits": len(only), "no_context": len(fp),
                                        "no_context_titles": sorted({text_idx.at[a, "project_title"][:60] for a in fp})[:6]})

    # geocoding
    us = proj["org_country"].fillna("").str.upper().eq("UNITED STATES")
    s["n_us"] = int(us.sum())
    s["n_nonus"] = int((~us).sum())
    s["n_zcta"] = int((proj["geocode_method"] == GEO_ZCTA).sum())
    s["n_reporter_point"] = int((proj["geocode_method"] == GEO_REPORTER).sum())
    s["n_geocoded"] = s["n_zcta"] + s["n_reporter_point"]
    s["n_us_geocoded"] = int((us & proj["geocode_method"].isin([GEO_ZCTA, GEO_REPORTER])).sum())
    s["fail_reasons"] = proj.loc[proj["geocode_method"] == "not_geocoded", "geocode_failure_reason"].value_counts().to_dict()
    s["state_mismatch"] = int(proj["geo_state_matches_org_state"].eq(False).sum())
    s["county_assignment"] = proj["county_assignment"].fillna("(none)").value_counts().to_dict()
    s["geocoded_without_county"] = int((proj["geocode_method"].isin([GEO_ZCTA, GEO_REPORTER])
                                        & proj["county_fips"].isna()).sum())
    # Connecticut: 2024 vintage = 9 planning regions (09110-09190); legacy counties 09001-09015
    ct = proj["county_fips"].fillna("").str.startswith("09")
    s["ct"] = {"rows_org_state_ct": int(proj["org_state"].eq("CT").sum()), "rows_with_ct_county": int(ct.sum()),
               "planning_region_rows": int(proj["county_fips"].fillna("").str.fullmatch(r"091[1-9]0").sum()),
               "legacy_county_rows": int(proj["county_fips"].fillna("").str.fullmatch(r"090(0[1-9]|1[0-5])").sum()),
               "regions": proj.loc[ct, "county_fips"].value_counts().sort_index().to_dict(),
               "by_assignment": proj.loc[ct, "county_assignment"].fillna("(none)").value_counts().to_dict()}
    fail = proj[proj["zip_geocode_failed"].astype(bool)].copy()
    fail["recovered_by_reporter_point"] = fail["geocode_method"] == GEO_REPORTER
    s["geocode_fail_us"] = (fail.groupby(["org_name", "org_city", "org_state", "org_zipcode", "recovered_by_reporter_point"],
                                         dropna=False).size().reset_index(name="projects")
                            .sort_values("projects", ascending=False))
    s["no_zip_orgs"] = (proj[proj["geocode_failure_reason"].fillna("").str.startswith("no_zip")]
                        .groupby("org_name", dropna=False).size().sort_values(ascending=False))
    z = proj["geocode_method"] == GEO_ZCTA
    has_both = z & proj["reporter_lat"].notna()
    if has_both.any():
        d = _haversine_km(proj.loc[has_both, "lat"], proj.loc[has_both, "lon"],
                          proj.loc[has_both, "reporter_lat"], proj.loc[has_both, "reporter_lon"])
        s["zcta_vs_reporter_km"] = {"n": int(has_both.sum()), "median": float(np.median(d)),
                                    "p90": float(np.percentile(d, 90)), "max": float(d.max()),
                                    "n_over_25km": int((d > 25).sum())}
        gap = proj.loc[has_both, ["org_name", "org_state", "zip5"]].assign(km=d)
        s["zcta_gap_top"] = (gap.sort_values("km", ascending=False).drop_duplicates("org_name").head(6)
                             .assign(km=lambda x: x["km"].round(0).astype(int)))
    # missingness
    key = ["project_title", "abstract_text", "terms", "pi_names", "contact_pi_name", "org_name", "org_city",
           "org_state", "org_zipcode", "org_country", "org_dept_type", "org_type_name", "award_amount",
           "direct_cost_amt", "indirect_cost_amt", "activity_code", "admin_ic_abbreviation", "funding_ics",
           "project_start_date", "project_end_date", "budget_start", "budget_end", "is_active",
           "spending_categories_desc", "reporter_lat", "lat", "county_fips"]
    miss = []
    for k in key:
        col = proj[k]
        n = int(col.isna().sum() + (col.astype(str).str.strip() == "").sum() * (col.dtype == object))
        miss.append({"variable": k, "missing": n, "pct": f"{100 * n / len(proj):.1f}%"})
    s["missingness"] = pd.DataFrame(miss)
    s["missing_award_by_agency"] = proj[proj["award_amount"].isna()]["agency_code"].fillna("(missing)").value_counts().to_dict()
    s["missing_award_subprojects"] = int((proj["award_amount"].isna() & proj["is_subproject"]).sum())
    return s


def _haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = (np.radians(np.asarray(x, dtype=float)) for x in (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def write_audit(s: dict, qmeta: dict, retrieved_at: str) -> Path:
    per = s["per_condition_main"]
    manifest = json.loads((raw_dir(SOURCE_ID) / "MANIFEST.json").read_text())
    files = "\n".join(f"| `{k}` | {v['bytes']:,} | `{v['sha256'][:16]}...` | {v['url']} |"
                      for k, v in sorted(manifest["files"].items()))
    fail = s["geocode_fail_us"]
    zc = s.get("zcta_vs_reporter_km", {})
    ps = s["parent_sub"]
    text = f"""# DATA AUDIT — NIH RePORTER (Project API v2)

| Field | Value |
|---|---|
| source_id | {SOURCE_ID} |
| Source (dataset/API name, exact files/endpoints) | NIH RePORTER Project API v2, `POST {API_URL}` (one request per page, limit {PAGE_LIMIT}, sorted by appl_id asc). Exact criteria JSON for all {s['n_queries']} queries in `queries.json`. |
| Publishing organization | National Institutes of Health, Office of Extramural Research (RePORT); records include non-NIH HHS/VA agencies that report to RePORTER |
| Retrieval date (UTC) | {qmeta.get('retrieved_at_first')} to {qmeta.get('retrieved_at_last')} (HTTP responses cached under data/_http_cache; reruns reuse them) |
| Source version / release | RePORTER API v2.0 (swagger `info.version`); no dataset release number. Latest `date_added` among retrieved records: {s['max_date_added']} |
| Source update date / cadence | RePORTER is refreshed weekly by NIH. FY2026 was still open at retrieval (fiscal year ends 2026-09-30), so FY2026 counts are incomplete. |
| License / access conditions | U.S. federal government data, public; no key. NIH asks for <= 1 request/second and large jobs off-hours; this module throttles to 1 request/second ({s['n_requests']} requests in total). |
| Unit of observation | One fiscal-year award record (`appl_id`): an application funded in one fiscal year. A multi-year project (`core_project_num`) has one appl_id per fiscal year; multi-component awards add subproject records with their own appl_ids. |
| Sample size (actual, as ingested) | {s['n_projects']:,} unique appl_ids ({s['n_core_projects']:,} core projects, {s['n_orgs']:,} organizations, {s['n_subprojects']:,} subproject records) from {s['n_condition_queries']} condition-term queries and {s['n_augmented_queries']} augmented queries, FY{FISCAL_YEARS[0]}-FY{FISCAL_YEARS[-1]} |
| Geography (resolution, vintage) | Awardee organization address (city/state/ZIP), U.S. only. ZIP5 -> 2024 Census ZCTA internal point (`zcta_centroid`) via `crosswalk.zip_to_geo`; if the ZIP has no ZCTA, RePORTER's own org coordinate verified to lie in the org's state (`point`). County = 2024 county containing the point. RePORTER's coordinate is always kept as `reporter_lat/lon`. |
| Person-level? | no (PI names are public award metadata about investigators, not study participants) |
| Geographic? | yes, facility/organization points (layer 4), not population burden |
| Omics? | no |
| Wearable? | no (grant text may mention wearables; no measurements) |
| True participant linkage across modalities? | no; not applicable (no participant data) |

## Files / endpoints retrieved
All files are in `data/raw/{SOURCE_ID}/`; hashes in `MANIFEST.json`.

| File | Bytes | sha256 | Source |
|---|---|---|---|
{files}

* `records.jsonl.gz`: every unique RePORTER record returned (full API JSON, one per appl_id).
* `query_hits.csv`: query_id x appl_id membership for every query.
* `queries.json`: for each query, the exact request body (`criteria` with `fiscal_years`, `advanced_text_search` {{operator: "{OPERATOR}", search_field: "{SEARCH_FIELD}", search_text}}), API total, RePORTER `search_id` and web URL, and a log of every page request.

Query design:
* Condition queries: one per `search_terms` entry in `configs/conditions.yaml`; search text is the exact quoted phrase, for example `"postural orthostatic tachycardia syndrome"`.
* Augmented queries (SPEC section D), for each demo-cluster condition ({', '.join(load_config('conditions')['demo_cluster']['members'])}): `("term 1" OR "term 2" ...) AND (<measurement variants>)`. Measurement variants: {'; '.join(f"{k} = {' OR '.join(v)}" for k, v in AUGMENT_TERMS.items())}. Plural variants are there because the API does no stemming (measured: "wearable" 5,902 hits, "wearables" 2,728, the OR of both 5,964; "biomarker" 135,407, "biomarkers" 82,679, both 143,201; FY2015-2026, all projects).
* Fiscal-year splitting: needed only when a criteria set returns more than {MAX_RETRIEVABLE:,} records (offset cap {MAX_OFFSET:,}). Split used in this run: {s['any_fy_split']}.
* Queries with zero API hits: {', '.join(s['zero_hit_queries']) or 'none'}.
* API total vs retrieved mismatches: {s['query_warnings'] or 'none'}.

## Key variables
| Variable | Meaning | Units | Downstream use |
|---|---|---|---|
| appl_id | RePORTER application id (one fiscal-year award record) | id | primary key of nih_projects |
| project_num / core_project_num | full award number / award number without type, support year and suffix | id | groups fiscal years of one project; parent/subproject grouping uses (project_num, fiscal_year) |
| subproject_id, is_subproject | component of a multi-project award | id / bool | totals count the parent's award_amount instead of its subprojects' (measured consistency below) |
| fiscal_year | federal fiscal year of the award | year | per-FY counts |
| project_title, abstract_text, public_health_relevance, terms, pref_terms | text; `terms` is RePORTER's thesaurus-expanded concept list | text | measurement text mining; match_tier |
| pi_names, pi_profile_ids, contact_pi_name, contact_pi_profile_id, principal_investigators (JSON incl. is_contact_pi) | investigators | text / id | investigator counts per organization |
| org_name, org_ipf_code, org_city, org_state, org_zipcode, org_country, org_dept_type, org_type_name | awardee organization | text | research-site aggregation |
| award_amount | total-cost award obligation of this appl_id; equals the sum of `ic_fundings[].total_cost` in {s['award_eq_ic_sum'][0]:,} of {s['award_eq_ic_sum'][1]:,} records that have both | USD | funding totals (de-duplicated) |
| direct_cost_amt, indirect_cost_amt | direct and indirect cost fields as published; they do **not** sum to award_amount in {s['award_ne_direct_indirect'][0]:,} of {s['award_ne_direct_indirect'][1]:,} records that have all three | USD | context only; never summed in place of award_amount |
| agency_code, admin_ic_abbreviation, funding_ics, ic_fundings (JSON) | awarding agency, administering IC, funding ICs with per-IC costs | code | agency breakdown |
| activity_code, funding_mechanism, award_type | mechanism (R01, U54, SBIR...) | code | capability type |
| project_start_date/end_date, budget_start/end | project and budget periods | date | recency |
| is_active | RePORTER active flag at retrieval | bool | active projects |
| spending_categories(_desc) | NIH RCDC spending categories | id / names | agreement check against text search |
| zip5, lat, lon, county_fips, state_fips, geocode_method | ZIP5 -> ZCTA internal point -> 2024 county (`zcta_centroid`); RePORTER's state-verified point when the ZIP has no ZCTA (`point`) | degrees / FIPS (zero-padded strings) | facility layer geography |
| reporter_lat, reporter_lon | RePORTER's own organization coordinate, always kept | degrees | canonical point only for the {s['n_reporter_point']:,} rows whose ZIP has no ZCTA; otherwise kept for comparison |
| matched_condition_ids, augmented_terms_hit | every condition / measurement term whose query returned the record, including terms-field-only and likely-false-positive links | text | convenience only; use `nih_project_conditions` or the helpers for the precision view |
| match_tier (links) | where the matched phrase occurs locally: title_abstract / terms_only / not_found_locally | category | precision filter |
| rcdc_category_assigned (links) | whether NIH's RCDC category for the condition is on the record | bool / null | precision check |

## Sample sizes

### Projects per condition
`projects (appl_id)` counts fiscal-year award records; `core projects` counts distinct multi-year projects.
`award obligations USD (dedup)` counts each appl_id once and excludes subprojects whose parent award is also in the condition set. `naive sum USD` is the plain sum and shows how much double counting that removes.
These are **award obligations** (total cost, summed over fiscal years 2015-2026), not expenditure, and not the portion of a grant spent on the condition: a grant that mentions the condition anywhere contributes its whole award.

{_md_table(per)}

**Precision view** (what `projects_for_condition` / `research_centers_for_condition` return by default): the condition phrase is in the title, abstract or narrative, and ambiguous-acronym hits without context are excluded.

{_md_table(s['per_condition_precision'])}

Why the precision view is the default, measured: among NIH-administered records of conditions that have an NIH RCDC spending category, the share carrying that category by match tier. Only fiscal years in which that category appears on at least one retrieved record are counted ({s['rcdc_years']}). RCDC categories are assigned after a fiscal year closes: {s['rcdc_fy2026'][0]} of {s['rcdc_fy2026'][1]} NIH FY2026 records have any category, and {s['rcdc_non_nih'][0]} of {s['rcdc_non_nih'][1]} non-NIH records do.

{_md_table(s['rcdc_by_tier'])}

Terms-field-only matches carry the NIH category less often than title/abstract matches for: {', '.join(s['rcdc_lower'])}; not for: {', '.join(s['rcdc_not_lower']) or 'none'} (ptlds is compared against the broader "Lyme Disease" category, so any Lyme grant passes). RCDC is NIH's own classification, not ground truth: a category is assigned when a project is judged relevant enough to report as spending, so even title/abstract matches do not all carry it.

* Unique appl_ids across all conditions (union): {s['n_projects']:,}. Union award obligations (dedup) = ${s['union_obligations']:,.0f} (naive sum ${s['union_naive_sum']:,.0f}; {s['union_sub_overlap']} subproject rows had their parent in the table).
* Parent vs subproject amounts, measured over the {ps['groups']} (project_num, fiscal year) groups in the table that hold a parent and at least one subproject: parent = sum of the retrieved subprojects in {ps['equal']}; parent larger in {ps['parent_larger']} (its other subprojects did not match any query); parent **smaller** in {ps['parent_smaller']} ({', '.join(ps['parent_smaller_groups'])}), where the subprojects sum to ${ps['parent_smaller_excess']:,.0f} more than their parents in total. RePORTER documents the parent amount as including all subprojects, and the parent amount always equals the award's IC funding, so the parent amount is what the totals count. A condition matched only by subprojects of such an award therefore carries the subprojects' larger sum.
* {s['multi_condition_appl_ids']:,} appl_ids matched more than one condition. The per-condition totals above therefore overlap: their sum (${s['sum_of_condition_obligations']:,.0f}) is **not** a total and must not be reported as one. Use the union figure.
* Active: {_pct(s['n_active'], s['n_projects'])} appl_ids carry `is_active = true`.
* Grouping condition `post_infectious_syndrome`: all 3 of its search phrases returned 0 records.

### Projects per condition and fiscal year (unique appl_ids)
{_md_table(s['per_condition_fy'])}

All matched projects by fiscal year: {', '.join(f'{k}: {v}' for k, v in s['overall_fy'].items())}.

### Condition-term queries (API total = records returned by RePORTER for that phrase, FY2015-2026)
{_md_table(s['per_term'])}

Note: several phrases return near-identical totals (e.g. "myalgic encephalomyelitis" vs "chronic fatigue syndrome", "long haul COVID" vs "post COVID syndrome"). The `terms` field carries thesaurus synonyms, so a phrase matches through the concept even when it is absent from the title and abstract. `match_tier` makes this visible; `terms-field only` rows are concept-level matches and lower precision.

### Measurement-augmented queries (unique appl_ids)
* `api_hits`: records RePORTER returned for `(condition terms) AND (measurement variants)`.
* `both_in_title_abstract`: the condition phrase and the measurement term both occur in the title, abstract or narrative, not only in the terms field.
* `local_regex_hits`: of the condition's own hits, how many contain a measurement variant anywhere (title, abstract, narrative or terms), by local regex. This is an independent check of the API's boolean logic.
* `api_only` / `local_only`: disagreements between the two.

A hit means a grant mentions the measurement term. It is not evidence that the measurement works for the condition.

{_md_table(s['augmented'])}

Augmented hits that were not also returned by a condition-term query for the same condition: {s['aug_not_in_condition_set']}.

### Ambiguous-acronym check
{_md_table(pd.DataFrame(s['acronym_checks'])) if s['acronym_checks'] else 'none'}

`no_context` rows are kept in `nih_project_conditions` with `likely_false_positive = true` and are excluded from the helpers by default.

## Missingness
Measured over the {s['n_projects']:,} nih_projects rows (null or empty string).

{_md_table(s['missingness'])}

* Missing award_amount by agency: {s['missing_award_by_agency']}; {s['missing_award_subprojects']} of the missing are subproject records.
* Agencies: {s['agency']}.
* Organization country: {s['country']}.

## Geocoding
Method, in order:
1. U.S. organizations only. Foreign postal codes are never read as ZIPs; a Nigerian 900246 would otherwise land on ZIP 90024 in Los Angeles.
2. ZIP5 to the Census ZCTA internal point with `crosswalk.zip_to_geo`: `geocode_method = {GEO_ZCTA}`, resolution `zcta_centroid`.
3. Only if the ZIP has no ZCTA (unique or PO-box ZIPs of large campuses): RePORTER's own organization coordinate, accepted only when it falls inside a 2024 county of the organization's stated state. `geocode_method = {GEO_REPORTER}`, resolution `point`.

Results:
* {_pct(s['n_geocoded'], s['n_projects'])} projects geocoded: {s['n_zcta']:,} by ZCTA and {s['n_reporter_point']:,} by the verified RePORTER point. That is {_pct(s['n_us_geocoded'], s['n_us'])} of U.S. organizations.
* Not geocoded, by reason: {s['fail_reasons']}.
* {s['n_nonus']} projects are at non-U.S. organizations.
* County assignment: {s['county_assignment']}. `zcta_county_rel2020_max_land_share` covers ZCTAs such as 60611 (Chicago), whose internal point lies in Lake Michigan outside the clipped 5m county polygons, so `zcta_centroids` has no county for them. Geocoded rows still without a county: {s['geocoded_without_county']}.
* ZCTA-geocoded rows whose ZCTA county lies in a different state than `org_state`: {s['state_mismatch']} (flag `geo_state_matches_org_state`).
* Connecticut (county vintage 2024 = 9 planning regions, 09110-09190): {s['ct']['rows_org_state_ct']} rows have org_state CT and {s['ct']['rows_with_ct_county']} carry a CT county FIPS. {s['ct']['planning_region_rows']} of those are 2024 planning regions and {s['ct']['legacy_county_rows']} are legacy counties (09001-09015). Regions: {s['ct']['regions']}; assignment: {s['ct']['by_assignment']}. The 2020 ZCTA-county fallback is never applied to CT, because its 2020 counties are not 2024 counties.
* Records without any address: {len(s['no_zip_orgs'])} organizations. These are mostly NIH intramural institutes, which RePORTER lists without city or ZIP. They are left ungeocoded rather than placed in Bethesda by assumption. Top: {dict(list(s['no_zip_orgs'].head(8).items()))}.
* U.S. organization ZIPs that failed `zip_to_geo` (the ZIP has no ZCTA), and whether the RePORTER point recovered them:

{_md_table(fail) if len(fail) else 'none'}

* Distance between the ZCTA internal point and RePORTER's own coordinate, for ZCTA-geocoded rows: n = {zc.get('n')}, median {zc.get('median', float('nan')):.1f} km, 90th percentile {zc.get('p90', float('nan')):.1f} km, max {zc.get('max', float('nan')):.1f} km, {zc.get('n_over_25km')} over 25 km. Largest gaps (organization, listed state and ZIP, km): {'; '.join(f"{r.org_name} ({r.org_state} {r.zip5}) {r.km}" for r in s['zcta_gap_top'].itertuples()) if 'zcta_gap_top' in s else 'n/a'}. In these cases RePORTER's point lies in another state from the ZIP on the record, which looks like a stale RePORTER coordinate, so the listed ZIP is primary.

## Linkage strategy
* `nih_projects.appl_id` joins `nih_project_conditions.appl_id` and `nih_project_query_hits.appl_id`.
* `condition_id` joins the ontology layer's condition registry (configs/conditions.yaml ids).
* `county_fips` / `state_fips` / `lat`,`lon` place the **organization** on the map (facility layer). Joining these to geographic burden rows is an ecological, place-level association: it does not mean the grant studies people from that county.
* `org_name` / `org_ipf_code` can be matched to ClinicalTrials.gov sites and NPPES organizations only by fuzzy name + location matching. No shared identifier exists; any such match must be labelled approximate.
* Not joinable: participants (no person-level data); PIs to NPPES NPIs (no shared id).

## Limitations and caveats
* Text-search recall and precision: condition assignment is RePORTER's phrase search over title, abstract and the thesaurus-expanded `terms` field. It includes projects that only mention the condition (e.g. as a comorbidity or exclusion), and misses projects that use other wording. `match_tier` and `rcdc_category_assigned` are measured quality signals; neither is ground truth.
* RCDC categories (measured): only NIH records carry them ({s['rcdc_non_nih'][0]} of {s['rcdc_non_nih'][1]} non-NIH records); none of the {s['rcdc_fy2026'][1]} NIH FY2026 records has one yet; category year ranges observed: {s['rcdc_years']}. dysautonomia, eds_hsd, mcas, ibs, gastroparesis and post_infectious_syndrome have no condition-specific category, so there is no agreement check for them. The `RCDC assigned` column counts over all fiscal years, including years without the category.
* Funding: `award_amount` is an obligation (total cost) for one fiscal year, not spending; a whole award counts toward every condition it mentions; parent/subproject double counting is removed by `obligation_total`; conditions overlap, so per-condition totals must not be added up.
* FY2026 is incomplete at retrieval.
* Non-NIH agencies that report to RePORTER are included and flagged by `agency_code` ({s['agency']}). VA records carry no award amount ({s['missing_award_by_agency'].get('VA', 0)} of {s['agency'].get('VA', 0)} missing), so VA-heavy organizations show unknown rather than $0 obligations.
* Organization location is the awardee's address. Multi-site studies and subaward sites are not represented.
* Research readiness is not patient burden: counts of grants reflect where research institutions are, which follows population and academic density.

## Processed outputs
| Table | Rows |
|---|---|
| nih_projects | {s['n_projects']:,} |
| nih_project_conditions | {s['n_links']:,} |
| nih_project_query_hits | {s['n_query_hits']:,} |

Helper functions for the MCP layer: `projects_for_condition(condition_id)` and `research_centers_for_condition(condition_id)` in `{MODULE}`.

## Reproduce
`uv run python -m {MODULE}` (uses the HTTP cache; `--refresh` re-queries the API; `--offline` rebuilds tables from the raw files only).
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


def write_registry(s: dict, retrieved_at: str, status: str = "ingested") -> Path:
    return write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "NIH RePORTER Project API v2",
        "publisher": "National Institutes of Health, Office of Extramural Research (RePORT)",
        "landing_url": LANDING_URL,
        "access_urls": [API_URL] + list(DOC_URLS.values()),
        "license": "U.S. federal government data (public domain)",
        "access_conditions": "open API, no key; NIH asks for <= 1 request/second and large jobs off-hours",
        "retrieved_at": retrieved_at,
        "source_version": f"RePORTER Project API v2.0; latest date_added {s['max_date_added']}",
        "update_date": "weekly (RePORTER refresh); FY2026 incomplete at retrieval",
        "data_layer": "facility",
        "unit_of_observation": "fiscal-year award record (appl_id) matched by text search to a target condition",
        "sample_size": {
            "projects_appl_id": s["n_projects"], "core_projects": s["n_core_projects"],
            "organizations": s["n_orgs"], "condition_links": s["n_links"], "augmented_query_hits": s["n_query_hits"],
            "queries": s["n_queries"], "fiscal_years": f"{FISCAL_YEARS[0]}-{FISCAL_YEARS[-1]}",
            "per_condition_projects": {r["condition_id"]: int(r["projects (appl_id)"])
                                       for _, r in s["per_condition"].iterrows()},
            "geocoded_projects": s["n_geocoded"],
        },
        "geographic_resolution": "zcta_centroid (organization ZIP5 as ZCTA internal point); point (RePORTER org "
                                 "coordinate, state-verified) when the ZIP has no ZCTA; organization city/state",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (award records, no participants)",
        "true_participant_linkage_across_modalities": False,
        "status": status,
        "processed_outputs": ["nih_projects", "nih_project_conditions", "nih_project_query_hits"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": MODULE,
        "limitations": [
            "Condition assignment is RePORTER phrase search over title/abstract/thesaurus-expanded terms; "
            "terms-only matches are concept-level (see match_tier)",
            "award_amount is a fiscal-year total-cost obligation, not spending; whole awards count toward every "
            "condition they mention; per-condition totals overlap and must not be summed",
            "Multi-project awards: obligation_total counts the parent amount (= IC funding) in place of its "
            "subprojects per (project_num, fiscal_year); in "
            f"{s['parent_sub']['parent_smaller']} of {s['parent_sub']['groups']} parent+subproject groups the "
            "subprojects sum to more than the parent (RePORTER inconsistency)",
            "FY2026 incomplete at retrieval",
            "Organization location = awardee address (ZIP5 -> ZCTA internal point, or RePORTER's state-verified "
            "point when the ZIP has no ZCTA); not study sites; NIH intramural records have no address",
            "Ambiguous acronym PASC also means pancreatic stellate cells; flagged likely_false_positive",
            "RCDC categories: NIH-only, not yet assigned for FY2026",
            "Research-capability signal only; not patient burden or evidence that a measurement works",
        ],
    })


# --------------------------------------------------------------------------- entry point
def run(refresh: bool = False, offline: bool = False) -> dict:
    if offline:
        records, hits, qmeta = load_raw()
        retrieved_at = qmeta["retrieved_at_last"]
    else:
        fetch_docs()
        records, hits, qlog = fetch_all(all_queries(), refresh=refresh)
        retrieved_at = save_raw(records, hits, qlog)
        _, _, qmeta = load_raw()
    proj, pc, qh = build_tables(records, hits, qmeta, retrieved_at)
    write_table(proj, "nih_projects", producer=MODULE,
                description="NIH RePORTER fiscal-year award records (one row per appl_id) matched to target conditions")
    write_table(pc, "nih_project_conditions", producer=MODULE,
                description="appl_id x condition_id x matched RePORTER condition-term query, with local match tier")
    write_table(qh, "nih_project_query_hits", producer=MODULE,
                description="appl_id x (demo condition AND measurement term) RePORTER query hits")
    s = compute_stats(proj, pc, qh, qmeta)
    s["n_links"], s["n_query_hits"] = int(len(pc)), int(len(qh))
    write_audit(s, qmeta, retrieved_at)
    write_registry(s, retrieved_at)
    print(f"nih_projects={len(proj)} nih_project_conditions={len(pc)} nih_project_query_hits={len(qh)}")
    return s


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refresh", action="store_true", help="bypass the HTTP cache and re-query RePORTER")
    ap.add_argument("--offline", action="store_true", help="rebuild tables from data/raw files only")
    a = ap.parse_args()
    run(refresh=a.refresh, offline=a.offline)
