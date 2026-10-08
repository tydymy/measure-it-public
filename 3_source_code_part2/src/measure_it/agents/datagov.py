"""Data.gov / U.S. open-data discovery agent (SPEC section H).

Data.gov is a METADATA catalog. It indexes dataset records harvested from
agency catalogs (data.json / DCAT-US) and points to the agency hosts that
actually serve the records. This module searches that catalog, normalises
each hit to a small routing record, and never treats Data.gov as the host of
the underlying data: every result carries ``retrieve_from``, the agency URL
the data should be retrieved from.

Current API (verified 2026-09-23; the legacy CKAN endpoint is gone)
-------------------------------------------------------------------
* Legacy CKAN ``https://catalog.data.gov/api/3/action/package_search`` answers
  HTTP 404 ``{"detail":{},"message":"Not Found"}``; the old gateway
  ``https://api.gsa.gov/technology/datagov/v3/action/package_search`` 301-redirects
  to that same 404.
* Replacement (Flask + OpenSearch app, DCAT-US metadata):
    - origin, keyless:  ``GET https://catalog.data.gov/search``
      (OpenAPI spec at ``https://catalog.data.gov/openapi.json``;
      robots.txt asks for ``Crawl-Delay: 10``, which this module honours)
    - api.data.gov gateway: ``GET https://api.gsa.gov/technology/datagov/v4/search``
      with header ``X-Api-Key`` (``DEMO_KEY`` works; documented limits are
      30 requests/IP/hour and 50/IP/day; the gateway reported
      ``X-RateLimit-Limit: 10`` on our DEMO_KEY calls). A personal key is read
      from $DATA_GOV_API_KEY / $API_DATA_GOV_KEY / $DATAGOV_API_KEY if set.
  Both return the same ordered results for the same query (the bytes differ only
  in OpenSearch ``_score`` values between replicas). Parameters used here:
  ``q`` (AND over terms; ``OR`` supported), ``per_page`` (1-1000), ``org_type``.
  The response has ``results``, ``after`` (cursor) and ``sort``; it does NOT
  report a total hit count, so ``n_results`` below is "results returned on the
  first page, capped at rows".

Outputs
-------
data/processed/data_gov_discovery_log.parquet   one row per executed query
data/processed/data_gov_search_results.parquet  one row per (query, result rank)
data/raw/data_gov_catalog/                      raw responses, endpoint probes,
                                                MANIFEST.json, DATA_AUDIT.md,
                                                registry_entry.yaml

Reproduce: ``uv run python -m measure_it.agents.datagov``
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import time
from collections import Counter
from urllib.parse import urlparse

import pandas as pd

from .. import http
from ..config import UNKNOWN, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..provenance import add_provenance
from ..registry import load_source_registry, write_registry_entry
from ..store import write_table

SOURCE_ID = "data_gov_catalog"
PRODUCER = "measure_it.agents.datagov"

CATALOG_HOST = "https://catalog.data.gov"
ORIGIN_SEARCH_URL = f"{CATALOG_HOST}/search"
GATEWAY_BASE = "https://api.gsa.gov/technology/datagov/v4"
GATEWAY_SEARCH_URL = f"{GATEWAY_BASE}/search"
OPENAPI_URL = f"{CATALOG_HOST}/openapi.json"
ORGANIZATIONS_URL = f"{CATALOG_HOST}/api/organizations"
STATS_URL = f"{CATALOG_HOST}/api/stats"
ROBOTS_URL = f"{CATALOG_HOST}/robots.txt"
DOCS_URL = "https://resources.data.gov/catalog-api/"
LEGACY_CKAN_URL = f"{CATALOG_HOST}/api/3/action/package_search"
LEGACY_GATEWAY_V3_URL = "https://api.gsa.gov/technology/datagov/v3/action/package_search"
DATASET_PAGE = f"{CATALOG_HOST}/dataset/{{slug}}"
BULK_DUMP_URL = "https://filestore.data.gov/gsa/catalog/jsonl/dataset.jsonl.gz"

API_KEY_ENV_VARS = ("DATA_GOV_API_KEY", "API_DATA_GOV_KEY", "DATAGOV_API_KEY")
DEMO_KEY = "DEMO_KEY"
ORIGIN_CRAWL_DELAY_S = 10.0   # catalog.data.gov/robots.txt: "Crawl-Delay: 10"
GATEWAY_MIN_INTERVAL_S = 1.0
DEFAULT_ROWS = 20
MAX_ROWS = 1000

# Hosts that are the catalog itself, never the agency that holds the records.
CATALOG_HOSTS = {"data.gov", "www.data.gov", "catalog.data.gov", "api.gsa.gov"}
# Secondary (department-level) catalogs whose dataset pages are often metadata stubs for records
# that another agency hosts (e.g. healthdata.gov/d/<id> for FDA 510(k)/PMA, NCBI PubChem patents,
# NIH RePORTER). A landing page there is only used for routing when no distribution URL points
# at a host outside the catalogs (HHS-hosted datasets keep their healthdata.gov URL).
SECONDARY_CATALOG_HOSTS = {"healthdata.gov", "www.healthdata.gov"}

DEMO_QUERIES = [
    "long COVID",
    "health center service delivery sites",
    "social vulnerability index",
    "Lyme disease",
    "PLACES local data for better health",
    "medical device",
]

RESULT_FIELDS = [
    "is_unknown", "result_status", "query", "rank", "title", "agency", "publisher",
    "publisher_hierarchy", "parent_organization", "organization_type", "organization_slug",
    "description", "access_level", "resources", "resource_urls", "resource_formats",
    "modified", "issued", "identifier", "landing_page", "retrieve_from", "retrieve_from_basis",
    "retrieve_from_host", "machine_readable_url", "license", "keywords", "datagov_url",
    "harvest_record", "last_harvested_date", "api_endpoint", "request_url", "http_status",
    "retrieved_at", "more_available", "note",
]

# ---------------------------------------------------------------------------
# Curated discovery rules for the planned sources in SOURCE_REGISTRY.yaml.
#
# `queries` are tried in order until one finds the planned dataset.
# A result is a RELATED match when its title matches `title` AND (when
# `publisher` is set) its publisher / publisher hierarchy / parent organisation
# matches `publisher` - i.e. a record of the same programme by the right agency.
# It is a PLANNED-dataset match when it is also matched by `planned` (defaults to
# `title`) and not by `exclude_planned`. Only a planned-dataset match counts as
# indexed = "yes"; related matches (derived indicators, locator tools, agency
# portals, blog-like records) are reported with match_scope = "related_dataset".
# Records whose title marks them as blog/news posts never match.
# Sources without a rule fall back to their registry name (see _rule_for).
# Regexes are case-insensitive unless they use a (?-i:...) group.
# ---------------------------------------------------------------------------
_CDC = r"Centers for Disease Control|\bCDC\b|National Center for Health Statistics|\bNCHS\b"
_NIH = r"National Institutes of Health|\bNIH\b"
_CMS = r"Centers for Medicare|\bCMS\b"
_HRSA = r"Health Resources and Services Administration|\bHRSA\b"
_CENSUS = r"Census"
_FDA = r"Food and Drug Administration|\bFDA\b"
NON_DATASET_TITLE = r"^\s*Blog\s*\|"
# `planned` for sources that are not a catalogued dataset (curated literature extractions): a topical hit is at most a
# related dataset, never "the dataset is indexed".
NEVER_PLANNED = r"(?!)"

DISCOVERY_RULES: dict[str, dict] = {
    "nhanes_2011_2014": dict(
        queries=["National Health and Nutrition Examination Survey", "NHANES 2013-2014", "NHANES 2011-2012"],
        title=r"NHANES|National Health and Nutrition Examination Survey", publisher=_CDC,
        # the survey itself or its 2011-2012 / 2013-2014 cycles, not indicators derived from it
        planned=r"^(National Health and Nutrition Examination Survey|NHANES)( \(NHANES\))?\s*$|(NHANES|Nutrition Examination Survey).*(2011|2013)"),
    "stanford_longcovid_wearables": dict(
        queries=["COVID-19 wearable", "wearable heart rate", "smartwatch"],
        title=r"wearable|smart ?watch|fitbit", publisher=r"Stanford"),
    "mapmecfs_nih_pi_mecfs": dict(
        queries=["myalgic encephalomyelitis", "chronic fatigue syndrome"],
        title=r"myalgic|chronic fatigue|ME/CFS|mapMECFS", publisher=_NIH + r"|NINDS|\bRTI\b"),
    "mondo_ontology": dict(
        queries=["Mondo Disease Ontology", "disease ontology"],
        title=r"\bMondo\b", publisher=None),
    "ebi_ols4": dict(
        queries=["Ontology Lookup Service"],
        title=r"Ontology Lookup Service|(?-i:\bOLS4?\b)", publisher=None),
    "open_targets": dict(
        queries=["Open Targets"],
        title=r"Open Targets", publisher=None),
    "gwas_catalog": dict(
        queries=["GWAS Catalog", "genome-wide association studies catalog"],
        title=r"GWAS Catalog|catalog of (published )?genome-wide association", publisher=None),
    "ncbi_geo_sra": dict(
        queries=["Gene Expression Omnibus", "Sequence Read Archive"],
        title=r"Gene Expression Omnibus|(?-i:\bGEO\b)|Sequence Read Archive|(?-i:\bSRA\b)",
        publisher=r"National Center for Biotechnology|\bNCBI\b|National Library of Medicine|\bNLM\b|" + _NIH,
        planned=r"^(GEO \()?Gene Expression Omnibus|^Sequence Read Archive|^(?-i:SRA)\b"),
    "clinicaltrials_gov": dict(
        queries=["ClinicalTrials.gov", "clinical trials registry"],
        title=r"ClinicalTrials\.gov", publisher=r"National Library of Medicine|\bNLM\b|" + _NIH,
        planned=r"^ClinicalTrials\.gov"),
    "nih_reporter": dict(
        queries=["NIH RePORTER", "NIH ExPORTER"],
        title=r"(?-i:RePORTER|ExPORTER)|Research Portfolio Online Reporting", publisher=_NIH),
    "openfda_device": dict(
        queries=["openFDA device", "510(k) premarket notification", "Product Classification medical device"],
        title=r"openFDA|510\(?k\)?|premarket|device classification|^Product Classification$|"
              r"registration (and|&) listing|medical device|\bGUDID\b|\bMAUDE\b|\bMDR\b",
        publisher=_FDA,
        # the device classification / 510(k) / PMA / registration & listing data behind openFDA device
        planned=r"openFDA.*device|510\(?k\)?|premarket (notification|approval)|(?-i:\bPMA\b)|"
                r"device classification|^Product Classification$|registration (and|&) listing"),
    "cdc_places": dict(
        queries=["PLACES Local Data for Better Health"],
        title=r"^PLACES\b|500 Cities", publisher=_CDC, planned=r"^PLACES\b"),
    "cdc_long_covid": dict(
        queries=["Long COVID Household Pulse Survey", "Post-COVID Conditions", "long COVID"],
        title=r"long covid|post-covid|post covid", publisher=_CDC + r"|" + _CENSUS),
    "cdc_lyme": dict(
        queries=["Lyme disease county", "Lyme disease"],
        title=r"Lyme", publisher=_CDC,
        planned=r"Lyme.*(geograph|county|counties)", exclude_planned=r"without geography"),
    "cdc_brfss": dict(
        queries=["Behavioral Risk Factor Surveillance System", "BRFSS"],
        title=r"Behavioral Risk Factor Surveillance|(?-i:\bBRFSS\b)", publisher=_CDC,
        planned=r"(Behavioral Risk Factor Surveillance|(?-i:\bBRFSS\b)).*(20(22|23)|annual|data)"),
    "cms_physician_service": dict(
        queries=["Medicare Physician & Other Practitioners - by Provider and Service"],
        title=r"Medicare Physician (&|and) Other Practitioners", publisher=_CMS,
        planned=r"Medicare Physician (&|and) Other Practitioners.*by Provider and Service"),
    "cms_partd_prescriber": dict(
        queries=["Medicare Part D Prescribers - by Provider and Drug"],
        title=r"Medicare Part D Prescribers", publisher=_CMS,
        planned=r"Medicare Part D Prescribers.*by Provider and Drug"),
    "cms_mmd": dict(
        queries=["Mapping Medicare Disparities", "Medicare disparities"],
        title=r"Mapping Medicare Disparities|(?-i:\bMMD\b)", publisher=_CMS),
    "cdc_svi": dict(
        queries=["Social Vulnerability Index", "CDC/ATSDR Social Vulnerability Index", "ATSDR SVI"],
        title=r"Social Vulnerability Index|(?-i:\bSVI\b)", publisher=_CDC + r"|ATSDR|Agency for Toxic Substances",
        # the index itself, not datasets stratified by it (e.g. COVID-19 deaths by county SVI)
        exclude_planned=r"covid|death|vaccin|dashboard|priority"),
    "census_acs": dict(
        queries=["American Community Survey 5-Year Estimates", "American Community Survey"],
        title=r"American Community Survey|(?-i:\bACS\b)", publisher=_CENSUS,
        planned=r"(American Community Survey|(?-i:\bACS\b)).*5-Year"),
    "census_geography": dict(
        queries=["cartographic boundary county",
                 "Cartographic Boundary Current County and Equivalent United States",
                 "Gazetteer counties", "ZCTA county relationship"],
        title=r"cartographic boundar|gazetteer|relationship file", publisher=_CENSUS,
        # national county boundaries, the Gazetteer, or the ZCTA-county relationship file
        # (not e.g. "119th Congressional District within Current County and Equivalent Entities")
        planned=r"Cartographic Boundary (File|Shapefile)[^,]*,\s*(Current )?County and Equivalent( Entities)? "
                r"for United States|Gazetteer.*Count|ZCTA.*County.*Relationship|Relationship.*ZCTA.*County"),
    "cms_icd10cm": dict(
        queries=["ICD-10-CM", "International Classification of Diseases Clinical Modification"],
        title=r"ICD-10-CM|ICD-10|International Classification of Diseases", publisher=r"Centers for Medicare|CMS|National Center for Health Statistics|NCHS|CDC"),
    "usda_rucc": dict(
        queries=["Rural-Urban Continuum Codes"],
        title=r"Rural[- ]Urban Continuum", publisher=r"Economic Research Service|\bERS\b|Agriculture|USDA"),
    "nppes": dict(
        queries=["National Plan and Provider Enumeration System", "National Provider Identifier", "NPPES"],
        title=r"NPPES|National Plan (and|&) Provider Enumeration|National Provider Identifier|(?-i:\bNPIs?\b)",
        publisher=_CMS,
        planned=r"NPPES|National Plan (and|&) Provider Enumeration|NPI (Registry|Downloadable|Files?)"),
    "nucc_taxonomy": dict(
        queries=["Health Care Provider Taxonomy", "provider taxonomy"],
        title=r"Health Care Provider Taxonomy|(?-i:\bNUCC\b)|Taxonomy Crosswalk",
        publisher=r"NUCC|National Uniform Claim|" + _CMS,
        planned=r"Health Care Provider Taxonomy|(?-i:\bNUCC\b)"),
    "hrsa_health_centers": dict(
        queries=["Health Center Service Delivery Sites", "HRSA health center", "health center look-alike sites"],
        title=r"health center|data\.hrsa\.gov", publisher=_HRSA,
        planned=r"service delivery|look-alike"),
    "hrsa_hpsa": dict(
        queries=["Health Professional Shortage Areas", "HPSA designation", "HPSA primary care"],
        title=r"Health Professional Shortage Area|(?-i:\bHPSAs?\b)|data\.hrsa\.gov", publisher=_HRSA,
        # the designation data themselves, not address-lookup tools built on them
        planned=r"Health Professional Shortage Area|(?-i:\bHPSAs?\b)",
        exclude_planned=r"^Find |by Address|Eligible|Bonus|data\.hrsa\.gov"),
    "data_gov_catalog": dict(
        queries=["Data.gov catalog API", "Data.gov CKAN API"],
        title=r"data\.gov", publisher=r"General Services Administration|\bGSA\b",
        planned=r"data\.gov.*(api|catalog)|catalog.*data\.gov"),
    # reference vocabularies added 2026-09-23 (measure_it.omics.reference_data); non-federal, not expected indexed
    "reactome": dict(
        queries=["Reactome"],
        title=r"\bReactome\b", publisher=None),
    "hgnc": dict(
        queries=["HGNC gene nomenclature", "HUGO Gene Nomenclature"],
        title=r"(?-i:\bHGNC\b)|HUGO Gene Nomenclature", publisher=None),
    # sources added 2026-09-24 (person-level replication cohort, GEO cohorts, literature extractions and the
    # journal-supplement datasets of the device and lab analyses). Only NHANES 2003-2006 and GEO are federal; the
    # journal/repository deposits and the curated literature tables are not expected to be indexed by Data.gov, and
    # their searches record that (as for reactome / hgnc). `title` catches the topic (a related dataset at most);
    # `planned` names the deposit itself, so a federal dataset on the same condition is never counted as "indexed".
    "nhanes_2003_2006": dict(
        queries=["NHANES 2003-2004", "NHANES 2005-2006"],
        title=r"NHANES|National Health and Nutrition Examination Survey", publisher=_CDC,
        planned=r"(NHANES|Nutrition Examination Survey).*(2003|2005)"),
    "geo_cohorts": dict(
        queries=["NCBI Gene Expression Omnibus", "long COVID gene expression"],
        title=r"Gene Expression Omnibus|(?-i:\bGEO\b)|gene expression|transcriptom",
        publisher=r"National Center for Biotechnology|\bNCBI\b|National Library of Medicine|\bNLM\b|" + _NIH,
        planned=r"^(GEO \()?Gene Expression Omnibus"),
    "published_device_evidence": dict(
        queries=["PubMed", "wearable device long COVID"],
        title=r"PubMed|Europe PMC|MEDLINE", publisher=r"National Library of Medicine|\bNLM\b|NCBI|" + _NIH,
        planned=NEVER_PLANNED),
    "published_lab_evidence": dict(
        queries=["MEDLINE", "biomarker chronic fatigue syndrome"],
        title=r"PubMed|Europe PMC|MEDLINE", publisher=r"National Library of Medicine|\bNLM\b|NCBI|" + _NIH,
        planned=NEVER_PLANNED),
    "charlton_lc_mecfs_cpet_source": dict(
        queries=["post-exertional malaise", "long COVID exercise"],
        title=r"MUSCLE-ME|post-exertional|exercise.*(long covid|ME/CFS|chronic fatigue)", publisher=None,
        planned=r"MUSCLE-ME"),
    "appelman_lc_pem_source": dict(
        queries=["long COVID metabolomics", "long COVID muscle"],
        title=r"(long covid|post-covid).*metabolom|metabolom.*(long covid|post-covid)|post-exertional",
        publisher=None, planned=r"(long covid|post-covid).*metabolom|metabolom.*(long covid|post-covid)"),
    "fm_thermography": dict(
        queries=["fibromyalgia thermography", "fibromyalgia"],
        title=r"fibromyalgia", publisher=None, planned=r"thermograph"),
    "endo_arg1_repod": dict(
        queries=["endometriosis"],
        title=r"endometriosis", publisher=None, planned=r"arginase"),
    "heds_hsd_olink_serum_cinquina2026": dict(
        queries=["Ehlers-Danlos syndrome", "hypermobility"],
        title=r"Ehlers[- ]Danlos|hypermobil", publisher=None, planned=r"Olink|proteom"),
    "klein2023_mylc_ml_table": dict(
        queries=["long COVID immune profiling", "long COVID cortisol"],
        title=r"(?-i:MY-LC)|long covid.*(immun|cortisol)", publisher=None, planned=r"(?-i:MY-LC)"),
}

# A query string belongs to one source only: the discovery log and result rows are keyed by (query_kind, query).
assert len({q for r in DISCOVERY_RULES.values() for q in r["queries"]}) == sum(
    len(r["queries"]) for r in DISCOVERY_RULES.values()), "duplicate Data.gov discovery query across sources"

# Best-first order of per-query index statuses (used for the per-source verdict).
STATUS_ORDER = {"dataset_indexed": 0, "related_dataset_only": 1, "publisher_indexed_dataset_not_matched": 2,
                "no_planned_source_matched": 3, "not_found": 4, "no_results": 5, "query_failed": 6}

_last_live_call: dict[str, float] = {}


# ---------------------------------------------------------------------------
# Low-level access
# ---------------------------------------------------------------------------
def api_key_from_env() -> str | None:
    for var in API_KEY_ENV_VARS:
        if os.environ.get(var):
            return os.environ[var]
    return None


def _polite_get(url: str, *, params: dict | None, headers: dict | None, min_interval: float,
                refresh: bool = False, **kw) -> http.CachedResponse:
    """measure_it.http.get plus a per-host delay between LIVE calls (cache hits are free)."""
    host = urlparse(url).netloc
    wait = _last_live_call.get(host, -1e12) + min_interval - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    resp = http.get(url, params=params, headers=headers, refresh=refresh, **kw)
    if not resp.from_cache:
        _last_live_call[host] = time.monotonic()
    return resp


def _endpoints(endpoint: str) -> list[dict]:
    key = api_key_from_env()
    origin = dict(name="catalog.data.gov/search (origin, keyless)", url=ORIGIN_SEARCH_URL,
                  headers=None, min_interval=ORIGIN_CRAWL_DELAY_S)
    gateway = dict(
        name=f"api.gsa.gov/technology/datagov/v4/search (gateway, {'personal key' if key else DEMO_KEY})",
        url=GATEWAY_SEARCH_URL, headers={"X-Api-Key": key or DEMO_KEY}, min_interval=GATEWAY_MIN_INTERVAL_S)
    if endpoint == "origin":
        return [origin]
    if endpoint == "gateway":
        return [gateway]
    if endpoint != "auto":
        raise ValueError(f"endpoint must be auto|origin|gateway, got {endpoint!r}")
    # A personal key makes the documented gateway the natural first choice; otherwise
    # the keyless origin (DEMO_KEY allows only a few dozen calls per hour).
    return [gateway, origin] if key else [origin, gateway]


def search_raw(query: str, rows: int = DEFAULT_ROWS, *, org_type: str | None = None,
               endpoint: str = "auto", refresh: bool = False) -> tuple[dict | None, dict]:
    """Call the catalog search. Returns (payload or None, call metadata incl. errors)."""
    if not query or not query.strip():
        raise ValueError("query must be a non-empty string")
    if not 1 <= int(rows) <= MAX_ROWS:
        raise ValueError(f"rows must be in [1, {MAX_ROWS}]")
    params = {"q": query.strip(), "per_page": int(rows)}
    if org_type:
        params["org_type"] = org_type
    errors = []
    for ep in _endpoints(endpoint):
        try:
            resp = _polite_get(ep["url"], params=params, headers=ep["headers"], min_interval=ep["min_interval"],
                               refresh=refresh, expect_json=True, max_retries=3, timeout=60)
            if resp.from_cache and resp.status >= 400:
                # The HTTP cache key ignores headers, and endpoint probes cache error responses
                # (e.g. the no-key 403). Never serve a cached error for a search: ask again live.
                resp = _polite_get(ep["url"], params=params, headers=ep["headers"],
                                   min_interval=ep["min_interval"], refresh=True, expect_json=True,
                                   max_retries=3, timeout=60)
        except Exception as exc:  # network failure, HTML error page, non-JSON body
            errors.append(f"{ep['name']}: {type(exc).__name__}: {str(exc)[:200]}")
            continue
        if resp.status >= 400:
            errors.append(f"{ep['name']}: HTTP {resp.status}: {resp.text[:200]}")
            continue
        payload = resp.json()
        if not isinstance(payload, dict) or "results" not in payload:
            errors.append(f"{ep['name']}: unexpected payload keys {list(payload)[:5] if isinstance(payload, dict) else type(payload)}")
            continue
        meta = dict(api_endpoint=ep["name"], request_url=resp.url, http_status=resp.status,
                    retrieved_at=resp.fetched_at, from_cache=resp.from_cache, errors=errors,
                    ratelimit_limit=resp.headers.get("x-ratelimit-limit") or resp.headers.get("X-RateLimit-Limit"),
                    params=params)
        return payload, meta
    return None, dict(api_endpoint=UNKNOWN, request_url=UNKNOWN, http_status=None,
                      retrieved_at=utc_now_iso(), from_cache=False, errors=errors, params=params)


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------
def _s(v) -> str:
    """A clean string or UNKNOWN."""
    if v is None:
        return UNKNOWN
    if isinstance(v, (list, tuple)):
        v = "; ".join(str(x) for x in v if x not in (None, ""))
    elif isinstance(v, dict):
        v = v.get("name") or v.get("@id") or json.dumps(v, sort_keys=True)
    v = str(v).strip()
    return v if v else UNKNOWN


def strip_html(text):
    """Agency descriptions sometimes carry HTML (<p>, <ul>); keep the text only."""
    if not isinstance(text, str):
        return text
    text = re.sub(r"<[^>]+>", " ", text)
    text = (text.replace("&nbsp;", " ").replace("&amp;", "&").replace("&lt;", "<")
            .replace("&gt;", ">").replace("&quot;", '"').replace("&#39;", "'"))
    return re.sub(r"\s+", " ", text).strip()


def _host(url: str) -> str:
    try:
        return urlparse(url).netloc.lower().split(":")[0]
    except Exception:
        return ""


def _is_http(url) -> bool:
    return isinstance(url, str) and url.lower().startswith(("http://", "https://"))


def _url_of(v):
    """A DCAT URL-valued field as a plain string. Some harvested records (e.g. data.va.gov) give
    landingPage as a Document object {"accessURL": ...} or a list; return the first http URL, else
    the value unchanged (a non-URL string or None)."""
    if isinstance(v, dict):
        for k in ("accessURL", "downloadURL", "@id", "url", "href"):
            if _is_http(v.get(k)):
                return v[k].strip()
        return None
    if isinstance(v, (list, tuple)):
        return next((u for u in (_url_of(x) for x in v) if _is_http(u)), None)
    return v.strip() if isinstance(v, str) else v


def is_catalog_url(url: str, own_hosts: frozenset | set = frozenset()) -> bool:
    """True for Data.gov catalog/gateway hosts, unless they are the publishing agency's own (GSA)."""
    return _host(url) in (CATALOG_HOSTS - set(own_hosts))


GSA_PUBLISHER = r"General Services Administration|\bGSA\b|^Data\.gov$"
GSA_OWN_HOSTS = frozenset({"data.gov", "www.data.gov"})


def _publisher_chain(pub) -> list[str]:
    """dcat publisher -> [publisher, its subOrganizationOf, ...] (DCAT-US nesting)."""
    chain, seen = [], 0
    while pub and seen < 10:
        if isinstance(pub, str):
            chain.append(pub.strip())
            break
        if not isinstance(pub, dict):
            break
        name = pub.get("name") or pub.get("foaf:name")
        if name:
            chain.append(str(name).strip())
        pub = pub.get("subOrganizationOf")
        seen += 1
    return [c for c in chain if c]


def normalise_resources(distribution) -> list[dict]:
    out = []
    for d in distribution or []:
        if not isinstance(d, dict):
            continue
        url, kind = d.get("downloadURL"), "downloadURL"
        if not _is_http(url):
            url, kind = d.get("accessURL"), "accessURL"
        if not _is_http(url):
            continue
        media = d.get("mediaType")
        out.append({
            "url": url.strip(),
            "format": _s(d.get("format") or media),
            "name": _s(d.get("title")),
            "media_type": _s(media),
            "url_kind": kind,
        })
    return out


_NON_MACHINE = {"text/html", UNKNOWN}
# Documents and images are not machine-readable data (a PDF user guide is not a data file).
_DOCUMENT_MEDIA = {"application/pdf", "application/msword",
                   "application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                   "application/vnd.ms-powerpoint",
                   "application/vnd.openxmlformats-officedocument.presentationml.presentation"}


def is_machine_readable_media(media_type: str) -> bool:
    mt = str(media_type or "").split(";")[0].strip().lower()
    return bool(mt) and mt not in _NON_MACHINE and mt != UNKNOWN.lower() and mt not in _DOCUMENT_MEDIA \
        and not mt.startswith(("image/", "video/", "audio/"))


def choose_retrieve_from(landing_page, resources: list[dict], identifier,
                         own_hosts: frozenset | set = frozenset()) -> tuple[str, str]:
    """The agency URL to retrieve the data from; never a Data.gov catalog page (for records that GSA
    itself publishes, data.gov/www.data.gov are the agency's own hosts and are allowed).

    A landing page on a secondary catalog (healthdata.gov) is a metadata stub when the record's
    distributions point at another agency host; the first such distribution is used instead."""
    landing_page = _url_of(landing_page)
    agency_res = [r for r in resources
                  if not is_catalog_url(r["url"], own_hosts) and _host(r["url"]) not in SECONDARY_CATALOG_HOSTS]
    if _is_http(landing_page) and not is_catalog_url(landing_page, own_hosts):
        if _host(landing_page) in SECONDARY_CATALOG_HOSTS and agency_res:
            r = agency_res[0]
            return r["url"], f"dcat.distribution.{r['url_kind']} (landingPage is a healthdata.gov catalog page)"
        return landing_page, "dcat.landingPage"
    for r in resources:
        if not is_catalog_url(r["url"], own_hosts):
            return r["url"], f"dcat.distribution.{r['url_kind']}"
    if _is_http(identifier) and not is_catalog_url(identifier, own_hosts):
        return identifier.strip(), "dcat.identifier"
    return UNKNOWN, "no agency URL in the catalog record"


def normalise_record(rec: dict, *, rank: int, query: str, meta: dict, more_available: bool) -> dict:
    dcat = rec.get("dcat") or {}
    org = rec.get("organization") or {}
    chain = _publisher_chain(dcat.get("publisher"))
    publisher = rec.get("publisher") or (chain[0] if chain else None)
    parent_org = org.get("name")
    distribution = dcat.get("distribution")
    if isinstance(distribution, dict):  # a single distribution object instead of a list
        distribution = [distribution]
    resources = normalise_resources(distribution)
    landing = _url_of(dcat.get("landingPage"))
    identifier = dcat.get("identifier") or rec.get("identifier")
    who = " | ".join([str(publisher or "")] + chain + [str(parent_org or "")])
    own = GSA_OWN_HOSTS if re.search(GSA_PUBLISHER, who, flags=re.I | re.M) else frozenset()
    retrieve_from, basis = choose_retrieve_from(landing, resources, identifier, own)
    machine = next((r["url"] for r in resources
                    if is_machine_readable_media(r["media_type"]) and not is_catalog_url(r["url"], own)), UNKNOWN)
    slug = rec.get("slug")
    hierarchy = chain + ([parent_org] if parent_org and parent_org not in chain else [])
    return {
        "is_unknown": False,
        "result_status": "ok",
        "query": query,
        "rank": rank,
        "title": _s(rec.get("title") or dcat.get("title")),
        "agency": _s(publisher or parent_org),
        "publisher": _s(publisher),
        "publisher_hierarchy": " > ".join(hierarchy) if hierarchy else UNKNOWN,
        "parent_organization": _s(parent_org),
        "organization_type": _s(org.get("organization_type")),
        "organization_slug": _s(org.get("slug")),
        "description": _s(strip_html(rec.get("description") or dcat.get("description"))),
        "access_level": _s(dcat.get("accessLevel")),
        "resources": resources,
        "resource_urls": [r["url"] for r in resources],
        "resource_formats": sorted({r["format"] for r in resources}),
        "modified": _s(dcat.get("modified")),
        "issued": _s(dcat.get("issued")),
        "identifier": _s(identifier),
        "landing_page": _s(landing),
        "retrieve_from": retrieve_from,
        "retrieve_from_basis": basis,
        "retrieve_from_host": _host(retrieve_from) or UNKNOWN,
        "machine_readable_url": machine,
        "license": _s(dcat.get("license")),
        "keywords": [str(k) for k in (rec.get("keyword") or dcat.get("keyword") or [])],
        "datagov_url": DATASET_PAGE.format(slug=slug) if slug else UNKNOWN,
        "harvest_record": _s(rec.get("harvest_record")),
        "last_harvested_date": _s(rec.get("last_harvested_date")),
        "api_endpoint": meta.get("api_endpoint", UNKNOWN),
        "request_url": meta.get("request_url", UNKNOWN),
        "http_status": meta.get("http_status"),
        "retrieved_at": meta.get("retrieved_at", UNKNOWN),
        "more_available": more_available,
        "note": "Catalog metadata only: Data.gov does not host the records; retrieve them from retrieve_from.",
    }


def unknown_result(query: str, status: str, note: str, meta: dict | None = None) -> dict:
    """Explicit 'nothing found / unreachable' result. Never a fabricated dataset."""
    meta = meta or {}
    out = {f: UNKNOWN for f in RESULT_FIELDS}
    out.update({
        "is_unknown": True,
        "result_status": status,
        "query": query,
        "rank": None,
        "resources": [],
        "resource_urls": [],
        "resource_formats": [],
        "keywords": [],
        "api_endpoint": meta.get("api_endpoint", UNKNOWN),
        "request_url": meta.get("request_url", UNKNOWN),
        "http_status": meta.get("http_status"),
        "retrieved_at": meta.get("retrieved_at", utc_now_iso()),
        "more_available": False,
        "note": note,
    })
    return out


def _search(query: str, rows: int = DEFAULT_ROWS, *, org_type: str | None = None, endpoint: str = "auto",
            refresh: bool = False) -> tuple[list[dict], dict | None, dict]:
    """(normalised results, raw payload or None, call metadata)."""
    payload, meta = search_raw(query, rows, org_type=org_type, endpoint=endpoint, refresh=refresh)
    if payload is None:
        return [unknown_result(query, "unreachable",
                               "Data.gov catalog API unreachable or errored: " + " | ".join(meta["errors"]),
                               meta)], None, meta
    results = payload.get("results") or []
    if not results:
        return [unknown_result(query, "no_results",
                               "Data.gov catalog returned 0 datasets for this query (AND over terms).",
                               meta)], payload, meta
    more = bool(payload.get("after"))
    return [normalise_record(r, rank=i + 1, query=query, meta=meta, more_available=more)
            for i, r in enumerate(results)], payload, meta


def search_us_open_data(query: str, rows: int = DEFAULT_ROWS, *, org_type: str | None = None,
                        endpoint: str = "auto", refresh: bool = False) -> list[dict]:
    """Search the Data.gov catalog and return routing records (metadata only).

    Each dict carries: title, agency/publisher, parent_organization, description,
    access_level, resources [{url, format, name, media_type, url_kind}], modified,
    identifier, landing_page, and retrieve_from (the agency URL that holds the data).
    When the API returns nothing or cannot be reached, the list holds exactly one
    record with is_unknown=True and every descriptive field = UNKNOWN.
    """
    results, _, _ = _search(query, rows, org_type=org_type, endpoint=endpoint, refresh=refresh)
    return results


# ---------------------------------------------------------------------------
# Source matching
# ---------------------------------------------------------------------------
_STOP = {"the", "and", "of", "for", "data", "dataset", "a", "an", "to", "in", "on", "by", "with"}


def _rule_for(source: dict) -> dict:
    sid = source["source_id"]
    if sid in DISCOVERY_RULES:
        return {**DISCOVERY_RULES[sid], "curated": True}
    name = re.sub(r"\(.*?\)", "", str(source.get("name", sid))).strip()
    words = [w for w in re.findall(r"[A-Za-z0-9.]+", name) if w.lower() not in _STOP][:3]
    title = "".join(f"(?=.*{re.escape(w)})" for w in words) or re.escape(name)
    return dict(queries=[name], title=title, publisher=None, curated=False)


def _who(result: dict) -> str:
    return " | ".join(str(result.get(k, "")) for k in ("publisher", "publisher_hierarchy", "parent_organization"))


def publisher_matches(result: dict, rule: dict) -> bool:
    if result.get("is_unknown") or not rule.get("publisher"):
        return False
    return bool(re.search(rule["publisher"], _who(result), flags=re.I))


def match_scope(result: dict, rule: dict) -> str | None:
    """'planned_dataset', 'related_dataset' or None for one result under one source rule."""
    if result.get("is_unknown"):
        return None
    title = str(result.get("title", ""))
    if re.search(NON_DATASET_TITLE, title, flags=re.I):
        return None
    if not re.search(rule["title"], title, flags=re.I):
        return None
    if rule.get("publisher") and not publisher_matches(result, rule):
        return None
    planned_ok = re.search(rule.get("planned") or rule["title"], title, flags=re.I)
    excluded = rule.get("exclude_planned") and re.search(rule["exclude_planned"], title, flags=re.I)
    return "planned_dataset" if planned_ok and not excluded else "related_dataset"


def matches_rule(result: dict, rule: dict) -> bool:
    """True when the result is the planned dataset of this source (counts as indexed)."""
    return match_scope(result, rule) == "planned_dataset"


def first_match(results: list[dict], rule: dict, scope: str = "planned_dataset") -> dict | None:
    return next((r for r in results if match_scope(r, rule) == scope), None)


def best_match(results: list[dict], rule: dict) -> tuple[dict | None, str | None]:
    """Best-ranked planned-dataset match, else best-ranked related match."""
    for scope in ("planned_dataset", "related_dataset"):
        m = first_match(results, rule, scope)
        if m is not None:
            return m, scope
    return None, None


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")[:60]


def _save_raw_response(kind: str, query: str, rows: int, payload: dict | None, meta: dict) -> str:
    """Store the query's raw payload + call metadata under data/raw (overwritten each run)."""
    d = raw_dir(SOURCE_ID) / "search_responses"
    d.mkdir(parents=True, exist_ok=True)
    body = json.dumps({"query": query, "rows": rows, "meta": meta, "payload": payload},
                      indent=1, sort_keys=True, default=str)
    name = f"{kind}__{_slug(query)}__rows{rows}.json"
    (d / name).write_text(body)
    return f"search_responses/{name}"


def _log_row(kind: str, query: str, attempt: int, results: list[dict], *, source_id=None,
             rule=None, raw_file="") -> dict:
    """One discovery-log row. matched_source_id = the planned source this row is about
    (for registry queries: the source being checked; for demo queries: the source whose
    rule the best-ranked result satisfies, else UNKNOWN). indexed = planned dataset found."""
    top = results[0]
    n = 0 if top["is_unknown"] else len(results)
    match, scope = best_match(results, rule) if rule else (None, None)
    pub_hit = any(publisher_matches(r, rule) for r in results) if rule else False
    if top["result_status"] == "unreachable":
        status = "query_failed"
    elif top["result_status"] == "no_results":
        status = "no_results"
    elif scope == "planned_dataset":
        status = "dataset_indexed"
    elif scope == "related_dataset":
        status = "related_dataset_only"
    elif pub_hit:
        status = "publisher_indexed_dataset_not_matched"
    else:
        status = "not_found"
    fmt_urls = lambda r: " | ".join(r["resource_urls"][:5]) if r and r["resource_urls"] else UNKNOWN  # noqa: E731
    return {
        "query_kind": kind,
        "query": query,
        "attempt": attempt,
        "n_results": n,
        "more_available": bool(top.get("more_available")),
        "result_status": top["result_status"],
        "top_title": top["title"],
        "top_publisher": top["publisher"],
        "top_parent_organization": top["parent_organization"],
        "top_resource_urls": fmt_urls(top),
        "top_retrieve_from": top["retrieve_from"],
        "matched_source_id": source_id if source_id else UNKNOWN,
        "match_rule_curated": bool(rule.get("curated", True)) if rule else False,
        "match_scope": scope or UNKNOWN,
        "match_rank": match["rank"] if match else None,
        "match_title": match["title"] if match else UNKNOWN,
        "match_publisher": match["publisher"] if match else UNKNOWN,
        "match_access_level": match["access_level"] if match else UNKNOWN,
        "match_modified": match["modified"] if match else UNKNOWN,
        "match_identifier": match["identifier"] if match else UNKNOWN,
        "match_landing_page": match["landing_page"] if match else UNKNOWN,
        "match_resource_urls": fmt_urls(match),
        "match_retrieve_from": match["retrieve_from"] if match else UNKNOWN,
        "match_retrieve_from_host": match["retrieve_from_host"] if match else UNKNOWN,
        "match_machine_readable_url": match["machine_readable_url"] if match else UNKNOWN,
        "indexed": "yes" if scope == "planned_dataset" else "no",
        "index_status": status,
        "api_endpoint": top["api_endpoint"],
        "request_url": top["request_url"],
        "raw_response_file": raw_file,
        "retrieved_at": top["retrieved_at"],
        "note": top["note"] if top["is_unknown"] else "",
    }


def _results_rows(kind: str, query: str, attempt: int, results: list[dict], source_id: str | None) -> list[dict]:
    rows = []
    for r in results:
        row = {k: r[k] for k in RESULT_FIELDS if k not in ("resources", "resource_urls", "resource_formats", "keywords")}
        row.update({
            "query_kind": kind,
            "attempt": attempt,
            "for_source_id": source_id or UNKNOWN,
            "resources_json": json.dumps(r["resources"]),
            "n_resources": len(r["resources"]),
            "resource_urls": " | ".join(r["resource_urls"]) if r["resource_urls"] else UNKNOWN,
            "resource_formats": "; ".join(r["resource_formats"]) if r["resource_formats"] else UNKNOWN,
            "keywords": "; ".join(r["keywords"]) if r["keywords"] else UNKNOWN,
            "record_key": f"{kind}|{query}|{r['rank'] if r['rank'] is not None else 'none'}",
        })
        rows.append(row)
    return rows


def _result_ids(content: bytes) -> list[str] | None:
    try:
        return [str((r.get("dcat") or {}).get("identifier") or r.get("identifier"))
                for r in json.loads(content).get("results", [])]
    except Exception:
        return None


def probe_endpoints(refresh: bool = False) -> dict:
    """Record the status of the legacy and current endpoints (cached, incl. error codes).

    measure_it.http keys its cache on method + URL + params (not headers), so probes that
    differ only by header use different params to get their own cache entries; the no-key
    probe adds ``sort``, which search_raw never sends, so its cached 403 cannot answer a search.
    """
    out = {"probed_at": utc_now_iso(), "probes": []}

    def probe(label, url, params=None, headers=None):
        try:
            interval = ORIGIN_CRAWL_DELAY_S if urlparse(url).netloc == "catalog.data.gov" else GATEWAY_MIN_INTERVAL_S
            r = _polite_get(url, params=params, headers=headers, min_interval=interval, cache_errors=True,
                            reject_html=False, max_retries=2, timeout=60, refresh=refresh)
            hdr = {k.lower(): v for k, v in r.headers.items()}
            entry = dict(label=label, url=url, params=params, sent_headers=sorted((headers or {}).keys()),
                         final_url=r.url, http_status=r.status,
                         content_type=hdr.get("content-type"), bytes=len(r.content),
                         sha256=hashlib.sha256(r.content).hexdigest(), body_head=r.text[:200],
                         x_ratelimit_limit=hdr.get("x-ratelimit-limit"), fetched_at=r.fetched_at,
                         result_ids=_result_ids(r.content) if r.status == 200 else None)
        except Exception as exc:
            entry = dict(label=label, url=url, params=params, error=f"{type(exc).__name__}: {exc}"[:300])
        out["probes"].append(entry)
        return entry

    probe("legacy CKAN package_search (catalog.data.gov)", LEGACY_CKAN_URL, {"q": "lyme", "rows": 1})
    probe("legacy CKAN via api.gsa.gov v3 gateway (DEMO_KEY)", LEGACY_GATEWAY_V3_URL, {"q": "lyme", "rows": 1},
          {"X-Api-Key": DEMO_KEY})
    o = probe("current search, origin (keyless)", ORIGIN_SEARCH_URL, {"q": "lyme", "per_page": 5})
    g = probe("current search, api.gsa.gov v4 gateway (DEMO_KEY)", GATEWAY_SEARCH_URL, {"q": "lyme", "per_page": 5},
              {"X-Api-Key": DEMO_KEY})
    probe("current search, api.gsa.gov v4 gateway (no key)", GATEWAY_SEARCH_URL,
          {"q": "lyme", "per_page": 5, "sort": "relevance"})
    probe("bulk metadata dump listed in the catalog's own 'Data.gov CKAN API' record", BULK_DUMP_URL)
    out["origin_and_gateway_same_results"] = (
        o.get("http_status") == 200 and g.get("http_status") == 200
        and bool(o.get("result_ids")) and o.get("result_ids") == g.get("result_ids"))
    out["origin_and_gateway_identical_bytes"] = bool(o.get("sha256")) and o.get("sha256") == g.get("sha256")
    return out


def _download_catalog_files() -> dict:
    files = {}
    for fname, url in [("openapi.json", OPENAPI_URL), ("organizations.json", ORGANIZATIONS_URL),
                       ("catalog_stats.json", STATS_URL), ("robots.txt", ROBOTS_URL)]:
        try:
            files[fname] = download_file(url, SOURCE_ID, fname)
        except Exception as exc:  # recorded, not fatal
            print(f"[datagov] could not download {url}: {exc}")
    return files


def run_discovery(rows: int = DEFAULT_ROWS) -> tuple[pd.DataFrame, pd.DataFrame]:
    reg = load_source_registry()
    log_rows, result_rows = [], []
    for src in reg.get("sources", []):
        sid = src["source_id"]
        rule = _rule_for(src)
        for attempt, q in enumerate(rule["queries"], start=1):
            res, payload, meta = _search(q, rows)
            raw_file = _save_raw_response("source", q, rows, payload, meta)
            log_rows.append(_log_row("registry_source", q, attempt, res, source_id=sid, rule=rule, raw_file=raw_file))
            result_rows += _results_rows("registry_source", q, attempt, res, sid)
            print(f"[datagov] {sid:32s} a{attempt} {q!r:55s} n={log_rows[-1]['n_results']:>3} "
                  f"-> {log_rows[-1]['index_status']}")
            if log_rows[-1]["indexed"] == "yes":
                break
    rules = {s["source_id"]: _rule_for(s) for s in reg.get("sources", [])}
    for q in DEMO_QUERIES:
        res, payload, meta = _search(q, rows)
        raw_file = _save_raw_response("demo", q, rows, payload, meta)
        # which planned source (if any) does the best-ranked matching result correspond to?
        # planned-dataset matches take precedence over related ones, then rank.
        hit_sid, hit_key = None, None
        for sid, rule in rules.items():
            m, scope = best_match(res, rule)
            if m is not None:
                key = (0 if scope == "planned_dataset" else 1, m["rank"])
                if hit_key is None or key < hit_key:
                    hit_sid, hit_key = sid, key
        row = _log_row("demo", q, 1, res, source_id=hit_sid, rule=rules.get(hit_sid) if hit_sid else None,
                       raw_file=raw_file)
        if hit_sid is None and row["result_status"] == "ok":
            row["index_status"] = "no_planned_source_matched"
        log_rows.append(row)
        result_rows += _results_rows("demo", q, 1, res, hit_sid)
        print(f"[datagov] demo {q!r:55s} n={row['n_results']:>3} matched={row['matched_source_id']}")
    log = pd.DataFrame(log_rows)
    # Source-level verdict over all attempts: best status wins (ties -> earliest attempt).
    src_mask = log["query_kind"] == "registry_source"
    rank = log["index_status"].map(STATUS_ORDER).fillna(len(STATUS_ORDER))
    log["_order"] = list(zip(rank, log["attempt"]))
    best_idx = log[src_mask].groupby("matched_source_id")["_order"].idxmin()
    log["is_best_attempt"] = False
    log.loc[best_idx.values, "is_best_attempt"] = True
    log.loc[~src_mask, "is_best_attempt"] = True
    best_status = log.loc[best_idx.values].set_index("matched_source_id")["index_status"]
    log["source_index_status"] = log["matched_source_id"].map(best_status).where(src_mask, log["index_status"])
    log["source_indexed"] = (log["source_index_status"] == "dataset_indexed").map({True: "yes", False: "no"})
    last_attempt = log[src_mask].groupby("matched_source_id")["attempt"].transform("max")
    log["is_final_attempt"] = True
    log.loc[src_mask, "is_final_attempt"] = (log.loc[src_mask, "attempt"] == last_attempt).astype(bool)
    log["is_final_attempt"] = log["is_final_attempt"].astype(bool)
    log = log.drop(columns="_order")
    return log, pd.DataFrame(result_rows)


def _prune_stale_raw_responses(keep: set[str]) -> list[str]:
    """Remove search_responses/*.json left by queries that are no longer run, so the raw directory
    holds exactly the responses behind the processed tables (the HTTP cache keeps everything)."""
    d = raw_dir(SOURCE_ID) / "search_responses"
    removed = []
    for p in sorted(d.glob("*.json")):
        rel = f"search_responses/{p.name}"
        if rel not in keep:
            p.unlink()
            removed.append(rel)
    if removed:
        print(f"[datagov] removed {len(removed)} stale raw response file(s): {removed}")
    return removed


def _version_string() -> str:
    p = raw_dir(SOURCE_ID) / "openapi.json"
    try:
        info = json.loads(p.read_text()).get("info", {})
        return f"{info.get('title', 'Datagov Catalog')} API {info.get('version', '?')} (openapi.json); DCAT-US metadata"
    except Exception:
        return "Datagov Catalog API (version unknown); DCAT-US metadata"


def build(rows: int = DEFAULT_ROWS, refresh_probes: bool = False) -> dict:
    rdir = raw_dir(SOURCE_ID)
    _download_catalog_files()
    probes = probe_endpoints(refresh=refresh_probes)
    (rdir / "endpoint_probes.json").write_text(json.dumps(probes, indent=2, default=str))

    log, results = run_discovery(rows)
    _prune_stale_raw_responses(set(log["raw_response_file"]))
    version = _version_string()
    notes = ("Data.gov catalog metadata (DCAT-US) returned by the catalog search API; Data.gov does not host "
             "the records. data_layer = metadata: rows describe datasets, not facilities or places.")
    retrieved = log["retrieved_at"].min()
    log_p = add_provenance(
        log, data_layer="metadata", source_name="Data.gov catalog search API (GSA)",
        source_version=version, retrieved_at=retrieved, evidence_type="metadata_catalog",
        source_record_id=lambda d: d["query_kind"] + "|" + d["query"],
        source_geographic_resolution="none", evidence_level="catalog_metadata", provenance_notes=notes)
    log_p["retrieved_at"] = log["retrieved_at"].values  # per-query fetch time from the HTTP cache
    res_p = add_provenance(
        results, data_layer="metadata", source_name="Data.gov catalog search API (GSA)",
        source_version=version, retrieved_at=retrieved, evidence_type="metadata_catalog",
        source_record_id=lambda d: d["identifier"].where(~d["is_unknown"].astype(bool), d["record_key"]),
        source_geographic_resolution="none", evidence_level="catalog_metadata", provenance_notes=notes)
    res_p["retrieved_at"] = results["retrieved_at"].values
    write_table(log_p, "data_gov_discovery_log", producer=PRODUCER,
                description="One row per Data.gov catalog query: registry sources (curated title/publisher rules) "
                            "and demo queries; indexed yes/no and the agency URL each match routes to.")
    write_table(res_p, "data_gov_search_results", producer=PRODUCER,
                description="One row per (query, result rank) of Data.gov catalog metadata; UNKNOWN rows for "
                            "queries with no results. retrieve_from = agency URL holding the data.")
    stats = _write_registry_and_audit(log_p, res_p, probes, version)
    return stats


# ---------------------------------------------------------------------------
# Registry entry + DATA_AUDIT.md (numbers computed from the outputs)
# ---------------------------------------------------------------------------
def _pct(n: int, d: int) -> str:
    return f"{n} ({100.0 * n / d:.1f}%)" if d else "0 (n/a)"


def _write_registry_and_audit(log: pd.DataFrame, res: pd.DataFrame, probes: dict, version: str) -> dict:
    rdir = raw_dir(SOURCE_ID)
    stats_path = rdir / "catalog_stats.json"
    catalog_total = UNKNOWN
    try:
        catalog_total = json.loads(stats_path.read_text())["results"]["datasets"]
    except Exception:
        pass
    src = log[log["query_kind"] == "registry_source"]
    final = src[src["is_best_attempt"]]
    n_sources = final["matched_source_id"].nunique()
    n_indexed = int((final["source_indexed"] == "yes").sum())
    status_counts = Counter(final["source_index_status"])
    real = res[~res["is_unknown"].astype(bool)]
    n_unknown = int(res["is_unknown"].astype(bool).sum())
    n_unique = real["identifier"].nunique()
    endpoints_used = sorted(log["api_endpoint"].unique())
    ratelimits = [p.get("x_ratelimit_limit") for p in probes["probes"] if p.get("x_ratelimit_limit")]
    manifest = load_manifest(SOURCE_ID)
    max_harvest = real["last_harvested_date"][real["last_harvested_date"] != UNKNOWN].max() if len(real) else UNKNOWN
    access_counts = Counter(real["access_level"])
    license_known = int((real["license"] != UNKNOWN).sum())
    gsa_rows = real["publisher_hierarchy"].str.contains(GSA_PUBLISHER, case=False, regex=True).fillna(False).astype(bool)
    # astype(bool): on an empty frame (every query failed) map() keeps the str dtype and '&' would raise
    retrieve_on_catalog = int((real["retrieve_from"].map(lambda u: is_catalog_url(u) if _is_http(u) else False)
                               .astype(bool) & ~gsa_rows).sum())
    # a Data.gov outage or rate limit is recorded (CONVENTIONS 3.4), never a crash: blocked when every query failed,
    # partial when some did
    n_failed = int((log["index_status"] == "query_failed").sum())
    status = "blocked" if len(log) and n_failed == len(log) else ("partial" if n_failed else "ingested")
    failed_queries = log.loc[log["index_status"] == "query_failed", "query"].astype(str).tolist()
    miss = {}
    for col in ["title", "publisher", "parent_organization", "description", "access_level", "modified",
                "identifier", "landing_page", "license", "retrieve_from", "machine_readable_url"]:
        miss[col] = int((real[col].astype(str) == UNKNOWN).sum())
    no_resources = int((real["n_resources"] == 0).sum())
    basis_counts = Counter(real["retrieve_from_basis"])
    host_counts = Counter(real["retrieve_from_host"]).most_common(12)
    org_type_counts = Counter(real["organization_type"])

    entry = {
        "source_id": SOURCE_ID,
        "name": "Data.gov catalog (metadata only)",
        "publisher": "U.S. General Services Administration (GSA), Data.gov",
        "landing_url": "https://catalog.data.gov/",
        "access_urls": [ORIGIN_SEARCH_URL, GATEWAY_SEARCH_URL, OPENAPI_URL, DOCS_URL],
        "license": ("Not stated by the API. Catalog metadata is published by a U.S. federal agency (GSA); each listed "
                    "dataset's own license (dcat.license, present on "
                    f"{license_known}/{len(real)} retrieved records) governs the records at the agency host."),
        "access_conditions": ("Keyless origin catalog.data.gov/search (robots.txt Crawl-Delay: 10, honoured). "
                              "Gateway api.gsa.gov/technology/datagov/v4 needs an api.data.gov key via X-Api-Key; "
                              "DEMO_KEY documented at 30 req/IP/hour and 50/IP/day (resources.data.gov/catalog-api); "
                              f"gateway reported X-RateLimit-Limit {ratelimits[0] if ratelimits else UNKNOWN} on DEMO_KEY. "
                              "Legacy CKAN package_search returns HTTP 404."),
        "retrieved_at": str(log["retrieved_at"].min()),
        "source_version": version,
        "update_date": f"continuous harvest; latest last_harvested_date among retrieved records: {max_harvest}",
        "data_layer": "metadata",
        "unit_of_observation": "catalog dataset record (DCAT-US metadata) per query and rank; not the data itself",
        "sample_size": {
            "queries_executed": int(len(log)),
            "registry_sources_checked": int(n_sources),
            "registry_sources_indexed": n_indexed,
            "registry_source_status_counts": dict(status_counts),
            "demo_queries": len(DEMO_QUERIES),
            "result_rows": int(len(res)),
            "result_rows_real": int(len(real)),
            "result_rows_unknown": n_unknown,
            "unique_dataset_identifiers": int(n_unique),
            "queries_failed": n_failed,
            "catalog_total_datasets_api_stats": catalog_total,
        },
        "geographic_resolution": "none (catalog metadata; spatial fields not used)",
        "person_level": False,
        "geographic": False,
        "omics": False,
        "wearable": False,
        "participant_linkage": "not applicable (metadata catalog; no participants)",
        "true_participant_linkage_across_modalities": False,
        "status": status,
        "processed_outputs": ["data_gov_discovery_log", "data_gov_search_results"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.agents.datagov",
        "limitations": [
            "Metadata only: Data.gov does not host records; retrieve_from routes to the agency host.",
            "Search is AND over terms with no total-hit count; n_results = first page, capped at rows.",
            "indexed yes/no = the planned dataset (curated title + planned + publisher regex) appears in the top "
            "rows results of a curated query; same-programme records (derived indicators, locator tools, portals) "
            "are reported as related_dataset_only, not as indexed. A source can exist in the catalog under wording "
            "the rules/queries miss (false 'no').",
            "Non-federal resources (Stanford, Monarch, EMBL-EBI, Open Targets, NUCC) are not expected in the catalog.",
            "Catalog relevance ranking and harvested contents change over time; results are a dated snapshot.",
        ] + ([f"{n_failed} of {len(log)} catalog queries failed (API unreachable, errored or rate-limited); their "
              "rows are explicit UNKNOWN results, not zero-result answers."] if n_failed else []),
    }
    write_registry_entry(entry)

    # ------------------------- DATA_AUDIT.md -------------------------
    probe_lines = []
    for p in probes["probes"]:
        if "error" in p:
            err = ("DNS: host name does not resolve (NameResolutionError)" if "NameResolution" in p["error"]
                   else p["error"][:200].replace("|", "/"))
            probe_lines.append(f"| {p['label']} | `{p['url']}` | error | {err} |")
        else:
            probe_lines.append(
                f"| {p['label']} | `{p['url']}` params={json.dumps(p['params'])} | HTTP {p['http_status']} "
                f"(final `{p['final_url'].split('?')[0]}`) | {p['content_type']}, {p['bytes']} B; "
                f"X-RateLimit-Limit={p.get('x_ratelimit_limit') or '-'}; head: `{(p['body_head'] or '')[:90].replace('|', '/').replace(chr(10), ' ')}` |")
    man_lines = [f"| `{k}` | {v['url']} | {v['bytes']} | `{v['sha256'][:16]}…` | {v['retrieved_at']} |"
                 for k, v in sorted(manifest["files"].items())]
    tries = src.groupby("matched_source_id")["query"].apply(lambda q: len(q))
    src_lines = []
    for _, r in final.sort_values("matched_source_id").iterrows():
        mt = r["match_title"] if r["match_title"] != UNKNOWN else "—"
        rf = r["match_retrieve_from"] if r["match_retrieve_from"] != UNKNOWN else "—"
        mr = r["match_machine_readable_url"] if r["match_machine_readable_url"] != UNKNOWN else "—"
        src_lines.append(
            f"| {r['matched_source_id']} | {r['query']} (attempt {r['attempt']} of {tries[r['matched_source_id']]}) | "
            f"{r['n_results']}{'+' if r['more_available'] else ''} | **{r['source_indexed']}** | "
            f"{r['source_index_status']} | {mt[:80].replace('|', '/')} | {rf} | {mr} |")
    # records worth flagging, computed from the retrieved results
    notable = []
    ii = real[real["title"].str.contains("Invisible Illness", case=False, na=False)].drop_duplicates("identifier")
    for _, r in ii.iterrows():
        notable.append(f"* **{r['title']}** — {r['publisher']} ({r['parent_organization']}); landing "
                       f"{r['landing_page']}; {r['n_resources']} distribution URLs; modified {r['modified']}. "
                       f"Returned for query '{r['query']}'. It is a challenge description page, not data.")
    dg_self = log[(log["matched_source_id"] == SOURCE_ID) & (log["indexed"] == "yes")]
    for _, r in dg_self.head(1).iterrows():
        rec = real[(real["query"] == r["query"]) & (real["title"] == r["match_title"])].head(1)
        if len(rec):
            x = rec.iloc[0]
            notable.append(f"* The catalog's own record for its API, **{x['title']}** (modified {x['modified']}), "
                           f"still describes the retired CKAN API: \"{x['description'][:160]}…\". Its "
                           f"distributions: {x['resource_urls']}.")
    dump = next((q for q in probes["probes"] if q["url"] == BULK_DUMP_URL), None)
    if dump:
        err = dump.get("error") or ""
        outcome = ("host name does not resolve in DNS (NameResolutionError)" if "NameResolution" in err
                   else err[:160] or f"HTTP {dump.get('http_status')}")
        notable.append(f"* Bulk metadata dump `{BULK_DUMP_URL}` listed in that record: {outcome}.")
    nonpub = real[real["access_level"].isin(["non-public", "restricted public"])].drop_duplicates("identifier")
    notable.append(f"* {len(nonpub)} unique retrieved records are flagged non-public/restricted public; "
                   f"{int((nonpub['n_resources'] > 0).sum())} of them still list distribution URLs "
                   f"(e.g. {'; '.join(nonpub['title'].head(3))}).")
    rerouted = real[real["retrieve_from_basis"].str.contains("healthdata.gov catalog page", regex=False)]
    on_hd = real[real["retrieve_from_host"].isin(SECONDARY_CATALOG_HOSTS)]
    notable.append(f"* {len(rerouted)} of {len(real)} result rows have a healthdata.gov (HHS catalog) landing page but a "
                   f"distribution on the agency host; `retrieve_from` uses that agency URL (e.g. "
                   f"{'; '.join((rerouted['title'].str[:45] + ' -> ' + rerouted['retrieve_from_host']).drop_duplicates().head(3))}). "
                   f"{len(on_hd)} rows still route to healthdata.gov: data hosted there by HHS, or records with no other URL.")
    demo_lines = []
    for _, r in log[log["query_kind"] == "demo"].iterrows():
        demo_lines.append(
            f"| {r['query']} | {r['n_results']}{'+' if r['more_available'] else ''} | {r['top_title'][:60]} | "
            f"{r['top_publisher']} | {r['matched_source_id']} | {r['index_status']} |")
    mr_web = real[real["machine_readable_url"].str.contains(r"(?:\.aspx?|\.cfm|\.html?|/)$", regex=True)]
    n_mr_webpage = len(mr_web)
    mr_webpage_eg = "; ".join((mr_web["title"].str[:40] + " -> " + mr_web["machine_readable_url"])
                              .drop_duplicates().head(2)) or "none"
    miss_lines = [f"| {k} | {_pct(v, len(real))} |" for k, v in miss.items()]
    miss_lines.append(f"| resources (0 distribution URLs) | {_pct(no_resources, len(real))} |")
    audit = f"""# DATA AUDIT — Data.gov catalog (metadata only)

Generated by `measure_it.agents.datagov` at {utc_now_iso()}. Every number below is computed from the
retrieved API responses in this directory.

| Field | Value |
|---|---|
| source_id | {SOURCE_ID} |
| Source (dataset/API name, exact files/endpoints) | Data.gov catalog search API: `GET {ORIGIN_SEARCH_URL}` (origin) and `GET {GATEWAY_SEARCH_URL}` (api.data.gov gateway); spec `{OPENAPI_URL}` |
| Publishing organization | U.S. General Services Administration (GSA), Data.gov program |
| Retrieval date (UTC) | {log['retrieved_at'].min()} to {log['retrieved_at'].max()} |
| Source version / release | {version} |
| Source update date / cadence | Continuous harvest from agency catalogs; latest `last_harvested_date` among retrieved records: {max_harvest} |
| License / access conditions | {entry['license']} {entry['access_conditions']} |
| Unit of observation | {entry['unit_of_observation']} |
| Sample size (actual, as ingested) | {len(log)} queries executed; {len(res)} result rows ({len(real)} dataset records, {n_unknown} explicit UNKNOWN rows); {n_unique} unique dataset identifiers. Catalog size per `/api/stats`: {catalog_total} datasets |
| Geography (resolution, vintage) | none (catalog metadata) |
| Person-level? | no |
| Geographic? | no |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — metadata catalog, no participants |
| Status | {status}{f": {n_failed} of {len(log)} queries failed (API unreachable, errored or rate-limited): {'; '.join(failed_queries[:10])}{' ...' if n_failed > 10 else ''}" if n_failed else " (every query answered)"} |

## Files / endpoints retrieved

### Endpoint status (measured; `endpoint_probes.json`)
| Probe | URL | Status | Detail |
|---|---|---|---|
{chr(10).join(probe_lines)}

Origin and gateway returned the same ordered result identifiers for the same query:
**{probes.get('origin_and_gateway_same_results')}** (byte-identical: {probes.get('origin_and_gateway_identical_bytes')};
bytes can differ only through OpenSearch `_score` jitter between replicas).

The CKAN action API (`/api/3/action/package_search`) that most client code and older docs
(including open.gsa.gov/api/datadotgov, which still shows the v3 CKAN example) describe is gone.
resources.data.gov/catalog-api states the prior endpoint "remains available in a read-only state";
that is not what we observe: it answers HTTP 404.

Endpoints used for the searches: {', '.join(endpoints_used)}.
API key: none configured (checked env vars {', '.join(API_KEY_ENV_VARS)}); no account was created.

### Bulk files (MANIFEST.json)
| File | URL | Bytes | sha256 | Retrieved |
|---|---|---|---|---|
{chr(10).join(man_lines)}

### Search responses
Raw JSON for every query is in `search_responses/` (query, parameters, endpoint, fetch time, full payload).
Query parameters: `q=<query>`, `per_page={DEFAULT_ROWS}`; the most results returned by any query was
{int(log['n_results'].max()) if len(log) else 0}. Responses are also cached under `data/_http_cache`.

## Key variables
| Output field | DCAT-US / API source | Meaning | Downstream use |
|---|---|---|---|
| title | `title` | dataset title | display, matching |
| agency / publisher | `publisher` (top level) / `dcat.publisher.name` | publishing agency or office | routing, matching |
| publisher_hierarchy | `dcat.publisher.subOrganizationOf` chain + `organization.name` | agency chain | routing |
| parent_organization | `organization.name` | Data.gov organisation that owns the harvest source (e.g. HHS) | grouping |
| description | `description` | abstract | display |
| access_level | `dcat.accessLevel` | public / restricted public / non-public | flag gated data |
| resources | `dcat.distribution[].downloadURL|accessURL, format|mediaType, title` | agency-hosted files/APIs | retrieval |
| modified / issued | `dcat.modified` / `dcat.issued` | agency-declared dates (strings, not normalised) | recency |
| identifier | `dcat.identifier` | agency identifier (often a URL) | dedup, provenance |
| landing_page | `dcat.landingPage` | agency landing page | retrieval |
| retrieve_from | landingPage (object-valued landingPage -> its accessURL), else first distribution URL, else identifier — never a catalog.data.gov/api.gsa.gov URL; a healthdata.gov landing page yields to a distribution on another agency host | where to get the data | routing |
| machine_readable_url | first distribution whose declared media type is data (not text/html, PDF/Word/PowerPoint documents, images, or missing) | direct API/file | retrieval |
| datagov_url | `https://catalog.data.gov/dataset/<slug>` | catalog page (metadata only) | citation |

## Missingness
Measured over the {len(real)} dataset records returned (all queries; records repeat across queries).

| Field | UNKNOWN / missing |
|---|---|
{chr(10).join(miss_lines)}

access_level distribution: {dict(access_counts)}.
Organisation types: {dict(org_type_counts)}.
retrieve_from basis: {dict(basis_counts)}. retrieve_from pointing at a Data.gov catalog host (records not
published by GSA itself): {retrieve_on_catalog}.
Most common retrieve_from hosts: {host_counts}.

## Registry sources: is each planned source indexed?
**{n_indexed} of {n_sources}** planned sources in SOURCE_REGISTRY.yaml have their planned dataset in the top
{DEFAULT_ROWS} results of at least one curated query. Per-source status counts: {dict(status_counts)}.

Statuses: `dataset_indexed` = a record of the planned dataset by the expected agency (title rule AND
`planned` rule AND publisher rule); `related_dataset_only` = only a same-programme record by the expected
agency (a derived indicator, a locator tool, an agency portal, or a different file of the same series);
`publisher_indexed_dataset_not_matched` = the agency appears in the results but not the programme;
`not_found` = neither. Blog/news records (title "Blog | …") never count. Queries are tried in order until
the planned dataset is found; the best attempt per source is shown (all attempts are rows in
`data_gov_discovery_log`). `n` is the first-page result count (`+` = the API reported more).

| source_id | query (best attempt) | n | indexed | status | matched title | routes to (retrieve_from) | machine-readable URL |
|---|---|---|---|---|---|---|---|
{chr(10).join(src_lines)}

## Notable catalog records
{chr(10).join(notable)}

## Demo queries
| query | n | top title | top publisher | matched planned source | status |
|---|---|---|---|---|---|
{chr(10).join(demo_lines)}

## Linkage strategy
None to people or places. `data_gov_discovery_log.matched_source_id` links a query to a planned
`source_id` in SOURCE_REGISTRY.yaml by rule (title regex AND publisher regex), not by identifier.
`data_gov_search_results.identifier` is the agency's own dataset identifier. Nothing here is joinable to
person, geographic or facility rows; it only says where an agency publishes a dataset.

## Limitations and caveats
* Data.gov is a metadata catalog. Records are harvested from agency catalogs; the data stay at the agency
  host (`retrieve_from`). Some landing pages are generic division pages (e.g. CDC PLACES points at
  www.cdc.gov/nccdphp/dph/); `machine_readable_url` then carries the specific agency file/API.
* The search is AND over terms (a long exact title can return 0 while its words OR'ed return many), and the
  response has no total count. `n_results` is first-page count, capped at rows={DEFAULT_ROWS}.
* "indexed = no" means our curated queries and rules found no record of the planned dataset in the top
  {DEFAULT_ROWS}; it is not proof of absence. The rules are regexes written after reading the returned titles
  (see `DISCOVERY_RULES` in the module); they encode a judgement of what counts as "the planned dataset".
* Non-federal resources (Stanford, Monarch/Mondo, EMBL-EBI OLS, Open Targets, NHGRI-EBI GWAS Catalog, NUCC)
  are not expected in a U.S. government catalog; their "no" is informative only as a negative control.
* Duplicate titles occur (the same dataset under several agency identifiers); rows are kept as returned.
* `modified`/`issued` are agency-declared strings; formats vary (dates, datetimes, ISO-8601 intervals).
* `access_level = non-public` records can still list distribution URLs; the flag is reported as-is.
* `machine_readable_url` trusts the agency-declared `mediaType`, which is sometimes wrong: {n_mr_webpage} rows'
  machine_readable_url is a web page (ends in .aspx/.cfm/.htm(l) or '/') declared as CSV/Excel
  (e.g. {mr_webpage_eg}).
* data_layer is `metadata` (catalog / inventory metadata about sources; it was `facility` before the `metadata`
  layer was added on 2026-09-23); the rows are dataset metadata, not facilities.

## Processed outputs
| Table | Rows |
|---|---|
| data_gov_discovery_log | {len(log)} |
| data_gov_search_results | {len(res)} |

## Reproduce
`uv run python -m measure_it.agents.datagov`  (set DATA_GOV_API_KEY to use the api.data.gov gateway first)
"""
    (rdir / "DATA_AUDIT.md").write_text(audit)
    return {"log_rows": len(log), "result_rows": len(res), "sources": n_sources, "indexed": n_indexed,
            "unknown_rows": n_unknown, "queries_failed": n_failed, "status": status}


def run() -> dict:
    stats = build()
    print(f"[datagov] done: {stats}")
    return stats


if __name__ == "__main__":
    run()
