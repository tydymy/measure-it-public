"""CMS Office of Minority Health — Mapping Medicare Disparities (MMD) Tool, Population View.

Where the data come from
------------------------
The MMD Population View (https://data.cms.gov/tools/mapping-medicare-disparities-by-population,
app at https://data.cms.gov/mmd-population/map.html) is not published as a catalogued bulk
dataset: it does not appear in https://data.cms.gov/data.json. The app itself reads
aggregate slices from the data.cms.gov data API

    https://data.cms.gov/data-api/v1/mmd-tool/?_source=<file>&<dimension>=<value>...

where <file> is chosen from the app's own lookup table
(/mmd-population/assets/codebook_crosswalk.csv) and each dimension filter accepts
`a|b|IS NULL` alternatives. This module repeats exactly those requests (same file lookup,
same filters, including the app's one-time 2019 ED-visit `elig` exception). It saves each
response as a raw JSON file with a MANIFEST entry, so the retrieval is fully reproducible.

What is ingested (Medicare fee-for-service population, all sexes/races/dual/eligibility):
    years       2012-2023 (2023 files carry the `_p` = preliminary claims suffix)
    geography   county, state/territory, nation
    prevalence  CCW 'Fibromyalgia, Chronic Pain and Fatigue' (MMD code 51),
                CCW 'Migraine and Other Chronic Headache' (59)            -> burden
                Disability as reason for entitlement (24), COVID-19 (134)  -> context
    context     hospitalization and ED-visit rates (all-cause 10, 51, 59, 24, COVID-19),
                average total / risk-adjusted total / principal cost for 51 and 59
    adjustment  unsmoothed actual (fltr 1) and unsmoothed age-standardized (fltr 2);
                empirical-Bayes smoothed rates (fltr 3/4) are NOT ingested
    age strata  all ages, plus <65 (the disability-entitled group) for prevalence

Outputs (data/processed)
    geo_condition_burden__cms_mmd      geography x year x target condition x adjustment x age
                                       stratum; prevalence % with burden_evidence_level
    geo_context__cms_mmd               utilization, spending and population-composition
                                       measures (long format); never burden estimates
    condition_burden_coverage__cms_mmd one row per target condition: which CMS/CCW measure,
                                       if any, covers it, the exact CCW algorithm and codes,
                                       and burden_evidence_level A/B/C/D with the rationale

Guardrails
    * Medicare FFS beneficiaries (mostly >= 65 plus disability-entitled) are not
      representative of working-age invisible-illness populations; coded-condition
      prevalence depends on diagnosis and coding practice.
    * Evidence levels: B for fibromyalgia (composite CCW algorithm that contains M79.7) and
      migraine (G43 + G44 headache syndromes); C for ME/CFS (the fibromyalgia/pain/fatigue
      composite contains R53.82 chronic fatigue but not G93.3x); D for everything else.
    * State rows keep source_geographic_resolution='state'; nothing is downscaled.
    * Legacy Connecticut county FIPS are matched to the flagged legacy rows in
      `geographies` and reported, never re-assigned to 2024 planning regions.
    * Two 1:1 county renames that CMS still publishes under their pre-2015 codes are bridged
      (46113 Shannon -> 46102 Oglala Lakota SD; 02270 Wade Hampton -> 02158 Kusilvak AK; FIPS_RENAMES):
      same boundary, so the value is the 2024 county's; `fips_as_published` / source_record_id keep
      the CMS code. Boundary changes (51515 Bedford city VA) are not bridged.
    * A published 0 % can be a suppression artefact (numerator 1-2 is set to 0 by CMS);
      such rows carry value_zero_possibly_suppressed=True.

Reproduce: uv run python -m measure_it.ingestion.cms_mmd
"""
from __future__ import annotations

import csv
import html
import io
import json
import re
import shutil
import subprocess
import time
from pathlib import Path
from urllib.parse import quote, urlencode

import numpy as np
import pandas as pd

from ..config import load_config, raw_dir
from ..download import download_file, load_manifest
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

SOURCE_ID = "cms_mmd"
PRODUCER = "ingestion.cms_mmd"
SOURCE_NAME = "CMS OMH Mapping Medicare Disparities (MMD) Tool - Population View (Medicare FFS)"
LANDING_URL = "https://data.cms.gov/tools/mapping-medicare-disparities-by-population"
APP_URL = "https://data.cms.gov/mmd-population/map.html"
API_BASE = "https://data.cms.gov/data-api/v1/mmd-tool/"
API_PAGE_SIZE = 500000  # the app requests _size=500000; a response this long could be truncated

REFERENCE_FILES = {
    # MMD app assets that define codes and the file lookup
    "app/map.html": APP_URL,
    "app/codebook_crosswalk.csv": "https://data.cms.gov/mmd-population/assets/codebook_crosswalk.csv",
    "app/strings.json": "https://data.cms.gov/mmd-population/assets/strings.json",
    "app/menus.js": "https://data.cms.gov/mmd-population/js/menus.js",
    "app/datadownload.js": "https://data.cms.gov/mmd-population/js/datadownload.js",
    # Documentation
    "docs/Mapping-Technical-Documentation.pdf":
        "https://www.cms.gov/About-CMS/Agency-Information/OMH/Downloads/Mapping-Technical-Documentation.pdf",
    "docs/MappingPublicFAQs.pdf":
        "https://www.cms.gov/About-CMS/Agency-Information/OMH/Downloads/MappingPublicFAQs.pdf",
    # Chronic Conditions Data Warehouse (CCW) condition algorithms used by MMD
    "ccw/condition-categories-other.html": "https://www2.ccwdata.org/web/guest/condition-categories-other",
    "ccw/condition-categories-chronic.html": "https://www2.ccwdata.org/web/guest/condition-categories-chronic",
    "ccw/other-condition-algorithms.pdf":
        "https://www2.ccwdata.org/documents/10280/19139421/other-condition-algorithms.pdf",
    "ccw/comparing-chronic-pain.pdf":
        "https://www2.ccwdata.org/documents/10280/19139421/comparing-chronic-pain.pdf",
    # The 30 CCW chronic conditions (MMD codes 1-23, 144-156) are published only as PDFs; the
    # condition-categories-chronic page has no code table. This combined PDF is text-searched.
    "ccw/chr-chronic-condition-algorithms.pdf":
        "https://www2.ccwdata.org/documents/10280/19139421/chr-chronic-condition-algorithms.pdf",
    # data.cms.gov DCAT catalog: used to document that MMD has no catalogued bulk dataset
    "catalog/data.cms.gov_data.json": "https://data.cms.gov/data.json",
}
HTML_REFERENCES = {"app/map.html", "ccw/condition-categories-other.html", "ccw/condition-categories-chronic.html"}
# run() cannot proceed without these; the others are documentation copies
REQUIRED_REFERENCES = {"app/codebook_crosswalk.csv", "app/strings.json", "app/menus.js",
                       "ccw/condition-categories-other.html", "ccw/condition-categories-chronic.html",
                       "catalog/data.cms.gov_data.json"}

# ---- MMD code books (decoded from app/menus.js, app/strings.json, app/datadownload.js) ----
YEAR_CODES = {"2": 2012, "3": 2013, "4": 2014, "5": 2015, "6": 2016, "7": 2017, "8": 2018,
              "9": 2019, "20": 2020, "21": 2021, "22": 2022, "23": 2023}
GEOGRAPHY_CODES = {"c": "county", "s": "state", "n": "national"}
FLTR_CODES = {"1": "unsmoothed_actual", "2": "unsmoothed_age_standardized",
              "3": "smoothed_actual", "4": "smoothed_age_standardized"}
AGE_CODES = {"all": "all", "0": "<65", "1": "65-74", "2": "75-84", "3": "85+", "4": "65+"}
DENCAT = {"1": ("11-499", 11, 499), "2": ("500-999", 500, 999), "3": ("1,000-4,999", 1000, 4999),
          "4": ("5,000-9,999", 5000, 9999), "5": ("10,000+", 10000, np.nan)}
MEASURES = {
    # code: (measure name, unit)  -- units from the app's tooltip suffixes
    "v": ("prevalence", "percent of beneficiaries"),
    "h": ("hospitalization_rate", "discharges per 1,000 beneficiaries"),
    "e": ("ed_visit_rate", "ED visits per 1,000 beneficiaries"),
    "t": ("average_total_cost", "USD per beneficiary with the condition"),
    "a": ("average_total_cost_risk_adjusted", "USD per beneficiary with the condition"),
    "p": ("average_principal_cost", "USD per beneficiary with the condition"),
}
# MMD condition codes used here -> label as it appears in the MMD menu (verified at run time)
MMD_CONDITIONS = {
    "51": "Fibromyalgia, Chronic Pain and Fatigue",
    "59": "Migraine and Other Chronic Headache",
    "24": "Disability (reason for Medicare eligibility)",
    "134": "COVID-19",
    "10": "All-cause (all beneficiaries)",
}
MMD_CONDITION_10_LABELS = {"h": "All-Cause Hospitalizations", "e": "All Emergency Department Visits"}

# What to request: measure -> (conditions, age strata, adjustments)
QUERY_PLAN = {
    "v": (["51", "59", "24", "134"], [".", "0"], ["1", "2"]),
    "h": (["10", "51", "59", "24", "134"], ["."], ["1", "2"]),
    "e": (["10", "51", "59", "24"], ["."], ["1", "2"]),
    "t": (["51", "59"], ["."], ["1", "2"]),
    "a": (["51", "59"], ["."], ["1", "2"]),
    "p": (["51", "59"], ["."], ["1", "2"]),
}
BURDEN_SOURCE_CONDITIONS = {"51", "59"}
CONTEXT_CATEGORY = {"v": "population_composition", "h": "utilization", "e": "utilization",
                    "t": "spending", "a": "spending", "p": "spending"}

# ---- target condition mapping (SPEC §E burden_evidence_level) ----
# Target ICD-10-CM codes that were searched for in every CCW algorithm on the downloaded CCW
# pages. These are search keys (prefix match), not an ontology mapping; the ontology module
# owns the authoritative code sets.
TARGET_SEARCH_CODES = {
    "long_covid": ["U09.9"],
    "me_cfs": ["G93.3", "G93.31", "G93.32", "G93.39", "R53.82"],
    "pots": ["G90.A"],
    "dysautonomia": ["G90.0", "G90.8", "G90.9", "I95.1"],
    "fibromyalgia": ["M79.7"],
    "eds_hsd": ["Q79.6", "M35.7"],
    "mcas": ["D89.4"],
    "lyme_disease": ["A69.2"],
    "ptlds": ["A69.2"],
    "post_infectious_syndrome": ["G93.3", "U09.9"],
    "migraine": ["G43", "G44"],
    "ibs": ["K58"],
    "gastroparesis": ["K31.84"],
    "endometriosis": ["N80"],
}
# (target condition_id, MMD code, level)
BURDEN_MAPPINGS = [
    ("fibromyalgia", "51", "B"),
    ("me_cfs", "51", "C"),
    ("migraine", "59", "B"),
]
CCW_ALGORITHM_FOR_CODE = {"51": "Fibromyalgia, Chronic Pain and Fatigue",
                          "59": "Migraine and Other Chronic Headache"}

LIMITATION = ("Medicare FFS beneficiaries (mostly >= 65 plus disability-entitled) - not representative of "
              "working-age invisible-illness populations; coded-condition prevalence depends on "
              "diagnosis/coding practice.")
ECOLOGICAL = "Ecological/population-level aggregate; not person-level; never join to layer-1 participants."
SUPPRESSION_NOTE = ("CMS suppression: cells with <11 beneficiaries are omitted; percentages with a numerator of "
                    "1-2 are published as 0 (value_zero_possibly_suppressed flags every published 0).")
# 1:1 county renames with a new FIPS code and an unchanged boundary (U.S. Census Bureau, "Substantial Changes to
# Counties and County Equivalent Entities: 2010-present"). MMD still publishes both counties under the pre-2015 code
# in every claims year (2012-2023). Bridged here, at ingestion, because it is a code-level identity of the source's
# county (no areal allocation, nothing inferred): the row takes the 2024 code in geo_id/county_fips and keeps the code
# CMS published in fips_as_published and source_record_id. Boundary changes are NOT bridged: 51515 Bedford city VA
# (merged into Bedford County 51019 in 2013) stays unmatched, like the CMS xx990 codes. The CDC Lyme and HRSA
# ingestions apply the same two renames (cdc_lyme.FIPS_RENAMES_2000S, hrsa.FIPS_RENAMES_1TO1).
FIPS_RENAMES = {
    "46113": ("46102", "Shannon County, South Dakota, renamed Oglala Lakota County (46102) in 2015"),
    "02270": ("02158", "Wade Hampton Census Area, Alaska, renamed Kusilvak Census Area (02158) in 2015"),
}


def fips_bridge_note(published: str) -> str | None:
    """Provenance sentence for a bridged 1:1 FIPS rename (None for any other code)."""
    if published not in FIPS_RENAMES:
        return None
    new, what = FIPS_RENAMES[published]
    return (f"FIPS bridge: CMS MMD publishes this county under its pre-2015 code {published}; {what}, a 1:1 rename "
            f"with an unchanged boundary, so the published value is carried to {new} unchanged (Census 'Substantial "
            "Changes to Counties').")


def bridge_fips_renames(geo_id: pd.Series, geo_level: pd.Series) -> tuple[pd.Series, pd.Series]:
    """(geo_id with the 1:1 renames applied to county rows, the code as published)."""
    published = geo_id.astype(str)
    cty = geo_level.eq("county")
    new = published.where(~(cty & published.isin(FIPS_RENAMES)),
                          published.map(lambda g: FIPS_RENAMES.get(g, (g,))[0]))
    return new, published


# =========================================================================================
# Pure helpers (unit-tested)
# =========================================================================================
def normalize_fips(fips, geography: str) -> str:
    """County FIPS -> 5 chars, state -> 2 chars, nation -> 'US'. MMD mixes '4' and '04'."""
    if geography == "n":
        return "US"
    s = "" if fips is None or (isinstance(fips, float) and np.isnan(fips)) else str(fips).strip()
    if s.endswith(".0"):
        s = s[:-2]
    if not s.isdigit():
        raise ValueError(f"non-numeric FIPS {fips!r} for geography {geography!r}")
    width = 5 if geography == "c" else 2
    if len(s) > width:
        raise ValueError(f"FIPS {fips!r} too long for geography {geography!r}")
    return s.zfill(width)


def normalize_stratum(value) -> str:
    """MMD encodes 'all' as '.', '' (older files) or null."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "all"
    s = str(value).strip()
    return "all" if s in ("", ".") else s


def decode_dencat(value) -> tuple[str | None, float, float]:
    """Return (label, min, max) for an MMD denominator category; max is NaN for 10,000+."""
    s = normalize_stratum(value)
    if s == "all" or s not in DENCAT:
        return None, np.nan, np.nan
    label, lo, hi = DENCAT[s]
    return label, float(lo), float(hi)


def load_crosswalk(text: str) -> list[dict]:
    return list(csv.DictReader(io.StringIO(text)))


_XWALK_DIMS = ["population", "measure", "year", "elig", "race_code", "sex_code", "adjust", "dual"]


def resolve_source_file(crosswalk: list[dict], options: dict) -> tuple[str | None, list[str]]:
    """Replicates newDataURL() in the MMD app: first crosswalk row whose non-blank dimension
    values all contain the selected option (JavaScript String.indexOf semantics)."""
    for cand in crosswalk:
        crit = [k for k in _XWALK_DIMS if (cand.get(k) or "") != ""]
        if all(str(options[k]) in cand[k] for k in crit):
            return cand["url"], crit
    return None, []


def _or(values: list[str]) -> str:
    """OR-list filter; the app encodes 'all' ('.') as '.|IS NULL' because older files leave it blank."""
    return "|".join(list(values) + (["IS NULL"] if "." in values else []))


def build_query(source_file: str, *, year_code: str, measure: str, conditions: list[str],
                agecats: list[str], fltrs: list[str], geographies=("c", "s", "n")) -> tuple[str, dict]:
    """Build the data-API URL the MMD app would request, with OR-lists for several selections.

    The app uses the `elig` column instead of `eligcat` for 2019 FFS ED-visit files
    (population f, year 9, measure e) -- a documented one-time exception in its source.
    """
    elig_key = "elig" if (year_code == "9" and measure == "e") else "eligcat"
    params = {
        "_source": source_file,
        "year": year_code,
        "geography": "|".join(geographies),
        "measure": measure,
        "condition": "|".join(conditions),
        "sexcat": _or(["."]),
        "agecat": _or(agecats),
        "dual": _or(["."]),
        elig_key: _or(["."]),
        "racecat": _or(["."]),
        "fltr": "|".join(fltrs),
        "_size": str(API_PAGE_SIZE),
    }
    return API_BASE + "?" + urlencode(params, safe="|.", quote_via=quote), params


def parse_ccw_algorithms(page_html: str) -> pd.DataFrame:
    """Parse a CCW condition-categories page (one <tr> per algorithm) into a table."""
    rows = []
    for tr in re.findall(r"<tr.*?</tr>", page_html, flags=re.S):
        cells = [re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", c))).strip()
                 for c in re.findall(r"<t[dh].*?</t[dh]>", tr, flags=re.S)]
        if len(cells) < 5 or cells[0].lower().startswith("algorithms") or not cells[0]:
            continue
        name = re.sub(r"\s+\d+\s*(,\s*including:)?$", "", cells[0]).strip()
        icd10 = cells[3]
        codes = re.findall(ICD10_TOKEN, icd10)
        rows.append({"ccw_algorithm_name": name, "reference_period": cells[1], "icd9_text": cells[2],
                     "icd10_text": icd10, "qualifying_claims": cells[4],
                     "icd10_codes": list(dict.fromkeys(codes))})
    return pd.DataFrame(rows)


ICD10_TOKEN = r"\b[A-Z]\d[0-9A-Z](?:\.[0-9A-Z]{1,4})?\b"


def icd10_codes_in_text(text: str) -> list[str]:
    """Distinct ICD-10-CM-shaped tokens in a text, in order of first appearance."""
    return list(dict.fromkeys(re.findall(ICD10_TOKEN, text)))


def pdf_text(path: Path) -> str | None:
    """Text of a PDF via the poppler `pdftotext` binary; None when it is not installed."""
    exe = shutil.which("pdftotext")
    if exe is None or not path.exists():
        return None
    out = subprocess.run([exe, "-layout", str(path), "-"], capture_output=True, check=True)
    return out.stdout.decode("utf-8", errors="replace")


def merge_search_codes(base: list[str], extra) -> list[str]:
    """base search prefixes + sourced codes not already covered by one of them (prefix match)."""
    out = list(base)
    for code in list(extra or []):
        if not any(code == s or code.startswith(s) for s in out):
            out.append(code)
    return out


def codes_matching(search_codes: list[str], algo_codes: list[str]) -> list[str]:
    """Algorithm codes equal to, or nested under, any search code (prefix match on ICD-10-CM)."""
    hits = []
    for code in algo_codes:
        if any(code == s or code.startswith(s) for s in search_codes):
            hits.append(code)
    return hits


def _norm_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", s.lower().replace("–", "-")).strip()


def evidence_level_for(condition_id: str) -> str:
    for cid, _code, level in BURDEN_MAPPINGS:
        if cid == condition_id:
            return level
    return "D"


# =========================================================================================
# Retrieval
# =========================================================================================
def fetch_reference_files() -> dict[str, Path]:
    raw = raw_dir(SOURCE_ID)
    out = {}
    for rel, url in REFERENCE_FILES.items():
        (raw / rel).parent.mkdir(parents=True, exist_ok=True)
        try:
            out[rel] = download_file(url, SOURCE_ID, rel, allow_html=rel in HTML_REFERENCES,
                                     headers={"Accept": "*/*"})
        except Exception as exc:  # recorded, not fatal for optional documents
            print(f"[cms_mmd] reference {rel} failed: {exc}")
            if rel in REQUIRED_REFERENCES:
                raise
    return out


def catalog_check(catalog_json: str) -> dict:
    """Search the data.cms.gov DCAT catalog for an MMD / chronic-condition prevalence dataset."""
    ds = json.loads(catalog_json).get("dataset", [])
    pat = re.compile(r"disparit|mapping medicare|chronic condition|fibromyalg|prevalence", re.I)
    hits = [d["title"] for d in ds if pat.search(d.get("title", "") + " " + (d.get("description") or ""))]
    gv = [d["title"] for d in ds if "geographic variation" in d.get("title", "").lower()]
    return {"n_datasets": len(ds), "mmd_or_prevalence_hits": hits, "geographic_variation_titles": gv}


def verify_condition_labels(menus_js: str, strings: dict) -> dict[str, str]:
    """Check that the MMD condition codes used here still carry the expected menu labels."""
    seg = menus_js[menus_js.index("conditions = ["):]
    rows = re.findall(r"get_string_by_id\(\"(\w+)\"\), 'val': '(\d*)', 'header': \w+, "
                      r"'domain': \w+, 'measure': '(\w+)'", seg)
    labels: dict[str, set] = {}
    for sid, val, meas in rows:
        labels.setdefault(val, set()).add(strings.get(sid, {}).get("english", "").strip())
    found = {}
    for code, expected in MMD_CONDITIONS.items():
        if code == "134":
            ok = "'val': \"134\"" in menus_js and "COVID-19" in menus_js
            found[code] = "COVID-19" if ok else ""
        elif code == "10":
            ok = set(MMD_CONDITION_10_LABELS.values()) <= labels.get("10", set())
            found[code] = "; ".join(sorted(labels.get("10", set())))
        else:
            ok = expected in labels.get(code, set())
            found[code] = "; ".join(sorted(labels.get(code, set())))
        if not ok:
            raise RuntimeError(f"MMD condition code {code} no longer labelled {expected!r}: {found[code]!r}")
    return found


def plan_requests(crosswalk: list[dict]) -> list[dict]:
    """One request per (claims year, measure, source file).

    The app resolves the file with the *selected* adjustment, and some measures are split into
    one file per adjustment (e.g. principal cost 2012-2018: `..._fltr_1_...` and `..._fltr_2_...`).
    Resolving once with adjust=1 and asking that file for fltr 1|2 silently returns fltr 1 only,
    so each adjustment is resolved separately and adjustments sharing a file share one request.
    """
    plan = []
    for yc in YEAR_CODES:
        for meas, (conds, ages, fltrs) in QUERY_PLAN.items():
            by_file: dict[tuple, list[str]] = {}
            for fltr in fltrs:
                opts = {"population": "f", "measure": meas, "year": yc, "elig": ".", "race_code": ".",
                        "sex_code": ".", "adjust": fltr, "dual": "."}
                src, crit = resolve_source_file(crosswalk, opts)
                by_file.setdefault((src, tuple(crit)), []).append(fltr)
            for (src, crit), file_fltrs in by_file.items():
                item = {"year_code": yc, "year": YEAR_CODES[yc], "measure": meas, "conditions": conds,
                        "agecats": ages, "fltrs": file_fltrs, "source_file": src, "crosswalk_keys": list(crit)}
                if src:
                    url, params = build_query(src, year_code=yc, measure=meas, conditions=conds,
                                              agecats=ages, fltrs=file_fltrs)
                    item.update(url=url, params=params,
                                raw_file=f"api/f_{YEAR_CODES[yc]}_{meas}__{src}.json")
                plan.append(item)
    return plan


def fetch_api_slices(plan: list[dict]) -> list[dict]:
    raw = raw_dir(SOURCE_ID)
    (raw / "api").mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(SOURCE_ID)["files"]
    for item in plan:
        if not item.get("url"):
            item["status"] = "no_source_file_in_crosswalk"
            continue
        cached = item["raw_file"] in manifest and (raw / item["raw_file"]).exists()
        # a cached slice fetched with a different query (e.g. other fltr list) is re-fetched so the
        # MANIFEST url always describes the bytes on disk
        stale = cached and manifest[item["raw_file"]].get("url") != item["url"]
        try:
            path = download_file(item["url"], SOURCE_ID, item["raw_file"], refresh=stale,
                                 headers={"Accept": "application/json"})
            data = json.loads(path.read_text())
            if not isinstance(data, list):
                raise ValueError(f"unexpected payload type {type(data).__name__}")
            item["status"] = "ok" if data else "empty"
            item["n_rows"] = len(data)
            if len(data) >= API_PAGE_SIZE:
                raise RuntimeError(f"{item['raw_file']}: {len(data)} rows = page size; response may be truncated")
            got = {str(r.get("fltr")) for r in data}
            item["fltrs_missing"] = [f for f in item["fltrs"] if f not in got]
            if item["fltrs_missing"]:
                print(f"[cms_mmd] WARNING {item['raw_file']}: requested fltr {item['fltrs_missing']} not in response")
        except Exception as exc:
            item["status"] = f"error: {exc}"
            print(f"[cms_mmd] {item['year']} {item['measure']} failed: {exc}")
        if not cached or stale:
            time.sleep(0.5)
        print(f"[cms_mmd] {item['year']} {item['measure']:>1} fltr={'|'.join(item['fltrs'])} {item['status']:<6} "
              f"{item.get('n_rows', 0):>7} rows  {item['source_file']}", flush=True)
    return plan


# =========================================================================================
# Parsing
# =========================================================================================
def parse_slices(plan: list[dict]) -> pd.DataFrame:
    raw = raw_dir(SOURCE_ID)
    manifest = load_manifest(SOURCE_ID)["files"]
    frames = []
    for item in plan:
        if item.get("status") != "ok":
            continue
        recs = json.loads((raw / item["raw_file"]).read_text())
        df = pd.DataFrame.from_records(recs)
        if "eligcat" not in df.columns and "elig" in df.columns:
            df["eligcat"] = df["elig"]
        for col in ("sexcat", "agecat", "dual", "racecat", "eligcat", "dencat", "rate", "fltr"):
            if col not in df.columns:
                df[col] = None
        df["raw_file"] = item["raw_file"]
        df["source_file"] = item["source_file"]
        df["retrieved_at"] = manifest[item["raw_file"]]["retrieved_at"]
        df["request_year_code"] = item["year_code"]
        frames.append(df[["year", "geography", "fips", "measure", "condition", "sexcat", "agecat", "dual",
                          "racecat", "eligcat", "fltr", "dencat", "rate", "raw_file", "source_file",
                          "retrieved_at", "request_year_code"]])
    df = pd.concat(frames, ignore_index=True)
    # guard: every row must belong to the requested year
    bad = df["year"].astype(str) != df["request_year_code"].astype(str)
    if bad.any():
        raise RuntimeError(f"{int(bad.sum())} rows carry a year different from the request")
    for col in ("sexcat", "agecat", "dual", "racecat", "eligcat", "fltr", "dencat"):
        df[col] = df[col].map(normalize_stratum)
    # keep only the all-sex/race/dual/eligibility cells (the filters requested exactly these)
    keep = (df["sexcat"] == "all") & (df["racecat"] == "all") & (df["dual"] == "all") & (df["eligcat"] == "all")
    df = df[keep].copy()
    df["geo_id"] = [normalize_fips(f, g) for f, g in zip(df["fips"], df["geography"])]
    df["geo_level"] = df["geography"].map(GEOGRAPHY_CODES)
    df["year"] = df["request_year_code"].map(YEAR_CODES).astype(int)
    df["value"] = pd.to_numeric(df["rate"], errors="coerce")
    df["adjustment"] = df["fltr"].map(FLTR_CODES)
    df["stratum_age"] = df["agecat"].map(AGE_CODES)
    dec = df["dencat"].map(decode_dencat)
    df["denominator_category"] = pd.to_numeric(df["dencat"].where(df["dencat"] != "all"), errors="coerce").astype("Int64")
    df["denominator_range"] = [d[0] for d in dec]
    df["denominator_min"] = [d[1] for d in dec]
    df["denominator_max"] = [d[2] for d in dec]
    key = ["geo_id", "geo_level", "year", "measure", "condition", "adjustment", "stratum_age"]
    dup = df.duplicated(key, keep=False)
    if dup.any():
        raise RuntimeError(f"{int(dup.sum())} duplicated MMD cells on {key}")
    return df.reset_index(drop=True)


def attach_geo_match(df: pd.DataFrame, geos: pd.DataFrame) -> pd.DataFrame:
    cty = geos[geos["geo_level"] == "county"][["geo_id", "ct_legacy"]]
    st = set(geos.loc[geos["geo_level"] == "state", "geo_id"])
    df = df.copy()
    df["geo_id"], df["fips_as_published"] = bridge_fips_renames(df["geo_id"], df["geo_level"])
    df["fips_bridge"] = df["fips_as_published"].map(
        lambda g: f"1:1 FIPS rename {g}->{FIPS_RENAMES[g][0]}" if g in FIPS_RENAMES else None)
    df.loc[df["geo_level"] != "county", "fips_bridge"] = None
    key = ["geo_id", "geo_level", "year", "measure", "condition", "adjustment", "stratum_age"]
    dup = df.duplicated(key, keep=False)
    if dup.any():   # CMS published both the old and the new code for the same cell: never merge them silently
        raise RuntimeError(f"{int(dup.sum())} MMD cells collide after the FIPS rename bridge: "
                           f"{sorted(df.loc[dup, 'fips_as_published'].unique())[:10]}")
    df = df.merge(cty.rename(columns={"ct_legacy": "_ct_legacy"}), on="geo_id", how="left")
    status = np.where(df["geo_level"] == "national", "national",
             np.where(df["geo_level"] == "state", np.where(df["geo_id"].isin(st), "matched_2024", "unmatched"),
             np.where(df["_ct_legacy"].isna(), "unmatched",
             np.where(df["_ct_legacy"].astype(bool), "matched_ct_legacy", "matched_2024"))))
    df["geo_match_status"] = status
    df["unmatched_reason"] = np.where(
        df["geo_match_status"] != "unmatched", None,
        np.where(df["geo_id"].str.fullmatch(r"\d{2}990"),
                 "cms_code_ending_990_not_a_census_county (undocumented in MMD tech doc)",
                 "fips_not_in_2024_census_vintage"))
    df = df.drop(columns="_ct_legacy")
    df["state_fips"] = np.where(df["geo_level"] == "national", None, df["geo_id"].str[:2])
    df["county_fips"] = np.where(df["geo_level"] == "county", df["geo_id"], None)
    return df


def _data_status(source_file: str) -> str:
    if source_file.endswith("_p"):
        return "preliminary (source file suffix _p)"
    return "final"


# =========================================================================================
# Output tables
# =========================================================================================
def build_burden(df: pd.DataFrame) -> pd.DataFrame:
    prev = df[(df["measure"] == "v") & df["condition"].isin(BURDEN_SOURCE_CONDITIONS)].copy()
    notes = _mapping_notes()
    parts = []
    for cid, code, level in BURDEN_MAPPINGS:
        p = prev[prev["condition"] == code].copy()
        p["condition_id"] = cid
        p["burden_evidence_level"] = level
        p["provenance_notes"] = f"{notes[cid]} | {LIMITATION} | {SUPPRESSION_NOTE} | {ECOLOGICAL}"
        parts.append(p)
    b = pd.concat(parts, ignore_index=True)
    b["provenance_notes"] = _with_bridge_note(b)
    b["source_condition_code"] = b["condition"]
    b["source_condition_name"] = b["condition"].map(MMD_CONDITIONS)
    b["ccw_algorithm_name"] = b["condition"].map(CCW_ALGORITHM_FOR_CODE)
    b["measure_code"] = "v"
    b["measure"] = "prevalence"
    b["unit"] = MEASURES["v"][1]
    b["prevalence_pct"] = b["value"]  # `value` is kept too, for union with other burden partitions
    b["value_zero_possibly_suppressed"] = b["value"].eq(0)
    b["population"] = "Medicare FFS beneficiaries continuously enrolled in Parts A and B, no MA enrollment"
    b["stratum_sex"] = "all"
    b["stratum_race_ethnicity"] = "all"
    b["stratum_dual"] = "all"
    b["stratum_eligibility"] = "all"
    b["default_view"] = (b["adjustment"] == "unsmoothed_actual") & (b["stratum_age"] == "all")
    b["claims_data_status"] = b["source_file"].map(_data_status)
    b["source_version"] = "MMD data API file " + b["source_file"]
    out = pd.concat([_with_prov(g, res) for res, g in b.groupby("geo_level")], ignore_index=True)
    cols = ["geo_id", "geo_level", "state_fips", "county_fips", "geo_match_status", "unmatched_reason",
            "fips_as_published", "fips_bridge", "year", "condition_id", "burden_evidence_level", "source_condition_code", "source_condition_name",
            "ccw_algorithm_name", "measure", "prevalence_pct", "value", "unit", "adjustment", "stratum_age", "stratum_sex",
            "stratum_race_ethnicity", "stratum_dual", "stratum_eligibility", "default_view",
            "denominator_category", "denominator_range", "denominator_min", "denominator_max",
            "value_zero_possibly_suppressed", "population", "claims_data_status", "raw_file"]
    return out[cols + [c for c in out.columns if c not in cols and c in _PROV]].sort_values(
        ["condition_id", "year", "geo_level", "geo_id", "adjustment", "stratum_age"]).reset_index(drop=True)


_PROV = ["data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
         "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes"]


def _record_id(d: pd.DataFrame) -> pd.Series:
    """The MMD cell as CMS published it (a bridged rename keeps the published FIPS here)."""
    cond = d["source_condition_code"] if "source_condition_code" in d else d["condition"]
    fips = d["fips_as_published"] if "fips_as_published" in d else d["geo_id"]
    return ("mmd:f:" + d["year"].astype(str) + ":" + d["measure_code"] + ":" + cond.astype(str) + ":"
            + d["geo_level"] + ":" + fips + ":" + d["adjustment"] + ":age=" + d["stratum_age"])


def _with_bridge_note(d: pd.DataFrame) -> pd.Series:
    """provenance_notes with the FIPS-bridge sentence prepended on bridged rows."""
    note = d["fips_as_published"].map(fips_bridge_note).where(d["geo_level"] == "county")
    return (note.fillna("") + np.where(note.notna(), " | ", "") + d["provenance_notes"]).astype(str)


def _with_prov(g: pd.DataFrame, geo_level: str) -> pd.DataFrame:
    """Attach provenance columns; g must carry measure_code, source_version, retrieved_at, provenance_notes."""
    g = g.copy()
    level = g["burden_evidence_level"] if "burden_evidence_level" in g else "context (not a burden estimate)"
    out = add_provenance(
        g.drop(columns=[c for c in ("retrieved_at", "source_version", "provenance_notes") if c in g]),
        data_layer="geographic", source_name=SOURCE_NAME, source_version="",
        retrieved_at="", evidence_type="administrative_claims", source_record_id=_record_id,
        source_geographic_resolution=geo_level, evidence_level=level, provenance_notes="")
    out["source_version"] = g["source_version"].values
    out["retrieved_at"] = g["retrieved_at"].values
    out["provenance_notes"] = g["provenance_notes"].values
    return out


def build_context(df: pd.DataFrame) -> pd.DataFrame:
    ctx = df[~((df["measure"] == "v") & df["condition"].isin(BURDEN_SOURCE_CONDITIONS))].copy()
    ctx["measure_code"] = ctx["measure"]
    ctx["measure"] = ctx["measure_code"].map(lambda m: MEASURES[m][0])
    ctx["unit"] = ctx["measure_code"].map(lambda m: MEASURES[m][1])
    ctx.loc[(ctx["condition"] == "24") & (ctx["measure_code"] == "v"), "unit"] = \
        "percent of beneficiaries entitled by disability (original or current reason)"
    ctx["context_category"] = ctx["measure_code"].map(CONTEXT_CATEGORY)
    ctx.loc[(ctx["condition"] == "134") & (ctx["measure_code"] == "v"), "context_category"] = "infection_exposure"
    ctx["source_condition_code"] = ctx["condition"]
    ctx["source_condition_name"] = ctx["condition"].map(MMD_CONDITIONS)
    ctx.loc[ctx["condition"] == "10", "source_condition_name"] = ctx["measure_code"].map(MMD_CONDITION_10_LABELS)
    ctx["related_condition_ids"] = ctx["condition"].map({"51": "fibromyalgia;me_cfs", "59": "migraine"})
    notes = {
        "51": "Condition-specific rate/cost for CCW 'Fibromyalgia, Chronic Pain and Fatigue'.",
        "59": "Condition-specific rate/cost for CCW 'Migraine and Other Chronic Headache'.",
        "24": "Beneficiaries whose original or current reason for entitlement is disability (MBSF); "
              "describes the <65 disability-entitled share of the FFS population, not a disease.",
        "134": "COVID-19 diagnosis on any claim (CMS Medicare COVID-19 snapshot methodology). Acute infection "
               "context ONLY: this is NOT Long COVID / post-COVID condition (U09.9) and must not be used as its "
               "burden.",
        "10": "All-cause rate among all FFS beneficiaries.",
    }
    measure_notes = {
        "h": "Hospitalization rate = inpatient discharges per 1,000 FFS beneficiaries (principal diagnosis for "
             "condition-specific rates; Part A continuous enrollment).",
        "e": "ED visit rate = ED visits per 1,000 FFS beneficiaries (inpatient + outpatient files).",
        "t": "Average total cost = mean annual Medicare spending, all claims, among beneficiaries with the condition.",
        "a": "Risk-adjusted total cost = average HCC risk score x fixed standard cost $9,276.26.",
        "p": "Principal cost = mean annual cost of claims whose principal diagnosis is the condition.",
        "v": "Prevalence.",
    }
    ctx["value_zero_possibly_suppressed"] = ctx["value"].eq(0) & ctx["measure_code"].isin(["v", "h", "e"])
    ctx["provenance_notes"] = (ctx["condition"].map(notes) + " " + ctx["measure_code"].map(measure_notes)
                               + " | " + LIMITATION + " | " + SUPPRESSION_NOTE + " | " + ECOLOGICAL)
    ctx["provenance_notes"] = _with_bridge_note(ctx)
    ctx["population"] = "Medicare FFS beneficiaries (MMD analysis population for the measure)"
    ctx["claims_data_status"] = ctx["source_file"].map(_data_status)
    ctx["source_version"] = "MMD data API file " + ctx["source_file"]
    out = pd.concat([_with_prov(g, res) for res, g in ctx.groupby("geo_level")], ignore_index=True)
    cols = ["geo_id", "geo_level", "state_fips", "county_fips", "geo_match_status", "unmatched_reason",
            "fips_as_published", "fips_bridge", "year", "context_category",
            "measure", "measure_code", "source_condition_code", "source_condition_name", "related_condition_ids",
            "value", "unit", "adjustment", "stratum_age", "denominator_category", "denominator_range",
            "denominator_min", "denominator_max", "value_zero_possibly_suppressed", "population",
            "claims_data_status", "raw_file"]
    return out[cols + _PROV].sort_values(["measure", "source_condition_code", "year", "geo_level", "geo_id",
                                          "adjustment"]).reset_index(drop=True)


def _mapping_notes() -> dict[str, str]:
    return {
        "fibromyalgia": (
            "burden_evidence_level B (closely matching coded condition): CCW 'Fibromyalgia, Chronic Pain and "
            "Fatigue' contains fibromyalgia M79.7 but is a composite that also counts chronic pain (G89.2x/G89.3/"
            "G89.4), radiculopathy (M54.1x), myositis (M60.8x/M60.9), myalgia (M79.1x), neuralgia (M79.2) and "
            "chronic fatigue (R53.82); it therefore overstates fibromyalgia alone. 2-year lookback; >=1 inpatient "
            "or >=2 other non-drug claims."),
        "me_cfs": (
            "burden_evidence_level C (symptom/comorbidity proxy): the same CCW 'Fibromyalgia, Chronic Pain and "
            "Fatigue' composite contains R53.82 (chronic fatigue, unspecified) but NOT G93.3x (G93.32 ME/CFS, "
            "G93.31 postviral fatigue syndrome); it is not an ME/CFS prevalence."),
        "migraine": (
            "burden_evidence_level B (closely matching coded condition): CCW 'Migraine and Other Chronic Headache' "
            "= migraine G43.x plus other headache syndromes G44.x (cluster, tension-type, post-traumatic, "
            "drug-induced, other); 2-year lookback; >=1 inpatient or >=2 other non-drug claims."),
    }


def build_coverage(conditions: list[dict], ccw: pd.DataFrame, mmd_names: set[str],
                   burden: pd.DataFrame, chronic30_codes: list[str] | None = None,
                   registry_codes: dict[str, list[str]] | None = None) -> pd.DataFrame:
    """chronic30_codes: ICD-10 tokens of the CCW 30-chronic-condition PDF (None = not searched).
    registry_codes: condition_id -> ICD-10-CM codes from the sourced condition_registry table,
    merged into the search keys so the search is at least as broad as the ontology layer."""
    exposed = ccw[ccw["ccw_algorithm_name"].map(_norm_name).isin({_norm_name(n) for n in mmd_names})]
    exposed_names = set(exposed["ccw_algorithm_name"])
    notes = _mapping_notes()
    registry_codes = registry_codes or {}
    searched_sets = (f"CCW 'other chronic or potentially disabling' algorithms ({len(ccw)} parsed from the CCW HTML "
                     f"table, {len(exposed_names)} exposed by MMD)")
    if chronic30_codes is not None:
        searched_sets += (f"; CCW 30 chronic-condition algorithms (the set MMD's primary chronic conditions are "
                          f"drawn from; current revision, MMD used the 27-condition set before 2021; "
                          f"{len(chronic30_codes)} distinct ICD-10 codes text-extracted from "
                          "ccw/chr-chronic-condition-algorithms.pdf, not attributable to a single algorithm)")
    else:
        searched_sets += ("; CCW 30 chronic-condition algorithms NOT searched (PDF only; pdftotext unavailable "
                          "or PDF not retrieved)")
    rows = []
    for c in conditions:
        cid = c["id"]
        level = evidence_level_for(cid)
        search = merge_search_codes(TARGET_SEARCH_CODES.get(cid, []), registry_codes.get(cid))
        hits_exposed, hits_other = [], []
        for _, a in ccw.iterrows():
            h = codes_matching(search, a["icd10_codes"])
            if h:
                tgt = hits_exposed if a["ccw_algorithm_name"] in exposed_names else hits_other
                tgt.append(f"{a['ccw_algorithm_name']}: {', '.join(h)}")
        hits30 = codes_matching(search, chronic30_codes) if chronic30_codes is not None and search else []
        code = next((code for i, code, _ in BURDEN_MAPPINGS if i == cid), None)
        algo = ccw[ccw["ccw_algorithm_name"] == CCW_ALGORITHM_FOR_CODE[code]].iloc[0] if code else None
        sub = burden[burden["condition_id"] == cid]
        scope = ("any CCW algorithm behind an MMD condition (other-chronic table and the 30 chronic conditions)"
                 if chronic30_codes is not None else
                 "any CCW 'other chronic' algorithm exposed by MMD (the 30 chronic-condition algorithms were not "
                 "searched)")
        if code:
            rationale = notes[cid]
        elif cid == "long_covid":
            rationale = (f"burden_evidence_level D: U09.9 (post COVID-19 condition) does not occur in {scope}. MMD "
                         "condition 'COVID-19' (134) counts any COVID-19 diagnosis (acute infection) and is kept only "
                         "in geo_context__cms_mmd as infection-exposure context, never as Long COVID burden.")
        elif c.get("grouping_only"):
            rationale = "burden_evidence_level D: grouping concept; no CMS/CCW measure."
        else:
            rationale = (f"burden_evidence_level D: none of the searched ICD-10-CM codes ({', '.join(search)}) occurs "
                         f"in {scope}; no usable CMS burden estimate.")
            if cid in ("pots", "dysautonomia", "long_covid", "eds_hsd"):
                rationale += (" The fibromyalgia/chronic-pain/fatigue composite is deliberately NOT used as a proxy: "
                              "it has no autonomic, post-infectious or connective-tissue component.")
        if hits30 and not code:
            rationale += (f" REVIEW: {', '.join(hits30)} occur(s) in the CCW 30 chronic-condition algorithms PDF; "
                          "check which algorithm before relying on level D.")
        rows.append({
            "condition_id": cid, "preferred_name": c["preferred_name"], "burden_evidence_level": level,
            "mmd_condition_code": code, "mmd_condition_name": MMD_CONDITIONS.get(code) if code else None,
            "ccw_algorithm_name": algo["ccw_algorithm_name"] if algo is not None else None,
            "ccw_reference_period": algo["reference_period"] if algo is not None else None,
            "ccw_qualifying_claims": algo["qualifying_claims"] if algo is not None else None,
            "ccw_icd10_codes": ", ".join(algo["icd10_codes"]) if algo is not None else None,
            "ccw_n_icd10_codes": len(algo["icd10_codes"]) if algo is not None else 0,
            "icd10_codes_searched": ", ".join(search),
            "search_hits_in_mmd_exposed_ccw_algorithms": "; ".join(hits_exposed) or None,
            "search_hits_in_ccw_algorithms_not_in_mmd": "; ".join(hits_other) or None,
            "search_hits_in_ccw_30_chronic_conditions": (", ".join(hits30) or "none")
            if chronic30_codes is not None else "not searched",
            "ccw_algorithm_sets_searched": searched_sets,
            "n_burden_rows": int(len(sub)),
            "years_available": (f"{sub['year'].min()}-{sub['year'].max()}" if len(sub) else None),
            "rationale": rationale,
        })
    cov = pd.DataFrame(rows)
    manifest = load_manifest(SOURCE_ID)["files"]
    return add_provenance(
        cov, data_layer="ontology", source_name="CCW condition algorithms (www2.ccwdata.org) x MMD condition menu",
        source_version="CCW condition-categories-other/chronic pages as retrieved; MMD menus.js as retrieved",
        retrieved_at=manifest["ccw/condition-categories-other.html"]["retrieved_at"],
        evidence_type="ontology_mapping", source_record_id="condition_id", source_geographic_resolution="none",
        evidence_level=cov["burden_evidence_level"],
        provenance_notes=cov["rationale"] + " | " + LIMITATION)


# =========================================================================================
# Audit statistics, registry entry and DATA_AUDIT.md
# =========================================================================================
def compute_stats(df: pd.DataFrame, burden: pd.DataFrame, context: pd.DataFrame, cov: pd.DataFrame,
                  plan: list[dict], geos: pd.DataFrame) -> dict:
    n2024 = int(((geos["geo_level"] == "county") & (~geos["ct_legacy"])).sum())
    cty2024 = set(geos.loc[(geos["geo_level"] == "county") & (~geos["ct_legacy"]), "geo_id"])
    s: dict = {"requests": [{k: it.get(k) for k in ("year", "measure", "fltrs", "source_file", "status", "n_rows",
                                                     "fltrs_missing", "raw_file")} for it in plan],
               "n_counties_2024": n2024}
    s["rows"] = {"parsed_cells": int(len(df)), "burden": int(len(burden)), "context": int(len(context)),
                 "coverage": int(len(cov))}
    # geography match per year (prevalence, condition 51, unsmoothed actual, all ages)
    per_year = []
    base = df[(df["measure"] == "v") & (df["condition"] == "51") & (df["adjustment"] == "unsmoothed_actual")
              & (df["stratum_age"] == "all")]
    for y, g in base.groupby("year"):
        c = g[g["geo_level"] == "county"]
        unmatched = sorted(c.loc[c["geo_match_status"] == "unmatched", "geo_id"])
        um990 = [u for u in unmatched if u.endswith("990")]
        per_year.append({
            "year": int(y), "county_rows": int(len(c)),
            "matched_2024": int((c["geo_match_status"] == "matched_2024").sum()),
            "matched_ct_legacy": int((c["geo_match_status"] == "matched_ct_legacy").sum()),
            "ct_planning_region_rows": int(c["geo_id"].str.match(r"^091[1-9]0$").sum()),
            "unmatched": len(unmatched), "unmatched_cms_990": len(um990),
            "unmatched_geo_ids": [u for u in unmatched if u not in um990],
            "bridged_renames": sorted(f"{a}->{b}" for a, b in c.loc[c["fips_bridge"].notna(),
                                                                    ["fips_as_published", "geo_id"]].values),
            "pct_rows_matched_2024": round(100 * (c["geo_match_status"] == "matched_2024").mean(), 2),
            "pct_2024_counties_covered": round(100 * len(cty2024 & set(c["geo_id"])) / n2024, 2),
            "state_rows": int((g["geo_level"] == "state").sum()),
            "state_unmatched": sorted(g.loc[(g["geo_level"] == "state") & (g["geo_match_status"] == "unmatched"),
                                            "geo_id"]),
            "national_rows": int((g["geo_level"] == "national").sum()),
        })
    s["geo_match_by_year"] = per_year
    latest = int(base["year"].max())
    have = set(base.loc[(base["year"] == latest) & (base["geo_level"] == "county"), "geo_id"])
    c24 = geos[(geos["geo_level"] == "county") & (~geos["ct_legacy"])]
    miss = c24[~c24["geo_id"].isin(have)]
    s["counties_2024_without_value_latest"] = {
        "year": latest, "n": int(len(miss)),
        "by_state": {k: sorted(v) for k, v in miss.groupby("state_abbr")["geo_id"]}}
    c990 = base[(base["geo_level"] == "county") & base["geo_id"].str.fullmatch(r"\d{2}990")]
    s["cms_990"] = {int(y): {"n": int(len(g)),
                             "denominator_ranges": {str(k): int(v) for k, v in
                                                    g["denominator_range"].value_counts(dropna=False).items()}}
                    for y, g in c990.groupby("year")}
    ak = base[(base["geo_level"] == "county") & (base["geo_id"].str[:2] == "02")]
    s["alaska"] = {"codes_by_year": {int(y): int(len(g)) for y, g in ak.groupby("year")},
                   "codes_identical_all_years": bool(ak.groupby("year")["geo_id"].apply(frozenset).nunique() == 1),
                   "n_2024_alaska": int((c24["state_fips"] == "02").sum())}
    # prevalence summaries per condition-year (unsmoothed actual, all ages)
    summ = []
    b = burden[(burden["adjustment"] == "unsmoothed_actual") & (burden["stratum_age"] == "all")
               & (burden["condition_id"].isin(["fibromyalgia", "migraine"]))]
    for (cid, y), g in b.groupby(["condition_id", "year"]):
        c = g[g["geo_level"] == "county"]
        nat = g.loc[g["geo_level"] == "national", "prevalence_pct"]
        st = g[g["geo_level"] == "state"]["prevalence_pct"]
        summ.append({
            "condition_id": cid, "source_condition": g["source_condition_name"].iloc[0], "year": int(y),
            "national_pct": float(nat.iloc[0]) if len(nat) else None,
            "state_min": float(st.min()), "state_max": float(st.max()),
            "county_n": int(len(c)), "county_value_missing": int(c["prevalence_pct"].isna().sum()),
            "county_zero": int(c["value_zero_possibly_suppressed"].sum()),
            "county_median": float(c["prevalence_pct"].median()),
            "county_p5": float(c["prevalence_pct"].quantile(0.05)),
            "county_p95": float(c["prevalence_pct"].quantile(0.95)),
            "counties_2024_without_value": int(n2024 - len(cty2024 & set(c.loc[c["prevalence_pct"].notna(), "geo_id"]))),
        })
    s["prevalence_summary"] = summ
    # <65 stratum and age-standardized coverage
    lt65 = burden[(burden["stratum_age"] == "<65") & (burden["geo_level"] == "county")
                  & (burden["condition_id"] == "fibromyalgia")]
    s["lt65_county_rows_by_year"] = {int(k): int(v) for k, v in lt65.groupby("year").size().items()}
    nat_lt65 = burden[(burden["stratum_age"] == "<65") & (burden["geo_level"] == "national")
                      & (burden["adjustment"] == "unsmoothed_actual")]
    s["national_lt65_pct"] = {f"{r.condition_id}:{r.year}": r.prevalence_pct for r in nat_lt65.itertuples()
                              if r.condition_id in ("fibromyalgia", "migraine")}
    # denominator categories (county, 51, actual, all ages)
    dc = burden[(burden["geo_level"] == "county") & (burden["condition_id"] == "fibromyalgia")
                & (burden["adjustment"] == "unsmoothed_actual") & (burden["stratum_age"] == "all")]
    s["denominator_category_counts_county_fibro"] = {
        int(y): {str(k): int(v) for k, v in g["denominator_range"].value_counts(dropna=False).items()}
        for y, g in dc.groupby("year")}
    # context summary
    cs = context[(context["adjustment"] == "unsmoothed_actual") & (context["stratum_age"] == "all")]
    ctx_sum = []
    for (m, cnd), g in cs.groupby(["measure", "source_condition_code"]):
        c = g[g["geo_level"] == "county"]
        nat = g[g["geo_level"] == "national"].sort_values("year")
        ctx_sum.append({"measure": m, "condition": cnd, "name": g["source_condition_name"].iloc[0],
                        "years": sorted(int(v) for v in g["year"].unique()),
                        "county_rows": int(len(c)), "county_missing_value": int(c["value"].isna().sum()),
                        "county_zero": int((c["value"] == 0).sum()),
                        "national_latest": (float(nat["value"].iloc[-1]) if len(nat) else None),
                        "national_latest_year": (int(nat["year"].iloc[-1]) if len(nat) else None)})
    s["context_summary"] = ctx_sum
    s["context_rows_by_measure"] = {str(k): int(v) for k, v in context.groupby("measure").size().items()}
    s["burden_rows_by_condition"] = {str(k): int(v) for k, v in burden.groupby("condition_id").size().items()}
    s["coverage"] = cov[["condition_id", "burden_evidence_level", "n_burden_rows"]].to_dict("records")
    return s


def registry_entry(stats: dict, retrieved: str) -> dict:
    ok = [r for r in stats["requests"] if r["status"] == "ok"]
    failed = [r for r in stats["requests"] if r["status"] not in ("ok",)]
    return {
        "source_id": SOURCE_ID,
        "name": "CMS Mapping Medicare Disparities (MMD) Tool - Population View",
        "publisher": "CMS Office of Minority Health (tool built and maintained by NORC for CMS OMH)",
        "landing_url": LANDING_URL,
        "access_urls": [API_BASE, APP_URL] + [u for u in REFERENCE_FILES.values() if u != APP_URL],
        "license": "U.S. federal government work (public domain); aggregate, CMS cell-suppressed public data",
        "access_conditions": ("Open, no registration or data-use agreement. MMD is not a catalogued data.cms.gov "
                              "dataset (absent from data.cms.gov/data.json at retrieval); aggregate slices are "
                              "served by the data API that backs the public MMD app, queried exactly as the app "
                              "does."),
        "retrieved_at": retrieved,
        "source_version": ("MMD Population View data API (Medicare FFS) files as resolved by the app's "
                           "codebook_crosswalk.csv; claims years 2012-2022 final, 2023 preliminary (_p files); "
                           "technical documentation revision 6/26/2024"),
        "update_date": "Periodic (roughly annual new claims year plus preliminary-to-final refresh)",
        "data_layer": "geographic",
        "unit_of_observation": ("geography (county / state-territory / nation) x claims year x MMD condition x "
                                "measure x adjustment x age stratum; aggregate of Medicare FFS beneficiaries"),
        "sample_size": {
            "api_slices_ok": len(ok), "api_slices_not_ok": len(failed),
            "parsed_cells": stats["rows"]["parsed_cells"],
            "geo_condition_burden_rows": stats["rows"]["burden"],
            "geo_context_rows": stats["rows"]["context"],
            "coverage_rows": stats["rows"]["coverage"],
            "county_rows_per_year_fibro_prevalence": {r["year"]: r["county_rows"] for r in stats["geo_match_by_year"]},
            "beneficiary_denominators": "published only as categories (11-499 ... 10,000+)",
        },
        "geographic_resolution": "county (FIPS as published by CMS; legacy CT counties), state/territory, national",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (aggregate claims statistics; no person identifiers)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested" if not failed else "partial",
        "processed_outputs": ["geo_condition_burden__cms_mmd", "geo_context__cms_mmd",
                              "condition_burden_coverage__cms_mmd"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.cms_mmd",
        "limitations": [
            LIMITATION,
            "Only fibromyalgia (B, via CCW 'Fibromyalgia, Chronic Pain and Fatigue'), migraine (B, via CCW 'Migraine "
            "and Other Chronic Headache') and ME/CFS (C, symptom proxy via the same fibromyalgia/pain/fatigue "
            "composite) have CMS burden rows; POTS, dysautonomia, Long COVID, EDS/HSD, MCAS, Lyme, PTLDS, IBS, "
            "gastroparesis and endometriosis are D (no usable CMS estimate).",
            "MMD 'COVID-19' is any COVID-19 diagnosis (acute infection), not Long COVID; kept as context only.",
            "Beneficiary denominators are published only as 5 size categories; prevalence is an integer percent.",
            "Suppression: cells <11 beneficiaries omitted; numerators of 1-2 published as 0 %.",
            "County FIPS follow CMS vintage (legacy Connecticut counties); unmatched rows are reported, not forced. "
            "The two 1:1 renames CMS still publishes under pre-2015 codes (46113 -> 46102 Oglala Lakota SD, 02270 -> "
            "02158 Kusilvak AK; same boundary) are bridged to the 2024 code; fips_as_published keeps the CMS code.",
            "Smoothed (empirical-Bayes) rates, sex/race/dual/eligibility strata and the Medicare Advantage view "
            "are available in MMD but not ingested.",
            "2023 claims are preliminary (_p files) and may change when CMS publishes final 2023.",
            "The data API is the MMD app's backend, not a documented public API; file names can change on refresh.",
        ],
    }


def write_audit(stats: dict, cov: pd.DataFrame, retrieved: str) -> Path:
    man = load_manifest(SOURCE_ID)["files"]
    api_files = {k: v for k, v in man.items() if k.startswith("api/")}
    ref_files = {k: v for k, v in man.items() if not k.startswith("api/")}
    total_api_bytes = sum(v["bytes"] for v in api_files.values())
    gy = stats["geo_match_by_year"]
    ps = stats["prevalence_summary"]
    lines = []
    A = lines.append
    A("# DATA AUDIT — CMS Mapping Medicare Disparities (MMD) Tool, Population View\n")
    A("| Field | Value |\n|---|---|")
    A(f"| source_id | {SOURCE_ID} |")
    A(f"| Source (dataset/API name, exact files/endpoints) | MMD Population View data API `{API_BASE}` "
      f"(the backend of {APP_URL}); {len(api_files)} JSON slices, file names resolved through the app's "
      "`codebook_crosswalk.csv`; see table below |")
    A("| Publishing organization | CMS Office of Minority Health (tool built by NORC for CMS OMH); condition algorithms "
      "from the CMS Chronic Conditions Data Warehouse (CCW) |")
    A(f"| Retrieval date (UTC) | {retrieved} (first API slice); see MANIFEST.json for each file |")
    A("| Source version / release | Medicare FFS claims 2012-2022 final, 2023 preliminary (`_p` file suffix); "
      "technical documentation revision dated 6/26/2024 (covers through 2022; the app now also serves 2023) |")
    A("| Source update date / cadence | Periodic: a new claims year roughly annually, preliminary then final "
      "(tech-doc revision history) |")
    A("| License / access conditions | U.S. federal government work, public domain; open access, no registration "
      "or DUA. MMD is not a catalogued data.cms.gov dataset (absent from https://data.cms.gov/data.json at "
      "retrieval); slices were requested exactly as the public app requests them |")
    A("| Unit of observation | Aggregate cell: geography x claims year x MMD condition x measure x adjustment x age "
      "stratum, for Medicare fee-for-service beneficiaries |")
    A(f"| Sample size (actual, as ingested) | {stats['rows']['parsed_cells']:,} parsed cells -> "
      f"{stats['rows']['burden']:,} burden rows, {stats['rows']['context']:,} context rows, "
      f"{stats['rows']['coverage']} coverage rows. Underlying FFS analysis population (tech doc Table A.1, national): "
      "34,066,594 (2012) ... 29,614,997 (2022) beneficiaries; per-cell denominators published only as categories |")
    A("| Geography (resolution, vintage) | County (CMS FIPS vintage, legacy Connecticut counties), state/territory, "
      "nation. Matched against `geographies` (2024 counties, CT planning regions; legacy CT flagged) |")
    A("| Person-level? | no |\n| Geographic? | yes |\n| Omics? | no |\n| Wearable? | no |")
    A("| True participant linkage across modalities? | no — aggregate claims statistics with no person identifiers |\n")

    A("## Why this source and route\n")
    A("* The MMD Population View is the CMS source named in SPEC §E. It offers no bulk download: its \"Download "
      "Data\" button exports the current map selection only. The app loads each selection from "
      f"`{API_BASE}?_source=<file>&year=..&geography=..&measure=..&condition=..&...`, choosing `<file>` from "
      "`/mmd-population/assets/codebook_crosswalk.csv`. `measure_it.ingestion.cms_mmd` replicates that lookup "
      "(`resolve_source_file`) and the filter syntax (`.|IS NULL` for 'all'; `a|b` alternatives), including the "
      "app's one-time 2019 ED-visit `elig` exception.")
    cat = stats["catalog"]
    A(f"* Alternatives checked at retrieval: the data.cms.gov catalog (`catalog/data.cms.gov_data.json`, "
      f"{cat['n_datasets']} datasets) has {len(cat['mmd_or_prevalence_hits'])} dataset(s) whose title/description "
      f"matches disparit/mapping medicare/chronic condition/fibromyalg/prevalence"
      f"{(': ' + '; '.join(cat['mmd_or_prevalence_hits'])) if cat['mmd_or_prevalence_hits'] else ''} (so no "
      "'Specific Chronic Conditions'-type dataset either). Medicare Geographic Variation ("
      f"{'; '.join(cat['geographic_variation_titles']) or 'not listed'}) covers utilization/spending, not these "
      "CCW 'other chronic' conditions. MMD is therefore the official CMS county-level source used, and the only one "
      "found that exposes CCW 'Fibromyalgia, Chronic Pain and Fatigue' and 'Migraine and Other Chronic Headache'.\n")

    A("## Files / endpoints retrieved\n")
    split = [r for r in stats["requests"] if len(r.get("fltrs") or []) == 1]
    split_desc = "; ".join(f"measure {m} years {min(ys)}-{max(ys)}" for m, ys in sorted(
        {m: [r["year"] for r in split if r["measure"] == m] for m in {r["measure"] for r in split}}.items()))
    missing = [r for r in stats["requests"] if r.get("fltrs_missing")]
    A(f"* {len(api_files)} API slices under `api/`, {total_api_bytes/1e6:.1f} MB total; each MANIFEST entry has the full "
      "query URL, bytes, sha256 and retrieved_at. Requests: one per (claims year x measure x source file), geography "
      "`c|s|n`, all sex/race/dual/eligibility. The file is resolved separately for each adjustment, as the app does: "
      f"where the crosswalk splits a measure into one file per adjustment ({len(split)} requests: "
      f"{split_desc or 'none'}) each adjustment is its own request; otherwise one request asks for `fltr=1|2`. Requested "
      f"adjustments absent from the response: {len(missing)} request(s)"
      + (": " + "; ".join(f"{r['year']} {r['measure']} fltr {r['fltrs_missing']}" for r in missing) if missing else "")
      + ".")
    A("\n| year | measure | fltr | MMD source file | status | rows |\n|---|---|---|---|---|---|")
    for r in stats["requests"]:
        A(f"| {r['year']} | {r['measure']} | {'&#124;'.join(r.get('fltrs') or [])} | `{r['source_file']}` | "
          f"{r['status']} | {r.get('n_rows') or 0:,} |")
    A("\nReference files (codes, definitions):\n")
    A("| file | url | bytes | sha256 (first 16) |\n|---|---|---|---|")
    for k, v in sorted(ref_files.items()):
        A(f"| `{k}` | {v['url']} | {v['bytes']:,} | {v['sha256'][:16]} |")

    A("\n## Code books (decoded from app/menus.js, app/strings.json, app/datadownload.js)\n")
    A("* year: 2..9 = 2012..2019, 20..23 = 2020..2023. geography: c county, s state/territory, n nation.")
    A("* measure: v prevalence (%), h hospitalization (per 1,000), e ED visits (per 1,000), t average total cost, "
      "a risk-adjusted total cost, p average principal cost ($ per beneficiary with the condition).")
    A("* fltr (adjustment): 1 unsmoothed actual, 2 unsmoothed age-standardized (ingested); 3/4 smoothed (not ingested).")
    A("* agecat: all, 0 <65, 1 65-74, 2 75-84, 3 85+, 4 65+. dencat: 1 11-499, 2 500-999, 3 1,000-4,999, "
      "4 5,000-9,999, 5 10,000+.")
    A("* MMD conditions used: 51 Fibromyalgia, Chronic Pain and Fatigue; 59 Migraine and Other Chronic Headache; "
      "24 Disability (reason for Medicare eligibility); 134 COVID-19; 10 all-cause (hospitalization/ED). Labels are "
      "re-verified against the downloaded menus.js on every run.\n")

    A("## Target conditions: what CMS/MMD covers and burden_evidence_level\n")
    cs = stats.get("ccw_search", {})
    pages = cs.get("algorithms_parsed_per_page", {})
    A("Target ICD-10-CM codes were searched (prefix match) in: "
      + "; ".join(f"`{k}`: {v} algorithms parsed from its code table" for k, v in pages.items())
      + " (the chronic page has no code table: it only links PDFs); and "
      + (f"the CCW 30 chronic-condition algorithms PDF (`ccw/chr-chronic-condition-algorithms.pdf`, text-extracted with "
         f"pdftotext, {cs.get('chronic30_pdf_distinct_icd10_codes')} distinct ICD-10 codes; hits are document-level, "
         "not attributed to one algorithm)." if cs.get("chronic30_pdf_searched") else
         "NOT the CCW 30 chronic-condition algorithms (PDF only; pdftotext unavailable), so level D is verified only "
         "against the 'other chronic' algorithms.")
      + " Search keys = module prefixes merged with the ICD-10-CM codes of the sourced `condition_registry` table.\n")
    A("| condition_id | level | MMD / CCW measure | searched codes | hits in MMD-exposed 'other chronic' CCW algorithms "
      "| hits in CCW 30 chronic-condition PDF |\n|---|---|---|---|---|---|")
    def _txt(v, empty="—"):
        return empty if v is None or (isinstance(v, float) and np.isnan(v)) or v == "" else v
    for r in cov.itertuples():
        A(f"| {r.condition_id} | **{r.burden_evidence_level}** | {_txt(r.ccw_algorithm_name)} | "
          f"{r.icd10_codes_searched} | {_txt(r.search_hits_in_mmd_exposed_ccw_algorithms, 'none')} | "
          f"{r.search_hits_in_ccw_30_chronic_conditions} |")
    fib = cov[cov["condition_id"] == "fibromyalgia"].iloc[0]
    mig = cov[cov["condition_id"] == "migraine"].iloc[0]
    A(f"\n* CCW **{fib.ccw_algorithm_name}** ({fib.ccw_reference_period}; {fib.ccw_qualifying_claims}): "
      f"{fib.ccw_n_icd10_codes} ICD-10 codes — {fib.ccw_icd10_codes}.")
    A(f"* CCW **{mig.ccw_algorithm_name}** ({mig.ccw_reference_period}; {mig.ccw_qualifying_claims}): "
      f"{mig.ccw_n_icd10_codes} ICD-10 codes (G43.x migraine + G44.x other headache syndromes).")
    A("* Rationale per level is stored in `condition_burden_coverage__cms_mmd.rationale` and in every burden row's "
      "`provenance_notes`. No condition is level A: MMD has no direct measure for any target condition.")
    A("* CCW also publishes a newer 'Chronic Pain' algorithm (2024) that contains M79.7 and G43/G44; MMD does not "
      "expose it (see ccw/comparing-chronic-pain.pdf).\n")

    A("## Key variables\n")
    A("`geo_condition_burden__cms_mmd`: `geo_id` (5-char county / 2-char state / `US`), `geo_level`, "
      "`geo_match_status` (matched_2024 / matched_ct_legacy / unmatched / national), `fips_as_published` (the code "
      "CMS published; differs from `geo_id` only for the two bridged 1:1 renames), `fips_bridge`, `year` (claims year), "
      "`condition_id` (target), `burden_evidence_level`, `source_condition_code/name`, `ccw_algorithm_name`, "
      "`prevalence_pct` (integer % of FFS beneficiaries, as published), `adjustment`, `stratum_age` (all or <65), "
      "`default_view` (unsmoothed actual, all ages = the MMD default map), `denominator_category/range/min/max`, "
      "`value_zero_possibly_suppressed`, `claims_data_status`. Use age-standardized rows to compare counties with "
      "different age structures; use actual rows to describe the coded-condition share.\n")
    A("`geo_context__cms_mmd`: same keys plus `context_category` (utilization / spending / population_composition "
      "/ infection_exposure), `measure`, `value`, `unit`, `related_condition_ids`. Never a burden estimate.\n")

    A("## Missingness (measured)\n")
    A(f"Reference: {stats['n_counties_2024']} non-legacy county rows in `geographies` (2024 vintage incl. PR and "
      "island areas).\n")
    A("| condition | year | national % | state range % | county rows | county median % (p5-p95) | county zeros "
      "(possibly suppressed) | 2024 counties without a value |\n|---|---|---|---|---|---|---|---|")
    for r in ps:
        A(f"| {r['condition_id']} | {r['year']} | {r['national_pct']} | {r['state_min']:.0f}-{r['state_max']:.0f} | "
          f"{r['county_n']} | {r['county_median']:.0f} ({r['county_p5']:.0f}-{r['county_p95']:.0f}) | "
          f"{r['county_zero']} | {r['counties_2024_without_value']} |")
    A("\nCounties never have a missing value inside a published row; missingness is whole cells omitted by CMS "
      "(<11 beneficiaries) or FIPS that do not exist in the 2024 vintage. <65 stratum county rows (fibromyalgia "
      "mapping) by year: " + ", ".join(f"{y}: {n}" for y, n in stats["lt65_county_rows_by_year"].items()) + ".")
    A("\nDenominator categories, county rows (fibromyalgia/pain/fatigue, actual, all ages):\n")
    for y, d in stats["denominator_category_counts_county_fibro"].items():
        A(f"* {y}: " + ", ".join(f"{k}: {v}" for k, v in d.items()))
    A("\nContext measures (unsmoothed actual, all ages):\n")
    A("| measure | MMD condition | years | county rows | county zeros | national (latest year) |\n|---|---|---|---|---|---|")
    for r in stats["context_summary"]:
        yrs = f"{r['years'][0]}-{r['years'][-1]} ({len(r['years'])})" if r["years"] else "—"
        A(f"| {r['measure']} | {r['condition']} {r['name']} | {yrs} | {r['county_rows']:,} | {r['county_zero']:,} | "
          f"{r['national_latest']} ({r['national_latest_year']}) |")

    A("\n## Geography match against `data/processed/geographies` (prevalence, condition 51, actual, all ages)\n")
    A("| year | county rows | matched 2024 | of which bridged 1:1 renames (published -> 2024) | matched legacy CT | "
      "CT planning-region rows | unmatched (all) | of which CMS xx990 codes | % rows matched 2024 | % of 2024 counties "
      "covered | other unmatched FIPS | states unmatched |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in gy:
        A(f"| {r['year']} | {r['county_rows']} | {r['matched_2024']} | {', '.join(r.get('bridged_renames') or []) or '—'} "
          f"| {r['matched_ct_legacy']} | "
          f"{r['ct_planning_region_rows']} | {r['unmatched']} | {r['unmatched_cms_990']} | "
          f"{r['pct_rows_matched_2024']} | {r['pct_2024_counties_covered']} | "
          f"{', '.join(r['unmatched_geo_ids']) or '—'} | {', '.join(r['state_unmatched']) or '—'} |")
    A("\n* **1:1 FIPS renames (bridged at ingestion).** CMS publishes two counties under codes retired in 2015 in every "
      "year: " + "; ".join(f"{k} -> {v[0]} ({v[1]})" for k, v in FIPS_RENAMES.items()) + ". Both are renames with an "
      "unchanged boundary (Census 'Substantial Changes to Counties and County Equivalent Entities'), so the published "
      "value belongs to the 2024 county unchanged: the rows carry the 2024 code in `geo_id` / `county_fips` "
      "(`geo_match_status='matched_2024'`), the code CMS published in `fips_as_published` and in `source_record_id`, "
      "`fips_bridge` names the rename, and `provenance_notes` starts with the bridge sentence. The CDC Lyme and HRSA "
      "ingestions apply the same two renames. Boundary changes are not bridged: 51515 (Bedford city VA, merged into "
      "Bedford County 51019 in 2013) stays unmatched.")
    c990 = stats["cms_990"]
    A("\n* `xx990` rows: CMS publishes county codes ending in 990 in several states (rows per year: "
      + ", ".join(f"{y}: {v['n']}" for y, v in c990.items()) + "; denominator classes across all years: "
      + ", ".join(f"{k}: {sum(v['denominator_ranges'].get(k, 0) for v in c990.values())}"
                  for k in sorted({k for v in c990.values() for k in v['denominator_ranges']}))
      + "). They are not Census county codes and are not described in the MMD technical documentation; rows are "
      "kept with `unmatched_reason` set and must not be mapped to any county.")
    ak = stats["alaska"]
    A(f"* Alaska: CMS publishes {'/'.join(str(v) for v in sorted(set(ak['codes_by_year'].values())))} Alaska "
      "county codes per year "
      f"(identical set every year: {ak['codes_identical_all_years']}); `geographies` has {ak['n_2024_alaska']} "
      "Alaska county-equivalents in 2024, so newer boroughs/census areas have no MMD value.")
    cw = stats["counties_2024_without_value_latest"]
    A(f"* 2024 county-equivalents without an MMD value in {cw['year']}: {cw['n']} — "
      + "; ".join(f"{k} ({len(v)}): {', '.join(v)}" for k, v in cw["by_state"].items()) + ".")
    legacy_years = [r["year"] for r in gy if r["matched_ct_legacy"] > 0]
    region_years = [r["year"] for r in gy if r["ct_planning_region_rows"] > 0]
    A(f"\nConnecticut: legacy CT county FIPS (09001-09015) appear in {len(legacy_years)} of {len(gy)} years "
      f"({', '.join(map(str, legacy_years)) or 'none'}); CT planning-region FIPS (09110-09190) appear in "
      f"{len(region_years)} years ({', '.join(map(str, region_years)) or 'none'}). Legacy rows match the flagged "
      "`ct_legacy=True` rows in `geographies` (`geo_match_status='matched_ct_legacy'`) and are NOT re-assigned to "
      "planning regions; planning regions only carry a value where CMS itself published one. Other unmatched FIPS "
      "are county codes that do not exist in the 2024 vintage (listed above); they are kept and reported, not "
      "forced.\n")

    A("## Linkage strategy\n")
    A("* `geo_id` / `county_fips` (CMS vintage, with the two 1:1 FIPS renames bridged to their 2024 codes; "
      "`fips_as_published` keeps the CMS code) -> `geographies.geo_id` (2024 counties + legacy CT). State rows "
      "-> state FIPS. All joins are ecological (population aggregates); no person-level join is possible or allowed.")
    A("* State rows keep `source_geographic_resolution='state'` and are never attached to counties as county values.")
    A("* `condition_id` -> `configs/conditions.yaml` ids; the mapping level is `burden_evidence_level`.\n")

    A("## Limitations and caveats\n")
    dis = next((r for r in stats["context_summary"] if r["measure"] == "prevalence" and r["condition"] == "24"), None)
    fib_nat = [r for r in ps if r["condition_id"] == "fibromyalgia" and r["national_pct"] is not None]
    fib_last = fib_nat[-1] if fib_nat else None
    A(f"* **{LIMITATION}** Measured: disability was the original or current reason for entitlement for "
      f"{dis['national_latest'] if dis else 'UNKNOWN'} % of FFS beneficiaries nationally in "
      f"{dis['national_latest_year'] if dis else 'UNKNOWN'} (context table, condition 24); the <65 stratum is "
      "provided for the burden measures.")
    A("* The fibromyalgia mapping is a broad composite (chronic pain, radiculopathy, myositis, myalgia, neuralgia, "
      f"chronic fatigue); its national prevalence was {fib_last['national_pct'] if fib_last else 'UNKNOWN'} % in "
      f"{fib_last['year'] if fib_last else 'UNKNOWN'} (and {fib_nat[0]['national_pct'] if fib_nat else 'UNKNOWN'} % "
      f"in {fib_nat[0]['year'] if fib_nat else 'UNKNOWN'}), far above fibromyalgia alone. For ME/CFS it is only a "
      "symptom proxy (level C): it lacks G93.3x.")
    A("* Migraine mapping includes all G44 headache syndromes (tension-type, cluster, post-traumatic, drug-induced).")
    fn = {r["year"]: r["national_pct"] for r in ps if r["condition_id"] == "fibromyalgia"}
    mn = {r["year"]: r["national_pct"] for r in ps if r["condition_id"] == "migraine"}
    A("* Coding-regime break: national fibromyalgia/pain/fatigue prevalence by year " +
      ", ".join(f"{y}: {v:g}" for y, v in sorted(fn.items())) + "; migraine/headache " +
      ", ".join(f"{y}: {v:g}" for y, v in sorted(mn.items())) + ". The largest step (2015 -> 2016) coincides with "
      "the ICD-9-CM -> ICD-10-CM switch (1 Oct 2015) and the 2-year lookback; treat 2012-2015 and 2016+ as different "
      "coding regimes and do not read the step as a change in illness.")
    A("* The tech doc's 27 -> 30 CCW indicator update (applied from 2021) concerns the primary chronic-condition "
      "indicators; conditions 51 and 59 come from CCW's separate 'Other Chronic or Potentially Disabling Conditions' "
      "algorithms, whose code lists also evolve with ICD-10-CM. Diagnosis/coding intensity varies by place and time, "
      "and the FFS population shrinks as enrolment shifts to Medicare Advantage (FFS analysis population 34,066,594 "
      "in 2012 vs 29,614,997 in 2022, tech doc Table A.1), so trends are not illness trends.")
    A("* MMD 'COVID-19' (condition 134, 2020+) counts any COVID-19 diagnosis on a claim; it is acute-infection "
      "context in `geo_context__cms_mmd` (context_category='infection_exposure') and is never used as Long COVID burden.")
    A(f"* {SUPPRESSION_NOTE}")
    A("* Prevalence is published as an integer percent; small differences between counties are below the "
      "published resolution. Denominators are published only as five size categories.")
    for r in stats["context_summary"]:
        if r["measure"] == "hospitalization_rate" and r["condition"] in ("51", "59") and r["county_rows"]:
            A(f"* Hospitalization rate, MMD {r['condition']} ({r['name']}): {r['county_zero']:,} of "
              f"{r['county_rows']:,} county rows ({100 * r['county_zero'] / r['county_rows']:.1f} %) are 0 per 1,000 "
              f"(principal-diagnosis discharges, integer-rounded); national latest {r['national_latest']}.")
    A("* 2023 files are preliminary (`_p`) claims. The data API is the MMD app backend, not a documented public API; "
      "file names can change when CMS refreshes the tool (the module re-resolves them from the crosswalk).")
    A("* Not ingested although available in MMD: smoothed rates, sex/race/dual/eligibility strata, Medicare "
      "Advantage view, other conditions, mortality/readmission/PQI/preventive-service measures.\n")

    A("## Processed outputs\n")
    A(f"* `data/processed/geo_condition_burden__cms_mmd.parquet` — {stats['rows']['burden']:,} rows "
      f"({', '.join(f'{k}: {v:,}' for k, v in stats['burden_rows_by_condition'].items())})")
    A(f"* `data/processed/geo_context__cms_mmd.parquet` — {stats['rows']['context']:,} rows "
      f"({', '.join(f'{k}: {v:,}' for k, v in stats['context_rows_by_measure'].items())})")
    A(f"* `data/processed/condition_burden_coverage__cms_mmd.parquet` — {stats['rows']['coverage']} rows "
      "(one per target condition in configs/conditions.yaml)\n")
    A("## Reproduce\n")
    A("`uv run python -m measure_it.ingestion.cms_mmd` (idempotent; raw slices are re-used from MANIFEST.json). "
      "Tests: `uv run pytest tests/test_cms_mmd.py`.")
    path = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    path.write_text("\n".join(lines) + "\n")
    return path


# =========================================================================================
# Entry point
# =========================================================================================
def run() -> dict:
    ref = fetch_reference_files()
    crosswalk = load_crosswalk(ref["app/codebook_crosswalk.csv"].read_text())
    strings = json.loads(ref["app/strings.json"].read_text())
    labels = verify_condition_labels(ref["app/menus.js"].read_text(), strings)
    ccw_pages = {rel: parse_ccw_algorithms(ref[rel].read_text(errors="replace"))
                 for rel in ("ccw/condition-categories-other.html", "ccw/condition-categories-chronic.html")}
    ccw = pd.concat(list(ccw_pages.values()), ignore_index=True)
    for code, name in CCW_ALGORITHM_FOR_CODE.items():
        if name not in set(ccw["ccw_algorithm_name"]):
            raise RuntimeError(f"CCW algorithm {name!r} (MMD {code}) not found on the CCW pages")
    # the chronic page only links PDFs; its 30 algorithms are searched in the combined PDF's text
    chr_pdf = ref.get("ccw/chr-chronic-condition-algorithms.pdf")
    chr_text = pdf_text(chr_pdf) if chr_pdf is not None else None
    chronic30_codes = icd10_codes_in_text(chr_text) if chr_text else None
    if chronic30_codes is not None and len(chronic30_codes) < 100:
        raise RuntimeError(f"only {len(chronic30_codes)} ICD-10 codes extracted from the CCW 30-condition PDF")
    try:
        reg = read_table("condition_registry", columns=["canonical_condition_id", "icd10cm_codes"])
        registry_codes = {r.canonical_condition_id: [str(x) for x in r.icd10cm_codes]
                          for r in reg.itertuples() if r.icd10cm_codes is not None}
    except FileNotFoundError:
        registry_codes = {}

    plan = fetch_api_slices(plan_requests(crosswalk))
    df = parse_slices(plan)
    geos = read_table("geographies")
    df = attach_geo_match(df, geos)

    burden = build_burden(df)
    context = build_context(df)
    js = ref["app/menus.js"].read_text()
    mmd_names = {strings[sid]["english"].strip() for sid in re.findall(r"get_string_by_id\(\"(\w+)\"\)", js)
                 if sid in strings}
    cov = build_coverage(load_config("conditions")["conditions"], ccw, mmd_names, burden,
                         chronic30_codes=chronic30_codes, registry_codes=registry_codes)

    write_table(burden, "geo_condition_burden__cms_mmd", producer=PRODUCER,
                description="CMS MMD (Medicare FFS) coded-condition prevalence mapped to target conditions with "
                            "burden_evidence_level; county/state/national x year")
    write_table(context, "geo_context__cms_mmd", producer=PRODUCER,
                description="CMS MMD (Medicare FFS) utilization, spending and population-composition context")
    write_table(cov, "condition_burden_coverage__cms_mmd", producer=PRODUCER,
                description="Per target condition: CMS/CCW measure coverage and burden_evidence_level A/B/C/D")

    stats = compute_stats(df, burden, context, cov, plan, geos)
    stats["condition_labels_verified"] = labels
    stats["ccw_search"] = {
        "algorithms_parsed_per_page": {k: int(len(v)) for k, v in ccw_pages.items()},
        "chronic30_pdf_searched": chronic30_codes is not None,
        "chronic30_pdf_distinct_icd10_codes": len(chronic30_codes) if chronic30_codes is not None else None,
        "registry_codes_merged": {k: v for k, v in registry_codes.items()},
    }
    stats["catalog"] = catalog_check(ref["catalog/data.cms.gov_data.json"].read_text())
    retrieved = min(v["retrieved_at"] for k, v in load_manifest(SOURCE_ID)["files"].items() if k.startswith("api/"))
    (raw_dir(SOURCE_ID) / "audit_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    write_registry_entry(registry_entry(stats, retrieved))
    write_audit(stats, cov, retrieved)
    print(f"[cms_mmd] burden={len(burden):,} context={len(context):,} coverage={len(cov)}")
    return stats


if __name__ == "__main__":
    run()
