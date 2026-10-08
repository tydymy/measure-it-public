"""openFDA device endpoints: which measurement classes have FDA-regulated device categories,
and which specific devices have 510(k) / De Novo / PMA records.

Pipeline
--------
1. Snapshot the full device-classification file (bulk download, MANIFEST.json) so every
   product code can be checked against FDA's own list.
2. DISCOVERY. For every measurement class in configs/measurements.yaml, query
   https://api.fda.gov/device/classification.json with keywords derived from that class,
   searching `device_name` OR `definition` (DISCOVERY_TERMS). Some devices sit under
   generically named product codes (e.g. the Zio patch under "Recorder, Magnetic Tape, Medical"),
   so a second recorded discovery path searches 510k.json by device/applicant name
   (NAMED_DEVICE_SEARCHES). Every call is written to data/raw/openfda_device/QUERY_LOG.json.
3. CURATION (human, not automated): configs/fda_product_code_map.yaml maps
   measurement_id -> product codes with a rationale and the exact recorded query that found
   each code. `validate_code_map` refuses codes that are not in FDA's classification file or
   that the named query did not return.
4. For each mapped product code: 510k.json (De Novo decisions live in the same endpoint with
   DEN numbers / decision_code DENG), pma.json, and registrationlisting.json counts.
   A code whose FDA name is generic (e.g. DQK "Computer, Diagnostic, Programmable") is only
   mapped with a `device_filter`, and all its counts are restricted to that filter.

Guardrail: a regulatory record is a deployment-readiness signal, not evidence that the
technology diagnoses the target illness (CONVENTIONS §3.7). Nothing here says a device is
cleared unless its K/DEN/P number came back from the API.

Outputs (data/processed): fda_device_classification, fda_510k, fda_pma,
fda_registration_listing_counts, measurement_regulatory_status, fda_named_device_findings,
fda_device_query_log.

Run: uv run python -m measure_it.ingestion.openfda [--discover-only]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

from ..config import CONFIGS, UNKNOWN, load_config, raw_dir, utc_now_iso
from ..download import download_file, load_manifest, sha256_file
from ..http import request
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table

SOURCE_ID = "openfda_device"
SOURCE_NAME = "openFDA device API (classification, 510k incl. De Novo, PMA, registration & listing)"
PRODUCER = "measure_it.ingestion.openfda"
API = "https://api.fda.gov/device"
CLASSIFICATION_BULK_URL = (
    "https://download.open.fda.gov/device/classification/device-classification-0001-of-0001.json.zip"
)
CODE_MAP_PATH = CONFIGS / "fda_product_code_map.yaml"
REGULATORY_NOTE = (
    "A regulatory record is a deployment-readiness signal, not evidence that the technology "
    "diagnoses the target illness"
)
NO_CODE = "no specific product code found"

# openFDA returns at most 1000 records per call. Up to that many we keep every record for a
# product code; above it we keep the 100 most recent decisions (sort=decision_date:desc) and
# take the counts from separate count queries. retrieval_mode records which happened.
API_MAX_LIMIT = 1000
RECENT_N = 100

# openFDA field doc (https://open.fda.gov/fields/deviceclass.yaml) documents 1-4 only.
SUBMISSION_TYPES = {"1": "510(k)", "2": "PMA", "3": "Contact ODE", "4": "510(k) exempt"}

# ---------------------------------------------------------------------------------------------
# Discovery vocabulary: keywords derived from each class's name, signals and text patterns in
# configs/measurements.yaml. Each term is searched in device_name OR definition. A trailing *
# is an openFDA prefix wildcard; multi-word terms are phrase searches.
# ---------------------------------------------------------------------------------------------
DISCOVERY_TERMS: dict[str, list[str]] = {
    "accelerometry": ["actigraph*", "acceleromet*", "sleep assessment", "activity monitor", "physical activity"],
    "wearable_heart_rate": ["heart rate", "cardiotachometer", "pulse rate", "wearable"],
    "hrv": ["heart rate variability", "variability", "autonomic"],
    "ecg_ambulatory": ["electrocardiograph*", "electrocardiogram*", "ambulatory", "holter", "arrhythmia",
                       "cardiac telemetry", "event recorder"],
    "ppg": ["photoplethysmograph*", "plethysmograph*", "pulse wave"],
    "posture_detection": ["posture", "body position", "upright", "inclinometer"],
    "continuous_spo2": ["oximeter", "oximetry", "oxygen saturation"],
    "respiratory_rate": ["breathing frequency", "respiratory rate", "breathing rate", "respiration rate",
                         "apnea monitor"],
    "sleep_objective": ["polysomnograph*", "sleep", "ventilatory effort"],
    "continuous_temperature": ["thermometer", "body temperature", "skin temperature", "temperature measurement"],
    "cpet": ["oxygen-uptake", "oxygen uptake", "spirometer", "pulmonary function", "exerciser", "ergometer",
             "cardiopulmonary exercise", "metabolic"],
    "autonomic_testing": ["tilt", "galvanic skin response", "sudomotor", "sweat", "valsalva", "autonomic",
                          "baroreflex"],
    "small_fiber_testing": ["nerve fiber", "temperature discrimination", "vibration threshold",
                            "evoked potential stimulator", "esthesiometer", "nerve conduction", "corneal",
                            "skin biopsy"],
    "endothelial_function": ["endothelial", "tonometry", "tonometer", "reactive hyperemia", "flow-mediated",
                             "plethysmograph*"],
    "vascular_imaging": ["pulsed doppler", "pulsed echo", "nuclear magnetic resonance", "echocardiograph*",
                         "computed emission", "vascular ultrasound"],
    "capillaroscopy": ["capillaroscop*", "capillary", "nailfold", "microscope", "dermatoscope",
                       "microcirculation"],
    "microvascular_function": ["laser doppler", "blood flow", "blood-flow", "flowmeter", "tissue saturation",
                               "near infrared", "microvascular"],
    "retinal_imaging": ["optical coherence", "ophthalmoscope", "camera, ophthalmic", "fundus", "retinal",
                        "retinopathy"],
    "digital_gait": ["gait", "goniometer", "kinematic", "inertial", "walking"],
    "digital_cognitive": ["cognitive", "cognitive assessment", "attention task", "neuropsychological"],
    "blood_biomarkers": ["c-reactive protein", "fibrin", "d-dimer", "interleukin", "cytokine", "cortisol",
                         "tryptase"],
    "metabolomics": ["metabolom*", "mass spectromet*", "metabolite"],
    "proteomics": ["proteom*", "aptamer", "multiplex protein", "protein panel"],
    "transcriptomics": ["gene expression", "transcriptom*", "mrna", "expression profiling"],
    "immune_assays": ["autoantibod*", "flow cytometr*", "lymphocyte", "immunophenotyp*", "t cell"],
    "remote_patient_monitoring": ["remote", "telemetry", "physiological signal", "physiological, patient",
                                  "transmitters and receivers"],
}

# Named-device discovery on 510k.json (device_name / applicant). `target` is what SPEC §D and
# the build brief asked us to look for specifically.
NAMED_DEVICE_SEARCHES: list[dict] = [
    # consumer wearable ECG / irregular-rhythm features
    *[{"target": "consumer_wearable_ecg_irregular_rhythm", "measurement_ids": ["ecg_ambulatory", "ppg", "wearable_heart_rate"],
       "search": s} for s in [
        'applicant:"apple inc"', "applicant:fitbit", "applicant:google", 'applicant:"samsung electronics"',
        "applicant:withings", "applicant:garmin", "applicant:whoop", "applicant:oura", "applicant:alivecor",
        'device_name:"irregular rhythm"']],
    # wearable patch / ambulatory ECG
    *[{"target": "wearable_patch_ecg", "measurement_ids": ["ecg_ambulatory"], "search": s} for s in [
        "applicant:irhythm", "device_name:zio", "applicant:bardy", "applicant:preventice", "device_name:holter"]],
    # peripheral arterial tonometry / endothelial function
    *[{"target": "peripheral_arterial_tonometry", "measurement_ids": ["endothelial_function"], "search": s} for s in [
        # not endopat*: that prefix also matches Ethicon "Endopath" laparoscopic instruments
        "applicant:itamar", "device_name:endopat", "device_name:endopatx", 'device_name:"endo pat"',
        'device_name:"peripheral arterial"', "device_name:vendys",
        'device_name:"endothelial function"']],
    # sudomotor testing
    *[{"target": "sudomotor_testing", "measurement_ids": ["autonomic_testing"], "search": s} for s in [
        "device_name:sudoscan", "applicant:impeto", "device_name:sudomotor", 'device_name:"quantitative sweat"',
        "device_name:qsart"]],
    # autonomic laboratory equipment
    *[{"target": "tilt_table_autonomic_lab", "measurement_ids": ["autonomic_testing"], "search": s} for s in [
        'device_name:"tilt table"', 'device_name:autonomic']],
    *[{"target": "heart_rate_variability_software", "measurement_ids": ["hrv"], "search": s} for s in [
        'device_name:"heart rate variability"', "device_name:hrv"]],
    # actigraphy
    *[{"target": "actigraphy", "measurement_ids": ["accelerometry", "sleep_objective"], "search": s} for s in [
        "device_name:actigraph*", "device_name:actiwatch*", 'device_name:"activity monitor"',
        "device_name:geneactiv", "device_name:motionwatch", "applicant:actigraph"]],
    # nailfold capillaroscopy / microcirculation imaging
    *[{"target": "nailfold_capillaroscopy", "measurement_ids": ["capillaroscopy"], "search": s} for s in [
        "device_name:capillaroscop*", "device_name:videocapillaroscop*", "device_name:nailfold",
        'device_name:"capillary microscope"', "device_name:cytocam", 'device_name:"dark field"']],
    # microvascular flow
    *[{"target": "laser_doppler_speckle", "measurement_ids": ["microvascular_function"], "search": s} for s in [
        'device_name:"laser doppler"', "device_name:speckle", "applicant:perimed"]],
    # corneal confocal microscopy (small-fiber)
    *[{"target": "corneal_confocal_microscopy", "measurement_ids": ["small_fiber_testing"], "search": s} for s in [
        "device_name:rostock", 'device_name:"corneal confocal"']],
    # quantitative sensory testing (small-fiber)
    *[{"target": "quantitative_sensory_testing", "measurement_ids": ["small_fiber_testing"], "search": s} for s in [
        "applicant:medoc", 'device_name:"quantitative sensory"']],
    # digital gait analysis
    *[{"target": "digital_gait_analysis", "measurement_ids": ["digital_gait"], "search": s} for s in [
        "device_name:gait", 'device_name:"gait analysis"', 'device_name:"inertial sensor"']],
    # CPET
    *[{"target": "cardiopulmonary_exercise_testing", "measurement_ids": ["cpet"], "search": s} for s in [
        'device_name:"cardiopulmonary exercise"', 'device_name:"metabolic cart"', 'device_name:"metabolic measurement"']],
    # D-dimer assays (blood biomarkers)
    {"target": "d_dimer_assay", "measurement_ids": ["blood_biomarkers"], "search": 'device_name:"d-dimer"'},
    # wearable vital-sign patches (remote monitoring)
    *[{"target": "wearable_vital_sign_patch", "measurement_ids": ["remote_patient_monitoring", "wearable_heart_rate"],
       "search": s} for s in ["applicant:vitalconnect", "applicant:biointellisense"]],
]


# ---------------------------------------------------------------------------------------------
# Pure helpers (unit-tested)
# ---------------------------------------------------------------------------------------------
def build_term(field: str, term: str) -> str:
    """openFDA field clause: prefix wildcard stays bare, anything with spaces/punctuation is a phrase."""
    term = term.strip()
    if term.endswith("*"):
        words = term.split()
        if len(words) == 1:
            return f"{field}:{term}"
        # multi-word prefix search: every word must match (a bare space would mean OR)
        return "(" + " AND ".join(f"{field}:{w}" for w in words) + ")"
    if re.search(r"[\s,\-/()]", term):
        return f'{field}:"{term}"'
    return f"{field}:{term}"


def classification_search(term: str, fields: tuple[str, ...] = ("device_name", "definition")) -> str:
    """device_name:<term> OR definition:<term> (a space between clauses is OR in openFDA)."""
    return " ".join(build_term(f, term) for f in fields)


def submission_kind(number: str) -> str:
    """Classify a 510k-endpoint record number: DEN... = De Novo, K... = 510(k)."""
    n = (number or "").upper()
    if n.startswith("DEN"):
        return "de_novo"
    if n.startswith("K"):
        return "510k"
    return "other"


def query_id(endpoint: str, params: dict) -> str:
    blob = json.dumps({"e": endpoint, "p": sorted(params.items())}, sort_keys=True)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]


def product_code_search(code: str, device_filter: str | None = None) -> str:
    base = f"product_code:{code}"
    return f"{base} AND {device_filter}" if device_filter else base


def example_devices(records: list[dict], n_recent: int = 5, max_total: int = 8) -> list[dict]:
    """Most recent decisions plus any De Novo grant (the category-creating decision)."""
    recs = sorted(records, key=lambda r: r.get("decision_date") or "", reverse=True)
    chosen = recs[:n_recent]
    for r in recs[n_recent:]:
        if len(chosen) >= max_total:
            break
        if submission_kind(r.get("k_number", "")) == "de_novo":
            chosen.append(r)
    return [
        {"clearance_id": r.get("k_number"), "device_name": r.get("device_name"),
         "applicant": r.get("applicant"), "decision_date": r.get("decision_date"),
         "submission_kind": submission_kind(r.get("k_number", ""))}
        for r in chosen
    ]


def pma_example_devices(records: list[dict], n: int = 5) -> list[dict]:
    originals = [r for r in records if not (r.get("supplement_number") or "").strip()]
    originals = sorted(originals, key=lambda r: r.get("decision_date") or "", reverse=True)[:n]
    return [
        {"clearance_id": r.get("pma_number"), "device_name": r.get("trade_name"),
         "applicant": r.get("applicant"), "decision_date": r.get("decision_date"),
         "submission_kind": "pma"}
        for r in originals
    ]


def regulatory_visibility(n_510k, n_denovo, n_pma, has_code: bool) -> str:
    if not has_code:
        return "no_product_code_found"
    total = sum(int(x) for x in (n_510k, n_denovo, n_pma) if x is not None and not pd.isna(x))
    return "decisions_on_record" if total > 0 else "category_exists_no_decisions_in_openfda"


# ---------------------------------------------------------------------------------------------
# API access (cached via measure_it.http; every call is logged)
# ---------------------------------------------------------------------------------------------
class QueryLog:
    def __init__(self) -> None:
        self.entries: dict[str, dict] = {}

    def add(self, entry: dict) -> None:
        """Insert or merge; the same query issued for several classes keeps all their ids."""
        old = self.entries.get(entry["query_id"], {})
        merged = {**old, **entry}
        mids = list(old.get("measurement_ids") or [])
        for m in ([entry["measurement_id"]] if entry.get("measurement_id") else []) + list(entry.get("measurement_ids") or []):
            if m not in mids:
                mids.append(m)
        if mids:
            merged["measurement_ids"] = mids
        self.entries[entry["query_id"]] = merged

    def find(self, endpoint: str, search: str) -> dict | None:
        for e in self.entries.values():
            if e["endpoint"] == endpoint and e["params"].get("search") == search:
                return e
        return None

    def latest_fetch(self, queries: list[str] | None = None) -> str | None:
        """Latest actual fetch time (from the HTTP cache metadata) over all calls, or over the
        given 'endpoint: search' strings. Used so provenance never stamps the rebuild time."""
        if queries is None:
            times = [e.get("fetched_at") for e in self.entries.values()]
        else:
            times = []
            for q in queries:
                m = re.match(r"^(classification|510k): (.+)$", str(q))
                e = self.find(m.group(1), m.group(2)) if m else None
                times.append(e.get("fetched_at") if e else None)
        times = [t for t in times if t]
        return max(times) if times else None


LOG = QueryLog()
_LAST_UPDATED: dict[str, str] = {}


def api_call(endpoint: str, params: dict, *, purpose: str, context: dict | None = None) -> dict:
    """GET api.fda.gov/device/<endpoint>.json. 404 NOT_FOUND = zero matches (cached, logged)."""
    url = f"{API}/{endpoint}.json"
    r = request("GET", url, params=params, expect_json=True, cache_errors=True)
    payload = json.loads(r.content) if r.content else {}
    if r.status >= 400:
        err = payload.get("error", {}) if isinstance(payload, dict) else {}
        if not (r.status == 404 and err.get("code") == "NOT_FOUND"):
            raise RuntimeError(f"openFDA {endpoint} {params}: HTTP {r.status} {r.text[:300]}")
        payload = {"meta": {}, "results": [], "_not_found": True}
    meta = payload.get("meta", {}) or {}
    if meta.get("last_updated"):
        _LAST_UPDATED[endpoint] = meta["last_updated"]
    results = payload.get("results", []) or []
    total = (meta.get("results") or {}).get("total")
    if payload.get("_not_found"):
        total = 0
    if "count" in params:
        total = None  # a count query returns buckets, not a record total
    qid = query_id(endpoint, params)
    LOG.add({
        "query_id": qid, "endpoint": endpoint, "url": url, "params": params, "http_status": r.status,
        "total": total, "n_returned": len(results), "fetched_at": r.fetched_at, "from_cache": r.from_cache,
        "last_updated": meta.get("last_updated"), "purpose": purpose, **(context or {}),
    })
    payload["_query_id"] = qid
    payload["_fetched_at"] = r.fetched_at
    return payload


def _save_json(path: Path, obj) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, sort_keys=False, default=str))
    return path


def _raw() -> Path:
    return raw_dir(SOURCE_ID)


# ---------------------------------------------------------------------------------------------
# Step 1: classification snapshot
# ---------------------------------------------------------------------------------------------
def load_classification_bulk() -> tuple[dict[str, dict], dict]:
    path = download_file(CLASSIFICATION_BULK_URL, SOURCE_ID, min_bytes=100_000)
    with zipfile.ZipFile(path) as zf:
        name = [n for n in zf.namelist() if n.endswith(".json")][0]
        data = json.loads(zf.read(name))
    return {r["product_code"]: r for r in data["results"]}, data.get("meta", {})


# ---------------------------------------------------------------------------------------------
# Step 2: discovery
# ---------------------------------------------------------------------------------------------
def run_discovery() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classification keyword queries + named-device 510k queries. Returns (candidates, named)."""
    cand_rows, raw_disc = [], {}
    for mid, terms in DISCOVERY_TERMS.items():
        for term in terms:
            search = classification_search(term)
            params = {"search": search, "limit": API_MAX_LIMIT}
            p = api_call("classification", params, purpose="discovery_classification",
                         context={"measurement_id": mid, "term": term})
            qid = p["_query_id"]
            codes = [r["product_code"] for r in p["results"]]
            LOG.entries[qid]["result_product_codes"] = codes
            raw_disc[qid] = {"measurement_id": mid, "term": term, "params": params,
                             "total": LOG.entries[qid]["total"], "results": p["results"]}
            for r in p["results"]:
                cand_rows.append({
                    "measurement_id": mid, "term": term, "query": f"classification: {search}",
                    "query_id": qid, "product_code": r["product_code"], "device_name": r.get("device_name"),
                    "device_class": r.get("device_class"), "regulation_number": r.get("regulation_number"),
                    "medical_specialty": r.get("medical_specialty_description"),
                    "definition": (r.get("definition") or "")[:300],
                })
    _save_json(_raw() / "api" / "classification_discovery.json", raw_disc)

    named_rows, raw_named = [], {}
    for spec in NAMED_DEVICE_SEARCHES:
        params = {"search": spec["search"], "limit": API_MAX_LIMIT}
        p = api_call("510k", params, purpose="discovery_named_device",
                     context={"target": spec["target"], "measurement_ids": spec["measurement_ids"]})
        qid = p["_query_id"]
        LOG.entries[qid]["result_product_codes"] = sorted({r.get("product_code") for r in p["results"]})
        raw_named[qid] = {"target": spec["target"], "params": params,
                          "total": LOG.entries[qid]["total"], "results": p["results"]}
        for r in p["results"]:
            named_rows.append({
                "target": spec["target"], "search": spec["search"], "query": f"510k: {spec['search']}",
                "query_id": qid, "target_measurement_ids": json.dumps(spec["measurement_ids"]),
                "clearance_id": r.get("k_number"), "submission_kind": submission_kind(r.get("k_number", "")),
                "device_name": r.get("device_name"), "applicant": r.get("applicant"),
                "decision_date": r.get("decision_date"), "decision_code": r.get("decision_code"),
                "decision_description": r.get("decision_description"),
                "product_code": r.get("product_code"),
                "product_code_device_name": (r.get("openfda") or {}).get("device_name"),
                "fetched_at": p["_fetched_at"],
            })
    _save_json(_raw() / "api" / "named_device_searches.json", raw_named)
    cands = pd.DataFrame(cand_rows)
    named = pd.DataFrame(named_rows)
    # human-readable candidate list for curation
    interim = raw_dir(SOURCE_ID) / "discovery"
    interim.mkdir(exist_ok=True)
    if len(cands):
        (cands.drop_duplicates(["measurement_id", "product_code"])
         .sort_values(["measurement_id", "product_code"])
         .to_csv(interim / "classification_candidates.csv", index=False))
    if len(named):
        named.to_csv(interim / "named_device_hits.csv", index=False)
    return cands, named


def endpoint_totals() -> dict[str, int]:
    """Whole-endpoint record counts (context for the audit; also shows De Novo coverage)."""
    out = {}
    for label, ep, search in [("classification_records", "classification", None),
                              ("510k_records", "510k", None),
                              ("510k_de_novo_records", "510k", "k_number:DEN*"),
                              ("pma_records", "pma", None),
                              ("registrationlisting_records", "registrationlisting", None)]:
        params = {"limit": 1, **({"search": search} if search else {})}
        p = api_call(ep, params, purpose="endpoint_total", context={"term": label})
        out[label] = LOG.entries[p["_query_id"]]["total"]
    _save_json(_raw() / "api" / "endpoint_totals.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# Step 3: curated map
# ---------------------------------------------------------------------------------------------
def load_code_map(path: Path = CODE_MAP_PATH) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def iter_mappings(code_map: dict):
    """Yield (measurement_id, entry) for every entry, including 'no specific product code' ones."""
    for mid, entries in (code_map.get("measurements") or {}).items():
        for e in entries or []:
            yield mid, e


def validate_code_map(code_map: dict, classification: dict[str, dict], log: QueryLog,
                      measurement_ids: set[str]) -> list[str]:
    """Return a list of problems (empty = valid)."""
    problems = []
    mapped_mids = set((code_map.get("measurements") or {}).keys())
    for mid in measurement_ids - mapped_mids:
        problems.append(f"{mid}: no entry (need codes or an explicit '{NO_CODE}')")
    for mid in mapped_mids - measurement_ids:
        problems.append(f"{mid}: not a measurement class in configs/measurements.yaml")
    for mid, e in iter_mappings(code_map):
        code = e.get("product_code")
        if not code:
            if e.get("status") != NO_CODE:
                problems.append(f"{mid}: entry without product_code must have status '{NO_CODE}'")
            if not e.get("queries_tried"):
                problems.append(f"{mid}: '{NO_CODE}' entry must list queries_tried")
            for q in e.get("queries_tried") or []:
                m = re.match(r"^(classification|510k): (.+)$", str(q))
                if not m or log.find(m.group(1), m.group(2)) is None:
                    problems.append(f"{mid}: queries_tried {q!r} is not in the query log")
            continue
        for k in ("device_name", "device_class", "regulation_number", "medical_specialty_description",
                  "mapping_rationale", "query_that_found_it", "mapping_confidence"):
            if k not in e:
                problems.append(f"{mid}/{code}: missing {k}")
        rec = classification.get(code)
        if rec is None:
            problems.append(f"{mid}/{code}: not in FDA classification file")
            continue
        if rec.get("device_name", "").strip().lower() != str(e.get("device_name", "")).strip().lower():
            problems.append(f"{mid}/{code}: device_name {e.get('device_name')!r} != FDA {rec.get('device_name')!r}")
        if e.get("mapping_confidence") not in {"high", "medium", "low"}:
            problems.append(f"{mid}/{code}: mapping_confidence must be high/medium/low")
        q = str(e.get("query_that_found_it", ""))
        m = re.match(r"^(classification|510k): (.+)$", q)
        if not m:
            problems.append(f"{mid}/{code}: query_that_found_it must start with 'classification: ' or '510k: '")
            continue
        entry = log.find(m.group(1), m.group(2))
        if entry is None:
            problems.append(f"{mid}/{code}: query {q!r} is not in the query log")
        elif code not in (entry.get("result_product_codes") or []):
            problems.append(f"{mid}/{code}: query {q!r} did not return {code}")
    return problems


# ---------------------------------------------------------------------------------------------
# Step 4: per-code records
# ---------------------------------------------------------------------------------------------
def _fname(code: str, device_filter: str | None) -> str:
    return code if not device_filter else f"{code}__f{hashlib.sha1(device_filter.encode()).hexdigest()[:8]}"


def fetch_510k(code: str, device_filter: str | None = None) -> dict:
    search = product_code_search(code, device_filter)
    ctx = {"product_code": code, "device_filter": device_filter}
    first = api_call("510k", {"search": search, "limit": API_MAX_LIMIT, "sort": "decision_date:desc"},
                     purpose="510k_records", context=ctx)
    total = LOG.entries[first["_query_id"]]["total"] or 0
    if total <= API_MAX_LIMIT:
        records, mode = first["results"], "all_records"
        n_denovo = sum(submission_kind(r.get("k_number", "")) == "de_novo" for r in records)
        n_510k = sum(submission_kind(r.get("k_number", "")) == "510k" for r in records)
    else:
        records, mode = first["results"][:RECENT_N], f"most_recent_{RECENT_N}_of_{total}"
        dn = api_call("510k", {"search": f"{search} AND k_number:DEN*", "limit": 1},
                      purpose="510k_denovo_count", context=ctx)
        n_denovo = LOG.entries[dn["_query_id"]]["total"] or 0
        n_510k = total - n_denovo
    out = {"product_code": code, "device_filter": device_filter, "search": search, "total": total,
           "n_510k": n_510k, "n_denovo": n_denovo, "retrieval_mode": mode, "records": records,
           "fetched_at": first["_fetched_at"], "query_id": first["_query_id"]}
    _save_json(_raw() / "api" / "510k" / f"{_fname(code, device_filter)}.json", out)
    return out


def fetch_pma(code: str, device_filter: str | None = None) -> dict:
    search = product_code_search(code, device_filter.replace("device_name:", "trade_name:")
                                 if device_filter else None)
    p = api_call("pma", {"search": search, "limit": API_MAX_LIMIT}, purpose="pma_records",
                 context={"product_code": code, "device_filter": device_filter})
    total = LOG.entries[p["_query_id"]]["total"] or 0
    records = p["results"]
    mode = "all_records" if total <= API_MAX_LIMIT else f"first_{API_MAX_LIMIT}_of_{total}"
    n_orig = len({r.get("pma_number") for r in records if not (r.get("supplement_number") or "").strip()})
    out = {"product_code": code, "device_filter": device_filter, "search": search, "total_records": total,
           "n_original_pma": n_orig, "n_pma_numbers": len({r.get("pma_number") for r in records}),
           "retrieval_mode": mode, "records": records, "fetched_at": p["_fetched_at"],
           "query_id": p["_query_id"]}
    _save_json(_raw() / "api" / "pma" / f"{_fname(code, device_filter)}.json", out)
    return out


def registration_search(code: str, k_numbers: list[str] | None = None) -> str:
    """Listings that carry this product code; for a filtered generic code, only the listings
    linked to the filtered clearances. registrationlisting records carry K/DEN numbers in
    `k_number` and PMA numbers in a separate `pma_number` field, so P numbers go there."""
    base = f"products.product_code:{code}"
    if k_numbers is None:
        return base
    ks = sorted(n for n in k_numbers if not n.upper().startswith("P"))
    ps = sorted(n for n in k_numbers if n.upper().startswith("P"))
    clauses = ([f"k_number:({' '.join(ks)})"] if ks else []) + ([f"pma_number:({' '.join(ps)})"] if ps else [])
    if len(clauses) == 1:
        return f"{base} AND {clauses[0]}"
    return f"{base} AND ({' OR '.join(clauses)})"


def fetch_registration_counts(code: str, device_filter: str | None = None,
                              k_numbers: list[str] | None = None) -> dict:
    """Device listings and distinct establishment registrations for a product code (all and US)."""
    ctx = {"product_code": code, "device_filter": device_filter}
    if device_filter is not None and not k_numbers:
        out = {"product_code": code, "device_filter": device_filter, "search": None, "n_listings": 0,
               "n_establishments": 0, "n_listings_us": 0, "n_establishments_us": 0,
               "establishments_truncated": False, "fetched_at": utc_now_iso(), "query_id": None,
               "note": "no clearances matched the device filter, so no listings can be linked"}
    else:
        search = registration_search(code, k_numbers if device_filter is not None else None)
        a = api_call("registrationlisting", {"search": search, "count": "registration.registration_number",
                                             "limit": API_MAX_LIMIT},
                     purpose="registration_listing_count", context=ctx)
        u = api_call("registrationlisting", {"search": f"{search} AND registration.iso_country_code:US",
                                             "count": "registration.registration_number", "limit": API_MAX_LIMIT},
                     purpose="registration_listing_count_us", context=ctx)
        out = {"product_code": code, "device_filter": device_filter, "search": search,
               "n_listings": int(sum(b["count"] for b in a["results"])),
               "n_establishments": len(a["results"]),
               "n_listings_us": int(sum(b["count"] for b in u["results"])),
               "n_establishments_us": len(u["results"]),
               "establishments_truncated": len(a["results"]) >= API_MAX_LIMIT,
               "fetched_at": a["_fetched_at"], "query_id": a["_query_id"], "note": ""}
    _save_json(_raw() / "api" / "registrationlisting" / f"{_fname(code, device_filter)}.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# Step 5: tables
# ---------------------------------------------------------------------------------------------
def _version(endpoint: str) -> str:
    return f"openFDA device/{endpoint} last_updated {_LAST_UPDATED.get(endpoint, UNKNOWN)}"


def build_tables(code_map: dict, classification: dict[str, dict], class_meta: dict,
                 named: pd.DataFrame) -> dict[str, pd.DataFrame]:
    measurements = {m["id"]: m for m in load_config("measurements")["measurement_classes"]}
    mappings = list(iter_mappings(code_map))
    coded = [(mid, e) for mid, e in mappings if e.get("product_code")]

    # unique (code, filter) fetch units
    units: dict[tuple[str, str | None], list[str]] = defaultdict(list)
    for mid, e in coded:
        units[(e["product_code"], e.get("device_filter"))].append(mid)
    codes = sorted({c for c, _ in units})

    k510, pma, reg = {}, {}, {}
    for (code, flt) in sorted(units, key=lambda u: (u[0], u[1] or "")):
        k510[(code, flt)] = fetch_510k(code, flt)
        pma[(code, flt)] = fetch_pma(code, flt)
        linked = None
        if flt:
            linked = [r["k_number"] for r in k510[(code, flt)]["records"]]
            linked += sorted({r["pma_number"] for r in pma[(code, flt)]["records"]})
        reg[(code, flt)] = fetch_registration_counts(code, flt, linked)

    # a re-curated map (changed code or filter) must not leave orphaned per-unit responses behind
    expected = {f"{_fname(code, flt)}.json" for code, flt in units}
    for sub in ("510k", "pma", "registrationlisting"):
        for p in (_raw() / "api" / sub).glob("*.json"):
            if p.name not in expected:
                p.unlink()

    mids_by_code = defaultdict(set)
    for (code, _), mids in units.items():
        mids_by_code[code].update(mids)

    # --- fda_device_classification: one row per mapped product code
    cls_rows = []
    for code in codes:
        r = classification[code]
        cls_rows.append({
            "product_code": code, "device_name": r.get("device_name"), "definition": r.get("definition") or "",
            "device_class": r.get("device_class"), "regulation_number": r.get("regulation_number") or "",
            "medical_specialty": r.get("medical_specialty"),
            "medical_specialty_description": r.get("medical_specialty_description"),
            "review_panel": r.get("review_panel"),
            "submission_type_id": r.get("submission_type_id"),
            "submission_type": SUBMISSION_TYPES.get(str(r.get("submission_type_id")),
                                                    f"undocumented value {r.get('submission_type_id')}"),
            "implant_flag": r.get("implant_flag"), "life_sustain_support_flag": r.get("life_sustain_support_flag"),
            "third_party_flag": r.get("third_party_flag"), "gmp_exempt_flag": r.get("gmp_exempt_flag"),
            "unclassified_reason": r.get("unclassified_reason") or "",
            "mapped_measurement_ids": json.dumps(sorted(mids_by_code[code])),
        })
    cls = pd.DataFrame(cls_rows)
    cls_retrieved = load_manifest(SOURCE_ID)["files"].get(CLASSIFICATION_BULK_URL.rsplit("/", 1)[-1], {}).get(
        "retrieved_at", utc_now_iso())
    cls = add_provenance(
        cls, data_layer="facility", source_name="openFDA device classification (bulk file)",
        source_version=f"device-classification-0001-of-0001 last_updated {class_meta.get('last_updated', UNKNOWN)}",
        retrieved_at=cls_retrieved, evidence_type="regulatory_record", source_record_id="product_code",
        evidence_level="fda_product_code",
        provenance_notes="FDA product-code category; curated to measurement classes in configs/fda_product_code_map.yaml. "
                         + REGULATORY_NOTE)

    # --- fda_510k: one row per K/DEN number; a record reached through several fetch units
    # (e.g. DPS unfiltered for ECG and DPS[HRV filter] for HRV) keeps every unit's filter and ids.
    by_k: dict[str, dict] = {}
    for (code, flt), res in k510.items():
        for r in res["records"]:
            k = r.get("k_number")
            if k in by_k:
                row = by_k[k]
                row["_filters"].add(flt or "")
                row["_mids"].update(units[(code, flt)])
                continue
            by_k[k] = {
                "clearance_id": k, "submission_kind": submission_kind(k), "product_code": r.get("product_code"),
                "device_name": r.get("device_name"), "applicant": r.get("applicant"),
                "decision_date": r.get("decision_date"), "date_received": r.get("date_received"),
                "decision_code": r.get("decision_code"), "decision_description": r.get("decision_description"),
                "clearance_type": r.get("clearance_type"), "statement_or_summary": r.get("statement_or_summary"),
                "third_party_flag": r.get("third_party_flag"), "expedited_review_flag": r.get("expedited_review_flag"),
                "applicant_city": r.get("city"), "applicant_state": r.get("state"),
                "applicant_country_code": r.get("country_code"), "advisory_committee": r.get("advisory_committee"),
                "retrieval_mode": res["retrieval_mode"], "retrieved_at": res["fetched_at"],
                "_filters": {flt or ""}, "_mids": set(units[(code, flt)]),
            }
    rows = []
    for row in by_k.values():
        row["device_filters"] = json.dumps(sorted(row.pop("_filters")))
        row["mapped_measurement_ids"] = json.dumps(sorted(row.pop("_mids")))
        rows.append(row)
    k510_df = pd.DataFrame(rows)
    k510_ret = k510_df.pop("retrieved_at") if len(k510_df) else pd.Series(dtype=str)
    k510_df = add_provenance(
        k510_df, data_layer="facility", source_name="openFDA device/510k (includes De Novo decisions)",
        source_version=_version("510k"), retrieved_at="", evidence_type="regulatory_record",
        source_record_id="clearance_id",
        evidence_level=k510_df["submission_kind"].map({"510k": "510k_substantially_equivalent_or_other_decision",
                                                       "de_novo": "de_novo_granted"}).fillna("other"),
        provenance_notes="Applicant address is the submitter's address, not a deployment site. " + REGULATORY_NOTE)
    k510_df["retrieved_at"] = k510_ret.values

    # --- fda_pma (same aggregation as fda_510k)
    by_p: dict[tuple, dict] = {}
    for (code, flt), res in pma.items():
        for r in res["records"]:
            key = (r.get("pma_number"), r.get("supplement_number") or "")
            if key in by_p:
                by_p[key]["_filters"].add(flt or "")
                by_p[key]["_mids"].update(units[(code, flt)])
                continue
            by_p[key] = {
                "pma_number": r.get("pma_number"), "supplement_number": r.get("supplement_number") or "",
                "is_original": not (r.get("supplement_number") or "").strip(),
                "product_code": r.get("product_code"), "trade_name": r.get("trade_name"),
                "generic_name": r.get("generic_name"), "applicant": r.get("applicant"),
                "decision_date": r.get("decision_date"), "decision_code": r.get("decision_code"),
                "supplement_type": r.get("supplement_type"), "supplement_reason": r.get("supplement_reason"),
                "applicant_city": r.get("city"), "applicant_state": r.get("state"),
                "retrieval_mode": res["retrieval_mode"], "retrieved_at": res["fetched_at"],
                "_filters": {flt or ""}, "_mids": set(units[(code, flt)]),
            }
    prow = []
    for row in by_p.values():
        row["device_filters"] = json.dumps(sorted(row.pop("_filters")))
        row["mapped_measurement_ids"] = json.dumps(sorted(row.pop("_mids")))
        prow.append(row)
    pma_cols = ["pma_number", "supplement_number", "is_original", "product_code", "trade_name", "generic_name",
                "applicant", "decision_date", "decision_code", "supplement_type", "supplement_reason",
                "applicant_city", "applicant_state", "retrieval_mode", "retrieved_at", "device_filters",
                "mapped_measurement_ids"]
    pma_df = pd.DataFrame(prow, columns=pma_cols)
    pma_ret = pma_df.pop("retrieved_at")
    pma_df["pma_record_id"] = pma_df["pma_number"].astype(str) + "/" + pma_df["supplement_number"].astype(str)
    pma_df = add_provenance(
        pma_df, data_layer="facility", source_name="openFDA device/pma", source_version=_version("pma"),
        retrieved_at="", evidence_type="regulatory_record", source_record_id="pma_record_id",
        evidence_level=pma_df["is_original"].map({True: "pma_original_approval", False: "pma_supplement"}),
        provenance_notes=REGULATORY_NOTE)
    pma_df["retrieved_at"] = pma_ret.values

    # --- fda_registration_listing_counts (one row per fetch unit: product code [+ device filter])
    rrows = []
    for (code, flt), res in reg.items():
        rrows.append({
            "product_code": code, "device_name": classification[code].get("device_name"),
            "device_filter": flt or "",
            "n_device_listings": res["n_listings"], "n_establishments": res["n_establishments"],
            "n_device_listings_us": res["n_listings_us"], "n_establishments_us": res["n_establishments_us"],
            "n_establishments_non_us": res["n_establishments"] - res["n_establishments_us"],
            "establishments_truncated_at_1000": res["establishments_truncated"],
            "query": (f"registrationlisting: search={res['search']} count=registration.registration_number"
                      if res["search"] else ""),
            "note": res["note"], "retrieved_at": res["fetched_at"],
        })
    reg_df = pd.DataFrame(rrows)
    reg_df["reg_record_id"] = reg_df["product_code"] + reg_df["device_filter"].map(lambda f: f"[{f}]" if f else "")
    reg_ret = reg_df.pop("retrieved_at")
    reg_df = add_provenance(
        reg_df, data_layer="facility", source_name="openFDA device/registrationlisting",
        source_version=_version("registrationlisting"), retrieved_at="", evidence_type="regulatory_record",
        source_record_id="reg_record_id", evidence_level="establishment_listing_count",
        provenance_notes="n_device_listings = registration & listing records carrying the product code; "
                         "n_establishments = distinct establishment registration numbers among them (manufacturers, "
                         "specification developers, relabelers, importers...). For rows with a device_filter the "
                         "listings are restricted to those linked to the filtered K/DEN/P numbers. Not a count of "
                         "U.S. deployment sites.")
    reg_df["retrieved_at"] = reg_ret.values

    # --- measurement_regulatory_status
    srows = []
    for mid, e in mappings:
        m = measurements.get(mid, {})
        code = e.get("product_code")
        base = {"measurement_id": mid, "measurement_name": m.get("name", UNKNOWN)}
        if not code:
            srows.append({**base, "product_code": None, "device_name": NO_CODE, "device_class": None,
                          "regulation_number": None, "medical_specialty": None, "submission_type": None,
                          "device_filter": "", "n_510k": pd.NA, "n_denovo": pd.NA, "n_pma": pd.NA,
                          "n_pma_supplements": pd.NA, "n_pma_numbers_any_record": pd.NA,
                          "n_510k_records_retrieved": pd.NA,
                          "record_retrieval_mode": "", "n_registered_establishments": pd.NA,
                          "n_registered_establishments_us": pd.NA,
                          "most_recent_decision_date": None, "example_devices": "[]",
                          "mapping_confidence": "none", "mapping_rationale": e.get("note", ""),
                          "query_that_found_it": json.dumps(e.get("queries_tried", [])),
                          "regulatory_visibility": regulatory_visibility(None, None, None, False),
                          "regulatory_note": REGULATORY_NOTE,
                          "retrieved_at": LOG.latest_fetch(e.get("queries_tried", [])) or utc_now_iso()})
            continue
        flt = e.get("device_filter")
        kr, pr, rr = k510[(code, flt)], pma[(code, flt)], reg[(code, flt)]
        dates = [r.get("decision_date") for r in kr["records"] if r.get("decision_date")]
        dates += [r.get("decision_date") for r in pr["records"] if r.get("decision_date")
                  and not (r.get("supplement_number") or "").strip()]
        ex = example_devices(kr["records"]) + pma_example_devices(pr["records"])
        cr = classification[code]
        srows.append({
            **base, "product_code": code, "device_name": cr.get("device_name"),
            "device_class": cr.get("device_class"), "regulation_number": cr.get("regulation_number") or "",
            "medical_specialty": cr.get("medical_specialty_description"),
            "submission_type": SUBMISSION_TYPES.get(str(cr.get("submission_type_id")),
                                                    f"undocumented value {cr.get('submission_type_id')}"),
            "device_filter": flt or "", "n_510k": kr["n_510k"], "n_denovo": kr["n_denovo"],
            "n_pma": pr["n_original_pma"], "n_pma_supplements": pr["total_records"] - pr["n_original_pma"],
            "n_pma_numbers_any_record": pr["n_pma_numbers"],
            "n_510k_records_retrieved": len(kr["records"]), "record_retrieval_mode": kr["retrieval_mode"],
            "n_registered_establishments": rr["n_establishments"],
            "n_registered_establishments_us": rr["n_establishments_us"],
            "most_recent_decision_date": max(dates) if dates else None,
            "example_devices": json.dumps(ex),
            "mapping_confidence": e["mapping_confidence"], "mapping_rationale": e["mapping_rationale"],
            "query_that_found_it": e["query_that_found_it"],
            "regulatory_visibility": regulatory_visibility(kr["n_510k"], kr["n_denovo"], pr["n_pma_numbers"], True),
            "regulatory_note": REGULATORY_NOTE, "retrieved_at": kr["fetched_at"],
        })
    st = pd.DataFrame(srows)
    for c in ("n_510k", "n_denovo", "n_pma", "n_pma_supplements", "n_pma_numbers_any_record",
              "n_510k_records_retrieved",
              "n_registered_establishments", "n_registered_establishments_us"):
        st[c] = st[c].astype("Int64")
    st_ret = st.pop("retrieved_at")
    st["status_record_id"] = st["measurement_id"] + ":" + st["product_code"].fillna("NONE") + (
        st["device_filter"].map(lambda f: f"[{f}]" if f else ""))
    st = add_provenance(
        st, data_layer="measurement", source_name="openFDA device API + curated configs/fda_product_code_map.yaml",
        source_version=f"{_version('510k')}; {_version('pma')}; {_version('registrationlisting')}",
        retrieved_at="", evidence_type="regulatory_record", source_record_id="status_record_id",
        evidence_level=st["regulatory_visibility"],
        provenance_notes="Product-code-to-measurement mapping is a human curation (mapping_confidence, "
                         "mapping_rationale); counts are measured from openFDA. Counts for rows with a "
                         "device_filter are restricted to records matching that filter. " + REGULATORY_NOTE)
    st["retrieved_at"] = st_ret.values

    # --- fda_named_device_findings
    mapped_codes = set(codes)
    nd = named.copy()
    if len(nd):
        # A record counts toward a measurement only if it falls inside a mapped fetch unit:
        # retrieved under that unit (so generic codes respect their device_filter), or - for a
        # unit where only the most recent records were kept - its code is mapped unfiltered.
        unfiltered = defaultdict(set)
        for (code, flt), mids in units.items():
            if not flt:
                unfiltered[code].update(mids)

        def _mids(row) -> list[str]:
            if row["clearance_id"] in by_k:
                return json.loads(by_k[row["clearance_id"]]["mapped_measurement_ids"])
            return sorted(unfiltered.get(row["product_code"], set()))

        nd["product_code_in_mapped_set"] = nd["product_code"].isin(mapped_codes)
        nd["mapped_measurement_ids"] = nd.apply(lambda r: json.dumps(_mids(r)), axis=1)
        nd["record_counted_in_status"] = nd["mapped_measurement_ids"] != "[]"
        nd_ret = nd.pop("fetched_at")
        nd["finding_id"] = nd["query_id"] + ":" + nd["clearance_id"].astype(str)
        nd = add_provenance(
            nd, data_layer="facility", source_name="openFDA device/510k (named-device searches)",
            source_version=_version("510k"), retrieved_at="", evidence_type="regulatory_record",
            source_record_id="finding_id", evidence_level=nd["submission_kind"],
            provenance_notes="Every record returned by a named device/applicant search, including unrelated "
                             "products of the same applicant (see record_counted_in_status / "
                             "product_code_in_mapped_set). " + REGULATORY_NOTE)
        nd["retrieved_at"] = nd_ret.values

    # --- fda_device_query_log
    ql = pd.DataFrame([
        {"query_id": e["query_id"], "endpoint": e["endpoint"], "url": e["url"],
         "params": json.dumps(e["params"], sort_keys=True), "search": e["params"].get("search"),
         "http_status": e["http_status"], "total": e.get("total"), "n_returned": e["n_returned"],
         "purpose": e["purpose"], "measurement_ids": json.dumps(e.get("measurement_ids") or []), "term": e.get("term"),
         "target": e.get("target"), "product_code": e.get("product_code"),
         "device_filter": e.get("device_filter") or "",
         "result_product_codes": json.dumps(e.get("result_product_codes")) if e.get("result_product_codes") is not None else "",
         "api_last_updated": e.get("last_updated"), "retrieved_at": e["fetched_at"]}
        for e in LOG.entries.values()
    ])
    ql["total"] = ql["total"].astype("Int64")
    ql_ret = ql.pop("retrieved_at")
    ql = add_provenance(
        ql, data_layer="facility", source_name="openFDA device API query log", source_version="api.fda.gov/device",
        retrieved_at="", evidence_type="metadata_catalog", source_record_id="query_id",
        evidence_level="query_record", provenance_notes="One row per distinct openFDA API call made by this module.")
    ql["retrieved_at"] = ql_ret.values

    return {
        "fda_device_classification": cls, "fda_510k": k510_df, "fda_pma": pma_df,
        "fda_registration_listing_counts": reg_df, "measurement_regulatory_status": st,
        "fda_named_device_findings": nd, "fda_device_query_log": ql,
    }


# ---------------------------------------------------------------------------------------------
# Audit + registry
# ---------------------------------------------------------------------------------------------
def write_query_log() -> Path:
    entries = sorted(LOG.entries.values(), key=lambda e: (e["purpose"], e["endpoint"], json.dumps(e["params"])))
    files = {}
    for p in sorted((_raw() / "api").rglob("*.json")):
        files[str(p.relative_to(_raw()))] = {"bytes": p.stat().st_size, "sha256": sha256_file(p)}
    return _save_json(_raw() / "QUERY_LOG.json", {
        "source_id": SOURCE_ID, "written_at": utc_now_iso(), "api_base": API,
        "api_last_updated": _LAST_UPDATED, "n_queries": len(entries), "queries": entries,
        "saved_api_responses": files,
    })


def _registry(tables: dict[str, pd.DataFrame], class_meta: dict, n_classification: int) -> None:
    st = tables["measurement_regulatory_status"]
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "openFDA device endpoints (classification, 510k incl. De Novo, PMA, registration & listing)",
        "publisher": "U.S. Food and Drug Administration (openFDA)",
        "landing_url": "https://open.fda.gov/apis/device/",
        "access_urls": [f"{API}/classification.json", f"{API}/510k.json", f"{API}/pma.json",
                        f"{API}/registrationlisting.json", CLASSIFICATION_BULK_URL],
        "license": "public domain, dedicated under Creative Commons CC0 1.0 Universal (https://open.fda.gov/license/)",
        "access_conditions": "open API, no key used; openFDA limits without a key: 240 requests/min and "
                             "1,000 requests/day per IP (https://open.fda.gov/apis/authentication/)",
        # latest actual API fetch (cached reruns keep the original retrieval time)
        "retrieved_at": LOG.latest_fetch() or utc_now_iso(),
        "source_version": "; ".join(f"{k} last_updated {v}" for k, v in sorted(_LAST_UPDATED.items()))
                          + f"; classification bulk last_updated {class_meta.get('last_updated')}",
        "update_date": "API meta.last_updated " + ", ".join(sorted(set(_LAST_UPDATED.values())))
                       + f"; classification bulk export last_updated {class_meta.get('last_updated')}",
        "data_layer": "facility",
        "unit_of_observation": "FDA product code; 510(k)/De Novo decision; PMA decision/supplement; "
                               "establishment registration count per product code",
        "sample_size": {
            "classification_product_codes_total": n_classification,
            "measurement_classes": int(st["measurement_id"].nunique()),
            "measurement_classes_with_mapped_code": int(st.loc[st["product_code"].notna(), "measurement_id"].nunique()),
            "mapped_product_codes": int(len(tables["fda_device_classification"])),
            "measurement_code_mappings": int(st["product_code"].notna().sum()),
            "fda_510k_records": int(len(tables["fda_510k"])),
            "fda_510k_de_novo_records": int((tables["fda_510k"]["submission_kind"] == "de_novo").sum()),
            "fda_pma_records": int(len(tables["fda_pma"])),
            "named_device_findings": int(len(tables["fda_named_device_findings"])),
            "api_queries": int(len(tables["fda_device_query_log"])),
        },
        "geographic_resolution": "none (applicant addresses are not deployment locations)",
        "person_level": False, "geographic": False, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (regulatory records, no people)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": list(tables.keys()),
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": PRODUCER,
        "limitations": [
            REGULATORY_NOTE,
            "Product-code-to-measurement mapping is a human curation of FDA category names; generic "
            "codes are mapped only with a device_filter and are low confidence.",
            "510(k)-exempt categories can have marketed devices with no 510(k) record.",
            "Consumer wellness wearables (e.g. step counters, sleep trackers without medical claims) are "
            "generally not FDA-cleared and are invisible here.",
            "De Novo decisions are taken from the 510k endpoint (DEN numbers); openFDA has no separate De Novo endpoint.",
            "Laboratory-developed tests (most omics assays) are not in these endpoints.",
            "Registration counts are establishment registrations whose current listings carry a code "
            "(manufacturers, relabelers, importers...), not distinct devices or deployment sites; older "
            "clearances no longer listed contribute 0.",
            "Applicant/establishment addresses are firm addresses, not clinics or deployment locations.",
        ],
    })


def run(discover_only: bool = False) -> dict[str, pd.DataFrame] | None:
    classification, class_meta = load_classification_bulk()
    totals = endpoint_totals()
    print("endpoint totals:", totals)
    _cands, named = run_discovery()
    if discover_only:
        write_query_log()
        print(f"discovery: {len(LOG.entries)} queries; candidates in {_raw() / 'discovery'}")
        return None
    code_map = load_code_map()
    mids = {m["id"] for m in load_config("measurements")["measurement_classes"]}
    problems = validate_code_map(code_map, classification, LOG, mids)
    if problems:
        raise ValueError("configs/fda_product_code_map.yaml invalid:\n  " + "\n  ".join(problems))
    tables = build_tables(code_map, classification, class_meta, named)
    for name, df in tables.items():
        write_table(df, name, producer=PRODUCER,
                    description=f"openFDA device data ({name}); {REGULATORY_NOTE}")
    write_query_log()
    _registry(tables, class_meta, len(classification))
    for name, df in tables.items():
        print(f"{name}: {len(df)} rows")
    return tables


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--discover-only", action="store_true",
                    help="run the classification/named-device discovery queries and stop (for curation)")
    run(discover_only=ap.parse_args().discover_only)
