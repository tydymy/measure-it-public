"""Condition ontology layer (SPEC section B).

Resolves every condition in configs/conditions.yaml to public ontology
identifiers. Nothing here is typed in from memory: every identifier comes from
a downloaded/queried public resource and carries its source, version and match
type.

Sources
  mondo_ontology  Mondo v2026-09-01 release (mondo.obo + mondo.sssom.tsv) from the
                  monarch-initiative/mondo GitHub release; Monarch KG API v3
                  disease-phenotype associations (Monarch publishes both).
  ebi_ols4        EMBL-EBI OLS4 API: EFO / MeSH / DOID / ORDO / HP / Mondo label
                  verification and exact-label searches.
  cms_icd10cm     CMS ICD-10-CM FY2026 code descriptions (April 1, 2026 update, in
                  force 2026-04-01..2026-09-30), plus FY2026 Oct-2025 and FY2027
                  files for change detection.

Method
  1. Each condition's names (preferred name and its parts, aliases, search
     terms) are matched against non-obsolete MONDO classes by exact label,
     exact synonym, or related/narrow/broad synonym (punctuation- and
     plural-insensitive). Acronym-only matches never choose a class (e.g.
     'PAIS' is an EXACT synonym of partial androgen insensitivity syndrome).
  2. The primary class is the best match on the highest-tier name; curated
     decisions in CURATION add secondary classes (component / narrower /
     broader_context / related_distinct) or make a manual choice, each with a
     rationale, and are checked against the downloaded release at build time.
  3. Mondo xrefs (with Mondo's source qualifiers and skos predicates) give
     ICD-10-CM, MeSH, EFO, UMLS, SNOMED CT, DOID, Orphanet, NCIT ... ids. The
     condition-level predicate composes the condition->class relation with the
     class->xref predicate.
  4. OLS4 verifies EFO/MeSH/DOID/ORDO labels and obsolescence, checks whether
     each Mondo class is imported into EFO, and adds exact-label EFO/MeSH hits.
  5. Every ICD-10-CM code is validated against the CMS FY2026 file. Clinically
     standard codes without a Mondo xref are added only when the CMS
     description text matches a condition name, with
     mapping_source='CMS ICD-10-CM FY2026 description match (curated)'.
  6. Monarch disease-phenotype associations (direct and descendant) give HPO
     terms with frequency/evidence/source; OLS4 HPO search resolves the
     wearable phenotype axes.

Outputs (data/processed): condition_registry, condition_ontology_mappings,
condition_phenotypes, phenotype_axes, ontology_icd10cm_codes.

Run: uv run python -m measure_it.ontology.build_registry
"""
from __future__ import annotations

import gzip
import json
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from ..config import UNKNOWN, load_config, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..http import get_json, request
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from .normalize import (
    PREDICATE_RANK, REGISTRY_PREDICATES, canonical_curie, compact_key, icd_dotted, icd_nodot,
    is_acronym, norm_text, singular_key, strip_parenthetical,
)

SOURCE_MONDO = "mondo_ontology"
SOURCE_OLS = "ebi_ols4"
SOURCE_CMS = "cms_icd10cm"
PRODUCER = "ontology.build_registry"

# Latest Mondo GitHub release on 2026-09-23 (GitHub releases/latest API). Pinned for
# reproducibility; build() reports if a newer release exists.
MONDO_TAG = "v2026-09-01"
MONDO_URLS = {
    "mondo.obo": f"https://github.com/monarch-initiative/mondo/releases/download/{MONDO_TAG}/mondo.obo",
    "mondo.sssom.tsv": f"https://raw.githubusercontent.com/monarch-initiative/mondo/{MONDO_TAG}/src/ontology/mappings/mondo.sssom.tsv",
    "source-versions.tsv": f"https://github.com/monarch-initiative/mondo/releases/download/{MONDO_TAG}/source-versions.tsv",
    f"mondo_release_{MONDO_TAG}.json": f"https://api.github.com/repos/monarch-initiative/mondo/releases/tags/{MONDO_TAG}",
}
MONDO_LATEST_API = "https://api.github.com/repos/monarch-initiative/mondo/releases/latest"
MONDO_OWL_VERSIONED = f"http://purl.obolibrary.org/obo/mondo/releases/{MONDO_TAG.lstrip('v')}/mondo.owl"
MONDO_OWL_HEADER_FILE = f"mondo_owl_header_first4KB_{MONDO_TAG}.xml"
MONDO_PURLS = ["http://purl.obolibrary.org/obo/mondo.obo", "http://purl.obolibrary.org/obo/mondo/mappings/mondo.sssom.tsv"]

CMS_URLS = {
    "fy2026_apr": "https://www.cms.gov/files/zip/april-1-2026-code-descriptions-tabular-order.zip",
    "fy2026_oct": "https://www.cms.gov/files/zip/2026-code-descriptions-tabular-order.zip",
    "fy2027": "https://www.cms.gov/files/zip/2027-code-descriptions-tabular-order.zip",
}
CMS_LANDING = "https://www.cms.gov/medicare/coding-billing/icd-10-codes"
CMS_VERSION = "ICD-10-CM FY2026 (April 1, 2026 update; valid 2026-04-01 to 2026-09-30)"
CMS_CURATED_SOURCE = "CMS ICD-10-CM FY2026 description match (curated)"

OLS = "https://www.ebi.ac.uk/ols4/api"
OLS_ONTOLOGIES = ["mondo", "efo", "mesh", "hp", "doid", "ordo"]
MONARCH = "https://api-v3.monarchinitiative.org/v3/api"

IRI_BASE = {
    "MONDO": "http://purl.obolibrary.org/obo/MONDO_",
    "HP": "http://purl.obolibrary.org/obo/HP_",
    "DOID": "http://purl.obolibrary.org/obo/DOID_",
    "EFO": "http://www.ebi.ac.uk/efo/EFO_",
    "ORPHANET": "http://www.orpha.net/ORDO/Orphanet_",
    "MESH": "http://id.nlm.nih.gov/mesh/",
}
OLS_ONTOLOGY_FOR = {"EFO": "efo", "MESH": "mesh", "DOID": "doid", "ORPHANET": "ordo", "HP": "hp", "MONDO": "mondo"}

FLAG_COLUMNS = ["infectious_or_post_infectious", "autoimmune", "autonomic", "vascular", "pain",
                "fatigue_pem", "gi", "neurologic"]

# --------------------------------------------------------------------------- curated decisions
# Relations are from the canonical condition to the target:
#   exact  = same concept; close = near-equivalent; narrow = target is narrower (a subtype or one
#   part of a composite condition); broad = target is broader; related = associated, not nested.
#
# NAME_RELATION: how a config name relates to the condition (default: preferred-name parts exact,
# aliases/search terms close). Used when a name, rather than a Mondo class, is matched in OLS.
NAME_RELATION = {
    ("dysautonomia", "orthostatic intolerance"): "narrow",
    ("dysautonomia", "autonomic neuropathy"): "narrow",
    ("eds_hsd", "Ehlers-Danlos syndromes"): "narrow",
    ("eds_hsd", "hypermobility spectrum disorders"): "narrow",
    ("eds_hsd", "Ehlers-Danlos syndrome"): "narrow",
    ("eds_hsd", "hypermobile Ehlers-Danlos syndrome"): "narrow",
    ("eds_hsd", "hypermobility spectrum disorder"): "narrow",
    ("eds_hsd", "joint hypermobility syndrome"): "narrow",
    ("eds_hsd", "Ehlers-Danlos"): "narrow",
    ("eds_hsd", "hypermobile Ehlers-Danlos"): "narrow",
    ("mcas", "idiopathic mast cell activation syndrome"): "narrow",
    ("migraine", "chronic migraine"): "narrow",
    ("migraine", "vestibular migraine"): "narrow",
    ("endometriosis", "endometriosis of uterus"): "narrow",
    ("long_covid", "post-acute sequelae of SARS-CoV-2 infection"): "exact",
    ("long_covid", "post-COVID-19 condition"): "exact",
    ("endometriosis", "adenomyosis"): "related",
    ("ptlds", "chronic Lyme disease"): "related",
    ("ptlds", "persistent symptoms after Lyme"): "close",
}

# CURATION: per condition
#   primary:  (mondo_id, expected_label, rationale) -> manual choice (only where no name matches)
#   relation: (relation, rationale) override for the primary class
#   classes:  [(mondo_id, expected_label, role, relation, rationale)] secondary classes
#   icd10cm:  [(code, relation, semantic_rationale_or_None)] CMS codes added without a Mondo xref.
#             semantic_rationale None => the CMS description must lexically match a condition name.
CURATION: dict[str, dict] = {
    "long_covid": {
        "classes": [("MONDO:0100320", "post-COVID-19 disorder", "broader_context", "broad",
                     "is_a parent of the primary class MONDO:0100233; groups all disorders arising after COVID-19")],
        "icd10cm": [("U09.9", "exact", None), ("U09", "exact", None)],
    },
    "me_cfs": {
        "icd10cm": [("G93.32", "exact", None)],
        "note": "Mondo v2026-09-01 spells the label 'myalgic encephalomeyelitis/chronic fatigue syndrome'; "
                "matched through exact synonyms. Mondo has no ICD-10-CM xref for this class.",
    },
    "pots": {
        "relation": ("broad", "Mondo v2026-09-01 has no stand-alone POTS class: 'postural orthostatic tachycardia syndrome' "
                              "and 'POTS' are EXACT synonyms of MONDO:0001315 'orthostatic intolerance', which also carries "
                              "'neurocirculatory asthenia', 'soldiers heart' and 'irritable heart' synonyms and a MeSH xref "
                              "to Neurocirculatory Asthenia. The class is treated as broader than POTS, so its xrefs are "
                              "broadMatch for this condition."),
        "classes": [("MONDO:0011479", "postural orthostatic tachycardia syndrome due to NET deficiency", "narrower", "narrow",
                     "monogenic norepinephrine-transporter-deficiency form; 'postural orthostatic tachycardia syndrome' is only a BROAD synonym")],
        "icd10cm": [("G90.A", "exact", None)],
    },
    "dysautonomia": {
        "classes": [
            ("MONDO:0001292", "autonomic nervous system disorder", "component", "exact",
             "named in the config preferred name 'Dysautonomia (autonomic nervous system disorder)'; parent of MONDO:0044872"),
            ("MONDO:0001315", "orthostatic intolerance", "narrower", "narrow", "config alias (exact label); is_a MONDO:0044872"),
            ("MONDO:0001300", "autonomic neuropathy", "narrower", "narrow", "config alias (exact label); is_a MONDO:0001292"),
            ("MONDO:0021809", "primary dysautonomia", "narrower", "narrow",
             "is_a MONDO:0044872; 'dysautonomia' is a RELATED synonym of this class"),
        ],
        "icd10cm": [("G90.9", "close", None), ("G90.89", "close", None)],
    },
    "eds_hsd": {
        "relation": ("narrow", "the condition is the union of the Ehlers-Danlos syndromes and hypermobility spectrum "
                               "disorders; each Mondo class covers one part"),
        "classes": [
            ("MONDO:0007523", "Ehlers-Danlos syndrome, hypermobility type", "component", "narrow",
             "hEDS; 'hypermobile Ehlers-Danlos syndrome' and 'hEDS' are RELATED synonyms of this class in Mondo"),
            ("MONDO:1040027", "hypermobility spectrum disorder", "component", "narrow",
             "exact label; the class has no xrefs in Mondo v2026-09-01"),
            ("MONDO:0001798", "hypermobility syndrome", "component", "narrow",
             "Mondo class carrying the exact ICD-10-CM xref M35.7 'Hypermobility syndrome' (and DOID:13781, whose OLS "
             "synonym is 'benign joint hypermobility'); its label is contained in the config alias 'joint hypermobility "
             "syndrome', and MONDO:1040027 lists 'benign joint hypermobility syndrome' as an exact synonym. Mondo does not "
             "place it under HSD (is_a arthropathy), so it is a separate component"),
        ],
        # M35.7 comes from the Mondo xref of MONDO:0001798, not from a CMS description match.
        "icd10cm": [("Q79.62", "narrow", None)],
    },
    "mcas": {
        "classes": [("MONDO:0100051", "idiopathic mast cell activation syndrome", "narrower", "narrow",
                     "config alias (exact label); is_a MONDO:0100004"),
                    ("MONDO:0033954", "monoclonal mast cell activation syndrome", "narrower", "narrow",
                     "is_a MONDO:0100004 (one of its three Mondo subclasses); exact ICD-10-CM xref D89.41"),
                    ("MONDO:0100006", "secondary mast cell activation syndrome", "narrower", "narrow",
                     "is_a MONDO:0100004 (one of its three Mondo subclasses); exact ICD-10-CM xref D89.43")],
        # D89.41 / D89.42 / D89.43 come from Mondo xrefs of the three MCAS subclasses.
        "icd10cm": [("D89.40", "close", None), ("D89.49", "narrow", None)],
        "note": "D89.41, D89.42 and D89.43 come from Mondo xrefs of the monoclonal, idiopathic and secondary MCAS "
                "subclasses. D89.44 'Hereditary alpha tryptasemia' sits under D89.4 but has no Mondo MCAS xref and its "
                "description does not match a condition name; not added.",
    },
    "post_infectious_syndrome": {
        "primary": ("MONDO:0021669", "post-infectious disorder",
                    "The lexical match MONDO:0021670 'post-infectious syndrome' is not used: in Mondo v2026-09-01 its "
                    "subclasses are congenital/acute infection syndromes (Zika virus congenital syndrome, TORCH syndrome, "
                    "KSHV inflammatory cytokine syndrome) and neither long COVID-19 nor PTLDS is classified under it. "
                    "MONDO:0021669 'post-infectious disorder' (synonym 'sequela of infectious disorder') is its parent and the "
                    "Mondo ancestor of both long COVID-19 (via post-viral disorder) and PTLDS (via post-bacterial disorder); "
                    "chosen manually as the closest grouping. 'PAIS' is an acronym homonym (partial androgen insensitivity "
                    "syndrome) and is rejected."),
        "relation": ("broad", "Mondo class covers all sequelae of infection (e.g. sequelae of tuberculosis), broader than "
                              "post-acute infection syndromes"),
        "classes": [("MONDO:0021674", "post-viral disorder", "narrower", "narrow",
                     "is_a MONDO:0021669; parent of MONDO:0100320 post-COVID-19 disorder"),
                    ("MONDO:0021670", "post-infectious syndrome", "related_distinct", "related",
                     "exact label match (plural-insensitive) to the config name, but its Mondo subclasses are congenital/"
                     "acute infection syndromes, not post-acute infection syndromes; kept for transparency, not used")],
        "icd10cm": [
            ("G93.31", "narrow", "CMS 'Postviral fatigue syndrome' names a post-infectious (post-viral) fatigue syndrome; "
                                 "lexical form differs from the config search term 'post-infectious fatigue'"),
            ("G93.39", "narrow", "CMS 'Other post infection and related fatigue syndromes' names post-infection fatigue "
                                 "syndromes; lexical form differs from 'post-infectious fatigue'"),
        ],
    },
    "endometriosis": {
        "classes": [("MONDO:0010888", "adenomyosis", "related_distinct", "related",
                     "config alias 'endometriosis of uterus' is an EXACT synonym of adenomyosis in Mondo; the config "
                     "marks adenomyosis as related but distinct")],
    },
}

ROLE_IN_REGISTRY = {"primary", "component", "narrower"}
ROLES_FOR_PHENOTYPES = {"primary", "component", "narrower"}

# Wearable phenotype axes. The HPO id is never typed in: it is resolved by OLS4 HPO search
# (exact label/synonym) from `queries`, in order. Signal columns are curated statements of what a
# wearable could observe (candidate measurable phenotype), not evidence that it does.
PHENOTYPE_AXES = [
    dict(axis_id="fatigue", axis_label="Fatigue", queries=["Fatigue"], accel=True, hr=False,
         signals=["accelerometry: daily activity volume", "accelerometry: activity fragmentation"],
         observability="indirect: fatigue is a symptom; wearables observe behavioural correlates only"),
    dict(axis_id="chronic_fatigue", axis_label="Chronic fatigue", queries=["Chronic fatigue"], accel=True, hr=False,
         signals=["accelerometry: sustained low daily activity over weeks", "accelerometry: day-to-day variability"],
         observability="indirect: requires symptom report to confirm"),
    dict(axis_id="exercise_intolerance", axis_label="Exercise intolerance", queries=["Exercise intolerance"], accel=True, hr=True,
         signals=["heart rate response to activity (PPG/ECG + accelerometry)", "accelerometry: peak and sustained activity bouts"],
         observability="partial: HR-activity coupling observable where HR is recorded"),
    dict(axis_id="post_exertional_malaise", axis_label="Post-exertional malaise",
         queries=["Post-exertional malaise", "Postexertional symptom exacerbation"], accel=True, hr=True,
         signals=["accelerometry: activity drop in the 1-3 days after high-activity days", "resting heart rate change after exertion"],
         observability="indirect: delayed activity crash is a candidate signature; symptom report required"),
    dict(axis_id="orthostatic_intolerance", axis_label="Orthostatic intolerance", queries=["Orthostatic intolerance"],
         closest_query="orthostatic", accel=True, hr=True,
         signals=["heart rate change at posture transitions (PPG/ECG + accelerometer posture)"],
         observability="partial: needs HR plus posture; blood pressure not observable by consumer wearables"),
    dict(axis_id="postural_tachycardia", axis_label="Postural tachycardia", queries=["Postural tachycardia"], accel=True, hr=True,
         signals=["heart rate increase within 10 minutes of standing (PPG/ECG + accelerometer posture)"],
         observability="partial: free-living posture transitions are not a controlled stand test"),
    dict(axis_id="tachycardia", axis_label="Tachycardia", queries=["Tachycardia"], accel=False, hr=True,
         signals=["resting and daytime heart rate (PPG/ECG)"], observability="direct where HR is recorded"),
    dict(axis_id="decreased_hrv", axis_label="Decreased heart rate variability", queries=["Decreased heart rate variability"],
         accel=False, hr=True, signals=["HRV from ECG or PPG inter-beat intervals (e.g. nocturnal RMSSD)"],
         observability="direct where beat-to-beat data are recorded"),
    dict(axis_id="orthostatic_hypotension", axis_label="Orthostatic hypotension", queries=["Orthostatic hypotension"],
         accel=False, hr=False, signals=[],
         observability="not observable by standard wrist wearables (requires blood pressure)"),
    dict(axis_id="autonomic_dysfunction", axis_label="Autonomic dysfunction",
         queries=["Autonomic dysfunction", "Dysautonomia"], accel=True, hr=True,
         signals=["HRV", "heart rate response to posture and activity"], observability="partial"),
    dict(axis_id="sleep_disturbance", axis_label="Sleep disturbance", queries=["Sleep disturbance"], accel=True, hr=False,
         signals=["actigraphy: sleep window, sleep fragmentation", "nocturnal heart rate where recorded"],
         observability="partial: actigraphy sleep estimates are approximate"),
    dict(axis_id="reduced_physical_activity", axis_label="Reduced physical activity",
         queries=["Reduced physical activity", "Physical inactivity"], accel=True, hr=False,
         signals=["accelerometry: total activity, sedentary fraction, steps"], observability="direct"),
    dict(axis_id="abnormal_circadian_rhythm", axis_label="Abnormal circadian rhythm",
         queries=["Abnormal circadian rhythm", "Sleep-wake cycle disturbance"], accel=True, hr=False,
         signals=["accelerometry rest-activity rhythm: interdaily stability, intradaily variability, relative amplitude"],
         observability="direct for rest-activity rhythm; not a measure of the molecular clock"),
]

# --------------------------------------------------------------------------- logging of API calls
_API_LOG: list[dict] = []


def _api_get(url: str, params: dict | None = None, *, tag: str) -> tuple[int, dict | None]:
    r = request("GET", url, params=params, cache_errors=True, reject_html=True)
    _API_LOG.append({"tag": tag, "url": url, "params": params, "status": r.status,
                     "fetched_at": r.fetched_at, "from_cache": r.from_cache})
    if r.status == 404:
        return 404, None
    r.raise_for_status()
    return r.status, r.json()


# --------------------------------------------------------------------------- downloads


def fetch_all() -> dict[str, Path]:
    files = {k: download_file(u, SOURCE_MONDO, k) for k, u in MONDO_URLS.items()}
    # First 4 KB of the release's OWL file: the owl:versionIRI lives in its header (mondo.obo carries only data-version).
    files["mondo_owl_header"] = download_file(MONDO_OWL_VERSIONED, SOURCE_MONDO, MONDO_OWL_HEADER_FILE,
                                              headers={"Range": "bytes=0-4095"})
    files["monarch_api_version.json"] = download_file(f"{MONARCH}/version", SOURCE_MONDO, "monarch_api_version.json")
    for key, url in CMS_URLS.items():
        files[f"cms_{key}"] = download_file(url, SOURCE_CMS)
    for ont in OLS_ONTOLOGIES:
        files[f"ols_{ont}"] = download_file(f"{OLS}/ontologies/{ont}", SOURCE_OLS, f"ols4_ontology_{ont}.json")
    return files


def _retrieved(source_id: str, filename: str) -> str:
    return load_manifest(source_id)["files"][filename]["retrieved_at"]


# --------------------------------------------------------------------------- Mondo


_SYN_RE = re.compile(r'^synonym: "((?:[^"\\]|\\.)*)" (EXACT|BROAD|NARROW|RELATED)(?: ([A-Za-z_:0-9]+))? \[(.*?)\]')
_XREF_RE = re.compile(r"^xref: (\S+)(?: \{(.*)\})?")


def parse_obo(path: Path) -> tuple[dict, dict]:
    """Parse [Term] stanzas of an OBO file. Returns (header, terms)."""
    header: dict = {}
    terms: dict = {}
    cur: dict | None = None
    in_header = True
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_header = False
                if cur and cur.get("id"):
                    terms[cur["id"]] = cur
                cur = {"syn": [], "xref": [], "isa": [], "skos": [], "subset": [], "obsolete": False,
                       "replaced_by": [], "consider": []} if line == "[Term]" else None
                continue
            if in_header:
                k, _, v = line.partition(": ")
                if k in ("data-version", "format-version", "ontology", "date"):
                    header[k] = v
                continue
            if cur is None or not line:
                continue
            k, _, v = line.partition(": ")
            if k == "id":
                cur["id"] = v
            elif k == "name":
                cur["name"] = v
            elif k == "synonym":
                m = _SYN_RE.match(line)
                if m:
                    cur["syn"].append((m.group(1).replace('\\"', '"'), m.group(2)))
            elif k == "xref":
                m = _XREF_RE.match(line)
                if m:
                    cur["xref"].append((m.group(1), re.findall(r'source="([^"]+)"', m.group(2) or "")))
            elif k == "is_a":
                cur["isa"].append(v.split()[0])
            elif k == "is_obsolete":
                cur["obsolete"] = v.strip() == "true"
            elif k == "replaced_by":
                cur["replaced_by"].append(v.split()[0])
            elif k == "consider":
                cur["consider"].append(v.split()[0])
            elif k == "subset":
                cur["subset"].append(v.split()[0])
            elif k == "property_value" and v.startswith("skos:"):
                parts = v.split()
                if len(parts) >= 2:
                    cur["skos"].append((parts[0], parts[1]))
    if cur and cur.get("id"):
        terms[cur["id"]] = cur
    return header, terms


SYN_MATCH = {"EXACT": ("exact synonym", "exact", 1), "RELATED": ("related synonym", "close", 3),
             "NARROW": ("narrow synonym", "broad", 3), "BROAD": ("broad synonym", "narrow", 3)}


def build_mondo_index(terms: dict) -> dict:
    """compact/singular key -> [(mondo_id, match_type, mondo_text, relation, rank)] for non-obsolete MONDO classes."""
    idx: dict = defaultdict(list)
    for tid, t in terms.items():
        if not tid.startswith("MONDO:") or t["obsolete"] or "name" not in t:
            continue
        entries = [(t["name"], "exact label", "exact", 0)] + [(s, *SYN_MATCH[sc]) for s, sc in t["syn"]]
        for text, mtype, rel, rank in entries:
            keys = {("c", compact_key(text))}
            if not is_acronym(text):  # never de-pluralise acronyms ('POTS' is not 'POT')
                keys.add(("s", singular_key(text)))
            for key in keys:
                idx[key].append((tid, mtype, text, rel, rank))
    return idx


def repair_flow_split(items) -> list[str]:
    """Re-join list items that a YAML flow sequence split at a comma inside parentheses.

    configs/conditions.yaml writes `aliases: [endometriosis of uterus, adenomyosis (related, distinct)]`, which YAML
    parses as ['endometriosis of uterus', 'adenomyosis (related', 'distinct)']. Without repair, 'distinct)' becomes an
    alias and normalize_condition('distinct') 'matches' endometriosis.
    """
    out: list[str] = []
    buf: str | None = None
    for it in list(items or []):
        s = str(it)
        if buf is not None:
            buf = f"{buf}, {s}"
            if buf.count("(") <= buf.count(")"):
                out.append(buf)
                buf = None
        elif s.count("(") > s.count(")"):
            buf = s
        else:
            out.append(s)
    if buf is not None:
        out.append(buf)
    return out


def name_variants(c: dict) -> list[dict]:
    """Condition names in priority order: tier 1 preferred-name parts, tier 2 aliases, tier 3 search terms."""
    out: list[dict] = []
    seen: set = set()

    def add(text: str, tier: int, kind: str):
        text = text.strip()
        key = compact_key(text)
        if not key or key in seen:
            return
        seen.add(key)
        default_rel = "exact" if tier == 1 else "close"
        rel = NAME_RELATION.get((c["id"], text), default_rel)
        out.append({"text": text, "tier": tier, "kind": kind, "position": len(out),
                    "acronym": is_acronym(text), "name_relation": rel})

    pn = c["preferred_name"]
    base = strip_parenthetical(pn)
    add(base, 1, "preferred_name")
    for part in re.split(r"\s*/\s*", base) if " / " in base else [base]:
        add(part, 1, "preferred_name_part")
    if "/" in base and " / " not in base:
        for part in base.split("/"):
            if len(part.split()) >= 2:
                add(part, 1, "preferred_name_part")
    for inner in re.findall(r"\(([^)]*)\)", pn):
        if len(inner.split()) >= 2:
            add(inner, 1, "preferred_name_part")
    for a in c.get("aliases", []):
        add(strip_parenthetical(a), 2, "alias")
    for s in c.get("search_terms", []):
        add(s, 3, "search_term")
    return out


def match_mondo(variants: list[dict], index: dict, terms: dict) -> list[dict]:
    rows = []
    for v in variants:
        hits = {}
        keys = [(("c", compact_key(v["text"])), "")]
        if not v["acronym"]:
            keys.append((("s", singular_key(v["text"])), " (plural-insensitive)"))
        for key, how in keys:
            for tid, mtype, mtext, rel, rank in index.get(key, []):
                if tid not in hits or rank < hits[tid]["mondo_rank"]:
                    hits[tid] = {"mondo_id": tid, "mondo_label": terms[tid]["name"], "match_type": mtype + how,
                                 "mondo_text": mtext, "mondo_relation": rel, "mondo_rank": rank}
        for h in hits.values():
            rows.append({**v, **h})
    return rows


def compose(a: str, b: str) -> str:
    """Relation condition->target from condition->class (a) and class->target (b)."""
    if a == "exact":
        return b
    if b == "exact":
        return a
    if a == b:
        return a
    if {a, b} == {"close", "narrow"}:
        return "narrow"
    if {a, b} == {"close", "broad"}:
        return "broad"
    return "related"


def description_match(desc: str, names: list[str]) -> str:
    """Lexical agreement between a CMS description and a list of names (longest name tried first).

    Match if the description (without a trailing ', unspecified') contains a name as whole words, or a
    name contains the description core. Plural- and punctuation-insensitive. Returns '' if no match.
    """
    core = re.sub(r",\s*unspecified$", "", str(desc), flags=re.I)
    ck = f" {singular_key(core)} "
    cand = sorted({x for x in names if x and not is_acronym(x) and singular_key(x)}, key=lambda x: (-len(norm_text(x)), x))
    for n in cand:
        if f" {singular_key(n)} " == ck:
            return f"CMS description equals name '{n}'"
    for n in cand:  # longest (most specific) name first
        if f" {singular_key(n)} " in ck:
            return f"CMS description contains name '{n}'"
    for n in reversed(cand):  # shortest name that still contains the description
        if ck.strip() and ck in f" {singular_key(n)} ":
            return f"name '{n}' contains CMS description '{core}'"
    return ""


SKOS_TO_REL = {"skos:exactMatch": "exact", "skos:closeMatch": "close", "skos:broadMatch": "broad",
               "skos:narrowMatch": "narrow", "skos:relatedMatch": "related"}


def mondo_xref_relation(qualifiers: list[str], skos: str | None) -> tuple[str, str]:
    """Return (relation, status) for a Mondo xref. status 'current' or 'obsolete'."""
    q = {x for x in qualifiers if x.startswith("MONDO:")}
    ql = {x.lower() for x in q}
    if any("obsolete" in x for x in ql):
        return ("exact" if any("equivalent" in x for x in ql) else "related"), "obsolete"
    if skos:
        return SKOS_TO_REL.get(skos, "related"), "current"
    if "mondo:equivalentto" in ql:
        return "exact", "current"
    if "mondo:mondoisnarrowerthansource" in ql:
        return "broad", "current"
    if "mondo:mondoisbroaderthansource" in ql:
        return "narrow", "current"
    return "related", "current"


def mondo_xrefs(term: dict) -> list[dict]:
    skos = {}
    for pred, obj in term["skos"]:
        skos[canonical_curie(obj)] = pred
    out = {}
    for val, quals in term["xref"]:
        cur = canonical_curie(val)
        rel, status = mondo_xref_relation(quals, skos.get(cur))
        out[cur] = {"target_id": cur, "mondo_predicate": skos.get(cur) or "", "mondo_xref_qualifiers": ";".join(quals),
                    "class_relation": rel, "target_status": "current" if status == "current" else "obsolete (per Mondo xref qualifier)"}
    for cur, pred in skos.items():
        if cur not in out:
            out[cur] = {"target_id": cur, "mondo_predicate": pred, "mondo_xref_qualifiers": "",
                        "class_relation": SKOS_TO_REL.get(pred, "related"), "target_status": "current"}
    return list(out.values())


def owl_version_iri(path: Path) -> str | None:
    """owl:versionIRI from the (partial) RDF/XML header of an OWL file."""
    m = re.search(r'<owl:versionIRI\s+rdf:resource="([^"]+)"', path.read_text(encoding="utf-8", errors="replace"))
    return m.group(1) if m else None


def mondo_icd10cm_xref_index(terms: dict) -> dict[str, list[tuple[str, str]]]:
    """ICD10CM CURIE -> [(mondo_id, label)] over current MONDO classes (to keep curated CMS codes honest)."""
    rev: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for tid, t in terms.items():
        if not tid.startswith("MONDO:") or t["obsolete"]:
            continue
        for v, _ in t["xref"]:
            cur = canonical_curie(v)
            if cur.startswith("ICD10CM:"):
                rev[cur].append((tid, t.get("name", "")))
    return rev


def load_sssom(path: Path) -> dict:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    return {(r.subject_id, canonical_curie(r.object_id)): r.predicate_id for r in df.itertuples(index=False)}


# --------------------------------------------------------------------------- CMS


def parse_cms_order(zip_path: Path) -> pd.DataFrame:
    """CMS icd10cm_order_YYYY.txt (fixed width): order, code, header flag, short, long description."""
    with zipfile.ZipFile(zip_path) as zf:
        name = next(n for n in zf.namelist() if re.search(r"icd10cm_order_\d{4}\.txt$", n))
        text = zf.read(name).decode("utf-8")
    rows = []
    for line in text.splitlines():
        if len(line) < 16:
            continue
        code = line[6:13].strip()
        rows.append((code, line[14] == "1", line[16:76].strip(), line[77:].strip()))
    df = pd.DataFrame(rows, columns=["code_nodot", "billable", "short_description", "description"])
    df["code"] = df["code_nodot"].map(lambda c: c[:3] + ("." + c[3:] if len(c) > 3 else ""))
    df.attrs["member"] = name
    return df


# --------------------------------------------------------------------------- OLS


def ols_term(curie: str) -> dict | None:
    """Look a CURIE up in its OLS4 ontology; None if absent. Label, obsolescence and replacement."""
    prefix, local = curie.split(":", 1)
    ont = OLS_ONTOLOGY_FOR.get(prefix)
    if not ont or prefix not in IRI_BASE:
        return None
    return ols_term_in(ont, IRI_BASE[prefix] + local)


def ols_term_in(ont: str, iri: str) -> dict | None:
    status, d = _api_get(f"{OLS}/ontologies/{ont}/terms", {"iri": iri}, tag=f"ols_term:{ont}")
    if status == 404 or not d:
        return None
    ts = d.get("_embedded", {}).get("terms", [])
    if not ts:
        return None
    t = ts[0]
    rb = t.get("term_replaced_by")
    return {"label": t.get("label"), "obo_id": t.get("obo_id"), "short_form": t.get("short_form"),
            "is_obsolete": bool(t.get("is_obsolete")), "replaced_by": canonical_curie(_iri_to_curie(rb)) if rb else None,
            "synonyms": t.get("synonyms") or [], "is_defining_ontology": t.get("is_defining_ontology")}


def _iri_to_curie(iri: str) -> str:
    for p, base in IRI_BASE.items():
        if iri.startswith(base):
            return f"{p}:{iri[len(base):]}"
    return iri


def ols_search_exact(ont: str, text: str, rows: int = 25) -> list[dict]:
    """OLS4 search restricted to hits whose label or a synonym equals `text` (punctuation/plural-insensitive).

    OLS `exact=true` still returns loose hits (e.g. 'tempol' for 'migraine'), so results are post-filtered.
    """
    _, d = _api_get(f"{OLS}/search", {"q": text, "ontology": ont, "exact": "true", "rows": rows,
                                      "fieldList": "iri,obo_id,short_form,label,synonym,is_obsolete,ontology_name,is_defining_ontology"},
                    tag=f"ols_search:{ont}")
    keys = {compact_key(text)} if is_acronym(text) else {compact_key(text), singular_key(text)}
    out = []
    for doc in (d or {}).get("response", {}).get("docs", []):
        if doc.get("is_obsolete"):
            continue
        label = doc.get("label") or ""
        how = None
        if {compact_key(label), singular_key(label)} & keys:
            how = "label"
        else:
            for s in doc.get("synonym") or []:
                if {compact_key(s), singular_key(s)} & keys:
                    how = "synonym"
                    break
        if how:
            out.append({"obo_id": canonical_curie(doc.get("obo_id") or ""), "label": label, "matched_on": how,
                        "iri": doc.get("iri"), "short_form": doc.get("short_form")})
    return out


def ols_versions() -> dict:
    out = {}
    for ont in OLS_ONTOLOGIES:
        d = json.loads((raw_dir(SOURCE_OLS) / f"ols4_ontology_{ont}.json").read_text())
        cfg = d.get("config", {})
        out[ont] = {"version": cfg.get("version"), "version_iri": cfg.get("versionIri"), "loaded": d.get("loaded"),
                    "number_of_terms": d.get("numberOfTerms"), "status": d.get("status")}
    return out


# --------------------------------------------------------------------------- Monarch


def monarch_phenotypes(mondo_id: str) -> list[dict]:
    """All Monarch disease->phenotype associations for a class and its descendants (paged, 500/page)."""
    items, offset = [], 0
    while True:
        _, d = _api_get(f"{MONARCH}/association", {
            "subject": mondo_id, "category": "biolink:DiseaseToPhenotypicFeatureAssociation",
            "direct": "false", "limit": 500, "offset": offset}, tag="monarch_association")
        items.extend(d.get("items", []))
        offset += 500
        if offset >= (d.get("total") or 0):
            return items


# --------------------------------------------------------------------------- build


def _latest_mondo_tag() -> str | None:
    try:
        return get_json(MONDO_LATEST_API).get("tag_name")
    except Exception as exc:  # network trouble must not break a cached rebuild
        return f"unavailable ({type(exc).__name__})"


def build() -> dict:
    _API_LOG.clear()
    files = fetch_all()
    cfg = load_config("conditions")
    # Copy (load_config is cached) and repair aliases/search terms split by YAML at a comma inside parentheses.
    conditions, alias_repairs = [], []
    for c0 in cfg["conditions"]:
        c = dict(c0)
        for key in ("aliases", "search_terms"):
            fixed = repair_flow_split(c0.get(key, []))
            if fixed != list(c0.get(key, []) or []):
                alias_repairs.append({"canonical_condition_id": c0["id"], "field": key,
                                      "as_parsed": list(c0.get(key, [])), "repaired": fixed})
            c[key] = fixed
        conditions.append(c)

    # ---- Mondo
    header, terms = parse_obo(files["mondo.obo"])
    data_version = header.get("data-version", "")
    if MONDO_TAG.lstrip("v") not in data_version:
        raise RuntimeError(f"mondo.obo data-version {data_version!r} does not match pinned tag {MONDO_TAG}")
    mondo_version = f"Mondo {MONDO_TAG} ({data_version})"
    mondo_version_iri = owl_version_iri(files["mondo_owl_header"])
    if not mondo_version_iri or data_version not in mondo_version_iri:
        raise RuntimeError(f"Mondo OWL versionIRI {mondo_version_iri!r} does not match OBO data-version {data_version!r}")
    latest = _latest_mondo_tag()
    mondo_index = build_mondo_index(terms)
    mondo_icd_rev = mondo_icd10cm_xref_index(terms)
    sssom = load_sssom(files["mondo.sssom.tsv"])
    mondo_retrieved = _retrieved(SOURCE_MONDO, "mondo.obo")

    # ---- CMS
    cms = parse_cms_order(files["cms_fy2026_apr"])
    cms_oct = parse_cms_order(files["cms_fy2026_oct"])
    cms27 = parse_cms_order(files["cms_fy2027"])
    cms_by = cms.set_index("code_nodot")
    cms27_by = cms27.set_index("code_nodot")
    cms_retrieved = _retrieved(SOURCE_CMS, CMS_URLS["fy2026_apr"].rsplit("/", 1)[-1])
    apr_vs_oct_identical = cms[["code_nodot", "billable", "description"]].reset_index(drop=True).equals(
        cms_oct[["code_nodot", "billable", "description"]].reset_index(drop=True))
    billable_nodot = cms.loc[cms["billable"], "code_nodot"].tolist()

    def icd_info(code: str) -> dict:
        """Validate one ICD-10-CM code or block against CMS FY2026 (and FY2027)."""
        if "-" in code:
            a, b = code.split("-")
            cats = sorted(c for c in cms["code_nodot"] if len(c) == 3 and a <= c <= b)
            label = "ICD-10-CM block " + code + ": " + "; ".join(
                f"{c} {cms_by.at[c, 'description']}" for c in cats) if cats else UNKNOWN
            return {"code_kind": "block", "in_cms_fy2026": bool(cats), "cms_description": label,
                    "cms_billable": False, "n_billable_codes_under": int(sum(1 for c in billable_nodot if a <= c[:3] <= b)),
                    "in_cms_fy2027": bool(cats), "cms_fy2027_description": UNKNOWN, "cms_description_changed_fy2027": False,
                    "validation_note": f"ICD-10-CM block, not a code; FY2026 categories inside: {','.join(cats)}"}
        nd = icd_nodot(code)
        in26 = nd in cms_by.index
        in27 = nd in cms27_by.index
        d26 = cms_by.at[nd, "description"] if in26 else UNKNOWN
        d27 = cms27_by.at[nd, "description"] if in27 else UNKNOWN
        return {"code_kind": "code", "in_cms_fy2026": in26, "cms_description": d26,
                "cms_billable": bool(cms_by.at[nd, "billable"]) if in26 else False,
                "n_billable_codes_under": int(sum(1 for c in billable_nodot if c.startswith(nd))),
                "in_cms_fy2027": in27, "cms_fy2027_description": d27,
                "cms_description_changed_fy2027": bool(in26 and in27 and d26 != d27),
                "validation_note": "" if in26 else "not in CMS FY2026 code set"}

    # ---- OLS versions + term cache
    olsv = ols_versions()
    ols_retrieved = _retrieved(SOURCE_OLS, "ols4_ontology_efo.json")
    ols_version = "; ".join(f"{k} {v['version']}" for k, v in olsv.items())
    ols_cache: dict[str, dict | None] = {}

    def ols(curie: str) -> dict | None:
        if curie not in ols_cache:
            ols_cache[curie] = ols_term(curie)
        return ols_cache[curie]

    monarch_version = json.loads(files["monarch_api_version.json"].read_text())
    monarch_kg = monarch_version.get("monarch_kg_version")

    registry_rows, mapping_rows, pheno_rows, candidate_rows = [], [], [], []
    curation_checks = []
    non_hpo_objects: list[dict] = []
    non_hpo_seen: set = set()
    freq_ids: set = set()

    for c in conditions:
        cid = c["id"]
        cur = CURATION.get(cid, {})
        variants = name_variants(c)
        cands = match_mondo(variants, mondo_index, terms)
        full = [x for x in cands if not x["acronym"]]
        full_ids = {x["mondo_id"] for x in full}

        # ---- primary class
        if "primary" in cur:
            pid, plabel, prat = cur["primary"]
            if full:  # never override an automatic match silently
                ab = min(full, key=lambda x: (x["tier"], x["position"], x["mondo_rank"]))
                if ab["mondo_id"] != pid:
                    prat += (f" [overrides automatic best match {ab['mondo_id']} '{ab['mondo_label']}' "
                             f"({ab['match_type']} on '{ab['text']}')]")
            primary = {"mondo_id": pid, "match_type": "manual choice", "matched_text": "", "mondo_text": "",
                       "name_relation": "exact", "mondo_relation": "exact", "rationale": prat}
            expected = [(pid, plabel)]
        elif full:
            best = min(full, key=lambda x: (x["tier"], x["position"], x["mondo_rank"]))
            primary = {"mondo_id": best["mondo_id"], "match_type": best["match_type"], "matched_text": best["text"],
                       "mondo_text": best["mondo_text"], "name_relation": best["name_relation"],
                       "mondo_relation": best["mondo_relation"],
                       "rationale": f"best match on {best['kind']} '{best['text']}' (tier {best['tier']})"}
            expected = []
        else:
            raise RuntimeError(f"{cid}: no Mondo match and no curated manual choice")
        rel_override = cur.get("relation")
        primary_relation = rel_override[0] if rel_override else compose(primary["name_relation"], primary["mondo_relation"])
        if rel_override:
            primary["rationale"] += f"; relation curated as {rel_override[0]}: {rel_override[1]}"

        classes = [{"mondo_id": primary["mondo_id"], "role": "primary", "relation": primary_relation,
                    "match_type": primary["match_type"], "matched_text": primary["matched_text"],
                    "mondo_text": primary["mondo_text"], "rationale": primary["rationale"]}]
        for mid, mlabel, role, rel, rat in cur.get("classes", []):
            expected.append((mid, mlabel))
            auto = [x for x in full if x["mondo_id"] == mid]
            if auto:
                a = min(auto, key=lambda x: (x["tier"], x["position"], x["mondo_rank"]))
                mt, mtext, mm = a["match_type"], a["text"], a["mondo_text"]
            else:
                mt, mtext, mm = "manual choice", "", ""
            classes.append({"mondo_id": mid, "role": role, "relation": rel, "match_type": mt, "matched_text": mtext,
                            "mondo_text": mm, "rationale": rat})

        # verify curated ids against this Mondo release
        for mid, mlabel in expected:
            t = terms.get(mid)
            ok = bool(t) and not t["obsolete"] and t.get("name") == mlabel
            curation_checks.append({"condition": cid, "mondo_id": mid, "expected_label": mlabel,
                                    "found_label": (t or {}).get("name"), "obsolete": (t or {}).get("obsolete"), "ok": ok})
            if not ok:
                raise RuntimeError(f"curation drift for {cid}: {mid} expected {mlabel!r}, found {(t or {}).get('name')!r}")

        used = {x["mondo_id"] for x in classes}
        for x in cands:
            status = "used" if x["mondo_id"] in used else (
                "rejected: acronym-only match (homonym risk)" if x["acronym"] and x["mondo_id"] not in full_ids
                else "not used (not reviewed as part of this condition)")
            candidate_rows.append({"canonical_condition_id": cid, "name": x["text"], "name_kind": x["kind"], "tier": x["tier"],
                                   "acronym": x["acronym"], "mondo_id": x["mondo_id"], "mondo_label": x["mondo_label"],
                                   "match_type": x["match_type"], "mondo_matched_text": x["mondo_text"], "status": status})

        # ---- mapping rows: Mondo classes, their xrefs, EFO import check
        def add_map(**kw):
            base = {"canonical_condition_id": cid, "mondo_id": None, "mondo_role": None, "target_ontology": None,
                    "target_id": None, "target_label": None, "target_label_source": None, "target_status": "current",
                    "predicate_condition": None, "class_relation": None, "mondo_predicate": "", "mondo_xref_qualifiers": "",
                    "in_mondo_sssom": None, "mapping_source": None, "source_version": None, "match_type": None,
                    "matched_text": "", "rationale": "", "code_kind": None, "in_cms_fy2026": None, "cms_description": None,
                    "cms_billable": None, "n_billable_codes_under": None, "in_cms_fy2027": None,
                    "cms_fy2027_description": None, "cms_description_changed_fy2027": None, "validation_note": "",
                    "description_match": ""}
            base.update(kw)
            mapping_rows.append(base)

        for cl in classes:
            t = terms[cl["mondo_id"]]
            o = ols(cl["mondo_id"])
            add_map(mondo_id=cl["mondo_id"], mondo_role=cl["role"], target_ontology="MONDO", target_id=cl["mondo_id"],
                    target_label=t["name"], target_label_source=f"{mondo_version} (OLS4 mondo {olsv['mondo']['version']}: "
                    f"{'label agrees' if o and o['label'] == t['name'] else 'label differs or absent'})",
                    predicate_condition=cl["relation"], class_relation="self", mapping_source=f"{mondo_version} label/synonym match"
                    if cl["match_type"] != "manual choice" else f"{mondo_version} manual choice (curated)",
                    source_version=mondo_version, match_type=cl["match_type"], matched_text=cl["matched_text"],
                    rationale=cl["rationale"])
            # EFO imports many Mondo classes; the Mondo IRI is then the EFO-usable id (Open Targets, GWAS Catalog).
            e = ols_term_in("efo", IRI_BASE["MONDO"] + cl["mondo_id"].split(":")[1])
            if e and not e["is_obsolete"]:
                add_map(mondo_id=cl["mondo_id"], mondo_role=cl["role"], target_ontology="EFO", target_id=cl["mondo_id"],
                        target_label=e["label"], target_label_source=f"OLS4 efo {olsv['efo']['version']}",
                        predicate_condition=cl["relation"], class_relation="exact (Mondo class imported into EFO)",
                        mapping_source=f"OLS4 EFO {olsv['efo']['version']}: Mondo class present in EFO",
                        source_version=f"EFO {olsv['efo']['version']}", match_type=cl["match_type"])
            for xr in mondo_xrefs(t):
                ont = xr["target_id"].split(":")[0]
                if ont == "ICD10CM" and not icd_dotted(xr["target_id"].split(":", 1)[1]):
                    continue
                row = dict(mondo_id=cl["mondo_id"], mondo_role=cl["role"], target_ontology=ont, target_id=xr["target_id"],
                           predicate_condition=compose(cl["relation"], xr["class_relation"]), class_relation=xr["class_relation"],
                           mondo_predicate=xr["mondo_predicate"], mondo_xref_qualifiers=xr["mondo_xref_qualifiers"],
                           in_mondo_sssom=(cl["mondo_id"], xr["target_id"]) in sssom,
                           mapping_source=f"{mondo_version} xref", source_version=mondo_version,
                           match_type=cl["match_type"], target_status=xr["target_status"])
                if ont == "ICD10CM":
                    info = icd_info(xr["target_id"].split(":", 1)[1])
                    row.update(info, target_label=info["cms_description"], target_label_source=CMS_VERSION)
                    if not info["in_cms_fy2026"]:
                        row["target_status"] = "not in CMS FY2026"
                    elif info["code_kind"] == "code":
                        dm = description_match(info["cms_description"], [t["name"]] + [x for x, _ in t["syn"]])
                        row["description_match"] = (f"Mondo xref; {dm}" if dm else
                                                    "Mondo xref; CMS description does not lexically match the Mondo label/synonyms (review)")
                elif ont in ("EFO", "MESH", "DOID", "ORPHANET"):
                    o = ols(xr["target_id"])
                    if o is None:
                        row.update(target_label=UNKNOWN, target_label_source=f"OLS4 {OLS_ONTOLOGY_FOR[ont]}: not found",
                                   target_status="not found in OLS4")
                    else:
                        row.update(target_label=o["label"], target_label_source=f"OLS4 {OLS_ONTOLOGY_FOR[ont]} "
                                   f"{olsv[OLS_ONTOLOGY_FOR[ont]]['version']}")
                        if o["is_obsolete"]:
                            row["target_status"] = "obsolete in OLS4" + (f" (replaced by {o['replaced_by']})" if o["replaced_by"] else "")
                else:
                    # UMLS / SNOMED CT / NCIT / MedGen / ... ids are Mondo's assertion only; no label is invented.
                    row.update(target_label_source=f"not verified: {ont} ids are not looked up (target_label left empty)")
                add_map(**row)

        # ---- OLS exact-label searches for EFO-native and MeSH ids by condition name
        # Added even when a Mondo xref already reached the id: the name-level evidence can be more specific
        # (e.g. MeSH 'Postural Orthostatic Tachycardia Syndrome' is exact for POTS although Mondo links it only
        # to the NET-deficiency subtype). Only repeated hits from later (lower-tier) names are skipped.
        seen_search = set()
        # An id that Mondo equates with a curated related-but-distinct class (e.g. DOID:288 'endometriosis of uterus',
        # whose DOID synonyms are 'adenomyosis') must not re-enter the condition as narrow/exact through a name search.
        related_distinct_targets = {r["target_id"]: r["mondo_id"] for r in mapping_rows
                                    if r["canonical_condition_id"] == cid and r["mondo_role"] == "related_distinct"
                                    and r["class_relation"] == "exact"}
        for v in variants:
            if v["acronym"]:
                continue
            for ont, prefix in (("efo", "EFO"), ("mesh", "MESH"), ("doid", "DOID")):
                for hit in ols_search_exact(ont, v["text"]):
                    if not hit["obo_id"].startswith(prefix + ":"):
                        continue  # Mondo/HP classes inside EFO are handled through the Mondo classes above
                    if (prefix, hit["obo_id"]) in seen_search:
                        continue
                    seen_search.add((prefix, hit["obo_id"]))
                    pred = compose(v["name_relation"], "exact" if hit["matched_on"] == "label" else "close")
                    rationale = ""
                    if hit["obo_id"] in related_distinct_targets:
                        pred = "related"
                        rationale = (f"name hit, but Mondo equates {hit['obo_id']} with the related-but-distinct class "
                                     f"{related_distinct_targets[hit['obo_id']]}; predicate set to related")
                    add_map(target_ontology=prefix, target_id=hit["obo_id"], target_label=hit["label"],
                            target_label_source=f"OLS4 {ont} {olsv[ont]['version']}",
                            predicate_condition=pred,
                            class_relation="", mapping_source=f"OLS4 {ont} {olsv[ont]['version']} exact {hit['matched_on']} search",
                            source_version=f"{ont} {olsv[ont]['version']}", match_type=f"exact {hit['matched_on']}",
                            matched_text=v["text"], rationale=rationale)

        # ---- curated CMS codes (no Mondo xref)
        names = [v["text"] for v in variants if not v["acronym"]]
        for cl in classes:
            if cl["role"] in ROLE_IN_REGISTRY:
                t = terms[cl["mondo_id"]]
                names += [t["name"]] + [s for s, sc in t["syn"] if sc == "EXACT" and not is_acronym(s)]
        for code, rel, semantic in cur.get("icd10cm", []):
            info = icd_info(code)
            if not info["in_cms_fy2026"]:
                raise RuntimeError(f"{cid}: curated ICD-10-CM {code} is not in the CMS FY2026 file")
            if mondo_icd_rev.get(f"ICD10CM:{code}"):
                # The brief allows CMS-description additions only where Mondo has no xref for the code.
                raise RuntimeError(f"{cid}: curated ICD-10-CM {code} is a Mondo xref of "
                                   f"{mondo_icd_rev[f'ICD10CM:{code}']}; take it from Mondo instead")
            desc = info["cms_description"]
            dm = description_match(desc, names)
            if not dm:
                if not semantic:
                    raise RuntimeError(f"{cid}: CMS description of {code} ({desc!r}) does not match any condition name")
                dm = f"curated semantic match (non-lexical): {semantic}"
            add_map(target_ontology="ICD10CM", target_id=f"ICD10CM:{code}", target_label=desc,
                    target_label_source=CMS_VERSION, predicate_condition=rel, class_relation="",
                    mapping_source=CMS_CURATED_SOURCE, source_version=CMS_VERSION, match_type="CMS description match",
                    description_match=dm, **info)

        # ---- Monarch phenotypes
        pheno_classes = [cl for cl in classes if cl["role"] in ROLES_FOR_PHENOTYPES]
        direct_ids = {cl["mondo_id"] for cl in pheno_classes if cl["role"] in ("primary", "component")}
        narrower_ids = {cl["mondo_id"] for cl in pheno_classes if cl["role"] == "narrower"}
        scope_rank = {"direct": 0, "narrower_class": 1, "descendant": 2}
        seen_assoc = {}
        for cl in pheno_classes:
            for it in monarch_phenotypes(cl["mondo_id"]):
                aid = it["id"]
                if not str(it.get("object", "")).startswith("HP:"):
                    if (cid, aid) in non_hpo_seen:
                        continue
                    non_hpo_seen.add((cid, aid))
                    non_hpo_objects.append({"canonical_condition_id": cid, "subject": it["subject"],
                                            "object": it.get("object"), "object_label": it.get("object_label"),
                                            "primary_knowledge_source": it.get("primary_knowledge_source")})
                    continue
                scope = ("direct" if it["subject"] in direct_ids else
                         "narrower_class" if it["subject"] in narrower_ids else "descendant")
                prev = seen_assoc.get(aid)
                if prev and scope_rank[prev["association_scope"]] <= scope_rank[scope]:
                    continue
                if it.get("frequency_qualifier"):
                    freq_ids.add(it["frequency_qualifier"])
                via = next((x for x in pheno_classes if x["mondo_id"] == it["subject"]), cl)
                seen_assoc[aid] = {
                    "canonical_condition_id": cid, "via_mondo_id": via["mondo_id"], "via_mondo_role": via["role"],
                    "annotated_disease_id": it["subject"], "annotated_disease_label": it.get("subject_label"),
                    "association_scope": scope, "hpo_id": it["object"], "hpo_label": it.get("object_label"),
                    "predicate": it.get("predicate"), "negated": bool(it.get("negated")),
                    "frequency_qualifier": it.get("frequency_qualifier"), "has_percentage": it.get("has_percentage"),
                    "has_count": it.get("has_count"), "has_total": it.get("has_total"),
                    "onset_qualifier": it.get("onset_qualifier"), "onset_qualifier_label": it.get("onset_qualifier_label"),
                    "sex_qualifier": it.get("sex_qualifier"),
                    "evidence_codes": ";".join(it.get("has_evidence") or []),
                    "publications": ";".join(it.get("publications") or []),
                    "primary_knowledge_source": it.get("primary_knowledge_source"),
                    "aggregator_knowledge_source": ";".join(it.get("aggregator_knowledge_source") or []),
                    "original_subject": it.get("original_subject"), "monarch_association_id": aid,
                    "_closure": set(it.get("object_closure") or []),
                }
        pheno_rows.extend(seen_assoc.values())

        registry_rows.append({"cfg": c, "classes": classes, "primary": primary, "primary_relation": primary_relation,
                              "note": cur.get("note", "")})

    # ---- frequency labels via OLS HPO
    freq_label = {f: (ols(f) or {}).get("label", UNKNOWN) for f in sorted(freq_ids)}

    # ---- phenotype axes
    axis_rows = []
    for ax in PHENOTYPE_AXES:
        hpo_id = hpo_label = None
        match_type = "no HPO term found by OLS4 exact search"
        matched_query = ""
        for i, q in enumerate(ax["queries"]):
            hits = [h for h in ols_search_exact("hp", q) if h["obo_id"].startswith("HP:")]
            if hits:
                h = sorted(hits, key=lambda x: x["matched_on"] != "label")[0]
                hpo_id, matched_query = h["obo_id"], q
                match_type = f"exact {h['matched_on']}" + ("" if i == 0 else f" (alternate query '{q}')")
                break
        closest = []
        if hpo_id is None:
            _, d = _api_get(f"{OLS}/search", {"q": ax.get("closest_query", ax["queries"][0]), "ontology": "hp", "rows": 5,
                                              "fieldList": "obo_id,label"}, tag="ols_search:hp_loose")
            closest = [f"{doc['obo_id']} {doc['label']}" for doc in (d or {}).get("response", {}).get("docs", [])
                       if str(doc.get("obo_id", "")).startswith("HP:")][:3]
        verified = ols(hpo_id) if hpo_id else None
        if verified:
            hpo_label = verified["label"]
        axis_rows.append({
            "axis_id": ax["axis_id"], "axis_label": ax["axis_label"], "hpo_id": hpo_id, "hpo_label": hpo_label,
            "hpo_match_type": match_type, "hpo_query": matched_query or ax["queries"][0],
            "hpo_obsolete": verified["is_obsolete"] if verified else None,
            "closest_hpo_terms_if_unmatched": closest,
            "wearable_signals": ax["signals"], "requires_accelerometry": ax["accel"], "requires_heart_rate": ax["hr"],
            "observability": ax["observability"],
            "wearable_signal_basis": "curated statement of candidate observability (not evidence of validity)",
        })

    # axis membership of each phenotype association (term or HPO descendant)
    axis_hpo = {a["hpo_id"]: a["axis_id"] for a in axis_rows if a["hpo_id"]}
    for p in pheno_rows:
        p["phenotype_axes"] = sorted(axis_hpo[h] for h in axis_hpo if h in p["_closure"] or h == p["hpo_id"])
        p["frequency_label"] = freq_label.get(p["frequency_qualifier"]) if p["frequency_qualifier"] else None
    for a in axis_rows:
        d = sorted({p["canonical_condition_id"] for p in pheno_rows
                    if a["axis_id"] in p["phenotype_axes"] and p["association_scope"] == "direct" and not p["negated"]})
        dd = sorted({p["canonical_condition_id"] for p in pheno_rows
                     if a["axis_id"] in p["phenotype_axes"] and p["association_scope"] != "direct" and not p["negated"]} - set(d))
        a["conditions_annotated_direct"] = d
        a["conditions_annotated_subclass_only"] = dd

    # ---- assemble tables
    maps = pd.DataFrame(mapping_rows)
    maps["in_registry_list"] = (
        maps["predicate_condition"].isin(REGISTRY_PREDICATES)
        & (maps["target_status"] == "current")
        & ((maps["target_ontology"] != "ICD10CM") | (maps["in_cms_fy2026"] == True))  # noqa: E712
        & (maps["mondo_role"].isna() | maps["mondo_role"].isin(ROLE_IN_REGISTRY))
    )
    maps["mapping_id"] = (maps["canonical_condition_id"] + "|" + maps["mondo_id"].fillna("-") + "|" + maps["target_id"]
                          + "|" + maps["mapping_source"].fillna(""))
    dup = maps["mapping_id"].duplicated()
    if dup.any():
        raise RuntimeError(f"duplicate mapping rows: {maps.loc[dup, 'mapping_id'].tolist()[:5]}")
    pheno = pd.DataFrame(pheno_rows).drop(columns=["_closure"]) if pheno_rows else pd.DataFrame()

    reg_out = []
    for rr in registry_rows:
        c, classes = rr["cfg"], rr["classes"]
        cid = c["id"]
        m = maps[(maps["canonical_condition_id"] == cid)]
        inl = m[m["in_registry_list"]]

        def ids(ont, frame=inl):
            sub = frame[frame["target_ontology"] == ont].copy()
            sub["r"] = sub["predicate_condition"].map(PREDICATE_RANK)
            return list(dict.fromkeys(sub.sort_values("r")["target_id"]))

        icd = inl[inl["target_ontology"] == "ICD10CM"]
        icd_codes = list(dict.fromkeys(icd.loc[icd["code_kind"] == "code", "target_id"].str.split(":").str[1]))
        icd_blocks = list(dict.fromkeys(icd.loc[icd["code_kind"] == "block", "target_id"].str.split(":").str[1]))
        broad = m[(m["predicate_condition"] == "broad") & (m["target_status"] == "current")]
        ph = pheno[pheno["canonical_condition_id"] == cid] if len(pheno) else pd.DataFrame(columns=["association_scope", "negated", "hpo_id"])
        hpo_direct = sorted(set(ph.loc[(ph["association_scope"] == "direct") & (~ph["negated"]), "hpo_id"]))
        reg_classes = [cl for cl in classes if cl["role"] in ROLE_IN_REGISTRY]
        # lexicon for search: labels/exact synonyms of primary + component classes only (not narrower subtypes)
        lex_classes = [cl for cl in classes if cl["role"] in ("primary", "component")]
        mondo_labels = [terms[cl["mondo_id"]]["name"] for cl in lex_classes]
        mondo_syn = sorted({s for cl in lex_classes for s, sc in terms[cl["mondo_id"]]["syn"] if sc == "EXACT"})
        efo = ids("EFO")
        prov = {
            "versions": {"mondo": mondo_version, "mondo_version_iri": mondo_version_iri,
                         "ols4": {k: v["version"] for k, v in olsv.items()},
                         "cms_icd10cm": CMS_VERSION, "monarch_kg": monarch_kg},
            "mondo_classes": [{k: cl[k] for k in ("mondo_id", "role", "relation", "match_type", "matched_text", "rationale")}
                              | {"label": terms[cl["mondo_id"]]["name"]} for cl in classes],
            "ids": [{"id": r.target_id, "ontology": r.target_ontology, "predicate": r.predicate_condition,
                     "source": r.mapping_source, "version": r.source_version, "via": r.mondo_id,
                     "match_type": r.match_type, "status": r.target_status}
                    for r in inl.itertuples(index=False)],
            "note": rr["note"],
        }
        reg_out.append({
            "canonical_condition_id": cid, "preferred_name": c["preferred_name"], "aliases": list(c.get("aliases", [])),
            "search_terms": list(c.get("search_terms", [])),
            "icd10cm_codes": icd_codes, "icd10cm_blocks": icd_blocks,
            "mondo_ids": [cl["mondo_id"] for cl in reg_classes],
            "primary_mondo_id": rr["primary"]["mondo_id"], "primary_mondo_label": terms[rr["primary"]["mondo_id"]]["name"],
            "primary_mondo_match_type": rr["primary"]["match_type"], "primary_mondo_relation": rr["primary_relation"],
            "related_mondo_ids": [f"{cl['mondo_id']}|{cl['role']}" for cl in classes if cl["role"] not in ROLE_IN_REGISTRY],
            "mondo_labels": mondo_labels, "mondo_exact_synonyms": mondo_syn,
            "efo_ids": efo, "efo_short_forms": [e.replace(":", "_") for e in efo],
            "mesh_ids": ids("MESH"), "doid_ids": ids("DOID"), "umls_ids": ids("UMLS"), "snomedct_ids": ids("SNOMEDCT"),
            "orphanet_ids": ids("ORPHANET"), "ncit_ids": ids("NCIT"),
            "broader_mapping_ids": list(dict.fromkeys(broad["target_id"])),
            "id_predicates": json.dumps({t: p for t, p in inl.sort_values(
                "predicate_condition", key=lambda x: x.map(PREDICATE_RANK), ascending=False)[
                ["target_id", "predicate_condition"]].itertuples(index=False)}, sort_keys=True),
            "hpo_terms": hpo_direct, "n_hpo_direct": len(hpo_direct),
            "n_hpo_narrower_class": int(len(set(ph.loc[ph["association_scope"] == "narrower_class", "hpo_id"]))),
            "n_hpo_descendant": int(len(set(ph.loc[ph["association_scope"] == "descendant", "hpo_id"]))),
            "hpo_annotation_status": ("direct Monarch annotations present" if hpo_direct else
                                      f"no direct Monarch KG {monarch_kg} disease-phenotype annotation"),
            "parent_condition": c.get("parent"), "grouping_only": bool(c.get("grouping_only", False)),
            **{f: bool(c["flags"][f]) for f in FLAG_COLUMNS},
            "config_notes": c.get("notes", ""),
            "ontology_provenance": json.dumps(prov, sort_keys=True),
        })
    reg = pd.DataFrame(reg_out)

    source_all = f"configs/conditions.yaml + {mondo_version} + OLS4 ({ols_version}) + CMS {CMS_VERSION} + Monarch KG {monarch_kg}"
    reg_t = add_provenance(
        reg, data_layer="ontology", source_name="configs/conditions.yaml + Mondo + EMBL-EBI OLS4 + CMS ICD-10-CM + Monarch KG",
        source_version=source_all, retrieved_at=mondo_retrieved, evidence_type="ontology_mapping",
        source_record_id="canonical_condition_id", evidence_level=reg["primary_mondo_relation"],
        provenance_notes="names, aliases, parent and flags are curated in configs/conditions.yaml (curated_config); "
                         "all identifiers come from Mondo/OLS4/CMS/Monarch with per-id provenance in ontology_provenance")
    write_table(reg_t, "condition_registry", producer=PRODUCER,
                description="One row per target condition with Mondo/EFO/MeSH/ICD-10-CM/HPO ids and per-id provenance")

    maps_t = add_provenance(
        maps, data_layer="ontology", source_name=maps["mapping_source"].fillna("unknown"),
        source_version=maps["source_version"].fillna(""), retrieved_at=mondo_retrieved,
        evidence_type="ontology_mapping", source_record_id="mapping_id", evidence_level=maps["predicate_condition"],
        provenance_notes="predicate_condition = relation from the canonical condition to target_id (exact/close/narrow/broad/related)")
    maps_t["source_name"] = maps["mapping_source"].fillna("unknown")
    maps_t["source_version"] = maps["source_version"].fillna("")
    maps_t.loc[maps_t["target_ontology"] == "ICD10CM", "retrieved_at"] = cms_retrieved
    maps_t.loc[maps_t["mapping_source"].str.startswith("OLS4", na=False), "retrieved_at"] = ols_retrieved
    write_table(maps_t, "condition_ontology_mappings", producer=PRODUCER,
                description="Long table: condition x ontology id with predicate, source, version, match type, validation")

    if len(pheno):
        pheno["phenotype_row_id"] = pheno["canonical_condition_id"] + "|" + pheno["monarch_association_id"]
        pheno_t = add_provenance(
            pheno, data_layer="ontology", source_name="Monarch KG API v3 (disease-phenotype associations)",
            source_version=f"Monarch KG {monarch_kg}; API {monarch_version.get('monarch_api_version')}",
            retrieved_at=_retrieved(SOURCE_MONDO, "monarch_api_version.json"), evidence_type="ontology_mapping",
            source_record_id="phenotype_row_id", evidence_level=pheno["association_scope"],
            provenance_notes="association_scope: direct = on the primary/component class; narrower_class = on a curated "
                             "narrower class; descendant = on another subclass (often a rare Mendelian subtype)")
        write_table(pheno_t, "condition_phenotypes", producer=PRODUCER,
                    description="Condition x HPO term from Monarch KG with frequency, evidence, source and scope")

    axes = pd.DataFrame(axis_rows)
    axes_t = add_provenance(
        axes, data_layer="ontology", source_name="EMBL-EBI OLS4 HPO search + curated wearable-signal statements",
        source_version=f"HP {olsv['hp']['version']}", retrieved_at=ols_retrieved, evidence_type="ontology_mapping",
        source_record_id="axis_id", evidence_level=axes["hpo_match_type"],
        provenance_notes="HPO ids resolved by OLS4 exact search; wearable_signals are curated candidate observability")
    write_table(axes_t, "phenotype_axes", producer=PRODUCER,
                description="Wearable-observable phenotype axes with verified HPO ids and candidate signals")

    codes = cms.merge(cms27[["code_nodot", "billable", "description"]].rename(
                          columns={"description": "description_fy2027", "billable": "billable_fy2027"}),
                      on="code_nodot", how="outer", indicator=True)
    codes["in_fy2026"] = codes["_merge"].isin(["both", "left_only"])
    codes["in_fy2027"] = codes["_merge"].isin(["both", "right_only"])
    codes["code"] = codes["code_nodot"].map(lambda c: c[:3] + ("." + c[3:] if len(c) > 3 else ""))
    codes["billable"] = codes["billable"].astype("boolean")
    codes["billable_fy2027"] = codes["billable_fy2027"].astype("boolean")
    codes = codes.drop(columns=["_merge"])
    codes_t = add_provenance(
        codes, data_layer="ontology", source_name="CMS ICD-10-CM code descriptions in tabular order",
        source_version=f"{CMS_VERSION}; FY2027 (effective 2026-10-01) for change detection",
        retrieved_at=cms_retrieved, evidence_type="ontology_mapping", source_record_id="code_nodot",
        provenance_notes="billable = CMS header flag 1 (FY2026; NA for FY2027-only codes); headers (0) are categories")
    write_table(codes_t, "ontology_icd10cm_codes", producer=PRODUCER,
                description="CMS ICD-10-CM FY2026 (Apr 2026) code set with FY2027 presence/descriptions")

    # ---- raw-side audit artefacts
    raw_m = raw_dir(SOURCE_MONDO)
    pd.DataFrame(candidate_rows).to_csv(raw_m / "condition_mondo_candidates.tsv", sep="\t", index=False)
    pd.DataFrame(curation_checks).to_csv(raw_m / "curation_checks.tsv", sep="\t", index=False)
    with gzip.open(raw_m / "monarch_disease_phenotype_associations.jsonl.gz", "wt") as fh:
        for p in pheno_rows:
            fh.write(json.dumps({k: (sorted(v) if isinstance(v, set) else v) for k, v in p.items() if k != "_closure"}) + "\n")
    (raw_dir(SOURCE_OLS) / "ols4_api_calls.json").write_text(json.dumps(
        [x for x in _API_LOG if x["tag"].startswith("ols")], indent=1, default=str))
    (raw_m / "monarch_api_calls.json").write_text(json.dumps(
        [x for x in _API_LOG if x["tag"].startswith("monarch")], indent=1, default=str))

    stats = _stats(terms, header, sssom, cms, cms_oct, cms27, apr_vs_oct_identical, reg, maps, pheno, axes,
                   candidate_rows, olsv, monarch_version, latest)
    stats["phenotypes"]["non_hpo_objects_excluded"] = non_hpo_objects
    stats["mondo"]["version_iri"] = mondo_version_iri
    stats["config_alias_repairs"] = alias_repairs
    _write_registry_entries(stats, mondo_retrieved, ols_retrieved, cms_retrieved)
    _write_audits(stats, files)
    return stats


# --------------------------------------------------------------------------- stats, registry entries, audits


def _stats(terms, header, sssom, cms, cms_oct, cms27, apr_eq, reg, maps, pheno, axes, cands, olsv, monarch_version, latest) -> dict:
    mondo_terms = {k: v for k, v in terms.items() if k.startswith("MONDO:")}
    n_obs = sum(t["obsolete"] for t in mondo_terms.values())
    xref_prefix = Counter(canonical_curie(v).split(":")[0] for t in mondo_terms.values() if not t["obsolete"] for v, _ in t["xref"])
    icd_maps = maps[maps["target_ontology"] == "ICD10CM"]
    miss = {}
    for col in ["icd10cm_codes", "efo_ids", "mesh_ids", "doid_ids", "umls_ids", "snomedct_ids", "orphanet_ids", "hpo_terms"]:
        n_empty = int(reg[col].map(len).eq(0).sum())
        miss[col] = {"empty": n_empty, "pct": round(100 * n_empty / len(reg), 1),
                     "conditions_empty": reg.loc[reg[col].map(len).eq(0), "canonical_condition_id"].tolist()}
    return {
        "mondo": {"data_version": header.get("data-version"), "n_terms_all_prefixes": len(terms),
                  "n_mondo_classes": len(mondo_terms), "n_mondo_obsolete": n_obs,
                  "n_mondo_current": len(mondo_terms) - n_obs, "xref_prefix_counts_current": dict(xref_prefix.most_common(20)),
                  "sssom_rows": len(sssom), "sssom_predicates": dict(Counter(sssom.values())), "latest_release_tag": latest},
        "cms": {"fy2026_apr_rows": len(cms), "fy2026_apr_billable": int(cms["billable"].sum()),
                "fy2026_apr_headers": int((~cms["billable"]).sum()), "fy2026_oct_rows": len(cms_oct),
                "apr_equals_oct_content": bool(apr_eq), "fy2027_rows": len(cms27),
                "fy2027_new_codes": int(len(set(cms27["code_nodot"]) - set(cms["code_nodot"]))),
                "fy2027_removed_codes": int(len(set(cms["code_nodot"]) - set(cms27["code_nodot"]))),
                "member": cms.attrs.get("member")},
        "registry": {"n_conditions": len(reg), "missingness": miss,
                     "primary": reg[["canonical_condition_id", "primary_mondo_id", "primary_mondo_label",
                                     "primary_mondo_match_type", "primary_mondo_relation"]].to_dict("records")},
        "mappings": {"n_rows": len(maps), "by_ontology": maps["target_ontology"].value_counts().to_dict(),
                     "by_source": maps["mapping_source"].value_counts().to_dict(),
                     "by_predicate": maps["predicate_condition"].value_counts().to_dict(),
                     "n_in_registry_list": int(maps["in_registry_list"].sum()),
                     "non_current": maps.loc[maps["target_status"] != "current",
                                             ["canonical_condition_id", "target_id", "target_label", "target_status"]].to_dict("records"),
                     "icd_rows": icd_maps[["canonical_condition_id", "target_id", "target_label", "mapping_source",
                                           "predicate_condition", "in_cms_fy2026", "cms_billable", "in_cms_fy2027",
                                           "cms_description_changed_fy2027", "description_match", "validation_note",
                                           "mondo_id"]].to_dict("records"),
                     "icd_not_in_cms": int((icd_maps["in_cms_fy2026"] == False).sum()),  # noqa: E712
                     "n_icd_mondo_conditions": int(icd_maps.loc[icd_maps["mapping_source"].str.contains("xref")
                                                                & icd_maps["in_registry_list"], "canonical_condition_id"].nunique()),
                     "n_icd_curated": int((icd_maps["mapping_source"] == CMS_CURATED_SOURCE).sum()),
                     "n_efo_xref": int((maps["target_id"].str.startswith("EFO:") & maps["mapping_source"].str.contains("xref")).sum()),
                     "n_efo_xref_obsolete": int((maps["target_id"].str.startswith("EFO:") & maps["mapping_source"].str.contains("xref")
                                                 & maps["target_status"].str.startswith("obsolete")).sum()),
                     "efo_rows": maps.loc[maps["target_ontology"] == "EFO", ["canonical_condition_id", "target_id", "target_label",
                                                                          "target_status", "mapping_source"]].to_dict("records")},
        "phenotypes": {"n_rows": int(len(pheno)),
                       "by_condition_scope": (pheno.groupby(["canonical_condition_id", "association_scope"]).size()
                                              .unstack(fill_value=0).to_dict("index") if len(pheno) else {}),
                       "by_source": pheno["primary_knowledge_source"].value_counts().to_dict() if len(pheno) else {},
                       "frequency_missing": int(pheno["frequency_qualifier"].isna().sum()) if len(pheno) else 0,
                       "percentage_missing": int(pheno["has_percentage"].isna().sum()) if len(pheno) else 0,
                       "n_negated": int(pheno["negated"].sum()) if len(pheno) else 0,
                       "n_distinct_hpo": int(pheno["hpo_id"].nunique()) if len(pheno) else 0,
                       "direct_small": (pheno[(pheno["association_scope"] == "direct") & pheno["canonical_condition_id"].isin(
                           ["me_cfs", "migraine", "gastroparesis", "dysautonomia"])][
                           ["canonical_condition_id", "hpo_id", "hpo_label", "primary_knowledge_source"]].to_dict("records")
                           if len(pheno) else [])},
        "axes": axes[["axis_id", "hpo_id", "hpo_label", "hpo_match_type", "closest_hpo_terms_if_unmatched",
                      "conditions_annotated_direct", "conditions_annotated_subclass_only"]].to_dict("records"),
        "candidates": {"n": len(cands), "by_status": dict(Counter(c["status"] for c in cands)),
                       "rejected": [c for c in cands if c["status"].startswith("rejected")]},
        "ols": olsv, "monarch": monarch_version,
        "api_calls": {"ols": sum(1 for x in _API_LOG if x["tag"].startswith("ols")),
                      "monarch": sum(1 for x in _API_LOG if x["tag"].startswith("monarch")),
                      "non_200": [x for x in _API_LOG if x["status"] != 200][:50]},
    }


def _write_registry_entries(s: dict, mondo_ret: str, ols_ret: str, cms_ret: str) -> None:
    common = {"person_level": False, "geographic": False, "omics": False, "wearable": False,
              "participant_linkage": "not applicable (concept-level vocabulary)",
              "true_participant_linkage_across_modalities": False, "geographic_resolution": "none",
              "ingestion_module": "measure_it.ontology.build_registry", "data_layer": "ontology"}
    outputs = ["condition_registry", "condition_ontology_mappings", "condition_phenotypes", "phenotype_axes",
               "ontology_icd10cm_codes"]
    ph = s["phenotypes"]
    write_registry_entry({
        "source_id": SOURCE_MONDO, "name": "Mondo Disease Ontology (+ Monarch KG disease-phenotype associations)",
        "publisher": "Monarch Initiative", "landing_url": "https://mondo.monarchinitiative.org/",
        "access_urls": list(MONDO_URLS.values()) + [f"{MONDO_OWL_VERSIONED} (first 4 KB, for owl:versionIRI)"]
                       + MONDO_PURLS + [f"{MONARCH}/association", f"{MONARCH}/version"],
        "license": "CC BY 4.0 (Mondo); Monarch KG aggregates sources under their own licenses (HPO annotations: HPO license)",
        "access_conditions": "open download and open API, no registration",
        "retrieved_at": mondo_ret, "source_version": f"Mondo {MONDO_TAG} ({s['mondo']['data_version']}; versionIRI "
                                                     f"{s['mondo']['version_iri']}); Monarch KG "
                                                     f"{s['monarch'].get('monarch_kg_version')}",
        "update_date": "Mondo: monthly releases; Monarch KG: rebuilt roughly monthly",
        "unit_of_observation": "disease class / disease-phenotype association",
        "sample_size": {"mondo_current_classes": s["mondo"]["n_mondo_current"], "mondo_obsolete_classes": s["mondo"]["n_mondo_obsolete"],
                        "registry_conditions": s["registry"]["n_conditions"], "mapping_rows": s["mappings"]["n_rows"],
                        "phenotype_rows": ph["n_rows"], "phenotype_distinct_hpo": ph["n_distinct_hpo"]},
        "status": "ingested", "processed_outputs": outputs, "audit": f"data/raw/{SOURCE_MONDO}/DATA_AUDIT.md",
        "limitations": [
            "No direct Monarch KG HPO annotation for: " + ", ".join(s["registry"]["missingness"]["hpo_terms"]["conditions_empty"]),
            "Descendant-scope phenotypes come mostly from rare Mendelian subtypes (e.g. EDS subtypes) and are flagged",
            "POTS has no stand-alone Mondo class; it is an exact synonym of 'orthostatic intolerance' (treated as broader)",
            f"Mondo ICD-10-CM xrefs reach {s['mappings']['n_icd_mondo_conditions']} of {s['registry']['n_conditions']} conditions; "
            f"{s['mappings']['n_icd_curated']} codes were added from CMS FY2026 descriptions (curated)",
            f"{s['mappings']['n_efo_xref_obsolete']} of {s['mappings']['n_efo_xref']} EFO xrefs carried by Mondo for these "
            f"conditions are obsolete in EFO {s['ols']['efo']['version']} (replaced by Mondo IRIs); EFO-usable ids are Mondo IRIs",
        ], **common})
    write_registry_entry({
        "source_id": SOURCE_OLS, "name": "EMBL-EBI Ontology Lookup Service (EFO, HPO, MeSH, DOID, ORDO, Mondo)",
        "publisher": "EMBL-EBI", "landing_url": "https://www.ebi.ac.uk/ols4/",
        "access_urls": [f"{OLS}/search", f"{OLS}/ontologies/{{ontology}}/terms"] + [f"{OLS}/ontologies/{o}" for o in OLS_ONTOLOGIES],
        "license": "OLS service: EMBL-EBI terms of use; each ontology under its own license (EFO/HPO/Mondo/DOID open; "
                   "MeSH NLM terms; ORDO CC BY 4.0)",
        "access_conditions": "open API, no registration",
        "retrieved_at": ols_ret, "source_version": "; ".join(f"{k} {v['version']}" for k, v in s["ols"].items()),
        "update_date": "OLS reloads ontologies nightly; loaded dates in ols4_ontology_*.json",
        "unit_of_observation": "ontology term (label, synonyms, obsolescence)",
        "sample_size": {"api_calls_logged": s["api_calls"]["ols"], "phenotype_axes": len(s["axes"]),
                        "efo_mapping_rows": len(s["mappings"]["efo_rows"])},
        "status": "ingested", "processed_outputs": ["condition_ontology_mappings", "condition_registry", "phenotype_axes"],
        "audit": f"data/raw/{SOURCE_OLS}/DATA_AUDIT.md",
        "limitations": ["OLS exact search returns loose hits; results are post-filtered to exact label/synonym equality",
                        "MeSH in OLS is the 2025 biopragmatics conversion", "No HPO term exists for 'orthostatic intolerance'"],
        **common})
    write_registry_entry({
        "source_id": SOURCE_CMS, "name": "CMS ICD-10-CM code descriptions in tabular order",
        "publisher": "CMS / CDC NCHS", "landing_url": CMS_LANDING, "access_urls": list(CMS_URLS.values()),
        "license": "U.S. federal government work (public domain)", "access_conditions": "open download, no registration",
        "retrieved_at": cms_ret, "source_version": CMS_VERSION, "update_date": "annual (Oct 1) with April 1 updates",
        "unit_of_observation": "ICD-10-CM code or category header",
        "sample_size": {"fy2026_rows": s["cms"]["fy2026_apr_rows"], "fy2026_billable": s["cms"]["fy2026_apr_billable"],
                        "fy2026_headers": s["cms"]["fy2026_apr_headers"], "fy2027_rows": s["cms"]["fy2027_rows"]},
        "status": "ingested", "processed_outputs": ["ontology_icd10cm_codes", "condition_ontology_mappings"],
        "audit": f"data/raw/{SOURCE_CMS}/DATA_AUDIT.md",
        "limitations": ["Code set only; no counts or prevalence", "FY2027 codes take effect 2026-10-01"], **common})


def _md_table(rows: list[dict], cols: list[str]) -> str:
    if not rows:
        return "_none_\n"
    out = "| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n"
    for r in rows:
        out += "| " + " | ".join(str(r.get(c, "")).replace("|", "/") for c in cols) + " |\n"
    return out


def _manifest_lines(source_id: str) -> str:
    m = load_manifest(source_id)["files"]
    return "".join(f"- `{k}` — {v['bytes']:,} bytes, sha256 `{v['sha256'][:16]}…`, retrieved {v['retrieved_at']}  \n  <{v['url']}>\n"
                   for k, v in sorted(m.items()))


def _write_audits(s: dict, files: dict) -> None:
    today = utc_now_iso()
    miss = s["registry"]["missingness"]
    miss_rows = [{"field": k, "empty_n": v["empty"], "empty_pct": v["pct"], "conditions_empty": ", ".join(v["conditions_empty"])}
                 for k, v in miss.items()]
    ph = s["phenotypes"]
    scope_rows = [{"condition": k, **v} for k, v in ph["by_condition_scope"].items()]
    mondo_md = f"""# DATA AUDIT — Mondo Disease Ontology (+ Monarch KG disease-phenotype associations)

_Generated by `measure_it.ontology.build_registry` at {today}. Every number below is computed from the retrieved files/API responses._

| Field | Value |
|---|---|
| source_id | {SOURCE_MONDO} |
| Source (dataset/API name, exact files/endpoints) | Mondo release {MONDO_TAG}: `mondo.obo`, `mondo.sssom.tsv` (tag-pinned), `source-versions.tsv`; Monarch KG API v3 `/association` (category biolink:DiseaseToPhenotypicFeatureAssociation, direct=false) and `/version` |
| Publishing organization | Monarch Initiative |
| Retrieval date (UTC) | {_retrieved(SOURCE_MONDO, 'mondo.obo')} |
| Source version / release | Mondo {MONDO_TAG} (OBO data-version `{s['mondo']['data_version']}`; owl:versionIRI `{s['mondo']['version_iri']}` read from the release OWL header); latest GitHub release at build time: {s['mondo']['latest_release_tag']}; Monarch KG {s['monarch'].get('monarch_kg_version')} (API {s['monarch'].get('monarch_api_version')}) |
| Source update date / cadence | Mondo monthly; Monarch KG rebuilt about monthly |
| License / access conditions | Mondo CC BY 4.0; Monarch KG aggregates HPO annotations (OMIM/Orphanet-derived) under their licenses; open, no registration |
| Unit of observation | disease class; disease-phenotype association |
| Sample size (actual, as ingested) | {s['mondo']['n_mondo_classes']:,} MONDO classes in mondo.obo ({s['mondo']['n_mondo_current']:,} current, {s['mondo']['n_mondo_obsolete']:,} obsolete; {s['mondo']['n_terms_all_prefixes']:,} [Term] stanzas incl. imported classes); SSSOM {s['mondo']['sssom_rows']:,} rows {s['mondo']['sssom_predicates']}; {s['registry']['n_conditions']} registry conditions; {ph['n_rows']:,} condition-phenotype rows ({ph['n_distinct_hpo']:,} distinct HPO terms) |
| Geography (resolution, vintage) | none (concept-level) |
| Person-level? | no |
| Geographic? | no |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — vocabulary, no participants |

## Files / endpoints retrieved
{_manifest_lines(SOURCE_MONDO)}
- `condition_mondo_candidates.tsv` — every Mondo class matched by any condition name, with status (used / rejected / not used)
- `curation_checks.tsv` — each curated Mondo id checked against this release (exists, non-obsolete, label as recorded)
- `monarch_disease_phenotype_associations.jsonl.gz` — the Monarch association rows used
- `monarch_api_calls.json` — {s['api_calls']['monarch']} Monarch API calls (URL, params, HTTP status, fetch time, cache hit)

Monarch queries: `GET {MONARCH}/association?subject=<MONDO id>&category=biolink:DiseaseToPhenotypicFeatureAssociation&direct=false&limit=500&offset=k`
for each Mondo class with role primary/component/narrower.

## Key variables
- Mondo `[Term]`: id, name, synonym (EXACT/RELATED/NARROW/BROAD), xref with `source=` qualifiers (MONDO:equivalentTo = exact,
  MONDO:relatedTo = related, MONDO:equivalentObsolete / obsoleteEquivalent = obsolete target, mondoIsNarrowerThanSource = broad,
  mondoIsBroaderThanSource = narrow), `property_value: skos:*Match` (used as the predicate where present), is_a, is_obsolete.
- Current-class xref prefix counts: {s['mondo']['xref_prefix_counts_current']}
- Monarch association: subject (annotated disease), object (HPO id), negated, frequency_qualifier (HPO frequency term; label
  resolved via OLS4), has_percentage / has_count / has_total, onset, has_evidence (ECO), publications, primary_knowledge_source.

## Matching results (primary Mondo class per condition)
{_md_table(s['registry']['primary'], ['canonical_condition_id', 'primary_mondo_id', 'primary_mondo_label', 'primary_mondo_match_type', 'primary_mondo_relation'])}
Candidate matches: {s['candidates']['n']} ({s['candidates']['by_status']}). Rejected acronym-only homonyms:
{_md_table(s['candidates']['rejected'], ['canonical_condition_id', 'name', 'mondo_id', 'mondo_label', 'match_type'])}

## Missingness
Registry list fields with no identifier (of {s['registry']['n_conditions']} conditions):
{_md_table(miss_rows, ['field', 'empty_n', 'empty_pct', 'conditions_empty'])}
Condition-phenotype rows by condition and scope (direct = annotated on the primary/component Mondo class; narrower_class =
annotated on a curated narrower class such as the monogenic NET-deficiency form of POTS; descendant = on any other subclass):
{_md_table(scope_rows, ['condition'] + sorted({k for r in scope_rows for k in r if k != 'condition'}))}
Conditions absent from this table have zero Monarch disease-phenotype associations. Frequency qualifier missing on
{ph['frequency_missing']:,} of {ph['n_rows']:,} rows; has_percentage missing on {ph['percentage_missing']:,}; negated rows: {ph['n_negated']}.
Primary knowledge sources: {ph['by_source']}.

Monarch returned {len(ph.get('non_hpo_objects_excluded', []))} disease-to-"phenotype" rows whose object is not an HPO term
(Mondo feature links such as PTLDS -> 'Lyme disease'); they are excluded:
{_md_table(ph.get('non_hpo_objects_excluded', []), ['canonical_condition_id', 'subject', 'object', 'object_label', 'primary_knowledge_source'])}

## Linkage strategy
Conditions join to Mondo classes by label/synonym match (plus documented manual choices); every other id is reached through
Mondo xrefs, OLS4 lookups, or CMS description matches and is recorded with its source in `condition_ontology_mappings`.
These are concept-level links; they never link people.

## Limitations and caveats
- HPO/Monarch disease-phenotype coverage is concentrated in rare Mendelian diseases. Several target conditions have no
  direct annotation at all (see table above). Direct annotations on ME/CFS, migraine, gastroparesis and dysautonomia:
  {'; '.join(f"{r['canonical_condition_id']}: {r['hpo_id']} {r['hpo_label']} ({r['primary_knowledge_source']})" for r in ph['direct_small'])}.
- Descendant-scope rows (e.g. {ph['by_condition_scope'].get('eds_hsd', {}).get('descendant', 0)} for EDS/HSD) mostly describe
  rare subtypes (classical, vascular, kyphoscoliotic EDS ...) and should not be read as the phenotype of the common condition.
- POTS: Mondo {MONDO_TAG} has no stand-alone class; 'postural orthostatic tachycardia syndrome' is an EXACT synonym of
  MONDO:0001315 'orthostatic intolerance', whose other synonyms (neurocirculatory asthenia, soldiers heart) and MeSH xref
  (Neurocirculatory Asthenia) are not POTS-specific. The class is curated as broader than POTS.
- ME/CFS: the Mondo label is spelled 'encephalomeyelitis' in this release, and the class sits under 'immunodeficiency disease'
  and 'hereditary neurological disease' in Mondo's hierarchy (affects Monarch closure-based queries).
- Mondo ICD-10-CM xrefs come from an older ICD-10-CM mirror (source-versions.tsv lists 2024ab) and are sparse; all were
  validated against CMS FY2026 (see cms_icd10cm audit).
- configs/conditions.yaml list items split by YAML at a comma inside parentheses were re-joined before matching
  (without this, 'distinct)' became an endometriosis alias): {'; '.join(f"{r['canonical_condition_id']}.{r['field']}: {r['as_parsed']} -> {r['repaired']}" for r in s.get('config_alias_repairs', [])) or 'none'}.
- Post-infectious grouping: the lexical match MONDO:0021670 'post-infectious syndrome' groups congenital/acute infection
  syndromes and contains neither long COVID-19 nor PTLDS; MONDO:0021669 'post-infectious disorder' is a documented manual
  choice (rationale in condition_ontology_mappings and condition_mondo_candidates.tsv).

## Processed outputs
condition_registry ({s['registry']['n_conditions']} rows), condition_ontology_mappings ({s['mappings']['n_rows']} rows),
condition_phenotypes ({ph['n_rows']} rows), phenotype_axes ({len(s['axes'])} rows), ontology_icd10cm_codes.

## Reproduce
`uv run python -m measure_it.ontology.build_registry`
"""
    (raw_dir(SOURCE_MONDO) / "DATA_AUDIT.md").write_text(mondo_md)

    efo_rows = s["mappings"]["efo_rows"]
    ols_md = f"""# DATA AUDIT — EMBL-EBI OLS4 (EFO, HPO, MeSH, DOID, ORDO, Mondo)

_Generated by `measure_it.ontology.build_registry` at {today}._

| Field | Value |
|---|---|
| source_id | {SOURCE_OLS} |
| Source (dataset/API name, exact files/endpoints) | OLS4 REST API: `{OLS}/ontologies/{{ont}}` (metadata), `{OLS}/ontologies/{{ont}}/terms?iri=` (term lookup), `{OLS}/search?q=&ontology=&exact=true` (search) |
| Publishing organization | EMBL-EBI |
| Retrieval date (UTC) | {_retrieved(SOURCE_OLS, 'ols4_ontology_efo.json')} |
| Source version / release | {'; '.join(f"{k} {v['version']} (loaded {v['loaded']})" for k, v in s['ols'].items())} |
| Source update date / cadence | OLS reloads ontologies nightly |
| License / access conditions | EMBL-EBI terms of use; ontologies under their own licenses; open API, no registration |
| Unit of observation | ontology term |
| Sample size (actual, as ingested) | {s['api_calls']['ols']} OLS API calls (logged); {len(efo_rows)} EFO mapping rows; {len(s['axes'])} phenotype axes; ontology sizes: {', '.join(f"{k} {v['number_of_terms']:,} terms" for k, v in s['ols'].items())} |
| Geography (resolution, vintage) | none |
| Person-level? | no |
| Geographic? | no |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — vocabulary |

## Files / endpoints retrieved
{_manifest_lines(SOURCE_OLS)}
- `ols4_api_calls.json` — every OLS call (URL, params, HTTP status, fetch time, cache hit). Non-200 responses: {len([x for x in s['api_calls']['non_200'] if 'ols' in x['tag']])} (404 = term not present in that ontology, e.g. a Mondo class not imported into EFO).

## Key variables
- term lookup: label, is_obsolete, term_replaced_by, synonyms → target_label / target_status in condition_ontology_mappings.
- search: label/synonym equality (post-filtered; OLS `exact=true` alone returns loose hits such as 'tempol' for 'migraine').

## EFO resolution
{_md_table(efo_rows, ['canonical_condition_id', 'target_id', 'target_label', 'target_status', 'mapping_source'])}
EFO identifiers for diseases are often Mondo IRIs imported into EFO (EFO obsoleted several native disease terms in favour of
Mondo; e.g. Mondo's xref EFO:0004540 for ME/CFS is obsolete in EFO, replaced by MONDO:0005404). `efo_short_forms` in
condition_registry gives the form used by Open Targets / GWAS Catalog (e.g. MONDO_0005404).

## Phenotype axes (HPO via OLS4)
{_md_table(s['axes'], ['axis_id', 'hpo_id', 'hpo_label', 'hpo_match_type', 'closest_hpo_terms_if_unmatched', 'conditions_annotated_direct', 'conditions_annotated_subclass_only'])}

## Missingness
Obsolete / not-found targets encountered: {len(s['mappings']['non_current'])} (excluded from registry id lists):
{_md_table(s['mappings']['non_current'], ['canonical_condition_id', 'target_id', 'target_label', 'target_status'])}
Axes without an HPO term: {sum(1 for a in s['axes'] if not isinstance(a['hpo_id'], str))} ({', '.join(a['axis_id'] for a in s['axes'] if not isinstance(a['hpo_id'], str)) or 'none'}).

## Linkage strategy
OLS verifies and enriches ids reached from Mondo; OLS search hits are added only when the label or a synonym equals a
condition name. Predicates: label hit on a preferred-name part = exact; synonym hit = close; hits on names curated as narrower
(e.g. 'orthostatic intolerance' for dysautonomia) inherit that relation.

## Limitations and caveats
- OLS is a live service; versions above are what was loaded on the retrieval date. HTTP responses are cached under data/_http_cache.
- MeSH in OLS is the biopragmatics 2025 conversion (obo_id prefix `mesh:`; normalized to `MESH:`).
- HPO has no 'orthostatic intolerance' term; 'reduced physical activity' and 'abnormal circadian rhythm' resolve only through
  alternate phrasings (see match types).
- Wearable signal columns in phenotype_axes are curated statements of candidate observability, not evidence.

## Processed outputs
phenotype_axes ({len(s['axes'])} rows); EFO/MeSH/DOID/ORDO rows in condition_ontology_mappings.

## Reproduce
`uv run python -m measure_it.ontology.build_registry`
"""
    (raw_dir(SOURCE_OLS) / "DATA_AUDIT.md").write_text(ols_md)

    c = s["cms"]
    icd_rows = s["mappings"]["icd_rows"]
    cms_md = f"""# DATA AUDIT — CMS ICD-10-CM code descriptions (FY2026)

_Generated by `measure_it.ontology.build_registry` at {today}._

| Field | Value |
|---|---|
| source_id | {SOURCE_CMS} |
| Source (dataset/API name, exact files/endpoints) | CMS "Code Descriptions in Tabular Order" ZIPs; member `{c['member']}` (fixed width: order, code, header flag, short, long description) |
| Publishing organization | CMS (code set maintained by CDC NCHS) |
| Retrieval date (UTC) | {_retrieved(SOURCE_CMS, CMS_URLS['fy2026_apr'].rsplit('/', 1)[-1])} |
| Source version / release | {CMS_VERSION}; FY2026 Oct-2025 file and FY2027 (effective 2026-10-01) file retrieved for change detection |
| Source update date / cadence | annual Oct 1 release, April 1 updates |
| License / access conditions | U.S. federal government work (public domain); open download |
| Unit of observation | ICD-10-CM code (billable) or category header |
| Sample size (actual, as ingested) | FY2026 Apr: {c['fy2026_apr_rows']:,} rows = {c['fy2026_apr_billable']:,} billable codes + {c['fy2026_apr_headers']:,} headers; FY2026 Oct: {c['fy2026_oct_rows']:,} rows (content identical to April: {c['apr_equals_oct_content']}); FY2027: {c['fy2027_rows']:,} rows ({c['fy2027_new_codes']} new, {c['fy2027_removed_codes']} removed vs FY2026) |
| Geography (resolution, vintage) | none |
| Person-level? | no |
| Geographic? | no |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no |

## Files / endpoints retrieved
Landing page: <{CMS_LANDING}>
{_manifest_lines(SOURCE_CMS)}

## Key variables
code (dotted), code_nodot, billable (header flag 1), short_description, description (long), description_fy2027, in_fy2026, in_fy2027.

## Validation of every ICD-10-CM id in the ontology layer
{_md_table(icd_rows, ['canonical_condition_id', 'target_id', 'target_label', 'mapping_source', 'predicate_condition', 'in_cms_fy2026', 'cms_billable', 'in_cms_fy2027', 'cms_description_changed_fy2027', 'description_match', 'validation_note'])}
ICD-10-CM ids not present in FY2026: {s['mappings']['icd_not_in_cms']}.

## Missingness
Conditions without any validated ICD-10-CM code: {', '.join(s['registry']['missingness']['icd10cm_codes']['conditions_empty']) or 'none'}.

## Linkage strategy
Mondo ICD10CM xrefs and curated additions are matched to the CMS file on the undotted code. Curated additions
(mapping_source '{CMS_CURATED_SOURCE}') are accepted only if no current Mondo class carries the code as an xref (the build
fails otherwise) and the CMS long description contains a condition name (or a condition name contains the description
core) — or, where marked 'curated semantic match (non-lexical)', with a stated rationale.

## Limitations and caveats
- Codes describe billing categories, not validated case definitions; code presence is not prevalence.
- Header categories (e.g. G43, N80, K58) include all subcodes; `n_billable_codes_under` gives the count.
- Mondo's adenomyosis xref ICD10CM:N80.0 is now the header 'Endometriosis of uterus'; FY2026 codes adenomyosis as N80.03.
- FY2027 changes take effect 2026-10-01; `in_cms_fy2027` and `cms_description_changed_fy2027` flag them.

## Processed outputs
ontology_icd10cm_codes ({c['fy2026_apr_rows']:,} FY2026 rows plus FY2027-only rows); ICD rows in condition_ontology_mappings ({len(icd_rows)}).

## Reproduce
`uv run python -m measure_it.ontology.build_registry`
"""
    (raw_dir(SOURCE_CMS) / "DATA_AUDIT.md").write_text(cms_md)


def run() -> dict:
    return build()


if __name__ == "__main__":
    st = build()
    print(json.dumps({"registry": st["registry"]["primary"], "mappings": st["mappings"]["by_ontology"],
                      "phenotypes": st["phenotypes"]["n_rows"], "axes": [(a["axis_id"], a["hpo_id"]) for a in st["axes"]]},
                     indent=1, default=str))
