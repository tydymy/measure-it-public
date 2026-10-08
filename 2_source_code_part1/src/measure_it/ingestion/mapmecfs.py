"""NIH intramural post-infectious ME/CFS deep phenotyping (Walitt et al., Nat Commun 2024).

Public deposits, participant-identifier linkage audit, linked omics subset and
condition-level molecular evidence for source_id ``mapmecfs_nih_pi_mecfs``.

What the module does
--------------------
1. Acquire every deposit named in the paper's Data Availability statement that
   is openly downloadable without an account: Nature supplementary + source
   data, the authors' GitHub repository, GEO GSE251872 / GSE245661 / GSE251790 /
   GSE254030 (processed files), SRA run + BioSample metadata for BioProject
   PRJNA954397, and Pennsieve dataset 356 metadata. mapMECFS data files answer
   anonymous requests with HTTP 403 "Login Required", so mapMECFS is catalogued
   from its public CKAN metadata only. No account is created, no login is
   attempted and the CKAN datastore API is not used to route around the gate.
2. Extract the participant identifiers of every modality actually obtained,
   classify each identifier namespace, and count exact overlaps.
3. Build the linked subset only where identifiers are shared (NIH study numbers
   1xx = healthy volunteer, 3xx = PI-ME/CFS across the GEO omics deposits and
   the ID-bearing rows of the CSF-metabolomics source data).
4. Write condition-level molecular evidence from Supplementary Data 14-21.
5. Write DATA_AUDIT.md, docs/MAP_MECFS_LINKAGE_AUDIT.md and registry_entry.yaml.

Wearable/actigraphy/HRV/orthostatic/CPET data are NOT linked to omics: in the
openly downloadable files they carry group labels only (HV / PI-ME/CFS), and
the per-participant files exist only behind the mapMECFS login. No
participant_wearable_features__mapmecfs table is written for that reason.

Reproduce: ``uv run python -m measure_it.ingestion.mapmecfs``
"""
from __future__ import annotations

import argparse
import gzip
import html
import io
import itertools
import json
import re
import tarfile
import warnings
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit

import numpy as np
import pandas as pd
import requests

from .. import download as _dl
from ..config import DOCS, PROCESSED, UNKNOWN, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..http import get, get_json, post_json
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

SOURCE_ID = "mapmecfs_nih_pi_mecfs"
DATASET_ID = "mapmecfs"
PRODUCER = "measure_it.ingestion.mapmecfs"
CONDITION_ID = "me_cfs"
CONDITION_LABEL = "post-infectious ME/CFS (PI-ME/CFS), NIH Clinical Center deep-phenotyping cohort"
LINKAGE_AUDIT_PATH = DOCS / "MAP_MECFS_LINKAGE_AUDIT.md"

PAPER = {
    "citation": "Walitt B, et al. Deep phenotyping of post-infectious myalgic encephalomyelitis/"
                "chronic fatigue syndrome. Nat Commun 15, 907 (2024)",
    "pmid": "38383456",
    "pmcid": "PMC10881493",
    "doi": "10.1038/s41467-024-45107-3",
    "url": "https://www.nature.com/articles/s41467-024-45107-3",
    "trial": "NCT02669212",
}
STUDY_POPULATION = ("NIH Clinical Center protocol 16-N-0058 (NCT02669212): adults with post-infectious "
                    "ME/CFS and healthy volunteers, Bethesda MD")

# --------------------------------------------------------------------------- endpoints
MAPMECFS_HOME = "https://www.mapmecfs.org/"
MAPMECFS_GROUP = "post-infectious-mecfs-at-the-nih"
MAPMECFS_GROUP_URL = f"https://www.mapmecfs.org/group/{MAPMECFS_GROUP}"
CKAN = "https://www.mapmecfs.org/api/3/action/"
ESM_URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-45107-3/"
           "MediaObjects/41467_2024_45107_MOESM{n}_ESM.{ext}")
ESM = {
    1: ("pdf", "Supplementary Information"),
    2: ("pdf", "Peer Review File"),
    3: ("pdf", "Description of Additional Supplementary Files"),
    4: ("zip", "Supplementary Data 1-24"),
    5: ("docx", "Reporting Summary"),
    6: ("zip", "Source Data"),
}
GITHUB_REPO = ("https://github.com/docwalitt/National-Institutes-of-Health-Myalgic-Encephalomyelitis-"
               "Chronic-Fatigue-Syndrome-Code-Repository")
GITHUB_ZIP = ("https://codeload.github.com/docwalitt/National-Institutes-of-Health-Myalgic-Encephalomyelitis-"
              "Chronic-Fatigue-Syndrome-Code-Repository/zip/refs/heads/main")
GITHUB_ROOT = "National-Institutes-of-Health-Myalgic-Encephalomyelitis-Chronic-Fatigue-Syndrome-Code-Repository-main"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/{sub}/{name}"
EPMC_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10881493/fullTextXML"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIOPROJECT = "PRJNA954397"
SRA_STUDY = "SRP467038"
PENNSIEVE_ID = 356
PENNSIEVE_API = f"https://api.pennsieve.io/discover/datasets/{PENNSIEVE_ID}"
PENNSIEVE_DOI = "10.26275/ile7-wrsk"
PENNSIEVE_SMALL = ["manifest.json", "metadata/records/file.csv", "metadata/schema.json", "files/Data_share.xlsx"]

GEO_SERIES = {
    "pbmc_rnaseq": {
        "gse": "GSE251872", "tissue": "PBMC", "assay": "bulk RNA-seq, per-sample gene counts",
        "matrices": ["GSE251872-GPL21290_series_matrix.txt.gz", "GSE251872-GPL24676_series_matrix.txt.gz"],
        "suppl": ["GSE251872_RAW.tar"], "unit": "read count",
    },
    "muscle_rnaseq": {
        "gse": "GSE245661", "tissue": "skeletal muscle (vastus lateralis)",
        "assay": "bulk RNA-seq, per-sample gene counts",
        "matrices": ["GSE245661_series_matrix.txt.gz"], "suppl": ["GSE245661_RAW.tar"], "unit": "read count",
    },
    "csf_somalogic": {
        "gse": "GSE251790", "tissue": "cerebrospinal fluid", "assay": "SOMAscan 1.3k aptamer proteomics",
        "matrices": ["GSE251790_series_matrix.txt.gz"],
        "suppl": ["GSE251790_CHI-19-027.Set_001.hybNorm.plateScale.medNorm.20200106.adat.txt.gz"],
        "unit": "RFU (hybNorm, plateScale, medNorm)",
    },
    "serum_somalogic": {
        "gse": "GSE254030", "tissue": "serum (GEO series title says plasma; sample source and paper say serum)",
        "assay": "SOMAscan 1.3k aptamer proteomics",
        "matrices": ["GSE254030_series_matrix.txt.gz"],
        "suppl": ["GSE254030_CHI-19-021.Plasma_for_CFS.adat.txt.gz"], "unit": "RFU",
    },
}
GEO_SUPERSERIES = "GSE251792"
EXTERNAL_GEO = {"GSE130353": "monocyte RNA-seq (external cohort; Supplementary Fig. S17)",
                "GSE156792": "DNA methylation array (external cohort)"}

GROUPS = ("HV", "PI-ME/CFS")
STUDY_NUMBER_NS = "nih_study_number"

# Sample sizes stated in "Description of Additional Supplementary Files" (MOESM3),
# as (n PI-ME/CFS, n HV). Used only when a table/figure file does not let us count.
STATED_N = {
    "14A": (17, 21), "15A": (15, 19), "15B": (10, 7), "15C": (16, 19), "15D": (9, 7),
    "16A": (14, 15), "16B": (6, 7), "16C": (8, 8), "17A": (15, 18), "17B": (15, 18),
    "18": (14, 15), "19A": (13, 12), "19B": (6, 7), "19C": (7, 5), "20": (15, 18), "21": (16, 17),
}
# Figure source-data files whose group-labelled sample columns/rows reproduce a table's samples.
FIGURE_FOR_TABLE = {"16A": "Figure 8A.xlsx", "16B": "Figure 8D.csv", "16C": "Figure 8F.csv",
                    "19A": "Figure 9A.csv", "19B": "Figure 9D.csv", "19C": "Figure 9F.csv"}
TABLE_META = {
    "14": dict(analyte_class="metabolite", tissue="cerebrospinal fluid",
               platform="Metabolon UPLC-MS/MS untargeted metabolomics", accession="source data only"),
    "15": dict(analyte_class="immune_cell_population", tissue=None, platform="flow cytometry (LSR II)",
               accession="source data only"),
    "16": dict(analyte_class="transcript", tissue="PBMC", platform="bulk RNA-seq (limma)", accession="GSE251872"),
    "17": dict(analyte_class="protein_aptamer", tissue=None, platform="SOMAscan 1.3k (SomaLogic)",
               accession=None),
    "18": dict(analyte_class="transposable_element", tissue="PBMC",
               platform="RNA-seq, TEtranscripts + DESeq2", accession="GSE251872"),
    "19": dict(analyte_class="transcript", tissue="skeletal muscle (vastus lateralis)",
               platform="bulk RNA-seq (limma)", accession="GSE245661"),
    "20": dict(analyte_class="lipid_class", tissue="blood (MOESM3: 'plasma ... from serum samples')",
               platform="lipidomics (lipid-class concentration, uM)", accession="source data only"),
    "21": dict(analyte_class="stool_metabolite", tissue="stool", platform="stool NMR metabolomics (uM)",
               accession="source data only"),
}

# Group-labelled figure source data (no participant identifiers) inventoried as modalities.
GROUP_ONLY_SPECS = [
    dict(modality_id="hrv_24h_time_domain", label="24-h ambulatory ECG HRV (avgHR, SDNNI, rMSSD, pNN50)",
         file="Figure 2A frequency domain.xlsx", layout="rows", data_kind="wearable_physiology", wearable=True),
    dict(modality_id="hrv_24h_heart_rate_5min", label="24-h ambulatory heart rate, 5-min intervals",
         file="Figure 2G.xlsx", layout="columns", data_kind="wearable_physiology", wearable=True),
    dict(modality_id="baroreflex_gain", label="Baroreflex cardiovagal gain (baroslope)",
         file="Figure 2H.xlsx", layout="rows", data_kind="autonomic_physiology", wearable=False),
    dict(modality_id="orthostatic_tilt_catecholamines", label="Head-up tilt plasma epinephrine",
         file="Figure S4A.csv", layout="rows", data_kind="autonomic_physiology", wearable=False),
    dict(modality_id="cpet_peak_vo2", label="Cardiopulmonary exercise test, peak VO2",
         file="Figure 5D.xlsx", layout="rows", data_kind="exercise_physiology", wearable=False),
    dict(modality_id="cpet_heart_rate_curve", label="CPET heart rate vs % test time",
         file="Figure 5G.xlsx", layout="columns", data_kind="exercise_physiology", wearable=False),
    dict(modality_id="grip_strength", label="Maximum grip force", file="Figure 3C.xlsx", layout="rows",
         data_kind="strength", wearable=False),
    dict(modality_id="body_composition_dexa", label="DEXA lean body mass", file="Figure S2A.xlsx",
         layout="rows", data_kind="body_composition", wearable=False),
    dict(modality_id="indirect_calorimetry", label="Whole-room calorimetry energy expenditure",
         file="Figure S9A.xlsx", layout="rows", data_kind="metabolic_physiology", wearable=False),
    dict(modality_id="csf_catecholamines", label="CSF DOPA (catechols)", file="Figure 6A.xlsx",
         layout="rows", data_kind="targeted_biochemistry", wearable=False),
    dict(modality_id="plasma_lipidomics_source_data", label="Plasma lipidomics (Fig. S19A)",
         file="Figure S19A.xlsx", layout="rows", data_kind="omics", wearable=False, omics=True),
    dict(modality_id="pbmc_rnaseq_source_data", label="PBMC RNA-seq counts (Fig. 8A source data)",
         file="Figure 8A.xlsx", layout="columns", data_kind="omics", wearable=False, omics=True),
    dict(modality_id="muscle_rnaseq_source_data", label="Muscle RNA-seq counts (Fig. 9A source data)",
         file="Figure 9A.csv", layout="columns", data_kind="omics", wearable=False, omics=True),
]
LETTER_SPECS = [
    dict(modality_id="eefrt_source_data", label="Effort-Expenditure for Rewards Task (Fig. 3A source data)",
         file="Figure 3A.xlsx", id_col="ID", data_kind="behavioural"),
    dict(modality_id="seahorse_post_cpet", label="PBMC Seahorse respiration after CPET (Fig. S9E)",
         file="Figure S9E.csv", id_col="Group", data_kind="cellular_bioenergetics"),
]
WEARABLE_KEYWORDS = ("accelerom", "actigraph", "heart rate variability", "hf power", "lf power",
                     "sd1 heart", "sd2 heart")

RAW = raw_dir(SOURCE_ID)
EXTRACT = RAW / "extracted"
SOURCE_DATA_DIR = EXTRACT / "source_data" / "Figures Source Data Files 12.30.23"
SUPP_DIR = EXTRACT / "supp_data" / "Individual Supplementary Data"
GITHUB_DIR = EXTRACT / "github" / GITHUB_ROOT


# =========================================================================== helpers
def _s(v) -> str:
    """Scalar -> stripped string; missing -> ''."""
    if v is None:
        return ""
    try:
        if pd.isna(v):
            return ""
    except (TypeError, ValueError):
        pass
    return str(v).strip()


_SN_PATTERNS = [re.compile(p, re.I) for p in (
    r"MECFS_(\d{3})\b",          # GEO RNA-seq titles: "... Baseline MECFS_101" / "[MECFS_101]"
    r"#\s*(\d{3})\b",            # GEO CSF SomaLogic titles: "Patient #312" / "Control #113"
    r"WALITT_?(\d{3})$",         # GEO count files: "Walitt101", "FibroWalitt101", "NT_FIBRO_WALITT_101"
    r"^S(\d{3})$",               # PBMC count-file column header "S101"
    r"^(\d{3})(?:\.0+)?$",       # Fig. 6J source data ID column "103"
)]


def study_number(label) -> str | None:
    """NIH study number (1xx healthy volunteer / 3xx PI-ME/CFS) from a label, else None.

    Only the label formats observed in the NIH deposits are recognised. Other
    namespaces (microbiome 'SID_340', mapMECFS 'map000288-00-02', SomaLogic
    'Fatigue-126', figure-local 'HV1'/'HV A') deliberately return None.
    """
    s = _s(label)
    for pat in _SN_PATTERNS:
        m = pat.search(s)
        if m and m.group(1)[0] in "13":
            return m.group(1)
    return None


def native_id(sn: str) -> str:
    return f"MECFS_{sn}"


def participant_id(sn: str) -> str:
    return f"{DATASET_ID}:{native_id(sn)}"


def normalise_group(text) -> str:
    """Map deposit group wording onto 'HV' / 'PI-ME/CFS' ('' when unknown)."""
    s = _s(text).lower()
    if not s:
        return ""
    if "control" in s or "healthy" in s or re.search(r"(^|[^a-z])hv([^a-z]|$)", s):
        return "HV"
    if "patient" in s or "me/cfs" in s or "mecfs" in s or s == "pi-me/cfs":
        return "PI-ME/CFS"
    return ""


def _strip_query(url: str) -> str:
    parts = urlsplit(url or "")
    return urlunsplit((parts.scheme, parts.netloc, parts.path, "", ""))


def _retrieved(rel: str) -> str:
    return load_manifest(SOURCE_ID)["files"].get(rel, {}).get("retrieved_at") or utc_now_iso()


def _write_json(rel: str, obj) -> Path:
    p = RAW / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str))
    return p


def _html_text(content: str, limit: int = 400) -> str:
    t = re.sub(r"<script.*?</script>|<style.*?</style>", " ", content, flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    return html.unescape(re.sub(r"\s+", " ", t)).strip()[:limit]


# =========================================================================== acquisition
@dataclass
class AccessRecord:
    deposit: str
    url: str
    method: str
    http_status: int | None
    outcome: str        # downloaded | cached | metadata_only | login_required | skipped | error
    detail: str
    checked_at: str = field(default_factory=utc_now_iso)


def _download(log: list[AccessRecord], deposit: str, url: str, rel: str, *, display_url: str | None = None,
              **kw) -> Path | None:
    shown = display_url or url
    dest = RAW / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and rel in load_manifest(SOURCE_ID)["files"]:
        log.append(AccessRecord(deposit, shown, "GET", None, "cached",
                                f"{rel} already in MANIFEST.json ({dest.stat().st_size} bytes, retrieved "
                                f"{_retrieved(rel)})"))
        return dest
    try:
        p = download_file(url, SOURCE_ID, rel, **kw)
    except requests.HTTPError as exc:
        status = exc.response.status_code if exc.response is not None else None
        body = _html_text(exc.response.text) if exc.response is not None else ""
        log.append(AccessRecord(deposit, shown, "GET", status, "error", body or str(exc)))
        return None
    except Exception as exc:  # network errors, HTML-instead-of-data, ...
        log.append(AccessRecord(deposit, shown, "GET", None, "error", f"{type(exc).__name__}: {exc}"[:400]))
        return None
    log.append(AccessRecord(deposit, shown, "GET", 200, "downloaded", f"{rel} ({p.stat().st_size} bytes)"))
    return p


def _extract_zip(zip_path: Path | None, dest: Path) -> None:
    if zip_path is None or not zip_path.exists():
        return
    if dest.exists() and any(dest.rglob("*")):
        return
    dest.mkdir(parents=True, exist_ok=True)
    root = dest.resolve()
    with zipfile.ZipFile(zip_path) as zf:
        for m in zf.infolist():
            target = (dest / m.filename).resolve()
            if not str(target).startswith(str(root)):
                raise ValueError(f"unsafe path in {zip_path.name}: {m.filename}")
        zf.extractall(dest)


def _probe(url: str, label: str) -> dict:
    """GET a URL anonymously and record status + what the page says (no credentials)."""
    try:
        r = get(url, reject_html=False, cache_errors=True, timeout=60, max_retries=2)
    except Exception as exc:
        return {"label": label, "url": url, "http_status": None, "error": f"{type(exc).__name__}: {exc}"[:300]}
    ctype = r.headers.get("Content-Type", r.headers.get("content-type", ""))
    text = r.text
    message = ""
    if "json" in ctype:
        try:
            message = json.dumps(r.json().get("error", {}))[:300]
        except Exception:
            message = text[:300]
    else:
        page = _html_text(text, 2000)
        m = re.search(r"Login Required[^.]*\.", page)
        message = m.group(0) if m else page[:300]
    return {"label": label, "url": url, "http_status": r.status, "content_type": ctype,
            "bytes": len(r.content), "login_required": ("Login Required" in text
                                                         or "Authentication is required" in text),
            "message": message, "fetched_at": r.fetched_at}


def _load_catalog(name: str) -> list[dict]:
    p = RAW / "mapmecfs" / name
    if not p.exists():
        return []
    try:
        return json.loads(p.read_text())["result"]["results"]
    except Exception:
        return []


def _acquire_mapmecfs(log: list[AccessRecord]) -> None:
    _download(log, "mapMECFS CKAN catalogue, NIH PI-ME/CFS group (public metadata)",
              f"{CKAN}package_search?fq=groups:{MAPMECFS_GROUP}&rows=1000", "mapmecfs/catalog_group_pi_mecfs.json")
    _download(log, "mapMECFS CKAN catalogue, all public datasets (metadata)",
              f"{CKAN}package_search?rows=1000", "mapmecfs/catalog_all_public.json")
    probes = [_probe(MAPMECFS_GROUP_URL, "group landing page named in the paper"),
              _probe(f"{CKAN}group_show?id={MAPMECFS_GROUP}", "CKAN group_show API"),
              _probe(f"{CKAN}package_list", "CKAN package_list API")]
    group = _load_catalog("catalog_group_pi_mecfs.json")
    for pkg in group:  # one uploaded file per NIH dataset
        res = [r for r in pkg.get("resources", []) if r.get("url_type") == "upload"]
        if res:
            probes.append(_probe(res[0]["url"], f"file download: {pkg['name']} / {res[0].get('name')}"))
    seen_orgs = {"nih-intramural"}
    for pkg in _load_catalog("catalog_all_public.json"):  # one file per other organisation
        org = (pkg.get("organization") or {}).get("name", "")
        res = [r for r in pkg.get("resources", []) if r.get("url_type") == "upload"]
        if org and org not in seen_orgs and res:
            seen_orgs.add(org)
            probes.append(_probe(res[0]["url"], f"file download (other organisation {org}): {pkg['name']}"))
    _write_json("mapmecfs/access_probe.json", probes)
    for p in probes:
        outcome = "login_required" if p.get("login_required") else ("error" if p.get("error") else "metadata_only")
        log.append(AccessRecord(f"mapMECFS {p['label']}", p["url"], "GET", p.get("http_status"), outcome,
                                p.get("message") or p.get("error", "")))


def _acquire_paper(log: list[AccessRecord]) -> None:
    _download(log, "Europe PMC full text (JATS XML) incl. Data Availability statement", EPMC_XML,
              "paper/PMC10881493_fullTextXML.xml")


def availability_statements() -> dict[str, str]:
    """Verbatim 'Data availability' / 'Code availability' text from the Europe PMC JATS XML."""
    p = RAW / "paper" / "PMC10881493_fullTextXML.xml"
    if not p.exists():
        return {}
    root = ET.fromstring(p.read_bytes())
    out = {}
    for sec in root.iter("sec"):
        title = re.sub(r"\s+", " ", "".join((sec.find("title").itertext()) if sec.find("title") is not None else ""))
        if title.strip().lower() in ("data availability", "code availability"):
            body = " ".join("".join(par.itertext()) for par in sec.findall("p"))
            links = [e.get("{http://www.w3.org/1999/xlink}href") for e in sec.iter("ext-link")]
            out[title.strip()] = re.sub(r"\s+", " ", body).strip()
            out[title.strip() + " links"] = " ".join(dict.fromkeys(x for x in links if x))
    return out


def _acquire_nature(log: list[AccessRecord]) -> None:
    for n, (ext, desc) in ESM.items():
        _download(log, f"Nature ESM MOESM{n} ({desc})", ESM_URL.format(n=n, ext=ext),
                  f"nature_esm/41467_2024_45107_MOESM{n}_ESM.{ext}")
    _extract_zip(RAW / "nature_esm/41467_2024_45107_MOESM4_ESM.zip", EXTRACT / "supp_data")
    _extract_zip(RAW / "nature_esm/41467_2024_45107_MOESM6_ESM.zip", EXTRACT / "source_data")


def _acquire_github(log: list[AccessRecord]) -> None:
    p = _download(log, "Authors' GitHub repository (code + source data), branch main", GITHUB_ZIP,
                  "github/docwalitt_code_repository_main.zip")
    _extract_zip(p, EXTRACT / "github")


def _acquire_geo(log: list[AccessRecord]) -> None:
    for spec in GEO_SERIES.values():
        gse = spec["gse"]
        stub = gse[:-3] + "nnn"
        for name in spec["matrices"]:
            _download(log, f"GEO {gse} series matrix", GEO_FTP.format(stub=stub, gse=gse, sub="matrix", name=name),
                      f"geo/{name}")
        for name in spec["suppl"]:
            _download(log, f"GEO {gse} supplementary file", GEO_FTP.format(stub=stub, gse=gse, sub="suppl", name=name),
                      f"geo/{name}")


def _acquire_sra(log: list[AccessRecord]) -> None:
    rel_runinfo = "sra/PRJNA954397_runinfo.csv"
    if not (RAW / rel_runinfo).exists():
        es = get_json(EUTILS + "esearch.fcgi", params={"db": "sra", "term": BIOPROJECT, "retmax": 1000,
                                                       "retmode": "json"})
        ids = es["esearchresult"]["idlist"]
        url = f"{EUTILS}efetch.fcgi?db=sra&id={','.join(ids)}&rettype=runinfo&retmode=text"
    else:
        url = f"{EUTILS}efetch.fcgi?db=sra&term={BIOPROJECT}&rettype=runinfo"
    p = _download(log, f"SRA run metadata for BioProject {BIOPROJECT}", url, rel_runinfo,
                  display_url=f"{EUTILS}efetch.fcgi?db=sra&rettype=runinfo (ids from esearch term={BIOPROJECT})")
    if p is None:
        return
    runinfo = pd.read_csv(p)
    bs = sorted(runinfo["BioSample"].dropna().astype(str).unique())
    _download(log, f"BioSample attributes for BioProject {BIOPROJECT}",
              f"{EUTILS}efetch.fcgi?db=biosample&id={','.join(bs)}&retmode=xml", "sra/PRJNA954397_biosamples.xml",
              display_url=f"{EUTILS}efetch.fcgi?db=biosample&retmode=xml ({len(bs)} BioSample accessions)")


def _acquire_pennsieve(log: list[AccessRecord], *, bulk: bool = False) -> None:
    meta = get_json(PENNSIEVE_API)
    _write_json("pennsieve/dataset_356.json", meta)
    ver = meta["version"]
    listing = {}
    for path in ("", "files", "files/EMG", "files/fMRI", "metadata", "metadata/records"):
        params = {"limit": 500}
        if path:
            params["path"] = path
        listing[path or "/"] = get_json(f"{PENNSIEVE_API}/versions/{ver}/files/browse", params=params)
    _write_json("pennsieve/files_listing.json", listing)
    log.append(AccessRecord("Pennsieve dataset 356 metadata", PENNSIEVE_API, "GET", 200, "metadata_only",
                            f"version {ver}; {meta.get('fileCount')} files; {meta.get('size')} bytes; "
                            f"license {meta.get('license')}"))
    targets = list(PENNSIEVE_SMALL)
    if bulk:
        for folder in ("files/EMG", "files/fMRI"):
            targets += [f["path"] for f in listing[folder]["files"] if f.get("type") == "File"]
    else:
        log.append(AccessRecord("Pennsieve dataset 356 EMG/fMRI signal files", f"{PENNSIEVE_API}/versions/{ver}",
                                "none", None, "skipped",
                                "openly downloadable (anonymous presigned URLs) but not mirrored by default "
                                "(~3.2 GB of .mat/.nii); identifiers taken from manifest.json. "
                                "Run with --pennsieve-bulk to mirror."))
    man = load_manifest(SOURCE_ID)["files"]
    todo = [p for p in targets if f"pennsieve/{p}" not in man or not (RAW / "pennsieve" / p).exists()]
    for p in targets:
        if p not in todo:
            log.append(AccessRecord("Pennsieve dataset 356", f"{PENNSIEVE_API}/versions/{ver}/files?path={p}",
                                    "GET", None, "cached", f"pennsieve/{p} already in MANIFEST.json"))
    if not todo:
        return
    resp = post_json(f"{PENNSIEVE_API}/versions/{ver}/files/download-manifest", json_body={"paths": todo},
                     refresh=True)
    for item in resp.get("data", []):
        path = "/".join(list(item.get("path", [])) + [item["fileName"]])
        rel = f"pennsieve/{path}"
        stable = f"{PENNSIEVE_API}/versions/{ver}/files?path={path}"
        if _download(log, "Pennsieve dataset 356", item["url"], rel, display_url=stable) is not None:
            m = load_manifest(SOURCE_ID)["files"].get(rel)
            if m and "X-Amz" in m.get("url", ""):
                m["url"] = stable
                m["final_url"] = _strip_query(m.get("final_url", ""))
                m["access_note"] = ("downloaded through a presigned S3 URL issued anonymously by POST "
                                    f"{PENNSIEVE_API}/versions/{ver}/files/download-manifest; query string removed")
                _dl._record(SOURCE_ID, rel, m)


def _acquire_future_sources(log: list[AccessRecord]) -> None:
    term = ('("chronic fatigue syndrome"[All Fields] OR "myalgic encephalomyelitis"[All Fields] OR '
            '"ME/CFS"[All Fields]) AND "gse"[Entry Type] AND "Homo sapiens"[Organism]')
    es = get_json(EUTILS + "esearch.fcgi", params={"db": "gds", "term": term, "retmax": 500, "retmode": "json"})
    ids = es["esearchresult"]["idlist"]
    rows = []
    if ids:
        summ = get_json(EUTILS + "esummary.fcgi", params={"db": "gds", "id": ",".join(ids), "retmode": "json"})
        for uid in summ["result"]["uids"]:
            r = summ["result"][uid]
            acc = r["accession"]
            probe = get(GEO_FTP.format(stub=acc[:-3] + "nnn", gse=acc, sub="matrix", name=""),
                        reject_html=False, cache_errors=True, timeout=60)
            n_matrix = len(re.findall(r'href="[^"]*series_matrix[^"]*\.gz"', probe.text)) if probe.status == 200 else 0
            title_mentions = bool(re.search(r"chronic fatigue|myalgic|ME/CFS|CFS", r["title"], re.I))
            rows.append({
                "accession": acc, "title": r["title"], "modality": r["gdstype"], "n_samples": int(r["n_samples"]),
                "platforms": r.get("gpl", ""), "public_date": r.get("pdat", ""),
                "series_matrix_http_status": probe.status, "n_series_matrix_files": n_matrix,
                "openly_downloadable": bool(probe.status == 200 and n_matrix > 0),
                "title_mentions_me_cfs": title_mentions,
                "this_nih_cohort": acc in {s["gse"] for s in GEO_SERIES.values()} | {GEO_SUPERSERIES},
                "query": term,
            })
    geo = pd.DataFrame(rows).sort_values("accession") if rows else pd.DataFrame()
    out = RAW / "future_sources"
    out.mkdir(parents=True, exist_ok=True)
    geo.to_csv(out / "geo_mecfs_series.csv", index=False)
    log.append(AccessRecord("GEO E-utilities search (future sources)", EUTILS + "esearch.fcgi?db=gds", "GET", 200,
                            "metadata_only", f"{len(geo)} human GEO series returned"))
    cat = _load_catalog("catalog_all_public.json")
    rows = []
    acc_re = re.compile(r"(GSE\d+|PRJ[NE][AB]\d+|SRP\d+|ERP\d+|PXD\d+|MTBLS\d+|ST\d{6}|phs\d+)")
    for pkg in cat:
        res = pkg.get("resources", [])
        ext = sorted({m.group(1) for r in res for m in [acc_re.search(r.get("url") or "")] if m})
        rows.append({
            "dataset_name": pkg.get("name"), "title": pkg.get("title"),
            "organization": (pkg.get("organization") or {}).get("title", ""),
            "data_type": pkg.get("data_type", ""), "assay": pkg.get("assay", ""),
            "sample": pkg.get("sample") or pkg.get("other_sample") or "",
            "n_uploaded_files": sum(1 for r in res if r.get("url_type") == "upload"),
            "uploaded_files_open_without_login": False if any(r.get("url_type") == "upload" for r in res) else None,
            "open_external_accessions": ";".join(ext),
            "license": pkg.get("license_id", ""), "metadata_modified": pkg.get("metadata_modified", ""),
            "url": f"https://www.mapmecfs.org/dataset/{pkg.get('name')}",
        })
    pd.DataFrame(rows).to_csv(out / "mapmecfs_catalog.csv", index=False)


def acquire(*, pennsieve_bulk: bool = False) -> list[AccessRecord]:
    """Fetch everything openly downloadable; never raises (errors are logged)."""
    log: list[AccessRecord] = []
    steps = [
        ("paper full text", _acquire_paper),
        ("mapMECFS", _acquire_mapmecfs),
        ("Nature supplementary + source data", _acquire_nature),
        ("GitHub repository", _acquire_github),
        ("GEO", _acquire_geo),
        ("SRA/BioSample", _acquire_sra),
        ("Pennsieve", lambda lg: _acquire_pennsieve(lg, bulk=pennsieve_bulk)),
        ("future-source catalogues", _acquire_future_sources),
    ]
    for name, fn in steps:
        try:
            fn(log)
        except Exception as exc:
            log.append(AccessRecord(name, "", "", None, "error", f"{type(exc).__name__}: {exc}"[:400]))
    _write_json("access_log.json", [asdict(a) for a in log])
    return log


# =========================================================================== parsers
def parse_series_matrix(path: Path) -> pd.DataFrame:
    """GEO series-matrix sample header -> one row per GSM."""
    fields: dict[str, list[str]] = {}
    chars: list[list[str]] = []
    desc: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key, vals = parts[0], [v.strip().strip('"') for v in parts[1:]]
            if key.startswith("!Sample_characteristics"):
                chars.append(vals)
            elif key == "!Sample_description":
                desc.append(vals)
            elif key not in fields:
                fields[key] = vals
    n = len(fields.get("!Sample_geo_accession", []))
    rows = []
    for i in range(n):
        row = {"gsm": fields["!Sample_geo_accession"][i], "title": fields.get("!Sample_title", [""] * n)[i],
               "source_name": fields.get("!Sample_source_name_ch1", [""] * n)[i],
               "description": " | ".join(d[i] for d in desc if i < len(d) and d[i])}
        for c in chars:
            if i < len(c) and ":" in c[i]:
                k, v = c[i].split(":", 1)
                row[f"char_{k.strip().lower()}"] = v.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def geo_samples(modality: str) -> pd.DataFrame:
    spec = GEO_SERIES[modality]
    frames = [parse_series_matrix(RAW / "geo" / m) for m in spec["matrices"] if (RAW / "geo" / m).exists()]
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True).drop_duplicates("gsm")
    df["modality"] = modality
    df["gse"] = spec["gse"]
    df["study_number"] = df["title"].map(study_number)
    df["fatigue_code"] = None
    if modality == "csf_somalogic":
        df["fatigue_code"] = df["description"].str.extract(r"(Fatigue-\d+)")[0]
    if modality == "serum_somalogic":
        df["fatigue_code"] = df["title"].str.extract(r"(Fatigue-\d+)")[0]
        df["study_number"] = None       # native IDs are SomaLogic codes; crosswalk applied later
    grp_src = df.get("char_group", df.get("char_subject status", df["title"]))
    df["group_label"] = [normalise_group(g) or normalise_group(t) for g, t in zip(grp_src, df["title"])]
    df["sex"] = df.get("char_sex", pd.Series([""] * len(df))).fillna("").str.title()
    df["age"] = pd.to_numeric(df.get("char_age", pd.Series([np.nan] * len(df))), errors="coerce")
    df["native_label"] = df["title"]
    return df


def parse_adat(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """SomaLogic ADAT -> (row metadata, SOMAmer column metadata, RFU matrix rows x SOMAmers)."""
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = fh.read().split("\n")
    start = next(i for i, ln in enumerate(lines) if ln.startswith("^TABLE_BEGIN"))
    rows = [ln.rstrip("\r").split("\t") for ln in lines[start + 1:] if ln.strip()]
    name_idx = next(j for j, c in enumerate(rows[0]) if c.strip())
    h = next(k for k, r in enumerate(rows) if r[0].strip())
    colmeta = pd.DataFrame({r[name_idx]: r[name_idx + 1:] for r in rows[:h]})
    nfeat = len(colmeta)
    header = [c.strip() for c in rows[h][:name_idx]]
    data = rows[h + 1:]
    meta = pd.DataFrame([r[:name_idx] for r in data], columns=header)
    vals = np.array([[float(x) if x.strip() else np.nan for x in r[name_idx + 1:name_idx + 1 + nfeat]]
                     for r in data], dtype=float)
    if vals.shape[1] != nfeat:
        raise ValueError(f"{path.name}: {vals.shape[1]} values per row but {nfeat} SOMAmers")
    return meta, colmeta, vals


def parse_count_tar(path: Path) -> pd.DataFrame:
    """GEO RAW.tar of per-sample count files -> long frame (gsm, in-file label, gene, count)."""
    out = []
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            raw = tf.extractfile(m).read()
            df = pd.read_csv(io.BytesIO(gzip.decompress(raw) if m.name.endswith(".gz") else raw), sep="\t")
            if df.shape[1] != 3:
                raise ValueError(f"{m.name}: expected 3 columns, found {list(df.columns)}")
            gid, gname, lab = df.columns
            out.append(pd.DataFrame({
                "gsm": m.name.split("_")[0], "file_member": m.name, "in_file_label": lab,
                "feature_id": df[gid].astype(str), "feature_name": df[gname].astype("string"),
                "value": pd.to_numeric(df[lab], errors="coerce"),
            }))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def read_figure(name: str) -> pd.DataFrame | None:
    p = SOURCE_DATA_DIR / name
    if not name or not p.is_file():
        return None
    if p.suffix == ".xlsx":
        return pd.read_excel(p, header=None)
    return pd.read_csv(p, header=None, low_memory=False)


def csf_metabolomics_rows() -> pd.DataFrame:
    """Fig. 6J source data: CSF metabolomics rows; ID column holds a study number on some rows only."""
    p = SOURCE_DATA_DIR / "Figure 6J.xlsx"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_excel(p).dropna(subset=["CLIENT_SAMPLE_ID"]).reset_index(drop=True).copy()
    df = df.assign(study_number=df["ID"].map(study_number), group_fig6K="", fig6K_row_matches=0)
    k = SOURCE_DATA_DIR / "Figure 6K.xlsx"
    if k.exists():
        # Fig. 6K holds the same values with the group label in CLIENT_SAMPLE_ID: recover each row's group by an
        # exact match on every shared metabolite column (NaN-aware).
        kk = pd.read_excel(k)
        feats = [c for c in df.columns if c in kk.columns and c not in ("CLIENT_SAMPLE_ID", "Birth Sex")]
        A = df[feats].apply(pd.to_numeric, errors="coerce").to_numpy()
        B = kk[feats].apply(pd.to_numeric, errors="coerce").to_numpy()
        groups, nmatch, sex6k = [], [], []
        for i in range(len(A)):
            m = np.where(np.isclose(B, A[i], rtol=1e-9, atol=1e-12, equal_nan=True).all(axis=1))[0]
            nmatch.append(len(m))
            groups.append(normalise_group(kk["CLIENT_SAMPLE_ID"].iloc[m[0]]) if len(m) == 1 else "")
            sex6k.append(kk["Birth Sex"].iloc[m[0]] if len(m) == 1 and "Birth Sex" in kk else np.nan)
        df["group_fig6K"], df["fig6K_row_matches"], df["birth_sex_fig6K"] = groups, nmatch, sex6k
    return df


ANDROGEN_RE = re.compile(r"androst|testosterone|dehydroepiandrosterone|\bdhea", re.I)


def birth_sex_value_alignment(metab: pd.DataFrame) -> dict | None:
    """Do the label-side columns (ID, Birth Sex) of Fig. 6J sit on the same rows as the values?

    Androgen sulfates differ strongly by sex, so if 'Birth Sex' is row-aligned with the metabolite values the
    androgen level separates the two Birth Sex codes within each group block. Returns AUC per block
    (P[value | Birth Sex 0] > value | Birth Sex 1]); ~1 = aligned, ~0.5 = unrelated, ~0 = reversed coding.
    """
    cols = [c for c in metab.columns if ANDROGEN_RE.search(str(c))]
    if not cols or "Birth Sex" not in metab:
        return None
    block = np.where(metab["study_number"].notna(),
                     np.where(metab["study_number"].astype(str).str[0] == "1", "HV", "PI-ME/CFS"),
                     metab["ID"].map(normalise_group))
    out = {"metabolite": cols[0]}
    v = pd.to_numeric(metab[cols[0]], errors="coerce")
    for b in ("HV", "PI-ME/CFS"):
        s0 = v[(block == b) & (metab["Birth Sex"] == 0)].dropna().to_numpy()
        s1 = v[(block == b) & (metab["Birth Sex"] == 1)].dropna().to_numpy()
        if len(s0) and len(s1):
            gt = (s0[:, None] > s1[None, :]).sum() + 0.5 * (s0[:, None] == s1[None, :]).sum()
            out[b] = {"n_sex0": int(len(s0)), "n_sex1": int(len(s1)), "auc": round(float(gt / (len(s0) * len(s1))), 3)}
    return out


# =========================================================================== modalities
@dataclass
class Modality:
    modality_id: str
    label: str
    deposit: str
    accession: str
    url: str
    obtained: bool
    open_download: bool | None
    file_format: str
    data_kind: str
    is_wearable: bool
    is_omics: bool
    tissue: str
    id_namespace: str
    native_ids: list = field(default_factory=list)
    study_numbers: list = field(default_factory=list)
    n_subjects: int | None = None
    n_rows_group_labelled_only: int | None = None
    linkage_method: str = ""
    blocked_reason: str = ""
    notes: str = ""

    @property
    def participant_level(self) -> bool:
        return bool(self.native_ids)


def _group_only_counts(df: pd.DataFrame, layout: str) -> tuple[int, dict]:
    """Count data rows (layout='rows') or sample columns (layout='columns') carrying only group labels."""
    if layout == "columns":
        best = {}
        for i in range(min(5, len(df))):
            row = [_s(x) for x in df.iloc[i].tolist()]
            c = {g: row.count(g) for g in GROUPS}
            if sum(c.values()) > sum(best.values() or [0]):
                best = c
        return sum(best.values()), best
    first = [_s(x) for x in df.iloc[0].tolist()]
    gcol = next((j for j, v in enumerate(first) if v.lower().startswith("group")), 0)
    g = df.iloc[1:, gcol].map(_s).replace("", np.nan).ffill()
    vals = df.iloc[1:, [j for j in range(df.shape[1]) if j != gcol]].apply(pd.to_numeric, errors="coerce")
    has = vals.notna().any(axis=1)
    counts = g[has].map(lambda x: "HV" if normalise_group(x) == "HV" else
                        ("PI-ME/CFS" if normalise_group(x) else "other")).value_counts().to_dict()
    return int(has.sum()), counts


def build_modalities(geo: dict[str, pd.DataFrame], crosswalk: dict[str, str], metab: pd.DataFrame
                     ) -> tuple[list[Modality], list[dict]]:
    """Every modality we obtained or could not obtain, with its identifiers."""
    mods: list[Modality] = []
    checks: list[dict] = []

    for mid, spec in GEO_SERIES.items():
        s = geo.get(mid)
        if s is None or s.empty:
            mods.append(Modality(mid, f"{spec['tissue']} {spec['assay']}", f"GEO {spec['gse']}", spec["gse"],
                                 f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={spec['gse']}", False, True,
                                 "series matrix + supplementary", "omics", False, True, spec["tissue"],
                                 STUDY_NUMBER_NS, blocked_reason="download failed (see access log)"))
            continue
        if mid == "serum_somalogic":
            sns = sorted({crosswalk[c] for c in s["fatigue_code"].dropna() if c in crosswalk})
            mods.append(Modality(
                mid, f"Serum SOMAscan 1.3k proteomics", f"GEO {spec['gse']}", spec["gse"],
                f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={spec['gse']}", True, True,
                "series matrix + ADAT (tab-delimited)", "omics", False, True, spec["tissue"],
                "somalogic_fatigue_code", native_ids=sorted(s["fatigue_code"].dropna().unique()),
                study_numbers=sns, n_subjects=int(s["fatigue_code"].nunique()),
                linkage_method="crosswalk: SomaLogic 'Fatigue-NNN' code -> study number, published per sample in "
                               "GSE251790 (title '#NNN' + description 'Fatigue-NNN')",
                notes="GEO sample titles carry only 'Fatigue-NNN'."))
            continue
        sns = sorted(s["study_number"].dropna().unique())
        mods.append(Modality(
            mid, f"{spec['tissue']} {spec['assay']}", f"GEO {spec['gse']}", spec["gse"],
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={spec['gse']}", True, True,
            "series matrix + " + ("per-sample count files (RAW.tar)" if "rnaseq" in mid else "ADAT (tab-delimited)"),
            "omics", False, True, spec["tissue"], STUDY_NUMBER_NS,
            native_ids=sorted(s["native_label"].unique()), study_numbers=sns, n_subjects=len(sns),
            linkage_method="direct: study number in GEO sample title",
            notes=("also carries SomaLogic 'Fatigue-NNN' codes in the sample description"
                   if mid == "csf_somalogic" else "")))
        prefix_ok = int(((s["study_number"].str[0] == "1") == (s["group_label"] == "HV")).sum())
        checks.append({"check": f"{mid}: study-number prefix (1xx HV / 3xx PI-ME/CFS) agrees with deposited group",
                       "n_checked": int(s["study_number"].notna().sum()), "n_agree": prefix_ok})

    if not metab.empty:
        ided = metab[metab["study_number"].notna()]
        mods.append(Modality(
            "csf_metabolomics", "CSF metabolomics (Metabolon UPLC-MS/MS), batch-normalised imputed values",
            "Nature Source Data, Figure 6J", "source data only", PAPER["url"], True, True, "xlsx", "omics", False,
            True, "cerebrospinal fluid", STUDY_NUMBER_NS,
            native_ids=sorted(ided["ID"].astype(str).map(lambda v: v.split(".")[0]).unique()),
            study_numbers=sorted(ided["study_number"].unique()), n_subjects=int(ided["study_number"].nunique()),
            n_rows_group_labelled_only=int(len(metab) - len(ided)),
            linkage_method="direct: study number in the source-data 'ID' column (only on rows with Birth Sex = 1)",
            notes=f"{len(metab)} rows; {len(metab) - len(ided)} rows carry only a group label in 'ID' and cannot "
                  "be linked; CLIENT_SAMPLE_ID 'NINR-xxxxx' is a lab sample code."))
        matched = ided[ided["group_fig6K"] != ""]
        checks.append({"check": "csf_metabolomics: Fig. 6J rows matched to exactly one Fig. 6K row (all metabolite "
                                "values equal)", "n_checked": int(len(metab)),
                       "n_agree": int((metab["fig6K_row_matches"] == 1).sum())})
        checks.append({"check": "csf_metabolomics: study-number prefix agrees with the Fig. 6K group label",
                       "n_checked": int(len(matched)),
                       "n_agree": int(((matched["study_number"].str[0] == "1") ==
                                       (matched["group_fig6K"] == "HV")).sum())})
        bs = ided["Birth Sex"].value_counts(dropna=False).to_dict()
        checks.append({"check": "csf_metabolomics: 'Birth Sex' of ID-bearing rows",
                       "n_checked": int(len(ided)), "n_agree": int(bs.get(1.0, 0)),
                       "detail": f"Birth Sex counts among ID-bearing rows {bs}; group-label-only rows "
                                 f"{metab.loc[metab['study_number'].isna(), 'Birth Sex'].value_counts().to_dict()}"})
        # corroboration that the 'ID' values are NIH study numbers: the same numbers carry sex in GEO
        geo_sex = {}
        for mid in ("pbmc_rnaseq", "muscle_rnaseq", "csf_somalogic"):
            s = geo.get(mid)
            if s is not None and not s.empty:
                for sn, sx in zip(s["study_number"], s["sex"]):
                    if sn and sx:
                        geo_sex.setdefault(sn, set()).add(sx)
        in_geo = [sn for sn in ided["study_number"] if sn in geo_sex]
        checks.append({"check": "csf_metabolomics: GEO sex of the ID-bearing study numbers (all have Birth Sex = 1)",
                       "n_checked": len(in_geo), "n_agree": sum(geo_sex[sn] == {"Female"} for sn in in_geo),
                       "detail": "n_agree = study numbers recorded Female in every GEO deposit; "
                                 f"{len(ided) - len(in_geo)} ID-bearing study number(s) not in GEO"})
        al = birth_sex_value_alignment(metab)
        if al:
            blocks = [b for b in ("HV", "PI-ME/CFS") if b in al]
            checks.append({"check": "csf_metabolomics: Fig. 6J label columns (ID, Birth Sex) are row-aligned with the "
                                    "metabolite values (androgen sulfate AUC by Birth Sex >= 0.8 within each group)",
                           "n_checked": len(blocks), "n_agree": sum(al[b]["auc"] >= 0.8 for b in blocks),
                           "detail": json.dumps(al)})
        if "birth_sex_fig6K" in metab:
            mt = metab[metab["fig6K_row_matches"] == 1]
            same = (pd.to_numeric(mt["birth_sex_fig6K"], errors="coerce") == mt["Birth Sex"])
            unl = int(((metab["fig6K_row_matches"] == 1) & (metab["group_fig6K"] == "")).sum())
            checks.append({"check": "csf_metabolomics: Fig. 6K 'Birth Sex' equals Fig. 6J 'Birth Sex' on the matched "
                                    "value row (published-file consistency)", "n_checked": int(len(mt)),
                           "n_agree": int(same.sum()),
                           "detail": f"disagreements by 6J block: {metab.loc[mt.index[~same.to_numpy()], 'ID'].astype(str).map(lambda x: 'HV' if normalise_group(x) == 'HV' or x.startswith('1') else 'PI-ME/CFS').value_counts().to_dict()}; "
                                     f"{unl} matched 6K value row(s) carry no group label. Fig. 6K 'Birth Sex' is not "
                                     "used; Fig. 6J is (it is row-aligned per the androgen check)"})

    # stool metagenomics: SRA BioSample + GitHub metadata use their own subject codes
    sids_sra, sids_gh = _microbiome_ids()
    if sids_sra or sids_gh:
        mods.append(Modality(
            "stool_metagenomics", "Stool shotgun metagenomics", f"SRA {SRA_STUDY} / BioProject {BIOPROJECT}; "
            "GitHub microbiome metadata", f"{SRA_STUDY}", f"https://www.ncbi.nlm.nih.gov/bioproject/{BIOPROJECT}",
            True, True, "SRA metadata (runinfo CSV, BioSample XML); xlsx", "omics", False, True, "stool",
            "microbiome_subject_sid", native_ids=sorted(set(sids_sra) | set(sids_gh)),
            n_subjects=len(set(sids_sra) | set(sids_gh)),
            notes=f"BioSample Subject_ID codes ({len(sids_sra)}) and GitHub metadata Subject_ID "
                  f"({len(sids_gh)}); overlap between the two {len(set(sids_sra) & set(sids_gh))}. The numbers in "
                  "'SID_NNN' are not study numbers (e.g. HV and ME/CFS codes both occur in the 1xx-9xx range). "
                  "Raw reads not mirrored."))

    # EEfRT with mapMECFS identifiers (GitHub)
    eef = GITHUB_DIR / "EEfRT" / "EEfRT_Full_Data_SAS_mapID.xlsx"
    if eef.exists():
        e = pd.read_excel(eef)
        ids = sorted(e["mapID"].dropna().astype(str).unique())
        seg = e.drop_duplicates("mapID")
        agree = int((seg["mapID"].str.split("-").str[1].astype(int) == seg["Status_MECFS_is_1"]).sum())
        checks.append({"check": "EEfRT mapID middle segment (-00-/-01-) equals Status_MECFS_is_1",
                       "n_checked": int(len(seg)), "n_agree": agree})
        mods.append(Modality(
            "eefrt_mapid", "Effort-Expenditure for Rewards Task trial data (with mapMECFS IDs)",
            "GitHub docwalitt repository, EEfRT/EEfRT_Full_Data_SAS_mapID.xlsx", "GitHub", GITHUB_REPO, True, True,
            "xlsx", "behavioural", False, False, "", "mapmecfs_mapid", native_ids=ids, n_subjects=len(ids),
            notes=f"{len(e)} trial rows; mapIDs 'map000NNN-GG-VV'; no public crosswalk from mapID to NIH study "
                  f"number; adjudicated flag present ({int(seg['Adjudicated_is_1'].sum())} of {len(seg)} "
                  "participants flagged adjudicated)."))

    # Pennsieve (TMS/EMG + fMRI)
    pm = RAW / "pennsieve" / "manifest.json"
    if pm.exists():
        files = json.loads(pm.read_text()).get("files", [])
        for folder, ns, lab, fmt in (("files/EMG/", "pennsieve_emg_index", "Repetitive grip EMG/TMS signals", ".mat"),
                                     ("files/fMRI/", "pennsieve_fmri_index", "Grip-task fMRI", ".nii/.1D")):
            labels = sorted({re.sub(r"^(HV|MECFS)0*(\d+)$", r"\1\2",
                                    re.sub(r"(_EMG\.mat|_4bl\.1D|\.nii)$", "", f["name"].strip()))
                             for f in files if f["path"].startswith(folder)})
            size = sum(int(f.get("size", 0)) for f in files if f["path"].startswith(folder))
            mods.append(Modality(
                f"pennsieve_{folder.split('/')[1].lower()}", lab, f"Pennsieve dataset {PENNSIEVE_ID} "
                f"(DOI {PENNSIEVE_DOI})", f"doi:{PENNSIEVE_DOI}", f"https://discover.pennsieve.io/datasets/{PENNSIEVE_ID}",
                True, True, fmt, "neurophysiology_imaging", False, False, "", ns, native_ids=labels,
                n_subjects=len(labels),
                notes=f"labels are dataset-local indices (HV1.. / MECFS1..); {size} bytes of signal files listed in "
                      "manifest.json (not mirrored by default)."))

    # figure-local letter codes
    for spec in LETTER_SPECS:
        df = read_figure(spec["file"])
        if df is None:
            continue
        hdr = [_s(x) for x in df.iloc[0].tolist()]
        col = hdr.index(spec["id_col"]) if spec["id_col"] in hdr else 0
        ids = sorted({v for v in df.iloc[1:, col].map(_s) if v})
        extra = ""
        if "Valid Data_is_1" in hdr:
            vcol = hdr.index("Valid Data_is_1")
            valid = df.iloc[1:].groupby(df.iloc[1:, col].map(_s))[vcol].apply(
                lambda s: pd.to_numeric(s, errors="coerce").max())
            extra = f" {int((valid == 1).sum())} codes have Valid Data = 1 ({', '.join(valid[valid != 1].index)} not)."
        mods.append(Modality(
            spec["modality_id"], spec["label"], f"Nature Source Data, {spec['file']}", "source data only",
            PAPER["url"], True, True, spec["file"].rsplit(".", 1)[-1], spec["data_kind"], False, False, "",
            f"figure_local_letter:{spec['file'].rsplit('.', 1)[0]}", native_ids=ids, n_subjects=len(ids),
            notes="pseudonymous letters ('HV A', 'PI-ME/CFS A', ...) assigned per figure; no crosswalk." + extra))

    # group-labelled source data (no identifiers)
    for spec in GROUP_ONLY_SPECS:
        df = read_figure(spec["file"])
        if df is None:
            continue
        n, counts = _group_only_counts(df, spec["layout"])
        mods.append(Modality(
            spec["modality_id"], spec["label"], f"Nature Source Data, {spec['file']}", "source data only",
            PAPER["url"], True, True, spec["file"].rsplit(".", 1)[-1], spec["data_kind"], spec["wearable"],
            spec.get("omics", False), "", "group_label_only", n_rows_group_labelled_only=n,
            notes=f"{n} per-participant {'columns' if spec['layout'] == 'columns' else 'rows'} labelled only by "
                  f"group {counts}; no participant identifier."))

    # mapMECFS gated files
    for pkg in _load_catalog("catalog_group_pi_mecfs.json"):
        if "dictionar" in pkg["name"] or "documentation" in pkg["name"]:
            continue
        for r in pkg.get("resources", []):
            if r.get("url_type") != "upload" or (r.get("name") or "").strip().lower() == "search terms":
                continue
            name = (r.get("name") or "").strip()
            low = name.lower()
            wear = any(k in low for k in WEARABLE_KEYWORDS)
            omics = any(k in pkg["name"] for k in ("rnaseq", "somalogic", "microbiome", "metabolomics", "lipidomics",
                                                   "nuclear-magnetic")) or "mitochondrial" in low
            kind = ("wearable_physiology" if wear else "omics" if omics else
                    "cognitive" if "neurocognitive" in low else
                    "targeted_biochemistry" if "csf" in low and "catechol" in low else
                    "autonomic_physiology" if any(k in low for k in ("tilt", "valsalva", "catechol")) else
                    "exercise_physiology" if "cpet" in low or "cardiopulmonary" in low else "clinical_or_other")
            mods.append(Modality(
                f"mapmecfs_{r['id'][:8]}", name, f"mapMECFS {pkg['name']}", r["id"],
                f"https://www.mapmecfs.org/dataset/{pkg['name']}/resource/{r['id']}", False, False,
                (r.get("format") or "").upper() or "UNKNOWN", kind, wear, omics, "", "not_obtained",
                blocked_reason="mapMECFS requires a registered account: anonymous download returns HTTP 403 "
                               "'Login Required. Please login to view this feature.'",
                notes=f"file size per CKAN metadata: {r.get('size')} bytes; datastore_active="
                      f"{r.get('datastore_active')} (not queried)"))
    return mods, checks


def _microbiome_ids() -> tuple[list[str], list[str]]:
    sra, gh = [], []
    xmlp = RAW / "sra" / "PRJNA954397_biosamples.xml"
    if xmlp.exists():
        root = ET.fromstring(xmlp.read_bytes())
        for b in root.findall("BioSample"):
            attrs = {a.get("attribute_name"): (a.text or "") for a in b.findall("./Attributes/Attribute")}
            if attrs.get("Subject_ID"):
                sra.append(attrs["Subject_ID"].strip())
    mp = GITHUB_DIR / "data" / "Walitt2023_Microbiome_metadata_and_PPM_relative_abundance.xlsx"
    if mp.exists():
        gh = pd.read_excel(mp, sheet_name="Metadata")["Subject_ID"].dropna().astype(str).str.strip().tolist()
    return sorted(set(sra)), sorted(set(gh))


def overlap_table(mods: list[Modality]) -> pd.DataFrame:
    """Exact identifier overlaps for every pair of modalities.

    Pairs are comparable only when both carry NIH study numbers (directly or via
    the published SomaLogic crosswalk) or both use the same identifier
    namespace. Everything else is reported as not comparable with the reason and
    n_overlap = NA: the overlap of participants is unknown there, not zero.
    """
    rows = []
    for a, b in itertools.combinations(mods, 2):
        rec = {"modality_a": a.modality_id, "modality_b": b.modality_id, "namespace_a": a.id_namespace,
               "namespace_b": b.id_namespace, "n_a": a.n_subjects, "n_b": b.n_subjects,
               "a_is_wearable": a.is_wearable, "b_is_wearable": b.is_wearable,
               "a_is_omics": a.is_omics, "b_is_omics": b.is_omics}
        if not a.obtained or not b.obtained:
            rec.update(comparable=False, n_overlap=None, overlap_ids="",
                       reason="not obtained: " + "; ".join(sorted({m.blocked_reason for m in (a, b) if not m.obtained})))
        elif not a.participant_level or not b.participant_level:
            rec.update(comparable=False, n_overlap=None, overlap_ids="",
                       reason="no participant identifiers in " + ", ".join(
                           m.modality_id for m in (a, b) if not m.participant_level) + " (group label only)")
        elif a.study_numbers and b.study_numbers:
            ov = sorted(set(a.study_numbers) & set(b.study_numbers))
            via = [m.modality_id for m in (a, b) if m.id_namespace != STUDY_NUMBER_NS]
            rec.update(comparable=True, n_overlap=len(ov), overlap_ids=",".join(ov),
                       reason="NIH study number" + (f" (after published crosswalk for {', '.join(via)})" if via else ""))
        elif a.id_namespace == b.id_namespace:
            ov = sorted(set(a.native_ids) & set(b.native_ids))
            rec.update(comparable=True, n_overlap=len(ov), overlap_ids=",".join(ov), reason=f"same namespace {a.id_namespace}")
        else:
            rec.update(comparable=False, n_overlap=None, overlap_ids="",
                       reason=f"different identifier namespaces ({a.id_namespace} vs {b.id_namespace}); "
                              "no public crosswalk")
        rows.append(rec)
    out = pd.DataFrame(rows)
    for c in ("n_a", "n_b", "n_overlap"):
        if c in out:
            out[c] = pd.array([None if pd.isna(v) else int(v) for v in out[c]], dtype="Int64")
    return out


# =========================================================================== linked subset
def build_linked(geo: dict[str, pd.DataFrame], crosswalk: dict[str, str], metab: pd.DataFrame
                 ) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """participants + long participant-linked omics values (study-number modalities only)."""
    checks: list[dict] = []
    frames = []
    samples = []   # (study_number, modality, gsm, group, sex, age, fatigue_code, linkage_method)

    for mid in ("pbmc_rnaseq", "muscle_rnaseq"):
        s = geo.get(mid)
        spec = GEO_SERIES[mid]
        tarp = RAW / "geo" / spec["suppl"][0]
        if s is None or s.empty or not tarp.exists():
            continue
        c = parse_count_tar(tarp)
        g2sn = dict(zip(s["gsm"], s["study_number"]))
        c["study_number"] = c["gsm"].map(g2sn)
        lab_sn = c.drop_duplicates("gsm").assign(
            label_sn=lambda d: d["in_file_label"].map(study_number),
            file_sn=lambda d: d["file_member"].str.replace(r"\.txt(\.gz)?$", "", regex=True).map(study_number))
        agree = int(((lab_sn["label_sn"] == lab_sn["study_number"]) & (lab_sn["file_sn"] == lab_sn["study_number"])).sum())
        checks.append({"check": f"{mid}: study number in GEO title = count-file name = in-file column header",
                       "n_checked": int(len(lab_sn)), "n_agree": agree})
        c = c[c["study_number"].notna()]
        frames.append(pd.DataFrame({
            "study_number": c["study_number"], "modality": mid, "tissue": spec["tissue"], "assay": spec["assay"],
            "feature_id": c["feature_id"], "feature_name": c["feature_name"], "value": c["value"],
            "value_unit": spec["unit"], "sample_accession": c["gsm"], "native_sample_label": c["in_file_label"],
            "linkage_method": "direct_study_number", "source_file": f"geo/{spec['suppl'][0]}::" + c["file_member"],
            "source_accession": spec["gse"],
        }))
        for _, r in s.iterrows():
            samples.append((r["study_number"], mid, r["gsm"], r["group_label"], r["sex"], np.nan, None,
                            "direct_study_number"))

    for mid in ("csf_somalogic", "serum_somalogic"):
        s = geo.get(mid)
        spec = GEO_SERIES[mid]
        adp = RAW / "geo" / spec["suppl"][0]
        if s is None or s.empty or not adp.exists():
            continue
        meta, colmeta, vals = parse_adat(adp)
        keep = (meta["SampleType"] == "Sample").to_numpy()
        meta, vals = meta[keep].reset_index(drop=True), vals[keep]
        fat2gsm = dict(zip(s["fatigue_code"], s["gsm"]))
        method = "direct_study_number" if mid == "csf_somalogic" else "crosswalk_fatigue_code_via_GSE251790"
        checks.append({"check": f"{mid}: ADAT SampleId (Fatigue-NNN) matched to a GEO sample",
                       "n_checked": int(len(meta)), "n_agree": int(meta["SampleId"].isin(fat2gsm).sum())})
        sn = meta["SampleId"].map(crosswalk)
        ok = sn.notna().to_numpy()
        nf = vals.shape[1]
        long = pd.DataFrame({
            "study_number": np.repeat(sn[ok].to_numpy(), nf),
            "feature_id": np.tile(colmeta["SeqId"].to_numpy(), int(ok.sum())),
            "feature_name": np.tile((colmeta["EntrezGeneSymbol"] + " | " + colmeta["Target"]).to_numpy(), int(ok.sum())),
            "value": vals[ok].ravel(),
            "sample_accession": np.repeat(meta.loc[ok, "SampleId"].map(fat2gsm).to_numpy(), nf),
            "native_sample_label": np.repeat(meta.loc[ok, "SampleId"].to_numpy(), nf),
        })
        long["modality"], long["tissue"], long["assay"] = mid, spec["tissue"], spec["assay"]
        long["value_unit"], long["linkage_method"] = spec["unit"], method
        long["source_file"], long["source_accession"] = f"geo/{spec['suppl'][0]}", spec["gse"]
        frames.append(long)
        for _, r in s.iterrows():
            snr = crosswalk.get(r["fatigue_code"])
            if snr:
                samples.append((snr, mid, r["gsm"], r["group_label"], r["sex"], r["age"], r["fatigue_code"], method))

    if not metab.empty:
        ided = metab[metab["study_number"].notna()]
        feats = [c for c in metab.columns if c not in ("CLIENT_SAMPLE_ID", "ID", "Birth Sex", "study_number",
                                                       "group_fig6K", "fig6K_row_matches", "birth_sex_fig6K")]
        m = ided.melt(id_vars=["study_number", "CLIENT_SAMPLE_ID"], value_vars=feats, var_name="feature_name",
                      value_name="value")
        frames.append(pd.DataFrame({
            "study_number": m["study_number"], "modality": "csf_metabolomics", "tissue": "cerebrospinal fluid",
            "assay": "Metabolon UPLC-MS/MS (batch-normalised, imputed)", "feature_id": m["feature_name"],
            "feature_name": m["feature_name"], "value": pd.to_numeric(m["value"], errors="coerce"),
            "value_unit": "batch-normalised relative abundance", "sample_accession": m["CLIENT_SAMPLE_ID"],
            "native_sample_label": m["study_number"], "linkage_method": "direct_study_number",
            "source_file": "extracted/source_data/Figure 6J.xlsx", "source_accession": "Nature Source Data Fig. 6J",
        }))
        for _, r in ided.iterrows():
            samples.append((r["study_number"], "csf_metabolomics", r["CLIENT_SAMPLE_ID"], r["group_fig6K"], "",
                            np.nan, None, "direct_study_number"))

    omics = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    smp = pd.DataFrame(samples, columns=["study_number", "modality", "sample", "group_label", "sex", "age",
                                         "fatigue_code", "linkage_method"])
    parts = _participants(smp, checks)
    return parts, omics, checks


def _participants(smp: pd.DataFrame, checks: list[dict]) -> pd.DataFrame:
    if smp.empty:
        return pd.DataFrame()
    mods = ["pbmc_rnaseq", "muscle_rnaseq", "csf_somalogic", "serum_somalogic", "csf_metabolomics"]
    rows = []
    for sn, g in smp.groupby("study_number"):
        groups = {m: x for m, x in zip(g["modality"], g["group_label"]) if x}
        sexes = {m: x for m, x in zip(g["modality"], g["sex"]) if x}
        ages = g["age"].dropna().unique()
        row = {
            "study_number": sn, "native_id": native_id(sn),
            "group_label": ("|".join(sorted(set(groups.values()))) if groups else UNKNOWN),
            "group_label_by_source": json.dumps(groups, sort_keys=True),
            "group_concordant_across_sources": len(set(groups.values())) <= 1,
            "group_from_study_number_prefix": "HV" if sn.startswith("1") else "PI-ME/CFS",
            "sex": ("|".join(sorted(set(sexes.values()))) if sexes else UNKNOWN),
            "sex_by_source": json.dumps(sexes, sort_keys=True),
            "sex_concordant_across_sources": len(set(sexes.values())) <= 1,
            "age_years_somalogic": float(ages[0]) if len(ages) == 1 else np.nan,
            "age_concordant_across_sources": len(ages) <= 1,
            "somalogic_fatigue_code": next((x for x in g["fatigue_code"] if isinstance(x, str) and x), None),
        }
        for m in mods:
            row[f"has_{m}"] = m in set(g["modality"])
        row["n_open_omics_modalities"] = int(sum(row[f"has_{m}"] for m in mods))
        row["has_wearable_or_hrv_data_open"] = False
        row["in_published_analytic_cohort"] = UNKNOWN
        rows.append(row)
    p = pd.DataFrame(rows).sort_values("study_number").reset_index(drop=True)
    checks.append({"check": "participants: group label concordant across open deposits",
                   "n_checked": int(len(p)), "n_agree": int(p["group_concordant_across_sources"].sum())})
    checks.append({"check": "participants: sex concordant across open deposits",
                   "n_checked": int(len(p)), "n_agree": int(p["sex_concordant_across_sources"].sum()),
                   "detail": "; ".join(f"{r.native_id}: {r.sex_by_source}" for r in p.itertuples()
                                       if not r.sex_concordant_across_sources)})
    return p


def somalogic_crosswalk(geo: dict[str, pd.DataFrame]) -> tuple[dict[str, str], list[dict]]:
    """Fatigue-NNN -> study number (from GSE251790) and its corroboration against GSE254030."""
    csf, ser = geo.get("csf_somalogic"), geo.get("serum_somalogic")
    if csf is None or csf.empty:
        return {}, []
    cw = {f: sn for f, sn in zip(csf["fatigue_code"], csf["study_number"]) if isinstance(f, str) and sn}
    checks = []
    if ser is not None and not ser.empty:
        m = csf.merge(ser, on="fatigue_code", suffixes=("_csf", "_ser"))
        checks.append({"check": "SomaLogic Fatigue-NNN codes present in both CSF (GSE251790) and serum (GSE254030)",
                       "n_checked": int(ser["fatigue_code"].nunique()), "n_agree": int(len(m))})
        for col in ("sex", "age", "group_label"):
            checks.append({"check": f"crosswalk corroboration: {col} agrees between CSF and serum for the same "
                                    "Fatigue code", "n_checked": int(len(m)),
                           "n_agree": int((m[f"{col}_csf"] == m[f"{col}_ser"]).sum())})
        adp = RAW / "geo" / GEO_SERIES["serum_somalogic"]["suppl"][0]
        if adp.exists():
            meta, _, _ = parse_adat(adp)
            notes = meta.loc[meta["SampleDescription"].fillna("").str.strip() != "", ["SampleId", "SampleDescription"]]
            checks.append({"check": "serum ADAT free-text SampleDescription notes", "n_checked": int(len(meta)),
                           "n_agree": int(len(notes)),
                           "detail": "; ".join(f"{a}: '{b}'" for a, b in notes.itertuples(index=False))})
    return cw, checks


# =========================================================================== condition evidence
def _stated(tid: str) -> tuple[int | None, int | None]:
    return STATED_N.get(tid, (None, None))


def _group_matrix(fig: str) -> tuple[pd.DataFrame, pd.Series, list[str]] | None:
    """Gene x sample matrix from a group-labelled figure file: (values, gene names, labels)."""
    df = read_figure(fig)
    if df is None:
        return None
    lab_row = next((i for i in range(min(5, len(df)))
                    if sum(_s(x) in GROUPS for x in df.iloc[i].tolist()) >= 3), None)
    if lab_row is None:
        return None
    labels = [_s(x) for x in df.iloc[lab_row].tolist()]
    cols = [j for j, x in enumerate(labels) if x in GROUPS]
    col0 = df.iloc[:, 0].map(_s)
    data = df[col0.str.startswith("ENSG")]
    vals = data.iloc[:, cols].apply(pd.to_numeric, errors="coerce")
    vals.index = data.iloc[:, 0].map(_s).str.split(".").str[0]
    names = pd.Series(data.iloc[:, 1].map(_s).to_numpy(), index=vals.index)
    return vals, names, [labels[j] for j in cols]


def _verify_logfc(tab: pd.DataFrame, fig: str) -> dict:
    """Sign agreement of published logFC with log2-CPM group means in the figure source data."""
    gm = _group_matrix(fig)
    if gm is None:
        return {"figure": fig, "n_joined": 0, "sign_concordance": None, "n_hv": None, "n_pi": None}
    vals, names, labs = gm
    cpm = np.log2(vals / vals.sum() * 1e6 + 1)
    hv = [i for i, x in enumerate(labs) if x == "HV"]
    pi = [i for i, x in enumerate(labs) if x == "PI-ME/CFS"]
    diff = cpm.iloc[:, pi].mean(axis=1) - cpm.iloc[:, hv].mean(axis=1)
    if tab["ensembl_gene_id"].notna().any():
        key, ref = tab["ensembl_gene_id"], diff[~diff.index.duplicated()]
    else:
        ref = pd.Series(diff.to_numpy(), index=names.to_numpy())
        ref = ref[~ref.index.duplicated()]
        key = tab["gene_symbol"]
    j = pd.DataFrame({"logfc": tab["effect_value"].to_numpy(), "ref": key.map(ref).to_numpy()}).dropna()
    j = j[(j["logfc"] != 0) & (j["ref"] != 0)]
    conc = float((np.sign(j["logfc"]) == np.sign(j["ref"])).mean()) if len(j) else None
    return {"figure": fig, "n_joined": int(len(j)), "sign_concordance": conc, "n_hv": len(hv), "n_pi": len(pi)}


def _blocks(header: list[str]) -> list[tuple[int, int]]:
    """Contiguous runs of non-empty header cells -> [(start, end_inclusive)]."""
    runs, start = [], None
    for j, v in enumerate(header + [""]):
        if v and start is None:
            start = j
        elif not v and start is not None:
            runs.append((start, j - 1))
            start = None
    return runs


def _sd(n: int) -> pd.DataFrame | None:
    p = SUPP_DIR / f"Supplementary Data {n}.xlsx"
    return pd.read_excel(p, header=None) if p.exists() else None


def _table_id(title: str, fallback: str) -> str:
    m = re.search(r"Supplementary\s+Data:?\s*(\d+[A-Z]?)", title)
    return m.group(1) if m else fallback


def _subgroup(text: str) -> str:
    t = text.lower()
    return "female" if re.search(r"\bfemale", t) else ("male" if re.search(r"\bmale", t) else "all")


def _num(series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def parse_de_tables() -> tuple[list[pd.DataFrame], list[dict]]:
    """Supplementary Data 16 and 19 (limma DE), with sign-convention verification."""
    out, checks = [], []
    for n in (16, 19):
        d = _sd(n)
        if d is None:
            continue
        hrow = next(i for i in range(len(d)) if "logFC" in [_s(x) for x in d.iloc[i].tolist()])
        header = [_s(x) for x in d.iloc[hrow].tolist()]
        titles = [(j, _s(d.iat[i, j])) for i in range(hrow) for j in range(d.shape[1])
                  if re.match(r"Supplementary\s+Data:?\s*\d+[A-Z]", _s(d.iat[i, j]))]
        for k, (a, b) in enumerate(_blocks(header)):
            tj = [t for j, t in titles if a - 1 <= j <= b]
            title = tj[0] if tj else f"Supplementary Data {n}"
            context = " ".join(_s(d.iat[i, j]) for i in range(hrow) for j in range(a, b + 1))
            tid = _table_id(title, f"{n}{'ABC'[k]}")
            blk = d.iloc[hrow + 1:, a:b + 1].copy()
            blk.columns = header[a:b + 1]
            blk = blk[_num(blk["p-value"]).notna()]
            ens = blk.get("ensembl_gene_id", blk.get("ensembl_gene_id_version"))
            tab = pd.DataFrame({
                "supplementary_table": f"Supplementary Data {tid}",
                "table_id": tid, "subgroup": _subgroup(context),
                "analyte_id": (ens if ens is not None else blk["external_gene_name"]).map(_s).to_numpy(),
                "analyte_id_type": "ensembl_gene_id" if ens is not None else "gene_symbol",
                "gene_symbol": blk["external_gene_name"].map(_s).to_numpy(),
                "ensembl_gene_id": (ens.map(_s).str.split(".").str[0].to_numpy() if ens is not None else None),
                "entrez_gene_id": blk["entrezgene_id"].map(lambda v: _s(v).split(".")[0]).to_numpy()
                if "entrezgene_id" in blk else None,
                "effect_measure": "limma logFC (PI-ME/CFS minus HV)",
                "effect_value": _num(blk["logFC"]).to_numpy(),
                "p_value": _num(blk["p-value"]).to_numpy(), "adj_p_value": _num(blk["adj p-value"]).to_numpy(),
                "source_row": (blk.index + 1).to_numpy(),
            })
            tab["analyte_name"] = tab["gene_symbol"]
            v = _verify_logfc(tab, FIGURE_FOR_TABLE.get(tid, ""))
            c = v["sign_concordance"]
            ok = c is not None and c >= 0.70     # a reversed convention would give c < 0.5
            tab["direction"] = np.where(not ok, UNKNOWN, np.where(tab["effect_value"] > 0, "up_in_PI-ME/CFS",
                                        np.where(tab["effect_value"] < 0, "down_in_PI-ME/CFS", "no_change")))
            tab["direction_check_concordance"] = c
            if c is not None and c >= 0.95:
                tab["direction_basis"] = (f"published logFC (PI-ME/CFS minus HV); orientation verified: sign agrees "
                                          f"with unadjusted log2-CPM group means in {v['figure']} for {c:.3f} of "
                                          f"{v['n_joined']} genes")
            elif ok:
                tab["direction_basis"] = (f"published logFC (PI-ME/CFS minus HV); orientation supported but only "
                                          f"{c:.3f} of {v['n_joined']} genes agree in sign with unadjusted log2-CPM "
                                          f"group means in {v['figure']} (published model likely covariate-adjusted); "
                                          "direction follows the published estimate")
            else:
                tab["direction_basis"] = f"published logFC; orientation NOT verified ({v})"
            if v["n_pi"]:
                tab["n_cases"], tab["n_controls"] = v["n_pi"], v["n_hv"]
                tab["n_basis"] = f"counted group-labelled sample columns in {v['figure']}"
            else:
                tab["n_cases"], tab["n_controls"] = _stated(tid)
                tab["n_basis"] = "stated in MOESM3"
            tab["statistical_test"] = ("limma moderated t (Figs. 8/9 legends); MOESM3 description says two-sided "
                                       "Wilcoxon")
            meta = TABLE_META[str(n)]
            tab["analyte_class"], tab["tissue"], tab["platform"] = meta["analyte_class"], meta["tissue"], meta["platform"]
            tab["underlying_data_accession"] = meta["accession"]
            out.append(tab)
            checks.append({"check": f"SD{tid}: logFC sign vs {v['figure']} group means",
                           "n_checked": v["n_joined"],
                           "n_agree": int(round((v["sign_concordance"] or 0) * v["n_joined"]))})
    return out, checks


def parse_metabolomics_table() -> tuple[list[pd.DataFrame], list[dict]]:
    d = _sd(14)
    if d is None:
        return [], []
    hrow = next(i for i in range(len(d)) if _s(d.iat[i, 0]) == "CHEM_ID")
    end = next((i for i in range(hrow + 1, len(d)) if _s(d.iat[i, 0]).startswith("Supplementary")), len(d))
    blk = d.iloc[hrow + 1:end, :9].copy()
    blk.columns = [_s(x) for x in d.iloc[hrow, :9].tolist()]
    blk = blk[blk["CHEM_ID"].map(_s).str.fullmatch(r"\d+")]
    fc, p, q = _num(blk["fold change"]), _num(blk["p-value"]), _num(blk["adjusted p-value"])
    # verify fold-change orientation (PI-ME/CFS / HV) against Fig. 6H group means
    conc, nj, n_pi, n_hv = None, 0, None, None
    h6 = SOURCE_DATA_DIR / "Figure 6H.xlsx"
    if h6.exists():
        h = pd.read_excel(h6)
        num = h.drop(columns=["ID"]).apply(pd.to_numeric, errors="coerce")
        ratio = num[h["ID"] == "PI-ME/CFS"].mean() / num[h["ID"] == "HV"].mean()
        ratio = ratio[~ratio.index.duplicated()]
        j = pd.DataFrame({"fc": fc.to_numpy(), "r": blk["CHEMICAL_NAME"].map(_s).map(ratio).to_numpy()}).dropna()
        j = j[(j["fc"] != 1) & (j["r"] != 1)]
        nj = len(j)
        conc = float(((j["fc"] > 1) == (j["r"] > 1)).mean()) if nj else None
        n_pi, n_hv = int((h["ID"] == "PI-ME/CFS").sum()), int((h["ID"] == "HV").sum())
    ok = conc is not None and conc >= 0.95
    tab = pd.DataFrame({
        "supplementary_table": "Supplementary Data 14A", "table_id": "14A", "subgroup": "all",
        "analyte_id": blk["CHEM_ID"].map(_s).to_numpy(), "analyte_id_type": "metabolon_chem_id",
        "analyte_name": blk["CHEMICAL_NAME"].map(_s).to_numpy(), "hmdb_id": blk["HMDB"].map(_s).to_numpy(),
        "kegg_id": blk["KEGG"].map(_s).to_numpy(),
        "pathway": (blk["SUPER_PATHWAY"].map(_s) + " / " + blk["SUB_PATHWAY"].map(_s)).to_numpy(),
        "effect_measure": "fold change (PI-ME/CFS / HV)", "effect_value": fc.to_numpy(),
        "p_value": p.to_numpy(), "adj_p_value": q.to_numpy(), "source_row": (blk.index + 1).to_numpy(),
    })
    tab["direction_check_concordance"] = conc
    tab["direction"] = np.where(not ok, UNKNOWN, np.where(tab["effect_value"] > 1, "up_in_PI-ME/CFS",
                                np.where(tab["effect_value"] < 1, "down_in_PI-ME/CFS", "no_change")))
    tab["direction_basis"] = (f"published fold change; orientation PI-ME/CFS/HV verified against Figure 6H "
                              f"group means: {conc:.3f} of {nj} metabolites concordant" if ok else
                              "published fold change; orientation NOT verified")
    tab["n_cases"], tab["n_controls"] = n_pi, n_hv
    tab["n_basis"] = "counted group labels in Figure 6H source data" if n_pi else "stated"
    tab["statistical_test"] = "two-sided Wilcoxon (MOESM3)"
    n_q_lt_p = int((q < p).sum())
    rho = float(p.rank().corr(q.rank())) if len(p) > 2 else float("nan")
    j6 = SOURCE_DATA_DIR / "Figure 6J.xlsx"
    n_measured = (len([c for c in pd.read_excel(j6, nrows=0).columns
                       if c not in ("CLIENT_SAMPLE_ID", "ID", "Birth Sex")]) if j6.exists() else None)
    tab["notes"] = (f"published 'adjusted p-value' is smaller than 'p-value' in {n_q_lt_p} of {len(tab)} rows and its "
                    f"rank correlation with 'p-value' is {rho:.2f}, so it cannot be a monotone adjustment of that "
                    "p-value column; fdr_significant is UNKNOWN for this table. The table is titled 'Differential "
                    f"metabolites' and lists {len(tab)} metabolites, all with published adjusted p < 0.05"
                    + (f", of {n_measured} metabolites in the Fig. 6J source data" if n_measured else "")
                    + " (a selected subset, not all tested analytes)")
    m = TABLE_META["14"]
    tab["analyte_class"], tab["tissue"], tab["platform"] = m["analyte_class"], m["tissue"], m["platform"]
    tab["underlying_data_accession"] = m["accession"]
    checks = [{"check": "SD14A fold-change orientation vs Figure 6H group means", "n_checked": nj,
               "n_agree": int(round((conc or 0) * nj))},
              {"check": "SD14A rows where published adjusted p < raw p (should be 0 for a p-value adjustment)",
               "n_checked": int(len(tab)), "n_agree": n_q_lt_p}]
    return [tab], checks


def parse_te_table() -> list[pd.DataFrame]:
    d = _sd(18)
    if d is None:
        return []
    hrow = next(i for i in range(len(d)) if _s(d.iat[i, 0]) == "GeneID")
    blk = d.iloc[hrow + 1:, :7].copy()
    blk.columns = [_s(x) for x in d.iloc[hrow, :7].tolist()]
    blk = blk[_num(blk["p-value"]).notna()]
    tab = pd.DataFrame({
        "supplementary_table": "Supplementary Data 18", "table_id": "18", "subgroup": "all",
        "analyte_id": blk["GeneID"].map(_s).to_numpy(), "analyte_id_type": "TEtranscripts feature",
        "analyte_name": blk["GeneID"].map(_s).to_numpy(), "effect_measure": "DESeq2 log2FoldChange",
        "effect_value": _num(blk["log2FoldChange"]).to_numpy(), "p_value": _num(blk["p-value"]).to_numpy(),
        "adj_p_value": _num(blk["adjusted p-value"]).to_numpy(), "source_row": (blk.index + 1).to_numpy(),
    })
    tab["direction"] = UNKNOWN
    tab["direction_basis"] = ("DESeq2 log2FoldChange published without the contrast orientation and without "
                              "per-sample TE counts, so the sign convention cannot be verified")
    tab["n_cases"], tab["n_controls"] = _stated("18")
    tab["n_basis"] = "stated in MOESM3"
    tab["statistical_test"] = "DESeq2 Wald test (MOESM3)"
    m = TABLE_META["18"]
    tab["analyte_class"], tab["tissue"], tab["platform"] = m["analyte_class"], m["tissue"], m["platform"]
    tab["underlying_data_accession"] = m["accession"]
    return [tab]


def somalogic_value_scale(seqids, values: pd.DataFrame, ref_seqids, ref_vals: np.ndarray | None,
                          max_features: int = 60) -> tuple[str, dict]:
    """Decide whether a table's SomaLogic values are raw RFU or log2(RFU).

    Each published value is compared with the nearest deposited ADAT RFU for the same SOMAmer, under both
    hypotheses (value = RFU, 2**value = RFU); the hypothesis with a median relative error below 5% wins.
    Without an ADAT, a table whose median value is below 25 is taken as log2 (flagged as heuristic).
    """
    x = values.apply(pd.to_numeric, errors="coerce").to_numpy(float)
    info: dict = {}
    if ref_vals is not None and ref_seqids is not None and len(ref_seqids):
        col = {s: j for j, s in enumerate(ref_seqids)}
        errs = {"RFU": [], "log2_RFU": []}
        used = 0
        for i, sid in enumerate(seqids):
            j = col.get(sid)
            if j is None:
                continue
            ref = ref_vals[:, j]
            ref = ref[np.isfinite(ref) & (ref > 0)]
            v = x[i][np.isfinite(x[i])]
            if not len(ref) or not len(v):
                continue
            for scale, y in (("RFU", v), ("log2_RFU", np.exp2(np.clip(v, -60, 60)))):
                k = np.abs(ref[None, :] - y[:, None]).argmin(axis=1)
                errs[scale].append(float(np.median(np.abs(ref[k] - y) / ref[k])))
            used += 1
            if used >= max_features:
                break
        if used:
            med = {s: float(np.median(e)) for s, e in errs.items()}
            info = {"n_features_compared": used, **{f"median_rel_error_{s}": round(v, 4) for s, v in med.items()}}
            best = min(med, key=med.get)
            return (best if med[best] < 0.05 else "unknown"), info
    finite = x[np.isfinite(x)]
    if finite.size:
        info = {"median_value": float(np.median(finite)), "basis": "heuristic (no ADAT to compare)"}
        return ("log2_RFU" if np.median(finite) < 25 else "RFU"), info
    return "unknown", info


def parse_individual_tables(seq_annot: pd.DataFrame | None,
                            adat_ref: dict[str, tuple] | None = None) -> tuple[list[pd.DataFrame], list[dict]]:
    """Supplementary Data 15, 17, 20, 21: group-labelled individual values + published p-values.

    adat_ref maps tissue -> (SeqId list, RFU matrix of the deposited ADAT) and is used to establish whether a
    Supplementary Data 17 block is on the RFU or log2(RFU) scale (17A serum is RFU, 17B CSF is log2 RFU).
    """
    out, checks = [], []
    for n in (15, 17, 20, 21):
        d = _sd(n)
        if d is None:
            continue
        cells = d.map(_s)
        hdr_rows = [i for i in range(len(d)) if int(cells.iloc[i].isin(GROUPS).sum()) >= 3]
        title_rows = [i for i in range(len(d)) if cells.iat[i, 0].startswith("Supplementary Data")]
        for k, h in enumerate(hdr_rows):
            header = cells.iloc[h].tolist()
            nxt = [i for i in title_rows + hdr_rows if i > h]
            end = min(nxt) if nxt else len(d)
            prev = [i for i in title_rows if i < h]
            title = cells.iat[prev[-1], 0] if prev else f"Supplementary Data {n}"
            tid = _table_id(title, str(n))
            hv = [j for j, x in enumerate(header) if x == "HV"]
            pi = [j for j, x in enumerate(header) if x == "PI-ME/CFS"]
            pcol = next((j for j, x in enumerate(header) if x.lower() in ("p-value", "pval")), None)
            if pcol is None:
                continue
            qcol = next((j for j, x in enumerate(header) if "fdr" in x.lower() or "adjusted" in x.lower()), None)
            blk = d.iloc[h + 1:end]
            names = blk.iloc[:, 0].map(_s)
            p = _num(blk.iloc[:, pcol])
            keep = (names != "") & p.notna()
            blk, names, p = blk[keep], names[keep], p[keep]
            vhv = blk.iloc[:, hv].apply(pd.to_numeric, errors="coerce")
            vpi = blk.iloc[:, pi].apply(pd.to_numeric, errors="coerce")
            mhv, mpi = vhv.median(axis=1), vpi.median(axis=1)
            meta = TABLE_META[str(n)]
            tissue = meta["tissue"]
            if tissue is None:
                tl = title.lower()
                tissue = "cerebrospinal fluid" if "cerebrospinal" in tl else ("serum" if "serum" in tl else "blood")
            scale = {15: "percent_or_ratio (as published)", 20: "uM (as published)", 21: "uM (as published)"}.get(n, "")
            if n == 17:
                seq = blk.iloc[:, 0].map(_s).str.replace(r"^SeqId\.", "", regex=True)
                ref = (adat_ref or {}).get(tissue, (None, None))
                scale, sinfo = somalogic_value_scale(seq.tolist(), blk.iloc[:, hv + pi], ref[0], ref[1])
                checks.append({"check": f"SD{tid}: value scale of the individual aptamer values vs the deposited "
                                        f"{tissue} ADAT RFUs", "n_checked": int(sinfo.get("n_features_compared", 0)),
                               "n_agree": int(scale != "unknown"), "detail": f"scale={scale}; {sinfo}"})
            with np.errstate(divide="ignore", invalid="ignore"):
                if scale == "log2_RFU":
                    l2 = (mpi - mhv).to_numpy()
                    measure = ("median PI-ME/CFS minus median HV of the published log2(RFU) values (log2 ratio "
                               "of group medians), computed from the table's individual values")
                else:
                    l2 = np.where((mhv > 0) & (mpi > 0), np.log2(mpi / mhv), np.nan)
                    measure = "log2(median PI-ME/CFS / median HV), computed from the table's individual values"
                    if scale == "unknown":
                        l2 = np.full(len(mhv), np.nan)
                        measure = "not computed: value scale of the published individual values is unknown"
            direction = np.where(mpi > mhv, "up_in_PI-ME/CFS", np.where(mpi < mhv, "down_in_PI-ME/CFS",
                                 "no_difference_in_medians"))
            tab = pd.DataFrame({
                "supplementary_table": f"Supplementary Data {tid}", "table_id": tid, "subgroup": "all",
                "analyte_name": names.to_numpy(), "effect_measure": measure, "value_scale": scale,
                "effect_value": l2, "p_value": p.to_numpy(),
                "adj_p_value": _num(blk.iloc[:, qcol]).to_numpy() if qcol is not None else np.nan,
                "direction": direction,
                "direction_basis": "computed: group medians of the individual values published in the same table",
                "median_pi_mecfs": mpi.to_numpy(), "median_hv": mhv.to_numpy(),
                "n_cases": len(pi), "n_controls": len(hv),
                "n_basis": "counted group-labelled value columns in the table",
                "source_row": (blk.index + 1).to_numpy(),
            })
            tab["analyte_class"], tab["platform"] = meta["analyte_class"], meta["platform"]
            tab["tissue"] = tissue
            tab["underlying_data_accession"] = meta["accession"] or (
                "GSE251790" if tissue == "cerebrospinal fluid" else "GSE254030")
            if n == 17:
                tab["analyte_id"], tab["analyte_id_type"] = seq.to_numpy(), "somalogic_seqid"
                tab["somalogic_seqid"] = seq.to_numpy()
                tab["analyte_name"] = blk.iloc[:, 2].map(_s).to_numpy()
                if seq_annot is not None:
                    a = seq_annot.set_index("SeqId")
                    tab["uniprot_id"] = seq.map(a["UniProt"]).to_numpy()
                    tab["gene_symbol"] = seq.map(a["EntrezGeneSymbol"]).to_numpy()
                    tab["entrez_gene_id"] = seq.map(a["EntrezGeneID"]).to_numpy()
                tab["statistical_test"] = "two-sided Wilcoxon (MOESM3)"
            else:
                tab["analyte_id"] = names.to_numpy()
                tab["analyte_id_type"] = {15: "flow_cytometry_population_label", 20: "lipid_class",
                                          21: "metabolite_name"}[n]
                tab["statistical_test"] = ("unadjusted two-sided t-test (MOESM3)" if n == 15
                                           else "two-sided Wilcoxon (MOESM3)")
            out.append(tab)
    return out, checks


def build_condition_evidence() -> tuple[pd.DataFrame, list[dict]]:
    annots = []
    adat_ref = {}
    for mid in ("csf_somalogic", "serum_somalogic"):   # the two runs use slightly different SOMAmer menus
        adp = RAW / "geo" / GEO_SERIES[mid]["suppl"][0]
        if adp.exists():
            meta, colmeta, vals = parse_adat(adp)
            annots.append(colmeta[["SeqId", "UniProt", "EntrezGeneSymbol", "EntrezGeneID"]])
            keep = (meta["SampleType"] == "Sample").to_numpy()
            adat_ref["cerebrospinal fluid" if mid == "csf_somalogic" else "serum"] = (
                colmeta["SeqId"].tolist(), vals[keep])
    seq_annot = pd.concat(annots).drop_duplicates("SeqId") if annots else None
    frames, checks = [], []
    f, c = parse_de_tables()
    frames += f
    checks += c
    f, c = parse_metabolomics_table()
    frames += f
    checks += c
    frames += parse_te_table()
    f, c = parse_individual_tables(seq_annot, adat_ref)
    frames += f
    checks += c
    if not frames:
        return pd.DataFrame(), checks
    ev = pd.concat(frames, ignore_index=True, sort=False)
    ev["n_cases_stated"] = ev["table_id"].map(lambda t: STATED_N.get(t, (np.nan, np.nan))[0])
    ev["n_controls_stated"] = ev["table_id"].map(lambda t: STATED_N.get(t, (np.nan, np.nan))[1])
    for col in ("gene_symbol", "ensembl_gene_id", "entrez_gene_id", "uniprot_id", "hmdb_id", "kegg_id",
                "somalogic_seqid", "pathway", "notes", "median_pi_mecfs", "median_hv", "direction_check_concordance",
                "value_scale"):
        if col not in ev:
            ev[col] = None
    ev["condition_id"] = CONDITION_ID
    ev["condition_label"] = CONDITION_LABEL
    ev["condition_ontology_id"], ev["condition_ontology_source"] = _condition_ontology_id()
    ev["comparison"] = "PI-ME/CFS vs healthy volunteers (HV)"
    ev["nominal_significant"] = ev["p_value"] < 0.05
    # A published 'adjusted p' is used for significance only when it behaves like a multiple-testing
    # adjustment of the table's own p column (>= p and non-decreasing in p). Supplementary Data 14A fails
    # this (adj < p in 70/142 rows, rank correlation ~0.2), so its FDR flag is unknown, not True.
    consistent = {t: adj_p_consistent_with_p(g["p_value"], g["adj_p_value"]) for t, g in ev.groupby("table_id")}
    ev["adj_p_consistent_with_p"] = pd.array([consistent[t] for t in ev["table_id"]], dtype="boolean")
    usable = ev["adj_p_consistent_with_p"].fillna(False).to_numpy(bool) & ev["adj_p_value"].notna().to_numpy()
    ev["fdr_significant"] = pd.array([bool(q < 0.05) if u else pd.NA for q, u in zip(ev["adj_p_value"], usable)],
                                     dtype="boolean")
    ev["adj_p_method"] = np.select(
        [ev["adj_p_value"].isna().to_numpy(), usable],
        ["no adjusted p-value published", "as published (FDR per MOESM3); consistent with the p-value column"],
        default="as published, but NOT a monotone adjustment of the table's p-value column; not used for "
                "fdr_significant (UNKNOWN)")
    # tables listing only selected analytes: every row nominally significant, or every row below 0.05 on the
    # published adjusted p (SD14A 'Differential metabolites': 142 of the metabolites measured)
    g = ev.groupby("table_id")
    selected = (g["p_value"].transform("max") <= 0.05) | (g["adj_p_value"].transform("max") < 0.05)
    ev["table_scope"] = np.where(selected, "significant_results_only", "all_tested_analytes")
    ev["pmid"], ev["doi"], ev["pmcid"] = PAPER["pmid"], PAPER["doi"], PAPER["pmcid"]
    ev["study_population"] = STUDY_POPULATION
    ev["source_file"] = "nature_esm/41467_2024_45107_MOESM4_ESM.zip::" + ev["supplementary_table"] + ".xlsx"
    ev["evidence_level_value"] = np.where(ev["fdr_significant"].fillna(False).to_numpy(bool),
                                          "single_cohort_fdr_significant",
                                          np.where(ev["nominal_significant"], "single_cohort_nominal_only",
                                                   "single_cohort_not_significant"))
    checks += [{"check": f"SD{t}: published adjusted p is >= p and non-decreasing in p (usable as FDR)",
                "n_checked": int(ev.loc[ev["table_id"] == t, "adj_p_value"].notna().sum()),
                "n_agree": int(bool(ok)), "detail": "consistent" if ok else "NOT consistent: fdr_significant set "
                "to UNKNOWN for this table"} for t, ok in consistent.items() if ok is not None]
    ev = add_common_evidence_columns(ev)
    return ev, checks


# Source label shared with omics.query's display view of this partition.
EVIDENCE_SOURCE_DATABASE = "mapMECFS / NIH PI-ME/CFS publication supplement"


def _registry_primary() -> dict:
    """me_cfs primary Mondo id / label / registry predicate from condition_registry (never typed in)."""
    try:
        reg = pd.read_parquet(PROCESSED / "condition_registry.parquet",
                              columns=["canonical_condition_id", "primary_mondo_id", "primary_mondo_label",
                                       "primary_mondo_relation"])
        r = reg[reg["canonical_condition_id"] == CONDITION_ID]
        if len(r) and _s(r["primary_mondo_id"].iloc[0]):
            return {"id": _s(r["primary_mondo_id"].iloc[0]), "label": _s(r["primary_mondo_label"].iloc[0]),
                    "relation": _s(r["primary_mondo_relation"].iloc[0]) or UNKNOWN}
    except Exception:
        pass
    return {"id": UNKNOWN, "label": UNKNOWN, "relation": UNKNOWN}


def add_common_evidence_columns(ev: pd.DataFrame) -> pd.DataFrame:
    """Give the published-analyte rows the common condition_molecular_evidence columns.

    Mirrors omics.query.EVIDENCE_COLUMNS so the union needs no per-source view:
    entity_type 'gene' for RNA-seq transcripts and 'analyte' for everything else (proteins keep
    their gene symbol in gene_symbol); entity_id = Ensembl gene id, else UniProt, HMDB, KEGG, else
    the published analyte id; entity_label = gene symbol, else analyte name; ontology_id = the
    registry's primary Mondo id for me_cfs; source_accession = the supplementary table of the paper
    (the deposit behind it stays in underlying_data_accession).
    """
    e = ev.copy()
    reg = _registry_primary()
    e["ontology_id"] = reg["id"]
    e["ontology_id_role"] = "primary" if reg["id"] != UNKNOWN else UNKNOWN
    e["ontology_match"] = reg["relation"]
    e["source_disease_id"] = ""   # a published cohort comparison, not an ontology-indexed source
    e["source_disease_label"] = CONDITION_LABEL
    e["source_evidence_category"] = e["supplementary_table"].astype(str) + ": " + e["analyte_class"].astype(str)
    e["entity_type"] = np.where(e["analyte_class"].astype(str) == "transcript", "gene", "analyte")
    eid = pd.Series([None] * len(e), index=e.index, dtype=object)
    for col in ("ensembl_gene_id", "uniprot_id", "hmdb_id", "kegg_id", "analyte_id"):
        if col in e:
            vals = e[col].map(_s)
            eid = eid.where(eid.notna(), vals.where(vals != "", None))
    e["entity_id"] = eid
    sym = e["gene_symbol"].map(_s) if "gene_symbol" in e else pd.Series("", index=e.index)
    e["entity_label"] = sym.where(sym != "", e["analyte_name"].map(_s))
    e["source_database"] = EVIDENCE_SOURCE_DATABASE
    e["source_accession"] = "PMID:" + PAPER["pmid"] + " " + e["supplementary_table"].astype(str)
    direction = e["direction"].map(_s)
    e["evidence_direction"] = direction.where(~direction.isin(["", UNKNOWN]), None)
    e["sample_size"] = (pd.to_numeric(e["n_cases"], errors="coerce")
                        + pd.to_numeric(e["n_controls"], errors="coerce")).astype("float64")
    e["sample_size_basis"] = "participants in the published comparison (n_cases + n_controls; see n_basis)"
    e["source_score"] = np.nan    # the study publishes statistics, not a score
    e["source_score_label"] = None
    e["date_retrieved"] = _retrieved("nature_esm/41467_2024_45107_MOESM4_ESM.zip")
    return e


def adj_p_consistent_with_p(p, q) -> bool | None:
    """True when q could be a multiple-testing adjustment of p: q >= p and q non-decreasing in p.

    None when no adjusted p-values are present.
    """
    d = pd.DataFrame({"p": pd.to_numeric(pd.Series(p), errors="coerce").to_numpy(),
                      "q": pd.to_numeric(pd.Series(q), errors="coerce").to_numpy()}).dropna()
    if d.empty:
        return None
    d = d.sort_values(["p", "q"], kind="mergesort")
    tol = 1e-9
    return bool((d["q"] >= d["p"] - tol).all() and (d["q"].diff().dropna() >= -tol).all())


def _condition_ontology_id() -> tuple[str, str]:
    """ME/CFS ontology id from the project's condition_registry (never typed from memory)."""
    p = PROCESSED / "condition_registry.parquet"
    try:
        reg = pd.read_parquet(p, columns=["canonical_condition_id", "primary_mondo_id"])
        v = reg.loc[reg["canonical_condition_id"] == CONDITION_ID, "primary_mondo_id"].map(_s)
        if len(v) and v.iloc[0]:
            return v.iloc[0], "condition_registry.primary_mondo_id"
    except Exception:
        pass
    return UNKNOWN, UNKNOWN


# =========================================================================== provenance + writing
def _prov(df: pd.DataFrame, *, layer: str, etype: str, rid, version: str, retrieved: str, notes: str,
          level=None) -> pd.DataFrame:
    return add_provenance(df, data_layer=layer, source_name=f"NIH intramural PI-ME/CFS deep phenotyping "
                          f"(Walitt et al. 2024) public deposits [{SOURCE_ID}]", source_version=version,
                          retrieved_at=retrieved, evidence_type=etype, source_record_id=rid,
                          source_geographic_resolution="none", evidence_level=level, provenance_notes=notes)


def _geo_version() -> str:
    dates = set()
    for spec in GEO_SERIES.values():
        for m in spec["matrices"]:
            p = RAW / "geo" / m
            if p.exists():
                with gzip.open(p, "rt", errors="replace") as fh:
                    for line in fh:
                        if line.startswith("!Series_last_update_date"):
                            dates.add(line.split("\t")[1].strip().strip('"'))
                            break
    return f"GEO SuperSeries {GEO_SUPERSERIES} SubSeries; last_update_date {', '.join(sorted(dates)) or UNKNOWN}"


def write_tables(parts, omics, mods, ovl, ev) -> dict[str, int]:
    written = {}
    geo_ret = _retrieved("geo/GSE251872_RAW.tar")
    if not parts.empty:
        p = parts.copy()
        p.insert(0, "participant_id", p["study_number"].map(participant_id))
        p["dataset_id"] = DATASET_ID
        p = _prov(p, layer="person", etype="person_derived_feature", rid="native_id", version=_geo_version(),
                  retrieved=geo_ret,
                  notes="one row per NIH study number seen in an openly downloadable omics deposit; native_id is the "
                        "NIH study number (MECFS_1xx healthy volunteer / MECFS_3xx PI-ME/CFS), NOT a mapMECFS "
                        "'map000NNN' ID; no wearable/actigraphy/HRV data are linked to these participants")
        write_table(p, "participants__mapmecfs", producer=PRODUCER,
                    description="NIH PI-ME/CFS participants identified by study number across open omics deposits")
        written["participants__mapmecfs"] = len(p)
    if not omics.empty:
        o = omics.copy()
        o.insert(0, "participant_id", o["study_number"].map(participant_id))
        o["source_record_id_tmp"] = o["sample_accession"].astype(str) + ":" + o["feature_id"].astype(str)
        # the PBMC count files list 25 Ensembl ids twice (X/Y PAR genes, e.g. CD99); keep both rows as provided
        # and make the record id unique with an occurrence suffix
        occ = o.groupby("source_record_id_tmp", sort=False).cumcount().to_numpy()
        o["source_record_id_tmp"] = np.where(occ > 0, o["source_record_id_tmp"] + "#" + (occ + 1).astype(str),
                                             o["source_record_id_tmp"])
        for col in ("modality", "tissue", "assay", "value_unit", "linkage_method", "source_file", "source_accession"):
            o[col] = o[col].astype("category")
        o = _prov(o, layer="person", etype="person_lab_measurement", rid="source_record_id_tmp",
                  version=_geo_version() + "; Nature Source Data Fig. 6J", retrieved=geo_ret,
                  notes="participant-linked omics: modalities joined only on NIH study numbers published in the "
                        "deposits (serum SomaLogic via the Fatigue-NNN crosswalk published in GSE251790); values as "
                        "provided by the depositors",
                  level=o["linkage_method"].astype(str))
        o = o.drop(columns=["source_record_id_tmp"])
        write_table(o, "participant_omics_linked__mapmecfs", producer=PRODUCER,
                    description="long-format participant-linked omics values (RNA-seq counts, SomaLogic RFU, "
                                "CSF metabolomics) for NIH PI-ME/CFS study numbers")
        written["participant_omics_linked__mapmecfs"] = len(o)
    if mods:
        inv = pd.DataFrame([{**{k: v for k, v in asdict(m).items() if k not in ("native_ids", "study_numbers")},
                             "n_native_ids": len(m.native_ids), "n_study_numbers": len(m.study_numbers),
                             "native_ids": ",".join(map(str, m.native_ids)),
                             "study_numbers": ",".join(m.study_numbers),
                             "participant_level_ids": m.participant_level} for m in mods])
        inv = _prov(inv, layer="metadata", etype="metadata_catalog", rid="modality_id",
                    version=f"audit built {utc_now_iso()}", retrieved=_retrieved("mapmecfs/catalog_group_pi_mecfs.json"),
                    notes="modality inventory for the NIH PI-ME/CFS study: identifiers per modality, open vs gated")
        write_table(inv, "modality_inventory__mapmecfs", producer=PRODUCER,
                    description="per-modality inventory: deposit, open/gated, identifier namespace, subjects")
        written["modality_inventory__mapmecfs"] = len(inv)
    if not ovl.empty:
        ov = ovl.copy()
        ov = _prov(ov, layer="metadata", etype="metadata_catalog",
                   rid=lambda d: d["modality_a"] + "|" + d["modality_b"], version=f"audit built {utc_now_iso()}",
                   retrieved=geo_ret, notes="exact participant-identifier overlap per modality pair")
        write_table(ov, "modality_overlap__mapmecfs", producer=PRODUCER,
                    description="pairwise exact identifier overlap between NIH PI-ME/CFS modalities")
        written["modality_overlap__mapmecfs"] = len(ov)
    if not ev.empty:
        e = ev.copy()
        e["source_record_id_tmp"] = e["supplementary_table"] + ":" + e["analyte_id"].astype(str) + ":row" + \
            e["source_row"].astype(str)
        e = _prov(e, layer="condition_molecular", etype="published_biomarker", rid="source_record_id_tmp",
                  version=f"{PAPER['citation']}; Supplementary Data 1-24 (MOESM4)",
                  retrieved=_retrieved("nature_esm/41467_2024_45107_MOESM4_ESM.zip"),
                  notes="condition-level molecular enrichment from one exploratory case-control cohort; not "
                        "participant-linked to any wearable data; exploratory single-study results",
                  level=e["evidence_level_value"])
        e = e.drop(columns=["source_record_id_tmp", "evidence_level_value"])
        write_table(e, "condition_molecular_evidence__mapmecfs", producer=PRODUCER,
                    description="published differential analytes/genes (Walitt 2024 Supplementary Data 14-21)")
        written["condition_molecular_evidence__mapmecfs"] = len(e)
    stale = PROCESSED / "participant_wearable_features__mapmecfs.parquet"
    if stale.exists():  # never leave a wearable partition behind when none can be linked
        stale.unlink()
        stale.with_suffix(".meta.json").unlink(missing_ok=True)
    return written


# =========================================================================== documents
def _md_table(df: pd.DataFrame, cols: list[str] | None = None, max_rows: int | None = None) -> str:
    if df is None or df.empty:
        return "_none_\n"
    d = df[cols] if cols else df
    if max_rows:
        d = d.head(max_rows)
    def esc(v) -> str:
        if isinstance(v, (float, np.floating)) and not pd.isna(v) and float(v).is_integer():
            v = int(v)
        return _s(v).replace("|", "\\|").replace("\n", " ")
    lines = ["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)]
    lines += ["| " + " | ".join(esc(v) for v in r) + " |" for r in d.itertuples(index=False)]
    return "\n".join(lines) + "\n"


def _overlap_matrix(mods: list[Modality], ovl: pd.DataFrame) -> pd.DataFrame:
    ids = [m.modality_id for m in mods]
    mat = pd.DataFrame("", index=ids, columns=ids)
    lookup = {(r.modality_a, r.modality_b): r for r in ovl.itertuples()}
    for m in mods:
        mat.loc[m.modality_id, m.modality_id] = (str(m.n_subjects) if m.participant_level else
                                                 ("no IDs" if m.obtained else "not obtained"))
    for a in ids:
        for b in ids:
            if a == b:
                continue
            r = lookup.get((a, b)) or lookup.get((b, a))
            if r is None:
                continue
            if r.comparable:
                mat.loc[a, b] = str(int(r.n_overlap))
            else:
                mat.loc[a, b] = "n/o" if str(r.reason).startswith("not obtained") else "n/c"
    mat.insert(0, "modality", ids)
    return mat


MATRIX_GROUP_ONLY = ("hrv_24h_time_domain", "orthostatic_tilt_catecholamines", "cpet_peak_vo2", "grip_strength")
MATRIX_GATED_LABELS = ("Free Living Accelerometry", "Heart Rate Variability Dataset", "Head-up Tilt Table Testing",
                       "Cardiopulmonary Exercise Tests (CPET)")


def _check_table(checks: list[dict]) -> pd.DataFrame:
    return pd.DataFrame([{"check": c["check"], "n_checked": c.get("n_checked"), "n_agree_or_count": c.get("n_agree"),
                          "detail": c.get("detail", "")} for c in checks])


def render_linkage_audit(ctx: dict) -> str:
    mods: list[Modality] = ctx["mods"]
    ovl: pd.DataFrame = ctx["ovl"]
    parts: pd.DataFrame = ctx["parts"]
    log: list[AccessRecord] = ctx["log"]
    obtained = [m for m in mods if m.obtained]
    pl = [m for m in obtained if m.participant_level]
    sn_mods = [m for m in pl if m.study_numbers]
    gated = [m for m in mods if not m.obtained]
    wear_obtained = [m for m in obtained if m.is_wearable]
    union = sorted(set().union(*[set(m.study_numbers) for m in sn_mods])) if sn_mods else []
    wear_omics_pairs = ovl[(ovl["a_is_wearable"] & ovl["b_is_omics"]) | (ovl["b_is_wearable"] & ovl["a_is_omics"])] \
        if not ovl.empty else ovl
    n_wear_omics_comparable = int(wear_omics_pairs["comparable"].sum()) if not wear_omics_pairs.empty else 0
    probes = []
    pp = RAW / "mapmecfs" / "access_probe.json"
    if pp.exists():
        probes = json.loads(pp.read_text())
    file_probes = [p for p in probes if p["label"].startswith("file download")]
    n_file_probes = len(file_probes)
    n_probe_403 = sum(1 for p in file_probes if p.get("http_status") == 403 and p.get("login_required"))
    api_403 = [p for p in probes if "API" in p["label"] and p.get("http_status") == 403]

    L = [f"# mapMECFS / NIH PI-ME/CFS participant-linkage audit\n",
         f"Generated by `{PRODUCER}` at {utc_now_iso()} from files retrieved into `data/raw/{SOURCE_ID}/`. "
         f"Every number below is computed from those files. Study: {PAPER['citation']} "
         f"(PMID {PAPER['pmid']}, doi:{PAPER['doi']}, {PAPER['trial']}).\n"]
    L.append("## Bottom line\n")
    if sn_mods:
        L.append(f"* **Wearable/actigraphy, HRV, orthostatic and CPET data cannot truthfully be linked to omics "
                 f"with openly downloadable data.** In every open file these measurements carry group labels only "
                 f"(HV / PI-ME/CFS) or figure-local letters, never a participant identifier "
                 f"({len(wear_obtained)} wearable-type modalities obtained, all without IDs; "
                 f"{n_wear_omics_comparable} wearable x omics pairs have comparable identifiers). The "
                 f"per-participant files (free-living accelerometry, HRV, tilt, CPET) exist only on mapMECFS, "
                 f"which requires an account.\n")
        L.append(f"* **Omics x omics linkage is real but limited to the NIH study number.** {len(sn_mods)} open omics "
                 f"modalities carry NIH study numbers (1xx = healthy volunteer, 3xx = PI-ME/CFS), directly or via the "
                 f"SomaLogic crosswalk published in GSE251790. Union: {len(union)} study numbers "
                 f"({sum(1 for s in union if s[0] == '1')} with 1xx, {sum(1 for s in union if s[0] == '3')} with 3xx). "
                 f"The paper's analytic cohort is 21 HV + 17 PI-ME/CFS, so the deposits include participants "
                 f"outside it, and no open file flags which.\n")
    else:
        L.append("* No participant-level identifiers could be obtained in this run (see access log). Molecular "
                 "results are condition-level evidence only.\n")
    L.append(f"* **mapMECFS**: the CKAN `package_search` catalogue metadata is public, but {n_probe_403} of "
             f"{n_file_probes} anonymous file-download probes returned HTTP 403 with the page text \"Login Required. "
             f"Please login to view this feature.\" ({len(api_403)} CKAN API actions also returned 403 "
             f"'Authentication is required'; the group landing page answers HTTP 200 but shows the same login "
             f"notice.) The paper states: \"Accessing Map ME/CFS data requires signing up for an account but "
             f"otherwise access to the data is unrestricted (Creative Commons BY 4.0).\" No account was created, so "
             f"{len(gated)} mapMECFS data files are catalogued as not automatable.\n")
    L.append("* Consequence for the engine: the molecular results from this study enter as **condition-level "
             "molecular enrichment** (`condition_molecular_evidence__mapmecfs`). The participant-linked omics table "
             "links omics to omics only; it is not \"patient multi-omics\" with wearable physiology.\n")

    L.append("## 1. Public deposits named in the paper and what was retrieved\n")
    stm = availability_statements()
    if stm:
        L.append("Verbatim Data Availability statement (Europe PMC JATS XML, `paper/PMC10881493_fullTextXML.xml`):\n\n"
                 f"> {stm.get('Data availability', UNKNOWN)}\n\nLinks in that section: "
                 f"{stm.get('Data availability links', '') or UNKNOWN}\n\n")
    L.append("Data Availability (paper) names: the mapMECFS group "
             f"{MAPMECFS_GROUP_URL}; GEO GSE251872 (PBMC RNA-seq), GSE245661 (muscle RNA-seq), GSE251790 and "
             f"GSE254030 (SomaLogic proteomics), all under SuperSeries {GEO_SUPERSERIES}; SRA {SRA_STUDY} / BioProject "
             f"{BIOPROJECT} (stool metagenomics); Pennsieve dataset {PENNSIEVE_ID} (doi:{PENNSIEVE_DOI}; TMS/EMG and "
             f"fMRI); Source Data with the paper; the authors' GitHub repository; and external comparison datasets "
             f"GSE130353 and GSE156792 (other cohorts, not ingested). The text once cites 'GSE13033' but links to "
             f"GSE130353.\n")
    dep = pd.DataFrame([{
        "modality": m.modality_id, "label": m.label, "deposit": m.deposit,
        "open without login": {True: "yes", False: "NO"}.get(m.open_download, "?"),
        "obtained": "yes" if m.obtained else "no", "format": m.file_format, "kind": m.data_kind,
        "subjects with IDs": m.n_subjects if m.participant_level else 0,
        "group-labelled rows/cols (no ID)": m.n_rows_group_labelled_only or "",
        "identifier scheme": m.id_namespace} for m in mods if m.obtained])
    L.append(_md_table(dep))
    L.append("\nAccess log (every request; presigned URLs shown by their stable API path):\n")
    L.append(_md_table(pd.DataFrame([asdict(a) for a in log]),
                       ["deposit", "url", "http_status", "outcome", "detail"]))

    L.append("\n## 2. Participant identifier schemes\n")
    ns = pd.DataFrame([
        {"namespace": STUDY_NUMBER_NS, "format": "MECFS_101 / 'Patient #312' / 'Walitt101' / 'S101' / '103'",
         "meaning": "NIH study number; 1xx healthy volunteer, 3xx PI-ME/CFS", "where": "GEO RNA-seq titles and "
         "count files, GEO CSF SomaLogic titles, Fig. 6J 'ID' column (some rows)"},
        {"namespace": "somalogic_fatigue_code", "format": "Fatigue-101 ... Fatigue-144",
         "meaning": "SomaLogic sample code shared by the CSF and serum runs", "where": "GEO GSE251790 description, "
         "GSE254030 titles, both ADAT files, Figs. S15/S16 source data"},
        {"namespace": "mapmecfs_mapid", "format": "map000288-00-02", "meaning": "mapMECFS participant ID "
         "(middle segment = ME/CFS status)", "where": "GitHub EEfRT file; mapMECFS files (gated, format unverified)"},
        {"namespace": "microbiome_subject_sid", "format": "SID_975", "meaning": "microbiome subject code",
         "where": "BioSample attributes (SRA), GitHub microbiome metadata"},
        {"namespace": "pennsieve_emg_index / pennsieve_fmri_index", "format": "HV1 / MECFS1 / HV01",
         "meaning": "dataset-local index", "where": "Pennsieve 356 file names, Fig. 4 source data"},
        {"namespace": "figure_local_letter:*", "format": "HV A / PI-ME/CFS A", "meaning": "per-figure pseudonym",
         "where": "Fig. 3A EEfRT, Fig. S9E-G Seahorse source data"},
        {"namespace": "group_label_only", "format": "HV / PI-ME/CFS", "meaning": "no participant identifier",
         "where": "all other source-data files and Supplementary Data 1-24"},
        {"namespace": "not_obtained", "format": UNKNOWN, "meaning": "file not accessible without login",
         "where": "mapMECFS data files"},
    ])
    L.append(_md_table(ns))
    L.append("\nNo public file maps one namespace onto another, except the per-sample pairing of 'Fatigue-NNN' "
             "with '#NNN' in GSE251790. Numeric parts of other namespaces (e.g. SID_340, map000288) are NOT study "
             "numbers and are never normalised into them.\n")

    L.append("## 3. Subjects per modality (participant-level identifiers)\n")
    sub = pd.DataFrame([{"modality": m.modality_id, "namespace": m.id_namespace, "n subjects": m.n_subjects,
                         "HV (1xx)": sum(1 for s in m.study_numbers if s[0] == "1") if m.study_numbers else "",
                         "PI-ME/CFS (3xx)": sum(1 for s in m.study_numbers if s[0] == "3") if m.study_numbers else "",
                         "linkage": m.linkage_method, "notes": m.notes} for m in pl])
    L.append(_md_table(sub))

    L.append("\n## 4. Exact identifier overlap between modalities\n")
    L.append("Diagonal = subjects with identifiers; off-diagonal = exact overlap; `n/c` = not comparable (different "
             "namespace with no crosswalk, or no identifiers at all); `n/o` = file not obtained (mapMECFS login). "
             "Full pairwise list: `modality_overlap__mapmecfs`.\n\n")
    core = (pl + [m for m in obtained if m.modality_id in MATRIX_GROUP_ONLY]
            + [m for m in gated if m.label in MATRIX_GATED_LABELS])
    if core and not ovl.empty:
        L.append(_md_table(_overlap_matrix(core, ovl)))
    if sn_mods and not ovl.empty:
        sn_ids = {m.modality_id for m in sn_mods}
        cmp_ = ovl[ovl["modality_a"].isin(sn_ids) & ovl["modality_b"].isin(sn_ids)]
        L.append("\nStudy-number overlaps (comparable pairs):\n\n")
        L.append(_md_table(cmp_, ["modality_a", "modality_b", "n_a", "n_b", "n_overlap", "reason", "overlap_ids"]))
        if len(sn_mods) >= 2:
            inter = set.intersection(*[set(m.study_numbers) for m in sn_mods])
            L.append(f"\nStudy numbers present in all {len(sn_mods)} study-number modalities: {len(inter)} "
                     f"({', '.join(sorted(inter))}).\n")
        if not parts.empty:
            dist = parts["n_open_omics_modalities"].value_counts().sort_index()
            L.append("\nParticipants by number of open omics modalities: " +
                     ", ".join(f"{k} modalities: {v}" for k, v in dist.items()) + ".\n")

    L.append("\n## 5. Consistency checks (measured)\n")
    L.append(_md_table(_check_table(ctx["checks"])))

    L.append("\n## 6. Can wearable/actigraphy and omics be linked?\n")
    L.append("**No, not with openly downloadable data.** Evidence:\n\n")
    for m in [m for m in obtained if m.is_wearable or m.data_kind in ("autonomic_physiology", "exercise_physiology")]:
        L.append(f"* `{m.modality_id}` ({m.deposit}): {m.notes}\n")
    wg = [m for m in gated if m.is_wearable or m.data_kind in ("autonomic_physiology", "exercise_physiology")]
    L.append(f"* {len(wg)} per-participant wearable/autonomic/exercise files are listed on mapMECFS and are gated: "
             + "; ".join(f"{m.label} ({m.file_format}, {re.search(r'(\d+|None) bytes', m.notes).group(0)})" for m in wg)
             + ".\n")
    L.append("* Their identifier scheme is unknown without access. The only mapMECFS-format IDs we could see "
             "(GitHub EEfRT file, 'map000NNN-GG-VV') have no public crosswalk to the NIH study numbers used by the GEO "
             "omics deposits. Even with a mapMECFS account, linkage to GEO would require that crosswalk; "
             "linkage within mapMECFS (e.g. accelerometry to mapMECFS-hosted RNA-seq/SomaLogic tables) might be "
             "possible but is UNKNOWN / NOT AVAILABLE until audited with access.\n")
    L.append("* Therefore `participant_wearable_features__mapmecfs` is not produced.\n")

    L.append("\n## 7. Paper source data and supplementary tables: identifiers\n")
    L.append("Per-participant values in the Nature Source Data count as person-level data only when they carry an "
             "identifier. Result of scanning every file (per-file list: `source_data_identifier_scan.csv`):\n\n")
    L.append(_md_table(pd.DataFrame(ctx.get("source_scan", []))))
    rscan = ctx.get("r_object_scan", [])
    if rscan:
        L.append("\nThe authors' GitHub repository also ships the figure source data as serialised R objects. A string "
                 "scan of the decompressed objects (character data is stored as plain bytes; numeric ID columns "
                 "would not be visible to it) finds:\n\n")
        L.append(_md_table(pd.DataFrame(rscan)))
        linking = [r["file"] for r in rscan if STUDY_NUMBER_NS in str(r["identifier_types"])
                   or "mapmecfs_mapid" in str(r["identifier_types"])]
        L.append("\n" + ("No NIH study number or mapMECFS ID string occurs in these objects, so they add no linkage "
                         "beyond the GEO SomaLogic 'Fatigue-NNN' codes and the SRA 'SID_NNN' codes.\n" if not linking
                         else f"Study-number or mapMECFS-ID strings occur in {linking}; these objects need a "
                              "column-level audit before any linkage claim.\n"))

    L.append("\n## 8. What is blocked and why\n")
    L.append(_md_table(pd.DataFrame(probes), ["label", "url", "http_status", "login_required", "message"])
             if probes else "_mapMECFS was not probed in this run_\n")
    L.append("\n* mapMECFS CKAN resources report `datastore_active=True`. The datastore API was not queried: the "
             "portal requires login for data and routing around that gate would not respect the access "
             "conditions.\n")
    L.append("* Pennsieve 356 EMG/fMRI signal files (~3.2 GB) are open but not mirrored by default; identifiers "
             "come from manifest.json. `--pennsieve-bulk` mirrors them.\n")
    L.append("* SRA raw reads (RNA-seq FASTQ and stool metagenomes) are open but not mirrored; run and BioSample "
             "metadata are.\n")
    L.append(f"* ClinicalTrials.gov {PAPER['trial']} has posted results, but they are group-level summaries with no "
             "participant identifiers.\n")

    L.append("\n## 9. Processed outputs\n")
    L.append(_md_table(pd.DataFrame([{"table": k, "rows": v} for k, v in ctx["written"].items()])))
    L.append("\n`participant_id = \"mapmecfs:MECFS_<study number>\"`. The native ID is the NIH study number, not a "
             "mapMECFS 'map000NNN' ID.\n")

    L.append("\n## 10. Other ME/CFS datasets (future sources, not ingested)\n")
    geo_csv = RAW / "future_sources" / "geo_mecfs_series.csv"
    if geo_csv.exists() and geo_csv.stat().st_size > 5:
        g = pd.read_csv(geo_csv)
        L.append(f"GEO human series matching the ME/CFS query ({len(g)}; openly downloadable = series-matrix "
                 f"directory answered HTTP 200 with at least one matrix file):\n\n")
        L.append(_md_table(g, ["accession", "title", "modality", "n_samples", "openly_downloadable",
                               "title_mentions_me_cfs", "this_nih_cohort"]))
    cat_csv = RAW / "future_sources" / "mapmecfs_catalog.csv"
    if cat_csv.exists() and cat_csv.stat().st_size > 5:
        c = pd.read_csv(cat_csv)
        L.append(f"\nmapMECFS public catalogue: {len(c)} datasets from {c['organization'].nunique()} organisations. "
                 "Uploaded files require login (verified with one anonymous probe per organisation, section 8); "
                 "entries with open external accessions can be fetched from those repositories:\n\n")
        L.append(_md_table(c.groupby("organization").agg(datasets=("dataset_name", "count"),
                                                         uploaded_files=("n_uploaded_files", "sum")).reset_index()))
        ext = c[c["open_external_accessions"].fillna("") != ""]
        L.append("\n" + _md_table(ext, ["title", "organization", "data_type", "open_external_accessions"]))
        L.append("\nFull list: `data/raw/mapmecfs_nih_pi_mecfs/future_sources/mapmecfs_catalog.csv`.\n")

    L.append("\n## 11. Reproduce\n\n`uv run python -m measure_it.ingestion.mapmecfs` (add `--pennsieve-bulk` to "
             "mirror the Pennsieve signal files). Tests: `uv run pytest tests/test_mapmecfs.py`.\n")
    return "".join(x if x.endswith("\n") else x + "\n" for x in L)


def scan_source_data_identifiers() -> list[dict]:
    """Classify every Nature source-data / supplementary file by the identifiers it carries."""
    out = []
    pat = {
        STUDY_NUMBER_NS: re.compile(r"^(MECFS_?\d{3}|[13]\d{2}(\.0)?)$"),
        "somalogic_fatigue_code": re.compile(r"^Fatigue-\d+$"),
        "figure_local_letter": re.compile(r"^(HV|PI-ME/CFS) [A-Z]$"),
        "figure_local_index": re.compile(r"^(HV|MECFS)\d+$"),
        "metabolon_client_id": re.compile(r"^NINR-\d+$"),
        "external_sra_run": re.compile(r"^SRR\d+$"),
    }
    files = sorted(SOURCE_DATA_DIR.glob("*")) + sorted(SUPP_DIR.glob("*.xlsx"))
    for f in files:
        if f.suffix not in (".xlsx", ".csv", ".txt"):
            out.append({"file": f.name, "identifier_types": f"not tabular ({f.suffix})", "n_distinct_ids": ""})
            continue
        try:
            if f.suffix == ".xlsx":
                df = pd.read_excel(f, header=None, nrows=60000)
            else:
                df = pd.read_csv(f, header=None, sep="\t" if f.suffix == ".txt" else ",", nrows=60000,
                                 low_memory=False, on_bad_lines="skip")
        except Exception as exc:
            out.append({"file": f.name, "identifier_types": f"unreadable: {type(exc).__name__}", "n_distinct_ids": ""})
            continue
        # identifiers live in the first two columns and the first three rows (headers)
        cand = set(df.iloc[:, :2].map(_s).to_numpy().ravel()) | set(df.iloc[:3].map(_s).to_numpy().ravel())
        found = {k: sorted(v for v in cand if p.match(v)) for k, p in pat.items()}
        # a 3-digit number in a value column is not an ID: only count study numbers in columns named ID
        if found[STUDY_NUMBER_NS]:
            idcols = [j for j in range(min(3, df.shape[1])) if _s(df.iat[0, j]).upper() == "ID"]
            found[STUDY_NUMBER_NS] = sorted({v for j in idcols for v in df.iloc[1:, j].map(_s) if pat[STUDY_NUMBER_NS].match(v)})
        kinds = {k: len(v) for k, v in found.items() if v}
        out.append({"file": f.name, "identifier_types": ", ".join(f"{k} ({n})" for k, n in kinds.items())
                    or "group label only / none", "n_distinct_ids": sum(kinds.values())})
    agg = pd.DataFrame(out)
    if agg.empty:
        return []
    agg.to_csv(RAW / "source_data_identifier_scan.csv", index=False)
    summary = agg.groupby("identifier_types")["file"].agg(
        lambda s: f"{len(s)} files: " + ", ".join(list(s)[:12]) + (", ..." if len(s) > 12 else "")).reset_index()
    return summary.rename(columns={"file": "files"}).to_dict("records")


R_ID_PATTERNS = {
    STUDY_NUMBER_NS: rb"MECFS_\d{3}\b|Walitt\d{3}\b|(?:Patient|Control) #\d{3}",
    "somalogic_fatigue_code": rb"Fatigue-\d+",
    "microbiome_subject_sid": rb"SID_\d+",
    "mapmecfs_mapid": rb"map\d{6}-\d{2}-\d{2}",
    "metabolon_client_id": rb"NINR-\d+",
}


def scan_github_r_objects() -> list[dict]:
    """Identifier strings inside the authors' serialised R objects (GitHub RData/rds).

    R serialisation stores character data as plain bytes, so a regex over the decompressed stream finds every
    string identifier. Numeric ID columns (doubles) are not visible to this scan.
    """
    rows = []
    for f in sorted(list(GITHUB_DIR.rglob("*.RData")) + list(GITHUB_DIR.rglob("*.rds"))):
        try:
            raw = f.read_bytes()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            elif raw[:3] == b"BZh":
                import bz2
                raw = bz2.decompress(raw)
            elif raw[:6] == b"\xfd7zXZ\x00":
                import lzma
                raw = lzma.decompress(raw)
        except Exception as exc:
            rows.append({"file": f"github: {f.relative_to(GITHUB_DIR)}", "identifier_types": f"unreadable: {exc}",
                         "n_distinct_ids": ""})
            continue
        found = {k: len(set(re.findall(p, raw))) for k, p in R_ID_PATTERNS.items()}
        kinds = {k: n for k, n in found.items() if n}
        rows.append({"file": f"github: {f.relative_to(GITHUB_DIR)}",
                     "identifier_types": ", ".join(f"{k} ({n})" for k, n in kinds.items()) or "none (string scan)",
                     "n_distinct_ids": sum(kinds.values())})
    return rows


def _missingness(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    for c in cols:
        if c not in df:
            continue
        s = df[c]
        miss = s.isna() | (s.astype(str).isin(["", UNKNOWN]))
        rows.append({"column": c, "n": len(s), "missing_or_unknown": int(miss.sum()),
                     "pct": round(100 * miss.mean(), 1) if len(s) else 0.0})
    return pd.DataFrame(rows)


def render_data_audit(ctx: dict) -> str:
    mods, parts, omics, ev = ctx["mods"], ctx["parts"], ctx["omics"], ctx["ev"]
    man = load_manifest(SOURCE_ID)["files"]
    pl = [m for m in mods if m.obtained and m.participant_level]
    gated = [m for m in mods if not m.obtained]
    status = ctx["status"]
    ss = ctx["sample_size"]
    L = [f"# DATA AUDIT — NIH intramural PI-ME/CFS deep phenotyping (mapMECFS + GEO/SRA/Pennsieve deposits)\n\n"]
    fields = [
        ("source_id", SOURCE_ID),
        ("Source (dataset/API name, exact files/endpoints)", "Walitt et al. 2024 Nature Communications supplementary "
         "+ source data; GitHub docwalitt code repository; GEO GSE251872, GSE245661, GSE251790, GSE254030 "
         f"(SuperSeries {GEO_SUPERSERIES}); SRA {SRA_STUDY}/BioProject {BIOPROJECT} metadata; Pennsieve dataset "
         f"{PENNSIEVE_ID} metadata; mapMECFS CKAN catalogue ({MAPMECFS_GROUP_URL}) metadata only"),
        ("Publishing organization", "NIH (NINDS Division of Intramural Research and collaborating institutes); "
         "mapMECFS hosted by RTI International; NCBI GEO/SRA; Pennsieve (Blackfynn Discover)"),
        ("Retrieval date (UTC)", ctx["retrieved"]),
        ("Source version / release", f"paper doi:{PAPER['doi']} (PMID {PAPER['pmid']}); {_geo_version()}; "
         f"GitHub main at commit {ctx.get('github_commit') or UNKNOWN}; Pennsieve 356 version "
         f"{ctx.get('pennsieve_version') or UNKNOWN}"),
        ("Source update date / cadence", "static publication deposits; mapMECFS group metadata_modified "
         f"{ctx.get('mapmecfs_modified') or UNKNOWN}"),
        ("License / access conditions", "Nature article and supplements CC BY 4.0; GEO/SRA public; Pennsieve CC BY; "
         "GitHub GPL-3.0; mapMECFS: CC BY 4.0 but files require a registered account (not created)"),
        ("Unit of observation", "participant x modality x feature (linked omics); analyte x comparison "
         "(condition evidence); modality (inventory)"),
        ("Sample size (actual, as ingested)", "; ".join(f"{k}: {v}" for k, v in ss.items())),
        ("Geography (resolution, vintage)", "none (single-site NIH Clinical Center study; no participant geography)"),
        ("Person-level?", "yes (omics linked by NIH study number)"),
        ("Geographic?", "no"),
        ("Omics?", "yes"),
        ("Wearable?", "no participant-level wearable data obtainable (group-labelled HRV source data only; "
         "actigraphy/HRV files gated on mapMECFS)"),
        ("True participant linkage across modalities?", "omics x omics: yes (NIH study numbers, overlaps counted); "
         "wearable/HRV/orthostatic/CPET x omics: NO (no identifiers in open files)"),
    ]
    L.append(_md_table(pd.DataFrame(fields, columns=["Field", "Value"])))
    L.append(f"\nStatus: **{status}**. Linkage details: `docs/MAP_MECFS_LINKAGE_AUDIT.md`.\n")
    L.append("\n## Files / endpoints retrieved\n\nFrom MANIFEST.json (sha256 there):\n\n")
    L.append(_md_table(pd.DataFrame([{"file": k, "bytes": v.get("bytes"), "retrieved_at": v.get("retrieved_at"),
                                      "url": v.get("url")} for k, v in sorted(man.items())])))
    L.append("\nOther records written: `access_log.json` (every request), `mapmecfs/access_probe.json` (anonymous "
             "mapMECFS probes), `pennsieve/dataset_356.json`, `pennsieve/files_listing.json`, "
             "`future_sources/*.csv`. Zips are extracted under `extracted/`.\n")
    L.append("\n## Key variables\n\n")
    kv = pd.DataFrame([
        ("participant_id", "mapmecfs:MECFS_<NIH study number>", "string", "person key; never joined to other datasets"),
        ("study_number / native_id", "NIH study number (1xx HV, 3xx PI-ME/CFS)", "string", "linkage key"),
        ("group_label", "HV / PI-ME/CFS as deposited", "string", "case-control label"),
        ("modality, feature_id, value", "omics modality, gene/SOMAmer/metabolite, measured value", "counts / RFU / "
         "normalised abundance", "participant-linked omics"),
        ("linkage_method", "direct_study_number / crosswalk_fatigue_code_via_GSE251790", "string",
         "how the sample was attached to a participant"),
        ("analyte_id, analyte_class, tissue", "published differential analyte", "", "condition evidence"),
        ("effect_value, direction, p_value, adj_p_value", "published effect + significance", "logFC / fold change "
         "/ log2 median ratio", "condition evidence"),
        ("n_cases, n_controls", "PI-ME/CFS and HV sample sizes", "count", "counted from files where possible"),
    ], columns=["variable", "meaning", "units", "use downstream"])
    L.append(_md_table(kv))
    L.append("\n## Missingness\n\nMeasured on the processed tables.\n\n")
    if not parts.empty:
        L.append("participants__mapmecfs:\n\n")
        L.append(_md_table(_missingness(parts, ["group_label", "sex", "age_years_somalogic", "somalogic_fatigue_code",
                                                "in_published_analytic_cohort"])))
        L.append("\nOpen omics modalities per participant: " + ", ".join(
            f"has_{m}: {int(parts[f'has_{m}'].sum())}" for m in ("pbmc_rnaseq", "muscle_rnaseq", "csf_somalogic",
                                                                  "serum_somalogic", "csf_metabolomics")) + ".\n")
    if not omics.empty:
        mm = omics.groupby("modality", observed=True).agg(rows=("value", "size"), missing=("value", lambda s: int(s.isna().sum())),
                                                          participants=("study_number", "nunique"),
                                                          features=("feature_id", "nunique")).reset_index()
        L.append("\nparticipant_omics_linked__mapmecfs (value missingness by modality):\n\n")
        L.append(_md_table(mm))
    if not ev.empty:
        em = ev.groupby("supplementary_table").agg(
            rows=("analyte_id", "size"), p_missing=("p_value", lambda s: int(s.isna().sum())),
            adj_p_missing=("adj_p_value", lambda s: int(s.isna().sum())),
            direction_unknown=("direction", lambda s: int((s == UNKNOWN).sum())),
            nominal_p_lt_05=("nominal_significant", "sum"), fdr_lt_05=("fdr_significant", lambda s: int(s.sum()) if s.notna().any() else UNKNOWN),
            n_cases=("n_cases", "first"), n_controls=("n_controls", "first"),
            n_cases_stated=("n_cases_stated", "first"), n_controls_stated=("n_controls_stated", "first"),
        ).reset_index()
        L.append("\ncondition_molecular_evidence__mapmecfs by supplementary table:\n\n")
        L.append(_md_table(em))
    L.append("\n## Linkage strategy\n\n")
    L.append(f"* Joinable: {', '.join(m.modality_id for m in pl if m.study_numbers)} on the NIH study number "
             "(serum SomaLogic through the Fatigue-NNN to study-number pairing published in GSE251790).\n"
             "* Not joinable: stool metagenomics (SID codes), EEfRT (mapMECFS IDs), Pennsieve TMS/fMRI "
             "(dataset-local indices), all group-labelled source data (HRV, CPET, tilt, grip, DEXA, ...), and "
             f"{len(gated)} mapMECFS files (login required).\n"
             "* No person in this source is joined to geographic or facility layers.\n")
    L.append("\n## Limitations and caveats\n\n")
    for lim in ctx["limitations"]:
        L.append(f"* {lim}\n")
    L.append("\n## Processed outputs\n\n")
    L.append(_md_table(pd.DataFrame([{"table": k, "rows": v} for k, v in ctx["written"].items()])))
    L.append("\n## Reproduce\n\n`uv run python -m measure_it.ingestion.mapmecfs`\n")
    return "".join(L)


# =========================================================================== run
def _github_commit() -> str | None:
    z = RAW / "github" / "docwalitt_code_repository_main.zip"
    if not z.exists():
        return None
    try:
        return zipfile.ZipFile(z).comment.decode().strip() or None
    except Exception:
        return None


def _load_access_log() -> list[AccessRecord]:
    p = RAW / "access_log.json"
    if not p.exists():
        return []
    try:
        return [AccessRecord(**r) for r in json.loads(p.read_text())]
    except Exception:
        return []


def run(*, fetch: bool = True, pennsieve_bulk: bool = False) -> dict:
    log = acquire(pennsieve_bulk=pennsieve_bulk) if fetch else _load_access_log()
    geo = {}
    for mid in GEO_SERIES:
        try:
            geo[mid] = geo_samples(mid)
        except Exception as exc:
            log.append(AccessRecord(f"parse {mid}", "", "", None, "error", f"{type(exc).__name__}: {exc}"[:300]))
            geo[mid] = pd.DataFrame()
    crosswalk, checks = somalogic_crosswalk(geo)
    metab = csf_metabolomics_rows()
    mods, c2 = build_modalities(geo, crosswalk, metab)
    checks += c2
    ovl = overlap_table(mods)
    parts, omics, c3 = build_linked(geo, crosswalk, metab)
    checks += c3
    ev, c4 = build_condition_evidence()
    checks += c4
    written = write_tables(parts, omics, mods, ovl, ev)

    pl = [m for m in mods if m.obtained and m.participant_level]
    sn_mods = [m for m in pl if m.study_numbers]
    union = sorted(set().union(*[set(m.study_numbers) for m in sn_mods])) if sn_mods else []
    gated = [m for m in mods if not m.obtained]
    sample_size = {
        "study_numbers_union_open_omics": len(union),
        **{f"subjects_{m.modality_id}": m.n_subjects for m in pl},
        "participant_omics_linked_rows": int(len(omics)),
        "condition_evidence_rows": int(len(ev)),
        "condition_evidence_nominal_p_lt_0.05": int(ev["nominal_significant"].sum()) if not ev.empty else 0,
        "mapmecfs_files_gated": len(gated),
    }
    status = "partial" if sn_mods or not ev.empty else "blocked"
    pen = RAW / "pennsieve" / "dataset_356.json"
    group_cat = _load_catalog("catalog_group_pi_mecfs.json")
    ctx = {
        "log": log, "mods": mods, "ovl": ovl, "parts": parts, "omics": omics, "ev": ev, "checks": checks,
        "written": written, "status": status, "sample_size": sample_size,
        "retrieved": min((v.get("retrieved_at") for v in load_manifest(SOURCE_ID)["files"].values()
                          if v.get("retrieved_at")), default=utc_now_iso()),
        "github_commit": _github_commit(),
        "pennsieve_version": json.loads(pen.read_text()).get("version") if pen.exists() else None,
        "mapmecfs_modified": max((p.get("metadata_modified", "") for p in group_cat), default=None),
        "source_scan": scan_source_data_identifiers() if SOURCE_DATA_DIR.exists() else [],
        "r_object_scan": scan_github_r_objects() if GITHUB_DIR.exists() else [],
        "limitations": [
            "Single exploratory case-control cohort (paper analytic n = 17 PI-ME/CFS + 21 HV); deposits include "
            "participants outside that cohort and no open file flags adjudication status.",
            "Wearable/actigraphy/HRV/orthostatic/CPET per-participant data are gated on mapMECFS; open source data "
            "for these carry group labels only, so they are not person-level data here.",
            "Participant 311 is labelled Male in both GEO RNA-seq deposits and Female in both SomaLogic deposits; "
            "flagged (sex_concordant_across_sources = False), not corrected." if not parts.empty and
            (~parts["sex_concordant_across_sources"]).any() else "Sex labels agree across deposits.",
            "Serum SomaLogic is attached to participants through the Fatigue-NNN pairing published only in the "
            "CSF deposit (GSE251790); corroborated by sex/age/group agreement, not by an ID in the serum deposit.",
            "CSF metabolomics source data carry study numbers only on rows with Birth Sex = 1; other rows are "
            "not linkable.",
            "Supplementary Data 14A 'adjusted p-value' is not a monotone adjustment of its 'p-value' column "
            "(see audit), so fdr_significant is UNKNOWN for its 142 rows (kept as published in adj_p_value); "
            "Supplementary Data 18 log2FoldChange orientation is unverifiable (direction UNKNOWN).",
            "Supplementary Data 17 value scales, established against the GEO ADAT RFUs: " + (", ".join(
                f"{t} {v}" for t, v in ev.loc[ev["table_id"].str.startswith("17")].groupby("table_id")["value_scale"]
                .first().items()) if not ev.empty and "value_scale" in ev else UNKNOWN) + ". effect_value is the "
            "log2 ratio of group medians for RFU values and the difference of group medians for log2(RFU) values "
            "(value_scale column).",
            "Condition evidence rows for SD15/17/20/21 take direction from group medians computed from the "
            "individual values published in the same table; p-values are as published.",
            "GEO PBMC deposit has 27 samples while the PBMC DE analysis used 29 (source data 15 HV + 14 PI); "
            "stated n in MOESM3 differs from source-data counts for SD16C (stated 8 HV, Fig. 8F has 9).",
        ],
    }
    ctx["written"] = written
    LINKAGE_AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    LINKAGE_AUDIT_PATH.write_text(render_linkage_audit(ctx))
    (RAW / "DATA_AUDIT.md").write_text(render_data_audit(ctx))

    access_urls = sorted({_strip_query(v.get("url", "")) if "X-Amz" in v.get("url", "") else v.get("url", "")
                          for v in load_manifest(SOURCE_ID)["files"].values()})
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "NIH intramural post-infectious ME/CFS deep phenotyping (Walitt et al. 2024; mapMECFS, GEO, SRA, "
                "Pennsieve deposits)",
        "publisher": "NIH NINDS Division of Intramural Research; mapMECFS (RTI International); NCBI GEO/SRA; Pennsieve",
        "landing_url": MAPMECFS_GROUP_URL,
        "access_urls": access_urls,
        "license": "CC BY 4.0 (article, supplements, mapMECFS, Pennsieve); GEO/SRA public; GitHub GPL-3.0",
        "access_conditions": "GEO/SRA/Nature/GitHub/Pennsieve: open, no registration. mapMECFS data files: "
                             "registered account required (anonymous requests return HTTP 403 'Login Required'); "
                             "no account created, so mapMECFS files are not automatable.",
        "retrieved_at": ctx["retrieved"],
        "source_version": f"doi:{PAPER['doi']}; {_geo_version()}; GitHub commit {ctx['github_commit'] or UNKNOWN}",
        "update_date": "static publication deposits (GEO last updated 2024-03-06 per series matrices)",
        "data_layer": "person",
        "unit_of_observation": "participant (NIH study number) x omics modality x feature; plus published "
                               "condition-level differential analytes",
        "sample_size": sample_size,
        "geographic_resolution": "none",
        "person_level": True,
        "geographic": False,
        "omics": True,
        "wearable": False,
        "participant_linkage": "omics x omics only, on NIH study numbers (1xx HV / 3xx PI-ME/CFS) published in GEO "
                               "titles/count files and Fig. 6J source data; serum SomaLogic via the Fatigue-NNN "
                               "crosswalk in GSE251790. Wearable/HRV/orthostatic/CPET data carry no participant IDs "
                               "in open files and are NOT linked.",
        "true_participant_linkage_across_modalities": bool(len(sn_mods) >= 2),
        "wearable_omics_linkage": False,
        "linked_modalities": [m.modality_id for m in sn_mods],
        "status": status,
        "processed_outputs": sorted(written),
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md; docs/MAP_MECFS_LINKAGE_AUDIT.md",
        "ingestion_module": PRODUCER,
        "limitations": ctx["limitations"],
    })
    return ctx


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-fetch", action="store_true", help="rebuild from files already in data/raw")
    ap.add_argument("--pennsieve-bulk", action="store_true", help="also mirror Pennsieve EMG/fMRI files (~3.2 GB)")
    a = ap.parse_args(argv)
    ctx = run(fetch=not a.no_fetch, pennsieve_bulk=a.pennsieve_bulk)
    print(json.dumps({"status": ctx["status"], "written": ctx["written"], "sample_size": ctx["sample_size"]},
                     indent=2, default=str))


if __name__ == "__main__":
    main()
