"""NHANES 2011-2012 (G) and 2013-2014 (H): the full public laboratory catalogue as person-level layers.

The core laboratory files (CBC, BIOPRO, GHB, GLU, TCHOL, HDL, TRIGLY, VID, VITB12, THYROD) are ingested by
``measure_it.ingestion.nhanes`` into ``participant_labs__nhanes``. This module inventories EVERY file on the
NHANES laboratory data listing pages for 2011 and 2013 (saved under ``data/raw/nhanes_2011_2014/docs/listings``),
ingests every remaining individually linkable biospecimen file, and assigns each analyte to a ``lab_layer``.

Sources (verified 2026-09-24):
  * listing pages   https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory&CycleBeginYear={2011|2013}
  * component files https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{2011|2013}/DataFiles/<FILE>.xpt (+ <FILE>.htm)
  * oral microbiome https://wwwn.cdc.gov/nchs/data/nhanes/omp/<file>.txt (NHANES Oral Microbiome Project, 2009-2012;
                    linked on SEQN; only the 2011-2012 participants fall inside this project's cycles)

Joins are made ONLY on SEQN; participant_id = "nhanes:<SEQN>" (the namespace of the core NHANES tables).

Detection limits: NHANES releases below-limit results already filled with LLOD/sqrt(2) and flags them in a
"comment code" variable (usually <prefix>D<analyte>LC; 0 = at or above the detection limit, 1 = below the
lower limit of detection, 2 = at or above the upper limit / other codes per codebook). This module keeps the
released value, attaches the matched comment-code variable (``lc_variable``) and its value (``lc_code``), and
sets ``below_lod`` = (lc_code == 1). Nothing is re-imputed.

Subsamples: each file's sample weight variable(s) (WTMEC2YR = full MEC sample; WTSA2YR / WTSB2YR / WTSC2YR =
one-third environmental subsamples A/B/C; WTSAF2YR = morning fasting subsample; WTSOG2YR = OGTT subsample;
WTFSM = special smoking sample; WTSSxxxx = surplus-serum subsamples ...) are detected from the file and
carried per row (``weight_variable``, ``weight_value``), so downstream analyses can restrict to the sampled
people instead of imputing an unsampled layer.

Outputs (data/processed):
  participant_labs_extended__nhanes      long: one row per participant x analyte (released value, LC flag,
                                         subsample weight, lab_layer)
  participant_layer_membership__nhanes   one row per NHANES 2011-2014 participant: 0/1 per data layer
                                         (lab layers from both lab tables, questionnaire, medications, exam,
                                         mortality, wearable)
Tables: results/tables/nhanes_lab_inventory.csv (every listed file with the ingest decision) and
        results/tables/nhanes_lab_layers.csv (layer x file x weight, measured counts).

Run: uv run python -m measure_it.ingestion.nhanes_labs_extended
"""
from __future__ import annotations

import argparse
import html as _html
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import TABLES, UNKNOWN, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..provenance import add_provenance
from ..store import read_table, write_table
from . import nhanes as N

SOURCE_ID = N.SOURCE_ID                      # same raw directory / registry entry as the core NHANES ingest
SOURCE_NAME = N.SOURCE_NAME
PRODUCER = "measure_it.ingestion.nhanes_labs_extended"
CYCLES = {"G": 2011, "H": 2013}
CYCLE_LABEL = {"G": "2011-2012", "H": "2013-2014"}
LISTING_URL = "https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Laboratory&CycleBeginYear={year}"
SITE = "https://wwwn.cdc.gov"
OMP_BASE = "https://wwwn.cdc.gov/nchs/data/nhanes/omp/{file}"
OMP_FILES = ["dada2rb-alpha.txt", "dada2rb-variablelist-alpha.txt", "dada2rb-genus-relative.txt",
             "dada2rb-taxonomy-annotate.txt", "OralMicrobiomeDataDocumentation-508.pdf"]
OMP_URL_DOC = "https://wwwn.cdc.gov/nchs/data/nhanes/omp/OralMicrobiomeDataDocumentation-508.pdf"

# Files already ingested by measure_it.ingestion.nhanes into participant_labs__nhanes (component base names).
CORE_COMPONENTS = {"CBC", "BIOPRO", "GHB", "GLU", "TCHOL", "HDL", "TRIGLY", "VID", "VITB12", "THYROD"}

# Listed files that are NOT ingested, with the reason (component base name -> reason).
EXCLUDED = {
    "FASTQX": "fasting questionnaire (hours since last meal, fasting status): a specimen-collection covariate, "
              "not a biomarker panel",
    "UCPREG": "urine pregnancy test used to screen eligibility for examination components; not a biomarker layer",
    "POOLTF": "pooled-sample technical support file: pools of sera, not linkable to one participant",
    "BFRPOL": "pooled samples (brominated flame retardants measured in pools of sera): not linkable to one SEQN",
    "PCBPOL": "pooled samples (PCBs measured in pools of sera): not linkable to one SEQN",
    "PSTPOL": "pooled samples (organochlorine pesticides measured in pools of sera): not linkable to one SEQN",
    "DOXPOL": "pooled samples (dioxins/furans measured in pools of sera): not linkable to one SEQN",
    "SSEVD": "withdrawn by NCHS (listing shows 'Withdrawn'; no data link)",
    "FLDEW": "fluoride in household tap water: an environmental water sample, not a participant biospecimen",
}

# Component base name -> lab_layer. Every listed individual-level file not in CORE/EXCLUDED must appear here
# (enforced in build_inventory). Special-sample files ("...S" variants) map to the same layer as the main file
# and are flagged special_sample=True.
LAYER_OF = {
    # metabolic / lipids / glycaemia (fasting and OGTT subsamples)
    "INS": "metabolic_fasting_subsample", "APOB": "metabolic_fasting_subsample",
    "OGTT": "metabolic_ogtt_subsample",
    # nutritional biomarkers
    "FOLATE": "nutrition_folate_b12_status", "FOLFMS": "nutrition_folate_b12_status",
    "MMA": "nutrition_folate_b12_status", "SSHOLO": "nutrition_holotranscobalamin_surplus",
    "CUSEZN": "nutrition_trace_elements_serum", "UIO": "nutrition_iodine_urine",
    "FAS": "nutrition_fatty_acids_serum", "CAFE": "diet_caffeine_metabolites_urine",
    # hormones
    "TST": "hormones_sex_steroids",
    # immune / autoimmune / aging-neuro (surplus sera)
    "SSANA2": "immune_autoantibodies_surplus", "SSDFS": "immune_autoantibodies_surplus",
    "TGEMA": "immune_celiac_serology",
    "SSKL": "aging_klotho_surplus", "SSSNFL": "neuro_neurofilament_light_surplus",
    # infectious serology / molecular
    "HEPA": "infectious_hepatitis", "HEPB_S": "infectious_hepatitis", "HEPBD": "infectious_hepatitis",
    "HEPC": "infectious_hepatitis", "HEPE": "infectious_hepatitis", "SSHCV": "infectious_hepatitis",
    "SSHEPC": "infectious_hepatitis",
    "HIV": "infectious_hiv", "HSV": "infectious_herpes_cmv", "CMV": "infectious_herpes_cmv",
    "CHLMDA": "infectious_sti_urine", "TRICH": "infectious_sti_urine", "SSCT": "infectious_sti_urine",
    "ORHPV": "infectious_hpv", "HPVSWR": "infectious_hpv", "HPVP": "infectious_hpv",
    "SSTOXO": "infectious_parasite_serology_surplus", "SSTOCA": "infectious_parasite_serology_surplus",
    "TB": "infectious_tuberculosis_igra",
    # urinary markers
    "ALB_CR": "urine_albumin_creatinine", "UCFLOW": "urine_flow_rate", "UCOSMO": "urine_osmolality",
    # environmental chemicals
    "PBCD": "env_metals_blood", "IHGEM": "env_mercury_species_blood",
    "UHM": "env_metals_urine", "UM": "env_metals_urine", "UHG": "env_metals_urine",
    "UAS": "env_arsenic_urine", "UTAS": "env_arsenic_urine",
    "PFC": "env_pfas_serum", "PFAS": "env_pfas_serum", "SSPFSU": "env_pfas_serum_surplus",
    "SSPFAS": "env_pfas_serum_surplus", "SSPFAC": "env_pfas_serum_surplus",
    "UVOC": "env_voc_metabolites_urine", "VOCWB": "env_vocs_blood", "VOCMWB": "env_vocs_blood",
    "PHTHTE": "env_phthalates_urine", "SSPHTE": "env_phthalates_urine_surplus",
    "EPH": "env_phenols_parabens_urine", "EPHPP": "env_phenols_parabens_urine",
    "PAH": "env_pahs_urine",
    "PP": "env_chlorophenols_herbicides_urine", "UPHOPM": "env_pesticides_urine", "OPD": "env_pesticides_urine",
    "DEET": "env_pesticides_urine", "SSGLYP": "env_glyphosate_urine_surplus",
    "PERNT": "env_perchlorate_nitrate_thiocyanate_urine",
    "SSFR": "env_flame_retardant_metabolites_urine_surplus", "SSFLRT": "env_flame_retardant_metabolites_urine_surplus",
    "COTNAL": "tobacco_cotinine", "COT": "tobacco_cotinine", "UCOT": "tobacco_nicotine_metabolites_urine",
    "TSNA": "tobacco_nitrosamines_urine", "VNA": "env_volatile_nitrosamines_urine",
    "AA": "env_aromatic_amines_urine", "HCAA": "env_heterocyclic_amines_urine",
    "AMDGYD": "env_acrylamide_hb_adducts", "AMDGDS": "env_acrylamide_hb_adducts",
    "ALD": "env_aldehydes_serum", "ETHOX": "env_ethylene_oxide_hb_adducts", "FORMAL": "env_formaldehyde_hb_adducts",
    "FLDEP": "env_fluoride_plasma", "SSTERP": "env_terpenes_serum_surplus",
}
# Special-sample twins (component -> main component), resolved before LAYER_OF lookup.
SPECIAL_SAMPLE_OF = {"UASS": "UAS", "UHMS": "UHM", "UMS": "UM", "UTASS": "UTAS", "PERNTS": "PERNT", "PAHS": "PAH",
                     "UVOCS": "UVOC", "VOCWBS": "VOCWB", "UCOTS": "UCOT", "AAS": "AA", "HCAAS": "HCAA",
                     "VNAS": "VNA", "ALDS": "ALD", "ETHOXS": "ETHOX", "FORMAS": "FORMAL", "AMDGDS": "AMDGYD"}

# Broad family for each layer (used in tables and the overlap figure).
LAYER_FAMILY_RULES = [
    ("core_", "core clinical chemistry / CBC"), ("lipids_glycaemia", "metabolic / lipids / glycaemia"),
    ("metabolic_", "metabolic / lipids / glycaemia"), ("nutrition_", "nutritional biomarkers"),
    ("diet_", "nutritional biomarkers"), ("hormones_", "hormones"), ("immune_", "immune / autoimmune"),
    ("aging_", "aging / neuro (surplus sera)"), ("neuro_", "aging / neuro (surplus sera)"),
    ("infectious_", "infectious serology / molecular"), ("urine_", "urinary markers"),
    ("tobacco_", "tobacco exposure"), ("env_", "environmental chemicals"), ("oral_microbiome", "oral microbiome (16S)"),
    ("thyroid", "hormones"),
]

# Layers of the core table (participant_labs__nhanes), by component.
CORE_LAYER_OF = {"CBC": "core_cbc", "BIOPRO": "core_biochemistry_profile",
                 "GHB": "lipids_glycaemia_full", "TCHOL": "lipids_glycaemia_full", "HDL": "lipids_glycaemia_full",
                 "VID": "nutrition_vitamin_d_b12", "VITB12": "nutrition_vitamin_d_b12",
                 "GLU": "metabolic_fasting_subsample", "TRIGLY": "metabolic_fasting_subsample",
                 "THYROD": "hormones_thyroid_subsample"}

WEIGHT_LABELS = {
    "WTMEC2YR": "full MEC examined sample", "WTSA2YR": "one-third subsample A", "WTSB2YR": "one-third subsample B",
    "WTSC2YR": "one-third subsample C", "WTSAF2YR": "morning fasting subsample", "WTSOG2YR": "OGTT subsample",
    "WTFSM": "special smoking sample", "WTSH2YR": "hepatitis / half-sample weight",
}


def layer_family(layer: str) -> str:
    for prefix, fam in LAYER_FAMILY_RULES:
        if layer.startswith(prefix):
            return fam
    return "other"


def component_of(file_stem: str) -> tuple[str, str]:
    """'UHMS_G' -> ('UHMS', 'G'); 'HEPB_S_G' -> ('HEPB_S', 'G'); 'SSHCV_E' -> ('SSHCV', 'E')."""
    base, _, suf = file_stem.rpartition("_")
    return base, suf


# --------------------------------------------------------------------------------------------
# Listing pages and inventory
# --------------------------------------------------------------------------------------------

def _clean(s: str) -> str:
    return re.sub(r"\s+", " ", _html.unescape(re.sub(r"<[^>]+>", "", s))).strip()


def fetch_listing(year: int) -> Path:
    (raw_dir(SOURCE_ID) / "docs" / "listings").mkdir(parents=True, exist_ok=True)
    return download_file(LISTING_URL.format(year=year), SOURCE_ID, f"docs/listings/laboratory_{year}.html",
                         allow_html=True, min_bytes=20000)


def parse_listing(path: Path, year: int) -> pd.DataFrame:
    """One row per listed data file: description, data/doc hrefs, publication date text."""
    t = path.read_text(encoding="utf-8", errors="replace")
    out = []
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", t, re.S):
        tds = re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)
        if len(tds) < 3:
            continue
        desc = _clean(tds[0])
        hrefs = re.findall(r"href=['\"]([^'\"]+)['\"]", row)
        data = next((h for h in hrefs if h.lower().endswith(".xpt")), "")
        doc = next((h for h in hrefs if h.lower().endswith(".htm")), "")
        other = next((h for h in hrefs if not h.lower().endswith((".xpt", ".htm"))), "")
        label = _clean(tds[2]) if len(tds) > 2 else ""
        stem = Path(data).stem if data else (label.split(" ")[0] if label else "")
        if "microbiome" in desc.lower():
            stem = "OMP_all_years"
        out.append({"cycle_begin_year": year, "file": stem, "description": desc,
                    "published": _clean(tds[-1]),
                    "data_url": (SITE + data) if data.startswith("/") else data,
                    "doc_url": (SITE + doc) if doc.startswith("/") else doc,
                    "other_url": (SITE + other) if other.startswith("/") else other})
    return pd.DataFrame(out)


def decide(file: str, description: str, published: str, data_url: str) -> tuple[str, str, str]:
    """(decision, lab_layer, reason) for one listed file."""
    if "microbiome" in description.lower():
        return "ingested", "oral_microbiome_16s", "NHANES Oral Microbiome Project text files (all years); SEQN-linked"
    comp, _ = component_of(file)
    if "withdrawn" in published.lower() or not data_url:
        return "excluded", "", ("withdrawn by NCHS (listing shows 'Withdrawn'; no data link)"
                                if "withdrawn" in published.lower() else EXCLUDED.get(comp, "no data link on the listing page"))
    if comp in CORE_COMPONENTS:
        return "core", CORE_LAYER_OF[comp], "already ingested by measure_it.ingestion.nhanes (participant_labs__nhanes)"
    if comp in EXCLUDED:
        return "excluded", "", EXCLUDED[comp]
    main = SPECIAL_SAMPLE_OF.get(comp, comp)
    if main in LAYER_OF:
        return "ingested", LAYER_OF[main], ("special smoking sample twin of " + main) if main != comp else ""
    raise KeyError(f"listed file {file} ({description}) has no layer assignment or exclusion reason")


def build_inventory(log=print) -> pd.DataFrame:
    frames = []
    for cyc, year in CYCLES.items():
        p = fetch_listing(year)
        df = parse_listing(p, year)
        df["cycle_code"] = cyc
        frames.append(df)
    inv = pd.concat(frames, ignore_index=True)
    dec = inv.apply(lambda r: decide(r["file"], r["description"], r["published"], r["data_url"]), axis=1)
    inv["decision"] = [d[0] for d in dec]
    inv["lab_layer"] = [d[1] for d in dec]
    inv["reason"] = [d[2] for d in dec]
    inv["layer_family"] = inv["lab_layer"].map(lambda x: layer_family(x) if x else "")
    log(f"[inventory] {len(inv)} listed files: " + ", ".join(f"{k}={v}" for k, v in inv["decision"].value_counts().items()))
    return inv


def fetch_files(inv: pd.DataFrame, log=print) -> None:
    (raw_dir(SOURCE_ID) / "docs").mkdir(exist_ok=True)
    todo = inv[(inv["decision"] == "ingested") & (inv["lab_layer"] != "oral_microbiome_16s")]
    for r in todo.itertuples():
        download_file(r.data_url, SOURCE_ID, f"{r.file}.xpt", min_bytes=500, max_retries=4)
        if r.doc_url:
            download_file(r.doc_url, SOURCE_ID, f"docs/{r.file}.htm", allow_html=True, min_bytes=3000, max_retries=4)
    (raw_dir(SOURCE_ID) / "omp").mkdir(exist_ok=True)
    for f in OMP_FILES:
        download_file(OMP_BASE.format(file=f), SOURCE_ID, f"omp/{f}", min_bytes=1000, max_retries=4)
    download_file("https://wwwn.cdc.gov/Nchs/Nhanes/Omp/Default.aspx", SOURCE_ID, "docs/listings/omp_default.html",
                  allow_html=True, min_bytes=5000)
    log(f"[fetch] {len(todo)} component files + {len(OMP_FILES)} oral-microbiome files present")


# --------------------------------------------------------------------------------------------
# Variable roles, comment-code (LC) matching and value types
# --------------------------------------------------------------------------------------------

LC_LABEL_RE = re.compile(r"comment|detection limit", re.I)
PROCESS_LABEL_RE = re.compile(
    r"\b(time|minutes?|hours?|volume of urine|food fast|administer|amount of glucose|drank|incomplete|retest|"
    r"reflex|retesting|no surplus|specimen|date|season|dilution|batch)\b|b/w", re.I)
LAB_PREFIX_RE = re.compile(r"^(LB|UR|SS|SD|OR)")
POSITIVE_RE = re.compile(r"^(positive|weakly positive|reactive|detected|present|positive, .*)$", re.I)
NEGATIVE_RE = re.compile(r"^(negative|non-?reactive|not detected|absent)$", re.I)
URINE_CREATININE_VARS = {"URXUCR"}


def _strip_unit(label: str) -> str:
    return re.sub(r"\s*\([^()]*\)\s*$", "", label or "").strip().lower()


def is_lc_variable(name: str, cb: dict) -> bool:
    """A detection-limit comment-code variable: codebook label says 'comment' or its codes mention a detection limit."""
    info = cb.get(name, {})
    codes = " ".join(info.get("codes", {}).values())
    if "detection limit" in codes.lower():
        return True
    return bool(re.search(r"comment", info.get("label", ""), re.I)) and name.upper().endswith(("LC", "L"))


def _norm_label(label: str) -> str:
    x = re.sub(r"\([^()]*\)", " ", label or "")
    x = re.sub(r"\b(comment\s*cod\w*|comment|cmt\.?(\s*code)?|code|above detect\w*\s*limit|below detect\w*\s*limit)\b",
               " ", x, flags=re.I)
    x = re.sub(r"\bbr\.?\s*iso( of)?\b", "branched iso of", x, flags=re.I)
    x = re.sub(r"\blin\.?\b", "linear", x, flags=re.I)
    return re.sub(r"[^a-z0-9]+", " ", x.lower()).strip()


def match_lc(lc_name: str, analytes: list[str], labels: dict[str, str] | None = None) -> str | None:
    """Map a detection-limit comment-code variable to its analyte.

    1. name rule: analyte <P>X<core> (LBX / URX / SSX) or <P>D<core>, comment code <P>D<core>LC (or ...L when the
       8-character SAS name limit bites);
    2. label rule: the comment code's label minus 'comment code' equals (or is a truncation of) the analyte label
       minus units (e.g. URDUA3LC 'Urinary Arsenous acid comment code' -> URXUAS3 'Urinary Arsenous acid (ug/L)');
    3. unique name-prefix fallback. Unmatched codes are reported (role comment_code_unmatched), never guessed.
    """
    labels = labels or {}
    if lc_name.startswith("SS"):  # surplus files: SS<core>L / SS<core>LC <-> SS<core> (no X/D letter)
        c2 = re.sub(r"(LC|L)$", "", lc_name[2:])
        hit = [a for a in analytes if a.startswith("SS") and a[2:] == c2]
        if len(hit) == 1:
            return hit[0]
    core = re.sub(r"(LC|L)$", "", lc_name[3:])
    pre = lc_name[:2]
    cands = [a for a in analytes if a[:2] == pre]
    exact = [a for a in cands if a[3:] == core] or [a for a in cands if a[3:] == core + "SI"]
    if exact:
        return exact[0]
    lab = _norm_label(labels.get(lc_name, ""))
    if lab:
        hits = []
        for a in analytes:
            al = _norm_label(labels.get(a, ""))
            if not al:
                continue
            abbrev = [_norm_label(x) for x in re.findall(r"\(([^()]*)\)", labels.get(a, ""))]
            if al == lab or (len(lab) >= 4 and (al.startswith(lab) or lab.startswith(al))) or lab in abbrev:
                hits.append(a)
        hits = [h for h in hits if not h.upper().endswith("SI")] or hits
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            eq = [h for h in hits if _norm_label(labels.get(h, "")) == lab]
            if len(eq) == 1:
                return eq[0]
    pref = [a for a in cands if a[3:].startswith(core) and len(core) >= 2]
    if len(pref) == 1:
        return pref[0]
    return None


def value_type(name: str, cb: dict, values: pd.Series) -> str:
    codes = cb.get(name, {}).get("codes", {})
    labels = [v.lower() for v in codes.values()]
    if any("range of values" in lab for lab in labels):
        return "continuous"
    if codes and len([k for k in codes if k not in (".", "< blank >")]) >= 1:
        return "categorical"
    return "continuous" if values.dropna().nunique() > 12 else "categorical"


def positive_code_map(name: str, cb: dict) -> dict[float, float]:
    """Categorical result code -> 1.0 (positive-like) / 0.0 (negative-like); other codes -> not mapped (NaN)."""
    out = {}
    for code, lab in cb.get(name, {}).get("codes", {}).items():
        try:
            c = float(code)
        except ValueError:
            continue
        if POSITIVE_RE.match(lab.strip()):
            out[c] = 1.0
        elif NEGATIVE_RE.match(lab.strip()):
            out[c] = 0.0
    return out if (1.0 in out.values() and 0.0 in out.values()) else {}


def _subsample_text(weight_var: str) -> str:
    if weight_var in WEIGHT_LABELS:
        return f"{WEIGHT_LABELS[weight_var]} (weight {weight_var})"
    if weight_var.startswith("WTSS"):
        return f"surplus-serum/urine subsample (weight {weight_var})"
    return f"component subsample (weight {weight_var})"


# --------------------------------------------------------------------------------------------
# Build the long extended-lab table
# --------------------------------------------------------------------------------------------

def _file_meta(file: str) -> dict:
    m = load_manifest(SOURCE_ID)["files"].get(f"{file}.xpt", {})
    return {"last_modified": m.get("last_modified", UNKNOWN), "retrieved_at": m.get("retrieved_at", UNKNOWN),
            "url": m.get("url", UNKNOWN)}


def extract_file(file: str, lab_layer: str, description: str, participants: pd.DataFrame) -> tuple[pd.DataFrame, list]:
    """Long rows for one component file + a variable dictionary (list of dicts)."""
    comp, suf = component_of(file)
    cb = N.parse_codebook(comp, suf)
    path = raw_dir(SOURCE_ID) / f"{file}.xpt"
    lay = N.xport_layout(path)
    df = N.read_xport(path)
    df["SEQN"] = df["SEQN"].astype("int64")
    df = df[df["SEQN"].isin(participants["seqn"])].reset_index(drop=True)
    names = [v["name"] for v in lay["variables"] if v["name"] != "SEQN"]
    numeric = {v["name"]: v["numeric"] for v in lay["variables"]}
    weights = [n for n in names if n.upper().startswith("WT")]
    lcs = [n for n in names if n not in weights and numeric[n] and is_lc_variable(n, cb)]
    # an 'above the upper detection limit' flag (e.g. SSNFLH) is not a lower-limit comment code: kept aside
    upper = [n for n in lcs if re.search(r"above (upper )?detect", cb.get(n, {}).get("label", ""), re.I)]
    lcs = [n for n in lcs if n not in upper]
    cand = [n for n in names if n not in weights and n not in lcs and n not in upper]
    vdict, analytes, seen_labels = [], [], {}
    for n in cand:
        info = cb.get(n, {})
        label = info.get("label", n)
        role = "analyte"
        if not numeric[n]:
            role = "character_variable_not_ingested"
        elif not LAB_PREFIX_RE.match(n):
            role = "questionnaire_item_not_ingested"   # e.g. VTQ* exposure questions bundled in the VOC files
        elif (not n.startswith(("LBX", "LBD")) and PROCESS_LABEL_RE.search(label)
              and n not in URINE_CREATININE_VARS):
            role = "process_variable"
        else:
            key = _strip_unit(label)
            if key in seen_labels and n.upper().endswith("SI"):
                role = "si_unit_duplicate_not_ingested"
            elif key in seen_labels and re.search(r"\((umol|nmol|mmol|pmol)/", label):
                role = "si_unit_duplicate_not_ingested"
            else:
                seen_labels.setdefault(key, n)
        if n in URINE_CREATININE_VARS:
            role = "urine_creatinine"
        vdict.append({"component_file": file, "lab_variable": n, "label": label, "unit": info.get("unit", ""),
                      "target": info.get("target", ""), "role": role})
        if role in ("analyte", "process_variable", "urine_creatinine"):
            analytes.append(n)
    lc_map = {}
    roles_now = {d["lab_variable"]: d["role"] for d in vdict}
    for lc in lcs:
        a = match_lc(lc, [x for x in analytes if roles_now.get(x) == "analyte"],
                     {v: cb.get(v, {}).get("label", "") for v in analytes + [lc]})
        vdict.append({"component_file": file, "lab_variable": lc, "label": cb.get(lc, {}).get("label", lc),
                      "unit": "", "target": cb.get(lc, {}).get("target", ""),
                      "role": "comment_code" + ("" if a else "_unmatched"), "matched_analyte": a or ""})
        if a:
            lc_map[a] = lc
    for u in upper:
        vdict.append({"component_file": file, "lab_variable": u, "label": cb.get(u, {}).get("label", u), "unit": "",
                      "target": cb.get(u, {}).get("target", ""), "role": "upper_limit_flag_not_ingested"})
    for w in weights:
        vdict.append({"component_file": file, "lab_variable": w, "label": cb.get(w, {}).get("label", w), "unit": "",
                      "target": cb.get(w, {}).get("target", ""), "role": "sample_weight"})
    weight_var = weights[0] if weights else ""
    if not analytes:
        return pd.DataFrame(), vdict
    long = df[["SEQN"] + analytes].melt(id_vars="SEQN", var_name="lab_variable", value_name="value")
    long = long.dropna(subset=["value"])
    vt = {a: value_type(a, cb, df[a]) for a in analytes}
    roles = {d["lab_variable"]: d["role"] for d in vdict}
    long["analyte_role"] = long["lab_variable"].map(roles)
    long["value_type"] = long["lab_variable"].map(vt)
    long["lab_name"] = long["lab_variable"].map(lambda v: re.sub(r"\s*\([^()]*\)\s*$", "", cb.get(v, {}).get("label", v)))
    long["unit"] = long["lab_variable"].map(lambda v: cb.get(v, {}).get("unit", ""))
    code_labels = {a: cb.get(a, {}).get("codes", {}) for a in analytes if vt[a] == "categorical"}

    def _code_label(row_var, val):
        codes = code_labels.get(row_var)
        if not codes:
            return ""
        key = str(int(val)) if float(val).is_integer() else str(val)
        return codes.get(key, "")

    long["value_label"] = [(_code_label(v, x) if vt[v] == "categorical" else "")
                           for v, x in zip(long["lab_variable"], long["value"])]
    pos = {a: positive_code_map(a, cb) for a in analytes if vt[a] == "categorical"}
    long["result_positive"] = [pos.get(v, {}).get(float(x), np.nan) if vt[v] == "categorical" else np.nan
                               for v, x in zip(long["lab_variable"], long["value"])]
    lc_vals = {}
    for a, lc in lc_map.items():
        lc_vals[a] = df.set_index("SEQN")[lc]
    long["lc_variable"] = long["lab_variable"].map(lambda v: lc_map.get(v, ""))
    lc_code = np.full(len(long), np.nan)
    for a, ser in lc_vals.items():
        m = (long["lab_variable"] == a).to_numpy()
        lc_code[m] = long.loc[m, "SEQN"].map(ser).to_numpy(dtype=float)
    long["lc_code"] = lc_code
    long["below_lod"] = pd.array(np.where(np.isnan(lc_code), None, lc_code == 1.0), dtype="boolean")
    if weight_var:
        wser = df.set_index("SEQN")[weight_var]
        long["weight_variable"] = weight_var
        long["weight_value"] = long["SEQN"].map(wser).astype(float)
        long["subsample"] = _subsample_text(weight_var)
    else:
        wmec = participants.set_index("seqn")["wtmec2yr"]
        long["weight_variable"] = "WTMEC2YR"
        long["weight_value"] = long["SEQN"].map(wmec).astype(float)
        long["subsample"] = ("full MEC examined sample within the component's eligible age/sex group "
                             "(no component weight; weight WTMEC2YR)")
    comp_main = SPECIAL_SAMPLE_OF.get(comp, comp)
    long["special_sample"] = comp_main != comp or weight_var == "WTFSM"
    long["component_file"] = file
    long["component_description"] = description
    long["lab_layer"] = lab_layer
    long["layer_family"] = layer_family(lab_layer)
    return long, vdict


def build_omp(participants: pd.DataFrame, log=print) -> pd.DataFrame:
    """Oral microbiome (16S rRNA V4, DADA2; rarefied) -> long rows for 2011-2012 participants.

    Kept: alpha diversity (observed ASVs, Faith PD, Shannon, inverse Simpson) at a rarefaction depth of 10,000
    reads, averaged over NCHS's 10 rarefactions; genus relative abundance (RB genus table) for genera detected in
    >= 5% of the 2011-2012 samples (zeros kept: they are measurements).
    """
    d = raw_dir(SOURCE_ID) / "omp"
    alpha = pd.read_csv(d / "dada2rb-alpha.txt", sep="\t")
    alpha["SEQN"] = alpha["SEQN"].astype("int64")
    alpha = alpha[alpha["SEQN"].isin(participants["seqn"])]
    rows = []
    metrics = sorted({re.sub(r"_\d+_\d+$", "", c) for c in alpha.columns if c.startswith("RB_")})
    for m in metrics:
        cols = [c for c in alpha.columns if re.fullmatch(rf"{m}_10000_\d+", c)]
        if not cols:
            continue
        v = alpha[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        rows.append(pd.DataFrame({"SEQN": alpha["SEQN"], "lab_variable": f"{m}_10000_mean", "value": v,
                                  "lab_name": f"oral microbiome alpha diversity: {m.replace('RB_', '')} "
                                              f"(rarefied to 10,000 reads, mean of 10 rarefactions)",
                                  "unit": "index"}))
    g = pd.read_csv(d / "dada2rb-genus-relative.txt", sep="\t")
    g = g[pd.to_numeric(g["SEQN"], errors="coerce").notna()].copy()
    g["SEQN"] = g["SEQN"].astype("int64")
    g = g[g["SEQN"].isin(participants["seqn"])]
    tax = pd.read_csv(d / "dada2rb-taxonomy-annotate.txt", sep="\t", dtype=str)
    tax = tax.apply(lambda c: c.str.strip())
    names = dict(zip(tax.iloc[:, 0], tax.iloc[:, 1]))  # Name -> "Taxonomy in SILVA v123" (Kingdom;Phylum;...;Genus)
    gcols = [c for c in g.columns if c.startswith("RB_genus")]
    vals = g[gcols].apply(pd.to_numeric, errors="coerce")
    prev = (vals > 0).mean()
    keep = [c for c in gcols if prev[c] >= 0.05]
    gl = pd.concat([g[["SEQN"]], vals[keep]], axis=1).melt(id_vars="SEQN", var_name="lab_variable", value_name="value")
    gl = gl.dropna(subset=["value"])
    gl["lab_name"] = gl["lab_variable"].map(lambda c: "oral microbiome genus relative abundance: "
                                            + str(names.get(c, c)).split(";")[-1] + f" [{names.get(c, c)}]")
    gl["unit"] = "relative abundance (0-1)"
    rows.append(gl)
    out = pd.concat(rows, ignore_index=True).dropna(subset=["value"])
    out["analyte_role"] = "analyte"
    out["value_type"] = "continuous"
    out["value_label"] = ""
    out["result_positive"] = np.nan
    out["lc_variable"] = ""
    out["lc_code"] = np.nan
    out["below_lod"] = pd.array([None] * len(out), dtype="boolean")
    out["weight_variable"] = "WTMEC2YR"
    out["weight_value"] = out["SEQN"].map(participants.set_index("seqn")["wtmec2yr"]).astype(float)
    out["subsample"] = ("NHANES Oral Microbiome Project: oral rinse, participants aged 14-69 with a sequenced "
                        "sample (no separate weight published; weight WTMEC2YR)")
    out["special_sample"] = False
    out["component_file"] = "OMP_dada2rb"
    out["component_description"] = "Oral Microbiome Project (DADA2, rarefied): alpha diversity + genus relative abundance"
    out["lab_layer"] = "oral_microbiome_16s"
    out["layer_family"] = layer_family("oral_microbiome_16s")
    log(f"[omp] {out['SEQN'].nunique()} participants, {len(keep)} genera kept (prevalence >= 5%) of {len(gcols)}")
    return out


def build_extended(inv: pd.DataFrame, log=print) -> tuple[pd.DataFrame, pd.DataFrame]:
    participants = read_table("participants__nhanes", columns=["participant_id", "seqn", "cycle_code", "wtmec2yr"])
    frames, vd = [], []
    todo = inv[(inv["decision"] == "ingested") & (inv["lab_layer"] != "oral_microbiome_16s")]
    for r in todo.itertuples():
        long, v = extract_file(r.file, r.lab_layer, r.description, participants)
        for x in v:
            x["lab_layer"] = r.lab_layer
        vd.extend(v)
        if len(long):
            frames.append(long)
    frames.append(build_omp(participants, log=log))
    ext = pd.concat(frames, ignore_index=True)
    ext = ext.merge(participants[["seqn", "participant_id", "cycle_code"]].rename(columns={"seqn": "SEQN"}),
                    on="SEQN", how="left", validate="many_to_one")
    ext["cycle"] = ext["cycle_code"].map(CYCLE_LABEL)
    ext["record_key"] = ext["component_file"] + ":SEQN=" + ext["SEQN"].astype(str) + ":" + ext["lab_variable"]
    ext = ext.rename(columns={"SEQN": "seqn"})
    order = ["participant_id", "seqn", "cycle", "cycle_code", "lab_layer", "layer_family", "component_file",
             "component_description", "lab_variable", "lab_name", "unit", "analyte_role", "value_type", "value",
             "value_label", "result_positive", "lc_variable", "lc_code", "below_lod", "subsample", "weight_variable",
             "weight_value", "special_sample", "record_key"]
    ext = ext[order]
    assert not ext.duplicated(["participant_id", "component_file", "lab_variable"]).any()
    log(f"[extended] {len(ext):,} rows, {ext['participant_id'].nunique():,} participants, "
        f"{ext['lab_layer'].nunique()} layers, {ext['component_file'].nunique()} files")
    return ext, pd.DataFrame(vd)


def write_extended(ext: pd.DataFrame) -> Path:
    parts = []
    for f, sub in ext.groupby("component_file", sort=True):
        if f.startswith("OMP"):
            man = load_manifest(SOURCE_ID)["files"].get("omp/dada2rb-genus-relative.txt", {})
            ver = (f"NHANES Oral Microbiome Project public files (dada2rb alpha/genus tables; last-modified "
                   f"{man.get('last_modified', UNKNOWN)})")
            retr = man.get("retrieved_at", UNKNOWN)
        else:
            meta = _file_meta(f)
            label = CYCLE_LABEL.get(component_of(f)[1], "2011-2014")
            ver = f"NHANES {label} public release; {f}.xpt last-modified {meta['last_modified']}"
            retr = meta["retrieved_at"]
        parts.append(add_provenance(
            sub, data_layer="person", source_name=SOURCE_NAME, source_version=ver, retrieved_at=retr,
            evidence_type="person_lab_measurement", source_record_id="record_key",
            source_geographic_resolution="none", evidence_level=sub["lab_layer"],
            provenance_notes=("Released values as published by NCHS; below-detection-limit results carry the NCHS "
                              "fill value (LLOD/sqrt(2)) and below_lod=True from the matched comment-code variable "
                              "(lc_variable). Subsample components carry their own weight (weight_variable). "
                              "evidence_level = lab_layer.")))
    out = pd.concat(parts, ignore_index=True).drop(columns=["record_key"])
    return write_table(out, "participant_labs_extended__nhanes", producer=PRODUCER,
                       description="NHANES 2011-2014 full laboratory catalogue beyond the core panels (long; lab_layer)")


# --------------------------------------------------------------------------------------------
# Layer membership (one row per participant) and layer summaries
# --------------------------------------------------------------------------------------------

NONLAB_LAYERS = {
    "questionnaire_phq9": ("questionnaire", "PHQ-9 complete (all nine DPQ010-DPQ090 items answered)"),
    "questionnaire_pfq": ("questionnaire", "physical-functioning composite answered (PFQ049/051/054/057/059; ages 20+)"),
    "questionnaire_self_rated_health": ("questionnaire", "HSD010 self-rated health answered"),
    "medications_rx_interview": ("medications", "prescription-medication question answered (RXDUSE yes/no)"),
    "exam_body_measures": ("exam", "BMI measured (BMX)"),
    "exam_blood_pressure": ("exam", "mean systolic BP available (BPX)"),
    "mortality_linked": ("mortality", "LMF-eligible with known vital status (follow-up to 2019-12-31)"),
    "wearable_valid": ("wearable", "passes the valid-wear rule (>= 4 valid days; wrist MIMS)"),
}


def build_membership(ext: pd.DataFrame, log=print) -> pd.DataFrame:
    parts = read_table("participants__nhanes", columns=["participant_id", "seqn", "cycle_code", "age_years",
                                                         "mec_examined", "pregnant_at_exam"])
    clin = read_table("participant_clinical_features__nhanes",
                      columns=["participant_id", "phq9_total", "any_functional_limitation_pfq",
                               "self_rated_health_mec_1excellent_5poor", "rx_any_past_month", "bmi", "sbp_mean_mmhg",
                               "mort_eligstat", "mort_status"])
    wear = read_table("participant_wearable_features__nhanes", columns=["participant_id", "passes_valid_wear_rule"])
    core = read_table("participant_labs__nhanes", columns=["participant_id", "component_file", "value"])
    m = parts.copy()
    ids = m["participant_id"]
    c = clin.set_index("participant_id").reindex(ids)
    m["questionnaire_phq9"] = c["phq9_total"].notna().to_numpy()
    m["questionnaire_pfq"] = c["any_functional_limitation_pfq"].notna().to_numpy()
    m["questionnaire_self_rated_health"] = c["self_rated_health_mec_1excellent_5poor"].notna().to_numpy()
    m["medications_rx_interview"] = c["rx_any_past_month"].notna().to_numpy()
    m["exam_body_measures"] = c["bmi"].notna().to_numpy()
    m["exam_blood_pressure"] = c["sbp_mean_mmhg"].notna().to_numpy()
    m["mortality_linked"] = ((c["mort_eligstat"] == "eligible") & c["mort_status"].notna()).to_numpy()
    wv = set(wear.loc[wear["passes_valid_wear_rule"], "participant_id"])
    m["wearable_valid"] = ids.isin(wv).to_numpy()
    core = core.dropna(subset=["value"])
    core["lab_layer"] = core["component_file"].str.rsplit("_", n=1).str[0].map(CORE_LAYER_OF)
    lab_sets = {k: set(g["participant_id"]) for k, g in core.groupby("lab_layer")}
    a = ext[ext["analyte_role"] == "analyte"]
    for k, g in a.groupby("lab_layer"):
        lab_sets[k] = lab_sets.get(k, set()) | set(g["participant_id"])
    for k in sorted(lab_sets):
        m[f"lab__{k}"] = ids.isin(lab_sets[k]).to_numpy()
    lab_cols = [c_ for c_ in m.columns if c_.startswith("lab__")]
    m["n_lab_layers"] = m[lab_cols].sum(axis=1)
    m["wearable_valid_adult"] = m["wearable_valid"] & (m["age_years"] >= 18)
    m["analysis_population"] = m["wearable_valid_adult"] & ~m["pregnant_at_exam"].fillna(False).astype(bool)
    for col in [c_ for c_ in m.columns if c_ in NONLAB_LAYERS or c_.startswith("lab__")]:
        m[col] = m[col].astype(int)
    log(f"[membership] {len(m)} participants; {len(lab_cols)} lab layers; "
        f"analysis population {int(m['analysis_population'].sum())}")
    return m


def write_membership(m: pd.DataFrame) -> Path:
    df = m.copy()
    df["record_key"] = "SEQN=" + df["seqn"].astype(str)
    parts = []
    for cyc, sub in df.groupby("cycle_code", sort=True):
        meta = N.cycle_meta(cyc)
        parts.append(add_provenance(
            sub, data_layer="person", source_name=SOURCE_NAME, source_version=meta["source_version"],
            retrieved_at=meta["retrieved_at"], evidence_type="person_derived_feature", source_record_id="record_key",
            source_geographic_resolution="none",
            provenance_notes="1 = the participant has >= 1 measured value in that layer (lab__* from "
                             "participant_labs__nhanes and participant_labs_extended__nhanes incl. special-sample "
                             "files; questionnaire/exam/medication/mortality/wearable definitions in "
                             "results/tables/nhanes_lab_layers.csv). Joined on SEQN only."))
    out = pd.concat(parts, ignore_index=True).drop(columns=["record_key"])
    return write_table(out, "participant_layer_membership__nhanes", producer=PRODUCER,
                       description="NHANES 2011-2014 per-participant data-layer membership (0/1 per layer)")


def layer_table(ext: pd.DataFrame, m: pd.DataFrame, inv: pd.DataFrame) -> pd.DataFrame:
    """One row per layer: files, weights, analytes, measured counts, LOD share, coverage in the wearable adults."""
    pop = m[m["analysis_population"]]
    rows = []
    a = ext[ext["analyte_role"] == "analyte"]
    core = read_table("participant_labs__nhanes", columns=["participant_id", "component_file", "lab_variable",
                                                           "subsample"])
    core["lab_layer"] = core["component_file"].str.rsplit("_", n=1).str[0].map(CORE_LAYER_OF)
    layers = sorted(set(a["lab_layer"]) | set(core["lab_layer"].dropna()))
    for lay in layers:
        d = a[a["lab_layer"] == lay]
        cd = core[core["lab_layer"] == lay]
        files = sorted(set(d["component_file"]) | set(cd["component_file"]))
        main = d[~d["special_sample"]]
        wvars = sorted(set(main["weight_variable"])) if len(main) else []
        if len(cd):
            wvars = sorted(set(wvars) | {re.search(r"weight (\w+)", s).group(1) for s in cd["subsample"].unique()})
        col = f"lab__{lay}"
        cont = d[d["value_type"] == "continuous"]
        lod = cont["below_lod"].dropna()
        rows.append({
            "lab_layer": lay, "layer_family": layer_family(lay),
            "source_table": ("participant_labs__nhanes" if len(cd) and not len(d) else
                             "participant_labs_extended__nhanes" if len(d) and not len(cd) else "both"),
            "files": ";".join(files), "weight_variables_main": ";".join(wvars),
            "n_analytes": int(d["lab_variable"].nunique() + cd["lab_variable"].nunique()),
            "n_participants_all_ages": int(m[col].sum()),
            "n_participants_G": int(m.loc[m["cycle_code"] == "G", col].sum()),
            "n_participants_H": int(m.loc[m["cycle_code"] == "H", col].sum()),
            "n_wearable_valid_adults": int(m.loc[m["wearable_valid_adult"], col].sum()),
            "n_analysis_population": int(pop[col].sum()),
            "coverage_analysis_population": float(pop[col].mean()),
            "coverage_analysis_G": float(pop.loc[pop["cycle_code"] == "G", col].mean()),
            "coverage_analysis_H": float(pop.loc[pop["cycle_code"] == "H", col].mean()),
            "n_participants_special_sample_only": int(len(set(d.loc[d["special_sample"], "participant_id"])
                                                          - set(main["participant_id"]))),
            "n_values": int(len(d) + len(cd)),
            "share_values_below_lod": float(lod.mean()) if len(lod) else np.nan,
            "n_values_with_lc_flag": int(len(lod)),
        })
    for k, (fam, definition) in NONLAB_LAYERS.items():
        rows.append({"lab_layer": k, "layer_family": fam, "source_table": "participant_clinical_features__nhanes / "
                     "participant_wearable_features__nhanes / participant_mortality__nhanes",
                     "files": definition, "n_participants_all_ages": int(m[k].sum()),
                     "n_participants_G": int(m.loc[m["cycle_code"] == "G", k].sum()),
                     "n_participants_H": int(m.loc[m["cycle_code"] == "H", k].sum()),
                     "n_wearable_valid_adults": int(m.loc[m["wearable_valid_adult"], k].sum()),
                     "n_analysis_population": int(pop[k].sum()),
                     "coverage_analysis_population": float(pop[k].mean()),
                     "coverage_analysis_G": float(pop.loc[pop["cycle_code"] == "G", k].mean()),
                     "coverage_analysis_H": float(pop.loc[pop["cycle_code"] == "H", k].mean())})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------
# Overlap map (who has which layers) among wearable-valid adults
# --------------------------------------------------------------------------------------------

FULL_SAMPLE_MIN_COVERAGE = 0.85   # full-sample layer: >= 85% of the analysis population in EACH cycle

FAMILY_ORDER = ["wearable", "questionnaire", "exam", "medications", "mortality", "core clinical chemistry / CBC",
                "metabolic / lipids / glycaemia", "nutritional biomarkers", "hormones", "immune / autoimmune",
                "infectious serology / molecular", "urinary markers", "tobacco exposure", "environmental chemicals",
                "aging / neuro (surplus sera)", "oral microbiome (16S)"]


def layer_families(m: pd.DataFrame) -> dict[str, list[str]]:
    fam: dict[str, list[str]] = {}
    for k, (f, _) in NONLAB_LAYERS.items():
        fam.setdefault(f, []).append(k)
    for c in m.columns:
        if c.startswith("lab__"):
            fam.setdefault(layer_family(c[5:]), []).append(c)
    return fam


def full_sample_layers(lt: pd.DataFrame) -> list[str]:
    lab = lt[lt["lab_layer"].isin([c for c in lt["lab_layer"] if c not in NONLAB_LAYERS])]
    ok = lab[(lab["coverage_analysis_G"] >= FULL_SAMPLE_MIN_COVERAGE) & (lab["coverage_analysis_H"] >= FULL_SAMPLE_MIN_COVERAGE)]
    return sorted(ok["lab_layer"])


def overlap_tables(m: pd.DataFrame, lt: pd.DataFrame) -> pd.DataFrame:
    """UpSet-style exact counts among wearable-valid adults (and the non-pregnant analysis population)."""
    rows = []
    base = {"wearable_valid_adults": m[m["wearable_valid_adult"]], "analysis_population": m[m["analysis_population"]]}
    # 1. each layer
    for r in lt.itertuples():
        col = r.lab_layer if r.lab_layer in NONLAB_LAYERS else f"lab__{r.lab_layer}"
        rows.append({"row_type": "layer", "set": r.lab_layer, "layer_family": r.layer_family,
                     "definition": "has the layer",
                     **{f"n_{k}": int(b[col].sum()) for k, b in base.items()},
                     **{f"pct_{k}": float(b[col].mean()) for k, b in base.items()}})
    # 2. family-level combination patterns (UpSet intersections; exact, mutually exclusive)
    fam = layer_families(m)
    fams = [f for f in FAMILY_ORDER if f in fam]
    for k, b in base.items():
        F = pd.DataFrame({f: b[fam[f]].max(axis=1).astype(int) for f in fams})
        pat = F.apply(lambda r: " & ".join(f for f in fams if r[f]) or "(none)", axis=1)
        vc = pat.value_counts()
        for p_, n in vc.items():
            rows.append({"row_type": f"family_combination__{k}", "set": p_, "layer_family": "",
                         "definition": "exact intersection: has every listed family and none of the others",
                         f"n_{k}": int(n), f"pct_{k}": float(n / len(b)), "n_families": p_.count("&") + 1})
    # 2b. analysis-group combination patterns (the UpSet figure's sets)
    for k, b in base.items():
        G = upset_matrix(b, lt)
        pat = G.apply(lambda r: " & ".join(g for g in G.columns if r[g]) or "(wearable only)", axis=1)
        for p_, n in pat.value_counts().items():
            rows.append({"row_type": f"group_combination__{k}", "set": "wearable & " + p_, "layer_family": "",
                         "definition": "exact intersection of the analysis groups (see UPSET_GROUPS)",
                         f"n_{k}": int(n), f"pct_{k}": float(n / len(b))})
    # 3. cumulative intersections in a fixed order (wearable first, then by decreasing coverage)
    full = full_sample_layers(lt)
    steps = [("wearable_valid", ["wearable_valid"]),
             ("+ questionnaires (PHQ-9, PFQ, self-rated health)",
              ["questionnaire_phq9", "questionnaire_pfq", "questionnaire_self_rated_health"]),
             ("+ exam (BMI, blood pressure)", ["exam_body_measures", "exam_blood_pressure"]),
             ("+ medications interview", ["medications_rx_interview"]),
             ("+ mortality linkage", ["mortality_linked"]),
             ("+ core labs (CBC, biochemistry, HbA1c/lipids, vitamin D/B12)",
              ["lab__core_cbc", "lab__core_biochemistry_profile", "lab__lipids_glycaemia_full", "lab__nutrition_vitamin_d_b12"]),
             ("+ every other full-sample lab layer (" + ", ".join(x for x in full if x not in (
                 "core_cbc", "core_biochemistry_profile", "lipids_glycaemia_full", "nutrition_vitamin_d_b12")) + ")",
              [f"lab__{x}" for x in full]),
             ("+ environmental subsample A core (urine metals, arsenic, perchlorate, VOC metabolites, PAHs, iodine, "
              "serum trace elements)",
              ["lab__env_metals_urine", "lab__env_arsenic_urine", "lab__env_perchlorate_nitrate_thiocyanate_urine",
               "lab__env_voc_metabolites_urine", "lab__env_pahs_urine", "lab__nutrition_iodine_urine",
               "lab__nutrition_trace_elements_serum"]),
             ("+ fasting subsample (glucose/insulin/TG/LDL/ApoB)", ["lab__metabolic_fasting_subsample"]),
             ("+ serum fatty acids", ["lab__nutrition_fatty_acids_serum"]),
             ("+ blood VOCs", ["lab__env_vocs_blood"])]
    for k, b in base.items():
        keep = pd.Series(True, index=b.index)
        for i, (name, cols) in enumerate(steps):
            cols = [c for c in cols if c in b.columns]
            keep &= b[cols].all(axis=1)
            rows.append({"row_type": f"cumulative__{k}", "set": name, "layer_family": "", "step": i,
                         "definition": "has every layer in this and all previous steps",
                         f"n_{k}": int(keep.sum()), f"pct_{k}": float(keep.mean())})
    # 4. named key overlaps
    key = {
        "wearable + questionnaires + exam + meds + mortality + core labs":
            [c for _, cols in steps[:6] for c in cols],
        "wearable + all full-sample lab layers": ["wearable_valid"] + [f"lab__{x}" for x in full],
        "wearable + oral microbiome (2011-2012 only)": ["wearable_valid", "lab__oral_microbiome_16s"],
        "wearable + fasting subsample": ["wearable_valid", "lab__metabolic_fasting_subsample"],
        "wearable + environmental subsample A core": ["wearable_valid"] + steps[7][1],
        "wearable + PFAS (serum)": ["wearable_valid", "lab__env_pfas_serum"],
        "wearable + phthalates + phenols/parabens": ["wearable_valid", "lab__env_phthalates_urine",
                                                      "lab__env_phenols_parabens_urine"],
        "wearable + fatty acids + fasting subsample": ["wearable_valid", "lab__nutrition_fatty_acids_serum",
                                                        "lab__metabolic_fasting_subsample"],
        "wearable + klotho (40-79 y surplus)": ["wearable_valid", "lab__aging_klotho_surplus"],
        "wearable + neurofilament light (2013-2014 surplus)": ["wearable_valid", "lab__neuro_neurofilament_light_surplus"],
        "wearable + autoantibodies (2011-2012 surplus)": ["wearable_valid", "lab__immune_autoantibodies_surplus"],
        "wearable + >= 20 lab layers": None, "wearable + >= 30 lab layers": None,
    }
    for name, cols in key.items():
        for k, b in base.items():
            if cols is None:
                thr = int(re.search(r">= (\d+)", name).group(1))
                sel = b["n_lab_layers"] >= thr
            else:
                sel = b[[c for c in cols if c in b.columns]].all(axis=1)
            rows.append({"row_type": f"key_overlap__{k}", "set": name, "layer_family": "",
                         "definition": "has all listed layers", f"n_{k}": int(sel.sum()), f"pct_{k}": float(sel.mean())})
    for k, b in base.items():
        for thr in (5, 10, 15, 20, 25, 30, 35, 40):
            rows.append({"row_type": f"n_lab_layers_at_least__{k}", "set": f">= {thr} lab layers", "layer_family": "",
                         "definition": "number of lab layers with >= 1 measured value",
                         f"n_{k}": int((b["n_lab_layers"] >= thr).sum()),
                         f"pct_{k}": float((b["n_lab_layers"] >= thr).mean())})
    out = pd.DataFrame(rows)
    return out


# Analysis-relevant groups for the UpSet view (each = has ALL listed layers, except "any" groups).
UPSET_GROUPS = {
    "questionnaires (PHQ-9 + PFQ + self-rated health)": ("all", ["questionnaire_phq9", "questionnaire_pfq",
                                                                  "questionnaire_self_rated_health"]),
    "exam + meds + mortality": ("all", ["exam_body_measures", "exam_blood_pressure", "medications_rx_interview",
                                        "mortality_linked"]),
    "core labs (CBC, biochemistry, HbA1c/lipids, vit D/B12)": ("all", ["lab__core_cbc", "lab__core_biochemistry_profile",
                                                                      "lab__lipids_glycaemia_full",
                                                                      "lab__nutrition_vitamin_d_b12"]),
    "other full-sample labs (7 layers)": ("all_full_other", []),
    "env. subsample A core (7 layers)": ("all", ["lab__env_metals_urine", "lab__env_arsenic_urine",
                                                 "lab__env_perchlorate_nitrate_thiocyanate_urine",
                                                 "lab__env_voc_metabolites_urine", "lab__env_pahs_urine",
                                                 "lab__nutrition_iodine_urine", "lab__nutrition_trace_elements_serum"]),
    "fasting subsample": ("all", ["lab__metabolic_fasting_subsample"]),
    "serum fatty acids": ("all", ["lab__nutrition_fatty_acids_serum"]),
    "blood VOCs": ("all", ["lab__env_vocs_blood"]),
    "oral microbiome 16S (2011-12)": ("all", ["lab__oral_microbiome_16s"]),
    "surplus sera (klotho / NfL / autoantibodies; any)": ("any", ["lab__aging_klotho_surplus",
                                                                  "lab__neuro_neurofilament_light_surplus",
                                                                  "lab__immune_autoantibodies_surplus"]),
}


def upset_matrix(b: pd.DataFrame, lt: pd.DataFrame) -> pd.DataFrame:
    full = [f"lab__{x}" for x in full_sample_layers(lt) if x not in (
        "core_cbc", "core_biochemistry_profile", "lipids_glycaemia_full", "nutrition_vitamin_d_b12")]
    out = {}
    for name, (how, cols) in UPSET_GROUPS.items():
        if how == "all_full_other":
            cols = full
        cols = [c for c in cols if c in b.columns]
        out[name] = (b[cols].all(axis=1) if how.startswith("all") else b[cols].any(axis=1)).astype(int)
    return pd.DataFrame(out)


def overlap_figure(m: pd.DataFrame, lt: pd.DataFrame, path: Path, top: int = 25) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    b = m[m["analysis_population"]]
    G = upset_matrix(b, lt)
    groups = list(G.columns)
    pat = G.apply(lambda r: tuple(int(r[g]) for g in groups), axis=1)
    vc = pat.value_counts().head(top)
    fig = plt.figure(figsize=(15, 15))
    gs = fig.add_gridspec(3, 2, width_ratios=[2.2, 4], height_ratios=[1.6, 1.8, 3.2], wspace=0.03, hspace=0.22)
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1])
    ax_set = fig.add_subplot(gs[1, 0])
    ax_lay = fig.add_subplot(gs[2, :])
    x = np.arange(len(vc))
    ax_bar.bar(x, vc.to_numpy(), color="#3b6ea5")
    for xi, v in zip(x, vc.to_numpy()):
        ax_bar.text(xi, v, f"{v:,}", ha="center", va="bottom", fontsize=6.5, rotation=90)
    ax_bar.set_xlim(-0.6, len(x) - 0.4)
    ax_bar.set_ylabel("participants\n(exact intersection)", fontsize=8)
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bar.set_title(f"NHANES 2011-2014 wearable-valid adults (analysis population, n = {len(b):,}): which data "
                     f"layers overlap (top {top} exact intersections; every row also has wrist accelerometry; draft)",
                     fontsize=10)
    for s_ in ("top", "right"):
        ax_bar.spines[s_].set_visible(False)
    for j in range(len(groups)):
        ax_mat.scatter(x, np.full(len(x), j), color="#dddddd", s=28, zorder=1)
    for xi, p_ in zip(x, vc.index):
        on = [j for j, v in enumerate(p_) if v]
        ax_mat.scatter(np.full(len(on), xi), on, color="#222222", s=28, zorder=2)
        if on:
            ax_mat.plot([xi, xi], [min(on), max(on)], color="#222222", lw=1, zorder=2)
    ax_mat.set_xlim(-0.6, len(x) - 0.4)
    ax_mat.set_ylim(len(groups) - 0.5, -0.5)
    ax_mat.set_yticks([])
    ax_mat.tick_params(axis="x", bottom=False, labelbottom=False)
    tot = G.sum().to_numpy()
    ax_set.barh(range(len(groups)), tot, color="#888888")
    ax_set.set_ylim(len(groups) - 0.5, -0.5)
    ax_set.set_yticks(range(len(groups)))
    ax_set.set_yticklabels([f"{g}  ({t:,})" for g, t in zip(groups, tot)], fontsize=7.5)
    ax_set.set_xlabel("participants with group", fontsize=8)
    lab = lt[~lt["lab_layer"].isin(list(NONLAB_LAYERS))].sort_values("n_analysis_population", ascending=False)
    fams = list(dict.fromkeys(lab["layer_family"]))
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(fams), 2)))
    col = {f: cmap[i] for i, f in enumerate(fams)}
    xx = np.arange(len(lab))
    ax_lay.bar(xx, lab["n_analysis_population"], color=[col[f] for f in lab["layer_family"]])
    ax_lay.set_xticks(xx)
    ax_lay.set_xticklabels(lab["lab_layer"], rotation=75, ha="right", fontsize=6.5)
    ax_lay.set_ylabel("wearable-valid adults with the layer", fontsize=8)
    ax_lay.axhline(len(b), color="grey", lw=0.7, ls="--")
    ax_lay.set_title("Each laboratory layer: participants in the analysis population with >= 1 measured value "
                     "(colour = family)", fontsize=9)
    from matplotlib.patches import Patch
    ax_lay.legend(handles=[Patch(color=col[f], label=f) for f in fams], fontsize=7, ncol=3, loc="upper right",
                  frameon=False)
    for s_ in ("top", "right"):
        ax_lay.spines[s_].set_visible(False)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------------------------
# DATA_AUDIT section, registry extension, entry point
# --------------------------------------------------------------------------------------------

AUDIT_SECTION_FILE = "labs_extended_audit_section.md"
AUDIT_BEGIN = "<!-- BEGIN labs_extended (generated by measure_it.ingestion.nhanes_labs_extended) -->"
AUDIT_END = "<!-- END labs_extended -->"


def _fmt(x) -> str:
    if isinstance(x, float):
        return "" if np.isnan(x) else f"{x:.3f}"
    return str(x)


def audit_section(inv: pd.DataFrame, vd: pd.DataFrame, lt: pd.DataFrame, ov: pd.DataFrame, ext: pd.DataFrame) -> str:
    L = [AUDIT_BEGIN, "", "## Extended laboratory catalogue (all 2011-2014 laboratory files) — layers",
         "",
         f"_Generated by `measure_it.ingestion.nhanes_labs_extended` on {utc_now_iso()}; every count is computed from "
         "the retrieved files. Tables: `participant_labs_extended__nhanes` (long), `participant_layer_membership__nhanes` "
         "(0/1 per layer), `results/tables/nhanes_lab_inventory.csv`, `nhanes_lab_layers.csv`, `nhanes_lab_variables.csv`, "
         "`nhanes_layer_overlap.csv`._", "",
         "Listing pages read (saved under `docs/listings/`): "
         + ", ".join(f"`{LISTING_URL.format(year=y)}`" for y in CYCLES.values())
         + "; the Oral Microbiome Project page `https://wwwn.cdc.gov/Nchs/Nhanes/Omp/Default.aspx` (linked from both "
         "listings) and its text files under `https://wwwn.cdc.gov/nchs/data/nhanes/omp/`.", ""]
    dc = inv["decision"].value_counts()
    L.append(f"Files listed: {len(inv)} ({dc.get('ingested', 0)} ingested here, {dc.get('core', 0)} already in "
             f"`participant_labs__nhanes`, {dc.get('excluded', 0)} excluded with a reason). Participants with >= 1 "
             f"extended value: {ext['participant_id'].nunique():,}; rows: {len(ext):,}.")
    L += ["", "### Every listed file and the decision", "",
          "| cycle | file | description | published | decision | lab_layer | reason |", "|---|---|---|---|---|---|---|"]
    for r in inv.sort_values(["cycle_code", "file"]).itertuples():
        L.append(f"| {r.cycle_code} | {r.file} | {r.description} | {r.published} | {r.decision} | {r.lab_layer} | "
                 f"{r.reason} |")
    L += ["", "### Layers: files, subsample weights and measured counts", "",
          "`n_analysis` = participants in the modelling population (wearable-valid, aged >= 18, not pregnant; n = "
          f"{int(ov.loc[(ov.row_type == 'layer') & (ov.set == 'wearable_valid'), 'n_analysis_population'].iloc[0]):,}) "
          "with >= 1 measured value in the layer; coverage G/H = share of that population per cycle. `special-only` = "
          "participants present only in the special smoking-sample file (weight WTFSM), not in the main subsample "
          "file. `% below LOD` = share of released continuous values whose comment code marks them below the lower "
          f"detection limit. Full-sample layers (coverage >= {FULL_SAMPLE_MIN_COVERAGE:.0%} in each cycle): "
          + ", ".join(f"`{x}`" for x in full_sample_layers(lt)) + ".", "",
          "| layer | family | files | main weight(s) | analytes | n all ages | n wearable-valid adults | n_analysis | "
          "coverage G | coverage H | special-only | % below LOD |", "|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in lt.sort_values("n_analysis_population", ascending=False).itertuples():
        L.append(f"| {r.lab_layer} | {r.layer_family} | {r.files} | {_fmt(r.weight_variables_main)} | "
                 f"{'' if pd.isna(r.n_analytes) else int(r.n_analytes)} | {r.n_participants_all_ages:,} | "
                 f"{r.n_wearable_valid_adults:,} | {r.n_analysis_population:,} | {r.coverage_analysis_G:.3f} | "
                 f"{r.coverage_analysis_H:.3f} | "
                 f"{'' if pd.isna(r.n_participants_special_sample_only) else int(r.n_participants_special_sample_only)} | "
                 f"{'' if pd.isna(r.share_values_below_lod) else f'{100 * r.share_values_below_lod:.1f}'} |")
    rc = vd["role"].value_counts()
    L += ["", "### Variable handling (detection limits, units, roles)", "",
          "* **Detection limits.** NCHS releases below-LLOD results already filled with LLOD/sqrt(2) and flags them "
          "with a comment-code variable (0 = at or above the detection limit, 1 = below the lower limit). The released "
          "value is kept unchanged in `value`; the matched comment code is in `lc_variable` / `lc_code`, and "
          "`below_lod` = (`lc_code` == 1). Comment codes are matched to analytes by the NHANES naming rule "
          "(`<P>X<core>` <-> `<P>D<core>LC`, truncated to 8 characters), then by label (comment-code label minus "
          "'comment code' = analyte label minus units, or the analyte's parenthesised abbreviation, e.g. URDMEALC "
          "'NMEA' -> URXNMEA 'N-Nitrosoethylmethylamine (NMEA)'), then by a unique name prefix; "
          f"{rc.get('comment_code', 0)} comment codes matched, {rc.get('comment_code_unmatched', 0)} unmatched. The "
          "serum neurofilament file flags the upper limit separately (SSNFLH, not ingested) and the lower limit in "
          "SSNFLL (used as the comment code).",
          "* **Units.** When a file carries both conventional and SI versions of the same analyte (label identical up "
          f"to the unit), only the first-listed (conventional) version is ingested ({rc.get('si_unit_duplicate_not_ingested', 0)} SI "
          "duplicates dropped); analytes published only in SI units (e.g. methylmalonic acid, serum folate forms) are kept.",
          f"* **Roles.** {rc.get('analyte', 0)} analyte variables; {rc.get('urine_creatinine', 0)} urinary creatinine "
          f"variables (kept, role `urine_creatinine`, for dilution adjustment); {rc.get('process_variable', 0)} "
          "specimen-process variables kept with role `process_variable` (urine collection volumes/times, the HDV "
          f"re-testing flag); {rc.get('questionnaire_item_not_ingested', 0)} questionnaire items bundled in the blood-VOC "
          f"files (VTQ*) and {rc.get('character_variable_not_ingested', 0)} character variable not ingested; "
          f"{rc.get('sample_weight', 0)} sample-weight variables carried per row as `weight_variable` / `weight_value`.",
          "* **Categorical results** (serology, STI, HPV, celiac, TB IGRA) keep the NCHS code in `value`, its label in "
          "`value_label`, and `result_positive` = 1 (Positive / Weakly positive / Reactive) or 0 (Negative); "
          "indeterminate / not evaluated codes stay NaN.",
          "* **Subsamples.** Each file's weight variable identifies the sampled people: WTSA2YR / WTSB2YR / WTSC2YR "
          "(one-third environmental subsamples A/B/C), WTSAF2YR (morning fasting), WTSOG2YR (OGTT), WTFAS2YR (fatty "
          "acids), WTSVOC2Y / WTSVS2YR (blood VOCs), WTSH2YR (2013-2014 blood metals, a half sample), WTALD2YR / "
          "WTALDS2Y (aldehydes), WTANA2YR (autoantibodies), WTFSM (special smoking sample), WTSS* (surplus-specimen "
          "projects). A participant not sampled for a component has no rows for it; nothing is imputed here. Files "
          "without a weight are measured on the full MEC sample within the component's age/sex eligibility (e.g. "
          "HSV 14-49, HIV 18-59, oral HPV 14-69, klotho 40-79 surplus); those restrictions are in each codebook's "
          "Target field (saved `docs/*.htm`) and in `nhanes_lab_variables.csv`.",
          "* **Special smoking sample** (`*S_G`/`*S_H` files, weight WTFSM): overlaps the main subsample-A files "
          "heavily (same participants measured once, released twice with different weights) and adds participants "
          "who are only in the special sample (column `special-only` above); rows are flagged `special_sample`.",
          "* **Oral microbiome (2009-2012; only the 2011-2012 participants are in this project).** Alpha diversity "
          "(observed ASVs, Faith's phylogenetic diversity, Shannon-Wiener, inverse Simpson) at a rarefaction depth "
          "of 10,000 reads averaged over NCHS's 10 rarefactions, and genus relative abundance (DADA2, rarefied, SILVA "
          "v123 annotation) for the genera present in >= 5% of the 2011-2012 samples (zeros are kept as measured "
          "zeros). Participants aged 14-69 with a sequenced oral rinse.",
          "* **Not individual-level / not ingested**: pooled-sample files (BFRPOL, PCBPOL, PSTPOL, DOXPOL, POOLTF), "
          "withdrawn files (SSEVD_G/H, PAHS_H), the fasting questionnaire, the urine pregnancy test and household "
          "tap-water fluoride (see the decision table).", "",
          "### Overlap with wearable data (measured)", "",
          "Cumulative intersection among wearable-valid adults (each step keeps people who have every layer so far):", "",
          "| step | n wearable-valid adults | n analysis population |", "|---|---:|---:|"]
    cw = ov[ov.row_type == "cumulative__wearable_valid_adults"].reset_index(drop=True)
    ca = ov[ov.row_type == "cumulative__analysis_population"].reset_index(drop=True)
    for i in range(len(cw)):
        L.append(f"| {cw.loc[i, 'set']} | {int(cw.loc[i, 'n_wearable_valid_adults']):,} | "
                 f"{int(ca.loc[i, 'n_analysis_population']):,} |")
    L += ["", "| key overlap | n wearable-valid adults | n analysis population |", "|---|---:|---:|"]
    kw = ov[ov.row_type == "key_overlap__wearable_valid_adults"].reset_index(drop=True)
    ka = ov[ov.row_type == "key_overlap__analysis_population"].reset_index(drop=True)
    for i in range(len(kw)):
        L.append(f"| {kw.loc[i, 'set']} | {int(kw.loc[i, 'n_wearable_valid_adults']):,} | "
                 f"{int(ka.loc[i, 'n_analysis_population']):,} |")
    L += ["", "Limitations: 2011-2014 has no CRP/hs-CRP and no ferritin (not on either listing); several layers exist "
          "in one cycle only (oral microbiome, autoantibodies, TB IGRA, urine osmolality: 2011-2012; tobacco "
          "nitrosamines, aromatic/heterocyclic amines, aldehydes, Hb adducts, NfL, surplus PFAS/phthalates/glyphosate/"
          "terpenes: 2013-2014); subsample layers must be analysed only among sampled people; values below the "
          "detection limit are substitution values, not measurements.", "", AUDIT_END, ""]
    return "\n".join(L)


def splice_audit(section: str) -> Path:
    d = raw_dir(SOURCE_ID)
    (d / AUDIT_SECTION_FILE).write_text(section)
    p = d / "DATA_AUDIT.md"
    text = p.read_text() if p.exists() else ""
    if AUDIT_BEGIN in text:
        pre = text.split(AUDIT_BEGIN)[0]
        post = text.split(AUDIT_END, 1)[1] if AUDIT_END in text else ""
        text = pre.rstrip() + "\n\n" + section + post.lstrip("\n")
    else:
        text = text.rstrip() + "\n\n" + section
    p.write_text(text)
    return p


def registry_extension(entry: dict) -> dict:
    """Merge the extended-lab facts into the nhanes_2011_2014 registry entry (idempotent)."""
    from ..store import processed_path
    meta_ext = processed_path("participant_labs_extended__nhanes").with_suffix(".meta.json")
    if not meta_ext.exists():
        return entry
    lt_path = TABLES / "nhanes_lab_layers.csv"
    lt = pd.read_csv(lt_path) if lt_path.exists() else pd.DataFrame()
    e = dict(entry)
    outs = list(e.get("processed_outputs", []))
    for t in ("participant_labs_extended__nhanes", "participant_layer_membership__nhanes"):
        if t not in outs:
            outs.append(t)
    e["processed_outputs"] = outs
    urls = list(e.get("access_urls", []))
    for u in ("https://wwwn.cdc.gov/nchs/data/nhanes/omp/",
              LISTING_URL.format(year=2011), LISTING_URL.format(year=2013)):
        if u not in urls:
            urls.append(u)
    e["access_urls"] = urls
    mod = str(e.get("ingestion_module", ""))
    if "nhanes_labs_extended" not in mod:
        e["ingestion_module"] = mod + " (+ measure_it.ingestion.nhanes_labs_extended: full laboratory catalogue)"
    ss = dict(e.get("sample_size", {}))
    if len(lt):
        lab = lt[~lt["lab_layer"].isin(list(NONLAB_LAYERS))]
        ss["lab_layers"] = int(len(lab))
        ss["lab_layer_participants_analysis_population"] = {r.lab_layer: int(r.n_analysis_population)
                                                             for r in lab.itertuples()}
    rows = json.loads(meta_ext.read_text())["rows"]
    ss["extended_lab_rows"] = int(rows)
    e["sample_size"] = ss
    lim = list(e.get("limitations", []))
    extra = ["Extended lab layers: many are one-third or surplus subsamples (own weights) or one cycle only; analyse "
             "only among sampled participants; below-LOD values are NCHS fill values (LLOD/sqrt(2))",
             "Oral microbiome (16S) only for 2011-2012 participants aged 14-69"]
    for x in extra:
        if x not in lim:
            lim.append(x)
    e["limitations"] = lim
    e["omics"] = e.get("omics", False)
    e["omics_note"] = ("oral microbiome 16S (NHANES Oral Microbiome Project, SEQN-linked) for 2011-2012; lab panels "
                       "are biomarker layers, not omics")
    return e


def update_registry() -> Path:
    import yaml
    from ..registry import write_registry_entry
    p = raw_dir(SOURCE_ID) / "registry_entry.yaml"
    entry = yaml.safe_load(p.read_text())
    return write_registry_entry(registry_extension(entry))


def run(log=print) -> dict:
    TABLES.mkdir(parents=True, exist_ok=True)
    inv = build_inventory(log=log)
    fetch_files(inv, log=log)
    ext, vd = build_extended(inv, log=log)
    write_extended(ext)
    m = build_membership(ext, log=log)
    write_membership(m)
    lt = layer_table(ext, m, inv)
    inv.to_csv(TABLES / "nhanes_lab_inventory.csv", index=False)
    lt.to_csv(TABLES / "nhanes_lab_layers.csv", index=False)
    vd.to_csv(TABLES / "nhanes_lab_variables.csv", index=False)
    ov = overlap_tables(m, lt)
    ov.to_csv(TABLES / "nhanes_layer_overlap.csv", index=False)
    from ..config import FIGURES
    (FIGURES / "drafts").mkdir(parents=True, exist_ok=True)
    overlap_figure(m, lt, FIGURES / "drafts" / "nhanes_layer_overlap.png")
    splice_audit(audit_section(inv, vd, lt, ov, ext))
    update_registry()
    log("[labs_extended] done")
    return {"inventory": inv, "extended": ext, "variables": vd, "membership": m, "layers": lt, "overlap": ov}


def _cli() -> None:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args()
    run()


if __name__ == "__main__":
    _cli()
