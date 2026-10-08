"""Ingest the fibromyalgia infrared-thermography dataset (PLOS ONE 2021, doi:10.1371/journal.pone.0253281, S1 Data).

Sempere-Rubio N, Aguilar-Rodriguez M, Ingles M, Izquierdo-Alventosa R, Serra-Ano P. "Thermal imaging ruled out as
a supplementary assessment in patients with fibromyalgia: A cross-sectional study." PLOS ONE 2021;16(6):e0253281.

One row per woman: 86 with fibromyalgia (ACR 2010, assessed by rheumatologists; IDs FM###, Group 1) and 92
age-matched women without FM symptoms (IDs CG###, Group 2); all post-menopausal (paper Methods). Resting skin
temperature (degC; minimum, maximum and average over each region of interest) from a FLIR E60bx camera at six
regions: neck, upper back, lower back, chest, medial knee (mean of both), lateral elbow (mean of both).

Outputs (person layer; no geography):
  participants__fm_thermography                 one row per participant (label, age, BMI, raw coded covariates)
  participant_device_features__fm_thermography  one row per participant, 18 temperature columns
Raw files, DATA_AUDIT.md and registry_entry.yaml under data/raw/fm_thermography/.

Run: uv run python -m measure_it.ingestion.fm_thermography
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from measure_it.config import UNKNOWN, raw_dir
from measure_it.download import download_file, load_manifest
from measure_it.provenance import add_provenance
from measure_it.registry import write_registry_entry
from measure_it.store import write_table

SOURCE_ID = "fm_thermography"
DATASET_ID = "fm_thermography"
DATA_URL = "https://ndownloader.figshare.com/files/28433660"
DATA_FILE = "pone.0253281.s001.xlsx"
ARTICLE_URL = "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0253281&type=manuscript"
ARTICLE_FILE = "pone.0253281.xml"
LANDING_URL = "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0253281"
DOI = "10.1371/journal.pone.0253281"
SOURCE_NAME = "PLOS ONE 0253281 S1 Data: infrared thermography in fibromyalgia vs controls"
CITATION = ("Sempere-Rubio N, Aguilar-Rodriguez M, Ingles M, Izquierdo-Alventosa R, Serra-Ano P. Thermal imaging "
            "ruled out as a supplementary assessment in patients with fibromyalgia: a cross-sectional study. "
            "PLOS ONE 2021;16(6):e0253281. doi:10.1371/journal.pone.0253281")
CONDITION_ID = "fibromyalgia"
MEASUREMENT_CLASS_ID = "continuous_temperature"
MEASUREMENT_CLASS_NOTE = ("nearest class in configs/measurements.yaml by SIGNAL (skin temperature); that class is a "
                          "wearable, this is a single clinic infrared-camera image at rest. Proposed new class: "
                          "infrared_thermography (not added; configs/measurements.yaml is shared).")
DEVICE = "FLIR E60bx infrared thermography camera (320x240, 0.045 degC sensitivity), 1 m, emissivity 0.98"
LABEL_BASIS = ("fibromyalgia = ACR 2010 criteria assessed by rheumatologists (paper Methods), recruited from FM "
               "associations and specialised units; control = age-matched women without FM symptoms recruited by "
               "advertisement/snowball sampling. Group 1 = FM (IDs FM###), Group 2 = control (IDs CG###), verified "
               "from the ID prefixes.")
PROTOCOL = ("resting, standing in underwear after 15 min acclimatisation; 8 m2 room held at 24 degC, 44% RH; "
            "3-6 pm; University of Valencia biomechanics lab (paper Methods). Per-session room temperature, date and "
            "season are not in the file.")

REGIONS = ["Neck", "Upper_back", "Lower_back", "Chest", "Knee", "Elbow"]
STATS = ["min", "max", "ave"]
RAW_TEMP_COLS = [f"{r}_{s}" for r in REGIONS for s in STATS]
AVE_COLS = [f"{r}_ave" for r in REGIONS]
CODED_COLS = ["Civil_status", "Employment_status", "Studies", "Tobacco", "Alcohol"]
EXPECTED_COLS = ["Participant", "Group", "Age", "BMI"] + RAW_TEMP_COLS + CODED_COLS
GROUP_MAP = {1: "fibromyalgia", 2: "control"}
# Code labels: paper Table 1 lists the categories in this order, and the per-group counts of codes 1..k in the file
# reproduce every Table 1 cell exactly (verified 2026-09-24; checked again in `verify_table1`).
CODE_LABELS = {
    "Civil_status": ["married/in union", "single", "widowed", "separated/divorced"],
    "Employment_status": ["active", "unemployed", "pensioner", "housekeeper"],
    "Studies": ["illiterate", "primary", "secondary", "university"],
    "Tobacco": ["daily", "not daily", "not at present", "never"],
    "Alcohol": ["daily", "some days per week", "some days per month", "some days in the last year",
                "not in the last year"],
}
# Paper Table 1 counts (FM, control) per category, in CODE_LABELS order.
TABLE1_COUNTS = {
    "Civil_status": [(68, 66), (6, 11), (1, 4), (11, 11)],
    "Employment_status": [(44, 60), (18, 14), (21, 14), (3, 4)],
    "Studies": [(1, 2), (37, 28), (36, 26), (12, 36)],
    "Tobacco": [(19, 14), (3, 7), (31, 30), (33, 41)],
    "Alcohol": [(6, 5), (25, 34), (15, 21), (16, 14), (24, 18)],
}


def feature_name(raw_col: str) -> str:
    """`Upper_back_ave` -> `skin_temp_upper_back_ave_c`."""
    return "skin_temp_" + raw_col.lower() + "_c"


FEATURE_COLS = [feature_name(c) for c in RAW_TEMP_COLS]
AVE_FEATURES = [feature_name(c) for c in AVE_COLS]


def participant_id(native_id: str) -> str:
    return f"{DATASET_ID}:{native_id}"


def fetch() -> tuple:
    data = download_file(DATA_URL, SOURCE_ID, DATA_FILE)
    article = download_file(ARTICLE_URL, SOURCE_ID, ARTICLE_FILE)
    return data, article


def read_raw(path=None) -> pd.DataFrame:
    path = path or raw_dir(SOURCE_ID) / DATA_FILE
    xl = pd.ExcelFile(path)
    if xl.sheet_names != ["Hoja5"]:
        raise ValueError(f"unexpected sheets {xl.sheet_names}")
    df = pd.read_excel(xl, "Hoja5")
    if list(df.columns) != EXPECTED_COLS:
        raise ValueError(f"unexpected columns: {list(df.columns)}")
    return df


def derive_label(df: pd.DataFrame) -> pd.Series:
    """Group code -> label, and fail if any ID prefix disagrees with the group code."""
    prefix = df["Participant"].astype(str).str.extract(r"^([A-Z]+)")[0]
    expected = df["Group"].map({1: "FM", 2: "CG"})
    if expected.isna().any():
        raise ValueError(f"unknown Group codes: {sorted(df.loc[expected.isna(), 'Group'].unique())}")
    if not (prefix == expected).all():
        raise ValueError("Group code disagrees with Participant ID prefix")
    return df["Group"].map(GROUP_MAP)


def verify_table1(df: pd.DataFrame) -> dict:
    """Per-group code counts must equal paper Table 1; this also checks the Group coding independently of the
    temperatures. Returns {variable: True}; raises on any mismatch."""
    label = derive_label(df)
    out = {}
    for var, counts in TABLE1_COUNTS.items():
        for code, (n_fm, n_c) in enumerate(counts, start=1):
            got = (int(((df[var] == code) & (label == "fibromyalgia")).sum()),
                   int(((df[var] == code) & (label == "control")).sum()))
            if got != (n_fm, n_c):
                raise ValueError(f"{var} code {code}: file {got} != paper Table 1 {(n_fm, n_c)}")
        out[var] = True
    return out


def article_text(path=None) -> str:
    path = path or raw_dir(SOURCE_ID) / ARTICLE_FILE
    raw = re.sub(r"<!DOCTYPE[^>]*>", "", path.read_text())
    return re.sub(r"\s+", " ", "".join(ET.fromstring(raw).itertext()))


def paper_facts(text: str) -> dict:
    """Facts checked in the article text (fail loudly if the text no longer says them)."""
    checks = {
        "post_menopausal": "postmenopausal" in text,
        "acr_2010": "ACR, 2010" in text,
        "room_24C": "kept at 24" in text,
        "camera_flir_e60bx": "FLIR E60BX" in text,
        "n_fm_86": "86 women with fibromyalgia" in text,
        "authors_conclusion_negative": "not an effective supplementary assessment tool" in text,
    }
    missing = [k for k, v in checks.items() if not v]
    if missing:
        raise ValueError(f"article text no longer states: {missing}")
    return checks


def build_frames(df: pd.DataFrame, retrieved_at: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    label = derive_label(df)
    nid = df["Participant"].astype(str)
    base = pd.DataFrame({
        "participant_id": nid.map(participant_id),
        "dataset_id": DATASET_ID,
        "native_id": nid,
    })
    n_missing = df[RAW_TEMP_COLS].isna().sum(axis=1)
    part = base.assign(
        group_label=label.values,
        source_group_code=df["Group"].astype(int).values,
        fibromyalgia=(label == "fibromyalgia").astype(int).values,
        condition_id=np.where(label == "fibromyalgia", CONDITION_ID, ""),
        label_basis=LABEL_BASIS,
        sex="female",
        sex_basis="paper Methods: all participants are women (no sex column in the file)",
        menopausal_status="post-menopausal",
        age_years=df["Age"].astype(float).values,
        bmi=df["BMI"].astype(float).values,
        civil_status_code=df["Civil_status"].astype(int).values,
        employment_status_code=df["Employment_status"].astype(int).values,
        studies_code=df["Studies"].astype(int).values,
        tobacco_code=df["Tobacco"].astype(int).values,
        alcohol_code=df["Alcohol"].astype(int).values,
        civil_status=df["Civil_status"].map(lambda c: CODE_LABELS["Civil_status"][c - 1]).values,
        employment_status=df["Employment_status"].map(lambda c: CODE_LABELS["Employment_status"][c - 1]).values,
        education=df["Studies"].map(lambda c: CODE_LABELS["Studies"][c - 1]).values,
        tobacco=df["Tobacco"].map(lambda c: CODE_LABELS["Tobacco"][c - 1]).values,
        alcohol=df["Alcohol"].map(lambda c: CODE_LABELS["Alcohol"][c - 1]).values,
        coded_covariate_labels=("code labels from paper Table 1 (category order); per-group counts reproduce every "
                                "Table 1 cell exactly"),
        n_temperature_values_missing=n_missing.values,
        complete_ave_temperatures=df[AVE_COLS].notna().all(axis=1).values,
    )
    part = add_provenance(
        part, data_layer="person", source_name=SOURCE_NAME, source_version=f"figshare file 28433660 ({DATA_FILE})",
        retrieved_at=retrieved_at, evidence_type="person_exam_measurement", source_record_id="native_id",
        evidence_level="clinical diagnosis (ACR 2010) vs self-described healthy controls",
        provenance_notes=f"{CITATION}. Label: {LABEL_BASIS}")
    feats = base.copy()
    for c in RAW_TEMP_COLS:
        feats[feature_name(c)] = df[c].astype(float).values
    feats = feats.assign(
        device=DEVICE,
        device_class="infrared_thermography_clinic",
        measurement_class_id=MEASUREMENT_CLASS_ID,
        measurement_class_note=MEASUREMENT_CLASS_NOTE,
        measurement_protocol=PROTOCOL,
        units="degC",
        n_temperature_values_missing=n_missing.values,
    )
    feats = add_provenance(
        feats, data_layer="person", source_name=SOURCE_NAME, source_version=f"figshare file 28433660 ({DATA_FILE})",
        retrieved_at=retrieved_at, evidence_type="person_device_measurement", source_record_id="native_id",
        evidence_level="device measurement as published (ROI summary statistics; images not released)",
        provenance_notes=f"{CITATION}. {PROTOCOL}")
    return part, feats


def missingness(df: pd.DataFrame) -> pd.DataFrame:
    label = derive_label(df)
    rows = []
    for c in ["Age", "BMI"] + RAW_TEMP_COLS + CODED_COLS:
        rows.append({"variable": c, "n_missing": int(df[c].isna().sum()),
                     "pct_missing": round(100 * df[c].isna().mean(), 2),
                     "n_missing_fm": int(df.loc[label == "fibromyalgia", c].isna().sum()),
                     "n_missing_control": int(df.loc[label == "control", c].isna().sum())})
    return pd.DataFrame(rows)


def write_audit(df: pd.DataFrame, part: pd.DataFrame, feats: pd.DataFrame, manifest: dict) -> None:
    miss = missingness(df)
    label = derive_label(df)
    n_fm, n_cg = int((label == "fibromyalgia").sum()), int((label == "control").sum())
    cc = part[part.complete_ave_temperatures]
    files = "\n".join(
        f"| `{k}` | {v.get('url')} | {v.get('bytes'):,} | `{v.get('sha256')}` | {v.get('retrieved_at')} |"
        for k, v in sorted(manifest.get("files", {}).items()))
    miss_nonzero = miss[miss.n_missing > 0]
    miss_md = "\n".join(f"| `{r.variable}` | {r.n_missing} ({r.pct_missing}%) | {r.n_missing_fm} | "
                        f"{r.n_missing_control} |" for r in miss_nonzero.itertuples())
    text = f"""# DATA AUDIT — Infrared thermography in fibromyalgia (PLOS ONE 0253281, S1 Data)

| Field | Value |
|---|---|
| source_id | `{SOURCE_ID}` |
| Source (dataset/API name, exact files/endpoints) | S1 Data `{DATA_FILE}` ({DATA_URL}); article JATS XML ({ARTICLE_URL}) |
| Publishing organization | PLOS ONE (data deposited on figshare by PLOS); authors at the University of Valencia, Department of Physiotherapy |
| Retrieval date (UTC) | {manifest['files'][DATA_FILE]['retrieved_at']} |
| Source version / release | figshare file 28433660, Last-Modified {manifest['files'][DATA_FILE].get('last_modified')}; article published 2021-06-24 |
| Source update date / cadence | static supplement (no updates) |
| License / access conditions | CC BY 4.0 (article licence covers the supplement); open download, no registration, no DUA |
| Unit of observation | participant (one resting thermography session) |
| Sample size (actual, as ingested) | {len(df)} women: {n_fm} fibromyalgia, {n_cg} controls; complete on the six average temperatures: {int(cc.fibromyalgia.sum())} FM, {int((1 - cc.fibromyalgia).sum())} controls |
| Geography (resolution, vintage) | none (single site, Valencia, Spain; no participant geography) |
| Person-level? | yes |
| Geographic? | no |
| Omics? | no |
| Wearable? | no — a clinic infrared camera (device measurement, not worn) |
| True participant linkage across modalities? | no — one modality (thermography ROI summaries + age/BMI + coded covariates). The same group published a force-plate FM set (PLOS ONE 0196575); a linkage between the two was **not** verified and is not assumed. |

## Files / endpoints retrieved

| file | url | bytes | sha256 | retrieved_at |
|---|---|---|---|---|
{files}

The xlsx has one sheet (`Hoja5`), 178 rows x 27 columns, no codebook. The data dictionary below is from the article
Methods (retrieved XML), not from the file.

## Key variables

| variable (raw) | meaning | units | downstream |
|---|---|---|---|
| `Participant` | native ID (`FM###` / `CG###`) | — | `participant_id = "{DATASET_ID}:<Participant>"` |
| `Group` | 1 = fibromyalgia, 2 = control (verified against ID prefix for all rows) | code | `fibromyalgia` (1/0), `group_label` |
| `Age` | age | years | covariate (model 1, S2) |
| `BMI` | body-mass index | kg/m2 | covariate (model 1, S2) |
| `<Region>_{{min,max,ave}}` | min / max / mean skin temperature within the region of interest; Region = Neck, Upper_back, Lower_back, Chest, Knee (medial, mean of both), Elbow (lateral, mean of both) | degC | `skin_temp_<region>_<stat>_c`; the six `_ave` are the primary features |
| `Civil_status`, `Employment_status`, `Studies`, `Tobacco`, `Alcohol` | sociodemographic / habit codes; labels are not in the file but are recovered from paper Table 1 (codes 1..k = the table's category order; the per-group counts of every code reproduce every Table 1 cell exactly) | code 1-4 (Alcohol 1-5) | carried as `*_code` plus decoded `civil_status`, `employment_status`, `education`, `tobacco`, `alcohol`; not analysed (not in the pre-specified plan) |

Sex is not a column: the paper states all participants are women, all post-menopausal (the authors excluded
pre-menopausal volunteers). Stored as `sex = "female"` with `sex_basis` recording that it comes from the paper.

Protocol (paper Methods): {PROTOCOL} Camera: {DEVICE}.

## Missingness

Measured on the ingested file. Every variable not listed has 0 missing of {len(df)}.

| variable | missing (all) | missing FM | missing controls |
|---|---|---|---|
{miss_md}

One FM participant ({', '.join(part.loc[part.n_temperature_values_missing > 0, 'native_id'])}) lacks all knee and elbow values (the paper
states the knees and elbows of one FM woman "could not be analyzed").

## Linkage strategy

One row per person; features and labels share the native ID in the same sheet. Nothing joins this dataset to any
other person-level dataset (`participant_id` is namespaced; there is no cross-dataset key by design). The dataset
has no geography and must not be joined to geographic or facility rows.

## Limitations and caveats

* **Published conclusion is negative:** "The infra-red thermography is not an effective supplementary assessment
  tool in women with fibromyalgia" (abstract).
* Single site, single research group, post-menopausal women only, age 38-70; no men, no pre-menopausal women.
* Case-control by recruitment channel (FM associations and specialised units vs advertisements/snowball); spectrum
  and selection effects apply to any discrimination estimate.
* Controls are "without FM symptoms"; no other chronic-pain or invisible-illness comparator, so the design cannot
  say whether a temperature difference is specific to FM.
* Room held at 24 degC by protocol, but per-session room temperature, date and season are not in the file.
* Only ROI summary statistics are released (no images), so no other thermal feature can be derived.
* The sociodemographic codes are not labelled in the file (decoded from paper Table 1, see above). Education
  differs sharply by group (university: 12/86 FM vs 36/92 controls, paper chi-square p < 0.05), consistent with
  the different recruitment channels; it is not adjusted for in the pre-specified analysis.
* **Paper text error (verified):** the Results text gives BMI as "CG 27.83 (4.75) ... FMG 25.97 (4.00)". In the
  file the FM group has mean BMI 27.83 and controls 25.97. The file's group coding is independently confirmed by
  Table 1 (all 23 category counts) and Table 2 (FM-column temperature means), so the paper text has the two BMI
  values swapped: FM women have the **higher** BMI.

## Processed outputs

| table | rows |
|---|---|
| `participants__fm_thermography` | {len(part)} |
| `participant_device_features__fm_thermography` | {len(feats)} |

The analysis (`measure_it.wearables.fm_thermography_analysis`) writes `phenotype_signatures__fm_thermography` and
`results/tables/fm_thermography_*.csv`; plan in `docs/ANALYSIS_PLAN_FM_THERMOGRAPHY.md`, results in
`results/FM_THERMOGRAPHY_RESULTS.md`.

## Reproduce

`uv run python -m measure_it.ingestion.fm_thermography`
"""
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text(text)


def registry(df: pd.DataFrame, part: pd.DataFrame, manifest: dict) -> None:
    label = derive_label(df)
    cc = part[part.complete_ave_temperatures]
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "Infrared thermography in fibromyalgia vs controls (PLOS ONE 2021, 0253281) S1 Data",
        "publisher": "PLOS ONE / figshare; University of Valencia (Sempere-Rubio et al.)",
        "landing_url": LANDING_URL,
        "access_urls": [DATA_URL, ARTICLE_URL],
        "license": "CC BY 4.0 (PLOS article licence)",
        "access_conditions": "Open download, no registration or DUA (verified 2026-09-24).",
        "retrieved_at": manifest["files"][DATA_FILE]["retrieved_at"],
        "source_version": f"figshare file 28433660 ({DATA_FILE}), sha256 {manifest['files'][DATA_FILE]['sha256']}",
        "update_date": "static supplement published 2021-06-24",
        "data_layer": "person",
        "unit_of_observation": "participant (one resting thermography session)",
        "sample_size": {"participants": int(len(df)), "fibromyalgia": int((label == "fibromyalgia").sum()),
                        "controls": int((label == "control").sum()),
                        "complete_six_average_temperatures": {"fibromyalgia": int(cc.fibromyalgia.sum()),
                                                              "controls": int((1 - cc.fibromyalgia).sum())}},
        "geographic_resolution": "none",
        "person_level": True,
        "geographic": False,
        "omics": False,
        "wearable": False,
        "participant_linkage": "labels and temperatures share the native ID in one sheet; no cross-dataset linkage",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": ["participants__fm_thermography", "participant_device_features__fm_thermography",
                              "phenotype_signatures__fm_thermography", "results/tables/fm_thermography_*.csv"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.fm_thermography",
        "limitations": [
            "Authors' published conclusion is negative (thermography not an effective supplementary assessment)",
            "Single site; post-menopausal women only; case-control by recruitment channel",
            "No non-FM pain or other invisible-illness comparator",
            "Per-session room temperature, date and season not recorded in the file (room held at 24 degC by protocol)",
            "ROI summary statistics only; images not released",
        ],
        "conditions": [CONDITION_ID],
        "measurement_class_ids": [MEASUREMENT_CLASS_ID],
        "measurement_class_note": MEASUREMENT_CLASS_NOTE,
        "citation": CITATION,
        "notes": "Discovered and ranked in docs/DEVICE_DATASET_DISCOVERY.md (rank 2); primary analysis locked there.",
    })


def run(skip_download: bool = False) -> dict:
    if not skip_download:
        fetch()
    manifest = load_manifest(SOURCE_ID)
    df = read_raw()
    facts = paper_facts(article_text())
    facts["table1_counts_reproduced"] = all(verify_table1(df).values())
    part, feats = build_frames(df, manifest["files"][DATA_FILE]["retrieved_at"])
    write_table(part, f"participants__{SOURCE_ID}", producer="measure_it.ingestion.fm_thermography",
                description="FM thermography participants (label, age, BMI, coded covariates)")
    write_table(feats, f"participant_device_features__{SOURCE_ID}", producer="measure_it.ingestion.fm_thermography",
                description="FM thermography: 18 regional resting skin temperatures (degC) per participant")
    write_audit(df, part, feats, manifest)
    registry(df, part, manifest)
    out = {"participants": len(part), "fibromyalgia": int(part.fibromyalgia.sum()),
           "controls": int((1 - part.fibromyalgia).sum()),
           "complete_ave": int(part.complete_ave_temperatures.sum()), "paper_facts": facts}
    print(out)
    return out


if __name__ == "__main__":
    run()
