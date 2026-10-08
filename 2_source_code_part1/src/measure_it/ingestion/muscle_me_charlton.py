"""MUSCLE-ME Source Data (Charlton et al., Nat Commun 2026) -> person-level steps + CPET + muscle tables.

Source: the article's Source Data workbook (41467_2026_75725_MOESM4_ESM.xlsx, CC BY 4.0), sheet ``Main``. Only
rows ``Group == "MUSCLE-ME"`` are ingested: 25 Long COVID (``Session`` LC), 26 ME/CFS (ME) and 30 age- and
sex-matched healthy controls (CON), one row per person. The 24 AGBRESA bed-rest volunteers (48 rows, two
sessions each) are a separate pre-pandemic cohort without step data and are counted in the audit but not ingested.

Outputs (data/processed):
* ``participants__charlton_lc_mecfs_cpet_source``  one row per person: label, sex, symptom duration, availability flags
* ``participant_wearable_features__charlton_lc_mecfs_cpet_source``  daily step count (hip ActiGraph, per Appelman 2024)
* ``participant_exam_features__charlton_lc_mecfs_cpet_source``  CPET, muscle histology/respirometry, blood lactate (wide)
* ``results/tables/charlton_lc_mecfs_cpet_source_variable_dictionary.csv``  column -> layer -> measurement class
* ``data/raw/charlton_lc_mecfs_cpet_source/{DATA_AUDIT.md, registry_entry.yaml}``

Run: uv run python -m measure_it.ingestion.muscle_me_charlton
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from measure_it.config import TABLES, UNKNOWN, raw_dir
from measure_it.download import download_file, load_manifest
from measure_it.ontology.normalize import normalize_condition
from measure_it.provenance import add_provenance
from measure_it.registry import write_registry_entry
from measure_it.store import write_table

SOURCE_ID = "charlton_lc_mecfs_cpet_source"
PRODUCER = "measure_it.ingestion.muscle_me_charlton"
URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-75725-y/MediaObjects/"
       "41467_2026_75725_MOESM4_ESM.xlsx")
FILENAME = "41467_2026_75725_MOESM4_ESM.xlsx"
DOI = "10.1038/s41467-026-75725-y"
LANDING = f"https://doi.org/{DOI}"
CITATION = ("Charlton BT, Slaghekke A, Appelman B, et al. Skeletal muscle properties in long COVID and ME/CFS differ "
            "from those induced by bed rest. Nat Commun 17 (2026). doi:10.1038/s41467-026-75725-y "
            "(PMID 42649155, PMC13518852; first published 2026-07-28)")
SOURCE_NAME = "MUSCLE-ME Source Data (Charlton et al., Nat Commun 2026), sheet Main"
SOURCE_VERSION = ("Source Data file 41467_2026_75725_MOESM4_ESM.xlsx (HTTP Last-Modified 2026-07-28); "
                  "rows Group == 'MUSCLE-ME'")
SHEET = "Main"
GROUP = "MUSCLE-ME"

SESSIONS = {
    "LC": {"cohort_group": "long_covid", "condition_query": "Long COVID", "is_patient": 1,
           "label_basis": "Long COVID >= 6 months after PCR/serology-proven SARS-CoV-2, PEM by DSQ-PEM, met the Canadian "
                          "Consensus Criteria for ME/CFS; assessed by two clinicians (Amsterdam UMC, NCT05225688)"},
    "ME": {"cohort_group": "me_cfs", "condition_query": "ME/CFS", "is_patient": 1,
           "label_basis": "ME/CFS diagnosed before 2020 (Canadian Consensus Criteria, PEM by DSQ-PEM, seen by an ME/CFS "
                          "specialist or post-COVID physician); all had a previous confirmed SARS-CoV-2 infection"},
    "CON": {"cohort_group": "healthy_control", "condition_query": None, "is_patient": 0,
            "label_basis": "age- and sex-matched healthy controls; per the paper none had residual symptoms after "
                           "SARS-CoV-2 infection and none were hospitalised within 6 months"},
}

STEPS_DEVICE = ("ActiGraph wGT3X-BT accelerometer, right hip, worn from early morning to bedtime during the study period; "
                "daily physical activity = steps per day (method stated in Appelman et al. 2024, same trial NCT05225688; "
                "Charlton 2026 reports 'Daily Steps (steps*day-1)' in Table 1 without restating the method)")

# Column -> (layer, measurement class id in configs/measurements.yaml or None, label)
CPET = "cpet"
VARIABLES: dict[str, tuple[str, str | None, str]] = {
    "Steps": ("wearable_accelerometry", "accelerometry", "daily step count (steps/day)"),
    "VO2_abs": ("cpet", CPET, "peak VO2, absolute"),
    "VO2_rel": ("cpet", CPET, "peak VO2 relative to body mass"),
    "VO2_pred": ("cpet", CPET, "predicted peak VO2"),
    "VO2_perc_pred": ("cpet", CPET, "peak VO2 % predicted"),
    "VO2_perc_pred_wasserman": ("cpet", CPET, "peak VO2 % predicted (Wasserman)"),
    "GET": ("cpet", CPET, "gas-exchange threshold"),
    "GET_rel": ("cpet", CPET, "gas-exchange threshold relative to body mass"),
    "GET_rel_perc_pred": ("cpet", CPET, "gas-exchange threshold, % predicted"),
    "VE_max": ("cpet", CPET, "maximal minute ventilation"),
    "EqCO2_max": ("cpet", CPET, "ventilatory equivalent for CO2 at max"),
    "Peak_power": ("cpet", CPET, "peak power output"),
    "GET_perc": ("cpet", CPET, "GET as % of peak VO2"),
    "EqO2_max": ("cpet", CPET, "ventilatory equivalent for O2 at max"),
    "VE_VCO2_slope": ("cpet", CPET, "VE/VCO2 slope"),
    "O2_pulse": ("cpet", CPET, "O2 pulse (VO2/HR)"),
    "HR_max": ("cpet", CPET, "maximal heart rate"),
    "HR_rest": ("cpet", CPET, "resting heart rate on the ergometer"),
    "AHRR": ("cpet", CPET, "adjusted heart-rate reserve"),
    "VO2_HR_slope": ("cpet", CPET, "VO2/HR slope"),
    "VE_VCO2_slope_submax": ("cpet", CPET, "submaximal VE/VCO2 slope"),
    "Rest_lactate": ("blood_lactate", "blood_biomarkers", "blood lactate at rest"),
    "Baseline_lactate": ("blood_lactate", "blood_biomarkers", "blood lactate at baseline cycling"),
    "Task_Failure_lactate": ("blood_lactate", "blood_biomarkers", "blood lactate at task failure"),
    "SDH": ("muscle_mitochondria", None, "succinate dehydrogenase activity (histochemistry)"),
    "Oxphos": ("muscle_mitochondria", None, "OXPHOS capacity (permeabilised-fibre respirometry)"),
    "Leak": ("muscle_mitochondria", None, "LEAK respiration"),
    "ADP": ("muscle_mitochondria", None, "ADP-stimulated respiration"),
    "N_linked": ("muscle_mitochondria", None, "N-linked respiration"),
    "Uncoupled": ("muscle_mitochondria", None, "uncoupled (ET) respiration"),
    "S_linked": ("muscle_mitochondria", None, "S-linked respiration"),
    "membrane_intact": ("muscle_mitochondria", None, "mitochondrial outer-membrane integrity (cytochrome c test)"),
}
for _c, _lab in [("FCSA", "mean fibre cross-sectional area"), ("Percent_I", "% type I fibres"),
                 ("Percent_I_IIa", "% hybrid I/IIa fibres"), ("Percent_IIa", "% type IIa fibres"),
                 ("Percent_IIa_IIx", "% hybrid IIa/IIx fibres"), ("Percent_IIx", "% type IIx fibres"),
                 ("Total_I_Positive", "count of type I-positive fibres"), ("Total_I_Negative", "count of type I-negative fibres"),
                 ("FT_Percent_I", "fibre-type area % type I"), ("FT_Percent_I_IIa", "fibre-type area % I/IIa"),
                 ("FT_Percent_IIa", "fibre-type area % IIa"), ("FT_Percent_IIa_IIx", "fibre-type area % IIa/IIx"),
                 ("FT_Percent_IIx", "fibre-type area % IIx"), ("TypeI_FCSA", "type I fibre CSA"),
                 ("TypeI_IIa_FCSA", "type I/IIa fibre CSA"), ("TypeIIa_FCSA", "type IIa fibre CSA"),
                 ("Type_IIa_IIx_FCSA", "type IIa/IIx fibre CSA"), ("TypeIIx_FCSA", "type IIx fibre CSA"),
                 ("Myoglobin", "intramyocyte myoglobin content")]:
    VARIABLES[_c] = ("muscle_histology", None, _lab)
for _c, _lab in [("CD", "capillary density"), ("CF", "capillary-to-fibre ratio"), ("Cap_size", "capillary size"),
                 ("TypeI_CF", "type I capillary-to-fibre ratio"), ("TypeII_CF", "type II capillary-to-fibre ratio")]:
    VARIABLES[_c] = ("muscle_capillarisation", None, _lab)
NON_FEATURE = ["Subject", "Session", "Group", "Sex", "Sx_duration"]
EXAM_COLUMNS = [c for c in VARIABLES if c != "Steps"]


def load_raw(path=None) -> pd.DataFrame:
    path = path or download_file(URL, SOURCE_ID, FILENAME)
    return pd.read_excel(path, sheet_name=SHEET)


def select_muscle_me(raw: pd.DataFrame) -> pd.DataFrame:
    m = raw[raw["Group"] == GROUP].copy()
    if m["Subject"].duplicated().any():
        raise ValueError("MUSCLE-ME Subject appears in more than one row")
    unknown = set(m["Session"]) - set(SESSIONS)
    if unknown:
        raise ValueError(f"unexpected Session values {unknown}")
    missing_cols = set(VARIABLES) - set(m.columns)
    if missing_cols:
        raise ValueError(f"Main sheet lacks {missing_cols}")
    return m.reset_index(drop=True)


def participant_id(native) -> str:
    return f"{SOURCE_ID}:{native}"


def condition_ids() -> dict[str, str | None]:
    out = {}
    for s, spec in SESSIONS.items():
        if spec["condition_query"] is None:
            out[s] = None
            continue
        r = normalize_condition(spec["condition_query"])
        if r["status"] != "matched" or not r["matches"]:
            raise ValueError(f"{spec['condition_query']!r} did not normalize: {r['status']}")
        out[s] = r["matches"][0]["canonical_condition_id"]
    return out


def build_participants(m: pd.DataFrame, retrieved_at: str) -> pd.DataFrame:
    cids = condition_ids()
    df = pd.DataFrame({
        "participant_id": m["Subject"].map(participant_id),
        "dataset_id": SOURCE_ID,
        "native_id": m["Subject"].astype(str),
        "session_code": m["Session"],
        "cohort_group": m["Session"].map(lambda s: SESSIONS[s]["cohort_group"]),
        "condition_id": m["Session"].map(cids),
        "is_patient": m["Session"].map(lambda s: SESSIONS[s]["is_patient"]).astype(int),
        "label_basis": m["Session"].map(lambda s: SESSIONS[s]["label_basis"]),
        "sex": m["Sex"].str.lower(),
        "age": UNKNOWN,  # not released in the Source Data (paper Table 1 gives group medians only)
        "symptom_duration_days": m["Sx_duration"].astype(float),
        "has_steps": m["Steps"].notna(),
        "has_vo2_rel": m["VO2_rel"].notna(),
        "n_exam_variables_nonmissing": m[EXAM_COLUMNS].notna().sum(axis=1).astype(int),
        "device": "ActiGraph wGT3X-BT (hip) for steps; electronically braked cycle ergometer (Lode) ramp CPET",
    })
    return add_provenance(
        df, data_layer="person", source_name=SOURCE_NAME, source_version=SOURCE_VERSION, retrieved_at=retrieved_at,
        evidence_type="person_exam_measurement", source_record_id="participant_id",
        provenance_notes=("Roster of MUSCLE-ME rows in the Source Data; label = Session (clinical case definition, not self-report). "
                          "Age not released. Probably overlaps Appelman et al. 2024 (paper: data from 21 controls and 25 Long "
                          "COVID patients previously published). No geography."))


def build_wearable(m: pd.DataFrame, retrieved_at: str) -> pd.DataFrame:
    df = pd.DataFrame({
        "participant_id": m["Subject"].map(participant_id),
        "dataset_id": SOURCE_ID,
        "native_id": m["Subject"].astype(str),
        "device_class": "actigraph_hip_daily_steps",
        "steps_mean_daily": m["Steps"].astype(float),
        "log10_steps_mean_daily": np.log10(m["Steps"].astype(float).where(m["Steps"] > 0)),
    })
    df = df[df["steps_mean_daily"].notna()].reset_index(drop=True)
    return add_provenance(
        df, data_layer="person", source_name=SOURCE_NAME, source_version=SOURCE_VERSION, retrieved_at=retrieved_at,
        evidence_type="person_device_measurement", source_record_id="participant_id",
        provenance_notes=f"steps_mean_daily = Source Data column 'Steps' as released. {STEPS_DEVICE}.")


def build_exam(m: pd.DataFrame, retrieved_at: str) -> pd.DataFrame:
    df = pd.concat([pd.DataFrame({"participant_id": m["Subject"].map(participant_id), "dataset_id": SOURCE_ID,
                                  "native_id": m["Subject"].astype(str)}),
                    m[EXAM_COLUMNS].astype(float).reset_index(drop=True)], axis=1)
    return add_provenance(
        df, data_layer="person", source_name=SOURCE_NAME, source_version=SOURCE_VERSION, retrieved_at=retrieved_at,
        evidence_type="person_exam_measurement", source_record_id="participant_id",
        provenance_notes=("Wide table; columns keep the Source Data names. CPET = ramp test on an electronically braked cycle "
                          "ergometer; muscle = vastus lateralis biopsy (histochemistry, immunohistochemistry, permeabilised-"
                          "fibre respirometry); lactate = blood lactate around the CPET. Units as released (not stated in the "
                          "file). Muscle variables are physiology/histology, not omics."))


def variable_dictionary(m: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col, (layer, mclass, label) in VARIABLES.items():
        row = {"column": col, "label": label, "layer": layer,
               "measurement_class": mclass if mclass else UNKNOWN,
               "measurement_class_note": "" if mclass else "no class in configs/measurements.yaml for muscle biopsy histology/respirometry",
               "n_nonmissing_total": int(m[col].notna().sum())}
        for s in SESSIONS:
            row[f"n_nonmissing_{s}"] = int(m.loc[m.Session == s, col].notna().sum())
        rows.append(row)
    return pd.DataFrame(rows)


def write_audit(raw: pd.DataFrame, m: pd.DataFrame, parts: pd.DataFrame, wear: pd.DataFrame, exam: pd.DataFrame,
                vd: pd.DataFrame, man: dict) -> None:
    f = man["files"][FILENAME]
    n = m.Session.value_counts().to_dict()
    sexct = pd.crosstab(m.Session, m.Sex)
    L: list[str] = []
    A = L.append
    A("# DATA AUDIT — MUSCLE-ME Source Data (Charlton et al., Nat Commun 2026)\n")
    A("| Field | Value |\n|---|---|")
    A(f"| source_id | {SOURCE_ID} |")
    A(f"| Source (dataset/API name, exact files/endpoints) | Source Data workbook `{FILENAME}` ({URL}), sheet `Main` |")
    A("| Publishing organization | Springer Nature (Nature Communications) on behalf of the authors: Amsterdam UMC / Vrije Universiteit Amsterdam (MUSCLE-ME, NCT05225688) and DLR/ESA/NASA (AGBRESA bed rest, DRKS00015677) |")
    A(f"| Retrieval date (UTC) | {f['retrieved_at']} |")
    A(f"| Source version / release | {SOURCE_VERSION}; sha256 {f['sha256']} |")
    A("| Source update date / cadence | article first published 2026-07-28; static (no updates expected) |")
    A("| License / access conditions | CC BY 4.0 (article licence covers the supplementary Source Data); anonymous download, no registration or DUA. Paper: processed data also on GitHub bcn95/LC_ME_CFS_BED_REST (Zenodo 10.5281/zenodo.20642303); raw data 'freely available from the corresponding author' (not requested). |")
    A("| Unit of observation | person (one row per MUSCLE-ME participant; cross-sectional) |")
    A(f"| Sample size (actual, as ingested) | {len(m)} MUSCLE-ME people: Long COVID {n.get('LC', 0)}, ME/CFS {n.get('ME', 0)}, healthy controls {n.get('CON', 0)}; with steps {int(parts.has_steps.sum())} ({int(parts[parts.is_patient == 1].has_steps.sum())} patients, {int(parts[parts.is_patient == 0].has_steps.sum())} controls); with VO2_rel {int(parts.has_vo2_rel.sum())}. Not ingested: {int((raw.Group != GROUP).sum())} AGBRESA bed-rest rows ({raw.loc[raw.Group != GROUP, 'Subject'].nunique()} people x 2 sessions). |")
    A("| Geography (resolution, vintage) | none (Amsterdam, the Netherlands, study site only; no participant geography) |")
    A("| Person-level? | yes |")
    A("| Geographic? | no |")
    A("| Omics? | no (muscle histology and respirometry are physiology, not omics) |")
    A("| Wearable? | yes (hip-worn ActiGraph daily step count; research-grade accelerometer, not a consumer wearable) |")
    A("| True participant linkage across modalities? | yes, within this file: steps, CPET, muscle biopsy and lactate share the `Subject` id on one row. No linkage to Appelman 2024 ids (different numbering) or to any other project dataset. |")
    A("\n## Files / endpoints retrieved\n")
    A(f"* `{FILENAME}` — {f['bytes']} bytes, sha256 `{f['sha256']}`, from `{URL}` (see MANIFEST.json).")
    A("* Paper full text read via Europe PMC (PMC13518852) for definitions; the step-measurement method is from Appelman et al. 2024 (PMC10766651), same trial.")
    A("* Other sheets in the workbook (`Lactate`, `FT_CF`, `DSQ_PEM`, `ME_pathogen_infection`, `Sx_percentages`, `SuppFig7`, `Fig2a`) are not ingested; `DSQ_PEM` (PEM questionnaire, patients only) could support a severity analysis later.")
    A("\n## Key variables\n")
    A("* `Subject` -> native id; `participant_id = charlton_lc_mecfs_cpet_source:<Subject>`.")
    A("* `Session` -> label: LC = Long COVID (Canadian Consensus Criteria), ME = ME/CFS (diagnosed before 2020), CON = healthy control. Normalised with `measure_it.ontology.normalize`: LC -> `long_covid`, ME -> `me_cfs`.")
    A("* `Sex` -> sex. **Age is not released** (paper Table 1 gives group medians only), so no age adjustment is possible.")
    A("* `Sx_duration` -> symptom duration in days (patients only; 45 of 51 non-missing).")
    A(f"* `Steps` -> `steps_mean_daily` (steps/day). {STEPS_DEVICE}. The averaging window and wear-time rule are not stated in the Source Data.")
    A("* 20 CPET columns (measurement class `cpet`), 3 blood lactate columns (`blood_biomarkers`), 32 muscle biopsy columns (histology, capillarisation, mitochondrial respiration; no class in `configs/measurements.yaml`). Full list with layer and class: `results/tables/charlton_lc_mecfs_cpet_source_variable_dictionary.csv`. Units are not stated in the file.")
    A("\n## Missingness\n")
    A(f"Sex by session: {json.dumps({s: sexct.loc[s].to_dict() for s in sexct.index})}.\n")
    A("| column | layer | non-missing LC | non-missing ME | non-missing CON | total (of 81) |")
    A("|---|---|---|---|---|---|")
    for _, r in vd.iterrows():
        A(f"| {r.column} | {r.layer} | {r.n_nonmissing_LC}/{n.get('LC', 0)} | {r.n_nonmissing_ME}/{n.get('ME', 0)} | {r.n_nonmissing_CON}/{n.get('CON', 0)} | {r.n_nonmissing_total} |")
    A("\nThe 8 controls without steps are the main gap: steps are missing only among controls, so a complete-case analysis could be biased if missingness relates to activity. The analysis reports worst/best-case bounds.")
    A("\n## Linkage strategy\n")
    A("All modalities for a person are on the same `Main` row keyed by `Subject`. There is no cross-dataset key: the same people probably appear in Appelman et al. 2024 (the paper says data from 21 controls and 25 Long COVID patients were previously published) under a different numbering, so the two sources are not independent evidence and must not be pooled. No join to geographic or facility layers.")
    A("\n## Limitations and caveats\n")
    A("* Healthy controls only: accuracy against healthy people overstates accuracy for the clinical question (patients vs other fatigued or deconditioned people).")
    A("* Incorporation bias: reduced activity / PEM is part of both case definitions.")
    A("* Patients were well enough to complete a maximal CPET and a muscle biopsy (the paper describes them as mildly affected).")
    A("* Controls: the paper says none had residual symptoms after SARS-CoV-2 infection, i.e. controls are post-COVID recovered people, which is a stronger comparison for Long COVID than never-infected controls.")
    A("* ME/CFS patients were enrolled about 20 months after the Long COVID patients (October 2023 vs February 2022 enrolment reference dates in Table 1).")
    A("* Step-count window, wear-time validity rule and whether the window overlapped the post-CPET PEM period are not stated for this release.")
    A("* Small: 51 vs 22 with steps.")
    A("\n## Processed outputs\n")
    A(f"* `participants__{SOURCE_ID}` — {len(parts)} rows")
    A(f"* `participant_wearable_features__{SOURCE_ID}` — {len(wear)} rows (people with steps)")
    A(f"* `participant_exam_features__{SOURCE_ID}` — {len(exam)} rows x {len(EXAM_COLUMNS)} exam columns")
    A(f"* `results/tables/{SOURCE_ID}_variable_dictionary.csv` — {len(vd)} rows")
    A(f"* analysis: `phenotype_signatures__{SOURCE_ID}` and `results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md` (module `measure_it.wearables.muscle_me_steps_analysis`; plan `docs/ANALYSIS_PLAN_CHARLTON_LC_MECFS_CPET_SOURCE.md`)")
    A("\n## Reproduce\n")
    A("`uv run python -m measure_it.ingestion.muscle_me_charlton` then `uv run python -m measure_it.wearables.muscle_me_steps_analysis`")
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text("\n".join(L) + "\n")


def write_registry(parts: pd.DataFrame, raw: pd.DataFrame, man: dict) -> None:
    f = man["files"][FILENAME]
    n = parts.session_code.value_counts().to_dict()
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "MUSCLE-ME: Long COVID and ME/CFS vs matched healthy controls, daily steps + CPET + muscle biopsy (Charlton et al., Nat Commun 2026, Source Data)",
        "publisher": "Springer Nature / Nature Communications (authors: Amsterdam UMC, Vrije Universiteit Amsterdam)",
        "landing_url": LANDING,
        "access_urls": [URL],
        "license": "CC BY 4.0 (article and supplementary Source Data)",
        "access_conditions": "Open anonymous download; no registration or DUA (verified 2026-09-24).",
        "retrieved_at": f["retrieved_at"],
        "source_version": f"{SOURCE_VERSION}; sha256 {f['sha256']}",
        "update_date": "2026-07-28 (article publication; static)",
        "data_layer": "person",
        "unit_of_observation": "person (cross-sectional)",
        "sample_size": {"people_ingested": int(len(parts)), "long_covid": int(n.get("LC", 0)), "me_cfs": int(n.get("ME", 0)),
                        "healthy_controls": int(n.get("CON", 0)), "with_steps": int(parts.has_steps.sum()),
                        "controls_with_steps": int(parts[parts.is_patient == 0].has_steps.sum()),
                        "with_vo2_rel": int(parts.has_vo2_rel.sum()),
                        "agbresa_rows_not_ingested": int((raw.Group != GROUP).sum())},
        "geographic_resolution": "none",
        "person_level": True, "geographic": False, "omics": False, "wearable": True,
        "participant_linkage": "steps, CPET, muscle biopsy and lactate on one row per Subject (within-file); no key to Appelman 2024 or other datasets",
        "true_participant_linkage_across_modalities": True,
        "status": "ingested",
        "processed_outputs": [f"participants__{SOURCE_ID}", f"participant_wearable_features__{SOURCE_ID}",
                              f"participant_exam_features__{SOURCE_ID}", f"phenotype_signatures__{SOURCE_ID}",
                              f"results/tables/{SOURCE_ID}_variable_dictionary.csv"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": PRODUCER,
        "limitations": [
            "healthy controls only (no disease or deconditioned comparator); incorporation bias (activity limitation is part of the case definitions)",
            "8 of 30 controls lack steps; age not released",
            "probably the same people as Appelman et al. 2024 (same trial NCT05225688): not independent evidence",
            "patients able to complete maximal CPET and biopsy (milder end)",
        ],
        "notes": f"{CITATION}. Steps device: {STEPS_DEVICE}.",
    })


def run(skip_download: bool = False) -> dict:
    path = raw_dir(SOURCE_ID) / FILENAME
    raw = load_raw(path if skip_download and path.exists() else None)
    man = load_manifest(SOURCE_ID)
    retrieved_at = man["files"][FILENAME]["retrieved_at"]
    m = select_muscle_me(raw)
    parts = build_participants(m, retrieved_at)
    wear = build_wearable(m, retrieved_at)
    exam = build_exam(m, retrieved_at)
    vd = variable_dictionary(m)
    write_table(parts, f"participants__{SOURCE_ID}", producer=PRODUCER,
                description="MUSCLE-ME roster (Long COVID, ME/CFS, healthy controls), one row per person")
    write_table(wear, f"participant_wearable_features__{SOURCE_ID}", producer=PRODUCER,
                description="Daily step count (hip ActiGraph) per person, as released")
    write_table(exam, f"participant_exam_features__{SOURCE_ID}", producer=PRODUCER,
                description="CPET, muscle biopsy and blood lactate per person (wide, Source Data column names)")
    TABLES.mkdir(parents=True, exist_ok=True)
    vd.to_csv(TABLES / f"{SOURCE_ID}_variable_dictionary.csv", index=False)
    write_audit(raw, m, parts, wear, exam, vd, man)
    write_registry(parts, raw, man)
    return {"participants": len(parts), "with_steps": len(wear), "exam_rows": len(exam),
            "sessions": parts.session_code.value_counts().to_dict()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--skip-download", action="store_true")
    a = ap.parse_args()
    print(json.dumps(run(skip_download=a.skip_download), indent=2, default=str))


if __name__ == "__main__":
    main()
