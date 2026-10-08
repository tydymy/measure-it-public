"""Launch agent part A (measure_it.agent.read / profile / mapper / llm) on SYNTHETIC messy inputs generated here.

No real data and no network: the LLM tests use a mock client that records what would have been sent.
"""
from __future__ import annotations

import json
import random
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml

from measure_it.agent import llm as L
from measure_it.agent import mapper as M
from measure_it.agent import profile as P
from measure_it.agent import read as R

SEED = 20260923
DICT_HEADER = ["Variable / Field Name", "Form Name", "Section Header", "Field Type", "Field Label",
               "Choices, Calculations, OR Slider Labels", "Field Note", "Text Validation Type OR Show Slider Number",
               "Text Validation Min", "Text Validation Max", "Identifier?", "Branching Logic (Show field only if...)",
               "Required Field?", "Custom Alignment", "Question Number (surveys only)", "Matrix Group Name",
               "Matrix Ranking?", "Field Annotation"]


# ------------------------------------------------------------------------------------------------------------------
# synthetic fixtures
# ------------------------------------------------------------------------------------------------------------------
def _dict_row(field, form, ftype, label, choices="", note="", validation="", identifier=""):
    r = dict.fromkeys(DICT_HEADER, "")
    r.update({"Variable / Field Name": field, "Form Name": form, "Field Type": ftype, "Field Label": label,
              "Choices, Calculations, OR Slider Labels": choices, "Field Note": note,
              "Text Validation Type OR Show Slider Number": validation, "Identifier?": identifier})
    return r


COHORTS = [("1", "LC", 40), ("2", "CLYME", 25), ("3", "ALYME", 8), ("4", "HC", 30)]   # ALYME n=8 < 11


def make_redcap(root: Path) -> Path:
    """A MAESTRO-like REDCap export: data CSV (with a repeating ring_daily instrument) + data dictionary CSV."""
    rng = np.random.default_rng(SEED)
    root.mkdir(parents=True, exist_ok=True)
    d = [
        _dict_row("record_id", "enrollment", "text", "Record ID"),
        _dict_row("cohort", "enrollment", "radio", "Study cohort",
                  "1, LC | 2, CLYME | 3, ALYME | 4, HC"),
        _dict_row("age", "enrollment", "text", "Age at enrollment (years)", validation="integer"),
        _dict_row("sex", "enrollment", "radio", "Sex assigned at birth", "1, Female | 2, Male"),
        _dict_row("zip", "enrollment", "text", "Home ZIP code", identifier="y"),
        _dict_row("visit_date", "enrollment", "text", "Visit date", validation="date_ymd"),
        _dict_row("notes", "enrollment", "notes", "Coordinator notes"),
        _dict_row("mh_conditions", "medical_history", "checkbox", "Which conditions have you been diagnosed with?",
                  "1, Raynaud's phenomenon | 2, Migraine | 3, POTS"),
        _dict_row("raynaud_hx", "medical_history", "yesno",
                  "Have you ever been diagnosed with Raynaud's phenomenon?"),
        _dict_row("compass31_q1", "compass31", "radio", "In the past year, have you ever felt faint?",
                  "0, No | 1, Yes"),
        _dict_row("compass31_vaso", "compass31", "text", "COMPASS-31 vasomotor score", validation="number"),
        _dict_row("compass31_total", "compass31", "text", "COMPASS-31 total weighted score", validation="number"),
        _dict_row("GlycA", "nightingale_nmr", "text", "Glycoprotein acetyls", note="mmol/L", validation="number"),
        _dict_row("Total_C", "nightingale_nmr", "text", "Total cholesterol", note="mmol/L", validation="number"),
        _dict_row("Glucose", "nightingale_nmr", "text", "Glucose", note="mmol/L", validation="number"),
        _dict_row("Creatinine", "nightingale_nmr", "text", "Creatinine", note="umol/L", validation="number"),
        _dict_row("XXL_VLDL_P", "nightingale_nmr", "text", "Concentration of chylomicrons and extremely large VLDL",
                  validation="number"),
        _dict_row("Ala", "nightingale_nmr", "text", "Alanine", note="mmol/L", validation="number"),
        _dict_row("olink_il6", "olink", "text", "IL-6 (NPX)", validation="number"),
        _dict_row("olink_tnf", "olink", "text", "TNF (NPX)", validation="number"),
        _dict_row("cap_density", "capillaroscopy", "text", "Capillary density (capillaries/mm)", validation="number"),
        _dict_row("cap_tortuosity_pct", "capillaroscopy", "text", "Tortuous capillaries (%)", validation="number"),
        _dict_row("cap_hemorrhages", "capillaroscopy", "text", "Number of microhemorrhages", validation="integer"),
        _dict_row("hba1c_pct", "labs_misc", "text", "HbA1c", validation="number"),
        _dict_row("mystery_17", "labs_misc", "text", "Score X", validation="number"),
        _dict_row("ring_date", "ring_daily", "text", "Day of ring summary", validation="date_ymd"),
        _dict_row("ring_rhr", "ring_daily", "text", "Resting heart rate (bpm)", validation="number"),
        _dict_row("ring_hrv_rmssd", "ring_daily", "text", "Nightly HRV RMSSD (ms)", validation="number"),
        _dict_row("ring_spo2", "ring_daily", "text", "Average SpO2 (%)", validation="number"),
        _dict_row("ring_sleep_hours", "ring_daily", "text", "Total sleep (hours)", validation="number"),
        _dict_row("ring_steps", "ring_daily", "text", "Steps", validation="integer"),
    ]
    pd.DataFrame(d, columns=DICT_HEADER).to_csv(root / "MAESTRO_DataDictionary_2026-10-07.csv", index=False)

    rows, i = [], 0
    for code, _, n in COHORTS:
        for _ in range(n):
            i += 1
            rid = f"MAE-{i:04d}"
            mh = rng.integers(0, 2, 3)
            rows.append({
                "record_id": rid, "redcap_event_name": "baseline_arm_1", "redcap_repeat_instrument": "",
                "redcap_repeat_instance": "", "cohort": code, "age": int(rng.integers(19, 75)),
                "sex": int(rng.integers(1, 3)), "zip": f"0{int(rng.integers(1000, 2999)):04d}",
                "visit_date": f"2025-{int(rng.integers(1, 13)):02d}-{int(rng.integers(1, 28)):02d}",
                "notes": rng.choice(["Participant reported feeling dizzy after the stand test today",
                                     "No issues; blood draw completed without complications at 9am",
                                     f"Called back on {int(rng.integers(1, 28))} March about scheduling"]),
                "mh_conditions___1": int(mh[0]), "mh_conditions___2": int(mh[1]), "mh_conditions___3": int(mh[2]),
                "raynaud_hx": int(rng.integers(0, 2)), "compass31_q1": int(rng.integers(0, 2)),
                "compass31_vaso": round(float(rng.uniform(0, 5)), 2),
                "compass31_total": round(float(rng.uniform(0, 80)), 1),
                "GlycA": round(float(rng.normal(0.8, 0.1)), 3), "Total_C": round(float(rng.normal(4.8, 0.8)), 2),
                "Glucose": round(float(rng.normal(5.2, 0.6)), 2), "Creatinine": round(float(rng.normal(75, 12)), 1),
                "XXL_VLDL_P": float(rng.normal(1e-4, 2e-5)), "Ala": round(float(rng.normal(0.4, 0.05)), 3),
                "olink_il6": round(float(rng.normal(3, 1)), 3), "olink_tnf": round(float(rng.normal(2, 0.5)), 3),
                "cap_density": round(float(rng.normal(8, 1.5)), 1),
                "cap_tortuosity_pct": round(float(rng.uniform(0, 40)), 1),
                "cap_hemorrhages": int(rng.integers(0, 4)), "hba1c_pct": round(float(rng.normal(5.4, 0.4)), 1),
                "mystery_17": round(float(rng.normal(50, 10)), 2),
                "enrollment_complete": 2, "compass31_complete": 2,
            })
            for day in range(1, 6):
                rows.append({"record_id": rid, "redcap_event_name": "baseline_arm_1",
                             "redcap_repeat_instrument": "ring_daily", "redcap_repeat_instance": day,
                             "ring_date": f"2025-06-{day:02d}", "ring_rhr": int(rng.integers(55, 90)),
                             "ring_hrv_rmssd": round(float(rng.normal(40, 10)), 1),
                             "ring_spo2": round(float(rng.normal(96, 1)), 1),
                             "ring_sleep_hours": round(float(rng.normal(7, 1)), 2),
                             "ring_steps": int(rng.integers(1000, 12000)), "ring_daily_complete": 2})
    pd.DataFrame(rows).to_csv(root / "MAESTRO_DATA_2026-10-07.csv", index=False)
    return root


def make_excel(path: Path) -> Path:
    rng = np.random.default_rng(SEED + 1)
    n = 40
    part = pd.DataFrame({"subject_id": [f"S{i:03d}" for i in range(n)],
                         "group": ["ME/CFS"] * 22 + ["Control"] * 18,
                         "age_band": rng.choice(["30-39", "40-49", "50-59"], n),
                         "gender": rng.choice(["F", "M"], n), "BMI": rng.normal(26, 4, n).round(1)})
    labs = []
    for sid in part["subject_id"]:
        for test, unit, mu in [("Glucose", "mg/dL", 95), ("Total cholesterol", "mg/dL", 190),
                               ("Creatinine", "mg/dL", 0.9), ("C-reactive protein", "mg/L", 2.0),
                               ("Lab panel Q", "U", 5.0)]:
            labs.append({"subject_id": sid, "test_name": test, "result": round(float(rng.normal(mu, mu / 10)), 2),
                         "units": unit, "collection_date": "2024-03-04"})
    with pd.ExcelWriter(path) as xw:
        part.to_excel(xw, sheet_name="Participants", index=False)
        pd.DataFrame(labs).to_excel(xw, sheet_name="Lab Results", index=False)
    return path


def make_xpt(path: Path) -> Path:
    import pyreadstat
    rng = np.random.default_rng(SEED + 2)
    n = 60
    df = pd.DataFrame({"SEQN": np.arange(62161, 62161 + n).astype(float), "RIAGENDR": rng.integers(1, 3, n).astype(float),
                       "RIDAGEYR": rng.integers(20, 80, n).astype(float), "LBXTC": rng.normal(190, 30, n).round(0),
                       "LBXGLU": rng.normal(100, 15, n).round(0)})
    pyreadstat.write_xport(df, str(path), file_format_version=5, table_name="LABDEMO",
                           column_labels=["Respondent sequence number", "Gender", "Age in years at screening",
                                          "Total Cholesterol (mg/dL)", "Fasting Glucose (mg/dL)"])
    return path


def make_fhir(path: Path) -> Path:
    random.seed(SEED)
    entries = []
    for i in range(15):
        pid = f"pat-{i}"
        entries.append({"fullUrl": f"urn:uuid:{pid}", "resource": {
            "resourceType": "Patient", "id": pid, "gender": random.choice(["female", "male"]),
            "birthDate": f"19{random.randint(50, 99)}-0{random.randint(1, 9)}-1{random.randint(0, 9)}",
            "address": [{"postalCode": f"0{random.randint(1000, 2999)}", "state": "MA"}]}})
        entries.append({"resource": {
            "resourceType": "Condition", "subject": {"reference": f"urn:uuid:{pid}"},
            "code": {"coding": [{"system": "http://hl7.org/fhir/sid/icd-10-cm", "code": "U09.9",
                                 "display": "Post COVID-19 condition, unspecified"},
                                {"system": "http://snomed.info/sct", "code": "1119303003"}]},
            "onsetDateTime": "2022-01-15"}})
        entries.append({"resource": {
            "resourceType": "Observation", "subject": {"reference": f"Patient/{pid}"},
            "code": {"coding": [{"system": "http://loinc.org", "code": "2093-3", "display": "Cholesterol"}]},
            "valueQuantity": {"value": round(random.gauss(190, 30)), "unit": "mg/dL"},
            "effectiveDateTime": "2023-05-02T09:30:00Z"}})
        entries.append({"resource": {
            "resourceType": "MedicationStatement", "subject": {"reference": f"Patient/{pid}"}, "status": "active",
            "medicationCodeableConcept": {"coding": [{"system": "http://www.nlm.nih.gov/research/umls/rxnorm",
                                                      "code": "6809", "display": "metformin"}]},
            "effectivePeriod": {"start": "2023-01-01"}}})
    entries.append({"resource": {"resourceType": "Encounter", "id": "enc-1"}})
    path.write_text(json.dumps({"resourceType": "Bundle", "type": "collection", "entry": entries}))
    return path


def make_omop(root: Path) -> Path:
    rng = np.random.default_rng(SEED + 3)
    root.mkdir(parents=True, exist_ok=True)
    n = 30
    pd.DataFrame({"person_id": np.arange(1, n + 1), "gender_concept_id": rng.choice([8507, 8532], n),
                  "year_of_birth": rng.integers(1950, 2000, n), "race_concept_id": 0,
                  "location_id": rng.integers(1, 5, n), "person_source_value": [f"MRN{i:07d}" for i in range(n)]}
                 ).to_csv(root / "person.csv", index=False)
    pd.DataFrame({"condition_occurrence_id": np.arange(1, 2 * n + 1), "person_id": np.repeat(np.arange(1, n + 1), 2),
                  "condition_concept_id": 0, "condition_start_date": "2021-03-01",
                  "condition_source_value": rng.choice(["U09.9", "G93.32", "I73.00", "R53.83"], 2 * n)}
                 ).to_csv(root / "condition_occurrence.csv", index=False)
    pd.DataFrame({"measurement_id": np.arange(1, 2 * n + 1), "person_id": np.repeat(np.arange(1, n + 1), 2),
                  "measurement_concept_id": 0, "measurement_date": "2021-04-01",
                  "measurement_source_value": rng.choice(["2093-3", "2345-7", "2160-0"], 2 * n),
                  "value_as_number": rng.normal(100, 20, 2 * n).round(1), "unit_source_value": "mg/dL"}
                 ).to_csv(root / "measurement.csv", index=False)
    return root


@pytest.fixture(scope="module")
def redcap_dir(tmp_path_factory):
    return make_redcap(tmp_path_factory.mktemp("redcap") / "maestro_export")


@pytest.fixture(scope="module")
def redcap_run(redcap_dir, tmp_path_factory):
    out = tmp_path_factory.mktemp("run_redcap")
    res = M.map_cmd(redcap_dir, out, llm=None, condition_hint="long_covid",
                    measurement_hint="nailfold_capillaroscopy", echo=lambda *a: None)
    return res, out


def _cols(mapping, table):
    return mapping["tables"][table]["columns"]


def _all_levels(profile):
    for t in profile["tables"]:
        for c in t["columns"]:
            for lv in c.get("levels") or []:
                yield t["name"], c["name"], lv


# ------------------------------------------------------------------------------------------------------------------
# REDCap (MAESTRO-like)
# ------------------------------------------------------------------------------------------------------------------
def test_redcap_read_layout(redcap_dir):
    ds = R.read_input(redcap_dir)
    assert ds.format == "redcap" and ds.container == "folder"
    assert ds.dictionary["n_fields"] == 31
    names = [t.name for t in ds.tables]
    assert names == ["records", "ring_daily"]
    rec, ring = ds.table("records"), ds.table("ring_daily")
    assert len(rec.df) == 103 and len(ring.df) == 515
    assert rec.labels["mh_conditions___1"] == "Which conditions have you been diagnosed with? (choice=Raynaud's phenomenon)"
    assert rec.choices["cohort"] == {"1": "LC", "2": "CLYME", "3": "ALYME", "4": "HC"}
    assert "zip" in rec.identifier_fields
    assert "ring_rhr" not in rec.df.columns and "ring_rhr" in ring.df.columns
    assert rec.df["zip"].map(lambda v: isinstance(v, str) and v.startswith("0")).all()   # leading zeros kept


def test_redcap_profile(redcap_run):
    prof = redcap_run[0]["profile"]
    rec = next(t for t in prof["tables"] if t["name"] == "records")
    ring = next(t for t in prof["tables"] if t["name"] == "ring_daily")
    assert rec["id_candidates"][0] == "record_id" and rec["shape"] == "wide"
    assert ring["shape"] == "per_day"
    cols = {c["name"]: c for c in rec["columns"]}
    assert "geography" in cols["zip"]["flags"] or "identifier" in cols["zip"]["flags"]
    assert cols["zip"]["levels"] is None and cols["zip"]["numeric"] is None
    assert cols["visit_date"]["dtype"] == "date" and "date" in cols["visit_date"]["flags"]
    assert "free_text" in cols["notes"]["flags"] and cols["notes"]["levels"] is None
    assert "exact_age" in cols["age"]["flags"] and cols["age"]["numeric"] is None
    assert cols["GlycA"]["unit"] == "mmol/L" and cols["Creatinine"]["unit"] == "umol/L"
    assert cols["cap_density"]["unit"] == "/mm"
    assert cols["olink_il6"]["unit"] == "NPX"
    assert cols["compass31_vaso"]["numeric"]["min"] >= 0
    lv = dict((k, v) for k, v in cols["cohort"]["levels"])
    assert lv["1"] == 40 and lv["4"] == 30 and "3" not in lv           # ALYME (n=8) never shown
    assert ["<11 other", None] in cols["cohort"]["levels"]
    assert json.loads(Path(redcap_run[1] / "profile.json").read_text())["format"] == "redcap"


def test_no_level_below_11_anywhere(redcap_run):
    for _, _, (v, n) in _all_levels(redcap_run[0]["profile"]):
        assert n is None or n >= P.SMALL_CELL


def test_redcap_mapping_key_columns(redcap_run):
    mapping = redcap_run[0]["mapping"]
    rec = _cols(mapping, "records")
    t = mapping["tables"]["records"]
    assert t["participant_id"] == "record_id" and t["day_index_from"] == "visit_date"
    assert t["event_column"] == "redcap_event_name"
    assert rec["record_id"] == {"role": "participant_id"}
    assert rec["zip"]["role"] == "drop" and rec["zip"]["reason"] in ("identifier", "geography")
    assert rec["notes"] == {"role": "drop", "reason": "free_text"}
    assert rec["visit_date"] == {"role": "date"}
    assert rec["age"] == {"role": "age", "transform": "age_band"}
    assert rec["sex"]["role"] == "sex" and rec["sex"]["value_map"] == {"1": "female", "2": "male"}
    lab = rec["cohort"]
    assert lab["role"] == "label" and lab["condition"] == "long_covid"
    assert lab["value_map"] == {"1": 1, "2": None, "3": None, "4": 0}
    assert lab["level_conditions"] == {"1": "long_covid", "2": "ptlds", "3": "lyme_disease", "4": "control"}
    assert rec["mh_conditions___1"]["role"] == "self_report_history" and rec["mh_conditions___1"]["icd10cm"] == "I73.0"
    assert rec["mh_conditions___1"]["value_map"] == {"0": 0, "1": 1}
    assert rec["mh_conditions___3"]["icd10cm"] == "G90.A"
    assert rec["raynaud_hx"]["role"] == "self_report_history" and rec["raynaud_hx"]["icd10cm"] == "I73.0"
    assert rec["compass31_vaso"] == {"role": "survey", "instrument": "COMPASS31", "item": "vaso"}
    assert rec["compass31_q1"]["role"] == "survey"                          # an instrument item, not a history item
    assert rec["GlycA"]["role"] == "lab" and rec["GlycA"]["loinc"] == "82730-3"
    assert rec["GlycA"]["platform"] == "nmr_nightingale" and rec["GlycA"]["unit"] == "mmol/L"
    assert rec["Total_C"]["loinc"] == "2093-3" and rec["Total_C"]["factor"] == 38.67
    assert rec["Glucose"]["loinc"] == "2345-7" and rec["Creatinine"]["loinc"] == "2160-0"
    assert rec["XXL_VLDL_P"] == {"role": "omics", "layer": "metabolomics", "platform": "nmr_nightingale"}
    assert rec["olink_il6"] == {"role": "omics", "layer": "proteomics", "platform": "olink"}
    assert rec["cap_density"]["role"] == "device" and rec["cap_density"]["measurement_class"] == "capillaroscopy"
    assert rec["enrollment_complete"]["role"] == "drop"
    assert rec["mystery_17"]["role"] == "unmapped" and "fuzzy" in rec["mystery_17"]["reason"]
    ring = _cols(mapping, "ring_daily")
    assert mapping["tables"]["ring_daily"]["participant_id"] == "record_id"
    assert mapping["tables"]["ring_daily"]["day_index_from"] == "ring_date"
    assert ring["ring_rhr"]["measurement_class"] == "wearable_heart_rate"
    assert ring["ring_hrv_rmssd"]["measurement_class"] == "hrv"
    assert ring["ring_spo2"]["measurement_class"] == "continuous_spo2"
    assert ring["ring_sleep_hours"]["measurement_class"] == "sleep_objective"
    assert ring["ring_steps"]["measurement_class"] == "accelerometry"
    assert {m["block"] for c, m in ring.items() if m["role"] == "device"} == {"ring_daily"}
    for tname, t in mapping["tables"].items():
        for c, m in t["columns"].items():
            assert m["role"] in L.ROLES
            assert t["proposed_by"][c]["method"] in ("rule", "dictionary", "fuzzy", "llm", "none")


def test_mapping_yaml_roundtrip(redcap_run):
    res, out = redcap_run
    text = (out / "mapping.yaml").read_text()
    assert text.startswith("# mapping.yaml - PROPOSED")
    m = yaml.safe_load(text)
    assert m["approved_by"] is None and m["approved_at"] is None
    assert m["tables"]["records"]["columns"]["sex"]["value_map"] == {"1": "female", "2": "male"}
    assert m["dataset_id"] == "maestro_export"
    assert (out / ".gitignore").read_text().strip() == "*"
    assert not (out / "llm_requests.jsonl").exists()                       # no model without llm="claude"


def test_identifier_values_never_written(redcap_run, redcap_dir):
    res, out = redcap_run
    data = pd.read_csv(redcap_dir / "MAESTRO_DATA_2026-10-07.csv", dtype=str)
    blob = (out / "profile.json").read_text() + (out / "mapping.yaml").read_text()
    for col in ("zip", "visit_date", "notes"):
        for v in data[col].dropna().unique():
            assert v not in blob, (col, v)


def test_zip_archive_of_redcap(redcap_dir, tmp_path):
    z = tmp_path / "export.zip"
    with zipfile.ZipFile(z, "w") as zf:
        for f in redcap_dir.iterdir():
            zf.write(f, f"maestro/{f.name}")
        zf.writestr("../evil.csv", "a,b\n1,2\n")                           # zip-slip member: ignored
    ds = R.read_input(z)
    assert ds.format == "redcap" and ds.container == "zip"
    assert [t.name for t in ds.tables] == ["records", "ring_daily"]
    assert all("evil" not in t.source_file for t in ds.tables)


# ------------------------------------------------------------------------------------------------------------------
# other formats
# ------------------------------------------------------------------------------------------------------------------
def test_excel_two_sheets(tmp_path):
    x = make_excel(tmp_path / "cohort_workbook.xlsx")
    res = M.map_cmd(x, tmp_path / "out", condition_hint="me_cfs", echo=lambda *a: None)
    prof, mapping = res["profile"], res["mapping"]
    assert prof["format"] == "excel" and [t["name"] for t in prof["tables"]] == ["participants", "lab_results"]
    labs_t = next(t for t in prof["tables"] if t["name"] == "lab_results")
    assert labs_t["shape"] == "long" and labs_t["id_candidates"][0] == "subject_id"
    pc = _cols(mapping, "participants")
    assert pc["subject_id"]["role"] == "participant_id"
    assert pc["group"]["role"] == "label" and pc["group"]["condition"] == "me_cfs"
    assert pc["group"]["value_map"] == {"Control": 0, "ME/CFS": 1}
    assert pc["age_band"] == {"role": "age_band"}
    assert pc["gender"]["value_map"] == {"F": "female", "M": "male"}
    assert pc["BMI"]["role"] == "lab" and pc["BMI"]["loinc"] == "39156-5" and pc["BMI"]["unit"] is None
    lc = _cols(mapping, "lab_results")
    t = mapping["tables"]["lab_results"]
    assert t["long_lab"] == {"code_column": "test_name", "value_column": "result", "unit_column": "units",
                             "code_system": "loinc"}
    assert lc["test_name"]["code_map"] == {"C-reactive protein": "1988-5", "Creatinine": "2160-0",
                                            "Glucose": "2345-7", "Total cholesterol": "2093-3"}
    assert lc["test_name"]["unmapped_levels"] == ["Lab panel Q"]
    assert lc["result"] == {"role": "lab", "part": "value"} and lc["units"] == {"role": "lab", "part": "unit"}
    assert lc["collection_date"] == {"role": "date"}


def test_xpt_nhanes_like(tmp_path):
    x = make_xpt(tmp_path / "labdemo.xpt")
    ds = R.read_input(x)
    assert ds.format == "xpt" and ds.tables[0].labels["LBXTC"] == "Total Cholesterol (mg/dL)"
    res = M.map_cmd(x, tmp_path / "out", echo=lambda *a: None)
    c = _cols(res["mapping"], "labdemo")
    assert c["SEQN"] == {"role": "participant_id"}
    assert c["RIDAGEYR"] == {"role": "age", "transform": "age_band"}
    assert c["LBXTC"]["loinc"] == "2093-3" and c["LBXTC"]["unit"] == "mg/dL"
    assert c["LBXGLU"]["loinc"] == "1558-6"
    assert c["RIAGENDR"]["role"] == "sex" and c["RIAGENDR"]["value_map"] is None     # no value labels: never guessed
    assert c["RIAGENDR"]["unmapped_values"] == ["1", "2"]


def test_fhir_bundle(tmp_path):
    f = make_fhir(tmp_path / "bundle.json")
    ds = R.read_input(f)
    assert ds.format == "fhir"
    assert {t.name for t in ds.tables} == {"patients", "conditions", "observations", "medications"}
    assert any("Encounter" in w for w in ds.warnings)
    obs = ds.table("observations").df
    assert set(obs["patient_id"]) == {f"pat-{i}" for i in range(15)}
    res = M.map_cmd(f, tmp_path / "out", echo=lambda *a: None)
    m = res["mapping"]
    p = _cols(m, "patients")
    assert p["patient_id"]["role"] == "participant_id"
    assert p["birth_date"]["role"] == "drop" and p["postal_code"] == {"role": "drop", "reason": "geography"}
    assert p["state"] == {"role": "drop", "reason": "geography"}
    assert p["gender"]["value_map"] == {"female": "female", "male": "male"}
    c = _cols(m, "conditions")
    assert c["icd10cm"] == {"role": "diagnosis_code"} and c["snomed"]["role"] == "unmapped"
    o = _cols(m, "observations")
    assert m["tables"]["observations"]["long_lab"]["code_column"] == "loinc"
    assert o["value"] == {"role": "lab", "part": "value"}
    assert _cols(m, "medications")["rxnorm"] == {"role": "medication_code"}


def test_omop_folder(tmp_path):
    root = make_omop(tmp_path / "omop_cdm")
    ds = R.read_input(root)
    assert ds.format == "omop" and {t.kind for t in ds.tables} == {"omop:person", "omop:condition_occurrence",
                                                                     "omop:measurement"}
    res = M.map_cmd(root, tmp_path / "out", echo=lambda *a: None)
    m = res["mapping"]
    p = _cols(m, "person")
    assert p["person_id"]["role"] == "participant_id"
    assert p["gender_concept_id"]["value_map"] == {"8507": "male", "8532": "female"}
    assert p["person_source_value"] == {"role": "drop", "reason": "identifier"}
    assert p["location_id"] == {"role": "drop", "reason": "geography"}
    assert p["year_of_birth"]["role"] == "drop"
    co = _cols(m, "condition_occurrence")
    assert co["condition_source_value"] == {"role": "diagnosis_code"}
    assert co["condition_occurrence_id"]["role"] == "drop"
    assert m["tables"]["measurement"]["long_lab"] == {"code_column": "measurement_source_value",
                                                      "value_column": "value_as_number",
                                                      "unit_column": "unit_source_value", "code_system": "loinc"}


def test_text_delimiters_and_encodings(tmp_path):
    p = tmp_path / "semi.csv"
    p.write_bytes("pid;Glucose (mg/dL);café\nA1;95;x\nA2;101;y\n".encode("cp1252"))
    ds = R.read_input(p)
    df = ds.tables[0].df
    assert list(df.columns) == ["pid", "Glucose (mg/dL)", "café"] and df["Glucose (mg/dL)"].tolist() == [95, 101]
    assert any("cp1252" in w for w in ds.warnings)
    t = tmp_path / "tab.txt"
    t.write_text("a\tb\n1\t2\n")
    assert list(R.read_input(t).tables[0].df.columns) == ["a", "b"]


def test_json_nested_and_jsonl(tmp_path):
    p = tmp_path / "nested.json"
    p.write_text(json.dumps({"study": "x", "participants": [{"id": "a", "vitals": {"hr": 60}, "tags": [1, 2]},
                                                             {"id": "b", "vitals": {"hr": 70}, "tags": []}]}))
    ds = R.read_input(p)
    assert ds.tables[0].name == "participants" and "vitals.hr" in ds.tables[0].df.columns
    q = tmp_path / "rows.jsonl"
    q.write_text('{"a": 1, "b": {"c": 2}}\n{"a": 3, "b": {"c": 4}}\n')
    ds = R.read_input(q)
    assert ds.format == "jsonl" and ds.tables[0].df["b.c"].tolist() == [2, 4]


def test_spss_and_stata_value_labels(tmp_path):
    import pyreadstat
    df = pd.DataFrame({"pid": ["p1", "p2", "p3"], "raynaud": [1.0, 0.0, 1.0]})
    pyreadstat.write_sav(df, str(tmp_path / "s.sav"), column_labels=["Participant", "Raynaud's phenomenon"],
                         variable_value_labels={"raynaud": {1: "Yes", 0: "No"}})
    pyreadstat.write_dta(df, str(tmp_path / "s.dta"), column_labels=["Participant", "Raynaud's phenomenon"],
                         variable_value_labels={"raynaud": {1: "Yes", 0: "No"}})
    for f in ("s.sav", "s.dta"):
        t = R.read_input(tmp_path / f).tables[0]
        assert t.labels["raynaud"] == "Raynaud's phenomenon" and t.choices["raynaud"] == {"1": "Yes", "0": "No"}
    prof = P.profile_dataset(R.read_input(tmp_path / "s.sav"))
    mp = M.propose_mapping(prof, R.read_input(tmp_path / "s.sav"), dataset_id="sav_test")
    assert mp["tables"]["s"]["columns"]["raynaud"]["icd10cm"] == "I73.0"


def test_unsupported_and_missing_dependency(tmp_path, monkeypatch):
    (tmp_path / "x.docx").write_text("no")
    with pytest.raises(R.ReadError, match="unsupported"):
        R.read_input(tmp_path / "x.docx")
    x = make_xpt(tmp_path / "a.xpt")
    monkeypatch.setitem(sys.modules, "pyreadstat", None)
    with pytest.raises(R.MissingDependency, match="pyreadstat"):
        R.read_input(x)


def test_check_digits_and_code_patterns():
    assert P.luhn_ok("2345-7") and P.luhn_ok("82730-3") and not P.luhn_ok("2345-8")
    assert P.snomed_ok("38341003") and not P.snomed_ok("38341004")
    assert P.parse_unit("Glucose (mg/dL)") == "mg/dL" and P.parse_unit("crea_umol_l") == "umol/L"
    t = R.Table(name="t", df=pd.DataFrame({"dx": ["U09.9", "G93.32", "I73.00", "R53.83"]}), source_file="t.csv")
    assert P.detect_code_pattern(t, "dx", t.df["dx"], None) == "icd10"
    t2 = R.Table(name="t", df=pd.DataFrame({"n": [1, 2, 3, 6809]}), source_file="t.csv")
    assert P.detect_code_pattern(t2, "n", t2.df["n"], None) == "none"           # bare integers are not RxNorm
    assert P.detect_code_pattern(t2, "rxcui", t2.df["n"], None) == "rxnorm"


# ------------------------------------------------------------------------------------------------------------------
# LLM step (mock client; never the network)
# ------------------------------------------------------------------------------------------------------------------
class _MockMessages:
    def __init__(self, answer):
        self.calls = []
        self.answer = answer

    def create(self, **kw):
        self.calls.append(kw)
        payload = json.loads(kw["messages"][0]["content"])
        props = [self.answer(c) for c in payload["columns"]]
        return SimpleNamespace(stop_reason="end_turn", stop_details=None,
                               content=[SimpleNamespace(type="text", text=json.dumps({"proposals": props}))])


def _mock_client(answer):
    msgs = _MockMessages(answer)
    return SimpleNamespace(beta=SimpleNamespace(messages=msgs)), msgs


def _answer(col):
    base = {"key": col["key"], "code": None, "instrument": None, "item": None, "measurement_class": None,
            "omics_layer": None, "platform": None, "confidence": 0.9, "rationale": "test"}
    if col["column"] == "hba1c_pct":
        return {**base, "role": "lab", "code_system": "loinc", "code": "4548-4"}
    if col["column"] == "mystery_17":
        return {**base, "role": "lab", "code_system": "loinc", "code": "99999-5"}     # valid check digit, not local
    return {**base, "role": "unmapped", "code_system": "none"}


def test_llm_payload_has_no_row_values_and_is_validated(redcap_dir, tmp_path):
    client, msgs = _mock_client(_answer)
    provider = L.ClaudeProvider(client=client)
    assert provider.model == L.DEFAULT_MODEL == "claude-fable-5-1"
    out = tmp_path / "run"
    res = M.map_cmd(redcap_dir, out, llm=provider, condition_hint="long_covid", echo=lambda *a: None)
    assert len(msgs.calls) == 1
    kw = msgs.calls[0]
    assert kw["model"] == "claude-fable-5-1" and kw["output_config"]["format"]["type"] == "json_schema"
    assert kw["fallbacks"] == "default" and kw["betas"] == [L.FALLBACK_BETA]
    assert "thinking" not in kw and "temperature" not in kw and "tool_choice" not in kw
    sent = kw["messages"][0]["content"]
    payload = json.loads(sent)
    sent_cols = {c["column"] for c in payload["columns"]}
    assert sent_cols == {"hba1c_pct", "mystery_17"}                          # only still-unmapped columns
    data = pd.read_csv(redcap_dir / "MAESTRO_DATA_2026-10-07.csv", dtype=str)
    for col in ("record_id", "zip", "visit_date", "notes", "hba1c_pct", "mystery_17", "ring_date"):
        for v in data[col].dropna().unique():
            if col in ("hba1c_pct", "mystery_17"):
                continue                                                     # min/median/max are allowed aggregates
            assert v not in sent, (col, v)
    for c in payload["columns"]:
        assert set(c) <= {"key", "table", "column", "label", "dtype", "unit", "share_missing", "numeric", "levels",
                          "choices", "code_pattern", "field_type", "form"}
        for v, n in c.get("levels") or []:
            assert n is None or n >= 11
    log = [json.loads(x) for x in (out / "llm_requests.jsonl").read_text().splitlines()]
    assert len(log) == 1 and log[0]["payload"] == payload and log[0]["model"] == "claude-fable-5-1"
    assert "api_key" not in json.dumps(log).lower()
    rec = _cols(res["mapping"], "records")
    assert rec["hba1c_pct"] == {"role": "lab", "loinc": "4548-4", "unit": "%"}       # unit parsed from '_pct'
    assert rec["mystery_17"]["role"] == "unmapped"
    assert res["mapping"]["tables"]["records"]["proposed_by"]["hba1c_pct"]["method"] == "llm"
    assert rec["mystery_17"]["role"] == "unmapped" and "not in the local list" in rec["mystery_17"]["reason"]
    assert res["mapping"]["llm"]["n_accepted"] == 1 and res["mapping"]["llm"]["n_rejected"] == 1


def test_llm_small_levels_never_sent():
    table = {"name": "t", "id_candidates": ["pid"]}
    col = {"name": "fav_colour", "dtype": "category", "label": "Colour", "flags": [], "share_missing": 0.0,
           "levels": [["red", 40], ["blue", 12], ["<11 other", 14]], "choices": None, "code_pattern": "none"}
    it = L.column_item("c1", table, col)
    assert it["levels"] == [["red", 40], ["blue", 12], ["<11 other", None]]
    assert L.column_item("c2", table, {**col, "name": "zip_code"}) is None
    assert L.column_item("c3", table, {**col, "flags": ["date"]}) is None
    assert L.column_item("c4", table, {**col, "name": "pid"}) is None


def test_llm_skipped_without_key(redcap_dir, tmp_path, monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.delenv("ANTHROPIC_AUTH_TOKEN", raising=False)
    msgs = []
    assert L.get_provider("claude", echo=msgs.append) is None
    assert msgs and ("ANTHROPIC_API_KEY" in msgs[0] or "anthropic" in msgs[0])
    res = M.map_cmd(redcap_dir, tmp_path / "o", llm="claude", echo=lambda *a: None)
    assert res["mapping"]["llm"] is None and not (tmp_path / "o" / "llm_requests.jsonl").exists()
    with pytest.raises(ValueError):
        L.get_provider("gpt")


def test_validate_proposal_rules(redcap_run):
    prof = redcap_run[0]["profile"]
    mp = M.Mapper(prof)
    num = {"name": "x", "dtype": "float", "unit": None, "flags": []}
    yn = {"name": "y", "dtype": "bool", "choices": {"1": "Yes", "0": "No"}, "flags": [], "_table": "none"}
    base = {"code": None, "code_system": "none", "confidence": 0.9}
    assert mp.validate_proposal({**base, "role": "participant_id"}, num)[0] is None
    assert mp.validate_proposal({**base, "role": "boss"}, num)[0] is None
    assert mp.validate_proposal({**base, "role": "lab", "code_system": "loinc", "code": "2345-8"}, num)[0] is None
    assert mp.validate_proposal({**base, "role": "lab", "code_system": "loinc", "code": "2345-7"}, num)[0]["loinc"] == "2345-7"
    bad_unit = mp.validate_proposal({**base, "role": "lab", "code_system": "loinc", "code": "2345-7"},
                                    {**num, "unit": "kg"})
    assert bad_unit[0] is None and "unit" in bad_unit[1]
    ok = mp.validate_proposal({**base, "role": "self_report_history", "code_system": "icd10cm", "code": "I73.0"}, yn)
    assert ok[0] == {"role": "self_report_history", "icd10cm": "I73.0", "value_map": {"1": 1, "0": 0}}
    assert mp.validate_proposal({**base, "role": "self_report_history", "code_system": "icd10cm", "code": "Q99.99X"},
                                yn)[0] is None
    assert mp.validate_proposal({**base, "role": "device", "measurement_class": "telepathy"}, num)[0] is None
    assert mp.validate_proposal({**base, "role": "device", "measurement_class": "hrv"}, num)[0]["measurement_class"] == "hrv"
    assert mp.validate_proposal({**base, "role": "medication_flag", "code_system": "rxnorm", "code": "abc"}, yn)[0] is None


def test_vocab_without_processed_tables(tmp_path):
    v = M.Vocab(processed=tmp_path)                                         # a fresh clone: no data/processed
    codes, src = v.icd_codes
    assert "I730" in codes and "U099" in codes and "configs" in src
    assert {"2093-3", "82730-3", "4548-4"} <= v.loinc_codes
    assert M.resolve_condition_text("HC", v)[0] == "control"
    assert M.resolve_condition_text("CLYME", v)[0] == "ptlds"
    assert M.resolve_condition_text("ME/CFS", v)[0] == "me_cfs"
