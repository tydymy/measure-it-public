"""Launch agent part B (measure_it.agent.apply / omop / omop_run): approval gate, de-identification transforms, the
BYOD folder passing its validator, the All of Us-shaped OMOP export and the first execution of a stage-D query pack.

Everything is SYNTHETIC and built here (tables + a mapping dict per docs/AGENT_CONTRACT.md §2), so the tests run on a
fresh clone without part A. The query-pack test needs sqlglot (`uv run --with sqlglot pytest ...`) and is skipped
with a reason otherwise.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from measure_it.agent import apply as A
from measure_it.agent import omop as O
from measure_it.agent import omop_run as R
from measure_it.byod import checks as C

HAS_SQLGLOT = importlib.util.find_spec("sqlglot") is not None
NAMES = ["Alice Smith", "Bob Jones", "Carol White", "Dan Brown", "Eve Black"]


def synthetic_tables(n_cases: int = 30, n_controls: int = 30, seed: int = 7) -> dict[str, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    n = n_cases + n_controls
    ids = [f"S{i + 1:03d}" for i in range(n)]
    case = np.r_[np.ones(n_cases), np.zeros(n_controls)].astype(bool)
    start = pd.Timestamp("2023-03-01") + pd.to_timedelta(rng.integers(0, 300, n), unit="D")
    records = pd.DataFrame({
        "record_id": ids,
        "name": [NAMES[i % 5] for i in range(n)],
        "zip": [f"0{2100 + i:04d}" for i in range(n)],
        "notes": ["called 617-555-0199 about visit" for _ in range(n)],
        "visit_date": (start + pd.to_timedelta(3, unit="D")).strftime("%Y-%m-%d"),
        "age": rng.integers(22, 95, n),
        "sex": rng.choice(["F", "M"], n),
        "race": rng.choice(["White", "Black", "Asian"], n),
        "cohort": np.where(case, "LC", "HC"),
        "raynaud_hx": np.where(rng.random(n) < np.where(case, 0.4, 0.1), "Yes", "No"),
        "glyca": np.round(0.8 + 0.05 * case + rng.normal(0, 0.06, n), 4),
        "compass31_vaso": np.round(np.clip(1 + 1.5 * case + rng.normal(0, 0.7, n), 0, 4.17), 2),
        "cap_abn": np.round(np.clip(1.5 + 0.8 * case + rng.normal(0, 0.5, n), 0, 4), 2),
        "cap_density": np.round(8.5 - 0.7 * case + rng.normal(0, 0.8, n), 2),
        "olink_il6": np.round(rng.normal(4 + 0.4 * case, 1, n), 3),
        "mystery_17": rng.normal(0, 1, n),
    })
    labs, dx, rx, ring = [], [], [], []
    for i, pid in enumerate(ids):
        labs.append({"record_id": pid, "lab_date": (start[i] + pd.Timedelta(days=10)).strftime("%m/%d/%Y"),
                     "loinc_code": "2093-3", "result": round(float(rng.normal(190, 30)), 1), "units": "mg/dL",
                     "test_name": "Cholesterol"})
        labs.append({"record_id": pid, "lab_date": (start[i] + pd.Timedelta(days=10)).strftime("%m/%d/%Y"),
                     "loinc_code": "not-a-code", "result": 1.0, "units": "x", "test_name": "junk"})
        if case[i]:
            dx.append({"record_id": pid, "dx_date": start[i].strftime("%Y-%m-%d"), "dx_code": "U099"})
            if rng.random() < 0.4:
                dx.append({"record_id": pid, "dx_date": start[i].strftime("%Y-%m-%d"), "dx_code": "I73.0"})
            if rng.random() < 0.5:
                rx.append({"record_id": pid, "rx_code": "6918"})
        else:
            dx.append({"record_id": pid, "dx_date": start[i].strftime("%Y-%m-%d"), "dx_code": "I10"})
        for d in range(5):
            ring.append({"record_id": pid, "day_date": (start[i] + pd.Timedelta(days=d)).strftime("%Y-%m-%d"),
                         "steps": int(rng.integers(2000, 12000)), "hr_mean_bpm": round(float(rng.normal(68, 5)), 1),
                         "resting_hr": round(float(rng.normal(58, 5)), 1),
                         "hrv_rmssd_ms": round(float(rng.normal(40, 8)), 1),
                         "sleep_minutes": int(rng.integers(300, 500))})
    return {"records": records, "labs": pd.DataFrame(labs), "dx": pd.DataFrame(dx), "rx": pd.DataFrame(rx),
            "ring": pd.DataFrame(ring)}


def synthetic_mapping(approved: bool = False) -> dict:
    return {
        "dataset_id": "agent_pytest_synth",
        "condition_hint": "long_covid",
        "measurement_hint": "nailfold_capillaroscopy",
        "approved_by": "pytest" if approved else None,
        "approved_at": "2026-10-07T00:00:00+00:00" if approved else None,
        "synthetic": True,
        "source": "SYNTHETIC tables built in tests/test_agent_apply_omop.py",
        "tables": {
            "records": {"participant_id": "record_id", "day_index_from": "visit_date", "columns": {
                "record_id": {"role": "participant_id"},
                "cohort": {"role": "label", "condition": "long_covid", "value_map": {"LC": 1, "HC": 0},
                           "label_basis": "clinical_case_definition"},
                "age": {"role": "age", "transform": "age_band"},
                "sex": {"role": "sex", "value_map": {"F": "female", "M": "male"}},
                "race": {"role": "demographic"},
                "raynaud_hx": {"role": "self_report_history", "icd10cm": "I73.0", "condition": "Raynaud's phenomenon",
                               "value_map": {"Yes": 1, "No": 0}},
                "glyca": {"role": "lab", "loinc": "82730-3", "unit": "mmol/L", "platform": "nmr_nightingale"},
                "compass31_vaso": {"role": "survey", "instrument": "compass31", "item": "vasomotor"},
                "cap_abn": {"role": "device", "block": "capillaroscopy", "feature": "abnormality_score",
                            "measurement_class": "capillaroscopy"},
                "cap_density": {"role": "device", "block": "capillaroscopy", "measurement_class": "capillaroscopy"},
                "olink_il6": {"role": "omics", "layer": "proteomics", "platform": "olink"},
                "zip": {"role": "drop", "reason": "geography"},
                "name": {"role": "drop", "reason": "identifier"},
                "notes": {"role": "drop", "reason": "free_text"},
                "visit_date": {"role": "date"},
                "mystery_17": {"role": "unmapped", "reason": "no rule, fuzzy score 0.41 < 0.80"}}},
            "labs": {"participant_id": "record_id", "day_index_from": "lab_date",
                     "long_labs": {"code": "loinc_code", "value": "result", "unit": "units", "name": "test_name"},
                     "columns": {"record_id": {"role": "participant_id"}, "lab_date": {"role": "date"}}},
            "dx": {"participant_id": "record_id", "day_index_from": "dx_date",
                   "columns": {"record_id": {"role": "participant_id"}, "dx_code": {"role": "diagnosis_code"}}},
            "rx": {"participant_id": "record_id", "columns": {"record_id": {"role": "participant_id"},
                                                              "rx_code": {"role": "medication_code"}}},
            "ring": {"participant_id": "record_id", "day_index_from": "day_date", "columns": {
                "record_id": {"role": "participant_id"},
                **{c: {"role": "device", "block": "ring", "measurement_class": "wearable_heart_rate"}
                   for c in ("steps", "hr_mean_bpm", "resting_hr", "hrv_rmssd_ms", "sleep_minutes")}}},
        },
        "subgroup_analysis": {"name": "capillaroscopy_abnormal", "within_label": "long_covid",
                              "device_positive": {"block": "capillaroscopy", "feature": "abnormality_score",
                                                  "threshold": 2.0, "direction": ">="},
                              "ehr_domains": ["condition", "measurement", "survey", "demographic"]},
    }


def synthetic_phenotype() -> dict:
    return {
        "phenotype_id": "byod_agent_pytest_synth:capillaroscopy_abnormal", "version": 1, "condition_id": "long_covid",
        "base_population": {"description": "Long COVID (U09.9)", "icd10cm": ["U09.9"], "label_column": "long_covid"},
        "intercept": -0.4, "decision_threshold": 0.4,
        "features": [
            {"feature_key": "ICD10CM:I73", "code_match": "category_prefix", "codes": ["I73.0"], "transform": "presence",
             "coefficient": 0.9, "coefficient_standardized": 0.4, "center": 0.3, "scale": 0.46},
            {"feature_key": "LOINC:2093-3", "codes": ["2093-3"], "transform": "standardized_value", "center": 190.0,
             "scale": 30.0, "unit": "mg/dL", "coefficient": 0.2, "coefficient_standardized": 0.2},
            {"feature_key": "RXNORM:6918", "codes": ["6918"], "transform": "presence", "coefficient": 0.3,
             "center": 0.2, "scale": 0.4},
            {"feature_key": "DEMOG:age_mid", "codes": ["age_mid"], "transform": "standardized_value", "center": 50.0,
             "scale": 15.0, "coefficient": 0.1},
            {"feature_key": "DEMOG:female", "codes": ["female"], "transform": "presence", "coefficient": 0.2,
             "center": 0.5, "scale": 0.5},
            {"feature_key": "SURVEY:COMPASS31:vasomotor", "codes": ["COMPASS31:vasomotor"],
             "transform": "standardized_value", "center": 2.0, "scale": 1.0, "coefficient": 0.6},
        ],
    }


@pytest.fixture(scope="module")
def applied(tmp_path_factory):
    out = tmp_path_factory.mktemp("agent_run")
    tables = synthetic_tables()
    report = A.apply_mapping(tables, synthetic_mapping(), out, approve=True, public=True)
    return out, tables, report


# ------------------------------------------------------------------------------------------------ approval gate

def test_unapproved_mapping_is_refused(tmp_path):
    with pytest.raises(A.ApplyRefused) as e:
        A.apply_mapping(synthetic_tables(), synthetic_mapping(), tmp_path / "run", approve=False)
    assert "not approved" in str(e.value)
    assert not (tmp_path / "run").exists()


def test_approve_stamps_and_preapproved_passes(applied, tmp_path):
    out, _, report = applied
    assert report["approved_by"] and report["approved_at"]
    rep2 = A.apply_mapping(synthetic_tables(), synthetic_mapping(approved=True), tmp_path / "r2", approve=False,
                           public=True, omop=False)
    assert rep2["approved_by"] == "pytest"


def test_identifying_ids_need_salt(tmp_path, monkeypatch):
    t = synthetic_tables()
    mrn = {f"S{i + 1:03d}": f"{1000000 + i}" for i in range(60)}
    for df in t.values():
        df["record_id"] = df["record_id"].map(mrn)
    monkeypatch.delenv("MEASURE_IT_BYOD_SALT", raising=False)
    with pytest.raises(A.ApplyRefused) as e:
        A.apply_mapping(t, synthetic_mapping(), tmp_path / "r", approve=True, public=True, omop=False)
    assert "MEASURE_IT_BYOD_SALT" in str(e.value)
    monkeypatch.setenv("MEASURE_IT_BYOD_SALT", "pytest-salt-0123456789")
    rep = A.apply_mapping(t, synthetic_mapping(), tmp_path / "r", approve=True, public=True, omop=False)
    assert rep["participant_id_transform"].startswith("salted_sha256")
    P = pd.read_csv(tmp_path / "r" / "byod" / "participants.csv", dtype=str)
    assert P["participant_id"].str.fullmatch(r"h[0-9a-f]{16}").all()
    assert not set(P["participant_id"]) & set(mrn.values())


def test_private_data_never_asserts_permission(tmp_path):
    rep = A.apply_mapping(synthetic_tables(), synthetic_mapping(), tmp_path / "p", approve=True, public=False,
                          omop=False)
    assert rep["data_use_confirmation_required"] is True
    assert not rep["ok"] and any("data_use" in e for e in rep["validation"]["errors"])
    assert rep["errors_other_than_data_use"] == []
    m = synthetic_mapping()
    m["data_use"] = {"deidentified": True, "use_permitted": True,
                     "statement": "SYNTHETIC test data; confirmed by the test as owner"}
    rep2 = A.apply_mapping(synthetic_tables(), m, tmp_path / "p2", approve=True, omop=False)
    assert rep2["ok"] and rep2["data_use_confirmation_required"] is False


# ------------------------------------------------------------------------------------------------ BYOD folder

def test_byod_validator_passes(applied):
    out, _, report = applied
    assert report["ok"], report["validation"]["errors"]
    assert C.check_dir(out / "byod").ok
    assert (out / ".gitignore").read_text().strip().endswith("*")
    files = set(report["files"])
    assert {"participants.csv", "device_capillaroscopy.csv", "device_ring.csv", "ehr_labs.csv", "ehr_diagnoses.csv",
            "ehr_medications.csv", "survey_compass31.csv", "omics_proteomics.csv", "omics_nightingale.csv",
            "self_report_history.csv"} <= files
    m = report["primary_analysis"]
    assert m["label"] == "long_covid" and m["feature_blocks"] == ["device_capillaroscopy"]
    assert report["measurement_class"] == "capillaroscopy"
    assert {d["column"].split(".")[-1] for d in report["dropped"]} == {"zip", "name", "notes"}
    assert any(u["column"] == "records.mystery_17" for u in report["unmapped"])


def _all_text(root: Path) -> str:
    parts = []
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.suffix in (".csv", ".yaml", ".md", ".json"):
            parts.append(p.read_text())
        elif p.is_file() and p.suffix == ".parquet":
            df = pd.read_parquet(p)
            parts.append(" ".join(map(str, df.columns)) + "\n" + df.astype(str).to_csv(index=False))
    return "\n".join(parts)


def test_identifiers_and_geography_absent(applied):
    out, tables, _ = applied
    rec = tables["records"]
    for sub in ("byod", "omop"):
        text = _all_text(out / sub)
        for v in list(rec["name"].unique()) + list(rec["zip"]) + ["617-555-0199"]:
            assert v not in text, (sub, "identifier/geography value written")
        for v in list(rec["visit_date"]) + list(tables["labs"]["lab_date"]) + list(tables["ring"]["day_date"]):
            assert v not in text, (sub, "calendar date written")
    for p in (out / "byod").glob("*.csv"):
        cols = set(pd.read_csv(p, nrows=0).columns)
        assert not cols & {"name", "zip", "notes", "visit_date", "age", "lab_date", "day_date", "record_id"}, p.name
        for c in cols:
            assert C.column_problem(c) is None, (p.name, c)


def test_dates_and_ages_converted(applied):
    out, tables, _ = applied
    P = pd.read_csv(out / "byod" / "participants.csv", dtype=str)
    assert P["age_band"].dropna().map(lambda v: C.age_band_problem(v) is None).all()
    assert set(P["sex"]) <= {"female", "male"}
    assert set(P["long_covid"]) == {"0", "1"}
    ring = pd.read_csv(out / "byod" / "device_ring.csv")
    assert ring.groupby("participant_id")["day_index"].min().eq(0).all()   # ring day 0 is each person's first date
    assert ring["day_index"].max() == 4
    labs = pd.read_csv(out / "byod" / "ehr_labs.csv")
    assert (labs["day_index"] == 10).all() and set(labs["loinc"]) == {"2093-3"}
    dx = pd.read_csv(out / "byod" / "ehr_diagnoses.csv")
    assert "U09.9" in set(dx["icd10cm"]) and (dx["day_index"] == 0).all()
    # exact ages >= 90 are top-coded
    age = tables["records"].set_index("record_id")["age"]
    old = age[age >= 90].index
    assert (P.set_index("participant_id").loc[old, "age_band"] == "90+").all()


def test_nmr_lab_routed_with_platform(applied):
    out, _, report = applied
    import yaml
    m = yaml.safe_load((out / "byod" / "manifest.yaml").read_text())
    assert m["omics"]["nightingale"] == {"platform": "nmr_nightingale", "units": "nightingale_standard",
                                         "description": "NMR measures routed from lab-role columns"}
    assert "GlycA" in pd.read_csv(out / "byod" / "omics_nightingale.csv", nrows=0).columns
    assert m["data_use"]["statement"].startswith("public dataset: ")
    assert m["comparator"]["type"] == "healthy" and report["caveats"]


# ------------------------------------------------------------------------------------------------ OMOP export

def test_omop_tables(applied):
    out, _, report = applied
    om = out / "omop"
    expected = {"person": O.PERSON_COLS, "observation_period": O.OBS_PERIOD_COLS,
                "condition_occurrence": O.CONDITION_COLS, "measurement": O.MEASUREMENT_COLS,
                "drug_exposure": O.DRUG_COLS, "observation": O.OBSERVATION_COLS, "concept": O.CONCEPT_COLS,
                "activity_summary": O.ACTIVITY_COLS, "heart_rate_summary": O.HR_SUMMARY_COLS,
                "sleep_daily_summary": O.SLEEP_COLS}
    for name, cols in expected.items():
        df = pd.read_parquet(om / f"{name}.parquet")
        assert list(df.columns) == cols, name
        assert len(df), name
    concept = pd.read_parquet(om / "concept.parquet")
    assert concept["concept_id"].min() >= 2_000_000_000 and concept["concept_id"].is_unique
    local = set(concept["concept_id"])
    for name in ("person", "condition_occurrence", "measurement", "drug_exposure", "observation",
                 "observation_period"):
        df = pd.read_parquet(om / f"{name}.parquet")
        for c in [c for c in df.columns if c.endswith("concept_id")]:
            v = set(df[c].dropna().astype(int)) - {0}
            assert v <= local, (name, c)
    person = pd.read_parquet(om / "person.parquet")
    assert person["year_of_birth"].isna().all() and person["age_band_source_value"].notna().any()
    co = pd.read_parquet(om / "condition_occurrence.parquet")
    assert pd.to_datetime(co["condition_start_date"]).min() >= pd.Timestamp("2000-01-01")
    m = pd.read_parquet(om / "measurement.parquet")
    assert m["measurement_source_value"].str.contains("platform=nmr_nightingale").any()
    obs = pd.read_parquet(om / "observation.parquet")
    assert "COMPASS31:vasomotor" in set(obs["observation_source_value"])
    assert "personal_medical_history" in set(obs["observation_source_value"])
    vocab = set(concept["vocabulary_id"])
    assert {"ICD10CM", "LOINC", "RxNorm", "UCUM", "Gender"} <= vocab
    assert (om / "README.md").exists()
    hr = pd.read_parquet(om / "heart_rate_summary.parquet")
    assert hr["mean_heart_rate_device"].notna().all() and hr["min_heart_rate"].notna().all()
    assert "hrv_rmssd_ms" in report["omop"]["fitbit_unmapped_columns"]["ring"]


def test_omop_from_example_byod_folder(tmp_path):
    spec = importlib.util.spec_from_file_location("byod_make_example", Path(__file__).resolve().parents[1]
                                                  / "examples" / "byod_synthetic" / "make_example.py")
    mk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mk)
    src = mk.make_maestro(tmp_path / "maestro", fast=True, dataset_id="agent_maestro")
    s = O.export_omop_from_byod(src, tmp_path / "omop")
    assert s["tables"]["person"] > 0 and s["tables"]["measurement"] > 0 and s["tables"]["observation"] > 0
    assert "activity_summary" in s["tables"] and "sleep_daily_summary" in s["tables"]


# ------------------------------------------------------------------------------------------------ query pack

def test_query_pack_check_end_to_end(applied):
    if not HAS_SQLGLOT:
        pytest.skip("sqlglot is not installed (not a project dependency); run "
                    "`uv run --with sqlglot pytest tests/test_agent_apply_omop.py`")
    out, _, _ = applied
    res = R.check_query_pack(out / "omop", synthetic_phenotype())
    assert (out / "omop" / "query_pack_check.json").exists()
    assert res["status"] == "executed", json.dumps(res.get("variants", res), default=str)[:3000]
    steps = {s["step"]: s["rows"] for s in res["variants"]["year_of_birth_from_age_band"]["steps"]}
    assert steps["base"] == 30                              # every study-label case carries U09.9
    assert 0 < steps["cohort"] <= 30
    assert steps["f1"] > 0                                  # LOINC 2093-3 in mg/dL found through the UCUM filter
    exported = {s["step"]: s["rows"] for s in res["variants"]["as_exported"]["steps"]}
    assert exported["base"] == 30 and exported["cohort"] == 0   # year_of_birth is null in the export
    assert res["cdm_tables_missing"] == []


def test_query_pack_check_reports_missing_sqlglot(applied, monkeypatch):
    import builtins
    real = builtins.__import__

    def fake(name, *a, **k):
        if name == "sqlglot":
            raise ImportError("no sqlglot")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake)
    out, _, _ = applied
    res = R.check_query_pack(out / "omop", synthetic_phenotype())
    assert res["status"] == "sqlglot_missing" and "uv run --with sqlglot" in res["error"]
