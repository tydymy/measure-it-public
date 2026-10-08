"""Launch agent part C (measure_it.agent.orchestrate / launch / report / cli): unit tests and an end-to-end run on a
SYNTHETIC, MAESTRO-like messy input (REDCap-style export + data dictionary; no EHR codes; cohorts Long COVID / chronic
Lyme / acute Lyme / healthy; capillaroscopy features with a device-positive column; questionnaires; Nightingale-like
NMR columns).

The end-to-end tests run offline with `llm=None`, mark the input synthetic (so the dataset is demo-tagged), and check
that the teardown leaves the processed tables as they were. They need data/processed (`@pytest.mark.data`).
Part A's mapper is used when `measure_it.agent.mapper` exists (followed by the human edit the workflow expects:
labels, device block and subgroup plan set by hand); otherwise a stub mapper local to this file writes the mapping.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from measure_it.agent import orchestrate as OR
from measure_it.agent import report as RP
from measure_it.config import PROCESSED

HAS_DATA = (PROCESSED / "deployment_opportunities.parquet").exists()
SEED = 20261007
COHORTS = {1: ("long_covid", 90), 2: ("chronic_lyme", 50), 3: ("acute_lyme", 40), 4: ("healthy", 80)}
# unique per test process: another process running these tests at the same time works on its own dataset id (a launch
# run refuses a dataset id that is already ingested and never tears down data it did not create)
DS_ID = f"launch_pytest_{os.getpid()}"


# ------------------------------------------------------------------------------------------------------------------
# the SYNTHETIC messy input
# ------------------------------------------------------------------------------------------------------------------
DICT_COLS = ["Variable / Field Name", "Form Name", "Section Header", "Field Type", "Field Label",
             "Choices, Calculations, OR Slider Labels", "Field Note", "Text Validation Type OR Show Slider Number",
             "Text Validation Min", "Text Validation Max", "Identifier?"]


def make_redcap_input(folder: Path, seed: int = SEED) -> Path:
    """SYNTHETIC REDCap export + data dictionary (nothing here is a real person or MAESTRO data)."""
    rng = np.random.default_rng(seed)
    folder.mkdir(parents=True, exist_ok=True)
    coh = np.concatenate([[k] * n for k, (_, n) in COHORTS.items()])
    coh = coh[rng.permutation(len(coh))]
    n = len(coh)
    lc = (coh == 1).astype(float)
    sick = np.isin(coh, [1, 2]).astype(float) + 0.5 * (coh == 3)
    m = rng.normal(0, 1, n) + np.select([coh == 1, coh == 2, coh == 4], [0.3, 0.1, -0.8], 0.0)   # microvascular
    female = rng.random(n) < np.where(coh == 4, 0.55, 0.68)
    age = np.clip(rng.normal(np.where(coh == 3, 47, 42), 12), 19, 84).round().astype(int)
    has_cap = rng.random(n) < 0.9
    abn = np.round(np.clip(1.5 + 0.9 * m + rng.normal(0, 0.45, n), 0, 4), 2)
    rec = pd.DataFrame({
        "record_id": [f"SYN-{i:04d}" for i in range(1, n + 1)],
        "redcap_event_name": "baseline_arm_1",
        "cohort": coh, "age": age, "sex": np.where(female, 1, 2),
        "zip_code": rng.choice(["02139", "21205", "10027"], n),
        "enroll_date": pd.to_datetime("2024-01-01") + pd.to_timedelta(rng.integers(0, 300, n), unit="D"),
        "initials": rng.choice(["AB", "CD", "EF"], n),
        "notes": rng.choice(["", "felt dizzy at visit", "rescheduled once"], n),
        "hx_raynaud": (rng.random(n) < 1 / (1 + np.exp(-(-2.6 + 1.6 * m * lc + 0.6 * sick)))).astype(int),
        "hx_covid": (rng.random(n) < np.where(lc == 1, 0.95, 0.45)).astype(int),
        "hx_lyme": (rng.random(n) < np.where(np.isin(coh, [2, 3]), 0.95, 0.03)).astype(int),
        "hx_pots": (rng.random(n) < 0.05 + 0.2 * sick).astype(int),
        "hx_migraine": (rng.random(n) < 0.12 + 0.15 * sick).astype(int),
        "cap_abnormality_score": np.where(has_cap, abn, np.nan),
        "cap_density": np.where(has_cap, np.round(8.5 - 0.7 * m + rng.normal(0, 0.8, n), 2), np.nan),
        "cap_giant_n": np.where(has_cap, rng.poisson(np.exp(-1.0 + 0.6 * m)), np.nan),
        "cap_hemorrhage_n": np.where(has_cap, rng.poisson(np.exp(-1.2 + 0.5 * m)), np.nan),
        "cap_tortuosity": np.where(has_cap, np.round(np.clip(1.2 + 0.5 * m + rng.normal(0, 0.6, n), 0, 3), 2),
                                   np.nan),
        "cap_positive": np.where(has_cap, (abn >= 2.0).astype(float), np.nan),
        "compass31_orthostatic": np.round(np.clip(14 * sick + rng.normal(4, 6, n), 0, 40), 1),
        "compass31_vasomotor": np.round(np.clip(0.8 + 0.9 * sick + 0.55 * lc * m + rng.normal(0, 0.6, n), 0, 4.17), 2),
        "fss9_total": np.round(np.clip(2.5 + 2.6 * sick + rng.normal(0, 1, n), 1, 7), 2),
        "dsqpem_frequency": np.round(np.clip(15 + 45 * sick + rng.normal(0, 15, n), 0, 100), 1),
        "dsqpem_severity": np.round(np.clip(12 + 40 * sick + rng.normal(0, 15, n), 0, 100), 1),
        "ng_total_c": np.round(rng.normal(5.0, 0.9, n), 3), "ng_hdl_c": np.round(rng.normal(1.5, 0.35, n), 3),
        "ng_tg": np.round(np.exp(rng.normal(0.2, 0.4, n)), 3), "ng_glucose": np.round(rng.normal(5.2, 0.5, n), 3),
        "ng_glyca": np.round(0.78 + 0.03 * sick + 0.07 * lc * m + rng.normal(0, 0.06, n), 4),
        "ng_creatinine": np.round(rng.normal(75, 14, n), 2),
        "q17_misc": rng.integers(0, 5, n),
        "baseline_complete": 2,
    })
    rec["enroll_date"] = rec["enroll_date"].dt.strftime("%Y-%m-%d")
    rec.to_csv(folder / "MAESTRO_SYNTHETIC_DATA_2026-10-07.csv", index=False)
    yn = ("yesno", "")
    fields = [
        ("record_id", "text", "Record ID", "", "", "", ""),
        ("cohort", "radio", "Study cohort", "1, Long COVID | 2, Chronic Lyme (persistent symptoms after treated "
                                           "Lyme) | 3, Acute Lyme | 4, Healthy volunteer", "", "", ""),
        ("age", "text", "Age at enrollment (years)", "", "", "integer", ""),
        ("sex", "radio", "Sex assigned at birth", "1, Female | 2, Male", "", "", ""),
        ("zip_code", "text", "Home ZIP code", "", "", "zipcode", "y"),
        ("enroll_date", "text", "Enrollment date", "", "", "date_ymd", "y"),
        ("initials", "text", "Participant initials", "", "", "", "y"),
        ("notes", "notes", "Coordinator notes", "", "", "", ""),
        ("hx_raynaud", *yn, "", "", "", ""), ("hx_covid", *yn, "", "", "", ""), ("hx_lyme", *yn, "", "", "", ""),
        ("hx_pots", *yn, "", "", "", ""), ("hx_migraine", *yn, "", "", "", ""),
        ("cap_abnormality_score", "text", "Nailfold capillaroscopy abnormality score (0-4)", "", "", "number", ""),
        ("cap_density", "text", "Capillary density", "", "per mm", "number", ""),
        ("cap_giant_n", "text", "Giant capillaries (count)", "", "", "integer", ""),
        ("cap_hemorrhage_n", "text", "Microhemorrhages (count)", "", "", "integer", ""),
        ("cap_tortuosity", "text", "Tortuosity score (0-3)", "", "", "number", ""),
        ("cap_positive", "yesno", "Abnormal capillaroscopy pattern (reader call)", "", "", "", ""),
        ("compass31_orthostatic", "text", "COMPASS-31 orthostatic intolerance domain", "", "", "number", ""),
        ("compass31_vasomotor", "text", "COMPASS-31 vasomotor domain", "", "", "number", ""),
        ("fss9_total", "text", "Fatigue Severity Scale (FSS-9) mean score", "", "", "number", ""),
        ("dsqpem_frequency", "text", "DSQ-PEM frequency score", "", "", "number", ""),
        ("dsqpem_severity", "text", "DSQ-PEM severity score", "", "", "number", ""),
        ("ng_total_c", "text", "Total cholesterol (NMR)", "", "mmol/L", "number", ""),
        ("ng_hdl_c", "text", "HDL cholesterol (NMR)", "", "mmol/L", "number", ""),
        ("ng_tg", "text", "Total triglycerides (NMR)", "", "mmol/L", "number", ""),
        ("ng_glucose", "text", "Glucose (NMR)", "", "mmol/L", "number", ""),
        ("ng_glyca", "text", "GlycA glycoprotein acetyls (NMR)", "", "mmol/L", "number", ""),
        ("ng_creatinine", "text", "Creatinine (NMR)", "", "umol/L", "number", ""),
        ("q17_misc", "text", "Q17", "", "", "integer", ""),
    ]
    labels = {"hx_raynaud": "History: Raynaud's phenomenon", "hx_covid": "History: COVID-19",
              "hx_lyme": "History: Lyme disease", "hx_pots": "History: POTS", "hx_migraine": "History: migraine"}
    rows = []
    for f in fields:
        name, ftype, label, choices, note, val, ident = f
        rows.append({"Variable / Field Name": name, "Form Name": "baseline", "Section Header": "",
                     "Field Type": ftype, "Field Label": labels.get(name, label),
                     "Choices, Calculations, OR Slider Labels": choices, "Field Note": note,
                     "Text Validation Type OR Show Slider Number": val, "Text Validation Min": "",
                     "Text Validation Max": "", "Identifier?": ident})
    pd.DataFrame(rows, columns=DICT_COLS).to_csv(folder / "MAESTRO_SYNTHETIC_DataDictionary.csv", index=False)
    (folder / "README_SYNTHETIC.txt").write_text("SYNTHETIC test input for tests/test_agent_launch.py; no real person.\n")
    return folder


HX = {"hx_raynaud": ("Raynaud's phenomenon", "I73.00"), "hx_covid": ("COVID-19", "U07.1"),
      "hx_lyme": ("Lyme disease", "A69.20"), "hx_pots": ("POTS", "G90.A"), "hx_migraine": ("migraine", "G43.909")}
NG = {"ng_total_c": ("2093-3", "mmol/L"), "ng_hdl_c": ("2085-9", "mmol/L"), "ng_tg": ("2571-8", "mmol/L"),
      "ng_glucose": ("2345-7", "mmol/L"), "ng_glyca": ("82730-3", "mmol/L"), "ng_creatinine": ("2160-0", "umol/L")}
CAP = {"cap_abnormality_score": "abnormality_score", "cap_density": "capillary_density_per_mm",
       "cap_giant_n": "giant_capillaries_n", "cap_hemorrhage_n": "microhemorrhages_n",
       "cap_tortuosity": "tortuosity_score", "cap_positive": "device_positive"}
SURVEY = {"compass31_orthostatic": ("compass31", "orthostatic"), "compass31_vasomotor": ("compass31", "vasomotor"),
          "fss9_total": ("fss9", "total"), "dsqpem_frequency": ("dsqpem", "frequency_score"),
          "dsqpem_severity": ("dsqpem", "severity_score")}


def reviewed_columns() -> dict:
    """The column specs a human reviewer settles for this export (contract section 2 shapes)."""
    cols = {
        "record_id": {"role": "participant_id"},
        "cohort": {"role": "label", "labels": [
            {"column": "long_covid", "condition": "long_covid", "value_map": {"1": 1, "4": 0},
             "label_basis": "clinical_diagnosis", "definition": "SYNTHETIC Long COVID cohort vs healthy volunteers"},
            {"column": "ptlds", "condition": "ptlds", "value_map": {"2": 1, "4": 0},
             "label_basis": "clinical_case_definition",
             "definition": "SYNTHETIC chronic Lyme (persistent symptoms after treated Lyme) vs healthy volunteers"},
            {"column": "lyme_disease", "condition": "lyme_disease", "value_map": {"3": 1, "4": 0},
             "label_basis": "clinical_diagnosis", "definition": "SYNTHETIC acute Lyme cohort vs healthy volunteers"}]},
        "age": {"role": "age", "transform": "age_band"},
        "sex": {"role": "sex", "value_map": {"1": "female", "2": "male"}},
        "redcap_event_name": {"role": "drop", "reason": "structural"},
        "baseline_complete": {"role": "drop", "reason": "structural"},
        "zip_code": {"role": "drop", "reason": "geography"},
        "enroll_date": {"role": "drop", "reason": "date"},
        "initials": {"role": "drop", "reason": "identifier"},
        "notes": {"role": "drop", "reason": "free_text"},
        "q17_misc": {"role": "unmapped", "reason": "no rule, no dictionary label"},
    }
    for c, (name, code) in HX.items():
        cols[c] = {"role": "self_report_history", "condition": name, "icd10cm": code, "value_map": {"1": 1, "0": 0}}
    for c, feat in CAP.items():
        cols[c] = {"role": "device", "block": "capillaroscopy", "feature": feat, "measurement_class": "capillaroscopy",
                   "description": "SYNTHETIC nailfold capillaroscopy features"}
    for c, (inst, item) in SURVEY.items():
        cols[c] = {"role": "survey", "instrument": inst, "item": item}
    for c, (loinc, unit) in NG.items():
        cols[c] = {"role": "lab", "loinc": loinc, "unit": unit, "platform": "nmr_nightingale"}
    return cols


SUBGROUP = {"name": "capillaroscopy_abnormal", "description": "SYNTHETIC: Long COVID cases with an abnormal nailfold "
                                                              "capillaroscopy pattern (reader call)",
            "within_label": "long_covid",
            "device_positive": {"block": "capillaroscopy", "feature": "device_positive", "threshold": 1,
                                "direction": ">="},
            "ehr_domains": ["condition", "measurement", "survey", "demographic"], "max_features": 8}


def stub_mapper(profile, dataset, condition, measurement, llm, out) -> dict:
    """Local stand-in for part A's mapper (used only when measure_it.agent.mapper is absent)."""
    t = profile["tables"][0]["name"]
    cols = reviewed_columns()
    return {"dataset_id": DS_ID, "condition_hint": condition, "measurement_hint": measurement,
            "approved_by": None, "approved_at": None,
            "tables": {t: {"participant_id": "record_id", "columns": cols,
                           "proposed_by": {c: {"method": "rule", "confidence": 1.0, "rationale": "test stub"}
                                           for c in cols}}},
            "comparator": {"type": "healthy", "description": "SYNTHETIC healthy volunteers"},
            "subgroup_analysis": SUBGROUP}


def human_edit(mapping: dict) -> dict:
    """The review step: dataset id, comparator, the key columns and the pre-specified subgroup plan."""
    mapping["dataset_id"] = DS_ID
    mapping["comparator"] = {"type": "healthy", "description": "SYNTHETIC healthy volunteers"}
    mapping["subgroup_analysis"] = SUBGROUP
    rc = reviewed_columns()
    for t in (mapping.get("tables") or {}).values():
        cols = t.setdefault("columns", {})
        if "record_id" in cols or "cohort" in cols:
            t["participant_id"] = "record_id"
            cols.update(rc)
    return mapping


def mapper_available() -> bool:
    try:
        return OR._module("mapper") is not None and any(
            callable(getattr(OR._module("mapper"), n, None)) for n in OR.MAPPER_FUNCS)
    except OR.AdapterMissing:
        return False


def prepare_run(inp: Path, out: Path, monkeypatch) -> str:
    """Profile + map (part A mapper or the stub), then the human edit of mapping.yaml. Returns which mapper ran."""
    which = "part A mapper" if mapper_available() else "test stub"
    if which == "test stub":
        monkeypatch.setitem(OR.PARTS, "mapper", stub_mapper)
    res = OR.run_launch(inp, "long_covid", "nailfold_capillaroscopy", out=out, stop_after="map", outreach=False,
                        echo=None)
    st = {s["step"]: s for s in res["plan"]["steps"]}
    assert st["profile"]["status"] == "done", st["profile"]
    assert st["map"]["status"] == "done", st["map"]
    mp = out / "mapping.yaml"
    m = yaml.safe_load(mp.read_text())
    mp.write_text(yaml.safe_dump(human_edit(m), sort_keys=False))
    return which


# ------------------------------------------------------------------------------------------------------------------
# unit tests (no processed data)
# ------------------------------------------------------------------------------------------------------------------
def test_data_use_problems_and_placeholders():
    assert OR.data_use_problems({"data_use": {"deidentified": True, "use_permitted": True,
                                              "statement": "IRB 2024-17 permits this local analysis"}}) == []
    probs = OR.data_use_problems({"data_use": {"deidentified": False, "use_permitted": False,
                                               "statement": "PLACEHOLDER: the data owner must confirm"}})
    assert len(probs) == 3 and any("deidentified" in p for p in probs) and any("statement" in p for p in probs)
    assert len(OR.data_use_problems({})) == 3


def test_prepare_mapping_records_its_edits():
    m = {"tables": {}, "subgroup_analysis": {"name": "x"}}
    edits = OR.prepare_mapping(m, synthetic=True, data_use=None, fast=True)
    assert m["synthetic"] and m["data_use"]["deidentified"] is True and "SYNTHETIC" in m["data_use"]["statement"]
    assert m["primary_analysis"]["n_permutations"] == OR.FAST_ANALYSIS["n_permutations"]
    assert m["subgroup_analysis"]["n_bootstrap"] == OR.FAST_ANALYSIS["n_bootstrap"]
    assert len(edits) == 3 and "data_use_confirmed_via" in m
    m2 = {"tables": {}}
    OR.prepare_mapping(m2, synthetic=False, data_use={"deidentified": True, "use_permitted": None,
                                                      "statement": "DUA 7"}, fast=False)
    assert m2["data_use"] == {"deidentified": True, "statement": "DUA 7"}


def test_mapping_summary_counts_methods_and_drops():
    m = stub_mapper({"tables": [{"name": "records"}]}, None, "long_covid", "nailfold_capillaroscopy", "none", None)
    s = OR.mapping_summary(m)
    assert s["roles"]["drop"] == 6 and s["roles"]["unmapped"] == 1 and s["methods"] == {"rule": s["n_columns"]}
    assert {d["column"] for d in s["dropped"]} >= {"zip_code", "initials", "notes", "enroll_date"}


def test_call_with_filters_and_names_missing_parameters():
    def f(profile, out=None, *, llm=None):
        return profile, out, llm
    assert OR.call_with(f, profile=1, out=2, llm=3, unused=4) == (1, 2, 3)
    with pytest.raises(TypeError, match="profile"):
        OR.call_with(f, out=2)


def test_keep_and_demo_rankings_are_exclusive(tmp_path):
    with pytest.raises(ValueError, match="keep"):
        OR.run_launch(tmp_path, "long_covid", "nailfold_capillaroscopy", out=tmp_path / "r", keep=True,
                      demo_rankings=True, echo=None)


def test_unapproved_mapping_stops_before_apply(tmp_path, monkeypatch):
    inp = make_redcap_input(tmp_path / "in")
    monkeypatch.setitem(OR.PARTS, "mapper", stub_mapper)
    applied = []
    monkeypatch.setitem(OR.PARTS, "apply", lambda *a, **k: applied.append(1))
    res = OR.run_launch(inp, "long_covid", "nailfold_capillaroscopy", out=tmp_path / "run", outreach=False,
                        echo=None)
    st = {s["step"]: s for s in res["plan"]["steps"]}
    assert st["approve"]["status"] == "stopped" and "launch apply" in st["approve"]["message"]
    assert st["apply"]["status"] == "skipped" and not applied
    assert st["byod_ingest"]["status"] == "skipped" and st["teardown"]["status"] == "skipped"
    assert st["report"]["status"] == "done" and res["plan"]["status"] == "stopped"
    assert (tmp_path / "run" / ".gitignore").read_text().strip().endswith("*")
    assert "The run stopped" in (tmp_path / "run" / "LAUNCH_REPORT.md").read_text()
    m = yaml.safe_load((tmp_path / "run" / "mapping.yaml").read_text())
    assert m["approved_by"] is None
    # the private-data rule: the profile never carries small categorical levels or identifier values
    prof = json.loads((tmp_path / "run" / "profile.json").read_text())
    for t in prof["tables"]:
        for c in t["columns"]:
            for lv in c.get("levels") or []:
                assert lv[0].startswith("<") or lv[1] >= 11
    assert "SYN-0001" not in json.dumps(prof)


def test_report_renders_a_minimal_result():
    r = {"generated_at": "t", "run": {"condition": "c", "measurement": "m", "out": "o", "input": "i",
                                      "status": "completed", "message": None},
         "steps": [{"step": "profile", "status": "done", "seconds": 1.0}],
         "input": {"profile": None, "mapping": None, "llm": "none"},
         "harmonization": {"byod_validation": None, "manifest": {}, "omop": None, "query_pack_check": None},
         "device_evidence": None, "subgroup": None, "similar": None, "launch": None, "outreach": None,
         "teardown": None, "caveats": ["x"]}
    md = RP.render(r)
    assert "# Launch report: c x m" in md and "no model was contacted" in md and "* x" in md


# ------------------------------------------------------------------------------------------------------------------
# end to end (needs data/processed)
# ------------------------------------------------------------------------------------------------------------------
def _snapshot() -> dict:
    from measure_it.byod import common as K
    from measure_it.store import read_table, table_exists
    snap = {}
    for t in (K.DATASETS_TABLE, K.RECORDS_TABLE, K.PHENOTYPES_TABLE, K.STRATA_TABLE, "participants",
              "geo_subgroup_burden", "similar_reference_prevalence", "similar_reference_profile",
              "similar_reference_distance"):
        snap[t] = len(read_table(t)) if table_exists(t) else 0
    snap["byod_partitions"] = sorted(p.name for p in PROCESSED.glob("*__byod_*"))
    d = read_table("deployment_opportunities", columns=["object_id", "rank_equal", "rank_evidence_weighted",
                                                         "measurement_performance_status"])
    snap["opp"] = d.sort_values("object_id").reset_index(drop=True)
    return snap


def _isolate(tmp_path, monkeypatch):
    from measure_it.byod import common as K
    from measure_it.byod import deploy as D
    from measure_it.similar import phenotype_io as PIO
    monkeypatch.setattr(K, "RESULTS_ROOT", tmp_path / "results_byod")
    monkeypatch.setattr(K, "LEDGER", tmp_path / "results_byod" / "LEDGER.jsonl")
    monkeypatch.setattr(D, "STATE", tmp_path / "results_byod" / "state.json")
    monkeypatch.setattr(PIO, "RESULTS_ROOT", tmp_path / "results_byod")
    monkeypatch.setenv(K.DEMO_ENV, "0")


def _check_teardown(before: dict, res: dict) -> None:
    from measure_it.byod import common as K
    ds = f"byod_{DS_ID}"
    assert not K.dataset_partitions(ds) and not K.raw_path(ds).exists() and not K.results_path(ds).exists()
    after = _snapshot()
    td = res["plan"]["teardown"]
    assert td["verified"] and td["partitions_left"] == [] and not any(td["rows_left"].values())
    if after["byod_partitions"] != before["byod_partitions"]:
        pytest.skip("another process ingested or removed a user dataset during this run: the processed tables are "
                    "shared, so only this run's own dataset could be checked (it is gone)")
    for k in before:
        if k == "opp":
            pd.testing.assert_frame_equal(before["opp"], after["opp"], check_dtype=False)
        else:
            assert before[k] == after[k], (k, before[k], after[k])


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_end_to_end_synthetic_launch(tmp_path, monkeypatch):
    """REDCap-style SYNTHETIC input -> profile/map/approve/apply/OMOP -> BYOD validate/ingest/evaluate/subgroup/
    similar -> subgroup burden -> launch ranking (engine + subgroup re-score) -> outreach -> teardown -> report."""
    _isolate(tmp_path, monkeypatch)
    inp = make_redcap_input(tmp_path / "maestro_redcap")
    out = tmp_path / "run"
    which = prepare_run(inp, out, monkeypatch)
    before = _snapshot()
    res = OR.run_launch(inp, "long_covid", "nailfold_capillaroscopy", out=out, llm=None, approve=True,
                        synthetic=True, fast=True, top=10, per_county=5, echo=None)
    st = {s["step"]: s for s in res["plan"]["steps"]}
    assert [s["step"] for s in res["plan"]["steps"]] == OR.STEPS
    for name in ("profile", "map", "approve", "apply", "omop_export", "byod_validate", "data_use_gate", "byod_ingest",
                 "byod_evaluate", "byod_subgroup", "byod_similar", "geo_subgroup_burden", "launch_ranking",
                 "outreach", "teardown", "report"):
        assert st[name]["status"] == "done", (name, st[name], which)
    # demo record: never in the default rankings
    assert st["byod_deploy"]["status"] == "skipped" and "demo" in st["byod_deploy"]["reason"]
    assert st["query_pack_check"]["status"] in ("done", "skipped")
    # BYOD validation passed; the manifest is synthetic + demo; OMOP export present
    val = json.loads((out / "byod_validation.json").read_text())
    assert val["ok"], val["errors"]
    man = yaml.safe_load((out / "byod" / "manifest.yaml").read_text())
    assert man["synthetic"] is True and man["demo"] is True and man["dataset_id"] == DS_ID
    for t in ("person", "observation", "concept"):
        assert (out / "omop" / f"{t}.parquet").exists(), t
    assert not any("zip" in p.name.lower() for p in (out / "byod").iterdir())
    # launch.json: top counties with every component, both rankings
    lj = json.loads((out / "launch.json").read_text())
    la = lj["launch"]
    assert la["primary"] == "subgroup" and la["engine"]["weight_set"] == "equal"
    eng, sub = la["engine"], la["subgroup"]
    assert len(eng["counties"]) == 10 and len(sub["counties"]) == 10
    for c in eng["counties"] + sub["counties"]:
        for k in ("burden_pct", "vulnerability_pct", "diagnostic_desert_pct", "clinic_capacity_pct",
                  "research_readiness_pct", "composite", "mc_rank_p05", "mc_rank_p95"):
            assert c[k] is not None, (k, c["geo_name"])
        assert {"billing_clinicians_n", "experience_sites_n", "candidate_facilities_n"} <= set(c["contacts"])
    assert eng["clinic_capacity_basis"] == "implementer_density"
    assert sub["replaced_members"] == ["long_covid"] and sub["counties"][0]["subgroup_estimate"] > 0
    assert -1 <= sub["effect"]["spearman_vs_engine"] <= 1
    # the engine ranking reproduces deployment_opportunities exactly
    assert [c["rank"] for c in eng["counties"]] == list(range(1, 11))
    assert lj["subgroup"]["phenotype_id"] == f"byod_{DS_ID}:capillaroscopy_abnormal"
    assert lj["device_evidence"]["demo"] in (True, 1)
    caveats = " ".join(lj["caveats"])
    for word in ("SYNTHETIC", "healthy controls", "Medicare fee-for-service", "small-area", "billing code"):
        assert word in caveats, word
    # outreach workbook under the run folder only
    wb = Path(lj["outreach"]["path"])
    assert wb.exists() and out in wb.parents
    # report
    md = (out / "LAUNCH_REPORT.md").read_text()
    for h in ("## 1. What went in", "## 2. Harmonization", "## 3. The device evidence", "## 4. The device subgroup",
              "## 5. Launch candidates", "## 6. Who to contact", "## 7. Engine state", "## 8. Caveats"):
        assert h in md
    assert "SYN-0" not in md and "SYN-0" not in json.dumps(lj)
    assert (out / "engine_results" / "computable_phenotype.json").exists()
    _check_teardown(before, res)
    print(json.dumps({"mapper": which,
                      "engine_top5": [(c["rank"], c["geo_name"]) for c in eng["counties"][:5]],
                      "subgroup_top5": [(c["rank"], c["geo_name"], c["rank_engine_same_weights"])
                                        for c in sub["counties"][:5]],
                      "effect": sub["effect"]}, indent=1))


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_end_to_end_demo_rankings_then_restored(tmp_path, monkeypatch):
    """--demo-rankings: the demo record enters the rankings for this run only (evidence_weighted), and the teardown
    restores deployment_opportunities exactly."""
    _isolate(tmp_path, monkeypatch)
    inp = make_redcap_input(tmp_path / "maestro_redcap")
    out = tmp_path / "run"
    prepare_run(inp, out, monkeypatch)
    before = _snapshot()
    res = OR.run_launch(inp, "long_covid", "nailfold_capillaroscopy", out=out, approve=True, synthetic=True,
                        fast=True, demo_rankings=True, outreach=False, top=5, echo=None)
    st = {s["step"]: s for s in res["plan"]["steps"]}
    assert st["byod_deploy"]["status"] == "done", st["byod_deploy"]
    la = res["launch"]
    assert la["engine"]["weight_set"] == "evidence_weighted" and la["engine"]["performance"]["status"] == "known"
    assert la["subgroup"]["weight_set"] == "evidence_weighted"
    assert all(c["measurement_evidence"] is not None for c in la["engine"]["counties"])
    _check_teardown(before, res)
    from measure_it.byod import deploy as D
    assert not D.demo_state()


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_private_data_without_data_use_stops_before_ingest(tmp_path, monkeypatch):
    _isolate(tmp_path, monkeypatch)
    inp = make_redcap_input(tmp_path / "maestro_redcap")
    out = tmp_path / "run"
    prepare_run(inp, out, monkeypatch)
    before = _snapshot()
    res = OR.run_launch(inp, "long_covid", "nailfold_capillaroscopy", out=out, approve=True, outreach=False,
                        echo=None)
    st = {s["step"]: s for s in res["plan"]["steps"]}
    assert st["data_use_gate"]["status"] == "stopped"
    msg = st["data_use_gate"]["message"]
    assert "data_use.deidentified" in msg and "data_use.use_permitted" in msg and "mapping.yaml" in msg
    assert st["byod_ingest"]["status"] == "skipped" and st["teardown"]["status"] == "skipped"
    assert res["plan"]["status"] == "stopped"
    now = _snapshot()
    for k in before:
        if k != "opp":
            assert before[k] == now[k], k


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_burden_only_fallback_outside_the_grid():
    from measure_it.agent import launch as L
    r = L.launch_ranking("ptlds", "nailfold_capillaroscopy", top=5, log=lambda *_: None)
    assert not r["scored_grid"] and r["primary"] == "engine" and r["engine"]["kind"] == "burden_only"
    assert "BURDEN-ONLY" in r["engine"]["why"] and len(r["engine"]["counties"]) == 5
    assert any("BURDEN-ONLY" in c for c in r["caveats"])


def test_device_positive_option_builds_a_locked_plan():
    """`--device-positive COLUMN[>=X]` -> subgroup_analysis for a device-role column; anything else is refused."""
    mapping = {"tables": {"records": {"columns": {
        "cap_positive": {"role": "device", "block": "capillaroscopy"},
        "cap_abnormality_score": {"role": "device", "block": "capillaroscopy"},
        "age": {"role": "age"}}}}}
    flag = OR.device_positive_plan(mapping, "cap_positive", "long_covid")
    assert flag["within_label"] == "long_covid"
    assert flag["device_positive"] == {"block": "device_capillaroscopy", "feature": "cap_positive",
                                       "threshold": 1.0, "direction": ">="}
    thr = OR.device_positive_plan(mapping, "cap_abnormality_score > 2.5", "long_covid")
    assert thr["device_positive"]["threshold"] == 2.5 and thr["device_positive"]["direction"] == ">"
    for bad in ("age", "nope", "cap_positive >= x"):
        with pytest.raises(ValueError):
            OR.device_positive_plan(mapping, bad, "long_covid")
