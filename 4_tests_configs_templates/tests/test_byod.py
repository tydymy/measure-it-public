"""Tests for bring-your-own-data (measure_it.byod): refusals, namespacing, provenance, plan locking, demo exclusion, the
no-geography guarantee and an end-to-end run of the SYNTHETIC example (ingest -> evaluate -> deploy -> remove).

Pure tests run on a fresh clone (the synthetic example is generated into a temporary folder and its label condition
is given as a canonical id, so no processed ontology table is needed). Tests marked `data` need data/processed; the
end-to-end test writes and then removes its own dataset (`byod_pytest_synth`) and restores the rankings.
"""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from typer.testing import CliRunner

from measure_it.byod import checks as C
from measure_it.byod import common as K
from measure_it.config import PROCESSED, PROJECT_ROOT, UNKNOWN

ROOT = PROJECT_ROOT
HAS_DATA = (PROCESSED / "deployment_opportunities.parquet").exists()
_spec = importlib.util.spec_from_file_location("byod_make_example", ROOT / "examples" / "byod_synthetic" / "make_example.py")
MAKE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(MAKE)


def _example(tmp_path: Path, dataset_id: str = "pytest_synth", n_cases: int = 40, n_controls: int = 50,
             condition: str = "me_cfs") -> Path:
    out = MAKE.make(tmp_path / dataset_id, n_cases=n_cases, n_controls=n_controls, fast=True, dataset_id=dataset_id)
    m = yaml.safe_load((out / "manifest.yaml").read_text())
    m["labels"][0]["condition"] = condition           # canonical id: resolves without the ontology tables
    (out / "manifest.yaml").write_text(yaml.safe_dump(m, sort_keys=False))
    return out


def _edit_csv(path: Path, fn) -> None:
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    fn(df).to_csv(path, index=False)


def _edit_manifest(path: Path, fn) -> None:
    m = yaml.safe_load((path / "manifest.yaml").read_text())
    fn(m)
    (path / "manifest.yaml").write_text(yaml.safe_dump(m, sort_keys=False))


# ------------------------------------------------------------------------------------------------ contract

def test_schema_template_matches_code():
    shipped = json.loads((ROOT / "templates" / "byod" / "manifest.schema.json").read_text())
    assert shipped == json.loads(json.dumps(C.MANIFEST_SCHEMA))
    assert (ROOT / "templates" / "byod" / "manifest.yaml").exists()


def test_synthetic_example_validates(tmp_path):
    rep = C.check_dir(_example(tmp_path))
    assert rep.ok, rep.errors
    assert rep.resolved["labels"]["me_cfs"]["condition_id"] == "me_cfs"
    assert rep.resolved["group_sizes"] == {"cases": 40, "controls": 50}
    assert any("participant_wearable_features__byod_pytest_synth" in w for w in rep.will_create)


def test_synthetic_example_is_marked_synthetic_everywhere(tmp_path):
    out = _example(tmp_path)
    m = yaml.safe_load((out / "manifest.yaml").read_text())
    assert m["synthetic"] is True and m["demo"] is True and "SYNTHETIC" in m["title"]
    assert "SYNTHETIC" in m["labels"][0]["definition"] and "SYNTHETIC" in m["comparator"]["description"]
    for f in out.glob("*.csv"):
        ids = pd.read_csv(f, dtype=str)["participant_id"]
        assert ids.str.startswith("SYNTHETIC-").all(), f.name
    assert (out / "README_SYNTHETIC.md").exists()


# ------------------------------------------------------------------------------------------------ refusals

@pytest.mark.parametrize("col", ["zip", "ZIP5", "postcode", "zip_code", "county", "county_fips", "latitude", "lon",
                                 "home_address", "street", "city", "state", "first_name", "patient_name", "name",
                                 "dob", "date_of_birth", "birth_date", "mrn", "phone", "email", "e_mail", "ssn",
                                 "visit_date", "admission_date", "age", "age_years", "gps_lat", "zcta"])
def test_identifying_or_geographic_column_names_are_refused(col):
    assert C.column_problem(col) is not None, col


@pytest.mark.parametrize("col", ["rmssd_ms", "steps_mean_daily", "sleep_state_fraction", "lab_name", "hip_accel_counts",
                                 "day_index", "age_band", "sex", "loinc", "PROT001", "upright_minutes"])
def test_ordinary_feature_names_pass(col):
    assert C.column_problem(col) is None, col


def test_identifying_values_are_detected_and_never_echoed():
    s = pd.Series(["ok", "jane.doe@example.org", "555-123-4567", "123-45-6789", "2024-03-15", "02139",
                   "12 Main Street", "34.0522, -118.2437", "3.4"])
    probs = C.value_problems(s, "notes", "x.csv")
    text = " ".join(probs)
    for kind in ("e-mail", "phone", "social security", "exact date", "ZIP", "street address", "coordinate"):
        assert kind in text, kind
    for secret in ("jane.doe", "555-123", "123-45-6789", "Main Street", "118.2437"):
        assert secret not in text          # the offending value is never printed
    assert C.value_problems(pd.Series(["12.5", "3", "14000"]), "steps", "x.csv") == []   # numbers are not scanned


def test_mrn_like_participant_ids_are_refused():
    for bad in (["1234567", "7654321"], ["MRN001234", "MRN001235"], ["123-45-6789"], ["John Smith"]):
        assert C.id_problems(pd.Series(bad), "participants.csv"), bad
    assert C.id_problems(pd.Series(["SYN-001", "P0002", "site3_17"]), "participants.csv") == []


def test_age_bands():
    assert C.age_band_problem("40-49") is None and C.age_band_problem("90+") is None
    assert C.age_band_problem("91-95") and C.age_band_problem("45") and C.age_band_problem("40-41")
    assert C.age_band_problem("95+")


@pytest.mark.parametrize("case", ["zip_column", "dob_column", "mrn_ids", "exact_dates", "email_in_lab_name",
                                  "latlon_device", "not_deidentified", "too_few_cases", "unknown_file",
                                  "synthetic_not_demo", "undeclared_participant_column", "orphan_ids",
                                  "bad_loinc", "hashing_without_salt"])
def test_validation_refuses(tmp_path, case, monkeypatch):
    out = _example(tmp_path)
    monkeypatch.delenv(K.SALT_ENV, raising=False)
    expect = {
        "zip_column": "geographic", "dob_column": "identifies a person", "mrn_ids": "MRN",
        "exact_dates": "calendar date", "email_in_lab_name": "e-mail", "latlon_device": "geographic",
        "not_deidentified": "de-identified", "too_few_cases": "at least 10 cases", "unknown_file": "not a file name",
        "synthetic_not_demo": "demo: true", "undeclared_participant_column": "undeclared column",
        "orphan_ids": "not in participants.csv", "bad_loinc": "LOINC", "hashing_without_salt": K.SALT_ENV}[case]
    P = out / "participants.csv"
    if case == "zip_column":
        _edit_csv(P, lambda d: d.assign(zip="02139"))
    elif case == "dob_column":
        _edit_csv(P, lambda d: d.assign(date_of_birth="1980"))
    elif case == "mrn_ids":
        for f in out.glob("*.csv"):
            _edit_csv(f, lambda d: d.assign(participant_id=d["participant_id"].str.replace("SYNTHETIC-", "88800",
                                                                                          regex=False)))
    elif case == "exact_dates":
        _edit_csv(out / "device_synthetic_wrist_actigraphy.csv", lambda d: d.assign(visit_date="2024-01-02"))
    elif case == "email_in_lab_name":
        _edit_csv(out / "ehr_labs.csv", lambda d: d.assign(lab_name="ordered by dr.who@example.org"))
    elif case == "latlon_device":
        _edit_csv(out / "device_synthetic_wrist_actigraphy.csv", lambda d: d.assign(latitude="34.05", longitude="-118.2"))
    elif case == "not_deidentified":
        _edit_manifest(out, lambda m: m["data_use"].update(deidentified=False))
    elif case == "too_few_cases":
        _edit_csv(P, lambda d: d.assign(me_cfs=np.where((d["me_cfs"] == "1").cumsum() > 8, "", d["me_cfs"])))
    elif case == "unknown_file":
        (out / "notes.xlsx").write_bytes(b"x")
    elif case == "synthetic_not_demo":
        _edit_manifest(out, lambda m: m.update(demo=False))
    elif case == "undeclared_participant_column":
        _edit_csv(P, lambda d: d.assign(bmi="24.1"))
    elif case == "orphan_ids":
        _edit_csv(out / "ehr_diagnoses.csv", lambda d: pd.concat([d, pd.DataFrame({"participant_id": ["SYNTHETIC-9999"],
                                                                                    "icd10cm": ["I10"]})]))
    elif case == "bad_loinc":
        _edit_csv(out / "ehr_labs.csv", lambda d: d.assign(loinc="CRP"))
    elif case == "hashing_without_salt":
        _edit_manifest(out, lambda m: m.update(id_hashing="salted_sha256"))
    rep = C.check_dir(out)
    assert not rep.ok
    assert any(expect.lower() in e.lower() for e in rep.errors), rep.errors


def test_cli_validate_exit_codes(tmp_path):
    from measure_it.cli import app
    out = _example(tmp_path)
    r = CliRunner().invoke(app, ["byod", "validate", str(out)])
    assert r.exit_code == 0 and "OK" in r.output
    _edit_csv(out / "participants.csv", lambda d: d.assign(postcode="SW1A 1AA"))
    r = CliRunner().invoke(app, ["byod", "validate", str(out)])
    assert r.exit_code == 1 and "REFUSED" in r.output and "postcode" in r.output
    assert "SW1A" not in r.output


# ------------------------------------------------------------------------------------------------ namespacing

@pytest.mark.data  # maps ICD codes through condition_registry (data/processed)
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_namespacing_and_salted_hashing(tmp_path, monkeypatch):
    from measure_it.byod.ingest import build_frames
    out = _example(tmp_path)
    rep = C.check_dir(out)
    fr = build_frames(rep)
    pid = fr["tables"]["participants"]["participant_id"]
    assert pid.str.startswith("byod_pytest_synth:SYNTHETIC-").all() and pid.is_unique
    assert fr["tables"]["participant_conditions"]["participant_id"].str.startswith("byod_pytest_synth:").all()
    # salted hashing: the partner's ids disappear, the mapping is stable for one salt and changes with the salt
    monkeypatch.setenv(K.SALT_ENV, "a-long-enough-test-salt-0001")
    _edit_manifest(out, lambda m: m.update(id_hashing="salted_sha256"))
    rep = C.check_dir(out)
    assert rep.ok, rep.errors
    a = build_frames(rep)["tables"]
    assert a["participants"]["participant_id"].str.fullmatch(r"byod_pytest_synth:h[0-9a-f]{16}").all()
    for t in a.values():
        assert not t.astype(str).apply(lambda s: s.str.contains("SYNTHETIC-", regex=False)).any().any()
    b = build_frames(C.check_dir(out))["tables"]["participants"]["participant_id"]
    assert (a["participants"]["participant_id"] == b).all()
    monkeypatch.setenv(K.SALT_ENV, "another-long-enough-salt-0002")
    c = build_frames(C.check_dir(out))["tables"]["participants"]["participant_id"]
    assert not (a["participants"]["participant_id"] == c).any()


def test_engine_ids():
    assert K.engine_id("abc") == "byod_abc" and K.engine_id("byod_abc") == "byod_abc"
    assert K.manifest_id_of("byod_abc") == "abc"
    assert K.partition_name("participants", "abc") == "participants__byod_abc"


# ------------------------------------------------------------------------------------------------ plan locking

def test_plan_lock_refuses_a_changed_plan(tmp_path, monkeypatch):
    from measure_it.byod import evaluate as E
    monkeypatch.setattr(K, "RESULTS_ROOT", tmp_path / "results_byod")
    monkeypatch.setattr(K, "LEDGER", tmp_path / "results_byod" / "LEDGER.jsonl")
    out = _example(tmp_path)
    m = yaml.safe_load((out / "manifest.yaml").read_text())
    meta = {"dataset_id": "byod_pytest_synth",
            "labels_json": json.dumps({"me_cfs": {"condition_id": "me_cfs", "definition": "d", "label_basis": "x"}})}
    a = E.lock_plan(m, meta)
    assert a["status"] == "locked" and re.fullmatch(r"[0-9a-f]{64}", a["sha256"])
    assert (tmp_path / "results_byod" / "pytest_synth" / "PLAN.md").exists()
    assert (tmp_path / "results_byod" / "pytest_synth" / ".gitignore").read_text().strip().endswith("*")
    assert E.lock_plan(m, meta)["status"] == "verified"
    m2 = json.loads(json.dumps(m))
    m2["primary_analysis"]["specificity"] = 0.95
    with pytest.raises(RuntimeError, match="locked"):
        E.lock_plan(m2, meta)
    b = E.lock_plan(m2, meta, amend=True)
    assert b["status"] == "amended" and b["amends"] == a["sha256"]
    ev = [e["event"] for e in E.ledger_events("byod_pytest_synth")]
    assert ev == ["plan_locked", "plan_amended"]


def test_plan_hash_does_not_depend_on_outcomes():
    from measure_it.byod.evaluate import plan_document
    m = {"primary_analysis": {"label": "y", "feature_blocks": ["device_a"]}, "comparator": {"type": "healthy"},
         "measurement": {"class": "hrv"}}
    meta = {"dataset_id": "byod_x", "labels_json": json.dumps({"y": {"condition_id": "me_cfs", "definition": "d",
                                                                      "label_basis": "proxy"}})}
    doc = plan_document(m, meta)          # built from the manifest and the ingest metadata only (no data read)
    assert set(doc) == {"dataset_id", "engine_version", "seed", "label", "condition_id", "label_definition",
                        "label_basis", "comparator", "measurement_class", "feature_blocks", "combined_blocks",
                        "covariates", "specificity", "metric_link_specificity", "cv", "model", "n_bootstrap",
                        "n_permutations_primary", "n_permutations_secondary", "secondary", "label_defining_exclusion",
                        "decision_rule", "record"}
    assert doc == plan_document(m, meta)
    assert doc["specificity"] == 0.90 and doc["cv"]["n_repeats"] == 10


# ------------------------------------------------------------------------------------------------ tier + demo

def test_user_tier_sits_between_own_and_published():
    from measure_it.scoring import metric_link as ML
    assert ML.tier_value(2.5) == 2.5 and ML.tier_value(1.0) == 1 and ML.tier_value(None) == 4
    assert ML.TIER_LABEL[2.5] == K.SOURCE_KIND and "2.5" in ML.tier_text(2.5)
    assert ML.TIER_FACTOR[2] > ML.TIER_FACTOR[2.5] > ML.TIER_FACTOR[3]
    assert ML.evidence_score(0.62, 2.5)[0] == pytest.approx(0.6 * 0.24)


def _recs_frame(ML, tiers):
    recs = [ML._rec(record_id=f"R{t}", source_kind="x", quality_tier=t, target_condition="me_cfs",
                    primary_class="hrv", auroc=0.75, auroc_ci_low=0.65, auroc_ci_high=0.85, op_sensitivity=0.4,
                    op_specificity=0.9, n_cases=50.0, n_controls=50.0) for t in tiers]
    df = pd.DataFrame(recs)
    for c in ("target_members", "source_object_ids", "source_tables", "caveats"):
        df[c] = df[c].map(json.dumps)
    return df


def test_selection_prefers_own_public_then_user_then_published():
    from measure_it.scoring import metric_link as ML
    assert ML.attach("me_cfs", "hrv", _recs_frame(ML, [3, 2.5, 1]))["selected_record_id"] == "R1"
    row = ML.attach("me_cfs", "hrv", _recs_frame(ML, [3, 2.5]))
    assert row["selected_record_id"] == "R2.5" and row["quality_tier"] == 2.5
    assert row["quality_tier_label"] == K.SOURCE_KIND
    assert ML.attach("me_cfs", "hrv", _recs_frame(ML, [3]))["selected_record_id"] == "R3"


def test_demo_records_excluded_unless_enabled(monkeypatch):
    from measure_it.byod import records as R
    rows = []
    for rid, demo in (("BYOD-A-PRIMARY", False), ("BYOD-DEMO-PRIMARY", True)):
        rows.append({"record_id": rid, "dataset_id": "byod_" + rid, "demo": demo, "target_condition": "me_cfs",
                     "primary_class": "hrv", "auroc": 0.7, "auroc_ci_low": 0.6, "auroc_ci_high": 0.8,
                     "measurement_only": True, "caveats": "[]", "source_object_ids": "[]", "source_tables": "[]",
                     "target_members": '["me_cfs"]', "user_supplied": True, "synthetic": demo, "plan_sha256": "x"})
    table = pd.DataFrame(rows)
    monkeypatch.setattr(R, "table_exists", lambda name: True)
    monkeypatch.setattr(R, "read_table", lambda name: table.copy())
    monkeypatch.setenv(K.DEMO_ENV, "0")
    assert [r["record_id"] for r in R.performance_records()] == ["BYOD-A-PRIMARY"]
    monkeypatch.setenv(K.DEMO_ENV, "1")
    assert {r["record_id"] for r in R.performance_records()} == {"BYOD-A-PRIMARY", "BYOD-DEMO-PRIMARY"}
    assert all(r["quality_tier"] == 2.5 for r in R.performance_records())


# ------------------------------------------------------------------------------------------------ no geography

READ_RE = re.compile(r"read_table\(\s*[\"']([A-Za-z_0-9]+)")


def test_scoring_code_reads_no_person_table():
    """The ranking side (scoring, facilities, geography) never reads a person-layer or BYOD person table: a user
    dataset reaches it only through the aggregate measurement_performance record."""
    person = ("participant", "participants", "phenotype_signatures", "byod_datasets")
    for pkg in ("scoring", "facilities", "geography"):
        for f in (ROOT / "src" / "measure_it" / pkg).glob("*.py"):
            names = READ_RE.findall(f.read_text())
            bad = [n for n in names if n.startswith(person)]
            assert not bad, (f.name, bad)
    rec_src = (ROOT / "src" / "measure_it" / "byod" / "records.py").read_text()
    assert set(READ_RE.findall(rec_src)) <= set()        # reads only via K.RECORDS_TABLE
    assert "read_table(K.RECORDS_TABLE)" in rec_src


def test_record_fields_are_aggregate_only():
    from measure_it.byod.records import EXTRA_FIELDS, RECORD_FIELDS
    fields = set(RECORD_FIELDS) | set(EXTRA_FIELDS)
    assert not {f for f in fields if "participant" in f or C.column_problem(f)}


# ------------------------------------------------------------------------------------------------ end to end

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_end_to_end_synthetic(tmp_path, monkeypatch):
    """ingest -> evaluate -> deploy (demo excluded) -> deploy --demo -> tools -> remove, on the SYNTHETIC example."""
    from measure_it.byod import deploy as D
    from measure_it.byod import tools as BT
    from measure_it.byod.evaluate import evaluate
    from measure_it.byod.ingest import ingest
    from measure_it.byod.remove import remove
    from measure_it.provenance import PROVENANCE_COLUMNS, check_provenance
    from measure_it.scoring import metric_link as ML
    from measure_it.store import read_table, table_exists
    from measure_it.validate import GEO_COLUMN_RE, check_person_geography
    monkeypatch.setattr(K, "RESULTS_ROOT", tmp_path / "results_byod")
    monkeypatch.setattr(K, "LEDGER", tmp_path / "results_byod" / "LEDGER.jsonl")
    monkeypatch.setattr(D, "STATE", tmp_path / "results_byod" / "state.json")
    monkeypatch.setenv(K.DEMO_ENV, "0")
    ds = "byod_pytest_synth"
    out = _example(tmp_path, condition="ME/CFS")
    before = read_table(ML.TABLE)
    cell = lambda t: t[(t["condition_id"] == "me_cfs") & (t["measurement_id"] == "autonomic_function_testing")]  # noqa
    assert cell(before)["performance_status"].iloc[0] == UNKNOWN, "baseline changed: rerun the pipeline"
    try:
        res = ingest(out, echo=lambda *_: None)
        # provenance + person layer + namespacing + no geography, on every partition written
        for name in res["partitions"]:
            t = read_table(name)
            check_provenance(t, name)
            assert set(t["data_layer"]) == {"person"} and t["evidence_type"].str.startswith("person_").all()
            assert t["provenance_notes"].str.startswith(K.PROVENANCE_NOTE).all()
            assert set(t["source_geographic_resolution"]) == {"none"}
            assert t["participant_id"].astype(str).str.startswith(f"{ds}:").all()
            assert not [c for c in t.columns if GEO_COLUMN_RE.search(c) and c not in PROVENANCE_COLUMNS]
        assert (K.raw_path(ds) / "DATA_AUDIT.md").exists() and (K.raw_path(ds) / ".gitignore").exists()
        entry = yaml.safe_load((K.raw_path(ds) / "byod_registry_entry.yaml").read_text())
        assert entry["user_supplied"] and entry["local_only"] and not (K.raw_path(ds) / "registry_entry.yaml").exists()
        dpv = read_table("digital_phenotype_datasets")
        assert (dpv["dataset_id"] == ds).sum() == 1                       # its own representation, not pooled
        assert check_person_geography()["ok"]
        ev = evaluate(out, echo=lambda *_: None)
        rec = ev["record"]
        assert rec["record_id"] == "BYOD-PYTEST_SYNTH-PRIMARY" and rec["demo"] and rec["measurement_only"]
        assert rec["quality_tier"] == 2.5 and rec["target_condition"] == "me_cfs" and rec["primary_class"] == "hrv"
        assert any("healthy controls" in c.lower() and "optimistic" in c.lower() for c in json.loads(rec["caveats"]))
        assert not any(f"{ds}:" in str(v) for v in rec.values())          # aggregate only
        sig = read_table(f"phenotype_signatures__{ds}")
        assert sig["object_id"].is_unique and sig["object_id"].str.startswith(f"signature:{ds}|").all()
        assert (sig["feature"] == "model_auroc:primary").sum() == 1
        # the label-defining ICD-10-CM code (G93.32 -> me_cfs) never enters a model
        feats = pd.read_csv(K.results_path(ds) / "evaluation_features.csv") if (
            K.results_path(ds) / "evaluation_features.csv").exists() else pd.DataFrame(columns=["analyte"])
        assert not feats["analyte"].astype(str).str.contains("G93").any()
        # deploy without --demo: the demo record stays out, nothing changes
        r0 = D.deploy(out, demo=False, echo=lambda *_: None)
        assert r0["changed_scored_combinations"] == []
        assert cell(read_table(ML.TABLE))["performance_status"].iloc[0] == UNKNOWN
        # deploy --demo: the (me_cfs, autonomic_function_testing) cell becomes known and is ranked
        r1 = D.deploy(out, demo=True, echo=lambda *_: None)
        after = cell(read_table(ML.TABLE)).iloc[0]
        assert after["selected_record_id"] == "BYOD-PYTEST_SYNTH-PRIMARY" and after["performance_status"] == "known"
        ch = [c for c in r1["ranking_changes"] if c["geo_level"] == "state"][0]
        assert ch["n_ranked_evidence_weighted_before"] == 0 and ch["n_ranked_evidence_weighted_after"] > 0
        assert ch["equal_weight_ranks_unchanged"]
        assert r1["guard"]["participant_ids_in_derived_tables"] == []
        w = [c for c in r1["relevant_cells"] if c["measurement_id"] == "wearable_autonomic_activity_monitoring"][0]
        assert not w["uses_user_record"]                                  # tier 1 MUSCLE-ME steps stays selected
        assert (K.results_path(ds) / "heatmap_measurement_evidence.csv").exists()
        assert (K.results_path(ds) / "embedding_dpv_pca.csv").exists()
        env = BT.list_user_datasets()
        assert env["status"] == "ok" and ds in [d["dataset_id"] for d in env["data"]["datasets"]]
        m = BT.get_user_dataset_metric("pytest_synth")
        assert m["status"] == "ok" and m["data"]["performance_record"]["record_id"] == "BYOD-PYTEST_SYNTH-PRIMARY"
        assert f"{ds}:" not in json.dumps(m)
    finally:
        remove(ds, echo=lambda *_: None)
        monkeypatch.setenv(K.DEMO_ENV, "0")
    assert not K.dataset_partitions(ds) and not K.raw_path(ds).exists()
    for t in (K.RECORDS_TABLE, K.DATASETS_TABLE):
        assert not table_exists(t) or ds not in set(read_table(t)["dataset_id"])
    restored = read_table(ML.TABLE)
    cols = ["condition_id", "measurement_id", "performance_status", "selected_record_id", "measurement_evidence"]
    pd.testing.assert_frame_equal(before[cols], restored[cols], check_dtype=False)
    assert not read_table("participants")["dataset_id"].astype(str).eq(ds).any()
    assert check_person_geography()["ok"]


# ------------------------------------------------------------------------------------------------ dashboard + MCP

def test_page6_source_has_limitations_panel_and_product_language():
    from measure_it.validate import forbidden_hits
    src = (ROOT / "dashboard" / "pages" / "6_Your_Data.py").read_text()
    assert "limitations_panel(" in src and not forbidden_hits(src)
    assert "pages/6_Your_Data.py" in (ROOT / "dashboard" / "app.py").read_text()


def test_page6_renders_without_user_data(monkeypatch):
    from streamlit.testing.v1 import AppTest

    from measure_it.byod import tools as BT
    empty = {"tool": "list_user_datasets", "status": UNKNOWN, "data": {"datasets": []}, "caveats": [],
             "reason": "no user dataset has been ingested (measure-it byod ingest <dir>)", "provenance": {},
             "truncation": [], "query": {}}
    monkeypatch.setattr(BT, "list_user_datasets", lambda: empty)
    at = AppTest.from_file(str(ROOT / "dashboard" / "pages" / "6_Your_Data.py"), default_timeout=120).run()
    assert not at.exception
    assert [t.value for t in at.title] == ["Your data"]
    assert any("make_example.py" in m.value for m in at.markdown)


def test_mcp_registers_user_dataset_tools_outside_the_spec_list():
    import asyncio

    from measure_it import tools as T
    from measure_it.mcp.server import BYOD_TOOLS, build_server
    names = {t.name for t in asyncio.run(build_server().list_tools())}
    assert set(BYOD_TOOLS) <= names and set(T.SPEC_TOOLS) <= names
    assert not set(BYOD_TOOLS) & set(T.TOOLS)          # the API tool list and its count are unchanged


# ------------------------------------------------------------------------------------------------ MAESTRO-like layers,
# ehr_medications, surveys, self-reported history, subgroup analysis (stage B+C)

def _maestro(tmp_path: Path, dataset_id: str = "pytest_maestro", cohorts: dict | None = None) -> Path:
    return MAKE.make_maestro(tmp_path / dataset_id, cohorts=cohorts, fast=True, dataset_id=dataset_id)


def test_ehr_example_has_medications(tmp_path):
    out = _example(tmp_path)
    rx = pd.read_csv(out / "ehr_medications.csv", dtype=str)
    assert set(rx.columns) == {"participant_id", "rxnorm"} and rx["rxnorm"].str.fullmatch(r"\d+").all()
    rep = C.check_dir(out)
    assert rep.ok, rep.errors
    assert "ehr_medications" in rep.info["blocks"]
    assert any("participant_medications__byod_pytest_synth" in w for w in rep.will_create)


def test_maestro_example_validates_without_ehr_files(tmp_path):
    out = _maestro(tmp_path)
    assert not any((out / f).exists() for f in ("ehr_labs.csv", "ehr_diagnoses.csv", "ehr_medications.csv"))
    m = yaml.safe_load((out / "manifest.yaml").read_text())
    assert m["synthetic"] is True and m["demo"] is True and "SYNTHETIC" in m["title"]
    assert {lab["condition"] for lab in m["labels"]} == {"long_covid", "ptlds", "lyme_disease"}
    for f in out.glob("*.csv"):
        assert pd.read_csv(f, dtype=str)["participant_id"].str.startswith("SYNTHETIC-M").all(), f.name
    rep = C.check_dir(out)
    assert rep.ok, rep.errors
    assert {"survey_compass31", "survey_fss9", "self_report_history", "omics_nightingale", "device_ring",
            "device_capillaroscopy"} <= set(rep.info["blocks"])
    assert any("unmapped" in w or "do not map" in w for w in rep.warnings)      # 'seasonal allergies' is counted
    assert any("person_concepts__byod_pytest_maestro" in w for w in rep.will_create)


@pytest.mark.parametrize("case", ["bad_rxnorm", "long_history_text", "bad_history_icd", "survey_text",
                                  "subgroup_unknown_label", "subgroup_unknown_feature", "subgroup_both_or_none",
                                  "dp_column_values", "omics_section_without_file"])
def test_validation_refuses_new_kinds(tmp_path, case):
    out = _maestro(tmp_path) if case != "bad_rxnorm" else _example(tmp_path)
    expect = {"bad_rxnorm": "RxNorm", "long_history_text": "80 characters", "bad_history_icd": "ICD-10-CM",
              "survey_text": "non-numeric", "subgroup_unknown_label": "within_label",
              "subgroup_unknown_feature": "device_positive.feature", "subgroup_both_or_none": "device_positive",
              "dp_column_values": "device_positive column", "omics_section_without_file": "omics.metabolon"}[case]
    if case == "bad_rxnorm":
        _edit_csv(out / "ehr_medications.csv", lambda d: d.assign(rxnorm="metoprolol"))
    elif case == "long_history_text":
        _edit_csv(out / "self_report_history.csv", lambda d: d.assign(condition="x" * 120))
    elif case == "bad_history_icd":
        _edit_csv(out / "self_report_history.csv", lambda d: d.assign(icd10cm="Raynaud"))
    elif case == "survey_text":
        _edit_csv(out / "survey_fss9.csv", lambda d: d.assign(total="quite tired"))
    elif case == "subgroup_unknown_label":
        _edit_manifest(out, lambda m: m["subgroup_analysis"].update(within_label="me_cfs"))
    elif case == "subgroup_unknown_feature":
        _edit_manifest(out, lambda m: m["subgroup_analysis"]["device_positive"].update(feature="no_such_score"))
    elif case == "subgroup_both_or_none":
        _edit_manifest(out, lambda m: m["subgroup_analysis"]["device_positive"].update(column="capi_pos"))
    elif case == "dp_column_values":
        _edit_manifest(out, lambda m: m["subgroup_analysis"].update(device_positive={"column": "capi_pos"}))
        _edit_csv(out / "participants.csv", lambda d: d.assign(capi_pos="maybe"))
    elif case == "omics_section_without_file":
        _edit_manifest(out, lambda m: m["omics"].update(metabolon={"platform": "other"}))
    rep = C.check_dir(out)
    assert not rep.ok
    assert any(expect.lower() in e.lower() for e in rep.errors), rep.errors


def test_device_positive_column_is_allowed(tmp_path):
    out = _maestro(tmp_path)
    _edit_manifest(out, lambda m: m["subgroup_analysis"].update(device_positive={"column": "capi_pos"}))
    _edit_csv(out / "participants.csv", lambda d: d.assign(capi_pos=np.where(d.index % 3 == 0, "1", "0")))
    rep = C.check_dir(out)
    assert rep.ok, rep.errors


@pytest.mark.data  # maps ICD codes through condition_registry (data/processed)
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_new_frames_maestro(tmp_path):
    from measure_it.byod.ingest import build_frames
    out = _maestro(tmp_path)
    fr = build_frames(C.check_dir(out))
    sv = fr["tables"]["participant_survey"]
    assert set(sv["instrument"]) == {"COMPASS31", "FSS9", "DSQPEM", "PROMIS_COG", "FUNCAP27"}
    assert sv["concept_code"].str.fullmatch(r"[A-Z0-9_]+:[a-z_]+").all()
    c = fr["tables"]["participant_conditions"]
    hx = c[c["condition_type"] == "self_report_history"]
    ray = hx[hx["reported_text"] == "Raynaud's phenomenon"]
    assert len(ray) and (ray["icd10cm_code"] == "I73.0").all() and (ray["mapping_confidence"] == "medium").all()
    assert (hx.loc[hx["reported_text"] == "seasonal allergies", "mapping_status"] == "unmapped").all()
    assert (hx.loc[hx["reported_text"] == "Long COVID", "canonical_condition_id"] == "long_covid").all()
    assert fr["blocks"]["survey_compass31"]["instrument"] == "COMPASS31"


def test_l1_model_respects_max_features_and_wilson():
    from measure_it.byod.subgroup import L1TopK, strata_fractions, wilson
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 30))
    y = (X[:, 0] + X[:, 1] + 0.5 * X[:, 2] + rng.normal(0, 1, 200) > 0).astype(int)
    X[rng.random(X.shape) < 0.05] = np.nan
    m = L1TopK(C=10.0, max_features=3).fit(X, y)
    assert 1 <= np.count_nonzero(m.model_.coef_) <= 3 and m.C_used_ < 10.0
    assert set(np.flatnonzero(m.model_.coef_[0])) <= {0, 1, 2, 3, 4, 5}
    lo, hi = wilson(14, 31)
    assert lo == pytest.approx(0.2916, abs=1e-3) and hi == pytest.approx(0.6223, abs=1e-3)   # by hand
    y_s = pd.Series([1, 0] * 10 + [1] * 5, index=[f"p{i}" for i in range(25)])
    demo = pd.DataFrame({"age_band": ["40-49"] * 20 + ["50-59"] * 5, "sex": ["female"] * 25}, index=y_s.index)
    rows = strata_fractions(y_s, demo)
    assert rows[0]["age_band"] == "all" and rows[0]["n_cases"] == 25 and rows[0]["n_positive"] == 15
    small = [r for r in rows if r["age_band"] == "50-59"][0]
    assert small["suppressed"] and small["fraction"] is None and small["n_positive"] is None


def test_subgroup_plan_document_is_outcome_free(tmp_path):
    from measure_it.byod.subgroup import plan_document
    out = _maestro(tmp_path)
    m = yaml.safe_load((out / "manifest.yaml").read_text())
    meta = {"dataset_id": "byod_pytest_maestro", "labels_json": json.dumps(
        {"long_covid": {"condition_id": "long_covid", "definition": "d", "label_basis": "clinical_diagnosis"}})}
    a = plan_document(m, meta)
    assert a == plan_document(m, meta) and a["command"] == "subgroup"
    assert a["device_positive"] == {"block": "device_capillaroscopy", "feature": "abnormality_score",
                                    "threshold": 2.0, "direction": ">="}
    assert "DEVICE and OMICS vocabularies never enter the model" in a["feature_source"]
    m["subgroup_analysis"]["max_features"] = 3
    assert K.sha256_json(plan_document(m, meta)) != K.sha256_json(a)


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_end_to_end_maestro_subgroup(tmp_path, monkeypatch):
    """ingest -> subgroup -> (changed plan refused) -> remove on the SYNTHETIC MAESTRO-like example."""
    from measure_it.byod import subgroup as SG
    from measure_it.byod.ingest import ingest
    from measure_it.byod.remove import remove
    from measure_it.harmonize import schema as HS
    from measure_it.harmonize.features import feature_matrix
    from measure_it.harmonize.phenotype import score_phenotype, validate_phenotype
    from measure_it.store import read_table, table_exists
    from measure_it.validate import check_person_geography
    monkeypatch.setattr(K, "RESULTS_ROOT", tmp_path / "results_byod")
    monkeypatch.setattr(K, "LEDGER", tmp_path / "results_byod" / "LEDGER.jsonl")
    monkeypatch.setenv(K.DEMO_ENV, "0")
    ds = "byod_pytest_maestro"
    out = _maestro(tmp_path)
    try:
        res = ingest(out, echo=lambda *_: None)
        assert f"person_concepts__{ds}" in res["partitions"]
        pc = read_table(f"person_concepts__{ds}")
        assert HS.check_concepts(pc) == []
        assert set(pc["vocabulary"]) == {"DEMOG", "ICD10CM", "LOINC", "SURVEY", "DEVICE", "OMICS"}
        ng = pc[pc["platform"].astype(str).str.startswith("nmr_nightingale;")]
        assert set(ng["concept_code"]) == {"2093-3", "2085-9", "2089-1", "2571-8", "2345-7", "2160-0", "1751-7",
                                           "1884-6", "1869-7", "82730-3"}
        assert ng.loc[ng["concept_code"] == "2093-3", "value_as_number"].median() == pytest.approx(5.0 * 38.67, rel=0.1)
        assert "nightingale:LDL_C" in set(pc["concept_code"]) and "nightingale:Total_C" not in set(pc["concept_code"])
        hx = pc[pc["mapping_method"] == "self_report_map"]
        assert len(hx) and set(hx["mapping_confidence"]) <= {"medium", "low"}
        assert check_person_geography()["ok"]
        r = SG.subgroup(out, echo=lambda *_: None)
        ph = r["phenotype"]
        assert validate_phenotype(ph) == []
        keys = [f["feature_key"] for f in ph["features"]]
        assert 1 <= len(keys) <= 8
        assert not [k for k in keys if k.split(":")[0] in ("DEVICE", "OMICS")]        # never the device itself
        assert "ICD10CM:U09" not in keys                                               # label-defining
        assert {"SURVEY:COMPASS31:vasomotor", "ICD10CM:I73", "LOINC:82730-3"} & set(keys)   # planted signal
        assert ph["performance"]["auroc"] > 0.6
        for f in ph["features"]:
            assert f["ehr_evaluable"] == (f["vocabulary"] in ("ICD10CM", "LOINC", "RXNORM", "DEMOG"))
        st = ph["strata_fractions"]
        assert st[0]["age_band"] == "all" and st[0]["n_cases"] == ph["base_population"][
            "n_cases_with_device_finding_defined"]
        assert all(s["fraction"] is None for s in st if s["suppressed"])
        # the phenotype scored on the dataset's own features reproduces the refit model (complete rows)
        fm = feature_matrix([ds], keys)
        X = r["feature_audit"]
        assert len(X) >= len(keys)
        p = score_phenotype(ph, fm)
        assert p.between(0, 1).all() and p.attrs["missing_features"] == []
        # aggregate tables: no participant ids, no geography
        for t in (K.PHENOTYPES_TABLE, K.STRATA_TABLE):
            tab = read_table(t)
            assert (tab["dataset_id"] == ds).any()
            assert not tab.astype(str).apply(lambda s: s.str.contains(f"{ds}:SYNTHETIC", regex=False)).any().any()
            assert not [c for c in tab.columns if C.column_problem(c) and c not in ("source_geographic_resolution",)]
        assert (K.results_path(ds) / "SUBGROUP.md").exists()
        assert json.loads((K.results_path(ds) / "computable_phenotype.json").read_text())["plan_sha256"] == \
            r["plan"]["sha256"]
        # the locked plan: unchanged -> verified, changed -> refused
        assert SG.subgroup(out, echo=lambda *_: None)["plan"]["status"] == "verified"
        _edit_manifest(out, lambda m: m["subgroup_analysis"].update(max_features=3))
        ingest(out, echo=lambda *_: None, rebuild=False)
        with pytest.raises(RuntimeError, match="locked"):
            SG.subgroup(out, echo=lambda *_: None)
    finally:
        remove(ds, rescore=False, echo=lambda *_: None)
    assert not K.dataset_partitions(ds)
    for t in (K.PHENOTYPES_TABLE, K.STRATA_TABLE):
        assert not table_exists(t) or ds not in set(read_table(t)["dataset_id"])
