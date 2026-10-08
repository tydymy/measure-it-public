"""Tests for measure_it.ingestion.openfda and measure_it.measurements.regulatory.

Unit tests need no network or data. Tests marked `data` read data/processed outputs
(run `uv run python -m measure_it.ingestion.openfda` first).
"""
from __future__ import annotations

import json
import re

import pandas as pd
import pytest

from measure_it.config import RAW, UNKNOWN, load_config
from measure_it.ingestion import openfda as of
from measure_it.measurements import regulatory as reg
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table, table_exists

MEASUREMENT_IDS = [m["id"] for m in load_config("measurements")["measurement_classes"]]
TABLES = ["fda_device_classification", "fda_510k", "fda_pma", "fda_registration_listing_counts",
          "measurement_regulatory_status", "fda_named_device_findings", "fda_device_query_log"]
GENERIC_CODES = {"DQK", "GWQ", "INQ", "KTB", "MXK"}  # FDA names that do not describe any measurement


# ------------------------------------------------------------------------------------------ unit
def test_build_term_forms():
    assert of.build_term("device_name", "oximeter") == "device_name:oximeter"
    assert of.build_term("device_name", "heart rate") == 'device_name:"heart rate"'
    assert of.build_term("definition", "c-reactive protein") == 'definition:"c-reactive protein"'
    assert of.build_term("device_name", "actigraph*") == "device_name:actigraph*"
    # a bare space is OR in openFDA, so multi-word prefix searches must AND their words
    assert of.build_term("device_name", "flow cytometr*") == "(device_name:flow AND device_name:cytometr*)"


def test_classification_search_ors_fields():
    assert of.classification_search("sleep assessment") == (
        'device_name:"sleep assessment" definition:"sleep assessment"')


@pytest.mark.parametrize("num,kind", [("K201525", "510k"), ("DEN180044", "de_novo"), ("den1", "de_novo"),
                                      ("P980016", "other"), ("", "other")])
def test_submission_kind(num, kind):
    assert of.submission_kind(num) == kind


def test_example_devices_keeps_de_novo_beyond_top_n():
    recs = [{"k_number": f"K2{i:05d}", "device_name": f"d{i}", "applicant": "a",
             "decision_date": f"2024-01-{i + 1:02d}"} for i in range(8)]
    recs.append({"k_number": "DEN180044", "device_name": "first", "applicant": "a", "decision_date": "2018-09-11"})
    ex = of.example_devices(recs, n_recent=3)
    ids = [e["clearance_id"] for e in ex]
    assert ids[:3] == ["K200007", "K200006", "K200005"]
    assert "DEN180044" in ids and ex[-1]["submission_kind"] == "de_novo"


def test_pma_examples_only_originals():
    recs = [{"pma_number": "P1", "supplement_number": "", "trade_name": "t", "applicant": "a", "decision_date": "2001-01-01"},
            {"pma_number": "P1", "supplement_number": "S001", "trade_name": "t", "applicant": "a", "decision_date": "2005-01-01"}]
    ex = of.pma_example_devices(recs)
    assert [e["clearance_id"] for e in ex] == ["P1"]


def test_regulatory_visibility():
    assert of.regulatory_visibility(None, None, None, False) == "no_product_code_found"
    assert of.regulatory_visibility(0, 0, 0, True) == "category_exists_no_decisions_in_openfda"
    assert of.regulatory_visibility(0, 1, 0, True) == "decisions_on_record"


def test_search_builders():
    assert of.product_code_search("DQK") == "product_code:DQK"
    assert of.product_code_search("DQK", "device_name:vendys") == "product_code:DQK AND device_name:vendys"
    assert of.registration_search("QDA") == "products.product_code:QDA"
    assert of.registration_search("DQK", ["K2", "K1"]) == "products.product_code:DQK AND k_number:(K1 K2)"
    # registrationlisting keeps PMA numbers in pma_number, not k_number
    assert of.registration_search("DQK", ["P880040"]) == "products.product_code:DQK AND pma_number:(P880040)"
    assert of.registration_search("DQK", ["K1", "DEN2", "P3"]) == (
        "products.product_code:DQK AND (k_number:(DEN2 K1) OR pma_number:(P3))")


def test_query_log_latest_fetch_uses_fetch_time():
    log = of.QueryLog()
    log.add({"query_id": "a", "endpoint": "510k", "params": {"search": "x"}, "fetched_at": "2026-09-23T20:35:00+00:00"})
    log.add({"query_id": "b", "endpoint": "510k", "params": {"search": "y"}, "fetched_at": "2026-09-23T20:45:00+00:00"})
    assert log.latest_fetch() == "2026-09-23T20:45:00+00:00"
    assert log.latest_fetch(["510k: x"]) == "2026-09-23T20:35:00+00:00"
    assert log.latest_fetch(["510k: never-run"]) is None


def _fake_log(search, codes, endpoint="classification"):
    log = of.QueryLog()
    log.add({"query_id": "q1", "endpoint": endpoint, "params": {"search": search},
             "result_product_codes": codes, "purpose": "t"})
    return log


def test_validate_code_map_catches_problems():
    cls = {"QDA": {"device_name": "Electrocardiograph Software For Over-The-Counter Use"}}
    search = of.classification_search("electrocardiograph*")
    good = {"product_code": "QDA", "device_name": "Electrocardiograph Software For Over-The-Counter Use",
            "device_class": "2", "regulation_number": "870.2345", "medical_specialty_description": "Cardiovascular",
            "mapping_rationale": "x", "mapping_confidence": "high",
            "query_that_found_it": f"classification: {search}"}
    ok_map = {"measurements": {"ecg_ambulatory": [good],
                               "capillaroscopy": [{"product_code": None, "status": of.NO_CODE,
                                                   "queries_tried": [f"classification: {search}"]}]}}
    mids = {"ecg_ambulatory", "capillaroscopy"}
    assert of.validate_code_map(ok_map, cls, _fake_log(search, ["QDA"]), mids) == []
    # a 'no code' claim must cite queries that were actually made
    unlogged = {"measurements": {**ok_map["measurements"], "capillaroscopy": [
        {"product_code": None, "status": of.NO_CODE, "queries_tried": ["classification: device_name:nailfold"]}]}}
    assert any("queries_tried" in p for p in of.validate_code_map(unlogged, cls, _fake_log(search, ["QDA"]), mids))
    # query did not return the code
    assert any("did not return" in p for p in of.validate_code_map(ok_map, cls, _fake_log(search, ["DPS"]), mids))
    # query never made
    assert any("not in the query log" in p for p in of.validate_code_map(ok_map, cls, _fake_log("x", ["QDA"]), mids))
    # name must equal FDA's
    bad = {"measurements": {**ok_map["measurements"], "ecg_ambulatory": [{**good, "device_name": "ECG"}]}}
    assert any("device_name" in p for p in of.validate_code_map(bad, cls, _fake_log(search, ["QDA"]), mids))
    # a class may not silently disappear, and an empty entry must say NO_CODE with queries
    missing = {"measurements": {"ecg_ambulatory": [good]}}
    assert any("capillaroscopy" in p for p in of.validate_code_map(missing, cls, _fake_log(search, ["QDA"]), mids))
    empty = {"measurements": {**ok_map["measurements"], "capillaroscopy": [{"product_code": None}]}}
    assert len(of.validate_code_map(empty, cls, _fake_log(search, ["QDA"]), mids)) == 2


def test_query_log_merges_measurement_ids():
    log = of.QueryLog()
    log.add({"query_id": "q", "endpoint": "classification", "params": {"search": "s"}, "measurement_id": "ppg"})
    log.add({"query_id": "q", "endpoint": "classification", "params": {"search": "s"},
             "measurement_id": "endothelial_function"})
    assert log.entries["q"]["measurement_ids"] == ["ppg", "endothelial_function"]


def test_code_map_config_is_complete():
    cmap = of.load_code_map()
    assert set(cmap["measurements"]) == set(MEASUREMENT_IDS)
    for mid, e in of.iter_mappings(cmap):
        if e.get("product_code") is None:
            assert e["status"] == of.NO_CODE and e["queries_tried"] and e["note"]
            continue
        assert e["mapping_confidence"] in {"high", "medium", "low"}
        assert re.match(r"^(classification|510k): ", e["query_that_found_it"])
        if e["product_code"] in GENERIC_CODES:
            assert e.get("device_filter"), f"{mid}/{e['product_code']} generic code needs a device_filter"
            assert e["mapping_confidence"] == "low"


def test_class_matching_patterns_and_fuzzy():
    ms = load_config("measurements")["measurement_classes"]
    ids = lambda q: [m["measurement_id"] for m in reg.match_measurement_classes(q, ms)]  # noqa: E731
    assert "ecg_ambulatory" in ids("Holter monitor")
    assert ids("nailfold capillaroscopy")[0] == "capillaroscopy"
    assert "autonomic_testing" in ids("wearable autonomic monitoring")
    assert "wearable_heart_rate" in ids("wearable autonomic monitoring")
    assert ids("hrv") == ["hrv"]
    # questionnaires are excluded from objective sleep measurement, in every matching layer
    assert "sleep_objective" not in ids("Pittsburgh Sleep Quality Index")
    assert ids("sleep questionnaire") == []
    assert ids("PSQI sleep score") == []
    assert ids("banana smoothie") == []


def test_class_matching_is_ordered_not_unioned():
    """A shared generic word must not attach other classes to an exact/pattern match."""
    ms = load_config("measurements")["measurement_classes"]
    ids = lambda q: [m["measurement_id"] for m in reg.match_measurement_classes(q, ms)]  # noqa: E731
    assert ids("respiratory rate") == ["respiratory_rate"]          # was + hrv/ecg/ppg/autonomic via 'rate'
    assert ids("heart rate variability") == ["hrv"]
    assert ids("pulse oximeter") == ["continuous_spo2"]             # was + ppg via 'pulse'
    assert ids("flow cytometry") == ["immune_assays"]               # was + capillaroscopy via 'flow'
    assert ids("gait speed") == ["digital_gait"]
    # a fuzzy match that explains an otherwise unexplained word is still kept
    assert set(ids("HRV and sleep")) == {"hrv", "sleep_objective"}
    # 4-letter tokens are not fuzzy-matched (heat~heart), but plurals of a stem are
    assert ids("heat therapy") == []
    assert "wearable_heart_rate" in ids("heart rates")
    assert reg.token_match("rate", "rates") and not reg.token_match("heat", "heart")
    assert not reg.token_match("sweat", "seat")


def test_device_fallback_matching():
    devices = pd.DataFrame([
        {"clearance_id": "K032519", "submission_kind": "510k", "device_name": "ENDO PAT 2000",
         "applicant": "Itamar Medical", "decision_date": "2003-11-12", "product_code": "DQK",
         "mapped_measurement_ids": '["endothelial_function"]', "source_name": "s", "retrieved_at": "t"},
        {"clearance_id": "K113862", "submission_kind": "510k", "device_name": "ZIO PATCH",
         "applicant": "iRhythm Technologies, Inc.", "decision_date": "2012-02-06", "product_code": "DSH",
         "mapped_measurement_ids": '["ecg_ambulatory"]', "source_name": "s", "retrieved_at": "t"},
    ])
    cls, rows = reg.match_devices("EndoPAT", devices)
    assert [c["measurement_id"] for c in cls] == ["endothelial_function"] and rows[0]["clearance_id"] == "K032519"
    cls, _ = reg.match_devices("iRhythm Zio", devices)
    assert [c["measurement_id"] for c in cls] == ["ecg_ambulatory"]
    assert reg.match_devices("Sudoscan", devices) == ([], [])
    # 'sweat' must not fuzzy-match 'seat' (was: Q-Sweat -> "The Heart Seat" -> remote monitoring)
    seat = pd.DataFrame([{"clearance_id": "K222330", "submission_kind": "510k", "device_name": "The Heart Seat",
                          "applicant": "Casana Care, Inc.", "decision_date": "2023-01-01", "product_code": "MWI",
                          "mapped_measurement_ids": '["remote_patient_monitoring"]', "source_name": "s",
                          "retrieved_at": "t"}])
    assert reg.match_devices("Q-Sweat", seat) == ([], [])


# ------------------------------------------------------------------------------------------ data
needs_data = pytest.mark.skipif(not table_exists("measurement_regulatory_status"),
                                reason="run uv run python -m measure_it.ingestion.openfda first")


@pytest.fixture(scope="module")
def t():
    return {name: read_table(name) for name in TABLES}


@pytest.mark.data
@needs_data
def test_tables_have_provenance(t):
    for name, df in t.items():
        assert len(df) > 0, name
        assert set(PROVENANCE_COLUMNS) <= set(df.columns), name
        assert (df["evidence_type"].isin({"regulatory_record", "metadata_catalog"})).all(), name
        assert df["retrieved_at"].astype(str).str.match(r"^\d{4}-\d\d-\d\dT").all(), name


@pytest.mark.data
@needs_data
def test_status_covers_every_measurement_class(t):
    st = t["measurement_regulatory_status"]
    assert set(st["measurement_id"]) == set(MEASUREMENT_IDS)
    assert (st["regulatory_note"] == of.REGULATORY_NOTE).all()
    nocode = st[st["product_code"].isna()]
    assert (nocode["device_name"] == of.NO_CODE).all()
    assert nocode["n_510k"].isna().all() and nocode["n_pma"].isna().all()  # unknown, not zero
    assert (nocode["regulatory_visibility"] == "no_product_code_found").all()
    coded = st[st["product_code"].notna()]
    for c in ("n_510k", "n_denovo", "n_pma", "n_registered_establishments"):
        assert coded[c].notna().all() and (coded[c] >= 0).all(), c
    assert (coded["n_registered_establishments_us"] <= coded["n_registered_establishments"]).all()


@pytest.mark.data
@needs_data
def test_example_devices_are_real_api_records(t):
    """Every example clearance id must be a record retrieved from openFDA (no invented clearances)."""
    known = set(t["fda_510k"]["clearance_id"]) | set(t["fda_pma"]["pma_number"])
    st = t["measurement_regulatory_status"]
    for _, r in st.iterrows():
        for ex in json.loads(r["example_devices"]):
            assert ex["clearance_id"] in known, (r["measurement_id"], ex)
            assert re.match(r"^(K\d{6}|DEN\d{6}|P\d{6})$", ex["clearance_id"]), ex


@pytest.mark.data
@needs_data
def test_510k_records_consistent(t):
    k = t["fda_510k"]
    assert k["clearance_id"].is_unique
    assert set(k["submission_kind"]) <= {"510k", "de_novo"}
    assert (k.loc[k["submission_kind"] == "de_novo", "clearance_id"].str.startswith("DEN")).all()
    assert k["decision_date"].str.match(r"^\d{4}-\d\d-\d\d$").all()
    # when every record was retrieved, the status counts equal the stored rows
    st = t["measurement_regulatory_status"]
    full = st[st["record_retrieval_mode"] == "all_records"]
    for _, r in full.iterrows():
        in_unit = k["device_filters"].map(lambda s, f=r["device_filter"]: f in json.loads(s))
        rows = k[(k["product_code"] == r["product_code"]) & in_unit]
        assert len(rows) == r["n_510k"] + r["n_denovo"], r["status_record_id"]
        assert (rows["mapped_measurement_ids"].map(lambda s, m=r["measurement_id"]: m in json.loads(s))).all()
    partial = st[st["record_retrieval_mode"].str.startswith("most_recent_")]
    for _, r in partial.iterrows():
        assert r["n_510k_records_retrieved"] == of.RECENT_N


@pytest.mark.data
@needs_data
def test_codes_match_fda_classification_and_queries(t):
    cls = t["fda_device_classification"].set_index("product_code")
    ql = t["fda_device_query_log"]
    for mid, e in of.iter_mappings(of.load_code_map()):
        code = e.get("product_code")
        if not code:
            continue
        assert cls.loc[code, "device_name"] == e["device_name"]
        ep, search = e["query_that_found_it"].split(": ", 1)
        q = ql[(ql["endpoint"] == ep) & (ql["search"] == search)]
        assert len(q) == 1, e["query_that_found_it"]
        assert code in json.loads(q.iloc[0]["result_product_codes"])


@pytest.mark.data
@needs_data
def test_filtered_generic_units_only_count_the_measurement(t):
    """Filtered generic codes must not sweep in unrelated devices of the same code."""
    k = t["fda_510k"]
    auto = k[k["mapped_measurement_ids"].map(lambda s: "autonomic_testing" in json.loads(s))]
    # cystic-fibrosis sweat-chloride systems (KTB) are not sudomotor/autonomic tests
    cf = auto["device_name"].str.contains(r"chloride|macroduct|nanoduct|cf indicator|sweat inducer|sweat collection",
                                         case=False, regex=True)
    assert not cf.any(), auto.loc[cf, ["clearance_id", "device_name"]].to_dict("records")
    ktb = auto[auto["product_code"] == "KTB"]
    assert set(ktb["clearance_id"]) == {"K992874"}  # Q-SWEAT
    inq = auto[auto["product_code"] == "INQ"]
    assert len(inq) > 0 and inq["device_name"].str.contains("tilt table", case=False).all()
    st = t["measurement_regulatory_status"]
    row = st[(st["measurement_id"] == "autonomic_testing") & (st["product_code"] == "KTB")].iloc[0]
    assert row["n_510k"] == 1


@pytest.mark.data
@needs_data
def test_registry_retrieved_at_is_fetch_time(t):
    import yaml
    entry = yaml.safe_load((RAW / of.SOURCE_ID / "registry_entry.yaml").read_text())
    assert str(entry["retrieved_at"]) == t["fda_device_query_log"]["retrieved_at"].max()
    st = t["measurement_regulatory_status"]
    nocode = st[st["product_code"].isna()]
    assert (nocode["retrieved_at"] <= t["fda_device_query_log"]["retrieved_at"].max()).all()


@pytest.mark.data
@needs_data
def test_registration_counts(t):
    r = t["fda_registration_listing_counts"]
    assert (r["n_establishments"] <= r["n_device_listings"]).all()
    assert (r["n_establishments_us"] <= r["n_establishments"]).all()
    assert (r["n_establishments_non_us"] >= 0).all()


@pytest.mark.data
@needs_data
def test_raw_audit_files_exist():
    d = RAW / of.SOURCE_ID
    for f in ("MANIFEST.json", "QUERY_LOG.json", "DATA_AUDIT.md", "registry_entry.yaml"):
        assert (d / f).exists(), f
    ql = json.loads((d / "QUERY_LOG.json").read_text())
    assert ql["n_queries"] == len(ql["queries"]) > 100


@pytest.mark.data
@needs_data
def test_regulatory_context_tool():
    reg._load.cache_clear()
    out = reg.get_regulatory_context("wearable ECG patch")
    assert out["status"] == "matched"
    assert "ecg_ambulatory" in [m["measurement_id"] for m in out["matched_measurements"]]
    assert out["regulatory_note"] == of.REGULATORY_NOTE and out["provenance"]["table"]
    codes = {r["product_code"] for r in out["regulatory_status"] if r["measurement_id"] == "ecg_ambulatory"}
    assert {"QDA", "MLO"} <= codes

    cap = reg.get_regulatory_context("nailfold capillaroscopy")
    rows = [r for r in cap["regulatory_status"] if r["measurement_id"] == "capillaroscopy"]
    assert rows and rows[0]["regulatory_visibility"] == "no_product_code_found"
    assert rows[0]["n_510k"] == UNKNOWN  # missing is UNKNOWN, never 0

    zio = reg.get_regulatory_context("Zio")  # brand name -> device-name fallback
    assert zio["status"] == "matched" and zio["matching_fda_records"]
    assert all(r["clearance_id"].startswith(("K", "DEN")) for r in zio["matching_fda_records"])

    none = reg.get_regulatory_context("quantum banana")
    assert none["status"] == UNKNOWN and none["regulatory_status"] == UNKNOWN
