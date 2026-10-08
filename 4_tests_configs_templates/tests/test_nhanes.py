"""Tests for measure_it.ingestion.nhanes and measure_it.wearables.nhanes_features.

Unit tests run without data; @pytest.mark.data tests check the processed outputs
(run `uv run python -m measure_it.ingestion.nhanes` first).
"""
import json

import numpy as np
import pandas as pd
import pytest
import yaml

from measure_it.config import PROCESSED, RAW
from measure_it.ingestion import nhanes as N
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.registry import REQUIRED_FIELDS
from measure_it.wearables import nhanes_features as NF


# ------------------------------------------------------------------ unit tests (no data needed)

def _ibm_bytes(hexstr: str) -> np.ndarray:
    return np.frombuffer(bytes.fromhex(hexstr), dtype=np.uint8).reshape(1, -1)


@pytest.mark.parametrize("hexstr,expected", [
    ("4110000000000000", 1.0), ("C110000000000000", -1.0), ("4080000000000000", 0.5),
    ("4264000000000000", 100.0), ("0000000000000000", 0.0), ("4210", 16.0),  # 2-byte truncated numeric
])
def test_ibm_to_float_known_values(hexstr, expected):
    assert N._ibm_to_float(_ibm_bytes(hexstr))[0] == pytest.approx(expected)


def test_ibm_to_float_sas_missing_is_nan():
    assert np.isnan(N._ibm_to_float(_ibm_bytes("2E00000000000000"))[0])   # '.'
    assert np.isnan(N._ibm_to_float(_ibm_bytes("4100000000000000"))[0])   # '.A' style special missing


def test_xport_reader_roundtrip(tmp_path):
    pyreadstat = pytest.importorskip("pyreadstat")
    df = pd.DataFrame({"SEQN": [62161.0, 62162.0, 62163.0], "VAL": [1.5, np.nan, -2.25],
                       "CODE": ["1", "", "9"], "NAME": ["ALPHA", "B", ""]})
    path = tmp_path / "TEST.xpt"
    pyreadstat.write_xport(df, str(path), file_format_version=5, table_name="TEST")
    lay = N.xport_layout(path)
    assert lay["n_rows"] == 3 and [v["name"] for v in lay["variables"]] == ["SEQN", "VAL", "CODE", "NAME"]
    out = N.read_xport(path)
    assert out["SEQN"].tolist() == [62161.0, 62162.0, 62163.0]
    assert out["VAL"].iloc[0] == 1.5 and np.isnan(out["VAL"].iloc[1]) and out["VAL"].iloc[2] == -2.25
    assert out["CODE"].iloc[0] == "1" and pd.isna(out["CODE"].iloc[1]) and pd.isna(out["NAME"].iloc[2])
    # chunked iteration and the uint8 fast path for 1-character columns
    chunks = list(N.iter_xport(path, ["SEQN", "CODE"], chunk_rows=2, char_as="uint8"))
    assert [len(c) for c in chunks] == [2, 1]
    assert chunks[0]["CODE"].dtype == np.uint8 and chunks[0]["CODE"].iloc[0] == ord("1")
    assert NF._digit(np.array([ord("1"), ord(" "), ord("9")], dtype=np.uint8)).tolist() == [1, -1, 9]


def test_xport_reader_strips_multiple_blank_padding_rows(tmp_path):
    # 2 numeric vars -> 16-byte rows; 1 row is padded with 64 blanks = 4 phantom rows (regression test)
    pyreadstat = pytest.importorskip("pyreadstat")
    path = tmp_path / "SHORT.xpt"
    pyreadstat.write_xport(pd.DataFrame({"SEQN": [62161.0], "LBXTC": [180.0]}), str(path),
                           file_format_version=5, table_name="SHORT")
    assert N.xport_layout(path)["n_rows"] == 1
    assert N.read_xport(path)["SEQN"].tolist() == [62161.0]


def test_recode_helpers():
    s = pd.Series([1.0, 2.0, 7.0, 9.0, np.nan])
    assert N.yes_no(s).tolist()[:2] == [1.0, 0.0] and N.yes_no(s).isna().sum() == 3
    r = N.in_range(pd.Series([0, 30, 77, 99, 5]), 0, 30)
    assert r.tolist()[:2] == [0, 30] and r.isna().sum() == 2 and r.iloc[4] == 5


def test_participant_id_namespacing():
    assert N.participant_id(pd.Series([62161.0, 73557])).tolist() == ["nhanes:62161", "nhanes:73557"]


def test_phq9_total_requires_all_items():
    rows = {c: [1, 1, 3] for c in N.PHQ9_ITEMS}
    rows["DPQ040"] = [2, 7, 3]           # 7 = refused -> missing item
    rows["SEQN"] = [1, 2, 3]
    out = N._phq9(pd.DataFrame(rows)).set_index("SEQN")
    assert out.loc[1, "phq9_total"] == 10 and out.loc[1, "phq9_ge10"] == 1.0
    assert np.isnan(out.loc[2, "phq9_total"]) and out.loc[2, "phq9_n_items"] == 8
    assert out.loc[3, "phq9_total"] == 27 and out.loc[1, "phq9_item4_tired_little_energy_0to3"] == 2


def test_rx_target_group_prefix_matching():
    assert N.rx_target_group("M79.7") == "fibromyalgia"
    assert N.rx_target_group("G43.P") == "migraine"          # NHANES 'P' = prevention use
    assert N.rx_target_group("K58.9") == "irritable bowel syndrome"
    assert N.rx_target_group("R53.83") == "fatigue / malaise"
    assert N.rx_target_group("M79.1") == "myalgia"           # not fibromyalgia
    assert N.rx_target_group("E11.2") == "" and N.rx_target_group(np.nan) == ""
    # G89 is 'Pain, not elsewhere classified': only G89.2x / G89.4 are chronic (regression: G89.18 acute
    # postprocedural pain and G89.3 neoplasm-related pain were labelled 'chronic pain')
    assert N.rx_target_group("G89.4") == "chronic pain" and N.rx_target_group("G89.29") == "chronic pain"
    assert N.rx_target_group("G89.18") == "" and N.rx_target_group("G89.3") == ""


def _synthetic_days(seqn=1, n_valid=5, weekend_valid=1):
    rows = []
    for d in range(1, 10):
        full = 2 <= d <= 8
        valid = full and (d - 2) < n_valid
        rows.append({"SEQN": seqn, "cycle_code": "G", "day_index": d, "day_of_week": 1 if d == 2 else 3,
                     "minutes_with_data": 1440 if full else 600, "wake_wear_min": 900 if valid else 100,
                     "sleep_wear_min": 450, "nonwear_min": 0, "unknown_min": 0,
                     "total_mims": 10000.0 + 1000 * d})
    df = pd.DataFrame(rows)
    df["weekend"] = df["day_of_week"].isin([1, 7])
    df["full_day"] = df["minutes_with_data"].eq(1440)
    df["valid_day"] = df["full_day"] & (df["wake_wear_min"] >= NF.VALID_DAY_MIN_WAKE_WEAR)
    return df


def test_day_level_features_rule_and_weekend_difference():
    f = NF.day_level_features(_synthetic_days(n_valid=5)).loc[1]
    assert f["n_valid_days"] == 5 and bool(f["passes_valid_wear_rule"])
    vals = np.array([12000, 13000, 14000, 15000, 16000.0])  # days 2..6
    assert f["mean_daily_mims"] == pytest.approx(vals.mean())
    assert f["cv_daily_mims"] == pytest.approx(vals.std(ddof=1) / vals.mean())
    assert f["weekend_minus_weekday_mims"] == pytest.approx(12000 - vals[1:].mean())  # day 2 is Sunday
    g = NF.day_level_features(_synthetic_days(seqn=2, n_valid=3)).loc[2]
    assert g["n_valid_days"] == 3 and not bool(g["passes_valid_wear_rule"])


def test_minute_features_on_synthetic_participant():
    n_days = 9
    minute = np.arange(n_days * 1440)
    tod = minute % 1440
    awake = (tod >= 7 * 60) & (tod < 23 * 60)
    mims = np.where(awake, np.where((tod // 60) % 2 == 0, 20.0, 5.0), 0.5)  # alternating active/inactive hours
    pred = np.where(awake, 1, 2).astype(np.int8)
    qf = np.zeros_like(minute, dtype=np.int16)
    f, daily = NF.participant_minute_features(minute, mims, pred, qf, valid_days=np.arange(2, 9), n_days=n_days)
    assert f["interdaily_stability"] == pytest.approx(1.0)          # identical days
    assert f["sedentary_fraction"] == pytest.approx(0.5, abs=0.02)  # half of waking hours below 10.558
    assert f["l5_mims"] == pytest.approx(0.5) and 0 < f["relative_amplitude"] < 1
    assert f["mean_active_bout_min"] == pytest.approx(60, abs=1) and f["astp"] == pytest.approx(1 / 60, rel=0.05)
    assert f["sleep_proxy_duration_h"] == pytest.approx(8.0)
    assert f["sleep_proxy_midpoint_h"] == pytest.approx(3.0, abs=0.05)
    assert f["sleep_proxy_onset_h"] == pytest.approx(23.0, abs=0.05)
    assert len(daily) == 7 and daily[0]["minute_sum_mims_qc_valid"] == pytest.approx(mims[1440:2880].sum())


def test_restore_mims_precision_keeps_cut_point_minutes_active():
    # regression: the interim parquet stores MIMS as float32 and float32(10.558) = 10.5579996 < 10.558,
    # which turned every minute recorded exactly at the cut-point into an 'inactive' minute
    f32 = np.array([10.558, 10.557, 0.0, -0.01, 608.252], dtype=np.float32)
    assert not float(f32[0]) >= NF.MIMS_ACTIVE_THRESHOLD
    r = NF.restore_mims_precision(f32)
    assert r.dtype == np.float64 and r.tolist() == [10.558, 10.557, 0.0, -0.01, 608.252]
    assert r[0] >= NF.MIMS_ACTIVE_THRESHOLD and r[1] < NF.MIMS_ACTIVE_THRESHOLD
    # end to end through the feature function: waking minutes exactly at the cut-point are active
    n_days = 9
    minute = np.arange(n_days * 1440)
    awake = (minute % 1440 >= 7 * 60) & (minute % 1440 < 23 * 60)
    mims32 = np.where(awake, 10.558, 0.5).astype(np.float32)
    pred = np.where(awake, 1, 2).astype(np.int8)
    f, _ = NF.participant_minute_features(minute, NF.restore_mims_precision(mims32), pred,
                                          np.zeros_like(minute, dtype=np.int16), valid_days=np.arange(2, 9),
                                          n_days=n_days)
    assert f["sedentary_fraction"] == 0.0 and f["active_min_per_day"] == 16 * 60


# ------------------------------------------------------------------ data tests (processed outputs)

def _table(name):
    path = PROCESSED / f"{name}.parquet"
    if not path.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(path)


@pytest.mark.data
@pytest.mark.parametrize("component,cycle", [("DEMO", "G"), ("MCQ", "H"), ("DPQ", "G"), ("SLQ", "H"),
                                             ("PFQ", "G"), ("RXQ_RX", "H"), ("PAXHD", "G"), ("DLQ", "H")])
def test_reader_reproduces_codebook_frequencies(component, cycle):
    """Every single-code value count printed in the NCHS codebook equals the count in the file we read."""
    path = N.xpt_path(component, cycle)
    if not path.exists():
        pytest.skip("raw file not downloaded")
    df = N.read_xport(path)
    cb = N.parse_codebook(component, cycle)
    checked = 0
    for var, meta in cb.items():
        if var not in df.columns or var == "SEQN":
            continue
        for code, n in meta.get("counts", {}).items():
            if not code.lstrip("-").isdigit():
                continue
            col = df[var]
            got = int((pd.to_numeric(col, errors="coerce") == int(code)).sum())
            assert got == n, (component, cycle, var, code, got, n)
            checked += 1
    assert checked > 5


@pytest.mark.data
def test_participants_table():
    p = _table("participants__nhanes")
    assert p["participant_id"].is_unique and p["participant_id"].str.startswith("nhanes:").all()
    assert p.groupby("cycle").size().to_dict() == {"2011-2012": 9756, "2013-2014": 10175}
    assert int(p["has_pam_data"].sum()) == 6917 + 7776  # PAXSTS = 1 counts stated in the PAX documentation
    for c in PROVENANCE_COLUMNS + ["sdmvpsu", "sdmvstra", "wtmec2yr", "wtmec4yr_pooled"]:
        assert c in p.columns
    assert (p["data_layer"] == "person").all() and (p["source_geographic_resolution"] == "none").all()
    assert np.allclose(p["wtmec4yr_pooled"], p["wtmec2yr"] / 2)


@pytest.mark.data
def test_clinical_features_one_row_per_participant_and_no_hr_stream():
    p = _table("participants__nhanes")
    c = _table("participant_clinical_features__nhanes")
    assert c["participant_id"].is_unique and set(c["participant_id"]) == set(p["participant_id"])
    hr_like = [x for x in c.columns if "hrv" in x.lower() or "heart_rate" in x.lower()]
    assert hr_like == []
    assert "exam_pulse_60s_bpm" in c.columns and c["exam_pulse_60s_bpm"].dropna().gt(0).all()
    assert c["phq9_total"].dropna().between(0, 27).all()
    assert c["healthcare_visits_past_year_cat_0none_1one_2two3_3four9_4ten12_5thirteenplus"].dropna().between(0, 5).all()
    # harmonised visits exist in both cycles
    assert c.groupby("cycle")["healthcare_visits_past_year_cat_0none_1one_2two3_3four9_4ten12_5thirteenplus"].count().min() > 9000


@pytest.mark.data
def test_conditions_labels_and_rx_reason_scope():
    c = _table("participant_conditions__nhanes")
    rx = c[c["evidence_basis"] == N.RX_BASIS]
    assert N.RX_BASIS == "condition inferred from prescription reason-for-use code"
    assert len(rx) > 0 and set(rx["cycle"]) == {"2013-2014"}  # RXDRSC exists only in RXQ_RX_H
    assert rx["icd10cm_code"].str.match(r"^[A-Z][0-9]").all()
    assert set(c["evidence_type"]) == {"person_self_report"}
    fib = rx[rx["rx_target_group"] == "fibromyalgia"]
    assert fib["icd10cm_code"].str.startswith("M79.7").all() and fib["participant_id"].nunique() > 0
    pain = rx[rx["rx_target_group"] == "chronic pain"]
    assert pain["icd10cm_code"].str.startswith(("G89.2", "G89.4")).all()
    assert (rx.loc[rx["icd10cm_code"] == "G89.18", "rx_target_group"] == "").all()


@pytest.mark.data
def test_labs_and_medications():
    labs = _table("participant_labs__nhanes")
    assert not labs.duplicated(["participant_id", "lab_variable"]).any()
    assert (labs["unit"] != "").mean() > 0.95
    assert set(labs.loc[labs["lab_variable"] == "LBXTSH1", "cycle"]) == {"2011-2012"}
    meds = _table("participant_medications__nhanes")
    assert meds.loc[meds["cycle"] == "2011-2012", "reason_icd10cm_1"].isna().all()
    assert meds.loc[meds["cycle"] == "2013-2014", "reason_icd10cm_1"].notna().mean() > 0.9
    # RXQSEEN code 3 ('Only pharmacy print out seen', 2013-2014) is kept, not dropped to missing
    printout = meds["container_seen_detail"] == "only pharmacy print out seen"
    assert printout.sum() > 0 and set(meds.loc[printout, "cycle"]) == {"2013-2014"}
    assert (meds.loc[printout, "container_seen"] == False).all()  # noqa: E712
    assert meds["container_seen"].isna().sum() == meds["container_seen_detail"].isna().sum()


@pytest.mark.data
def test_mortality_under18_not_released():
    m = _table("participant_mortality__nhanes")
    c = _table("participant_clinical_features__nhanes").set_index("participant_id")
    kids = m[m["participant_id"].map(c["age_years"]) < 18]
    assert (kids["mort_eligstat"] == "under age 18, not released").all()
    assert (m["mort_followup_end"] == "2019-12-31").all()


@pytest.mark.data
def test_wearable_features_rule_ranges_and_no_hr():
    f = _table("participant_wearable_features__nhanes")
    p = _table("participants__nhanes")
    assert f["participant_id"].is_unique and set(f["participant_id"]) <= set(p["participant_id"])
    assert (f["n_valid_days"] >= NF.MIN_VALID_DAYS).all()
    assert not [x for x in f.columns if "hr" in x.lower().split("_") or "hrv" in x.lower() or "heart" in x.lower()]
    for col, lo, hi in [("sedentary_fraction", 0, 1), ("interdaily_stability", 0, 1), ("relative_amplitude", 0, 1),
                        ("m10_onset_h", 0, 24), ("l5_onset_h", 0, 24), ("sleep_proxy_midpoint_h", 0, 24)]:
        if col in f.columns and f[col].notna().any():
            assert f[col].dropna().between(lo, hi).all(), col
    if "interdaily_stability" in f.columns:
        assert f["interdaily_stability"].notna().mean() > 0.95
        assert f["intradaily_variability"].dropna().between(0, 4).all()


@pytest.mark.data
def test_wearable_features_match_independent_raw_recomputation():
    """Values recomputed during the 2026-09-23 audit straight from PAXMIN_H.xpt (read with pyreadstat, own
    alignment from PAXFTIME + PAXSSNMP/80, own masks). SEQN 79583 has one waking minute at exactly
    MIMS = 10.558, so it also guards the float32 cut-point regression (the old output was 319.2857 / 0.34635)."""
    f = _table("participant_wearable_features__nhanes").set_index("participant_id")
    if "nhanes:79583" not in f.index:
        pytest.skip("participant not in table")
    r = f.loc["nhanes:79583"]
    assert r["n_valid_days"] == 7
    assert r["sedentary_min_per_day"] == pytest.approx(2234 / 7, abs=1e-6)
    assert r["sedentary_fraction"] == pytest.approx(0.34620, abs=5e-6)
    assert r["interdaily_stability"] == pytest.approx(0.72538, abs=5e-6)
    assert r["relative_amplitude"] == pytest.approx(0.92004, abs=5e-6)
    assert r["mean_daily_mims"] == pytest.approx(16340.81214, abs=1e-3)
    assert 0 <= f["sleep_proxy_onset_h"].dropna().max() < 24  # clock hours in [0, 24)


@pytest.mark.data
def test_wearable_daily_valid_rule_and_minute_day_sum_agreement():
    d = _table("participant_wearable_daily__nhanes")
    v = d[d["valid_day"]]
    assert v["full_day"].all() and (v["wake_wear_min"] >= NF.VALID_DAY_MIN_WAKE_WEAR).all()
    if "minute_sum_mims_qc_valid" in d.columns:
        chk = v.dropna(subset=["minute_sum_mims_qc_valid"])
        assert len(chk) > 0
        assert (chk["minute_sum_mims_qc_valid"] - chk["total_mims"]).abs().max() < 0.05


@pytest.mark.data
def test_registry_entry_and_audit():
    entry_path = RAW / N.SOURCE_ID / "registry_entry.yaml"
    audit = RAW / N.SOURCE_ID / "DATA_AUDIT.md"
    if not entry_path.exists() or not audit.exists():
        pytest.skip("registry/audit not written")
    entry = yaml.safe_load(entry_path.read_text())
    assert all(k in entry for k in REQUIRED_FIELDS)
    f = _table("participant_wearable_features__nhanes")
    assert entry["sample_size"]["passing_valid_wear_rule"] == len(f)
    text = audit.read_text()
    assert f"| pass valid-wear rule (>= 4 valid days) |" in text and str(len(f)) in text
    manifest = json.loads((RAW / N.SOURCE_ID / "MANIFEST.json").read_text())
    assert "DEMO_G.xpt" in manifest["files"] and len(manifest["files"]["DEMO_G.xpt"]["sha256"]) == 64
