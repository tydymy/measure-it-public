"""Tests for NHANES 2003-2006 ingestion (measure_it.ingestion.nhanes_2003_2006), hip accelerometry features
(measure_it.wearables.nhanes0306_features) and the replication models (measure_it.wearables.nhanes0306_models)."""
import numpy as np
import pandas as pd
import pytest

from measure_it.config import PROCESSED, TABLES
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.validate import GEO_COLUMN_RE
from measure_it.wearables import nhanes0306_features as H

# ------------------------------------------------------------------ unit tests


def test_nonwear_troiano_allows_two_minute_spikes():
    c = np.r_[np.full(30, 500), np.zeros(40), [50, 60], np.zeros(30), np.full(10, 300)]
    m = H.nonwear_mask(c)
    assert m[30:102].all() and not m[:30].any() and not m[102:].any()


def test_nonwear_three_spikes_break_the_interval():
    c = np.r_[np.zeros(40), [50, 60, 70], np.zeros(40)]
    assert not H.nonwear_mask(c).any()           # neither side reaches 60 zero minutes


def test_nonwear_big_count_breaks_and_short_runs_are_wear():
    c = np.r_[np.zeros(59), [5000], np.zeros(59)]
    assert not H.nonwear_mask(c).any()
    assert H.nonwear_mask(np.zeros(60)).all()


def test_bout_minutes():
    a = np.r_[np.ones(5), [0, 0], np.ones(4), [0, 0, 0], np.ones(12)].astype(bool)
    assert H.bout_minutes(a, np.ones(len(a), bool)) == 11 + 12
    v = np.ones(len(a), bool)
    v[20] = False                                # an invalid minute splits the 12-minute run
    assert H.bout_minutes(a, v) == 11


def test_participant_features_synthetic_week():
    rng = np.random.default_rng(1)
    counts = np.zeros(7 * 1440)
    for d in range(7):
        day = counts[d * 1440:(d + 1) * 1440]
        day[7 * 60:21 * 60] = rng.integers(0, 3000, 14 * 60)     # 14 h of waking wear
    f, daily = H.participant_features(counts, None, np.arange(1, 8, dtype=float))
    assert f["n_valid_days"] == 7 and f["passes_valid_wear_rule"]
    assert 800 <= f["mean_wear_min"] <= 14 * 60
    assert 0 < f["sedentary_fraction"] < 0.1
    assert f["mvpa_min_per_day"] > 0 and np.isnan(f["steps_per_day"])
    assert 0 <= f["interdaily_stability"] <= 1
    assert len(daily) == 7


def test_invalid_counts_not_wear():
    counts = np.zeros(7 * 1440)
    counts[:] = 30000                           # implausible everywhere
    f, _ = H.participant_features(counts, None, np.arange(1, 8, dtype=float))
    assert f["n_valid_days"] == 0 and not f["passes_valid_wear_rule"]


# ------------------------------------------------------------------ data tests

needs_data = pytest.mark.data


def _t(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p)


@needs_data
@pytest.mark.parametrize("name", ["participants__nhanes0306", "participant_clinical_features__nhanes0306",
                                  "participant_labs__nhanes0306", "participant_mortality__nhanes0306",
                                  "participant_medications__nhanes0306", "participant_wearable_features__nhanes0306",
                                  "participant_wearable_daily__nhanes0306"])
def test_tables_contract(name):
    df = _t(name)
    assert set(PROVENANCE_COLUMNS) <= set(df.columns)
    assert df["participant_id"].str.match(r"^nhanes0306:\d+$").all()
    assert (df["data_layer"] == "person").all()
    assert not [c for c in df.columns if GEO_COLUMN_RE.search(c) and c not in PROVENANCE_COLUMNS]


@needs_data
def test_namespaces_do_not_collide_with_2011_2014():
    a = _t("participants__nhanes0306")
    b = _t("participants__nhanes")
    assert set(a["cycle_code"]) == {"C", "D"}
    assert not set(a["seqn"]) & set(b["seqn"])
    assert a["participant_id"].is_unique


@needs_data
def test_hip_features_have_no_mims_or_heart_rate():
    w = _t("participant_wearable_features__nhanes0306")
    bad = [c for c in w.columns if "mims" in c.lower() or "heart" in c.lower() or c.lower().startswith("hr")]
    assert not bad
    ok = w[w["passes_valid_wear_rule"].fillna(False).astype(bool)]
    assert (ok["n_valid_days"] >= 4).all() and ok["monitor_reliable"].all()
    assert ok["sedentary_fraction"].between(0, 1).all()
    assert ok.loc[ok["cycle_code"] == "C", "steps_per_day"].isna().all()      # steps exist only in 2005-2006
    assert ok.loc[ok["cycle_code"] == "D", "steps_per_day"].notna().mean() > 0.99
    assert ok["pctl_within_cohort__mean_daily_counts"].between(0, 1).all()


@needs_data
def test_paxraw_rows_fully_processed():
    import json
    from measure_it.config import interim_dir
    qa = json.loads((interim_dir("nhanes_2003_2006") / "paxraw_qa.json").read_text())["qa"]
    for c in ("C", "D"):
        assert qa[c]["rows_in_file"] == qa[c]["rows_processed"]


@needs_data
def test_crp_present_and_positive():
    c = _t("participant_clinical_features__nhanes0306")
    crp = c.loc[c["age_years"] >= 18, "lab_crp_mg_dl"]
    assert crp.notna().mean() > 0.8 and (crp.dropna() >= 0).all()


@needs_data
def test_replication_results_if_built():
    p = TABLES / "nhanes0306_replication.csv"
    if not p.exists():
        pytest.skip("replication not run")
    r = pd.read_csv(p).dropna(subset=["delta_2003_2006"])
    assert {"functional_limitation", "fair_poor_health", "mortality"} <= set(r["target"])
    assert r["replicated"].isin([True, False]).all()
    perf = pd.read_csv(TABLES / "nhanes0306_model_performance.csv")
    m = perf[(perf.target == "mortality") & (perf.metric == "auroc") & (perf.model == "pen_lr")
             & (perf.feature_set == "wearable") & (perf.split == "cv_5x5")]
    assert m["ci_low"].iloc[0] > 0.5          # positive control visible
