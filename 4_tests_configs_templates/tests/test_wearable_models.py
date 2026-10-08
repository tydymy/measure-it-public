"""Tests for NHANES wearable phenotype modelling (models, signatures, unsupervised, cohort)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from measure_it.config import PROCESSED, TABLES, UNKNOWN
from measure_it.wearables import nhanes_cohort as C
from measure_it.wearables import signatures as S
from measure_it.wearables.models import (
    FEATURE_SETS, bootstrap_weights, calibration, make_model, metric_arrays, rank_metrics,
)

# --------------------------------------------------------------------------- pure functions


def _toy(n=2000, seed=0):
    rng = np.random.default_rng(seed)
    y = (rng.random(n) < 0.2).astype(int)
    s = np.round(rng.random(n) * 0.6 + 0.25 * y, 2)  # rounded -> ties
    return s, y


def test_rank_metrics_match_sklearn_with_ties():
    from sklearn.metrics import average_precision_score, roc_auc_score
    s, y = _toy()
    auc, ap = rank_metrics(s, y)
    assert auc[0] == pytest.approx(roc_auc_score(y, s), abs=1e-12)
    assert ap[0] == pytest.approx(average_precision_score(y, s), abs=1e-12)


def test_bootstrap_weights_equal_expanded_resample():
    from sklearn.metrics import roc_auc_score
    s, y = _toy()
    W = bootstrap_weights(len(y), 2, seed=3)
    assert W.sum(1).tolist() == [len(y), len(y)]
    idx = np.repeat(np.arange(len(y)), W[1].astype(int))
    auc, _ = rank_metrics(s, y, W)
    assert auc[1] == pytest.approx(roc_auc_score(y[idx], s[idx]), abs=1e-12)


def test_calibration_matches_statsmodels():
    import statsmodels.api as sm
    from scipy.special import logit
    s, y = _toy()
    p = np.clip(0.05 + 0.6 * s, 0.01, 0.99)
    a0, b = calibration(p, y)
    lp = logit(p)
    slope = sm.GLM(y, sm.add_constant(lp), family=sm.families.Binomial()).fit().params[1]
    citl = sm.GLM(y, np.ones((len(y), 1)), offset=lp, family=sm.families.Binomial()).fit().params[0]
    assert b[0] == pytest.approx(slope, abs=1e-6)
    assert a0[0] == pytest.approx(citl, abs=1e-6)


def test_perfectly_calibrated_predictions_have_unit_slope():
    rng = np.random.default_rng(1)
    p = rng.uniform(0.02, 0.6, 50_000)
    y = (rng.random(len(p)) < p).astype(int)
    a0, b = calibration(p, y)
    assert abs(a0[0]) < 0.05 and abs(b[0] - 1) < 0.05


def test_metric_arrays_sens_spec_from_calls():
    y = np.array([1, 1, 0, 0, 0])
    pred = np.array([[0.9, 0.2, 0.8, 0.1, 0.1]])
    call = np.array([[1.0, 0.0, 1.0, 0.0, 0.0]])
    m = metric_arrays(pred, call, y, None)
    assert m["sensitivity"][0] == pytest.approx(0.5)
    assert m["specificity"][0] == pytest.approx(2 / 3)


def test_preprocessing_inside_pipeline():
    for name in ("lr", "pen_lr", "hgb"):
        steps = [s for s, _ in make_model(name, 0).steps]
        assert steps == ["impute", "scale", "model"]
        assert make_model(name, 0).named_steps["impute"].add_indicator


def test_feature_sets_exclude_target_defining_items():
    forbidden = set(C.PFQ_ITEMS + C.PHQ_ITEMS + C.SRH_ITEMS) | {
        c for c in C.NEVER_PREDICTORS if not (" " in c)}
    for fs in FEATURE_SETS:
        for t in C.TARGETS:
            cols = set(C.feature_columns(fs, t))
            assert not cols & forbidden, (fs, t, cols & forbidden)
            assert not any(c.startswith(("pfq", "phq", "dlq", "y_")) for c in cols)
    # the proxy is defined by exclusion diagnoses -> the comorbidity count must not predict it
    assert "comorbidity_count" not in C.feature_columns("clinical_wearable", "mecfs_like_proxy")
    assert "comorbidity_count" in C.feature_columns("clinical_wearable", "fatigue")


def test_wearable_set_has_no_heart_rate():
    wear = C.feature_columns("wearable")
    assert "exam_pulse_60s_bpm" not in wear
    assert not any("hr" in c.split("_") or "hrv" in c for c in wear)
    assert "exam_pulse_60s_bpm" in C.feature_columns("clinical")


def test_clock_transform_wraps_at_cut():
    df = pd.DataFrame({f: [1.0, 1.0] for f in C.WEARABLE_FEATURES})
    df["sleep_proxy_midpoint_h"] = [3.0, 11.5]
    df["m10_onset_h"] = [9.0, 19.0]
    out = C.transform_wearable(df)
    assert out["sleep_proxy_midpoint_h"].tolist() == [15.0, 23.5]
    assert out["m10_onset_h"].tolist() == [15.0, 1.0]


def test_signature_estimator_recovers_known_shift():
    rng = np.random.default_rng(7)
    n = 4000
    age = rng.uniform(20, 80, n)
    female = (rng.random(n) < 0.5).astype(float)
    case = (rng.random(n) < 0.2).astype(float)
    x = 0.05 * age + 0.5 * case + rng.normal(0, 1, n)
    r = S.adjusted_std_difference(x, case, age, female)
    sd = x.std(ddof=1)
    assert r["ci_low"] < 0.5 / sd < r["ci_high"]
    assert r["n_cases"] == int(case.sum())


def test_bh_fdr_monotone():
    q = S.bh_fdr(np.array([0.001, 0.01, 0.04, 0.5, np.nan]))
    assert np.isnan(q[-1]) and (np.diff(q[:4]) >= 0).all() and q[0] <= 0.005


# --------------------------------------------------------------------------- query function (synthetic partitions)


def _fake_partition(tmp_path, name, df):
    p = tmp_path / f"{name}.parquet"
    df.to_parquet(p)
    return p


def test_signature_query_unknown_when_no_partitions(monkeypatch):
    monkeypatch.setattr("measure_it.store.partitions", lambda table: [])
    r = S.get_patient_phenotype_signature("ME/CFS")
    assert r["status"] == UNKNOWN and "reason" in r


def test_signature_query_reads_all_partitions_and_alt_schema(monkeypatch, tmp_path):
    a = pd.DataFrame({
        "object_id": ["signature:nhanes|mecfs_like_proxy|mean_daily_mims", "signature:nhanes|mecfs_like_proxy|l5_mims"],
        "dataset_id": "nhanes", "phenotype_id": "mecfs_like_proxy",
        "phenotype_label": C.PROXY_LABEL, "phenotype_definition": "def", "label_basis": "constructed_proxy",
        "condition_id": "me_cfs", "is_proxy": True, "feature": ["mean_daily_mims", "l5_mims"],
        "feature_label": ["a", "b"], "effect_size": [-0.4, 0.05], "ci_low": [-0.6, -0.1], "ci_high": [-0.2, 0.2],
        "p_value": [1e-4, 0.6], "q_value_bh": [2e-4, 0.6], "n_cases": 100, "n_controls": 5000,
        "effect_measure": "x", "caveats": "proxy", "evidence_level": "proxy"})
    b = pd.DataFrame({  # another producer, different column names
        "phenotype": "long_covid_status", "canonical_condition_id": "long_covid", "feature_name": ["rhr_mean"],
        "std_diff": [0.3], "ci_lower": [0.1], "ci_upper": [0.5], "q_value": [0.01], "n_case": [40]})
    pa = _fake_partition(tmp_path, "phenotype_signatures__nhanes", a)
    pb = _fake_partition(tmp_path, "phenotype_signatures__stanford_covid", b)
    monkeypatch.setattr("measure_it.store.partitions", lambda table: [pa, pb])
    r = S.get_patient_phenotype_signature("ME/CFS")
    assert r["status"] == "found"
    e = r["datasets"][0]
    assert e["dataset_id"] == "nhanes" and e["is_proxy"] and "PROXY" in e["phenotype_label"]
    assert [f["feature"] for f in e["top_features"]] == ["mean_daily_mims"]
    assert e["null_results"]["features"] == ["l5_mims"]
    r2 = S.get_patient_phenotype_signature("Long COVID")
    assert r2["status"] == "found" and r2["datasets"][0]["dataset_id"] == "stanford_covid"
    assert r2["datasets"][0]["top_features"][0]["object_id"] == "signature:stanford_covid|long_covid_status|rhr_mean"
    r3 = S.get_patient_phenotype_signature("gastroparesis")
    assert r3["status"] == UNKNOWN and "reason" in r3


def test_signature_query_text_match_ignores_exclusion_variants(monkeypatch, tmp_path):
    """A variant that EXCLUDES arthritis must not answer a query for arthritis; 'health' must not hit 'healthy'."""
    base = dict(dataset_id="nhanes", phenotype_definition="d", label_basis="constructed_proxy", condition_id="me_cfs",
                is_proxy=True, feature="l5_mims", feature_label="x", effect_size=0.3, ci_low=0.1, ci_high=0.5,
                p_value=1e-3, q_value_bh=1e-3, n_cases=80, n_controls=5000, effect_measure="smd", null_value=0.0,
                caveats="c", evidence_level="e")
    rows = [
        {**base, "phenotype_id": "mecfs_like_proxy_excl_arthritis",
         "phenotype_label": C.PROXY_LABEL + " - sensitivity: also excluding any arthritis"},
        {**base, "phenotype_id": "mecfs_like_proxy_vs_healthy_reference",
         "phenotype_label": C.PROXY_LABEL + " - contrast vs healthy reference"},
        {**base, "phenotype_id": "fair_poor_health", "phenotype_label": "Fair or poor self-rated health",
         "condition_id": None, "is_proxy": False},
    ]
    df = pd.DataFrame(rows)
    df["object_id"] = "signature:nhanes|" + df["phenotype_id"] + "|" + df["feature"]
    p = _fake_partition(tmp_path, "phenotype_signatures__nhanes", df)
    monkeypatch.setattr("measure_it.store.partitions", lambda table: [p])
    assert S.get_patient_phenotype_signature("arthritis")["status"] == UNKNOWN
    r = S.get_patient_phenotype_signature("health")
    assert r["status"] == "found" and [e["phenotype_id"] for e in r["datasets"]] == ["fair_poor_health"]


def test_signature_query_separates_model_rows_and_untested_rows(monkeypatch, tmp_path):
    """Model-level AUROC rows are not wearable-feature signatures; rows without q are not null results;
    a missing null value is UNKNOWN, never assumed to be 0."""
    df = pd.DataFrame({
        "phenotype_id": "long_covid_self_report", "canonical_condition_id": "long_covid",
        "feature_name": ["rhr_mean", "model_auroc:demographics:l2_logistic", "recovery_days_km_median"],
        "std_diff": [0.1, 0.70, 6.0], "ci_lower": [-0.2, 0.6, 3.0], "ci_upper": [0.4, 0.8, 9.0],
        "q_value": [0.6, 0.01, np.nan], "effect_unit": ["Hedges g", "AUROC (CV)", "days"],
        "null_value": [np.nan, 0.5, np.nan], "n_case": 31})
    p = _fake_partition(tmp_path, "phenotype_signatures__stanford", df)
    monkeypatch.setattr("measure_it.store.partitions", lambda table: [p])
    r = S.get_patient_phenotype_signature("Long COVID")
    e = r["datasets"][0]
    assert e["top_features"] == []                                  # the demographics AUROC is not a top feature
    assert e["null_results"]["phenotype_level_null"] is True
    assert e["null_results"]["features"] == ["rhr_mean"]
    assert [x["feature"] for x in e["model_level_results"]["rows"]] == ["model_auroc:demographics:l2_logistic"]
    assert [x["feature"] for x in e["descriptive_rows_not_tested"]] == ["recovery_days_km_median"]
    assert e["descriptive_rows_not_tested"][0]["null_value"] == UNKNOWN


def test_holm_boot_adjustment():
    from measure_it.wearables.nhanes_report import SUP, holm_boot
    p = [0.0, 0.05, 0.712, 0.034, 0.096, 0.24]
    inc = pd.DataFrame({"target": SUP, "comparison": "test2_clinical_wearable_vs_clinical", "model": "pen_lr",
                        "metric": "auroc", "split": "cv_5x5", "p_boot_two_sided": p, "n_boot": 1000})
    h = holm_boot(inc, "test2_clinical_wearable_vs_clinical")
    assert h["functional_limitation"] == (pytest.approx(0.002), pytest.approx(0.012))   # floored at 2/B, x6
    assert h["depression_phq9"][1] == pytest.approx(max(0.012, 5 * 0.034))
    assert h["fatigue"][1] == pytest.approx(max(0.17, 4 * 0.05))
    vals = [h[t][1] for t in sorted(h, key=lambda t: h[t][0])]
    assert all(np.diff(vals) >= -1e-12) and max(vals) <= 1.0


# --------------------------------------------------------------------------- data tests (processed outputs)

needs_data = pytest.mark.data


@needs_data
def test_cohort_counts_match_audit():
    co = C.load_cohort()
    steps = {s["step"]: s["n"] for s in co.flow}
    assert steps["pass valid-wear rule (>= 4 valid days)"] == 12955
    assert steps["aged >= 18 at screening"] == 8954
    f = co.frame
    assert f["participant_id"].is_unique and (f["age_years"] >= 18).all()
    assert f["participant_id"].str.startswith("nhanes:").all()
    # rx groups are 2013-2014 only
    assert f.loc[f["cycle_code"] == "G", "y_rx_fibromyalgia"].isna().all()
    assert int(f["y_rx_fibromyalgia"].sum()) > 0
    # the proxy is a subset of fatigue AND limitation
    prox = f["y_mecfs_like_proxy"] == 1
    assert (f.loc[prox, "y_fatigue"] == 1).all() and (f.loc[prox, "y_functional_limitation"] == 1).all()
    assert (f.loc[prox, "proxy_any_exclusion_dx"] == 0).all()


@needs_data
def test_phenotype_signatures_table_contract():
    p = PROCESSED / "phenotype_signatures__nhanes.parquet"
    if not p.exists():
        pytest.skip("phenotype_signatures__nhanes not built")
    from measure_it.provenance import PROVENANCE_COLUMNS, check_provenance
    s = pd.read_parquet(p)
    check_provenance(s, "phenotype_signatures__nhanes")
    assert set(PROVENANCE_COLUMNS) <= set(s.columns)
    assert s["object_id"].is_unique
    assert s["object_id"].str.match(r"^signature:nhanes\|[a-z0-9_]+\|[a-z0-9_]+$").all()
    # condition ids only where allowed; me_cfs only on labelled proxy rows
    allowed = {"me_cfs", "fibromyalgia", "migraine", "ibs"}
    assert set(s["condition_id"].dropna()) <= allowed
    me = s[s["condition_id"] == "me_cfs"]
    assert len(me) and me["phenotype_label"].str.contains("ME/CFS-LIKE PROXY", case=False).all()
    assert me["is_proxy"].all()
    assert not s["phenotype_label"].str.contains("diagnosed", case=False).any()
    assert (s["ci_low"] <= s["effect_size"]).all() and (s["effect_size"] <= s["ci_high"]).all()
    assert (s["null_value"] == 0.0).all()
    # provenance cites only the cycles a phenotype used: rx groups are 2013-2014 only
    rx = s[s["phenotype_id"].str.startswith("rx_")]
    assert rx["source_version"].str.contains("2013-2014").all()
    assert not rx["source_version"].str.contains("2011-2012").any()
    both = s[s["phenotype_id"] == "fatigue"]
    assert both["source_version"].str.contains("2011-2012").all() and both["source_version"].str.contains("2013-2014").all()


@needs_data
def test_cluster_assignments_contract():
    p = PROCESSED / "participant_cluster_assignments__nhanes.parquet"
    if not p.exists():
        pytest.skip("participant_cluster_assignments__nhanes not built")
    from measure_it.provenance import check_provenance
    a = pd.read_parquet(p)
    check_provenance(a, "participant_cluster_assignments__nhanes")
    assert a.groupby("method")["participant_id"].nunique().min() == a.groupby("method").size().min()
    assert a["cluster_note"].str.contains("not a").all()
    # each row cites its own cycle's wearable release
    assert a.loc[a["cycle_code"] == "H", "source_version"].str.contains("2013-2014").all()
    assert a.loc[a["cycle_code"] == "G", "source_version"].str.contains("2011-2012").all()


@needs_data
def test_positive_control_mortality_is_visible():
    p = TABLES / "nhanes_model_performance.csv"
    if not p.exists():
        pytest.skip("model results not built")
    perf = pd.read_csv(p)
    r = perf[(perf.target == "mortality") & (perf.feature_set == "wearable") & (perf.model == "pen_lr")
             & (perf.metric == "auroc") & (perf.split == "cv_5x5")]
    assert len(r) == 1 and r["ci_low"].iloc[0] > 0.6


@needs_data
def test_permutation_null_centred_on_chance():
    p = TABLES / "nhanes_permutation_null.csv"
    if not p.exists():
        pytest.skip("permutation results not built")
    pn = pd.read_csv(p)
    assert (pn["n_permutations"] >= 200).all()
    assert pn["null_mean"].between(0.45, 0.55).all()


@needs_data
def test_query_returns_proxy_for_mecfs():
    if not (PROCESSED / "phenotype_signatures__nhanes.parquet").exists():
        pytest.skip("phenotype_signatures__nhanes not built")
    r = S.get_patient_phenotype_signature("ME/CFS")
    assert r["status"] == "found"
    nh = [e for e in r["datasets"] if e["dataset_id"] == "nhanes"]
    assert nh and all(e["is_proxy"] for e in nh)
    assert all("PROXY" in e["phenotype_label"] for e in nh)


@pytest.mark.data
def test_sleep_self_report_signatures_present():
    """Integration item A8: SLQ050 / SLQ060 signatures exist, are signature-only and flow into Phase 3 signal rows."""
    from measure_it.store import read_table, table_exists
    sig = read_table("phenotype_signatures__nhanes")
    for pid in ("slq050_told_trouble_sleeping", "slq060_told_sleep_disorder"):
        g = sig[sig["phenotype_id"] == pid]
        assert len(g) == sig["feature"].nunique() and (g["n_cases"] > 100).all()
        assert (g["label_basis"] == "self_report_told_by_health_professional").all()
    if table_exists("measurement_phenotype_signal"):
        s = read_table("measurement_phenotype_signal")
        assert {"slq050_told_trouble_sleeping", "slq060_told_sleep_disorder"} <= set(s["label_id"])
        assert (s.loc[s["label_id"].str.startswith("slq"), "label_class"] == "self_reported_symptom_or_function").all()
