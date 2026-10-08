"""Tests for measure_it.wearables.stanford_models (Stanford wearable analyses).

Unit tests use synthetic data; tests marked `data` check the outputs written by
`uv run python -m measure_it.wearables.stanford_models`.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

from measure_it.config import RESULTS, TABLES
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table, table_exists
from measure_it.wearables import stanford_models as sm


# ---------------------------------------------------------------------------
# pure statistics
# ---------------------------------------------------------------------------
def test_auroc_rows_matches_sklearn_with_ties():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 60)
    S = np.round(rng.normal(size=(3, 60)), 1)  # rounding creates ties
    got = sm.auroc_rows(y, S)
    want = [roc_auc_score(y, s) for s in S]
    assert np.allclose(got, want)


def test_auroc_rows_single_class_is_nan():
    assert np.isnan(sm.auroc_rows(np.ones(5), np.arange(5.0))[0])


def test_average_precision_matches_sklearn_with_ties():
    rng = np.random.default_rng(1)
    for _ in range(20):
        y = rng.integers(0, 2, 40)
        if y.sum() == 0:
            continue
        s = np.round(rng.random(40), 1)
        assert sm.average_precision(y, s) == pytest.approx(average_precision_score(y, s))


def test_hedges_g_known_value():
    x1 = np.array([2.0, 4.0, 6.0])
    x0 = np.array([1.0, 2.0, 3.0])
    sp = np.sqrt((2 * 4.0 + 2 * 1.0) / 4)
    j = 1 - 3 / (4 * 6 - 9)
    assert sm.hedges_g(x1, x0) == pytest.approx((4 - 2) / sp * j)


def test_bh_qvalues_skips_nan_and_is_monotone():
    p = pd.Series([0.01, np.nan, 0.04, 0.03, 0.5])
    q = sm.bh_qvalues(p)
    assert np.isnan(q.iloc[1])
    assert q.iloc[0] == pytest.approx(0.04)          # 0.01 * 4 / 1
    assert (q.dropna() >= p.dropna()).all()


def test_wilson_ci_brackets_estimate():
    lo, hi = sm.wilson_ci(7, 10)
    assert 0 < lo < 0.7 < hi < 1
    assert np.isnan(sm.wilson_ci(0, 0)[0])


def test_signflip_pvalue_behaviour():
    assert sm.signflip_pvalue(np.full(30, 1.0), n_perm=2000) < 0.01
    d = np.r_[np.full(15, 1.0), np.full(15, -1.0)]
    assert sm.signflip_pvalue(d, n_perm=2000) > 0.5


def test_signature_evidence_level_rules():
    assert sm.signature_evidence_level(np.nan, np.nan, 0.1, 0.2, 0.0) == "descriptive"
    assert sm.signature_evidence_level(0.001, 0.01, 0.1, 0.4, 0.0) == "supported_single_dataset"
    assert sm.signature_evidence_level(0.001, 0.01, -0.1, 0.4, 0.0) == "nominal_single_dataset"
    assert sm.signature_evidence_level(0.04, 0.2, 0.1, 0.4, 0.0) == "nominal_single_dataset"
    assert sm.signature_evidence_level(0.3, 0.5, 0.45, 0.6, 0.5) == "null_single_dataset"


def test_welch_t_matches_scipy():
    from scipy import stats
    rng = np.random.default_rng(3)
    x = np.r_[rng.normal(1, 3, 20), rng.normal(0, 1, 60)]
    g = np.r_[np.ones(20, bool), np.zeros(60, bool)]
    want = stats.ttest_ind(x[g], x[~g], equal_var=False).statistic
    assert sm.welch_t(x, g[None, :])[0] == pytest.approx(want)
    # one row per permutation
    G = np.array([rng.permutation(g) for _ in range(5)])
    assert sm.welch_t(x, G).shape == (5,)


def test_km_summary_drops_missing_rows_instead_of_casting_nan():
    t = np.array([1.0, 2.0, 3.0, np.nan, 5.0])
    e = np.array([1.0, 1.0, 0.0, np.nan, 1.0])
    s, tab = sm.km_summary(t, e)
    assert s["n"] == 4 and s["n_dropped_no_followup"] == 1
    assert s["events"] == 3 and s["censored"] == 1  # a NaN cast to int would make this sum absurd


def test_km_weighted_median_equals_unweighted_with_equal_weights():
    rng = np.random.default_rng(4)
    t = rng.integers(0, 15, 60).astype(float)
    e = (rng.random(60) < 0.7).astype(float)
    s, _ = sm.km_summary(t, e)
    assert sm.km_weighted_median(t, e, np.ones(60)) == pytest.approx(s["median"])
    assert sm.km_weighted_median(t, e, np.full(60, 0.25)) == pytest.approx(s["median"])


def test_km_participant_weighted_gives_each_participant_weight_one():
    # participant a: 9 onsets recovering at day 1; participant b: 1 onset at day 10; c: 1 onset at day 12.
    # Unweighted, the median is day 1; participant-weighted (a, b, c each weight 1) it is day 10.
    df = pd.DataFrame({"participant_id": ["a"] * 9 + ["b", "c"], "t": [1.0] * 9 + [10.0, 12.0], "e": 1.0})
    s = sm.km_participant_weighted(df, "t", "e", n_boot=50, seed=0)
    assert s["n"] == 3 and s["n_onsets"] == 11
    assert s["median"] == 10.0


# ---------------------------------------------------------------------------
# cross-validation machinery
# ---------------------------------------------------------------------------
def test_cv_splits_do_not_depend_on_features_and_cover_everyone():
    y = np.r_[np.ones(20), np.zeros(60)].astype(int)
    a = sm.cv_splits(y, 3)
    b = sm.cv_splits(y, 3)
    assert all((x[1] == z[1]).all() for x, z in zip(a, b))
    test_idx = np.sort(np.concatenate([te for _, te in a]))
    assert (test_idx == np.arange(80)).all()


def test_grouped_splits_never_share_a_participant():
    y = np.tile([1, 0], 40)
    groups = np.repeat(np.arange(40), 2)
    assert sm.grouped_fold_overlap(y, groups, n_repeats=3) == 0


def test_swap_labels_within_groups_keeps_one_positive_per_pair():
    y = np.tile([1, 0], 50)
    groups = np.repeat(np.arange(50), 2)
    yp = sm.swap_labels_within_groups(y, groups, np.random.default_rng(0))
    assert (pd.Series(yp).groupby(groups).sum() == 1).all()
    assert 0 < (yp != y).sum() < len(y)


def test_cv_oof_random_labels_near_chance_and_signal_detected():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(150, 3))
    y_null = (rng.random(150) < 0.3).astype(int)
    P, B = sm.cv_oof(X, y_null, "l2_logistic_fixed", n_repeats=2)
    assert np.isfinite(P).all() and set(np.unique(B)) <= {0.0, 1.0}
    assert abs(sm.auroc_rows(y_null, P).mean() - 0.5) < 0.15
    y_sig = (X[:, 0] + 0.5 * rng.normal(size=150) > 0.5).astype(int)
    P, _ = sm.cv_oof(X, y_sig, "l2_logistic_fixed", n_repeats=2)
    assert sm.auroc_rows(y_sig, P).mean() > 0.8


# ---------------------------------------------------------------------------
# placebo onsets on synthetic daily data
# ---------------------------------------------------------------------------
def _synthetic_daily(first: int = -120, last: int = 60, spike=(0, 3), seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    days = np.arange(first, last + 1).astype(float)
    rhr = 60 + rng.normal(0, 1.5, days.size)
    steps = 8000 + rng.normal(0, 800, days.size)
    ill = (days >= spike[0]) & (days <= spike[1])
    rhr[ill] += 8
    steps[ill] -= 4000
    ref = pd.Timestamp("2024-06-01")
    return pd.DataFrame({"participant_id": "x:1", "date": ref + pd.to_timedelta(days, unit="D"), "days_from_onset": days,
                         "rhr_mean": rhr, "steps": steps, "valid_wear_day": True})


def test_deviations_at_uses_only_days_before_the_shifted_baseline_end():
    d = _synthetic_daily()
    dev, base = sm.deviations_at(d, -50)
    assert base.ok
    bdays = dev.loc[dev.is_baseline_day, "days_from_onset"]
    assert bdays.max() == -15 and bdays.min() == -42
    # baseline relative to the true reference is days -92..-65
    assert (dev.loc[dev.is_baseline_day, "days_from_onset"] - 50).max() == -65


def test_placebo_windows_never_reach_the_true_window_and_respect_blocked_dates():
    d = _synthetic_daily()
    res = sm.process_acute_participant("x:1", d, np.array([-60.0]))
    assert res["baseline_ok"] and res["true"]["rhr_acute_elevation"]
    offs = np.array([st["placebo_offset"] for st in res["placebo"]])
    assert offs.size > 0
    assert offs.max() <= sm.PLACEBO_MAX_OFFSET
    assert all(st["max_day_used_rel_true"] <= -15 for st in res["placebo"])
    # blocked date -60 must not fall inside any placebo acute window [P-7, P+14]
    assert not np.any((offs - 7 <= -60) & (offs + 14 >= -60))
    # the synthetic spike is at the true onset only: placebo windows should mostly be quiet
    mean_true = res["true"]["rhr_mean_z_post0_14"]
    mean_plc = np.nanmean([st["rhr_mean_z_post0_14"] for st in res["placebo"]])
    assert mean_true > mean_plc + 0.5


def test_wear_stats_counts_only_days_inside_the_record_span():
    d = _synthetic_daily(first=-30)  # record starts inside the [-42, -15] baseline window
    d["hr_minutes"] = 900.0
    dev, _ = sm.deviations_at(d, 0)
    w = sm.wear_stats(dev)
    assert w["wear_valid_day_frac_base"] == pytest.approx(1.0)  # 16 recorded days of 16 in-span days, not 16/28
    assert w["wear_valid_day_frac_change"] == pytest.approx(0.0)
    assert w["hr_minutes_per_valid_day_change"] == pytest.approx(0.0)


def test_placebo_requires_its_own_baseline():
    d = _synthetic_daily(first=-45)  # too short for any placebo baseline + P <= -29
    res = sm.process_acute_participant("x:1", d, None)
    assert res["placebo"] == []


# ---------------------------------------------------------------------------
# processed outputs
# ---------------------------------------------------------------------------
REQUIRED_SIG_COLS = ["dataset_id", "phenotype_id", "phenotype_label", "condition_id", "feature", "effect_size", "effect_unit",
                     "ci_low", "ci_high", "p_value", "q_value", "n_cases", "n_controls", "label_basis", "object_id"]
# columns read by measure_it.wearables.signatures.get_patient_phenotype_signature (UNKNOWN when absent)
QUERY_SIG_COLS = ["feature_label", "effect_measure", "phenotype_definition", "is_proxy", "caveats"]


@pytest.mark.data
def test_signature_table_schema_and_ids():
    if not table_exists(sm.SIGNATURE_TABLE):
        pytest.skip("run the stanford_models pipeline first")
    sig = read_table(sm.SIGNATURE_TABLE)
    for c in REQUIRED_SIG_COLS + PROVENANCE_COLUMNS:
        assert c in sig.columns, c
    assert sig.object_id.is_unique
    pat = re.compile(r"^signature:[^|]+\|[^|]+\|[^|]+$")
    assert sig.object_id.map(lambda s: bool(pat.match(s))).all()
    assert (sig.object_id == "signature:" + sig.dataset_id + "|" + sig.phenotype_id + "|" + sig.feature).all()
    uw = sig[sig.dataset_id == sm.UWAKWE]
    ac = sig[sig.dataset_id == sm.ACUTE_POOLED]
    assert len(uw) > 0 and len(ac) > 0 and len(uw) + len(ac) == len(sig)
    assert (uw.condition_id == "long_covid").all()
    assert ac.condition_id.isna().all()
    assert (ac.phenotype_label == "post-infection recovery trajectory (acute COVID-19 cohort; not a Long COVID label)").all()
    assert not sig.feature.str.contains("symptom", case=False).any()  # leakage reference never becomes a signature
    assert set(sig.evidence_level) <= {"supported_single_dataset", "nominal_single_dataset", "null_single_dataset", "descriptive"}
    assert (sig.data_layer == "person").all() and (sig.source_geographic_resolution == "none").all()
    ok = sig.ci_low.notna()
    assert (sig.loc[ok, "ci_low"] <= sig.loc[ok, "ci_high"]).all()
    for c in QUERY_SIG_COLS:
        assert c in sig.columns, c
    assert sig.caveats.str.len().gt(0).all() and sig.feature_label.notna().all() and sig.effect_measure.notna().all()
    assert ac.caveats.str.contains("NOT a Long COVID label").all()
    # only wearable-derived quantities are signature rows
    assert not sig.feature.str.startswith(("model_auroc:demographics:", "model_auroc:demographics_plus_wearable:")).any()


@pytest.mark.data
def test_signature_query_does_not_surface_non_wearable_rows_for_long_covid():
    if not table_exists(sm.SIGNATURE_TABLE):
        pytest.skip("run the stanford_models pipeline first")
    from measure_it.wearables.signatures import get_patient_phenotype_signature
    res = get_patient_phenotype_signature("long_covid")
    if res.get("status") != "found":
        pytest.skip(f"signature query unavailable: {res.get('reason')}")
    uw = [e for e in res["datasets"] if e["dataset_id"] == sm.UWAKWE]
    assert uw, "Uwakwe partition not returned for long_covid"
    for e in uw:
        feats = [t["feature"] for t in e["top_features"]]
        assert not any("demographics" in f and "increment" not in f for f in feats), feats
        assert e["caveats"] != "UNKNOWN / NOT AVAILABLE"
    for e in res["datasets"]:
        if e["dataset_id"] == sm.ACUTE_POOLED:  # may be returned by phenotype-text match; must carry the caveat
            assert "NOT a Long COVID label" in str(e["caveats"])


@pytest.mark.data
def test_recovery_tables_include_peak_timed_and_weighted_placebo_rows():
    p = TABLES / "stanford_acute_recovery_km.csv"
    if not p.exists():
        pytest.skip("run the stanford_models pipeline first")
    km = pd.read_csv(p)
    assert "time_origin" in km and (km.time_origin == "acute peak day").sum() >= 3
    assert km.analysis.str.contains("participant-weighted").sum() == 2
    assert (km.n_dropped_no_followup.fillna(0) >= 0).all()
    assert (km.events + km.censored == km.n_onsets.fillna(km.n)).all()


@pytest.mark.data
def test_leakage_and_consistency_checks_all_pass():
    p = TABLES / "stanford_temporal_leakage_checks.csv"
    if not p.exists():
        pytest.skip("run the stanford_models pipeline first")
    leak = pd.read_csv(p)
    assert len(leak) >= 8
    assert (leak.violations == 0).all(), leak[leak.violations > 0]


@pytest.mark.data
def test_recomputed_true_onset_z_equals_ingestion_table():
    if not table_exists("participant_wearable_daily__stanford_covid"):
        pytest.skip("ingestion outputs missing")
    part, daily, _ = sm.load_acute()
    pids = part[part.baseline_ok.astype(bool)].participant_id.head(12)
    for pid in pids:
        dp = daily[daily.participant_id == pid]
        dev, _ = sm.deviations_at(dp, 0)
        a, b = dev.rhr_z.values, dp.rhr_z.values
        assert np.allclose(a, b, equal_nan=True), pid
        assert np.allclose(dev.steps_z.values, dp.steps_z.values, equal_nan=True), pid


@pytest.mark.data
def test_uwakwe_frame_matches_published_cohort():
    if not table_exists("participants__stanford_covid"):
        pytest.skip("ingestion outputs missing")
    df = sm.build_uwakwe_frame()
    assert len(df) == 126 and int(df.y.sum()) == 31
    assert int(df.loc[df.y == 1, "female"].sum()) == 23 and int(df.loc[df.y == 0, "female"].sum()) == 42
    for fs, cols in sm.FEATURE_SETS.items():
        assert set(cols) <= set(df.columns), fs
    assert not set(sm.SYMPTOM_FEATURES) & set(sm.FEATURE_SETS["demographics_plus_wearable"])


@pytest.mark.data
def test_model_tables_and_permutation_nulls():
    p = TABLES / "stanford_uwakwe_model_performance.csv"
    if not p.exists():
        pytest.skip("run the stanford_models pipeline first")
    perf = pd.read_csv(p)
    assert len(perf) == len(sm.FEATURE_SETS) * len(sm.MODELS)
    prim = perf[perf.feature_set.isin(sm.PRIMARY_SETS)]
    assert prim.label_perm_p.notna().all()
    assert ((perf.auroc_ci_low <= perf.auroc) & (perf.auroc <= perf.auroc_ci_high)).all()
    nulls = pd.read_csv(TABLES / "stanford_uwakwe_permutation_nulls.csv")
    lab = nulls[nulls.null_type == "label_permutation"]
    assert (lab.n_perm >= 1000).all()
    assert lab.null_mean.between(0.4, 0.6).all()  # shuffled labels -> no structure
    blk = nulls[nulls.null_type == "wearable_block_permutation"]
    assert len(blk) == len(sm.MODELS) and (blk.n_perm >= 1000).all()


@pytest.mark.data
def test_results_markdown_language():
    p = RESULTS / "WEARABLE_STANFORD_RESULTS.md"
    if not p.exists():
        pytest.skip("results not rendered yet")
    text = p.read_text().lower()
    for banned in ("proves", "diagnosed by ai", "patient has", "best clinic", "definitive biomarker", "confirms mechanism"):
        assert banned not in text, banned
    assert "not a long covid label" in text or "not long covid" in text
