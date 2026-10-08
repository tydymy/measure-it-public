"""Tests for the Stanford Snyder Lab COVID-19 / Long COVID wearable ingestion.

Unit tests exercise the pure feature functions on synthetic data; tests marked
`data` check the processed outputs written by
`uv run python -m measure_it.ingestion.stanford_wearables`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from measure_it.config import TABLES
from measure_it.ingestion import stanford_wearables as sw
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table, table_exists
from measure_it.wearables import stanford_features as sf


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _minutes(start: str, n: int) -> pd.DatetimeIndex:
    return pd.date_range(start, periods=n, freq="min")


def _daily_from_levels(levels: dict[int, float], ref="2024-03-01") -> pd.DataFrame:
    ref = pd.Timestamp(ref)
    rows = [{"date": ref + pd.Timedelta(days=d), "rhr_mean": v, "valid_wear_day": True, "steps": 5000.0,
             "hr_minutes": 1200, "rhr_minutes": 600, "rhr_night_mean": v, "hr_mean_all": v + 5}
            for d, v in sorted(levels.items())]
    return sf.add_reference_days(pd.DataFrame(rows), ref)


# ---------------------------------------------------------------------------
# minute-level preprocessing
# ---------------------------------------------------------------------------
def test_minute_median_hr_takes_median_and_drops_header_rows():
    t = pd.Series(["2020-01-01 00:00:05", "datetime", "2020-01-01 00:00:35", "2020-01-01 00:00:50", "2020-01-01 00:01:10"])
    h = pd.Series(["60", "heartrate", "90", "70", "0"])  # 0 bpm is not a valid sample
    out = sf.minute_median_hr(t, h)
    assert list(out.index) == [pd.Timestamp("2020-01-01 00:00")]
    assert out.iloc[0] == 70.0


def test_resting_mask_excludes_step_minute_and_next_ten_within_day():
    idx = _minutes("2020-01-01 09:55", 30)
    hr = pd.Series(60.0, index=idx)
    steps = pd.Series([12.0], index=[pd.Timestamp("2020-01-01 10:00")])
    m = sf.resting_mask(hr, steps)
    assert m[pd.Timestamp("2020-01-01 09:59")]
    for k in range(11):
        assert not m[pd.Timestamp("2020-01-01 10:00") + pd.Timedelta(minutes=k)]
    assert m[pd.Timestamp("2020-01-01 10:11")]


def test_resting_mask_does_not_cross_midnight():
    idx = _minutes("2020-01-01 23:50", 20)
    hr = pd.Series(60.0, index=idx)
    steps = pd.Series([5.0], index=[pd.Timestamp("2020-01-01 23:58")])
    m = sf.resting_mask(hr, steps)
    assert not m[pd.Timestamp("2020-01-01 23:59")]
    assert m[pd.Timestamp("2020-01-02 00:01")]  # the R implementation masks within a day column


def test_step_intervals_are_spread_over_covered_minutes():
    out = sf.step_minutes_from_intervals(pd.Series(["2020-01-01 00:00:00", "start_datetime"]),
                                         pd.Series(["2020-01-01 00:02:59", "end_datetime"]), pd.Series(["30", "steps"]))
    assert len(out) == 3 and np.isclose(out.sum(), 30)


def test_daily_summary_wear_threshold_and_steps_on_worn_days():
    worn = pd.Series(70.0, index=_minutes("2020-01-01 00:00", 700))
    unworn = pd.Series(70.0, index=_minutes("2020-01-02 00:00", 100))
    hr = pd.concat([worn, unworn])
    steps = pd.Series([10.0], index=[pd.Timestamp("2020-01-01 10:00")])
    d = sf.daily_summary(hr, steps, min_wear_minutes=600, min_rhr_minutes=60).set_index("date")
    assert bool(d.loc["2020-01-01", "valid_wear_day"]) and not bool(d.loc["2020-01-02", "valid_wear_day"])
    assert d.loc["2020-01-01", "steps"] == 10
    assert np.isnan(d.loc["2020-01-02", "steps"])
    assert d.loc["2020-01-01", "rhr_mean"] == 70.0 and d.loc["2020-01-01", "rhr_minutes"] == 700 - 11


def test_worn_day_without_any_recorded_step_gets_no_rhr_and_no_steps():
    """A worn day whose step stream is absent or all zeros is a dropout: the step filter
    would call every (including active) minute resting. Regression for AJWW3IY days -8..-4
    (rhr_mean == mean of all HR, rhr_z ~ +3, steps 0)."""
    idx = _minutes("2020-01-01 00:00", 3 * 1440)
    hr = pd.Series(np.where(idx.hour.isin([12, 13]), 110.0, 60.0), index=idx)  # active midday hours
    # day 1: steps at noon; day 2: explicit zero rows only (Mishra style); day 3: no step rows at all
    steps = pd.concat([
        pd.Series(50.0, index=_minutes("2020-01-01 12:00", 120)),
        pd.Series(0.0, index=_minutes("2020-01-02 00:00", 1440)),
    ])
    d = sf.daily_summary(hr, steps, min_wear_minutes=600, min_rhr_minutes=60).set_index("date")
    assert d.loc["2020-01-01", "step_stream_present"] and abs(d.loc["2020-01-01", "rhr_mean"] - 60.0) < 1e-9
    for day in ("2020-01-02", "2020-01-03"):
        assert not d.loc[day, "step_stream_present"]
        assert np.isnan(d.loc[day, "rhr_mean"]) and np.isnan(d.loc[day, "rhr_night_mean"]) and np.isnan(d.loc[day, "steps"])
        assert d.loc[day, "rhr_minutes"] == 0 and bool(d.loc[day, "valid_wear_day"])
    # idempotent, so it can be re-applied to cached daily tables
    again = sf.mask_stepless_days(d.reset_index()).set_index("date")
    pd.testing.assert_frame_equal(again, d)
    # an unmasked cached core is corrected by derive_features; pre-masked HR is left alone
    raw = d.reset_index().copy()
    raw.loc[raw.date == pd.Timestamp("2020-01-02"), ["rhr_mean", "steps"]] = [80.0, 0.0]
    core = {"daily": raw.drop(columns="step_stream_present"), "alarms": pd.DataFrame(columns=["alarm_time"]), "cusum_evaluable": False}
    dd, _ = sw.derive_features(core, None, "fitbit_minute")
    assert dd.set_index("date").loc["2020-01-02", ["rhr_mean", "steps"]].isna().all()
    dp, _ = sw.derive_features(core, None, "preprocessed_resting_hr")
    assert dp.set_index("date").loc["2020-01-02", "rhr_mean"] == 80.0


def test_no_step_data_means_no_rhr():
    hr = pd.Series(70.0, index=_minutes("2020-01-01 00:00", 1440 * 30))
    core = sw.wearable_core(hr, None, "fitbit_minute")
    assert not core["daily"].rhr_computable.any()
    assert core["daily"].rhr_mean.isna().all() and core["daily"].hr_mean_all.notna().all()
    assert not core["cusum_evaluable"]
    pre = sw.wearable_core(hr, None, "preprocessed_resting_hr", run_cusum=False)  # authors' pre-masked RHR
    assert pre["daily"].rhr_mean.notna().all()


# ---------------------------------------------------------------------------
# baseline / deviations / recovery
# ---------------------------------------------------------------------------
def test_baseline_window_must_precede_reference():
    with pytest.raises(ValueError):
        sf.baseline_stats(pd.DataFrame({"days_from_onset": [0]}), window=(-10, 0))


def test_baseline_ignores_post_reference_data():
    levels = {d: 60.0 + (d % 2) for d in range(-42, -14)}
    levels.update({d: 200.0 for d in range(-14, 30)})  # post-window data must not leak into the baseline
    d = _daily_from_levels(levels)
    b = sf.baseline_stats(d)
    assert b.ok and b.n_days == 28
    assert abs(b.rhr_mean - 60.5) < 1e-9


def test_recovery_event_and_censoring():
    base = {d: 60.0 + (1 if d % 2 else -1) for d in range(-42, -14)}  # mean 60, sd ~1.02
    elevated = {d: 70.0 for d in range(0, 5)}
    back = {d: 60.0 for d in range(5, 12)}
    d = sf.add_deviations(_daily_from_levels({**base, **elevated, **back}), sf.baseline_stats(_daily_from_levels(base)))
    r = sf.recovery_trajectory(d)
    assert r["rhr_acute_elevation"] and r["rhr_acute_peak_day"] == 0
    assert r["rhr_recovery_event"] == 1 and r["days_to_rhr_recovery"] == 5 and r["rhr_recovery_days_from_peak"] == 5
    assert r["rhr_elevated_days_post"] == 5
    d2 = sf.add_deviations(_daily_from_levels({**base, **{k: 70.0 for k in range(0, 20)}}), sf.baseline_stats(_daily_from_levels(base)))
    r2 = sf.recovery_trajectory(d2)
    assert r2["rhr_recovery_event"] == 0 and r2["days_to_rhr_recovery"] == r2["last_followup_day"] == 19


def test_recovery_not_applicable_without_acute_elevation():
    base = {d: 60.0 + (1 if d % 2 else -1) for d in range(-42, -14)}
    flat = {d: 60.5 for d in range(-7, 30)}
    d = sf.add_deviations(_daily_from_levels({**base, **flat}), sf.baseline_stats(_daily_from_levels(base)))
    r = sf.recovery_trajectory(d)
    assert r["rhr_acute_elevation"] is False
    assert np.isnan(r["days_to_rhr_recovery"]) and np.isnan(r["rhr_recovery_event"])
    assert abs(r["rhr_mean_z_post0_14"] - 0.5 / d.baseline_rhr_sd.iloc[0]) < 1e-9


def test_recovery_undefined_without_baseline():
    d = _daily_from_levels({d: 60.0 for d in range(-5, 10)})
    d = sf.add_deviations(d, sf.baseline_stats(d))
    assert np.isnan(sf.recovery_trajectory(d)["days_to_rhr_recovery"])


# ---------------------------------------------------------------------------
# CuSum port
# ---------------------------------------------------------------------------
def _synthetic_hr(days: int, elevate_from: int | None, seed: int = 1) -> pd.Series:
    rng = np.random.default_rng(seed)
    idx = _minutes("2021-01-01", days * 1440)
    hod = idx.hour.to_numpy()
    base = 60 + 5 * np.sin(2 * np.pi * hod / 24) + rng.normal(0, 2, idx.size)
    if elevate_from is not None:
        base = base + np.where(np.arange(idx.size) >= elevate_from * 1440, 12.0, 0.0)
    return pd.Series(base, index=idx)


def test_cusum_needs_more_than_28_days():
    assert sf.cusum_online_stats(_synthetic_hr(28, None), None) is None


def test_cusum_alarms_on_sustained_elevation():
    hr = _synthetic_hr(40, elevate_from=33)
    al = sf.cusum_alarms(sf.cusum_online_stats(hr, None))
    assert len(al) >= 1
    first = al.alarm_time.min()
    assert pd.Timestamp("2021-02-03") <= first <= pd.Timestamp("2021-02-05")


def test_select_infection_alarm_window_and_nearest():
    al = pd.DataFrame({"alarm_time": pd.to_datetime(["2024-01-01 10:00", "2024-01-20 12:00", "2024-01-29 02:00"])})
    ref = pd.Timestamp("2024-01-30")
    assert sf.select_infection_alarm(al, ref) == -1.0          # nearest of the two in [-14, +7]
    assert np.isnan(sf.select_infection_alarm(al.iloc[:1], ref))  # -29 d is outside the window -> miss


# ---------------------------------------------------------------------------
# label helpers
# ---------------------------------------------------------------------------
def test_parse_timestamp_list():
    v = "[Timestamp('2023-08-29 00:00:00'), Timestamp('2022-11-14 00:00:00')]"
    assert sw.parse_timestamp_list(v) == [pd.Timestamp("2023-08-29"), pd.Timestamp("2022-11-14")]
    assert sw.parse_timestamp_list("[NaT]") == []


def test_mishra_reference_day_rules():
    T = pd.Timestamp
    # earliest symptom within 30 d before the earliest diagnosis (older, unrelated symptom dates ignored)
    assert sw.mishra_reference_day([T("2028-01-21"), T("2028-06-20")], [T("2028-06-21")]) == (T("2028-06-20"), "symptom_onset")
    assert sw.mishra_reference_day([T("2023-12-26"), T("2023-12-21")], [T("2023-12-26"), T("2023-12-21")])[0] == T("2023-12-21")
    assert sw.mishra_reference_day([], [T("2025-03-06")]) == (T("2025-03-06"), "diagnosis")
    assert sw.mishra_reference_day([T("2028-01-16")], []) == (T("2028-01-16"), "earliest_symptom_date")
    # symptom lists are unordered: onset never depends on list position
    assert sw.mishra_reference_day([T("2025-09-03"), T("2025-09-02")], []) == (T("2025-09-02"), "earliest_symptom_date")
    # AIFDJZB: only symptom is 43 d before diagnosis -> that symptom (paper Supp. Table 27 day 0)
    assert sw.mishra_reference_day([T("2023-11-07")], [T("2023-12-20")]) == (T("2023-11-07"), "symptom_onset_gt30d_before_diagnosis")


# ---------------------------------------------------------------------------
# processed outputs
# ---------------------------------------------------------------------------
TABLES_ = ["participants__stanford_covid", "participant_conditions__stanford_covid",
           "participant_wearable_features__stanford_covid", "participant_wearable_daily__stanford_covid"]
PREFIXES = ("stanford_covid_mishra2020:", "stanford_covid_alavi2022:", "stanford_longcovid_uwakwe2025:")


def _need(name):
    if not table_exists(name):
        pytest.skip(f"{name} not built; run the ingestion first")
    return read_table(name)


@pytest.mark.data
@pytest.mark.parametrize("name", TABLES_)
def test_tables_have_provenance_and_prefixed_ids(name):
    df = _need(name)
    assert len(df) > 0
    assert set(PROVENANCE_COLUMNS) <= set(df.columns)
    assert (df.data_layer == "person").all()
    assert (df.source_geographic_resolution == "none").all()
    assert df.participant_id.str.startswith(PREFIXES).all()


@pytest.mark.data
def test_participant_counts():
    p = _need("participants__stanford_covid")
    assert p.participant_id.is_unique
    n = p.groupby("dataset_id").size()
    assert n["stanford_covid_mishra2020"] == 118
    assert n["stanford_covid_alavi2022"] == 2123
    assert n["stanford_longcovid_uwakwe2025"] == 126
    m = p[p.dataset_id == "stanford_covid_mishra2020"].cohort_group.value_counts()
    assert m["covid19_positive"] == 32 and m["potential_healthy_no_illness_reported"] == 73
    a = p[(p.dataset_id == "stanford_covid_alavi2022") & (p.cohort_group == "covid19_positive")]
    assert len(a) == 84 and a.wearable_processed.all()


@pytest.mark.data
def test_long_covid_label_only_in_uwakwe():
    c = _need("participant_conditions__stanford_covid")
    lc = c[c.is_long_covid_label]
    assert set(lc.dataset_id) == {"stanford_longcovid_uwakwe2025"}
    assert len(lc) == 126 and int(lc.label_value.sum()) == 31
    acute = c[c.dataset_id != "stanford_longcovid_uwakwe2025"]
    assert not acute.condition_label.str.contains("long covid", case=False).any()
    covid = c[(c.dataset_id == "stanford_covid_mishra2020") & (c.condition_type == "infection_episode")]
    assert covid.native_id.nunique() == 32 and covid.reference_date.notna().all()


@pytest.mark.data
def test_daily_rows_unique_and_days_from_onset_per_participant():
    d = _need("participant_wearable_daily__stanford_covid")
    f = _need("participant_wearable_features__stanford_covid")
    assert not d.duplicated(["participant_id", "date"]).any()
    ref = f.set_index("participant_id").reference_date
    x = d[d.days_from_onset.notna()].copy()
    x["ref"] = x.participant_id.map(ref)
    assert ((x.date - x.ref).dt.days == x.days_from_onset).all()
    # the Long COVID subset has no released infection date -> never aligned
    assert d.loc[d.dataset_id == "stanford_longcovid_uwakwe2025", "days_from_onset"].isna().all()


@pytest.mark.data
def test_baseline_uses_only_pre_reference_days():
    d = _need("participant_wearable_daily__stanford_covid")
    f = _need("participant_wearable_features__stanford_covid")
    b = d[d.is_baseline_day.fillna(False)]
    assert (b.days_from_onset <= sf.BASELINE_WINDOW[1]).all() and (b.days_from_onset >= sf.BASELINE_WINDOW[0]).all()
    ok = f[f.baseline_ok.fillna(False)].set_index("participant_id")
    assert len(ok) > 0
    recomputed = b[b.rhr_mean.notna()].groupby("participant_id").rhr_mean.mean()
    np.testing.assert_allclose(recomputed.reindex(ok.index).values, ok.baseline_rhr_mean.values, rtol=1e-9)


@pytest.mark.data
def test_rhr_only_where_step_filter_applies():
    d = _need("participant_wearable_daily__stanford_covid")
    assert d.loc[~d.rhr_computable.astype(bool), "rhr_mean"].isna().all()
    assert d.loc[d.dataset_id == "stanford_covid_mishra2020", "rhr_computable"].astype(bool).all()


@pytest.mark.data
def test_recovery_censoring_consistent():
    f = _need("participant_wearable_features__stanford_covid")
    cens = f[f.rhr_recovery_event == 0]
    assert (cens.days_to_rhr_recovery == cens.last_followup_day).all()
    ev = f[f.rhr_recovery_event == 1]
    assert (ev.days_to_rhr_recovery >= 0).all()
    # recovery is only defined after an acute elevation
    assert f.loc[f.rhr_recovery_event.notna(), "rhr_acute_elevation"].fillna(False).all()
    assert f.trajectory_label[f.dataset_id != "stanford_longcovid_uwakwe2025"].str.contains("not a Long COVID label").all()
    # no published infection/illness date -> never labelled as a post-infection trajectory
    noref = f[(f.dataset_id != "stanford_longcovid_uwakwe2025") & f.reference_date.isna()]
    assert len(noref) > 0 and not noref.trajectory_label.str.startswith("post-infection").any()
    assert noref.days_to_rhr_recovery.isna().all() and noref.baseline_rhr_mean.isna().all()


@pytest.mark.data
def test_no_rhr_on_days_without_any_recorded_step():
    d = _need("participant_wearable_daily__stanford_covid")
    x = d[d.device_class != "preprocessed_resting_hr"]
    stepless = x.steps_recorded.isna() | (x.steps_recorded == 0)
    assert stepless.sum() > 0
    assert x.loc[stepless, ["rhr_mean", "rhr_night_mean", "steps", "rhr_z"]].isna().all().all()
    assert not x.loc[stepless, "step_stream_present"].astype(bool).any()
    # the authors' pre-masked Long COVID HR is untouched
    u = d[d.device_class == "preprocessed_resting_hr"]
    assert u.rhr_mean.notna().mean() > 0.9


@pytest.mark.data
def test_reproduction_benchmark_written():
    p = TABLES / "stanford_reproduction_benchmark.csv"
    if not p.exists():
        pytest.skip("benchmark not built")
    b = pd.read_csv(p).set_index("benchmark_id")
    port = b.loc["mishra2020_cusum_alarm_port"]
    # 220/223 published alarms reproduce exactly; the 3 misses all belong to one truncated public extract
    assert int(port.reproduced_value) >= 0.98 * int(port.paper_value)
    assert "only-in-ours 0" in port.details
    assert b.loc["uwakwe2025_cohort_size", "agreement"] == "agree"
    assert b.loc["uwakwe2025_table1", "agreement"] == "agree"
    assert b.loc["uwakwe2025_rfhrm_rfcm", "agreement"] == "not_reproducible"
    assert "mishra2020_cusum_presymptomatic_fraction" in b.index
