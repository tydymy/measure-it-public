"""Unit tests for measure_it.wearables.circadian (pure functions, hand-computable cases)."""
import math

import numpy as np
import pytest

from measure_it.wearables import circadian as C


def test_bin_mean_basic_and_min_valid():
    x = np.array([1, 3, np.nan, np.nan, 5, np.nan], dtype=float)
    out = C.bin_mean(x, 2)  # default min_valid = 1 for bin of 2
    assert out[0] == 2 and np.isnan(out[1]) and out[2] == 5
    out2 = C.bin_mean(x, 2, min_valid=2)
    assert out2[0] == 2 and np.isnan(out2[2])
    with pytest.raises(ValueError):
        C.bin_mean(np.ones(5), 2)


def test_is_perfectly_repeating_profile_is_one():
    day = np.array([0, 0, 1, 5, 9, 4, 2, 0], dtype=float)
    x = np.tile(day, 7)
    assert C.interdaily_stability(x, period=8) == pytest.approx(1.0)


def test_is_hand_computed_two_days():
    # two days, period 2: day1 = [0, 2], day2 = [2, 4]
    # grand mean 2; total var = (4+0+0+4)/4 = 2; profile = [1, 3] -> between var = (1+1)/2 = 1
    x = np.array([0, 2, 2, 4], dtype=float)
    assert C.interdaily_stability(x, 2) == pytest.approx(0.5)
    # textbook formula n*sum_h(xh-x)^2 / (p*sum_i(xi-x)^2) = 4*2/(2*8) = 0.5
    assert C.interdaily_stability(x, 2) == pytest.approx(4 * 2 / (2 * 8))


def test_is_white_noise_is_small_and_iv_near_two():
    rng = np.random.default_rng(0)
    x = rng.normal(size=24 * 200)
    assert C.interdaily_stability(x, 24) < 0.05
    assert C.intradaily_variability(x) == pytest.approx(2.0, abs=0.1)


def test_iv_hand_computed_and_smooth_signal():
    x = np.array([1, 2, 3, 4], dtype=float)
    # textbook: n*sum(diff^2) / ((n-1)*sum(dev^2)) = 4*3 / (3*5) = 0.8
    assert C.intradaily_variability(x) == pytest.approx(0.8)
    t = np.arange(24 * 7)
    smooth = np.sin(2 * np.pi * t / 24)
    assert C.intradaily_variability(smooth) < 0.1
    alternating = np.tile([1.0, -1.0], 50)
    assert C.intradaily_variability(alternating) == pytest.approx(4.0)


def test_iv_ignores_pairs_with_missing():
    x = np.array([1, 2, np.nan, 3, 4], dtype=float)
    # valid pairs: (1,2),(3,4) -> mean diff^2 = 1; valid mean 2.5, var = (2.25+0.25+0.25+2.25)/4 = 1.25
    assert C.intradaily_variability(x) == pytest.approx(1 / 1.25)


def test_is_missing_data_equals_complete_case_when_balanced():
    day = np.array([0, 1, 4, 1], dtype=float)
    x = np.tile(day, 3)
    y = x.copy()
    y[4:8] = np.nan  # drop the whole second day
    assert C.interdaily_stability(y, 4) == pytest.approx(1.0)
    assert np.isnan(C.interdaily_stability(np.full(8, np.nan), 4))


def test_m10_l5_hourly_profile_with_wrap():
    # hourly profile: active 08:00-18:00 (value 10), rest 00:00-05:00 and 23:00 (value 0), else 2
    prof = np.full(24, 2.0)
    prof[8:18] = 10.0
    prof[0:5] = 0.0
    r = C.m10_l5(prof, period=24, epoch_minutes=60, profile=True)
    assert r["m10"] == pytest.approx(10.0) and r["m10_onset_h"] == 8
    assert r["l5"] == pytest.approx(0.0) and r["l5_onset_h"] == 0
    assert r["ra"] == pytest.approx(1.0)
    # least-active block wrapping midnight: 22:00-03:00
    prof2 = np.full(24, 5.0)
    prof2[[22, 23, 0, 1, 2]] = 1.0
    r2 = C.m10_l5(prof2, period=24, epoch_minutes=60, profile=True)
    assert r2["l5_onset_h"] == 22 and r2["l5"] == pytest.approx(1.0)


def test_m10_l5_from_minute_series_and_ra():
    minute_day = np.full(1440, 1.0)
    minute_day[9 * 60:19 * 60] = 21.0  # M10 = 21 starting 09:00
    minute_day[1 * 60:6 * 60] = 0.0    # L5 = 0 starting 01:00
    x = np.tile(minute_day, 3)
    r = C.m10_l5(x, period=1440, epoch_minutes=1)
    assert r["m10"] == pytest.approx(21.0) and r["m10_onset_h"] == pytest.approx(9.0)
    assert r["l5"] == pytest.approx(0.0) and r["l5_onset_h"] == pytest.approx(1.0)
    assert C.relative_amplitude(3.0, 1.0) == pytest.approx(0.5)
    assert np.isnan(C.relative_amplitude(0.0, 0.0))


def test_run_lengths_and_transition_probabilities():
    state = np.array([1, 1, 0, 0, 0, 1, 0, 1, 1, 1], dtype=bool)
    act, sed = C.run_lengths(state)
    assert list(act) == [2, 1, 3] and list(sed) == [3, 1]
    tp = C.transition_probabilities(state)
    assert tp["mean_active_bout"] == pytest.approx(2.0) and tp["astp"] == pytest.approx(0.5)
    assert tp["mean_sedentary_bout"] == pytest.approx(2.0) and tp["satp"] == pytest.approx(0.5)
    # an invalid epoch splits a bout
    valid = np.ones(10, dtype=bool)
    valid[8] = False
    act2, sed2 = C.run_lengths(state, valid)
    assert list(act2) == [2, 1, 1, 1] and list(sed2) == [3, 1]


def test_transition_probabilities_all_sedentary():
    tp = C.transition_probabilities(np.zeros(20, dtype=bool))
    assert np.isnan(tp["astp"]) and tp["satp"] == pytest.approx(1 / 20)


def test_main_rest_window_bridging():
    r = np.zeros(30, dtype=bool)
    r[5:10] = True
    r[12:20] = True   # gap of 2 between 10 and 12
    r[25:27] = True
    assert C.main_rest_window(r) == (12, 20)
    assert C.main_rest_window(r, gap_tolerance=2) == (5, 20)
    assert C.main_rest_window(np.zeros(5, dtype=bool)) is None
    assert C.main_rest_window(r, min_length=100) is None


def test_circular_mean_and_sd_hours():
    # clock hours must lie in [0, 24): a mean at midnight is 0.0, never 24.0 (regression: 24.0 was returned)
    for hours in ([23.0, 1.0], [23.5, 0.5], [22.0, 2.0], [0.0, 0.0]):
        m = C.circular_mean_hours(hours)
        assert 0.0 <= m < 24.0 and m == pytest.approx(0.0, abs=1e-9), (hours, m)
    assert C.circular_mean_hours([23.0, 23.5]) == pytest.approx(23.25)
    assert C.circular_mean_hours([2.0, 4.0]) == pytest.approx(3.0)
    assert C.circular_sd_hours([3.0, 3.0, 3.0]) == pytest.approx(0.0, abs=1e-5)
    assert C.circular_sd_hours([2.0, 4.0]) > 0.9
    assert math.isnan(C.circular_mean_hours([np.nan]))


def test_rolling_mean_centred():
    x = np.array([0, 0, 3, 0, 0], dtype=float)
    out = C.rolling_mean(x, 3, min_valid=1)
    assert out[2] == pytest.approx(1.0) and out[1] == pytest.approx(1.0) and out[0] == pytest.approx(0.0)
