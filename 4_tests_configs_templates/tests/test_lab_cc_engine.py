"""Unit tests for measure_it.labs.lab_cc_engine (synthetic data only)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from measure_it.labs import lab_cc_engine as E
from measure_it.wearables.stanford_models import auroc_rows


def test_bootstrap_keeps_class_and_stratum_counts():
    y = np.array([1] * 7 + [0] * 13)
    st = np.array(["a", "b"] * 10)
    idx = E.stratified_boot_idx(y, n_boot=40, seed=1, strata=st)
    assert idx.shape == (40, 20)
    assert (y[idx].sum(axis=1) == 7).all()
    assert ((st[idx] == "a").sum(axis=1) == 10).all()


def test_cv_imputes_inside_folds_and_is_not_optimistic_on_noise():
    rng = np.random.default_rng(3)
    y = np.array([1] * 60 + [0] * 60)
    X = rng.normal(size=(120, 30))
    X[rng.random(X.shape) < 0.1] = np.nan
    P = E.cv_oof(X, y, C=0.1, n_repeats=3)
    assert np.isfinite(P).all()
    assert auroc_rows(y, P).mean() < 0.62
    Xs = np.nan_to_num(X) + 1.2 * y[:, None] * (np.arange(30) < 3)
    assert auroc_rows(y, E.cv_oof(Xs, y, C=0.1, n_repeats=3)).mean() > 0.8


def test_fitfree_direction_and_permutation():
    rng = np.random.default_rng(0)
    y = np.array([1] * 40 + [0] * 40)
    x = rng.normal(size=80) + 2.0 * y
    hi = E.fitfree_test(x, y, n_boot=200, n_perm=500)
    lo = E.fitfree_test(x, y, n_boot=200, n_perm=500, higher_is_case=False)
    assert hi["auroc"] > 0.8 and abs(hi["auroc"] + lo["auroc"] - 1) < 1e-9
    assert hi["perm_p"] < 0.01 and lo["perm_p"] > 0.5
    assert hi["auroc_ci_low"] <= hi["auroc"] <= hi["auroc_ci_high"]


def test_fitfree_strata_drop_single_class_strata():
    y = np.array([1, 1, 0, 0, 1, 1])
    x = np.array([2.0, 3.0, 1.0, 0.0, 9.0, 9.0])
    st = np.array(["p1", "p1", "p1", "p1", "p2", "p2"])      # p2 has no control -> ignored
    r = E.fitfree_test(x, y, n_boot=50, n_perm=50, strata=st)
    assert r["auroc"] == 1.0


def test_per_analyte_effects_bh_and_signs():
    rng = np.random.default_rng(1)
    df = pd.DataFrame({"y": [1] * 50 + [0] * 50, "up": np.r_[rng.normal(1, 1, 50), rng.normal(0, 1, 50)],
                       "flat": rng.normal(size=100)})
    e = E.per_analyte_effects(df, "y", ["up", "flat"], n_boot=100).set_index("analyte")
    assert e.loc["up", "hedges_g"] > 0.5 and e.loc["up", "q_value"] < 0.01
    assert e.loc["up", "q_value"] >= e.loc["up", "p_value"]


def test_delta_is_paired():
    rng = np.random.default_rng(2)
    y = np.array([1] * 50 + [0] * 50)
    X = rng.normal(size=(100, 2))
    b = E.stratified_boot_idx(y, n_boot=100)
    a1 = E.cv_block(X, y, 1.0, b)
    d = E.delta(a1, a1, y)
    assert d["delta_auroc"] == 0 and d["ci_low"] == 0 and d["ci_high"] == 0
