"""Tests for the labs-vs-diagnosis analyses (measure_it.labs.*) under docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md.

Unit tests use synthetic data; tests marked `data` check the outputs of
`uv run python -m measure_it.labs.labs_vs_diagnosis` and `uv run python -m measure_it.labs.appelman_metabolomics`.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from measure_it.config import TABLES
from measure_it.labs import appelman_metabolomics as ap
from measure_it.labs import labs_vs_diagnosis as lv
from measure_it.labs import plan
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table, table_exists


# --------------------------------------------------------------------------------------------- plan lock
def test_plan_is_locked_and_unchanged():
    assert plan.PLAN_PATH.exists()
    assert hashlib.sha256(plan.PLAN_PATH.read_text().encode()).hexdigest() == plan.plan_sha256()
    assert plan.check_plan_locked() == plan.plan_sha256()
    lock = json.loads(plan.LOCK_PATH.read_text())
    assert lock["sha256"] == plan.plan_sha256()


def test_plan_lists_every_label():
    for lab in lv.LABELS:
        if lab.group != "pairwise":
            assert f"`{lab.label_id}`" in plan.PLAN_MD, lab.label_id


# --------------------------------------------------------------------------------------------- estimators
def test_ols_hc3_matches_statsmodels():
    import statsmodels.api as sm
    rng = np.random.default_rng(0)
    n = 300
    X = np.column_stack([np.ones(n), rng.integers(0, 2, n), rng.normal(size=n), rng.integers(0, 2, n)])
    z = X @ np.array([0.1, 0.4, -0.2, 0.3]) + rng.standard_t(4, n)
    b, se, p = lv.ols_hc3(z, X)
    fit = sm.OLS(z, X).fit(cov_type="HC3")
    assert b == pytest.approx(fit.params[1])
    assert se == pytest.approx(fit.bse[1])
    assert p == pytest.approx(fit.pvalues[1], rel=1e-6)


def test_adjusted_effect_matches_project_signature_estimator():
    from measure_it.wearables.signatures import adjusted_std_difference
    rng = np.random.default_rng(1)
    n = 500
    age, female = rng.uniform(20, 80, n), rng.integers(0, 2, n).astype(float)
    case = (rng.uniform(size=n) < 0.2).astype(float)
    x = 0.5 * case + 0.02 * age + rng.normal(size=n)
    a, b = lv.adjusted_effect(x, case, age, female), adjusted_std_difference(x, case, age, female)
    assert a["effect_size"] == pytest.approx(b["effect_size"])
    assert a["se"] == pytest.approx(b["se"])


def test_adjusted_effect_needs_min_cases():
    rng = np.random.default_rng(2)
    case = np.zeros(200)
    case[:5] = 1
    r = lv.adjusted_effect(rng.normal(size=200), case, rng.uniform(20, 80, 200), rng.integers(0, 2, 200) * 1.0)
    assert np.isnan(r["effect_size"]) and r["n_cases"] == 5


def test_hedges_g_known_value():
    x1, x0 = np.array([2.0, 3.0, 4.0]), np.array([1.0, 2.0, 3.0])
    # pooled SD 1, d = 1, J = 1 - 3/(4*6-9) = 0.8
    assert ap.hedges_g(x1, x0) == pytest.approx(0.8)


def test_skew_transform_is_label_blind_and_logs_skewed():
    s = pd.Series(np.exp(np.random.default_rng(3).normal(size=500)))
    t, name = lv._skew_transform(s)
    assert name == "log" and np.allclose(t, np.log(s))
    s2 = pd.Series(np.random.default_rng(4).normal(size=500))
    assert lv._skew_transform(s2)[1] == "none"


def test_stratified_bootstrap_keeps_group_sizes():
    y = np.r_[np.ones(7), np.zeros(5)].astype(int)
    W = ap.strat_boot_weights(y, 50, 0)
    assert np.allclose(W[:, y == 1].sum(1), 7) and np.allclose(W[:, y == 0].sum(1), 5)


def test_top10_selection_is_inside_the_fold():
    """Predictions for a test fold must not change when the test fold's labels are permuted."""
    rng = np.random.default_rng(5)
    n, p = 40, 30
    X = rng.normal(size=(n, p))
    y = np.r_[np.ones(20), np.zeros(20)].astype(int)
    tr, te = np.arange(0, 32), np.arange(32, 40)
    a = ap.fit_predict(X, None, y, tr, te, "top10", 0)
    y2 = y.copy()
    y2[te] = 1 - y2[te]
    b = ap.fit_predict(X, None, y2, tr, te, "top10", 0)
    assert np.allclose(a, b)


def test_label_definitions():
    ids = [x.label_id for x in lv.LABELS]
    assert len(ids) == len(set(ids))
    dia = lv.LABEL_BY_ID["diabetes_no_glycaemic_labs"]
    assert {"LBXGH", "LBXSGL", "LBXGLU", "LBXGLT"} <= set(dia.drop_analytes)
    assert not lv.LABEL_BY_ID["diabetes_with_glycaemic_labs"].drop_analytes
    assert lv.LABEL_BY_ID["mecfs_like_proxy"].is_proxy
    assert all(x.cycles == "H" for x in lv.LABELS if x.label_id.startswith("rx_"))


# --------------------------------------------------------------------------------------------- data tests
@pytest.mark.data
def test_population_and_pairwise_labels_are_disjoint():
    frame, flow = lv.load_population()
    assert flow[-1]["n"] == len(frame) and (frame["age_years"] >= 20).all()
    assert frame["participant_id"].is_unique and frame["participant_id"].str.startswith("nhanes:").all()
    for pid, a, b, _ in lv.PAIRS:
        y = frame[f"y_{pid}"]
        assert ((y == 1) <= (frame[f"y_{a}"] == 1)).all()
    rx = frame[frame["y_rx_fibromyalgia"].notna()]
    assert (rx["cycle_code"] == "H").all()


@pytest.mark.data
def test_analyte_dictionary():
    dic = pd.read_csv(TABLES / "labs_dx_nhanes_analytes.csv")
    assert dic["analyte"].is_unique
    prim = dic[(dic["tier"] == "primary") & dic["included"]]
    assert len(prim) >= 50
    assert "LBXSCH" not in set(prim["analyte"])
    # NHANES reuses LBXHCT (hematocrit in CBC, hydroxycotinine in COT_H): the latter must be namespaced
    assert (dic.loc[dic["analyte"] == "LBXHCT", "lab_name"] == "Hematocrit").all()
    assert dic["analyte"].str.startswith("LBXHCT__").any()
    assert not dic.loc[dic["included"], "layer"].eq("urine_flow_rate").any()


@pytest.mark.data
def test_effects_respect_label_defining_exclusions():
    eff = pd.read_csv(TABLES / "labs_dx_nhanes_effects.csv")
    d = eff[eff["label_id"] == "diabetes_no_glycaemic_labs"]
    assert not d["analyte"].isin(lv.GLYCAEMIC).any()
    assert eff.loc[eff["label_id"] == "diabetes_with_glycaemic_labs", "analyte"].eq("LBXGH").any()
    ok = eff["q_value_bh"].notna()
    assert ((eff.loc[ok, "q_value_bh"] >= eff.loc[ok, "p_value"] - 1e-12)).all()


@pytest.mark.data
def test_nhanes_signature_partition():
    assert table_exists("phenotype_signatures__nhanes_labs")
    s = read_table("phenotype_signatures__nhanes_labs")
    assert set(PROVENANCE_COLUMNS) <= set(s.columns)
    assert set(lv.SIG_COLUMNS) <= set(s.columns)
    assert s["object_id"].is_unique and s["object_id"].str.startswith("signature:nhanes_labs|").all()
    assert (s["data_layer"] == "person").all()
    assert s.loc[s["phenotype_id"] == "mecfs_like_proxy", "is_proxy"].all()


@pytest.mark.data
def test_models_and_positive_control_gate_recorded():
    meta = json.loads((TABLES / "labs_dx_nhanes_run_metadata.json").read_text())
    assert meta["plan_sha256"] == plan.plan_sha256()
    gate = meta["positive_control_gate"]
    assert set(gate["required"]) == {"diabetes_no_glycaemic_labs", "weak_failing_kidneys", "gout"}
    perf = pd.read_csv(TABLES / "labs_dx_nhanes_model_performance.csv")
    assert perf["auroc"].between(0, 1).all()
    triv = perf[(perf.label_id == "diabetes_with_glycaemic_labs") & (perf.feature_set == "labs")]
    assert triv["auroc"].iloc[0] > 0.85  # HbA1c alone should nearly separate self-reported diabetes
    perm = pd.read_csv(TABLES / "labs_dx_nhanes_permutation_null.csv")
    assert perm["null_mean"].between(0.4, 0.6).all()


@pytest.mark.data
def test_appelman_outputs():
    m = pd.read_csv(TABLES / "labs_dx_appelman_metabolites.csv")
    b = m[(m["tissue"] == "blood") & (m["analysis"] == "baseline")]
    assert len(b) == 83 and (b["n_lc"] == 25).all() and (b["n_healthy"] == 21).all()
    mu = m[(m["tissue"] == "muscle") & (m["analysis"] == "baseline")]
    assert len(mu) == 116 and (mu["n_healthy"] == 19).all()
    mods = pd.read_csv(TABLES / "labs_dx_appelman_models.csv")
    prim = mods[mods["role"] == "primary"]
    assert set(prim["tissue"]) == {"blood", "muscle"} and prim["perm_p"].notna().all()
    s = read_table("phenotype_signatures__appelman_metabolomics")
    assert s["object_id"].is_unique and set(PROVENANCE_COLUMNS) <= set(s.columns)
    d = pd.read_csv(TABLES / "labs_dx_appelman_paper_direction.csv")
    assert len(d) == len(ap.CLAIMS)
