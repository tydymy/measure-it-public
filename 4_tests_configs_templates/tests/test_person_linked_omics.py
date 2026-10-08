"""Tests for measure_it.omics.person_linked (mapMECFS participant-linked omics, person-level modelling).

Unit tests check the statistics helpers and that every data-dependent preprocessing step uses the training fold
only. Tests marked `data` check the outputs of `uv run python -m measure_it.omics.person_linked`.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from sklearn.metrics import roc_auc_score

from measure_it import config
from measure_it.omics import person_linked as pl
from measure_it.provenance import check_provenance

T = config.TABLES


def _layer(n=30, p=200, rna=False, seed=0):
    rng = np.random.default_rng(seed)
    y = np.r_[np.zeros(n // 2, int), np.ones(n - n // 2, int)]
    X = rng.normal(size=(n, p))
    cov = np.column_stack([rng.integers(0, 2, n).astype(float), rng.normal(45, 12, n)])
    expressed = rng.random((n, p)) > 0.3 if rna else None
    return pl.Layer("synthetic", np.array([f"p{i}" for i in range(n)]), y, X, np.array([f"f{j}" for j in range(p)]),
                    np.array([f"g{j}" for j in range(p)]), cov, ("female", "age"), expressed)


# --------------------------------------------------------------------------- unit tests
def test_welch_t_matches_scipy():
    L = _layer()
    ours = pl.welch_t(L.X, L.y)
    ref = stats.ttest_ind(L.X[L.y == 1], L.X[L.y == 0], axis=0, equal_var=False).statistic
    np.testing.assert_allclose(ours, ref, rtol=1e-10)


def test_welch_t_zero_variance_is_zero():
    X = np.ones((6, 2))
    assert np.all(pl.welch_t(X, np.array([0, 0, 0, 1, 1, 1])) == 0)


def test_auroc_matches_sklearn_with_ties():
    rng = np.random.default_rng(1)
    y = rng.integers(0, 2, 40)
    s = rng.integers(0, 5, 40).astype(float)
    assert pl.auroc(y, s) == pytest.approx(roc_auc_score(y, s))
    S = np.vstack([s, -s])
    np.testing.assert_allclose(pl.auroc_rows(y, S), [roc_auc_score(y, s), roc_auc_score(y, -s)])


def test_bh_matches_scipy_and_keeps_nan():
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.5])
    q = pl.bh(p)
    assert np.isnan(q[2])
    np.testing.assert_allclose(q[[0, 1, 3, 4]], stats.false_discovery_control(p[[0, 1, 3, 4]]))


def test_perm_p_formula():
    assert pl._perm_p(0.7, np.array([0.5, 0.6, 0.7, 0.8])) == pytest.approx(3 / 5)
    assert pl._perm_p(0.9, np.zeros(999)) == pytest.approx(1 / 1000)


def test_rna_log_transform_is_per_sample_cpm():
    wide = pd.DataFrame([[10, 90], [1, 999]], columns=["a", "b"])
    X, expressed = pl.log_transform("pbmc_rnaseq", wide)
    np.testing.assert_allclose(X[0], np.log2(np.array([1e5, 9e5]) + 1))
    assert expressed.all()
    with pytest.raises(ValueError):
        pl.log_transform("csf_metabolomics", pd.DataFrame([[0.0, 1.0]]))


@pytest.mark.parametrize("spec_id", ["l2_top50", "l2_top50_resid", "enet_top50"])
def test_selection_and_scaling_use_training_fold_only(spec_id):
    """Changing test-fold values must not change which features are selected or the test rows' neighbours' scores."""
    L = _layer(rna=True)
    tr, te = np.arange(0, 24), np.arange(24, 30)
    s1, i1 = pl.fit_predict(L, L.y, tr, te, pl.MODELS[spec_id])
    L2 = pl.Layer(L.name, L.ids, L.y, L.X.copy(), L.features, L.feature_names, L.cov.copy(), L.cov_names,
                  L.expressed.copy())
    L2.X[te] = L2.X[te] * 50 + 7              # wreck the test rows
    L2.expressed[te] = ~L2.expressed[te]
    L2.cov[te] = L2.cov[te] + 100
    _, i2 = pl.fit_predict(L2, L2.y, tr, te, pl.MODELS[spec_id])
    np.testing.assert_array_equal(i1["selected"], i2["selected"])
    # and the test labels are never read: flipping them leaves the scores unchanged
    y_flip = L.y.copy()
    y_flip[te] = 1 - y_flip[te]
    s3, _ = pl.fit_predict(L, y_flip, tr, te, pl.MODELS[spec_id])
    np.testing.assert_allclose(s1, s3)


def test_no_signal_cv_auroc_is_not_inflated():
    """Pure noise with 10,000 features and in-fold selection must not look predictive."""
    L = _layer(n=30, p=10000, seed=3)
    S, _ = pl.cv_scores([L], L.y, pl.MODELS["l2_top50"], seed=config.SEED)
    assert pl.auroc_rows(L.y, S).mean() < 0.7


def test_strong_signal_is_found():
    L = _layer(n=30, p=500, seed=4)
    L.X[L.y == 1, :20] += 2.0
    S, _ = pl.cv_scores([L], L.y, pl.MODELS["l2_top50"], seed=config.SEED)
    assert pl.auroc_rows(L.y, S).mean() > 0.9


def test_splits_are_stratified_and_complete():
    y = np.r_[np.zeros(11, int), np.ones(10, int)]
    for folds in pl.make_splits(y, config.SEED):
        te_all = np.concatenate([te for _, te in folds])
        assert sorted(te_all) == list(range(len(y)))
        for _, te in folds:
            assert 1 <= y[te].sum() <= 3


# --------------------------------------------------------------------------- data tests
needs_outputs = pytest.mark.skipif(not (T / "linked_omics_model_performance.csv").exists(),
                                   reason="run `uv run python -m measure_it.omics.person_linked` first")


@pytest.mark.data
@needs_outputs
def test_model_table_covers_every_layer_and_model_with_full_permutations():
    perf = pd.read_csv(T / "linked_omics_model_performance.csv")
    assert set(perf["layer"]) == set(pl.LAYERS)
    for lay in pl.LAYERS:
        assert set(perf.loc[perf.layer == lay, "model"]) == set(pl.MODELS)
    prim = perf[perf.model.isin(pl.PERMUTED_MODELS)]
    assert (prim["n_perm"] == pl.N_PERM).all()
    assert prim["perm_p"].between(1 / (pl.N_PERM + 1) - 1e-12, 1).all()
    assert (perf["ci_low"] <= perf["auroc"] + 1e-9).all() and (perf["auroc"] <= perf["ci_high"] + 1e-9).all()


@pytest.mark.data
@needs_outputs
def test_layer_sizes_match_linkage_audit():
    perf = pd.read_csv(T / "linked_omics_model_performance.csv")
    n = perf[perf.model == "l2_top50"].set_index("layer")[["n_cases", "n_controls"]]
    expected = {"csf_somalogic": (21, 21), "serum_somalogic": (21, 21), "pbmc_rnaseq": (12, 15),
                "muscle_rnaseq": (13, 12), "csf_metabolomics": (10, 11)}
    for lay, (c, h) in expected.items():
        assert tuple(n.loc[lay]) == (c, h)


@pytest.mark.data
@needs_outputs
def test_top_features_only_for_passing_layers():
    perf = pd.read_csv(T / "linked_omics_model_performance.csv")
    passing = set(perf[(perf.model == "l2_top50") & (perf.perm_p < pl.PASS_ALPHA)]["layer"])
    top = pd.read_csv(T / "linked_omics_top_features.csv")
    assert set(top["layer"].dropna()) <= passing


@pytest.mark.data
@needs_outputs
def test_fusion_n_matches_audit_overlaps():
    f = pd.read_csv(T / "linked_omics_fusion.csv").set_index("fusion")
    assert f.loc["F1", "n"] == 42 and f.loc["F2", "n"] == 27 and f.loc["F3", "n"] == 21
    assert f.loc["F4", "n"] == 17 and f.loc["F5", "n"] == 9
    assert not bool(f.loc["F5", "analysed"])


@pytest.mark.data
@needs_outputs
def test_run_metadata_records_plan_and_seed():
    meta = json.loads((T / "linked_omics_run_metadata.json").read_text())
    assert meta["seed"] == config.SEED and meta["n_perm"] == pl.N_PERM and meta["plan"] == pl.PLAN


@pytest.mark.data
def test_signature_partition_schema_and_labels():
    path = config.PROCESSED / f"{pl.SIGNATURE_TABLE}.parquet"
    if not path.exists():
        pytest.skip("signature partition not built")
    sig = pd.read_parquet(path)
    check_provenance(sig, pl.SIGNATURE_TABLE)
    canon = config.PROCESSED / "phenotype_signatures__nhanes.parquet"
    if canon.exists():
        other = set(pd.read_parquet(canon).columns)
        assert other <= set(sig.columns)
    assert (sig["condition_id"] == "me_cfs").all()
    assert (sig["label_basis"] == pl.LABEL_BASIS).all()
    assert sig["object_id"].is_unique
    assert sig["object_id"].str.startswith(f"signature:{pl.DATASET_ID}|pi_mecfs_vs_hv__").all()
    assert (sig["data_layer"] == "person").all()
    model = sig[sig["feature"].str.startswith("model_auroc:")]
    assert set(model["phenotype_id"]) == {f"pi_mecfs_vs_hv__{lay}" for lay in pl.LAYERS}
    assert model["q_value_bh"].isna().all()          # model rows are never counted as FDR feature hits
    assert "wearable" in sig["caveats"].iloc[0]


@pytest.mark.data
def test_signature_tool_reports_the_partition():
    if not (config.PROCESSED / f"{pl.SIGNATURE_TABLE}.parquet").exists():
        pytest.skip("signature partition not built")
    from measure_it.wearables.signatures import get_patient_phenotype_signature
    r = get_patient_phenotype_signature("ME/CFS")
    assert r["status"] == "found"
    ours = [d for d in r["datasets"] if d["dataset_id"] == pl.DATASET_ID]
    assert {d["phenotype_id"] for d in ours} == {f"pi_mecfs_vs_hv__{lay}" for lay in pl.LAYERS}
    for d in ours:
        assert d["model_level_results"]["rows"]
