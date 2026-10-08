"""Tests for the fibromyalgia infrared-thermography source (measure_it.ingestion.fm_thermography) and its
pre-specified analysis (measure_it.wearables.fm_thermography_analysis).

Unit tests use synthetic data. Data tests (marked `data`) check the processed partitions, the result tables and the
guardrails (namespaced ids, no geography, provenance, the locked decision rule applied as written).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from measure_it.config import DOCS, PROCESSED, RESULTS, TABLES
from measure_it.ingestion import fm_thermography as ING
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.wearables import fm_thermography_analysis as A


# --------------------------------------------------------------------------- unit
def test_feature_and_id_naming():
    assert ING.feature_name("Upper_back_ave") == "skin_temp_upper_back_ave_c"
    assert ING.participant_id("FM001") == "fm_thermography:FM001"
    assert len(ING.FEATURE_COLS) == 18 and len(ING.AVE_FEATURES) == 6
    assert all(f.endswith("_ave_c") for f in ING.AVE_FEATURES)


def test_derive_label_checks_prefix_against_group():
    ok = pd.DataFrame({"Participant": ["FM001", "CG001"], "Group": [1, 2]})
    assert ING.derive_label(ok).tolist() == ["fibromyalgia", "control"]
    with pytest.raises(ValueError):
        ING.derive_label(pd.DataFrame({"Participant": ["FM001"], "Group": [2]}))
    with pytest.raises(ValueError):
        ING.derive_label(pd.DataFrame({"Participant": ["FM001"], "Group": [3]}))


def test_sens_at_spec_perfect_and_random():
    y = np.array([1] * 10 + [0] * 10)
    perfect = np.r_[np.linspace(0.6, 1, 10), np.linspace(0, 0.4, 10)]
    assert A.sens_at_spec(y, perfect)[0] == 1.0
    tied = np.full(20, 0.5)
    assert A.sens_at_spec(y, tied)[0] == 0.0       # ties never count as a positive call above the threshold


def test_stratified_bootstrap_keeps_class_counts():
    y = np.array([1] * 7 + [0] * 13)
    idx = A.stratified_boot_idx(y, n_boot=50, seed=1)
    assert idx.shape == (50, 20)
    assert (y[idx].sum(axis=1) == 7).all()


def test_folds_depend_only_on_labels_and_cv_is_paired():
    rng = np.random.default_rng(0)
    y = np.array([1] * 30 + [0] * 30)
    f1, f2 = A.folds(y, 3), A.folds(y, 3)
    assert all((a[1] == b[1]).all() for a, b in zip(f1, f2))
    X = rng.normal(size=(60, 2))
    P = A.cv_oof(X, y, n_repeats=2)
    assert P.shape == (2, 60) and np.isfinite(P).all()


def test_cv_is_not_optimistic_on_noise_and_finds_real_signal():
    """Out-of-fold AUROC must not reward noise (pooled OOF AUROC on noise is, if anything, below 0.5) and must
    recover a real shift."""
    from measure_it.wearables.stanford_models import auroc_rows
    rng = np.random.default_rng(5)
    y = np.array([1] * 80 + [0] * 80)
    noise = rng.normal(size=(160, 6))
    assert auroc_rows(y, A.cv_oof(noise, y, n_repeats=3)).mean() < 0.6
    signal = noise + 1.0 * y[:, None]
    assert auroc_rows(y, A.cv_oof(signal, y, n_repeats=3)).mean() > 0.85


# --------------------------------------------------------------------------- data


def _t(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p)


@pytest.mark.data
def test_participants_partition():
    p = _t("participants__fm_thermography")
    assert len(p) == 178 and p.participant_id.is_unique
    assert p.participant_id.str.match(r"^fm_thermography:(FM|CG)\d{3}$").all()
    assert int(p.fibromyalgia.sum()) == 86 and int((1 - p.fibromyalgia).sum()) == 92
    assert (p.loc[p.fibromyalgia == 1, "native_id"].str.startswith("FM")).all()
    assert int(p.complete_ave_temperatures.sum()) == 177
    assert set(PROVENANCE_COLUMNS) <= set(p.columns)
    assert (p.data_layer == "person").all() and (p.source_geographic_resolution == "none").all()
    assert not any(c in p.columns for c in ("state", "county", "zip", "fips", "zcta", "latitude", "longitude"))


@pytest.mark.data
def test_device_features_partition():
    f = _t("participant_device_features__fm_thermography")
    assert len(f) == 178
    assert set(ING.FEATURE_COLS) <= set(f.columns)
    vals = f[ING.FEATURE_COLS].to_numpy(float)
    assert np.nanmin(vals) > 20 and np.nanmax(vals) < 40          # degC skin temperatures
    assert int(f[ING.FEATURE_COLS].isna().sum().sum()) == 6         # one FM woman lacks knee + elbow (3 stats each)
    assert (f.evidence_type == "person_device_measurement").all()


@pytest.mark.data
def test_primary_decision_applies_locked_rule():
    d = pd.read_csv(TABLES / "fm_thermography_primary_decision.csv").iloc[0]
    assert d.n_cases == 85 and d.n_controls == 92 and d.n_perm == 1000
    met = (d.delta_ci_low > 0) and (d.permutation_p_device_only < 0.05)
    assert bool(d.criterion_delta_ci_low_gt_0) == (d.delta_ci_low > 0)
    assert bool(d.criterion_perm_p_lt_0_05) == (d.permutation_p_device_only < 0.05)
    assert d.decision.startswith("thermography adds information") == met
    assert d.delta_ci_low <= d.delta_auroc <= d.delta_ci_high


@pytest.mark.data
def test_model_table_and_null():
    perf = pd.read_csv(TABLES / "fm_thermography_model_performance.csv")
    assert set(perf.feature_set) == {"demo", "device", "demo_device", "demo_device18"}
    assert ((perf.auroc_ci_low <= perf.auroc) & (perf.auroc <= perf.auroc_ci_high)).all()
    null = pd.read_csv(TABLES / "fm_thermography_permutation_null.csv")
    assert len(null) == 1000 and abs(null.device_only_auroc_null.mean() - 0.5) < 0.03


@pytest.mark.data
def test_reproduces_paper_table2():
    c = pd.read_csv(TABLES / "fm_thermography_paper_concordance.csv")
    assert len(c) == 18
    assert c.means_concordant_0_06C.all()
    assert c.significance_concordant.all()


@pytest.mark.data
def test_signatures_partition():
    s = _t("phenotype_signatures__fm_thermography")
    assert s.object_id.is_unique
    assert s.object_id.str.match(r"^signature:fm_thermography\|fibromyalgia\|\S+$").all()
    assert (s.condition_id == "fibromyalgia").all() and (~s.is_proxy).all()
    prim = s[s.effect_measure == "delta_auroc"]
    assert len(prim) == 1
    d = pd.read_csv(TABLES / "fm_thermography_primary_decision.csv").iloc[0]
    expected = ("supported_single_dataset" if d.criterion_delta_ci_low_gt_0 and d.criterion_perm_p_lt_0_05
                else "null_single_dataset")
    assert prim.evidence_level.iloc[0] == expected
    assert set(PROVENANCE_COLUMNS) <= set(s.columns)


@pytest.mark.data
def test_plan_and_report_exist_and_state_the_rule():
    plan = (DOCS / "ANALYSIS_PLAN_FM_THERMOGRAPHY.md").read_text()
    assert "before any case-control comparison was computed" in plan
    assert "CI lower bound" in plan and "1,000 label shuffles" in plan
    rep = (RESULTS / "FM_THERMOGRAPHY_RESULTS.md").read_text()
    assert "pre-specified primary test is" in rep
    meta = json.loads((TABLES / "fm_thermography_run_metadata.json").read_text())
    assert meta["n_perm"] == 1000 and meta["n_boot"] == 2000 and meta["n_repeats"] == 20


@pytest.mark.data   # the normaliser reads data/processed/condition_registry.parquet
def test_condition_id_matches_ontology_normalizer():
    from measure_it.ontology.normalize import normalize_condition
    m = normalize_condition("fibromyalgia")
    assert m["status"] == "matched"
    assert m["matches"][0]["canonical_condition_id"] == ING.CONDITION_ID
