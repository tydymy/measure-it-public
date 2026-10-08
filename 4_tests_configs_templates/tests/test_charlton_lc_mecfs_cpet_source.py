"""Tests for the MUSCLE-ME ingestion (measure_it.ingestion.muscle_me_charlton) and the locked steps analysis
(measure_it.wearables.muscle_me_steps_analysis).

Unit tests use synthetic data; tests marked `data` check the outputs written by
`uv run python -m measure_it.ingestion.muscle_me_charlton` and `... wearables.muscle_me_steps_analysis`.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from measure_it.config import TABLES
from measure_it.ingestion import muscle_me_charlton as ing
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table, table_exists
from measure_it.wearables import muscle_me_steps_analysis as an

SID = ing.SOURCE_ID


# --------------------------------------------------------------------------------------------- pure functions
def test_auroc_matches_sklearn_with_ties():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 70)
    s = np.round(rng.normal(size=70), 1)
    assert an.auroc(y, s) == pytest.approx(roc_auc_score(y, s))


def test_auroc_pairs_matches_scalar():
    rng = np.random.default_rng(1)
    s1, s0 = rng.normal(size=(5, 12)), rng.normal(size=(5, 9))
    got = an.auroc_pairs(s1, s0)
    want = [an.auroc(np.r_[np.ones(12), np.zeros(9)], np.r_[s1[b], s0[b]]) for b in range(5)]
    assert np.allclose(got, want)


def test_fixed_direction_lower_is_case():
    y = np.r_[np.ones(10), np.zeros(10)].astype(int)
    x = np.r_[np.arange(10), np.arange(10) + 100.0]  # cases all lower
    r = an.fixed_direction_auroc(y, x, n_boot=200)
    assert r["auroc"] == 1.0 and r["decision_rule_result"].startswith("wearable captures")


def test_decision_rule():
    assert an.decide(0.72, 0.9).startswith("wearable captures")
    assert an.decide(0.45, 0.8) == "null"
    assert an.decide(0.6, 0.8) == "weak"
    assert an.decide(0.2, 0.45).startswith("reverse")


def test_strat_boot_keeps_group_sizes():
    y = np.r_[np.ones(7), np.zeros(4)].astype(int)
    B = an.strat_boot(y, 50)
    assert B.shape == (50, 11)
    assert (y[B[:, :7]] == 1).all() and (y[B[:, 7:]] == 0).all()


def test_hedges_g_sign():
    assert an.hedges_g(np.array([2.0, 3, 4]), np.array([0.0, 1, 2])) > 0


def test_participant_id_namespaced():
    assert ing.participant_id("P25") == f"{SID}:P25"


def test_select_muscle_me_drops_bed_rest_and_rejects_duplicates():
    cols = {c: [1.0, 2.0, 3.0] for c in ing.VARIABLES}
    raw = pd.DataFrame({"Subject": ["P1", "P25", "P80"], "Session": ["BDC", "LC", "CON"],
                        "Group": ["AGBRESA", "MUSCLE-ME", "MUSCLE-ME"], "Sex": ["Male"] * 3, "Sx_duration": np.nan, **cols})
    m = ing.select_muscle_me(raw)
    assert list(m.Subject) == ["P25", "P80"]
    dup = pd.concat([raw, raw.iloc[[1]]])
    with pytest.raises(ValueError):
        ing.select_muscle_me(dup)


def test_plan_is_locked_and_unchanged():
    assert an.PLAN_PATH.exists()
    assert hashlib.sha256(an.PLAN_PATH.read_text().encode()).hexdigest() == an.plan_sha256()


# --------------------------------------------------------------------------------------------- processed outputs
needs = pytest.mark.skipif(not table_exists(f"phenotype_signatures__{SID}"), reason="run the MUSCLE-ME modules first")


@pytest.mark.data
@needs
def test_participants_counts_and_labels():
    p = read_table(f"participants__{SID}")
    assert len(p) == 81 and p.participant_id.is_unique
    assert p.participant_id.str.startswith(f"{SID}:").all()
    assert p.session_code.value_counts().to_dict() == {"CON": 30, "ME": 26, "LC": 25}
    assert set(p.loc[p.session_code == "LC", "condition_id"]) == {"long_covid"}
    assert set(p.loc[p.session_code == "ME", "condition_id"]) == {"me_cfs"}
    assert p.loc[p.session_code == "CON", "condition_id"].isna().all()
    assert int(p.has_steps.sum()) == 73 and int(p[p.is_patient == 0].has_steps.sum()) == 22
    assert set(PROVENANCE_COLUMNS) <= set(p.columns)
    assert (p.data_layer == "person").all() and (p.source_geographic_resolution == "none").all()


@pytest.mark.data
@needs
def test_wearable_and_exam_tables():
    w = read_table(f"participant_wearable_features__{SID}")
    e = read_table(f"participant_exam_features__{SID}")
    assert len(w) == 73 and w.steps_mean_daily.notna().all()
    assert len(e) == 81 and set(ing.EXAM_COLUMNS) <= set(e.columns)
    assert (w.evidence_type == "person_device_measurement").all()


@pytest.mark.data
@needs
def test_group_medians_reproduce_published_table1():
    g = pd.read_csv(TABLES / f"{SID}_group_medians.csv").set_index("session_code")
    for s, pub in {"CON": 7153, "LC": 4718, "ME": 3704}.items():
        assert abs(g.loc[s, "Steps_median"] - pub) <= 1


@pytest.mark.data
@needs
def test_primary_has_complete_case_and_both_bounds():
    pr = pd.read_csv(TABLES / f"{SID}_primary.csv")
    assert list(pr.analysis) == ["primary", "primary bound", "primary bound"]
    assert pr.iloc[0].n_cases == 51 and pr.iloc[0].n_controls == 22
    assert (pr.iloc[1:].n_controls == 30).all()
    assert pr.iloc[1].auroc <= pr.iloc[0].auroc <= pr.iloc[2].auroc  # worst <= complete case <= best
    assert pr.iloc[0].n_permutations == an.N_PERM_PRIMARY


@pytest.mark.data
@needs
def test_all_locked_secondaries_reported():
    s = pd.read_csv(TABLES / f"{SID}_secondary.csv")
    assert sorted(s.analysis) == ["(a)", "(a)", "(b)", "(c)"]
    assert (TABLES / f"{SID}_head_to_head.csv").exists()
    d = pd.read_csv(TABLES / f"{SID}_cv_delta_auroc.csv")
    assert d.role.str.startswith("pre-specified").sum() == 1 and d.role.str.startswith("POST HOC").sum() == 1


@pytest.mark.data
@needs
def test_signatures_schema_and_ids():
    sig = read_table(f"phenotype_signatures__{SID}")
    assert sig.object_id.is_unique
    assert sig.object_id.str.match(rf"^signature:{SID}\|[a-z_]+\|.+").all()
    assert set(PROVENANCE_COLUMNS) <= set(sig.columns)
    prim = sig[sig.method.str.startswith("PRIMARY")]
    assert len(prim) == 1 and prim.iloc[0].effect_measure == "auroc"
    assert sig.evidence_level.isin({"supported_single_dataset", "nominal_single_dataset", "null_single_dataset",
                                    "descriptive"}).all()
    g = sig[sig.effect_measure == "hedges_g"]
    assert g.method.str.contains("EXPLORATORY").all()


@pytest.mark.data
@needs
def test_posthoc_missing_controls_diagnostic_is_labelled():
    d = pd.read_csv(TABLES / f"{SID}_posthoc_missing_controls.csv")
    assert d.analysis.str.startswith("POST HOC").all()
    fem = d.set_index("variable").loc["female (count)"]
    assert int(fem.missing_steps_n) == 8 and int(fem.missing_steps_value) == 7 and int(fem.with_steps_value) == 8


def test_report_does_not_claim_complete_case_sex_matching():
    src = an.write_report.__code__.co_consts
    text = " ".join(c for c in src if isinstance(c, str))
    assert "cannot carry" not in text and "sex carries none" not in text
