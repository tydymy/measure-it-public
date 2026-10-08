"""Tests for the temporal-holdout validation (Test 8; measure_it.scoring.temporal_holdout).

Plan: docs/ANALYSIS_PLAN_TEMPORAL_HOLDOUT.md. Pure tests run anywhere; tests marked `data` need the processed tables
(and, for the results tests, a completed `scoring_temporal_holdout` step).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from measure_it.config import RESULTS, TABLES
from measure_it.scoring import opportunity as O
from measure_it.scoring import temporal_holdout as TH
from measure_it.geography import features as F


# ------------------------------------------------------------------------------------------------ pure

def test_ref_percentile_matches_engine_when_reference_is_everything():
    rng = np.random.default_rng(0)
    x = rng.integers(0, 20, 200).astype(float)
    x[[3, 17]] = np.nan
    ref = np.ones(len(x), bool)
    np.testing.assert_allclose(TH.ref_percentile(x, ref), O.percentile_rank(x), equal_nan=True)


def test_ref_percentile_out_of_reference_row_is_ranked_as_if_added():
    x = np.array([1.0, 2.0, 2.0, 3.0, 2.0, 10.0])
    ref = np.array([True, True, True, True, False, False])
    out = TH.ref_percentile(x, ref)
    # reference rows: engine definition within the reference (1, 2, 2, 3)
    np.testing.assert_allclose(out[:4], O.percentile_rank(x[:4]))
    # a non-reference row: its percentile among reference + itself
    for i in (4, 5):
        pooled = np.r_[x[:4], x[i]]
        assert out[i] == pytest.approx(O.percentile_rank(pooled)[-1])


def test_fast_auroc_equals_sklearn_with_ties():
    rng = np.random.default_rng(1)
    for _ in range(5):
        s = rng.integers(0, 6, 300).astype(float)
        y = (rng.random(300) < 0.3).astype(int)
        assert TH.auroc(s, y) == pytest.approx(roc_auc_score(y, s))
    y = np.zeros(10, int)
    assert np.isnan(TH.auroc(np.arange(10.0), y))


def test_point_metrics_on_a_perfect_and_an_inverted_score():
    y = np.r_[np.ones(20), np.zeros(180)].astype(np.int8)
    s = np.r_[np.linspace(2, 3, 20), np.linspace(0, 1, 180)]
    ids = np.array([f"{i:05d}" for i in range(200)])
    m = TH.point_metrics(s, y, y.astype(float), ids)
    assert m["auroc"] == 1.0 and m["precision_at_10"] == 1.0 and m["top_decile_hit_rate"] == 1.0
    assert m["top_decile_lift"] == pytest.approx(10.0)
    m = TH.point_metrics(-s, y, y.astype(float), ids)
    assert m["auroc"] == 0.0 and m["precision_at_10"] == 0.0


def test_ties_in_the_ranking_are_broken_by_geo_id():
    s = np.array([1.0, 1.0, 1.0, 0.0])
    ids = np.array(["03", "01", "02", "04"])
    pos = TH._order_positions(s, ids)
    assert list(pos) == [2, 0, 1, 3]


def test_within_state_permutation_keeps_each_states_positives():
    rng = np.random.default_rng(2)
    n = 400
    groups = rng.choice(["01", "02", "06", "48"], n)
    y = (rng.random(n) < 0.25).astype(np.int8)
    s = rng.random(n)
    null = TH.perm_auroc(s, y, groups, 50, 3)
    assert null.shape == (50,) and np.all((null > 0) & (null < 1))
    # re-derive one permutation the same way and check group totals are preserved
    key = np.random.default_rng(9).random(n)
    g = pd.factorize(groups)[0]
    order, base = np.lexsort((key, g)), np.lexsort((np.arange(n), g))
    yp = np.empty(n, dtype=y.dtype)
    yp[base] = y[order]
    for st in np.unique(groups):
        assert yp[groups == st].sum() == y[groups == st].sum()


def test_bootstrap_ci_brackets_the_point_estimate():
    rng = np.random.default_rng(4)
    n = 500
    y = (rng.random(n) < 0.3).astype(np.int8)
    s = y + rng.normal(0, 1, n)
    ids = np.array([f"{i:05d}" for i in range(n)])
    idx = TH.boot_index(n, 300, 5)
    b = TH.boot_metrics(s, y, y.astype(float), ids, idx)
    lo, hi = TH._ci(b["auroc"])
    assert lo < TH.auroc(s, y) < hi


def test_engine_view_restores_the_engine_even_after_an_error():
    saved = (O.spatial, O._features, O._npi_state_groups)
    fake = TH.View(T=None, variant="x", S=None, npi_groups=pd.DataFrame(), features=pd.DataFrame())
    with pytest.raises(RuntimeError):
        with TH.engine_view(fake):
            assert O.spatial() is None
            raise RuntimeError("boom")
    assert (O.spatial, O._features, O._npi_state_groups) == saved


def test_reference_normalisation_restores_percentile_functions():
    saved = (O.percentile_rank, O.pct_rank_rows, F.percentile_rank)
    ref = np.array([True, True, False])
    with TH.reference_normalisation(ref):
        out = O.percentile_rank(np.array([1.0, 2.0, 3.0]))
        assert out[2] == pytest.approx(1.0)     # 3 among (1, 2, 3) -> 3/3
        assert out[0] == pytest.approx(0.5)     # reference row: 1 among (1, 2) -> 1/2
        # arrays of another length use the engine's function
        np.testing.assert_allclose(O.percentile_rank(np.array([1.0, 2.0])), [0.5, 1.0])
    assert (O.percentile_rank, O.pct_rank_rows, F.percentile_rank) == saved


def test_expected_directions_are_declared_for_every_component():
    for c in TH.COMPONENT_SCORES:
        assert c in TH.EXPECTED
    assert TH.EXPECTED["diagnostic_desert_pct"] == "negative"
    assert TH.EXPECTED["research_readiness_pct"] == "positive"


# ------------------------------------------------------------------------------------------------ data

@pytest.mark.data
def test_trial_partition_is_disjoint_and_dated():
    pre, post = TH.trials_before(2021), TH.trials_after(2021)
    assert not (pre & post)
    d = TH.trial_dates().set_index("nct_id")
    assert (d.loc[list(pre), "first_post"] <= pd.Timestamp("2021-12-31")).all()
    assert (d.loc[list(post), "first_post"] > pd.Timestamp("2021-12-31")).all()
    assert d["start"].notna().mean() > 0.99      # month-precision start dates parse ("2003-12")


@pytest.mark.data
def test_view_none_reproduces_the_stored_engine_output():
    chk = TH.view_check()
    assert chk["nan_pattern_equal"].all()
    assert (chk["max_abs_diff"].fillna(0) == 0).all(), chk


@pytest.mark.data
def test_trial_counter_with_everything_pre_period_reproduces_the_geography_module():
    c, s = TH.pool_counts(2100, "pre", ("me_cfs",), "me_cfs")
    f = O._features()
    fc = f[(f["condition_id"] == "me_cfs") & (f["geo_level"] == "county")].set_index("geo_id")
    got = c["trials_in_geo_or_within_50km_n"].reindex(fc.index).fillna(0)
    assert (got.to_numpy() == fc["trials_in_geo_or_within_50km_n"].to_numpy()).all()


@pytest.mark.data
def test_dated_view_removes_every_post_cutoff_input():
    v = TH.make_view(2021, "A_clean")
    d = TH.trial_dates().set_index("nct_id")["first_post"]
    assert (v.S.base.ft["nct_id"].astype(str).map(d) <= pd.Timestamp("2021-12-31")).all()
    assert (v.S.base.fn["fiscal_year"].astype(float) <= 2021).all()
    dates = TH._spatial_npi_dates()
    assert (v.S.npi_bits[dates > np.datetime64("2021-12-31")] == 0).all()
    assert (v.npi_groups["npi"].isin(set(TH._npi_state_groups_dated().query("enumeration_date > '2021-12-31'")["npi"]))
            .sum() == 0)
    exists = TH.facility_exists(2021, v.S.base.fac["facility_id"])
    assert not v.S.base.G[~exists].any()
    f = v.features
    lc = f[f["condition_id"] == "long_covid"]
    assert lc["burden_value"].isna().all() and (lc["burden_defined_level"] == "D").all()
    me = f[(f["condition_id"] == "me_cfs") & (f["geo_level"] == "county")].set_index("geo_id")
    cms = TH.cms_burden()
    cms = cms[(cms["geo_level"] == "county") & (cms["year"] == 2021.0)].set_index("geo_id")["value"]
    common = me.index[me["burden_value"].notna()]
    np.testing.assert_allclose(me.loc[common, "burden_value"], cms.reindex(common))
    # the engine really runs on the view: no Long COVID burden enters the pre-period set burden
    fr = TH.score_view(v)
    ranked = fr[fr["rank_eligible"]]
    assert (ranked["burden_members_contributing"] == "me_cfs").all()
    assert ranked["composite_equal"].notna().all()


@pytest.mark.data
def test_processed_rankings_table():
    from measure_it.store import read_table
    try:
        t = read_table(TH.TABLE)
    except FileNotFoundError:
        pytest.skip("scoring_temporal_holdout has not run")
    assert t["object_id"].is_unique and t["object_id"].str.startswith("temporal_holdout:").all()
    for c in ("source_name", "source_version", "retrieved_at", "evidence_type", "provenance_notes"):
        assert c in t and t[c].notna().all()
    assert {"current", "2021", "2020", "2019"} <= set(t["holdout_cutoff"])
    for c in TH.COMPONENT_SCORES + ["composite_equal", "rank_equal"]:
        assert c in t
    a = t[(t["holdout_cutoff"] == "2021") & (t["holdout_variant"] == "A_clean") & (t["holdout_query"] == TH.QUERY)
          & (t["geo_level"] == "county")]
    assert "post_o1_trials_pool_n" in a and a["post_o1_trials_pool_n"].notna().all()


@pytest.mark.data
def test_results_tables_and_report():
    p = TABLES / "temporal_holdout_metrics.csv"
    if not p.exists():
        pytest.skip("scoring_temporal_holdout has not run")
    m = pd.read_csv(p)
    prim = m[(m["T"] == TH.PRIMARY_T) & (m["variant"] == "A_clean") & (m["query"] == TH.QUERY)
             & (m["outcome"] == "O1")].set_index("score")
    for s in ("composite_equal", "composite_evidence_weighted", "research_readiness_pct", "diagnostic_desert_pct",
              "baseline_population_within_50km", "composite_burden_only"):
        r = prim.loc[s]
        assert r["auroc_ci_low"] <= r["auroc"] <= r["auroc_ci_high"]
        assert 0 < r["base_rate"] < 1
    assert "delta_auroc_vs_population_within_50km" in prim
    assert (prim["expected_direction"].notna()).all()
    md = (RESULTS / "TEMPORAL_HOLDOUT_RESULTS.md").read_text()
    assert "proxy outcome" in md.lower() and "not a clean holdout" in md.lower()
    for name in ("loso_summary", "rolling", "leakage", "population_adjusted", "state", "burden_persistence"):
        assert (TABLES / f"temporal_holdout_{name}.csv").exists(), name
