"""Tests for the metric -> translation link (measure_it.scoring.metric_link and its wiring into the scoring).

Plan: docs/ANALYSIS_PLAN_METRIC_LINK.md. Pure tests run anywhere; tests marked `data` need the processed tables.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from measure_it.config import RESULTS, TABLES, UNKNOWN
from measure_it.scoring import metric_link as ML
from measure_it.scoring import opportunity as O

BUNDLES = O.MEASUREMENTS


# ------------------------------------------------------------------------------------------------ pure

def test_sens_at_spec_threshold_rule():
    # 10 controls 0..9 (higher = case), cases 5..14: threshold = 9th smallest control (ceil(0.9 * 10) = 9) = 8
    y = np.r_[np.zeros(10), np.ones(10)]
    s = np.r_[np.arange(10), np.arange(5, 15)]
    se, sp, c = ML.sens_at_spec(y, s, 0.90)
    assert c == 8 and sp == pytest.approx(0.9) and se == pytest.approx(6 / 10)   # cases 9..14 exceed 8
    # 22 controls: ceil(19.8) = 20 -> 2 controls may exceed the threshold, specificity 20/22 >= 0.90
    y = np.r_[np.zeros(22), np.ones(5)]
    s = np.r_[np.arange(22), np.full(5, 100)]
    se, sp, _ = ML.sens_at_spec(y, s, 0.90)
    assert sp == pytest.approx(20 / 22) and sp >= 0.9 and se == 1.0


def test_sens_at_spec_bootstrap_is_seeded():
    rng = np.random.default_rng(1)
    y = np.r_[np.zeros(30), np.ones(30)]
    s = np.r_[rng.normal(0, 1, 30), rng.normal(1, 1, 30)]
    a, b = ML.sens_at_spec_boot(y, s, n_boot=200), ML.sens_at_spec_boot(y, s, n_boot=200)
    assert a == b
    assert a["op_sensitivity_ci_low"] <= a["op_sensitivity"] <= a["op_sensitivity_ci_high"]
    assert a["op_specificity"] >= 0.9


def test_wilson_hanley_and_evidence_score():
    lo, hi = ML.wilson(0.9, 100)
    assert 0.8 < lo < 0.9 < hi < 0.96
    assert all(np.isnan(ML.wilson(0.9, np.nan)))
    lo, hi = ML.hanley_mcneil_ci(0.9, 66, 20)
    assert lo < 0.9 < hi
    assert ML.evidence_score(0.75, 1) == (pytest.approx(0.5), pytest.approx(0.5))
    assert ML.evidence_score(0.75, 3) == (pytest.approx(0.25), pytest.approx(0.5))
    assert ML.evidence_score(0.40, 1) == (0.0, 0.0)            # below chance -> 0, never negative
    assert all(np.isnan(ML.evidence_score(np.nan, 1)))        # no AUROC -> UNKNOWN, not 0
    assert all(np.isnan(ML.evidence_score(0.9, 4)))


def test_draw_logit_stays_in_unit_interval_and_centres():
    x = ML.draw_logit(0.55, 0.41, 0.76, np.random.default_rng(0), 5000)
    assert ((x > 0) & (x < 1)).all() and abs(np.median(x) - 0.55) < 0.02
    assert np.isnan(ML.draw_logit(np.nan, 0.1, 0.2, np.random.default_rng(0), 3)).all()


def test_evidence_weighted_set_is_declared():
    ws = O.weight_sets()
    assert O.EVIDENCE_WEIGHT_SET in ws and "equal" in ws
    ew = ws[O.EVIDENCE_WEIGHT_SET]
    assert set(ew) == set(O.COMPONENTS) | set(O.EVIDENCE_COMPONENTS)
    assert sum(ew.values()) == pytest.approx(1.0)
    assert O.uses_evidence(ew) and not O.uses_evidence(ws["equal"])
    assert O.EVIDENCE_WEIGHT_SET not in O.weight_sets(include_extra=False)
    assert O.default_weights() == ws["equal"]               # the default is unchanged (continuity)


def test_published_map_rows_are_diagnostic_accuracy_rows_of_grid_conditions():
    grid = set(ML.SINGLE_CONDITIONS)
    for pid, (target, cls, _) in ML.PUBLISHED_MAP.items():
        assert target in grid and cls


# ------------------------------------------------------------------------------------------------ data

@pytest.fixture(scope="module")
def perf():
    from measure_it.store import read_table, table_exists
    if not table_exists(ML.TABLE):
        pytest.skip("measurement_performance not built")
    return read_table(ML.TABLE)


@pytest.fixture(scope="module")
def opp():
    from measure_it.store import read_table, table_exists
    if not table_exists(O.TABLE):
        pytest.skip("deployment_opportunities not built")
    return read_table(O.TABLE)


@pytest.mark.data
def test_performance_table_grid_ids_and_provenance(perf):
    from measure_it.provenance import PROVENANCE_COLUMNS
    assert set(PROVENANCE_COLUMNS) <= set(perf.columns)
    assert perf["object_id"].is_unique
    assert (perf["object_id"] == "measurement_performance:" + perf["condition_id"] + "|" + perf["measurement_id"]).all()
    assert len(perf) == len(ML.grid())
    assert set(perf["performance_status"]) <= {"known", "partial", UNKNOWN}
    unk = perf[perf["performance_status"] == UNKNOWN]
    assert unk["measurement_evidence"].isna().all() and unk["auroc"].isna().all() and (unk["quality_tier"] == 4).all()


@pytest.mark.data
def test_nothing_is_imputed_across_conditions_or_measurements(perf):
    rec = ML.records()
    by_id = rec.set_index("record_id")
    for r in perf[perf["performance_status"] != UNKNOWN].itertuples():
        src = by_id.loc[r.selected_record_id]
        assert sorted(json.loads(src["target_members"])) == sorted(ML._members(r.condition_id))
        assert src["primary_class"] in r.measurement_members.split(";")
        assert src["auroc"] == pytest.approx(r.auroc)
    # capillaroscopy and the demo cluster have no record at all
    assert (perf.loc[perf["measurement_id"] == "nailfold_capillaroscopy", "performance_status"] == UNKNOWN).all()
    assert (perf.loc[perf["condition_id"] == "autonomic_activity_invisible_illness", "performance_status"]
            == UNKNOWN).all()
    assert (perf.loc[perf["condition_id"] == "dysautonomia", "performance_status"] == UNKNOWN).all()


@pytest.mark.data
def test_scored_records_of_the_primary_query(perf):
    p = perf.set_index(["condition_id", "measurement_id"])
    w = p.loc[(O.QUERY_SET_ID, "wearable_autonomic_activity_monitoring")]
    assert w["selected_record_id"] == "OWN-MM-STEPS-POOLED" and w["quality_tier"] == 1
    assert w["auroc"] == pytest.approx(0.8307, abs=1e-3) and w["performance_status"] == "known"
    assert w["op_specificity"] >= ML.TARGET_SPECIFICITY
    lc = p.loc[("long_covid", "wearable_autonomic_activity_monitoring")]
    assert lc["selected_record_id"] == "OWN-MM-STEPS-LC" and bool(lc["evidence_conflict"])
    assert "OWN-UW-RHR" in lc["evidence_conflict_note"]
    uw = p.loc[("long_covid", "wearable_heart_rate")]
    assert uw["selected_record_id"] == "OWN-UW-RHR" and uw["quality_tier"] == 2
    assert uw["measurement_evidence"] == 0.0          # AUROC lower bound below 0.5: no demonstrated discrimination
    pots = p.loc[("pots", "wearable_autonomic_activity_monitoring")]
    assert pots["selected_record_id"] == "PDE-007" and pots["quality_tier"] == 3   # published, flagged
    assert p.loc[("pots", "autonomic_function_testing"), "performance_status"] == "partial"
    assert p.loc[("me_cfs", "autonomic_function_testing"), "performance_status"] == UNKNOWN


@pytest.mark.data
def test_existing_muscle_me_aurocs_are_reproduced():
    from measure_it.wearables.muscle_me_steps_analysis import fixed_direction_auroc, load_frame
    mm = load_frame()
    prim = pd.read_csv(TABLES / "charlton_lc_mecfs_cpet_source_primary.csv")
    d = mm[mm["Steps"].notna()]
    r = fixed_direction_auroc(d["is_patient"].to_numpy(), d["Steps"].to_numpy(float))
    p = prim[prim["analysis"] == "primary"].iloc[0]
    assert r["auroc"] == pytest.approx(p["auroc"]) and r["ci_low"] == pytest.approx(p["ci_low"])


@pytest.mark.data
def test_unknown_performance_is_never_ranked_under_evidence_weighted(opp):
    ew = f"composite_{O.EVIDENCE_WEIGHT_SET}"
    not_known = opp["measurement_performance_status"] != "known"
    assert opp.loc[not_known, ew].isna().all() and opp.loc[not_known, f"rank_{O.EVIDENCE_WEIGHT_SET}"].isna().all()
    assert opp.loc[not_known, "evidence_weighted_status"].str.startswith("not ranked").all()
    # never scored 0: the components are NaN, not 0
    assert opp.loc[not_known, "measurement_evidence"].isna().all()
    assert opp.loc[not_known, "expected_yield"].isna().all()
    # the equal-weight ranking is unchanged in kind for them (continuity)
    assert opp.loc[not_known & opp["rank_eligible"].astype(bool), "rank_equal"].notna().all()
    cap = opp[(opp["measurement_id"] == "nailfold_capillaroscopy")]
    assert cap[f"rank_{O.EVIDENCE_WEIGHT_SET}_joint"].isna().all() and cap["rank_equal_joint"].notna().any()


@pytest.mark.data
def test_evidence_weighted_composite_formula(opp):
    ws = O.weight_sets()[O.EVIDENCE_WEIGHT_SET]
    sub = opp[opp[f"composite_{O.EVIDENCE_WEIGHT_SET}"].notna()]
    cols = {k: (O.PCT[k] if k in O.PCT else k) for k in ws}
    X = sub[[cols[k] for k in ws]].to_numpy(float)
    W = np.array([ws[k] for k in ws])
    avail = ~np.isnan(X)
    comp = (np.where(avail, X, 0) * W).sum(1) / (avail * W).sum(1)
    assert np.allclose(comp, sub[f"composite_{O.EVIDENCE_WEIGHT_SET}"].to_numpy(float))
    # measurement_evidence is constant within a combination
    assert (sub.groupby(["condition_id", "measurement_id", "geo_level"])["measurement_evidence"].nunique() == 1).all()


@pytest.mark.data
def test_expected_yield_factors(opp):
    known = opp[(opp["measurement_performance_status"] == "known") & opp["expected_yield"].notna()]
    pc = known[known["expected_yield_basis"] == "percentile"]
    assert np.allclose(pc["expected_yield"], pc["burden_pct"] * pc["op_sensitivity"] * pc["reach"])
    assert np.allclose(pc["expected_yield_index"], pc["expected_yield"])
    assert np.allclose(pc["expected_false_positives_upper"],
                       (1 - pc["op_specificity"]) * pc["adults_18plus"] * pc["reach"])
    assert pc["expected_detectable_cases"].isna().all()        # no count where no defensible count exists
    cnt = opp[opp["expected_yield_basis"] == "count"]
    assert set(cnt["condition_id"]) == {"long_covid"} and set(cnt["geo_level"]) == {"state"}
    k = cnt[cnt["measurement_performance_status"] == "known"]
    bc = k["burden_value"] / 100 * k["adults_18plus"]
    assert np.allclose(k["burden_count"], bc)
    assert np.allclose(k["expected_detectable_cases"], bc * k["op_sensitivity"] * k["reach"])
    assert np.allclose(k["expected_false_positives"], (1 - k["op_specificity"]) * (k["adults_18plus"] - bc) * k["reach"])
    assert (k["expected_ppv"] < 0.5).all()           # low prevalence: most flags are false positives
    # county long-COVID burden is inherited or modelled (BRFSS small-area): never a county count
    lc = opp[(opp["condition_id"] == "long_covid") & (opp["geo_level"] == "county")]
    assert (lc["expected_yield_basis"] == "percentile").all() and lc["burden_count"].isna().all()
    # intervals bracket the observed value where both exist
    m = known["expected_yield_mc_p05"].notna()
    obs = np.where(known["expected_yield_basis"] == "count", known["expected_detectable_cases"],
                   known["expected_yield_index"])[m.to_numpy()]
    assert (known.loc[m, "expected_yield_mc_p05"].to_numpy() <= obs * 1.0001 + 1e-9).mean() > 0.95


@pytest.mark.data
def test_level_d_with_known_performance_drops_yield_weight(opp):
    p = opp[(opp["condition_id"] == "pots") & (opp["measurement_id"] == "wearable_autonomic_activity_monitoring")]
    assert (p["expected_yield_basis"] == "level_D").all() and p["expected_yield"].isna().all()
    assert p["expected_yield_excluded_level_D"].all()
    assert p.loc[p["rank_eligible"].astype(bool), f"rank_{O.EVIDENCE_WEIGHT_SET}"].notna().all()


@pytest.mark.data
def test_uncertainties_name_the_performance_tier(opp):
    row = opp[(opp["condition_id"] == O.QUERY_SET_ID) & (opp["measurement_id"] == O.PRIMARY["measurement"])].iloc[0]
    u = " ".join(json.loads(row["uncertainties"]))
    assert "tier 1" in u and "healthy" in u and "not predictions of diagnoses" in u
    cap = opp[(opp["measurement_id"] == "nailfold_capillaroscopy")].iloc[0]
    assert "UNKNOWN" in " ".join(json.loads(cap["uncertainties"]))


@pytest.mark.data
def test_recommendation_carries_performance_and_yield():
    from measure_it.scoring.recommend import rank_deployment_opportunities
    r = rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=2,
                                      weight_set=O.EVIDENCE_WEIGHT_SET)
    assert r["status"] == "ok" and r["measurement_performance"]["selected_record_id"] == "OWN-MM-STEPS-POOLED"
    rec = r["recommendations"][0]
    tech = rec["technology"]
    assert tech["measurement_performance"]["performance_status"] == "known"
    ey = rec["geography"]["expected_yield"]
    assert ey["basis"] == "percentile" and ey["expected_false_positives_upper"] > 0
    assert "not predictions of diagnoses" in ey["language"]
    assert "EVIDENCE_WEIGHTED" in rec["geography"]["monte_carlo_basis"]
    assert rec["geography"]["rank"] == 1
    assert any("tier 1" in u for u in rec["uncertainties"])
    assert "measurement_performance:" + O.QUERY_SET_ID + "|" + O.PRIMARY["measurement"] in \
        rec["provenance"]["object_ids"]
    cap = rank_deployment_opportunities("Long COVID or ME/CFS", "capillaroscopy", top_n=2,
                                        weight_set=O.EVIDENCE_WEIGHT_SET, with_context=False)
    assert cap["status"] == UNKNOWN and cap["measurement_performance"]["performance_status"] == UNKNOWN
    eq = rank_deployment_opportunities("Long COVID or ME/CFS", "capillaroscopy", top_n=1, with_context=False)
    assert eq["status"] == "ok"
    assert "carries no measurement-performance information" in eq["recommendations"][0]["recommended_next_step"]


@pytest.mark.data
def test_trace_resolves_performance_ids(perf):
    from measure_it.mcp.trace import trace_evidence
    oid = f"measurement_performance:{O.QUERY_SET_ID}|wearable_autonomic_activity_monitoring"
    r = trace_evidence(oid)
    assert r["status"] == "ok" and r["object_id"] == oid


@pytest.mark.data
def test_result_tables_and_reports():
    t5 = pd.read_csv(TABLES / "test5_contrasts.csv")
    t5e = t5[t5["contrast"] == O.EVIDENCE_WEIGHT_SET]
    assert len(t5e) and set(t5e["status"].str.split(":").str[0]) <= {"ranked", "not ranked under evidence_weighted"}
    s5 = pd.read_csv(TABLES / "test7_sensitivity_s5_evidence.csv")
    cap = s5[(s5["measurement_id_b"] == "nailfold_capillaroscopy") & (s5["weight_set"] == O.EVIDENCE_WEIGHT_SET)]
    assert cap["status"].str.startswith("not ranked").all()
    jc = pd.read_csv(TABLES / "scoring_joint_ranking_composition.csv")
    ew = jc[(jc["weight_set"] == O.EVIDENCE_WEIGHT_SET)]
    assert (ew["n_nailfold_capillaroscopy"] == 0).all()
    txt = (RESULTS / "SCORING_RESULTS.md").read_text()
    assert "## 12. The metric -> translation link" in txt and "not predictions of diagnoses" in txt
    ex = json.loads((RESULTS / "example_deployment_recommendation.json").read_text())
    assert "measurement_performance" in ex["recommendations"][0]["technology"]
    assert "expected_yield" in ex["recommendations"][0]["geography"]
