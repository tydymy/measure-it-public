"""Tests for measure_it.scoring (opportunity, controls, sensitivity, recommend).

Pure-function tests run anywhere; tests marked `data` need data/processed (run
`uv run python -m measure_it.scoring.report` first).
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import pytest

from measure_it.config import RESULTS, TABLES, UNKNOWN
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table
from measure_it.scoring import controls as C
from measure_it.scoring import opportunity as O
from measure_it.scoring import sensitivity as SENS

BANNED = [r"\bproves?\b", r"diagnosed by ai", r"\bpatient has\b", r"\bbest clinic\b", r"definitive biomarker",
          r"confirms mechanism", r"optimal treatment"]
SPEC_KEYS = {"condition", "phenotype", "measurement", "technology", "geography", "candidate_sites",
             "research_evidence", "molecular_context", "uncertainties", "provenance", "recommended_next_step"}


# ------------------------------------------------------------------------------------------------------------------
# pure functions
# ------------------------------------------------------------------------------------------------------------------

def test_percentile_rank_matches_pandas_and_keeps_nan():
    x = np.array([3.0, 1.0, np.nan, 2.0, 2.0])
    r = O.percentile_rank(x)
    assert np.isnan(r[2])
    assert np.allclose(r[[0, 1, 3, 4]], [1.0, 0.25, 0.625, 0.625])


def test_pct_rank_rows_equals_rowwise_percentile_rank():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(4, 30))
    X[:, 5] = np.nan
    R = O.pct_rank_rows(X)
    for i in range(4):
        assert np.allclose(R[i], O.percentile_rank(X[i]), equal_nan=True)


def test_weighted_composite_renormalises_missing_component():
    comps = {"burden": np.array([np.nan, 1.0]), "vulnerability": np.array([0.5, 0.0])}
    c = O.weighted_composite(comps, {"burden": 0.5, "vulnerability": 0.5}, ["burden", "vulnerability"])[0]
    assert c[0] == pytest.approx(0.5)            # burden missing -> vulnerability alone
    assert c[1] == pytest.approx(0.5)
    bo = O.weighted_composite(comps, {"burden": 1.0}, ["burden"])[0]
    assert np.isnan(bo[0])                         # burden_only undefined without burden
    W = np.array([[1.0, 0.0], [0.0, 1.0]])
    cc = O.weighted_composite(comps, W, ["burden", "vulnerability"])
    assert cc.shape == (2, 2) and cc[1, 0] == pytest.approx(0.5) and cc[0, 1] == pytest.approx(1.0)


def test_rank_rows_ties_by_column_order_and_eligibility():
    c = np.array([0.5, 0.9, 0.5, np.nan, 0.95])
    elig = np.array([True, True, True, True, False])
    r = O.rank_rows(c, elig)[0]
    assert r[1] == 1 and r[0] == 2 and r[2] == 3
    assert np.isnan(r[3]) and np.isnan(r[4])


def test_kendalls_w_bounds():
    R = np.tile(np.arange(1, 11), (5, 1)).astype(float)
    assert SENS.kendalls_w(R) == pytest.approx(1.0)
    R2 = np.vstack([np.arange(1, 11), np.arange(10, 0, -1)]).astype(float)
    assert SENS.kendalls_w(R2) == pytest.approx(0.0)


def test_size_class_and_cms_sd():
    assert O.size_class_n("Medicare FFS beneficiaries; size class 1,000-4,999") == pytest.approx(2236.0)
    assert O.size_class_n("Medicare FFS beneficiaries; size class 10,000+") == pytest.approx(10000.0)
    assert np.isnan(O.size_class_n(None))
    sd = O.cms_sd_pp(23.0, 10000.0)
    assert sd == pytest.approx(np.sqrt(0.23 * 0.77 / 10000 * 1e4 + 1 / 12))


def test_between_state_tau_method_of_moments():
    v = np.array([1.0, 2.0, 3.0, 4.0])
    assert O.between_state_tau(v, np.zeros(4)) == pytest.approx(np.std(v, ddof=1))
    assert O.between_state_tau(v, np.full(4, 10.0)) == 0.0


def test_tie_averaged_overlap_and_jaccard():
    s = np.array([3.0, 2.0, 2.0, 2.0, 1.0])
    elig = np.ones(5, bool)
    ov, tied = C.tie_averaged_overlap(s, elig, np.array([0, 1]), 2)
    assert tied == 3 and ov == pytest.approx(1 + 1 / 3)
    assert C.jaccard([1, 2], [2, 3]) == pytest.approx(1 / 3)


def test_within_group_permutation_stays_in_group():
    rng = np.random.default_rng(1)
    g = np.array(list("aabbbcc") * 20)
    p = C.permutation(len(g), g, rng)
    assert sorted(p) == list(range(len(g)))
    assert (g[p] == g).all()


def test_weight_sets_sum_to_one():
    for name, w in O.weight_sets().items():
        assert sum(w.values()) == pytest.approx(1.0), name


def test_saturation_adjusted_lives_in_config_and_inherited_multiplier_is_explicit():
    from measure_it.config import load_config
    cfg = load_config("scoring")
    assert "saturation_adjusted" in cfg["weight_sets"] and not hasattr(O, "EXTRA_WEIGHT_SETS")
    assert all(v == 1 / 6 for v in O.weight_sets()["saturation_adjusted"].values())   # exact float64 1/6
    assert "saturation_adjusted" not in O.weight_sets(include_extra=False)
    assert O.inherited_multiplier() == cfg["inherited_burden_uncertainty_multiplier"] == O.multiplier("C")


def test_condition_specs():
    s = O.condition_spec(O.QUERY_SET_ID)
    assert s["kind"] == "condition_set" and s["members"] == ["long_covid", "me_cfs"]
    d = O.condition_spec("autonomic_activity_invisible_illness")
    assert set(d["members"]) == {"long_covid", "me_cfs", "pots", "dysautonomia"}
    assert O.condition_spec("pots")["members"] == ["pots"]
    with pytest.raises(KeyError):
        O.condition_spec("not_a_condition")


# ------------------------------------------------------------------------------------------------------------------
# processed outputs
# ------------------------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def opp():
    from measure_it.store import read_table, table_exists
    if not table_exists(O.TABLE):
        pytest.skip("deployment_opportunities not built")
    return read_table(O.TABLE)


@pytest.mark.data
def test_table_shape_ids_and_provenance(opp):
    assert set(PROVENANCE_COLUMNS) <= set(opp.columns)
    assert opp["object_id"].is_unique
    assert opp["object_id"].str.match(r"^opportunity:[a-z_+]+\|[a-z_]+\|\d{2,5}$").all()
    assert (opp["data_layer"] == "derived").all()
    n_cond = len(O.all_condition_ids())
    assert len(opp) == n_cond * len(O.MEASUREMENTS) * (3144 + 51)
    raw = {"burden": "burden_value", "vulnerability": "svi_overall"}
    for c in O.COMPONENTS:
        assert raw.get(c, c) in opp.columns
        assert O.PCT[c] in opp.columns
        v = opp[O.PCT[c]].dropna()
        assert ((v > 0) & (v <= 1)).all()
    for c in ["burden_value", "svi_overall", "diagnostic_desert", "clinic_capacity", "research_readiness",
              "technology_saturation", "technology_saturation_pct", "burden_evidence_level", "burden_inherited",
              "burden_source_resolution", "rank_mc_p05", "rank_mc_p95", "uncertainties"]:
        assert c in opp.columns, c
    for ws in O.weight_sets():
        assert f"composite_{ws}" in opp.columns and f"rank_{ws}" in opp.columns


@pytest.mark.data
def test_composite_is_weighted_mean_of_components(opp):
    w = O.default_weights()
    X = opp[[O.PCT[k] for k in O.COMPONENTS]].to_numpy(float)
    W = np.array([w[k] for k in O.COMPONENTS])
    avail = ~np.isnan(X)
    comp = (np.where(avail, X, 0) * W).sum(1) / (avail * W).sum(1)
    comp[opp["burden_incomplete"].to_numpy(bool)] = np.nan      # incomplete burden: no composite (deviation 7)
    assert np.allclose(comp, opp["composite_equal"].to_numpy(float), equal_nan=True)


@pytest.mark.data
def test_level_d_excluded_and_flagged(opp):
    d = opp[opp["condition_id"].isin(["pots", "dysautonomia"])]
    assert d["burden_excluded_level_D"].all() and d["burden_pct"].isna().all()
    assert d["rank_burden_only"].isna().all()
    assert (d["burden_evidence_level"] == "D").all()
    assert d["diagnostic_desert_variant"].str.startswith("access_only").all()
    assert d["uncertainties"].str.contains("level D").all()


@pytest.mark.data
def test_inherited_burden_is_labelled(opp):
    lc = opp[(opp["condition_id"] == "long_covid") & (opp["geo_level"] == "county")]
    st = opp[(opp["condition_id"] == "long_covid") & (opp["geo_level"] == "state")]
    assert not st["burden_inherited"].any()
    assert (st["burden_inherited_uncertainty_multiplier"] == 1.0).all()
    assert (opp["mc_inherited_multiplier"] == O.inherited_multiplier()).all()
    if lc["burden_inherited"].any():   # long_covid_burden: hps_inherited (the rule before 2026-10-07)
        assert lc["burden_inherited"].all()
        assert (lc["source_geographic_resolution"] == "state").all()
        assert lc["uncertainties"].str.contains("inherited").all()
        assert (lc.groupby(["measurement_id", "state_fips"])["burden_value"].nunique() == 1).all()
        # the explicit inherited-burden multiplier is applied (and exposed) only on inherited rows
        assert (lc["burden_inherited_uncertainty_multiplier"] == O.inherited_multiplier()).all()
        assert lc["burden_mc_sd_basis"].str.contains("inherited-burden multiplier").all()
    else:                               # brfss_sae: a modelled county estimate, labelled as such in every row
        assert lc["burden_source_name"].str.contains("small-area").all()
        assert lc["uncertainties"].str.contains("modelled small-area estimate").all()
        assert (lc["burden_inherited_uncertainty_multiplier"] == 1.0).all()
        assert (lc.groupby(["measurement_id", "state_fips"])["burden_value"].nunique() > 1).any()


@pytest.mark.data
def test_monte_carlo_intervals(opp):
    assert (opp["mc_n_draws"] >= 500).all()
    e = opp[opp["rank_eligible"]]
    assert (e["rank_mc_p05"] <= e["rank_mc_p50"]).all() and (e["rank_mc_p50"] <= e["rank_mc_p95"]).all()
    small = opp[(opp["geo_level"] == "county") & ~opp["rank_eligible"]]
    assert small["rank_equal"].isna().all() and small["rank_mc_p05"].isna().all()
    inc = opp["burden_incomplete"].astype(bool)
    assert small.loc[~inc.loc[small.index], "rank_equal_incl_small"].notna().all()
    assert opp.loc[inc, "rank_equal_incl_small"].isna().all()


@pytest.mark.data
def test_set_burden_is_renormalised_mean_of_members_and_capacity_condition_free(opp):
    key = ["measurement_id", "geo_level", "geo_id"]
    s = opp[opp["condition_id"] == O.QUERY_SET_ID].set_index(key)
    lc = opp[opp["condition_id"] == "long_covid"].set_index(key)["burden_pct"]
    me = opp[opp["condition_id"] == "me_cfs"].set_index(key)["burden_pct"]
    # raw set burden = mean of member percentiles, only where BOTH defined members have a value (deviation 7:
    # no set burden rests on one member); normalised = pr() of it over the complete regions, per measurement x level
    exp = pd.concat([lc, me], axis=1).mean(axis=1, skipna=False)
    assert np.allclose(s["burden_set_member_mean_pct"].to_numpy(float), exp.reindex(s.index).to_numpy(float),
                       equal_nan=True)
    renorm = s.groupby(level=["measurement_id", "geo_level"])["burden_set_member_mean_pct"].rank(
        pct=True, method="average")
    assert np.allclose(s["burden_pct"].to_numpy(float), renorm.to_numpy(float), equal_nan=True)
    # the normalised set burden has the full spread of a percentile (like every other component)
    e = s[s["rank_eligible"]].xs(("wearable_autonomic_activity_monitoring", "county"), level=[0, 1])
    assert e["burden_pct"].std() > 0.26
    single = opp[opp["condition_kind"] == "condition"]
    assert single["burden_set_member_mean_pct"].isna().all()
    cc = opp.groupby(key)["clinic_capacity"].nunique(dropna=False)
    assert (cc == 1).all()


def test_set_burden_renormalisation_leaves_single_conditions_unchanged():
    x = np.array([[5.0, 1.0, np.nan, 3.0, 3.0, 9.0]])
    p = O.pct_rank_rows(x)
    assert np.allclose(O.pct_rank_rows(p), p, equal_nan=True)  # pr(pr(x)) == pr(x)


@pytest.mark.data
def test_rank_deployment_argument_validation():
    from measure_it.scoring.recommend import rank_deployment_opportunities
    for bad in (0, -3, 2.5, "ten", None):
        r = rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=bad,
                                          with_context=False)
        assert r["status"] == UNKNOWN and r["recommendations"] == [], bad
    r = rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=5,
                                      weight_set="burden_led", include_small_population=True, with_context=False)
    assert r["status"] == "ok" and r["rank_column"] == "rank_burden_led_incl_small"
    assert [x["geography"]["rank"] for x in r["recommendations"]] == [1, 2, 3, 4, 5]
    r = rank_deployment_opportunities("POTS or dysautonomia", "wearable autonomic monitoring", top_n=2,
                                      with_context=False)
    assert r["status"] == "ok"
    for rec in r["recommendations"]:
        assert rec["geography"]["burden_evidence_level"] == "D"
        assert not any("rests on the other member" in u for u in rec["uncertainties"])
        assert any("no member of the set has a usable burden" in u for u in rec["uncertainties"])


@pytest.mark.data
def test_on_the_fly_ranking_cites_only_stored_rows():
    """An ad hoc set is scored on the fly: its opportunity id has no stored row, so provenance must not cite it; it
    cites traceable component ids (geo, conditions, measurement, the members' stored opportunity rows) instead."""
    from measure_it.mcp.trace import trace_evidence
    from measure_it.scoring.recommend import rank_deployment_opportunities
    r = rank_deployment_opportunities("POTS or dysautonomia", "wearable autonomic monitoring", top_n=2,
                                      with_context=False)
    assert r["status"] == "ok" and r["opportunity_rows_stored"] is False
    for rec in r["recommendations"]:
        row = rec["provenance"]["opportunity_row"]
        assert row["stored"] is False and "not a stored row" in row["note"]
        ids = rec["provenance"]["object_ids"]
        assert row["object_id"] not in ids
        gid = rec["geography"]["fips"]
        assert {f"geo:{gid}", "condition:pots", "condition:dysautonomia",
                "measurement:wearable_autonomic_activity_monitoring"} <= set(ids)
        opp_ids = [i for i in ids if i.startswith("opportunity:")]
        assert sorted(opp_ids) == sorted(f"opportunity:{c}|wearable_autonomic_activity_monitoring|{gid}"
                                         for c in ("dysautonomia", "pots"))
        for i in opp_ids + [f"geo:{gid}"]:
            assert trace_evidence(i)["status"] == "ok", i
        assert trace_evidence(row["object_id"])["status"] == UNKNOWN
    assert all(e["object_id"] is None for e in r["excluded_incomplete_burden"])
    # the tool facade's envelope provenance collects object ids from the whole answer: none may be an unstored row
    from measure_it import tools as T
    env = T.rank_deployment_opportunities("POTS or dysautonomia", "wearable autonomic monitoring", top_n=2)
    assert env["status"] == "ok"
    ids = T.collect_object_ids(env["data"])
    assert not any(i.startswith("opportunity:dysautonomia+pots|") for i in ids)
    assert not any(i.startswith("opportunity:dysautonomia+pots|") for i in env["provenance"]["object_ids"])
    s = rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=1,
                                      with_context=False)
    row = s["recommendations"][0]["provenance"]["opportunity_row"]
    assert s["opportunity_rows_stored"] and row["stored"] and s["recommendations"][0]["provenance"]["object_ids"][0] \
        == row["object_id"] and trace_evidence(row["object_id"])["status"] == "ok"


@pytest.mark.data
def test_next_step_says_which_rank_the_interval_belongs_to():
    """The Monte Carlo interval is drawn around the equal weights among population-eligible regions; a rank under
    another weight set (or among counties incl. small ones) is never presented as the interval's rank."""
    from measure_it.scoring.recommend import rank_deployment_opportunities
    eq = rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=2,
                                       with_context=False)["recommendations"][0]
    assert "under the 'equal' weights; 5th-95th percentile rank interval" in eq["recommended_next_step"]
    assert "equal-weight rank" not in eq["recommended_next_step"]
    bl = rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=2,
                                       weight_set="burden_led", with_context=False)["recommendations"]
    for rec in bl:
        g, txt = rec["geography"], rec["recommended_next_step"]
        assert f"rank {g['rank']:.0f} of" in txt and "under the 'burden_led' weights" in txt
        if g["rank_interval_5_95"][0] is not None:
            assert "belongs to the equal-weight rank" in txt and "not an interval for the 'burden_led'-weight rank" in txt
            eqr = g["ranks_under_weight_sets"]["equal"]
            assert f"({eqr:.0f})" in txt
        assert "not the rank under 'burden_led'" in g["monte_carlo_basis"]
    sm = rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=3,
                                       include_small_population=True, with_context=False)["recommendations"]
    assert all("among all counties incl. small ones" in x["recommended_next_step"] for x in sm)


@pytest.mark.data
def test_candidate_sites_carry_coordinates():
    p = RESULTS / "example_deployment_recommendation.json"
    if not p.exists():
        pytest.skip("demo not run")
    res = json.loads(p.read_text())
    sites = [s for rec in res["recommendations"] for s in (rec["candidate_sites"] if isinstance(rec["candidate_sites"],
                                                                                            list) else [])]
    assert sites and all(isinstance(s["lat"], float) and isinstance(s["lon"], float) and s["geocode_precision"]
                         for s in sites)


@pytest.mark.data
def test_desert_matches_geography_module(opp):
    from measure_it.store import read_table
    f = read_table("geo_condition_features", columns=["geo_id", "geo_level", "condition_id", "diagnostic_desert",
                                                      "diagnostic_desert_access_only", "in_analysis_universe"])
    f = f[f["in_analysis_universe"]]
    for cid in ("long_covid", "me_cfs", "pots"):
        o = opp[(opp["condition_id"] == cid) & (opp["measurement_id"] == O.PRIMARY["measurement"])]
        m = o.merge(f[f["condition_id"] == cid], on=["geo_id", "geo_level"], suffixes=("", "_geo"))
        ref = np.where(m["diagnostic_desert_geo"].isna(), m["diagnostic_desert_access_only_geo"],
                       m["diagnostic_desert_geo"])
        assert np.allclose(m["diagnostic_desert"].to_numpy(float), ref.astype(float), equal_nan=True), cid


@pytest.mark.data
def test_pools_match_find_candidate_clinics():
    from measure_it.facilities import matching as M
    cb = O.make_combo("long_covid", O.PRIMARY["measurement"], "county")
    for gid in ("06073", "21095", "37165"):
        i = int(np.flatnonzero(cb.geo["geo_id"].to_numpy() == gid)[0])
        r = M.find_candidate_clinics(gid, "long_covid", O.PRIMARY["measurement"])
        assert int(cb.fm["pool_facilities_n"][i]) == r["n_facilities_in_pool"]


@pytest.mark.data
def test_resolve_condition_aliases():
    from measure_it.scoring.recommend import resolve_condition_query
    assert resolve_condition_query("Long COVID or ME/CFS")["condition_id"] == O.QUERY_SET_ID
    assert resolve_condition_query("ME/CFS or Long COVID")["condition_id"] == O.QUERY_SET_ID
    assert resolve_condition_query("POTS")["condition_id"] == "pots"
    assert resolve_condition_query("not a real illness xyz")["status"] == UNKNOWN


@pytest.mark.data
def test_rank_deployment_unknowns():
    from measure_it.scoring.recommend import rank_deployment_opportunities
    r = rank_deployment_opportunities("Long COVID", "wearable autonomic monitoring", weight_set="no_such_set")
    assert r["status"] == UNKNOWN and r["recommendations"] == []
    r = rank_deployment_opportunities("POTS", "autonomic_function_testing", weight_set="burden_only",
                                      with_context=False)
    assert r["status"] == UNKNOWN  # level D: no burden-only ranking


@pytest.mark.data
def test_example_recommendation_schema_and_language():
    p = RESULTS / "example_deployment_recommendation.json"
    if not p.exists():
        pytest.skip("demo not run")
    res = json.loads(p.read_text())
    assert res["status"] == "ok" and len(res["recommendations"]) == 10
    assert res["condition_resolution"]["condition_id"] == O.QUERY_SET_ID
    for rec in res["recommendations"]:
        assert SPEC_KEYS <= set(rec)
        assert rec["technology"].keys() >= {"name", "regulatory_context", "technology_evidence"}
        assert rec["geography"].keys() >= {"name", "fips", "burden", "burden_evidence_level", "vulnerability",
                                           "diagnostic_desert"}
        assert "candidate deployment opportunity" in rec["recommended_next_step"].lower()
        assert "not a validated diagnostic pathway" in rec["recommended_next_step"].lower()
        assert any("inherited" in u for u in rec["uncertainties"])
        assert rec["provenance"]["object_ids"][0].startswith("opportunity:")
        for m in rec["molecular_context"]:
            assert "condition-level molecular enrichment" in m["label"]
    text = p.read_text().lower()
    for b in BANNED:
        assert not re.search(b, text), b


@pytest.mark.data
def test_deployment_candidates_and_results_files():
    from measure_it.store import read_table, table_exists
    if not table_exists("deployment_candidates"):
        pytest.skip("demo not run")
    dc = read_table("deployment_candidates")
    assert len(dc) == 10 and set(PROVENANCE_COLUMNS) <= set(dc.columns)
    assert dc["recommendation_json"].map(lambda s: SPEC_KEYS <= set(json.loads(s))).all()
    # inherited (state) burden keeps source_geographic_resolution = state (CONVENTIONS 3.5)
    assert (dc.loc[dc["burden_inherited"], "source_geographic_resolution"] == "state").all()
    for s in dc["recommendation_json"]:
        rec = json.loads(s)
        ids = set(rec["provenance"]["object_ids"])
        rev = rec["research_evidence"]
        for k in ("condition_trials", "technology_experience_trials"):
            assert set(rev[k]["object_ids"]) <= ids, k
        if rec["geography"]["burden"]["inherited"]:
            assert "inherited" in (rec["geography"]["burden_evidence_level_note"] or "")
    for f in ("test5_contrasts", "test5_movements", "test6_scoring_label_shuffle", "test6_scoring_clinic_shuffle",
              "test6_scoring_phenotype_label_nulls", "test7_summary", "test7_unstable_regions",
              "demo_top10_regions"):
        assert (TABLES / f"{f}.csv").exists(), f
    md = RESULTS / "SCORING_RESULTS.md"
    if md.exists():
        text = md.read_text().lower()
        for b in BANNED:
            assert not re.search(b, text), b


# ------------------------------------------------------------------------------------------------------------------
# burden completeness rule (plan deviation 7, 2026-09-24) and sensitivity S6 specialist_only (deviation 8)
# ------------------------------------------------------------------------------------------------------------------

def test_set_burden_member_mean_rules():
    B = np.array([[[0.9, 0.5, 0.2]], [[np.nan, 0.7, 0.4]], [[np.nan, np.nan, np.nan]]])  # members x draws x n
    defined = [True, True, False]   # the third member has no defined burden (level D everywhere)
    new = O.set_burden_member_mean(B, defined, "complete")[0]
    assert np.isnan(new[0])                                   # a defined member is missing: incomplete, no burden
    assert new[1] == pytest.approx(0.6) and new[2] == pytest.approx(0.3)
    old = O.set_burden_member_mean(B, defined, "available")[0]
    assert old[0] == pytest.approx(0.9)                       # the old rule rested on the other member alone
    assert np.isnan(O.set_burden_member_mean(B, [False, False, False], "complete")).all()
    with pytest.raises(ValueError):
        O.set_burden_member_mean(B, defined, "nope")


def test_permutation_among_keeps_rows_outside_the_mask():
    mask = np.array([True, False, True, True, False, True])
    p = C.permutation_among(mask, None, np.random.default_rng(3))
    assert sorted(p) == list(range(6)) and p[1] == 1 and p[4] == 4
    assert set(p[mask]) == set(np.flatnonzero(mask))
    full = np.ones(6, bool)
    assert (C.permutation_among(full, None, np.random.default_rng(5)) ==
            C.permutation(6, None, np.random.default_rng(5))).all()


@pytest.mark.data
def test_incomplete_burden_is_not_ranked_and_listed(opp):
    ranks = [c for c in opp.columns if c.startswith("rank_") and not c.startswith("rank_eligible")]
    comps = [c for c in opp.columns if c.startswith("composite")]
    inc = opp[opp["burden_incomplete"].astype(bool)]
    assert len(inc) and inc[ranks + comps + ["burden_pct"]].isna().all().all()
    assert (~inc["rank_eligible"]).all() and inc["burden_incomplete_reason"].notna().all()
    assert (inc["burden_evidence_level"] == "D").all() and not inc["burden_excluded_level_D"].any()
    assert set(inc["condition_id"]) == {"me_cfs", O.QUERY_SET_ID, "autonomic_activity_invisible_illness"}
    assert (inc["geo_level"] == "county").all()
    # exactly the counties where ME/CFS (a member with a defined burden) has no value
    me = opp[(opp["condition_id"] == "me_cfs") & (opp["geo_level"] == "county")]
    miss = set(me.loc[me["burden_value"].isna(), "geo_id"])
    for cid in (O.QUERY_SET_ID, "autonomic_activity_invisible_illness", "me_cfs"):
        c = inc[inc["condition_id"] == cid]
        assert c.groupby("measurement_id")["geo_id"].apply(set).map(lambda x: x == miss).all(), cid
        assert c["burden_members_missing"].eq("me_cfs").all()
    assert "09120" in miss                                    # Greater Bridgeport Planning Region CT
    # no ranked row of a set rests on part of its defined members
    ranked = opp[opp["rank_equal"].notna()]
    assert ranked["burden_members_missing"].fillna("").eq("").all()
    s = ranked[ranked["condition_id"] == O.QUERY_SET_ID]
    assert s["burden_members_contributing"].map(lambda x: set(str(x).split("|")) == {"long_covid", "me_cfs"}).all()
    # level-D-only conditions are never incomplete (no defined burden)
    assert not opp[opp["condition_id"].isin(["pots", "dysautonomia"])]["burden_incomplete"].any()
    lst = pd.read_csv(TABLES / "scoring_incomplete_burden.csv", dtype={"geo_id": str})
    assert len(lst) == len(inc) and lst["burden_incomplete_reason"].notna().all()


@pytest.mark.data
def test_burden_rule_before_after_tables():
    b = pd.read_csv(TABLES / "scoring_set_burden_rule_before_after.csv")
    p = b[(b["condition_id"] == O.QUERY_SET_ID) & (b["measurement_id"] == O.PRIMARY["measurement"]) &
          (b["geo_level"] == "county")].iloc[0]
    assert p["n_ranked_before"] - p["n_ranked_after"] == p["n_incomplete_population_eligible"] > 0
    assert 0 <= p["top10_overlap_before_after"] <= 10 and -1 <= p["kendall_tau_b_common"] <= 1
    # combinations without an incomplete region are unchanged by the rule
    z = b[b["n_incomplete"] == 0]
    assert (z["top25_overlap_before_after"] == 25).all() and (z["top10_overlap_before_after"] == 10).all()
    assert (z["n_ranked_before"] == z["n_ranked_after"]).all()
    t = pd.read_csv(TABLES / "scoring_set_burden_rule_top10.csv")
    left = t[t["status"].str.startswith("left: incomplete")]
    assert left["burden_incomplete"].all() and left["rank_after"].isna().all()


@pytest.mark.data
def test_specialist_provider_term_matches_geography_module():
    from measure_it.store import read_table
    f = read_table("geo_condition_features", columns=["geo_id", "geo_level", "condition_id", "relevant_specialists_n",
                                                      "diagnostic_desert_sens_specialists", "in_analysis_universe"])
    f = f[f["in_analysis_universe"]]
    for level in O.LEVELS:
        geo = O.spatial().levels[level].geo
        for cid in ("long_covid", "me_cfs"):
            g = f[(f["condition_id"] == cid) & (f["geo_level"] == level)].set_index("geo_id").reindex(geo["geo_id"])
            n, _ = O.relevant_provider_counts(cid, level, exclude_groups=("primary_care",))
            assert np.array_equal(n, g["relevant_specialists_n"].to_numpy(float)), (cid, level)
            m = O.member_data(cid, level, geo, provider_variant="specialists_only")
            assert np.allclose(m.desert, g["diagnostic_desert_sens_specialists"].to_numpy(float), equal_nan=True)
    with pytest.raises(ValueError):
        O.make_combo("long_covid", O.PRIMARY["measurement"], "county", provider_variant="nope")


@pytest.mark.data
def test_specialist_only_sensitivity_tables():
    s = pd.read_csv(TABLES / "test7_sensitivity_specialist_only.csv")
    p = s[(s["condition_id"] == O.QUERY_SET_ID) & (s["measurement_id"] == O.PRIMARY["measurement"]) &
          (s["geo_level"] == "county") & (s["variant"] == "specialist_only")]
    assert len(p) == 1 and bool(p["declared"].iloc[0])
    assert set(s["variant"]) == {"specialist_only", "desert_provider_term_only", "clinic_capacity_only"}
    assert s["kendall_tau_b"].between(-1, 1).all() and (s["top10_overlap"] <= 10).all()
    assert (s["top25_overlap"] <= 25).all()
    assert len(s[(s["geo_level"] == "county") & (s["variant"] == "specialist_only")]) == \
        len(O.all_condition_ids()) * len(O.MEASUREMENTS)
    qa = pd.read_csv(TABLES / "test7_sensitivity_specialist_only_qa.csv")
    chk = qa[qa["check"].str.contains("this module")]
    assert len(chk) and (chk["n_identical"] == chk["n_geographies"]).all()
    sh = qa["share_primary_care"].dropna()
    assert len(sh) and sh.between(0.8, 1.0).all()


@pytest.mark.data
def test_fips_renames_are_bridged_not_incomplete():
    """Oglala Lakota SD (46102, renamed from Shannon County 46113 in 2015) and Kusilvak AK (02158, from Wade Hampton
    02270) are published by CMS under the legacy FIPS; the CMS ingestion bridges them, so they carry a ME/CFS value,
    are not incomplete, and no legacy code is left in the burden table."""
    from measure_it.scoring import report as R
    r = {x["legacy_fips"]: x for x in R.bridged_fips_renames()}
    assert set(r) == {"46113", "02270"}
    for x in r.values():
        assert x["n_cms_rows_bridged"] > 0 and x["has_primary_value"] and not x["legacy_code_left_in_burden"]
    d = read_table("deployment_opportunities", columns=["geo_id", "condition_id", "measurement_id", "geo_level",
                                                        "burden_incomplete", "rank_equal"])
    p = d[(d["condition_id"] == "long_covid_or_me_cfs") & (d["measurement_id"] == "wearable_autonomic_activity_monitoring")
          & (d["geo_level"] == "county")].set_index("geo_id")
    assert not p.loc["46102", "burden_incomplete"] and p.loc["46102", "rank_equal"] == p.loc["46102", "rank_equal"]
    assert not p.loc["02158", "burden_incomplete"]
    md = (RESULTS / "SCORING_RESULTS.md").read_text()
    assert "1:1 FIPS renames are bridged" in md and "fips_as_published = 46113" in md
    assert "No published value for this geography in the primary measure (suppressed" not in md
    # the research-readiness sentence must not attribute a floor-level percentile to technology-experience trials
    assert "yet their research_readiness percentile is" not in md


@pytest.mark.data
def test_recommendation_lists_incomplete_regions():
    p = RESULTS / "example_deployment_recommendation.json"
    if not p.exists():
        pytest.skip("demo not run")
    res = json.loads(p.read_text())
    ex = res["excluded_incomplete_burden"]
    fips = {e["fips"] for e in ex}
    assert "09120" in fips and res["excluded_incomplete_burden_note"]
    assert not fips & {r["geography"]["fips"] for r in res["recommendations"]}
    for e in ex:
        assert e["rank"] == UNKNOWN and e["burden"] == UNKNOWN and "me_cfs" in e["reason"]
        assert not e["small_population"]
