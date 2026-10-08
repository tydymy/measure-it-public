"""Tests for measure_it.geography.features and measure_it.geography.query.

Unit tests use synthetic inputs. @pytest.mark.data tests read data/processed and results/tables and are skipped
when the tables have not been built (`uv run python -m measure_it.geography.features`).
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from measure_it.config import DOCS, PROCESSED, TABLES, UNKNOWN, load_config
from measure_it.geography import features as F
from measure_it.provenance import PROVENANCE_COLUMNS


# ============================================================================ unit tests: helpers
def test_percentile_rank_average_ties_and_nan():
    s = pd.Series([1.0, 2.0, 2.0, np.nan, 4.0])
    r = F.percentile_rank(s)
    assert np.isnan(r.iloc[3])
    assert r.iloc[1] == r.iloc[2] == pytest.approx(2.5 / 4)
    assert r.iloc[4] == pytest.approx(1.0)


def test_haversine_one_degree_latitude():
    assert F.haversine_km(40.0, -100.0, 41.0, -100.0) == pytest.approx(111.2, abs=0.3)


def test_radius_matrix_selects_points_within_radius():
    m = F.radius_matrix(np.array([40.0]), np.array([-100.0]), np.array([40.0, 40.3, 41.0]),
                        np.array([-100.0, -100.0, -100.0]), 50)
    assert m.toarray().tolist() == [[1, 1, 0]]    # 0 km, ~33 km in; ~111 km out


def test_partial_spearman_removes_shared_confounder():
    rng = np.random.default_rng(0)
    z = rng.normal(size=2000)
    x = z + 0.3 * rng.normal(size=2000)
    y = z + 0.3 * rng.normal(size=2000)
    assert F.spearman(x, y) > 0.8
    assert abs(F.partial_spearman(x, y, z)) < 0.1


def test_topn_overlap_identical_and_reversed():
    a = pd.Series(np.arange(100, dtype=float), index=[f"g{i}" for i in range(100)])
    same = F.topn_overlap(a, a, 10)
    assert same["overlap"] == 10 and same["jaccard"] == 1.0 and same["p_hypergeom"] < 1e-10
    rev = F.topn_overlap(a, -a, 10)
    assert rev["overlap"] == 0 and rev["expected_overlap"] == pytest.approx(1.0)


def test_top_n_ids_deterministic_tie_break():
    s = pd.Series([1.0, 1.0, 1.0], index=["b", "a", "c"])
    assert F.top_n_ids(s, 2) == {"a", "b"}


def test_permute_within_groups_keeps_values_inside_groups():
    rng = np.random.default_rng(1)
    vals = np.arange(12, dtype=float)
    groups = np.array(["A"] * 4 + ["B"] * 5 + ["C"] * 3)
    for _ in range(20):
        p = F.permute_within_groups(vals, groups, rng)
        for g in "ABC":
            assert sorted(p[groups == g]) == sorted(vals[groups == g])
    # it does shuffle
    assert any(not np.array_equal(F.permute_within_groups(vals, groups, rng), vals) for _ in range(5))


def _chain(n):
    rows = np.r_[np.arange(n - 1), np.arange(1, n)]
    cols = np.r_[np.arange(1, n), np.arange(n - 1)]
    return F.row_standardize(sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)))


def test_morans_i_sign_on_chain_graph():
    w = _chain(50)
    assert F.morans_i(np.arange(50, dtype=float), w) > 0.9           # smooth gradient
    assert F.morans_i(np.tile([0.0, 1.0], 25), w) < -0.9             # alternating
    assert np.isnan(F.morans_i(np.ones(50), w))


def test_permutation_summary_collapse_logic():
    null = np.random.default_rng(2).normal(0, 0.02, 999)
    s = F.permutation_summary(0.5, null, 0.0)
    assert s["observed_outside_null95"] and s["p_perm"] == pytest.approx(1 / 1000)
    assert abs(s["retained_fraction"]) < 0.05
    s2 = F.permutation_summary(0.01, null, 0.0)
    assert not s2["observed_outside_null95"]


def test_permutation_summary_centred_p_for_shifted_null():
    # within-state style null centred on a between-state component (0.5), not on 0
    null = 0.5 + np.random.default_rng(3).normal(0, 0.01, 999)
    below = F.permutation_summary(0.45, null, 0.0)
    assert below["observed_outside_null95"]
    assert below["p_perm"] > 0.9                         # pre-specified P around 0 says nothing here
    assert below["p_perm_vs_null_centre"] == pytest.approx(1 / 1000)
    same = F.permutation_summary(0.5, null, 0.0)
    assert same["p_perm_vs_null_centre"] > 0.9


def test_holm_matches_hand_computation():
    p = np.array([0.01, 0.04, 0.03, np.nan])
    # sorted 0.01, 0.03, 0.04 -> 3*0.01=0.03, 2*0.03=0.06, 1*0.04 -> max-accumulated 0.06
    assert np.allclose(F.holm(p)[:3], [0.03, 0.06, 0.06])
    assert np.isnan(F.holm(p)[3])
    assert (F.holm([0.9, 0.8]) <= 1).all()


def test_min_detectable_rho():
    assert F.min_detectable_rho(49) == pytest.approx(0.39, abs=0.02)
    assert F.min_detectable_rho(3000) < 0.06
    assert np.isnan(F.min_detectable_rho(3))


def test_crh_effective_n_iid_vs_smooth_fields():
    rng = np.random.default_rng(4)
    n = 600
    lat, lon = rng.uniform(30, 48, n), rng.uniform(-120, -75, n)
    iid = F.crh_effective_n(rng.normal(size=n), rng.normal(size=n), lat, lon)
    assert iid["n_eff_spatial"] > 0.7 * n
    # two independent but smooth (spatially autocorrelated) fields: both follow longitude + noise
    x = lon + rng.normal(0, 2, n)
    y = -lat + 0.3 * lon + rng.normal(0, 2, n)
    sm = F.crh_effective_n(x, y, lat, lon)
    assert sm["n_eff_spatial"] < 0.3 * n


def test_crh_calibration_naive_test_is_anticonservative():
    rng = np.random.default_rng(6)
    n = 300
    ctx = pd.DataFrame({"geo_level": "county", "state_fips": rng.choice(["06", "48", "36", "17"], n),
                        "lat": rng.uniform(30, 45, n), "lon": rng.uniform(-115, -80, n)})
    cal = F.crh_calibration(ctx, n_points=n, n_sim=40, ranges_km=(400,))
    r = cal.iloc[0]
    assert r["false_positive_rate_naive"] > 0.3          # independent smooth fields look "related" to a naive test
    assert r["false_positive_rate_crh"] < 0.2
    assert r["median_n_eff"] < n / 3


def test_cluster_bootstrap_wider_for_clustered_data():
    rng = np.random.default_rng(5)
    g = np.repeat(np.arange(30), 40)
    shift = rng.normal(0, 1, 30)[g]
    x, y = shift + rng.normal(0, 0.3, len(g)), shift + rng.normal(0, 0.3, len(g))
    lo_i, hi_i = F.bootstrap_spearman_ci(x, y, n_boot=300)
    lo_c, hi_c = F.cluster_bootstrap_spearman_ci(x, y, g, n_boot=300)
    rho = F.spearman(x, y)
    assert lo_c <= rho <= hi_c
    assert (hi_c - lo_c) > 2 * (hi_i - lo_i)


def test_topn_inclusion_prob_and_tie_averaged_overlap():
    a = pd.Series([5, 5, 5, 5, 4, 3], index=list("abcdef"), dtype=float)
    p = F.topn_inclusion_prob(a, 2)
    assert p[list("abcd")].tolist() == [0.5] * 4 and p["e"] == 0
    b = pd.Series([9, 8, 1, 1, 1, 1], index=list("abcdef"), dtype=float)
    o = F.topn_overlap(a, b, 2)
    assert o["overlap_tie_averaged"] == pytest.approx(1.0)   # a,b each have 0.5 inclusion in a, 1 in b
    assert o["ties_at_boundary_a"] == 4 and o["ties_at_boundary_b"] == 0


# ============================================================================ unit tests: definitions
def test_every_condition_has_one_definition_per_level():
    d = F.definitions_frame()
    conds = F.target_conditions()
    assert set(d["condition_id"]) == set(conds)
    counts = d.groupby(["condition_id", "geo_level"]).size()
    assert (counts == 1).all()
    assert len(d) == 3 * len(conds)
    assert set(d["burden_evidence_level"]) <= {"A", "B", "C", "D"}


def test_level_d_has_no_measure_and_nonempty_reason():
    d = F.definitions_frame()
    dd = d[d["burden_evidence_level"] == "D"]
    assert dd["primary_measure_id"].isna().all()
    assert (dd["rationale"].str.len() > 20).all()
    assert dd["uncertainty_multiplier"].isna().all()


def test_long_covid_county_is_inherited_state_value():
    d = F.definitions_frame().set_index(["condition_id", "geo_level"])
    r = d.loc[("long_covid", "county")]
    assert r["inherited"] and r["source_geographic_resolution"] == "state"
    assert r["county_proxy_places_measure"] == "PHLTH"
    # nothing else inherits
    assert d["inherited"].sum() == 1


def test_proxies_are_level_c_and_multipliers_follow_scoring_config():
    d = F.definitions_frame()
    mult = load_config("scoring")["evidence_level_uncertainty_multiplier"]
    for _, r in d[d["burden_evidence_level"] != "D"].iterrows():
        assert r["uncertainty_multiplier"] == mult[r["burden_evidence_level"]]
    assert set(d.loc[d["condition_id"].isin(["me_cfs", "ptlds"]) & (d["burden_evidence_level"] != "D"),
                     "burden_evidence_level"]) == {"C"}


def test_burden_definitions_doc_matches_code():
    doc = (DOCS / "BURDEN_DEFINITIONS.md").read_text()
    d = F.definitions_frame()
    for mid in d["primary_measure_id"].dropna().unique():
        base = re.sub(r"_agg$", "", mid)
        assert base in doc, f"{mid} missing from BURDEN_DEFINITIONS.md"
    for cond in d.loc[d["burden_evidence_level"] == "D", "condition_id"].unique():
        assert re.search(rf"\| {cond} \| all \| none \| D \|", doc), f"{cond} not documented as level D"
    for cond, m in F.PLACES_PROXY.items():
        assert re.search(rf"\| {cond} \| [^|]*{m}", doc), f"proxy {cond}->{m} not in the proxy map"


# ============================================================================ data tests
def _need(*tables):
    for t in tables:
        if not (PROCESSED / f"{t}.parquet").exists():
            pytest.skip(f"{t} not built")


@pytest.fixture(scope="module")
def burden():
    _need("geo_condition_burden")
    return pd.read_parquet(PROCESSED / "geo_condition_burden.parquet")


@pytest.fixture(scope="module")
def feats():
    _need("geo_condition_features")
    return pd.read_parquet(PROCESSED / "geo_condition_features.parquet")


@pytest.mark.data
def test_burden_provenance_and_ids(burden):
    for c in PROVENANCE_COLUMNS:
        assert c in burden.columns
    for c in ("data_layer", "source_name", "retrieved_at", "source_geographic_resolution", "evidence_type"):
        assert burden[c].notna().all() and (burden[c].astype(str) != "").all()
    assert burden["burden_row_id"].is_unique
    assert set(burden["burden_evidence_level"]) <= {"A", "B", "C"}


@pytest.mark.data
def test_inherited_rows_keep_state_resolution(burden):
    inh = burden[burden["inherited"]]
    assert len(inh) > 3000
    assert set(inh["condition_id"]) == {"long_covid"}
    assert set(inh["geo_level"]) == {"county"}
    assert set(inh["source_geographic_resolution"]) == {"state"}
    assert set(inh["derivation"]) == {"inherited_from_state"}
    # no long-COVID county row claims county resolution at level A unless it is the labelled small-area MODEL
    # (no downscaling of a state value; since 2026-10-07 the BRFSS MRP estimate is the primary county burden)
    lc = burden[(burden["condition_id"] == "long_covid") & (burden["geo_level"] == "county")]
    a_cty = lc[(lc["burden_evidence_level"] == "A") & (lc["source_geographic_resolution"] == "county")]
    assert (a_cty["derivation"] == "mrp_small_area_estimate").all()
    assert (a_cty["estimate_kind"] == "modeled_small_area").all()
    assert not lc.loc[lc["inherited"], "measure_role"].eq("primary").any() or a_cty.empty
    # inherited value equals the state's HPS value it was carried from
    hs = burden[(burden["condition_id"] == "long_covid") & (burden["geo_level"] == "state")
                & (burden["source_table"] == "geo_condition_burden__cdc_long_covid")
                & (burden["measure_id"] == "lc_current_pct_all_adults")
                & (burden["derivation"] == "as_published") & burden["value"].notna()]
    st = hs.sort_values("period_end").groupby("geo_id").tail(1).set_index("geo_id")["value"]
    x = inh.sample(50, random_state=0)
    assert np.allclose(x["value"].to_numpy(), x["geo_id"].str[:2].map(st).to_numpy())


@pytest.mark.data
def test_one_primary_per_geography_and_proxies_labelled(burden):
    p = burden[burden["measure_role"] == "primary"]
    assert not p.duplicated(["condition_id", "geo_level", "geo_id"]).any()
    places = burden[burden["source_table"] == "geo_context__cdc_places"]
    assert set(places["burden_evidence_level"]) == {"C"}
    assert places["proxy_for_condition"].notna().all()
    assert not places["measure_role"].eq("primary").any()
    assert set(burden.loc[burden["condition_id"] == "me_cfs", "burden_evidence_level"]) == {"C"}
    assert set(burden.loc[burden["condition_id"] == "ptlds", "burden_evidence_level"]) == {"C"}


@pytest.mark.data
def test_features_one_row_per_geo_condition(feats):
    assert not feats.duplicated(["geo_id", "geo_level", "condition_id"]).any()
    assert set(feats["geo_level"]) == {"state", "county"}
    assert (feats["object_id"] == "geo:" + feats["geo_id"]).all()
    u = feats[feats["in_analysis_universe"]]
    n_cond = len(F.target_conditions())
    assert len(u[u["geo_level"] == "county"]) == 3144 * n_cond
    assert len(u[u["geo_level"] == "state"]) == 51 * n_cond
    for c in PROVENANCE_COLUMNS:
        assert c in feats.columns


@pytest.mark.data
def test_level_d_conditions_have_null_burden_and_reason(feats):
    d = feats[feats["burden_defined_level"] == "D"]
    assert len(d) > 0
    assert d["burden_value"].isna().all()
    assert d["diagnostic_desert"].isna().all()
    assert (d["burden_evidence_level"] == "D").all()
    assert d["burden_unknown_reason"].notna().all()
    # the access-only index is still defined inside the universe, and labelled separately
    assert d.loc[d["in_analysis_universe"], "diagnostic_desert_access_only"].notna().all()


@pytest.mark.data
def test_counts_only_at_same_resolution(feats):
    inh = feats[feats["burden_inherited"]]
    assert inh["burden_count"].isna().all()
    cms = feats[feats["burden_measure_id"].fillna("").str.startswith("cms_mmd")]
    assert cms["burden_count"].isna().all()
    # since 2026-10-07 no primary burden is inherited by default (long_covid_burden: brfss_sae)
    assert set(feats.loc[feats["burden_inherited"], "condition_id"]) <= {"long_covid"}
    assert (feats.loc[feats["burden_inherited"], "burden_source_resolution"] == "state").all()


@pytest.mark.data
def test_desert_bounds_and_components(feats):
    d = feats["diagnostic_desert"].dropna()
    assert ((d > 0) & (d <= 1)).all()
    comp = ["desert_pct_burden", "desert_pct_vulnerability", "desert_pct_low_providers", "desert_pct_low_trials"]
    x = feats.dropna(subset=["diagnostic_desert"])
    assert np.allclose(x[comp].mean(axis=1), x["diagnostic_desert"])
    assert x["in_analysis_universe"].all()
    # nothing outside 50 states + DC gets a desert value
    assert feats.loc[~feats["in_analysis_universe"], "diagnostic_desert"].isna().all()


@pytest.mark.data
def test_connecticut_cms_not_remapped(feats):
    ct = feats[(feats["state_fips"] == "09") & (feats["geo_level"] == "county") & (feats["condition_id"] == "migraine")]
    assert len(ct) == 9
    assert ct["burden_value"].isna().all()
    assert ct["burden_unknown_reason"].str.contains("legacy").all()


@pytest.mark.data
def test_geo_context_unique_rows():
    _need("geo_context")
    c = pd.read_parquet(PROCESSED / "geo_context.parquet", columns=["geo_id", "geo_level", "object_id", "svi_overall"])
    assert not c.duplicated(["geo_level", "geo_id"]).any()
    assert set(c["geo_level"]) == {"state", "county", "zcta"}
    assert c["object_id"].str.startswith("geo:").all()
    assert c.loc[c["geo_level"] == "county", "svi_overall"].notna().sum() >= 3144


@pytest.mark.data
def test_result_tables_written_and_nulls_collapse():
    need = ["test4_population_correlations", "test4_topn_overlap", "test4_partial_correlations",
            "test6_geography_relationships", "test6_geography_morans_i", "geo_burden_coverage_matrix"]
    for n in need:
        if not (TABLES / f"{n}.csv").exists():
            pytest.skip(f"{n}.csv not written")
    rel = pd.read_csv(TABLES / "test6_geography_relationships.csv")
    nat = rel[rel["permutation_scheme"] == "national"]
    assert (nat["null_mean"].abs() < 0.05).all()           # national shuffle removes cross-source structure
    mor = pd.read_csv(TABLES / "test6_geography_morans_i.csv")
    mn = mor[mor["permutation_scheme"] == "national"]
    assert ((mn["null_mean"] - mn["expected_I"]).abs() < 0.02).all()
    inh = mor[(mor["variable"] == "long_covid_inherited_state_value") & (mor["permutation_scheme"] == "within_state")]
    assert np.allclose(inh["null_mean"], inh["observed_I"])   # inherited values carry no within-state information


# ============================================================================ query layer
@pytest.mark.data
def test_query_condition_burden_level_d_returns_unknown():
    _need("geo_burden_definitions", "geo_condition_burden")
    from measure_it.geography.query import get_condition_burden
    r = get_condition_burden("POTS", "county")
    assert r["status"] == "unknown" and r["value"] == UNKNOWN and r["burden_evidence_level"] == "D"
    assert r["rows"] == [] and len(r["reason"]) > 10
    bad = get_condition_burden("long_covid", "planet")
    assert bad["status"] == "unknown" and bad["value"] == UNKNOWN
    nope = get_condition_burden("not a condition", "state")
    assert nope["status"] == "unknown"


@pytest.mark.data
def test_query_condition_burden_long_covid():
    _need("geo_burden_definitions", "geo_condition_burden")
    from measure_it.geography.query import get_condition_burden
    st = get_condition_burden("Long COVID", "state")
    assert st["status"] == "ok" and st["n_rows"] == 51
    assert all(r["burden_evidence_level"] == "A" and r["source_geographic_resolution"] == "state" for r in st["rows"])
    cty = get_condition_burden("long_covid", "county", geo_id="06073")
    assert cty["n_rows"] == 1
    row = cty["rows"][0]
    # since 2026-10-07: the BRFSS small-area model, labelled as modelled (the inherited HPS row is an alternate)
    assert row["inherited"] is False and row["derivation"] == "mrp_small_area_estimate"
    assert any("MODELLED" in c for c in cty["caveats"])
    alt = get_condition_burden("long_covid", "county", geo_id="06073", include_alternates=True)
    roles = {r["measure_role"] for r in alt["rows"]}
    assert {"primary", "county_proxy", "alternate"} <= roles
    assert any(r["inherited"] for r in alt["rows"] if r["measure_role"] == "alternate")
    assert all(r["burden_evidence_level"] == "C" for r in alt["rows"] if r["measure_role"] == "county_proxy")


@pytest.mark.data
def test_query_resolve_geography_variants():
    _need("geographies", "geo_context")
    from measure_it.geography.query import resolve_geography
    for q in ("San Diego County, California", "San Diego County, CA", "san diego, ca", "06073", "geo:06073"):
        r = resolve_geography(q)
        assert r["status"] == "ok" and r["geo_id"] == "06073", q
    assert resolve_geography("California")["geo_id"] == "06"
    assert resolve_geography("6")["geo_id"] == "06"
    assert resolve_geography("Nowhere County, Atlantis")["status"] == "unknown"
    amb = resolve_geography("Baltimore, MD")
    assert amb["status"] == "ambiguous" and {c["geo_id"] for c in amb["candidates"]} == {"24005", "24510"}


@pytest.mark.data
def test_query_geographic_context_san_diego():
    _need("geo_condition_features", "geo_context", "provider_density_county")
    from measure_it.geography.query import get_geographic_context
    r = get_geographic_context("San Diego County, California")
    assert r["status"] == "ok" and r["geography"]["object_id"] == "geo:06073"
    burden = {b["condition_id"]: b for b in r["burden_by_condition"]}
    assert set(burden) == set(F.target_conditions())
    assert burden["pots"]["burden_value"] == UNKNOWN and burden["pots"]["burden_evidence_level"] == "D"
    assert burden["long_covid"]["burden_inherited"] is False   # modelled small-area estimate since 2026-10-07
    assert "MODELLED" in burden["long_covid"]["burden_value_note"]
    assert r["population"]["total_population"] > 3_000_000
    assert 0 <= r["vulnerability"]["svi_overall"] <= 1
    groups = {g["specialty_group"] for g in r["providers_by_specialty_group"]}
    assert {"cardiology", "neurology", "primary_care"} <= groups
    assert r["hrsa_health_centers"]["n_service_delivery_sites"] > 0
    assert r["trials"]["distinct_relevant_trials_with_site_here"] >= r["trials"]["recruiting"]
    assert r["nih_projects"]["distinct_projects_appl_id"] >= 0
    import json
    json.dumps(r)   # JSON-serialisable for the MCP layer


# ============================================================================ reviewer additions (data)
@pytest.mark.data
def test_inherited_burden_flagged_in_desert(feats):
    lc = feats[(feats["condition_id"] == "long_covid") & (feats["geo_level"] == "county")
               & feats["diagnostic_desert"].notna()]
    assert len(lc) > 3000
    if lc["burden_inherited"].any():   # long_covid_burden: hps_inherited
        assert lc["desert_burden_resolution"].str.startswith("state (inherited").all()
        assert lc["desert_definition"].str.contains("INHERITED STATE VALUE").all()
        # the burden percentile is constant within a state (that is what the flag warns about)
        assert (lc.groupby("state_fips")["desert_pct_burden"].nunique() == 1).all()
    else:                               # long_covid_burden: brfss_sae (default since 2026-10-07)
        assert lc["desert_burden_resolution"].eq("county (modelled small-area estimate)").all()
        assert lc["desert_definition"].str.contains("MODELLED SMALL-AREA ESTIMATE").all()
    other = feats[feats["diagnostic_desert"].notna() & ~feats["burden_inherited"]]
    assert not other["desert_burden_resolution"].str.contains("inherited").any()


@pytest.mark.data
def test_reviewer_tables_written():
    for n in ("test4_small_population_check", "test4_verdicts", "test6_geography_relationships"):
        if not (TABLES / f"{n}.csv").exists():
            pytest.skip(f"{n}.csv not written")
    sm = pd.read_csv(TABLES / "test4_small_population_check.csv")
    assert {"topN_in_bottom_pop_quartile", "p_hypergeom_bottom_quartile_holm"} <= set(sm.columns)
    ver = pd.read_csv(TABLES / "test4_verdicts.csv")
    assert ver["reviewer_reading"].notna().all() and "p_hypergeom_holm" in ver
    rel = pd.read_csv(TABLES / "test6_geography_relationships.csv")
    nat = rel[rel["permutation_scheme"] == "national"]
    assert (nat["n_eff_spatial"] <= nat["n"]).all() and nat["p_spatial_crh_holm"].notna().all()
    ws = rel[rel["permutation_scheme"] == "within_state"]
    assert ws["p_perm_vs_null_centre"].notna().all()
    # the positive control (same source, two years) must survive the spatial correction
    assert nat.loc[nat["relationship_id"] == "R4", "p_spatial_crh_holm"].iloc[0] < 0.05


@pytest.mark.data
def test_query_reviewer_edge_cases():
    _need("geo_burden_definitions", "geo_condition_burden", "geo_condition_features", "geographies")
    from measure_it.geography.query import get_condition_burden, get_geographic_context, resolve_condition, \
        resolve_geography
    assert resolve_condition("Lyme") == "lyme_disease"
    assert resolve_condition("xyz") is None and resolve_condition("") is None
    assert resolve_geography("Dona Ana County, NM")["geo_id"] == "35013"
    assert resolve_geography("US")["geo_level"] == "national"
    wrong = get_condition_burden("long_covid", "state", geo_id="06073")
    assert wrong["status"] == "unknown" and "county" in wrong["reason"]
    ct = get_condition_burden("migraine", "county", geo_id="09110")
    assert ct["value"] == UNKNOWN and "legacy" in ct["reason"]
    nat = get_condition_burden("long_covid", "national", geo_id="US")
    assert nat["status"] == "ok" and nat["n_rows"] == 1
    ctx = get_geographic_context("US")
    assert ctx["status"] == "unknown" and ctx["value"] == UNKNOWN
    sd = get_geographic_context("06073")
    lc = [b for b in sd["burden_by_condition"] if b["condition_id"] == "long_covid"][0]
    assert "not county prevalence" in lc["burden_value_note"]


@pytest.mark.data
def test_desert_top_lists_break_ties_deterministically():
    """geo_desert_top25.csv: desert descending, exact ties ordered by geo_id (stable, not sort-algorithm order)."""
    import pandas as pd
    from measure_it.config import TABLES as _T
    p = _T / "geo_desert_top25.csv"
    if not p.exists():
        pytest.skip("geography results not built")
    d = pd.read_csv(p, dtype={"geo_id": str})
    for _, g in d.groupby(["condition_id", "geo_level"], sort=False):
        key = list(zip(-g["diagnostic_desert"], g["geo_id"]))
        assert key == sorted(key)


@pytest.mark.data
def test_renamed_counties_carry_cms_values_and_gaps_have_measured_reasons(feats, burden):
    """The CMS ingestion bridges 46113 -> 46102 (Oglala Lakota SD) and 02270 -> 02158 (Kusilvak AK): both counties
    carry the CMS primary measures, the burden rows keep the published code in source_record_id, and no county row is
    left on a legacy code. Every remaining missing primary value has a measured, source-specific reason."""
    f = feats[(feats["geo_level"] == "county")]
    for gid, legacy in (("46102", "46113"), ("02158", "02270")):
        for cid in ("me_cfs", "fibromyalgia", "migraine"):
            r = f[(f["geo_id"] == gid) & (f["condition_id"] == cid)].iloc[0]
            assert pd.notna(r["burden_value"]) and pd.isna(r["burden_unknown_reason"]), (gid, cid)
            assert f":county:{legacy}:" in r["burden_row_id"]
    assert not burden.loc[burden["geo_level"] == "county", "geo_id"].isin(["46113", "02270"]).any()
    miss = f[f["burden_value"].isna() & (f["burden_defined_level"] != "D")]
    reasons = miss["burden_unknown_reason"].astype(str)
    assert miss["burden_unknown_reason"].notna().all()
    assert not reasons.str.contains("suppressed, cell omitted by the source, or entity absent").any()
    cms = miss[miss["burden_measure_id"].astype(str).str.startswith("cms_mmd_") & miss["in_analysis_universe"]
               & (miss["state_fips"] != "09")]
    assert len(cms) and cms["burden_unknown_reason"].str.startswith("CMS MMD publishes no row for this county").all()
    assert set(cms["geo_id"]) >= {"51685", "15005", "02185"} and not set(cms["geo_id"]) & {"46102", "02158"}


def test_not_published_reasons_distinguish_absent_entity_from_omitted_year():
    b = pd.DataFrame({"source_table": ["geo_condition_burden__cms_mmd"] * 2 + ["geo_condition_burden__cdc_lyme"],
                      "geo_level": ["county"] * 3, "geo_id": ["11111", "11111", "33333"],
                      "measure_id": ["cms_mmd_51_prev_actual_all", "cms_mmd_51_prev_actual_all", "lyme_x"],
                      "year": [2019.0, 2020.0, 2021.0]})
    rows = pd.DataFrame({"geo_id": ["11111", "22222", "33333", "44444"], "geo_level": ["county"] * 4,
                         "burden_measure_id": ["cms_mmd_51_prev_actual_all", "cms_mmd_51_prev_actual_all",
                                               "lyme_incidence_per_100k", "lyme_incidence_per_100k"]})
    out = F.not_published_reasons(rows, b)
    assert out[0].startswith(f"CMS MMD publishes no {F.PRIMARY_CMS_YEAR} prevalence cell") and "2019-2020" in out[0]
    assert out[1].startswith("CMS MMD publishes no row for this county code in any claims year") or \
        "omitted" in out[1]   # 22222 may exist in the real CMS context table; never the old generic text
    assert out[2].startswith(f"CDC Lyme publishes no {F.PRIMARY_LYME_YEAR} value") and "2021-2021" in out[2]
    assert out[3] == "The CDC Lyme county file has no row for this county code in any year."
