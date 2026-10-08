"""Device-defined subgroup burden by geography (measure_it.geography.subgroup), built against the contract with a
synthetic `subgroup_strata_fractions` table (HARMONIZATION_CONTRACT section 2)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from measure_it.geography import subgroup as G
from measure_it.store import read_table, table_exists


def synthetic_fractions(pid="byod_test:lc_capillaroscopy_abnormal", cond="long_covid") -> pd.DataFrame:
    rows = [dict(phenotype_id=pid, condition_id=cond, age_band="all", sex="all", n_cases=60, n_positive=24,
                 fraction=0.4, ci_low=0.28, ci_high=0.53, suppressed=False)]
    for band in ["20-29", "30-39", "40-49", "50-59", "60-69"]:
        for sex in ["female", "male"]:
            sup = band == "60-69" and sex == "male"
            rows.append(dict(phenotype_id=pid, condition_id=cond, age_band=band, sex=sex, n_cases=5 if sup else 12,
                             n_positive=None if sup else 6, fraction=None if sup else 0.5, ci_low=None, ci_high=None,
                             suppressed=sup))
    return pd.DataFrame(rows)


def test_parse_age_band():
    assert G.parse_age_band("40-49") == (40, 49)
    assert G.parse_age_band("90+") == (90, 200)
    assert G.parse_age_band("all") is None
    with pytest.raises(ValueError):
        G.parse_age_band("forties")


def test_stratum_key_fallbacks():
    avail = {("all", "all"): 1, ("40-49", "female"): 1, ("50-59", "all"): 1, ("all", "male"): 1}
    assert G.stratum_key(45, "female", avail) == (("40-49", "female"), "age_band_x_sex")
    assert G.stratum_key(52, "female", avail) == (("50-59", "all"), "age_band")
    assert G.stratum_key(45, "male", avail) == (("all", "male"), "sex")
    assert G.stratum_key(18, "female", avail) == (("all", "all"), "pooled")
    # ACS 85+ cell goes to the band containing its lower bound
    assert G.stratum_key(85, "female", {("all", "all"): 1, ("80-89", "female"): 1, ("90+", "female"): 1})[0] == \
        ("80-89", "female")


def test_fraction_draws_jeffreys_and_suppression():
    fr = synthetic_fractions()
    d = G.fraction_draws(fr, 20000, np.random.default_rng(0))
    assert ("60-69", "male") not in d                       # suppressed -> not usable
    assert d[("all", "all")].mean() == pytest.approx((24 + 0.5) / 61, abs=0.005)
    assert d[("40-49", "female")].mean() == pytest.approx(0.5, abs=0.01)
    with pytest.raises(ValueError):
        G.fraction_draws(fr[fr.age_band != "all"], 10, np.random.default_rng(0))


def test_fraction_draws_from_ci_when_counts_missing():
    fr = pd.DataFrame([dict(age_band="all", sex="all", n_cases=None, n_positive=None, fraction=0.3, ci_low=0.2,
                            ci_high=0.4, suppressed=False)])
    d = G.fraction_draws(fr, 40000, np.random.default_rng(1))[("all", "all")]
    assert d.mean() == pytest.approx(0.3, abs=0.005)
    assert np.quantile(d, 0.975) - np.quantile(d, 0.025) == pytest.approx(0.2, abs=0.02)


def test_subgroup_counts_hand_check():
    N = np.array([[100.0, 200.0]])
    prev = np.array([[[0.10], [0.20]]])
    frac = np.array([[0.5], [0.25]])
    assert G.subgroup_counts(N, prev, frac)[0, 0] == pytest.approx(100 * 0.1 * 0.5 + 200 * 0.2 * 0.25)


def test_rows_and_aggregate_consistency():
    rng = np.random.default_rng(2)
    base = rng.uniform(10, 20, (3, 50))
    sub = base * 0.3
    pop = np.array([1000.0, 2000.0, 500.0])
    keys, B, S, P = G._aggregate(np.array(["01001", "01003", "02001"]), base, sub, pop, np.array(["01", "01", "02"]))
    assert list(keys) == ["01", "02"] and P.tolist() == [3000.0, 500.0]
    assert np.allclose(S[0], sub[0] + sub[1])
    r = G._rows(keys, "state", B, S, P, {"burden_basis": "x"})
    assert (r["subgroup_ci_low"] <= r["subgroup_estimate"]).all()
    assert r["subgroup_rate_per_100k"].iloc[1] == pytest.approx(sub[2].mean() / 500 * 1e5)


@pytest.mark.data
def test_build_long_covid_and_lyme_with_synthetic_fractions():
    if not table_exists("acs_county_adult_poststrat_cells"):
        pytest.skip("run measure_it.geography.sae")
    fr = pd.concat([synthetic_fractions(), synthetic_fractions("byod_test:lyme_x", "lyme_disease"),
                    synthetic_fractions("byod_test:ptlds_x", "ptlds")], ignore_index=True)
    out, notes = G.build(fr)
    assert not notes
    lc = out[out.phenotype_id == "byod_test:lc_capillaroscopy_abnormal"]
    assert (lc["burden_basis"] == "brfss_2023_mrp_sae").all()
    cty, st = lc[lc.geo_level == "county"], lc[lc.geo_level == "state"]
    assert len(cty) == 3144 and len(st) == 51
    assert (cty["subgroup_estimate"] < cty["base_burden"]).all()
    assert st["subgroup_estimate"].sum() == pytest.approx(cty["subgroup_estimate"].sum(), rel=1e-6)
    # fractions near 0.4-0.5 of ~7% of 261M adults
    nat = lc[lc.geo_level == "national"].iloc[0]
    assert 5e6 < nat["subgroup_estimate"] < 12e6 and nat["subgroup_ci_low"] < nat["subgroup_estimate"]
    ly = out[out.phenotype_id == "byod_test:lyme_x"]
    nat_ly = ly[ly.geo_level == "national"].iloc[0]
    cases = read_table("geo_condition_burden__cdc_lyme")
    c23 = cases[(cases.year == 2023) & (cases.measure_id == "lyme_reported_cases") & (cases.geo_level == "county")
                & cases.in_canonical_geographies.astype(bool)]["value"].sum()
    assert nat_ly["base_burden"] == pytest.approx(c23, rel=0.01)
    pt = out[(out.phenotype_id == "byod_test:ptlds_x") & (out.geo_level == "national")].iloc[0]
    assert pt["base_burden"] == pytest.approx(c23 * 0.137, rel=0.05)
    assert pt["burden_basis"] == "derived_ptlds_from_lyme_cases_aucott2022"


@pytest.mark.data
def test_generic_fallback_never_counts_inherited_state_values():
    if not table_exists("geo_condition_features"):
        pytest.skip("geo_condition_features missing")
    fd = G.fraction_draws(synthetic_fractions(), 100, np.random.default_rng(0))
    t, _ = G.generic_rows("long_covid", fd)
    cty = t[t.geo_level == "county"]
    assert cty["subgroup_estimate"].isna().all()            # no county count from an inherited state value
    f = read_table("geo_condition_features").set_index(["geo_id", "geo_level", "condition_id"])
    inh = f.loc[[(g, "county", "long_covid") for g in cty["geo_id"]], "burden_inherited"].fillna(False).to_numpy()
    assert (cty.loc[inh.astype(bool), "burden_basis"] == "inherited_state_value_no_county_count").all()


@pytest.mark.data
def test_real_subgroup_table_if_present():
    if not table_exists(G.TABLE):
        pytest.skip("geo_subgroup_burden not built (needs subgroup_strata_fractions from stage B+C)")
    t = read_table(G.TABLE)
    assert {"phenotype_id", "geo_id", "geo_level", "subgroup_estimate", "subgroup_ci_low", "subgroup_ci_high",
            "burden_basis"} <= set(t.columns)
    assert (t["data_layer"] == "derived").all()
