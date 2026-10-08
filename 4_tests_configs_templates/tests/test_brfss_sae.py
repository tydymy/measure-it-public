"""BRFSS ingestion (measure_it.ingestion.brfss) and the county small-area estimate (measure_it.geography.sae)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

from measure_it.geography import sae as S
from measure_it.ingestion import brfss as B
from measure_it.store import read_table, table_exists


# ------------------------------------------------------------------------------------------------ pure: ingestion
def test_weighted_prop_equal_weights_one_stratum():
    y = np.array([1, 0, 0, 1, 0, 0, 0, 0, 1, 0], float)
    r = B.weighted_prop(y, np.ones(10), np.ones(10), np.arange(10))
    assert r["p"] == pytest.approx(0.3)
    # every record its own PSU in one stratum: var = n/(n-1) * sum((y-p)/n)^2 = p(1-p)/(n-1)
    assert r["se"] == pytest.approx(np.sqrt(0.3 * 0.7 / 9))
    assert r["n"] == 10 and r["cases"] == 3


def test_weighted_prop_weights_nan_and_domain():
    y = np.array([1, 0, np.nan, 1], float)
    w = np.array([2, 1, 5, 1], float)
    r = B.weighted_prop(y, w, np.zeros(4), np.arange(4))
    assert r["p"] == pytest.approx(3 / 4)
    dom = np.array([True, True, True, False])
    r2 = B.weighted_prop(y, w, np.zeros(4), np.arange(4), domain=dom)
    assert r2["p"] == pytest.approx(2 / 3) and r2["n"] == 2


def test_logit_ci_contains_point_and_is_asymmetric():
    lo, hi = B.logit_ci(0.05, 0.01)
    assert lo < 0.05 < hi and (hi - 0.05) > (0.05 - lo)
    assert np.isnan(B.logit_ci(0.0, 0.01)[0])


def _raw23(rows):
    base = {"_STATE": 1.0, "_STSTR": 1011.0, "_PSU": 1.0, "_LLCPWT": 1.0, "DISPCODE": 1100.0, "_AGEG5YR": 5.0,
            "_SEX": 2.0, "_IMPRACE": 1.0, "_EDUCAG": 4.0, "_METSTAT": 1.0, "_URBSTAT": 1.0}
    return pd.DataFrame([{**base, **r} for r in rows])


def test_derive_2023_outcome_coding():
    rows = [
        {"COVIDPO1": 1, "COVIDSM1": 1, "COVIDACT": 1},       # current LC, significant
        {"COVIDPO1": 1, "COVIDSM1": 1, "COVIDACT": 2},       # current LC, a little
        {"COVIDPO1": 1, "COVIDSM1": 1, "COVIDACT": 7},       # current LC, act unknown
        {"COVIDPO1": 1, "COVIDSM1": 2, "COVIDACT": np.nan},  # positive, no LC
        {"COVIDPO1": 2, "COVIDSM1": np.nan, "COVIDACT": np.nan},  # never positive
        {"COVIDPO1": 7, "COVIDSM1": np.nan, "COVIDACT": np.nan},  # don't know
        {"COVIDPO1": np.nan, "COVIDSM1": np.nan, "COVIDACT": np.nan, "DISPCODE": 1200},  # break-off
        {"COVIDPO1": 1, "COVIDSM1": 7, "COVIDACT": np.nan},  # LC unknown
    ]
    d = B.derive(_raw23(rows), 2023)
    cc = d["lc_current_pct_all_adults"].tolist()
    assert cc[:5] == [1, 1, 1, 0, 0] and all(np.isnan(cc[5:]))
    assert d["lc_current_pct_all_adults_mmwr_def"].tolist() == [1, 1, 1, 0, 0, 0, 0, 0]
    sig = d["lc_significant_limitation_pct_all_adults"].tolist()
    assert sig[:2] == [1, 0] and np.isnan(sig[2]) and sig[3:5] == [0, 0] and all(np.isnan(sig[5:]))
    among = d["lc_significant_limitation_pct_current_lc_mmwr_def"].tolist()
    assert among[:3] == [1, 0, 0] and all(np.isnan(among[3:]))
    assert d["partial_interview"].tolist()[6]


def test_derive_covariate_categories():
    d = B.derive(_raw23([{"COVIDPO1": 2, "COVIDSM1": np.nan, "COVIDACT": np.nan, "_AGEG5YR": a, "_IMPRACE": r,
                          "_METSTAT": m, "_URBSTAT": u, "_EDUCAG": e}
                         for a, r, m, u, e in [(1, 1, 1, 1, 1), (13, 5, 2, 1, 4), (14, 4, 2, 2, 9), (6, 2, np.nan, np.nan, 2)]]),
                 2023)
    assert d["age7"].tolist()[:2] == [0, 6] and pd.isna(d["age7"].iloc[2]) and d["age7"].iloc[3] == 3
    assert d["race4"].tolist() == [0, 2, 3, 1]
    assert d["urban3"].tolist()[:3] == [0, 1, 2] and pd.isna(d["urban3"].iloc[3])
    assert d["edu4"].iloc[0] == 0 and d["edu4"].iloc[1] == 3 and pd.isna(d["edu4"].iloc[2])


def test_parse_mmwr_table_lines():
    text = """ Jurisdiction   No.   % (95% CI)   No.  % (95% CI)
 National                          27,074           6.4 (6.3, 6.6)          5,722           19.8 (18.9, 20.8)
 Alabama                             382            9.1 (8.0, 10.3)           90            19.7 (15.2, 25.1)
 Kentucky                             —                                       —
 District of Columbia                115            3.8 (3.1, 4.7)            17             12.8 (7.8, 20.3)
"""
    t = B.parse_mmwr_table(text)
    assert t["jurisdiction"].tolist() == ["National", "Alabama", "Kentucky", "District of Columbia"]
    assert t.loc[0, "n_cases_current"] == 27074 and t.loc[1, "mmwr_current_hi"] == 10.3
    assert t.loc[3, "mmwr_sig_among_lc_pct"] == 12.8 and pd.isna(t.loc[2, "mmwr_current_pct"])


def test_codebook_and_layout_parsers():
    cb = ("<p>Label: Do you currently have Covid-19 symptoms? Section Name: Long-term COVID Effects Core Section Number:"
          " 16 Question Number: 2 Column: 258 Type of Variable: Num SAS Variable Name: COVIDSM1 Question Prologue: "
          "Question: Do you currently have symptoms lasting 3 months or longer? Value Value Label Frequency</p>")
    e = B.codebook_entries(cb)
    assert e["COVIDSM1"]["section"] == "Long-term COVID Effects Core"
    assert "3 months" in e["COVIDSM1"]["question"]
    lay = "<tr><td>257</td><td>COVIDPO1</td><td>1</td></tr><tr><td>1</td><td>_STATE</td></tr>"
    assert B.layout_variables(lay) == ["COVIDPO1", "_STATE"]


# ------------------------------------------------------------------------------------------------ pure: SAE
def test_age7_mapping_and_design_shape():
    assert [S.age7_of(x) for x in (18, 24, 25, 44, 45, 64, 65, 74, 75, 85)] == [0, 0, 1, 2, 3, 4, 5, 5, 6, 6]
    g = S.demo_cell_grid()
    X = S.demo_design(g["age7"], g["sex"], g["race4"], g["edu4"])
    assert X.shape == (S.N_DEMO_CELLS, len(S.DEMO_NAMES))
    assert (X[:, 0] == 1).all()


def test_kish_scale_equal_weights_is_one():
    w = pd.Series([2.0, 2.0, 2.0, 5.0, 5.0])
    g = pd.Series(["a", "a", "a", "b", "b"])
    f = S.kish_scale(w, g)
    assert (w * f).groupby(g).sum().tolist() == pytest.approx([3.0, 2.0])


def test_fit_glmm_recovers_simulated_parameters():
    rng = np.random.default_rng(1)
    S_, n_cells = 40, 200
    sigma, beta = 0.3, np.array([-2.5, 0.8])
    u = rng.normal(0, sigma, S_)
    g = np.repeat(np.arange(S_), n_cells)
    x = rng.integers(0, 2, S_ * n_cells).astype(float)
    W = np.full(S_ * n_cells, 50.0)
    p = expit(beta[0] + beta[1] * x + u[g])
    Y = rng.binomial(50, p).astype(float)
    X = np.column_stack([np.ones_like(x), x])
    f = S.fit_glmm(X, Y, W, g, [str(i) for i in range(S_)], ["b0", "b1"])
    assert f.converged
    assert f.beta[1] == pytest.approx(0.8, abs=0.06)
    assert np.sqrt(f.sigma2) == pytest.approx(sigma, abs=0.1)
    assert np.corrcoef(f.u, u)[0, 1] > 0.9


def test_poststratification_hand_check():
    """One county, two populated cells: the county rate is the population-weighted mean of the cell rates."""
    names = S.DEMO_NAMES + ["urban_micropolitan", "urban_noncore"]
    beta = np.zeros(len(names))
    beta[0] = -3.0                                    # intercept
    beta[names.index("female")] = 1.0
    beta[names.index("urban_noncore")] = 0.5
    fit = S.Fit(beta=beta, u=np.array([0.2]), sigma2=0.04, cov=np.zeros((len(beta) + 1,) * 2), names=names,
                groups=["01"], loglik_laplace=0.0, n_cells=0, converged=True)
    N = np.zeros((1, S.N_AGE, S.N_SEX, S.N_RACE, S.N_EDU))
    N[0, 0, 0, 0, 0] = 300.0      # 18-24 male white NH < HS
    N[0, 0, 1, 0, 0] = 100.0      # 18-24 female white NH < HS
    counties = pd.DataFrame({"state_fips": ["01"], "urban3": [2], "adults": [400.0], "urbanicity": ["noncore"]},
                            index=pd.Index(["01001"], name="geo_id"))
    geo = S.Geo(counties, N, pd.DataFrame(index=counties.index), pd.DataFrame())
    bd, ud = beta[:, None].repeat(3, 1), np.array([[0.2, 0.2, 0.2]])
    dr, cells = S.predict_draws(fit, bd, ud, geo, False, True, np.random.default_rng(0), want_cells=True)
    pm, pf = expit(-3.0 + 0.5 + 0.2), expit(-3.0 + 1.0 + 0.5 + 0.2)
    assert dr[0, 0] == pytest.approx((300 * pm + 100 * pf) / 400, rel=1e-5)
    assert cells[0, 0, 0] == pytest.approx(pm, rel=1e-5) and cells[0, 1, 0] == pytest.approx(pf, rel=1e-5)


def test_absent_state_effect_is_drawn_from_prior():
    names = S.DEMO_NAMES + ["urban_micropolitan", "urban_noncore"]
    beta = np.zeros(len(names))
    beta[0] = -2.0
    fit = S.Fit(beta=beta, u=np.array([0.0]), sigma2=0.09, cov=np.eye(len(beta) + 1), names=names, groups=["01"],
                loglik_laplace=0.0, n_cells=0, converged=True)
    N = np.zeros((1, S.N_AGE, S.N_SEX, S.N_RACE, S.N_EDU))
    N[0, 1, 0, 0, 0] = 10.0
    counties = pd.DataFrame({"state_fips": ["21"], "urban3": [0], "adults": [10.0]}, index=pd.Index(["21001"]))
    geo = S.Geo(counties, N, pd.DataFrame(index=counties.index), pd.DataFrame())
    D = 4000
    dr, _ = S.predict_draws(fit, beta[:, None].repeat(D, 1), np.zeros((1, D)), geo, False, True,
                            np.random.default_rng(3))
    lg = np.log(dr[0] / (1 - dr[0]))
    assert lg.mean() == pytest.approx(-2.0, abs=0.02) and lg.std() == pytest.approx(0.3, abs=0.02)


# ------------------------------------------------------------------------------------------------ data
@pytest.mark.data
def test_direct_estimates_and_mmwr_agreement():
    d = read_table("brfss_long_covid_direct_estimates")
    nat = d[(d.geo_id == "US") & (d.measure_id == "lc_current_pct_all_adults_mmwr_def") & (d.estimate_type == "crude")]
    assert nat["value"].iloc[0] == pytest.approx(6.40, abs=0.05)
    cmp_ = pd.read_csv(S.TABLES / "brfss_mmwr_2023_comparison.csv", dtype={"geo_id": str})
    ok = cmp_.dropna(subset=["mmwr_current_pct"])
    st = ok[ok.geo_id != "US50DC"]          # MMWR 'National' includes the territories; ours is 50 states + DC
    assert len(st) == 52 and (st["n_cases_current"] == st["ours_current_n_cases"]).all()
    assert (st["ours_current_std"] - st["mmwr_current_pct"]).abs().max() < 0.5
    assert set(d.loc[d.geo_level == "state", "geo_id"]) >= {"01", "06", "48"}
    assert not {"21", "42"} & set(d.loc[(d.geo_level == "state") & (d.year == 2023), "geo_id"])


@pytest.mark.data
def test_sae_table_contract():
    t = read_table(S.TABLE)
    assert set(t["evidence_type"]) == {"modeled_small_area_estimate"}
    c = t[t.geo_level == "county"]
    for m in S.OUTCOMES:
        cm = c[c.measure_id == m]
        assert len(cm) == 3144 and cm["geo_id"].is_unique
        assert (cm["ci_low"] <= cm["value"]).all() and (cm["value"] <= cm["ci_high"]).all()
        assert (cm["numerator"] <= cm["denominator"]).all()
    assert set(c["source_geographic_resolution"]) == {"county"}
    assert set(t.loc[t.geo_level == "state", "source_geographic_resolution"]) == {"state"}
    assert not c["inherited"].any()
    ky = c[(c.state_fips == "21") & (c.measure_id == S.PRIMARY)]
    assert (~ky["state_in_brfss_2023"]).all() and ky["state_effect_source"].str.startswith("drawn").all()
    # composition-only variation exists within a state but is bounded
    tx = c[(c.state_fips == "48") & (c.measure_id == S.PRIMARY)]["value"]
    assert 0 < tx.max() - tx.min() < 15


@pytest.mark.data
def test_sae_validation_tables():
    sv = pd.read_csv(S.TABLES / f"brfss_sae_state_validation_{S.PRIMARY}.csv", dtype={"geo_id": str})
    have = sv.dropna(subset=["value"])
    assert len(have) == 49 and (have["diff"].abs().mean()) < 0.5
    lo = pd.read_csv(S.TABLES / "brfss_sae_loso.csv", dtype={"state_fips": str})
    assert set(lo["model"]) == {"M0_composition_urbanicity", "M1_plus_acs_covariates"}
    assert lo.groupby("model").size().eq(49).all()
    mm = pd.read_csv(S.TABLES / "brfss_sae_mmsa_validation.csv", dtype={"mmsa": str})
    assert (mm["status"] == "ok").sum() >= 100


@pytest.mark.data
def test_poststrat_cells_table():
    if not table_exists(S.CELLS_TABLE):
        pytest.skip("run measure_it.geography.sae")
    cells = read_table(S.CELLS_TABLE)
    tot = cells.groupby("geo_id")["population"].sum()
    assert len(tot) == 3144 and tot.sum() == pytest.approx(261.4e6, rel=0.02)
