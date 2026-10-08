"""BRFSS small-area estimates (multilevel regression and post-stratification) of adult Long COVID by county.

Why a model: the Household Pulse Survey ended in 2024-09 and was state-level only; the public BRFSS 2023 file asks
current Long COVID but carries no county identifier. A county value can therefore only come from a model. This module
fits one on BRFSS 2023 respondents and post-stratifies it to ACS 2020-2024 county population cells, in the spirit of
CDC PLACES, with the restriction that the public file has no county: see "What the model can and cannot represent".

Model (per outcome; pre-specified, see results/BRFSS_SAE_RESULTS.md):
  logit P(y_i = 1) = age7 x sex + race/ethnicity (4) + education (4) + urbanicity (3: NCHS metro / micropolitan /
                     noncore of the respondent's county, BRFSS _METSTAT x _URBSTAT) [+ county covariates] + u_state,
  u_state ~ N(0, sigma^2). Fitted by survey-weighted pseudo-likelihood on cells (state x urbanicity x age x sex x race
  x education) with each state's weights scaled to its Kish effective sample size, and a Laplace approximation for
  sigma^2 (maximised) and for the joint mode / covariance of the fixed and state effects.
  Candidate county covariates (ACS 2020-2024 % below poverty, % uninsured; none is a Long COVID estimate and PLACES is
  not used because PLACES has no Kentucky / Pennsylvania county values in its BRFSS-2023 release) enter the fit as the
  adult-population-weighted mean over the counties of the respondent's state x urbanicity stratum, and the prediction
  as the county's own value (a cross-level assumption, reported as such).
Post-stratification: ACS county adult population by age7 x sex x race/ethnicity (B01001, B01001H/I/B with B03002 for
  non-Hispanic Black) x education (B15001 shares within county x age band x sex; race and education are combined
  assuming conditional independence within county x age x sex, the IPF solution with these two margins).
Uncertainty: 1,000 draws from the Laplace (Gaussian) approximation of the joint fixed + state-effect posterior; states
  absent from BRFSS 2023 (Kentucky, Pennsylvania) draw their effect from N(0, sigma^2). Sigma^2 is fixed at its estimate.

What the model can and cannot represent: within a state, county estimates differ ONLY through (1) the county's
demographic composition (age, sex, race/ethnicity, education), (2) its urbanicity class and (3) the covariate term.
The respondent's county is never observed, so a county whose residents differ from similar residents elsewhere in
the state for reasons not in the model (a local outbreak history, local vaccination, a local clinic) is invisible.
The state effect is shared by every county of the state.

Outputs: geo_condition_burden__brfss_long_covid_sae (county + state-aggregate + national rows, evidence_type
"modeled_small_area_estimate"), acs_county_adult_poststrat_cells (ACS fine age x sex adult cells, used by
geography.subgroup), draws for geography.subgroup in data/interim/cdc_brfss/, results/BRFSS_SAE_RESULTS.md and
validation tables under results/tables/brfss_sae_*.csv; an SAE section in data/raw/cdc_brfss/DATA_AUDIT.md.

Run: uv run python -m measure_it.geography.sae      (needs measure_it.ingestion.brfss first)
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import pearsonr
from scipy.special import expit

from .. import config
from ..config import RAW, RESULTS, TABLES, interim_dir, raw_dir, utc_now_iso
from ..geography.crosswalk import STATE_FIPS
from ..ingestion import brfss as B
from ..provenance import add_provenance
from ..store import read_table, write_table

PRODUCER = "measure_it.geography.sae"
TABLE = "geo_condition_burden__brfss_long_covid_sae"
CELLS_TABLE = "acs_county_adult_poststrat_cells"
SOURCE_NAME = "Measure-It MRP small-area estimate from CDC BRFSS 2023 + ACS 2020-2024"
N_DRAWS = 1000
N_DRAWS_SAVED = 500
LOSO_DRAWS = 400
PRIMARY = "lc_current_pct_all_adults"
OUTCOMES = {
    "lc_current_pct_all_adults": ("currently_has_long_covid", "Current long COVID (BRFSS 2023 COVIDSM1), modelled",
                                  "adults 18+ (ACS 2020-2024 county population)"),
    "lc_significant_limitation_pct_all_adults": ("activity_limitation_significant",
                                                 "Current long COVID with activity reduced 'a lot', modelled",
                                                 "adults 18+ (ACS 2020-2024 county population)"),
    "lc_current_pct_all_adults_mmwr_def": ("currently_has_long_covid",
                                           "Current long COVID, MMWR definition (break-offs counted as no), modelled",
                                           "adults 18+ (ACS 2020-2024 county population)"),
}
COVARIATES = ["pct_below_poverty", "pct_uninsured"]
SAE_START = "<!-- SAE SECTION START (written by measure_it.geography.sae) -->"
SAE_END = "<!-- SAE SECTION END -->"

N_AGE, N_SEX, N_RACE, N_EDU, N_URB = len(B.AGE7), 2, len(B.RACE4), len(B.EDU4), len(B.URBAN3)
N_DEMO_CELLS = N_AGE * N_SEX * N_RACE * N_EDU

# ACS cell bounds -------------------------------------------------------------------------------------------------
B01001_LOWS = [18, 20, 21, 22, 25, 30, 35, 40, 45, 50, 55, 60, 62, 65, 67, 70, 75, 80, 85]
B01001_HIGHS = [19, 20, 21, 24, 29, 34, 39, 44, 49, 54, 59, 61, 64, 66, 69, 74, 79, 84, None]
RACE_TABLE_LOWS = [18, 20, 25, 30, 35, 45, 55, 65, 75, 85]   # B01001H/I/B adult cells (indent-2 rows 007.. per sex)
EDU_BANDS = {"18 to 24 years": [0], "25 to 34 years": [1], "35 to 44 years": [2], "45 to 64 years": [3, 4],
             "65 years and over": [5, 6]}
EDU_LABEL_TO_EDU4 = {"Less than 9th grade": 0, "9th to 12th grade, no diploma": 0,
                     "High school graduate (includes equivalency)": 1, "Some college, no degree": 2,
                     "Associate's degree": 2, "Bachelor's degree": 3, "Graduate or professional degree": 3}


def age7_of(lo: int) -> int:
    for i, (a, b) in enumerate([(18, 24), (25, 34), (35, 44), (45, 54), (55, 64), (65, 74), (75, 200)]):
        if a <= lo <= b:
            return i
    raise ValueError(lo)


# =================================================================================================== ACS cells
def _read_sf(path: Path) -> pd.DataFrame:
    con = duckdb.connect()
    df = con.execute("SELECT * FROM read_csv(?, delim='|', header=true, all_varchar=true) "
                     "WHERE left(GEO_ID, 9) = '0500000US'", [str(path)]).df()
    con.close()
    df["geo_id"] = df["GEO_ID"].str[9:]
    return df.set_index("geo_id")


def _shells() -> pd.DataFrame:
    return pd.read_csv(RAW / "census_acs" / "ACS20245YR_Table_Shells.txt", sep="|", dtype=str, keep_default_na=False)


def _est(df: pd.DataFrame, uid: str) -> pd.Series:
    t, n = uid.split("_")
    v = pd.to_numeric(df[f"{t}_E{n}"], errors="coerce")
    return v.where(v >= 0)


def _sex_age_cells(shells: pd.DataFrame, table: str) -> list[tuple[str, int, int]]:
    """(unique id, sex 0/1, lower age) for adult sex x age cells of a B01001-type table (labels parsed, not assumed)."""
    rows = shells[shells["Table ID"] == table]
    out, sex = [], None
    for uid, label, indent in zip(rows["Unique ID"], rows["Label"], rows["Indent"]):
        if indent == "1":
            sex = {"Male:": 0, "Female:": 1}[label]
        elif indent == "2":
            m = re.match(r"(Under )?(\d+)", label)
            lo = 0 if m.group(1) else int(m.group(2))
            if lo >= 18:
                out.append((uid, sex, lo))
    return out


def build_poststrat(counties: pd.Index) -> tuple[np.ndarray, pd.DataFrame, dict]:
    """County x (age7, sex, race4, edu4) adult population array + fine B01001 age x sex cells + QA counts."""
    shells = _shells()
    acs = RAW / "census_acs"
    ps = RAW / B.SOURCE_ID / "acs_poststrat"
    tot = _read_sf(acs / "acsdt5y2024-b01001.dat").reindex(counties)
    fine_rows = []
    total7 = np.zeros((len(counties), N_AGE, N_SEX))
    for uid, sex, lo in _sex_age_cells(shells, "B01001"):
        v = _est(tot, uid).to_numpy(float)
        i = B01001_LOWS.index(lo)
        fine_rows.append(pd.DataFrame({"geo_id": counties, "sex": B.SEX2[sex], "age_low": lo,
                                       "age_high": B01001_HIGHS[i] if B01001_HIGHS[i] else 200,
                                       "age7": B.AGE7[age7_of(lo)], "population": v}))
        total7[:, age7_of(lo), sex] += np.nan_to_num(v)
    fine = pd.concat(fine_rows, ignore_index=True)

    race = np.zeros((len(counties), N_AGE, N_SEX, N_RACE))
    for table, r in (("B01001H", 0), ("B01001B", 1), ("B01001I", 2)):
        df = _read_sf(ps / f"acsdt5y2024-{table.lower()}.dat").reindex(counties)
        cells = _sex_age_cells(shells, table)
        assert sorted({c[2] for c in cells}) == RACE_TABLE_LOWS, table
        for uid, sex, lo in cells:
            race[:, age7_of(lo), sex, r] += np.nan_to_num(_est(df, uid).to_numpy(float))
    b3 = _read_sf(ps / "acsdt5y2024-b03002.dat").reindex(counties)
    nh_black, h_black = _est(b3, "B03002_004").to_numpy(float), _est(b3, "B03002_014").to_numpy(float)
    denom = nh_black + h_black
    ratio = np.where(denom > 0, nh_black / np.where(denom > 0, denom, 1), 1.0)
    race[..., 1] *= ratio[:, None, None]
    other = total7 - race[..., 0] - race[..., 1] - race[..., 2]
    qa = {"other_nh_negative_cells": int((other < -0.5).sum()), "other_nh_negative_sum": float(other[other < 0].sum())}
    race[..., 3] = np.clip(other, 0, None)

    # education shares within county x B15001 age band x sex
    e = _read_sf(ps / "acsdt5y2024-b15001.dat").reindex(counties)
    rows = shells[shells["Table ID"] == "B15001"]
    edu_counts = np.zeros((len(counties), N_AGE, N_SEX, N_EDU))
    sex = band = None
    for uid, label, indent in zip(rows["Unique ID"], rows["Label"], rows["Indent"]):
        if indent == "1":
            sex = {"Male:": 0, "Female:": 1}[label]
        elif indent == "2":
            band = EDU_BANDS[label.rstrip(":")]
        elif indent == "3":
            v = np.nan_to_num(_est(e, uid).to_numpy(float))
            for a in band:
                edu_counts[:, a, sex, EDU_LABEL_TO_EDU4[label]] += v
    esum = edu_counts.sum(-1, keepdims=True)
    national = edu_counts.sum(0, keepdims=True) / edu_counts.sum(0, keepdims=True).sum(-1, keepdims=True)
    shares = np.where(esum > 0, edu_counts / np.where(esum > 0, esum, 1), national)
    qa["edu_share_fallback_cells"] = int((esum[..., 0] == 0).sum())
    N = race[..., :, None] * shares[:, :, :, None, :]           # (C, age, sex, race, edu)
    qa["poststrat_total_vs_b01001_adults_max_abs_diff"] = float(np.abs(N.sum((1, 2, 3, 4)) - total7.sum((1, 2))).max())
    return N, fine, qa


# =================================================================================================== urbanicity
def county_urbanicity(counties: pd.Index) -> pd.DataFrame:
    """NCHS 2013 urban-rural code -> metro (1-4) / micropolitan (5) / noncore (6), as in BRFSS _METSTAT/_URBSTAT.

    Counties absent from the 2013 file (2024-vintage Connecticut planning regions, Alaska's 2019 split) take their
    OMB July 2023 delineation status (metropolitan / micropolitan / neither)."""
    raw = raw_dir(B.SOURCE_ID) / "reference"
    n = pd.read_excel(raw / "NCHSURCodes2013.xlsx")
    n["geo_id"] = n["FIPS code"].astype(int).astype(str).str.zfill(5)
    code = n.set_index("geo_id")["2013 code"].astype(int)
    out = pd.DataFrame(index=counties)
    out["nchs_2013_code"] = code.reindex(counties)
    out["urban3"] = out["nchs_2013_code"].map(lambda c: np.nan if pd.isna(c) else (0 if c <= 4 else (1 if c == 5 else 2)))
    out["urbanicity_source"] = np.where(out["urban3"].notna(), "NCHS 2013 urban-rural code", "")
    omb = pd.read_excel(raw / "omb_delineation_list1_2023.xlsx", header=2)
    omb = omb.dropna(subset=["FIPS State Code", "FIPS County Code"])
    omb["geo_id"] = (omb["FIPS State Code"].astype(int).astype(str).str.zfill(2)
                     + omb["FIPS County Code"].astype(int).astype(str).str.zfill(3))
    status = omb.set_index("geo_id")["Metropolitan/Micropolitan Statistical Area"]
    miss = out["urban3"].isna()
    st = status.reindex(out.index[miss])
    out.loc[miss, "urban3"] = st.map({"Metropolitan Statistical Area": 0, "Micropolitan Statistical Area": 1}).fillna(2).to_numpy()
    out.loc[miss, "urbanicity_source"] = "OMB 2023 delineation (not in NCHS 2013 file)"
    out["urban3"] = out["urban3"].astype(int)
    out["urbanicity"] = out["urban3"].map(dict(enumerate(B.URBAN3)))
    return out


# =================================================================================================== design
def demo_design(age7, sex, race4, edu4) -> np.ndarray:
    """Intercept, age (6), female, female x age (6), race (3), education (3)."""
    age7, sex, race4, edu4 = (np.asarray(x, int) for x in (age7, sex, race4, edu4))
    n = len(age7)
    cols = [np.ones(n)]
    cols += [(age7 == a).astype(float) for a in range(1, N_AGE)]
    cols += [sex.astype(float)]
    cols += [((age7 == a) & (sex == 1)).astype(float) for a in range(1, N_AGE)]
    cols += [(race4 == r).astype(float) for r in range(1, N_RACE)]
    cols += [(edu4 == e).astype(float) for e in range(1, N_EDU)]
    return np.column_stack(cols)


DEMO_NAMES = (["intercept"] + [f"age_{B.AGE7[a]}" for a in range(1, N_AGE)] + ["female"]
              + [f"female_x_age_{B.AGE7[a]}" for a in range(1, N_AGE)] + [f"race_{B.RACE4[r]}" for r in range(1, N_RACE)]
              + [f"edu_{B.EDU4[e]}" for e in range(1, N_EDU)])


def urban_design(urban3) -> np.ndarray:
    u = np.asarray(urban3, int)
    return np.column_stack([(u == k).astype(float) for k in range(1, N_URB)])


def demo_cell_grid() -> pd.DataFrame:
    idx = pd.MultiIndex.from_product([range(N_AGE), range(N_SEX), range(N_RACE), range(N_EDU)],
                                     names=["age7", "sex", "race4", "edu4"])
    return idx.to_frame(index=False)


# =================================================================================================== GLMM
@dataclass
class Fit:
    beta: np.ndarray
    u: np.ndarray
    sigma2: float
    cov: np.ndarray
    names: list[str]
    groups: list[str]
    loglik_laplace: float
    n_cells: int
    converged: bool
    extra: dict = field(default_factory=dict)


def _mode(X, Y, W, g, S, sig2, theta0, max_iter=200):
    p = X.shape[1]
    theta = theta0.copy()

    def objective(th):
        eta = X @ th[:p] + th[p:][g]
        return -(Y * eta - W * np.logaddexp(0, eta)).sum() + (th[p:] ** 2).sum() / (2 * sig2)

    f = objective(theta)
    converged = False
    for _ in range(max_iter):
        beta, u = theta[:p], theta[p:]
        eta = X @ beta + u[g]
        mu = expit(eta)
        r = Y - W * mu
        v = W * mu * (1 - mu)
        grad = np.r_[X.T @ r, np.bincount(g, r, S) - u / sig2]
        H = _hessian(X, v, g, S, sig2)
        step = np.linalg.solve(H, grad)
        t = 1.0
        while True:
            cand = theta + t * step
            fc = objective(cand)
            if fc <= f + 1e-10 or t < 1e-6:
                break
            t /= 2
        theta, fold, f = cand, f, fc
        if np.max(np.abs(t * step)) < 1e-9 or abs(fold - f) < 1e-12 * max(1, abs(f)):
            converged = True
            break
    eta = X @ theta[:p] + theta[p:][g]
    mu = expit(eta)
    v = W * mu * (1 - mu)
    return theta, f, _hessian(X, v, g, S, sig2), np.bincount(g, v, S), converged


def _hessian(X, v, g, S, sig2):
    p = X.shape[1]
    Hbb = X.T @ (X * v[:, None])
    Hbu = np.column_stack([np.bincount(g, X[:, j] * v, S) for j in range(p)]).T
    Huu = np.diag(np.bincount(g, v, S) + 1.0 / sig2)
    return np.block([[Hbb, Hbu], [Hbu.T, Huu]])


def fit_glmm(X, Y, W, g, groups: list[str], names: list[str], sigma2: float | None = None) -> Fit:
    """Random-intercept logistic model on aggregated (pseudo-)binomial cells; Laplace approximation."""
    X, Y, W = (np.asarray(a, float) for a in (X, Y, W))
    g = np.asarray(g, int)
    S, p = len(groups), X.shape[1]
    # start: fixed-effects-only logistic (one Newton run with tiny group variance)
    theta0 = np.zeros(p + S)
    theta0[0] = np.log(Y.sum() / (W.sum() - Y.sum()))
    cache = {}

    def neg_marginal(log_s2):
        s2 = float(np.exp(log_s2))
        th0 = cache.get("theta", theta0)
        theta, f, H, vsum, conv = _mode(X, Y, W, g, S, s2, th0)
        cache["theta"] = theta
        return f + 0.5 * S * np.log(s2) + 0.5 * np.log(vsum + 1.0 / s2).sum()

    if sigma2 is None:
        res = minimize_scalar(neg_marginal, bounds=(np.log(1e-5), np.log(4.0)), method="bounded",
                              options={"xatol": 1e-4})
        sigma2 = float(np.exp(res.x))
    theta, f, H, vsum, conv = _mode(X, Y, W, g, S, sigma2, cache.get("theta", theta0))
    ll = -(f + 0.5 * S * np.log(sigma2) + 0.5 * np.log(vsum + 1.0 / sigma2).sum())
    cov = np.linalg.inv(H)
    return Fit(theta[:p], theta[p:], sigma2, cov, list(names), list(groups), float(ll), len(Y), conv)


def draw_params(fit: Fit, n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    theta = np.r_[fit.beta, fit.u]
    L = np.linalg.cholesky(fit.cov + 1e-12 * np.eye(len(theta)))
    z = rng.standard_normal((len(theta), n))
    d = theta[:, None] + L @ z
    p = len(fit.beta)
    return d[:p], d[p:]          # (p, n), (S, n)


# =================================================================================================== data prep
@dataclass
class Prepared:
    frame: pd.DataFrame          # respondent-level analytic frame (50 states + DC)
    cells: pd.DataFrame          # aggregated cells
    exclusions: dict


def kish_scale(w: pd.Series, by: pd.Series) -> pd.Series:
    """Per-group factor so that scaled weights sum to the Kish effective sample size of the group."""
    g = pd.DataFrame({"w": w, "g": by})
    s = g.groupby("g")["w"].agg(lambda x: (x.sum() ** 2 / (x ** 2).sum()) / x.sum())
    return by.map(s)


def prepare(frame: pd.DataFrame, outcome: str, universe: set[str]) -> Prepared:
    d = frame[frame["state_fips"].isin(universe)].copy()
    excl = {"records_50_states_dc": len(d)}
    need = [outcome, "age7", "sex", "race4", "edu4", "urban3"]
    for c in need:
        excl[f"missing_{c}"] = int(d[c].isna().sum())
    d = d.dropna(subset=need)
    excl["analytic_records"] = len(d)
    for c in ["age7", "sex", "race4", "edu4", "urban3"]:
        d[c] = d[c].astype(int)
    d["y"] = d[outcome].astype(float)
    d["ws"] = d["weight"] * kish_scale(d["weight"], d["state_fips"])
    keys = ["state_fips", "urban3", "age7", "sex", "race4", "edu4"]
    cells = d.assign(wy=d["ws"] * d["y"]).groupby(keys, as_index=False).agg(W=("ws", "sum"), Y=("wy", "sum"),
                                                                          n=("y", "size"), cases=("y", "sum"))
    return Prepared(d, cells, excl)


def covariate_frames(counties: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """County z-scored covariates and their state x urbanicity adult-weighted means."""
    acs = read_table("geo_context__acs")
    acs = acs[acs["geo_level"] == "county"].set_index("geo_id")
    z = pd.DataFrame(index=counties.index)
    stats = {}
    w = counties["adults"]
    for c in COVARIATES:
        x = acs[c].reindex(counties.index).astype(float)
        stats[f"{c}_missing_counties"] = int(x.isna().sum())
        x = x.fillna(x.groupby(counties["state_fips"]).transform("median"))
        m = np.average(x, weights=w)
        s = np.sqrt(np.average((x - m) ** 2, weights=w))
        z[c] = (x - m) / s
        stats[c] = {"mean": float(m), "sd": float(s)}
    return z, stats


def group_means(z: pd.DataFrame, counties: pd.DataFrame) -> pd.DataFrame:
    df = z.join(counties[["state_fips", "urban3", "adults"]])
    out = {}
    for c in z.columns:
        num = (df[c] * df["adults"]).groupby([df["state_fips"], df["urban3"]]).sum()
        den = df["adults"].groupby([df["state_fips"], df["urban3"]]).sum()
        out[c] = num / den
    gm = pd.DataFrame(out)
    st = {c: (df[c] * df["adults"]).groupby(df["state_fips"]).sum() / df["adults"].groupby(df["state_fips"]).sum()
          for c in z.columns}
    gm.attrs["state_means"] = pd.DataFrame(st)
    return gm


def cell_design(cells: pd.DataFrame, use_cov: bool, gm: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    X = [demo_design(cells["age7"], cells["sex"], cells["race4"], cells["edu4"]), urban_design(cells["urban3"])]
    names = DEMO_NAMES + [f"urban_{B.URBAN3[k]}" for k in range(1, N_URB)]
    if use_cov:
        key = pd.MultiIndex.from_arrays([cells["state_fips"], cells["urban3"]])
        zz = gm.reindex(key)
        sm = gm.attrs["state_means"].reindex(cells["state_fips"])
        zz = zz.fillna(pd.DataFrame(sm.to_numpy(), index=zz.index, columns=sm.columns))
        X.append(zz.to_numpy(float))
        names += [f"cov_{c}_z" for c in COVARIATES]
    return np.column_stack(X), names


# =================================================================================================== prediction
@dataclass
class Geo:
    counties: pd.DataFrame       # index geo_id; state_fips, urban3, adults
    N: np.ndarray                # (C, age, sex, race, edu)
    z: pd.DataFrame              # county z covariates (raw, unclamped)
    gm: pd.DataFrame             # state x urban covariate means
    zclip: pd.DataFrame | None = None   # county z clamped to the support of the fitting group means


def predict_draws(fit: Fit, beta_d: np.ndarray, u_d: np.ndarray, geo: Geo, use_cov: bool, county_cov: bool,
                  rng: np.random.Generator, idx: np.ndarray | None = None, chunk: int = 150,
                  want_cells: bool = False) -> tuple[np.ndarray, np.ndarray | None]:
    """County prevalence draws (C, D) [and county x age7 x sex draws (C, 14, D)] by post-stratification."""
    C_all = geo.counties
    idx = np.arange(len(C_all)) if idx is None else idx
    D = beta_d.shape[1]
    grid = demo_cell_grid()
    Xd = demo_design(grid["age7"], grid["sex"], grid["race4"], grid["edu4"])
    pd_ = Xd.shape[1]
    a = (Xd @ beta_d[:pd_]).astype(np.float32)                         # (K, D)
    b_urban = np.vstack([np.zeros(D), beta_d[pd_:pd_ + N_URB - 1]])    # (3, D)
    # state effects, absent states drawn from N(0, sigma^2)
    gpos = {s: i for i, s in enumerate(fit.groups)}
    states = C_all["state_fips"].to_numpy()
    uniq = np.unique(states[idx])
    u_state = {}
    for s in uniq:
        u_state[s] = u_d[gpos[s]] if s in gpos else rng.normal(0.0, np.sqrt(fit.sigma2), D)
    out = np.empty((len(idx), D), np.float32)
    cells_out = np.empty((len(idx), N_AGE * N_SEX, D), np.float32) if want_cells else None
    for start in range(0, len(idx), chunk):
        ii = idx[start:start + chunk]
        sub = C_all.iloc[ii]
        bc = b_urban[sub["urban3"].to_numpy()] + np.vstack([u_state[s] for s in sub["state_fips"]])
        if use_cov:
            if county_cov:
                zc = (geo.zclip if geo.zclip is not None else geo.z).iloc[ii][COVARIATES].to_numpy(float)
            else:
                key = pd.MultiIndex.from_arrays([sub["state_fips"], sub["urban3"]])
                zc = geo.gm.reindex(key)[COVARIATES].to_numpy(float)
            bc = bc + zc @ beta_d[pd_ + N_URB - 1:pd_ + N_URB - 1 + len(COVARIATES)]
        p = expit(a[None, :, :] + bc[:, None, :].astype(np.float32))   # (c, K, D)
        Nk = geo.N[ii].reshape(len(ii), -1).astype(np.float32)          # (c, K)
        tot = Nk.sum(1)
        out[start:start + len(ii)] = np.einsum("ck,ckd->cd", Nk, p) / np.where(tot > 0, tot, 1)[:, None]
        if want_cells:
            pk = (Nk[:, :, None] * p).reshape(len(ii), N_AGE * N_SEX, N_RACE * N_EDU, D).sum(2)
            nk = Nk.reshape(len(ii), N_AGE * N_SEX, N_RACE * N_EDU).sum(2)
            cells_out[start:start + len(ii)] = pk / np.where(nk > 0, nk, 1)[:, :, None]
    return out, cells_out


def aggregate(draws: np.ndarray, pop: np.ndarray, groups: np.ndarray) -> pd.DataFrame:
    """Population-weighted aggregate of county draws by group -> mean, 2.5%, 97.5%."""
    rows = []
    for gid in pd.unique(groups):
        m = groups == gid
        agg = (pop[m, None] * draws[m]).sum(0) / pop[m].sum()
        rows.append({"geo_id": gid, "mean": float(agg.mean()), "lo": float(np.quantile(agg, 0.025)),
                     "hi": float(np.quantile(agg, 0.975)), "pop": float(pop[m].sum())})
    return pd.DataFrame(rows)


# =================================================================================================== validation
def loso(prep: Prepared, geo: Geo, direct: pd.Series, use_cov: bool, seed: int) -> pd.DataFrame:
    """Leave-one-state-out: refit without state s (sigma re-estimated), predict s by post-stratification."""
    rng = np.random.default_rng(seed)
    rows = []
    cells = prep.cells
    states = sorted(cells["state_fips"].unique())
    pop = geo.counties["adults"].to_numpy(float)
    st_arr = geo.counties["state_fips"].to_numpy()
    nat = prep.frame
    for s in states:
        c = cells[cells["state_fips"] != s]
        groups = sorted(c["state_fips"].unique())
        gidx = c["state_fips"].map({g: i for i, g in enumerate(groups)}).to_numpy()
        X, names = cell_design(c, use_cov, geo.gm)
        f = fit_glmm(X, c["Y"], c["W"], gidx, groups, names)
        bd, ud = draw_params(f, LOSO_DRAWS, rng)
        idx = np.where(st_arr == s)[0]
        dr, _ = predict_draws(f, bd, ud, geo, use_cov, True, rng, idx=idx)
        agg = (pop[idx, None] * dr).sum(0) / pop[idx].sum()
        rest = nat[nat["state_fips"] != s]
        base = float(np.average(rest["y"], weights=rest["weight"]))
        rows.append({"state_fips": s, "direct": float(direct.get(s, np.nan)), "loso_pred": 100 * float(agg.mean()),
                     "loso_lo": 100 * float(np.quantile(agg, 0.025)), "loso_hi": 100 * float(np.quantile(agg, 0.975)),
                     "baseline_national_rest": 100 * base, "sigma2": f.sigma2})
    out = pd.DataFrame(rows)
    out["err_model"] = out["loso_pred"] - out["direct"]
    out["err_baseline"] = out["baseline_national_rest"] - out["direct"]
    out["covered"] = (out["direct"] >= out["loso_lo"]) & (out["direct"] <= out["loso_hi"])
    return out


def loso_summary(t: pd.DataFrame) -> dict:
    return {"n_states": int(len(t)), "mae_model": float(t["err_model"].abs().mean()),
            "rmse_model": float(np.sqrt((t["err_model"] ** 2).mean())),
            "mae_baseline": float(t["err_baseline"].abs().mean()),
            "rmse_baseline": float(np.sqrt((t["err_baseline"] ** 2).mean())),
            "pearson_model_vs_direct": float(np.corrcoef(t["loso_pred"], t["direct"])[0, 1]),
            "coverage_95": float(t["covered"].mean()), "mean_bias": float(t["err_model"].mean())}


def mmsa_county_map(counties: pd.Index) -> pd.DataFrame:
    """MMSA code -> member counties from the OMB March 2020 delineation (SMART BRFSS 2023 vintage)."""
    raw = raw_dir(B.SOURCE_ID) / "reference" / "omb_delineation_list1_2020.xls"
    o = pd.read_excel(raw, header=2).dropna(subset=["FIPS State Code", "FIPS County Code"])
    o["geo_id"] = (o["FIPS State Code"].astype(int).astype(str).str.zfill(2)
                   + o["FIPS County Code"].astype(int).astype(str).str.zfill(3))
    rows = []
    for _, r in o.iterrows():
        rows.append({"mmsa": str(int(r["CBSA Code"])), "geo_id": r["geo_id"]})
        if pd.notna(r["Metropolitan Division Code"]):
            rows.append({"mmsa": str(int(r["Metropolitan Division Code"])), "geo_id": r["geo_id"]})
    m = pd.DataFrame(rows).drop_duplicates()
    m["in_2024_counties"] = m["geo_id"].isin(counties)
    return m


def mmsa_validation(county_est: pd.DataFrame, comp_est: pd.DataFrame, geo: Geo, state_direct: pd.Series,
                    state_sae: pd.Series) -> tuple[pd.DataFrame, dict]:
    mm = pd.read_csv(TABLES / "brfss_mmsa_2023_direct.csv", dtype={"mmsa": str})
    cmap = mmsa_county_map(geo.counties.index)
    pop = geo.counties["adults"]
    rows = []
    for _, r in mm.iterrows():
        mem = cmap[cmap["mmsa"] == r["mmsa"]]
        if mem.empty:
            rows.append({**r.to_dict(), "status": "no counties in OMB 2020 delineation"})
            continue
        if not mem["in_2024_counties"].all():
            rows.append({**r.to_dict(), "status": "member counties not in 2024 vintage (e.g. legacy CT)",
                         "n_counties": len(mem)})
            continue
        ids = mem["geo_id"].to_numpy()
        w = pop.reindex(ids).to_numpy(float)
        sts = geo.counties.loc[ids, "state_fips"]
        rows.append({**r.to_dict(), "status": "ok", "n_counties": len(ids),
                     "states": ",".join(sorted(sts.unique())),
                     "sae": float(np.average(county_est.reindex(ids), weights=w)),
                     "sae_composition_only": float(np.average(comp_est.reindex(ids), weights=w)),
                     "inherited_state_direct": float(np.average(state_direct.reindex(sts).to_numpy(float), weights=w)),
                     "state_sae": float(np.average(state_sae.reindex(sts).to_numpy(float), weights=w))})
    t = pd.DataFrame(rows)
    ok = t[t["status"] == "ok"].dropna(subset=["value", "inherited_state_direct"]).copy()
    ok["dev_direct"] = ok["value"] - ok["inherited_state_direct"]
    ok["dev_sae"] = ok["sae"] - ok["state_sae"]
    ok["dev_comp"] = ok["sae_composition_only"] - ok["state_sae"]
    signal_var = float(ok["dev_direct"].var() - (ok["se"] ** 2).mean())
    s = {
        "n_mmsa_total": int(len(t)), "n_mmsa_matched": int(len(ok)),
        "unmatched": t.loc[t["status"] != "ok", ["mmsa_name", "status"]].to_dict("records"),
        "mae_sae": float((ok["sae"] - ok["value"]).abs().mean()),
        "mae_inherited_state": float((ok["inherited_state_direct"] - ok["value"]).abs().mean()),
        "mae_sae_composition_only": float((ok["sae_composition_only"] - ok["value"]).abs().mean()),
        "pearson_dev_sae_vs_direct": float(np.corrcoef(ok["dev_sae"], ok["dev_direct"])[0, 1]),
        "p_dev_sae_vs_direct": float(pearsonr(ok["dev_sae"], ok["dev_direct"]).pvalue),
        "pearson_dev_comp_vs_direct": float(np.corrcoef(ok["dev_comp"], ok["dev_direct"])[0, 1]),
        "sd_dev_direct": float(ok["dev_direct"].std()), "mean_se_direct": float(np.sqrt((ok["se"] ** 2).mean())),
        "signal_sd_dev_direct": float(np.sqrt(max(signal_var, 0.0))), "signal_var_raw": signal_var,
        "sd_dev_sae": float(ok["dev_sae"].std()),
        "slope_dev_direct_on_dev_sae": float(np.polyfit(ok["dev_sae"], ok["dev_direct"], 1)[0]),
        "mean_dev_direct": float(ok["dev_direct"].mean()), "mean_dev_sae": float(ok["dev_sae"].mean()),
    }
    return t, s


# =================================================================================================== within-state
def within_state_variation(est: pd.Series, geo: Geo) -> tuple[pd.DataFrame, dict]:
    df = geo.counties[["state_fips", "adults"]].assign(v=est)
    rows = []
    for s, g in df.groupby("state_fips"):
        m = np.average(g["v"], weights=g["adults"])
        rows.append({"state_fips": s, "state_abbr": STATE_FIPS.get(s, s), "n_counties": len(g),
                     "pop_weighted_mean": m, "min": g["v"].min(), "max": g["v"].max(),
                     "range_pp": g["v"].max() - g["v"].min(),
                     "pop_weighted_sd_pp": float(np.sqrt(np.average((g["v"] - m) ** 2, weights=g["adults"])))})
    t = pd.DataFrame(rows)
    sm = df.groupby("state_fips").apply(lambda g: np.average(g["v"], weights=g["adults"]), include_groups=False)
    within = df["v"] - df["state_fips"].map(sm)
    tot_var = df["v"].var()
    s = {"county_sd_pp": float(df["v"].std()), "within_state_sd_pp": float(within.std()),
         "share_variance_within_state": float(within.var() / tot_var),
         "median_state_range_pp": float(t["range_pp"].median()), "max_state_range_pp": float(t["range_pp"].max()),
         "median_within_state_sd_pp": float(t["pop_weighted_sd_pp"].median())}
    return t, s


# =================================================================================================== outputs
def burden_rows(measure: str, county: pd.DataFrame, state: pd.DataFrame, national: pd.DataFrame, geo: Geo,
                fit: Fit, model_id: str, retrieved: str, version: str, absent: set[str]) -> pd.DataFrame:
    family, label, denom = OUTCOMES[measure]
    g = read_table("geographies")
    names = g.set_index("geo_id")["name"]
    st = g[g["geo_level"] == "state"].set_index("state_fips")["name"]
    parts = []
    for level, t in (("county", county), ("state", state), ("national", national)):
        t = t.copy()
        t["geo_level"] = level
        parts.append(t)
    df = pd.concat(parts, ignore_index=True)
    df["state_fips"] = np.where(df["geo_level"] == "county", df["geo_id"].str[:2],
                                np.where(df["geo_level"] == "state", df["geo_id"], None))
    out = pd.DataFrame({
        "geo_id": df["geo_id"], "geo_level": df["geo_level"],
        "geo_name": df["geo_id"].map(names).fillna(df["geo_id"].map({"US": "United States (50 states + DC)"})),
        "condition_id": "long_covid", "measure_id": measure, "measure_label": label, "indicator_family": family,
        "denominator_population": denom, "metric_type": "modeled_prevalence_pct",
        "value": 100 * df["mean"], "value_unit": "percent", "ci_low": 100 * df["lo"], "ci_high": 100 * df["hi"],
        "ci_level": "95% (model; Laplace-approximation posterior draws, sigma fixed)",
        "state_fips": df["state_fips"], "state_abbr": df["state_fips"].map(STATE_FIPS),
        "state_name": df["state_fips"].map(st),
        "county_fips": np.where(df["geo_level"] == "county", df["geo_id"], None),
        "geo_vintage": np.where(df["geo_level"] == "national", "national", "2024"),
        "in_canonical_geographies": df["geo_id"].isin(set(g["geo_id"])),
        "numerator": (df["mean"] * df["pop"]).round(0), "numerator_ci_low": (df["lo"] * df["pop"]).round(0),
        "numerator_ci_high": (df["hi"] * df["pop"]).round(0), "denominator": df["pop"].round(0),
        "denominator_source": "ACS 2020-2024 5-year, B01001 population 18+ (incl. group quarters)",
        "primary_burden_measure": measure == PRIMARY,
        "burden_evidence_level": "A", "year": 2023, "period_start": "2023-01-01", "period_end": "2023-12-31",
        "period_label": "BRFSS 2023", "suppressed": False, "estimate_kind": "modeled", "inherited": False,
        "derivation": "mrp_small_area_estimate", "model_id": model_id,
        "state_in_brfss_2023": ~df["state_fips"].isin(absent) & df["geo_level"].ne("national"),
        "state_effect_source": np.where(df["state_fips"].isin(absent),
                                        "drawn from N(0, sigma^2): state absent from the BRFSS 2023 public file",
                                        np.where(df["geo_level"] == "national", "all states", "estimated")),
        "within_state_variation_basis": "county demographic composition (age, sex, race/ethnicity, education), "
                                        "urbanicity class and ACS covariate term only; no county-specific survey data",
        "urbanicity": df["geo_id"].map(geo.counties["urbanicity"]) if "urbanicity" in geo.counties else None,
        "survey": "CDC BRFSS 2023 (public LLCP file) + ACS 2020-2024",
    })
    out.loc[out["geo_level"] != "county", "urbanicity"] = None
    res = np.where(out["geo_level"] == "county", "county", np.where(out["geo_level"] == "state", "state", "national"))
    out = add_provenance(
        out, data_layer="geographic", source_name=SOURCE_NAME, source_version=version, retrieved_at=retrieved,
        evidence_type="modeled_small_area_estimate",
        source_record_id=lambda d: "brfss_sae:" + d["measure_id"] + ":" + d["geo_level"] + ":" + d["geo_id"],
        source_geographic_resolution="county", evidence_level="A",
        provenance_notes=f"Multilevel regression and post-stratification ({model_id}); modelled, not observed county "
                         "prevalence. Within-state variation comes only from composition, urbanicity and covariates. "
                         "Method and validation: data/raw/cdc_brfss/DATA_AUDIT.md, results/BRFSS_SAE_RESULTS.md.")
    out["source_geographic_resolution"] = res
    return out


def save_draws(measure: str, geo: Geo, cells: np.ndarray, county: np.ndarray, rng) -> Path:
    keep = rng.choice(cells.shape[2], size=min(N_DRAWS_SAVED, cells.shape[2]), replace=False)
    p = interim_dir(B.SOURCE_ID) / f"sae_draws_{measure}.npz"
    labels = [f"{B.AGE7[a]}|{B.SEX2[s]}" for a in range(N_AGE) for s in range(N_SEX)]
    np.savez(p, geo_id=np.array([str(x) for x in geo.counties.index], dtype="U5"), cell_labels=np.array(labels, dtype="U16"),
             cell_draws=cells[:, :, keep].astype(np.float32), county_draws=county[:, keep].astype(np.float32))
    return p


def load_draws(measure: str = PRIMARY) -> dict | None:
    p = interim_dir(B.SOURCE_ID) / f"sae_draws_{measure}.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return {k: z[k] for k in z.files}


# =================================================================================================== run
def _fmt(x, nd=2):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def build_geo() -> tuple[Geo, pd.DataFrame, dict]:
    g = read_table("geographies")
    universe = {f for f in g.loc[g["geo_level"] == "state", "geo_id"] if int(f) <= 56}
    cty = g[(g["geo_level"] == "county") & (~g["ct_legacy"].fillna(False).astype(bool))
            & g["state_fips"].isin(universe)].sort_values("geo_id")
    idx = pd.Index(cty["geo_id"].to_numpy(), name="geo_id")
    N, fine, qa = build_poststrat(idx)
    counties = pd.DataFrame(index=idx)
    counties["state_fips"] = cty.set_index("geo_id")["state_fips"].reindex(idx)
    counties["adults"] = N.sum((1, 2, 3, 4))
    urb = county_urbanicity(idx)
    counties = counties.join(urb)
    z, cstats = covariate_frames(counties)
    gm = group_means(z, counties)
    qa.update(cstats)
    # No extrapolation beyond the fitted support: the covariate coefficients are estimated from state x urbanicity
    # means, so a county value outside the range of those means is clamped to it (deviation added 2026-10-07 after
    # the first run produced 21% for Oglala Lakota and 14% for Amish-majority LaGrange/Holmes counties, whose ACS
    # 'uninsured' share (religious non-participation) lies far outside the fitted range).
    lo, hi = gm.min(), gm.max()
    zclip = z.clip(lower=lo, upper=hi, axis=1)
    qa["covariate_support"] = {c: [float(lo[c]), float(hi[c])] for c in z.columns}
    qa["counties_clamped"] = {c: int((z[c] != zclip[c]).sum()) for c in z.columns}
    qa["counties"] = len(idx)
    qa["counties_zero_adults"] = int((counties["adults"] <= 0).sum())
    qa["urbanicity_counts"] = counties["urbanicity"].value_counts().to_dict()
    qa["urbanicity_from_omb"] = counties.loc[counties["urbanicity_source"].str.startswith("OMB"), "urbanicity"].to_dict()
    return Geo(counties, N, z, gm, zclip), fine, qa


def write_cells_table(fine: pd.DataFrame, retrieved: str) -> None:
    f = fine.copy()
    f["state_fips"] = f["geo_id"].str[:2]
    f["geo_level"] = "county"
    out = add_provenance(
        f, data_layer="geographic", source_name="U.S. Census Bureau ACS 2020-2024 5-year, B01001 (Summary File)",
        source_version="ACS 2020-2024 5-year", retrieved_at=retrieved, evidence_type="census_estimate",
        source_record_id=lambda d: "B01001:" + d["geo_id"] + ":" + d["sex"] + ":" + d["age_low"].astype(str),
        source_geographic_resolution="county",
        provenance_notes="Adult sex x age cells (ACS bands) used for post-stratification and subgroup burden.")
    write_table(out, CELLS_TABLE, producer=PRODUCER, description="ACS county adult population by sex x ACS age band")


def run() -> dict:
    rng = np.random.default_rng(config.SEED)
    t0 = utc_now_iso()
    frame = B.analysis_frame(2023)
    geo, fine, qa = build_geo()
    universe = set(geo.counties["state_fips"].unique())
    absent = universe - set(frame["state_fips"].unique())
    acs_manifest = json.loads((RAW / "census_acs" / "MANIFEST.json").read_text())["files"]
    write_cells_table(fine, acs_manifest["acsdt5y2024-b01001.dat"]["retrieved_at"])
    direct = read_table("brfss_long_covid_direct_estimates")
    direct = direct[(direct["year"] == 2023) & (direct["estimate_type"] == "crude")]
    bmanifest = json.loads((RAW / B.SOURCE_ID / "MANIFEST.json").read_text())["files"]
    retrieved = bmanifest["2023/LLCP2023XPT.zip"]["retrieved_at"]
    pop = geo.counties["adults"].to_numpy(float)
    st_arr = geo.counties["state_fips"].to_numpy()

    results, frames, val_tables = {}, [], {}
    # ---------------- model selection on the primary outcome (pre-specified rule: M1 unless LOSO RMSE worse than M0)
    prep = prepare(frame, PRIMARY, universe)
    d_primary = direct[direct["measure_id"] == PRIMARY].set_index("geo_id")["value"]
    loso_tabs = {}
    for mid, use_cov in (("M0_composition_urbanicity", False), ("M1_plus_acs_covariates", True)):
        loso_tabs[mid] = loso(prep, geo, d_primary, use_cov, config.SEED)
    loso_sum = {k: loso_summary(v) for k, v in loso_tabs.items()}
    chosen = ("M1_plus_acs_covariates" if loso_sum["M1_plus_acs_covariates"]["rmse_model"]
              <= loso_sum["M0_composition_urbanicity"]["rmse_model"] else "M0_composition_urbanicity")
    use_cov = chosen.startswith("M1")
    lt = pd.concat([v.assign(model=k) for k, v in loso_tabs.items()], ignore_index=True)
    lt["state_abbr"] = lt["state_fips"].map(STATE_FIPS)
    lt.to_csv(TABLES / "brfss_sae_loso.csv", index=False)

    version = (f"BRFSS 2023 LLCP public file + ACS 2020-2024 5-year; model {chosen}; run {t0}")
    coef_rows = []
    for measure in OUTCOMES:
        prep_m = prep if measure == PRIMARY else prepare(frame, measure, universe)
        groups = sorted(prep_m.cells["state_fips"].unique())
        gidx = prep_m.cells["state_fips"].map({g: i for i, g in enumerate(groups)}).to_numpy()
        X, names = cell_design(prep_m.cells, use_cov, geo.gm)
        fit = fit_glmm(X, prep_m.cells["Y"], prep_m.cells["W"], gidx, groups, names)
        se = np.sqrt(np.diag(fit.cov))
        coef_rows += [{"measure_id": measure, "term": n, "estimate": b, "se": s, "odds_ratio": np.exp(b)}
                      for n, b, s in zip(names, fit.beta, se[:len(names)])]
        coef_rows += [{"measure_id": measure, "term": f"state_effect_{g}", "estimate": u, "se": s}
                      for g, u, s in zip(groups, fit.u, se[len(names):])]
        coef_rows.append({"measure_id": measure, "term": "sigma_state", "estimate": float(np.sqrt(fit.sigma2))})
        bd, ud = draw_params(fit, N_DRAWS, rng)
        cdraw, celldraw = predict_draws(fit, bd, ud, geo, use_cov, True, rng, want_cells=True)
        county = pd.DataFrame({"geo_id": geo.counties.index, "mean": cdraw.mean(1),
                               "lo": np.quantile(cdraw, 0.025, axis=1), "hi": np.quantile(cdraw, 0.975, axis=1),
                               "pop": pop})
        state = aggregate(cdraw, pop, st_arr)
        national = aggregate(cdraw, pop, np.array(["US"] * len(pop)))
        frames.append(burden_rows(measure, county, state, national, geo, fit, chosen, retrieved, version, absent))
        save_draws(measure, geo, celldraw, cdraw, rng)
        # composition-only variant (covariate at its state x urbanicity mean) for the within-state analysis
        if use_cov:
            comp, _ = predict_draws(fit, bd[:, :300], ud[:, :300], geo, use_cov, False, rng)
            comp_mean = pd.Series(comp.mean(1), index=geo.counties.index)
        else:
            comp_mean = pd.Series(cdraw.mean(1), index=geo.counties.index)
        cest = pd.Series(cdraw.mean(1), index=geo.counties.index)
        dm = direct[direct["measure_id"] == measure].set_index("geo_id")
        sv = state.set_index("geo_id").join(dm[["value", "ci_low", "ci_high", "n_unweighted"]], how="left")
        sv["sae"] = 100 * sv["mean"]
        sv["sae_lo"], sv["sae_hi"] = 100 * sv["lo"], 100 * sv["hi"]
        sv["diff"] = sv["sae"] - sv["value"]
        sv["direct_in_sae_interval"] = (sv["value"] >= sv["sae_lo"]) & (sv["value"] <= sv["sae_hi"])
        sv["state_abbr"] = sv.index.map(STATE_FIPS)
        sv["in_brfss_2023"] = ~sv.index.isin(absent)
        sv = sv.reset_index()
        sv.to_csv(TABLES / f"brfss_sae_state_validation_{measure}.csv", index=False)
        have = sv.dropna(subset=["value"])
        wsv, wss = within_state_variation(100 * cest, geo)
        wsv_c, wss_c = within_state_variation(100 * comp_mean, geo)
        res = {"fit": {"sigma_state": float(np.sqrt(fit.sigma2)), "n_cells": fit.n_cells, "converged": fit.converged,
                       "states_with_data": len(groups), "exclusions": prep_m.exclusions},
               "state_validation": {"n": int(len(have)), "mae_pp": float(have["diff"].abs().mean()),
                                    "max_abs_pp": float(have["diff"].abs().max()),
                                    "pearson": float(np.corrcoef(have["sae"], have["value"])[0, 1]),
                                    "coverage": float(have["direct_in_sae_interval"].mean()),
                                    "sd_direct": float(have["value"].std()), "sd_sae": float(have["sae"].std())},
               "national": national.iloc[0].to_dict(),
               "within_state": wss, "within_state_composition_only": wss_c,
               "county_range": [float(100 * cest.min()), float(100 * cest.max())],
               "ci_width_median_pp": float(100 * (county["hi"] - county["lo"]).median())}
        if measure == PRIMARY:
            wsv.to_csv(TABLES / "brfss_sae_within_state_variation.csv", index=False)
            mt, ms = mmsa_validation(cest * 100, comp_mean * 100, geo,
                                     d_primary, state.set_index("geo_id")["mean"] * 100)
            mt.to_csv(TABLES / "brfss_sae_mmsa_validation.csv", index=False)
            res["mmsa"] = ms
            res["state_table"] = sv
            res["county_table"] = county
        results[measure] = res
    pd.DataFrame(coef_rows).to_csv(TABLES / "brfss_sae_model_coefficients.csv", index=False)
    out = pd.concat(frames, ignore_index=True)
    write_table(out, TABLE, producer=PRODUCER,
                description=f"BRFSS 2023 MRP county Long COVID estimates ({chosen}); modelled, not observed")
    stats = {"national_direct": float(direct.loc[(direct["measure_id"] == PRIMARY) & (direct["geo_id"] == "US"),
                                                  "value"].iloc[0]),
             "chosen": chosen, "loso": loso_sum, "loso_table": lt, "qa": qa, "absent": sorted(absent),
             "results": results, "rows": len(out), "run_at": t0}
    write_results(stats)
    write_sae_audit(stats)
    summary = {"chosen": chosen, "loso": loso_sum, "absent": sorted(absent), "rows": len(out),
               **{m: {k: v for k, v in r.items() if k not in ("state_table", "county_table")}
                  for m, r in results.items()}}
    print(json.dumps(summary, indent=2, default=str))
    return summary


# =================================================================================================== reporting
def write_results(s: dict) -> Path:
    r = s["results"][PRIMARY]
    sig = s["results"]["lc_significant_limitation_pct_all_adults"]
    mmwr = s["results"]["lc_current_pct_all_adults_mmwr_def"]
    sv, ws, wc, mm = r["state_validation"], r["within_state"], r["within_state_composition_only"], r["mmsa"]
    L0, L1 = s["loso"]["M0_composition_urbanicity"], s["loso"]["M1_plus_acs_covariates"]
    st = r["state_table"].sort_values("diff")
    lt = s["loso_table"][s["loso_table"]["model"] == s["chosen"]].sort_values("err_model")
    worst_loso = "; ".join(f"{x.state_abbr} {x.loso_pred:.2f} vs {x.direct:.2f}" for x in
                           pd.concat([lt.head(3), lt.tail(3)]).itertuples())
    wsv = pd.read_csv(TABLES / "brfss_sae_within_state_variation.csv", dtype={"state_fips": str})
    top_range = "; ".join(f"{x.state_abbr} {x.min:.2f}-{x.max:.2f}" for x in
                          wsv.sort_values("range_pp", ascending=False).head(5).itertuples())
    coefs = pd.read_csv(TABLES / "brfss_sae_model_coefficients.csv")
    cp = coefs[(coefs["measure_id"] == PRIMARY) & ~coefs["term"].str.startswith("state_effect")]
    coef_md = "\n".join(f"| {x.term} | {x.estimate:.3f} | {_fmt(x.se, 3)} | {_fmt(x.odds_ratio, 2)} |"
                        for x in cp.itertuples())
    absent = ", ".join(STATE_FIPS.get(a, a) for a in s["absent"])
    nat = r["national"]
    text = f"""# BRFSS small-area estimates of adult Long COVID (county), 2023

Generated by `measure_it.geography.sae` at {s['run_at']}. Inputs: BRFSS 2023 LLCP public file (`measure_it.ingestion.brfss`,
audit `data/raw/cdc_brfss/DATA_AUDIT.md`), ACS 2020-2024 5-year Summary File tables, NCHS 2013 urban-rural codes,
ACS covariates from `geo_context__acs`. Output: `data/processed/{TABLE}.parquet`. Every number below is read from the
tables this run wrote (`results/tables/brfss_sae_*.csv`).

## What this is, and what it is not
A multilevel regression and post-stratification (MRP) model. The public BRFSS file has **no county identifier**, so no
county is observed: a county estimate is the model's prediction for that county's adult population. Within a state,
counties differ **only** through (1) their demographic composition (age x sex, race/ethnicity, education), (2) their
urbanicity class (metro / micropolitan / noncore; observed for respondents through `_METSTAT`/`_URBSTAT`) and (3) the
ACS covariate term (if selected). Every county of a state shares that state's random effect. A local excess not
explained by these terms (outbreak history, vaccination, local care) is invisible to the model. These are modelled
estimates, not observed county prevalence, and they are labelled `evidence_type = modeled_small_area_estimate`.
{absent} are absent from the 2023 public file: their counties use a state effect drawn from N(0, sigma^2), so their
values are fixed-effect predictions with wider intervals (`state_in_brfss_2023 = False`).

## Model (pre-specified before fitting)
* Outcome (primary): current Long COVID, complete case (`lc_current_pct_all_adults`; COVIDPO1 = 1 and COVIDSM1 = 1;
  no = COVIDPO1 = 2 or COVIDSM1 = 2; don't know / refused / break-off excluded). Also fitted: significant activity
  limitation (`lc_significant_limitation_pct_all_adults`) and the MMWR-definition outcome (`..._mmwr_def`, break-offs = no).
* Fixed effects: age (7 bands) x sex, race/ethnicity (`_IMPRACE`: white NH, Black NH, Hispanic, other NH), education
  (`_EDUCAG`, 4), urbanicity (3). Random effect: state intercept. Income is not used (no county joint cells; 20% unknown).
* Candidate covariates (M1): ACS % below poverty and % uninsured (z-scored), entered at the state x urbanicity mean for
  respondents and at the county value for prediction. Not used: any Long COVID estimate; PLACES measures (missing for all
  Kentucky and Pennsylvania counties in the BRFSS-2023 PLACES release, and symptom measures such as frequent physical
  distress or disability would be partly caused by Long COVID); SVI (overlaps the ACS covariates).
* Selection rule (fixed in advance): use M1 unless its leave-one-state-out RMSE is larger than M0's. **Chosen: {s['chosen']}.**
* Estimation: survey-weighted pseudo-likelihood on {r['fit']['n_cells']:,} cells, weights scaled per state to the Kish
  effective sample size, Laplace approximation (sigma_state = {r['fit']['sigma_state']:.3f} on the logit scale); 1,000 draws from the
  Gaussian approximation of the joint fixed + state-effect posterior (sigma fixed). Analytic records: {r['fit']['exclusions']['analytic_records']:,}
  of {r['fit']['exclusions']['records_50_states_dc']:,} in the 50 states + DC (exclusions: {', '.join(f'{k} {v:,}' for k, v in r['fit']['exclusions'].items() if k.startswith('missing'))}).
* Post-stratification: {s['qa']['counties']:,} counties (2024 vintage, CT planning regions) x 224 cells; ACS adults 18+
  ({nat['pop']:,.0f}); non-Hispanic Black cells = B01001B x county NH share of Black (B03002); other NH = total - white NH -
  Hispanic - Black NH ({s['qa']['other_nh_negative_cells']} cells clipped at 0); education shares from B15001 assumed
  independent of race within county x age x sex. Urbanicity: NCHS 2013 codes, OMB 2023 status for {len(s['qa']['urbanicity_from_omb'])}
  counties not in the 2013 file.

Deviation from the pre-specified model (2026-10-07, after the first run): county covariate values are clamped to the
range of the state x urbanicity means the coefficients were estimated on ({', '.join(f"{k}: {v} counties" for k, v in s['qa']['counties_clamped'].items())}
clamped). Without the clamp the first run extrapolated to 21.2% (Oglala Lakota, SD; 57.6% poverty, 44.6% uninsured)
and 14.2% / 13.9% (LaGrange, IN / Holmes, OH: Amish-majority counties whose ACS 'uninsured' share of ~44% reflects
religious non-participation in insurance), and the within-state SD was 0.871 pp (49% of county variance within states).

Primary-model coefficients (log-odds; `results/tables/brfss_sae_model_coefficients.csv`):

| term | estimate | se | odds ratio |
|---|---|---|---|
{coef_md}

## Results
* National (50 states + DC, ACS adult population): {100 * nat['mean']:.2f}% ({100 * nat['lo']:.2f}-{100 * nat['hi']:.2f}); BRFSS direct
  complete-case national {s['national_direct']:.2f}%. MMWR-definition model: {100 * mmwr['national']['mean']:.2f}%. Significant limitation:
  {100 * sig['national']['mean']:.2f}% of adults.
* County estimates range {r['county_range'][0]:.2f}% to {r['county_range'][1]:.2f}%; median 95% interval width {r['ci_width_median_pp']:.2f} pp.

## Validation (a): county estimates aggregated back to the state vs the direct state estimate
In-sample (the state effect is fitted to the same respondents), so this checks the post-stratification and the
shrinkage, not out-of-sample skill. Population bases differ (ACS adults incl. group quarters vs BRFSS weights).

| measure | states | MAE (pp) | max abs (pp) | Pearson r | direct inside SAE 95% interval | SD direct | SD SAE |
|---|---|---|---|---|---|---|---|
{chr(10).join(f"| {m} | {x['state_validation']['n']} | {x['state_validation']['mae_pp']:.3f} | {x['state_validation']['max_abs_pp']:.2f} | {x['state_validation']['pearson']:.3f} | {100 * x['state_validation']['coverage']:.0f}% | {x['state_validation']['sd_direct']:.2f} | {x['state_validation']['sd_sae']:.2f} |" for m, x in s['results'].items())}

Largest differences (primary): {'; '.join(f"{x.state_abbr} SAE {x.sae:.2f} vs direct {x.value:.2f}" for x in pd.concat([st.dropna(subset=['value']).head(3), st.dropna(subset=['value']).tail(3)]).itertuples())}.
Table: `results/tables/brfss_sae_state_validation_<measure>.csv`.

## Validation (b): leave-one-state-out prediction of state prevalence (primary outcome)
Each state's respondents are removed, the model is refitted (sigma re-estimated), and the state's prevalence is
predicted from its counties' composition, urbanicity and covariates with a state effect drawn from N(0, sigma^2).
Baseline: the weighted national prevalence of the other states.

| model | states | MAE (pp) | RMSE (pp) | Pearson r (pred vs direct) | 95% interval coverage | mean bias (pp) | baseline MAE | baseline RMSE |
|---|---|---|---|---|---|---|---|---|
| M0 composition + urbanicity | {L0['n_states']} | {L0['mae_model']:.3f} | {L0['rmse_model']:.3f} | {L0['pearson_model_vs_direct']:.3f} | {100 * L0['coverage_95']:.0f}% | {L0['mean_bias']:.3f} | {L0['mae_baseline']:.3f} | {L0['rmse_baseline']:.3f} |
| M1 + ACS covariates | {L1['n_states']} | {L1['mae_model']:.3f} | {L1['rmse_model']:.3f} | {L1['pearson_model_vs_direct']:.3f} | {100 * L1['coverage_95']:.0f}% | {L1['mean_bias']:.3f} | {L1['mae_baseline']:.3f} | {L1['rmse_baseline']:.3f} |

Largest LOSO errors ({s['chosen']}): {worst_loso}. Table: `results/tables/brfss_sae_loso.csv`.
Reading: the LOSO error is the error to expect for a state the survey did not sample (Kentucky, Pennsylvania), and
an upper bound on how well composition + urbanicity + covariates explain between-area differences.

## Validation (c): sub-state check against SMART BRFSS MMSA direct estimates (primary outcome)
The 2023 SMART MMSA file gives direct current-Long-COVID estimates for {mm['n_mmsa_total']} metropolitan areas/divisions; {mm['n_mmsa_matched']}
map to 2024-vintage counties through the OMB March 2020 delineation (unmatched: {', '.join(x['mmsa_name'] for x in mm['unmatched']) or 'none'}).
The model never sees MMSA membership, so the within-state contrast (MMSA vs its state) is out-of-sample geographically,
although the MMSA respondents are part of the fitted sample.

* MAE vs MMSA direct: SAE {mm['mae_sae']:.3f} pp; inherited state direct value {mm['mae_inherited_state']:.3f} pp; SAE with covariates at
  their state x urbanicity mean (composition + urbanicity only) {mm['mae_sae_composition_only']:.3f} pp.
* Within-state deviation (MMSA minus its state): observed SD {mm['sd_dev_direct']:.3f} pp, of which sampling noise accounts for
  RMS SE {mm['mean_se_direct']:.3f} pp, leaving a signal SD of about {mm['signal_sd_dev_direct']:.3f} pp; the SAE predicts deviations with SD
  {mm['sd_dev_sae']:.3f} pp. Correlation of predicted with observed deviation: r = {mm['pearson_dev_sae_vs_direct']:.3f} (composition-only
  variant r = {mm['pearson_dev_comp_vs_direct']:.3f}); regression slope of observed on predicted deviation {mm['slope_dev_direct_on_dev_sae']:.2f}.
  Mean deviation observed {mm['mean_dev_direct']:.3f} pp vs predicted {mm['mean_dev_sae']:.3f} pp.
* Table: `results/tables/brfss_sae_mmsa_validation.csv`.

## How much within-state variation the SAE shows, and why
* Across all counties the SD of the estimate is {ws['county_sd_pp']:.3f} pp; the within-state SD is {ws['within_state_sd_pp']:.3f} pp
  ({100 * ws['share_variance_within_state']:.0f}% of the county variance is within states). Median within-state range {ws['median_state_range_pp']:.2f} pp
  (max {ws['max_state_range_pp']:.2f} pp); median population-weighted within-state SD {ws['median_within_state_sd_pp']:.3f} pp.
  Widest ranges: {top_range}.
* With the covariates held at their state x urbanicity mean (composition + urbanicity only) the within-state SD is
  {wc['within_state_sd_pp']:.3f} pp ({100 * wc['share_variance_within_state']:.0f}% of variance within states).
* Where it comes from: the state effect is common to all counties and no respondent's county is known, so every
  within-state difference is produced by composition (age x sex, race/ethnicity, education), urbanicity and the
  poverty/insurance term; about half of the within-state SD comes from the covariate term (compare the two SDs above).
* How far to trust it: the covariate term improves the between-state LOSO prediction (RMSE {L1['rmse_model']:.3f} vs {L0['rmse_model']:.3f} pp
  without it), but the sub-state MMSA check is weak: predicted and observed MMSA-vs-state deviations correlate at
  r = {mm['pearson_dev_sae_vs_direct']:.3f} (P = {mm['p_dev_sae_vs_direct']:.2f}, n = {mm['n_mmsa_matched']}), and the SAE is no closer to MMSA direct estimates than the inherited state value
  (MAE {mm['mae_sae']:.3f} vs {mm['mae_inherited_state']:.3f} pp). The MMSA data cannot confirm the within-state pattern (their sampling noise is
  {mm['mean_se_direct'] / max(mm['signal_sd_dev_direct'], 1e-9):.1f} times the signal SD), and they do not contradict it. The within-state ordering of counties is therefore an
  ordering by composition, urbanicity and poverty/insurance, not by observed Long COVID. Use the county values as
  population-scaled expectations (adults x modelled rate) whose state totals are validated, not as evidence that one
  county of a state has more Long COVID than its neighbour.
* LOSO 95% intervals covered {100 * L1['coverage_95']:.0f}% (M1) of held-out states: intervals for states without data (Kentucky,
  Pennsylvania) are somewhat too narrow (sigma is fixed at its estimate).

## Use downstream
* `numerator` = modelled adults with the outcome (value x ACS adults) with `numerator_ci_*`; `burden_evidence_level = A`
  refers to the measure (a direct self-reported Long COVID item); the modelled nature is carried by
  `evidence_type = modeled_small_area_estimate`, `estimate_kind = modeled`, `within_state_variation_basis` and
  `state_effect_source`.
* Draws for subgroup burden: `data/interim/cdc_brfss/sae_draws_<measure>.npz` (county x age7 x sex, 500 draws).
* `measure_it.geography.subgroup` (table `geo_subgroup_burden`) multiplies these draws, cell by cell, by a device-defined
  subgroup's cohort fractions (`subgroup_strata_fractions`, pooled fraction for suppressed cells). For `lyme_disease` it
  uses 2023 CDC county reported cases (all ages, 2022 case definition, under-ascertained) and for `ptlds` a labelled
  derived estimate (reported cases x 13.7% PTLD among ideally treated early Lyme, Aucott et al. 2022, PMID 35066160);
  `burden_basis` names the base in every row.
* Reproduce: `uv run python -m measure_it.ingestion.brfss && uv run python -m measure_it.geography.sae`.
"""
    p = RESULTS / "BRFSS_SAE_RESULTS.md"
    p.write_text(text)
    return p


def write_sae_audit(s: dict) -> Path:
    r = s["results"][PRIMARY]
    L = s["loso"][s["chosen"]]
    sec = f"""{SAE_START}
## Small-area estimate built on this source (`measure_it.geography.sae`, {s['run_at']})
* Output `geo_condition_burden__brfss_long_covid_sae`: county (3,144), state-aggregate and national rows for
  `lc_current_pct_all_adults` (primary), `lc_significant_limitation_pct_all_adults`, `lc_current_pct_all_adults_mmwr_def`;
  `evidence_type = modeled_small_area_estimate`, `source_geographic_resolution = county` for county rows.
* Method: MRP, random state intercept, age x sex + race/ethnicity + education + urbanicity (+ ACS poverty/uninsured if
  selected; chosen: {s['chosen']}); post-stratified to ACS 2020-2024 county cells. Full description and every number:
  `results/BRFSS_SAE_RESULTS.md`.
* burden_evidence_level = **A**: the outcome is the direct self-reported Long COVID item (same as the HPS level-A
  measure); the county value is modelled, which is carried by `evidence_type`, `estimate_kind = modeled`,
  `within_state_variation_basis` and `state_effect_source` rather than by downgrading the measure level. Downstream
  users who rank counties WITHIN a state should know that the within-state variation is composition/urbanicity/
  covariate-driven only.
* Validation: state aggregate vs direct MAE {r['state_validation']['mae_pp']:.3f} pp (in-sample); leave-one-state-out MAE
  {L['mae_model']:.3f} pp (RMSE {L['rmse_model']:.3f}; national-mean baseline MAE {L['mae_baseline']:.3f}); MMSA sub-state check
  r = {r['mmsa']['pearson_dev_sae_vs_direct']:.3f} between predicted and observed within-state deviations.
* Kentucky and Pennsylvania (absent from the 2023 file): state effect drawn from N(0, sigma^2).
{SAE_END}
"""
    p = raw_dir(B.SOURCE_ID) / "DATA_AUDIT.md"
    text = p.read_text() if p.exists() else ""
    if SAE_START in text:
        text = text.split(SAE_START)[0] + sec + text.split(SAE_END, 1)[1].lstrip("\n")
    else:
        text = text.rstrip("\n") + "\n\n" + sec
    p.write_text(text)
    return p


if __name__ == "__main__":
    run()
