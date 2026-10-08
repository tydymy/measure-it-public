"""Design-based estimation for complex surveys (NHANES): weighted domain proportions and means with
Taylor-linearised standard errors, and the NCHS reliability rules for proportions.

Estimator (domain mean / proportion, ratio form):

    p_D = sum_i d_i w_i y_i / sum_i d_i w_i

Linearised variable (all sample units, zero outside the domain or where y is missing):

    z_i = d_i w_i (y_i - p_D) / sum_i d_i w_i

With-replacement PSU variance within strata (the NHANES masked-variance design, SDMVSTRA / SDMVPSU):

    v(p_D) = sum_h n_h / (n_h - 1) * sum_j (z_hj - zbar_h)^2,   z_hj = sum of z_i over PSU j in stratum h

This is the estimator of R `survey::svymean(~y, subset(design, D))` and SUDAAN's DESCRIPT with SUBPOPN (domain
estimation keeps the full design; the domain is never subset before the variance is computed). Strata with a single
PSU contribute 0 by default (`lonely_psu="remove"`, as R `options(survey.lonely.psu="remove")`), or are centred on
the grand mean of PSU totals (`"adjust"`).

Confidence intervals: Wald on the Taylor SE with a t quantile on the design degrees of freedom (number of PSUs minus
number of strata among those holding domain units), clipped to [0, 1] for proportions; plus the Korn-Graubard interval
used by the NCHS Data Presentation Standards for Proportions (Parker et al. 2017, Vital Health Stat 2(175)), from
which the `nchs_reliability` flag is derived.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class Estimate:
    estimate: float
    se: float
    ci_low: float
    ci_high: float
    df: int
    n_unweighted: int          # domain units with a non-missing y and a positive weight
    sum_weights: float         # estimated population size of the domain
    kg_ci_low: float = np.nan  # Korn-Graubard (proportions only)
    kg_ci_high: float = np.nan
    effective_n: float = np.nan
    nchs_reliability: str = ""  # "reliable", "suppress: ...", "review: ..." (proportions only)

    def as_dict(self) -> dict:
        return asdict(self)


def _psu_totals(z: np.ndarray, strata: np.ndarray, psu: np.ndarray) -> pd.DataFrame:
    df = pd.DataFrame({"h": strata, "j": psu, "z": z})
    return df.groupby(["h", "j"], sort=True, observed=True)["z"].sum().reset_index()


def taylor_variance(z: np.ndarray, strata: np.ndarray, psu: np.ndarray, lonely_psu: str = "remove") -> float:
    """Variance of a total of linearised values z under the stratified with-replacement PSU design."""
    t = _psu_totals(z, strata, psu)
    nh = t.groupby("h")["j"].transform("size").to_numpy()
    zbar = t.groupby("h")["z"].transform("mean").to_numpy()
    zz = t["z"].to_numpy()
    multi = nh > 1
    v = float(np.sum((nh[multi] / (nh[multi] - 1.0)) * (zz[multi] - zbar[multi]) ** 2))
    if lonely_psu == "adjust" and (~multi).any():
        grand = zz.mean()
        v += float(np.sum((zz[~multi] - grand) ** 2))
    elif lonely_psu not in ("remove", "adjust"):
        raise ValueError("lonely_psu must be 'remove' or 'adjust'")
    return v


def design_df(strata: np.ndarray, psu: np.ndarray, mask: np.ndarray | None = None) -> int:
    """Number of PSUs minus number of strata (restricted to strata/PSUs holding units in `mask` if given)."""
    s, p = np.asarray(strata), np.asarray(psu)
    if mask is not None:
        s, p = s[mask], p[mask]
    if len(s) == 0:
        return 0
    pairs = pd.DataFrame({"h": s, "j": p}).drop_duplicates()
    return int(len(pairs) - pairs["h"].nunique())


def korn_graubard(p: float, se: float, n: int, df: int, level: float = 0.95) -> tuple[float, float, float]:
    """Korn-Graubard CI for a survey proportion (Korn & Graubard 1998; NCHS 2017 standards).

    Returns (low, high, effective_n). effective_n = p(1-p)/var, capped at n (n when var == 0), then multiplied by
    (t_{n-1} / t_{df})^2 for the interval."""
    a = 1.0 - level
    if n <= 0 or not np.isfinite(p):
        return np.nan, np.nan, np.nan
    var = se ** 2
    n_eff = n if (var <= 0 or p in (0.0, 1.0)) else min(n, p * (1 - p) / var)
    if df > 0 and n > 1:
        n_eff_adj = n_eff * (stats.t.ppf(1 - a / 2, n - 1) / stats.t.ppf(1 - a / 2, df)) ** 2
    else:
        n_eff_adj = n_eff
    x = p * n_eff_adj
    low = 0.0 if x <= 0 else float(stats.beta.ppf(a / 2, x, n_eff_adj - x + 1))
    high = 1.0 if x >= n_eff_adj else float(stats.beta.ppf(1 - a / 2, x + 1, n_eff_adj - x))
    return low, high, float(n_eff)


def nchs_reliability(p: float, kg_low: float, kg_high: float, n_eff: float, df: int) -> str:
    """NCHS Data Presentation Standards for Proportions (2017), simplified to the published decision rules."""
    if not np.isfinite(p):
        return "suppress: no estimate"
    if not np.isfinite(n_eff) or n_eff < 30:
        return "suppress: effective sample size < 30"
    width = kg_high - kg_low
    if width >= 0.30:
        return "suppress: Korn-Graubard CI width >= 0.30"
    if 0.05 < width < 0.30 and p > 0 and width / p > 1.30:
        return "suppress: relative CI width > 130%"
    if df < 8:
        return "review: design degrees of freedom < 8"
    if p == 0 or p == 1:
        return "review: estimate of 0 or 1"
    return "reliable"


def domain_estimate(y, w, strata, psu, domain=None, *, proportion: bool = True, level: float = 0.95,
                    lonely_psu: str = "remove") -> Estimate:
    """Weighted mean (or proportion) of y in a domain, with a Taylor-linearised SE over the FULL design.

    y, w, strata, psu, domain are aligned arrays over every sample unit in the design (do not pre-subset to the
    domain). Units with missing y or non-positive/missing weight are treated as outside the domain."""
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    strata = np.asarray(strata)
    psu = np.asarray(psu)
    d = np.ones(len(y), dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    d = d & np.isfinite(y) & np.isfinite(w) & (w > 0)
    n = int(d.sum())
    W = float(w[d].sum())
    if n == 0 or W <= 0:
        return Estimate(np.nan, np.nan, np.nan, np.nan, 0, n, W, nchs_reliability=("suppress: no units"
                                                                                   if proportion else ""))
    est = float(np.sum(w[d] * y[d]) / W)
    z = np.zeros(len(y))
    z[d] = w[d] * (y[d] - est) / W
    var = taylor_variance(z, strata, psu, lonely_psu=lonely_psu)
    se = float(np.sqrt(max(var, 0.0)))
    df = design_df(strata, psu, d)
    q = stats.t.ppf(1 - (1 - level) / 2, df) if df > 0 else np.inf
    lo, hi = est - q * se, est + q * se
    if proportion:
        lo, hi = max(0.0, lo), min(1.0, hi)
        kg_lo, kg_hi, n_eff = korn_graubard(est, se, n, df, level)
        rel = nchs_reliability(est, kg_lo, kg_hi, n_eff, df)
        return Estimate(est, se, float(lo), float(hi), df, n, W, kg_lo, kg_hi, n_eff, rel)
    return Estimate(est, se, float(lo), float(hi), df, n, W)


def weighted_quantiles(x, w, qs) -> np.ndarray:
    """Weighted quantiles (type: inverse of the weighted empirical CDF at the midpoint of each weight step)."""
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[ok], w[ok]
    if len(x) == 0:
        return np.full(len(qs), np.nan)
    o = np.argsort(x)
    x, w = x[o], w[o]
    cdf = (np.cumsum(w) - 0.5 * w) / w.sum()
    return np.interp(qs, cdf, x)
