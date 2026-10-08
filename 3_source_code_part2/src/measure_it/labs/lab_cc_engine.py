"""Shared case-control statistics for the per-dataset lab analyses in measure_it.labs (lab biomarker vs diagnosis).

Used by heds_hsd_olink_serum_cinquina2026, klein2023_mylc_ml_table and endo_arg1_repod. Every function is
deterministic given SEED. The conventions follow measure_it.wearables.fm_thermography_analysis:

* cross-validated AUROC = AUROC of the pooled out-of-fold predictions of one repeat, averaged over repeats;
  folds depend only on (y, seed, repeat), so models compared on the same people are paired;
* bootstrap CIs resample PEOPLE within each class, holding the out-of-fold predictions fixed (so they include
  person sampling, not refit variability);
* permutation nulls shuffle the labels and re-run the whole cross-validation per shuffle;
* fit-free single-analyte AUROC has the direction fixed in advance (no fitting, so no CV is needed).
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import stats

from measure_it.config import SEED
from measure_it.wearables.stanford_models import auroc_rows, bh_qvalues, hedges_g, percentile_ci

N_SPLITS, N_REPEATS, N_BOOT = 5, 20, 2000


def make_model(C: float = 1.0):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         LogisticRegression(C=C, max_iter=20000))


def folds(y: np.ndarray, repeat: int, seed: int = SEED) -> list:
    from sklearn.model_selection import StratifiedKFold
    return list(StratifiedKFold(N_SPLITS, shuffle=True, random_state=seed + repeat).split(np.zeros((len(y), 1)), y))


def cv_oof(X: np.ndarray, y: np.ndarray, C: float = 1.0, n_repeats: int = N_REPEATS, seed: int = SEED) -> np.ndarray:
    """Out-of-fold P(case), shape (n_repeats, n). Imputation and scaling are fitted inside each training fold."""
    X, y = np.asarray(X, float), np.asarray(y).astype(int)
    if X.ndim == 1:
        X = X[:, None]
    P = np.full((n_repeats, len(y)), np.nan)
    for r in range(n_repeats):
        for tr, te in folds(y, r, seed):
            m = make_model(C)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m.fit(X[tr], y[tr])
            P[r, te] = m.predict_proba(X[te])[:, 1]
    return P


def stratified_boot_idx(y: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED, strata: np.ndarray | None = None):
    """(n_boot, n) person indices resampled with replacement within each class (and within each stratum if given)."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    keys = y if strata is None else np.array([f"{a}|{b}" for a, b in zip(y, strata)])
    blocks = [np.where(keys == k)[0] for k in pd.unique(keys)]
    return np.hstack([rng.choice(b, (n_boot, b.size)) for b in blocks])


def sens_at_spec(y: np.ndarray, P: np.ndarray, spec: float = 0.90) -> np.ndarray:
    y = np.asarray(y).astype(bool)
    P = np.atleast_2d(P)
    thr = np.quantile(P[:, ~y], spec, axis=1, method="higher")
    return (P[:, y] > thr[:, None]).mean(axis=1)


def summarize(y: np.ndarray, P: np.ndarray, boot_idx: np.ndarray) -> tuple[dict, np.ndarray]:
    """AUROC (mean over repeats), CI from the person bootstrap, sensitivity at 90% specificity."""
    y = np.asarray(y).astype(int)
    P = np.atleast_2d(P)
    auc = auroc_rows(y, P)
    sas = sens_at_spec(y, P)
    b_auc = np.array([auroc_rows(y[i], P[:, i]).mean() for i in boot_idx])
    b_sas = np.array([sens_at_spec(y[i], P[:, i]).mean() for i in boot_idx])
    lo, hi = percentile_ci(b_auc)
    slo, shi = percentile_ci(b_sas)
    return ({"auroc": float(auc.mean()), "auroc_ci_low": lo, "auroc_ci_high": hi,
             "auroc_repeat_min": float(auc.min()), "auroc_repeat_max": float(auc.max()),
             "sens_at_90spec": float(sas.mean()), "sens_at_90spec_ci_low": slo, "sens_at_90spec_ci_high": shi,
             "n": int(len(y)), "n_cases": int(y.sum()), "n_controls": int((1 - y).sum())}, b_auc)


def _perm_task(X, y, C, i, seed):
    yp = np.random.default_rng([seed, i]).permutation(y)
    return float(auroc_rows(yp, cv_oof(X, yp, C)).mean())


def permutation_null(X: np.ndarray, y: np.ndarray, C: float = 1.0, n_perm: int = 1000, jobs: int = 8,
                     seed: int = SEED) -> np.ndarray:
    from joblib import Parallel, delayed
    return np.array(Parallel(n_jobs=jobs)(delayed(_perm_task)(X, y, C, i, seed) for i in range(n_perm)))


def perm_p(null: np.ndarray, observed: float) -> float:
    return float((1 + np.sum(null >= observed - 1e-12)) / (1 + len(null)))


def fitfree_auroc(x_case: np.ndarray, x_ctrl: np.ndarray) -> float:
    """P(case > control) + 0.5 P(tie)."""
    x1, x0 = np.asarray(x_case, float), np.asarray(x_ctrl, float)
    return float(stats.mannwhitneyu(x1, x0, alternative="two-sided").statistic / (x1.size * x0.size))


def fitfree_test(x: np.ndarray, y: np.ndarray, n_boot: int = N_BOOT, n_perm: int = 10000, seed: int = SEED,
                 higher_is_case: bool = True, strata: np.ndarray | None = None) -> dict:
    """Fit-free single-analyte AUROC with the direction fixed in advance; stratified person-bootstrap CI and a
    one-sided label-permutation p (larger AUROC = more extreme). With strata, the AUROC is the pair-weighted mean of
    within-stratum AUROCs (strata lacking either class are dropped) and both resampling and shuffling stay within
    stratum."""
    x, y = np.asarray(x, float), np.asarray(y).astype(int)
    s = -x if not higher_is_case else x

    def stat(yy, ss, st):
        if st is None:
            return auroc_rows(yy, ss[None, :])[0]
        num = den = 0.0
        for k in pd.unique(st):
            m = st == k
            n1, n0 = int(yy[m].sum()), int((1 - yy[m]).sum())
            if n1 and n0:
                num += auroc_rows(yy[m], ss[m][None, :])[0] * n1 * n0
                den += n1 * n0
        return num / den if den else np.nan

    obs = stat(y, s, strata)
    bidx = stratified_boot_idx(y, n_boot, seed, strata)
    boots = np.array([stat(y[i], s[i], None if strata is None else strata[i]) for i in bidx])
    rng = np.random.default_rng(seed + 1)
    null = np.empty(n_perm)
    for j in range(n_perm):
        if strata is None:
            yp = rng.permutation(y)
        else:
            yp = y.copy()
            for k in pd.unique(strata):
                m = np.where(strata == k)[0]
                yp[m] = rng.permutation(y[m])
        null[j] = stat(yp, s, strata)
    lo, hi = percentile_ci(boots)
    ctrl = s[y == 0]
    thr = np.quantile(ctrl, 0.90, method="higher")
    return {"auroc": float(obs), "auroc_ci_low": lo, "auroc_ci_high": hi, "perm_p": perm_p(null, obs),
            "n_perm": n_perm, "null_q95": float(np.quantile(null, 0.95)), "n_cases": int(y.sum()),
            "n_controls": int((1 - y).sum()), "sens_at_90spec": float((s[y == 1] > thr).mean()),
            "direction": "higher = case" if higher_is_case else "lower = case"}


def per_analyte_effects(df: pd.DataFrame, label: str, analytes: list[str], transform=None,
                        n_boot: int = N_BOOT, seed: int = SEED) -> pd.DataFrame:
    """Per analyte: n, medians, Hedges g (case minus control, on transform(values)), stratified bootstrap CI,
    fit-free AUROC P(case > control), two-sided Mann-Whitney p; BH q across the analytes given."""
    rng = np.random.default_rng(seed)
    rows = []
    for a in analytes:
        v = pd.to_numeric(df[a], errors="coerce")
        ok = v.notna()
        x1 = v[ok & (df[label] == 1)].to_numpy(float)
        x0 = v[ok & (df[label] == 0)].to_numpy(float)
        if x1.size < 3 or x0.size < 3:
            rows.append({"analyte": a, "n_cases": x1.size, "n_controls": x0.size})
            continue
        t1, t0 = (transform(x1), transform(x0)) if transform else (x1, x0)
        g = hedges_g(t1, t0)
        gb = [hedges_g(rng.choice(t1, t1.size), rng.choice(t0, t0.size)) for _ in range(n_boot)]
        lo, hi = percentile_ci(np.array(gb))
        u = stats.mannwhitneyu(x1, x0, alternative="two-sided")
        rows.append({"analyte": a, "n_cases": x1.size, "n_controls": x0.size, "median_cases": float(np.median(x1)),
                     "median_controls": float(np.median(x0)), "hedges_g": g, "g_ci_low": lo, "g_ci_high": hi,
                     "auroc_case_higher": float(u.statistic / (x1.size * x0.size)), "p_value": float(u.pvalue)})
    out = pd.DataFrame(rows)
    out["q_value"] = bh_qvalues(out["p_value"]) if "p_value" in out else np.nan
    return out


def cv_block(X: np.ndarray, y: np.ndarray, C: float, boot_idx: np.ndarray, n_perm: int = 0, jobs: int = 8) -> dict:
    """CV AUROC + bootstrap CI (+ permutation null when n_perm > 0). Returns summary, OOF matrix, boots, null."""
    P = cv_oof(X, y, C)
    summ, boots = summarize(y, P, boot_idx)
    null = permutation_null(X, y, C, n_perm, jobs) if n_perm else np.array([])
    if n_perm:
        summ.update({"perm_p": perm_p(null, summ["auroc"]), "n_perm": n_perm, "null_mean": float(null.mean()),
                     "null_q95": float(np.quantile(null, 0.95))})
    return {"summary": summ, "P": P, "boots": boots, "null": null}


def delta(a: dict, b: dict, y: np.ndarray) -> dict:
    """Paired Delta AUROC (a minus b) with the paired bootstrap CI (both blocks must share boot_idx and folds)."""
    d = auroc_rows(y, a["P"]) - auroc_rows(y, b["P"])
    db = a["boots"] - b["boots"]
    lo, hi = percentile_ci(db)
    return {"delta_auroc": float(d.mean()), "ci_low": lo, "ci_high": hi, "delta_repeat_min": float(d.min()),
            "delta_repeat_max": float(d.max()), "bootstrap_fraction_le_0": float(np.mean(db <= 0))}


def pfmt(p: float, n_perm: int) -> str:
    """Permutation p for reports: at the floor 1/(n+1) say so instead of printing a spuriously exact value."""
    if p is None or not np.isfinite(p):
        return "not run"
    floor = 1 / (n_perm + 1)
    return f"<= {floor:.4f} (no shuffle reached the observed value)" if p <= floor + 1e-12 else f"{p:.4f}"
