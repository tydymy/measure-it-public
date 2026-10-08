"""Person-level modelling of the mapMECFS participant-linked omics (NIH PI-ME/CFS deep phenotyping deposits).

Pre-registered plan: ``docs/ANALYSIS_PLAN_LINKED_OMICS.md`` (written before any of this was run). This is the one
dataset in the project where several omics layers belong to the SAME people (NIH study number). None of it links to
wearable, actigraphy, HRV, orthostatic or CPET data (0 identifier pairs in open files); the cohort is small,
single-site, and adjudication status is not in the open files.

Per layer (csf_somalogic, serum_somalogic, pbmc_rnaseq, muscle_rnaseq, csf_metabolomics):
  * penalised logistic classifiers (primary L2, secondary elastic net) with ALL preprocessing and top-k (k = 50)
    Welch-t feature selection inside the training folds; 10 x 5 repeated stratified CV; AUROC with stratified
    participant-bootstrap CI; 1,000-permutation label null re-running the whole procedure; age/sex baseline;
  * sensitivity S1 (features residualised on sex/age inside the fold) and S2 (k = 10, 200);
  * late fusion (mean per-layer log-odds) on the people who have every fused layer;
  * univariate Welch t / Hedges' g with BH-FDR, OLS adjusted for sex/age, and agreement with the paper's published
    directions (a reproduction check on overlapping people, not replication);
  * CSF vs serum SomaScan consistency (same 42 people) with a joint label-permutation null.

Outputs: results/tables/linked_omics_*.csv, results/figures/drafts/linked_omics_*.png,
results/LINKED_OMICS_RESULTS.md, processed partition ``phenotype_signatures__mapmecfs_omics``.

Run: ``uv run python -m measure_it.omics.person_linked [--workers 16] [--n-perm 1000]``;
``--report-only`` rebuilds the markdown report, figures and signature partition from the saved tables.
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .. import config
from ..config import UNKNOWN

PRODUCER = "measure_it.omics.person_linked"
ANALYSIS_VERSION = "2026-09-24.1"
DATASET_ID = "mapmecfs_nih_pi_mecfs"
SIGNATURE_TABLE = "phenotype_signatures__mapmecfs_omics"
PLAN = "docs/ANALYSIS_PLAN_LINKED_OMICS.md"
LAYERS = ("csf_somalogic", "serum_somalogic", "pbmc_rnaseq", "muscle_rnaseq", "csf_metabolomics")
RNA_LAYERS = frozenset({"pbmc_rnaseq", "muscle_rnaseq"})
LAYER_LABEL = {
    "csf_somalogic": "CSF SomaScan 1.3k proteomics",
    "serum_somalogic": "serum SomaScan 1.3k proteomics",
    "pbmc_rnaseq": "PBMC bulk RNA-seq",
    "muscle_rnaseq": "skeletal-muscle bulk RNA-seq",
    "csf_metabolomics": "CSF untargeted metabolomics (Metabolon)",
}
LAYER_UNITS = {
    "csf_somalogic": "log2 RFU", "serum_somalogic": "log2 RFU",
    "pbmc_rnaseq": "log2(CPM + 1)", "muscle_rnaseq": "log2(CPM + 1)",
    "csf_metabolomics": "log2 batch-normalised relative abundance",
}
# which covariates exist per layer: CSF metabolomics carries only Birth Sex = 1 rows, so sex is constant there
LAYER_COVARIATES = {lay: ("female", "age") for lay in LAYERS} | {"csf_metabolomics": ("age",)}

K_PRIMARY = 50
K_SENSITIVITY = (10, 200)
N_REPEATS = 10
N_FOLDS = 5
N_PERM = 1000
N_BOOT = 2000
FDR_ALPHA = 0.05
PASS_ALPHA = 0.05 / len(LAYERS)          # Bonferroni over the 5 layers, primary L2 model
EXPRESSED_CPM = 1.0
EXPRESSED_FRACTION = 0.5
MIN_CLASS_FOR_CV = 6

FUSIONS = {
    "F1": ("csf_somalogic", "serum_somalogic"),
    "F2": ("csf_somalogic", "serum_somalogic", "pbmc_rnaseq"),
    "F3": ("csf_somalogic", "serum_somalogic", "muscle_rnaseq"),
    "F4": ("pbmc_rnaseq", "muscle_rnaseq"),
    "F5": LAYERS,
}
FUSION_ROLE = {"F1": "primary fusion", "F2": "secondary", "F3": "secondary", "F4": "secondary",
               "F5": "n only (pre-specified: not analysed)"}

LABEL_BASIS = "NIH intramural adjudicated PI-ME/CFS vs healthy volunteer (per paper)"
PAPER_TABLE = {"pbmc_rnaseq": "16A", "muscle_rnaseq": "19A", "serum_somalogic": "17A", "csf_somalogic": "17B",
               "csf_metabolomics": "14A"}

COHORT_CAVEAT = (
    "Small single-site case-control cohort (NIH Clinical Center protocol 16-N-0058). Adjudication status is not in "
    "the open files: the deposits hold 21 PI-ME/CFS + 21 HV SomaScan samples while the paper's analytic cohort is "
    "17 PI-ME/CFS + 21 HV, so some deposited PI-ME/CFS participants may be outside it. No wearable, actigraphy, HRV, "
    "orthostatic or CPET measurement is linked to these people in open data (0 identifier pairs): this is omics-only "
    "person-level signal, not a device evaluation and not patient multi-omics with wearables.")

MODELS = {
    "l2_top50": {"kind": "omics", "clf": "l2", "k": K_PRIMARY, "resid": False,
                 "label": "L2 logistic, top-50 Welch-t features (primary)"},
    "enet_top50": {"kind": "omics", "clf": "enet", "k": K_PRIMARY, "resid": False,
                   "label": "elastic-net logistic (l1_ratio 0.5), top-50 (secondary)"},
    "age_sex": {"kind": "baseline", "clf": "l2", "k": None, "resid": False,
                "label": "age/sex-only L2 logistic (baseline)"},
    "l2_top50_resid": {"kind": "omics", "clf": "l2", "k": K_PRIMARY, "resid": True,
                       "label": "S1: L2 top-50 on sex/age-residualised features"},
    "l2_top10": {"kind": "omics", "clf": "l2", "k": 10, "resid": False, "label": "S2: L2 top-10"},
    "l2_top200": {"kind": "omics", "clf": "l2", "k": 200, "resid": False, "label": "S2: L2 top-200"},
}
PERMUTED_MODELS = ("l2_top50", "enet_top50", "age_sex", "l2_top50_resid")


# ============================================================================================ data


@dataclass
class Layer:
    name: str
    ids: np.ndarray            # participant ids (sorted)
    y: np.ndarray              # 1 = PI-ME/CFS, 0 = HV
    X: np.ndarray              # log-scale values (people x features)
    features: np.ndarray       # feature ids
    feature_names: np.ndarray
    cov: np.ndarray            # covariates (people x c), NaN = missing
    cov_names: tuple
    expressed: np.ndarray | None = None   # RNA-seq: CPM >= 1 (per-sample quantity)

    def subset(self, ids) -> "Layer":
        pos = pd.Index(self.ids).get_indexer(ids)
        if (pos < 0).any():
            raise KeyError(f"{self.name}: ids not in layer")
        return Layer(self.name, self.ids[pos], self.y[pos], self.X[pos], self.features, self.feature_names,
                     self.cov[pos], self.cov_names, None if self.expressed is None else self.expressed[pos])


def _sex_for_layer(sex_by_source: str, layer: str) -> float:
    try:
        d = json.loads(sex_by_source) if isinstance(sex_by_source, str) else {}
    except ValueError:
        d = {}
    s = d.get(layer)
    return {"Female": 1.0, "Male": 0.0}.get(s, np.nan)


def log_transform(layer: str, wide: pd.DataFrame) -> tuple[np.ndarray, np.ndarray | None]:
    """Per-sample log transform (no cross-sample information, so safe outside the CV folds).
    RNA-seq: log2(CPM + 1) with each sample's own library size, plus the CPM >= 1 indicator for the in-fold filter."""
    v = wide.to_numpy(dtype=float)
    if layer in RNA_LAYERS:
        cpm = v / v.sum(axis=1, keepdims=True) * 1e6
        return np.log2(cpm + 1.0), cpm >= EXPRESSED_CPM
    if (v <= 0).any():
        raise ValueError(f"{layer}: non-positive values; the plan's log2 transform assumes none")
    return np.log2(v), None


def load_layers() -> dict[str, Layer]:
    from ..store import read_table

    omics = read_table("participant_omics_linked__mapmecfs",
                       columns=["participant_id", "modality", "feature_id", "feature_name", "value"])
    people = read_table("participants__mapmecfs").set_index("participant_id")
    y_map = {"PI-ME/CFS": 1, "HV": 0}
    out = {}
    for lay in LAYERS:
        s = omics[omics["modality"] == lay].drop_duplicates(["participant_id", "feature_id"], keep="first")
        wide = s.pivot(index="participant_id", columns="feature_id", values="value").sort_index().sort_index(axis=1)
        if wide.isna().any().any():
            raise ValueError(f"{lay}: missing cells in the linked omics matrix")
        names = s.drop_duplicates("feature_id").set_index("feature_id")["feature_name"].reindex(wide.columns)
        X, expressed = log_transform(lay, wide)
        ids = wide.index.to_numpy()
        pp = people.loc[ids]
        y = pp["group_label"].map(y_map)
        if y.isna().any():
            raise ValueError(f"{lay}: unlabelled participants {list(y[y.isna()].index)}")
        cols = []
        for c in LAYER_COVARIATES[lay]:
            if c == "female":
                cols.append([_sex_for_layer(v, lay) for v in pp["sex_by_source"]])
            else:
                cols.append(pd.to_numeric(pp["age_years_somalogic"], errors="coerce").to_list())
        out[lay] = Layer(lay, ids, y.to_numpy(int), X, wide.columns.to_numpy(), names.fillna("").to_numpy(),
                         np.asarray(cols, dtype=float).T, LAYER_COVARIATES[lay], expressed)
    return out


# ============================================================================================ statistics helpers


def welch_t(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Vectorised Welch t statistic per column (group 1 minus group 0); 0 where both variances are 0."""
    a, b = X[y == 1], X[y == 0]
    va, vb = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
    se = np.sqrt(va / len(a) + vb / len(b))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (a.mean(axis=0) - b.mean(axis=0)) / se
    return np.where(se > 0, t, 0.0)


def auroc(y: np.ndarray, s: np.ndarray) -> float:
    """Mann-Whitney AUROC with mid-ranks for ties."""
    y = np.asarray(y)
    n1 = int(y.sum())
    n0 = len(y) - n1
    r = stats.rankdata(s)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def auroc_rows(y: np.ndarray, S: np.ndarray) -> np.ndarray:
    """AUROC of every row of S (repeats x people) against the same labels."""
    n1 = int(y.sum())
    n0 = len(y) - n1
    R = stats.rankdata(S, axis=1)
    return (R[:, y == 1].sum(axis=1) - n1 * (n1 + 1) / 2) / (n1 * n0)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        from scipy.stats import false_discovery_control
        q[ok] = false_discovery_control(p[ok], method="bh")
    return q


def _impute(Ctr: np.ndarray, Cte: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    med = np.nanmedian(Ctr, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    return np.where(np.isnan(Ctr), med, Ctr), np.where(np.isnan(Cte), med, Cte)


def _make_clf(kind: str):
    from sklearn.linear_model import LogisticRegression
    if kind == "l2":
        return LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
    if kind == "enet":
        return LogisticRegression(C=1.0, l1_ratio=0.5, solver="saga", class_weight="balanced", max_iter=10000,
                                  random_state=config.SEED)   # saga shuffles; fixed seed for reproducibility
    raise ValueError(kind)


def fit_predict(L: Layer, y: np.ndarray, tr: np.ndarray, te: np.ndarray, spec: dict) -> tuple[np.ndarray, dict]:
    """Fit one pipeline on the training people only and return test log-odds. Every data-dependent step
    (expression filter, residualisation, imputation, feature selection, scaling) uses the training fold only."""
    info: dict = {}
    if spec["kind"] == "baseline":
        Xtr, Xte = _impute(L.cov[tr], L.cov[te])
    else:
        Xtr, Xte = L.X[tr], L.X[te]
        if L.expressed is not None:
            keep = L.expressed[tr].mean(axis=0) >= EXPRESSED_FRACTION
            Xtr, Xte = Xtr[:, keep], Xte[:, keep]
            kept_idx = np.flatnonzero(keep)
        else:
            kept_idx = np.arange(Xtr.shape[1])
        if spec["resid"]:
            Ctr, Cte = _impute(L.cov[tr], L.cov[te])
            Dtr = np.column_stack([np.ones(len(tr)), Ctr])
            Dte = np.column_stack([np.ones(len(te)), Cte])
            B = np.linalg.solve(Dtr.T @ Dtr, Dtr.T @ Xtr)      # OLS on the (people x 3) design, training fold only
            Xtr, Xte = Xtr - Dtr @ B, Xte - Dte @ B
        t = welch_t(Xtr, y[tr])
        sel = np.argsort(-np.abs(t), kind="stable")[: spec["k"]]
        Xtr, Xte = Xtr[:, sel], Xte[:, sel]
        info["selected"] = kept_idx[sel]
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0, ddof=1)
    sd = np.where(sd > 0, sd, 1.0)
    Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
    clf = _make_clf(spec["clf"])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        clf.fit(Xtr, y[tr])
    info["n_warnings"] = sum(1 for x in w if "converge" in str(x.message).lower())
    info["n_nonzero"] = int((clf.coef_ != 0).sum())
    return clf.decision_function(Xte), info


def make_splits(y: np.ndarray, seed: int) -> list[list[tuple[np.ndarray, np.ndarray]]]:
    from sklearn.model_selection import StratifiedKFold
    return [list(StratifiedKFold(N_FOLDS, shuffle=True, random_state=seed + r).split(np.zeros(len(y)), y))
            for r in range(N_REPEATS)]


def cv_scores(layers: list[Layer], y: np.ndarray, spec: dict, seed: int, collect: bool = False):
    """Repeated stratified CV. With several layers (same people, same order) the per-layer log-odds are averaged
    (late fusion). Returns the OOF score matrix (repeats x people) and, if collect, per-fold diagnostics."""
    splits = make_splits(y, seed)
    S = np.zeros((N_REPEATS, len(y)))
    diag = []
    for r, folds in enumerate(splits):
        for tr, te in folds:
            acc = np.zeros(len(te))
            for L in layers:
                s, info = fit_predict(L, y, tr, te, spec)
                acc += s
                if collect:
                    diag.append({"repeat": r, "layer": L.name, **info})
            S[r, te] = acc / len(layers)
    return S, diag


def _perm_worker(layers: list[Layer], y: np.ndarray, spec: dict, perm_ids: list[int]) -> list[float]:
    out = []
    for i in perm_ids:
        rng = np.random.default_rng([config.SEED, 7, i])
        yp = rng.permutation(y)
        S, _ = cv_scores(layers, yp, spec, seed=config.SEED + 10_000 + 97 * i)
        out.append(float(auroc_rows(yp, S).mean()))
    return out


def permutation_null(layers: list[Layer], y: np.ndarray, spec: dict, n_perm: int, workers: int) -> np.ndarray:
    from joblib import Parallel, delayed
    chunks = [list(c) for c in np.array_split(np.arange(n_perm), max(1, min(workers * 4, n_perm)))]
    res = Parallel(n_jobs=workers)(delayed(_perm_worker)(layers, y, spec, c) for c in chunks if len(c))
    return np.array([v for r in res for v in r])


def boot_auroc(y: np.ndarray, S: np.ndarray, n_boot: int = N_BOOT, seed: int = config.SEED,
               S_ref: np.ndarray | None = None) -> np.ndarray:
    """Stratified participant bootstrap of the mean-over-repeats OOF AUROC (or of S minus S_ref, paired)."""
    rng = np.random.default_rng([seed, 11])
    cases, ctrls = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(cases, len(cases)), rng.choice(ctrls, len(ctrls))])
        yb = y[idx]
        v = auroc_rows(yb, S[:, idx]).mean()
        if S_ref is not None:
            v -= auroc_rows(yb, S_ref[:, idx]).mean()
        out[b] = v
    return out


def _perm_p(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= obs - 1e-12)) / (1 + len(null)))


# ============================================================================================ analyses


def cohort_table(layers: dict[str, Layer]) -> pd.DataFrame:
    rows = []
    for lay, L in layers.items():
        for g, lab in ((1, "PI-ME/CFS"), (0, "HV")):
            m = L.y == g
            fem = L.cov[m, L.cov_names.index("female")] if "female" in L.cov_names else np.full(m.sum(), np.nan)
            age = L.cov[m, L.cov_names.index("age")]
            rows.append({"layer": lay, "group": lab, "n": int(m.sum()),
                         "n_female": int(np.nansum(fem)) if "female" in L.cov_names else UNKNOWN,
                         "n_sex_known": int(np.isfinite(fem).sum()) if "female" in L.cov_names else 0,
                         "age_mean": round(float(np.nanmean(age)), 1) if np.isfinite(age).any() else np.nan,
                         "age_sd": round(float(np.nanstd(age, ddof=1)), 1) if np.isfinite(age).sum() > 1 else np.nan,
                         "n_age_missing": int(np.isnan(age).sum()),
                         "n_features": L.X.shape[1]})
    return pd.DataFrame(rows)


def run_models(layers: dict[str, Layer], n_perm: int, workers: int, log=print):
    perf, nulls, sel_rows = [], [], []
    oof = {}
    for lay, L in layers.items():
        for mid, spec in MODELS.items():
            t0 = time.time()
            S, diag = cv_scores([L], L.y, spec, seed=config.SEED, collect=True)
            oof[(lay, mid)] = S
            per_rep = auroc_rows(L.y, S)
            obs = float(per_rep.mean())
            bt = boot_auroc(L.y, S)
            row = {"layer": lay, "model": mid, "model_label": spec["label"], "n": len(L.y),
                   "n_cases": int(L.y.sum()), "n_controls": int((L.y == 0).sum()),
                   "auroc": obs, "ci_low": float(np.percentile(bt, 2.5)), "ci_high": float(np.percentile(bt, 97.5)),
                   "repeat_min": float(per_rep.min()), "repeat_max": float(per_rep.max()),
                   "n_convergence_warnings": int(sum(d["n_warnings"] for d in diag)),
                   "median_nonzero_coefs": float(np.median([d["n_nonzero"] for d in diag])),
                   "folds_all_zero_coefs": int(sum(d["n_nonzero"] == 0 for d in diag)),
                   "covariates": ",".join(L.cov_names) if spec["kind"] == "baseline" or spec["resid"] else ""}
            if mid in PERMUTED_MODELS and n_perm > 0:
                null = permutation_null([L], L.y, spec, n_perm, workers)
                row.update({"perm_p": _perm_p(obs, null), "n_perm": len(null), "null_mean": float(null.mean()),
                            "null_q95": float(np.quantile(null, 0.95))})
                nulls.append(pd.DataFrame({"layer": lay, "model": mid, "perm_index": np.arange(len(null)),
                                           "null_auroc": null}))
            else:
                row.update({"perm_p": np.nan, "n_perm": 0, "null_mean": np.nan, "null_q95": np.nan})
            if mid == "l2_top50":
                cnt = pd.Series(np.concatenate([d["selected"] for d in diag])).value_counts()
                sel_rows.append(pd.DataFrame({"layer": lay, "feature_id": L.features[cnt.index.to_numpy()],
                                              "feature_name": L.feature_names[cnt.index.to_numpy()],
                                              "n_folds_selected": cnt.to_numpy(),
                                              "n_folds": N_REPEATS * N_FOLDS}))
            perf.append(row)
            log(f"[linked_omics] {lay:17s} {mid:15s} AUROC {obs:.3f} [{row['ci_low']:.3f}, {row['ci_high']:.3f}] "
                f"perm p {row['perm_p']:.4f} ({time.time() - t0:.0f}s)")
        # omics vs baseline, paired
        for mid in ("l2_top50", "enet_top50"):
            d = boot_auroc(L.y, oof[(lay, mid)], S_ref=oof[(lay, "age_sex")])
            obs = auroc_rows(L.y, oof[(lay, mid)]).mean() - auroc_rows(L.y, oof[(lay, "age_sex")]).mean()
            perf.append({"layer": lay, "model": f"delta:{mid}-age_sex", "model_label": f"{mid} minus age/sex (paired)",
                         "n": len(L.y), "n_cases": int(L.y.sum()), "n_controls": int((L.y == 0).sum()),
                         "auroc": float(obs), "ci_low": float(np.percentile(d, 2.5)),
                         "ci_high": float(np.percentile(d, 97.5))})
    perf = pd.DataFrame(perf)
    perf["passes_primary"] = np.where(perf["model"] == "l2_top50", perf["perm_p"] < PASS_ALPHA, pd.NA)
    return perf, (pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame()), pd.concat(sel_rows), oof


def run_fusion(layers: dict[str, Layer], n_perm: int, workers: int, log=print) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, nulls = [], []
    spec = MODELS["l2_top50"]
    for fid, lays in FUSIONS.items():
        common = sorted(set.intersection(*(set(layers[x].ids) for x in lays)))
        subs = [layers[x].subset(common) for x in lays]
        y = subs[0].y
        base = {"fusion": fid, "layers": "+".join(lays), "role": FUSION_ROLE[fid], "n": len(y),
                "n_cases": int(y.sum()), "n_controls": int((y == 0).sum())}
        if fid == "F5" or min(y.sum(), (y == 0).sum()) < MIN_CLASS_FOR_CV:
            rows.append({**base, "analysed": False,
                         "reason": "pre-specified n-only" if fid == "F5" else f"a class has < {MIN_CLASS_FOR_CV}"})
            continue
        t0 = time.time()
        S, _ = cv_scores(subs, y, spec, seed=config.SEED)
        obs = float(auroc_rows(y, S).mean())
        bt = boot_auroc(y, S)
        singles = {}
        for L in subs:
            Ss, _ = cv_scores([L], y, spec, seed=config.SEED)
            singles[L.name] = (float(auroc_rows(y, Ss).mean()), Ss)
        best = max(singles, key=lambda k: singles[k][0])
        dlt = boot_auroc(y, S, S_ref=singles[best][1])
        null = permutation_null(subs, y, spec, n_perm, workers) if n_perm > 0 else np.array([])
        rows.append({**base, "analysed": True, "reason": "", "auroc": obs,
                     "ci_low": float(np.percentile(bt, 2.5)), "ci_high": float(np.percentile(bt, 97.5)),
                     "perm_p": _perm_p(obs, null) if len(null) else np.nan, "n_perm": len(null),
                     "null_q95": float(np.quantile(null, 0.95)) if len(null) else np.nan,
                     **{f"single_{k}": v[0] for k, v in singles.items()},
                     "best_single_layer": best, "best_single_auroc": singles[best][0],
                     "delta_vs_best_single": obs - singles[best][0],
                     "delta_ci_low": float(np.percentile(dlt, 2.5)), "delta_ci_high": float(np.percentile(dlt, 97.5))})
        if len(null):
            nulls.append(pd.DataFrame({"layer": fid, "model": "fusion_l2_top50", "perm_index": np.arange(len(null)),
                                       "null_auroc": null}))
        log(f"[linked_omics] fusion {fid} n={len(y)} AUROC {obs:.3f} perm p {rows[-1]['perm_p']:.4f} "
            f"({time.time() - t0:.0f}s)")
    return pd.DataFrame(rows), (pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame())


def _ols_group_p(X: np.ndarray, y: np.ndarray, cov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    C, _ = _impute(cov, cov[:1])
    D = np.column_stack([np.ones(len(y)), y, C])
    B, *_ = np.linalg.lstsq(D, X, rcond=None)
    resid = X - D @ B
    df = len(y) - D.shape[1]
    s2 = (resid ** 2).sum(axis=0) / df
    XtX_inv = np.linalg.inv(D.T @ D)
    se = np.sqrt(s2 * XtX_inv[1, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = B[1] / se
    return B[1], 2 * stats.t.sf(np.abs(t), df)


def univariate(layers: dict[str, Layer]) -> pd.DataFrame:
    frames = []
    for lay, L in layers.items():
        X = L.X
        keep = np.ones(X.shape[1], bool) if L.expressed is None else L.expressed.mean(axis=0) >= EXPRESSED_FRACTION
        X = X[:, keep]
        a, b = X[L.y == 1], X[L.y == 0]
        n1, n0 = len(a), len(b)
        t, p = stats.ttest_ind(a, b, axis=0, equal_var=False)
        diff = a.mean(axis=0) - b.mean(axis=0)
        sp = np.sqrt(((n1 - 1) * a.var(axis=0, ddof=1) + (n0 - 1) * b.var(axis=0, ddof=1)) / (n1 + n0 - 2))
        J = 1 - 3 / (4 * (n1 + n0) - 9)
        with np.errstate(divide="ignore", invalid="ignore"):
            g = J * diff / sp
        se = np.sqrt((n1 + n0) / (n1 * n0) + g ** 2 / (2 * (n1 + n0)))
        beta, p_adj = _ols_group_p(X, L.y.astype(float), L.cov)
        frames.append(pd.DataFrame({
            "layer": lay, "feature_id": L.features[keep], "feature_name": L.feature_names[keep],
            "n_cases": n1, "n_controls": n0, "mean_cases": a.mean(axis=0), "mean_controls": b.mean(axis=0),
            "diff_log2": diff, "hedges_g": g, "g_se": se, "g_ci_low": g - 1.96 * se, "g_ci_high": g + 1.96 * se,
            "welch_t": t, "p_value": p, "q_value_bh": bh(p),
            "ols_group_beta": beta, "ols_p_value": p_adj, "ols_q_value_bh": bh(p_adj),
            "ols_covariates": "+".join(L.cov_names)}))
    u = pd.concat(frames, ignore_index=True)
    u["fdr_significant"] = u["q_value_bh"] < FDR_ALPHA
    u["ols_fdr_significant"] = u["ols_q_value_bh"] < FDR_ALPHA
    return u


def univariate_summary(u: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lay, g in u.groupby("layer", sort=False):
        p = g["p_value"].to_numpy()
        rows.append({"layer": lay, "n_features_tested": len(g), "n_cases": int(g["n_cases"].iloc[0]),
                     "n_controls": int(g["n_controls"].iloc[0]),
                     "n_nominal_p_lt_0.05": int((p < 0.05).sum()),
                     "expected_nominal_under_null": round(0.05 * len(g), 1),
                     "n_fdr_significant": int(g["fdr_significant"].sum()),
                     "n_ols_fdr_significant": int(g["ols_fdr_significant"].sum()),
                     "min_q_value": float(np.nanmin(g["q_value_bh"])),
                     "pi0_storey_lambda0.5": float(min(1.0, (p > 0.5).mean() / 0.5)),
                     "median_abs_hedges_g": float(np.nanmedian(np.abs(g["hedges_g"])))})
    return pd.DataFrame(rows)


def _norm_name(s) -> str:
    return " ".join(str(s).lower().replace("*", "").split())


def paper_concordance(u: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    from ..store import read_table
    ev = read_table("condition_molecular_evidence__mapmecfs")
    rows, joined_all = [], []
    for lay, tid in PAPER_TABLE.items():
        pub = ev[ev["table_id"] == tid].copy()
        ours = u[u["layer"] == lay].copy()
        if lay == "pbmc_rnaseq":
            pub["key"], ours["key"] = pub["ensembl_gene_id"].astype(str), ours["feature_id"].astype(str)
        elif lay == "muscle_rnaseq":
            pub["key"], ours["key"] = pub["gene_symbol"].astype(str), ours["feature_name"].astype(str)
        elif lay.endswith("somalogic"):
            pub["key"], ours["key"] = pub["somalogic_seqid"].astype(str), ours["feature_id"].astype(str)
        else:
            pub["key"], ours["key"] = pub["analyte_name"].map(_norm_name), ours["feature_id"].map(_norm_name)
        pub = pub[pub["direction"].isin(["up_in_PI-ME/CFS", "down_in_PI-ME/CFS"])].drop_duplicates("key")
        ours = ours.drop_duplicates("key")
        j = pub[["key", "direction", "effect_value", "p_value", "adj_p_value", "n_cases", "n_controls"]].rename(
            columns={"p_value": "pub_p", "adj_p_value": "pub_adj_p", "effect_value": "pub_effect",
                     "n_cases": "pub_n_cases", "n_controls": "pub_n_controls"}).merge(
            ours[["key", "feature_id", "feature_name", "diff_log2", "hedges_g", "p_value", "q_value_bh"]], on="key")
        j["pub_sign"] = np.where(j["direction"] == "up_in_PI-ME/CFS", 1, -1)
        j["our_sign"] = np.sign(j["diff_log2"])
        j["agree"] = j["pub_sign"] == j["our_sign"]
        j.insert(0, "layer", lay)
        j.insert(1, "paper_table", f"Supplementary Data {tid}")
        joined_all.append(j)
        nom = j[j["pub_p"] < 0.05]
        k, n = int(nom["agree"].sum()), len(nom)
        rho = stats.spearmanr(j["pub_effect"], j["diff_log2"]).statistic if lay.endswith("somalogic") else np.nan
        rows.append({"layer": lay, "paper_table": f"Supplementary Data {tid}",
                     "paper_scope": pub["table_scope"].iloc[0] if len(pub) else UNKNOWN,
                     "paper_n_cases": int(pub["n_cases"].iloc[0]) if len(pub) else np.nan,
                     "paper_n_controls": int(pub["n_controls"].iloc[0]) if len(pub) else np.nan,
                     "n_published_with_direction": len(pub), "n_matched": len(j),
                     "n_paper_nominal_matched": n, "n_sign_agree": k,
                     "sign_agreement": k / n if n else np.nan,
                     "binomial_p": float(stats.binomtest(k, n, 0.5).pvalue) if n else np.nan,
                     "n_paper_nominal_also_our_p_lt_0.05": int((nom["p_value"] < 0.05).sum()),
                     "spearman_effect_all_matched": float(rho) if np.isfinite(rho) else np.nan,
                     "interpretation": "reproduction check on overlapping people (not independent replication)"})
    return pd.DataFrame(rows), pd.concat(joined_all, ignore_index=True)


def csf_serum_consistency(layers: dict[str, Layer], n_perm: int, top: int = 50) -> tuple[pd.DataFrame, pd.DataFrame]:
    A, B = layers["csf_somalogic"], layers["serum_somalogic"]
    common = sorted(set(A.ids) & set(B.ids))
    A, B = A.subset(common), B.subset(common)
    # exact SomaLogic SeqId match (a few aptamers carry a different SeqId version suffix in the two runs; excluded)
    shared = np.array(sorted(set(A.features) & set(B.features)))
    n_unmatched = len(set(A.features) ^ set(B.features))

    def _cols(L):
        pos = pd.Index(L.features).get_indexer(shared)
        return Layer(L.name, L.ids, L.y, L.X[:, pos], L.features[pos], L.feature_names[pos], L.cov, L.cov_names)

    A, B = _cols(A), _cols(B)
    y = A.y

    def stat(yy):
        ta, tb = welch_t(A.X, yy), welch_t(B.X, yy)
        rho = stats.spearmanr(ta, tb).statistic
        ov = len(set(np.argsort(-np.abs(ta))[:top]) & set(np.argsort(-np.abs(tb))[:top]))
        return float(rho), ov, ta, tb

    rho, ov, ta, tb = stat(y)
    rng = np.random.default_rng([config.SEED, 23])
    null = np.array([stat(rng.permutation(y))[:2] for _ in range(n_perm)]) if n_perm else np.empty((0, 2))
    pa = stats.ttest_ind(A.X[y == 1], A.X[y == 0], axis=0, equal_var=False).pvalue
    pb = stats.ttest_ind(B.X[y == 1], B.X[y == 0], axis=0, equal_var=False).pvalue
    qa, qb = bh(pa), bh(pb)
    sa, sb = qa < FDR_ALPHA, qb < FDR_ALPHA
    either = sa | sb
    coupling = np.array([stats.spearmanr(A.X[:, j], B.X[:, j]).statistic for j in range(A.X.shape[1])])
    summary = pd.DataFrame([
        {"quantity": "people (same in both fluids)", "value": len(y)},
        {"quantity": "aptamers (exact SeqId match in both runs)", "value": A.X.shape[1]},
        {"quantity": "aptamer ids present in only one run (SeqId version differs; excluded)", "value": n_unmatched},
        {"quantity": "spearman_rho_of_welch_t_csf_vs_serum", "value": rho},
        {"quantity": "joint_permutation_null_mean_rho", "value": float(null[:, 0].mean()) if len(null) else np.nan},
        {"quantity": "joint_permutation_null_q95_rho", "value": float(np.quantile(null[:, 0], 0.95)) if len(null) else np.nan},
        {"quantity": "joint_permutation_p_rho (one-sided)", "value": _perm_p(rho, null[:, 0]) if len(null) else np.nan},
        {"quantity": f"top{top}_overlap_count", "value": ov},
        {"quantity": f"top{top}_overlap_expected_random", "value": top * top / A.X.shape[1]},
        {"quantity": f"top{top}_overlap_joint_permutation_mean", "value": float(null[:, 1].mean()) if len(null) else np.nan},
        {"quantity": f"top{top}_overlap_joint_permutation_p (one-sided)", "value": _perm_p(ov, null[:, 1]) if len(null) else np.nan},
        {"quantity": "n_fdr_significant_csf", "value": int(sa.sum())},
        {"quantity": "n_fdr_significant_serum", "value": int(sb.sum())},
        {"quantity": "n_fdr_significant_both", "value": int((sa & sb).sum())},
        {"quantity": "sign_agreement_among_fdr_significant_in_either",
         "value": float((np.sign(ta[either]) == np.sign(tb[either])).mean()) if either.any() else np.nan},
        {"quantity": "n_permutations", "value": len(null)},
        {"quantity": "coupling_median_spearman_csf_vs_serum_level_across_people", "value": float(np.nanmedian(coupling))},
        {"quantity": "coupling_fraction_aptamers_rho_gt_0.5", "value": float(np.nanmean(coupling > 0.5))},
        {"quantity": "coupling_fraction_aptamers_rho_gt_0.3", "value": float(np.nanmean(coupling > 0.3))},
    ])
    per = pd.DataFrame({"feature_id": A.features, "feature_name": A.feature_names, "t_csf": ta, "t_serum": tb,
                        "q_csf": qa, "q_serum": qb, "csf_serum_level_spearman": coupling})
    return summary, per


# ============================================================================================ signature partition


def _provenance_source() -> dict:
    from ..store import read_table
    o = read_table("participant_omics_linked__mapmecfs",
                   columns=["modality", "source_name", "source_version", "retrieved_at"]).drop_duplicates("modality")
    return {r.modality: r for r in o.itertuples(index=False)}


def build_signatures(u: pd.DataFrame, perf: pd.DataFrame) -> pd.DataFrame:
    from ..provenance import add_provenance

    src = _provenance_source()
    frames = []
    for lay, g in u.groupby("layer", sort=False):
        pid = f"pi_mecfs_vs_hv__{lay}"
        passes = bool(perf.loc[(perf.layer == lay) & (perf.model == "l2_top50"), "passes_primary"].fillna(False).iloc[0])
        ev_level = ("single_cohort_person_level_omics_classifier_passes_permutation" if passes
                    else "single_cohort_person_level_omics_classifier_null")
        n1, n0 = int(g["n_cases"].iloc[0]), int(g["n_controls"].iloc[0])
        common = {
            "dataset_id": DATASET_ID, "phenotype_id": pid,
            "phenotype_label": f"PI-ME/CFS vs healthy volunteer: {LAYER_LABEL[lay]} (NIH intramural cohort)",
            "condition_id": "me_cfs", "is_proxy": False, "label_basis": LABEL_BASIS,
            "phenotype_definition": ("group_label as deposited (1xx study number = healthy volunteer, 3xx = PI-ME/CFS; "
                                     "concordant across deposits); cohort per Walitt et al. 2024 Nat Commun 15:907"),
            "population": f"NIH Clinical Center PI-ME/CFS study participants with open {lay} data "
                          f"({n1} PI-ME/CFS, {n0} HV)",
            "cycles": "not applicable", "feature_units": LAYER_UNITS[lay],
            "n_cases": n1, "n_controls": n0, "underpowered": True, "caveats": COHORT_CAVEAT,
            "evidence_level": ev_level,
        }
        f = pd.DataFrame({
            **common,
            "feature": g["feature_id"].astype(str).str.replace("|", "/", regex=False).to_numpy(),
            "feature_label": g["feature_name"].astype(str).to_numpy(),
            "effect_measure": "hedges_g", "effect_unit": "Hedges' g (PI-ME/CFS minus HV), unadjusted",
            "effect_size": g["hedges_g"].to_numpy(), "ci_low": g["g_ci_low"].to_numpy(),
            "ci_high": g["g_ci_high"].to_numpy(), "null_value": 0.0, "se": g["g_se"].to_numpy(),
            "p_value": g["p_value"].to_numpy(), "q_value_bh": g["q_value_bh"].to_numpy(),
            "q_value": g["q_value_bh"].to_numpy(), "fdr_significant": g["fdr_significant"].to_numpy(bool),
            "direction": np.where(g["hedges_g"] > 0, "higher_in_cases", "lower_in_cases"),
            "feature_transform": (f"{LAYER_UNITS[lay]}" + ("; genes with CPM >= 1 in >= 50% of samples"
                                                              if lay in RNA_LAYERS else "")),
            "covariates": "none (Welch t); OLS adjusted for " + "+".join(LAYER_COVARIATES[lay])
                          + " in results/tables/linked_omics_univariate.csv",
            "mean_cases": g["mean_cases"].to_numpy(), "mean_controls": g["mean_controls"].to_numpy(),
            "method": "Welch t-test on log-scale values; BH-FDR within layer; Hedges' g with normal-approximation CI",
        })
        frames.append(f)
        pr = perf[(perf.layer == lay) & perf.model.isin(list(MODELS))]
        m = pd.DataFrame({
            **common,
            "feature": "model_auroc:" + lay + ":" + pr["model"].astype(str),
            "feature_label": pr["model_label"].to_numpy(),
            "effect_measure": "cross-validated AUROC (10 x 5 repeated stratified CV, mean over repeats)",
            "effect_unit": "AUROC", "effect_size": pr["auroc"].to_numpy(), "ci_low": pr["ci_low"].to_numpy(),
            "ci_high": pr["ci_high"].to_numpy(), "null_value": 0.5, "se": np.nan,
            "p_value": pr["perm_p"].to_numpy(), "q_value_bh": np.nan, "q_value": np.nan, "fdr_significant": False,
            "direction": np.where(pr["auroc"] > 0.5, "higher_in_cases", "lower_in_cases"),
            "feature_transform": "in-fold preprocessing and top-k selection (docs/ANALYSIS_PLAN_LINKED_OMICS.md)",
            "covariates": pr["covariates"].fillna("").to_numpy(),
            "mean_cases": np.nan, "mean_controls": np.nan,
            "method": "p_value = one-sided label-permutation p (1,000 permutations of the full CV) where run; "
                      "CI = stratified participant bootstrap (2,000)",
        })
        frames.append(m)
    sig = pd.concat(frames, ignore_index=True)
    sig["object_id"] = "signature:" + sig["dataset_id"] + "|" + sig["phenotype_id"] + "|" + sig["feature"]
    if sig["object_id"].duplicated().any():
        raise ValueError("duplicate signature object ids")
    parts = []
    for lay, g in sig.groupby(sig["phenotype_id"].str.split("__").str[1], sort=False):
        r = src[lay]
        parts.append(add_provenance(
            g, data_layer="person", source_name=r.source_name, source_version=f"{r.source_version}; analysis "
            f"{PRODUCER} {ANALYSIS_VERSION}", retrieved_at=r.retrieved_at, evidence_type="person_derived_feature",
            source_record_id=lambda d: d["object_id"].str.replace("signature:", "", regex=False),
            evidence_level=g["evidence_level"],
            provenance_notes=(f"Cohort-level contrast computed by {PRODUCER} from person-level linked omics "
                              "(participant_omics_linked__mapmecfs, NIH study number); aggregate row, not a person "
                              f"row; no geography; not linked to any wearable data; plan {PLAN}")))
    out = pd.concat(parts, ignore_index=True)
    cols = ["object_id", "dataset_id", "phenotype_id", "phenotype_label", "condition_id", "is_proxy", "label_basis",
            "feature", "feature_label", "effect_measure", "effect_size", "ci_low", "ci_high", "null_value", "se",
            "p_value", "q_value_bh", "fdr_significant", "direction", "n_cases", "n_controls", "underpowered",
            "phenotype_definition", "population", "cycles", "feature_units", "feature_transform", "covariates",
            "mean_cases", "mean_controls", "caveats", "data_layer", "source_name", "source_record_id",
            "source_version", "retrieved_at", "source_geographic_resolution", "evidence_type", "evidence_level",
            "provenance_notes", "effect_unit", "q_value", "method"]
    out = out[cols]
    out["condition_id"] = out["condition_id"].astype(object)
    out["n_cases"] = out["n_cases"].astype("int64")
    out["n_controls"] = out["n_controls"].astype("int64")
    return out


# ============================================================================================ figures


def figures(perf: pd.DataFrame, fusion: pd.DataFrame, per_apt: pd.DataFrame, cons: pd.DataFrame) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = config.FIGURES / "drafts"
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    show = ["l2_top50", "enet_top50", "l2_top50_resid", "age_sex"]
    colors = {"l2_top50": "#1f5fa8", "enet_top50": "#6f9bd1", "l2_top50_resid": "#8a6fb3", "age_sex": "#8c8c8c"}
    rows = []
    for lay in LAYERS:
        for mid in show:
            r = perf[(perf.layer == lay) & (perf.model == mid)]
            if len(r):
                rows.append((f"{lay} (n={int(r.n.iloc[0])})", mid, r.iloc[0]))
    for _, f in fusion[fusion.get("analysed", False) == True].iterrows():  # noqa: E712
        rows.append((f"{f['fusion']} {f['layers']} (n={int(f['n'])})", "fusion", f))
    fig, ax = plt.subplots(figsize=(10, 0.32 * len(rows) + 1.9))
    ylabels = []
    for i, (lab, mid, r) in enumerate(rows):
        yv = len(rows) - 1 - i
        c = colors.get(mid, "#c0392b")
        ax.plot([r["ci_low"], r["ci_high"]], [yv, yv], color=c, lw=2)
        ax.plot(r["auroc"], yv, "o", color=c, ms=5)
        if np.isfinite(r.get("null_q95", np.nan)):
            ax.plot(r["null_q95"], yv, "|", color="black", ms=9)
        p = r.get("perm_p", np.nan)
        ax.text(1.005, yv, f"p={p:.3f}" if np.isfinite(p) else "", va="center", fontsize=7,
                transform=ax.get_yaxis_transform())
        ylabels.append(f"{lab}  {mid}")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(ylabels[::-1], fontsize=7)
    ax.axvline(0.5, color="#999999", ls="--", lw=1)
    ax.set_xlim(0.1, 1.0)
    ax.set_xlabel("cross-validated AUROC (mean of 10 x 5 CV) with 95% participant-bootstrap CI\n"
                  "black tick = 95th percentile of the 1,000-permutation null; p = one-sided permutation p", fontsize=8)
    ax.set_title("mapMECFS linked omics, PI-ME/CFS vs healthy volunteer\n(person level, omics only; no wearable data "
                 "linked; pre-specified pass threshold p < 0.01)", fontsize=9)
    fig.tight_layout()
    p1 = out_dir / "linked_omics_auroc_forest.png"
    fig.savefig(p1, dpi=160)
    plt.close(fig)
    written.append(str(p1.relative_to(config.PROJECT_ROOT)))

    cv = dict(zip(cons["quantity"], cons["value"]))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    ax = axes[0]
    ax.scatter(per_apt["t_csf"], per_apt["t_serum"], s=6, alpha=0.5, color="#1f5fa8", edgecolor="none")
    ax.axhline(0, color="#999", lw=0.8)
    ax.axvline(0, color="#999", lw=0.8)
    ax.set_xlabel("Welch t, CSF SomaScan (PI-ME/CFS minus HV)")
    ax.set_ylabel("Welch t, serum SomaScan (same people)")
    ax.set_title(f"per-aptamer t: Spearman rho {cv['spearman_rho_of_welch_t_csf_vs_serum']:.3f}, "
                 f"joint-perm p {cv['joint_permutation_p_rho (one-sided)']:.3f}", fontsize=9)
    ax = axes[1]
    ax.hist(per_apt["csf_serum_level_spearman"].dropna(), bins=40, color="#8a6fb3")
    ax.set_xlabel("per-aptamer Spearman rho, CSF vs serum level across the 42 people")
    ax.set_ylabel("aptamers")
    ax.set_title("does a person's CSF level track their serum level?", fontsize=9)
    fig.tight_layout()
    p2 = out_dir / "linked_omics_csf_vs_serum.png"
    fig.savefig(p2, dpi=160)
    plt.close(fig)
    written.append(str(p2.relative_to(config.PROJECT_ROOT)))
    return written


# ============================================================================================ entry points


def _tables_dir():
    config.TABLES.mkdir(parents=True, exist_ok=True)
    return config.TABLES


def run(n_perm: int = N_PERM, workers: int = 16, log=print) -> dict:
    t0 = time.time()
    T = _tables_dir()
    layers = load_layers()
    cohort = cohort_table(layers)
    cohort.to_csv(T / "linked_omics_cohort.csv", index=False)
    perf, nulls, sel, _ = run_models(layers, n_perm, workers, log=log)
    fusion, fnull = run_fusion(layers, n_perm, workers, log=log)
    u = univariate(layers)
    usum = univariate_summary(u)
    conc, conc_rows = paper_concordance(u)
    cons, per_apt = csf_serum_consistency(layers, n_perm)
    perf["analysis_version"] = ANALYSIS_VERSION
    primary = perf[perf.model.isin(list(MODELS))]
    primary.to_csv(T / "linked_omics_model_performance.csv", index=False)
    perf[perf.model.str.startswith("delta:")].to_csv(T / "linked_omics_vs_baseline.csv", index=False)
    perf[perf.model.isin(["l2_top50", "l2_top50_resid", "l2_top10", "l2_top200"])].to_csv(
        T / "linked_omics_sensitivity.csv", index=False)
    pd.concat([nulls, fnull], ignore_index=True).to_csv(T / "linked_omics_permutation_null.csv", index=False)
    fusion.to_csv(T / "linked_omics_fusion.csv", index=False)
    u.to_csv(T / "linked_omics_univariate.csv", index=False)
    usum.to_csv(T / "linked_omics_univariate_summary.csv", index=False)
    conc.to_csv(T / "linked_omics_paper_concordance.csv", index=False)
    conc_rows.to_csv(T / "linked_omics_paper_concordance_features.csv", index=False)
    cons.to_csv(T / "linked_omics_csf_serum_consistency.csv", index=False)
    per_apt.to_csv(T / "linked_omics_csf_serum_per_aptamer.csv", index=False)
    sel.to_csv(T / "linked_omics_selection_frequency.csv", index=False)
    meta = {"analysis_version": ANALYSIS_VERSION, "seed": config.SEED, "n_perm": n_perm, "n_boot": N_BOOT,
            "cv": f"{N_REPEATS}x{N_FOLDS} repeated stratified", "k_primary": K_PRIMARY, "pass_alpha": PASS_ALPHA,
            "workers": workers, "runtime_s": round(time.time() - t0, 1),
            "finished_at": pd.Timestamp.now(tz="UTC").isoformat(), "plan": PLAN}
    (T / "linked_omics_run_metadata.json").write_text(json.dumps(meta, indent=2))
    finish(log=log)
    return meta


def finish(log=print) -> None:
    """From the saved tables: top features of passing layers, signature partition, figures, report."""
    from ..store import write_table

    T = config.TABLES
    perf = pd.read_csv(T / "linked_omics_model_performance.csv")
    u = pd.read_csv(T / "linked_omics_univariate.csv")
    fusion = pd.read_csv(T / "linked_omics_fusion.csv")
    cons = pd.read_csv(T / "linked_omics_csf_serum_consistency.csv")
    per_apt = pd.read_csv(T / "linked_omics_csf_serum_per_aptamer.csv")
    sel = pd.read_csv(T / "linked_omics_selection_frequency.csv")
    passing = perf[(perf.model == "l2_top50") & (perf.perm_p < PASS_ALPHA)]["layer"].tolist()
    tops = []
    for lay in passing:
        g = u[u.layer == lay].merge(sel[sel.layer == lay][["feature_id", "n_folds_selected"]], on="feature_id",
                                    how="left").fillna({"n_folds_selected": 0})
        tops.append(g.sort_values("p_value").head(25))
    top_cols = ["layer", "feature_id", "feature_name", "hedges_g", "g_ci_low", "g_ci_high", "diff_log2", "p_value",
                "q_value_bh", "ols_p_value", "n_folds_selected"]
    (pd.concat(tops)[top_cols] if tops else pd.DataFrame(columns=top_cols + ["note"])).assign(
        note="reported only because this layer passed its pre-specified permutation test" if tops else
        "no layer passed its pre-specified permutation test; top features deliberately not reported").to_csv(
        T / "linked_omics_top_features.csv", index=False)
    sig = build_signatures(u, perf.assign(passes_primary=perf["perm_p"].lt(PASS_ALPHA).where(perf.model == "l2_top50")))
    write_table(sig, SIGNATURE_TABLE, producer=PRODUCER,
                description="PI-ME/CFS vs HV person-level omics contrasts (mapMECFS linked omics, NIH study number): "
                            "per-feature Hedges' g + BH q per layer and cross-validated AUROC rows; not linked to "
                            "wearables; plan docs/ANALYSIS_PLAN_LINKED_OMICS.md")
    figs = figures(perf, fusion, per_apt, cons)
    write_report()
    log(f"[linked_omics] signatures {len(sig)} rows; figures {figs}; passing layers {passing or 'none'}")


def _f(v, d=3):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "n/a"
    return "n/a" if not np.isfinite(v) else f"{v:.{d}f}"


def _pfmt(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "n/a"
    if not np.isfinite(v):
        return "n/a"
    return f"{v:.4f}" if v >= 0.001 else f"{v:.1e}"


def write_report() -> None:
    """results/LINKED_OMICS_RESULTS.md, generated from the saved tables (numbers are never typed by hand)."""
    T = config.TABLES
    perf = pd.read_csv(T / "linked_omics_model_performance.csv")
    delta = pd.read_csv(T / "linked_omics_vs_baseline.csv")
    fusion = pd.read_csv(T / "linked_omics_fusion.csv")
    usum = pd.read_csv(T / "linked_omics_univariate_summary.csv")
    conc = pd.read_csv(T / "linked_omics_paper_concordance.csv")
    cons = pd.read_csv(T / "linked_omics_csf_serum_consistency.csv")
    cohort = pd.read_csv(T / "linked_omics_cohort.csv")
    top = pd.read_csv(T / "linked_omics_top_features.csv")
    nulls = pd.read_csv(T / "linked_omics_permutation_null.csv")
    meta = json.loads((T / "linked_omics_run_metadata.json").read_text())
    cv = dict(zip(cons["quantity"], cons["value"]))
    prim = perf[perf.model == "l2_top50"].set_index("layer")
    passing = [lay for lay in LAYERS if prim.loc[lay, "perm_p"] < PASS_ALPHA]
    nominal = [lay for lay in LAYERS if PASS_ALPHA <= prim.loc[lay, "perm_p"] < 0.05]

    def row(lay, mid):
        r = perf[(perf.layer == lay) & (perf.model == mid)]
        return r.iloc[0] if len(r) else None

    L = []
    w = L.append
    w("# Person-level linked omics (mapMECFS / NIH PI-ME/CFS): results")
    w("")
    w(f"Generated by `{PRODUCER}` (analysis {ANALYSIS_VERSION}) from `results/tables/linked_omics_*.csv`; run finished "
      f"{meta.get('finished_at', UNKNOWN)} ({meta.get('n_perm')} label permutations per permuted model, "
      f"{meta.get('n_boot', N_BOOT)} bootstrap resamples, seed {meta.get('seed')}). Pre-registered plan: `{PLAN}` "
      "(written before any analysis ran; deviations listed in section 9).")
    w("")
    w("## Bottom line")
    w("")
    w("* **None of this is wearable or device evidence.** No wearable, actigraphy, HRV, orthostatic or CPET measurement "
      "is linked to these people in open data (0 identifier pairs; `docs/MAP_MECFS_LINKAGE_AUDIT.md`). The omics layers "
      "are linked to each other only, by NIH study number.")
    if passing:
        w(f"* **Layers passing the pre-specified test** (primary L2 model, permutation p < {PASS_ALPHA:.2f}, Bonferroni "
          f"over 5 layers): {', '.join(passing)}.")
    else:
        w(f"* **No layer passed the pre-specified test** (primary L2 model, one-sided permutation p < {PASS_ALPHA:.2f}, "
          "Bonferroni over 5 layers).")
    if nominal:
        w(f"* Nominal only (0.01 <= p < 0.05, not passing): {', '.join(nominal)}.")
    w("* Per layer (primary L2 model, 10 x 5 CV, 1,000 permutations):")
    for lay in LAYERS:
        r = prim.loc[lay]
        w(f"  * {lay}: {int(r.n_cases)} PI-ME/CFS / {int(r.n_controls)} HV; AUROC {_f(r.auroc)} "
          f"(95% CI {_f(r.ci_low)}-{_f(r.ci_high)}); permutation-null mean {_f(r.null_mean)}; p = {_pfmt(r.perm_p)}; "
          f"FDR-significant features {int(usum.set_index('layer').loc[lay, 'n_fdr_significant'])}.")
    f1 = fusion[fusion.fusion == "F1"].iloc[0]
    w(f"* Fusing CSF + serum SomaScan (same 42 people): AUROC {_f(f1.auroc)} (CI {_f(f1.ci_low)}-{_f(f1.ci_high)}), "
      f"permutation p = {_pfmt(f1.perm_p)}; vs the better single fluid {_f(f1.delta_vs_best_single, 3)} "
      f"(CI {_f(f1.delta_ci_low)} to {_f(f1.delta_ci_high)}).")
    w(f"* CSF and serum SomaScan do **not** agree on which proteins differ: Spearman rho of per-aptamer t statistics "
      f"{_f(cv['spearman_rho_of_welch_t_csf_vs_serum'])} against a joint-permutation null mean of "
      f"{_f(cv['joint_permutation_null_mean_rho'])} (one-sided p {_pfmt(cv['joint_permutation_p_rho (one-sided)'])}); "
      f"{int(cv['top50_overlap_count'])} of the top-50 aptamers shared (null mean "
      f"{_f(cv['top50_overlap_joint_permutation_mean'], 1)}).")
    w("* Agreement with the paper's published directions is high, but it is a **reproduction check on largely the same "
      "people**, not independent replication (section 6).")
    w("* Small (n = 21-42), single-site, adjudication status not in open files. A layer that passes is a candidate "
      "for replication in an independent cohort, not a biomarker and not a diagnostic.")
    w("")
    w("## Interpretation")
    w("")
    r = {lay: prim.loc[lay] for lay in LAYERS}
    s1 = perf[perf.model == "l2_top50_resid"].set_index("layer")
    us = usum.set_index("layer")
    nom_layers = [lay for lay in LAYERS if prim.loc[lay, "perm_p"] < 0.05]
    w(f"* **The pre-specified answer is null**: no layer reached p < {PASS_ALPHA:.2f}. "
      f"{len(nom_layers)} of 5 layers have nominal p < 0.05 ({', '.join(f'{x} p = {_pfmt(prim.loc[x, 'perm_p'])}' for x in nom_layers)}). "
      "Five nominal tests would give 0.25 such layers on average if everything were null, so this is more than "
      "chance would usually produce, but that count was **not pre-specified**, the layers share people "
      "(CSF and serum SomaScan are the same 42), and no p-value is attached to it. It is a reason to test these "
      "layers in an independent cohort, nothing more.")
    w("* **Monte Carlo error around the thresholds.** A permutation p from 1,000 permutations has a Monte Carlo "
      "standard error of about sqrt(p(1 - p) / 1000): about 0.004 at p = 0.013 and 0.007 at p = 0.049. So CSF "
      f"metabolomics (p = {_pfmt(prim.loc['csf_metabolomics', 'perm_p'])}) is within Monte Carlo error of the "
      "0.01 pass line, and PBMC RNA-seq "
      f"(p = {_pfmt(prim.loc['pbmc_rnaseq', 'perm_p'])}) is within it of the 0.05 nominal line. An independent "
      "re-implementation at review, with a different permutation seed and 400 permutations, gave p = 0.010 for "
      "CSF metabolomics. The pre-specified decision is the one computed here, and it is 'not passing'. Neither "
      "the 'fails' label nor the 'nominal' label is sharp, though. CSF metabolomics would not become evidence even "
      "if it crossed the line: its signal weakens once age is removed (below).")
    w(f"* **CSF metabolomics** has the highest AUROC ({_f(r['csf_metabolomics'].auroc)}) and the only univariate excess "
      f"that survives FDR ({int(us.loc['csf_metabolomics', 'n_fdr_significant'])} metabolite; "
      f"{int(us.loc['csf_metabolomics', 'n_nominal_p_lt_0.05'])} nominal vs "
      f"{us.loc['csf_metabolomics', 'expected_nominal_under_null']} expected; Storey pi0 "
      f"{_f(us.loc['csf_metabolomics', 'pi0_storey_lambda0.5'], 2)}). But the PI-ME/CFS participants in this layer are "
      "about 11 years younger, and after removing age inside the folds (S1) the AUROC falls to "
      f"{_f(s1.loc['csf_metabolomics', 'auroc'])} (p = {_pfmt(s1.loc['csf_metabolomics', 'perm_p'])}); with the OLS "
      f"age adjustment {int(us.loc['csf_metabolomics', 'n_ols_fdr_significant'])} metabolites pass FDR. The layer is "
      "also single-sex, n = 21, from figure source data with an unpublished batch structure. Part of its signal is "
      "age.")
    w(f"* **Serum SomaScan** (AUROC {_f(r['serum_somalogic'].auroc)}, p = {_pfmt(r['serum_somalogic'].perm_p)}) weakens "
      f"after sex/age residualisation (S1 {_f(s1.loc['serum_somalogic', 'auroc'])}, p = "
      f"{_pfmt(s1.loc['serum_somalogic', 'perm_p'])}), has no FDR-significant aptamer and no nominal excess "
      f"({int(us.loc['serum_somalogic', 'n_nominal_p_lt_0.05'])} vs {us.loc['serum_somalogic', 'expected_nominal_under_null']} "
      "expected), and depends on k (S2: k = 10 and k = 200 are lower).")
    w(f"* **PBMC RNA-seq** (AUROC {_f(r['pbmc_rnaseq'].auroc)}, p = {_pfmt(r['pbmc_rnaseq'].perm_p)}) strengthens after "
      f"residualisation (S1 {_f(s1.loc['pbmc_rnaseq', 'auroc'])}, p = {_pfmt(s1.loc['pbmc_rnaseq', 'perm_p'])}), so it "
      "is not a sex-gene artefact, yet the univariate test shows no excess of small p-values "
      f"({int(us.loc['pbmc_rnaseq', 'n_nominal_p_lt_0.05'])} nominal vs {us.loc['pbmc_rnaseq', 'expected_nominal_under_null']} "
      "expected, pi0 " + _f(us.loc['pbmc_rnaseq', 'pi0_storey_lambda0.5'], 2) + "): whatever separates the groups is "
      "diffuse and borderline at n = 27.")
    w(f"* **Muscle RNA-seq** is the mirror image: a univariate excess ({int(us.loc['muscle_rnaseq', 'n_nominal_p_lt_0.05'])} "
      f"nominal vs {us.loc['muscle_rnaseq', 'expected_nominal_under_null']} expected, pi0 "
      f"{_f(us.loc['muscle_rnaseq', 'pi0_storey_lambda0.5'], 2)}, 0 FDR hits) without a classifier that beats its null "
      f"(p = {_pfmt(r['muscle_rnaseq'].perm_p)}). **CSF SomaScan** is null on every test.")
    w("* **Fusion never helped**: every fused model scored below the best single layer on the same people (section 4). "
      "Adding layers adds noise faster than signal at these sample sizes.")
    w("* **Person-level omics beyond condition-level enrichment**: in this cohort, at most weak and unconfirmed "
      "person-level signal (serum proteins, PBMC transcripts, CSF metabolites), none of which passes the "
      "pre-specified bar. It does not upgrade any molecular layer to 'evidence-supported' person-level "
      "discrimination.")
    w("")
    w("## 1. What was tested")
    w("")
    w("Per layer: log transform (RNA-seq log2(CPM + 1); SomaScan log2 RFU; metabolomics log2), then **inside each "
      "training fold** the RNA expression filter, top-50 Welch-t feature selection and standardisation, then a "
      "penalised logistic classifier (primary L2, C = 1; secondary elastic net, l1_ratio 0.5). 10 x 5 repeated "
      "stratified CV; AUROC = mean over repeats of the pooled out-of-fold AUROC; 95% CI from a stratified participant "
      "bootstrap (2,000); one-sided p from 1,000 label permutations each re-running the whole procedure. Baseline: "
      "age + sex only. Sensitivity: S1 features residualised on sex/age inside the fold; S2 k = 10 / 200.")
    w("")
    nm = perf[perf.model.isin(["l2_top50", "enet_top50", "l2_top50_resid"])]["null_mean"]
    w(f"**Reading the AUROCs.** Each AUROC is compared with its own permutation null (the whole procedure re-run on "
      f"shuffled labels). For the omics models the null means are {_f(nm.min(), 2)}-{_f(nm.max(), 2)} and their 95th "
      "percentiles are 0.66-0.75: with 21-42 people, an AUROC near 0.70 is what chance alone produces about one time "
      "in twenty. The two-feature age/sex baselines fall below 0.5 (as low as "
      f"{_f(perf[perf.model == 'age_sex']['auroc'].min(), 2)}): that is small-sample cross-validation anti-learning "
      "with a nearly uninformative model (each baseline is inside its permutation null), not an inverse association.")
    w("")
    w("## 2. Cohort")
    w("")
    w("| layer | group | n | female | age mean (SD) | age missing | features |")
    w("|---|---|---|---|---|---|---|")
    for r in cohort.itertuples(index=False):
        w(f"| {r.layer} | {r.group} | {r.n} | {r.n_female} | {_f(r.age_mean, 1)} ({_f(r.age_sd, 1)}) | "
          f"{r.n_age_missing} | {r.n_features} |")
    w("")
    w("CSF metabolomics: only the 21 rows that carry a study number are linkable, and all of them have Birth Sex = 1 "
      "(recorded Female in GEO for the 20 also in GEO), so the layer is single-sex and its baseline is age only. "
      "The PI-ME/CFS participants in this layer are on average ~11 years younger than the HV.")
    w("")
    w("## 3. Classifiers per layer")
    w("")
    w("| layer | model | cases / controls | AUROC (95% CI) | repeat range | null mean | null 95th pct | "
      "permutation p |")
    w("|---|---|---|---|---|---|---|---|")
    for lay in LAYERS:
        for mid in ("l2_top50", "enet_top50", "age_sex", "l2_top50_resid"):
            r = row(lay, mid)
            if r is None:
                continue
            w(f"| {lay} | {mid} | {int(r.n_cases)} / {int(r.n_controls)} | {_f(r.auroc)} ({_f(r.ci_low)}-"
              f"{_f(r.ci_high)}) | {_f(r.repeat_min)}-{_f(r.repeat_max)} | {_f(r.null_mean)} | {_f(r.null_q95)} | "
              f"{_pfmt(r.perm_p)} |")
    w("")
    w("Omics minus age/sex baseline, same people (paired participant bootstrap; no permutation):")
    w("")
    w("| layer | contrast | ΔAUROC (95% CI) |")
    w("|---|---|---|")
    for r in delta.itertuples(index=False):
        w(f"| {r.layer} | {r.model.replace('delta:', '')} | {_f(r.auroc)} ({_f(r.ci_low)} to {_f(r.ci_high)}) |")
    w("")
    w("Because the baseline AUROCs are themselves pushed below 0.5 by the small-sample CV bias, these deltas "
      "overstate the omics increment; the permutation p in the table above is the pre-specified test.")
    w("")
    w("S2 (feature count; AUROC + bootstrap CI only, no permutation):")
    w("")
    w("| layer | k = 10 | k = 50 (primary) | k = 200 |")
    w("|---|---|---|---|")
    for lay in LAYERS:
        cells = []
        for mid in ("l2_top10", "l2_top50", "l2_top200"):
            r = row(lay, mid)
            cells.append(f"{_f(r.auroc)} ({_f(r.ci_low)}-{_f(r.ci_high)})")
        w(f"| {lay} | " + " | ".join(cells) + " |")
    w("")
    w("## 4. Late fusion (people with every fused layer)")
    w("")
    w("| fusion | layers | n (PI / HV) | AUROC (95% CI) | permutation p | best single layer on the same people | "
      "Δ vs best single (95% CI) |")
    w("|---|---|---|---|---|---|---|")
    for r in fusion.itertuples(index=False):
        if not bool(r.analysed):
            w(f"| {r.fusion} | {r.layers} | {r.n} ({r.n_cases} / {r.n_controls}) | not analysed ({r.reason}) | | | |")
            continue
        w(f"| {r.fusion} ({r.role}) | {r.layers} | {r.n} ({r.n_cases} / {r.n_controls}) | {_f(r.auroc)} "
          f"({_f(r.ci_low)}-{_f(r.ci_high)}) | {_pfmt(r.perm_p)} | {r.best_single_layer} {_f(r.best_single_auroc)} | "
          f"{_f(r.delta_vs_best_single)} ({_f(r.delta_ci_low)} to {_f(r.delta_ci_high)}) |")
    w("")
    w("The 'best single layer' is picked after seeing the single-layer AUROCs on the fusion cohort, which favours the "
      "single layer; a fusion gain would have to overcome that.")
    w("")
    w("## 5. Univariate differential analysis (all people in each layer; BH-FDR within layer)")
    w("")
    w("| layer | cases / controls | features tested | nominal p < 0.05 | expected under null | FDR q < 0.05 "
      "(Welch) | FDR q < 0.05 (OLS + sex/age) | min q | Storey pi0 (lambda 0.5) | median abs(g) |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for _, r in usum.iterrows():
        w(f"| {r['layer']} | {r['n_cases']} / {r['n_controls']} | {r['n_features_tested']} | "
          f"{r['n_nominal_p_lt_0.05']} | {r['expected_nominal_under_null']} | {r['n_fdr_significant']} | "
          f"{r['n_ols_fdr_significant']} | {_pfmt(r['min_q_value'])} | {_f(r['pi0_storey_lambda0.5'], 2)} | "
          f"{_f(r['median_abs_hedges_g'], 2)} |")
    w("")
    w("A nominal count near the null expectation, pi0 near 1 and zero FDR hits mean the layer shows no detectable "
      "per-feature difference at this n. The median |g| of a null layer is not zero: with 10-21 per group, "
      "sampling noise alone gives |g| around 0.2-0.3.")
    w("")
    w("## 6. Agreement with the paper's published direction (reproduction check)")
    w("")
    w("| layer | paper table (scope; paper n PI / HV) | matched | paper-nominal matched | sign agreement "
      "(binomial p) | also p < 0.05 here | Spearman of effects (all matched) |")
    w("|---|---|---|---|---|---|---|")
    for _, r in conc.iterrows():
        w(f"| {r['layer']} | {r['paper_table']} ({r['paper_scope']}; {r['paper_n_cases']} / {r['paper_n_controls']}) | "
          f"{r['n_matched']} | {r['n_paper_nominal_matched']} | {_f(r['sign_agreement'])} ({_pfmt(r['binomial_p'])}) | "
          f"{r['n_paper_nominal_also_our_p_lt_0.05']} | {_f(r['spearman_effect_all_matched'])} |")
    w("")
    w("This agreement is expected and says little: the paper's analyses used largely the same people and the same "
      "deposits, and for PBMC (16A), muscle (19A) and CSF metabolomics (14A) the paper published only its own "
      "significant results, so matching their signs on overlapping people is close to circular. It confirms that the "
      "linked files, labels and orientation are read correctly; it is not replication. The SomaScan tables (17A/17B) "
      "list every aptamer, so their effect correlations are the more informative check.")
    w("")
    w("## 7. CSF vs serum SomaScan (same 42 people, same aptamers)")
    w("")
    w("| quantity | value |")
    w("|---|---|")
    for q, v in cv.items():
        w(f"| {q} | {_f(v, 4) if not float(v).is_integer() else int(v)} |")
    w("")
    w("The joint permutation (same shuffled labels in both fluids) keeps each person's CSF-serum coupling, so its null "
      "for the t-statistic correlation is centred above zero; the observed value is compared with that null, not "
      "with 0. `coupling_*` rows describe how far a person's CSF level of an aptamer tracks their serum level at all "
      "(it mostly does not), which is why the two fluids can disagree even in the same people.")
    w("")
    w("## 8. Top features")
    w("")
    if passing:
        w(f"Reported only for layers that passed ({', '.join(passing)}); `results/tables/linked_omics_top_features.csv` "
          "(top 25 by Welch p, with the number of the 50 CV folds in which the feature was selected).")
        w("")
        w("| layer | feature | name | Hedges g (95% CI) | Welch p | BH q | folds selected (of 50) |")
        w("|---|---|---|---|---|---|---|")
        for r in top.head(15).itertuples(index=False):
            w(f"| {r.layer} | {r.feature_id} | {r.feature_name} | {_f(r.hedges_g, 2)} ({_f(r.g_ci_low, 2)} to "
              f"{_f(r.g_ci_high, 2)}) | {_pfmt(r.p_value)} | {_pfmt(r.q_value_bh)} | {int(r.n_folds_selected)} |")
    else:
        w("No layer passed its pre-specified permutation test, so no feature list is reported (by design; "
          "`linked_omics_top_features.csv` is intentionally empty).")
    w("")
    w("## 9. Deviations from the plan and errata")
    w("")
    w("* Plan erratum: PBMC has 18,369 distinct genes *after* de-duplicating the 25 pseudo-autosomal genes (the plan "
      "said 18,344).")
    w(f"* CSF and serum SomaScan share {int(cv['aptamers (exact SeqId match in both runs)'])} aptamers by exact "
      f"SomaLogic SeqId; {int(cv['aptamer ids present in only one run (SeqId version differs; excluded)'])} ids (6 per "
      "run) differ only in the SeqId version suffix and were excluded from the cross-fluid comparison (the plan said "
      "'same 1,317 aptamers').")
    w("* Before the full run the pipeline was executed twice with 2-3 permutations to test the plumbing. Between those "
      "runs and the full run only two things changed: the CSF/serum aptamer matching above (the step crashed) and the "
      "in-fold residualisation was rewritten from `lstsq` to the equivalent normal equations for speed (identical "
      "AUROC). After the first full run a unit test showed that the secondary elastic net (saga solver) was not "
      "bit-reproducible because saga shuffles with an unset random_state; it was given random_state = the project seed "
      "and the whole analysis was re-run (this report). Every AUROC and permutation p (elastic net included) was "
      "identical to 6 decimals to the first full run; one bootstrap CI bound moved in the third decimal (csf_somalogic "
      "elastic net lower bound 0.399 in the first run, 0.398 here). No model, k, hyper-parameter, fold scheme or "
      "threshold was changed after seeing any result.")
    w("* Plan erratum: the plan said the 25 duplicated PBMC pseudo-autosomal genes carry identical values. In 36 of the "
      "675 duplicated person x gene pairs the two copies differ by a few read counts (< 1%); the first copy is kept, "
      "as the plan specifies. This cannot move any result materially.")
    w("* Independent review (2026-09-24). The plan's timestamp and the author's session log show the plan was written "
      "before the first model run. A full re-run reproduced every `linked_omics_*.csv` byte for byte. An "
      "independent re-implementation matched the cohort sizes, the primary L2 AUROC of all five layers (to 4 "
      "decimals), the univariate nominal and FDR counts, the CSF-serum t-statistic rho (0.036), the top-50 overlap "
      "(5) and the SomaScan paper-effect Spearman values (0.560 and 0.478).")
    w("* Five CSF metabolites have an undefined Welch p (constant, or constant within a group, after the depositors' "
      "imputation). BH-FDR is applied to the 440 defined p-values; the 'expected under null' column uses all 445.")
    w("* In the CSF-metabolomics layer sex is constant (Birth Sex = 1 only), so its baseline and S1 use age only, as "
      "the plan states.")
    w("")
    w("## 10. What this means for the device question")
    w("")
    w("Nothing directly: these are blood/CSF/tissue assays, not devices, and they cannot be linked to the study's "
      "wearable or autonomic measurements in open data. The study's per-person HRV, tilt, actigraphy and CPET files "
      "sit behind the mapMECFS login (`docs/MAP_MECFS_LINKAGE_AUDIT.md` section 8); with an account, whether they "
      "carry the NIH study numbers used here is untested. Until then the only public ME/CFS device data from this "
      "cohort are group-labelled figure source data (e.g. 24-h HRV, 33 rows, no ids), which support group "
      "comparisons only.")
    w("")
    w("## 11. Limitations")
    w("")
    w("* n = 21-42 per layer, 17-42 per fusion; a single site; wide CIs; AUROCs this small-sample are unstable across "
      "CV repeats (see repeat range).")
    w("* Adjudication status is not in the open files: 21 PI-ME/CFS are deposited in SomaScan versus 17 in the "
      "paper's analytic cohort, so the label may include people the paper excluded. Label noise would bias towards "
      "the null.")
    w("* Serum SomaScan is linked through the Fatigue-NNN crosswalk published only in the CSF deposit (corroborated by "
      "sex/age/group, not by an id in the serum deposit); Fatigue-124's serum and CSF came from separate visits.")
    w("* CSF metabolomics: single-sex, batch-normalised and imputed values from figure source data; the batch structure "
      "is not published, so batch-group confounding cannot be excluded; the PI-ME/CFS participants are younger.")
    w("* MECFS_311 is Female in SomaScan and Male in RNA-seq; sex is used as recorded in each layer.")
    w("* Univariate agreement with the paper is not independent (same people).")
    w("")
    w("## 12. Files")
    w("")
    w("* Tables: `results/tables/linked_omics_{cohort,model_performance,vs_baseline,sensitivity,fusion,"
      "permutation_null,univariate,univariate_summary,paper_concordance,paper_concordance_features,"
      "csf_serum_consistency,csf_serum_per_aptamer,selection_frequency,top_features}.csv`, "
      "`linked_omics_run_metadata.json`.")
    w("* Figures: `results/figures/drafts/linked_omics_auroc_forest.png`, `linked_omics_csf_vs_serum.png`.")
    w(f"* Processed: `data/processed/{SIGNATURE_TABLE}.parquet` (object ids "
      f"`signature:{DATASET_ID}|pi_mecfs_vs_hv__<layer>|<feature>`; model rows `model_auroc:<layer>:<model>`), read by "
      "`get_patient_phenotype_signature`.")
    w(f"* Reproduce: `uv run python -m {PRODUCER}` (about {int(meta.get('runtime_s', 0) // 60)} min with "
      f"{meta.get('workers')} workers); `--report-only` rebuilds this file, the figures and the partition from the "
      "tables.")
    w("")
    (config.RESULTS / "LINKED_OMICS_RESULTS.md").write_text("\n".join(L))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args(argv)
    if a.report_only:
        finish()
    else:
        run(n_perm=a.n_perm, workers=a.workers)


if __name__ == "__main__":
    main()
