"""Part 2 — Appelman 2024 Source Data: plasma and muscle metabolomics, Long COVID vs healthy.

Follows docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md (Part 2), locked by `measure_it.labs.plan`; refuses to run if the
plan was edited. Reads the file registered by `measure_it.ingestion.appelman_lc_pem`.

    uv run python -m measure_it.labs.appelman_metabolomics
"""
from __future__ import annotations

import json
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from ..config import SEED, TABLES, utc_now_iso
from ..ingestion import appelman_lc_pem as src
from .labs_vs_diagnosis import bh_fdr, ols_hc3
from .plan import check_plan_locked

PRODUCER = "measure_it.labs.appelman_metabolomics"
DATASET_ID = "appelman_lc_pem_source"
PREFIX = "labs_dx_appelman"
N_SPLITS, N_REPEATS = 5, 50
N_BOOT, N_PERM = 2000, 1000
TOP_K = 10
C_FIXED = 1.0
CS_INNER = (0.01, 0.1, 1.0)
FDR_ALPHA = 0.05

# Published claims (baseline, LC vs healthy), recorded in the plan before computation. sign: -1 lower, +1 higher,
# 0 "no/few differences".
CLAIMS = [
    ("muscle", "metabolite", "Glutamate", -1, "TCA-cycle metabolites incl. glutamate lower"),
    ("muscle", "metabolite", "FAD", -1, "TCA-cycle metabolites incl. FAD+ lower"),
    ("muscle", "metabolite", "Alpha-Ketoglutarate", -1, "TCA-cycle metabolites incl. alpha-ketoglutarate lower"),
    ("muscle", "metabolite", "Citric acid", -1, "TCA-cycle metabolites incl. citric acid lower"),
    ("muscle", "domain", "TCA", -1, "key metabolites of the TCA cycle lower"),
    ("muscle", "domain", "Glycolysis", 0, "glycolytic metabolites displayed few group differences"),
    ("muscle", "metabolite", "Creatine", -1, "creatine concentrations lower"),
    ("muscle", "metabolite", "S-adenosyl methionine", -1, "lower S-adenosylmethionine (SAM)"),
    ("muscle", "metabolite", "Hydroxyphenyllactic acid", -1,
     "lower hydroxyphenyl acetic acid (NAME MISMATCH: file has hydroxyphenyllactic acid)"),
    ("muscle", "ratio", "Citric acid/Lactate", -1, "ratio of citric acid to lactate lower"),
    ("muscle", "domain", "Aminoacid", 0, "many amino acids not different at rest"),
    ("blood", "domain", "Glycolysis", 1, "glycolytic metabolites in venous blood significantly higher"),
    ("blood", "metabolite", "Pyruvate", -1, "pyruvate lower"),
    ("blood", "domain", "TCA", -1, "other TCA-cycle metabolites lower"),
    ("blood", "domain", "Purines", -1, "various purine-pathway metabolites lower"),
]


# ======================================================================================================== data
def baseline(tissue: str) -> tuple[pd.DataFrame, list[str]]:
    d, mets = src.load_metabolomics(tissue)
    b = d[d["Time"] == "Baseline"].copy()
    demo = src.load_demographics()[["ID", "female", "bmi", "long_covid"]]
    b = b.merge(demo, on="ID", how="left", validate="one_to_one")
    assert (b["long_covid"] == b["Group"].map(src.GROUPS)).all()
    b[mets] = np.log2(b[mets].astype(float))
    return b.sort_values("ID").reset_index(drop=True), mets


def pem_change(tissue: str) -> tuple[pd.DataFrame, list[str]]:
    d, mets = src.load_metabolomics(tissue)
    b = d[d["Time"] == "Baseline"].set_index("ID")
    p = d[d["Time"] == "1-day after PEM"].set_index("ID")
    ids = sorted(set(b.index) & set(p.index))
    ch = np.log2(p.loc[ids, mets].astype(float)) - np.log2(b.loc[ids, mets].astype(float))
    ch["long_covid"] = b.loc[ids, "Group"].map(src.GROUPS)
    return ch.reset_index(), mets


# ======================================================================================================== stats
def hedges_g(x1: np.ndarray, x0: np.ndarray) -> float:
    n1, n0 = len(x1), len(x0)
    sp = np.sqrt(((n1 - 1) * x1.var(ddof=1) + (n0 - 1) * x0.var(ddof=1)) / (n1 + n0 - 2))
    if not np.isfinite(sp) or sp == 0:
        return np.nan
    j = 1 - 3 / (4 * (n1 + n0) - 9)
    return float(j * (x1.mean() - x0.mean()) / sp)


def strat_boot_idx(y: np.ndarray, n_boot: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    i1, i0 = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    return np.hstack([rng.choice(i1, (n_boot, len(i1))), rng.choice(i0, (n_boot, len(i0)))])


def strat_boot_weights(y: np.ndarray, n_boot: int, seed: int) -> np.ndarray:
    idx = strat_boot_idx(y, n_boot, seed)
    W = np.zeros((n_boot, len(y)), dtype=float)
    for b in range(n_boot):
        np.add.at(W[b], idx[b], 1.0)
    return W


def per_metabolite(df: pd.DataFrame, mets: list[str], y: np.ndarray, domains: dict, seed: int,
                   covariates: bool = True) -> pd.DataFrame:
    idx = strat_boot_idx(y, N_BOOT, seed)
    rows = []
    for m in mets:
        x = df[m].to_numpy(float)
        g = hedges_g(x[y == 1], x[y == 0])
        gb = np.array([hedges_g(x[ii][y[ii] == 1], x[ii][y[ii] == 0]) for ii in idx])
        t = stats.ttest_ind(x[y == 1], x[y == 0], equal_var=False)
        mw = stats.mannwhitneyu(x[y == 1], x[y == 0], alternative="two-sided")
        r = {"metabolite": m, "domain": domains.get(m, "unassigned"), "n_lc": int((y == 1).sum()),
             "n_healthy": int((y == 0).sum()), "mean_log2_lc": float(x[y == 1].mean()),
             "mean_log2_healthy": float(x[y == 0].mean()), "log2_fold_change": float(x[y == 1].mean() - x[y == 0].mean()),
             "hedges_g": g, "g_ci_low": float(np.nanquantile(gb, 0.025)), "g_ci_high": float(np.nanquantile(gb, 0.975)),
             "welch_p": float(t.pvalue), "mannwhitney_p": float(mw.pvalue)}
        if covariates:
            z = (x - x.mean()) / x.std(ddof=1)
            female = df["female"].to_numpy(float)
            b, se, p = ols_hc3(z, np.column_stack([np.ones_like(z), y, female]))
            r.update({"sex_adj_std_diff": b, "sex_adj_ci_low": b - 1.959964 * se, "sex_adj_ci_high": b + 1.959964 * se,
                      "sex_adj_p": p})
            bmi = df["bmi"].to_numpy(float)
            b2, se2, p2 = ols_hc3(z, np.column_stack([np.ones_like(z), y, female, bmi]))
            r.update({"sex_bmi_adj_std_diff": b2, "sex_bmi_adj_ci_low": b2 - 1.959964 * se2,
                      "sex_bmi_adj_ci_high": b2 + 1.959964 * se2, "sex_bmi_adj_p": p2})
        rows.append(r)
    out = pd.DataFrame(rows)
    out["q_welch_bh_tissue"] = bh_fdr(out["welch_p"].to_numpy())
    out["q_mannwhitney_bh_tissue"] = bh_fdr(out["mannwhitney_p"].to_numpy())
    if covariates:
        out["q_sex_adj_bh_tissue"] = bh_fdr(out["sex_adj_p"].to_numpy())
    out["q_welch_bh_within_domain_reproduction"] = np.nan
    for _, g in out.groupby("domain"):
        out.loc[g.index, "q_welch_bh_within_domain_reproduction"] = bh_fdr(g["welch_p"].to_numpy())
    out["fdr_significant_tissue"] = out["q_welch_bh_tissue"] < FDR_ALPHA
    return out


# ======================================================================================================== models
def cv_folds(y: np.ndarray, repeat: int):
    from sklearn.model_selection import StratifiedKFold
    return list(StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED + repeat).split(np.zeros(len(y)), y))


def _welch_t(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    a, b = X[y == 1], X[y == 0]
    va, vb = a.var(0, ddof=1) / len(a), b.var(0, ddof=1) / len(b)
    return (a.mean(0) - b.mean(0)) / np.sqrt(va + vb + 1e-12)


def fit_predict(Xm: np.ndarray | None, Xf: np.ndarray | None, y: np.ndarray, tr: np.ndarray, te: np.ndarray,
                method: str, seed: int) -> np.ndarray:
    """Xm: metabolite block (selection applies), Xf: forced covariates (sex / VO2max; never selected out).
    method: 'top10' (in-fold |Welch t| top-10 + L2 C=1) or 'all_l2cv' (all metabolites, C by inner 3-fold CV)."""
    from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
    from sklearn.model_selection import StratifiedKFold
    blocks_tr, blocks_te = [], []
    if Xf is not None:
        blocks_tr.append(Xf[tr])
        blocks_te.append(Xf[te])
    if Xm is not None:
        mu, sd = Xm[tr].mean(0), Xm[tr].std(0, ddof=1)
        sd[sd == 0] = 1.0
        Ztr, Zte = (Xm[tr] - mu) / sd, (Xm[te] - mu) / sd
        if method == "top10":
            keep = np.argsort(-np.abs(_welch_t(Ztr, y[tr])))[:TOP_K]
            Ztr, Zte = Ztr[:, keep], Zte[:, keep]
        blocks_tr.append(Ztr)
        blocks_te.append(Zte)
    A, B = np.column_stack(blocks_tr), np.column_stack(blocks_te)
    if Xf is not None:  # standardise forced covariates too (fitted on train)
        k = Xf.shape[1]
        mu, sd = A[:, :k].mean(0), A[:, :k].std(0, ddof=1)
        sd[sd == 0] = 1.0
        A[:, :k], B[:, :k] = (A[:, :k] - mu) / sd, (B[:, :k] - mu) / sd
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if method == "all_l2cv" and Xm is not None:
            m = LogisticRegressionCV(Cs=list(CS_INNER), cv=StratifiedKFold(3, shuffle=True, random_state=seed),
                                     scoring="neg_log_loss", solver="lbfgs", max_iter=5000)
        else:
            m = LogisticRegression(C=C_FIXED, solver="lbfgs", max_iter=5000)
        m.fit(A, y[tr])
    return m.predict_proba(B)[:, 1]


def oof(Xm, Xf, y, repeat: int, method: str) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    for k, (tr, te) in enumerate(cv_folds(y, repeat)):
        pred[te] = fit_predict(Xm, Xf, y, tr, te, method, SEED + 100 * repeat + k)
    return pred


def run_model(name: str, Xm, Xf, y, method: str = "top10", n_perm: int = 0, n_jobs: int = 16) -> dict:
    from joblib import Parallel, delayed

    from ..wearables.models import rank_metrics
    P = np.vstack(Parallel(n_jobs=n_jobs)(delayed(oof)(Xm, Xf, y, r, method) for r in range(N_REPEATS)))
    W = strat_boot_weights(y, N_BOOT, SEED + 17)
    per_rep = np.array([rank_metrics(P[r], y)[0][0] for r in range(N_REPEATS)])
    ap_rep = np.array([rank_metrics(P[r], y)[1][0] for r in range(N_REPEATS)])
    boot = np.mean([rank_metrics(P[r], y, W)[0] for r in range(N_REPEATS)], axis=0)
    res = {"model": name, "method": method, "n_lc": int(y.sum()), "n_healthy": int((1 - y).sum()),
           "auroc": float(per_rep.mean()), "auroc_ci_low": float(np.quantile(boot, 0.025)),
           "auroc_ci_high": float(np.quantile(boot, 0.975)), "auroc_repeat_p2_5": float(np.quantile(per_rep, 0.025)),
           "auroc_repeat_p97_5": float(np.quantile(per_rep, 0.975)), "auprc": float(ap_rep.mean()),
           "auprc_no_skill": float(y.mean()), "_P": P, "_boot": boot}
    if n_perm:
        def _perm(i):
            yp = np.random.default_rng(SEED + 7919 * (i + 1)).permutation(y)
            return float(rank_metrics(oof(Xm, Xf, yp, 0, method), yp)[0][0])
        null = np.array(Parallel(n_jobs=n_jobs)(delayed(_perm)(i) for i in range(n_perm)))
        obs = float(rank_metrics(P[0], y)[0][0])
        res.update({"perm_observed_repeat1": obs, "perm_null_mean": float(null.mean()),
                    "perm_null_p95": float(np.quantile(null, 0.95)),
                    "perm_p": float((1 + (null >= obs).sum()) / (1 + n_perm)), "n_permutations": n_perm})
    return res


def decision(r: dict) -> str:
    if r.get("perm_p") is None:
        return ""
    if r["auroc_ci_low"] > 0.5 and r["perm_p"] <= 0.05:
        return "strong separation" if r["auroc_ci_low"] >= 0.70 else "metabolites separate LC from healthy"
    return "null"


def delta(a: dict, b: dict) -> dict:
    d = a["_boot"] - b["_boot"]
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return {"comparison": f"{a['model']} minus {b['model']}", "delta_auroc": a["auroc"] - b["auroc"],
            "ci_low": float(np.quantile(d, 0.025)), "ci_high": float(np.quantile(d, 0.975)),
            "boot_p": float(max(p, 1 / len(d)))}


# ======================================================================================================== paper check
def direction_check(tissue: str, df: pd.DataFrame, mets: list[str], y: np.ndarray, eff: pd.DataFrame,
                    domains: dict) -> list[dict]:
    rows = []
    idx = strat_boot_idx(y, N_BOOT, SEED + 3)
    Z = (df[mets] - df[mets].mean()) / df[mets].std(ddof=1)
    for t, kind, name, sign, claim in CLAIMS:
        if t != tissue:
            continue
        r = {"tissue": tissue, "kind": kind, "item": name, "claimed_sign": sign, "published_claim": claim}
        if kind == "metabolite":
            e = eff.set_index("metabolite").loc[name]
            g, lo, hi, q = e.hedges_g, e.g_ci_low, e.g_ci_high, e.q_welch_bh_tissue
            r.update({"n_in_domain": 1, "share_claimed_sign": float(np.sign(g) == sign) if sign else np.nan})
        else:
            if kind == "ratio":
                a, b = name.split("/")
                score = (df[a] - df[b]).to_numpy(float)
                members = [a, b]
            else:
                members = [m for m in mets if domains.get(m) == name]
                score = Z[members].mean(1).to_numpy(float)
            g = hedges_g(score[y == 1], score[y == 0])
            gb = np.array([hedges_g(score[ii][y[ii] == 1], score[ii][y[ii] == 0]) for ii in idx])
            lo, hi = float(np.nanquantile(gb, 0.025)), float(np.nanquantile(gb, 0.975))
            q = np.nan
            sub = eff[eff["metabolite"].isin(members)]
            r.update({"n_in_domain": len(members),
                      "share_claimed_sign": float((np.sign(sub["hedges_g"]) == sign).mean()) if sign and kind == "domain"
                      else np.nan,
                      "n_domain_fdr_tissue": int(sub["fdr_significant_tissue"].sum()) if kind == "domain" else np.nan,
                      "n_domain_q_within_domain_lt_0_05": int((sub["q_welch_bh_within_domain_reproduction"] < 0.05).sum())
                      if kind == "domain" else np.nan})
        r.update({"our_hedges_g": g, "ci_low": lo, "ci_high": hi, "q_tissue": q})
        if sign == 0:
            r["concordance"] = "concordant (CI includes 0)" if lo <= 0 <= hi else "discordant (difference found)"
        else:
            same = np.sign(g) == sign
            ci_ex = (hi < 0) if sign < 0 else (lo > 0)
            opp = (lo > 0) if sign < 0 else (hi < 0)
            r["concordance"] = ("same direction, CI excludes 0" if same and ci_ex else
                                "same direction, CI includes 0" if same else
                                "opposite direction, CI excludes 0" if opp else "opposite direction, CI includes 0")
        rows.append(r)
    return rows


# ======================================================================================================== signatures
def build_signatures(eff_all: pd.DataFrame, models: pd.DataFrame, meta: dict) -> pd.DataFrame:
    from ..provenance import add_provenance
    from ..store import write_table
    from .labs_vs_diagnosis import SIG_COLUMNS
    caveat = ("Appelman 2024 Source Data (NCT05225688): 25 Long COVID vs healthy people recovered from mild COVID-19; "
              "healthy controls only; age not released (sex/BMI adjustment only); relative normalised LC-MS "
              "intensities; same trial as Charlton 2026 (not independent). Group difference, not a diagnostic test.")
    definition = ("Long COVID (>= 6 months, PEM by DSQ-PEM, aged 18-65) vs age- and sex-matched healthy controls who "
                  "recovered from mild COVID-19 (per paper)")
    rows = []
    for _, e in eff_all.iterrows():
        pid = f"long_covid_vs_healthy__{e['tissue']}_{e['analysis']}"
        rows.append({
            "phenotype_id": pid,
            "phenotype_label": f"Long COVID vs healthy: {e['tissue']} metabolomics ({e['analysis'].replace('_', ' ')})",
            "condition_id": "long_covid", "is_proxy": False, "label_basis": "clinical case definition (study enrolment)",
            "feature": e["metabolite"], "feature_label": f"{e['metabolite']} [{e['domain']}]",
            "effect_measure": "hedges_g", "effect_size": e["hedges_g"], "ci_low": e["g_ci_low"],
            "ci_high": e["g_ci_high"], "null_value": 0.0, "se": np.nan, "p_value": e["welch_p"],
            "q_value_bh": e["q_welch_bh_tissue"], "fdr_significant": bool(e["q_welch_bh_tissue"] < FDR_ALPHA),
            "direction": "higher_in_cases" if e["hedges_g"] > 0 else "lower_in_cases",
            "n_cases": e["n_lc"], "n_controls": e["n_healthy"], "underpowered": True,
            "phenotype_definition": definition, "cycles": "", "feature_units": "log2 relative intensity" +
            (" change (1 day after PEM minus baseline)" if e["analysis"] == "pem_change" else ""),
            "feature_transform": "log2", "covariates": "none (sex-adjusted estimate in labs_dx_appelman_metabolites.csv)",
            "mean_cases": e["mean_log2_lc"], "mean_controls": e["mean_log2_healthy"], "caveats": caveat,
            "effect_unit": "Hedges g (LC minus healthy)", "q_value": e["q_welch_bh_tissue"],
            "method": "Hedges g with stratified bootstrap CI (2,000); Welch t p; BH across the tissue's metabolites"})
    for _, m in models.iterrows():
        pid = f"long_covid_vs_healthy__{m['tissue']}_baseline"
        rows.append({
            "phenotype_id": pid, "phenotype_label": f"Long COVID vs healthy: {m['tissue']} metabolomics (baseline)",
            "condition_id": "long_covid", "is_proxy": False, "label_basis": "clinical case definition (study enrolment)",
            "feature": f"model_auroc:{m['model']}", "feature_label": f"cross-validated AUROC, {m['model']} ({m['method']})",
            "effect_measure": "cv_auroc", "effect_size": m["auroc"], "ci_low": m["auroc_ci_low"],
            "ci_high": m["auroc_ci_high"], "null_value": 0.5, "se": np.nan, "p_value": m.get("perm_p", np.nan),
            "q_value_bh": np.nan, "fdr_significant": False, "direction": str(m.get("decision", "") or ""),
            "n_cases": m["n_lc"], "n_controls": m["n_healthy"], "underpowered": True,
            "phenotype_definition": definition, "cycles": "", "feature_units": "AUROC", "feature_transform": "",
            "covariates": "", "mean_cases": np.nan, "mean_controls": np.nan, "caveats": caveat, "effect_unit": "AUROC",
            "q_value": np.nan,
            "method": "stratified 5-fold x 50 CV, in-fold selection; stratified person bootstrap CI (2,000); "
                      "p = label-permutation null (1,000) where run"})
    sig = pd.DataFrame(rows)
    sig["dataset_id"] = DATASET_ID
    sig["population"] = "Appelman 2024 PEM study participants with the sample at the stated timepoint"
    sig["object_id"] = "signature:" + DATASET_ID + "|" + sig["phenotype_id"] + "|" + sig["feature"]
    sig = sig[SIG_COLUMNS]
    sig = add_provenance(
        sig, data_layer="person", source_name="Appelman et al. 2024 Nat Commun 15:17 Source Data",
        source_version=meta["source_version"], retrieved_at=meta["retrieved_at"],
        evidence_type="person_lab_measurement", source_record_id="object_id",
        evidence_level="study case definition; participant-linked metabolomics (within-file ids)",
        provenance_notes=("Our own computation from the released Source Data (plan docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS"
                          ".md, sha256 " + meta["plan_sha256"][:16] + "); published claims are compared separately "
                          "in labs_dx_appelman_paper_direction.csv."))
    assert sig["object_id"].is_unique
    write_table(sig, "phenotype_signatures__appelman_metabolomics", producer=PRODUCER,
                description="Appelman 2024 plasma/muscle metabolomics, Long COVID vs healthy: per-metabolite Hedges g "
                            "and cross-validated model rows")
    return sig


# ======================================================================================================== run
def run(n_jobs: int = 16, log=print) -> dict:
    plan_sha = check_plan_locked()
    t0 = time.time()
    man = src.load_manifest(src.SOURCE_ID)["files"][src.FILENAME]
    meta = {"computed_at": utc_now_iso(), "plan_sha256": plan_sha, "seed": SEED,
            "source_version": f"{src.FILENAME} sha256 {man['sha256']} (Last-Modified {man.get('last_modified', '')})",
            "retrieved_at": man["retrieved_at"], "n_splits": N_SPLITS, "n_repeats": N_REPEATS, "n_boot": N_BOOT,
            "n_perm": N_PERM, "top_k": TOP_K}
    effs, pem_effs, models, deltas, direction = [], [], [], [], []
    data = {}
    for tissue in ("blood", "muscle"):
        df, mets = baseline(tissue)
        dom = src.load_domains(tissue).set_index("Metabolite")["Domain"].to_dict()
        y = df["long_covid"].to_numpy(int)
        data[tissue] = (df, mets, y)
        eff = per_metabolite(df, mets, y, dom, SEED + (1 if tissue == "blood" else 2))
        eff.insert(0, "tissue", tissue)
        eff.insert(1, "analysis", "baseline")
        effs.append(eff)
        direction += direction_check(tissue, df, mets, y, eff, dom)
        log(f"{tissue}: per-metabolite done ({time.time() - t0:.0f}s)")
        Xm, Xs = df[mets].to_numpy(float), df[["female"]].to_numpy(float)
        prim = run_model(f"{tissue}_top10", Xm, None, y, "top10", N_PERM, n_jobs)
        prim["decision"] = decision(prim)
        sexm = run_model(f"{tissue}_sex_only", None, Xs, y, "top10", 0, n_jobs)
        sexmet = run_model(f"{tissue}_sex_plus_top10", Xm, Xs, y, "top10", 0, n_jobs)
        allm = run_model(f"{tissue}_all_l2cv", Xm, None, y, "all_l2cv", 0, n_jobs)
        for r in (prim, sexm, sexmet, allm):
            r["tissue"] = tissue
            r["role"] = "primary" if r is prim else "secondary"
            models.append(r)
        deltas.append({"tissue": tissue, **delta(sexmet, sexm)})
        log(f"{tissue}: models done ({time.time() - t0:.0f}s) primary AUROC {prim['auroc']:.3f}")
        ch, cmets = pem_change(tissue)
        yc = ch["long_covid"].to_numpy(int)
        pe = per_metabolite(ch, cmets, yc, dom, SEED + 5, covariates=False)
        pe.insert(0, "tissue", tissue)
        pe.insert(1, "analysis", "pem_change")
        pem_effs.append(pe)
    # (d) blood + muscle combined
    bdf, bm, _ = data["blood"]
    mdf, mm, _ = data["muscle"]
    comb = bdf[["ID", "long_covid", "female"] + bm].rename(columns={m: f"blood:{m}" for m in bm}).merge(
        mdf[["ID"] + mm].rename(columns={m: f"muscle:{m}" for m in mm}), on="ID")
    cm = [c for c in comb.columns if ":" in c]
    yb = comb["long_covid"].to_numpy(int)
    r = run_model("blood_plus_muscle_top10", comb[cm].to_numpy(float), None, yb, "top10", N_PERM, n_jobs)
    r.update({"tissue": "blood+muscle", "role": "secondary"})
    r["decision"] = decision(r)
    models.append(r)
    # (e) VO2max alone vs VO2max + blood metabolites
    vo = src.load_vo2().dropna(subset=["VO2max"])
    vb = bdf.merge(vo, on="ID")
    yv = vb["long_covid"].to_numpy(int)
    rv = run_model("vo2max_only", None, vb[["VO2max"]].to_numpy(float), yv, "top10", 0, n_jobs)
    rvb = run_model("vo2max_plus_blood_top10", vb[bm].to_numpy(float), vb[["VO2max"]].to_numpy(float), yv, "top10",
                    0, n_jobs)
    for x in (rv, rvb):
        x.update({"tissue": "blood (people with VO2max)", "role": "secondary"})
        models.append(x)
    deltas.append({"tissue": "blood (people with VO2max)", **delta(rvb, rv)})
    log(f"combined/VO2 done ({time.time() - t0:.0f}s)")

    eff_all = pd.concat(effs + pem_effs, ignore_index=True)
    eff_all.to_csv(TABLES / f"{PREFIX}_metabolites.csv", index=False)
    mdf_out = pd.DataFrame([{k: v for k, v in m.items() if not k.startswith("_")} for m in models])
    mdf_out.to_csv(TABLES / f"{PREFIX}_models.csv", index=False)
    pd.DataFrame(deltas).to_csv(TABLES / f"{PREFIX}_deltas.csv", index=False)
    pd.DataFrame(direction).to_csv(TABLES / f"{PREFIX}_paper_direction.csv", index=False)
    build_signatures(eff_all, mdf_out[mdf_out["tissue"].isin(["blood", "muscle"])], meta)
    meta["seconds"] = round(time.time() - t0)
    (TABLES / f"{PREFIX}_run_metadata.json").write_text(json.dumps(meta, indent=2, default=str))
    log(mdf_out[["model", "auroc", "auroc_ci_low", "auroc_ci_high", "perm_p"]].to_string())
    return meta


if __name__ == "__main__":
    run()
