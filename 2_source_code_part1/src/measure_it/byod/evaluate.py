"""`measure-it byod evaluate <dir>`: lock the manifest's primary analysis, then run it and nothing else.

Order (enforced):
 1. the plan (label, condition, comparator, feature blocks, covariates, operating point, CV / bootstrap / permutation
    settings, model, decision rule, engine version) is written to results/byod/<id>/PLAN.md and its sha256 to
    plan_lock.json BEFORE any outcome is computed (as measure_it.labs.plan does for the lab analyses). A later run
    with a different plan is refused unless --amend is given; an amendment is recorded, and the record says so.
    Every lock and evaluation is appended to results/byod/LEDGER.jsonl, which `byod remove` keeps, so a dataset that
    is removed and evaluated again under another plan shows its earlier evaluations.
 2. models (participant-level repeated stratified k-fold CV; L2 logistic regression; median imputation and scaling
    fitted inside each training fold; lab_cc_engine conventions):
      primary       the pre-specified feature blocks          -> the performance record
      combined      the pre-specified combined blocks (optional)
      <block>       each supplied block alone                 (secondary, descriptive)
      covariates    age band + sex only (when supplied)       (baseline)
      primary+covariates vs covariates: paired Delta AUROC   (secondary)
    AUROC = mean over repeats of the pooled out-of-fold AUROC; 95% CI from a participant bootstrap (resampling people
    within class, out-of-fold predictions held fixed); AUPRC of the repeat-averaged out-of-fold probabilities (same
    bootstrap); sensitivity at the manifest's specificity (threshold on the controls' repeat-averaged out-of-fold
    probabilities, re-chosen in every bootstrap resample; metric_link.sens_at_spec_boot); label-permutation null
    (labels shuffled, CV re-run with one repeat per shuffle; one-sided p = (1 + #null >= observed) / (1 + n)).
 3. outputs: phenotype_signatures__byod_<id> (model rows + per-feature effects), one record in
    byod_performance_records (tier user_supplied_own_computation), results/byod/<id>/EVALUATION.md.

Decision rule (fixed): 'supported' = AUROC 95% CI lower bound > 0.5 AND permutation p <= 0.05; otherwise 'null'.
Against healthy controls every number is optimistic (the clinical question is the illness versus look-alikes).
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import SEED, UNKNOWN, load_config, utc_now_iso
from ..provenance import add_provenance
from ..store import read_table, table_exists, write_table
from . import common as K
from . import person as PER

PRODUCER = "measure_it.byod.evaluate"
EVAL_VERSION = "byod_evaluate 2026-09-28.1"
N_PERM_SECONDARY = 100
METRIC_LINK_SPEC = 0.90     # the metric link's operating point (docs/ANALYSIS_PLAN_METRIC_LINK.md section 2.4)
DEVICE_FAMILIES = {"wearable", "clinical_test", "imaging", "digital"}


# ------------------------------------------------------------------------------------------------------------------
# plan + lock
# ------------------------------------------------------------------------------------------------------------------

def plan_document(manifest: dict, meta: dict) -> dict:
    pa = K.analysis_defaults(manifest["primary_analysis"])
    labels = json.loads(meta["labels_json"])
    lab = labels[pa["label"]]
    return {
        "dataset_id": meta["dataset_id"], "engine_version": EVAL_VERSION, "seed": SEED,
        "label": pa["label"], "condition_id": lab["condition_id"], "label_definition": lab["definition"],
        "label_basis": lab["label_basis"], "comparator": manifest["comparator"],
        "measurement_class": manifest["measurement"]["class"], "feature_blocks": list(pa["feature_blocks"]),
        "combined_blocks": list(pa["combined_blocks"]), "covariates": list(pa["covariates"]),
        "specificity": float(pa["specificity"]), "metric_link_specificity": METRIC_LINK_SPEC,
        "cv": {"scheme": "participant-level stratified k-fold, shuffled, seed = SEED + repeat",
               "n_splits": int(pa["n_splits"]), "n_repeats": int(pa["n_repeats"])},
        "model": f"L2 logistic regression, C = {pa['C']}, median imputation + standardisation fitted inside each "
                 "training fold (measure_it.labs.lab_cc_engine.make_model)",
        "n_bootstrap": int(pa["n_bootstrap"]), "n_permutations_primary": int(pa["n_permutations"]),
        "n_permutations_secondary": min(N_PERM_SECONDARY, int(pa["n_permutations"])),
        "secondary": ["each supplied block alone (participants with data in that block)",
                      "covariate baseline (age band midpoint + female indicator) when supplied",
                      "primary + covariates vs covariates: paired Delta AUROC (same folds and resamples)"],
        "label_defining_exclusion": "ICD-10-CM categories that map to the label's condition are removed from models",
        "decision_rule": "supported = AUROC 95% CI lower bound > 0.5 AND permutation p <= 0.05; otherwise null",
        "record": "the primary model only enters byod_performance_records (tier user_supplied_own_computation)",
    }


def plan_markdown(plan: dict, sha: str, locked_at: str, amends: str | None) -> str:
    cmd = plan.get("command", "evaluate")
    lines = [f"# Locked {'subgroup ' if cmd == 'subgroup' else ''}analysis plan: {plan['dataset_id']}", "",
             f"Locked {locked_at} (sha256 `{sha}`) by `measure-it byod {cmd}`, before any outcome was computed."
             + (f" Amends plan `{amends}` (recorded in results/byod/LEDGER.jsonl)." if amends else ""), "",
             "```json", json.dumps(plan, indent=1, sort_keys=True), "```", ""]
    return "\n".join(lines)


def _ledger(event: dict) -> None:
    K.local_only_dir(K.RESULTS_ROOT)
    with open(K.LEDGER, "a") as fh:
        fh.write(json.dumps({**event, "at": utc_now_iso()}) + "\n")


def ledger_events(dataset_id: str) -> list[dict]:
    if not K.LEDGER.exists():
        return []
    ds = K.engine_id(dataset_id)
    out = []
    for line in K.LEDGER.read_text().splitlines():
        try:
            e = json.loads(line)
        except ValueError:
            continue
        if e.get("dataset_id") == ds:
            out.append(e)
    return out


def lock_plan(manifest: dict, meta: dict, amend: bool = False) -> dict:
    return lock_document(meta["dataset_id"], plan_document(manifest, meta), amend=amend)


def lock_document(ds: str, plan: dict, amend: bool = False, *, lock_name: str = "plan_lock.json",
                  md_name: str = "PLAN.md", event: str = "plan", what: str = "analysis plan") -> dict:
    """Hash `plan` (sha256) and write it before any outcome is computed; refuse a different plan unless amend.
    Used by `evaluate` (plan_lock.json / PLAN.md) and `subgroup` (subgroup_plan_lock.json / SUBGROUP_PLAN.md)."""
    out_dir = K.local_only_dir(K.results_path(ds))
    sha = K.sha256_json(plan)
    lock_path = out_dir / lock_name
    if lock_path.exists():
        old = json.loads(lock_path.read_text())
        if old["sha256"] == sha:
            return {**old, "status": "verified"}
        if not amend:
            raise RuntimeError(
                f"the {what} of {ds} was locked on {old['locked_at']} (sha256 {old['sha256'][:12]}) and the "
                "manifest now describes a different plan. Run under the locked plan, or re-run with --amend to "
                "record an amendment (the outputs will say the plan was amended).")
        stem = lock_name.rsplit(".", 1)[0]
        (out_dir / f"{stem}.superseded.{old['sha256'][:12]}.json").write_text(json.dumps(old, indent=1))
        amends = old["sha256"]
    else:
        amends = None
    now = utc_now_iso()
    lock = {"dataset_id": ds, "sha256": sha, "locked_at": now, "plan": plan, "amends": amends}
    lock_path.write_text(json.dumps(lock, indent=1))
    (out_dir / md_name).write_text(plan_markdown(plan, sha, now, amends))
    _ledger({"event": f"{event}_amended" if amends else f"{event}_locked", "dataset_id": ds, "plan_sha256": sha,
             "amends": amends})
    return {**lock, "status": "amended" if amends else "locked"}


# ------------------------------------------------------------------------------------------------------------------
# statistics (lab_cc_engine conventions, with the plan's k and repeats)
# ------------------------------------------------------------------------------------------------------------------

def cv_oof(X: np.ndarray, y: np.ndarray, C: float, n_splits: int, n_repeats: int, seed: int = SEED,
           model_factory=None) -> np.ndarray:
    """Out-of-fold P(y=1), shape (n_repeats, n). `model_factory(C)` builds the model (default: the L2 logistic
    pipeline of measure_it.labs.lab_cc_engine.make_model; `byod subgroup` passes its L1 model)."""
    from sklearn.model_selection import StratifiedKFold

    from ..labs.lab_cc_engine import make_model
    make = model_factory or make_model
    X, y = np.asarray(X, float), np.asarray(y).astype(int)
    P = np.full((n_repeats, len(y)), np.nan)
    for r in range(n_repeats):
        for tr, te in StratifiedKFold(n_splits, shuffle=True, random_state=seed + r).split(np.zeros((len(y), 1)), y):
            m = make(C)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m.fit(X[tr], y[tr])
            P[r, te] = m.predict_proba(X[te])[:, 1]
    return P


def _perm_one(X, y, C, n_splits, i, seed, model_factory=None):
    from ..wearables.stanford_models import auroc_rows
    yp = np.random.default_rng([seed, i]).permutation(y)
    return float(auroc_rows(yp, cv_oof(X, yp, C, n_splits, 1, seed, model_factory)).mean())


def permutation_null(X, y, C, n_splits, n_perm, jobs, seed=SEED, model_factory=None) -> np.ndarray:
    if n_perm <= 0:
        return np.array([])
    from joblib import Parallel, delayed
    return np.array(Parallel(n_jobs=jobs)(delayed(_perm_one)(X, y, C, n_splits, i, seed, model_factory)
                                          for i in range(n_perm)))


def summarize(y: np.ndarray, P: np.ndarray, boot_idx: np.ndarray, spec: float, n_boot: int) -> tuple[dict, np.ndarray]:
    from sklearn.metrics import average_precision_score

    from ..scoring.metric_link import sens_at_spec_boot
    from ..wearables.stanford_models import auroc_rows, percentile_ci
    y = np.asarray(y).astype(int)
    auc = auroc_rows(y, P)
    b_auc = np.array([auroc_rows(y[i], P[:, i]).mean() for i in boot_idx])
    pm = P.mean(axis=0)
    ap = float(average_precision_score(y, pm))
    b_ap = np.array([average_precision_score(y[i], pm[i]) for i in boot_idx])
    lo, hi = percentile_ci(b_auc)
    alo, ahi = percentile_ci(b_ap)
    op = sens_at_spec_boot(y, pm, spec, n_boot=n_boot)
    out = {"auroc": float(auc.mean()), "auroc_ci_low": lo, "auroc_ci_high": hi, "auroc_repeat_min": float(auc.min()),
           "auroc_repeat_max": float(auc.max()), "auprc": ap, "auprc_ci_low": alo, "auprc_ci_high": ahi,
           "auprc_prevalence_baseline": float(y.mean()), "specificity_target": spec, **op,
           "n": int(len(y)), "n_cases": int(y.sum()), "n_controls": int((1 - y).sum())}
    if spec < METRIC_LINK_SPEC:
        op90 = sens_at_spec_boot(y, pm, METRIC_LINK_SPEC, n_boot=n_boot)
        out.update({f"ml90_{k}": v for k, v in op90.items()})
    return out, b_auc


def decide(s: dict) -> str:
    p = s.get("perm_p", np.nan)
    return "supported" if (s["auroc_ci_low"] > 0.5 and np.isfinite(p) and p <= 0.05) else "null"


# ------------------------------------------------------------------------------------------------------------------
# evaluation
# ------------------------------------------------------------------------------------------------------------------

def _design(mats: dict, blocks: list[str], idx: pd.Index) -> pd.DataFrame:
    parts = [mats[b].reindex(idx) for b in blocks if b in mats]
    if not parts:
        return pd.DataFrame(index=idx)
    X = pd.concat(parts, axis=1)
    X = X.loc[:, X.notna().any() & (X.std(ddof=0).fillna(0) > 0)]
    return X


def _has_data(mats: dict, blocks: list[str], idx: pd.Index) -> pd.Series:
    has = pd.Series(False, index=idx)
    for b in blocks:
        if b in mats:
            has |= mats[b].reindex(idx).notna().any(axis=1)
    return has


def measurement_only(blocks: list[str], block_info: dict, measurement_class: str) -> tuple[bool, str | None]:
    classes = {c["id"]: c for c in load_config("measurements")["measurement_classes"]}
    fam = classes.get(measurement_class, {}).get("modality_family")
    kinds = {block_info.get(b, {}).get("kind") for b in blocks}
    dev_classes = {block_info[b].get("measurement_class") for b in blocks if block_info.get(b, {}).get("kind") == "device"}
    if kinds == {"device"} and fam in DEVICE_FAMILIES and dev_classes <= {measurement_class}:
        return True, None
    if kinds == {"ehr_labs"} and fam == "lab":
        return True, None
    if kinds == {"omics"} and fam == "omics":
        return True, None
    return False, (f"the primary model combines blocks {sorted(k for k in kinds if k)} that are not all measurements "
                   f"of class {measurement_class!r}: not a measurement-only metric, kept but not scored")


def evaluate(root: str | Path, *, amend: bool = False, jobs: int | None = None, echo=print) -> dict:
    root = Path(root).resolve()
    inp = K.discover(root)
    manifest = K.load_manifest(inp.manifest_path)
    ds = K.engine_id(manifest["dataset_id"])
    meta = PER.dataset_row(ds)
    if meta is None:
        raise RuntimeError(f"{ds} is not ingested: run `measure-it byod ingest {root}` first")
    if K.sha256_file(inp.manifest_path) != meta["manifest_sha256"]:
        raise RuntimeError(f"manifest.yaml changed since {ds} was ingested: re-run `measure-it byod ingest` first")
    # ---- 1. lock before anything outcome-related is computed
    lock = lock_plan(manifest, meta, amend=amend)
    plan = lock["plan"]
    echo(f"byod evaluate {ds}: plan {lock['status']} (sha256 {lock['sha256'][:12]}, locked {lock['locked_at']})")
    prior = [e for e in ledger_events(ds) if e.get("event") == "evaluated"]
    jobs = jobs or max(1, min(8, (os.cpu_count() or 2) - 1))
    # ---- 2. data
    P = read_table(K.partition_name("participants", ds)).set_index("participant_id")
    y_all = pd.to_numeric(P[plan["label"]], errors="coerce")
    idx_lab = y_all.index[y_all.notna()]
    mats = PER.block_matrices(ds, plan["condition_id"])
    block_info = json.loads(meta["blocks_json"])
    prim = plan["feature_blocks"]
    has_prim = _has_data(mats, prim, idx_lab)
    pop = idx_lab[has_prim.to_numpy()]
    y = y_all.reindex(pop).astype(int).to_numpy()
    n_boot, spec, C = plan["n_bootstrap"], plan["specificity"], float(K.analysis_defaults(
        manifest["primary_analysis"])["C"])
    ks, reps = plan["cv"]["n_splits"], plan["cv"]["n_repeats"]
    from ..labs.lab_cc_engine import delta as paired_delta
    from ..labs.lab_cc_engine import perm_p, stratified_boot_idx
    boot = stratified_boot_idx(y, n_boot, SEED)
    models, store = [], {}

    def run(name: str, blocks: list[str], idx: pd.Index, yy: np.ndarray, bidx, n_perm: int, role: str,
            extra: pd.DataFrame | None = None):
        X = _design(mats, blocks, idx)
        if extra is not None:
            X = pd.concat([X, extra.reindex(idx)], axis=1)
        if X.shape[1] == 0 or len(yy) == 0 or yy.sum() < ks or (1 - yy).sum() < ks:
            models.append({"model": name, "role": role, "blocks": ";".join(blocks), "status": UNKNOWN,
                           "reason": "no usable features or too few cases / controls"})
            return None
        Pm = cv_oof(X.to_numpy(float), yy, C, ks, reps)
        s, b = summarize(yy, Pm, bidx, spec, n_boot)
        null = permutation_null(X.to_numpy(float), yy, C, ks, n_perm, jobs)
        s.update({"perm_p": perm_p(null, s["auroc"]) if len(null) else np.nan, "n_perm": int(len(null)),
                  "null_q95": float(np.quantile(null, 0.95)) if len(null) else np.nan})
        s["decision"] = decide(s)
        models.append({"model": name, "role": role,
                       "blocks": ";".join(list(blocks) + (["covariates"] if extra is not None else [])),
                       "n_features": int(X.shape[1]), "status": "ok", **s})
        store[name] = {"P": Pm, "boots": b, "y": yy, "X": X}
        echo(f"  {name:28s} AUROC {s['auroc']:.3f} ({s['auroc_ci_low']:.3f}-{s['auroc_ci_high']:.3f}), "
             f"sens {s['op_sensitivity']:.2f} at spec {s['op_specificity']:.2f}, perm p "
             f"{s['perm_p'] if np.isfinite(s['perm_p']) else float('nan'):.4f} [{s['decision']}]")
        return s

    echo(f"  population: {int(y.sum())} cases, {int((1 - y).sum())} controls with data in {prim}")
    run("primary", prim, pop, y, boot, plan["n_permutations_primary"], "primary (pre-specified)")
    n2 = plan["n_permutations_secondary"]
    if plan["combined_blocks"]:
        has_c = _has_data(mats, plan["combined_blocks"], idx_lab)
        pc = idx_lab[has_c.to_numpy()]
        yc = y_all.reindex(pc).astype(int).to_numpy()
        run("combined", plan["combined_blocks"], pc, yc, stratified_boot_idx(yc, n_boot, SEED), n2,
            "combined (pre-specified)")
    for b in sorted(mats):
        if [b] == list(prim):
            continue
        hb = _has_data(mats, [b], idx_lab)
        pb = idx_lab[hb.to_numpy()]
        yb = y_all.reindex(pb).astype(int).to_numpy()
        run(f"block:{b}", [b], pb, yb, stratified_boot_idx(yb, n_boot, SEED), n2, "secondary (single block)")
    cov = PER.covariate_matrix(ds, plan["covariates"])
    deltas = []
    if cov is not None and cov.shape[1]:
        run("covariates", [], pop, y, boot, n2, "baseline (age band + sex)", extra=cov)
        run("primary+covariates", prim, pop, y, boot, n2, "secondary (increment)", extra=cov)
        if "covariates" in store and "primary+covariates" in store:
            d = paired_delta(store["primary+covariates"], store["covariates"], y)
            deltas.append({"contrast": "primary+covariates minus covariates", **d})
            echo(f"  Delta AUROC (primary+covariates - covariates) {d['delta_auroc']:+.3f} "
                 f"({d['ci_low']:+.3f} to {d['ci_high']:+.3f})")
    mdf = pd.DataFrame(models)
    prim_row = mdf[mdf["model"] == "primary"].iloc[0].to_dict()
    if prim_row.get("status") != "ok":
        raise RuntimeError(f"primary model could not be fitted: {prim_row.get('reason')}")
    # ---- per-feature effects (descriptive; BH within block)
    feats = feature_effects(mats, block_info, y_all, idx_lab)
    # ---- 3. outputs
    now = utc_now_iso()
    rec = performance_record(manifest, meta, plan, lock, prim_row, prior)
    sig = signatures(ds, manifest, meta, plan, mdf, feats, lock, now)
    write_outputs(ds, manifest, meta, plan, lock, mdf, pd.DataFrame(deltas), feats, rec, sig, now)
    _ledger({"event": "evaluated", "dataset_id": ds, "plan_sha256": lock["sha256"], "auroc": prim_row["auroc"],
             "auroc_ci_low": prim_row["auroc_ci_low"], "auroc_ci_high": prim_row["auroc_ci_high"]})
    echo(f"  record {rec['record_id']}: tier user_supplied_own_computation, class {rec['primary_class']}, "
         f"measurement-only {rec['measurement_only']}" + (" [demo: excluded from default rankings]"
                                                          if rec["demo"] else ""))
    return {"dataset_id": ds, "plan": lock, "models": mdf, "deltas": deltas, "record": rec,
            "n_prior_evaluations": len(prior)}


def feature_effects(mats: dict, block_info: dict, y_all: pd.Series, idx_lab: pd.Index) -> pd.DataFrame:
    from ..labs.lab_cc_engine import per_analyte_effects
    rows = []
    for b, M in mats.items():
        if b in K.BINARY_BLOCK_KINDS:
            continue
        d = M.reindex(idx_lab).copy()
        d["_y"] = y_all.reindex(idx_lab).astype(int).to_numpy()
        cols = [c for c in M.columns if d[c].notna().sum() >= 6]
        if not cols:
            continue
        e = per_analyte_effects(d, "_y", cols, n_boot=1000)
        e.insert(0, "block", b)
        rows.append(e)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _caveats(manifest: dict, meta: dict, lock: dict, prior: list, prim: dict) -> list[str]:
    ct = manifest["comparator"]["type"]
    lab = json.loads(meta["labels_json"])[lock["plan"]["label"]]
    out = []
    if manifest.get("synthetic"):
        out.append("SYNTHETIC tutorial data (generated, not real people): demonstrates the mechanics only.")
    out.append("User-supplied, local-only data: not public and not re-analysable by others; de-identification, labels "
               "and permission attested by the data owner.")
    if ct == "healthy":
        out.append("Comparator: healthy controls. Performance against healthy controls is optimistic: it overstates "
                   "real-world performance, where the question is this illness versus other causes of the same "
                   "symptoms.")
    elif ct == "general_population":
        out.append("Comparator: general population (not clinical look-alikes); prevalence and case mix differ from a "
                   "clinic.")
    else:
        out.append("Comparator: clinical look-alikes (the clinically relevant contrast), as defined by the owner: "
                   f"{manifest['comparator']['description']}.")
    if lab["label_basis"] in ("self_report", "proxy"):
        out.append(f"Label basis {lab['label_basis']}: not a clinical case definition.")
    if lock.get("amends"):
        out.append(f"The analysis plan was amended after an earlier lock ({lock['amends'][:12]}): read as post hoc.")
    if prior:
        out.append(f"This dataset id was evaluated {len(prior)} time(s) before (results/byod/LEDGER.jsonl).")
    if prim.get("decision") == "null":
        out.append("NULL by the locked decision rule (CI includes 0.5 or permutation p > 0.05).")
    out.append("Operating point chosen on the same controls (in-sample threshold on out-of-fold probabilities); the "
               "bootstrap re-chooses it per resample, which only partly corrects the optimism.")
    return out


def performance_record(manifest: dict, meta: dict, plan: dict, lock: dict, prim: dict, prior: list) -> dict:
    ds = meta["dataset_id"]
    block_info = json.loads(meta["blocks_json"])
    mo, why = measurement_only(plan["feature_blocks"], block_info, plan["measurement_class"])
    lab = json.loads(meta["labels_json"])[plan["label"]]
    spec = plan["specificity"]
    use90 = spec < METRIC_LINK_SPEC
    g = (lambda k: prim[f"ml90_{k}"]) if use90 else (lambda k: prim[k])
    from ..scoring.opportunity import condition_sets
    sets = condition_sets()
    members = sorted(sets[lab["condition_id"]]["members"]) if lab["condition_id"] in sets else [lab["condition_id"]]
    return {
        "record_id": f"BYOD-{K.manifest_id_of(ds).upper()}-PRIMARY", "source_kind": K.SOURCE_KIND,
        "quality_tier": 2.5, "quality_tier_label": K.SOURCE_KIND, "target_condition": lab["condition_id"],
        "target_members": json.dumps(members), "primary_class": plan["measurement_class"], "dataset_id": ds,
        "analysis": (f"{'SYNTHETIC ' if manifest.get('synthetic') else ''}user-supplied '{manifest['title']}': blocks "
                     f"{', '.join(plan['feature_blocks'])}; L2 logistic, {plan['cv']['n_splits']}-fold x "
                     f"{plan['cv']['n_repeats']} repeated CV; locked plan {lock['sha256'][:12]}"),
        "label_basis": f"{lab['label_basis']}: {lab['definition']}",
        "comparator": manifest["comparator"]["description"], "comparator_kind": meta["comparator_kind"],
        "n_cases": float(prim["n_cases"]), "n_controls": float(prim["n_controls"]), "auroc": float(prim["auroc"]),
        "auroc_ci_low": float(prim["auroc_ci_low"]), "auroc_ci_high": float(prim["auroc_ci_high"]),
        "auroc_basis": "engine CV AUROC (mean over repeats) with participant-bootstrap 95% CI, under the locked plan",
        "op_sensitivity": float(g("op_sensitivity")), "op_sensitivity_ci_low": float(g("op_sensitivity_ci_low")),
        "op_sensitivity_ci_high": float(g("op_sensitivity_ci_high")), "op_specificity": float(g("op_specificity")),
        "op_specificity_ci_low": float(g("op_specificity_ci_low")),
        "op_specificity_ci_high": float(g("op_specificity_ci_high")), "op_threshold": float(g("op_threshold")),
        "op_basis": (f"threshold at >= {METRIC_LINK_SPEC if use90 else spec:.2f} specificity on the controls' "
                     "repeat-averaged out-of-fold probabilities, stratified bootstrap CI with the threshold re-chosen "
                     "per resample" + (f" (the manifest's {spec:.2f} is below the metric link's 0.90; reported in "
                                       "EVALUATION.md)" if use90 else "")),
        "measurement_only": bool(mo), "exclusion_reason": why,
        "source_object_ids": json.dumps([f"signature:{ds}|{plan['label']}|model_auroc:primary"]),
        "source_tables": json.dumps([f"data/processed/phenotype_signatures__{ds}",
                                     f"results/byod/{K.manifest_id_of(ds)}/EVALUATION.md"]),
        "caveats": json.dumps(_caveats(manifest, meta, lock, prior, prim)),
        "computed_here": "everything: engine run on local user-supplied data under the locked plan",
        "auprc": float(prim["auprc"]), "auprc_ci_low": float(prim["auprc_ci_low"]),
        "auprc_ci_high": float(prim["auprc_ci_high"]), "perm_p": float(prim["perm_p"]),
        "n_perm": int(prim["n_perm"]), "decision": prim["decision"],
        "user_supplied": True, "demo": bool(manifest.get("demo")), "synthetic": bool(manifest.get("synthetic")),
        "plan_sha256": lock["sha256"], "plan_amended": bool(lock.get("amends")), "n_prior_evaluations": len(prior),
    }


def signatures(ds: str, manifest: dict, meta: dict, plan: dict, mdf: pd.DataFrame, feats: pd.DataFrame, lock: dict,
               now: str) -> pd.DataFrame:
    lab = json.loads(meta["labels_json"])[plan["label"]]
    cav = ("SYNTHETIC tutorial data. " if manifest.get("synthetic") else "") + (
        "Healthy controls only: optimistic. " if manifest["comparator"]["type"] == "healthy" else "") + \
        "User-supplied, local-only data (not public)."
    base = {"dataset_id": ds, "phenotype_id": plan["label"],
            "phenotype_label": f"{lab['definition']} vs {manifest['comparator']['description']}",
            "condition_id": lab["condition_id"], "is_proxy": lab["label_basis"] in ("proxy", "self_report"),
            "label_basis": lab["label_basis"], "phenotype_definition": lab["definition"], "caveats": cav,
            "measurement_class_id": plan["measurement_class"], "null_value": 0.5}
    rows = []
    for r in mdf.to_dict("records"):
        if r.get("status") != "ok":
            continue
        rows.append({**base, "object_id": f"signature:{ds}|{plan['label']}|model_auroc:{r['model']}",
                     "feature": f"model_auroc:{r['model']}", "feature_label": f"{r['model']} model ({r['blocks']})",
                     "effect_measure": "auroc", "effect_size": r["auroc"], "effect_unit": "cross-validated AUROC",
                     "ci_low": r["auroc_ci_low"], "ci_high": r["auroc_ci_high"], "p_value": r.get("perm_p"),
                     "q_value": np.nan, "n_cases": r["n_cases"], "n_controls": r["n_controls"],
                     "method": (f"{r['role']}; L2 logistic, {plan['cv']['n_splits']}-fold x {plan['cv']['n_repeats']} "
                                f"CV; participant bootstrap CI ({plan['n_bootstrap']}); permutation p ({r['n_perm']}); "
                                f"plan {lock['sha256'][:12]}"),
                     "evidence_level": f"{r['decision']}_single_dataset"})
    for r in feats.to_dict("records") if len(feats) else []:
        if not np.isfinite(r.get("hedges_g", np.nan)):
            continue
        rows.append({**base, "object_id": f"signature:{ds}|{plan['label']}|{r['analyte']}",
                     "feature": r["analyte"], "feature_label": f"{r['analyte']} ({r['block']})",
                     "effect_measure": "hedges_g", "effect_size": r["hedges_g"],
                     "effect_unit": "Hedges' g (cases minus controls)", "ci_low": r["g_ci_low"], "ci_high": r["g_ci_high"],
                     "p_value": r["p_value"], "q_value": r["q_value"], "n_cases": r["n_cases"],
                     "n_controls": r["n_controls"], "null_value": 0.0,
                     "method": "per-feature descriptive contrast: stratified bootstrap CI (1,000), Mann-Whitney p, BH q "
                               "within block (secondary, not pre-specified as a claim)",
                     "evidence_level": "descriptive_single_dataset"})
    sig = pd.DataFrame(rows)
    return sig


def write_outputs(ds, manifest, meta, plan, lock, mdf, deltas, feats, rec, sig, now) -> None:
    from .ingest import union_family, upsert_dataset_row
    src = f"user-supplied dataset '{manifest['title']}' ({manifest['owner']}) [{ds}]"
    ver = f"{EVAL_VERSION}; plan sha256 {lock['sha256'][:16]}"
    notes = (f"{K.PROVENANCE_NOTE}; local-only; " + ("SYNTHETIC tutorial data; " if manifest.get("synthetic") else "")
             + "cohort-level results computed from person-level rows under a locked plan")
    s = add_provenance(sig, data_layer="person", source_name=src, source_version=ver, retrieved_at=now,
                       evidence_type="person_derived_feature", source_record_id="object_id",
                       evidence_level=sig["evidence_level"].astype(str).to_numpy(), provenance_notes=notes)
    write_table(s, K.partition_name("phenotype_signatures", ds), producer=PRODUCER,
                description=f"signatures of user-supplied dataset {ds} ({K.PROVENANCE_NOTE})")
    union_family("phenotype_signatures")
    # the performance record (aggregate only: no participant id, no person row)
    r = pd.DataFrame([rec])
    old = read_table(K.RECORDS_TABLE) if table_exists(K.RECORDS_TABLE) else None
    if old is not None:
        old = old[old["dataset_id"] != ds].drop(columns=[c for c in old.columns if c in (
            "data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
            "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes")])
        r = pd.concat([old, r], ignore_index=True)
    write_records_table(r)
    meta = dict(meta)
    for k in ("data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
              "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes"):
        meta.pop(k, None)
    meta.update({"status": "evaluated", "plan_sha256": lock["sha256"], "evaluated_at": now})
    upsert_dataset_row(meta)
    out = K.local_only_dir(K.results_path(ds))
    mdf.to_csv(out / "evaluation_models.csv", index=False)
    if len(deltas):
        deltas.to_csv(out / "evaluation_deltas.csv", index=False)
    if len(feats):
        feats.to_csv(out / "evaluation_features.csv", index=False)
    (out / "performance_record.json").write_text(json.dumps(rec, indent=1, default=str))
    (out / "EVALUATION.md").write_text(evaluation_markdown(ds, manifest, plan, lock, mdf, deltas, feats, rec))


def write_records_table(r: pd.DataFrame) -> None:
    from ..config import PROCESSED
    if r.empty:
        for suf in (".parquet", ".meta.json"):
            (PROCESSED / f"{K.RECORDS_TABLE}{suf}").unlink(missing_ok=True)
        return
    r = add_provenance(r.reset_index(drop=True), data_layer="measurement",
                       source_name="measure_it.byod.evaluate (user-supplied datasets, locked plans)",
                       source_version=EVAL_VERSION, retrieved_at=utc_now_iso(),
                       evidence_type="derived_evidence_summary", source_record_id="record_id",
                       evidence_level=K.SOURCE_KIND,
                       provenance_notes=("aggregate performance records only (no person rows, no participant ids, no "
                                         "geography); read by measure_it.scoring.metric_link; demo records excluded "
                                         "unless MEASURE_IT_BYOD_DEMO=1"))
    write_table(r, K.RECORDS_TABLE, producer=PRODUCER,
                description="Performance records of user-supplied datasets (tier user_supplied_own_computation)")


def _f(v, nd=3) -> str:
    try:
        return f"{float(v):.{nd}f}" if np.isfinite(float(v)) else "n/a"
    except (TypeError, ValueError):
        return "n/a"


def evaluation_markdown(ds, manifest, plan, lock, mdf, deltas, feats, rec) -> str:
    lines = [f"# Evaluation: {manifest['title']} ({ds})", "",
             ("**SYNTHETIC tutorial dataset: generated data, not real people.**\n\n" if manifest.get("synthetic")
              else "") + f"Plan locked {lock['locked_at']} (sha256 `{lock['sha256'][:16]}`; PLAN.md). Label "
             f"`{plan['label']}` -> `{plan['condition_id']}` ({plan['label_basis']}); comparator "
             f"{manifest['comparator']['type']}: {manifest['comparator']['description']}.", "",
             "| model | role | n cases / controls | AUROC (95% CI) | AUPRC (95% CI; prevalence) | sensitivity at "
             "specificity | permutation p | decision |", "|---|---|---|---|---|---|---|---|"]
    for r in mdf.to_dict("records"):
        if r.get("status") != "ok":
            lines.append(f"| {r['model']} | {r['role']} | | {UNKNOWN}: {r.get('reason')} | | | | |")
            continue
        lines.append(f"| {r['model']} | {r['role']} | {r['n_cases']} / {r['n_controls']} | {_f(r['auroc'])} "
                     f"({_f(r['auroc_ci_low'])}-{_f(r['auroc_ci_high'])}) | {_f(r['auprc'])} ({_f(r['auprc_ci_low'])}-"
                     f"{_f(r['auprc_ci_high'])}; {_f(r['auprc_prevalence_baseline'], 2)}) | {_f(r['op_sensitivity'], 2)} "
                     f"({_f(r['op_sensitivity_ci_low'], 2)}-{_f(r['op_sensitivity_ci_high'], 2)}) at "
                     f"{_f(r['op_specificity'], 2)} | {_f(r['perm_p'], 4)} ({r['n_perm']}) | {r['decision']} |")
    if len(deltas):
        for d in deltas.to_dict("records"):
            lines.append(f"\nDelta AUROC, {d['contrast']}: {d['delta_auroc']:+.3f} ({d['ci_low']:+.3f} to "
                         f"{d['ci_high']:+.3f}).")
    lines += ["", "## Performance record (what crosses to the ranking)", "",
              f"`{rec['record_id']}`: condition `{rec['target_condition']}`, class `{rec['primary_class']}`, tier "
              f"user_supplied_own_computation (2.5), AUROC {_f(rec['auroc'])} ({_f(rec['auroc_ci_low'])}-"
              f"{_f(rec['auroc_ci_high'])}), sensitivity {_f(rec['op_sensitivity'], 2)} at specificity "
              f"{_f(rec['op_specificity'], 2)}, n {int(rec['n_cases'])} + {int(rec['n_controls'])}; measurement-only "
              f"{rec['measurement_only']}" + (f" ({rec['exclusion_reason']})" if rec["exclusion_reason"] else "")
              + (". Demo record: excluded from rankings unless `byod deploy --demo`." if rec["demo"] else "."), "",
              "Caveats:", ""] + [f"* {c}" for c in json.loads(rec["caveats"])]
    if len(feats):
        top = feats.dropna(subset=["q_value"]).sort_values("q_value").head(10)
        lines += ["", "## Per-feature contrasts (descriptive, top 10 by q)", "",
                  "| block | feature | Hedges g (95% CI) | q |", "|---|---|---|---|"]
        lines += [f"| {r['block']} | {r['analyte']} | {_f(r['hedges_g'], 2)} ({_f(r['g_ci_low'], 2)} to "
                  f"{_f(r['g_ci_high'], 2)}) | {_f(r['q_value'], 4)} |" for r in top.to_dict("records")]
    return "\n".join(lines) + "\n"
