"""`measure-it byod subgroup <dir>`: a device-defined subgroup inside a condition and its computable phenotype.

The manifest's `subgroup_analysis` section (schema: templates/byod/manifest.schema.json) says, BEFORE any outcome is
computed: inside the cases of label `within_label`, who is device-positive (a participants.csv column of 1/0, or a
device feature against a threshold), which harmonized domains may describe them (`ehr_domains`: condition,
measurement, drug, demographic, survey; never the device or omics features), the maximum number of features and the
CV / bootstrap / permutation settings. Order (enforced):

 1. lock: the subgroup plan is hashed (sha256) and written to results/byod/<id>/SUBGROUP_PLAN.md and
    subgroup_plan_lock.json (byod.evaluate.lock_document; a changed plan is refused unless --amend; LEDGER.jsonl).
 2. features: person_concepts__byod_<id> -> measure_it.harmonize.features.feature_matrix (ICD-10-CM categories,
    LOINC medians, RxNorm presence, SURVEY items, DEMOG age midpoint / female), restricted to the cases with a defined
    device finding. Removed: label-defining ICD-10-CM categories (any code that maps to the label's condition, and the
    condition's own ICD-10-CM codes from the ontology), the plan's `exclude_feature_keys`, binary features with fewer
    than `min_feature_participants` participants on either side, numeric features observed in under half of the
    cases or constant. Label-free filters only.
 3. model: L1 logistic regression (median imputation + standardisation fitted inside each training fold; the largest
    C on the grid C x 2^-k, k = 0..13, with at most `max_features` non-zero coefficients), participant-level repeated
    stratified k-fold CV (byod.evaluate.cv_oof), participant bootstrap CI and label-permutation null
    (byod.evaluate.summarize / permutation_null; lab_cc_engine.stratified_boot_idx / perm_p); the same decision rule
    as `evaluate` (supported = CI lower bound > 0.5 and permutation p <= 0.05).
 4. computable phenotype: refit on all cases -> results/byod/<id>/computable_phenotype.json (contract section 2;
    every feature says whether a clinic EHR can evaluate it), strata fractions (age band x sex, n < 10 suppressed,
    pooled 'all' row with a Wilson 95% CI), aggregate tables `computable_phenotypes` and `subgroup_strata_fractions`
    (no participant id, no geography), results/byod/<id>/SUBGROUP.md.
"""
from __future__ import annotations

import json
import os
import warnings
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin

from ..config import SEED, utc_now_iso
from ..provenance import add_provenance
from ..store import read_table, table_exists, write_table
from . import common as K
from . import evaluate as E
from . import person as PER

PRODUCER = "measure_it.byod.subgroup"
SUBGROUP_VERSION = "byod_subgroup 2026-10-07.1"
SUPPRESS_BELOW = 10
C_GRID_STEPS = 14
THRESHOLD_RULE = "Youden J on the repeat-averaged out-of-fold probabilities (in-sample choice: optimistic)"
_PROV = ("data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
         "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes")


# ------------------------------------------------------------------------------------------------------------------
# model
# ------------------------------------------------------------------------------------------------------------------

class L1TopK(ClassifierMixin, BaseEstimator):
    """Median imputation + standardisation + L1 logistic regression with at most `max_features` non-zero
    coefficients (the largest C on the grid C * 2^-k giving that many). Deterministic."""

    def __init__(self, C: float = 1.0, max_features: int = 10, n_grid: int = C_GRID_STEPS):
        self.C = C
        self.max_features = max_features
        self.n_grid = n_grid

    def fit(self, X, y):
        from sklearn.impute import SimpleImputer
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        X = np.asarray(X, float)
        self.imputer_ = SimpleImputer(strategy="median", keep_empty_features=True).fit(X)
        Z0 = self.imputer_.transform(X)
        self.scaler_ = StandardScaler().fit(Z0)
        Z = self.scaler_.transform(Z0)
        m = None
        for k in range(int(self.n_grid)):
            c = float(self.C) * 0.5 ** k
            m = LogisticRegression(penalty="l1", solver="liblinear", C=c, max_iter=10000, random_state=SEED)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m.fit(Z, y)
            self.C_used_ = c
            if np.count_nonzero(m.coef_) <= int(self.max_features):
                break
        self.model_ = m
        self.classes_ = m.classes_
        return self

    def predict_proba(self, X):
        return self.model_.predict_proba(self.scaler_.transform(self.imputer_.transform(np.asarray(X, float))))

    def decision_function(self, X):
        return self.model_.decision_function(self.scaler_.transform(self.imputer_.transform(np.asarray(X, float))))


def make_l1_model(C: float, max_features: int = 10) -> L1TopK:
    return L1TopK(C=C, max_features=max_features)


# ------------------------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------------------------

def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n <= 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (float(max(0.0, c - h)), float(min(1.0, c + h)))


def condition_icd_codes(condition_id: str) -> list[str]:
    """ICD-10-CM codes the ontology maps to a condition (and its set members); [] when none (e.g. ptlds)."""
    if not table_exists("condition_ontology_mappings"):
        return []
    m = read_table("condition_ontology_mappings", columns=["canonical_condition_id", "target_ontology", "target_id"])
    members = PER._members(condition_id)
    sel = m[(m["canonical_condition_id"].isin(members)) & (m["target_ontology"] == "ICD10CM")]
    return sorted({str(t).split(":", 1)[1] for t in sel["target_id"].dropna()})


def plan_document(manifest: dict, meta: dict) -> dict:
    sa = K.subgroup_defaults(manifest["subgroup_analysis"])
    labels = json.loads(meta["labels_json"])
    lab = labels[sa["within_label"]]
    return {
        "command": "subgroup", "dataset_id": meta["dataset_id"], "engine_version": SUBGROUP_VERSION, "seed": SEED,
        "name": sa["name"], "description": sa["description"], "within_label": sa["within_label"],
        "condition_id": lab["condition_id"], "label_definition": lab["definition"], "label_basis": lab["label_basis"],
        "base_population_icd10cm": condition_icd_codes(lab["condition_id"]),
        "device_positive": sa["device_positive"],
        "ehr_domains": sorted(sa["ehr_domains"]),
        "feature_source": "person_concepts (measure_it.harmonize.features.feature_matrix, condition_level=category, "
                          "absent=auto); DEVICE and OMICS vocabularies never enter the model",
        "feature_filters": {"min_feature_participants": int(sa["min_feature_participants"]),
                            "numeric_min_observed_fraction": 0.5,
                            "exclude_feature_keys": sorted(sa["exclude_feature_keys"]),
                            "label_defining_exclusion": "ICD-10-CM categories containing a code that maps to the "
                                                        "label condition (ontology normaliser) or one of the "
                                                        "condition's own ICD-10-CM codes"},
        "model": f"L1 logistic (liblinear), median imputation + standardisation inside each training fold; largest C "
                 f"on the grid {float(sa['C'])} x 2^-k (k = 0..{C_GRID_STEPS - 1}) with <= {int(sa['max_features'])} "
                 "non-zero coefficients",
        "max_features": int(sa["max_features"]), "C": float(sa["C"]),
        "cv": {"scheme": "participant-level stratified k-fold, shuffled, seed = SEED + repeat",
               "n_splits": int(sa["n_splits"]), "n_repeats": int(sa["n_repeats"])},
        "n_bootstrap": int(sa["n_bootstrap"]), "n_permutations": int(sa["n_permutations"]),
        "decision_rule": "supported = AUROC 95% CI lower bound > 0.5 AND permutation p <= 0.05; otherwise null",
        "decision_threshold_rule": THRESHOLD_RULE,
        "strata": f"age band (10 years) x sex among the cases with a defined device finding; cells with fewer than "
                  f"{SUPPRESS_BELOW} cases suppressed; pooled row age_band=all, sex=all with a Wilson 95% CI",
    }


def device_positive(ds: str, dp: dict, idx: pd.Index) -> tuple[pd.Series, str]:
    """1 / 0 / NaN per participant in idx, and a human description."""
    if dp.get("column"):
        P = read_table(K.partition_name("participants", ds)).set_index("participant_id")
        if dp["column"] not in P.columns:
            raise RuntimeError(f"device_positive column {dp['column']!r} is not in participants__{ds}: re-ingest")
        y = pd.to_numeric(P[dp["column"]], errors="coerce").reindex(idx).astype(float)
        return y, f"device-positive = participants.csv column `{dp['column']}` == 1 (owner-supplied device finding)"
    mats = PER.device_matrices(ds)
    block = dp["block"]
    if block not in mats:
        raise RuntimeError(f"device block {block!r} has no features in {ds}")
    M = mats[block]
    dev = block[len("device_"):]
    feat = str(dp["feature"])
    cands = [f"{dev}__{feat}", feat, f"{dev}__{feat}__mean"]
    col = next((c for c in cands if c in M.columns), None)
    if col is None:
        raise RuntimeError(f"device feature {feat!r} not among {block} features {list(M.columns)[:12]}")
    x = M[col].reindex(idx).astype(float)
    thr, op = float(dp["threshold"]), dp.get("direction", ">=")
    cmp = {">=": x >= thr, ">": x > thr, "<=": x <= thr, "<": x < thr}[op]
    y = cmp.astype(float).where(x.notna())
    return y, f"device-positive = {block}.{col.split('__', 1)[-1]} {op} {thr:g} (locked plan)"


def label_defining_categories(concepts: pd.DataFrame, condition_id: str, base_codes: list[str]) -> set[str]:
    from ..harmonize.features import icd10_category
    members = PER._members(condition_id)
    icd = concepts.loc[concepts["vocabulary"] == "ICD10CM", "concept_code"].astype(str).unique()
    out = {icd10_category(c) for c in icd if PER.icd_condition(c) in members}
    out |= {icd10_category(c) for c in base_codes}
    return {c for c in out if c}


def select_features(fm: pd.DataFrame, plan: dict, defining: set[str]) -> tuple[pd.DataFrame, list[dict]]:
    """Label-free feature filters; returns the design matrix and one audit row per candidate feature."""
    rules = plan["feature_filters"]
    mn = int(rules["min_feature_participants"])
    excl = set(rules["exclude_feature_keys"])
    keep, audit = [], []
    for c in fm.columns:
        x = fm[c]
        vocab = c.split(":", 1)[0]
        why = None
        if vocab in ("DEVICE", "OMICS"):
            why = "device / omics features never enter the subgroup model"
        elif c in excl:
            why = "excluded by the locked plan"
        elif vocab == "ICD10CM" and c.split(":", 1)[1] in defining:
            why = "label-defining ICD-10-CM category"
        else:
            obs = x.dropna()
            binary = len(obs) > 0 and set(np.unique(obs)) <= {0.0, 1.0}
            if binary:
                if (obs == 1).sum() < mn or (obs == 0).sum() < mn:
                    why = f"fewer than {mn} participants on one side"
            elif len(obs) < max(mn, 0.5 * len(x)):
                why = "observed in fewer than half of the cases"
            elif obs.std(ddof=0) == 0:
                why = "constant"
        audit.append({"feature_key": c, "included": why is None, "reason": why or "",
                      "n_observed": int(x.notna().sum())})
        if why is None:
            keep.append(c)
    return fm[keep], audit


def phenotype_features(model: L1TopK, X: pd.DataFrame, catalog: pd.DataFrame) -> tuple[list[dict], float]:
    """Refit model -> contract feature list (non-zero coefficients) and the intercept on the transformed scale."""
    from ..harmonize.self_report import crosswalk_for
    beta = model.model_.coef_[0]
    intercept = float(model.model_.intercept_[0])
    mu, sd = model.scaler_.mean_, model.scaler_.scale_
    cat = catalog.set_index("feature_key")
    feats = []
    for j, key in enumerate(X.columns):
        b = float(beta[j])
        if b == 0.0:
            continue
        info = {k: (None if not isinstance(v, (list, tuple, np.ndarray)) and pd.isna(v) else v)
                for k, v in (cat.loc[key].to_dict() if key in cat.index else {}).items()}
        vocab = key.split(":", 1)[0]
        obs = X[key].dropna()
        binary = set(np.unique(obs)) <= {0.0, 1.0}
        f = {"feature_key": key, "domain": info.get("domain", "demographic" if vocab == "DEMOG" else ""),
             "vocabulary": vocab, "codes": list(info.get("codes") or []),
             "code_match": ("category_prefix" if vocab == "ICD10CM" else "derived" if vocab == "DEMOG" else "exact"),
             "label": str(info.get("label", key)), "center": float(mu[j]), "scale": float(sd[j]),
             "coefficient_standardized": b, "unit": info.get("unit"),
             "mapping_methods": list(info.get("mapping_methods") or []),
             "mapping_confidences": list(info.get("mapping_confidences") or []),
             "platform": info.get("platform")}
        if binary and vocab in ("ICD10CM", "RXNORM", "DEMOG"):
            f["transform"] = "presence"
            f["coefficient"] = b / float(sd[j])
            intercept -= b * float(mu[j]) / float(sd[j])
        else:
            f["transform"] = "standardized_value"
            f["coefficient"] = b
        f["ehr_evaluable"], f["ehr_note"] = ehr_evaluability(f)
        if vocab == "SURVEY":
            f["reference_crosswalk"] = crosswalk_for(key.split(":", 1)[1])
        feats.append(f)
    feats.sort(key=lambda f: -abs(f["coefficient_standardized"]))
    return feats, intercept


def ehr_evaluability(f: dict) -> tuple[bool, str]:
    vocab = f["vocabulary"]
    methods = set(f.get("mapping_methods") or [])
    if vocab == "ICD10CM":
        note = "any code in the ICD-10-CM category (problem list / encounter diagnoses)"
        if "self_report_map" in methods:
            note += "; in this dataset the feature came from self-reported history, so an EHR code is a stricter proxy"
        return True, note
    if vocab == "LOINC":
        note = f"LOINC-coded result (median; unit {f.get('unit')})"
        if f.get("platform") and "nmr_nightingale" in str(f.get("platform")):
            note += ("; derived here from Nightingale NMR: routine clinical-chemistry values are not interchangeable "
                     "without calibration")
        return True, note
    if vocab == "RXNORM":
        return True, "medication list mapped to the RxNorm ingredient"
    if vocab == "DEMOG":
        return True, "date of birth / administrative sex in the EHR (age as the band midpoint)"
    if vocab == "SURVEY":
        return False, "patient-reported instrument not stored in routine EHR data; administer the instrument"
    return False, "not observable outside the study"


def strata_fractions(y: pd.Series, demo: pd.DataFrame) -> list[dict]:
    d = pd.DataFrame({"y": y}).join(demo, how="left")
    d["age_band"] = d["age_band"].fillna("unknown")
    d["sex"] = d["sex"].fillna("unknown")
    rows = []
    n, k = int(len(d)), int(d["y"].sum())
    lo, hi = wilson(k, n)
    rows.append({"age_band": "all", "sex": "all", "n_cases": n, "n_positive": k, "fraction": k / n if n else None,
                 "ci_low": lo, "ci_high": hi, "suppressed": False})
    for (a, s), g in d.groupby(["age_band", "sex"], sort=True):
        n, k = int(len(g)), int(g["y"].sum())
        if n < SUPPRESS_BELOW:
            rows.append({"age_band": str(a), "sex": str(s), "n_cases": None, "n_positive": None, "fraction": None,
                         "ci_low": None, "ci_high": None, "suppressed": True})
            continue
        lo, hi = wilson(k, n)
        rows.append({"age_band": str(a), "sex": str(s), "n_cases": n, "n_positive": k, "fraction": k / n,
                     "ci_low": lo, "ci_high": hi, "suppressed": False})
    return rows


def youden_threshold(y: np.ndarray, p: np.ndarray) -> tuple[float, float, float]:
    from sklearn.metrics import roc_curve
    fpr, tpr, thr = roc_curve(y, p)
    j = int(np.argmax(tpr - fpr))
    t = float(min(1.0, max(0.0, thr[j])))
    return t, float(tpr[j]), float(1 - fpr[j])


# ------------------------------------------------------------------------------------------------------------------
# the command
# ------------------------------------------------------------------------------------------------------------------

def subgroup(root: str | Path, *, amend: bool = False, jobs: int | None = None, echo=print) -> dict:
    from ..harmonize import schema as HS
    from ..harmonize.features import feature_catalog, feature_matrix
    from ..labs.lab_cc_engine import perm_p, stratified_boot_idx
    root = Path(root).resolve()
    inp = K.discover(root)
    manifest = K.load_manifest(inp.manifest_path)
    ds = K.engine_id(manifest["dataset_id"])
    if not manifest.get("subgroup_analysis"):
        raise RuntimeError("manifest.yaml has no subgroup_analysis section (docs/BRING_YOUR_OWN_DATA.md)")
    meta = PER.dataset_row(ds)
    if meta is None:
        raise RuntimeError(f"{ds} is not ingested: run `measure-it byod ingest {root}` first")
    if K.sha256_file(inp.manifest_path) != meta["manifest_sha256"]:
        raise RuntimeError(f"manifest.yaml changed since {ds} was ingested: re-run `measure-it byod ingest` first")
    if not table_exists(f"{HS.TABLE}__{ds}"):
        from ..harmonize.byod import build as hbuild
        hbuild(ds, manifest=manifest, echo=echo)
    # ---- 1. lock before anything outcome-related is computed
    plan = plan_document(manifest, meta)
    lock = E.lock_document(ds, plan, amend=amend, lock_name="subgroup_plan_lock.json", md_name="SUBGROUP_PLAN.md",
                           event="subgroup_plan", what="subgroup plan")
    echo(f"byod subgroup {ds}: plan {lock['status']} (sha256 {lock['sha256'][:12]}, locked {lock['locked_at']})")
    jobs = jobs or max(1, min(8, (os.cpu_count() or 2) - 1))
    # ---- 2. population + features
    P = read_table(K.partition_name("participants", ds)).set_index("participant_id")
    lab = pd.to_numeric(P[plan["within_label"]], errors="coerce")
    cases = lab.index[lab == 1]
    y_all, sub_def = device_positive(ds, plan["device_positive"], cases)
    pop = y_all.index[y_all.notna()]
    y = y_all.reindex(pop).astype(int).to_numpy()
    n1, n0 = int(y.sum()), int(len(y) - y.sum())
    echo(f"  {len(cases)} cases of {plan['within_label']}; {len(pop)} with a device finding: {n1} positive, "
         f"{n0} negative")
    ks = plan["cv"]["n_splits"]
    if n1 < ks or n0 < ks or min(n1, n0) < 5:
        raise RuntimeError(f"too few device-positive ({n1}) or device-negative ({n0}) cases for {ks}-fold CV")
    concepts = HS.read_concepts(ds)
    fm = feature_matrix(None, concepts=concepts, domains=plan["ehr_domains"], condition_level="category")
    fm = fm.reindex(pop)
    defining = label_defining_categories(concepts, plan["condition_id"], plan["base_population_icd10cm"])
    X, feat_audit = select_features(fm, plan, defining)
    if X.shape[1] == 0:
        raise RuntimeError("no usable harmonized feature in the plan's domains")
    echo(f"  {X.shape[1]} candidate features ({fm.shape[1]} before filters; label-defining categories excluded: "
         f"{sorted(defining) or 'none'})")
    # ---- 3. CV, bootstrap, permutation null
    factory = partial(make_l1_model, max_features=plan["max_features"])
    Xn = X.to_numpy(float)
    C, reps, n_boot = plan["C"], plan["cv"]["n_repeats"], plan["n_bootstrap"]
    Pm = E.cv_oof(Xn, y, C, ks, reps, model_factory=factory)
    boot = stratified_boot_idx(y, n_boot, SEED)
    s, _ = E.summarize(y, Pm, boot, 0.90, n_boot)
    null = E.permutation_null(Xn, y, C, ks, plan["n_permutations"], jobs, model_factory=factory)
    s.update({"perm_p": perm_p(null, s["auroc"]) if len(null) else np.nan, "n_perm": int(len(null)),
              "null_q95": float(np.quantile(null, 0.95)) if len(null) else np.nan})
    s["decision"] = E.decide(s)
    echo(f"  L1 model AUROC {s['auroc']:.3f} ({s['auroc_ci_low']:.3f}-{s['auroc_ci_high']:.3f}), perm p "
         f"{s['perm_p']:.4f} [{s['decision']}]")
    thr, t_sens, t_spec = youden_threshold(y, Pm.mean(axis=0))
    # ---- 4. refit on all cases -> computable phenotype
    model = make_l1_model(C, plan["max_features"]).fit(Xn, y)
    catalog = feature_catalog(concepts[concepts["participant_id"].isin(pop)])
    feats, intercept = phenotype_features(model, X, catalog)
    demo = concepts[concepts["vocabulary"] == "DEMOG"].pivot_table(
        index="participant_id", columns="concept_code", values="value_as_string", aggfunc="first")
    demo = demo.rename(columns={"age_band": "age_band", "sex": "sex"}).reindex(columns=["age_band", "sex"])
    strata = strata_fractions(pd.Series(y, index=pop), demo)
    prior = [e for e in E.ledger_events(ds) if e.get("event") == "subgroup_plan_amended"]
    labels = json.loads(meta["labels_json"])
    labinfo = labels[plan["within_label"]]
    n_ehr = sum(f["ehr_evaluable"] for f in feats)
    mass = sum(abs(f["coefficient_standardized"]) for f in feats) or 1.0
    ehr_mass = sum(abs(f["coefficient_standardized"]) for f in feats if f["ehr_evaluable"]) / mass
    caveats = subgroup_caveats(manifest, lock, s, feats, n1, n0, plan)
    phen = {
        "phenotype_id": f"{ds}:{plan['name']}", "version": 1 + len(prior), "dataset_id": ds,
        "condition_id": plan["condition_id"],
        "base_population": {"description": f"cases of label `{plan['within_label']}` ({labinfo['definition']}) with a "
                                           "defined device finding", "icd10cm": plan["base_population_icd10cm"],
                            "label_column": plan["within_label"], "label_basis": labinfo["label_basis"],
                            "n_cases": int(len(cases)), "n_cases_with_device_finding_defined": int(len(pop))},
        "subgroup_definition": sub_def, "model": "l1_logistic", "intercept": float(intercept), "features": feats,
        "decision_threshold": thr, "decision_threshold_rule": THRESHOLD_RULE,
        "performance": {"auroc": s["auroc"], "auroc_ci_low": s["auroc_ci_low"], "auroc_ci_high": s["auroc_ci_high"],
                        "n_positive": n1, "n_negative": n0,
                        "comparator": "condition cases without the device finding",
                        "perm_p": float(s["perm_p"]) if np.isfinite(s["perm_p"]) else None, "n_perm": s["n_perm"],
                        "decision": s["decision"], "auprc": s["auprc"],
                        "auprc_prevalence_baseline": s["auprc_prevalence_baseline"],
                        "threshold_sensitivity_oof": t_sens, "threshold_specificity_oof": t_spec,
                        "cv": f"{ks}-fold x {reps} repeats", "C_used_refit": float(model.C_used_),
                        "n_candidate_features": int(X.shape[1])},
        "strata_fractions": strata,
        "ehr_evaluable_summary": {"n_features": len(feats), "n_ehr_evaluable": int(n_ehr),
                                  "ehr_evaluable_coefficient_share": float(ehr_mass),
                                  "note": "share of sum |standardized coefficient| carried by features a clinic EHR "
                                          "can evaluate (ICD-10-CM, LOINC, RxNorm, demographics)"},
        "plan_sha256": lock["sha256"], "caveats": caveats,
        "synthetic": bool(manifest.get("synthetic")), "demo": bool(manifest.get("demo")),
    }
    from ..harmonize.phenotype import validate_phenotype
    errs = validate_phenotype(phen)
    if errs:
        raise RuntimeError(f"computable phenotype fails its schema: {errs[:3]}")
    text = json.dumps(phen, default=float)
    leaked = [p for p in P.index.astype(str) if p in text]
    if leaked:
        raise RuntimeError("computable phenotype would carry participant ids; refusing to write it")
    write_outputs(ds, manifest, plan, lock, phen, feat_audit, s, echo)
    E._ledger({"event": "subgroup_evaluated", "dataset_id": ds, "plan_sha256": lock["sha256"],
               "phenotype_id": phen["phenotype_id"], "auroc": s["auroc"], "auroc_ci_low": s["auroc_ci_low"],
               "auroc_ci_high": s["auroc_ci_high"]})
    return {"dataset_id": ds, "plan": lock, "phenotype": phen, "summary": s, "feature_audit": feat_audit}


def subgroup_caveats(manifest, lock, s, feats, n1, n0, plan) -> list[str]:
    out = []
    if manifest.get("synthetic"):
        out.append("SYNTHETIC tutorial data (generated, not real people): demonstrates the mechanics only.")
    out.append("User-supplied, local-only data; single dataset; no external validation. The phenotype describes who "
               "resembles the device-positive cases in this cohort, not a diagnosis.")
    out.append(f"Small groups ({n1} device-positive, {n0} device-negative cases): wide intervals; coefficients of an "
               "L1 model are unstable across resamples.")
    out.append("Decision threshold chosen in-sample (Youden J on out-of-fold probabilities): optimistic.")
    if any(f["vocabulary"] == "SURVEY" for f in feats):
        out.append("Some features are patient-reported instruments that a clinic EHR does not hold "
                   "(ehr_evaluable=false); scoring elsewhere uses the feature center for them and reports coverage.")
    if any("self_report_map" in (f.get("mapping_methods") or []) for f in feats):
        out.append("ICD-10-CM features derived from self-reported history (mapping confidence <= medium); EHR codes "
                   "are a different, usually stricter, ascertainment.")
    if any("nmr_nightingale" in str(f.get("platform") or "") for f in feats):
        out.append("LOINC features derived from Nightingale NMR (explicit unit conversion); routine laboratory values "
                   "are not interchangeable without calibration.")
    if lock.get("amends"):
        out.append(f"The subgroup plan was amended after an earlier lock ({lock['amends'][:12]}): read as post hoc.")
    if s.get("decision") == "null":
        out.append("NULL by the locked decision rule (CI includes 0.5 or permutation p > 0.05): the harmonized "
                   "features do not separate the device subgroup reliably.")
    return out


def write_outputs(ds, manifest, plan, lock, phen, feat_audit, s, echo) -> None:
    out = K.local_only_dir(K.results_path(ds))
    (out / "computable_phenotype.json").write_text(json.dumps(phen, indent=1, default=float))
    pd.DataFrame(feat_audit).to_csv(out / "subgroup_feature_audit.csv", index=False)
    pd.DataFrame(phen["features"]).to_csv(out / "subgroup_features.csv", index=False)
    (out / "SUBGROUP.md").write_text(subgroup_markdown(ds, manifest, plan, lock, phen, s))
    now = utc_now_iso()
    src = f"measure_it.byod.subgroup (user-supplied dataset {ds}; locked plan {lock['sha256'][:12]})"
    row = {"phenotype_id": phen["phenotype_id"], "dataset_id": ds, "condition_id": phen["condition_id"],
           "subgroup_name": plan["name"], "version": phen["version"], "model": phen["model"],
           "subgroup_definition": phen["subgroup_definition"], "n_features": len(phen["features"]),
           "n_features_ehr_evaluable": phen["ehr_evaluable_summary"]["n_ehr_evaluable"],
           "ehr_evaluable_coefficient_share": phen["ehr_evaluable_summary"]["ehr_evaluable_coefficient_share"],
           "feature_keys": json.dumps([f["feature_key"] for f in phen["features"]]),
           "auroc": s["auroc"], "auroc_ci_low": s["auroc_ci_low"], "auroc_ci_high": s["auroc_ci_high"],
           "perm_p": phen["performance"]["perm_p"], "decision": s["decision"],
           "n_positive": phen["performance"]["n_positive"], "n_negative": phen["performance"]["n_negative"],
           "positive_fraction": phen["strata_fractions"][0]["fraction"],
           "decision_threshold": phen["decision_threshold"], "plan_sha256": lock["sha256"],
           "plan_amended": bool(lock.get("amends")), "synthetic": phen["synthetic"], "demo": phen["demo"],
           "phenotype_json": json.dumps(phen, default=float), "caveats": json.dumps(phen["caveats"])}
    _upsert(K.PHENOTYPES_TABLE, pd.DataFrame([row]), ds, src, now, "phenotype_id",
            "Computable phenotypes of device-defined subgroups (aggregate; docs/HARMONIZATION_CONTRACT.md section 2)")
    st = pd.DataFrame(phen["strata_fractions"])
    st.insert(0, "phenotype_id", phen["phenotype_id"])
    st.insert(1, "condition_id", phen["condition_id"])
    st.insert(2, "dataset_id", ds)
    st["demo"] = bool(phen["demo"])   # consumers (geography.subgroup) leave demo phenotypes out by default
    st["stratum_id"] = st["phenotype_id"] + "|" + st["age_band"] + "|" + st["sex"]
    for c in ("n_cases", "n_positive"):
        st[c] = st[c].astype("Int64")
    for c in ("fraction", "ci_low", "ci_high"):
        st[c] = pd.to_numeric(st[c], errors="coerce").astype(float)
    _upsert(K.STRATA_TABLE, st, ds, src, now, "stratum_id",
            "Device-subgroup fractions among condition cases by age band x sex (aggregate; n < 10 suppressed)")
    echo(f"  wrote results/byod/{K.manifest_id_of(ds)}/computable_phenotype.json, SUBGROUP.md; tables "
         f"{K.PHENOTYPES_TABLE}, {K.STRATA_TABLE}")


def _upsert(table: str, new: pd.DataFrame, ds: str, src: str, now: str, rid: str, desc: str) -> None:
    old = read_table(table) if table_exists(table) else None
    if old is not None:
        old = old[old["dataset_id"] != ds].drop(columns=[c for c in old.columns if c in _PROV])
        new = pd.concat([old, new], ignore_index=True)
    write_aggregate_table(table, new, src, now, rid, desc)


def write_aggregate_table(table: str, df: pd.DataFrame, src: str, now: str, rid: str, desc: str) -> None:
    from ..config import PROCESSED
    if df.empty:
        for suf in (".parquet", ".meta.json"):
            (PROCESSED / f"{table}{suf}").unlink(missing_ok=True)
        return
    df = add_provenance(df.reset_index(drop=True), data_layer="measurement", source_name=src,
                        source_version=SUBGROUP_VERSION, retrieved_at=now, evidence_type="derived_evidence_summary",
                        source_record_id=rid, evidence_level=K.SOURCE_KIND,
                        provenance_notes=("aggregate only (no person rows, no participant ids, no geography); "
                                          f"{K.PROVENANCE_NOTE}; computed under a locked subgroup plan"))
    write_table(df, table, producer=PRODUCER, description=desc)


def drop_dataset(ds: str) -> None:
    """Remove a dataset's rows from the aggregate subgroup tables (used by `byod remove`)."""
    for table in (K.PHENOTYPES_TABLE, K.STRATA_TABLE):
        if not table_exists(table):
            continue
        t = read_table(table)
        rid = "phenotype_id" if table == K.PHENOTYPES_TABLE else "stratum_id"
        write_aggregate_table(table, t[t["dataset_id"] != ds].drop(columns=[c for c in t.columns if c in _PROV]),
                              f"measure_it.byod.remove ({ds} removed)", utc_now_iso(), rid,
                              f"{table} without the removed user dataset {ds}")


def _f(v, nd=3) -> str:
    try:
        return f"{float(v):.{nd}f}" if v is not None and np.isfinite(float(v)) else "n/a"
    except (TypeError, ValueError):
        return "n/a"


def subgroup_markdown(ds, manifest, plan, lock, phen, s) -> str:
    perf = phen["performance"]
    bp = phen["base_population"]
    lines = [f"# Device-defined subgroup: {plan['name']} ({ds})", "",
             ("**SYNTHETIC tutorial dataset: generated data, not real people.**\n\n" if manifest.get("synthetic")
              else "")
             + f"Subgroup plan locked {lock['locked_at']} (sha256 `{lock['sha256'][:16]}`; SUBGROUP_PLAN.md). "
               f"Base population: {bp['description']}; condition `{phen['condition_id']}`"
             + (f" (ICD-10-CM {', '.join(bp['icd10cm'])})" if bp["icd10cm"] else " (no ICD-10-CM code)") + ".",
             "", f"* {phen['subgroup_definition']}",
             f"* cases: {bp['n_cases']}; with a device finding: {bp['n_cases_with_device_finding_defined']} "
             f"({perf['n_positive']} positive, {perf['n_negative']} negative; pooled fraction "
             f"{_f(phen['strata_fractions'][0]['fraction'], 2)} [{_f(phen['strata_fractions'][0]['ci_low'], 2)}-"
             f"{_f(phen['strata_fractions'][0]['ci_high'], 2)}])",
             f"* harmonized domains: {', '.join(plan['ehr_domains'])}; {perf['n_candidate_features']} candidate "
             f"features after the locked filters", "",
             "## Can harmonized features find the device subgroup?", "",
             f"L1 logistic, {perf['cv']} CV: AUROC {_f(perf['auroc'])} ({_f(perf['auroc_ci_low'])}-"
             f"{_f(perf['auroc_ci_high'])}), AUPRC {_f(perf['auprc'])} (prevalence {_f(perf['auprc_prevalence_baseline'], 2)}), "
             f"permutation p {_f(perf['perm_p'], 4)} ({perf['n_perm']} shuffles): **{perf['decision']}**. "
             f"Decision threshold {_f(phen['decision_threshold'])} ({phen['decision_threshold_rule']}; sensitivity "
             f"{_f(perf['threshold_sensitivity_oof'], 2)}, specificity {_f(perf['threshold_specificity_oof'], 2)}).",
             "", "## Computable phenotype (refit on all cases)", "",
             f"Intercept {_f(phen['intercept'])}. {phen['ehr_evaluable_summary']['n_ehr_evaluable']} of "
             f"{phen['ehr_evaluable_summary']['n_features']} features can be evaluated in a clinic EHR "
             f"({_f(100 * phen['ehr_evaluable_summary']['ehr_evaluable_coefficient_share'], 0)}% of the coefficient "
             "mass).", "",
             "| feature | label | transform | coefficient (standardized) | EHR-evaluable | mapping |",
             "|---|---|---|---|---|---|"]
    for f in phen["features"]:
        lines.append(f"| `{f['feature_key']}` | {f['label']} | {f['transform']} | {_f(f['coefficient'])} "
                     f"({_f(f['coefficient_standardized'])}) | {'yes' if f['ehr_evaluable'] else 'no'} | "
                     f"{', '.join(f.get('mapping_methods') or [])} |")
    lines += ["", "## Fraction device-positive by stratum (aggregate; n < 10 suppressed)", "",
              "| age band | sex | n cases | n positive | fraction (95% CI) |", "|---|---|---|---|---|"]
    for r in phen["strata_fractions"]:
        if r["suppressed"]:
            lines.append(f"| {r['age_band']} | {r['sex']} | <{SUPPRESS_BELOW} | suppressed | suppressed |")
        else:
            lines.append(f"| {r['age_band']} | {r['sex']} | {r['n_cases']} | {r['n_positive']} | "
                         f"{_f(r['fraction'], 2)} ({_f(r['ci_low'], 2)}-{_f(r['ci_high'], 2)}) |")
    lines += ["", "Caveats:", ""] + [f"* {c}" for c in phen["caveats"]]
    lines += ["", "Files: `computable_phenotype.json` (contract section 2; scored elsewhere with "
              "`measure_it.harmonize.phenotype.score_phenotype`), `subgroup_features.csv`, "
              "`subgroup_feature_audit.csv` (every candidate feature and why it was kept or dropped).", ""]
    text = "\n".join(lines)
    from ..validate import forbidden_hits
    if forbidden_hits(text):
        import re
        text = re.sub(r"(?i)\b(proves|cures|patient has|best clinic|definitive biomarker|optimal treatment|"
                      r"confirms mechanism|diagnosed by AI)\b", "[term withheld]", text)
    return text
