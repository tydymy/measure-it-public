"""The metric -> translation link: measurement performance, implementation reach, expected yield.

Pre-specified in docs/ANALYSIS_PLAN_METRIC_LINK.md (written before any number here was computed); method reference
docs/SCORING.md section 9.

    build()                             -> data/processed/measurement_performance (+ measurement_performance_records),
                                           results/tables/measurement_performance*.csv, measurement_reach*.csv
    performance_for(condition, measurement) -> the scored performance record of one (condition, measurement)
    reach_for(level, implementer_groups)    -> share of each region's population within the default radius of an
                                               implementer-group facility

A performance **record** is one discrimination result (target condition, comparator, measurement class, AUROC and CI,
sensitivity at the pre-specified 0.90 specificity where it exists, n, label basis, quality tier). Records come only
from this project's own result tables (tiers 1-2) and, separately flagged, from `published_device_evidence` (tier 3,
never reproduced here). No value is ever carried from one condition, condition set, class or dataset to another: a
(condition, measurement) without a record is UNKNOWN.

Everything here describes measurements and places. Performance against healthy controls overstates real-world
performance, and an expected yield is a planning estimate for a candidate pilot, not a prediction of diagnoses.
"""
from __future__ import annotations

import json
import math
import sys
import time
from functools import lru_cache

import numpy as np
import pandas as pd

from ..config import SEED, TABLES, UNKNOWN, load_config, utc_now_iso
from ..provenance import add_provenance
from ..store import read_table, table_exists, write_table

PRODUCER = "measure_it.scoring.metric_link"
TABLE = "measurement_performance"
RECORDS_TABLE = "measurement_performance_records"
PLAN = "docs/ANALYSIS_PLAN_METRIC_LINK.md"
TARGET_SPECIFICITY = 0.90
N_BOOT = 2000
TIER_LABEL = {1: "own_computed_case_definition", 2: "own_computed_proxy_label", 3: "published_claim", 4: "none"}
TIER_TEXT = {1: "tier 1: this project's own computation on a clinical / study case definition",
             2: "tier 2: this project's own computation on a proxy or self-reported label",
             3: "tier 3: published claim, not reproduced by this project",
             4: "tier 4: no performance record (UNKNOWN)"}
TIER_FACTOR = {1: 1.0, 2: 0.75, 3: 0.5}
# Tier 2.5 (docs/BRING_YOUR_OWN_DATA.md): computed by this engine under a locked plan on user-supplied, local-only data
# (measure_it.byod). Ranked after tiers 1-2 (public data anyone can re-analyse) and before tier 3 (published claims never
# re-analysed); evidence factor 0.6. Present only when a user dataset has been evaluated (no effect otherwise).
USER_TIER = 2.5
TIER_LABEL[USER_TIER] = "user_supplied_own_computation"
TIER_TEXT[USER_TIER] = ("tier 2.5: computed by this engine under a locked plan on user-supplied, local-only data "
                        "(not public; not re-analysable by others)")
TIER_FACTOR[USER_TIER] = 0.6


def tier_value(q):
    """A stored quality tier as its dictionary key: int for 1-4, 2.5 for user-supplied records."""
    q = float(q if q is not None and q == q else 4)
    return int(q) if q.is_integer() else q


def tier_text(q) -> str:
    return TIER_TEXT.get(tier_value(q), TIER_TEXT[4])
STATUS_ORDER = {"known": 0, "partial": 1, UNKNOWN: 2}
EVIDENCE_COMPONENTS = ("measurement_evidence", "expected_yield")
QUERY_SET_ID = "long_covid_or_me_cfs"
SINGLE_CONDITIONS = ["long_covid", "me_cfs", "pots", "dysautonomia"]
BUNDLES = ["wearable_autonomic_activity_monitoring", "nailfold_capillaroscopy", "autonomic_function_testing",
           "exercise_capacity_testing"]
EXTRA_CLASSES = ["blood_biomarkers", "metabolomics"]
HEALTHY_CAVEAT = ("Performance measured against healthy (or recovered) controls overstates real-world performance: the "
                  "clinical question is this illness versus other causes of the same symptoms.")
YIELD_LANGUAGE = ("Expected yields are planning estimates for a candidate pilot, not predictions of diagnoses; where "
                  "prevalence is low, most flagged people are false positives.")

# ------------------------------------------------------------------------------------------------------------------
# small statistics (pure)
# ------------------------------------------------------------------------------------------------------------------


def wilson(p: float, n: float, z: float = 1.959964) -> tuple[float, float]:
    """Wilson score 95% interval of a proportion p observed in n (NaN if n is missing)."""
    if p is None or n is None or not np.isfinite(p) or not np.isfinite(n) or n <= 0:
        return (np.nan, np.nan)
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def hanley_mcneil_ci(auc: float, n1: float, n0: float, z: float = 1.959964) -> tuple[float, float]:
    """Hanley & McNeil (1982) 95% CI of an AUROC from its value and the case / control counts."""
    if not all(np.isfinite(x) for x in (auc, n1, n0)) or n1 <= 0 or n0 <= 0:
        return (np.nan, np.nan)
    q1, q2 = auc / (2 - auc), 2 * auc * auc / (1 + auc)
    var = (auc * (1 - auc) + (n1 - 1) * (q1 - auc ** 2) + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0)
    se = math.sqrt(max(var, 0.0))
    return (max(0.0, auc - z * se), min(1.0, auc + z * se))


def sens_at_spec(y: np.ndarray, score: np.ndarray, spec: float = TARGET_SPECIFICITY) -> tuple[float, float, float]:
    """(sensitivity, achieved specificity, threshold) at >= `spec` specificity for a score where higher = case.

    Threshold c = the ceil(spec * n_controls)-th smallest control score; flagged = score > c (plan section 2.4)."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(score, dtype=float)
    s0, s1 = np.sort(s[~y]), s[y]
    if len(s0) == 0 or len(s1) == 0:
        return (np.nan, np.nan, np.nan)
    k = int(math.ceil(spec * len(s0) - 1e-12))
    c = s0[min(max(k, 1), len(s0)) - 1]
    return (float((s1 > c).mean()), float((s0 <= c).mean()), float(c))


def sens_at_spec_boot(y: np.ndarray, score: np.ndarray, spec: float = TARGET_SPECIFICITY, n_boot: int = N_BOOT,
                      seed: int = SEED) -> dict:
    """Point estimate + stratified bootstrap 95% CI (threshold re-chosen in every resample)."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(score, dtype=float)
    se, sp, c = sens_at_spec(y, s, spec)
    rng = np.random.default_rng(seed)
    i1, i0 = np.flatnonzero(y), np.flatnonzero(~y)
    bs, bp = np.empty(n_boot), np.empty(n_boot)
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))])
        bs[b], bp[b], _ = sens_at_spec(y[idx], s[idx], spec)
    return {"op_sensitivity": se, "op_sensitivity_ci_low": float(np.percentile(bs, 2.5)),
            "op_sensitivity_ci_high": float(np.percentile(bs, 97.5)), "op_specificity": sp,
            "op_specificity_ci_low": float(np.percentile(bp, 2.5)),
            "op_specificity_ci_high": float(np.percentile(bp, 97.5)), "op_threshold": c}


def evidence_score(auroc_lower: float, tier: int) -> tuple[float, float]:
    """(tier-discounted, undiscounted) measurement_evidence = clip((AUROC_lower - 0.5) / 0.5, 0, 1) x tier factor."""
    if auroc_lower is None or not np.isfinite(auroc_lower) or tier not in TIER_FACTOR:
        return (np.nan, np.nan)
    raw = float(np.clip((auroc_lower - 0.5) / 0.5, 0.0, 1.0))
    return (raw * TIER_FACTOR[tier], raw)


def logit_sd(lo: float, hi: float) -> float:
    """SD on the logit scale implied by a 95% CI (bounds clipped to [0.001, 0.999])."""
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return 0.0
    lo, hi = np.clip([lo, hi], 0.001, 0.999)
    return float((np.log(hi / (1 - hi)) - np.log(lo / (1 - lo))) / (2 * 1.959964))


def draw_logit(p: float, lo: float, hi: float, rng: np.random.Generator, n: int) -> np.ndarray:
    """n draws of a proportion around p on the logit scale with the SD of its 95% CI."""
    if not np.isfinite(p):
        return np.full(n, np.nan)
    pc = float(np.clip(p, 0.001, 0.999))
    z = np.log(pc / (1 - pc)) + rng.standard_normal(n) * logit_sd(lo, hi)
    return 1.0 / (1.0 + np.exp(-z))


# ------------------------------------------------------------------------------------------------------------------
# records
# ------------------------------------------------------------------------------------------------------------------

def _rec(**kw) -> dict:
    base = {"record_id": None, "source_kind": None, "quality_tier": 4, "target_condition": None,
            "target_members": None, "primary_class": None, "dataset_id": None, "analysis": None, "label_basis": None,
            "comparator": None, "comparator_kind": None, "n_cases": np.nan, "n_controls": np.nan, "auroc": np.nan,
            "auroc_ci_low": np.nan, "auroc_ci_high": np.nan, "auroc_basis": UNKNOWN, "op_sensitivity": np.nan,
            "op_sensitivity_ci_low": np.nan, "op_sensitivity_ci_high": np.nan, "op_specificity": np.nan,
            "op_specificity_ci_low": np.nan, "op_specificity_ci_high": np.nan, "op_threshold": np.nan,
            "op_basis": UNKNOWN, "measurement_only": True, "exclusion_reason": None, "source_object_ids": [],
            "source_tables": [], "caveats": [], "computed_here": ""}
    base.update(kw)
    if base["target_members"] is None and base["target_condition"]:
        base["target_members"] = _members(base["target_condition"])
    lo = base["auroc_ci_low"]
    base["measurement_evidence"], base["measurement_evidence_undiscounted"] = evidence_score(lo, base["quality_tier"])
    has_auc = np.isfinite(lo)
    has_op = np.isfinite(base["op_sensitivity"]) and np.isfinite(base["op_specificity"])
    base["performance_status"] = "known" if (has_auc and has_op) else ("partial" if has_auc else UNKNOWN)
    base["quality_tier_label"] = TIER_LABEL[base["quality_tier"]]
    return base


def _members(cid: str) -> list[str]:
    if cid == QUERY_SET_ID:
        return ["long_covid", "me_cfs"]
    dc = load_config("conditions")["demo_cluster"]
    if cid == dc["id"]:
        return list(dc["members"])
    return sorted(set(cid.split("+"))) if "+" in cid else [cid]


def _csv(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / f"{name}.csv")


MM_LABEL = ("clinical / study case definition (MUSCLE-ME: Long COVID meeting the Canadian Consensus Criteria for "
            "ME/CFS; ME/CFS diagnosed before 2019)")
MM_COMPARATOR = "healthy controls (age- and sex-matched; recovered from SARS-CoV-2 without residual symptoms)"
MM_CAVEATS = ["Single-site case-control study (Amsterdam UMC); healthy controls only.",
              "Activity limitation / PEM is part of both case definitions (incorporation bias for steps).",
              "The Long COVID group meets the Canadian Consensus Criteria for ME/CFS: a subset of Long COVID.",
              "Research-grade hip ActiGraph, not a consumer wrist wearable.",
              "Probably the same people as Appelman 2024 (not an independent cohort)."]
MM_STEPS_SIG = "signature:charlton_lc_mecfs_cpet_source|{}|steps_mean_daily"


def _mm_frame() -> pd.DataFrame:
    from ..wearables.muscle_me_steps_analysis import load_frame
    return load_frame()


def own_records(compute: bool = True) -> list[dict]:
    """Tier 1-2 records from this project's result tables (plan section 2.1). compute=True runs the NEW computations
    (sensitivity at 0.90 specificity; the two per-condition VO2_rel AUROCs) on the stored person-level tables."""
    from ..wearables.muscle_me_steps_analysis import fixed_direction_auroc
    out = []
    prim = _csv("charlton_lc_mecfs_cpet_source_primary.csv".removesuffix(".csv"))
    sec = _csv("charlton_lc_mecfs_cpet_source_secondary")
    mm = _mm_frame() if compute else None
    worst = prim[prim["population"].astype(str).str.startswith("worst")].iloc[0]
    best = prim[prim["population"].astype(str).str.startswith("best")].iloc[0]

    def mm_rec(rid, target, cls, feature, sessions, prod_row, sig, extra_caveats=(), new=False):
        d = None
        kw = {}
        if mm is not None:
            d = mm[mm.session_code.isin(list(sessions) + ["CON"]) & mm[feature].notna()]
            y = d.is_patient.to_numpy()
            score = -d[feature].to_numpy(float)  # lower = case (fixed direction in the producer's plan)
            kw = sens_at_spec_boot(y, score)
            if new:
                r = fixed_direction_auroc(y, d[feature].to_numpy(float))
                prod_row = {"auroc": r["auroc"], "ci_low": r["ci_low"], "ci_high": r["ci_high"],
                            "n_cases": r["n_cases"], "n_controls": r["n_controls"]}
        feat_label = "daily steps (hip accelerometer)" if feature == "Steps" else "relative VO2peak (CPET)"
        return _rec(record_id=rid, source_kind="own_computed", quality_tier=1, target_condition=target,
                    primary_class=cls, dataset_id="charlton_lc_mecfs_cpet_source",
                    analysis=f"MUSCLE-ME {feat_label}, fixed direction (lower = case), nothing fitted",
                    label_basis=MM_LABEL, comparator=MM_COMPARATOR, comparator_kind="healthy",
                    n_cases=float(prod_row["n_cases"]), n_controls=float(prod_row["n_controls"]),
                    auroc=float(prod_row["auroc"]), auroc_ci_low=float(prod_row["ci_low"]),
                    auroc_ci_high=float(prod_row["ci_high"]),
                    auroc_basis=("computed here with the producer's function (NEW; plan section 2.1)" if new else
                                 "producer's 95% stratified bootstrap CI"),
                    op_basis=("threshold at >= 90% specificity on the controls (in-sample), stratified bootstrap CI "
                              "with the threshold re-chosen per resample (computed here)") if kw else UNKNOWN,
                    source_object_ids=[s for s in [sig] if s],
                    source_tables=["results/tables/charlton_lc_mecfs_cpet_source_primary.csv",
                                   "results/tables/charlton_lc_mecfs_cpet_source_secondary.csv",
                                   "data/processed/participant_exam_features__charlton_lc_mecfs_cpet_source",
                                   "data/processed/participant_wearable_features__charlton_lc_mecfs_cpet_source"],
                    caveats=MM_CAVEATS + list(extra_caveats),
                    computed_here=("AUROC and CI (NEW); " if new else "") + ("sensitivity at 0.90 specificity" if kw
                                                                            else ""),
                    **kw)

    p = prim[prim["analysis"] == "primary"].iloc[0]
    miss = (f"8 of 30 controls have no step data: worst case (missing controls = lowest observed steps) AUROC "
            f"{worst['auroc']:.3f} ({worst['ci_low']:.3f}-{worst['ci_high']:.3f}), CI includes 0.5; best case "
            f"{best['auroc']:.3f}; the missing controls are not missing completely at random (7 of 8 women).")
    out.append(mm_rec("OWN-MM-STEPS-POOLED", QUERY_SET_ID, "accelerometry", "Steps", ("LC", "ME"), p,
                      MM_STEPS_SIG.format("long_covid_or_me_cfs_vs_healthy"), [miss]))
    a = sec[sec["analysis"] == "(a)"]
    out.append(mm_rec("OWN-MM-STEPS-LC", "long_covid", "accelerometry", "Steps", ("LC",),
                      a[a["contrast"].str.startswith("Long COVID")].iloc[0],
                      MM_STEPS_SIG.format("long_covid_ccc_vs_healthy"), [miss.replace("30 controls", "30 controls "
                                                                                      "(pooled analysis)")]))
    out.append(mm_rec("OWN-MM-STEPS-ME", "me_cfs", "accelerometry", "Steps", ("ME",),
                      a[a["contrast"].str.startswith("ME/CFS")].iloc[0],
                      MM_STEPS_SIG.format("me_cfs_vs_healthy"), [miss.replace("30 controls", "30 controls "
                                                                               "(pooled analysis)")]))
    b = sec[sec["analysis"] == "(b)"].iloc[0]
    vo2_cav = ["Maximal CPET is burdensome and can provoke post-exertional malaise; deconditioning lowers VO2peak "
               "in any sedentary group."]
    out.append(mm_rec("OWN-MM-VO2-POOLED", QUERY_SET_ID, "cpet", "VO2_rel", ("LC", "ME"), b,
                      "signature:charlton_lc_mecfs_cpet_source|long_covid_or_me_cfs_vs_healthy|vo2_rel", vo2_cav))
    if mm is not None:
        out.append(mm_rec("OWN-MM-VO2-LC", "long_covid", "cpet", "VO2_rel", ("LC",), None,
                          "signature:charlton_lc_mecfs_cpet_source|long_covid_or_me_cfs_vs_healthy|vo2_rel", vo2_cav,
                          new=True))
        out.append(mm_rec("OWN-MM-VO2-ME", "me_cfs", "cpet", "VO2_rel", ("ME",), None,
                          "signature:charlton_lc_mecfs_cpet_source|long_covid_or_me_cfs_vs_healthy|vo2_rel", vo2_cav,
                          new=True))
    # ---- Uwakwe 2025: whole-record resting HR, L2 logistic (tier 2: self-reported label)
    perf = _csv("stanford_uwakwe_model_performance")
    u = perf[(perf["feature_set"] == "wearable") & (perf["model"] == "l2_logistic")].iloc[0]
    kw = {}
    if compute:
        oof = _csv("stanford_uwakwe_oof_predictions")
        o = oof[(oof["feature_set"] == "wearable") & (oof["model"] == "l2_logistic")]
        kw = sens_at_spec_boot(o["y"].to_numpy(), o["oof_prob_mean_over_repeats"].to_numpy(float))
    out.append(_rec(record_id="OWN-UW-RHR", source_kind="own_computed", quality_tier=2, target_condition="long_covid",
                    primary_class="wearable_heart_rate", dataset_id="stanford_longcovid_uwakwe2025",
                    analysis="Uwakwe 2025 whole-record resting-HR features (7), L2 logistic, 10x5 repeated CV",
                    label_basis="self-reported Long COVID (symptoms >= 12 weeks after confirmed COVID-19)",
                    comparator="participants without Long COVID after confirmed infection (recovered)",
                    comparator_kind="recovered after infection", n_cases=float(u["n_cases"]),
                    n_controls=float(u["n_controls"]), auroc=float(u["auroc"]), auroc_ci_low=float(u["auroc_ci_low"]),
                    auroc_ci_high=float(u["auroc_ci_high"]),
                    auroc_basis="producer's mean CV AUROC with participant-bootstrap 95% CI",
                    op_basis=("threshold at >= 90% specificity on the stored repeat-averaged out-of-fold probabilities "
                              "(in-sample threshold), stratified bootstrap CI (computed here)") if kw else UNKNOWN,
                    source_object_ids=["signature:stanford_longcovid_uwakwe2025|long_covid_self_report|"
                                       "model_auroc:wearable:l2_logistic"],
                    source_tables=["results/tables/stanford_uwakwe_model_performance.csv",
                                   "results/tables/stanford_uwakwe_oof_predictions.csv"],
                    caveats=["NULL result: whole-record resting HR does not separate self-reported Long COVID "
                             "(label-permutation p 0.30).", "n = 126 (31 cases); not aligned to infection; no steps; "
                             "not posture-linked HR and not beat-to-beat HRV."],
                    computed_here="sensitivity at 0.90 specificity" if kw else "", **kw))
    # ---- NHANES ME/CFS-like proxy, wearable-only pen_lr (tier 2: constructed proxy)
    nh = _csv("nhanes_model_performance")
    n = nh[(nh["target"] == "mecfs_like_proxy") & (nh["split"] == "cv_5x5") & (nh["feature_set"] == "wearable")
           & (nh["model"] == "pen_lr") & (nh["metric"] == "auroc")].iloc[0]
    out.append(_rec(record_id="OWN-NH-MECFS-PROXY", source_kind="own_computed", quality_tier=2,
                    target_condition="me_cfs", primary_class="accelerometry", dataset_id="nhanes",
                    analysis="NHANES 2011-2014 wrist accelerometry + sleep-proxy features (26), wearable-only "
                             "penalised logistic regression, 5x5 CV",
                    label_basis="constructed ME/CFS-LIKE PROXY (unexplained fatigue with functional limitation; not "
                                "ME/CFS; partly definitional for activity features)",
                    comparator="the rest of the wearable-valid NHANES adult sample (general population)",
                    comparator_kind="general population", n_cases=float(n["n_cases"]),
                    n_controls=float(n["n"] - n["n_cases"]), auroc=float(n["estimate"]),
                    auroc_ci_low=float(n["ci_low"]), auroc_ci_high=float(n["ci_high"]),
                    auroc_basis="producer's mean CV AUROC with participant-bootstrap 95% CI",
                    op_basis="UNKNOWN: no stored person-level predictions (not refitted here)",
                    source_object_ids=["measurement_signal:nhanes|mecfs_like_proxy|accelerometry|class_all"],
                    source_tables=["results/tables/nhanes_model_performance.csv"],
                    caveats=["Proxy label, not ME/CFS; the wearable adds nothing over demographics "
                             "(Delta AUROC +0.014, CI includes 0)."]))
    # ---- Appelman 2024 (tier 1: clinical Long COVID cohort vs healthy)
    ap = _csv("labs_dx_appelman_models")
    for rid, model, cls, lab, cav in (
            ("OWN-APP-VO2", "vo2max_only", "cpet", "Appelman 2024 VO2max-only CV model",
             ["Probably the same people as MUSCLE-ME (not independent).", "CV model; no stored predictions."]),
            ("OWN-APP-METAB", "blood_top10", "metabolomics", "Appelman 2024 plasma metabolomics, in-fold top-10 CV "
             "model (locked primary)", ["NULL by its locked rule (permutation p 0.063).", "25 vs 21 people."])):
        r = ap[ap["model"] == model].iloc[0]
        out.append(_rec(record_id=rid, source_kind="own_computed", quality_tier=1, target_condition="long_covid",
                        primary_class=cls, dataset_id="appelman_lc_pem_source", analysis=lab,
                        label_basis="clinical Long COVID cohort (Appelman 2024)", comparator="healthy controls",
                        comparator_kind="healthy", n_cases=float(r["n_lc"]), n_controls=float(r["n_healthy"]),
                        auroc=float(r["auroc"]), auroc_ci_low=float(r["auroc_ci_low"]),
                        auroc_ci_high=float(r["auroc_ci_high"]),
                        auroc_basis="producer's CV AUROC with stratified bootstrap 95% CI",
                        op_basis="UNKNOWN: no stored person-level predictions (not refitted here)",
                        source_object_ids=([f"signature:appelman_lc_pem_source|long_covid_vs_healthy__blood_baseline|"
                                            f"model_auroc:{model}"] if model == "blood_top10" else []),
                        source_tables=["results/tables/labs_dx_appelman_models.csv"], caveats=cav))
    # ---- MY-LC cortisol (tier 1)
    kl = _csv("klein2023_mylc_ml_table_cortisol_fitfree")
    for rid, key, comp, ck in (("OWN-KL-CORT", "S6", "healthy + convalescent controls (pooled)", "mixed"),
                               ("OWN-KL-CORT-CC", "S2 LC vs CC", "convalescent controls (recovered after COVID-19, "
                                                                 "no Long COVID)", "recovered after infection")):
        r = kl[kl["contrast"].astype(str).str.startswith(key)].iloc[0]
        se = float(r["sens_at_90spec"])
        lo, hi = wilson(se, float(r["n_cases"]))
        out.append(_rec(record_id=rid, source_kind="own_computed", quality_tier=1, target_condition="long_covid",
                        primary_class="blood_biomarkers", dataset_id="klein2023_mylc_ml_table",
                        analysis="MY-LC serum cortisol, fit-free (lower = case)",
                        label_basis="MY-LC research-cohort Long COVID definition", comparator=comp,
                        comparator_kind=ck, n_cases=float(r["n_cases"]), n_controls=float(r["n_controls"]),
                        auroc=float(r["auroc"]), auroc_ci_low=float(r["auroc_ci_low"]),
                        auroc_ci_high=float(r["auroc_ci_high"]), auroc_basis="producer's bootstrap 95% CI",
                        op_sensitivity=se, op_sensitivity_ci_low=lo, op_sensitivity_ci_high=hi,
                        op_specificity=TARGET_SPECIFICITY, op_specificity_ci_low=np.nan, op_specificity_ci_high=np.nan,
                        op_basis=("producer's sensitivity at 90% specificity; Wilson CI from that proportion and "
                                  "n_cases computed here (threshold-selection uncertainty not included)"),
                        source_object_ids=["signature:klein2023_mylc_ml_table|long_covid|model_auroc:log_cortisol:"
                                           "l2_logistic"],
                        source_tables=["results/tables/klein2023_mylc_ml_table_cortisol_fitfree.csv"],
                        caveats=["Draw time and demographics differ between groups (covariates alone AUROC 0.70).",
                                 "Cortisol is a blood test, not part of any scored bundle."],
                        computed_here="Wilson CI of the stored sensitivity"))
    return out


# published rows: (target condition, primary class, bundle attachment curated in plan section 2.2, notes)
PUBLISHED_MAP = {
    "PDE-006": ("pots", "ecg_ambulatory", None),
    "PDE-007": ("pots", "ecg_ambulatory", {"spec_n": 55, "comparator_kind": "look-alike (Long COVID without POTS) + "
                                                                     "healthy"}),
    "PDE-008": ("pots", "ecg_ambulatory", None),
    "PDE-002": ("pots", "autonomic_testing", None),
    "PDE-003": ("pots", "autonomic_testing", None),
    "PDE-022": ("me_cfs", "autonomic_testing", None),
    "PDE-014": ("me_cfs", "cpet", None),
    "PDE-040": ("long_covid", "wearable_heart_rate", {"exclude": "includes symptom features: not a measurement-only "
                                                                 "metric (and not reproducible from the public release)",
                                                      "comparator_kind": "recovered after infection"}),
    "PDE-043": ("long_covid", "hrv", None),
}


def published_records() -> list[dict]:
    """Tier 3 records from published_device_evidence (plan section 2.2); nothing here is reproduced."""
    if not table_exists("published_device_evidence"):
        return []
    pde = read_table("published_device_evidence").set_index("evidence_id")
    out = []
    for pid, (target, cls, opt) in PUBLISHED_MAP.items():
        if pid not in pde.index:
            continue
        r = pde.loc[pid]
        opt = opt or {}
        n1 = float(r["n_cases"]) if pd.notna(r["n_cases"]) else np.nan
        n0 = float(opt.get("spec_n", r["n_controls"])) if pd.notna(r["n_controls"]) else np.nan
        auc = float(r["auroc"]) if pd.notna(r["auroc"]) else np.nan
        se = float(r["sensitivity"]) if pd.notna(r["sensitivity"]) else np.nan
        sp = float(r["specificity"]) if pd.notna(r["specificity"]) else np.nan
        se_lo, se_hi = wilson(se, n1)
        sp_lo, sp_hi = wilson(sp, n0)
        if np.isfinite(auc):
            lo, hi = hanley_mcneil_ci(auc, n1, n0)
            basis = "published AUROC; 95% CI computed here (Hanley-McNeil from the published AUROC and counts)"
        elif np.isfinite(se) and np.isfinite(sp):
            auc, lo, hi = (se + sp) / 2, (se_lo + sp_lo) / 2, (se_hi + sp_hi) / 2
            basis = ("single-operating-point lower bound AUROC >= (se + sp) / 2; CI from Wilson intervals of the "
                     "published counts (computed here)")
        else:
            lo = hi = np.nan
            basis = UNKNOWN
        op = {}
        if np.isfinite(se) and np.isfinite(sp) and sp >= TARGET_SPECIFICITY:
            op = {"op_sensitivity": se, "op_sensitivity_ci_low": se_lo, "op_sensitivity_ci_high": se_hi,
                  "op_specificity": sp, "op_specificity_ci_low": sp_lo, "op_specificity_ci_high": sp_hi,
                  "op_basis": "published operating point (specificity >= 0.90; sensitivity at exactly 0.90 can only be "
                              "equal or higher); Wilson CIs from the published counts (computed here)"}
        elif np.isfinite(se):
            op = {"op_basis": (f"UNKNOWN: published specificity {sp:.2f} < 0.90" if np.isfinite(sp) else
                               "UNKNOWN: no published specificity") + " (never extrapolated)"}
        out.append(_rec(record_id=pid, source_kind="published_claim", quality_tier=3, target_condition=target,
                        primary_class=cls, dataset_id=f"{r['first_author']} {r['year']} (PMID {r['pmid']})",
                        analysis=str(r["device_or_test"]), label_basis=str(r["condition_definition_in_study"]),
                        comparator=str(r["control_group"]),
                        comparator_kind=opt.get("comparator_kind", "healthy"),
                        n_cases=n1, n_controls=float(r["n_controls"]) if pd.notna(r["n_controls"]) else np.nan,
                        auroc=auc, auroc_ci_low=lo, auroc_ci_high=hi, auroc_basis=basis,
                        measurement_only=not bool(opt.get("exclude")), exclusion_reason=opt.get("exclude"),
                        source_object_ids=[f"published_evidence:{pid}"],
                        source_tables=["data/processed/published_device_evidence"],
                        caveats=[str(r["claim_label"]), str(r["risk_of_bias_notes"])],
                        computed_here="CI / bound from published numbers" if basis != UNKNOWN else "", **op))
    return out


def user_supplied_records() -> list[dict]:
    """Tier 2.5 records of evaluated user datasets (measure_it.byod; demo records only with MEASURE_IT_BYOD_DEMO=1)."""
    from ..byod.records import performance_records
    return [_rec(**r) for r in performance_records()]


def build_records(compute: bool = True) -> pd.DataFrame:
    recs = own_records(compute=compute) + user_supplied_records() + published_records()
    df = pd.DataFrame(recs)
    for c in ("target_members", "source_object_ids", "source_tables", "caveats"):
        df[c] = df[c].map(json.dumps)
    return df


@lru_cache(maxsize=1)
def records() -> pd.DataFrame:
    """All performance records: the stored table when it exists, else computed now (same code)."""
    if table_exists(RECORDS_TABLE):
        return read_table(RECORDS_TABLE)
    return build_records()


# ------------------------------------------------------------------------------------------------------------------
# attachment and the scored record of one (condition, measurement)
# ------------------------------------------------------------------------------------------------------------------

def _measurement(mid: str) -> dict:
    from ..facilities import matching as M
    m = M.resolve_measurement(mid)
    if m.get("status") != "matched":
        raise KeyError(f"measurement {mid!r}: {m.get('reason')}")
    return m


def attach(condition_id: str, measurement_id: str, recs: pd.DataFrame | None = None) -> dict:
    """The performance row of (condition or set, measurement class or bundle); plan section 2.5."""
    recs = records() if recs is None else recs
    m = _measurement(measurement_id)
    members = sorted(_members(condition_id))
    cls = set(m["members"])
    tm = recs["target_members"].map(lambda s: sorted(json.loads(s)) if isinstance(s, str) else sorted(s))
    att = recs[(tm.map(lambda x: x == members)) & recs["primary_class"].isin(cls)].copy()
    usable = att[att["measurement_only"].astype(bool)].copy()
    no_metric = usable[usable["performance_status"] == UNKNOWN]
    usable = usable[usable["performance_status"] != UNKNOWN]
    row = {"condition_id": condition_id, "measurement_id": m["measurement_id"], "measurement_kind": m["kind"],
           "measurement_members": ";".join(m["members"]),
           "object_id": f"measurement_performance:{condition_id}|{m['measurement_id']}",
           "n_records": int(len(att)), "n_records_usable": int(len(usable)),
           "records": json.dumps([{k: (None if isinstance(v, float) and not np.isfinite(v) else v)
                                   for k, v in r.items() if k in ("record_id", "quality_tier", "primary_class",
                                                                   "auroc", "auroc_ci_low", "auroc_ci_high",
                                                                   "op_sensitivity", "op_specificity",
                                                                   "performance_status", "comparator_kind",
                                                                   "measurement_only", "exclusion_reason")}
                                  for r in att.to_dict("records")])}
    if usable.empty:
        row.update({"performance_status": UNKNOWN, "quality_tier": 4, "quality_tier_label": TIER_LABEL[4],
                    "selected_record_id": None, "measurement_evidence": np.nan,
                    "measurement_evidence_undiscounted": np.nan, "evidence_conflict": False, "evidence_conflict_note": None,
                    "performance_note": ("UNKNOWN / NOT AVAILABLE: no performance record of this project or of the "
                                         "curated published evidence for this condition and measurement"
                                         + (" with a usable metric (attached records without an AUROC or an "
                                            "operating point: " + "; ".join(no_metric["record_id"]) + ")"
                                            if len(no_metric) else "")
                                         + (" (records excluded as not measurement-only: "
                                            + "; ".join(att.loc[~att["measurement_only"].astype(bool), "record_id"])
                                            + ")" if (~att["measurement_only"].astype(bool)).any() else "")
                                         + ". Nothing is imputed from another condition or measurement.")})
        for c in ("dataset_id", "analysis", "label_basis", "comparator", "comparator_kind", "auroc_basis", "op_basis",
                  "source_kind"):
            row[c] = None
        for c in ("auroc", "auroc_ci_low", "auroc_ci_high", "op_sensitivity", "op_sensitivity_ci_low",
                  "op_sensitivity_ci_high", "op_specificity", "op_specificity_ci_low", "op_specificity_ci_high",
                  "n_cases", "n_controls"):
            row[c] = np.nan
        row["source_object_ids"] = "[]"
        row["source_tables"] = "[]"
        row["caveats"] = "[]"
        return row
    usable["_st"] = usable["performance_status"].map(STATUS_ORDER)
    usable["_me"] = usable["measurement_evidence"].fillna(np.inf)
    usable["_n"] = -usable["n_cases"].fillna(0)
    sel = usable.sort_values(["quality_tier", "_st", "_me", "_n", "record_id"]).iloc[0]
    conflict = False
    note = None
    if np.isfinite(sel["auroc_ci_low"]) and sel["auroc_ci_low"] > 0.5:
        other = usable[(usable["record_id"] != sel["record_id"]) & (usable["quality_tier"] <= USER_TIER)
                       & (usable["auroc_ci_low"] <= 0.5) & (usable["auroc_ci_high"] >= 0.5)]
        if len(other):
            conflict = True
            note = ("conflicting own-computed evidence: " + "; ".join(
                f"{o['record_id']} ({o['primary_class']}, {o['label_basis']}) AUROC {o['auroc']:.3f} "
                f"({o['auroc_ci_low']:.3f}-{o['auroc_ci_high']:.3f}) includes 0.5" for _, o in other.iterrows()))
    keep = ("dataset_id", "analysis", "label_basis", "comparator", "comparator_kind", "auroc", "auroc_ci_low",
            "auroc_ci_high", "auroc_basis", "op_sensitivity", "op_sensitivity_ci_low", "op_sensitivity_ci_high",
            "op_specificity", "op_specificity_ci_low", "op_specificity_ci_high", "op_basis", "n_cases", "n_controls",
            "quality_tier", "quality_tier_label", "measurement_evidence", "measurement_evidence_undiscounted",
            "performance_status", "source_object_ids", "source_tables", "caveats", "source_kind")
    row.update({k: sel[k] for k in keep})
    row["quality_tier"] = tier_value(row["quality_tier"])
    row["selected_record_id"] = sel["record_id"]
    row["evidence_conflict"] = conflict
    row["evidence_conflict_note"] = note
    st = row["performance_status"]
    row["performance_note"] = (
        f"{TIER_TEXT[row['quality_tier']]}; {sel['record_id']}: AUROC {sel['auroc']:.3f} "
        f"({sel['auroc_ci_low']:.3f}-{sel['auroc_ci_high']:.3f}) vs {sel['comparator']}"
        + (f"; sensitivity {sel['op_sensitivity']:.2f} at specificity {sel['op_specificity']:.2f}"
           if st == "known" else "; no operating point at >= 0.90 specificity (expected yield UNKNOWN)"
           if st == "partial" else "")
        + ". " + HEALTHY_CAVEAT)
    return row


def grid() -> list[tuple[str, str]]:
    from .opportunity import condition_sets
    conds = SINGLE_CONDITIONS + list(condition_sets())
    meas = list(BUNDLES)
    for b in BUNDLES:
        for c in _measurement(b)["members"]:
            if c not in meas:
                meas.append(c)
    meas += [c for c in EXTRA_CLASSES if c not in meas]
    cells = [(c, m) for c in conds for m in meas]
    from ..byod.records import extra_grid_cells   # user-dataset cells outside the grid (none without user data)
    return cells + extra_grid_cells(set(cells))


@lru_cache(maxsize=256)
def _performance_cached(condition_id: str, measurement_id: str, stamp: float) -> dict:
    if table_exists(TABLE):
        t = read_table(TABLE)
        hit = t[(t["condition_id"] == condition_id) & (t["measurement_id"] == measurement_id)]
        if len(hit):
            return hit.iloc[0].to_dict()
    return attach(condition_id, measurement_id)


def performance_for(condition_id: str, measurement_id: str) -> dict:
    """The scored performance row (stored table first; attached on the fly otherwise, same rules)."""
    from ..store import processed_path
    p = processed_path(TABLE)
    stamp = p.stat().st_mtime if p.exists() else 0.0
    mid = _measurement(measurement_id)["measurement_id"]
    return dict(_performance_cached(condition_id, mid, stamp))


# ------------------------------------------------------------------------------------------------------------------
# implementation reach
# ------------------------------------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def zcta_population() -> pd.DataFrame:
    z = read_table("zcta_centroids", columns=["zcta", "lat", "lon"])
    a = read_table("geo_context_zcta__acs", columns=["zcta", "state_fips", "county_fips_2024_of_internal_point",
                                                     "total_population"])
    d = a.merge(z, on="zcta", how="inner")
    d = d[d["lat"].notna() & d["total_population"].fillna(0).gt(0)].reset_index(drop=True)
    return d.rename(columns={"county_fips_2024_of_internal_point": "county_fips"})


@lru_cache(maxsize=32)
def _zcta_covered(groups: tuple, radius: float) -> np.ndarray:
    """Boolean per ZCTA: internal point within `radius` km of >= 1 geocoded clinic-candidate implementer facility."""
    from sklearn.neighbors import BallTree
    from . import opportunity as O
    S = O.spatial()
    fac_impl, _ = O._group_mask(S, list(groups))
    locs = np.unique(S.fac_loc[fac_impl])
    z = zcta_population()
    if len(locs) == 0:
        return np.zeros(len(z), bool)
    pts = np.radians(S.loc.loc[locs, ["lat", "lon"]].to_numpy(float))
    tree = BallTree(pts, metric="haversine")
    d, _ = tree.query(np.radians(z[["lat", "lon"]].to_numpy(float)), k=1)
    return (d[:, 0] * 6371.0088) <= radius


def reach_for(level: str, groups, geo: pd.DataFrame | None = None) -> pd.DataFrame:
    """reach, reach_population, zcta_population_total, reach_zcta_n per region (aligned to the level's universe)."""
    from . import opportunity as O
    geo = O.spatial().levels[level].geo if geo is None else geo
    z = zcta_population()
    cov = _zcta_covered(tuple(sorted(groups)), O.radius_km())
    key = "county_fips" if level == "county" else "state_fips"
    g = pd.DataFrame({"k": z[key].astype(str), "pop": z["total_population"].astype(float), "cov": cov})
    agg = g.groupby("k").agg(pop_total=("pop", "sum"), n=("pop", "size"),
                             pop_cov=("pop", lambda s: float(s[g.loc[s.index, "cov"]].sum())))
    ids = geo["geo_id"].astype(str) if level == "county" else geo["state_fips"].astype(str)
    a = agg.reindex(ids.to_numpy())
    with np.errstate(invalid="ignore", divide="ignore"):
        reach = (a["pop_cov"] / a["pop_total"]).to_numpy(float)
    return pd.DataFrame({"geo_id": geo["geo_id"].to_numpy(), "reach": reach,
                         "reach_population": a["pop_cov"].to_numpy(float),
                         "reach_zcta_population_total": a["pop_total"].to_numpy(float),
                         "reach_zcta_n": a["n"].fillna(0).to_numpy(float)})


# ------------------------------------------------------------------------------------------------------------------
# build
# ------------------------------------------------------------------------------------------------------------------

def build(write: bool = True, write_csv: bool = True) -> dict:
    """write_csv=False (measure-it byod deploy / remove) leaves the committed results/tables CSVs untouched."""
    t0 = time.time()
    rec = build_records()
    rows = [attach(c, m, rec) for c, m in grid()]
    perf = pd.DataFrame(rows)
    built = utc_now_iso()
    perf = add_provenance(
        perf, data_layer="measurement",
        source_name=f"{PRODUCER} (this project's result tables: MUSCLE-ME, Stanford Uwakwe, NHANES, Appelman, MY-LC; "
                    "published_device_evidence flagged as tier 3)",
        source_version=f"metric link {built}; plan {PLAN}", retrieved_at=built,
        evidence_type="derived_evidence_summary", source_record_id="object_id", source_geographic_resolution="none",
        evidence_level=perf["quality_tier_label"].astype(str).to_numpy(),
        provenance_notes=("Performance of a measurement for a condition, drawn only from this project's own computed "
                          "results (tiers 1-2) and, flagged, from published claims (tier 3). UNKNOWN where nothing "
                          "exists; never imputed from a neighbouring condition or measurement. " + HEALTHY_CAVEAT))
    rec_t = add_provenance(
        rec, data_layer="measurement", source_name=f"{PRODUCER} (performance records)",
        source_version=f"metric link {built}; plan {PLAN}", retrieved_at=built,
        evidence_type="derived_evidence_summary", source_record_id="record_id", source_geographic_resolution="none",
        evidence_level=rec["quality_tier_label"].astype(str).to_numpy(),
        provenance_notes="One discrimination result per row; tier 3 rows are published claims not reproduced here.")
    # reach per scored bundle x level (context table; the scoring computes the same function)
    from . import opportunity as O
    reach_rows = []
    for level in O.LEVELS:
        for b in BUNDLES:
            r = reach_for(level, _measurement(b)["implementer_groups"])
            r.insert(0, "measurement_id", b)
            r.insert(1, "geo_level", level)
            reach_rows.append(r)
    reach = pd.concat(reach_rows, ignore_index=True)
    if write:
        assert perf["object_id"].is_unique
        write_table(perf, TABLE, producer=PRODUCER,
                    description="Measurement performance per (condition or set, measurement class or bundle): the "
                                "scored record, its tier, AUROC and CI, sensitivity at 0.90 specificity, "
                                "measurement_evidence; UNKNOWN where no record exists. docs/ANALYSIS_PLAN_METRIC_LINK.md")
        write_table(rec_t, RECORDS_TABLE, producer=PRODUCER,
                    description="All performance records (own-computed tiers 1-2; published claims tier 3).")
        if write_csv:
            TABLES.mkdir(parents=True, exist_ok=True)
            perf.drop(columns=[c for c in perf.columns if c in ("records",)]).to_csv(
                TABLES / "measurement_performance.csv", index=False)
            rec.to_csv(TABLES / "measurement_performance_records.csv", index=False)
            reach.to_csv(TABLES / "measurement_reach.csv", index=False)
        records.cache_clear()
        _performance_cached.cache_clear()
    print(f"[metric_link] {len(rec)} records, {len(perf)} condition x measurement rows "
          f"({(perf['performance_status'] == 'known').sum()} known, {(perf['performance_status'] == 'partial').sum()} "
          f"partial, {(perf['performance_status'] == UNKNOWN).sum()} UNKNOWN) ({time.time() - t0:.0f}s)",
          file=sys.stderr, flush=True)
    return {"records": rec, "performance": perf, "reach": reach}


def run() -> dict:
    return build()


if __name__ == "__main__":
    run()
