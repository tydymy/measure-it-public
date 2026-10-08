"""Part 1 — NHANES 2011-2014 laboratory analytes vs diagnosis-like labels.

Follows docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md (locked by `measure_it.labs.plan`; this module refuses to run if the
plan was edited). Person-level only: every join is on SEQN (participant_id = nhanes:<SEQN>).

Stages (each writes results/tables/labs_dx_nhanes_*.csv):
  effects   per-analyte age/sex(/BMI)-adjusted standardised differences (BH-FDR within label x tier)
  models    penalised logistic CV (demographics vs + labs), bootstrap CIs, label-permutation null, lab-block shuffle
  signatures  phenotype_signatures__nhanes_labs

    uv run python -m measure_it.labs.labs_vs_diagnosis [--stages effects models signatures] [--jobs 16]
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
import zlib
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from ..config import SEED, TABLES, utc_now_iso
from .plan import check_plan_locked

PRODUCER = "measure_it.labs.labs_vs_diagnosis"
DATASET_ID = "nhanes_labs"
SOURCE_ID = "nhanes_2011_2014"
PREFIX = "labs_dx_nhanes"
FDR_ALPHA = 0.05
UNDERPOWERED_CASES = 50
MIN_CASES_EFFECT = 10
MIN_CASES_MODEL = 25
MIN_POS_BINARY = 5
N_SPLITS, N_REPEATS = 5, 5
N_BOOT = 1000
N_PERM = 200
N_SHUFFLE = 20
CS = (1e-3, 1e-2, 1e-1, 1.0, 10.0)
PRIMARY_MIN_NONMISSING = 0.70
MAX_BELOW_LOD = 0.50

CORE_FULL_SAMPLE_DROP = {"LBXSCH": "duplicate of total cholesterol (LBXTC)"}
PRIMARY_EXTENDED_LAYERS = ("urine_albumin_creatinine", "nutrition_folate_b12_status", "hormones_sex_steroids",
                           "tobacco_cotinine")
GLYCAEMIC = {"LBXGH", "LBXSGL", "LBXGLU", "LBXIN", "LBXIN__metabolic_fasting_subsample", "LBXGLT"}
URINE_CREATININE = "URXUCR"

FEATURE_SETS = ("demographics", "demographics_bmi", "labs", "demographics_labs", "demographics_bmi_labs")
COMPARISONS = {  # name -> (with, without)
    "demographics_labs_minus_demographics": ("demographics_labs", "demographics"),
    "demographics_bmi_labs_minus_demographics_bmi": ("demographics_bmi_labs", "demographics_bmi"),
    "labs_minus_demographics": ("labs", "demographics"),
}
PRIMARY_COMPARISON = "demographics_labs_minus_demographics"


# ======================================================================================================== labels
@dataclass(frozen=True)
class Label:
    label_id: str
    group: str            # umbrella / comparator / positive_control / pairwise
    title: str
    definition: str
    label_basis: str
    is_proxy: bool = False
    condition_id: str | None = None
    drop_analytes: frozenset = frozenset()
    cycles: str = "G,H"
    caveat: str = ""


RX_CAVEAT = ("Rx reason-code label: participants taking a prescription for this reason (2013-2014 only); lab "
             "differences may reflect the drugs as much as the condition.")
LABELS: list[Label] = [
    Label("rx_fibromyalgia", "umbrella", "Fibromyalgia (Rx reason code M79.7)",
          "any prescription with reason-for-use ICD-10-CM M79.7", "rx_reason_code", condition_id="fibromyalgia",
          cycles="H", caveat=RX_CAVEAT),
    Label("rx_migraine", "umbrella", "Migraine (Rx reason code G43)", "any prescription with reason-for-use G43",
          "rx_reason_code", condition_id="migraine", cycles="H", caveat=RX_CAVEAT),
    Label("rx_ibs", "umbrella", "Irritable bowel syndrome (Rx reason code K58)",
          "any prescription with reason-for-use K58", "rx_reason_code", condition_id="ibs", cycles="H",
          caveat=RX_CAVEAT),
    Label("rx_insomnia", "umbrella", "Insomnia (Rx reason code G47.0/F51.0)",
          "any prescription with reason-for-use G47.0 or F51.0", "rx_reason_code", cycles="H", caveat=RX_CAVEAT),
    Label("rx_myalgia", "umbrella", "Myalgia (Rx reason code M79.1)", "any prescription with reason-for-use M79.1",
          "rx_reason_code", cycles="H",
          caveat=RX_CAVEAT + " Statin-associated muscle pain is a likely contributor."),
    Label("mecfs_like_proxy", "umbrella", "Unexplained fatigue with functional limitation (ME/CFS-LIKE PROXY; not ME/CFS)",
          "DPQ040 >= 2 AND PFQ functional-limitation composite AND no affirmative answer to any exclusion diagnosis "
          "(nhanes_cohort.PROXY_EXCLUSION_ITEMS); controls: everyone else with DPQ040 and PFQ answered",
          "constructed_proxy", is_proxy=True, condition_id="me_cfs",
          caveat="PROXY, not an ME/CFS case definition (no PEM, duration or clinician assessment). The exclusion list "
                 "removes people reporting diabetes, thyroid problems, anaemia treatment, liver or kidney-adjacent "
                 "conditions, so labs linked to those are partly defined out of the cases."),
    Label("fatigue_symptom", "umbrella", "Fatigue symptom (PHQ-9 item 4 >= 2)",
          "DPQ040 ('tired or little energy') more than half the days", "self_report_symptom_item", is_proxy=True),
    Label("told_sleep_disorder", "umbrella", "Told by a doctor: sleep disorder (SLQ060)", "SLQ060 = yes vs no",
          "self_report_told_by_health_professional"),
    Label("depression_phq9", "umbrella", "Depressive symptoms (PHQ-9 >= 10)", "PHQ-9 total >= 10 (complete PHQ-9)",
          "questionnaire_screen", is_proxy=True),
    Label("celiac_disease", "comparator", "Celiac disease (self-reported diagnosis)", "MCQ082 = yes vs no",
          "self_reported_diagnosis"),
    Label("thyroid_problem", "comparator", "Thyroid problem (self-reported diagnosis)", "MCQ160M = yes vs no",
          "self_reported_diagnosis"),
    Label("rheumatoid_arthritis", "comparator", "Rheumatoid arthritis (self-reported type)",
          "MCQ195 = rheumatoid arthritis; controls MCQ160A = no arthritis", "self_reported_diagnosis"),
    Label("osteoarthritis", "comparator", "Osteoarthritis (self-reported type)",
          "MCQ195 = osteoarthritis; controls MCQ160A = no arthritis", "self_reported_diagnosis"),
    Label("psoriasis", "comparator", "Psoriasis (self-reported diagnosis)", "MCQ070 = yes vs no",
          "self_reported_diagnosis"),
    Label("asthma", "comparator", "Asthma (self-reported diagnosis)", "MCQ010 = yes vs no", "self_reported_diagnosis"),
    Label("diabetes_no_glycaemic_labs", "positive_control", "Diabetes, glycaemic labs removed (positive control)",
          "DIQ010 = yes vs no (borderline excluded); HbA1c, serum/fasting glucose, insulin, OGTT removed",
          "self_reported_diagnosis", drop_analytes=frozenset(GLYCAEMIC)),
    Label("diabetes_with_glycaemic_labs", "positive_control", "Diabetes, all labs (trivial sanity check)",
          "DIQ010 = yes vs no (borderline excluded); all analytes", "self_reported_diagnosis"),
    Label("weak_failing_kidneys", "positive_control", "Weak/failing kidneys (KIQ022; positive control)",
          "KIQ022 = yes vs no", "self_reported_diagnosis"),
    Label("gout", "positive_control", "Gout (self-reported diagnosis; positive control)", "MCQ160N = yes vs no",
          "self_reported_diagnosis"),
    Label("congestive_heart_failure", "positive_control", "Congestive heart failure (positive control)",
          "MCQ160B = yes vs no", "self_reported_diagnosis"),
    Label("anemia_treatment", "positive_control", "Anaemia treatment in past 3 months (positive control)",
          "MCQ053 = yes vs no", "self_reported_item"),
    Label("liver_condition", "positive_control", "Liver condition (self-reported; positive control)",
          "MCQ160L = yes vs no", "self_reported_diagnosis"),
]
PAIRS = [  # (label_id, case source, control source, title)
    ("pair_rx_fibromyalgia_vs_rx_migraine", "rx_fibromyalgia", "rx_migraine",
     "Rx fibromyalgia (case) vs Rx migraine (reference)"),
    ("pair_rx_fibromyalgia_vs_rx_insomnia", "rx_fibromyalgia", "rx_insomnia",
     "Rx fibromyalgia (case) vs Rx insomnia (reference)"),
    ("pair_rx_migraine_vs_rx_insomnia", "rx_migraine", "rx_insomnia", "Rx migraine (case) vs Rx insomnia (reference)"),
    ("pair_rx_fibromyalgia_vs_rheumatoid_arthritis", "rx_fibromyalgia", "rheumatoid_arthritis",
     "Rx fibromyalgia (case) vs rheumatoid arthritis (reference; 2013-2014 only)"),
    ("pair_mecfs_like_proxy_vs_depression_non_proxy", "mecfs_like_proxy", "depression_phq9",
     "ME/CFS-like PROXY (case) vs PHQ-9 >= 10 non-proxy (reference)"),
]
for pid, a, b, title in PAIRS:
    la = next(x for x in LABELS if x.label_id == a)
    LABELS.append(Label(pid, "pairwise", title, f"cases of {a} vs cases of {b}; people in both dropped",
                        "pairwise_contrast", is_proxy=la.is_proxy or a == "mecfs_like_proxy",
                        cycles="H" if a.startswith("rx_") else "G,H", caveat=la.caveat))
LABEL_BY_ID = {x.label_id: x for x in LABELS}


def _yes_no(s: pd.Series) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    return v.where(v.isin([0, 1]))


def _kiq022() -> pd.DataFrame:
    frames = []
    from ..download import download_file
    for c, year in (("G", 2011), ("H", 2013)):
        # KIQ_U (kidney conditions questionnaire) is not fetched by measure_it.ingestion.nhanes; downloaded here into
        # the same raw directory (MANIFEST.json records it; offline runs read the recorded file)
        p = download_file(f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/KIQ_U_{c}.xpt", SOURCE_ID,
                          min_bytes=10_000)
        d = pd.read_sas(p)[["SEQN", "KIQ022"]]
        frames.append(d)
    k = pd.concat(frames)
    k["participant_id"] = "nhanes:" + k["SEQN"].astype(int).astype(str)
    k["weak_failing_kidneys"] = k["KIQ022"].map({1.0: 1.0, 2.0: 0.0})
    return k[["participant_id", "weak_failing_kidneys"]]


def load_population() -> tuple[pd.DataFrame, list[dict]]:
    """Adults >= 20, not pregnant, with the biochemistry profile; label columns y_<label_id>."""
    from ..store import read_table
    from ..wearables.nhanes_cohort import EDUCATION_ORDER, PROXY_EXCLUSION_ITEMS, RACE_LEVELS

    clin = read_table("participant_clinical_features__nhanes")
    parts = read_table("participants__nhanes", columns=["participant_id", "pregnant_at_exam"])
    df = clin.merge(parts, on="participant_id", how="left", validate="one_to_one")
    flow = [{"step": "NHANES 2011-2014 participants (DEMO)", "n": int(len(df))}]
    df = df[df["age_years"] >= 20]
    flow.append({"step": "aged >= 20", "n": int(len(df))})
    preg = df["pregnant_at_exam"].fillna(False).astype(bool)
    df = df[~preg]
    flow.append({"step": "not pregnant at MEC exam", "n": int(len(df))})
    df = df[df["lab_albumin_g_dl"].notna()].copy()
    flow.append({"step": "POPULATION: biochemistry profile measured", "n": int(len(df))})
    df = df.merge(_kiq022(), on="participant_id", how="left", validate="one_to_one")

    df["female"] = (df["sex"] == "female").astype(float)
    for i, lvl in enumerate(RACE_LEVELS):
        df[f"race_{i}"] = (df["race_ethnicity"] == lvl).astype(float)
    df["education_ord"] = df["education_adult20plus"].map(EDUCATION_ORDER).astype(float)
    is_h = df["cycle_code"] == "H"

    # rx reason-code groups (H only)
    cond = read_table("participant_conditions__nhanes", columns=["participant_id", "rx_target_group"])
    cond = cond[cond["rx_target_group"].fillna("") != ""]
    rx_ok = is_h & df["rx_any_past_month"].notna()
    for grp, lid in (("fibromyalgia", "rx_fibromyalgia"), ("migraine", "rx_migraine"),
                     ("irritable bowel syndrome", "rx_ibs"), ("insomnia", "rx_insomnia"), ("myalgia", "rx_myalgia")):
        ids = set(cond.loc[cond["rx_target_group"] == grp, "participant_id"])
        df[f"y_{lid}"] = np.where(rx_ok, df["participant_id"].isin(ids).astype(float), np.nan)

    # ME/CFS-like proxy (same definition as nhanes_cohort, in this population)
    lim = df["any_functional_limitation_pfq"]
    fat = df["phq9_item4_tired_little_energy_0to3"]
    df["ra_or_psoriatic_arthritis"] = df["arthritis_type"].fillna("").str.contains(
        "heumatoid|soriatic", regex=True).astype(float)
    excl = pd.Series(0.0, index=df.index)
    for c in PROXY_EXCLUSION_ITEMS:
        excl = np.maximum(excl, (df[c] == 1).astype(float))
    ok = fat.notna() & lim.notna()
    case = (fat >= 2) & (lim == 1) & (excl == 0)
    df["y_mecfs_like_proxy"] = np.where(ok, case.astype(float), np.nan)
    df["proxy_any_exclusion_dx"] = excl
    df["y_fatigue_symptom"] = np.where(fat.notna(), (fat >= 2).astype(float), np.nan)
    df["y_told_sleep_disorder"] = _yes_no(df["told_sleep_disorder"])
    df["y_depression_phq9"] = _yes_no(df["phq9_ge10"])
    df["y_celiac_disease"] = _yes_no(df["dx_celiac_disease"])
    df["y_thyroid_problem"] = _yes_no(df["dx_thyroid_problem"])
    no_arth = df["dx_arthritis"] == 0
    at = df["arthritis_type"].fillna("")
    df["y_rheumatoid_arthritis"] = np.where(at == "rheumatoid arthritis", 1.0, np.where(no_arth, 0.0, np.nan))
    df["y_osteoarthritis"] = np.where(at == "osteoarthritis", 1.0, np.where(no_arth, 0.0, np.nan))
    df["y_psoriasis"] = _yes_no(df["dx_psoriasis"])
    df["y_asthma"] = _yes_no(df["dx_asthma"])
    dia = _yes_no(df["dx_diabetes"]).where(~((df["dx_diabetes"] == 0) & (df["dx_borderline_diabetes"] == 1)))
    df["y_diabetes_no_glycaemic_labs"] = dia
    df["y_diabetes_with_glycaemic_labs"] = dia
    df["y_weak_failing_kidneys"] = df["weak_failing_kidneys"]
    df["y_gout"] = _yes_no(df["dx_gout"])
    df["y_congestive_heart_failure"] = _yes_no(df["dx_congestive_heart_failure"])
    df["y_anemia_treatment"] = _yes_no(df["anemia_treatment_past_3_months"])
    df["y_liver_condition"] = _yes_no(df["dx_liver_condition"])

    # pairwise contrasts
    for pid, a, b, _ in PAIRS:
        ya, yb = df[f"y_{a}"] == 1, df[f"y_{b}"] == 1
        if pid == "pair_mecfs_like_proxy_vs_depression_non_proxy":
            yb = (df["y_depression_phq9"] == 1) & (df["y_mecfs_like_proxy"] == 0)
        if pid == "pair_rx_fibromyalgia_vs_rheumatoid_arthritis":
            yb = yb & is_h
        both = ya & yb
        df[f"y_{pid}"] = np.where(ya & ~both, 1.0, np.where(yb & ~both, 0.0, np.nan))
    return df.reset_index(drop=True), flow


# ======================================================================================================== analytes
def _skew_transform(x: pd.Series) -> tuple[pd.Series, str]:
    v = x.dropna()
    if len(v) < 20 or v.min() < 0 or stats.skew(v) <= 1:
        return x, "none"
    pos = v[v > 0]
    c = 0.0 if v.min() > 0 else (0.5 * pos.min() if len(pos) else 1.0)
    return np.log(x + c), ("log" if c == 0 else f"log(x+{c:.4g})")


def load_analytes(pop_ids: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(continuous wide matrix [transformed], binary wide matrix, analyte dictionary) for the population."""
    from ..store import read_table

    ids = set(pop_ids)
    core = read_table("participant_labs__nhanes",
                      columns=["participant_id", "lab_variable", "value", "lab_name", "unit", "subsample",
                               "component_file", "lab_note"])
    core = core[core["participant_id"].isin(ids)].copy()
    core.loc[core["lab_variable"].isin(["LBXB12", "LBDB12"]), "lab_variable"] = "VITB12"
    core_meta = (core.drop_duplicates("lab_variable")
                 .assign(layer=lambda d: np.where(d["subsample"].str.startswith("full MEC"), "core_full_sample",
                                                  "core_subsample"),
                         below_lod_share=np.nan, component=lambda d: d["component_file"].str[:-2])
                 [["lab_variable", "lab_name", "unit", "layer", "subsample", "component"]])
    core_meta.loc[core_meta["lab_variable"] == "VITB12", "lab_name"] = (
        "Vitamin B12 (2011-12 LBXB12; 2013-14 Deming-adjusted LBDB12)")
    core_w = core.pivot_table(index="participant_id", columns="lab_variable", values="value", aggfunc="mean")

    ext = read_table("participant_labs_extended__nhanes",
                     columns=["participant_id", "lab_layer", "component_file", "component_description",
                              "lab_variable", "lab_name", "unit", "analyte_role", "value_type", "value",
                              "result_positive", "below_lod", "special_sample", "subsample"])
    ext = ext[ext["participant_id"].isin(ids) & ~ext["special_sample"]]
    ucr = (ext[(ext["lab_variable"] == URINE_CREATININE) & (ext["lab_layer"] == "urine_albumin_creatinine")]
           .groupby("participant_id")["value"].mean().rename("_urine_creatinine"))
    ext = ext[ext["analyte_role"] == "analyte"].copy()
    # NHANES reuses variable names across files (LBXHCT = hematocrit in CBC, hydroxycotinine in COT_H): extended
    # variables that collide with a core variable are namespaced <var>__<layer>
    clash = ext["lab_variable"].isin(set(core["lab_variable"]))
    ext.loc[clash, "lab_variable"] = ext.loc[clash, "lab_variable"] + "__" + ext.loc[clash, "lab_layer"]
    cont = ext[ext["value_type"] == "continuous"]
    binr = ext[(ext["value_type"] == "categorical") & ext["result_positive"].notna()]
    # a variable released as categorical in one cycle and numeric in the other (LBDRPI) is analysed as binary only
    cont = cont[~cont["lab_variable"].isin(set(binr["lab_variable"]))]
    lod = cont.groupby("lab_variable")["below_lod"].apply(lambda s: float(s.fillna(False).astype(bool).mean()))
    emeta = (cont.drop_duplicates("lab_variable")
             [["lab_variable", "lab_name", "unit", "lab_layer", "subsample", "component_file", "component_description"]]
             .rename(columns={"lab_layer": "layer", "component_file": "component"}))
    emeta["component"] = emeta["component"].str.replace(r"_[GH]$", "", regex=True)
    emeta["below_lod_share"] = emeta["lab_variable"].map(lod)
    emeta["urinary"] = (emeta["component_description"].str.contains("Urine", case=False)
                        & ~emeta["lab_name"].str.contains("ratio|creatinine", case=False))
    ext_w = cont.pivot_table(index="participant_id", columns="lab_variable", values="value", aggfunc="mean")
    bmeta = (binr.drop_duplicates("lab_variable")
             [["lab_variable", "lab_name", "unit", "lab_layer", "subsample", "component_file"]]
             .rename(columns={"lab_layer": "layer", "component_file": "component"}))
    bin_w = binr.pivot_table(index="participant_id", columns="lab_variable", values="result_positive",
                             aggfunc="max")

    wide = pd.DataFrame(index=pd.Index(pop_ids, name="participant_id"))
    assert not set(ext_w.columns) & set(core_w.columns)
    wide = pd.concat([wide, core_w.reindex(wide.index), ext_w.reindex(wide.index)], axis=1)
    wide = wide.join(ucr)
    bin_w = pd.DataFrame(index=wide.index).join(bin_w)

    rows = []
    seen_names: dict[str, list[str]] = {}
    n = len(wide)
    for _, m in pd.concat([core_meta.assign(urinary=False), emeta], ignore_index=True).iterrows():
        v = m["lab_variable"]
        if v not in wide.columns:
            continue
        nonmiss = float(wide[v].notna().mean())
        lodshare = m.get("below_lod_share", np.nan)
        tier = "primary" if (m["layer"] == "core_full_sample" or m["layer"] in PRIMARY_EXTENDED_LAYERS) else "extended"
        note = ""
        if tier == "primary" and nonmiss < PRIMARY_MIN_NONMISSING:
            # failing the primary-panel coverage rule: screened in the extended tier instead (deviation, documented)
            tier, note = "extended", f"moved from the primary panel: < {PRIMARY_MIN_NONMISSING:.0%} non-missing"
        dup = unit_duplicate_of(wide, v, m["lab_name"], seen_names)
        reason = ""
        if v in CORE_FULL_SAMPLE_DROP:
            reason = CORE_FULL_SAMPLE_DROP[v]
        elif m["layer"] == "urine_flow_rate":
            reason = "urine flow-rate process variable, not an analyte (plan)"
        elif dup:
            reason = f"unit duplicate of {dup} (same analyte label, rank correlation > 0.999)"
        elif pd.notna(lodshare) and lodshare >= MAX_BELOW_LOD:
            reason = f">= {MAX_BELOW_LOD:.0%} below the detection limit"
        elif wide[v].notna().sum() < 50 or wide[v].nunique() < 5:
            reason = "fewer than 50 measured or 5 distinct values in the population"
        wide[v], tr = _skew_transform(wide[v]) if not reason else (wide[v], "none")
        rows.append({"analyte": v, "lab_name": m["lab_name"], "unit": m["unit"], "layer": m["layer"],
                     "component": m["component"], "subsample": m["subsample"], "tier": tier,
                     "value_type": "continuous", "urinary_creatinine_adjusted": bool(m.get("urinary", False)),
                     "n_measured": int(wide[v].notna().sum()), "share_nonmissing": round(nonmiss, 4),
                     "share_below_lod": lodshare, "transform": tr, "included": not reason,
                     "exclusion_reason": reason, "note": note, "n_population": n})
        if not reason:
            seen_names.setdefault(str(m["lab_name"]).strip().lower(), []).append(v)
    for _, m in bmeta.iterrows():
        v = m["lab_variable"]
        pos = bin_w[v].sum()
        reason = "" if pos >= 2 * MIN_POS_BINARY else "fewer than 10 positives in the population"
        rows.append({"analyte": v, "lab_name": m["lab_name"], "unit": m["unit"], "layer": m["layer"],
                     "component": m["component"], "subsample": m["subsample"], "tier": "extended_binary",
                     "value_type": "binary", "urinary_creatinine_adjusted": False,
                     "n_measured": int(bin_w[v].notna().sum()), "share_nonmissing": round(float(bin_w[v].notna().mean()), 4),
                     "share_below_lod": np.nan, "transform": "result_positive (0/1)", "included": not reason,
                     "exclusion_reason": reason, "n_population": n})
    dic = pd.DataFrame(rows)
    wide = wide.copy()
    wide["_log_urine_creatinine"] = np.log(wide.pop("_urine_creatinine"))
    return wide, bin_w, dic


def unit_duplicate_of(wide: pd.DataFrame, v: str, name: str, seen: dict) -> str:
    """An already-kept analyte with the same label whose values correlate > 0.999 (one analyte released in two
    units, e.g. URXUMA/URXUMS, LBDFOT/LBDFOTSI), else ''."""
    for other in seen.get(str(name).strip().lower(), []):
        both = wide[[v, other]].dropna()
        if both[v].nunique() < 3 or both[other].nunique() < 3:
            continue
        # rank correlation: the kept twin may already be log-transformed
        if len(both) > 50 and abs(both[v].corr(both[other], method="spearman")) > 0.999:
            return other
    return ""


# ======================================================================================================== effects
def ols_hc3(z: np.ndarray, X: np.ndarray, j: int = 1) -> tuple[float, float, float]:
    """(beta_j, HC3 se_j, two-sided normal p) — same estimator as statsmodels OLS(...).fit(cov_type='HC3')."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ z
    e = z - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
    w = (e / (1 - h)) ** 2
    cov = XtX_inv @ (X.T * w) @ X @ XtX_inv
    se = float(np.sqrt(cov[j, j]))
    p = float(2 * stats.norm.sf(abs(b[j] / se))) if se > 0 else np.nan
    return float(b[j]), se, p


def adjusted_effect(x, case, age, female, extra=None) -> dict:
    ok = np.isfinite(x) & np.isfinite(case) & np.isfinite(age) & np.isfinite(female)
    if extra is not None:
        extra = np.atleast_2d(extra.T).T if extra.ndim == 1 else extra
        ok &= np.isfinite(extra).all(1)
    x, case, age, female = x[ok], case[ok], age[ok], female[ok]
    n1, n0 = int(case.sum()), int((1 - case).sum())
    out = {"n_cases": n1, "n_controls": n0,
           "mean_cases": float(x[case == 1].mean()) if n1 else np.nan,
           "mean_controls": float(x[case == 0].mean()) if n0 else np.nan}
    sd = x.std(ddof=1) if len(x) > 1 else np.nan
    if n1 < MIN_CASES_EFFECT or n0 < MIN_CASES_EFFECT or not np.isfinite(sd) or sd == 0:
        return {**out, "effect_size": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan}
    z = (x - x.mean()) / sd
    ac = age - 50.0
    cols = [np.ones_like(z), case, ac, ac ** 2, female]
    if extra is not None:
        cols += [extra[ok][:, k] for k in range(extra.shape[1])]
    b, se, p = ols_hc3(z, np.column_stack(cols))
    return {**out, "effect_size": b, "se": se, "ci_low": b - 1.959964 * se, "ci_high": b + 1.959964 * se,
            "p_value": p}


def adjusted_logodds(x, case, age, female) -> dict:
    import statsmodels.api as sm
    ok = np.isfinite(x) & np.isfinite(case) & np.isfinite(age) & np.isfinite(female)
    x, case, age, female = x[ok], case[ok], age[ok], female[ok]
    n1, n0 = int(case.sum()), int((1 - case).sum())
    out = {"n_cases": n1, "n_controls": n0, "mean_cases": float(x[case == 1].mean()) if n1 else np.nan,
           "mean_controls": float(x[case == 0].mean()) if n0 else np.nan,
           "effect_size": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan}
    if (n1 < MIN_CASES_EFFECT or n0 < MIN_CASES_EFFECT or x[case == 1].sum() + x[case == 0].sum() < 2 * MIN_POS_BINARY
            or x.sum() < MIN_POS_BINARY or (len(x) - x.sum()) < MIN_POS_BINARY):
        return out
    ac = age - 50.0
    X = np.column_stack([np.ones_like(x), case, ac, ac ** 2, female])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = sm.Logit(x, X).fit(disp=0, maxiter=200)
        b, se = float(fit.params[1]), float(fit.bse[1])
        if not (np.isfinite(b) and np.isfinite(se)) or abs(b) > 15:
            return out
        return {**out, "effect_size": b, "se": se, "ci_low": b - 1.959964 * se, "ci_high": b + 1.959964 * se,
                "p_value": float(fit.pvalues[1])}
    except Exception:  # noqa: BLE001 - separation / singular fits are reported as not estimable
        return out


def bh_fdr(p: np.ndarray) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def _effects_for_label(lab: Label, frame: pd.DataFrame, wide: pd.DataFrame, bin_w: pd.DataFrame,
                       dic: pd.DataFrame) -> list[dict]:
    y = frame[f"y_{lab.label_id}"].to_numpy(float)
    sel = np.isfinite(y)
    d = frame[sel]
    case = y[sel]
    age, female, bmi = (d["age_years"].to_numpy(float), d["female"].to_numpy(float), d["bmi"].to_numpy(float))
    W = wide.loc[d["participant_id"]]
    B = bin_w.loc[d["participant_id"]]
    lucr = W["_log_urine_creatinine"].to_numpy(float)
    rows = []
    for _, a in dic[dic["included"]].iterrows():
        v = a["analyte"]
        if v in lab.drop_analytes:
            continue
        if a["value_type"] == "binary":
            r = adjusted_logodds(B[v].to_numpy(float), case, age, female)
            rows.append({"label_id": lab.label_id, "analyte": v, "tier": a["tier"],
                         "effect_measure": "age_sex_adjusted_log_odds_ratio", **r})
            continue
        x = W[v].to_numpy(float)
        extra = lucr if a["urinary_creatinine_adjusted"] else None
        r = adjusted_effect(x, case, age, female, extra=extra)
        extra_b = bmi if extra is None else np.column_stack([bmi, lucr])
        rb = adjusted_effect(x, case, age, female, extra=extra_b)
        rows.append({"label_id": lab.label_id, "analyte": v, "tier": a["tier"],
                     "effect_measure": "age_sex_adjusted_standardized_difference", **r,
                     "effect_size_bmi_adj": rb["effect_size"], "ci_low_bmi_adj": rb["ci_low"],
                     "ci_high_bmi_adj": rb["ci_high"], "p_value_bmi_adj": rb["p_value"]})
    return rows


def stage_effects(frame, wide, bin_w, dic, n_jobs: int = 16, log=print) -> pd.DataFrame:
    from joblib import Parallel, delayed
    t0 = time.time()
    res = Parallel(n_jobs=n_jobs)(delayed(_effects_for_label)(lab, frame, wide, bin_w, dic) for lab in LABELS)
    eff = pd.DataFrame([r for block in res for r in block])
    eff["q_value_bh"] = np.nan
    for (_, _), g in eff.groupby(["label_id", "tier"]):
        eff.loc[g.index, "q_value_bh"] = bh_fdr(g["p_value"].to_numpy())
        eff.loc[g.index, "q_value_bh_bmi_adj"] = bh_fdr(g["p_value_bmi_adj"].to_numpy()) \
            if "p_value_bmi_adj" in g else np.nan
    eff["fdr_significant"] = eff["q_value_bh"] < FDR_ALPHA
    eff["direction"] = np.where(eff["effect_size"].isna(), "not_estimable",
                                np.where(eff["effect_size"] > 0, "higher_in_cases", "lower_in_cases"))
    eff["underpowered"] = eff["n_cases"] < UNDERPOWERED_CASES
    meta = dic.set_index("analyte")[["lab_name", "unit", "layer", "transform"]]
    eff = eff.join(meta, on="analyte")
    eff.to_csv(TABLES / f"{PREFIX}_effects.csv", index=False)
    log(f"effects: {len(eff)} rows in {time.time() - t0:.0f}s")
    return eff


# ======================================================================================================== models
def _demog_cols() -> list[str]:
    from ..wearables.nhanes_cohort import DEMOGRAPHIC_FEATURES
    return list(DEMOGRAPHIC_FEATURES)


def primary_panel(dic: pd.DataFrame, lab: Label) -> list[str]:
    p = dic[(dic["tier"] == "primary") & dic["included"]]["analyte"].tolist()
    return [a for a in p if a not in lab.drop_analytes]


def feature_matrix(frame: pd.DataFrame, wide: pd.DataFrame, labs: list[str], fset: str) -> np.ndarray:
    cols = []
    if fset.startswith("demographics"):
        cols.append(frame[_demog_cols()].to_numpy(float))
    if "bmi" in fset:
        cols.append(frame[["bmi"]].to_numpy(float))
    if fset.endswith("labs"):
        cols.append(wide.loc[frame["participant_id"], labs].to_numpy(float))
    return np.column_stack(cols)


def make_model(seed: int, C: float | None = None):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    if C is None:
        est = LogisticRegressionCV(Cs=list(CS), cv=StratifiedKFold(3, shuffle=True, random_state=seed),
                                   scoring="neg_log_loss", solver="newton-cholesky", max_iter=300)
    else:
        est = LogisticRegression(C=C, solver="newton-cholesky", max_iter=300)
    return Pipeline([("impute", SimpleImputer(strategy="median", add_indicator=True)),
                     ("scale", StandardScaler()), ("model", est)])


def cv_folds(y: np.ndarray, repeat: int, n_splits: int = N_SPLITS):
    from sklearn.model_selection import StratifiedKFold
    return list(StratifiedKFold(n_splits, shuffle=True, random_state=SEED + repeat).split(np.zeros(len(y)), y))


def oof_predict(X: np.ndarray, y: np.ndarray, folds, seed: int, C: float | None = None) -> tuple[np.ndarray, list]:
    from sklearn.exceptions import ConvergenceWarning
    pred = np.full(len(y), np.nan)
    cs = []
    for k, (tr, te) in enumerate(folds):
        m = make_model(seed + k, C)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.simplefilter("ignore", RuntimeWarning)
            m.fit(X[tr], y[tr])
        pred[te] = m.predict_proba(X[te])[:, 1]
        est = m.named_steps["model"]
        cs.append(float(np.ravel(getattr(est, "C_", [C]))[0]))
    return pred, cs


def _cv_task(label_id: str, fset: str, repeat: int, X: np.ndarray, y: np.ndarray) -> dict:
    pred, cs = oof_predict(X, y, cv_folds(y, repeat), SEED + 100 * repeat)
    return {"label_id": label_id, "fset": fset, "repeat": repeat, "pred": pred, "Cs": cs}


def _perm_task(label_id: str, i: int, X: np.ndarray, y: np.ndarray, C: float) -> dict:
    from ..wearables.models import auroc
    rng = np.random.default_rng(SEED + 7919 * (i + 1))
    yp = rng.permutation(y)
    pred, _ = oof_predict(X, yp, cv_folds(yp, 0), SEED, C=C)
    return {"label_id": label_id, "perm": i, "auroc": auroc(pred, yp)}


def _shuffle_task(label_id: str, s: int, Xd: np.ndarray, Xl: np.ndarray, y: np.ndarray, base_pred: np.ndarray
                  ) -> dict:
    from ..wearables.models import auroc
    rng = np.random.default_rng(SEED + 104729 * (s + 1))
    Xs = np.column_stack([Xd, Xl[rng.permutation(len(y))]])
    pred, _ = oof_predict(Xs, y, cv_folds(y, 0), SEED)
    return {"label_id": label_id, "shuffle": s, "delta": auroc(pred, y) - auroc(base_pred, y)}


def boot_metrics(preds: np.ndarray, y: np.ndarray, Wb: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
    """AUROC/AUPRC averaged over repeats (preds: R x n) for each bootstrap weight row."""
    from ..wearables.models import rank_metrics
    au, ap = [], []
    for r in range(preds.shape[0]):
        a, p = rank_metrics(preds[r], y, Wb)
        au.append(a)
        ap.append(p)
    return np.mean(au, 0), np.mean(ap, 0)


def _ci(a):
    a = a[np.isfinite(a)]
    return float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975))


def _boot_p(d: np.ndarray) -> float:
    d = d[np.isfinite(d)]
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return float(max(p, 1.0 / len(d)))


def model_labels(frame: pd.DataFrame) -> list[Label]:
    out = []
    for lab in LABELS:
        y = frame[f"y_{lab.label_id}"]
        if (y == 1).sum() >= MIN_CASES_MODEL and (y == 0).sum() >= MIN_CASES_MODEL:
            out.append(lab)
    return out


def stage_models(frame, wide, dic, n_jobs: int = 16, log=print) -> dict:
    from joblib import Parallel, delayed

    from ..wearables.models import bootstrap_weights
    t0 = time.time()
    labs_m = model_labels(frame)
    data = {}
    for lab in labs_m:
        d = frame[frame[f"y_{lab.label_id}"].notna()].reset_index(drop=True)
        y = d[f"y_{lab.label_id}"].to_numpy(int)
        panel = primary_panel(dic, lab)
        fsets = FEATURE_SETS if lab.group != "pairwise" else ("demographics", "labs", "demographics_labs")
        data[lab.label_id] = {"y": y, "panel": panel, "fsets": fsets,
                              "X": {f: feature_matrix(d, wide, panel, f) for f in fsets}}
    tasks = [(lid, f, r) for lid, v in data.items() for f in v["fsets"] for r in range(N_REPEATS)]
    log(f"models: {len(labs_m)} labels, {len(tasks)} CV tasks")
    res = Parallel(n_jobs=n_jobs)(delayed(_cv_task)(lid, f, r, data[lid]["X"][f], data[lid]["y"])
                                  for lid, f, r in tasks)
    preds, cs_used = {}, {}
    for r in res:
        preds.setdefault((r["label_id"], r["fset"]), {})[r["repeat"]] = r["pred"]
        cs_used.setdefault((r["label_id"], r["fset"]), []).extend(r["Cs"])
    log(f"  CV done {time.time() - t0:.0f}s")

    perm_tasks = []
    for lid, v in data.items():
        C = float(np.median(cs_used[(lid, "labs")]))
        perm_tasks += [(lid, i, v["X"]["labs"], v["y"], C) for i in range(N_PERM)]
    perm = pd.DataFrame(Parallel(n_jobs=n_jobs)(delayed(_perm_task)(*t) for t in perm_tasks))
    log(f"  permutations done {time.time() - t0:.0f}s")
    sh_tasks = [(lid, s, v["X"]["demographics"], v["X"]["labs"], v["y"], preds[(lid, "demographics")][0])
                for lid, v in data.items() for s in range(N_SHUFFLE)]
    shuf = pd.DataFrame(Parallel(n_jobs=n_jobs)(delayed(_shuffle_task)(*t) for t in sh_tasks))
    log(f"  shuffles done {time.time() - t0:.0f}s")

    from ..wearables.models import auroc
    perf_rows, inc_rows, perm_rows, sh_rows, oof_rows = [], [], [], [], []
    for lid, v in data.items():
        y = v["y"]
        Wb = bootstrap_weights(len(y), N_BOOT, SEED + zlib.crc32(lid.encode()) % 10_000)
        P = {f: np.vstack([preds[(lid, f)][r] for r in range(N_REPEATS)]) for f in v["fsets"]}
        boot, point = {}, {}
        for f, pr in P.items():
            point[f] = boot_metrics(pr, y, None)
            boot[f] = boot_metrics(pr, y, Wb)
            per_rep = [auroc(pr[r], y) for r in range(N_REPEATS)]
            perf_rows.append({"label_id": lid, "feature_set": f, "n": len(y), "n_cases": int(y.sum()),
                              "auroc": float(point[f][0][0]), "auroc_ci_low": _ci(boot[f][0])[0],
                              "auroc_ci_high": _ci(boot[f][0])[1], "auprc": float(point[f][1][0]),
                              "auprc_ci_low": _ci(boot[f][1])[0], "auprc_ci_high": _ci(boot[f][1])[1],
                              "auprc_no_skill": float(y.mean()), "auroc_repeat_min": min(per_rep),
                              "auroc_repeat_max": max(per_rep),
                              "median_C": float(np.median(cs_used[(lid, f)])),
                              "n_lab_features": len(v["panel"]) if f.endswith("labs") else 0})
        for comp, (a, b) in COMPARISONS.items():
            if a not in P or b not in P:
                continue
            d_pt = float(point[a][0][0] - point[b][0][0])
            d_bt = boot[a][0] - boot[b][0]
            dp_bt = boot[a][1] - boot[b][1]
            inc_rows.append({"label_id": lid, "comparison": comp, "metric": "auroc", "delta": d_pt,
                             "ci_low": _ci(d_bt)[0], "ci_high": _ci(d_bt)[1], "boot_p": _boot_p(d_bt)})
            inc_rows.append({"label_id": lid, "comparison": comp, "metric": "auprc",
                             "delta": float(point[a][1][0] - point[b][1][0]),
                             "ci_low": _ci(dp_bt)[0], "ci_high": _ci(dp_bt)[1], "boot_p": _boot_p(dp_bt)})
        obs = auroc(P["labs"][0], y)
        pn = perm.loc[perm["label_id"] == lid, "auroc"].to_numpy()
        perm_rows.append({"label_id": lid, "feature_set": "labs", "observed_auroc_repeat1": obs,
                          "n_permutations": len(pn), "fixed_C": float(np.median(cs_used[(lid, "labs")])),
                          "null_mean": float(pn.mean()), "null_p95": float(np.quantile(pn, 0.95)),
                          "null_max": float(pn.max()), "empirical_p": float((1 + (pn >= obs).sum()) / (1 + len(pn)))})
        obs_d = auroc(P["demographics_labs"][0], y) - auroc(P["demographics"][0], y)
        sd = shuf.loc[shuf["label_id"] == lid, "delta"].to_numpy()
        sh_rows.append({"label_id": lid, "observed_delta_repeat1": obs_d, "n_shuffles": len(sd),
                        "shuffled_delta_mean": float(sd.mean()), "shuffled_delta_max": float(sd.max()),
                        "empirical_p": float((1 + (sd >= obs_d).sum()) / (1 + len(sd)))})
    perf, inc = pd.DataFrame(perf_rows), pd.DataFrame(inc_rows)
    permdf, shdf = pd.DataFrame(perm_rows), pd.DataFrame(sh_rows)
    perf.to_csv(TABLES / f"{PREFIX}_model_performance.csv", index=False)
    inc.to_csv(TABLES / f"{PREFIX}_increments.csv", index=False)
    permdf.to_csv(TABLES / f"{PREFIX}_permutation_null.csv", index=False)
    shdf.to_csv(TABLES / f"{PREFIX}_lab_block_shuffle.csv", index=False)
    perm.to_csv(TABLES / f"{PREFIX}_permutation_draws.csv", index=False)
    coef = full_data_coefficients(frame, wide, dic, labs_m)
    coef.to_csv(TABLES / f"{PREFIX}_coefficients_descriptive.csv", index=False)
    log(f"models done in {time.time() - t0:.0f}s")
    return {"perf": perf, "inc": inc, "perm": permdf, "shuffle": shdf}


def full_data_coefficients(frame, wide, dic, labs_m) -> pd.DataFrame:
    """Descriptive only (not a pre-specified test): standardised coefficients of the demographics+labs model fitted
    on all rows of each label (C re-tuned by 5-fold CV), to show which analytes carry the model."""
    rows = []
    for lab in labs_m:
        d = frame[frame[f"y_{lab.label_id}"].notna()].reset_index(drop=True)
        y = d[f"y_{lab.label_id}"].to_numpy(int)
        panel = primary_panel(dic, lab)
        X = feature_matrix(d, wide, panel, "demographics_labs")
        m = make_model(SEED)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(X, y)
        names = _demog_cols() + panel
        coef = m.named_steps["model"].coef_.ravel()[:len(names)]
        for nm, c in zip(names, coef):
            rows.append({"label_id": lab.label_id, "feature": nm, "std_coefficient": float(c),
                         "is_lab": nm in panel, "C": float(np.ravel(m.named_steps["model"].C_)[0])})
    return pd.DataFrame(rows)


# ======================================================================================================== summary
def decision_table(frame: pd.DataFrame, eff: pd.DataFrame, models: dict | None) -> pd.DataFrame:
    rows = []
    perf = models["perf"].set_index(["label_id", "feature_set"]) if models else None
    inc = models["inc"] if models else None
    perm = models["perm"].set_index("label_id") if models else None
    sh = models["shuffle"].set_index("label_id") if models else None
    for lab in LABELS:
        y = frame[f"y_{lab.label_id}"]
        e = eff[(eff["label_id"] == lab.label_id)]
        ep = e[e["tier"] == "primary"]
        ex = e[e["tier"] != "primary"]
        top = ep[ep["fdr_significant"]].reindex(ep[ep["fdr_significant"]]["effect_size"].abs()
                                               .sort_values(ascending=False).index).head(5)
        r = {"label_id": lab.label_id, "group": lab.group, "title": lab.title, "is_proxy": lab.is_proxy,
             "cycles": lab.cycles, "n_cases": int((y == 1).sum()), "n_controls": int((y == 0).sum()),
             "n_primary_analytes_tested": int(ep["p_value"].notna().sum()),
             "n_primary_fdr": int(ep["fdr_significant"].sum()),
             "n_extended_tested": int(ex["p_value"].notna().sum()), "n_extended_fdr": int(ex["fdr_significant"].sum()),
             "top_primary_fdr": "; ".join(f"{a} ({n}) {s:+.2f}" for a, n, s in
                                          zip(top["analyte"], top["lab_name"], top["effect_size"]))}
        if perf is not None and (lab.label_id, "demographics") in perf.index:
            pdm, plb = perf.loc[(lab.label_id, "demographics")], perf.loc[(lab.label_id, "labs")]
            pdl = perf.loc[(lab.label_id, "demographics_labs")]
            i = inc[(inc.label_id == lab.label_id) & (inc.comparison == PRIMARY_COMPARISON) & (inc.metric == "auroc")].iloc[0]
            ib = inc[(inc.label_id == lab.label_id) & (inc.comparison == "demographics_bmi_labs_minus_demographics_bmi")
                     & (inc.metric == "auroc")]
            add = bool(i.ci_low > 0 and sh.loc[lab.label_id, "empirical_p"] <= 0.05)
            alone = bool(plb.auroc_ci_low > 0.5 and perm.loc[lab.label_id, "empirical_p"] <= 0.05)
            r.update({
                "auroc_demographics": pdm.auroc, "auroc_demographics_ci": f"{pdm.auroc_ci_low:.3f}-{pdm.auroc_ci_high:.3f}",
                "auroc_labs": plb.auroc, "auroc_labs_ci": f"{plb.auroc_ci_low:.3f}-{plb.auroc_ci_high:.3f}",
                "auroc_demographics_labs": pdl.auroc,
                "auroc_demographics_labs_ci": f"{pdl.auroc_ci_low:.3f}-{pdl.auroc_ci_high:.3f}",
                "delta_auroc": i.delta, "delta_ci_low": i.ci_low, "delta_ci_high": i.ci_high, "delta_boot_p": i.boot_p,
                "delta_auroc_bmi_baseline": ib.delta.iloc[0] if len(ib) else np.nan,
                "delta_bmi_ci": f"{ib.ci_low.iloc[0]:+.3f} to {ib.ci_high.iloc[0]:+.3f}" if len(ib) else "",
                "perm_p_labs": perm.loc[lab.label_id, "empirical_p"],
                "shuffle_p": sh.loc[lab.label_id, "empirical_p"],
                "labs_add_to_demographics": add, "labs_alone_discriminate": alone,
                "reading": ("NULL: labs do not add to demographics" if not add else
                            ("labs add, but NEGLIGIBLE (Delta < 0.02)" if i.delta < 0.02 else "labs add to demographics")),
            })
        else:
            r["reading"] = (f"no model (fewer than {MIN_CASES_MODEL} cases or controls); per-analyte contrasts only")
        rows.append(r)
    out = pd.DataFrame(rows)
    um = out[(out.group == "umbrella") & out["delta_boot_p"].notna()] if "delta_boot_p" in out else out.iloc[0:0]
    out["delta_boot_p_bh_umbrella"] = np.nan
    if len(um):
        out.loc[um.index, "delta_boot_p_bh_umbrella"] = bh_fdr(um["delta_boot_p"].to_numpy())
    out.to_csv(TABLES / f"{PREFIX}_summary.csv", index=False)
    return out


def positive_control_gate(summary: pd.DataFrame) -> dict:
    need = ["diabetes_no_glycaemic_labs", "weak_failing_kidneys", "gout"]
    s = summary.set_index("label_id")
    res = {k: bool(s.loc[k, "labs_add_to_demographics"]) if k in s.index and "labs_add_to_demographics" in s
           and pd.notna(s.loc[k, "labs_add_to_demographics"]) else False for k in need}
    return {"required": need, "passed": res, "gate_passed": all(res.values())}


def posthoc_cotinine_adjusted(frame, wide, dic, log=print) -> pd.DataFrame:
    """POST HOC (added after the lock, after seeing that serum cotinine / urinary NNAL top several symptom labels):
    primary-panel contrasts additionally adjusted for log serum cotinine (LBXCOT), tobacco analytes themselves
    omitted. Descriptive sensitivity only; not part of the pre-specified decision rules."""
    rows = []
    cot = wide["LBXCOT"]
    prim = dic[(dic["tier"] == "primary") & dic["included"] & (dic["layer"] != "tobacco_cotinine")]
    for lab in LABELS:
        y = frame[f"y_{lab.label_id}"].to_numpy(float)
        sel = np.isfinite(y)
        d = frame[sel]
        case = y[sel]
        age, female = d["age_years"].to_numpy(float), d["female"].to_numpy(float)
        W = wide.loc[d["participant_id"]]
        c = cot.loc[d["participant_id"]].to_numpy(float)
        lucr = W["_log_urine_creatinine"].to_numpy(float)
        block = []
        for _, a in prim.iterrows():
            v = a["analyte"]
            if v in lab.drop_analytes:
                continue
            extra = np.column_stack([c, lucr]) if a["urinary_creatinine_adjusted"] else c
            r = adjusted_effect(W[v].to_numpy(float), case, age, female, extra=extra)
            block.append({"label_id": lab.label_id, "analyte": v, "lab_name": a["lab_name"],
                          "effect_size_cotinine_adj": r["effect_size"], "ci_low": r["ci_low"], "ci_high": r["ci_high"],
                          "p_value": r["p_value"], "n_cases": r["n_cases"], "n_controls": r["n_controls"]})
        q = bh_fdr(np.array([b["p_value"] for b in block]))
        for b, qq in zip(block, q):
            b["q_value_bh"] = qq
        rows += block
    out = pd.DataFrame(rows)
    out["analysis"] = "POST HOC: age/sex + log serum cotinine adjusted (not pre-specified)"
    out.to_csv(TABLES / f"{PREFIX}_posthoc_cotinine_adjusted.csv", index=False)
    log(f"post hoc cotinine-adjusted contrasts: {len(out)} rows")
    return out


POSTHOC_MODEL_LABELS = ("mecfs_like_proxy", "fatigue_symptom", "told_sleep_disorder", "depression_phq9",
                        "pair_mecfs_like_proxy_vs_depression_non_proxy", "asthma", "rheumatoid_arthritis")


def posthoc_models(frame, wide, dic, n_jobs: int = 16, log=print) -> pd.DataFrame:
    """POST HOC (after seeing the locked results; not pre-specified, descriptive):
    (i) smoking: demographics + log serum cotinine [+ BMI] vs the same + primary-panel labs without tobacco analytes;
    (ii) the ME/CFS-like proxy against controls with none of the proxy's exclusion diagnoses (so labs cannot win by
        recognising the excluded diseases, e.g. diabetes via HbA1c)."""
    from joblib import Parallel, delayed

    from ..wearables.models import bootstrap_weights
    demog = _demog_cols()
    nt = [a for a in dic[(dic["tier"] == "primary") & dic["included"] & (dic["layer"] != "tobacco_cotinine")]["analyte"]]
    specs = []
    for lid in POSTHOC_MODEL_LABELS:
        d = frame[frame[f"y_{lid}"].notna()].reset_index(drop=True)
        specs.append((lid, "smoking_adjusted", d))
    d = frame[frame["y_mecfs_like_proxy"].notna()]
    d = d[(d["y_mecfs_like_proxy"] == 1) | (d["proxy_any_exclusion_dx"] == 0)].reset_index(drop=True)
    specs.append(("mecfs_like_proxy", "controls_without_exclusion_dx", d))
    sets = {
        "smoking_adjusted": {"base": demog + ["_cot"], "full": demog + ["_cot"] + nt,
                             "base_bmi": demog + ["bmi", "_cot"], "full_bmi": demog + ["bmi", "_cot"] + nt},
        "controls_without_exclusion_dx": {"base": demog, "full": demog + nt, "base_bmi": demog + ["bmi"],
                                          "full_bmi": demog + ["bmi"] + nt},
    }
    jobs, mats = [], {}
    for lid, kind, d in specs:
        W = wide.loc[d["participant_id"]].reset_index(drop=True)
        Z = pd.concat([d[demog + ["bmi"]].reset_index(drop=True), W[nt], W[["LBXCOT"]].rename(columns={"LBXCOT": "_cot"})],
                      axis=1)
        y = d[f"y_{lid}"].to_numpy(int)
        for fs, cols in sets[kind].items():
            mats[(lid, kind, fs)] = (Z[cols].to_numpy(float), y)
            jobs += [(lid, kind, fs, r) for r in range(N_REPEATS)]
    res = Parallel(n_jobs=n_jobs)(delayed(_cv_task)(f"{l}|{k}", fs, r, *mats[(l, k, fs)]) for l, k, fs, r in jobs)
    preds = {}
    for r in res:
        preds.setdefault((r["label_id"], r["fset"]), {})[r["repeat"]] = r["pred"]
    rows = []
    for lid, kind, _ in specs:
        y = mats[(lid, kind, "base")][1]
        Wb = bootstrap_weights(len(y), N_BOOT, SEED + zlib.crc32(f"{lid}|{kind}".encode()) % 10_000)
        P = {fs: np.vstack([preds[(f"{lid}|{kind}", fs)][r] for r in range(N_REPEATS)]) for fs in sets[kind]}
        pt = {fs: boot_metrics(P[fs], y, None)[0][0] for fs in P}
        bt = {fs: boot_metrics(P[fs], y, Wb)[0] for fs in P}
        for full, base in (("full", "base"), ("full_bmi", "base_bmi")):
            dd = bt[full] - bt[base]
            rows.append({"label_id": lid, "analysis": kind, "baseline": "+".join(sets[kind][base][len(demog):]) or
                         "demographics", "n": len(y), "n_cases": int(y.sum()), "auroc_base": pt[base],
                         "auroc_full": pt[full], "delta_auroc": pt[full] - pt[base], "ci_low": _ci(dd)[0],
                         "ci_high": _ci(dd)[1], "boot_p": _boot_p(dd),
                         "note": "POST HOC (not pre-specified); labs exclude tobacco analytes"})
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / f"{PREFIX}_posthoc_models.csv", index=False)
    log(out.to_string())
    return out


# ======================================================================================================== signatures
SIG_COLUMNS = ["object_id", "dataset_id", "phenotype_id", "phenotype_label", "condition_id", "is_proxy", "label_basis",
               "feature", "feature_label", "effect_measure", "effect_size", "ci_low", "ci_high", "null_value", "se",
               "p_value", "q_value_bh", "fdr_significant", "direction", "n_cases", "n_controls", "underpowered",
               "phenotype_definition", "population", "cycles", "feature_units", "feature_transform", "covariates",
               "mean_cases", "mean_controls", "caveats", "effect_unit", "q_value", "method"]
GENERIC_CAVEATS = ("NHANES 2011-2014, adults >= 20 with the biochemistry profile; unweighted internal comparison (not "
                   "a US prevalence or population estimate); self-reported labels; general-population controls; no "
                   "CRP, ferritin or ESR in these cycles. A per-analyte difference is a group difference, not "
                   "individual-level discrimination (see the model rows).")


def build_signatures(frame, eff, dic, models, source_meta: dict) -> pd.DataFrame:
    from ..provenance import add_provenance
    from ..store import write_table
    rows = []
    unit = dic.set_index("analyte")["unit"].to_dict()
    for _, e in eff.iterrows():
        lab = LABEL_BY_ID[e["label_id"]]
        cov = "age (centred at 50), age^2, sex" + ("; log urinary creatinine" if dic.set_index("analyte").loc[
            e["analyte"], "urinary_creatinine_adjusted"] else "")
        rows.append({
            "phenotype_id": lab.label_id, "phenotype_label": lab.title, "condition_id": lab.condition_id,
            "is_proxy": lab.is_proxy, "label_basis": lab.label_basis, "feature": e["analyte"],
            "feature_label": f"{e['lab_name']} [{e['layer']}; tier {e['tier']}]",
            "effect_measure": e["effect_measure"], "effect_size": e["effect_size"], "ci_low": e["ci_low"],
            "ci_high": e["ci_high"], "null_value": 0.0, "se": e["se"], "p_value": e["p_value"],
            "q_value_bh": e["q_value_bh"], "fdr_significant": bool(e["fdr_significant"]), "direction": e["direction"],
            "n_cases": e["n_cases"], "n_controls": e["n_controls"], "underpowered": bool(e["underpowered"]),
            "phenotype_definition": lab.definition, "cycles": lab.cycles, "feature_units": unit.get(e["analyte"], ""),
            "feature_transform": e["transform"], "covariates": cov, "mean_cases": e["mean_cases"],
            "mean_controls": e["mean_controls"],
            "caveats": " ".join(x for x in (lab.caveat, GENERIC_CAVEATS) if x),
            "effect_unit": ("SD units of the transformed analyte" if e["effect_measure"].endswith("difference")
                            else "log odds ratio"),
            "q_value": e["q_value_bh"],
            "method": (f"OLS (HC3) of z-scored analyte on case + covariates; BH-FDR within label x tier "
                       f"({e['tier']}); plan docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md"),
            "evidence_level": lab.label_basis,
        })
    if models:
        perf = models["perf"].set_index(["label_id", "feature_set"])
        inc = models["inc"]
        perm = models["perm"].set_index("label_id")
        for (lid, f), p in perf.iterrows():
            lab = LABEL_BY_ID[lid]
            rows.append({
                "phenotype_id": lid, "phenotype_label": lab.title, "condition_id": lab.condition_id,
                "is_proxy": lab.is_proxy, "label_basis": lab.label_basis, "feature": f"model_auroc:{f}",
                "feature_label": f"cross-validated AUROC, feature set {f} (penalised logistic, 5x5 CV)",
                "effect_measure": "cv_auroc", "effect_size": p.auroc, "ci_low": p.auroc_ci_low,
                "ci_high": p.auroc_ci_high, "null_value": 0.5, "se": np.nan,
                "p_value": perm.loc[lid, "empirical_p"] if f == "labs" else np.nan, "q_value_bh": np.nan,
                "fdr_significant": False, "direction": "", "n_cases": p.n_cases, "n_controls": p.n - p.n_cases,
                "underpowered": p.n_cases < UNDERPOWERED_CASES, "phenotype_definition": lab.definition,
                "cycles": lab.cycles, "feature_units": "AUROC", "feature_transform": "", "covariates": "",
                "mean_cases": np.nan, "mean_controls": np.nan, "caveats": " ".join(x for x in (lab.caveat, GENERIC_CAVEATS) if x),
                "effect_unit": "AUROC", "q_value": np.nan,
                "method": "participant bootstrap CI (1,000); p (labs only) = label-permutation null (200)",
                "evidence_level": lab.label_basis})
        for _, i in inc[inc.metric == "auroc"].iterrows():
            lab = LABEL_BY_ID[i.label_id]
            p = perf.loc[(i.label_id, "demographics")]
            rows.append({
                "phenotype_id": i.label_id, "phenotype_label": lab.title, "condition_id": lab.condition_id,
                "is_proxy": lab.is_proxy, "label_basis": lab.label_basis, "feature": f"model_delta_auroc:{i.comparison}",
                "feature_label": f"Delta AUROC {i.comparison}", "effect_measure": "cv_delta_auroc",
                "effect_size": i.delta, "ci_low": i.ci_low, "ci_high": i.ci_high, "null_value": 0.0, "se": np.nan,
                "p_value": i.boot_p, "q_value_bh": np.nan, "fdr_significant": False,
                "direction": "labs_add" if i.ci_low > 0 else ("labs_worse" if i.ci_high < 0 else "no_detectable_increment"),
                "n_cases": p.n_cases, "n_controls": p.n - p.n_cases, "underpowered": p.n_cases < UNDERPOWERED_CASES,
                "phenotype_definition": lab.definition, "cycles": lab.cycles, "feature_units": "AUROC",
                "feature_transform": "", "covariates": "", "mean_cases": np.nan, "mean_controls": np.nan,
                "caveats": " ".join(x for x in (lab.caveat, GENERIC_CAVEATS) if x), "effect_unit": "AUROC difference",
                "q_value": np.nan, "method": "paired participant bootstrap (1,000) on repeat-averaged CV predictions",
                "evidence_level": lab.label_basis})
    sig = pd.DataFrame(rows)
    sig["dataset_id"] = DATASET_ID
    sig["population"] = "NHANES 2011-2014 adults >= 20, not pregnant, biochemistry profile measured"
    sig["object_id"] = "signature:" + DATASET_ID + "|" + sig["phenotype_id"] + "|" + sig["feature"]
    sig["condition_id"] = sig["condition_id"].astype(object).where(sig["condition_id"].notna(), None)
    ev = sig.pop("evidence_level")
    sig = sig[SIG_COLUMNS]
    sig = add_provenance(
        sig, data_layer="person", source_name=source_meta["source_name"],
        source_version=source_meta["source_version"], retrieved_at=source_meta["retrieved_at"],
        evidence_type="person_lab_measurement", source_record_id="object_id", evidence_level=ev,
        provenance_notes=("Derived per-label contrasts of NHANES laboratory analytes (participant_labs__nhanes, "
                          "participant_labs_extended__nhanes, KIQ_U) joined on SEQN; plan sha256 "
                          + source_meta["plan_sha256"][:16] + "."))
    assert sig["object_id"].is_unique
    write_table(sig, "phenotype_signatures__nhanes_labs", producer=PRODUCER,
                description="NHANES 2011-2014 lab analytes vs diagnosis-like labels: adjusted standardised "
                            "differences + cross-validated model rows (docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md)")
    return sig


def _source_meta(frame: pd.DataFrame, plan_sha: str) -> dict:
    vers = "; ".join(sorted(set(frame["source_version"].astype(str))))
    return {"source_name": str(frame["source_name"].iloc[0]),
            "source_version": vers + "; KIQ_U_G/KIQ_U_H (kidney questionnaire) retrieved 2026-09-24",
            "retrieved_at": max(frame["retrieved_at"].astype(str)), "plan_sha256": plan_sha}


# ======================================================================================================== run
def run(stages=("effects", "models", "signatures", "posthoc", "posthoc_models"), n_jobs: int = 16, log=print) -> dict:
    plan_sha = check_plan_locked()
    t0 = time.time()
    frame, flow = load_population()
    wide, bin_w, dic = load_analytes(frame["participant_id"])
    TABLES.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(flow).to_csv(TABLES / f"{PREFIX}_population_flow.csv", index=False)
    dic.to_csv(TABLES / f"{PREFIX}_analytes.csv", index=False)
    lab_rows = []
    for lab in LABELS:
        y = frame[f"y_{lab.label_id}"]
        lab_rows.append({"label_id": lab.label_id, "group": lab.group, "title": lab.title,
                         "definition": lab.definition, "label_basis": lab.label_basis, "is_proxy": lab.is_proxy,
                         "cycles": lab.cycles, "n_cases": int((y == 1).sum()), "n_controls": int((y == 0).sum()),
                         "dropped_analytes": ";".join(sorted(lab.drop_analytes)), "caveat": lab.caveat})
    pd.DataFrame(lab_rows).to_csv(TABLES / f"{PREFIX}_labels.csv", index=False)
    log(f"population {len(frame)}; analytes included {int(dic.included.sum())} "
        f"(primary {int(((dic.tier == 'primary') & dic.included).sum())})")
    eff = stage_effects(frame, wide, bin_w, dic, n_jobs=n_jobs, log=log) if "effects" in stages else \
        pd.read_csv(TABLES / f"{PREFIX}_effects.csv")
    if "models" in stages:
        models = stage_models(frame, wide, dic, n_jobs=n_jobs, log=log)
    else:
        try:
            models = {"perf": pd.read_csv(TABLES / f"{PREFIX}_model_performance.csv"),
                      "inc": pd.read_csv(TABLES / f"{PREFIX}_increments.csv"),
                      "perm": pd.read_csv(TABLES / f"{PREFIX}_permutation_null.csv"),
                      "shuffle": pd.read_csv(TABLES / f"{PREFIX}_lab_block_shuffle.csv")}
        except FileNotFoundError:
            models = None
    summary = decision_table(frame, eff, models)
    gate = positive_control_gate(summary) if models else {}
    meta = {"computed_at": utc_now_iso(), "plan_sha256": plan_sha, "seed": SEED, "n_population": len(frame),
            "n_splits": N_SPLITS, "n_repeats": N_REPEATS, "n_boot": N_BOOT, "n_perm": N_PERM,
            "n_shuffle": N_SHUFFLE, "Cs": CS, "positive_control_gate": gate, "seconds": round(time.time() - t0)}
    if "posthoc" in stages:
        posthoc_cotinine_adjusted(frame, wide, dic, log=log)
    if "posthoc_models" in stages:
        posthoc_models(frame, wide, dic, n_jobs=n_jobs, log=log)
    if "signatures" in stages:
        sig = build_signatures(frame, eff, dic, models, _source_meta(frame, plan_sha))
        meta["n_signature_rows"] = len(sig)
    (TABLES / f"{PREFIX}_run_metadata.json").write_text(json.dumps(meta, indent=2, default=str))
    log(json.dumps(gate))
    return meta


def _cli() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stages", nargs="+", default=["effects", "models", "signatures", "posthoc", "posthoc_models"])
    ap.add_argument("--jobs", type=int, default=16)
    a = ap.parse_args()
    run(tuple(a.stages), n_jobs=a.jobs)


if __name__ == "__main__":
    _cli()
