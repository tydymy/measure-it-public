"""MY-LC (Mount Sinai-Yale Long COVID): serum cortisol and 144 plasma mediators vs Long COVID status.

Klein J, Wood J, Jaycox JR, et al. "Distinguishing features of long COVID identified through immune profiling."
Nature 2023;623:139-148. doi:10.1038/s41586-023-06651-y (PMC10620090). Supplementary Table 3 (MOESM4) is the
authors' integrated machine-learning table: one row per participant (x0_LC_ID) with demographics, draw time,
exclusion flags, symptom surveys and 7,259 x1_* immune features (144 x1_Cytokines analytes incl. cortisol).

Cohort codes (x0_Censor_Cohort_ID): 1 = healthy uninfected (HC), 2 = convalescent without persistent symptoms (CC),
3 = Long COVID (LC; persistent symptoms > 6 weeks after WHO-confirmed/probable COVID-19, recruited from Mount Sinai LC
clinics). x0_Censor_Complete == 0 leaves 40 HC, 39 CC, 99 LC, the paper's analysed groups.

Pre-specified primary analysis: docs/ANALYSIS_PLAN_KLEIN2023_MYLC_ML_TABLE.md.

Run: uv run python -m measure_it.labs.klein2023_mylc_ml_table [--jobs 16] [--n-perm 1000]
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from measure_it.config import FIGURES, PROJECT_ROOT, SEED, TABLES, raw_dir, utc_now_iso
from measure_it.download import download_file, load_manifest
from measure_it.labs import lab_cc_engine as E
from measure_it.provenance import add_provenance
from measure_it.registry import write_registry_entry
from measure_it.store import write_table

SOURCE_ID = "klein2023_mylc_ml_table"
DATA_URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06651-y/MediaObjects/"
            "41586_2023_6651_MOESM4_ESM.xlsx")
DATA_FILE = "41586_2023_6651_MOESM4_ESM.xlsx"
LANDING_URL = "https://doi.org/10.1038/s41586-023-06651-y"
SOURCE_NAME = "Klein et al. 2023 Nature (MY-LC) Supplementary Table 3: integrated per-participant feature table"
CITATION = ("Klein J, Wood J, Jaycox JR, et al. Distinguishing features of long COVID identified through immune "
            "profiling. Nature 2023;623:139-148. doi:10.1038/s41586-023-06651-y")
CONDITION_ID = "long_covid"
MEASUREMENT_CLASS_ID = "blood_biomarkers"
PLAN = "docs/ANALYSIS_PLAN_KLEIN2023_MYLC_ML_TABLE.md"
REPORT = "results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md"
ANALYSIS_VERSION = f"{SOURCE_ID} 2026-09-24.1"
BUA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"}
CORT = "x1_Cytokines_Cortisol_Obs_Conc_ng/mL_GG21.1_ML"
COHORT = {1: "HC", 2: "CC", 3: "LC"}
COVS = ["draw_time_min", "age_years", "sex_code", "bmi"]
LABEL_BASIS = ("MY-LC study groups: LC = age >= 18, previous confirmed or probable COVID-19 (WHO) and persistent "
               "symptoms > 6 weeks after infection, recruited from Mount Sinai LC clinics; HC = no previous SARS-CoV-2 "
               "infection and no active symptoms on screening; CC = previous infection without persistent symptoms. "
               "Coded x0_Censor_Cohort_ID 3/1/2 (mapping checked against the infection-test columns).")
CAVEAT = ("Case-control by recruitment channel (LC clinics vs advertisement); one blood draw per person; cortisol is "
          "diurnal (draw time adjusted linearly) and measured on a multiplex immunoassay; values and labels are the "
          "authors' released rows, so this re-analyses the same data, it is not an independent replication.")


def fetch():
    return download_file(DATA_URL, SOURCE_ID, DATA_FILE, headers=BUA, max_retries=3)


def read_raw() -> pd.DataFrame:
    d = pd.read_excel(raw_dir(SOURCE_ID) / DATA_FILE, sheet_name="Sheet1")
    if d.shape != (185, 7431) or CORT not in d.columns:
        raise ValueError(f"unexpected table shape {d.shape}")
    return d


def analytes(d: pd.DataFrame) -> list[str]:
    return [c for c in d.columns if c.startswith("x1_Cytokines_")]


def _demojibake(t: str) -> str:
    try:
        return t.encode("latin-1").decode("utf-8")
    except (UnicodeEncodeError, UnicodeDecodeError):
        return t


def analyte_name(col: str) -> str:
    return _demojibake(col.replace("x1_Cytokines_", "").split("_Obs_Conc_")[0])


def analyte_unit(col: str) -> str:
    return _demojibake(col.split("_Obs_Conc_")[1].rsplit("_GG", 1)[0]) if "_Obs_Conc_" in col else ""


def participant_id(nid: str) -> str:
    return f"{SOURCE_ID}:{nid}"


def build_frames(d: pd.DataFrame, retrieved_at: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    nid = d["x0_LC_ID"].astype(str)
    grp = d["x0_Censor_Cohort_ID"].map(COHORT)
    part = pd.DataFrame({
        "participant_id": nid.map(participant_id), "dataset_id": SOURCE_ID, "native_id": nid, "group_label": grp,
        "long_covid": (grp == "LC").astype(int), "condition_id": np.where(grp == "LC", CONDITION_ID, ""),
        "label_basis": LABEL_BASIS, "in_analysis_set": d["x0_Censor_Complete"].eq(0).values,
        "censor_complete_code": d["x0_Censor_Complete"].values,
        "oral_steroid_flag": d["x0_Censor_Oral_Steroid"].values, "pit_adrenal_flag": d["x0_Censor_Pit_Adre_Dysfunction"].values,
        "draw_time_min": pd.to_numeric(d["x0_Sample_Time_Min"], errors="coerce").values,
        "age_years": pd.to_numeric(d["x0_Demographics_Age"], errors="coerce").values,
        "sex_code": pd.to_numeric(d["x0_Demographics_Sex"], errors="coerce").values,
        "bmi": pd.to_numeric(d["x0_Demographics_BMI"], errors="coerce").values,
        "hospitalised_acute": pd.to_numeric(d["x0_Acute_COVID19_Hosp"], errors="coerce").values,
        "days_from_acute_onset": pd.to_numeric(d["x0_Acute_COVID19_DFSO"], errors="coerce").values,
        "lc_symptom_total": pd.to_numeric(d["x0_LC_Symptom_totalsympt"], errors="coerce").values,
    })
    part = add_provenance(part, data_layer="person", source_name=SOURCE_NAME, source_version=f"{DATA_FILE} (Sheet1)",
                          retrieved_at=retrieved_at, evidence_type="person_lab_measurement", source_record_id="native_id",
                          provenance_notes="one row per MY-LC participant; sex_code as released (1/2, meaning not stated)")
    cols = analytes(d)
    long = d[["x0_LC_ID"] + cols].melt(id_vars="x0_LC_ID", var_name="lab_variable", value_name="value")
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    long = long.assign(participant_id=long.x0_LC_ID.astype(str).map(participant_id), dataset_id=SOURCE_ID,
                       lab_name=long.lab_variable.map(analyte_name),
                       unit="z-score (batch-integrated; column name states " + long.lab_variable.map(analyte_unit) + ")",
                       specimen="plasma/serum (multiplex immunoassay)", harmonized_name=long.lab_variable.map(analyte_name)
                       .str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True),
                       lab_note=np.where(long.value.isna(), "missing", ""))
    long["source_record_id"] = long.x0_LC_ID.astype(str) + ":" + long.lab_variable
    long = long.drop(columns=["x0_LC_ID"])
    long = add_provenance(long, data_layer="person", source_name=SOURCE_NAME, source_version=DATA_FILE,
                          retrieved_at=retrieved_at, evidence_type="person_lab_measurement",
                          source_record_id="source_record_id", provenance_notes="x1_Cytokines_* values as released")
    return part, long


def ingest() -> dict:
    fetch()
    ret = load_manifest(SOURCE_ID)["files"][DATA_FILE]["retrieved_at"]
    d = read_raw()
    part, labs = build_frames(d, ret)
    write_table(part, f"participants__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="MY-LC participants: group, analysis-set flag, draw time, demographics")
    write_table(labs, f"participant_labs__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="MY-LC 144 plasma analytes incl. cortisol, long format")
    return {"raw": d, "part": part, "retrieved_at": ret}


def _floor(v) -> float:
    """Half the smallest positive value (log offset for zeros); 1e-9 if none is positive."""
    v = pd.to_numeric(pd.Series(v), errors="coerce")
    pos = v[v > 0]
    return float(pos.min()) / 2 if len(pos) else 1e-9


def log10p(x: np.ndarray, floor: float) -> np.ndarray:
    return np.log10(np.asarray(x, float) + floor)


def analyse(dd: dict, jobs: int = 8, n_perm: int = 1000, log=print) -> dict:
    d, part = dd["raw"], dd["part"]
    cols = analytes(d)
    df = pd.concat([part[["native_id", "group_label", "long_covid", "in_analysis_set", "oral_steroid_flag",
                          "pit_adrenal_flag"] + COVS].reset_index(drop=True), d[cols].reset_index(drop=True)], axis=1)
    # Deviation 1 (plan): the released x1_Cytokines_* values are batch-integrated z-scores (negative values
    # included), so they enter untransformed; the column keeps its planned name for traceability.
    df["log_cortisol"] = pd.to_numeric(df[CORT], errors="coerce")
    a = df[df.in_analysis_set].copy()
    prim = a.dropna(subset=["log_cortisol"] + COVS).reset_index(drop=True)
    y = prim.long_covid.to_numpy(int)
    log(f"primary sample: {int(y.sum())} LC vs {int((1 - y).sum())} controls "
        f"(HC {int((prim.group_label == 'HC').sum())}, CC {int((prim.group_label == 'CC').sum())})")
    bidx = E.stratified_boot_idx(y)
    m1 = E.cv_block(prim[COVS].to_numpy(float), y, 1.0, bidx)
    m2 = E.cv_block(prim[["log_cortisol"]].to_numpy(float), y, 1.0, bidx, n_perm=n_perm, jobs=jobs)
    m3 = E.cv_block(prim[COVS + ["log_cortisol"]].to_numpy(float), y, 1.0, bidx)
    dl = E.delta(m3, m1, y)
    rows = [{"model": "1 covariates (draw time, age, sex, BMI)", **m1["summary"]},
            {"model": "2 cortisol alone (released z-score)", **m2["summary"]},
            {"model": "3 covariates + cortisol", **m3["summary"]}]
    # S3 multi-analyte panel on the analysis set (analytes only, median-imputed in fold)
    pan = a.reset_index(drop=True)
    Xp = np.column_stack([pd.to_numeric(pan[c], errors="coerce").to_numpy(float) for c in cols])
    yp = pan.long_covid.to_numpy(int)
    s3 = E.cv_block(Xp, yp, 0.1, E.stratified_boot_idx(yp), n_perm=n_perm, jobs=jobs)
    rows.append({"model": "S3 144-analyte panel (released z-scores, L2 C=0.1)", **s3["summary"]})
    # S4 exclusions
    s4d = prim[(prim.oral_steroid_flag != 1) & (prim.pit_adrenal_flag != 1)].reset_index(drop=True)
    y4 = s4d.long_covid.to_numpy(int)
    b4 = E.stratified_boot_idx(y4)
    s4a = E.cv_block(s4d[COVS].to_numpy(float), y4, 1.0, b4)
    s4b = E.cv_block(s4d[COVS + ["log_cortisol"]].to_numpy(float), y4, 1.0, b4)
    dl4 = E.delta(s4b, s4a, y4)
    rows.append({"model": "S4 covariates + cortisol, steroid / pituitary-adrenal flags excluded", **s4b["summary"]})
    perf = pd.DataFrame(rows)
    perf["cv"] = f"stratified {E.N_SPLITS}-fold x {E.N_REPEATS} (seed {SEED}+r); L2 logistic, standardised in fold"
    deltas = pd.DataFrame([{"comparison": "PRIMARY model 3 minus model 1", **dl, "n_cases": int(y.sum()),
                            "n_controls": int((1 - y).sum())},
                           {"comparison": "S4 (exclusions) covariates+cortisol minus covariates", **dl4,
                            "n_cases": int(y4.sum()), "n_controls": int((1 - y4).sum())}])
    # S2 + S6 fit-free cortisol (lower = LC)
    ff = []
    for name, ctrl in (("S6 LC vs HC+CC (primary sample)", ["HC", "CC"]), ("S2 LC vs HC", ["HC"]), ("S2 LC vs CC", ["CC"])):
        s = a[a.group_label.isin(["LC"] + ctrl) & a.log_cortisol.notna()]
        r = E.fitfree_test(s.log_cortisol.to_numpy(), s.long_covid.to_numpy(), n_perm=n_perm, higher_is_case=False)
        ff.append({"contrast": name, **r})
    fitfree = pd.DataFrame(ff)
    # S1 per-analyte effects
    eff_df = a.copy()
    for c in cols:
        v = pd.to_numeric(eff_df[c], errors="coerce")
        eff_df[c] = v
    eff = E.per_analyte_effects(eff_df, "long_covid", cols)     # released z-scores (Deviation 1)
    eff["analyte_name"] = eff.analyte.map(analyte_name)
    eff["unit"] = eff.analyte.map(analyte_unit)
    # E1 (POST HOC, exploratory; added after the primary was seen): covariates alone reach AUROC ~0.70, so check the
    # fit-free cortisol AUROC within draw-time tertiles (strata-preserving bootstrap and permutation).
    tert = pd.qcut(prim.draw_time_min, 3, labels=["early", "middle", "late"]).astype(str).to_numpy()
    e1 = E.fitfree_test(prim.log_cortisol.to_numpy(), y, n_perm=n_perm, higher_is_case=False, strata=tert)
    tcomp = pd.crosstab(tert, prim.group_label).reset_index().rename(columns={"row_0": "draw_time_tertile"})
    tcomp["tertile_minutes"] = tcomp.iloc[:, 0].map(
        prim.groupby(tert).draw_time_min.agg(lambda s: f"{int(s.min())}-{int(s.max())}"))
    fitfree = pd.concat([fitfree, pd.DataFrame([{"contrast": "E1 POST HOC exploratory: LC vs HC+CC within draw-time "
                                                             "tertiles (pair-weighted)", **e1}])], ignore_index=True)
    ok = bool(dl["ci_low"] > 0 and m2["summary"]["perm_p"] < 0.05)
    decision = pd.DataFrame([{
        "question": "Does serum cortisol add information about Long COVID beyond draw time, age, sex and BMI?",
        "delta_auroc": dl["delta_auroc"], "delta_ci_low": dl["ci_low"], "delta_ci_high": dl["ci_high"],
        "cortisol_only_auroc": m2["summary"]["auroc"], "cortisol_only_ci_low": m2["summary"]["auroc_ci_low"],
        "cortisol_only_ci_high": m2["summary"]["auroc_ci_high"], "cortisol_only_perm_p": m2["summary"]["perm_p"],
        "n_perm": n_perm, "null_q95": m2["summary"]["null_q95"], "n_cases": int(y.sum()), "n_controls": int((1 - y).sum()),
        "rule": "adds information only if Delta CI low > 0 AND cortisol-only permutation p < 0.05",
        "decision": ("cortisol adds information beyond draw time and demographics (pre-specified rule met)" if ok
                     else "not shown (pre-specified rule not met)"),
        "caveat": CAVEAT, "plan": PLAN, "analysis_version": ANALYSIS_VERSION, "computed_at": utc_now_iso()}])
    log(decision[["delta_auroc", "delta_ci_low", "delta_ci_high", "cortisol_only_auroc", "cortisol_only_perm_p",
                  "decision"]].to_string())
    return {"perf": perf, "deltas": deltas, "fitfree": fitfree, "effects": eff, "decision": decision, "tertiles": tcomp,
            "null": m2["null"], "panel_null": s3["null"], "y": y}


def _level(p, q, lo, hi, null):
    excl = (lo > null) or (hi < null)
    if np.isfinite(q) and q < 0.05 and excl:
        return "supported_single_dataset"
    if np.isfinite(p) and p < 0.05:
        return "nominal_single_dataset"
    return "null_single_dataset"


def signatures(res: dict, retrieved_at: str) -> pd.DataFrame:
    rows = []
    for r in res["effects"].dropna(subset=["p_value"]).itertuples():
        rows.append({"feature": "plasma_" + r.analyte_name.lower().replace(" ", "_"),
                     "feature_label": f"{r.analyte_name} (released z-score)", "effect_measure": "hedges_g",
                     "effect_size": r.hedges_g, "effect_unit": "Hedges' g on released z-scores (LC minus HC+CC)",
                     "ci_low": r.g_ci_low, "ci_high": r.g_ci_high, "p_value": r.p_value, "q_value": r.q_value,
                     "n_cases": r.n_cases, "n_controls": r.n_controls, "null_value": 0.0,
                     "method": "S1: stratified bootstrap CI (2,000); two-sided Mann-Whitney p; BH q over 144 analytes (z-scores, Deviation 1)"})
    dec = res["decision"].iloc[0]
    m2 = res["perf"].iloc[1]
    rows.append({"feature": "model_auroc:log_cortisol:l2_logistic", "feature_label": "serum cortisol alone (CV)",
                 "effect_measure": "auroc", "effect_size": m2.auroc, "effect_unit": "cross-validated AUROC",
                 "ci_low": m2.auroc_ci_low, "ci_high": m2.auroc_ci_high, "p_value": m2.perm_p, "q_value": np.nan,
                 "n_cases": m2.n_cases, "n_controls": m2.n_controls, "null_value": 0.5,
                 "method": "primary model 2: 5-fold x 20 CV; person bootstrap CI; label-permutation p"})
    rows.append({"feature": "model_delta_auroc:cortisol_over_drawtime_age_sex_bmi:l2_logistic",
                 "feature_label": "cortisol added to draw time + age + sex + BMI", "effect_measure": "delta_auroc",
                 "effect_size": dec.delta_auroc, "effect_unit": "AUROC(model 3) - AUROC(model 1)",
                 "ci_low": dec.delta_ci_low, "ci_high": dec.delta_ci_high, "p_value": np.nan, "q_value": np.nan,
                 "n_cases": dec.n_cases, "n_controls": dec.n_controls, "null_value": 0.0,
                 "method": f"PRIMARY (locked): paired person bootstrap CI (2,000). Decision: {dec.decision}"})
    sig = pd.DataFrame(rows)
    sig.insert(0, "object_id", "signature:" + SOURCE_ID + "|long_covid|" + sig.feature)
    sig.insert(1, "dataset_id", SOURCE_ID)
    sig.insert(2, "phenotype_id", "long_covid")
    sig.insert(3, "phenotype_label", "Long COVID (MY-LC clinic cohort) vs healthy + convalescent controls")
    sig.insert(4, "condition_id", CONDITION_ID)
    sig.insert(5, "is_proxy", False)
    sig["label_basis"] = LABEL_BASIS
    sig["phenotype_definition"] = "MY-LC analysed set: 99 LC vs 40 HC + 39 CC; one blood draw per person"
    sig["caveats"] = CAVEAT
    sig["measurement_class_id"] = MEASUREMENT_CLASS_ID
    sig = sig.drop_duplicates("object_id")
    sig["n_cases"] = sig.n_cases.astype("Int64")
    sig["n_controls"] = sig.n_controls.astype("Int64")
    lv = []
    for r in sig.itertuples():
        if r.effect_measure == "delta_auroc":
            lv.append("supported_single_dataset" if dec.decision.startswith("cortisol adds") else "null_single_dataset")
        elif r.effect_measure == "auroc":
            lv.append("supported_single_dataset" if (r.ci_low > 0.5 and r.p_value < 0.05) else "null_single_dataset")
        else:
            lv.append(_level(r.p_value, r.q_value, r.ci_low, r.ci_high, r.null_value))
    return add_provenance(sig.reset_index(drop=True), data_layer="person", source_name=SOURCE_NAME,
                          source_version=f"{DATA_FILE}; analysis {ANALYSIS_VERSION}", retrieved_at=retrieved_at,
                          evidence_type="person_derived_feature", source_record_id="object_id",
                          evidence_level=pd.Series(lv),
                          provenance_notes=f"Cohort-level contrast from person-level MY-LC rows; plan {PLAN}")


def write_tables(res: dict) -> None:
    T = TABLES
    res["perf"].to_csv(T / f"{SOURCE_ID}_model_performance.csv", index=False)
    res["deltas"].to_csv(T / f"{SOURCE_ID}_delta_auroc.csv", index=False)
    res["fitfree"].to_csv(T / f"{SOURCE_ID}_cortisol_fitfree.csv", index=False)
    res["tertiles"].to_csv(T / f"{SOURCE_ID}_drawtime_tertile_composition.csv", index=False)
    res["decision"].to_csv(T / f"{SOURCE_ID}_primary_decision.csv", index=False)
    res["effects"].sort_values("p_value").to_csv(T / f"{SOURCE_ID}_analyte_effects.csv", index=False)
    pd.DataFrame({"permutation": np.arange(len(res["null"])), "cortisol_only_auroc_null": res["null"],
                  "panel_auroc_null": res["panel_null"]}).to_csv(T / f"{SOURCE_ID}_permutation_null.csv", index=False)


def figures(res: dict) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = FIGURES / "drafts"
    out.mkdir(parents=True, exist_ok=True)
    perf = res["perf"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for i, r in enumerate(perf.itertuples()):
        ax[0].errorbar(r.auroc, i, xerr=[[r.auroc - r.auroc_ci_low], [r.auroc_ci_high - r.auroc]], fmt="o", capsize=3,
                       color="#c0392b" if i == 1 else "#444")
    ax[0].axvline(0.5, ls="--", color="grey")
    ax[0].set_yticks(range(len(perf)), [m[:45] for m in perf.model])
    ax[0].set_xlabel("cross-validated AUROC (95% CI)")
    ax[1].hist(res["null"], bins=40, color="#bbb")
    ax[1].axvline(perf.auroc.iloc[1], color="#c0392b", lw=2)
    ax[1].set_xlabel("cortisol-only AUROC under label shuffling")
    fig.tight_layout()
    p = out / f"{SOURCE_ID}_auroc_and_null.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return [p]


def write_registry_and_audit(dd: dict) -> None:
    man = load_manifest(SOURCE_ID)["files"][DATA_FILE]
    part = dd["part"]
    a = part[part.in_analysis_set]
    d = dd["raw"]
    cort_n = d.loc[d.x0_Censor_Complete.eq(0)].groupby("x0_Censor_Cohort_ID")[CORT].apply(lambda s: s.notna().sum())
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "MY-LC integrated per-participant table (Klein et al. 2023 Nature, Suppl. Table 3)",
        "publisher": "Nature / Springer Nature; Yale University and Icahn School of Medicine at Mount Sinai",
        "landing_url": LANDING_URL, "access_urls": [DATA_URL], "license": "CC BY 4.0 (article supplementary information)",
        "access_conditions": "Open download (Springer static content; browser User-Agent), no registration (2026-09-24).",
        "retrieved_at": man["retrieved_at"], "source_version": f"{DATA_FILE} sha256 {man['sha256']}",
        "update_date": "static supplement (2023)", "data_layer": "person", "unit_of_observation": "participant",
        "sample_size": {"rows": int(len(part)), "analysis_set": int(len(a)),
                        **{f"analysis_set_{k}": int(v) for k, v in a.group_label.value_counts().items()},
                        **{f"cortisol_{COHORT[k]}": int(v) for k, v in cort_n.items()}},
        "geographic_resolution": "none", "person_level": True, "geographic": False, "omics": True, "wearable": False,
        "participant_linkage": "all 7,431 columns share x0_LC_ID (one row per person)",
        "true_participant_linkage_across_modalities": True, "status": "ingested",
        "processed_outputs": [f"participants__{SOURCE_ID}", f"participant_labs__{SOURCE_ID}",
                              f"phenotype_signatures__{SOURCE_ID}", f"results/tables/{SOURCE_ID}_*.csv"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": f"measure_it.labs.{SOURCE_ID}",
        "limitations": ["recruitment-channel case-control", "single draw; cortisol diurnal", "no codebook in the file"],
        "conditions": ["long_covid"], "measurement_class_ids": [MEASUREMENT_CLASS_ID, "immune_assays"], "citation": CITATION,
        "notes": "Ranked (tied 2nd/3rd) in docs/LAB_DATASET_DISCOVERY.md; plan " + PLAN})
    miss = {c: int(part.loc[part.in_analysis_set, c].isna().sum()) for c in COVS}
    text = f"""# DATA AUDIT — {SOURCE_NAME}

| Field | Value |
|---|---|
| source_id | {SOURCE_ID} |
| Source (dataset/API name, exact files/endpoints) | {DATA_URL} (Supplementary Table 3, sheet `Sheet1`) |
| Publishing organization | Nature (Springer Nature); Yale / Mount Sinai MY-LC investigators |
| Retrieval date (UTC) | {man['retrieved_at']} |
| Source version / release | {DATA_FILE}, {man['bytes']} bytes, sha256 {man['sha256']} |
| Source update date / cadence | static supplement |
| License / access conditions | CC BY 4.0; anonymous download (browser User-Agent) |
| Unit of observation | participant (one blood draw) |
| Sample size (actual, as ingested) | 185 rows; analysis set (x0_Censor_Complete == 0) {len(a)}: {a.group_label.value_counts().to_dict()} |
| Geography (resolution, vintage) | none |
| Person-level? | yes |
| Geographic? | no |
| Omics? | yes (immune profiling: cytokines/hormones, flow, antibody reactivity) |
| Wearable? | no |
| True participant linkage across modalities? | yes: every column is on the same row, keyed by x0_LC_ID |

## Files / endpoints retrieved

| file | bytes | sha256 | url |
|---|---|---|---|
| {DATA_FILE} | {man['bytes']} | {man['sha256']} | {DATA_URL} |

## Key variables

* `x0_LC_ID` person id; `x0_Censor_Cohort_ID` 1 = HC, 2 = CC, 3 = LC (inferred: cohort 1 has no infection-test fields,
  35/42 of cohort 2 had a positive test; the paper's analysed n 40/39/99 is reproduced by `x0_Censor_Complete == 0`).
* `x0_Sample_Time_Min` draw time (minutes, clock time assumed), `x0_Demographics_Age/Sex/BMI`.
* 144 `x1_Cytokines_*` analytes with units in the column name, incl. `{CORT}` (ng/mL).
* Exclusion flags `x0_Censor_Oral_Steroid`, `x0_Censor_Pit_Adre_Dysfunction` (used in S4).

## Missingness (analysis set, n = {len(a)})

Cortisol non-missing by cohort: {json.dumps({COHORT[k]: int(v) for k, v in cort_n.items()})}. Covariates missing:
{json.dumps(miss)}.

## Linkage strategy

Single table; participant id `{SOURCE_ID}:<x0_LC_ID>`. Not linkable to any other dataset.

## Limitations and caveats

{CAVEAT} No codebook is shipped with the table; the cohort-code mapping and the draw-time unit are inferred.

## Processed outputs

`participants__{SOURCE_ID}` (185), `participant_labs__{SOURCE_ID}` ({185 * 144}), `phenotype_signatures__{SOURCE_ID}`.

## Reproduce

`uv run python -m measure_it.labs.{SOURCE_ID}`
"""
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text(text)


def _f(x, d=3):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def write_report() -> str:
    T = lambda n: pd.read_csv(TABLES / f"{SOURCE_ID}_{n}.csv")  # noqa: E731
    perf, dl, ff, dec, eff = T("model_performance"), T("delta_auroc"), T("cortisol_fitfree"), \
        T("primary_decision").iloc[0], T("analyte_effects")
    meta = json.loads((TABLES / f"{SOURCE_ID}_run_metadata.json").read_text())
    pr = dl.iloc[0]

    def ci(r, a="auroc"):
        return f"{_f(r[a])} ({_f(r[a + '_ci_low'])} to {_f(r[a + '_ci_high'])})"
    prow = "\n".join(f"| {r.model} | {int(r.n_cases)} / {int(r.n_controls)} | {ci(r._asdict())} | "
                     f"{E.pfmt(r.perm_p, r.n_perm) if pd.notna(r.perm_p) else 'not run'} | {_f(r.sens_at_90spec, 2)} |"
                     for r in perf.itertuples())
    frow = "\n".join(f"| {r.contrast} | {int(r.n_cases)} / {int(r.n_controls)} | {_f(r.auroc)} ({_f(r.auroc_ci_low)} to "
                     f"{_f(r.auroc_ci_high)}) | {E.pfmt(r.perm_p, r.n_perm)} | {_f(r.sens_at_90spec, 2)} |" for r in ff.itertuples())
    nq = int((eff.q_value < 0.05).sum())
    top = eff.sort_values("p_value").head(15)
    trow = "\n".join(f"| {r.analyte_name} | {r.unit} | {int(r.n_cases)} / {int(r.n_controls)} | {_f(r.median_cases, 2)} / "
                     f"{_f(r.median_controls, 2)} | {_f(r.hedges_g, 2)} ({_f(r.g_ci_low, 2)} to {_f(r.g_ci_high, 2)}) | "
                     f"{r.p_value:.2e} | {_f(r.q_value, 4)} |" for r in top.itertuples())
    s4 = dl.iloc[1]
    s6 = ff.iloc[0]
    e1r = ff[ff.contrast.str.startswith("E1")].iloc[0]
    cg = eff[eff.analyte_name == "Cortisol"].iloc[0]
    text = f"""# Serum cortisol and plasma mediators in Long COVID (MY-LC): pre-specified re-analysis

Dataset `{SOURCE_ID}`: Klein et al., Nature 2023, Supplementary Table 3, one row per participant. Analysis set
(`x0_Censor_Complete == 0`): 99 LC, 40 HC, 39 CC. Plan: `{PLAN}` (written before any analyte was compared). Analysis
`{meta['analysis_version']}`, computed {meta['computed_at']}.

## Answer first

**Primary: adding serum cortisol to draw time, age, sex and BMI changes cross-validated AUROC by
{_f(pr.delta_auroc)} (95% CI {_f(pr.ci_low)} to {_f(pr.ci_high)}); cortisol alone: AUROC
{_f(dec.cortisol_only_auroc)} ({_f(dec.cortisol_only_ci_low)} to {_f(dec.cortisol_only_ci_high)}), label-permutation
p {E.pfmt(dec.cortisol_only_perm_p, int(dec.n_perm))} ({int(dec.n_perm)} shuffles, null 95th percentile {_f(dec.null_q95)}).** Sample:
{int(dec.n_cases)} LC vs {int(dec.n_controls)} controls with complete cortisol and covariates.
Decision by the locked rule: **{dec.decision}**.

* Fit-free cortisol AUROC (lower = LC) on the analysis set (S6): {_f(s6.auroc)} ({_f(s6.auroc_ci_low)} to
  {_f(s6.auroc_ci_high)}). The paper reports 0.96 (0.93-0.99) on its Gale-Shapley-matched subset.
* Excluding steroid / pituitary-adrenal flags (S4): Delta {_f(s4.delta_auroc)} ({_f(s4.ci_low)} to {_f(s4.ci_high)}).
* {nq} of 144 plasma analytes differ at BH q < 0.05 (S1).

## How to read this

* The released rows **reproduce the published cortisol result**: fit-free AUROC {_f(s6.auroc)} and CV AUROC
  {_f(dec.cortisol_only_auroc)} against the paper's 0.96. The primary increment over draw time and demographics is
  large, because cortisol alone already carries almost all the separation (model 3 {_f(perf.auroc.iloc[2])} vs model 2
  {_f(perf.auroc.iloc[1])}).
* **The size of the effect is itself a caveat.** Cortisol in LC vs controls has Hedges g {_f(cg.hedges_g, 2)} (on the
  released z-scores). That is a very large separation for a single draw of a diurnal, pulsatile hormone. The table
  releases **batch-integrated z-scores**, not concentrations, and no assay-batch or plate column. So a group-by-batch
  or recruitment-channel difference can be **neither checked nor excluded** from the public file.
* The covariates alone separate the groups (AUROC {_f(perf.auroc.iloc[0])}), so draw time and demographics differ by
  group. **Post hoc, exploratory (E1, not part of the decision):** within draw-time tertiles the fit-free cortisol
  AUROC is {_f(e1r.auroc)} ({_f(e1r.auroc_ci_low)} to {_f(e1r.auroc_ci_high)}). Group counts per tertile are in
  `results/tables/{SOURCE_ID}_drawtime_tertile_composition.csv`.
* This is the same released data re-analysed, not an independent replication. Whether hypocortisolaemia marks Long
  COVID outside this recruitment design (clinic patients vs advertised volunteers) needs other cohorts. The open
  candidates for that check are listed in `docs/LAB_DATASET_DISCOVERY.md`: Appelman 2024, with cortisol and minutes
  since waking; and a PLOS ONE cohort with 10 vs 7 people.

## Published claim vs our computation

| | authors (published) | this analysis (same released rows) |
|---|---|---|
| Cortisol | lower in LC, similar in HC and CC; "cortisol alone achieved an AUC of 0.96 (95% CI 0.93-0.99)" in a matched subset; significant after adjusting for demographics and sample-collection time | fit-free AUROC {_f(s6.auroc)} on the full analysis set; CV cortisol-only {_f(dec.cortisol_only_auroc)}; increment over draw time + demographics {_f(pr.delta_auroc)} ({_f(pr.ci_low)} to {_f(pr.ci_high)}) |
| Other mediators | C4b, CCL19, CCL20, galectin-1, CCL4, APRIL, LH higher; IL-5 lower (Kruskal-Wallis across 3 groups) | {nq} of 144 at q < 0.05, LC vs pooled controls (table below) |

## Models (CV, 5-fold x 20)

| model | cases / controls | AUROC (95% CI) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
{prow}

## Fit-free cortisol AUROC by control type (S2, S6)

| contrast | LC / controls | AUROC (lower cortisol = LC) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
{frow}

## Top 15 plasma analytes (S1, sorted by p)

| analyte | unit named in column | n LC / controls | median z LC / controls | Hedges g on z (95% CI) | Mann-Whitney p | BH q |
|---|---|---|---|---|---|---|
{trow}

## Caveats carried with every number

{CAVEAT} The cohort-code mapping (1 = HC, 2 = CC, 3 = LC) and draw-time unit are inferred from the file (no codebook);
the analysed n reproduces the paper's. Bootstrap CIs hold out-of-fold predictions fixed.

## Outputs

`participants__{SOURCE_ID}`, `participant_labs__{SOURCE_ID}`, `phenotype_signatures__{SOURCE_ID}` ({meta['signature_rows']} rows);
`results/tables/{SOURCE_ID}_{{model_performance, delta_auroc, cortisol_fitfree, primary_decision, analyte_effects,
permutation_null}}.csv`; figure {', '.join(meta['figures'])}.

## Reproduce

```bash
uv run python -m measure_it.labs.{SOURCE_ID} --jobs 16     # deterministic, seed {SEED}
uv run pytest tests/test_{SOURCE_ID}.py
```
"""
    path = PROJECT_ROOT / REPORT
    path.write_text(text)
    return str(path)


def run(jobs: int = 8, n_perm: int = 1000, log=print) -> dict:
    dd = ingest()
    res = analyse(dd, jobs, n_perm, log)
    write_tables(res)
    figs = figures(res)
    sig = signatures(res, dd["retrieved_at"])
    write_table(sig, f"phenotype_signatures__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="MY-LC: 144 analyte effects, cortisol CV AUROC and primary Delta")
    write_registry_and_audit(dd)
    meta = {"analysis_version": ANALYSIS_VERSION, "plan": PLAN, "computed_at": utc_now_iso(), "seed": SEED,
            "n_perm": n_perm, "signature_rows": len(sig), "figures": [str(p.relative_to(PROJECT_ROOT)) for p in figs]}
    (TABLES / f"{SOURCE_ID}_run_metadata.json").write_text(json.dumps(meta, indent=2))
    log(f"report: {write_report()}")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args()
    print(write_report()) if a.report_only else run(a.jobs, a.n_perm)
