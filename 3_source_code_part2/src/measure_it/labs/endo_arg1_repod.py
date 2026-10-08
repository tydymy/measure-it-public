"""Serum arginase-1 in endometriosis vs surgical and healthy controls (RepOD doi:10.18150/PY1P9X).

Per-person file behind Pliszkiewicz M, Czystowska-Kuzmicz M, Soroczynska K, Siekierski BP, Safranow K, et al.
"Determination of serum arginase-1 concentrations and serum arginase activity for the non-invasive diagnosis of
endometriosis." J Clin Med 2024;13:1489 (PMC10933979). One tab-separated table: Group (Patient = surgically confirmed
endometriosis; Ctrl = benign gynaecological surgery without endometriosis; Ctrl2 = uterine myomas; Ctrl3 = healthy
women), rASRM stage, pre/post-operative ARG1 and ARG2 (ng/mL), arginase activity (U/L), glycaemia, creatinine,
AST, ALT, ELISA kit / plate numbers.

Pre-specified primary analysis: docs/ANALYSIS_PLAN_ENDO_ARG1_REPOD.md.

Run: uv run python -m measure_it.labs.endo_arg1_repod
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats

from measure_it.config import PROJECT_ROOT, SEED, TABLES, raw_dir, utc_now_iso
from measure_it.download import download_file, load_manifest
from measure_it.labs import lab_cc_engine as E
from measure_it.provenance import add_provenance
from measure_it.registry import write_registry_entry
from measure_it.store import write_table

SOURCE_ID = "endo_arg1_repod"
DATA_URL = "https://repod.icm.edu.pl/api/access/datafile/48824"
DATA_FILE = "endo_arg1.tab"
LANDING_URL = "https://repod.icm.edu.pl/dataset.xhtml?persistentId=doi:10.18150/PY1P9X"
SOURCE_NAME = "RepOD PY1P9X: serum arginase in endometriosis (Pliszkiewicz et al. 2024 J Clin Med data)"
CITATION = ("Pliszkiewicz M, Czystowska-Kuzmicz M, Soroczynska K, Siekierski BP, Safranow K, et al. Determination of "
            "serum arginase-1 concentrations and serum arginase activity for the non-invasive diagnosis of "
            "endometriosis. J Clin Med 2024;13:1489. Data: RepOD doi:10.18150/PY1P9X")
CONDITION_ID = "endometriosis"
MEASUREMENT_CLASS_ID = "blood_biomarkers"
PLAN = "docs/ANALYSIS_PLAN_ENDO_ARG1_REPOD.md"
REPORT = "results/ENDO_ARG1_REPOD_RESULTS.md"
ANALYSIS_VERSION = f"{SOURCE_ID} 2026-09-24.1"
BUA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"}
ARG1 = "Arg1PREng/ml"
PANEL = {"Arg1PREng/ml": ("ARG1, pre-operative", "ng/mL"), "Argact.PREU/L": ("arginase activity, pre-operative", "U/L"),
         "Arg2PREng/ml": ("ARG2, pre-operative (authors: not reproducible)", "ng/mL"),
         "glycemia": ("glycaemia", "mg/dL"), "creatinine": ("creatinine", "mg/dL"), "AspAT": ("AST", "U/L"),
         "AlAT": ("ALT", "U/L")}
POST = {"Arg1POSTng/ml": ("ARG1, post-operative", "ng/mL"), "Argact.POSTU/L": ("arginase activity, post-operative", "U/L"),
        "Arg2POSTng/ml": ("ARG2, post-operative", "ng/mL")}
LABEL_BASIS = ("Patient = endometriosis confirmed at laparoscopy with rASRM stage I-IV; Ctrl = women operated on for other "
               "benign gynaecological conditions (file 'diagnosis' column: myoma, cysts, pain, diagnostic laparoscopy, "
               "...); Ctrl2 = uterine myomas; Ctrl3 = healthy women without features of endometriosis (51 'healthy' + 3 "
               "other). Paper groups: GB (endometriosis), K1 (surgical controls), K2 (healthy).")
CAVEAT = ("Single-centre (Warsaw) case-control; the healthy controls (Ctrl3) have no ELISA plate recorded, so a batch "
          "difference cannot be excluded for patient-vs-healthy contrasts; the released file has more people than the "
          "paper (127/34/54 vs 105/22/53); surgical controls are few.")


def fetch():
    return download_file(DATA_URL, SOURCE_ID, DATA_FILE, headers=BUA, max_retries=3)


def read_raw() -> pd.DataFrame:
    d = pd.read_csv(raw_dir(SOURCE_ID) / DATA_FILE, sep="\t")
    need = {"Group", "Nr", "age", "Grading", "diagnosis", "platenr", "ELISAkitnr", *PANEL, *POST}
    if not need <= set(d.columns):
        raise ValueError(f"missing columns: {need - set(d.columns)}")
    return d


def participant_id(g: str, nr) -> str:
    return f"{SOURCE_ID}:{g}-{nr}"


def num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.replace(",", ".").str.strip().replace({"n/a": None, "N/A": None}),
                         errors="coerce")


def build_frames(d: pd.DataFrame, retrieved_at: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = d[d.Group.notna()].copy()
    d["native_id"] = d.Group.astype(str) + "-" + d.Nr.astype(str)
    if d.native_id.duplicated().any():
        raise ValueError("duplicate Group-Nr ids")
    grp = d.Group.astype(str)
    part = pd.DataFrame({
        "participant_id": d.native_id.map(lambda x: f"{SOURCE_ID}:{x}"), "dataset_id": SOURCE_ID, "native_id": d.native_id,
        "group_label": grp.map({"Patient": "endometriosis", "Ctrl": "surgical_control", "Ctrl2": "surgical_control_myoma",
                                "Ctrl3": "healthy_control"}),
        "source_group": grp, "endometriosis": (grp == "Patient").astype(int),
        "condition_id": np.where(grp == "Patient", CONDITION_ID, ""), "label_basis": LABEL_BASIS,
        "diagnosis_text": d.diagnosis.astype(str), "rasrm_stage": d.Grading.astype(str),
        "age_years": num(d.age), "elisa_plate": num(d.platenr), "elisa_kit": num(d.ELISAkitnr),
        "sex": "female"})
    part = add_provenance(part.reset_index(drop=True), data_layer="person", source_name=SOURCE_NAME,
                          source_version=DATA_FILE, retrieved_at=retrieved_at, evidence_type="person_lab_measurement",
                          source_record_id="native_id", provenance_notes="one row per woman; ids = Group-Nr")
    rows = []
    for col, (name, unit) in {**PANEL, **POST}.items():
        v = num(d[col])
        rows.append(pd.DataFrame({"native_id": d.native_id.values, "lab_variable": col, "value": v.values,
                                  "lab_name": name, "unit": unit}))
    labs = pd.concat(rows, ignore_index=True)
    labs = labs.assign(participant_id=labs.native_id.map(lambda x: f"{SOURCE_ID}:{x}"), dataset_id=SOURCE_ID,
                       specimen="serum", harmonized_name=labs.lab_variable.str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True),
                       lab_note=np.where(labs.value.isna(), "missing / n/a", ""),
                       source_record_id=labs.native_id + ":" + labs.lab_variable)
    labs = add_provenance(labs, data_layer="person", source_name=SOURCE_NAME, source_version=DATA_FILE,
                          retrieved_at=retrieved_at, evidence_type="person_lab_measurement",
                          source_record_id="source_record_id", provenance_notes="values as released ('n/a' -> missing)")
    return part, labs


def ingest() -> dict:
    fetch()
    ret = load_manifest(SOURCE_ID)["files"][DATA_FILE]["retrieved_at"]
    d = read_raw()
    part, labs = build_frames(d, ret)
    write_table(part, f"participants__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="endometriosis / control women: group, stage, age, ELISA plate")
    write_table(labs, f"participant_labs__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="arginase and routine labs, long format")
    wide = labs.pivot(index="native_id", columns="lab_variable", values="value")
    return {"raw": d, "part": part, "wide": part.set_index("native_id").join(wide), "retrieved_at": ret}


def analyse(dd: dict, n_perm: int = 10000, log=print) -> dict:
    w = dd["wide"].copy()
    surg = w.source_group.isin(["Patient", "Ctrl", "Ctrl2"]) & w[ARG1].notna()
    ps = w[surg].copy()
    ps["y"] = ps.source_group.eq("Patient").astype(int)
    y = ps.y.to_numpy()
    log(f"primary: {int(y.sum())} endometriosis vs {int((1 - y).sum())} surgical controls")
    prim = E.fitfree_test(ps[ARG1].to_numpy(), y, n_perm=n_perm)
    # S1 plate-stratified
    ps_pl = ps[ps.elisa_plate.notna()]
    strat = E.fitfree_test(ps_pl[ARG1].to_numpy(), ps_pl.y.to_numpy(), n_perm=n_perm,
                           strata=ps_pl.elisa_plate.astype(int).astype(str).to_numpy())
    plates = ps_pl.groupby("elisa_plate").y.agg(["sum", "count"]).rename(columns={"sum": "cases", "count": "n"})
    plates["controls"] = plates.n - plates.cases
    # S2 healthy
    hs = w[w.source_group.isin(["Patient", "Ctrl3"]) & w[ARG1].notna()].copy()
    hs["y"] = hs.source_group.eq("Patient").astype(int)
    healthy = E.fitfree_test(hs[ARG1].to_numpy(), hs.y.to_numpy(), n_perm=n_perm)
    fit = pd.DataFrame([{"analysis": "PRIMARY endometriosis vs surgical controls (Ctrl + Ctrl2), pooled", **prim},
                        {"analysis": "S1 same, ELISA-plate-stratified (plates holding both groups)", **strat},
                        {"analysis": "S2 endometriosis vs healthy (Ctrl3; batch unknown)", **healthy}])
    # S3 age adjustment (CV)
    pa = ps.dropna(subset=["age_years"]).reset_index(drop=True)
    ya = pa.y.to_numpy()
    b = E.stratified_boot_idx(ya)
    m_age = E.cv_block(pa[["age_years"]].to_numpy(float), ya, 1.0, b)
    m_both = E.cv_block(np.column_stack([pa.age_years, np.log10(pa[ARG1].clip(lower=1e-3))]), ya, 1.0, b)
    dl = E.delta(m_both, m_age, ya)
    cv = pd.DataFrame([{"model": "age", **m_age["summary"]}, {"model": "age + log10 ARG1", **m_both["summary"]}])
    delta = pd.DataFrame([{"comparison": "S3 (age + log ARG1) minus age", **dl}])
    # S4 panel
    pan = ps.copy()
    pan["y"] = pan.y
    full = w[w.source_group.isin(["Patient", "Ctrl", "Ctrl2"])].copy()
    full["y"] = full.source_group.eq("Patient").astype(int)
    eff = E.per_analyte_effects(full, "y", list(PANEL), transform=lambda x: np.log10(np.clip(x, 1e-3, None)))
    eff["analyte_name"] = eff.analyte.map(lambda c: PANEL[c][0])
    eff["unit"] = eff.analyte.map(lambda c: PANEL[c][1])
    # S5 stage
    st = w[(w.source_group == "Patient") & w[ARG1].notna() & w.rasrm_stage.isin(["I", "II", "III", "IV"])]
    rho, sp = stats.spearmanr(st.rasrm_stage.map({"I": 1, "II": 2, "III": 3, "IV": 4}), st[ARG1])
    stage = pd.DataFrame([{"n": len(st), "spearman_rho": rho, "p_value": sp,
                           "stage_counts": json.dumps(st.rasrm_stage.value_counts().to_dict())}])
    ok = bool(prim["auroc_ci_low"] > 0.5 and prim["perm_p"] < 0.05)
    batch_sensitive = bool(prim["auroc_ci_low"] > 0.5 and not strat["auroc_ci_low"] > 0.5)
    decision = pd.DataFrame([{
        "question": "Does pre-operative serum ARG1 separate endometriosis from surgical non-endometriosis controls?",
        "auroc": prim["auroc"], "ci_low": prim["auroc_ci_low"], "ci_high": prim["auroc_ci_high"],
        "perm_p": prim["perm_p"], "n_perm": n_perm, "n_cases": prim["n_cases"], "n_controls": prim["n_controls"],
        "plate_stratified_auroc": strat["auroc"], "plate_stratified_ci_low": strat["auroc_ci_low"],
        "plate_stratified_ci_high": strat["auroc_ci_high"], "batch_sensitive": batch_sensitive,
        "published_auc_gb_vs_k1": 0.848, "published_ci": "0.769-0.926",
        "rule": "separates if CI low > 0.5 AND permutation p < 0.05; batch-sensitive if the plate-stratified CI low <= 0.5",
        "decision": ("ARG1 separates endometriosis from surgical controls (pre-specified rule met)" if ok
                     else "not shown (pre-specified rule not met)") + ("; BATCH-SENSITIVE" if batch_sensitive else ""),
        "caveat": CAVEAT, "plan": PLAN, "analysis_version": ANALYSIS_VERSION, "computed_at": utc_now_iso()}])
    log(decision[["auroc", "ci_low", "ci_high", "perm_p", "plate_stratified_auroc", "decision"]].to_string())
    return {"fitfree": fit, "plates": plates.reset_index(), "cv": cv, "delta": delta, "effects": eff, "stage": stage,
            "decision": decision}


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
        rows.append({"feature": "serum_" + r.analyte.lower().replace(".", "_").replace("/", "_per_"),
                     "feature_label": f"{r.analyte_name} ({r.unit}), log10", "effect_measure": "hedges_g",
                     "effect_size": r.hedges_g, "effect_unit": "Hedges' g on log10 (endometriosis minus surgical controls)",
                     "ci_low": r.g_ci_low, "ci_high": r.g_ci_high, "p_value": r.p_value, "q_value": r.q_value,
                     "n_cases": r.n_cases, "n_controls": r.n_controls, "null_value": 0.0,
                     "method": "S4: stratified bootstrap CI (2,000); two-sided Mann-Whitney p; BH q over 7 analytes"})
    for r in res["fitfree"].itertuples():
        tag = "primary" if r.analysis.startswith("PRIMARY") else r.analysis.split()[0].lower()
        rows.append({"feature": f"fitfree_auroc:arg1_pre:{tag}", "feature_label": r.analysis,
                     "effect_measure": "auroc", "effect_size": r.auroc, "effect_unit": "fit-free AUROC, higher ARG1 = case",
                     "ci_low": r.auroc_ci_low, "ci_high": r.auroc_ci_high, "p_value": r.perm_p, "q_value": np.nan,
                     "n_cases": r.n_cases, "n_controls": r.n_controls, "null_value": 0.5,
                     "method": "stratified person bootstrap CI (2,000); one-sided label-permutation p"})
    sig = pd.DataFrame(rows)
    sig.insert(0, "object_id", "signature:" + SOURCE_ID + "|endometriosis|" + sig.feature)
    sig.insert(1, "dataset_id", SOURCE_ID)
    sig.insert(2, "phenotype_id", "endometriosis")
    sig.insert(3, "phenotype_label", "Surgically confirmed endometriosis vs benign-surgery controls (Warsaw)")
    sig.insert(4, "condition_id", CONDITION_ID)
    sig.insert(5, "is_proxy", False)
    sig["label_basis"] = LABEL_BASIS
    sig["phenotype_definition"] = "laparoscopy-confirmed endometriosis (rASRM I-IV) vs women operated for other benign conditions"
    sig["caveats"] = CAVEAT
    sig["measurement_class_id"] = MEASUREMENT_CLASS_ID
    sig["n_cases"] = sig.n_cases.astype("Int64")
    sig["n_controls"] = sig.n_controls.astype("Int64")
    lv = [("supported_single_dataset" if (r.ci_low > 0.5 and r.p_value < 0.05) else "null_single_dataset")
          if r.effect_measure == "auroc" else _level(r.p_value, r.q_value, r.ci_low, r.ci_high, r.null_value)
          for r in sig.itertuples()]
    return add_provenance(sig, data_layer="person", source_name=SOURCE_NAME,
                          source_version=f"{DATA_FILE}; analysis {ANALYSIS_VERSION}", retrieved_at=retrieved_at,
                          evidence_type="person_derived_feature", source_record_id="object_id",
                          evidence_level=pd.Series(lv), provenance_notes=f"Cohort-level contrast; plan {PLAN}")


def write_tables(res: dict) -> None:
    T = TABLES
    res["fitfree"].to_csv(T / f"{SOURCE_ID}_fitfree_auroc.csv", index=False)
    res["plates"].to_csv(T / f"{SOURCE_ID}_plate_composition.csv", index=False)
    res["cv"].to_csv(T / f"{SOURCE_ID}_cv_age_models.csv", index=False)
    res["delta"].to_csv(T / f"{SOURCE_ID}_delta_auroc.csv", index=False)
    res["effects"].to_csv(T / f"{SOURCE_ID}_analyte_effects.csv", index=False)
    res["stage"].to_csv(T / f"{SOURCE_ID}_stage_correlation.csv", index=False)
    res["decision"].to_csv(T / f"{SOURCE_ID}_primary_decision.csv", index=False)


def write_registry_and_audit(dd: dict) -> None:
    man = load_manifest(SOURCE_ID)["files"][DATA_FILE]
    part = dd["part"]
    w = dd["wide"]
    counts = part.source_group.value_counts().to_dict()
    arg_n = w.groupby("source_group")[ARG1].apply(lambda s: int(s.notna().sum())).to_dict()
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "Serum arginase in endometriosis vs controls (RepOD PY1P9X)",
        "publisher": "RepOD (ICM University of Warsaw Dataverse); Medical University of Warsaw",
        "landing_url": LANDING_URL, "access_urls": [DATA_URL],
        "license": "not stated in the file API response; RepOD default terms (check landing page)",
        "access_conditions": "Open anonymous download via the Dataverse access API (verified 2026-09-24).",
        "retrieved_at": man["retrieved_at"], "source_version": f"{DATA_FILE} sha256 {man['sha256']}",
        "update_date": "static deposit", "data_layer": "person", "unit_of_observation": "participant",
        "sample_size": {"rows": int(len(part)), **{str(k): int(v) for k, v in counts.items()},
                        **{f"arg1_pre_{k}": v for k, v in arg_n.items()}},
        "geographic_resolution": "none", "person_level": True, "geographic": False, "omics": False, "wearable": False,
        "participant_linkage": "single table", "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": [f"participants__{SOURCE_ID}", f"participant_labs__{SOURCE_ID}",
                              f"phenotype_signatures__{SOURCE_ID}", f"results/tables/{SOURCE_ID}_*.csv"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": f"measure_it.labs.{SOURCE_ID}",
        "limitations": ["healthy controls have no ELISA plate recorded", "more people than the paper", "few surgical controls",
                        "licence not stated in the API"],
        "conditions": ["endometriosis"], "measurement_class_ids": [MEASUREMENT_CLASS_ID], "citation": CITATION,
        "notes": "Tied rank 2 in docs/LAB_DATASET_DISCOVERY.md; plan " + PLAN})
    miss = {c: {g: int(w.loc[w.source_group == g, c].isna().sum()) for g in counts} for c in PANEL}
    text = f"""# DATA AUDIT — {SOURCE_NAME}

| Field | Value |
|---|---|
| source_id | {SOURCE_ID} |
| Source (dataset/API name, exact files/endpoints) | {DATA_URL} (RepOD Dataverse datafile 48824, tab-separated) |
| Publishing organization | RepOD (University of Warsaw ICM); Medical University of Warsaw |
| Retrieval date (UTC) | {man['retrieved_at']} |
| Source version / release | {DATA_FILE}, {man['bytes']} bytes, sha256 {man['sha256']} |
| Source update date / cadence | static deposit |
| License / access conditions | licence not stated in the API response; anonymous download |
| Unit of observation | participant (woman) |
| Sample size (actual, as ingested) | {len(part)} women with a Group: {json.dumps(counts)}; with pre-operative ARG1: {json.dumps(arg_n)} (2 rows without Group dropped) |
| Geography (resolution, vintage) | none |
| Person-level? | yes |
| Geographic? | no |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | single table |

## Files / endpoints retrieved

| file | bytes | sha256 | url |
|---|---|---|---|
| {DATA_FILE} | {man['bytes']} | {man['sha256']} | {DATA_URL} |

## Key variables

`Group` (Patient / Ctrl / Ctrl2 / Ctrl3), `diagnosis` (free text; Ctrl = other benign gynaecological surgery, Ctrl2 =
myomas, Ctrl3 = 'healthy'), `Grading` (rASRM 0/I-IV), `age`, `Arg1PREng/ml` (primary; the `...old` column is an
earlier assay run and is not used), `Arg1POSTng/ml`, `Arg2PRE/POST`, `Argact.PRE/POSTU/L`, `glycemia`, `creatinine`,
`AspAT`, `AlAT`, `ELISAkitnr`, `platenr`.

## Missingness

Missing values per analyte and group: {json.dumps(miss)}. ELISA plate is missing for all 54 Ctrl3 (healthy) women.

## Linkage strategy

Participant id `{SOURCE_ID}:<Group>-<Nr>`. No link to other datasets.

## Limitations and caveats

{CAVEAT}

## Processed outputs

`participants__{SOURCE_ID}` ({len(part)}), `participant_labs__{SOURCE_ID}`, `phenotype_signatures__{SOURCE_ID}`.

## Reproduce

`uv run python -m measure_it.labs.{SOURCE_ID}`
"""
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text(text)


def _f(x, d=3):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def write_report() -> str:
    T = lambda n: pd.read_csv(TABLES / f"{SOURCE_ID}_{n}.csv")  # noqa: E731
    ff, pl, cv, dl, eff, stg, dec = (T("fitfree_auroc"), T("plate_composition"), T("cv_age_models"), T("delta_auroc").iloc[0],
                                     T("analyte_effects"), T("stage_correlation").iloc[0], T("primary_decision").iloc[0])
    meta = json.loads((TABLES / f"{SOURCE_ID}_run_metadata.json").read_text())
    frow = "\n".join(f"| {r.analysis} | {int(r.n_cases)} / {int(r.n_controls)} | {_f(r.auroc)} ({_f(r.auroc_ci_low)} to "
                     f"{_f(r.auroc_ci_high)}) | {E.pfmt(r.perm_p, r.n_perm)} | {_f(r.sens_at_90spec, 2)} |" for r in ff.itertuples())
    erow = "\n".join(f"| {r.analyte_name} | {int(r.n_cases)} / {int(r.n_controls)} | {_f(r.median_cases, 2)} / "
                     f"{_f(r.median_controls, 2)} | {_f(r.auroc_case_higher)} | {_f(r.hedges_g, 2)} ({_f(r.g_ci_low, 2)} to "
                     f"{_f(r.g_ci_high, 2)}) | {r.p_value:.2e} | {_f(r.q_value, 4)} |" for r in eff.dropna(subset=['p_value']).itertuples())
    prow = "\n".join(f"| {int(r.elisa_plate)} | {int(r.cases)} | {int(r.controls)} |" for r in pl.itertuples())
    s1, s2 = ff.iloc[1], ff.iloc[2]
    text = f"""# Serum arginase-1 in endometriosis: pre-specified re-analysis of the released rows

Dataset `{SOURCE_ID}` (RepOD doi:10.18150/PY1P9X; Pliszkiewicz et al., J Clin Med 2024). Plan: `{PLAN}` (written
before any arginase value was compared). Analysis `{meta['analysis_version']}`, computed {meta['computed_at']}.

## Answer first

**Primary: pre-operative serum ARG1, endometriosis ({int(dec.n_cases)}) vs women operated on for other benign
gynaecological conditions ({int(dec.n_controls)}): fit-free AUROC {_f(dec.auroc)} (95% CI {_f(dec.ci_low)} to
{_f(dec.ci_high)}), label-permutation p {E.pfmt(dec.perm_p, int(dec.n_perm))} ({int(dec.n_perm)} shuffles).** Decision by the locked
rule: **{dec.decision}**. Published (paper, GB vs K1, 105 vs 22): AUC 0.848 (0.769-0.926).

* **Batch check (S1):** within ELISA plates that hold both groups, AUROC {_f(s1.auroc)} ({_f(s1.auroc_ci_low)} to
  {_f(s1.auroc_ci_high)}), permutation p {E.pfmt(s1.perm_p, int(s1.n_perm))}. The pooled result is therefore
  not explained by plate: surgical controls and patients were assayed on shared plates.
* **Healthy controls (S2, batch unknown):** AUROC {_f(s2.auroc)} ({_f(s2.auroc_ci_low)} to {_f(s2.auroc_ci_high)}); the
  paper reports 0.912. No ELISA plate is recorded for any healthy control, so this contrast cannot be checked for batch.
* **Age (S3, CV):** age alone AUROC {_f(cv.auroc.iloc[0])}; age + log ARG1 {_f(cv.auroc.iloc[1])}; Delta
  {_f(dl.delta_auroc)} ({_f(dl.ci_low)} to {_f(dl.ci_high)}). The Delta is large because age alone carries almost no
  information here.
* **This reproduces the published estimate in a larger released sample.** The file has 120 vs 32 people with the
  primary ARG1 value; the paper had 105 vs 22. We get 0.846 against the paper's 0.848, and 0.918 against 0.912 vs
  healthy.
* **Stage:** ARG1 does not track rASRM stage, which is noted here as a descriptive finding.
* **Stage (S5):** Spearman rho {_f(stg.spearman_rho, 2)} (p = {_f(stg.p_value, 3)}, n = {int(stg.n)}).

## Published claim vs our computation

| | authors (published) | this analysis (released rows) |
|---|---|---|
| Groups | 105 GB / 22 K1 / 53 K2 | {int(dec.n_cases)} endometriosis / {int(dec.n_controls)} surgical controls (Ctrl + Ctrl2) with ARG1 / {int(s2.n_controls)} healthy |
| ARG1 vs surgical controls | AUC 0.848 (0.769-0.926) | {_f(dec.auroc)} ({_f(dec.ci_low)} to {_f(dec.ci_high)}); plate-stratified {_f(s1.auroc)} |
| ARG1 vs healthy | AUC 0.912 (0.859-0.965) | {_f(s2.auroc)} ({_f(s2.auroc_ci_low)} to {_f(s2.auroc_ci_high)}), batch unknown |
| Arginase activity | AUC < 0.7, "not clinically useful" | see panel |

## Fit-free ARG1 AUROC

| analysis | cases / controls | AUROC (95% CI) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
{frow}

ELISA plate composition of the primary sample (plates with both groups enter S1):

| plate | endometriosis | surgical controls |
|---|---|---|
{prow}

## Panel (S4): endometriosis vs surgical controls

| analyte | n | median case / control | AUROC (higher = case) | Hedges g on log10 (95% CI) | Mann-Whitney p | BH q (7) |
|---|---|---|---|---|---|---|
{erow}

## Caveats carried with every number

{CAVEAT} Controls operated on for myomas, cysts or pain are the clinically closer comparator, but there are only
{int(dec.n_controls)}. Nothing here evaluates ARG1 in women with pelvic pain referred for suspected endometriosis,
which is the intended-use population.

## Outputs

`participants__{SOURCE_ID}`, `participant_labs__{SOURCE_ID}`, `phenotype_signatures__{SOURCE_ID}` ({meta['signature_rows']}
rows); `results/tables/{SOURCE_ID}_*.csv`.

## Reproduce

```bash
uv run python -m measure_it.labs.{SOURCE_ID}      # deterministic, seed {SEED}
uv run pytest tests/test_{SOURCE_ID}.py
```
"""
    path = PROJECT_ROOT / REPORT
    path.write_text(text)
    return str(path)


def run(n_perm: int = 10000, log=print) -> dict:
    dd = ingest()
    res = analyse(dd, n_perm, log)
    write_tables(res)
    sig = signatures(res, dd["retrieved_at"])
    write_table(sig, f"phenotype_signatures__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="endometriosis arginase: fit-free AUROCs and panel effects")
    write_registry_and_audit(dd)
    meta = {"analysis_version": ANALYSIS_VERSION, "plan": PLAN, "computed_at": utc_now_iso(), "seed": SEED,
            "n_perm": n_perm, "signature_rows": len(sig)}
    (TABLES / f"{SOURCE_ID}_run_metadata.json").write_text(json.dumps(meta, indent=2))
    log(f"report: {write_report()}")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args()
    print(write_report()) if a.report_only else run(a.n_perm)
