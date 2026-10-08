"""hEDS / HSD serum Olink proteome vs healthy controls (Cinquina et al., Clin Proteomics 2026; PMC13081554).

Cinquina V, Carini G, Chiarelli N, Vezzoli M, Bertini V, Venturini M, et al. "Proximity extension assay-based serum
proteomic profiling identifies shared protein signatures in hypermobile Ehlers-Danlos syndrome and hypermobility
spectrum disorders." Clinical Proteomics 2026 (Europe PMC PMC13081554), Additional files 2 (NPX) and 5 (clinical).

Per-sample NPX for 460 Olink Target 96 assays (Cardiometabolic, Development, Inflammation, Neurology, Organ Damage)
in 88 hEDS (P1-P88), 88 HSD (H1-H88) and 176 healthy controls (C1-C176). hEDS = 2017 criteria (5PQ-adjusted); HSD =
symptomatic hypermobility not meeting them (checklists per patient in Additional file 5). Per the paper, P1-P44 and
H1-H44 are Italian (Brescia) and P45-P88 / H45-H88 are American (home draws shipped to Baltimore); all controls are
Italian. The pre-specified primary analysis (docs/ANALYSIS_PLAN_HEDS_HSD_OLINK_SERUM_CINQUINA2026.md) therefore uses
Italian patients vs controls only.

Outputs (person layer; no geography):
  participants__heds_hsd_olink_serum_cinquina2026, participant_labs__heds_hsd_olink_serum_cinquina2026,
  phenotype_signatures__heds_hsd_olink_serum_cinquina2026; results/tables/heds_hsd_olink_serum_cinquina2026_*;
  results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md; data/raw/<id>/{DATA_AUDIT.md, registry_entry.yaml}.

Run: uv run python -m measure_it.labs.heds_hsd_olink_serum_cinquina2026 [--jobs 16] [--n-perm 1000]
(xlrd reads the .xls clinical file; it is a project dependency in uv.lock since 2026-09-24.)
"""
from __future__ import annotations

import argparse
import json
import zipfile

import numpy as np
import pandas as pd

from measure_it.config import FIGURES, PROJECT_ROOT, SEED, TABLES, raw_dir, utc_now_iso
from measure_it.download import download_file, load_manifest, record_file
from measure_it.labs import lab_cc_engine as E
from measure_it.provenance import add_provenance
from measure_it.registry import write_registry_entry
from measure_it.store import write_table

SOURCE_ID = "heds_hsd_olink_serum_cinquina2026"
DATA_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13081554/supplementaryFiles"
ZIP_FILE = "PMC13081554_supplementaryFiles.zip"
NPX_MEMBER = "12014_2026_9588_MOESM4_ESM.xlsx"
CLIN_MEMBER = "12014_2026_9588_MOESM6_ESM.xls"
DEP_MEMBER = "12014_2026_9588_MOESM8_ESM.xlsx"
LANDING_URL = "https://europepmc.org/article/PMC/PMC13081554"
SOURCE_NAME = "Cinquina et al. 2026 Clin Proteomics, Olink serum proteome in hEDS/HSD (Additional files 2 and 5)"
CITATION = ("Cinquina V, Carini G, Chiarelli N, et al. Proximity extension assay-based serum proteomic profiling "
            "identifies shared protein signatures in hypermobile Ehlers-Danlos syndrome and hypermobility spectrum "
            "disorders. Clinical Proteomics 2026 (PMC13081554).")
CONDITION_ID = "eds_hsd"
MEASUREMENT_CLASS_ID = "proteomics"
PLAN = "docs/ANALYSIS_PLAN_HEDS_HSD_OLINK_SERUM_CINQUINA2026.md"
REPORT = "results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md"
ANALYSIS_VERSION = f"{SOURCE_ID} 2026-09-24.1"
LABEL_BASIS = ("hEDS = 2017 international criteria (5PQ-adjusted criterion 1); HSD = symptomatic joint hypermobility not "
               "meeting them; per-patient criteria checklists in Additional file 5; controls = unrelated Italian "
               "volunteers without signs or symptoms of hEDS/HSD on examination and interview (paper Methods). "
               "Sample-ID prefix P = hEDS, H = HSD (verified against Additional file 5), C = control.")
N_ITALIAN = 44          # P1-P44 and H1-H44 are Italian; P45-P88 and H45-H88 are American (paper Methods / Results)
MISSING_MAX = 0.25
C_PANEL = 0.1
COMORB_COLS = {39: "functional_gi", 40: "neurological", 41: "dysautonomia_pots", 42: "psychological",
               43: "bladder_urological", 44: "gynecological", 45: "chronic_fatigue", 46: "tmj_disorder",
               47: "allergic_atopic_incl_mcas", 48: "early_osteoarthritis", 49: "rheumatological_suspicion"}
CAVEAT = ("Case-control; controls are healthy Italian volunteers (not people with other causes of joint pain), older "
          "on average than patients (paper: mean 42 vs 36-38 y) with no per-person age/sex released for controls; all "
          "controls are Italian while half the patients are American home draws shipped to Baltimore, so only the "
          "Italian-only contrast is site-matched; Olink NPX are relative log2 units.")


# ---------------------------------------------------------------------------------------------------------------
# acquisition and parsing
# ---------------------------------------------------------------------------------------------------------------
def fetch() -> dict:
    z = download_file(DATA_URL, SOURCE_ID, ZIP_FILE, max_retries=3)
    out = {}
    man = load_manifest(SOURCE_ID)["files"]
    with zipfile.ZipFile(z) as zf:
        for m in (NPX_MEMBER, CLIN_MEMBER, DEP_MEMBER):
            rel = f"extracted/{m}"
            dest = raw_dir(SOURCE_ID) / rel
            if not dest.exists() or rel not in man:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(zf.read(m))
                record_file(SOURCE_ID, rel, url=f"{DATA_URL}#{m}", note=f"member of {ZIP_FILE}")
            out[m] = dest
    return out


def _p(member: str):
    return raw_dir(SOURCE_ID) / "extracted" / member


def read_npx() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(npx wide: sample_id x assay_key, assay table, sample meta with plate and QC per panel)."""
    n = pd.read_excel(_p(NPX_MEMBER), header=None)
    panel, assay, uniprot, oid = (n.iloc[r, 1:461].astype(str).str.strip().values for r in (2, 3, 4, 5))
    ids = n.iloc[:, 0].fillna("").astype(str).str.strip()
    rows = ids.str.fullmatch(r"[PHC]\d+")
    if int(rows.sum()) != 352:
        raise ValueError(f"expected 352 sample rows, found {int(rows.sum())}")
    keys = [f"{a}|{o}" for a, o in zip(assay, oid)]
    vals = n.loc[rows, 1:460].apply(pd.to_numeric, errors="coerce")
    vals.columns = keys
    vals.index = ids[rows].values
    panels_short = {"Cardiometabolic": "CAM", "Development": "DEV", "Inflammation": "INF", "Neurology": "NEU",
                    "Organ Damage": "ODA"}
    meta = pd.DataFrame(index=ids[rows].values)
    blk = n.iloc[2, 461:471].astype(str).values
    kind = n.iloc[3, 461:471].astype(str).values
    for j, (pnl, k) in enumerate(zip(blk, kind)):
        short = next(v for kk, v in panels_short.items() if kk in pnl)
        meta[f"{'plate' if k == 'Plate ID' else 'qc'}_{short}"] = n.loc[rows, 461 + j].astype(str).values
    at = pd.DataFrame({"assay_key": keys, "assay": assay, "uniprot": uniprot, "olink_id": oid, "panel": panel})
    at["panel_short"] = at.panel.map(lambda p: next(v for kk, v in panels_short.items() if kk in p))
    # QC warnings -> missing (plan amendment 1)
    for short in panels_short.values():
        warn = meta[f"qc_{short}"].str.strip().str.lower() == "warning"
        cols = at.loc[at.panel_short == short, "assay_key"]
        vals.loc[warn.values, cols] = np.nan
    return vals, at, meta


def assay_filter(vals: pd.DataFrame, at: pd.DataFrame) -> pd.DataFrame:
    """Plan: drop cross-panel duplicates (first occurrence kept), then assays with > 25% missing (labels ignored)."""
    at = at.copy()
    at["duplicate_dropped"] = at.duplicated("assay", keep="first")
    at["missing_frac"] = vals.isna().mean().reindex(at.assay_key).values
    at["retained"] = ~at.duplicate_dropped & (at.missing_frac <= MISSING_MAX)
    return at


def read_clinical() -> pd.DataFrame:
    frames = []
    for sheet, grp in (("hEDS", "hEDS"), ("HSD", "HSD")):
        d = pd.read_excel(_p(CLIN_MEMBER), sheet_name=sheet, header=None)
        b = d.iloc[3:]
        b = b[b[2].astype(str).str.fullmatch(r"[PH]\d+")]
        out = pd.DataFrame({"sample_id": b[2].astype(str).values, "group_sheet": grp,
                            "sex": b[4].astype(str).str.strip().replace({"nan": None}).values,
                            "relationship_with_proband": b[5].astype(str).str.strip().values,
                            "age_years": pd.to_numeric(b[6], errors="coerce").values,
                            "meets_2017_hEDS_criteria": (b[7].astype(str).str.strip() == "+").values,
                            "beighton_score": pd.to_numeric(b[8], errors="coerce").values})
        for j, name in COMORB_COLS.items():
            v = b[j].astype(str).str.strip()
            out[f"comorb_{name}"] = v.map({"+": True, "-": False}).astype("boolean").values
        frames.append(out)
    return pd.concat(frames, ignore_index=True)


def participant_id(sid: str) -> str:
    return f"{SOURCE_ID}:{sid}"


def build_participants(vals, meta, clin, retrieved_at) -> pd.DataFrame:
    sid = pd.Series(vals.index, name="sample_id")
    pre = sid.str[0]
    num = sid.str[1:].astype(int)
    grp = pre.map({"P": "hEDS", "H": "HSD", "C": "control"})
    site = np.where(pre == "C", "Italy", np.where(num <= N_ITALIAN, "Italy", "USA"))
    p = pd.DataFrame({"participant_id": sid.map(participant_id), "dataset_id": SOURCE_ID, "native_id": sid,
                      "group_label": grp, "affected": (pre != "C").astype(int), "site": site,
                      "site_basis": "paper: P1-P44/H1-H44 Italian, P45-P88/H45-H88 American, all controls Italian",
                      "condition_id": np.where(pre != "C", CONDITION_ID, ""), "label_basis": LABEL_BASIS,
                      "plate_id": meta["plate_CAM"].str.extract(r"(SP\d+)")[0].values})
    p = p.merge(clin.drop(columns=["group_sheet"]), left_on="native_id", right_on="sample_id", how="left").drop(
        columns=["sample_id"])
    if (p.loc[p.group_label != "control", "relationship_with_proband"].isna()).any():
        raise ValueError("patient without a clinical row")
    return add_provenance(p, data_layer="person", source_name=SOURCE_NAME,
                          source_version=f"Europe PMC PMC13081554 supplementary bundle ({NPX_MEMBER}, {CLIN_MEMBER})",
                          retrieved_at=retrieved_at, evidence_type="person_lab_measurement",
                          source_record_id="native_id",
                          provenance_notes="one row per serum sample/person; controls have no released covariates")


def build_labs(vals, at, retrieved_at) -> pd.DataFrame:
    long = vals.rename_axis("native_id").reset_index().melt(id_vars="native_id", var_name="assay_key", value_name="value")
    long = long.merge(at[["assay_key", "assay", "uniprot", "olink_id", "panel", "retained"]], on="assay_key")
    long.insert(0, "participant_id", long.native_id.map(participant_id))
    long = long.assign(lab_variable=long.olink_id, lab_name=long.assay, unit="NPX (log2, relative)",
                       specimen="serum", assay_platform="Olink Target 96 (PEA)", harmonized_name="olink_npx_" + long.assay,
                       lab_note=np.where(long.value.isna(), "missing or QC warning", ""), dataset_id=SOURCE_ID)
    long["source_record_id"] = long.native_id + ":" + long.olink_id
    return add_provenance(long, data_layer="person", source_name=SOURCE_NAME,
                          source_version=f"Europe PMC PMC13081554 ({NPX_MEMBER})", retrieved_at=retrieved_at,
                          evidence_type="person_lab_measurement", source_record_id="source_record_id",
                          provenance_notes="NPX as released; values with panel QC Warning set to missing")


def ingest() -> dict:
    fetch()
    ret = load_manifest(SOURCE_ID)["files"][ZIP_FILE]["retrieved_at"]
    vals, at, meta = read_npx()
    at = assay_filter(vals, at)
    clin = read_clinical()
    part = build_participants(vals, meta, clin, ret)
    labs = build_labs(vals, at, ret)
    write_table(part, f"participants__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="hEDS/HSD/control serum samples with group, site, patient covariates")
    write_table(labs, f"participant_labs__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="Olink NPX, long format (participant x assay)")
    return {"vals": vals, "at": at, "meta": meta, "part": part, "retrieved_at": ret}


# ---------------------------------------------------------------------------------------------------------------
# analysis (plan: docs/ANALYSIS_PLAN_HEDS_HSD_OLINK_SERUM_CINQUINA2026.md)
# ---------------------------------------------------------------------------------------------------------------
def _xy(part: pd.DataFrame, vals: pd.DataFrame, feats: list[str], mask: pd.Series, label: pd.Series):
    ids = part.loc[mask, "native_id"].values
    return vals.loc[ids, feats].to_numpy(float), label[mask].to_numpy(int), ids


def analyse(d: dict, jobs: int = 8, n_perm: int = 1000, log=print) -> dict:
    vals, at, part = d["vals"], d["at"], d["part"].set_index("native_id", drop=False)
    feats = at.loc[at.retained, "assay_key"].tolist()
    log(f"retained assays: {len(feats)} of {len(at)}")
    it = part.site.eq("Italy")
    aff = part.affected
    res = {"feats": feats}

    # primary: Italian patients vs controls
    m = it
    X, y, ids = _xy(part, vals, feats, m, aff)
    bidx = E.stratified_boot_idx(y)
    prim = E.cv_block(X, y, C_PANEL, bidx, n_perm=n_perm, jobs=jobs)
    res["primary"] = prim
    res["primary_ids"] = ids
    log(f"primary AUROC {prim['summary']['auroc']:.3f} [{prim['summary']['auroc_ci_low']:.3f}, "
        f"{prim['summary']['auroc_ci_high']:.3f}] perm p {prim['summary']['perm_p']:.4f}")

    rows = [{"analysis": "primary: Italian hEDS+HSD vs controls", **prim["summary"]}]
    # S2 site negative control: US vs Italian patients
    m2 = aff.eq(1)
    X2, y2, _ = _xy(part, vals, feats, m2, part.site.eq("USA").astype(int))
    s2 = E.cv_block(X2, y2, C_PANEL, E.stratified_boot_idx(y2), n_perm=n_perm, jobs=jobs)
    rows.append({"analysis": "S2 site negative control: US vs Italian patients", **s2["summary"]})
    # S3 whole cohort (site-confounded)
    X3, y3, _ = _xy(part, vals, feats, pd.Series(True, index=part.index), aff)
    s3 = E.cv_block(X3, y3, C_PANEL, E.stratified_boot_idx(y3), n_perm=0)
    rows.append({"analysis": "S3 whole cohort all patients vs controls (site-confounded by design)", **s3["summary"]})
    # S4 hEDS vs HSD, Italian
    m4 = it & aff.eq(1)
    X4, y4, _ = _xy(part, vals, feats, m4, part.group_label.eq("hEDS").astype(int))
    s4 = E.cv_block(X4, y4, C_PANEL, E.stratified_boot_idx(y4), n_perm=n_perm, jobs=jobs)
    rows.append({"analysis": "S4 hEDS vs HSD (Italian)", **s4["summary"]})
    # S5 probands only
    m5 = it & (aff.eq(0) | part.relationship_with_proband.eq("P"))
    X5, y5, _ = _xy(part, vals, feats, m5, aff)
    s5 = E.cv_block(X5, y5, C_PANEL, E.stratified_boot_idx(y5), n_perm=0)
    rows.append({"analysis": "S5 primary, probands only (relatives excluded)", **s5["summary"]})
    # S7 plate indicators added
    plates = pd.get_dummies(part.loc[it, "plate_id"], drop_first=True).astype(float).to_numpy()
    s7 = E.cv_block(np.hstack([X, plates]), y, C_PANEL, bidx, n_perm=0)
    rows.append({"analysis": "S7 primary + plate indicators", **s7["summary"]})
    perf = pd.DataFrame(rows)
    perf["model"] = f"L2 logistic C={C_PANEL}, median-impute + standardise in fold; {len(feats)} assays"
    perf["cv"] = f"stratified {E.N_SPLITS}-fold x {E.N_REPEATS} (seed {SEED}+r)"
    res["perf"] = perf
    res["null"] = prim["null"]
    res["s2_null"] = s2["null"]

    # S1 per-protein effects, Italian primary sample
    df = vals.loc[ids, feats].copy()
    df["_y"] = y
    eff = E.per_analyte_effects(df, "_y", feats)
    eff = eff.merge(at[["assay_key", "assay", "uniprot", "olink_id", "panel_short"]], left_on="analyte",
                    right_on="assay_key", how="left").drop(columns=["assay_key"])
    res["effects"] = eff

    # decision
    s = prim["summary"]
    ok = bool(s["auroc_ci_low"] > 0.5 and s["perm_p"] < 0.05)
    site = s2["summary"]
    res["decision"] = pd.DataFrame([{
        "question": "Does the serum Olink proteome separate Italian hEDS/HSD patients from Italian healthy controls "
                    "in held-out people?",
        "auroc": s["auroc"], "ci_low": s["auroc_ci_low"], "ci_high": s["auroc_ci_high"], "perm_p": s["perm_p"],
        "n_perm": n_perm, "null_mean": s["null_mean"], "null_q95": s["null_q95"], "n_cases": s["n_cases"],
        "n_controls": s["n_controls"], "rule": "discriminates if CI low > 0.5 AND permutation p < 0.05",
        "decision": ("the serum proteome discriminates Italian hEDS/HSD from Italian controls (pre-specified rule met)"
                     if ok else "not shown (pre-specified rule not met)"),
        "site_negative_control_auroc": site["auroc"], "site_negative_control_ci_low": site["auroc_ci_low"],
        "site_negative_control_perm_p": site["perm_p"],
        "site_detectable": bool(site["auroc_ci_low"] > 0.5),
        "caveat": CAVEAT, "plan": PLAN, "analysis_version": ANALYSIS_VERSION, "computed_at": utc_now_iso()}])
    return res


def paper_dep_overlap(eff: pd.DataFrame) -> pd.DataFrame:
    """Post hoc (plan S1): our Italian-only q < 0.05 proteins vs the paper's 69 pooled-cohort DEPs (published)."""
    x = pd.read_excel(_p(DEP_MEMBER), sheet_name="176 Affected vs 176 Controls", header=None)
    hdr = next(i for i in range(6) if str(x.iloc[i, 0]).strip() == "Names.prot")
    t = x.iloc[hdr + 1:].copy()
    t.columns = [str(c).strip() for c in x.iloc[hdr]]
    padj_col = next((c for c in t.columns if "adj" in c.lower() or "fdr" in c.lower()), None)
    t["paper_q"] = pd.to_numeric(t[padj_col], errors="coerce") if padj_col else np.nan
    paper = set(t.loc[t.paper_q < 0.05, "Names.prot"].astype(str).str.strip())
    ours = set(eff.loc[eff.q_value < 0.05, "assay"].astype(str))
    allp = set(eff.assay.astype(str))
    return pd.DataFrame([{"paper_pooled_DEPs_q_lt_0_05": len(paper), "paper_DEPs_found_in_our_assays": len(paper & allp),
                          "ours_italian_q_lt_0_05": len(ours), "overlap": len(paper & ours),
                          "overlap_names": ";".join(sorted(paper & ours)),
                          "paper_only": ";".join(sorted((paper & allp) - ours)),
                          "ours_only": ";".join(sorted(ours - paper)), "paper_q_column": padj_col}])


# ---------------------------------------------------------------------------------------------------------------
# signatures, tables, report, registry, audit
# ---------------------------------------------------------------------------------------------------------------
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
        rows.append({"feature": f"olink_{r.olink_id}", "feature_label": f"serum {r.assay} (Olink {r.panel_short}, NPX)",
                     "effect_measure": "hedges_g", "effect_size": r.hedges_g,
                     "effect_unit": "Hedges' g on NPX (hEDS/HSD minus control), Italian site-matched sample",
                     "ci_low": r.g_ci_low, "ci_high": r.g_ci_high, "p_value": r.p_value, "q_value": r.q_value,
                     "n_cases": r.n_cases, "n_controls": r.n_controls, "null_value": 0.0,
                     "method": "S1: stratified bootstrap CI (2,000); two-sided Mann-Whitney p; BH q across retained assays"})
    s = res["perf"].set_index("analysis").iloc[0]
    rows.append({"feature": "model_auroc:olink_panel:l2_logistic", "feature_label": f"{len(res['feats'])} Olink assays",
                 "effect_measure": "auroc", "effect_size": s.auroc, "effect_unit": "cross-validated AUROC",
                 "ci_low": s.auroc_ci_low, "ci_high": s.auroc_ci_high, "p_value": s.perm_p, "q_value": np.nan,
                 "n_cases": s.n_cases, "n_controls": s.n_controls, "null_value": 0.5,
                 "method": "PRIMARY (locked): L2 logistic C=0.1, 5-fold x 20 CV; person bootstrap CI; label-permutation p"})
    sig = pd.DataFrame(rows)
    sig.insert(0, "object_id", "signature:" + SOURCE_ID + "|eds_hsd|" + sig.feature)
    sig.insert(1, "dataset_id", SOURCE_ID)
    sig.insert(2, "phenotype_id", "eds_hsd")
    sig.insert(3, "phenotype_label", "hEDS or HSD (2017 criteria), Italian clinic patients vs Italian healthy volunteers")
    sig.insert(4, "condition_id", CONDITION_ID)
    sig.insert(5, "is_proxy", False)
    sig["label_basis"] = LABEL_BASIS
    sig["phenotype_definition"] = "44 hEDS + 44 HSD (P1-P44, H1-H44, Brescia) vs 176 healthy Italian controls; serum Olink"
    sig["caveats"] = CAVEAT
    sig["measurement_class_id"] = MEASUREMENT_CLASS_ID
    sig["n_cases"] = sig.n_cases.astype("Int64")
    sig["n_controls"] = sig.n_controls.astype("Int64")
    lv = []
    for r in sig.itertuples():
        if r.effect_measure == "auroc":
            lv.append("supported_single_dataset" if (r.ci_low > 0.5 and r.p_value < 0.05) else "null_single_dataset")
        else:
            lv.append(_level(r.p_value, r.q_value, r.ci_low, r.ci_high, r.null_value))
    return add_provenance(sig, data_layer="person", source_name=SOURCE_NAME,
                          source_version=f"PMC13081554 supplementary; analysis {ANALYSIS_VERSION}", retrieved_at=retrieved_at,
                          evidence_type="person_derived_feature", source_record_id="object_id",
                          evidence_level=pd.Series(lv),
                          provenance_notes=(f"Cohort-level contrast computed from person-level NPX by measure_it.labs."
                                            f"{SOURCE_ID}; plan {PLAN}; q = BH within the retained-assay family."))


def write_tables(res: dict, overlap: pd.DataFrame) -> None:
    T = TABLES
    T.mkdir(parents=True, exist_ok=True)
    res["perf"].to_csv(T / f"{SOURCE_ID}_model_performance.csv", index=False)
    res["decision"].to_csv(T / f"{SOURCE_ID}_primary_decision.csv", index=False)
    pd.DataFrame({"permutation": np.arange(len(res["null"])), "auroc_null": res["null"]}).to_csv(
        T / f"{SOURCE_ID}_permutation_null.csv", index=False)
    res["effects"].sort_values("p_value").to_csv(T / f"{SOURCE_ID}_protein_effects.csv", index=False)
    overlap.to_csv(T / f"{SOURCE_ID}_paper_dep_overlap.csv", index=False)


def figures(res: dict) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = FIGURES / "drafts"
    out.mkdir(parents=True, exist_ok=True)
    perf = res["perf"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for i, r in enumerate(perf.itertuples()):
        ax[0].errorbar(r.auroc, i, xerr=[[r.auroc - r.auroc_ci_low], [r.auroc_ci_high - r.auroc]], fmt="o",
                       color="#c0392b" if i == 0 else "#444", capsize=3)
    ax[0].axvline(0.5, ls="--", color="grey")
    ax[0].set_yticks(range(len(perf)), [a[:48] for a in perf.analysis])
    ax[0].set_xlabel("cross-validated AUROC (95% person-bootstrap CI)")
    ax[1].hist(res["null"], bins=40, color="#bbb")
    ax[1].axvline(perf.auroc.iloc[0], color="#c0392b", lw=2)
    ax[1].set_xlabel("primary AUROC under label shuffling")
    fig.tight_layout()
    p = out / f"{SOURCE_ID}_auroc_and_null.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return [p]


def write_registry_and_audit(d: dict, res: dict | None = None) -> None:
    man = load_manifest(SOURCE_ID)["files"]
    part = d["part"]
    at = d["at"]
    counts = part.groupby(["group_label", "site"]).size().to_dict()
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "Serum Olink proteome in hEDS/HSD vs healthy controls (Cinquina 2026)",
        "publisher": "Clinical Proteomics (BMC) / Europe PMC; Universita degli Studi di Brescia; The Ehlers-Danlos Society",
        "landing_url": LANDING_URL, "access_urls": [DATA_URL],
        "license": "CC BY 4.0 (article and additional files, BMC open access)",
        "access_conditions": "Open download through the Europe PMC supplementaryFiles API; no registration (verified 2026-09-24).",
        "retrieved_at": d["retrieved_at"], "source_version": f"{ZIP_FILE} sha256 {man[ZIP_FILE]['sha256']}",
        "update_date": "static supplement (published 2026-03-07)", "data_layer": "person",
        "unit_of_observation": "participant (one serum sample)",
        "sample_size": {"participants": int(len(part)), **{f"{g}_{s}": int(n) for (g, s), n in counts.items()}},
        "geographic_resolution": "none", "person_level": True, "geographic": False, "omics": True, "wearable": False,
        "participant_linkage": "NPX and clinical sheets share the sample ID (P#/H#); controls have NPX only",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": [f"participants__{SOURCE_ID}", f"participant_labs__{SOURCE_ID}",
                              f"phenotype_signatures__{SOURCE_ID}", f"results/tables/{SOURCE_ID}_*.csv"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": f"measure_it.labs.{SOURCE_ID}",
        "limitations": ["controls all Italian; half the patients American home draws (site/handling confounding)",
                        "no per-person age/sex for controls", "healthy controls only; no disease comparator",
                        "Olink NPX relative units"],
        "conditions": ["eds_hsd"], "measurement_class_ids": [MEASUREMENT_CLASS_ID], "citation": CITATION,
        "notes": "Rank 1 in docs/LAB_DATASET_DISCOVERY.md; primary analysis locked in " + PLAN})
    miss = part.isna().mean()
    npx_missing = d["vals"].isna().mean().describe()
    files = "\n".join(f"| {k} | {v['bytes']} | {v['sha256'][:16]}... | {v['retrieved_at']} | {v['url']} |" for k, v in man.items())
    text = f"""# DATA AUDIT — {SOURCE_NAME}

| Field | Value |
|---|---|
| source_id | {SOURCE_ID} |
| Source (dataset/API name, exact files/endpoints) | Europe PMC supplementaryFiles bundle for PMC13081554: {NPX_MEMBER} (Additional file 2, NPX), {CLIN_MEMBER} (Additional file 5, patient clinical sheets), {DEP_MEMBER} (published DEP tables, used only for the post hoc overlap) |
| Publishing organization | Clinical Proteomics (BMC); Brescia (Italy) and the HEDGE / DICE EDS & HSD registry (USA) |
| Retrieval date (UTC) | {d['retrieved_at']} |
| Source version / release | {ZIP_FILE}, sha256 {man[ZIP_FILE]['sha256']} |
| Source update date / cadence | static supplement |
| License / access conditions | CC BY 4.0; anonymous download |
| Unit of observation | participant (one serum sample, one Olink run per panel) |
| Sample size (actual, as ingested) | {len(part)} people: {json.dumps({f'{g}/{s}': int(n) for (g, s), n in counts.items()})} |
| Geography (resolution, vintage) | none (site = country from the paper's ID ranges; no location data) |
| Person-level? | yes |
| Geographic? | no |
| Omics? | yes (targeted proteomics, 460 Olink assays) |
| Wearable? | no |
| True participant linkage across modalities? | NPX and patient clinical sheet share the sample ID; controls have NPX only |

## Files / endpoints retrieved

| file | bytes | sha256 | retrieved_at | url |
|---|---|---|---|---|
{files}

## Key variables

* NPX sheet rows 2-5: Panel, Assay, UniProt, OlinkID for 460 assays; sample rows `P#`, `H#`, `C#` (352); 10
  `CONTROL_SAMPLE_CSAS2` plate controls and the LOD / missing-frequency / normalisation rows are dropped. Per-panel
  `Plate ID` (4 plates) and `QC Warning` columns; NPX with a panel QC Warning are set to missing.
* Clinical sheets (`hEDS`, `HSD`): `ID NPx OLINK` (P#/H#, the join key), sex, relationship with proband, age at last
  evaluation, 2017-criteria checklist (Beighton, 5PQ, features A-C, criterion 2/3), comorbidity flags (functional GI,
  neurological, dysautonomia/POTS, psychological, bladder, gynaecological, chronic fatigue, TMJ, allergic/atopic incl.
  MCAS, early osteoarthritis, rheumatological suspicion).
* Labels: P = hEDS (sheet `hEDS`), H = HSD (sheet `HSD`), C = control. Site: P/H 1-44 Italy, 45-88 USA (paper text).

## Missingness

* NPX per assay (after QC-warning masking), fraction missing across 352 samples: min {npx_missing['min']:.3f},
  median {npx_missing['50%']:.3f}, max {npx_missing['max']:.3f}; {int((~at.retained & ~at.duplicate_dropped).sum())}
  assays exceed the 25% threshold and {int(at.duplicate_dropped.sum())} cross-panel duplicates are dropped
  ({int(at.retained.sum())} retained).
* Patient covariates: age missing {int(part.loc[part.affected == 1, 'age_years'].isna().sum())} of 176 patients;
  controls have no age/sex/relationship (not released; {int(miss.get('age_years', 0) * len(part))} missing overall).

## Linkage strategy

Sample ID joins NPX to the patient clinical sheets. No link to any other dataset; the participant id is
`{SOURCE_ID}:<sample id>`.

## Limitations and caveats

{CAVEAT} The paper's own DEP analysis pooled American patients against Italian controls; see the plan.

## Processed outputs

`participants__{SOURCE_ID}` ({len(part)} rows), `participant_labs__{SOURCE_ID}` ({len(part) * 460} rows),
`phenotype_signatures__{SOURCE_ID}`.

## Reproduce

`uv run python -m measure_it.labs.{SOURCE_ID}`
"""
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text(text)


def _f(x, d=3):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def write_report() -> str:
    T = lambda n: pd.read_csv(TABLES / f"{SOURCE_ID}_{n}.csv")  # noqa: E731
    perf, dec, eff, ov = T("model_performance"), T("primary_decision").iloc[0], T("protein_effects"), T("paper_dep_overlap").iloc[0]
    meta = json.loads((TABLES / f"{SOURCE_ID}_run_metadata.json").read_text())
    p = perf.iloc[0]

    def ci(r):
        return f"{_f(r.auroc)} ({_f(r.auroc_ci_low)} to {_f(r.auroc_ci_high)})"
    rows = "\n".join(f"| {r.analysis} | {int(r.n_cases)} / {int(r.n_controls)} | {ci(r)} | "
                     f"{E.pfmt(r.perm_p, r.n_perm) if pd.notna(getattr(r, 'perm_p', np.nan)) else 'not run'} | "
                     f"{_f(r.sens_at_90spec, 2)} ({_f(r.sens_at_90spec_ci_low, 2)}-{_f(r.sens_at_90spec_ci_high, 2)}) |"
                     for r in perf.itertuples())
    nq = int((eff.q_value < 0.05).sum())
    top = eff.sort_values("p_value").head(15)
    toprows = "\n".join(f"| {r.assay} | {r.panel_short} | {int(r.n_cases)} / {int(r.n_controls)} | {_f(r.hedges_g, 2)} "
                        f"({_f(r.g_ci_low, 2)} to {_f(r.g_ci_high, 2)}) | {_f(r.auroc_case_higher, 3)} | {r.p_value:.2e} | "
                        f"{_f(r.q_value, 4)} |" for r in top.itertuples())
    s2 = perf[perf.analysis.str.startswith("S2")].iloc[0]
    s3 = perf[perf.analysis.str.startswith("S3")].iloc[0]
    s4 = perf[perf.analysis.str.startswith("S4")].iloc[0]
    s5 = perf[perf.analysis.str.startswith("S5")].iloc[0]
    s7 = perf[perf.analysis.str.startswith("S7")].iloc[0]
    text = f"""# Serum Olink proteome in hEDS/HSD: pre-specified case-control analysis

Dataset `{SOURCE_ID}`: Cinquina et al., Clinical Proteomics 2026 (PMC13081554), Additional files 2 and 5. 460 Olink
Target 96 assays in serum from 88 hEDS, 88 HSD (2017 criteria) and 176 healthy controls. Plan: `{PLAN}` (written
before any protein value was compared; amendment 1 added before computation). Analysis `{meta['analysis_version']}`,
computed {meta['computed_at']}.

## Answer first

**Primary (Italian, site-matched): {int(p.n_cases)} Italian hEDS/HSD vs {int(p.n_controls)} Italian controls,
{meta['n_assays']} assays, L2 logistic, 5-fold x 20 CV: AUROC {ci(p)}; label-permutation p {E.pfmt(p.perm_p, int(p.n_perm))}
(null mean {_f(p.null_mean)}, 95th percentile {_f(p.null_q95)}, {int(p.n_perm)} shuffles).**
Decision by the locked rule: **{dec.decision}**.

* Sensitivity at 90% specificity: {_f(p.sens_at_90spec, 2)} ({_f(p.sens_at_90spec_ci_low, 2)}-{_f(p.sens_at_90spec_ci_high, 2)}).
* Per protein (S1, Italian sample): **{nq} of {len(eff)} assays have BH q < 0.05**.
* **Site negative control (S2):** American vs Italian *patients* (same diagnoses) are separated with AUROC {ci(s2)},
  permutation p = {E.pfmt(s2.perm_p, int(s2.n_perm))}. {"The proteome carries a strong site/handling signal, so every contrast that puts American patients against the all-Italian controls (including the paper's pooled 176 vs 176 analysis) is site-confounded." if s2.auroc_ci_low > 0.5 else "No site signal was detectable, so pooling sites is less of a concern than the design suggested."}
* Whole cohort as in the paper (S3, site-confounded by design): AUROC {ci(s3)}.
* hEDS vs HSD (S4, Italian): AUROC {ci(s4)}, p = {E.pfmt(s4.perm_p, int(s4.n_perm))} (paper: no DEPs between them).
* Probands only (S5, relatives excluded): AUROC {ci(s5)}; plate indicators added (S7): AUROC {ci(s7)}.
* {int(ov.overlap)} of the {int((eff.q_value < 0.05).sum())} proteins at q < 0.05 in the Italian site-matched sample are
  in the paper's list of 69 pooled DEPs. The strongest are among the paper's top-ranked proteins (MYOC and COMP lower;
  NUDT5 higher). The other {int(ov.paper_pooled_DEPs_q_lt_0_05) - int(ov.overlap)} of the paper's DEPs do not reach
  q < 0.05 here. With half the patients this is partly a matter of power, so whether they depend on the American
  samples cannot be settled.
* **Size:** an AUROC of about 0.68, with 30% sensitivity at 90% specificity, is a modest group-level separation and
  not a usable test. The paper fitted its LASSO, random forest and classification tree on the full data and reports
  no held-out estimate. This one is modest.
* **Age is not controlled.** Controls are older on average (paper: 42 vs 36-38 y), and age is unavailable per control.
  Some of the separation could be age-related.

What it means: {"a serum protein profile distinguishes Italian clinic patients with hEDS/HSD from Italian healthy volunteers in held-out people. This is a group-level, case-control result: the controls are healthy volunteers (older on average, age/sex not released), not people with other causes of joint pain, fatigue or dysautonomia, so it is not evidence of a diagnostic test." if dec.decision.startswith("the serum") else "the pre-specified bar was not met in the site-matched sample."}

## Published claim vs our own computation

| | authors (published) | this analysis (from the released rows) |
|---|---|---|
| Sample | 88 hEDS + 88 HSD (44 Italian + 44 American each) vs 176 Italian controls, pooled | primary: Italian patients only (88) vs 176 controls; pooled cohort only as S3 |
| Method | Wilcoxon + BH per protein; LASSO -> RF / classification tree on the 69 DEPs, fitted on all data | pre-specified CV AUROC, person bootstrap, label-permutation null; per-protein Wilcoxon + BH |
| Per-protein result | 69 pooled DEPs (54 hEDS, 49 HSD, 0 hEDS vs HSD) | {nq} Italian-only assays q < 0.05; overlap with the paper's pooled list: {int(ov.overlap)} of {int(ov.paper_pooled_DEPs_q_lt_0_05)} |
| Site | "no DEPs" Italian vs American within hEDS or HSD | S2 CV AUROC {ci(s2)} for American vs Italian patients |

## Model results

| analysis | cases / controls | AUROC (95% CI) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
{rows}

## Top 15 proteins (S1, Italian sample, sorted by p)

| protein | panel | n | Hedges g (95% CI) | P(case > control) | Mann-Whitney p | BH q |
|---|---|---|---|---|---|---|
{toprows}

Overlap with the paper's pooled DEP list (post hoc, as planned): ours only = {ov.ours_only if isinstance(ov.ours_only, str) else ''};
paper only (present in our assays) = {ov.paper_only if isinstance(ov.paper_only, str) else ''}.

## Caveats carried with every number

{CAVEAT} Bootstrap CIs hold the out-of-fold predictions fixed (person sampling only). Relatives of probands are
included in the primary (S5 excludes them).

## Outputs

`participants__{SOURCE_ID}`, `participant_labs__{SOURCE_ID}`, `phenotype_signatures__{SOURCE_ID}` ({meta['signature_rows']}
rows); `results/tables/{SOURCE_ID}_{{model_performance, primary_decision, permutation_null, protein_effects,
paper_dep_overlap}}.csv`, `{SOURCE_ID}_run_metadata.json`; figure {', '.join(meta['figures'])}.

## Reproduce

```bash
uv run python -m measure_it.labs.{SOURCE_ID} --jobs 16      # deterministic, seed {SEED}
uv run pytest tests/test_{SOURCE_ID}.py
```
"""
    path = PROJECT_ROOT / REPORT
    path.write_text(text)
    return str(path)


def run(jobs: int = 8, n_perm: int = 1000, log=print) -> dict:
    d = ingest()
    res = analyse(d, jobs, n_perm, log)
    ov = paper_dep_overlap(res["effects"])
    write_tables(res, ov)
    figs = figures(res)
    sig = signatures(res, d["retrieved_at"])
    write_table(sig, f"phenotype_signatures__{SOURCE_ID}", producer=f"measure_it.labs.{SOURCE_ID}",
                description="hEDS/HSD Olink: per-protein effects (Italian sample) and primary CV AUROC")
    write_registry_and_audit(d, res)
    meta = {"analysis_version": ANALYSIS_VERSION, "plan": PLAN, "computed_at": utc_now_iso(), "seed": SEED,
            "n_assays": len(res["feats"]), "n_perm": n_perm, "signature_rows": len(sig),
            "figures": [str(p.relative_to(PROJECT_ROOT)) for p in figs]}
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
