"""Measured harmonization coverage -> docs/HARMONIZATION_AUDIT.md (+ results/tables/harmonization_coverage.csv).

    uv run python -m measure_it.harmonize.audit                    # NHANES + any ingested user dataset
    uv run python -m measure_it.harmonize.audit --with-synthetic   # also ingests the two SYNTHETIC BYOD examples into
                                                                   # a temporary folder, measures them, removes them
    uv run python -m measure_it.harmonize.audit --verify-codes     # re-checks every curated LOINC code against the
                                                                   # NLM Clinical Tables index (cached HTTP)
Everything in the document is computed from the processed tables at run time (no number is typed in by hand).
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import pandas as pd

from ..config import PROJECT_ROOT, utc_now_iso
from ..store import partitions, read_table, table_exists
from . import schema as S
from .features import feature_catalog

DOC = PROJECT_ROOT / "docs" / "HARMONIZATION_AUDIT.md"
VARS = PROJECT_ROOT / "results" / "tables" / "harmonization_variables.csv"
COVERAGE = PROJECT_ROOT / "results" / "tables" / "harmonization_coverage.csv"
LOINC_SEARCH = "https://clinicaltables.nlm.nih.gov/api/loinc_items/v3/search"


def verify_loinc_codes() -> pd.DataFrame:
    """Every LOINC code in the curated configs vs the NLM Clinical Tables LONG_COMMON_NAME (exact string match)."""
    from .. import http
    rows = []
    for cfg in ("harmonize_nhanes_loinc", "harmonize_nightingale_loinc"):
        for r in S.config(cfg)["mappings"]:
            try:
                res = http.get_json(LOINC_SEARCH, params={"terms": r["loinc"], "sf": "LOINC_NUM",
                                                          "df": "LOINC_NUM,LONG_COMMON_NAME", "maxList": 5})
                hit = [x for x in res[3] if x[0] == r["loinc"]]
                name = hit[0][1] if hit else None
                status = "verified" if name == r["loinc_name"] else ("name differs" if name else "not found")
            except Exception as exc:  # noqa: BLE001 - recorded, not raised
                name, status = None, f"unavailable ({type(exc).__name__})"
            rows.append({"config": cfg, "variable": r.get("variable", r.get("measure")), "loinc": r["loinc"],
                         "config_name": r["loinc_name"], "nlm_name": name, "status": status})
    return pd.DataFrame(rows)


def verify_icd_codes() -> pd.DataFrame:
    if not table_exists("ontology_icd10cm_codes"):
        return pd.DataFrame()
    known = read_table("ontology_icd10cm_codes", columns=["code", "description"])
    d = dict(zip(known["code"], known["description"]))
    rows = []
    for cfg in ("harmonize_nhanes_conditions", "harmonize_self_report_icd10"):
        for r in S.config(cfg)["mappings"]:
            rows.append({"config": cfg, "icd10cm": r["icd10cm"], "label": r["label"],
                         "cms_description": d.get(r["icd10cm"]),
                         "status": "in CMS FY2026 list" if r["icd10cm"] in d else "NOT FOUND"})
    return pd.DataFrame(rows).drop_duplicates(["config", "icd10cm"])


def _num(s):
    return pd.to_numeric(s, errors="coerce").fillna(0)


def nhanes_coverage(ds: str, v: pd.DataFrame) -> dict:
    g = v[v["dataset_id"] == ds].copy()
    for c in ("n_rows", "n_mapped_rows", "n_yes"):
        g[c] = _num(g.get(c))
    out = {}
    sr = g[(g["domain"] == "condition") & (g["source"] == "clinical_features") & (g["status"] != "not_a_condition")]
    out["condition_selfreport_items"] = f"{int((sr['status'] == 'mapped').sum())}/{len(sr)}"
    yes_m, yes_t = sr.loc[sr["status"] == "mapped", "n_yes"].sum(), sr["n_yes"].sum()
    out["condition_selfreport_yes_mapped_pct"] = 100 * yes_m / yes_t if yes_t else None
    rx = g[(g["domain"] == "condition") & (g["source"] != "clinical_features")]
    out["condition_rx_reason_rows"] = int(rx["n_rows"].sum())
    out["condition_rx_reason_mapped_pct"] = (100 * rx["n_mapped_rows"].sum() / rx["n_rows"].sum()
                                             if rx["n_rows"].sum() else None)
    for src, key in (("labs", "lab_core"), ("labs_extended", "lab_extended"), ("clinical_features", "exam")):
        m = g[(g["domain"] == "measurement") & (g["source"] == src) & (g["status"] != "not_a_numeric_result")]
        if not len(m):
            continue
        out[f"{key}_variables"] = f"{m.loc[m['status'] == 'mapped', 'source_variable'].nunique()}/" \
                                  f"{m['source_variable'].nunique()}"
        out[f"{key}_rows_mapped_pct"] = 100 * m["n_mapped_rows"].sum() / m["n_rows"].sum()
        dropped = _num(m.get("n_dropped_unit")).sum()
        out[f"{key}_rows_dropped_unit"] = int(dropped)
    d = g[g["domain"] == "drug"]
    out["drug_names"] = f"{int((d['status'] == 'mapped').sum())}/{len(d)}"
    out["drug_records_mapped_pct"] = 100 * d["n_mapped_rows"].sum() / d["n_rows"].sum() if d["n_rows"].sum() else None
    s = g[g["domain"] == "survey"]
    out["survey_items"] = int(len(s))
    return out


def dataset_summary(ds: str) -> tuple[dict, pd.DataFrame]:
    t = read_table(f"{S.TABLE}__{ds}", columns=["participant_id", "dataset_id", "domain", "vocabulary",
                                                  "concept_code", "concept_label", "value_as_number",
                                                  "value_as_string", "unit", "mapping_method",
                                                  "mapping_confidence", "platform"])
    cat = feature_catalog(t)
    row = {"dataset_id": ds, "rows": int(len(t)), "participants": int(t["participant_id"].nunique())}
    for v, n in t["vocabulary"].astype(str).value_counts().items():
        row[f"rows_{v}"] = int(n)
    for v, n in cat["vocabulary"].value_counts().items():
        row[f"keys_{v}"] = int(n)
    row["methods"] = json.dumps(t["mapping_method"].astype(str).value_counts().to_dict())
    return row, cat


def overlap_table(cats: dict[str, pd.DataFrame], byod: list[str]) -> list[dict]:
    keys = {ds: set(c["feature_key"]) for ds, c in cats.items()}
    voc = {ds: c.set_index("feature_key")["vocabulary"].to_dict() for ds, c in cats.items()}
    rows = []
    nh = keys.get("nhanes", set()) & keys.get("nhanes0306", set())
    rows.append({"comparison": "nhanes AND nhanes0306", "n_keys": len(nh),
                 "by_vocabulary": pd.Series([voc["nhanes"][k] for k in nh]).value_counts().to_dict() if nh else {}})
    for b in byod:
        bk = {k for k in keys[b] if not k.startswith(("DEVICE:", "OMICS:"))}
        for ref in ("nhanes", "nhanes0306"):
            if ref not in keys:
                continue
            ov = bk & keys[ref]
            rows.append({"comparison": f"{b} AND {ref}", "n_keys": len(ov), "n_byod_keys_non_device": len(bk),
                         "by_vocabulary": pd.Series([voc[b][k] for k in ov]).value_counts().to_dict() if ov else {},
                         "keys": sorted(ov)})
        ov3 = bk & nh
        rows.append({"comparison": f"{b} AND nhanes AND nhanes0306", "n_keys": len(ov3),
                     "n_byod_keys_non_device": len(bk), "keys": sorted(ov3)})
    return rows


def _pct(x) -> str:
    return "n/a" if x is None or (isinstance(x, float) and pd.isna(x)) else f"{x:.1f}%"


def write(verify_codes: bool = False, byod_ids: list[str] | None = None, synthetic_note: str = "") -> dict:
    v = pd.read_csv(VARS, dtype=str) if VARS.exists() else pd.DataFrame(columns=["dataset_id"])
    nh_ids = [d for d in ("nhanes", "nhanes0306") if table_exists(f"{S.TABLE}__{d}")]
    if byod_ids is None:
        byod_ids = sorted(p.stem.split("__", 1)[1] for p in partitions(S.TABLE)
                          if p.stem.split("__", 1)[1].startswith("byod_"))
    rows, cats = [], {}
    for ds in nh_ids + byod_ids:
        r, cats[ds] = dataset_summary(ds)
        if ds in nh_ids:
            r.update(nhanes_coverage(ds, v))
        rows.append(r)
    cov = pd.DataFrame(rows)
    COVERAGE.parent.mkdir(parents=True, exist_ok=True)
    cov[cov["dataset_id"].isin(nh_ids)].to_csv(COVERAGE, index=False)     # user datasets stay out of tracked files
    ov = overlap_table(cats, byod_ids)
    lv = verify_loinc_codes() if verify_codes else None
    iv = verify_icd_codes()
    DOC.write_text(markdown(cov, ov, v, lv, iv, nh_ids, byod_ids, synthetic_note))
    return {"coverage": cov, "overlap": ov}


def markdown(cov, ov, v, lv, iv, nh_ids, byod_ids, synthetic_note) -> str:
    L = ["# Harmonization audit (person_concepts)", "",
         f"Generated {utc_now_iso()} by `uv run python -m measure_it.harmonize.audit"
         + (" --with-synthetic" if synthetic_note else "") + (" --verify-codes" if lv is not None else "")
         + "`. Every number below is measured from the processed tables; contract: "
           "[HARMONIZATION_CONTRACT.md](HARMONIZATION_CONTRACT.md). Per-variable detail (aggregate, NHANES): "
           "`results/tables/harmonization_variables.csv`; per-dataset coverage: "
           "`results/tables/harmonization_coverage.csv`.", "",
         "## 1. Rows and participants", "",
         "| dataset | rows | participants | DEMOG | ICD10CM | LOINC | RXNORM | SURVEY | DEVICE | OMICS |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    def _i(x) -> int:
        return 0 if x is None or (isinstance(x, float) and pd.isna(x)) else int(x)
    for r in cov.to_dict("records"):
        L.append(f"| {r['dataset_id']} | {r['rows']:,} | {r['participants']:,} | "
                 + " | ".join(f"{_i(r.get(f'rows_{x}')):,} ({_i(r.get(f'keys_{x}'))} keys)"
                              for x in ("DEMOG", "ICD10CM", "LOINC", "RXNORM", "SURVEY", "DEVICE", "OMICS")) + " |")
    L += ["", "Keys are feature keys at the level `feature_matrix` uses (ICD-10-CM 3-character categories).", "",
          "## 2. Mapping coverage, NHANES", "",
          "| measure | " + " | ".join(nh_ids) + " |", "|---|" + "---|" * len(nh_ids)]
    nhc = cov.set_index("dataset_id")
    items = [("Self-report condition items mapped to an ICD-10-CM category", "condition_selfreport_items", str),
             ("'Yes' answers on mapped items (share of all 'yes' answers)", "condition_selfreport_yes_mapped_pct", _pct),
             ("Prescription reason-for-use ICD-10-CM rows (native)", "condition_rx_reason_rows", lambda x: f"{int(x):,}"),
             ("  of which well-formed ICD-10-CM (kept)", "condition_rx_reason_mapped_pct", _pct),
             ("Core lab variables mapped to LOINC (participant_labs)", "lab_core_variables", str),
             ("Core lab result rows mapped", "lab_core_rows_mapped_pct", _pct),
             ("Core lab rows dropped for a unit outside the config", "lab_core_rows_dropped_unit", lambda x: f"{int(x):,}"),
             ("Extended lab variables mapped (participant_labs_extended)", "lab_extended_variables", str),
             ("Extended lab result rows mapped", "lab_extended_rows_mapped_pct", _pct),
             ("Exam measures mapped (BP, pulse, anthropometry)", "exam_variables", str),
             ("Distinct drug names (components) mapped to an RxNorm ingredient", "drug_names", str),
             ("Medication records (components) mapped", "drug_records_mapped_pct", _pct),
             ("SURVEY items (PHQ-9, CDC Healthy Days, PFQ / DLQ / SLQ)", "survey_items", lambda x: str(int(x)))]
    for label, key, fmt in items:
        vals = []
        for ds in nh_ids:
            x = nhc.loc[ds].get(key) if key in nhc.columns else None
            vals.append("n/a" if x is None or (isinstance(x, float) and pd.isna(x)) else fmt(x))
        L.append(f"| {label} | " + " | ".join(vals) + " |")
    if len(v):
        L += ["", "Unmapped by design (counted, not guessed):", ""]
        g = v[(v["status"] == "unmapped")].copy()
        g["n_rows"] = _num(g["n_rows"])
        for (ds, dom), gg in g.groupby(["dataset_id", "domain"]):
            top = gg.sort_values("n_rows", ascending=False).head(6)
            L.append(f"* {ds} {dom}: {len(gg)} variables / names, {int(gg['n_rows'].sum()):,} rows; largest: "
                     + "; ".join(f"`{r['source_variable'] if dom != 'drug' else r['label']}` "
                                 f"({int(r['n_rows']):,}: {str(r.get('reason') or '')[:70]})"
                                 for r in top.to_dict("records")))
    L += ["", "## 3. User (BYOD) datasets", ""]
    if not byod_ids:
        L.append("No user dataset was harmonized when this audit ran.")
    else:
        if synthetic_note:
            L.append(synthetic_note)
            L.append("")
        for ds in byod_ids:
            r = nhc.loc[ds]
            L.append(f"* `{ds}`: {int(r['rows']):,} rows, {int(r['participants'])} participants; mapping methods "
                     f"{r['methods']}")
    L += ["", "## 4. Overlap of feature keys (what a phenotype built in one dataset can be evaluated on)", "",
          "| comparison | shared keys | by vocabulary | keys (user datasets) |", "|---|---|---|---|"]
    for r in ov:
        ks = r.get("keys")
        L.append(f"| {r['comparison']} | {r['n_keys']}"
                 + (f" of {r['n_byod_keys_non_device']} non-device keys" if "n_byod_keys_non_device" in r else "")
                 + f" | {r.get('by_vocabulary', '')} | "
                 + (", ".join(f"`{k}`" for k in ks[:40]) + (" ..." if len(ks) > 40 else "") if ks else "") + " |")
    L += ["", "## 5. Code verification", ""]
    if lv is not None:
        L.append(f"LOINC (NLM Clinical Tables, exact LONG_COMMON_NAME match): {int((lv['status'] == 'verified').sum())}"
                 f"/{len(lv)} verified; " + "; ".join(f"{r['loinc']} {r['status']}" for r in
                                                   lv[lv["status"] != "verified"].to_dict("records")))
    else:
        L.append("LOINC: not re-checked in this run (`--verify-codes`); every config row carries the LONG_COMMON_NAME "
                 "looked up on 2026-10-07, and tests verify the LOINC check digits.")
    if len(iv):
        L.append(f"ICD-10-CM: {int((iv['status'] != 'NOT FOUND').sum())}/{len(iv)} curated codes found in the CMS "
                 "FY2026 code list (data/processed/ontology_icd10cm_codes).")
    L += ["", "## 6. Known limits", "",
          "* Self-reported conditions map to ICD-10-CM categories at most at confidence `medium`; several NHANES items "
          "(any thyroid problem, any liver condition, cancer, arthritis of unspecified type) have no single category "
          "and stay unmapped.",
          "* Feature keys roll ICD-10-CM up to 3-character categories (contract), which can join unrelated codes: "
          "e.g. `ICD10CM:Z86` holds both Z86.16 (history of COVID-19, self-reported in a user dataset) and Z86.73 "
          "(history of stroke, NHANES self-report). Consumers that need precision use condition_level='code'.",
          "* NHANES top-codes age (80 in 2011-2014, 85 in 2003-2006): their top band is `80+`, wider than the "
          "contract's 10-year bands.",
          "* Prescription reason-for-use codes exist only in NHANES 2011-2014; categories that come only from them "
          "are missing (NaN), not absent, in 2003-2006.",
          "* Environmental-chemical, infectious-serology (categorical) and microbiome laboratory layers are not "
          "mapped to LOINC: no clinical EHR counterpart is needed for the phenotype use case, and guessing codes for "
          "them is not allowed.",
          "* Survey instruments of user datasets (COMPASS-31, DSQ-PEM, FUNCAP27, ...) have no NHANES counterpart except "
          "the weak crosswalks in configs/harmonize_survey_crosswalk.yaml; SURVEY overlap with NHANES is therefore "
          "small by construction.",
          "* NMR-derived LOINC values (Nightingale) and routine clinical-chemistry values share codes after an "
          "explicit unit conversion but are not calibrated against each other.", ""]
    return "\n".join(L)


def run(verify_codes: bool = False, with_synthetic: bool = False, log=print) -> dict:
    if not with_synthetic:
        return write(verify_codes)
    import importlib.util
    from ..byod.ingest import ingest
    from ..byod.remove import remove
    spec = importlib.util.spec_from_file_location("byod_make_example",
                                                  PROJECT_ROOT / "examples" / "byod_synthetic" / "make_example.py")
    mk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mk)
    ids = ["byod_audit_synth_ehr", "byod_audit_synth_maestro"]
    with tempfile.TemporaryDirectory() as tmp:
        a = mk.make(Path(tmp) / "ehr", fast=True, dataset_id="audit_synth_ehr")
        b = mk.make_maestro(Path(tmp) / "maestro", fast=True, dataset_id="audit_synth_maestro")
        try:
            ingest(a, echo=log, rebuild=False)
            ingest(b, echo=log, rebuild=False)
            note = ("The two user datasets below are the SYNTHETIC examples (`examples/byod_synthetic/make_example.py`: "
                    "`ehr` = ME/CFS with EHR labs, diagnoses and medications; `maestro` = MAESTRO-like layers without "
                    "an EHR extract), ingested into a temporary folder for this audit and removed afterwards. They "
                    "show the mechanics, not the coverage of a real study.")
            return write(verify_codes, byod_ids=ids, synthetic_note=note)
        finally:
            for ds in ids:
                try:
                    remove(ds, rescore=False, purge_ledger=True, echo=log)
                except RuntimeError:
                    pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--verify-codes", action="store_true")
    ap.add_argument("--with-synthetic", action="store_true")
    a = ap.parse_args()
    res = run(a.verify_codes, a.with_synthetic)
    print(res["coverage"].to_string())
    print(json.dumps(res["overlap"], indent=1, default=str))
