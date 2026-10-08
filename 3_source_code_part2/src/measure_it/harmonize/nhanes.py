"""NHANES 2011-2014 (`nhanes`) and 2003-2006 (`nhanes0306`) -> person_concepts__<dataset_id>.

    uv run python -m measure_it.harmonize.nhanes            # both cycles
    uv run python -m measure_it.harmonize.nhanes nhanes0306

Inputs (processed tables, unchanged): participants__<ds>, participant_clinical_features__<ds>,
participant_medications__<ds>, participant_labs__<ds>, participant_labs_extended__nhanes (2011-2014 only) and the
prescription reason-for-use rows of participant_conditions__nhanes (2011-2014 only).

  demographic  DEMOG:age_band (10-year bands; NHANES top-codes age at 80 (2011-14) and 85 (2003-06) -> '80+'),
               DEMOG:sex
  condition    ICD10CM: self-report questions via configs/harmonize_nhanes_conditions.yaml (curated_map; yes = 1.0,
               an explicit no = 0.0) + prescription reason-for-use codes (native, 2011-2014)
  measurement  LOINC: labs and exam measures via configs/harmonize_nhanes_loinc.yaml (curated_map; rows whose unit
               differs from the config unit are dropped and counted, never converted)
  drug         RXNORM ingredients via RxNav (rxnav_api; combination products split into their ingredients)
  survey       SURVEY: PHQ-9, CDC Healthy Days, PFQ / DLQ / SLQ items (configs/harmonize_survey_crosswalk.yaml)
Unmapped variables / names are counted in results/tables/harmonization_unmapped.csv and docs/HARMONIZATION_AUDIT.md.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from ..config import PROJECT_ROOT
from ..store import read_table, table_exists
from . import rxnorm
from . import schema as S

DATASETS = {"nhanes": {"topcode_col": "age_topcoded_80", "topcode": 80.0, "label": "NHANES 2011-2014"},
            "nhanes0306": {"topcode_col": "age_topcoded_85", "topcode": 85.0, "label": "NHANES 2003-2006"}}
RESULTS_TABLES = PROJECT_ROOT / "results" / "tables"


def _norm_unit(u) -> str:
    return str(u or "").strip()


def loinc_lookup(cfg: dict | None = None) -> dict:
    """{(source, variable): [config rows]} with source in labs | labs_extended | clinical_features."""
    cfg = cfg or S.config("harmonize_nhanes_loinc")
    out: dict = {}
    for r in cfg["mappings"]:
        srcs = [r["source"]] if r.get("source") else ["labs", "labs_extended"]
        for s in srcs:
            out.setdefault((s, r["variable"]), []).append(r)
    return out


def _pick(rows: list[dict], component_file: str) -> dict | None:
    pref = [r for r in rows if r.get("component_file_prefix") and str(component_file).startswith(r["component_file_prefix"])]
    if pref:
        return pref[0]
    plain = [r for r in rows if not r.get("component_file_prefix")]
    return plain[0] if plain else None


def lab_concepts(labs: pd.DataFrame, ds: str, source: str, lookup: dict) -> tuple[pd.DataFrame, list[dict]]:
    """Long lab rows (participant_id, lab_variable, value, unit, component_file, lab_name) -> LOINC concept rows."""
    labs = labs.copy()
    labs["component_file"] = labs.get("component_file", pd.Series("", index=labs.index)).astype(str)
    keys = labs[["lab_variable", "component_file", "unit", "lab_name"]].drop_duplicates()
    rows, audit = [], []
    for (var, cf), g in keys.groupby(["lab_variable", "component_file"], sort=True):
        cand = lookup.get((source, var), [])
        r = _pick(cand, cf) if cand else None
        sel = (labs["lab_variable"] == var) & (labs["component_file"] == cf)
        sub = labs.loc[sel]
        n = int(len(sub))
        name = str(g["lab_name"].iloc[0])
        if r is None:
            audit.append({"dataset_id": ds, "domain": "measurement", "source": source, "source_variable": var,
                          "component_file": cf, "label": name, "n_rows": n, "status": "unmapped",
                          "reason": "no curated LOINC mapping (configs/harmonize_nhanes_loinc.yaml)"})
            continue
        ok_units = {r["unit"], *r.get("accept_units", [])}
        unit = sub["unit"].map(_norm_unit)
        val = pd.to_numeric(sub["value"], errors="coerce")
        keep = unit.isin(ok_units) & val.notna()
        n_unit = int((~unit.isin(ok_units)).sum())
        n_null = int((unit.isin(ok_units) & val.isna()).sum())
        audit.append({"dataset_id": ds, "domain": "measurement", "source": source, "source_variable": var,
                      "component_file": cf, "label": name, "n_rows": n, "status": "mapped", "code": r["loinc"],
                      "n_mapped_rows": int(keep.sum()), "n_dropped_unit": n_unit, "n_dropped_missing_value": n_null,
                      "confidence": r["confidence"]})
        if keep.any():
            k = sub[keep]
            rows.append(S.concepts_frame(
                participant_id=k["participant_id"], dataset_id=ds, domain="measurement", vocabulary="LOINC",
                concept_code=r["loinc"], concept_label=r["loinc_name"], value_as_number=val[keep].to_numpy(),
                unit=r["unit"], mapping_method="curated_map", mapping_confidence=r["confidence"],
                source_variable=f"{cf}.{var}" if cf else var).assign(_evidence_type="person_lab_measurement"))
    out = pd.concat(rows, ignore_index=True) if rows else S.coerce(pd.DataFrame(columns=S.COLUMNS))
    return out, audit


def exam_concepts(cf: pd.DataFrame, ds: str, lookup: dict) -> tuple[pd.DataFrame, list[dict]]:
    rows, audit = [], []
    for (src, var), cand in sorted(lookup.items()):
        if src != "clinical_features" or var not in cf.columns:
            continue
        r = cand[0]
        v = pd.to_numeric(cf[var], errors="coerce")
        m = v.notna()
        audit.append({"dataset_id": ds, "domain": "measurement", "source": "clinical_features", "source_variable": var,
                      "label": r["loinc_name"], "n_rows": int(m.sum()), "status": "mapped", "code": r["loinc"],
                      "n_mapped_rows": int(m.sum()), "confidence": r["confidence"]})
        if m.any():
            rows.append(S.concepts_frame(
                participant_id=cf.loc[m, "participant_id"], dataset_id=ds, domain="measurement", vocabulary="LOINC",
                concept_code=r["loinc"], concept_label=r["loinc_name"], value_as_number=v[m].to_numpy(),
                unit=r["unit"], mapping_method="curated_map", mapping_confidence=r["confidence"],
                source_variable=f"participant_clinical_features.{var}").assign(
                _evidence_type="person_exam_measurement"))
    out = pd.concat(rows, ignore_index=True) if rows else S.coerce(pd.DataFrame(columns=S.COLUMNS))
    return out, audit


def condition_concepts(cf: pd.DataFrame, ds: str, cfg: dict | None = None) -> tuple[pd.DataFrame, list[dict]]:
    """Self-report columns -> ICD-10-CM rows (1.0 yes, 0.0 explicit no)."""
    cfg = cfg or S.config("harmonize_nhanes_conditions")
    rows, audit = [], []
    for r in cfg["mappings"]:
        col = r["column"]
        if col not in cf.columns:
            continue
        if "value" in r:            # categorical column (arthritis_type): this type vs any other answer / no arthritis
            t = cf[col].astype("string")
            known = t.notna() | (pd.to_numeric(cf.get("dx_arthritis"), errors="coerce") == 0).fillna(False)
            v = pd.Series(np.where((t == r["value"]).fillna(False).to_numpy(bool), 1.0, 0.0),
                          index=cf.index).where(known.to_numpy(bool))
            qual = f"{col}={r['value']}"
        else:
            v = pd.to_numeric(cf[col], errors="coerce").where(lambda s: s.isin([0, 1]))
            qual = col
        m = v.notna()
        audit.append({"dataset_id": ds, "domain": "condition", "source": "clinical_features", "source_variable": qual,
                      "label": r["label"], "n_rows": int(m.sum()), "n_yes": int((v == 1).sum()), "status": "mapped",
                      "code": r["icd10cm"], "n_mapped_rows": int(m.sum()), "confidence": r["confidence"]})
        if m.any():
            rows.append(S.concepts_frame(
                participant_id=cf.loc[m, "participant_id"], dataset_id=ds, domain="condition", vocabulary="ICD10CM",
                concept_code=r["icd10cm"], concept_label=r["label"], value_as_number=v[m].to_numpy(),
                mapping_method="curated_map", mapping_confidence=r["confidence"],
                source_variable=f"participant_clinical_features.{qual}").assign(_evidence_type="person_self_report"))
    for u in cfg.get("unmapped", []):
        col = u["column"]
        if col not in cf.columns:
            continue
        if "value" in u:
            n_yes = int((cf[col].astype("string") == u["value"]).sum())
            n = n_yes
            qual = f"{col}={u['value']}"
        else:
            v = pd.to_numeric(cf[col], errors="coerce")
            n, n_yes, qual = int(v.notna().sum()), int((v == 1).sum()), col
        audit.append({"dataset_id": ds, "domain": "condition", "source": "clinical_features", "source_variable": qual,
                      "label": col, "n_rows": n, "n_yes": n_yes,
                      "status": "not_a_condition" if u.get("kind") == "not_a_condition" else "unmapped",
                      "reason": u["reason"]})
    out = pd.concat(rows, ignore_index=True) if rows else S.coerce(pd.DataFrame(columns=S.COLUMNS))
    return out, audit


def rx_reason_concepts(ds: str) -> tuple[pd.DataFrame, list[dict]]:
    name = f"participant_conditions__{ds}"
    if not table_exists(name):
        return S.coerce(pd.DataFrame(columns=S.COLUMNS)), []
    c = read_table(name)
    if "source_variable" not in c.columns:
        return S.coerce(pd.DataFrame(columns=S.COLUMNS)), []
    r = c[c["source_variable"].astype(str).str.startswith("RXQ_RX.") & c["icd10cm_code"].astype(str).ne("")]
    ok = r["icd10cm_code"].astype(str).str.fullmatch(S.ICD10CM_RE)
    audit = [{"dataset_id": ds, "domain": "condition", "source": "participant_conditions (RXQ_RX reason codes)",
              "source_variable": "RXQ_RX.RXDRSC1-3", "label": "prescription reason-for-use ICD-10-CM codes",
              "n_rows": int(len(r)), "n_yes": int(len(r)), "status": "mapped" if ok.all() else "partly mapped",
              "code": f"{r.loc[ok, 'icd10cm_code'].nunique()} distinct codes (native)",
              "n_mapped_rows": int(ok.sum()), "confidence": "high"}]
    r = r[ok].drop_duplicates(["participant_id", "icd10cm_code"])
    out = S.concepts_frame(participant_id=r["participant_id"], dataset_id=ds, domain="condition", vocabulary="ICD10CM",
                           concept_code=r["icd10cm_code"], concept_label=r["condition_label"], value_as_number=1.0,
                           mapping_method="native", mapping_confidence="high",
                           source_variable=r["source_variable"]).assign(_evidence_type="person_rx_reason_code")
    return out, audit


def drug_concepts(ds: str, echo=print) -> tuple[pd.DataFrame, list[dict]]:
    name = f"participant_medications__{ds}"
    if not table_exists(name):
        return S.coerce(pd.DataFrame(columns=S.COLUMNS)), []
    m = read_table(name)
    m = m[m["rx_generic_name"].notna() & m["rx_generic_name"].astype(str).str.strip().ne("")]
    if "rx_name_status" in m.columns:
        m = m[m["rx_name_status"].astype(str).eq("recorded")]
    long = m[["participant_id", "rx_generic_name"]].assign(
        component=m["rx_generic_name"].map(rxnorm.split_components)).explode("component")
    long = long[long["component"].notna()]
    res = rxnorm.map_names(long["rx_generic_name"].unique(), echo=echo)
    long["status"] = long["component"].map(lambda c: res[c]["status"])
    audit = []
    for comp, g in long.groupby("component"):
        rr = res[comp]
        audit.append({"dataset_id": ds, "domain": "drug", "source": "participant_medications",
                      "source_variable": "rx_generic_name", "label": comp, "n_rows": int(len(g)), "status": rr["status"],
                      "code": rr.get("rxcui"), "n_mapped_rows": int(len(g)) if rr["status"] == "mapped" else 0,
                      "confidence": ("high" if "via" not in rr else "medium") if rr["status"] == "mapped" else None,
                      "reason": rr.get("reason")})
    ok = long[long["status"] == "mapped"].copy()
    ok["rxcui"] = ok["component"].map(lambda c: res[c]["rxcui"])
    ok["iname"] = ok["component"].map(lambda c: res[c]["name"])
    ok["conf"] = ok["component"].map(lambda c: "medium" if "via" in res[c] else "high")
    ok = ok.drop_duplicates(["participant_id", "rxcui"])
    out = S.concepts_frame(participant_id=ok["participant_id"], dataset_id=ds, domain="drug", vocabulary="RXNORM",
                           concept_code=ok["rxcui"], concept_label=ok["iname"], value_as_number=1.0,
                           mapping_method="rxnav_api", mapping_confidence=ok["conf"],
                           source_variable="RXQ_RX.RXDDRUG:" + ok["component"].astype(str)).assign(
        _evidence_type="person_self_report")
    return out, audit


def survey_concepts(cf: pd.DataFrame, ds: str) -> tuple[pd.DataFrame, list[dict]]:
    cfg = S.config("harmonize_survey_crosswalk")
    rows, audit = [], []
    for r in cfg["reference_items"]:
        col = r["column"]
        if ds not in r.get("datasets", [ds]) or col not in cf.columns:
            continue
        v = pd.to_numeric(cf[col], errors="coerce")
        m = v.notna()
        audit.append({"dataset_id": ds, "domain": "survey", "source": "clinical_features", "source_variable": col,
                      "label": r["label"], "n_rows": int(m.sum()), "status": "mapped", "code": r["code"],
                      "n_mapped_rows": int(m.sum()), "confidence": "high"})
        if m.any():
            rows.append(S.concepts_frame(
                participant_id=cf.loc[m, "participant_id"], dataset_id=ds, domain="survey", vocabulary="SURVEY",
                concept_code=r["code"], concept_label=r["label"], value_as_number=v[m].to_numpy(),
                mapping_method="curated_map", mapping_confidence="high",
                source_variable=f"participant_clinical_features.{col}").assign(_evidence_type="person_self_report"))
    out = pd.concat(rows, ignore_index=True) if rows else S.coerce(pd.DataFrame(columns=S.COLUMNS))
    return out, audit


def build(ds: str = "nhanes", *, echo=print, write: bool = True) -> dict:
    if ds not in DATASETS:
        raise ValueError(f"unknown NHANES dataset {ds!r} (one of {sorted(DATASETS)})")
    info = DATASETS[ds]
    P = read_table(f"participants__{ds}")
    cf = read_table(f"participant_clinical_features__{ds}")
    lookup = loinc_lookup()
    top = P[info["topcode_col"]].astype("boolean").fillna(False) if info["topcode_col"] in P.columns else None
    age = pd.to_numeric(P["age_years"], errors="coerce")
    band = S.age_band_from_years(age, info["topcode"])
    if top is not None:     # a top-coded person is 'X+' even if the recorded value is below the decade boundary
        band = band.where(~top.to_numpy(), f"{int(info['topcode'] // 10 * 10)}+")
    parts = [S.demographic_rows(P["participant_id"], ds, band, P["sex"], band_var="participants.age_years",
                                sex_var="participants.sex").assign(_evidence_type="person_self_report")]
    audit: list[dict] = []
    c1, a = condition_concepts(cf, ds)
    parts.append(c1)
    audit += a
    c2, a = rx_reason_concepts(ds)
    parts.append(c2)
    audit += a
    labs = read_table(f"participant_labs__{ds}",
                      columns=["participant_id", "lab_variable", "value", "unit", "component_file", "lab_name"])
    l1, a = lab_concepts(labs, ds, "labs", lookup)
    parts.append(l1)
    audit += a
    if table_exists(f"participant_labs_extended__{ds}"):
        ext = read_table(f"participant_labs_extended__{ds}",
                         columns=["participant_id", "lab_variable", "value", "unit", "component_file", "lab_name",
                                  "analyte_role", "value_type"])
        proc = ext[ext["analyte_role"] == "process_variable"]
        cat = ext[(ext["analyte_role"] != "process_variable") & (ext["value_type"] != "continuous")]
        ext = ext[(ext["analyte_role"] != "process_variable") & (ext["value_type"] == "continuous")]
        l2, a = lab_concepts(ext, ds, "labs_extended", lookup)
        parts.append(l2)
        audit += a
        for label, sub in (("process variable (not a lab result)", proc), ("categorical result", cat)):
            for var, g in sub.groupby("lab_variable"):
                audit.append({"dataset_id": ds, "domain": "measurement", "source": "labs_extended",
                              "source_variable": var, "label": str(g["lab_name"].iloc[0]), "n_rows": int(len(g)),
                              "status": "not_a_numeric_result" if label.startswith("process") else "unmapped",
                              "reason": label})
    l3, a = exam_concepts(cf, ds, lookup)
    parts.append(l3)
    audit += a
    d, a = drug_concepts(ds, echo=echo)
    parts.append(d)
    audit += a
    s, a = survey_concepts(cf, ds)
    parts.append(s)
    audit += a
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    res = {"dataset_id": ds, "rows": int(len(df)), "participants": int(df["participant_id"].nunique()),
           "audit": audit}
    if write:
        name = S.write_concepts(df, ds, source_name=f"{info['label']} processed tables (measure_it.ingestion)",
                                source_version=f"{S.VERSION}; configs harmonize_nhanes_conditions / _loinc / "
                                               "_survey_crosswalk; RxNav REST",
                                notes=("harmonized from NHANES public-use processed tables; self-report conditions "
                                       "mapped to ICD-10-CM categories (curated); labs to LOINC (curated, unit-checked); "
                                       "drugs to RxNorm ingredients (RxNav)"))
        res["table"] = name
        write_audit_rows(ds, audit)
    echo(f"harmonize {ds}: {res['rows']} rows, {res['participants']} participants")
    return res


def write_audit_rows(ds: str, audit: list[dict]) -> None:
    """results/tables/harmonization_unmapped.csv: one row per source variable / name and its status (aggregate)."""
    RESULTS_TABLES.mkdir(parents=True, exist_ok=True)
    path = RESULTS_TABLES / "harmonization_variables.csv"
    new = pd.DataFrame(audit)
    if path.exists():
        old = pd.read_csv(path, dtype=str)
        old = old[old["dataset_id"] != ds]
        new = pd.concat([old, new.astype(str)], ignore_index=True)
    new.to_csv(path, index=False)


def run(datasets: list[str] | None = None, echo=print) -> dict:
    out = {}
    for ds in datasets or list(DATASETS):
        r = build(ds, echo=echo)
        out[ds] = {k: v for k, v in r.items() if k != "audit"}
    return out


if __name__ == "__main__":
    print(json.dumps(run(sys.argv[1:] or None), indent=1))
