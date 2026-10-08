"""A user (BYOD) dataset's partitions -> person_concepts__byod_<id> (called by `measure-it byod ingest`).

    participants__byod_<id>                 DEMOG age_band (supplied band -> 10-year band), sex
    participant_conditions__byod_<id>       ehr_diagnosis rows -> ICD10CM (native, high);
                                            self_report_history rows -> ICD10CM (self_report_map, <= medium);
                                            study-label rows are outcomes and are NOT harmonized
    participant_labs__byod_<id>             LOINC (native; one unit per code, enforced at ingest)
    participant_medications__byod_<id>      RXNORM (native; the owner supplies ingredient CUIs, not re-checked online:
                                            user-supplied codes never leave the machine)
    participant_survey__byod_<id>           SURVEY <INSTRUMENT>:<item> (native)
    participant_wearable_/device_features   DEVICE <device>:<feature> (native)
    participant_omics_linked__byod_<id>     OMICS <layer>:<feature>; for a layer declared
                                            `omics: {<layer>: {platform: nmr_nightingale, units: nightingale_standard}}`
                                            the clinical-chemistry-equivalent measures become LOINC rows instead
                                            (curated_map, configs/harmonize_nightingale_loinc.yaml, explicit unit
                                            factor recorded in `platform`)
Everything stays local (results/byod/<id>/harmonization_variables.csv lists what was and was not mapped).
"""
from __future__ import annotations

import json

import pandas as pd

from ..store import read_table, table_exists
from . import schema as S
from .self_report import instrument_label


def _read(name: str) -> pd.DataFrame | None:
    return read_table(name) if table_exists(name) else None


def _norm(s: str) -> str:
    return "".join(ch for ch in str(s).lower() if ch.isalnum())


def nightingale_lookup() -> dict:
    cfg = S.config("harmonize_nightingale_loinc")
    out = {}
    for r in cfg["mappings"]:
        for a in [r["measure"], *r.get("aliases", [])]:
            out[_norm(a)] = r
    return out


def build(dataset_id: str, *, manifest: dict | None = None, echo=print, write: bool = True) -> dict:
    from ..byod import common as K
    from ..byod.person import dataset_row
    ds = K.engine_id(dataset_id)
    meta = dataset_row(ds) or {}
    blocks = json.loads(meta.get("blocks_json") or "{}")
    if manifest is None:
        manifest = {}
    omics_spec = manifest.get("omics") or {}
    parts, audit, notes = [], [], []
    # ---- demographics
    P = read_table(K.partition_name("participants", ds))
    band = S.age_band_from_supplied(P["age_range"]) if "age_range" in P.columns else pd.Series([None] * len(P))
    sex = P["sex"] if "sex" in P.columns else pd.Series(["unknown"] * len(P))
    parts.append(S.demographic_rows(P["participant_id"], ds, band, sex, band_var="participants.age_band",
                                    sex_var="participants.sex").assign(_evidence_type="person_self_report"))
    # ---- conditions
    c = _read(K.partition_name("participant_conditions", ds))
    if c is not None and len(c):
        dx = c[(c["condition_type"] == "ehr_diagnosis") & c["icd10cm_code"].astype(str).str.len().gt(0)]
        if len(dx):
            parts.append(S.concepts_frame(
                participant_id=dx["participant_id"], dataset_id=ds, domain="condition", vocabulary="ICD10CM",
                concept_code=dx["icd10cm_code"], concept_label=dx["condition_label"], value_as_number=1.0,
                day_index=dx["day_index"], mapping_method="native", mapping_confidence="high",
                source_variable="ehr_diagnoses.icd10cm").assign(_evidence_type="person_derived_feature"))
            audit.append({"dataset_id": ds, "domain": "condition", "source_variable": "ehr_diagnoses.icd10cm",
                          "n_rows": int(len(dx)), "n_mapped_rows": int(len(dx)), "status": "mapped",
                          "code": f"{dx['icd10cm_code'].nunique()} codes (native)"})
        hx = c[c["condition_type"] == "self_report_history"]
        if len(hx):
            ok = hx[hx["mapping_status"] == "mapped"]
            parts.append(S.concepts_frame(
                participant_id=ok["participant_id"], dataset_id=ds, domain="condition", vocabulary="ICD10CM",
                concept_code=ok["icd10cm_code"], concept_label=ok["condition_label"], value_as_number=1.0,
                day_index=ok["day_index"], mapping_method="self_report_map",
                mapping_confidence=ok["mapping_confidence"].where(ok["mapping_confidence"].isin(["medium", "low"]),
                                                                  "medium"),
                source_variable="self_report_history.condition").assign(_evidence_type="person_self_report"))
            for txt, g in hx.groupby(hx["reported_text"].astype(str).str.lower()):
                st = g["mapping_status"].iloc[0]
                audit.append({"dataset_id": ds, "domain": "condition", "source_variable": "self_report_history",
                              "label": txt if st == "mapped" else "(unmapped text withheld)", "n_rows": int(len(g)),
                              "n_mapped_rows": int(len(g)) if st == "mapped" else 0, "status": st,
                              "code": g["icd10cm_code"].iloc[0] if st == "mapped" else None,
                              "confidence": g["mapping_confidence"].iloc[0] if st == "mapped" else None})
            if (hx["mapping_status"] != "mapped").any():
                notes.append(f"self-reported history: {int((hx['mapping_status'] != 'mapped').sum())} of {len(hx)} "
                             "reports not mapped to ICD-10-CM (kept in participant_conditions, not harmonized)")
    # ---- labs
    lab = _read(K.partition_name("participant_labs", ds))
    if lab is not None and len(lab):
        parts.append(S.concepts_frame(
            participant_id=lab["participant_id"], dataset_id=ds, domain="measurement", vocabulary="LOINC",
            concept_code=lab["loinc_code"], concept_label=lab["lab_name"].where(lab["lab_name"].astype(str) != "",
                                                                                lab["loinc_code"]),
            value_as_number=lab["value"], unit=lab["unit"], day_index=lab["day_index"], mapping_method="native",
            mapping_confidence="high", source_variable="ehr_labs.loinc").assign(
            _evidence_type="person_lab_measurement"))
        audit.append({"dataset_id": ds, "domain": "measurement", "source_variable": "ehr_labs.loinc",
                      "n_rows": int(len(lab)), "n_mapped_rows": int(len(lab)), "status": "mapped",
                      "code": f"{lab['loinc_code'].nunique()} codes (native)"})
    # ---- medications
    med = _read(K.partition_name("participant_medications", ds))
    if med is not None and len(med):
        med = med.drop_duplicates(["participant_id", "rxnorm_cui", "day_index"])
        parts.append(S.concepts_frame(
            participant_id=med["participant_id"], dataset_id=ds, domain="drug", vocabulary="RXNORM",
            concept_code=med["rxnorm_cui"], concept_label="RxNorm " + med["rxnorm_cui"].astype(str),
            value_as_number=1.0, day_index=med["day_index"], mapping_method="native", mapping_confidence="high",
            source_variable="ehr_medications.rxnorm").assign(_evidence_type="person_derived_feature"))
        audit.append({"dataset_id": ds, "domain": "drug", "source_variable": "ehr_medications.rxnorm",
                      "n_rows": int(len(med)), "n_mapped_rows": int(len(med)), "status": "mapped",
                      "code": f"{med['rxnorm_cui'].nunique()} CUIs (native; owner-attested ingredients)"})
    # ---- surveys
    sv = _read(K.partition_name("participant_survey", ds))
    if sv is not None and len(sv):
        lab_ = sv["instrument"].map(instrument_label) + ": " + sv["item"].astype(str)
        parts.append(S.concepts_frame(
            participant_id=sv["participant_id"], dataset_id=ds, domain="survey", vocabulary="SURVEY",
            concept_code=sv["concept_code"], concept_label=lab_, value_as_number=sv["value"],
            day_index=sv["day_index"] if "day_index" in sv.columns else None, mapping_method="native",
            mapping_confidence="high", source_variable=sv["instrument_file"] + ":" + sv["item"].astype(str)).assign(
            _evidence_type="person_self_report"))
        audit.append({"dataset_id": ds, "domain": "survey", "source_variable": "survey_*",
                      "n_rows": int(len(sv)), "n_mapped_rows": int(len(sv)), "status": "mapped",
                      "code": f"{sv['concept_code'].nunique()} items (native)"})
    # ---- devices
    for table in ("participant_wearable_features", "participant_device_features"):
        t = _read(K.partition_name(table, ds))
        if t is None:
            continue
        for b, info in blocks.items():
            if info.get("kind") != "device" or info.get("table") != table:
                continue
            dev = b[len("device_"):]
            cols = [x for x in info.get("features", []) if x in t.columns]
            if not cols:
                continue
            long = t[["participant_id", *cols]].melt(id_vars="participant_id", var_name="col", value_name="v")
            long["v"] = pd.to_numeric(long["v"], errors="coerce")
            long = long[long["v"].notna()]
            code = dev + ":" + long["col"].str.split("__", n=1).str[1]
            parts.append(S.concepts_frame(
                participant_id=long["participant_id"], dataset_id=ds, domain="device", vocabulary="DEVICE",
                concept_code=code, concept_label=code + f" ({info.get('measurement_class')})",
                value_as_number=long["v"], mapping_method="native", mapping_confidence="high",
                source_variable=f"{table}:" + long["col"]).assign(_evidence_type="person_device_measurement"))
    # ---- omics (+ Nightingale -> LOINC)
    om = _read(K.partition_name("participant_omics_linked", ds))
    if om is not None and len(om):
        ng = nightingale_lookup()
        flag = S.config("harmonize_nightingale_loinc")["units_flag"]
        for layer, g in om.groupby("modality"):
            spec = omics_spec.get(layer) or {}
            is_ng = spec.get("platform") == "nmr_nightingale" and spec.get("units") == flag
            hit = g["feature_id"].map(lambda f: ng.get(_norm(f))) if is_ng else pd.Series([None] * len(g), index=g.index)
            mapped = hit.notna()
            rest = g[~mapped]
            parts.append(S.concepts_frame(
                participant_id=rest["participant_id"], dataset_id=ds, domain="omics", vocabulary="OMICS",
                concept_code=f"{layer}:" + rest["feature_id"].astype(str),
                concept_label=f"{layer}:" + rest["feature_id"].astype(str), value_as_number=rest["value"],
                mapping_method="native", mapping_confidence="high",
                platform=spec.get("platform") or None,
                source_variable=f"omics_{layer}:" + rest["feature_id"].astype(str)).assign(
                _evidence_type="person_lab_measurement"))
            if is_ng:
                for feat, gg in g[mapped].groupby("feature_id"):
                    r = ng[_norm(feat)]
                    parts.append(S.concepts_frame(
                        participant_id=gg["participant_id"], dataset_id=ds, domain="measurement", vocabulary="LOINC",
                        concept_code=r["loinc"], concept_label=r["loinc_name"],
                        value_as_number=gg["value"].astype(float) * float(r["factor"]), unit=r["unit"],
                        mapping_method="curated_map", mapping_confidence=r["confidence"],
                        platform=(f"nmr_nightingale; {r['measure']} {r['source_unit']} x {r['factor']} -> {r['unit']}"),
                        source_variable=f"omics_{layer}:{feat}").assign(_evidence_type="person_lab_measurement"))
                    audit.append({"dataset_id": ds, "domain": "measurement", "source_variable": f"omics_{layer}:{feat}",
                                  "label": r["measure"], "n_rows": int(len(gg)), "n_mapped_rows": int(len(gg)),
                                  "status": "mapped", "code": r["loinc"], "confidence": r["confidence"],
                                  "platform": "nmr_nightingale"})
                notes.append(f"omics_{layer} (nmr_nightingale): {g.loc[mapped, 'feature_id'].nunique()} of "
                             f"{g['feature_id'].nunique()} measures mapped to LOINC (explicit unit factors); the rest "
                             "kept as OMICS")
            elif spec.get("platform") == "nmr_nightingale":
                notes.append(f"omics_{layer}: platform nmr_nightingale without units: {flag}; no LOINC mapping")
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    res = {"dataset_id": ds, "rows": int(len(df)), "participants": int(df["participant_id"].nunique()),
           "by_vocabulary": df["vocabulary"].astype(str).value_counts().to_dict(), "notes": notes, "audit": audit}
    if write:
        syn = "SYNTHETIC tutorial data; " if meta.get("synthetic") else ""
        res["table"] = S.write_concepts(
            df, ds, source_name=f"user-supplied dataset '{meta.get('title', ds)}' ({meta.get('owner', '')}) [{ds}]",
            source_version=f"{S.VERSION}; manifest sha256 {str(meta.get('manifest_sha256', ''))[:16]}",
            notes=f"{K.PROVENANCE_NOTE}; local-only; {syn}harmonized vocabulary (docs/HARMONIZATION_CONTRACT.md)",
            description=f"harmonized person concepts of user-supplied dataset {ds} ({K.PROVENANCE_NOTE}; local-only)")
        out = K.local_only_dir(K.results_path(ds))
        pd.DataFrame(audit).to_csv(out / "harmonization_variables.csv", index=False)
    echo(f"  harmonized {ds}: {res['rows']} person_concepts rows ({', '.join(f'{k} {v}' for k, v in sorted(res['by_vocabulary'].items()))})")
    return res
