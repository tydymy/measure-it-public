"""`measure-it byod ingest <dir>`: a validated user folder -> namespaced person-layer partitions with provenance.

Writes (dataset id `byod_<id>`; participant ids `byod_<id>:<native id>` or, with `id_hashing: salted_sha256`,
`byod_<id>:h<first 16 hex of sha256(salt ":" native id)>`, the salt read from MEASURE_IT_BYOD_SALT and never stored):

    participants__byod_<id>                 labels, age band, sex
    participant_conditions__byod_<id>       one row per label value + one per ICD-10-CM diagnosis code
    participant_labs__byod_<id>             LOINC-coded labs
    participant_wearable_features__byod_<id> / participant_device_features__byod_<id>   adapter embeddings
    participant_wearable_daily__byod_<id>   per-day wearable rows (when supplied per day)
    participant_adapter_scores__byod_<id>   adapter phenotype scores (measurements.adapters.scores_to_long)
    participant_omics_linked__byod_<id>     omics layers (long)
    participant_medications__byod_<id>      EHR medications (RxNorm CUIs), when ehr_medications.csv is supplied
    participant_survey__byod_<id>           PRO instruments (long; survey_<instrument>.csv)
    person_concepts__byod_<id>              the harmonized vocabulary (measure_it.harmonize.byod)
    byod_datasets                           one metadata row per user dataset

Every person-layer row carries the provenance columns with data_layer = person, evidence_type person_*,
source_geographic_resolution = none and provenance_notes starting 'user-supplied; not redistributed'. Then the
canonical unions and this dataset's Digital Phenotype Vector are rebuilt (never pooled with another dataset), and
data/raw/byod_<id>/ gets byod_registry_entry.yaml, DATA_AUDIT.md and input_manifest.json (local-only).
"""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from ..config import utc_now_iso
from ..provenance import add_provenance
from ..registry import REQUIRED_FIELDS
from ..store import partitions, read_table, table_exists, write_table
from ..validate import GEO_COLUMN_RE
from . import checks
from . import common as K

PRODUCER = "measure_it.byod.ingest"


_PROV = ("data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
         "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes")


class Refused(RuntimeError):
    """Validation refused the folder; `report` holds every reason."""

    def __init__(self, report):
        super().__init__("\n".join(report.errors))
        self.report = report


def id_mapper(manifest: dict):
    """native participant id -> the token used after the dataset prefix."""
    if manifest.get("id_hashing") == "salted_sha256":
        salt = os.environ.get(K.SALT_ENV, "")
        if len(salt) < 16:
            raise ValueError(f"{K.SALT_ENV} must hold a salt of >= 16 characters")
        return lambda v: "h" + hashlib.sha256(f"{salt}:{v}".encode()).hexdigest()[:16]
    return lambda v: str(v).strip()


class _Prov:
    def __init__(self, ds: str, manifest: dict, msha: str, files: list[str], now: str):
        self.ds, self.m, self.now = ds, manifest, now
        syn = "SYNTHETIC tutorial data (not real people); " if manifest.get("synthetic") else ""
        self.source_name = f"user-supplied dataset '{manifest['title']}' ({manifest['owner']}) [{ds}]"
        self.source_version = f"manifest sha256 {msha[:16]}; files {', '.join(files)}"
        self.notes = f"{K.PROVENANCE_NOTE}; local-only; {syn}licence / data use: {manifest['licence']}"

    def person(self, df: pd.DataFrame, evidence_type: str, record_id, evidence_level: str = "", note: str = ""):
        return add_provenance(df, data_layer="person", source_name=self.source_name,
                              source_version=self.source_version, retrieved_at=self.now, evidence_type=evidence_type,
                              source_record_id=record_id, source_geographic_resolution="none",
                              evidence_level=evidence_level, provenance_notes=self.notes + (f"; {note}" if note else ""))


def _write(df: pd.DataFrame, table: str, ds: str, description: str) -> str:
    geo = [c for c in df.columns if GEO_COLUMN_RE.search(c) and c not in (
        "source_geographic_resolution",)]
    if geo:   # defence in depth: validation already refused these
        raise ValueError(f"{table}: geographic column(s) {geo} in a person-layer table")
    name = K.partition_name(table, ds)
    write_table(df, name, producer=PRODUCER, description=description + f" ({K.PROVENANCE_NOTE}; local-only)")
    return name


def build_frames(rep: checks.Report) -> dict:
    """Validated text frames -> the partition frames (no provenance yet) + metadata. Pure given the report."""
    m = rep.manifest
    ds = K.engine_id(m["dataset_id"])
    tok = id_mapper(m)
    pid = lambda s: ds + ":" + s.astype(str).str.strip().map(tok)  # noqa: E731
    P = rep.frames["participants.csv"].copy()
    labels = {lab["column"]: lab for lab in m["labels"]}
    pa = K.analysis_defaults(m["primary_analysis"])
    primary = labels[pa["label"]]
    out: dict = {"ds": ds, "tables": {}, "blocks": {}, "daily": []}
    # ---- participants
    part = pd.DataFrame({"participant_id": pid(P["participant_id"]), "dataset_id": ds})
    part["native_id"] = part["participant_id"].str.split(":", n=1).str[1]
    if "age_band" in P.columns:
        part["age_range"] = P["age_band"].astype("string").str.strip().to_numpy()
    if "sex" in P.columns:
        part["sex"] = P["sex"].astype(str).str.strip().str.lower().map(checks.SEX_VALUES).fillna("unknown").to_numpy()
    for c in labels:
        part[c] = pd.to_numeric(P[c], errors="coerce").astype("Int64").to_numpy()
    dp_col = ((m.get("subgroup_analysis") or {}).get("device_positive") or {}).get("column")
    if dp_col and dp_col in P.columns:      # owner-supplied device finding (1/0/blank) for `byod subgroup`
        part[dp_col] = pd.to_numeric(P[dp_col], errors="coerce").astype("Int64").to_numpy()
    part["condition_id"] = rep.resolved["labels"][pa["label"]]["condition_id"]
    part["label_basis"] = primary["label_basis"]
    part["user_supplied"] = True
    part["synthetic"] = bool(m.get("synthetic"))
    out["tables"]["participants"] = part
    idmap = dict(zip(P["participant_id"].astype(str).str.strip(), part["participant_id"]))
    # ---- conditions: label rows (+ diagnosis rows below)
    cond_rows = []
    for c, lab in labels.items():
        v = part[["participant_id", "native_id", c]].dropna()
        cond_rows.append(pd.DataFrame({
            "participant_id": v["participant_id"], "dataset_id": ds, "native_id": v["native_id"],
            "condition_type": "study_label", "condition_label": lab["definition"],
            "canonical_condition_id": rep.resolved["labels"][c]["condition_id"], "label_column": c,
            "label_value": v[c].astype(int), "label_basis": lab["label_basis"], "icd10cm_code": "",
            "day_index": pd.NA, "evidence_basis": f"label supplied by the data owner ({lab['label_basis']})"}))
    # ---- data files
    for f in rep.inputs.files:
        if f.kind == "participants":
            continue
        df = rep.frames[f.path.name].copy()
        df["participant_id"] = df["participant_id"].astype(str).str.strip().map(idmap)
        if f.kind == "ehr_labs":
            lab = pd.DataFrame({"participant_id": df["participant_id"], "dataset_id": ds,
                                "native_id": df["participant_id"].str.split(":", n=1).str[1],
                                "loinc_code": df["loinc"].astype(str).str.strip(),
                                "value": pd.to_numeric(df["value"], errors="coerce").astype(float),
                                "unit": df["unit"].astype(str).str.strip()})
            lab["lab_variable"] = "loinc:" + lab["loinc_code"]
            lab["harmonized_name"] = "loinc_" + lab["loinc_code"]
            lab["lab_name"] = df["lab_name"].astype(str) if "lab_name" in df.columns else ""
            lab["day_index"] = pd.to_numeric(df["day_index"], errors="coerce").astype("Int64") \
                if "day_index" in df.columns else pd.array([pd.NA] * len(df), dtype="Int64")
            # one unit per LOINC code: rows in a minority unit are dropped (never converted), and counted
            unit_mode = lab.groupby("loinc_code")["unit"].agg(lambda s: s.value_counts().index[0])
            keep = lab["unit"].eq(lab["loinc_code"].map(unit_mode))
            out["labs_unit_dropped"] = int((~keep).sum())
            lab = lab[keep].reset_index(drop=True)
            lab["record_key"] = lab["participant_id"] + "|" + lab["loinc_code"] + "|" + lab.index.astype(str)
            out["tables"]["participant_labs"] = lab
            out["blocks"]["ehr_labs"] = {"kind": "ehr_labs", "table": "participant_labs",
                                         "features": sorted("loinc_" + lab["loinc_code"].unique())}
        elif f.kind == "ehr_diagnoses":
            codes = df["icd10cm"].astype(str).str.strip().str.upper()
            codes = codes.where(codes.str.contains(r"\.") | (codes.str.len() <= 3),
                                codes.str[:3] + "." + codes.str[3:])
            from .person import icd_condition
            cond_rows.append(pd.DataFrame({
                "participant_id": df["participant_id"], "dataset_id": ds,
                "native_id": df["participant_id"].str.split(":", n=1).str[1], "condition_type": "ehr_diagnosis",
                "condition_label": codes, "canonical_condition_id": codes.map(lambda c: icd_condition(c) or ""),
                "label_column": "", "label_value": 1, "label_basis": "ehr_codes", "icd10cm_code": codes,
                "day_index": pd.to_numeric(df["day_index"], errors="coerce").astype("Int64")
                if "day_index" in df.columns else pd.NA,
                "evidence_basis": "EHR ICD-10-CM diagnosis code (user-supplied)"}))
            out["blocks"]["ehr_diagnoses"] = {"kind": "ehr_diagnoses", "table": "participant_conditions",
                                              "features": sorted(codes.unique())}
        elif f.kind == "ehr_medications":
            med = pd.DataFrame({"participant_id": df["participant_id"], "dataset_id": ds,
                                "native_id": df["participant_id"].str.split(":", n=1).str[1],
                                "rxnorm_cui": df["rxnorm"].astype(str).str.strip(),
                                "day_index": pd.to_numeric(df["day_index"], errors="coerce").astype("Int64")
                                if "day_index" in df.columns else pd.array([pd.NA] * len(df), dtype="Int64"),
                                "evidence_basis": "EHR medication (RxNorm CUI, user-supplied)"})
            med = med.reset_index(drop=True)
            med["record_key"] = med["participant_id"] + "|rx|" + med.index.astype(str)
            out["tables"]["participant_medications"] = med
            out["blocks"]["ehr_medications"] = {"kind": "ehr_medications", "table": "participant_medications",
                                                "features": sorted(med["rxnorm_cui"].unique())}
        elif f.kind == "self_report_history":
            from ..harmonize.self_report import map_history
            from .person import icd_condition
            icd_in = df["icd10cm"] if "icd10cm" in df.columns else pd.Series([None] * len(df), index=df.index)
            maps = [map_history(t, c) for t, c in zip(df["condition"], icd_in)]
            codes = pd.Series([mm["icd10cm"] or "" for mm in maps], index=df.index)
            cond_rows.append(pd.DataFrame({
                "participant_id": df["participant_id"], "dataset_id": ds,
                "native_id": df["participant_id"].str.split(":", n=1).str[1], "condition_type": "self_report_history",
                "condition_label": [mm["label"] for mm in maps],
                "canonical_condition_id": codes.map(lambda c: (icd_condition(c) or "") if c else ""),
                "label_column": "", "label_value": 1, "label_basis": "self_report", "icd10cm_code": codes,
                "day_index": pd.to_numeric(df["day_index"], errors="coerce").astype("Int64")
                if "day_index" in df.columns else pd.NA,
                "evidence_basis": "self-reported medical history (configs/harmonize_self_report_icd10.yaml)",
                "reported_text": df["condition"].astype(str).str.strip().str.slice(0, 80),
                "mapping_status": [mm["status"] for mm in maps],
                "mapping_confidence": [mm["confidence"] or "" for mm in maps]}))
            out["history_unmapped"] = int(sum(mm["status"] != "mapped" for mm in maps))
            out["blocks"]["self_report_history"] = {"kind": "self_report_history", "table": "participant_conditions",
                                                    "features": sorted(c for c in codes.unique() if c)}
        elif f.kind == "survey":
            from ..harmonize.self_report import instrument_code
            code = instrument_code(f.name)
            keys = ["participant_id"] + (["day_index"] if "day_index" in df.columns else [])
            items = [c for c in df.columns if c not in keys]
            long = df.melt(id_vars=keys, value_vars=items, var_name="item", value_name="value")
            long["value"] = pd.to_numeric(long["value"], errors="coerce").astype(float)
            long = long[long["value"].notna()]
            long.insert(1, "dataset_id", ds)
            long["native_id"] = long["participant_id"].str.split(":", n=1).str[1]
            long["instrument"] = code
            long["instrument_file"] = f.path.name
            long["concept_code"] = code + ":" + long["item"].astype(str)
            if "day_index" in long.columns:
                long["day_index"] = pd.to_numeric(long["day_index"], errors="coerce").astype("Int64")
            out.setdefault("survey_long", []).append(long)
            out["blocks"][f.block] = {"kind": "survey", "table": "participant_survey", "instrument": code,
                                      "features": [f"{code}:{c}" for c in items]}
        elif f.kind == "omics":
            feats = [c for c in df.columns if c != "participant_id"]
            long = df.melt(id_vars="participant_id", value_vars=feats, var_name="feature_id", value_name="value")
            long["value"] = pd.to_numeric(long["value"], errors="coerce").astype(float)
            long = long[long["value"].notna()]
            long.insert(1, "dataset_id", ds)
            long["modality"] = f.name
            long["feature_name"] = long["feature_id"]
            long["value_unit"] = "as supplied by the data owner"
            long["linkage_method"] = "owner-supplied participant_id shared with the other files (overlap counted)"
            long["source_file"] = f.path.name
            out.setdefault("omics_long", []).append(long)
            out["blocks"][f.block] = {"kind": "omics", "table": "participant_omics_linked", "features": feats,
                                      "layer": f.name}
    cond = pd.concat(cond_rows, ignore_index=True)
    cond["record_key"] = cond["participant_id"] + "|" + cond["condition_type"] + "|" + cond.index.astype(str)
    out["tables"]["participant_conditions"] = cond
    if out.get("survey_long"):
        sv = pd.concat(out.pop("survey_long"), ignore_index=True, sort=False)
        sv["record_key"] = sv["participant_id"] + "|" + sv["concept_code"] + "|" + sv.index.astype(str)
        out["tables"]["participant_survey"] = sv
    if out.get("omics_long"):
        om = pd.concat(out.pop("omics_long"), ignore_index=True)
        om["record_key"] = om["participant_id"] + "|" + om["modality"] + "|" + om["feature_id"].astype(str)
        out["tables"]["participant_omics_linked"] = om
    return out


def _device_tables(rep: checks.Report, out: dict, echo) -> dict:
    """Run each device file through its adapter; returns {table: wide frame}, daily rows and adapter scores."""
    ds = out["ds"]
    idmap = dict(zip(rep.frames["participants.csv"]["participant_id"].astype(str).str.strip(),
                     out["tables"]["participants"]["participant_id"]))
    classes = {c["id"]: c for c in checks.load_config("measurements")["measurement_classes"]}
    wide: dict[str, pd.DataFrame] = {}
    scores, daily, adapters, registered = [], [], {}, []
    try:
        _run_devices(rep, out, ds, idmap, classes, wide, scores, daily, adapters, registered, echo)
    finally:
        # the adapter registry is per process and is needed only while ingesting: a long-lived process (the test
        # suite, a notebook) must not keep a partner's adapter registered after the ingest
        from ..measurements.adapters import unregister_adapter
        for name in registered:
            unregister_adapter(name)
    return {"wide": wide, "scores": scores, "daily": daily, "adapters": adapters}


def _run_devices(rep, out, ds, idmap, classes, wide, scores, daily, adapters, registered, echo) -> None:
    from ..measurements.adapters import scores_to_long
    from .adapters import adapter_for_device, run_adapter
    m = rep.manifest
    for f in rep.inputs.by_kind("device"):
        raw = rep.frames[f.path.name].copy()
        raw["participant_id"] = raw["participant_id"].astype(str).str.strip().map(idmap)
        adapter, reg_name = adapter_for_device(rep.root, m, f.name, ds)
        registered.append(reg_name)
        emb, sc, desc = run_adapter(adapter, raw)
        cls = desc.measurement_ids[0]
        fam = classes.get(cls, {}).get("modality_family")
        table = "participant_wearable_features" if fam == "wearable" else "participant_device_features"
        feats = [c for c in emb.columns if c != "participant_id"]
        emb = emb.rename(columns={c: f"{f.name}__{c}" for c in feats})
        pref = [f"{f.name}__{c}" for c in feats]
        w = wide.get(table)
        wide[table] = emb if w is None else w.merge(emb, on="participant_id", how="outer")
        out["blocks"][f.block] = {"kind": "device", "table": table, "features": pref, "measurement_class": cls,
                                  "adapter_id": desc.adapter_id, "adapter_registry_name": reg_name,
                                  "adapter_class": type(adapter).__name__, "per_day": "day_index" in raw.columns,
                                  "modality": desc.modality}
        adapters[f.block] = desc.to_dict()
        long = scores_to_long(desc, sc, ds)
        scores.append(long)
        if "day_index" in raw.columns and fam == "wearable":
            d = raw.copy()
            fcols = [c for c in d.columns if c not in ("participant_id", "day_index")]
            for c in fcols:
                d[c] = pd.to_numeric(d[c], errors="coerce").astype(float)
            d = d.rename(columns={c: f"{f.name}__{c}" for c in fcols})
            d.insert(1, "dataset_id", ds)
            d["device"] = f.name
            d["day_index"] = pd.to_numeric(d["day_index"]).astype(int)
            daily.append(d)
        echo(f"  device {f.name}: adapter {desc.adapter_id} ({type(adapter).__name__}), class {cls}, "
             f"{len(emb)} participants x {len(pref)} features -> {table}")


def ingest(root: str | Path, *, min_cases: int | None = None, min_controls: int | None = None, rebuild: bool = True,
           echo=print) -> dict:
    root = Path(root).resolve()
    rep = checks.check_dir(root, min_cases, min_controls)
    if not rep.ok:
        raise Refused(rep)
    m = rep.manifest
    ds = K.engine_id(m["dataset_id"])
    now = utc_now_iso()
    msha = K.sha256_file(rep.inputs.manifest_path)
    files = sorted(f.path.name for f in rep.inputs.files)
    prov = _Prov(ds, m, msha, files, now)
    echo(f"byod ingest {ds}: {len(files)} file(s), manifest sha256 {msha[:12]}")
    # remove partitions of an earlier ingest of the same dataset (re-ingest replaces, never appends); an earlier
    # evaluation described the old data, so its record is discarded too (the plan lock and the ledger are kept)
    for p in K.dataset_partitions(ds):
        p.unlink()
    if table_exists(K.RECORDS_TABLE):
        r = read_table(K.RECORDS_TABLE)
        if (r["dataset_id"] == ds).any():
            from .evaluate import write_records_table
            write_records_table(r[r["dataset_id"] != ds].drop(columns=[c for c in r.columns if c in _PROV]))
            echo("  an earlier evaluation of this dataset was discarded: run `byod evaluate` (and `byod deploy`) again")
    out = build_frames(rep)
    dev = _device_tables(rep, out, echo)
    written = []
    T = out["tables"]
    written.append(_write(prov.person(T["participants"], "person_derived_feature", "native_id",
                                      evidence_level=T["participants"]["label_basis"].iloc[0],
                                      note="labels as supplied by the owner; age as a band"),
                          "participants", ds, "participants of a user-supplied dataset"))
    c = T["participant_conditions"]
    written.append(_write(prov.person(c, "person_derived_feature", "record_key",
                                      evidence_level=np.select(
                                          [c["condition_type"] == "ehr_diagnosis",
                                           c["condition_type"] == "self_report_history"],
                                          ["ehr_icd10cm_code", "self_report_history"], "owner_label"),
                                      note="study labels, EHR ICD-10-CM codes and self-reported history")
                              .drop(columns=["record_key"]),
                          "participant_conditions", ds, "labels, EHR diagnosis codes, self-reported history"))
    if "participant_labs" in T:
        lab = T["participant_labs"]
        written.append(_write(prov.person(lab, "person_lab_measurement", "record_key", note="LOINC-coded EHR labs")
                              .drop(columns=["record_key"]), "participant_labs", ds, "EHR labs (LOINC)"))
    if "participant_omics_linked" in T:
        om = T["participant_omics_linked"]
        written.append(_write(prov.person(om, "person_lab_measurement", "record_key",
                                          evidence_level="owner_declared_participant_id",
                                          note="participant-linked omics (owner-supplied shared participant_id)")
                              .drop(columns=["record_key"]), "participant_omics_linked", ds, "omics layers (long)"))
    if "participant_medications" in T:
        med = T["participant_medications"]
        written.append(_write(prov.person(med, "person_derived_feature", "record_key", note="EHR medications (RxNorm)")
                              .drop(columns=["record_key"]), "participant_medications", ds, "EHR medications (RxNorm)"))
    if "participant_survey" in T:
        sv = T["participant_survey"]
        written.append(_write(prov.person(sv, "person_self_report", "record_key",
                                          note="patient-reported outcome instruments (items / scores as supplied)")
                              .drop(columns=["record_key"]), "participant_survey", ds, "PRO instruments (long)"))
    for table, w in dev["wide"].items():
        w = w.copy()
        w.insert(1, "dataset_id", ds)
        w.insert(2, "native_id", w["participant_id"].str.split(":", n=1).str[1])
        w["device_blocks"] = ";".join(b for b, i in out["blocks"].items() if i.get("table") == table)
        written.append(_write(prov.person(w, "person_device_measurement", "native_id",
                                          note="features from the dataset's measurement adapter(s)"),
                              table, ds, "user-supplied device features"))
    if dev["daily"]:
        d = pd.concat(dev["daily"], ignore_index=True, sort=False)
        d["record_key"] = d["participant_id"] + "|" + d["device"] + "|" + d["day_index"].astype(str)
        written.append(_write(prov.person(d, "person_device_measurement", "record_key", note="per-day rows")
                              .drop(columns=["record_key"]), "participant_wearable_daily", ds, "per-day wearable rows"))
    if dev["scores"]:
        s = pd.concat(dev["scores"], ignore_index=True)
        written.append(_write(prov.person(s, "person_derived_feature", "score_record_id",
                                          note="adapter phenotype scores: descriptive, not diagnoses"),
                              "participant_adapter_scores", ds, "measurement-adapter scores"))
    # ---- metadata row, registry entry, audit
    pa = K.analysis_defaults(m["primary_analysis"])
    gs = rep.resolved.get("group_sizes", {})
    meta = {
        "dataset_id": ds, "manifest_dataset_id": m["dataset_id"], "title": m["title"], "owner": m["owner"],
        "description": m.get("description", ""), "licence": m["licence"],
        "data_use_statement": m["data_use"]["statement"],
        "irb_or_dua_reference": m["data_use"].get("irb_or_dua_reference", ""),
        "synthetic": bool(m.get("synthetic")), "demo": bool(m.get("demo")),
        "comparator_type": m["comparator"]["type"], "comparator_kind": K.COMPARATOR_TYPES[m["comparator"]["type"]],
        "comparator_description": m["comparator"]["description"],
        "labels_json": json.dumps(rep.resolved["labels"]), "primary_label": pa["label"],
        "primary_condition_id": rep.resolved["labels"][pa["label"]]["condition_id"],
        "measurement_class": m["measurement"]["class"], "measurement_bundle": m["measurement"].get("bundle", ""),
        "bundles_containing_class": ";".join(rep.resolved.get("bundles_containing_class", [])),
        "blocks_json": json.dumps(out["blocks"]), "adapters_json": json.dumps(dev["adapters"], default=str),
        "n_participants": int(len(T["participants"])), "n_cases_primary": int(gs.get("cases", 0)),
        "n_controls_primary": int(gs.get("controls", 0)), "id_hashing": m.get("id_hashing", "none"),
        "manifest_sha256": msha, "ingested_at": now, "status": "ingested", "plan_sha256": "", "evaluated_at": "",
        "deployed_at": "", "deployed_with_demo": False, "user_supplied": True, "local_only": True,
        "partitions": ";".join(written)}
    upsert_dataset_row(meta)
    raw = K.local_only_dir(K.raw_path(ds))
    inputs = []
    for f in rep.inputs.files:
        df = rep.frames[f.path.name]
        inputs.append({"file": f.path.name, "kind": f.kind, "block": f.block, "bytes": f.path.stat().st_size,
                       "sha256": K.sha256_file(f.path), "rows": int(len(df)), "columns": int(df.shape[1])})
    (raw / "input_manifest.json").write_text(json.dumps({"dataset_id": ds, "manifest_sha256": msha, "files": inputs,
                                                         "ingested_at": now, "note": "file names, sizes and hashes "
                                                         "only; the data stay where the owner keeps them"}, indent=1))
    # the shared vocabulary (person_concepts__byod_<id>; docs/HARMONIZATION_CONTRACT.md) from the partitions above
    from ..harmonize.byod import build as harmonize_build
    hz = harmonize_build(ds, manifest=m, echo=echo)
    written.append(hz["table"])
    meta["partitions"] = ";".join(written)
    upsert_dataset_row(meta)
    entry = registry_entry(meta, rep, written, inputs)
    (raw / "byod_registry_entry.yaml").write_text(yaml.safe_dump(entry, sort_keys=False, allow_unicode=True))
    (raw / "DATA_AUDIT.md").write_text(data_audit(meta, rep, out, dev, written, inputs, hz))
    echo(f"  wrote {len(written)} partition(s): {', '.join(written)}")
    res = {"dataset_id": ds, "partitions": written, "meta": meta, "warnings": rep.warnings}
    if rebuild:
        res["rebuild"] = rebuild_person_layer(echo=echo)
    return res


def upsert_dataset_row(meta: dict) -> None:
    old = read_table(K.DATASETS_TABLE) if table_exists(K.DATASETS_TABLE) else None
    row = pd.DataFrame([meta])
    if old is not None:
        old = old.drop(columns=[c for c in old.columns if c in (
            "data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
            "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes")])
        row = pd.concat([old[old["dataset_id"] != meta["dataset_id"]], row], ignore_index=True)
    write_datasets_table(row)


def write_datasets_table(df: pd.DataFrame) -> None:
    from ..config import PROCESSED
    if df.empty:
        for suf in (".parquet", ".meta.json"):
            (PROCESSED / f"{K.DATASETS_TABLE}{suf}").unlink(missing_ok=True)
        return
    df = add_provenance(df.reset_index(drop=True), data_layer="metadata",
                        source_name="measure_it.byod (user dataset manifests)",
                        source_version="byod manifest schema v1", retrieved_at=utc_now_iso(),
                        evidence_type="metadata_catalog", source_record_id="dataset_id",
                        provenance_notes=f"one row per user-supplied dataset on this machine; {K.PROVENANCE_NOTE}; "
                                         "local-only metadata (no person rows)")
    write_table(df, K.DATASETS_TABLE, producer=PRODUCER,
                description="User-supplied (BYOD) datasets on this machine: manifest, plan hash, status")


def registry_entry(meta: dict, rep: checks.Report, written: list[str], inputs: list[dict]) -> dict:
    blocks = json.loads(meta["blocks_json"])
    e = {
        "source_id": meta["dataset_id"], "name": meta["title"], "publisher": meta["owner"],
        "landing_url": "not applicable (user-supplied, local-only)", "access_urls": [],
        "license": meta["licence"], "access_conditions": meta["data_use_statement"],
        "retrieved_at": meta["ingested_at"], "source_version": f"manifest sha256 {meta['manifest_sha256'][:16]}",
        "update_date": "as supplied", "data_layer": "person",
        "unit_of_observation": "participant (" + ", ".join(sorted({b.get("kind") for b in blocks.values()})) + ")",
        "sample_size": {"participants": meta["n_participants"], "primary_cases": meta["n_cases_primary"],
                        "primary_controls": meta["n_controls_primary"]},
        "geographic_resolution": "none", "person_level": True, "geographic": False,
        "omics": any(b.get("kind") == "omics" for b in blocks.values()),
        "wearable": any(b.get("table") == "participant_wearable_features" for b in blocks.values()),
        "participant_linkage": f"namespaced ids {meta['dataset_id']}:<native id>; never linked to another dataset",
        "true_participant_linkage_across_modalities": (
            "owner-declared: the same participant_id in every file (overlap counted in DATA_AUDIT.md)"),
        "status": "user_supplied", "processed_outputs": written, "audit": f"data/raw/{meta['dataset_id']}/DATA_AUDIT.md",
        "ingestion_module": PRODUCER,
        "limitations": [K.PROVENANCE_NOTE, "local-only: never merged into SOURCE_REGISTRY.yaml and never sent to any "
                        "external service", "de-identification, labels and permission attested by the data owner"]
                       + (["SYNTHETIC tutorial data"] if meta["synthetic"] else []),
        "user_supplied": True, "local_only": True, "demo": meta["demo"], "synthetic": meta["synthetic"],
        "input_files": inputs,
    }
    missing = [k for k in REQUIRED_FIELDS if k not in e]
    if missing:
        raise ValueError(f"registry entry misses {missing}")
    return e


def _missing_line(df: pd.DataFrame, cols: list[str]) -> str:
    if not cols:
        return "no feature columns"
    miss = df[cols].isna().mean()
    return (f"{len(cols)} features; missing per feature median {miss.median():.1%}, max {miss.max():.1%} "
            f"({miss.idxmax()})")


def data_audit(meta: dict, rep: checks.Report, out: dict, dev: dict, written: list[str], inputs: list[dict],
               harmonized: dict | None = None) -> str:
    T = out["tables"]
    ds = meta["dataset_id"]
    blocks = json.loads(meta["blocks_json"])
    labels = json.loads(meta["labels_json"])
    lines = [f"# DATA AUDIT — {meta['title']} ({ds})", "",
             ("**SYNTHETIC tutorial dataset: generated data, not real people.**\n" if meta["synthetic"] else "")
             + f"User-supplied, local-only ({K.PROVENANCE_NOTE}). Generated by `measure-it byod ingest` at "
               f"{meta['ingested_at']}; nothing here is published or redistributed.", "",
             "| Field | Value |", "|---|---|",
             f"| source_id | {ds} |",
             f"| Source (files) | {', '.join(i['file'] for i in inputs)} |",
             f"| Publishing organization | {meta['owner']} (data owner) |",
             f"| Retrieval date (UTC) | {meta['ingested_at']} (ingestion) |",
             f"| Source version / release | manifest sha256 {meta['manifest_sha256'][:16]} |",
             "| Source update date / cadence | as supplied |",
             f"| License / access conditions | {meta['licence']}; data use: {meta['data_use_statement']} |",
             "| Unit of observation | participant |",
             f"| Sample size (actual, as ingested) | {meta['n_participants']} participants; primary label "
             f"{meta['primary_label']}: {meta['n_cases_primary']} cases, {meta['n_controls_primary']} controls |",
             "| Geography (resolution, vintage) | none (geographic columns refused at validation) |",
             "| Person-level? | yes |", "| Geographic? | no |",
             f"| Omics? | {'yes' if any(b['kind'] == 'omics' for b in blocks.values()) else 'no'} |",
             f"| Wearable? | {'yes' if any(b.get('table') == 'participant_wearable_features' for b in blocks.values()) else 'no'} |",
             "| True participant linkage across modalities? | owner-declared shared participant_id; overlap below |",
             "", "## Files", "", "| file | kind | rows | columns | bytes | sha256 |", "|---|---|---|---|---|---|"]
    lines += [f"| {i['file']} | {i['kind']} | {i['rows']} | {i['columns']} | {i['bytes']} | {i['sha256'][:16]} |"
              for i in inputs]
    lines += ["", "## Labels and comparator", ""]
    for c, v in labels.items():
        y = T["participants"][c]
        lines.append(f"* `{c}` -> condition `{v['condition_id']}` ({v['mapping']}); basis {v['label_basis']}; "
                     f"{int((y == 1).sum())} cases, {int((y == 0).sum())} controls, {int(y.isna().sum())} blank. "
                     f"Definition: {v['definition']}")
    lines.append(f"* Comparator: {meta['comparator_kind']} ({meta['comparator_description']})"
                 + ("; performance against healthy controls overstates real-world performance"
                    if meta["comparator_type"] == "healthy" else ""))
    lines += ["", "## Key variables and missingness (measured)", ""]
    part = T["participants"]
    for c in ("age_range", "sex"):
        if c in part.columns:
            lines.append(f"* participants.{c}: {part[c].isna().mean():.1%} missing")
    for b, info in blocks.items():
        if info["kind"] == "device":
            w = dev["wide"][info["table"]]
            cols = [c for c in info["features"] if c in w.columns]
            wb = w[w[cols].notna().any(axis=1)] if cols else w
            lines.append(f"* {b} ({info['measurement_class']}, adapter {info['adapter_id']}): {len(wb)} participants; "
                         + _missing_line(wb, cols))
        elif info["kind"] == "ehr_labs":
            lab = T["participant_labs"]
            lines.append(f"* ehr_labs: {lab['participant_id'].nunique()} participants, {lab['loinc_code'].nunique()} "
                         f"LOINC codes, {len(lab)} rows; rows dropped for a minority unit: "
                         f"{out.get('labs_unit_dropped', 0)}")
        elif info["kind"] == "ehr_diagnoses":
            dx = T["participant_conditions"]
            dx = dx[dx["condition_type"] == "ehr_diagnosis"]
            lines.append(f"* ehr_diagnoses: {dx['participant_id'].nunique()} participants with >= 1 code, "
                         f"{dx['icd10cm_code'].nunique()} distinct ICD-10-CM codes")
        elif info["kind"] == "ehr_medications":
            med = T["participant_medications"]
            lines.append(f"* ehr_medications: {med['participant_id'].nunique()} participants, "
                         f"{med['rxnorm_cui'].nunique()} distinct RxNorm CUIs, {len(med)} rows")
        elif info["kind"] == "survey":
            sv = T["participant_survey"]
            g = sv[sv["instrument"] == info["instrument"]]
            lines.append(f"* {b} ({info['instrument']}): {g['participant_id'].nunique()} participants x "
                         f"{len(info['features'])} items; {len(g)} non-missing values")
        elif info["kind"] == "self_report_history":
            h = T["participant_conditions"]
            h = h[h["condition_type"] == "self_report_history"]
            lines.append(f"* self_report_history: {h['participant_id'].nunique()} participants, {len(h)} reports, "
                         f"{int((h['mapping_status'] == 'mapped').sum())} mapped to ICD-10-CM "
                         f"({h.loc[h['mapping_status'] == 'mapped', 'icd10cm_code'].nunique()} codes), "
                         f"{int((h['mapping_status'] != 'mapped').sum())} unmapped (kept, not modelled)")
        elif info["kind"] == "omics":
            om = T["participant_omics_linked"]
            g = om[om["modality"] == info["layer"]]
            n_feat = len(info["features"])
            n_p = g["participant_id"].nunique()
            lines.append(f"* {b}: {n_p} participants x {n_feat} features; {1 - len(g) / max(1, n_p * n_feat):.1%} "
                         "missing cells")
    if harmonized:
        lines += ["", "## Harmonized vocabulary (person_concepts)", "",
                  f"`person_concepts__{ds}`: {harmonized['rows']} rows, {harmonized['participants']} participants."]
        for k, v in sorted(harmonized.get("by_vocabulary", {}).items()):
            lines.append(f"* {k}: {v} rows")
        for line in harmonized.get("notes", []):
            lines.append(f"* {line}")
    lines += ["", "## Linkage strategy", "",
              f"Every file is joined on the owner's participant_id only, namespaced `{ds}:<id>`"
              + (" after salted hashing (salt held by the owner)" if meta["id_hashing"] == "salted_sha256" else "")
              + ". No row is linked to any other dataset, place, facility or registry: there is no cross-dataset "
                "participant key, and person rows carry no geographic column."]
    for k, v in rep.info.items():
        if "overlap" in k or "per block" in k:
            lines.append(f"* {k}: {v}")
    lines += ["", "## Limitations and caveats", "",
              "* User-supplied data: de-identification, label definitions and the permission to use the data here are "
              "attested by the data owner; the engine's validation is a safety net, not a certification.",
              "* Not public and not reproducible by others; results carry the tier `user_supplied_own_computation`.",
              "* Performance against the stated comparator only"
              + (" (healthy controls: optimistic)." if meta["comparator_type"] == "healthy" else "."),
              *[f"* Validation warning: {w}" for w in rep.warnings]]
    lines += ["", "## Processed outputs", ""]
    for name in written:
        n = len(read_table(name, columns=["participant_id"])) if table_exists(name) else 0
        lines.append(f"* `data/processed/{name}.parquet`: {n} rows")
    lines += ["", "## Reproduce", "", "```bash", "uv run measure-it byod validate <dataset folder>",
              "uv run measure-it byod ingest <dataset folder>", "```", ""]
    text = "\n".join(lines)
    from ..validate import forbidden_hits
    if forbidden_hits(text):   # owner text could contain a forbidden product term; keep the audit clean
        text = re.sub(r"(?i)\b(proves|cures|patient has|best clinic|definitive biomarker|optimal treatment|"
                      r"confirms mechanism|diagnosed by AI)\b", "[term withheld]", text)
    return text


# ------------------------------------------------------------------------------------------------------------------
# rebuilding the canonical unions and the digital person
# ------------------------------------------------------------------------------------------------------------------

def union_family(table: str) -> str:
    """Rebuild one non-owned canonical table from its partitions (pipeline.canonical_unions rules)."""
    from ..pipeline import _common_dtype
    parts = partitions(table)
    if not parts:
        return "no partitions"
    frames = [pd.read_parquet(p) for p in parts]
    for c in sorted({c for f in frames for c in f.columns}):
        dtypes = [f[c].dtype for f in frames if c in f.columns]
        if len({str(d) for d in dtypes}) > 1:
            target = _common_dtype(dtypes)
            for f in frames:
                if c in f.columns:
                    f[c] = f[c].astype(target)
    df = pd.concat(frames, ignore_index=True, sort=False)
    write_table(df, table, producer="measure_it.pipeline.canonical_unions",
                description=f"union of {[p.stem for p in parts]}")
    return f"{len(df)} rows from {len(parts)} partitions"


def rebuild_person_layer(echo=print, families: list[str] | None = None) -> dict:
    """Canonical unions touched by user datasets + digital person unions and DPV (no tracked result files)."""
    from ..pipeline import OWNED_CANONICALS
    from ..wearables import digital_person as DP
    done = {}
    fams = families or ["phenotype_signatures", "participant_device_features"]
    for t in fams:
        if t not in OWNED_CANONICALS:
            done[t] = union_family(t)
    if partitions("participant_adapter_scores"):
        done["participant_adapter_scores"] = union_family("participant_adapter_scores")
    if partitions("participant_survey"):
        done["participant_survey"] = union_family("participant_survey")
    # person_concepts: refreshed only when the pipeline has built the canonical union (it holds NHANES too)
    if partitions("person_concepts") and table_exists("person_concepts"):
        done["person_concepts"] = union_family("person_concepts")
    for t in DP.UNION_TABLES:
        if not partitions(t):
            continue
        df, info = DP.build_union(t)
        write_table(df, t, producer=DP.PRODUCER,
                    description=(f"canonical union of {info['partitions']}; dataset_id + namespaced participant_id "
                                 f"kept; no cross-dataset participant key. Columns cast to string because of "
                                 f"incompatible dtypes across partitions: {info['dtype_casts_to_string'] or 'none'}"))
        done[t] = f"{info['rows_written']} rows ({info['datasets']})"
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = DP.build_all(log=lambda s: None)
        DP.write_all(out)
    ds = out["digital_phenotype_datasets"]
    user = ds[ds["dataset_id"].astype(str).str.startswith(K.PREFIX)]
    done["digital_person"] = f"{len(ds)} dataset representations ({len(user)} user-supplied)"
    echo(f"  rebuilt canonical unions ({', '.join(done)}) and the Digital Phenotype Vectors")
    return done
