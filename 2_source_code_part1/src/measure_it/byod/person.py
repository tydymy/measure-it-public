"""A user dataset's person-level feature blocks and its Digital Phenotype Vector specification.

block_matrices(dataset_id, label_conditions)  -> {block: participant x feature frame} used by `byod evaluate`
digital_person_specs()                         -> one DatasetSpec per ingested user dataset, appended by
                                                  measure_it.wearables.digital_person.build_all (per dataset, within-
                                                  dataset z-scores and PCA: never pooled with any other dataset)

Diagnosis codes that map (measure_it.ontology.normalize, exact ICD-10-CM match) to the condition being analysed are
label-defining: they are left out of every model and flagged (and left out of the PCA fit) in the DPV.
"""
from __future__ import annotations

import json
from functools import lru_cache

import numpy as np
import pandas as pd

from ..store import partitions, read_table, table_exists
from . import common as K

MIN_DX_POSITIVES = 5
OMICS_MAX_FEATURES = 2000
OMICS_MAX_PCS = 5


def user_dataset_ids() -> list[str]:
    """Engine ids of every ingested user dataset (from the participants partitions on disk)."""
    return sorted(p.stem.split("__", 1)[1] for p in partitions("participants")
                  if p.stem.split("__", 1)[1].startswith(K.PREFIX))


def dataset_row(dataset_id: str) -> dict | None:
    ds = K.engine_id(dataset_id)
    if not table_exists(K.DATASETS_TABLE):
        return None
    d = read_table(K.DATASETS_TABLE)
    r = d[d["dataset_id"] == ds]
    return r.iloc[0].to_dict() if len(r) else None


def _members(condition_id: str) -> set[str]:
    try:
        from ..scoring.opportunity import condition_sets
        s = condition_sets().get(condition_id)
        if s:
            return set(s["members"])
    except Exception:  # noqa: BLE001
        pass
    return {condition_id}


@lru_cache(maxsize=4096)
def icd_condition(code: str) -> str | None:
    """Canonical condition of an ICD-10-CM code when the normaliser matches it exactly (else None)."""
    from ..ontology.normalize import normalize_condition
    r = normalize_condition(code)
    ids = {m.get("canonical_condition_id") for m in r.get("matches", []) or [] if m.get("canonical_condition_id")}
    return next(iter(ids)) if r.get("status") == "matched" and len(ids) == 1 else None


def _read(name: str) -> pd.DataFrame | None:
    return read_table(name) if table_exists(name) else None


def device_matrices(ds: str) -> dict[str, pd.DataFrame]:
    meta = dataset_row(ds) or {}
    blocks = json.loads(meta.get("blocks_json") or "{}")
    out = {}
    for table in ("participant_wearable_features", "participant_device_features"):
        t = _read(K.partition_name(table, ds))
        if t is None:
            continue
        t = t.set_index("participant_id")
        for b, info in blocks.items():
            if info.get("kind") != "device" or info.get("table") != table:
                continue
            cols = [c for c in info.get("features", []) if c in t.columns]
            if cols:
                out[b] = t[cols].apply(pd.to_numeric, errors="coerce").astype(float)
    return out


def labs_matrix(ds: str) -> pd.DataFrame | None:
    t = _read(K.partition_name("participant_labs", ds))
    if t is None or t.empty:
        return None
    t = t.assign(value=pd.to_numeric(t["value"], errors="coerce"))
    w = t.pivot_table(index="participant_id", columns="loinc_code", values="value", aggfunc="median")
    w.columns = [f"loinc_{c}" for c in w.columns]
    return w.astype(float)


def diagnosis_matrix(ds: str, participants: pd.Index, label_conditions: set[str] | None = None,
                     keep_label_defining: bool = False, condition_type: str = "ehr_diagnosis",
                     prefix: str = "icd3") -> tuple[pd.DataFrame | None, list[str]]:
    """ICD-10-CM 3-character category indicators (>= MIN_DX_POSITIVES participants). Returns (matrix, label-defining
    categories). Categories containing a code that maps to one of `label_conditions` are dropped unless
    keep_label_defining (the DPV keeps them, flagged). condition_type 'self_report_history' gives the same matrix for
    self-reported history mapped to ICD-10-CM (columns `hx3=<category>`)."""
    t = _read(K.partition_name("participant_conditions", ds))
    if t is None or "condition_type" not in t.columns:
        return None, []
    dx = t[(t["condition_type"] == condition_type) & t["icd10cm_code"].astype(str).str.len().gt(0)]
    if dx.empty:
        return None, []
    dx = dx.assign(cat=dx["icd10cm_code"].astype(str).str.replace(".", "", regex=False).str[:3])
    label_conditions = label_conditions or set()
    defining = sorted({c for c, code in zip(dx["cat"], dx["icd10cm_code"].astype(str))
                       if label_conditions and (icd_condition(code) in label_conditions)})
    counts = dx.drop_duplicates(["participant_id", "cat"])["cat"].value_counts()
    cats = [c for c in sorted(counts.index) if counts[c] >= MIN_DX_POSITIVES
            and (keep_label_defining or c not in defining)]
    if not cats:
        return None, defining
    has = dx.drop_duplicates(["participant_id", "cat"])
    M = pd.DataFrame(0.0, index=participants, columns=[f"{prefix}={c}" for c in cats])
    for c in cats:
        ids = has.loc[has["cat"] == c, "participant_id"]
        M.loc[M.index.isin(ids), f"{prefix}={c}"] = 1.0
    return M, defining


def medications_matrix(ds: str, participants: pd.Index) -> pd.DataFrame | None:
    """RxNorm CUI presence indicators (>= MIN_DX_POSITIVES participants; absence of a row = not recorded = 0)."""
    t = _read(K.partition_name("participant_medications", ds))
    if t is None or t.empty or "rxnorm_cui" not in t.columns:
        return None
    has = t.drop_duplicates(["participant_id", "rxnorm_cui"])
    counts = has["rxnorm_cui"].value_counts()
    cuis = sorted(c for c in counts.index if counts[c] >= MIN_DX_POSITIVES)
    if not cuis:
        return None
    M = pd.DataFrame(0.0, index=participants, columns=[f"rxcui={c}" for c in cuis])
    for c in cuis:
        M.loc[M.index.isin(has.loc[has["rxnorm_cui"] == c, "participant_id"]), f"rxcui={c}"] = 1.0
    return M


def survey_matrices(ds: str) -> dict[str, pd.DataFrame]:
    """{survey_<file name>: participant x '<INSTRUMENT>:<item>' (median over days when repeated)}."""
    t = _read(K.partition_name("participant_survey", ds))
    if t is None or t.empty:
        return {}
    out = {}
    for fname, g in t.groupby("instrument_file"):
        block = "survey_" + str(fname).rsplit(".", 1)[0][len("survey_"):]
        w = g.pivot_table(index="participant_id", columns="concept_code", values="value", aggfunc="median")
        w.columns = [str(c) for c in w.columns]
        out[block] = w.astype(float)
    return out


def omics_matrices(ds: str) -> dict[str, pd.DataFrame]:
    t = _read(K.partition_name("participant_omics_linked", ds))
    if t is None or t.empty:
        return {}
    out = {}
    for layer, g in t.groupby("modality"):
        w = g.pivot_table(index="participant_id", columns="feature_id", values="value", aggfunc="first")
        w.columns = [f"{layer}__{c}" for c in w.columns]
        out[f"omics_{layer}"] = w.astype(float)
    return out


def block_matrices(dataset_id: str, label_condition: str | None = None) -> dict[str, pd.DataFrame]:
    """{block name as in the manifest: participant_id x features} for evaluation (label-defining codes removed)."""
    ds = K.engine_id(dataset_id)
    P = read_table(K.partition_name("participants", ds))
    idx = pd.Index(P["participant_id"].astype(str))
    out = dict(device_matrices(ds))
    labs = labs_matrix(ds)
    if labs is not None:
        out["ehr_labs"] = labs
    conds = _members(label_condition) if label_condition else set()
    dxm, _ = diagnosis_matrix(ds, idx, conds, keep_label_defining=False)
    if dxm is not None:
        out["ehr_diagnoses"] = dxm
    hxm, _ = diagnosis_matrix(ds, idx, conds, keep_label_defining=False, condition_type="self_report_history",
                              prefix="hx3")
    if hxm is not None:
        out["self_report_history"] = hxm
    rxm = medications_matrix(ds, idx)
    if rxm is not None:
        out["ehr_medications"] = rxm
    out.update(survey_matrices(ds))
    out.update(omics_matrices(ds))
    return out


def covariate_matrix(dataset_id: str, covariates: list[str]) -> pd.DataFrame | None:
    """age_band -> band midpoint (years; 90+ -> 92.5); sex -> female indicator (other / unknown -> missing)."""
    if not covariates:
        return None
    P = read_table(K.partition_name("participants", K.engine_id(dataset_id))).set_index("participant_id")
    cols = {}
    if "age_band" in covariates and "age_range" in P.columns:
        def mid(v):
            v = str(v)
            if v.endswith("+"):
                return float(v[:-1]) + 2.5
            if "-" in v:
                a, b = v.split("-", 1)
                try:
                    return (float(a) + float(b)) / 2
                except ValueError:
                    return np.nan
            return np.nan
        cols["age_band_midpoint"] = P["age_range"].map(mid).astype(float)
    if "sex" in covariates and "sex" in P.columns:
        cols["sex_female"] = P["sex"].map({"female": 1.0, "male": 0.0}).astype(float)
    return pd.DataFrame(cols, index=P.index) if cols else None


# ------------------------------------------------------------------------------------------------------------------
# Digital Phenotype Vector specification (read by measure_it.wearables.digital_person.build_all)
# ------------------------------------------------------------------------------------------------------------------

def _omics_pca(w: pd.DataFrame, layer: str, ds: str) -> tuple[pd.DataFrame, list[dict]]:
    from sklearn.decomposition import PCA

    from ..config import SEED
    X = w.copy()
    X = X.loc[:, X.notna().mean() >= 0.5]
    if X.shape[1] > OMICS_MAX_FEATURES:
        X = X[X.var().sort_values(ascending=False).index[:OMICS_MAX_FEATURES]]
    Z = (X - X.mean()) / X.std(ddof=0).replace(0, np.nan)
    Z = Z.loc[:, Z.notna().any()].fillna(0.0)
    k = min(OMICS_MAX_PCS, Z.shape[1], max(Z.shape[0] - 1, 0))
    if k < 1:
        return pd.DataFrame(index=w.index), []
    pca = PCA(n_components=k, random_state=SEED, svd_solver="full").fit(Z.to_numpy())
    coords = pd.DataFrame(pca.transform(Z.to_numpy()), index=Z.index, columns=[f"pc{i + 1}" for i in range(k)])
    cum = np.cumsum(pca.explained_variance_ratio_)
    rows = [{"kind": "pca_row", "dataset_id": ds, "block": f"omics_{layer}", "component": i + 1, "n_fit": int(len(Z)),
             "n_features": int(Z.shape[1]), "explained_variance_ratio": float(pca.explained_variance_ratio_[i]),
             "cumulative_explained_variance": float(cum[i]),
             "top_loadings_json": json.dumps({c: round(float(v), 4) for c, v in
                                              pd.Series(pca.components_[i], index=Z.columns)
                                              .pipe(lambda s: s.reindex(s.abs().sort_values(ascending=False).index[:8]))
                                              .items()}),
             "weakly_structured_block": bool(pca.explained_variance_ratio_[0] < 0.10),
             "n_label_defining_excluded_from_fit": 0,
             "note": "owner-declared participant-linked omics; within-dataset z-scores; descriptive coordinates"}
            for i in range(k)]
    return coords, rows


def dataset_spec(ds: str):
    from ..wearables.digital_person import DatasetSpec, Feat
    meta = dataset_row(ds) or {}
    P = read_table(K.partition_name("participants", ds)).set_index("participant_id")
    labels = json.loads(meta.get("labels_json") or "{}")
    label_conditions = set().union(*[_members(v.get("condition_id")) for v in labels.values()
                                     if v.get("condition_id")] or [set()])
    frame = pd.DataFrame(index=P.index)
    feats: list = []
    part = K.partition_name("participants", ds)
    if "age_range" in P.columns and P["age_range"].notna().any():
        frame["age_range"] = P["age_range"]
        feats.append(Feat("age_range", "demographics", "categorical", part, note="age band as supplied"))
    if "sex" in P.columns:
        frame["sex"] = P["sex"].where(P["sex"].isin(["female", "male"]))
        feats.append(Feat("sex", "demographics", "binary", part))
    unavailable = {}
    blocks = json.loads(meta.get("blocks_json") or "{}")
    for b, M in device_matrices(ds).items():
        table = blocks.get(b, {}).get("table", "")
        dpv_block = "wearable" if table == "participant_wearable_features" else "clinical_exam"
        for c in M.columns:
            frame[c] = M[c].reindex(frame.index)
            feats.append(Feat(c, dpv_block, "numeric", K.partition_name(table, ds), note=f"{b} (user-supplied)"))
    labs = labs_matrix(ds)
    if labs is not None:
        for c in labs.columns:
            frame[c] = labs[c].reindex(frame.index)
            feats.append(Feat(c, "labs", "numeric", K.partition_name("participant_labs", ds), note="LOINC median"))
    else:
        unavailable["labs"] = "no ehr_labs.csv supplied"
    dxm, defining = diagnosis_matrix(ds, frame.index, label_conditions, keep_label_defining=True)
    if dxm is not None:
        for c in dxm.columns:
            frame[c] = dxm[c]
            feats.append(Feat(c, "symptoms_conditions", "binary", K.partition_name("participant_conditions", ds),
                              label_defining=c.split("=", 1)[1] in defining,
                              note="ICD-10-CM category from EHR codes" + (
                                  "; label-defining (maps to a label condition): flagged, left out of the PCA fit"
                                  if c.split("=", 1)[1] in defining else "")))
    else:
        unavailable["symptoms_conditions"] = "no ehr_diagnoses.csv supplied (or no category with >= 5 participants)"
    hxm, hx_def = diagnosis_matrix(ds, frame.index, label_conditions, keep_label_defining=True,
                                   condition_type="self_report_history", prefix="hx3")
    if hxm is not None:
        for c in hxm.columns:
            frame[c] = hxm[c]
            feats.append(Feat(c, "symptoms_conditions", "binary", K.partition_name("participant_conditions", ds),
                              label_defining=c.split("=", 1)[1] in hx_def,
                              note="self-reported history mapped to an ICD-10-CM category" + (
                                  "; label-defining: flagged, left out of the PCA fit"
                                  if c.split("=", 1)[1] in hx_def else "")))
    rxm = medications_matrix(ds, frame.index)
    if rxm is not None:
        for c in rxm.columns:
            frame[c] = rxm[c]
            feats.append(Feat(c, "medications", "binary", K.partition_name("participant_medications", ds),
                              note="RxNorm CUI (EHR medication list)"))
    else:
        unavailable.setdefault("medications", "no ehr_medications.csv supplied (or no CUI in >= 5 participants)")
    sv = survey_matrices(ds)
    for b, w in sv.items():
        for c in w.columns:
            frame[c] = w[c].reindex(frame.index)
            feats.append(Feat(c, "functional_status", "numeric", K.partition_name("participant_survey", ds),
                              note=f"{b} (patient-reported, user-supplied)"))
    if not sv:
        unavailable.setdefault("functional_status", "no survey_<instrument> file supplied")
    extra_blocks, extra_meta = {}, []
    om = omics_matrices(ds)
    for b, w in om.items():
        coords, rows = _omics_pca(w, b[len("omics_"):], ds)
        if len(coords.columns):
            extra_blocks[b] = coords
            extra_meta += rows
            extra_meta.append({"kind": "feature_row", "dataset_id": ds, "block": b,
                               "feature": f"(omics layer: {w.shape[1]} features -> {len(coords.columns)} PCs)",
                               "source_table": K.partition_name("participant_omics_linked", ds), "label_defining": False,
                               "note": "participant-linked omics (owner-declared); precomputed within-dataset PCA",
                               "n_population": int(len(frame)), "included": True, "exclusion_reason": "",
                               "n_observed": int(w.notna().any(axis=1).sum()),
                               "missing_fraction": round(1 - float(w.notna().any(axis=1).sum()) / max(1, len(frame)), 6),
                               "transform": "z-score within dataset -> PCA",
                               "dimensions": ";".join(f"pca:{b}:{c}" for c in coords.columns),
                               "has_missing_indicator": False})
    linked = bool(om)
    title = meta.get("title", ds)
    synthetic = bool(meta.get("synthetic"))
    statement = ("Omics here are participant-linked because the data owner supplied the same participant_id in the "
                 "omics and the other files (overlap counted at ingestion; the owner attests the linkage). "
                 if linked else
                 "No omics were supplied: molecular context is condition-level molecular enrichment reached through "
                 "the ontology, not individual multi-omics. ")
    return DatasetSpec(
        dataset_id=ds, frame=frame, feats=feats,
        population=f"all participants of user-supplied dataset {ds} (n = {len(frame)})",
        source_tables=[part] + sorted({f.source_table for f in feats if f.source_table != part}),
        unavailable_blocks=unavailable,
        labels_excluded=[f"{k} (label: outcome, never a dimension)" for k in labels] + ["participant ids"],
        representation_type="multimodal_person_representation" if linked else "digital_phenotype_vector",
        molecular_context_type=("participant_linked_omics_owner_declared" if linked
                                else "condition_level_molecular_enrichment_via_ontology"),
        molecular_context_statement=statement,
        representation_label=(("SYNTHETIC tutorial data. " if synthetic else "")
                              + f"User-supplied dataset '{title}' ({K.PROVENANCE_NOTE}; local-only). "
                              "Digital Phenotype Vector built within this dataset only; never pooled. " + statement),
        extra_blocks=extra_blocks, extra_meta=extra_meta)


def digital_person_specs() -> list:
    """DatasetSpec of every ingested user dataset ([] when there is none)."""
    out = []
    for ds in user_dataset_ids():
        if dataset_row(ds) is None:
            continue
        out.append(dataset_spec(ds))
    return out
