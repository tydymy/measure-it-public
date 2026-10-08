"""Wide participant x feature frames from person_concepts (docs/HARMONIZATION_CONTRACT.md section 1).

    feature_matrix(["nhanes", "nhanes0306"], ["ICD10CM:I10", "LOINC:2093-3", "DEMOG:age_mid"])

Feature keys `<VOCAB>:<code>`:
    ICD10CM:I10      presence (max over rows) of any code in the 3-character category (condition_level="category")
                     or of the exact code (condition_level="code")
    RXNORM:6809      presence of the ingredient
    LOINC:2093-3     median value per participant (one unit per code and dataset, by contract)
    SURVEY:<instrument>:<item>, DEVICE:<device>:<feature>, OMICS:<layer>:<feature>   median value
    DEMOG:age_mid    midpoint of the 10-year age band;   DEMOG:female  1 female / 0 male / NaN other, unknown

Absent presence features (ICD10CM / RXNORM), `absent="auto"`: within a dataset that has the key, a participant without
a row is 0 ("not recorded") unless the dataset also has explicit 0 rows for that key (survey questions answered "no",
e.g. NHANES self-report), in which case a participant without a row is missing (not asked / not answered). A key that
a dataset never records is missing (NaN) for all its participants: absence there is not evidence. `absent="zero"`
fills 0 everywhere a dataset has the key; `absent="nan"` never fills.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import schema as S

PRESENCE_VOCABS = {"ICD10CM", "RXNORM"}
VALUE_VOCABS = {"LOINC", "SURVEY", "DEVICE", "OMICS"}


def icd10_category(code) -> str | None:
    """ICD-10-CM code -> 3-character category ('I73.00' -> 'I73', 'G9332' -> 'G93'); None for empty input."""
    if code is None or (isinstance(code, float) and np.isnan(code)):
        return None
    s = str(code).strip().upper().replace(".", "")
    return s[:3] if len(s) >= 3 else None


def feature_key(vocab: str, code: str, condition_level: str = "category") -> str:
    if vocab == "ICD10CM" and condition_level == "category":
        return f"ICD10CM:{icd10_category(code)}"
    return f"{vocab}:{code}"


def keyed(concepts: pd.DataFrame, condition_level: str = "category") -> pd.DataFrame:
    """person_concepts rows -> (participant_id, dataset_id, key, value, vocabulary) with DEMOG turned numeric."""
    c = concepts
    out = []
    demo = c[c["vocabulary"] == "DEMOG"]
    if len(demo):
        a = demo[demo["concept_code"] == "age_band"]
        out.append(pd.DataFrame({"participant_id": a["participant_id"].to_numpy(), "dataset_id": a["dataset_id"].to_numpy(),
                                 "key": "DEMOG:age_mid", "value": S.age_mid(a["value_as_string"]).to_numpy(),
                                 "vocabulary": "DEMOG"}))
        s = demo[demo["concept_code"] == "sex"]
        out.append(pd.DataFrame({"participant_id": s["participant_id"].to_numpy(), "dataset_id": s["dataset_id"].to_numpy(),
                                 "key": "DEMOG:female",
                                 "value": s["value_as_string"].map({"female": 1.0, "male": 0.0}).astype(float).to_numpy(),
                                 "vocabulary": "DEMOG"}))
    rest = c[c["vocabulary"] != "DEMOG"]
    if len(rest):
        codes = rest["concept_code"].astype(str)
        if condition_level == "category":
            icd = rest["vocabulary"] == "ICD10CM"
            codes = codes.where(~icd, codes.map(icd10_category))
        out.append(pd.DataFrame({"participant_id": rest["participant_id"].to_numpy(),
                                 "dataset_id": rest["dataset_id"].to_numpy(),
                                 "key": (rest["vocabulary"].astype(str) + ":" + codes).to_numpy(),
                                 "value": rest["value_as_number"].astype(float).to_numpy(),
                                 "vocabulary": rest["vocabulary"].astype(str).to_numpy()}))
    if not out:
        return pd.DataFrame(columns=["participant_id", "dataset_id", "key", "value", "vocabulary"])
    return pd.concat(out, ignore_index=True)


def feature_matrix(dataset_ids, feature_keys: list[str] | None = None, *, condition_level: str = "category",
                   domains: list[str] | None = None, absent: str = "auto",
                   concepts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Wide participant x feature frame (index participant_id; float columns named by feature key).

    dataset_ids: one id or a list (each read from person_concepts__<id>), or None when `concepts` is given.
    domains: restrict to these person_concepts domains (e.g. ['condition', 'measurement', 'survey']).
    """
    if condition_level not in ("category", "code"):
        raise ValueError("condition_level must be 'category' or 'code'")
    if absent not in ("auto", "zero", "nan"):
        raise ValueError("absent must be 'auto', 'zero' or 'nan'")
    c = concepts if concepts is not None else S.read_concepts(
        dataset_ids, columns=["participant_id", "dataset_id", "domain", "vocabulary", "concept_code",
                              "value_as_number", "value_as_string"])
    if isinstance(dataset_ids, str):
        dataset_ids = [dataset_ids]
    if dataset_ids is not None and concepts is not None:
        c = c[c["dataset_id"].isin(dataset_ids)]
    universe = c[["participant_id", "dataset_id"]].drop_duplicates("participant_id")
    if domains is not None:
        c = c[c["domain"].isin(domains)]
    k = keyed(c, condition_level)
    if feature_keys is not None:
        k = k[k["key"].isin(set(feature_keys))]
    k = k[k["value"].notna()]
    pres = k[k["vocabulary"].isin(PRESENCE_VOCABS)]
    vals = k[~k["vocabulary"].isin(PRESENCE_VOCABS)]
    frames = []
    if len(vals):
        frames.append(vals.groupby(["participant_id", "key"])["value"].median().unstack("key"))
    if len(pres):
        w = pres.groupby(["participant_id", "key"])["value"].max().unstack("key")
        if absent != "nan":
            w = w.reindex(universe["participant_id"])
            ds_of = universe.set_index("participant_id")["dataset_id"].reindex(w.index)
            has = pres.groupby(["dataset_id", "key"])["value"].agg(lambda s: bool((s == 0).any()))
            for (ds, key), explicit in has.items():
                if absent == "auto" and explicit:
                    continue
                rows = (ds_of == ds).to_numpy()
                col = w[key]
                w.loc[rows, key] = col[rows].fillna(0.0)
        frames.append(w)
    if frames:
        fm = pd.concat(frames, axis=1)
    else:
        fm = pd.DataFrame(index=pd.Index([], name="participant_id"))
    fm = fm.reindex(universe["participant_id"])
    fm.index.name = "participant_id"
    fm.columns.name = None
    if feature_keys is not None:
        fm = fm.reindex(columns=list(feature_keys))
    else:
        fm = fm[sorted(fm.columns)]
    return fm.astype(float)


def feature_catalog(concepts: pd.DataFrame, condition_level: str = "category") -> pd.DataFrame:
    """One row per feature key: vocabulary, domain, label, codes (list), unit, platform, mapping methods and
    confidences, n participants, ehr_evaluable (contract section 0)."""
    c = concepts.copy()
    c["key"] = [feature_key(v, cc, condition_level) for v, cc in zip(c["vocabulary"].astype(str),
                                                                     c["concept_code"].astype(str))]
    demo = c["vocabulary"] == "DEMOG"
    c.loc[demo & (c["concept_code"] == "age_band"), "key"] = "DEMOG:age_mid"
    c.loc[demo & (c["concept_code"] == "sex"), "key"] = "DEMOG:female"
    rows = []
    for key, g in c.groupby("key", sort=True):
        vocab = str(g["vocabulary"].iloc[0])
        labels = g["concept_label"].dropna().astype(str).value_counts()
        label = ((labels.index[0] if len(labels) else key.split(":", 1)[1]) if vocab != "DEMOG"
                 else {"DEMOG:age_mid": "Age (band midpoint, years)", "DEMOG:female": "Female sex"}.get(key, key))
        rows.append({"feature_key": key, "vocabulary": vocab, "domain": str(g["domain"].iloc[0]), "label": label,
                     "codes": sorted(g["concept_code"].astype(str).unique().tolist()) if vocab != "DEMOG" else [],
                     "unit": (g["unit"].dropna().astype(str).iloc[0] if g["unit"].notna().any() else None),
                     "platform": (g["platform"].dropna().astype(str).iloc[0]
                                  if "platform" in g and g["platform"].notna().any() else None),
                     "mapping_methods": sorted(g["mapping_method"].dropna().astype(str).unique().tolist()),
                     "mapping_confidences": sorted(g["mapping_confidence"].dropna().astype(str).unique().tolist()),
                     "n_participants": int(g["participant_id"].nunique()),
                     "ehr_evaluable": bool(S.EHR_EVALUABLE.get(vocab, False))})
    return pd.DataFrame(rows)
