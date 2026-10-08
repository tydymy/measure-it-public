"""The `person_concepts` schema (docs/HARMONIZATION_CONTRACT.md section 1), shared helpers and the writer.

One row per participant x concept (x day when known). Partitions `person_concepts__<dataset_id>`; the canonical
`person_concepts` is the union of the partitions (pipeline `canonical_unions`, or `store.union_partitions`).
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from ..config import load_config, utc_now_iso
from ..provenance import add_provenance
from ..store import read_table, table_exists, write_table
from ..validate import GEO_COLUMN_RE

TABLE = "person_concepts"
PRODUCER = "measure_it.harmonize"
VERSION = "harmonize 2026-10-07.1"

COLUMNS = ["participant_id", "dataset_id", "domain", "vocabulary", "concept_code", "concept_label", "value_as_number",
           "unit", "value_as_string", "day_index", "mapping_method", "mapping_confidence", "source_variable",
           "platform"]
DOMAINS = {"demographic", "condition", "measurement", "drug", "device", "omics", "survey"}
VOCABULARIES = {"ICD10CM", "LOINC", "RXNORM", "DEMOG", "DEVICE", "OMICS", "SURVEY"}
DOMAIN_VOCAB = {"demographic": "DEMOG", "condition": "ICD10CM", "measurement": "LOINC", "drug": "RXNORM",
                "device": "DEVICE", "omics": "OMICS", "survey": "SURVEY"}
METHODS = {"native", "curated_map", "self_report_map", "rxnav_api", "rule"}
CONFIDENCE = {"high", "medium", "low"}
SEX_VALUES = {"female", "male", "other", "unknown"}
ICD10CM_RE = re.compile(r"^[A-Z]\d[0-9A-Z](?:\.[0-9A-Z]{1,4})?$")

# Which vocabularies a clinic EHR can evaluate (contract section 0: the exported phenotype says so per feature).
EHR_EVALUABLE = {"ICD10CM": True, "LOINC": True, "RXNORM": True, "DEMOG": True,
                 "SURVEY": False, "DEVICE": False, "OMICS": False}


def concepts_frame(**cols) -> pd.DataFrame:
    """Build a person_concepts frame from column arrays/scalars (missing columns -> null) with the contract dtypes."""
    n = max((len(v) for v in cols.values() if hasattr(v, "__len__") and not isinstance(v, str)), default=0)
    out = pd.DataFrame(index=range(n))
    for c in COLUMNS:
        v = cols.get(c)
        if isinstance(v, (pd.Series, pd.Index)):
            v = v.to_numpy()
        out[c] = v if v is not None else None
    return coerce(out)


def coerce(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = None
    for c in ("participant_id", "dataset_id", "domain", "vocabulary", "concept_code", "concept_label", "unit",
              "value_as_string", "mapping_method", "mapping_confidence", "source_variable", "platform"):
        df[c] = df[c].astype("string")
    df["value_as_number"] = pd.to_numeric(df["value_as_number"], errors="coerce").astype(float)
    df["day_index"] = pd.to_numeric(df["day_index"], errors="coerce").astype("Int64")
    return df[COLUMNS + [c for c in df.columns if c not in COLUMNS]]


def check_concepts(df: pd.DataFrame) -> list[str]:
    """Contract violations (empty list = ok)."""
    errs = []
    for c in COLUMNS:
        if c not in df.columns:
            errs.append(f"missing column {c}")
    if errs:
        return errs
    bad = set(df["domain"].dropna().unique()) - DOMAINS
    if bad:
        errs.append(f"unknown domain(s) {sorted(bad)}")
    bad = set(df["vocabulary"].dropna().unique()) - VOCABULARIES
    if bad:
        errs.append(f"unknown vocabulary(ies) {sorted(bad)}")
    pair = df[["domain", "vocabulary"]].drop_duplicates()
    wrong = pair[pair["domain"].map(DOMAIN_VOCAB) != pair["vocabulary"]]
    if len(wrong):
        errs.append(f"domain/vocabulary pairs not in the contract: {wrong.values.tolist()}")
    bad = set(df["mapping_method"].dropna().unique()) - METHODS
    if bad:
        errs.append(f"unknown mapping_method(s) {sorted(bad)}")
    bad = set(df["mapping_confidence"].dropna().unique()) - CONFIDENCE
    if bad:
        errs.append(f"unknown mapping_confidence(s) {sorted(bad)}")
    sr = df[(df["mapping_method"] == "self_report_map") & (df["mapping_confidence"] == "high")]
    if len(sr):
        errs.append(f"{len(sr)} self_report_map row(s) with confidence high (capped at medium)")
    icd = df.loc[df["vocabulary"] == "ICD10CM", "concept_code"].dropna().astype(str)
    badc = icd[~icd.str.fullmatch(ICD10CM_RE)]
    if len(badc):
        errs.append(f"{len(badc)} malformed ICD-10-CM code(s), e.g. {badc.iloc[0]!r}")
    if df["participant_id"].isna().any():
        errs.append("null participant_id")
    if not (df["participant_id"].astype(str).str.split(":", n=1).str[0] == df["dataset_id"].astype(str)).all():
        errs.append("participant_id not namespaced with its dataset_id")
    geo = [c for c in df.columns if GEO_COLUMN_RE.search(c) and c != "source_geographic_resolution"]
    if geo:
        errs.append(f"geographic column(s) {geo} (guardrail 2)")
    # one unit per LOINC code per dataset
    lab = df[df["vocabulary"] == "LOINC"]
    if len(lab):
        nu = lab.groupby(["dataset_id", "concept_code"])["unit"].nunique()
        if (nu > 1).any():
            errs.append(f"LOINC code(s) with more than one unit in a dataset: {nu[nu > 1].index.tolist()[:5]}")
    return errs


def write_concepts(df: pd.DataFrame, dataset_id: str, *, source_name: str, source_version: str, notes: str,
                   evidence_col: str = "_evidence_type", description: str = "") -> str:
    """Add provenance (per-row evidence type from `evidence_col`), check the contract, write the partition."""
    now = utc_now_iso()
    parts = []
    for et, g in df.groupby(evidence_col, sort=True):
        g = g.drop(columns=[evidence_col])
        parts.append(add_provenance(g, data_layer="person", source_name=source_name, source_version=source_version,
                                    retrieved_at=now, evidence_type=str(et), source_record_id=None,
                                    source_geographic_resolution="none",
                                    evidence_level=g["mapping_confidence"].astype(str).to_numpy(),
                                    provenance_notes=notes))
    out = pd.concat(parts, ignore_index=True) if parts else coerce(pd.DataFrame(columns=COLUMNS))
    out = out.sort_values(["participant_id", "domain", "concept_code"], kind="stable").reset_index(drop=True)
    out["source_record_id"] = (out["participant_id"].astype(str) + "|" + out["vocabulary"].astype(str) + "|"
                               + out["concept_code"].astype(str) + "|" + out.index.astype(str))
    errs = check_concepts(out)
    if errs:
        raise ValueError(f"person_concepts__{dataset_id}: " + "; ".join(errs))
    name = f"{TABLE}__{dataset_id}"
    write_table(out, name, producer=PRODUCER,
                description=description or f"harmonized person concepts of {dataset_id} (contract section 1)")
    return name


def read_concepts(dataset_ids: list[str] | str, columns: list[str] | None = None) -> pd.DataFrame:
    if isinstance(dataset_ids, str):
        dataset_ids = [dataset_ids]
    frames = []
    for ds in dataset_ids:
        name = f"{TABLE}__{ds}"
        if not table_exists(name):
            raise FileNotFoundError(f"{name} not found: run `python -m measure_it.harmonize.run` (NHANES) or "
                                    "`measure-it byod ingest` (user datasets)")
        frames.append(read_table(name, columns=columns))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns or COLUMNS)


# ------------------------------------------------------------------------------------------------------------------
# demographics
# ------------------------------------------------------------------------------------------------------------------

def age_band_from_years(age: pd.Series, topcode: float | None = None) -> pd.Series:
    """Exact (or top-coded) age -> 10-year band ('40-49'); >= 90 -> '90+'. A source top-coded below 90 (NHANES 80 or
    85) gets '<decade of the top code>+' for every age in that decade or above (e.g. '80+'), because the true band of
    a top-coded person is unknown."""
    a = pd.to_numeric(age, errors="coerce")
    lo = (np.floor(a / 10) * 10)
    band = lo.map(lambda v: f"{int(v)}-{int(v) + 9}" if pd.notna(v) else None)
    band = band.where(a < 90, "90+")
    if topcode is not None and topcode < 90:
        top = int(np.floor(topcode / 10) * 10)
        band = band.where(a < top, f"{top}+")
    band = band.astype(object)
    band[a.isna().to_numpy()] = None
    return band


def age_band_from_supplied(band: pd.Series) -> pd.Series:
    """A supplied band (>= 5 years, '90+' top-coded) -> the 10-year band containing it ('40-44' -> '40-49'); a band
    that straddles two decades is kept as supplied."""
    def one(v):
        if v is None or (isinstance(v, float) and np.isnan(v)) or pd.isna(v):
            return None
        s = str(v).strip()
        m = re.fullmatch(r"(\d{1,2})\+", s)
        if m:
            lo = int(m.group(1))
            return "90+" if lo >= 90 else f"{lo // 10 * 10}+"
        m = re.fullmatch(r"(\d{1,2})-(\d{1,2})", s)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            return f"{lo // 10 * 10}-{lo // 10 * 10 + 9}" if lo // 10 == hi // 10 else s
        return None
    return pd.Series([one(v) for v in band], index=band.index, dtype=object)


def age_mid(band: pd.Series) -> pd.Series:
    """Band midpoint in years ('40-49' -> 44.5; 'X+' -> X + 2.5, as measure_it.byod.person.covariate_matrix)."""
    def one(v):
        s = str(v)
        if s.endswith("+"):
            try:
                return float(s[:-1]) + 2.5
            except ValueError:
                return np.nan
        if "-" in s:
            a, b = s.split("-", 1)
            try:
                return (float(a) + float(b)) / 2
            except ValueError:
                return np.nan
        return np.nan
    return band.map(one).astype(float)


def demographic_rows(pid: pd.Series, ds: str, band: pd.Series, sex: pd.Series, *, band_var: str, sex_var: str,
                     band_method: str = "rule") -> pd.DataFrame:
    sex = sex.astype("string").str.lower().where(lambda s: s.isin(list(SEX_VALUES)), "unknown").fillna("unknown")
    a = concepts_frame(participant_id=pid, dataset_id=ds, domain="demographic", vocabulary="DEMOG",
                       concept_code="age_band", concept_label="Age band (10 years)", value_as_string=band,
                       mapping_method=band_method, mapping_confidence="high", source_variable=band_var)
    a = a[a["value_as_string"].notna()]
    s = concepts_frame(participant_id=pid, dataset_id=ds, domain="demographic", vocabulary="DEMOG",
                       concept_code="sex", concept_label="Sex", value_as_string=sex.to_numpy(),
                       mapping_method="native", mapping_confidence="high", source_variable=sex_var)
    return pd.concat([a, s], ignore_index=True)


def config(name: str) -> dict:
    return load_config(name)
