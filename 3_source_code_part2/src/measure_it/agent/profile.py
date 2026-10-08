"""Profile a `read.Dataset` into `profile.json` (docs/AGENT_CONTRACT.md section 1).

    prof = profile_dataset(read_input(path))      # dict, JSON-serialisable
    write_profile(prof, out_dir)                  # <out>/profile.json

Per column: inferred type, missingness, numeric min / median / max, categorical levels (only levels seen in >= 11 rows;
smaller ones collapse into ``<11 other``), unit (parsed from the label / header / REDCap field note, e.g.
``glucose (mg/dL)``), choices (REDCap / SPSS / Stata value labels), code pattern (ICD-10-CM, LOINC with its check
digit, RxNorm, SNOMED CT with its Verhoeff check digit) and privacy flags from `measure_it.byod.checks`
(identifier, geography, date, exact age, free text). Flagged columns never carry levels or numeric summaries. Per
table: participant-id candidates and the shape (wide / long / per_day).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..byod import checks as C
from .read import Dataset, Table, _code_str

SMALL_CELL = 11                  # a level is shown only when seen in at least this many rows (contract rule 2)
MAX_LEVELS = 20
OTHER_SMALL = f"<{SMALL_CELL} other"
OTHER_LARGE = "<other levels>"
SENSITIVE_FLAGS = {"identifier", "geography", "date", "exact_age", "free_text", "identifier?", "geography?"}

ID_NAME_RE = re.compile(r"^(record_id|participant_id|participant|participantid|subject_id|subjectid|subject|subjid|"
                        r"usubjid|study_id|studyid|person_id|patient_id|pid|ptid|seqn|id|record|redcap_record_id|"
                        r"eid|case_id|study_code|participant_code|subject_code|record_number_study)$")
DAY_NAME_RE = re.compile(r"^(day|day_index|study_day|visit_day|day_number|dayno|day_no|ady|lbdy)$")
CODE_NAME_RE = re.compile(r"(^|_)(code|loinc|icd|icd10|icd10cm|rxnorm|rxcui|snomed|sct|testcd|lbtestcd|analyte|test|"
                          r"parameter|concept_id|source_value)(_|$)")
VALUE_NAME_RE = re.compile(r"(^|_)(value|result|value_as_number|lborres|lbstresn|aval|measurement_value)(_|$)")
DATE_VALUE_RE = re.compile(r"^(\d{4}-\d{1,2}-\d{1,2}|\d{1,2}/\d{1,2}/\d{2,4})([ T]\d{1,2}:\d{2}(:\d{2})?(\.\d+)?(Z|[+-]\d{2}:?\d{2})?)?$")
TIME_PART_RE = re.compile(r"[ T]\d{1,2}:\d{2}")
BOOL_STRINGS = {"true", "false", "yes", "no", "y", "n", "t", "f", "checked", "unchecked"}
EXTRA_IDENTIFIER_NAMES = {"person_source_value", "provider_id", "provider_source_value", "care_site_source_value",
                          "redcap_survey_identifier"}
EXTRA_GEO_NAMES = {"location_id", "care_site_id", "redcap_data_access_group"}

LOINC_RE = re.compile(r"^\d{1,7}-\d$")
ICD10_RE = re.compile(r"^[A-Z]\d[0-9A-Z](?:\.[0-9A-Z]{1,4}|[0-9A-Z]{1,4})?$")
RXNORM_RE = re.compile(r"^\d{1,9}$")
SNOMED_RE = re.compile(r"^\d{6,18}$")

UNIT_CANON = {"mg/dl": "mg/dL", "mmol/l": "mmol/L", "umol/l": "umol/L", "µmol/l": "umol/L", "μmol/l": "umol/L",
              "nmol/l": "nmol/L", "pmol/l": "pmol/L", "g/l": "g/L", "g/dl": "g/dL", "mg/l": "mg/L", "ng/ml": "ng/mL",
              "pg/ml": "pg/mL", "ug/ml": "ug/mL", "ug/dl": "ug/dL", "ug/l": "ug/L", "iu/l": "IU/L", "u/l": "U/L",
              "miu/l": "mIU/L", "uiu/ml": "uIU/mL", "%": "%", "pct": "%", "percent": "%", "bpm": "/min",
              "beats/min": "/min", "/min": "/min", "breaths/min": "/min", "ms": "ms", "msec": "ms", "kg": "kg",
              "cm": "cm", "mm": "mm", "m": "m", "mmhg": "mm[Hg]", "mm[hg]": "mm[Hg]", "kg/m2": "kg/m2",
              "kg/m^2": "kg/m2", "npx": "NPX", "mmol/mol": "mmol/mol", "10^9/l": "10*9/L", "10*9/l": "10*9/L",
              "10^3/ul": "10*3/uL", "cells/ul": "/uL", "°c": "Cel", "degc": "Cel", "celsius": "Cel", "h": "h",
              "hours": "h", "hrs": "h", "min": "min", "minutes": "min", "steps": "steps", "years": "a", "yrs": "a",
              "ml/min/1.73m2": "mL/min/{1.73_m2}", "ratio": "ratio", "a.u.": "a.u.", "au": "a.u.",
              "capillaries/mm": "/mm", "/mm": "/mm", "um": "um", "µm": "um", "μm": "um"}
UNIT_SUFFIX_RE = re.compile(r"_(mg_?dl|mmol_?l|umol_?l|nmol_?l|g_?l|g_?dl|mg_?l|ng_?ml|pg_?ml|pct|percent|bpm|ms|kg|"
                            r"cm|mmhg|npx|kg_?m2)$", re.IGNORECASE)


# ------------------------------------------------------------------------------------------------------------------
# check digits
# ------------------------------------------------------------------------------------------------------------------
def luhn_ok(number: str) -> bool:
    """LOINC check digit (mod 10, Luhn): '2345-7' -> True."""
    if not LOINC_RE.match(number):
        return False
    body, check = number.split("-")
    total = 0
    for i, ch in enumerate(reversed(body)):
        d = int(ch)
        if i % 2 == 0:
            d *= 2
            if d > 9:
                d -= 9
        total += d
    return (10 - total % 10) % 10 == int(check)


_VD = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 2, 3, 4, 0, 6, 7, 8, 9, 5], [2, 3, 4, 0, 1, 7, 8, 9, 5, 6],
       [3, 4, 0, 1, 2, 8, 9, 5, 6, 7], [4, 0, 1, 2, 3, 9, 5, 6, 7, 8], [5, 9, 8, 7, 6, 0, 4, 3, 2, 1],
       [6, 5, 9, 8, 7, 1, 0, 4, 3, 2], [7, 6, 5, 9, 8, 2, 1, 0, 4, 3], [8, 7, 6, 5, 9, 3, 2, 1, 0, 4],
       [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]]
_VP = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 5, 7, 6, 2, 8, 3, 0, 9, 4], [5, 8, 0, 3, 7, 9, 6, 1, 4, 2],
       [8, 9, 1, 6, 0, 4, 3, 5, 2, 7], [9, 4, 5, 3, 1, 2, 6, 8, 7, 0], [4, 2, 8, 6, 5, 7, 3, 9, 0, 1],
       [2, 7, 9, 3, 8, 0, 6, 4, 1, 5], [7, 0, 4, 6, 9, 1, 3, 2, 5, 8]]


def verhoeff_ok(number: str) -> bool:
    """SNOMED CT identifiers end in a Verhoeff check digit."""
    if not number.isdigit():
        return False
    c = 0
    for i, ch in enumerate(reversed(number)):
        c = _VD[c][_VP[i % 8][int(ch)]]
    return c == 0


def snomed_ok(code: str) -> bool:
    return bool(SNOMED_RE.match(code)) and verhoeff_ok(code) and code[-3:-1] in ("00", "10")   # concept partition


def icd10_ok(code: str) -> bool:
    return bool(ICD10_RE.match(str(code).strip().upper()))


# ------------------------------------------------------------------------------------------------------------------
# column helpers
# ------------------------------------------------------------------------------------------------------------------
def _num(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    v = float(v)
    return int(v) if v.is_integer() and abs(v) < 1e15 else float(f"{v:.6g}")


def parse_unit(*texts: str | None) -> str | None:
    """'Glucose (mg/dL)' -> 'mg/dL'; 'glucose_mmol_l' -> 'mmol/L'; a bare REDCap note 'mmol/L' -> 'mmol/L'."""
    for t in texts:
        if not t:
            continue
        t = str(t).strip()
        for m in re.finditer(r"[\(\[]([^\(\)\[\]]{1,20})[\)\]]", t):
            cand = m.group(1).strip().lower().replace(" ", "")
            if cand in UNIT_CANON:
                return UNIT_CANON[cand]
            if re.fullmatch(r"(?:[a-zµμ]{1,5}|10[\^*]\d)/[a-z0-9.\[\]]{1,8}", cand):
                return m.group(1).strip()
        low = t.lower().replace(" ", "")
        if low in UNIT_CANON:
            return UNIT_CANON[low]
        m = UNIT_SUFFIX_RE.search(t)
        if m:
            key = m.group(1).lower()
            key = {"mgdl": "mg/dl", "mg_dl": "mg/dl", "mmoll": "mmol/l", "mmol_l": "mmol/l", "umoll": "umol/l",
                   "umol_l": "umol/l", "nmoll": "nmol/l", "nmol_l": "nmol/l", "gl": "g/l", "g_l": "g/l",
                   "gdl": "g/dl", "g_dl": "g/dl", "mgl": "mg/l", "mg_l": "mg/l", "ngml": "ng/ml", "ng_ml": "ng/ml",
                   "pgml": "pg/ml", "pg_ml": "pg/ml", "kgm2": "kg/m2", "kg_m2": "kg/m2"}.get(key, key)
            if key in UNIT_CANON:
                return UNIT_CANON[key]
    return None


def _str_values(s: pd.Series) -> pd.Series:
    return s.dropna().map(_code_str)


def detect_code_pattern(t: Table, col: str, s: pd.Series, label: str | None) -> str:
    sys_ = t.code_systems.get(col)
    if sys_:
        return {"icd10cm": "icd10", "icd10": "icd10"}.get(sys_, sys_)
    vals = _str_values(s)
    if vals.empty:
        return "none"
    uniq = pd.Series(vals.unique()[:500])
    hint = f"{col} {label or ''}".lower()
    share = lambda f: float(uniq.map(f).mean())  # noqa: E731
    if share(luhn_ok) >= 0.9 and (len(uniq) >= 3 or "loinc" in hint):
        return "loinc"
    if share(lambda v: icd10_ok(v)) >= 0.9 and (re.search(r"icd|diag|dx", hint) or (
            len(uniq) >= 3 and uniq.map(lambda v: "." in v or len(v) >= 4).any())):
        return "icd10"
    if re.search(r"rxnorm|rxcui|rx_?code|drug_?code|med_?code|drug_source|medication_code", hint) and share(lambda v: bool(RXNORM_RE.match(v))) >= 0.9:
        return "rxnorm"
    if re.search(r"snomed|sct", hint) and share(snomed_ok) >= 0.9:
        return "snomed"
    return "none"


def _value_flags(s: pd.Series, col: str, file: str) -> list[str]:
    out = []
    for msg in C.value_problems(s, col, file):
        m = re.search(r"look like: (.+?) \(e\.g\.", msg)
        kind = m.group(1) if m else ""
        if kind in ("ZIP code", "coordinate pair"):
            out.append("geography")
        elif kind == "exact date":
            out.append("date")
        elif kind:
            out.append("identifier")
    return out


def name_flags(col: str, label: str | None) -> list[str]:
    flags = []
    low = re.sub(r"[^a-z0-9]+", "_", str(col).lower()).strip("_")
    prob = C.column_problem(col)
    if low in EXTRA_IDENTIFIER_NAMES:
        flags.append("identifier")
    elif low in EXTRA_GEO_NAMES:
        flags.append("geography")
    elif prob:
        flags.append("geography" if "geographic" in prob else "identifier" if "identifies" in prob else
                     "date" if "calendar date" in prob else "exact_age" if "exact age" in prob else "identifier")
    elif re.fullmatch(r"age(_at_\w+)?|age_?(yrs?|years?)|ridageyr|age_years_\w+", low):
        flags.append("exact_age")
    if label and not flags:
        lab = re.sub(r"\(choice=[^)]*\)", "", label)
        lab = re.sub(r"\b(assigned|recorded)?\s*at birth\b", " ", lab, flags=re.IGNORECASE)   # 'sex assigned at birth'
        lp = C.column_problem(lab)
        if lp and "geographic" in lp:
            flags.append("geography?")
        elif lp and "identifies" in lp:
            flags.append("identifier?")
        elif lp and "calendar date" in lp:
            flags.append("date?")
        elif re.match(r"^\s*age\b", label, re.IGNORECASE) and re.search(r"\byears?\b", label, re.IGNORECASE) \
                and not re.search(r"band|group|range|onset|diagnos", label, re.IGNORECASE):
            flags.append("exact_age")
    return flags


def infer_dtype(s: pd.Series, choices: dict | None, validation: str | None, field_type: str | None) -> str:
    nn = s.dropna()
    if pd.api.types.is_datetime64_any_dtype(s):
        if nn.empty:
            return "date"
        t = pd.to_datetime(nn)
        return "datetime" if ((t.dt.hour != 0) | (t.dt.minute != 0) | (t.dt.second != 0)).any() else "date"
    if validation and validation.startswith("date"):
        return "datetime" if "datetime" in validation else "date"
    if pd.api.types.is_bool_dtype(s):
        return "bool"
    if field_type in ("notes",):
        return "text"
    if pd.api.types.is_numeric_dtype(s):
        if nn.empty:
            return "float"
        u = set(pd.unique(nn))
        if u <= {0, 1}:
            return "bool"
        if choices:
            return "category"
        return "int" if np.all(np.mod(nn.astype(float), 1) == 0) else "float"
    vals = nn.astype(str)
    if vals.empty:
        return "text"
    low = set(vals.str.lower().unique())
    if low <= BOOL_STRINGS and len(low) <= 3:
        return "bool"
    dshare = vals.map(lambda v: bool(DATE_VALUE_RE.match(v))).mean()
    if dshare >= 0.9:
        return "datetime" if vals.map(lambda v: bool(TIME_PART_RE.search(v))).mean() > 0.5 else "date"
    if choices:
        return "category"
    n_u = vals.nunique()
    mean_len = vals.str.len().mean()
    if n_u <= max(25, 0.05 * len(vals)) and mean_len < 40:
        return "category"
    return "text"


def small_cell_levels(s: pd.Series) -> list[list]:
    """Top levels with counts; counts below SMALL_CELL are never shown (collapsed into '<11 other')."""
    vc = s.dropna().map(_code_str).value_counts()
    shown, small, large = [], 0, 0
    for v, c in vc.items():
        c = int(c)
        if c < SMALL_CELL:
            small += c
        elif len(shown) < MAX_LEVELS:
            shown.append([str(v), c])
        else:
            large += c
    if large:
        shown.append([OTHER_LARGE, large])
    if small:
        shown.append([OTHER_SMALL, small if small >= SMALL_CELL else None])
    return shown


def is_free_text(s: pd.Series, dtype: str, col: str) -> bool:
    if dtype != "text":
        return False
    vals = s.dropna().astype(str)
    if vals.empty:
        return False
    if re.search(r"(^|_)(notes?|comments?|remarks?|free_?text|other_?specify|specify|describe|description|narrative)(_|$)",
                 col.lower()):
        return True
    return vals.str.len().mean() > 30 or (vals.str.contains(r"\s").mean() > 0.5 and vals.nunique() > 0.5 * len(vals))


def profile_column(t: Table, col: str) -> dict:
    s = t.df[col]
    n = len(s)
    label = t.labels.get(col)
    choices = t.choices.get(col) or None
    ftype = t.field_types.get(col)
    dtype = infer_dtype(s, choices, t.validation.get(col), ftype)
    flags = name_flags(col, label)
    if col in t.identifier_fields or ftype == "file":
        flags.append("identifier")
    if col in t.structural:
        flags.append("structural")
    if dtype in ("category", "text", "bool") or (dtype in ("int",) and s.astype("string").str.match(r"^0\d{4}$").any()):
        flags += _value_flags(s.astype("string") if not pd.api.types.is_string_dtype(s) else s, col, t.source_file)
    if dtype in ("date", "datetime") and "date" not in flags:
        flags.append("date")
    if is_free_text(s, dtype, col):
        flags.append("free_text")
    flags = list(dict.fromkeys(flags))
    sensitive = bool(set(flags) & SENSITIVE_FLAGS)
    nn = s.dropna()
    out = {"name": col, "label": label, "dtype": dtype, "n_missing": int(s.isna().sum()),
           "share_missing": round(float(s.isna().mean()), 4) if n else None,
           "n_unique": int(nn.map(_code_str).nunique()) if len(nn) else 0,
           "numeric": None, "levels": None,
           "unit": parse_unit(label, col, t.notes.get(col)),
           "choices": choices, "code_pattern": "none", "flags": flags}
    if ftype:
        out["field_type"] = ftype
    if t.forms.get(col):
        out["form"] = t.forms[col]
    if not sensitive:
        if dtype in ("int", "float") and len(nn):
            x = pd.to_numeric(nn, errors="coerce").dropna()
            if len(x):
                out["numeric"] = {"min": _num(x.min()), "median": _num(x.median()), "max": _num(x.max())}
        if dtype in ("bool", "category"):
            out["levels"] = small_cell_levels(s)
        out["code_pattern"] = detect_code_pattern(t, col, s, label)
    return out


def id_candidates(t: Table, cols: list[dict]) -> list[str]:
    bad = {c["name"] for c in cols if {"identifier", "geography", "date", "free_text", "structural"} & set(c["flags"])}
    by = {c["name"]: c for c in cols}
    n = len(t.df)
    exact = [c for c in t.df.columns if c not in bad and ID_NAME_RE.match(re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_"))]
    order = {"record_id": 0, "participant_id": 0, "subject_id": 0, "usubjid": 0, "person_id": 0, "patient_id": 1,
             "seqn": 1, "study_id": 1}
    exact.sort(key=lambda c: order.get(c.lower(), 2))
    if t.kind == "redcap_records" or (t.kind or "").startswith("redcap_repeating"):
        first = t.df.columns[0]
        if first not in bad and first not in exact and not first.startswith("redcap_"):
            exact.insert(0, first)
    other = []
    for c in t.df.columns:
        if c in bad or c in exact or by[c]["dtype"] in ("float", "date", "datetime", "bool", "text"):
            continue
        low = c.lower()
        if (low.endswith("_id") or low.endswith("id")) and by[c]["n_unique"] >= 0.5 * max(1, n - by[c]["n_missing"]) \
                and not re.search(r"(occurrence|exposure|measurement|observation|visit|concept|sample|specimen|"
                                  r"event|row|record_line)_?id$", low):
            other.append(c)
    return exact + other


def detect_shape(t: Table, cols: list[dict], ids: list[str]) -> str:
    if (t.kind or "").startswith(("fhir:Condition", "fhir:Observation", "fhir:Medication")) or t.kind in (
            "omop:condition_occurrence", "omop:measurement", "omop:drug_exposure", "omop:observation",
            "omop:procedure_occurrence", "omop:device_exposure"):
        return "long"
    if not ids:
        return "wide"
    pid = t.df[ids[0]]
    if not pid.dropna().duplicated().any():
        return "wide"
    names = {c["name"]: c for c in cols}
    has_code = any(c["code_pattern"] != "none" or CODE_NAME_RE.search(c["name"].lower()) for c in cols)
    has_value = any(VALUE_NAME_RE.search(c["name"].lower()) for c in cols)
    if has_code and has_value:
        return "long"
    if any(DAY_NAME_RE.match(n.lower()) or names[n]["dtype"] in ("date", "datetime") for n in names):
        return "per_day"
    return "long"


def profile_table(t: Table) -> dict:
    cols = [profile_column(t, c) for c in t.df.columns]
    ids = id_candidates(t, cols)
    for c in cols:
        if c["name"] in ids[:1]:
            probs = C.id_problems(t.df[c["name"]].dropna().map(_code_str), t.source_file)
            if probs:
                c["flags"].append("id_values_look_identifying")
    shape = detect_shape(t, cols, ids)
    labels_from = ("redcap_data_dictionary" if t.field_types else "variable_labels" if t.labels else "none")
    return {"name": t.name, "source_file": t.source_file, "sheet": t.sheet, "kind": t.kind, "n_rows": int(len(t.df)),
            "n_columns": int(t.df.shape[1]), "shape": shape, "id_candidates": ids, "labels_from": labels_from,
            "columns": cols}


def profile_dataset(ds: Dataset) -> dict:
    return {"input": ds.input, "format": ds.format, "container": ds.container,
            "tables": [profile_table(t) for t in ds.tables],
            "dictionary": ({k: v for k, v in ds.dictionary.items()} if ds.dictionary else None),
            "warnings": list(ds.warnings), "small_cell_threshold": SMALL_CELL,
            "generated_by": "measure_it.agent.profile"}


def write_profile(profile: dict, out: str | Path) -> Path:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    p = out / "profile.json"
    p.write_text(json.dumps(profile, indent=2, default=str) + "\n")
    return p
