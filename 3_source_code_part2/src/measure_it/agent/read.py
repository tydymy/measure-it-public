"""Load (almost) any input path into named pandas tables (launch agent part A, docs/AGENT_CONTRACT.md).

    ds = read_input("export.zip")          # -> Dataset(format, container, tables=[Table, ...], dictionary, warnings)
    for t in ds.tables: t.name, t.df, t.labels, t.choices

Supported: CSV / TSV / TXT (delimiter sniffing; utf-8, utf-8-sig, cp1252, latin-1), Excel (.xlsx / .xlsm / .xls, every
sheet), Parquet, JSON (records, nested -> flattened, one table per list of records) and JSONL / NDJSON, SAS XPT and
sas7bdat, SPSS .sav, Stata .dta (pyreadstat; variable labels and value labels kept), ZIP archives of these and folders of
these. Special layouts:

* REDCap export (a data CSV + the data dictionary CSV, UI or API column names): field labels, choice maps (radio,
  dropdown, checkbox ``field___code`` columns, yesno, truefalse), field types, form names and the dictionary's
  ``Identifier?`` column are attached; repeating instruments become their own tables; REDCap structural columns
  (``redcap_event_name``, ``redcap_repeat_*``, ``<form>_complete``, data access group) are marked.
* OMOP CDM tables (person, condition_occurrence, measurement, ...): format ``omop``, table kind ``omop:<table>``.
* FHIR R4 Bundle JSON / NDJSON (Patient, Condition, Observation, MedicationStatement / MedicationRequest): long tables
  with one column per code system (ICD-10-CM, SNOMED CT, LOINC, RxNorm).

Nothing here sends data anywhere. A format whose reader package is missing raises `MissingDependency` naming it.
"""
from __future__ import annotations

import csv
import io
import json
import re
import tempfile
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

TEXT_SUFFIXES = {".csv", ".tsv", ".txt", ".tab"}
EXCEL_SUFFIXES = {".xlsx", ".xlsm", ".xls"}
STAT_SUFFIXES = {".xpt", ".sas7bdat", ".sav", ".zsav", ".dta"}
JSON_SUFFIXES = {".json", ".jsonl", ".ndjson"}
SUPPORTED_SUFFIXES = TEXT_SUFFIXES | EXCEL_SUFFIXES | STAT_SUFFIXES | JSON_SUFFIXES | {".parquet", ".pq", ".zip"}
FORMAT_OF_SUFFIX = {".csv": "csv", ".tsv": "tsv", ".tab": "tsv", ".txt": "csv", ".xlsx": "excel", ".xlsm": "excel",
                    ".xls": "excel", ".parquet": "parquet", ".pq": "parquet", ".json": "json", ".jsonl": "jsonl",
                    ".ndjson": "jsonl", ".xpt": "xpt", ".sas7bdat": "sas7bdat", ".sav": "sav", ".zsav": "sav",
                    ".dta": "dta"}
ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")
NA_VALUES = ["", "NA", "N/A", "n/a", "na", "NaN", "nan", "NULL", "null", "None", "#N/A", ".", "-"]
NUMERIC_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")
LEADING_ZERO_RE = re.compile(r"^[+-]?0\d")

OMOP_TABLES = {"person", "observation_period", "visit_occurrence", "visit_detail", "condition_occurrence",
               "drug_exposure", "procedure_occurrence", "device_exposure", "measurement", "observation", "death",
               "note", "note_nlp", "specimen", "fact_relationship", "location", "care_site", "provider",
               "payer_plan_period", "cost", "drug_era", "dose_era", "condition_era", "episode", "concept",
               "vocabulary", "concept_relationship", "concept_ancestor", "cdm_source"}
REDCAP_STRUCTURAL = {"redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance",
                     "redcap_data_access_group", "redcap_survey_identifier", "redcap_arm"}
FHIR_SYSTEMS = {
    "http://hl7.org/fhir/sid/icd-10-cm": "icd10cm", "http://hl7.org/fhir/sid/icd-10": "icd10",
    "urn:oid:2.16.840.1.113883.6.90": "icd10cm", "http://snomed.info/sct": "snomed", "http://loinc.org": "loinc",
    "http://www.nlm.nih.gov/research/umls/rxnorm": "rxnorm", "urn:oid:2.16.840.1.113883.6.88": "rxnorm",
}


class ReadError(ValueError):
    """The input cannot be read (unsupported, empty, corrupt)."""


class MissingDependency(ImportError):
    """A reader package is not installed. The message names it."""


@dataclass
class Table:
    name: str
    df: pd.DataFrame
    source_file: str
    sheet: str | None = None
    format: str = "csv"
    kind: str | None = None                       # redcap_records | redcap_repeating:<form> | omop:<t> | fhir:<r>
    labels: dict[str, str] = field(default_factory=dict)
    choices: dict[str, dict[str, str]] = field(default_factory=dict)
    field_types: dict[str, str] = field(default_factory=dict)      # REDCap field type per column
    forms: dict[str, str] = field(default_factory=dict)            # REDCap form (instrument) per column
    notes: dict[str, str] = field(default_factory=dict)            # REDCap field note (often the unit)
    validation: dict[str, str] = field(default_factory=dict)       # REDCap text validation (date_ymd, number, ...)
    identifier_fields: set[str] = field(default_factory=set)       # REDCap 'Identifier?' = y
    structural: set[str] = field(default_factory=set)              # REDCap event / repeat / form-status columns
    code_systems: dict[str, str] = field(default_factory=dict)     # column -> icd10cm|snomed|loinc|rxnorm (FHIR)


@dataclass
class Dataset:
    input: str
    format: str
    container: str | None = None                 # zip | folder | None
    tables: list[Table] = field(default_factory=list)
    dictionary: dict | None = None
    warnings: list[str] = field(default_factory=list)

    def table(self, name: str) -> Table:
        for t in self.tables:
            if t.name == name:
                return t
        raise KeyError(name)


# ------------------------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------------------------
def slug(text: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")
    return s or "table"


def _code_str(v) -> str:
    """Value-label / choice code as a string: 1.0 -> '1'."""
    if isinstance(v, (float, np.floating)) and float(v).is_integer():
        return str(int(v))
    return str(v).strip()


def coerce_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Strip text, turn blank strings into missing and convert all-numeric text columns to numbers. Columns with a
    leading-zero value ('02139', '007') stay text: they are codes, not quantities."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    df = df.loc[:, [c for c in df.columns if not (c.startswith("Unnamed:") and df[c].isna().all())]]
    for c in df.columns:
        s = df[c]
        if s.dtype != object and not pd.api.types.is_string_dtype(s):
            continue
        st = s.map(lambda v: v.strip() if isinstance(v, str) else v)
        st = st.map(lambda v: np.nan if (isinstance(v, str) and v in NA_VALUES) else v)
        nn = st.dropna()
        if len(nn) and nn.map(lambda v: isinstance(v, str)).all():
            if nn.map(lambda v: bool(NUMERIC_RE.match(v))).all() and not nn.map(lambda v: bool(LEADING_ZERO_RE.match(v)) and not v.lstrip("+-").startswith("0.")).any():
                st = pd.to_numeric(st, errors="coerce")
        df[c] = st
    return df.dropna(how="all")


def _decode(raw: bytes) -> tuple[str, str]:
    for enc in ENCODINGS:
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace"), "latin-1"


def _sniff_delimiter(text: str, suffix: str) -> str:
    sample = text[:65536]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;|").delimiter
    except csv.Error:
        if suffix in (".tsv", ".tab"):
            return "\t"
        first = sample.splitlines()[0] if sample else ""
        return max(",\t;|", key=first.count) if first else ","


def _text_table(raw: bytes, suffix: str) -> tuple[pd.DataFrame, dict]:
    text, enc = _decode(raw)
    delim = _sniff_delimiter(text, suffix)
    df = pd.read_csv(io.StringIO(text), sep=delim, dtype=str, keep_default_na=False, na_values=NA_VALUES,
                     engine="python")
    return df, {"encoding": enc, "delimiter": delim}


def _stat_reader(suffix: str):
    try:
        import pyreadstat
    except ImportError as e:  # pragma: no cover - pyreadstat is a project dependency
        raise MissingDependency(f"reading {suffix} files needs the package 'pyreadstat' (pip install pyreadstat)") from e
    return {".xpt": pyreadstat.read_xport, ".sas7bdat": pyreadstat.read_sas7bdat, ".sav": pyreadstat.read_sav,
            ".zsav": pyreadstat.read_sav, ".dta": pyreadstat.read_dta}[suffix]


def _read_stat(path: Path, suffix: str, label: str) -> Table:
    reader = _stat_reader(suffix)
    try:
        df, meta = reader(str(path))
    except Exception as e:  # noqa: BLE001 - pyreadstat raises its own error types
        raise ReadError(f"{label}: cannot read as {suffix}: {e}") from e
    labels = {k: v for k, v in (meta.column_names_to_labels or {}).items() if v}
    choices = {col: {_code_str(k): str(v) for k, v in vl.items()}
               for col, vl in (getattr(meta, "variable_value_labels", None) or {}).items() if vl}
    return Table(name=slug(path.stem), df=coerce_frame(df), source_file=label, format=FORMAT_OF_SUFFIX[suffix],
                 labels=labels, choices=choices)


def _read_excel(path: Path, label: str) -> list[Table]:
    engine = "xlrd" if path.suffix.lower() == ".xls" else "openpyxl"
    try:
        __import__(engine)
    except ImportError as e:  # pragma: no cover - both are project dependencies
        raise MissingDependency(f"reading {path.suffix} files needs the package '{engine}'") from e
    try:
        sheets = pd.read_excel(path, sheet_name=None, dtype=object, engine=engine)
    except Exception as e:  # noqa: BLE001
        raise ReadError(f"{label}: cannot read as Excel: {e}") from e
    out = []
    for sheet, df in sheets.items():
        df = df.dropna(how="all").dropna(axis=1, how="all")
        if df.empty:
            continue
        out.append(Table(name=slug(sheet), df=coerce_frame(df), source_file=label, sheet=str(sheet), format="excel"))
    return out


def _flatten_records(records: list) -> pd.DataFrame:
    df = pd.json_normalize(records, sep=".")
    for c in df.columns:
        if df[c].map(lambda v: isinstance(v, (list, dict))).any():
            df[c] = df[c].map(lambda v: json.dumps(v, sort_keys=True) if isinstance(v, (list, dict)) else v)
    return df


def _json_tables(obj, base: str, label: str, fmt: str) -> list[Table]:
    if isinstance(obj, list):
        recs = [r for r in obj if isinstance(r, dict)]
        if not recs:
            return [Table(name=base, df=pd.DataFrame({"value": obj}), source_file=label, format=fmt)]
        return [Table(name=base, df=coerce_frame(_flatten_records(recs)), source_file=label, format=fmt)]
    if isinstance(obj, dict):
        lists = {k: v for k, v in obj.items() if isinstance(v, list) and v and all(isinstance(r, dict) for r in v)}
        if lists:
            return [Table(name=slug(k), df=coerce_frame(_flatten_records(v)), source_file=label, format=fmt)
                    for k, v in lists.items()]
        return [Table(name=base, df=coerce_frame(_flatten_records([obj])), source_file=label, format=fmt)]
    raise ReadError(f"{label}: JSON holds neither records nor an object")


def _is_fhir(obj) -> bool:
    if isinstance(obj, dict):
        return obj.get("resourceType") == "Bundle" or (isinstance(obj.get("resourceType"), str) and "id" in obj)
    if isinstance(obj, list):
        return bool(obj) and all(isinstance(r, dict) and isinstance(r.get("resourceType"), str) for r in obj[:50])
    return False


# ------------------------------------------------------------------------------------------------------------------
# FHIR R4
# ------------------------------------------------------------------------------------------------------------------
def _fhir_resources(objs: list) -> list[dict]:
    out = []
    for o in objs:
        if isinstance(o, list):
            out.extend(_fhir_resources(o))
        elif isinstance(o, dict) and o.get("resourceType") == "Bundle":
            for e in o.get("entry", []) or []:
                r = e.get("resource")
                if isinstance(r, dict):
                    r = dict(r)
                    r.setdefault("_fullUrl", e.get("fullUrl"))
                    out.extend(_fhir_resources([r]) if r.get("resourceType") == "Bundle" else [r])
        elif isinstance(o, dict) and o.get("resourceType"):
            out.append(o)
    return out


def _codings(cc: dict | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if not isinstance(cc, dict):
        return out
    for c in cc.get("coding", []) or []:
        sys_ = FHIR_SYSTEMS.get(str(c.get("system", "")).rstrip("/"))
        if sys_ and c.get("code") and sys_ not in out:
            out[sys_] = str(c["code"]).strip()
    return out


def _cc_text(cc: dict | None) -> str | None:
    if not isinstance(cc, dict):
        return None
    if cc.get("text"):
        return str(cc["text"])
    for c in cc.get("coding", []) or []:
        if c.get("display"):
            return str(c["display"])
    return None


def fhir_tables(objs: list, label: str) -> tuple[list[Table], list[str]]:
    """FHIR R4 resources (Bundles or bare resources) -> long tables with one column per code system."""
    res = _fhir_resources(objs)
    if not res:
        raise ReadError(f"{label}: no FHIR resources found")
    url_to_pid = {}
    for r in res:
        if r.get("resourceType") == "Patient":
            for k in (r.get("_fullUrl"), f"Patient/{r.get('id')}"):
                if k:
                    url_to_pid[k] = str(r.get("id"))
    meds = {f"Medication/{r.get('id')}": r for r in res if r.get("resourceType") == "Medication"}
    meds.update({r.get("_fullUrl"): r for r in res if r.get("resourceType") == "Medication" and r.get("_fullUrl")})

    def pid(r):
        ref = (r.get("subject") or r.get("patient") or {}).get("reference")
        if not ref:
            return None
        return url_to_pid.get(ref) or ref.split("/")[-1].replace("urn:uuid:", "")

    pats, conds, obs, rx = [], [], [], []
    other: dict[str, int] = {}
    for r in res:
        t = r.get("resourceType")
        if t == "Patient":
            addr = (r.get("address") or [{}])[0]
            pats.append({"patient_id": str(r.get("id")), "gender": r.get("gender"), "birth_date": r.get("birthDate"),
                         "deceased": r.get("deceasedBoolean", r.get("deceasedDateTime")),
                         "postal_code": addr.get("postalCode"), "state": addr.get("state")})
        elif t == "Condition":
            c = _codings(r.get("code"))
            conds.append({"patient_id": pid(r), "icd10cm": c.get("icd10cm") or c.get("icd10"), "snomed": c.get("snomed"),
                          "condition_text": _cc_text(r.get("code")),
                          "clinical_status": ((r.get("clinicalStatus") or {}).get("coding") or [{}])[0].get("code"),
                          "onset_date": r.get("onsetDateTime"), "recorded_date": r.get("recordedDate")})
        elif t == "Observation":
            c = _codings(r.get("code"))
            vq = r.get("valueQuantity") or {}
            obs.append({"patient_id": pid(r), "loinc": c.get("loinc"), "snomed": c.get("snomed"),
                        "observation_text": _cc_text(r.get("code")), "value": vq.get("value"),
                        "unit": vq.get("unit") or vq.get("code"),
                        "value_code": _cc_text(r.get("valueCodeableConcept")) or r.get("valueString"),
                        "category": _cc_text((r.get("category") or [None])[0]),
                        "effective_date": r.get("effectiveDateTime") or (r.get("effectivePeriod") or {}).get("start")})
        elif t in ("MedicationStatement", "MedicationRequest"):
            cc = r.get("medicationCodeableConcept")
            if cc is None and (r.get("medicationReference") or {}).get("reference"):
                cc = (meds.get(r["medicationReference"]["reference"]) or {}).get("code")
            c = _codings(cc)
            rx.append({"patient_id": pid(r), "rxnorm": c.get("rxnorm"), "medication_text": _cc_text(cc),
                       "resource_type": t, "status": r.get("status"),
                       "start_date": (r.get("effectiveDateTime") or (r.get("effectivePeriod") or {}).get("start")
                                      or r.get("authoredOn"))})
        elif t != "Medication":
            other[t] = other.get(t, 0) + 1
    tables, warns = [], []
    spec = [("patients", pats, "fhir:Patient", {}), ("conditions", conds, "fhir:Condition",
                                                     {"icd10cm": "icd10cm", "snomed": "snomed"}),
            ("observations", obs, "fhir:Observation", {"loinc": "loinc", "snomed": "snomed"}),
            ("medications", rx, "fhir:Medication", {"rxnorm": "rxnorm"})]
    labels = {"icd10cm": "Condition.code coding (ICD-10-CM)", "snomed": "code coding (SNOMED CT)",
              "loinc": "Observation.code coding (LOINC)", "rxnorm": "medication coding (RxNorm)",
              "value": "Observation.valueQuantity.value", "unit": "Observation.valueQuantity.unit",
              "patient_id": "FHIR Patient.id"}
    for name, rows, kind, systems in spec:
        if not rows:
            continue
        df = pd.DataFrame(rows).dropna(axis=1, how="all")
        tables.append(Table(name=name, df=coerce_frame(df), source_file=label, format="fhir", kind=kind,
                            labels={k: v for k, v in labels.items() if k in df.columns},
                            code_systems={k: v for k, v in systems.items() if k in df.columns}))
    if other:
        warns.append(f"{label}: FHIR resource types not converted: "
                     + ", ".join(f"{k} ({v})" for k, v in sorted(other.items())))
    return tables, warns


# ------------------------------------------------------------------------------------------------------------------
# REDCap
# ------------------------------------------------------------------------------------------------------------------
_DICT_COLS = {"variable_field_name": "field", "field_name": "field", "form_name": "form", "field_type": "type",
              "field_label": "label", "choices_calculations_or_slider_labels": "choices",
              "select_choices_or_calculations": "choices", "field_note": "note",
              "text_validation_type_or_show_slider_number": "validation", "identifier": "identifier"}


def is_redcap_dictionary(df: pd.DataFrame) -> bool:
    cols = {slug(c) for c in df.columns}
    return ("variable_field_name" in cols or "field_name" in cols) and "field_label" in cols and (
        "field_type" in cols) and ("choices_calculations_or_slider_labels" in cols or "select_choices_or_calculations" in cols)


def _strip_html(s) -> str:
    s = re.sub(r"<[^>]+>", " ", str(s or ""))
    return re.sub(r"\s+", " ", s).strip()


def parse_choices(text) -> dict[str, str]:
    """'1, Yes | 0, No' -> {'1': 'Yes', '0': 'No'}."""
    out = {}
    if not isinstance(text, str) or not text.strip():
        return out
    for part in text.split("|"):
        if "," in part:
            code, lab = part.split(",", 1)
            out[code.strip()] = _strip_html(lab)
    return out


def parse_redcap_dictionary(df: pd.DataFrame) -> list[dict]:
    ren = {c: _DICT_COLS[slug(c)] for c in df.columns if slug(c) in _DICT_COLS}
    d = df.rename(columns=ren)
    fields = []
    for _, r in d.iterrows():
        name = str(r.get("field", "") or "").strip()
        if not name or name == "nan":
            continue
        ftype = str(r.get("type", "") or "").strip().lower()
        if ftype == "descriptive":
            continue
        ch = {}
        if ftype in ("radio", "dropdown", "checkbox"):
            ch = parse_choices(r.get("choices"))
        elif ftype == "yesno":
            ch = {"1": "Yes", "0": "No"}
        elif ftype == "truefalse":
            ch = {"1": "True", "0": "False"}
        note = r.get("note")
        val = r.get("validation")
        fields.append({"field": name, "form": str(r.get("form", "") or "").strip(), "type": ftype,
                       "label": _strip_html(r.get("label")), "choices": ch,
                       "note": _strip_html(note) if isinstance(note, str) else "",
                       "validation": str(val).strip() if isinstance(val, str) else "",
                       "identifier": str(r.get("identifier", "") or "").strip().lower() in ("y", "yes", "1")})
    return fields


def _annotate_redcap(t: Table, fields: list[dict]) -> None:
    by_name = {f["field"]: f for f in fields}
    forms = {f["form"] for f in fields if f["form"]}
    for c in t.df.columns:
        f = by_name.get(c)
        if f is None and "___" in c:
            base, code = c.split("___", 1)
            fb = by_name.get(base)
            if fb is not None and fb["type"] == "checkbox":
                lab = fb["choices"].get(code) or fb["choices"].get(code.replace("_", "-")) or code
                t.labels[c] = f"{fb['label']} (choice={lab})"
                t.choices[c] = {"0": "Unchecked", "1": "Checked"}
                t.field_types[c] = "checkbox"
                t.forms[c] = fb["form"]
                if fb["identifier"]:
                    t.identifier_fields.add(c)
                continue
        if f is not None:
            if f["label"]:
                t.labels[c] = f["label"]
            if f["choices"]:
                t.choices[c] = f["choices"]
            t.field_types[c] = f["type"]
            t.forms[c] = f["form"]
            if f["note"]:
                t.notes[c] = f["note"]
            if f["validation"]:
                t.validation[c] = f["validation"]
            if f["identifier"]:
                t.identifier_fields.add(c)
        elif c.endswith("_complete") and c[:-9] in forms:
            t.labels[c] = f"Form status: {c[:-9]}"
            t.choices[c] = {"0": "Incomplete", "1": "Unverified", "2": "Complete"}
            t.structural.add(c)
            t.forms[c] = c[:-9]
        if c in REDCAP_STRUCTURAL:
            t.structural.add(c)


def is_redcap_data(df: pd.DataFrame) -> bool:
    return any(c in df.columns for c in ("redcap_event_name", "redcap_repeat_instrument", "redcap_repeat_instance"))


def split_redcap(t: Table, fields: list[dict] | None) -> list[Table]:
    """A REDCap data table -> base records table + one table per repeating instrument."""
    df = t.df
    if fields:
        _annotate_redcap(t, fields)
    for c in df.columns:
        if c in REDCAP_STRUCTURAL:
            t.structural.add(c)
    t.kind = "redcap_records"
    if "redcap_repeat_instrument" not in df.columns or df["redcap_repeat_instrument"].isna().all():
        return [t]
    id_col = (fields[0]["field"] if fields else None) or df.columns[0]
    if id_col not in df.columns:
        id_col = df.columns[0]
    keys = [c for c in (id_col, "redcap_event_name", "redcap_repeat_instance") if c in df.columns]
    rep = df["redcap_repeat_instrument"]
    base = df[rep.isna()]
    out = []
    form_cols: dict[str, list[str]] = {}
    for inst in sorted(rep.dropna().astype(str).unique()):
        rows = df[rep.astype(str) == inst]
        if t.forms:
            cols = [c for c in df.columns if t.forms.get(c) == inst and c not in keys]
        else:
            cols = [c for c in df.columns if c not in keys and c not in REDCAP_STRUCTURAL and rows[c].notna().any()
                    and base[c].isna().all()]
        form_cols[inst] = cols
        sub = rows[keys + cols].reset_index(drop=True)
        out.append(Table(name=slug(inst), df=sub, source_file=t.source_file, format="redcap",
                         kind=f"redcap_repeating:{inst}",
                         labels={c: t.labels[c] for c in sub.columns if c in t.labels},
                         choices={c: t.choices[c] for c in sub.columns if c in t.choices},
                         field_types={c: t.field_types[c] for c in sub.columns if c in t.field_types},
                         forms={c: t.forms[c] for c in sub.columns if c in t.forms},
                         notes={c: t.notes[c] for c in sub.columns if c in t.notes},
                         validation={c: t.validation[c] for c in sub.columns if c in t.validation},
                         identifier_fields={c for c in t.identifier_fields if c in sub.columns},
                         structural={c for c in t.structural if c in sub.columns} | {"redcap_repeat_instance"} & set(sub.columns)))
    rep_only = {c for cols in form_cols.values() for c in cols if base[c].isna().all()}
    base_cols = [c for c in df.columns if c not in rep_only and c not in ("redcap_repeat_instrument", "redcap_repeat_instance")]
    t.df = base[base_cols].reset_index(drop=True)
    t.name = "records"
    for d in (t.labels, t.choices, t.field_types, t.forms, t.notes, t.validation):
        for c in list(d):
            if c not in base_cols:
                d.pop(c)
    t.identifier_fields &= set(base_cols)
    t.structural &= set(base_cols)
    return [t, *out]


# ------------------------------------------------------------------------------------------------------------------
# files, folders, archives
# ------------------------------------------------------------------------------------------------------------------
def read_file(path: Path, label: str | None = None) -> tuple[list[Table], list[str], list]:
    """One file -> (tables, warnings, FHIR objects found). FHIR JSON is returned raw so a folder of bundles can be
    combined (references across files)."""
    label = label or path.name
    suf = path.suffix.lower()
    warns: list[str] = []
    if suf in TEXT_SUFFIXES:
        raw = path.read_bytes()
        if not raw.strip():
            return [], [f"{label}: empty file"], []
        df, meta = _text_table(raw, suf)
        if meta["encoding"] not in ("utf-8-sig",):
            warns.append(f"{label}: decoded as {meta['encoding']}")
        fmt = "tsv" if meta["delimiter"] == "\t" else "csv"
        return [Table(name=slug(path.stem), df=coerce_frame(df), source_file=label, format=fmt)], warns, []
    if suf in EXCEL_SUFFIXES:
        return _read_excel(path, label), warns, []
    if suf in (".parquet", ".pq"):
        try:
            df = pd.read_parquet(path)
        except ImportError as e:  # pragma: no cover
            raise MissingDependency("reading Parquet needs the package 'pyarrow'") from e
        return [Table(name=slug(path.stem), df=coerce_frame(df), source_file=label, format="parquet")], warns, []
    if suf in STAT_SUFFIXES:
        return [_read_stat(path, suf, label)], warns, []
    if suf in JSON_SUFFIXES:
        text, _ = _decode(path.read_bytes())
        if suf in (".jsonl", ".ndjson"):
            objs = [json.loads(ln) for ln in text.splitlines() if ln.strip()]
            fmt = "jsonl"
        else:
            try:
                objs = json.loads(text)
                fmt = "json"
            except json.JSONDecodeError:
                objs = [json.loads(ln) for ln in text.splitlines() if ln.strip()]
                fmt = "jsonl"
        if _is_fhir(objs):
            return [], warns, objs if isinstance(objs, list) else [objs]
        return _json_tables(objs, slug(path.stem), label, fmt), warns, []
    raise ReadError(f"{label}: unsupported file type {suf or '(none)'}; supported: "
                    + ", ".join(sorted(SUPPORTED_SUFFIXES)))


def _unique_names(tables: list[Table], multi_file: bool) -> None:
    seen: dict[str, int] = {}
    for t in tables:
        if multi_file and t.sheet is not None:
            t.name = f"{slug(Path(t.source_file).stem)}_{slug(t.sheet)}"
        n = t.name
        if n in seen:
            seen[n] += 1
            t.name = f"{n}_{seen[n]}"
        else:
            seen[n] = 1


def _read_many(files: list[tuple[Path, str]], input_label: str, container: str | None) -> Dataset:
    tables: list[Table] = []
    warns: list[str] = []
    fhir_objs: list = []
    fhir_files: list[str] = []
    for p, lab in files:
        if p.suffix.lower() == ".zip":
            warns.append(f"{lab}: nested ZIP archives are not opened")
            continue
        if p.suffix.lower() not in SUPPORTED_SUFFIXES:
            warns.append(f"{lab}: skipped (unsupported file type {p.suffix or '(none)'})")
            continue
        try:
            t, w, f = read_file(p, lab)
        except ReadError as e:
            warns.append(str(e))
            continue
        tables += t
        warns += w
        if f:
            fhir_objs += f
            fhir_files.append(lab)
    if fhir_objs:
        ft, fw = fhir_tables(fhir_objs, ", ".join(fhir_files))
        tables += ft
        warns += fw
    if not tables:
        raise ReadError(f"{input_label}: no readable table found" + (f" ({'; '.join(warns[:5])})" if warns else ""))
    return assemble(tables, input_label, container, warns, multi_file=len(files) > 1)


def assemble(tables: list[Table], input_label: str, container: str | None, warns: list[str],
             multi_file: bool = False) -> Dataset:
    """Detect special layouts (REDCap, OMOP, FHIR) across the tables read and decide the dataset format."""
    dictionary = None
    dicts = [t for t in tables if is_redcap_dictionary(t.df)]
    fields = None
    if dicts:
        d = dicts[0]
        fields = parse_redcap_dictionary(d.df)
        dictionary = {"source": d.source_file, "n_fields": len(fields), "kind": "redcap_data_dictionary",
                      "forms": sorted({f["form"] for f in fields if f["form"]})}
        tables = [t for t in tables if t not in dicts]
        if len(dicts) > 1:
            warns.append(f"{len(dicts)} REDCap data dictionaries found; used {d.source_file}")
    fmt = None
    if fields is not None or any(is_redcap_data(t.df) for t in tables):
        fmt = "redcap"
        names = {f["field"] for f in fields or []}
        out = []
        for t in tables:
            overlap = len(set(t.df.columns) & names) if names else 0
            if is_redcap_data(t.df) or (names and overlap >= max(2, 0.3 * len(t.df.columns))):
                out += split_redcap(t, fields)
            else:
                out.append(t)
        tables = out
    stems = {Path(t.source_file).stem.lower() for t in tables}
    omop_hits = stems & OMOP_TABLES
    if fmt is None and ("person" in omop_hits and len(omop_hits) >= 2 or len(omop_hits) >= 3):
        fmt = "omop"
        for t in tables:
            st = Path(t.source_file).stem.lower()
            if st in OMOP_TABLES and t.sheet is None:
                t.kind = f"omop:{st}"
                t.name = st
    if fmt is None and tables and all(t.format == "fhir" for t in tables):
        fmt = "fhir"
    elif fmt is None and any(t.format == "fhir" for t in tables):
        fmt = "mixed"
    if fmt is None:
        kinds = {t.format for t in tables}
        fmt = kinds.pop() if len(kinds) == 1 else "mixed"
    _unique_names(tables, multi_file and fmt not in ("redcap", "omop", "fhir"))
    return Dataset(input=input_label, format=fmt, container=container, tables=tables, dictionary=dictionary,
                   warnings=warns)


def _safe_members(zf: zipfile.ZipFile) -> list[zipfile.ZipInfo]:
    out = []
    for info in zf.infolist():
        n = info.filename
        if info.is_dir() or n.startswith("/") or ".." in Path(n).parts or "__MACOSX" in n or Path(n).name.startswith("."):
            continue
        out.append(info)
    return out


def read_input(path: str | Path) -> Dataset:
    """Read a file, a ZIP archive or a folder into a `Dataset` of named tables (module docstring)."""
    p = Path(path)
    if not p.exists():
        raise ReadError(f"{p}: no such file or folder")
    if p.is_dir():
        files = sorted(f for f in p.rglob("*") if f.is_file() and not any(part.startswith(".") for part in f.relative_to(p).parts)
                       and "__MACOSX" not in f.parts)
        return _read_many([(f, str(f.relative_to(p))) for f in files], str(p), "folder")
    if p.suffix.lower() == ".zip" or zipfile.is_zipfile(p) and p.suffix.lower() not in EXCEL_SUFFIXES:
        try:
            zf = zipfile.ZipFile(p)
        except zipfile.BadZipFile as e:
            raise ReadError(f"{p}: not a valid ZIP archive") from e
        with zf, tempfile.TemporaryDirectory(prefix="measure_it_agent_") as tmp:
            files = []
            for info in _safe_members(zf):
                dst = Path(tmp) / info.filename
                dst.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info) as src, open(dst, "wb") as fh:
                    fh.write(src.read())
                files.append((dst, f"{p.name}!{info.filename}"))
            return _read_many(sorted(files, key=lambda x: x[1]), str(p), "zip")
    tables, warns, fhir_objs = read_file(p)
    if fhir_objs:
        tables, fw = fhir_tables(fhir_objs, p.name)
        warns += fw
    if not tables:
        raise ReadError(f"{p}: no readable table found" + (f" ({'; '.join(warns)})" if warns else ""))
    return assemble(tables, str(p), None, warns)
