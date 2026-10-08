"""`measure-it byod validate <dir>`: schema checks and the refusals that keep identifying or geographic content out.

The engine refuses (non-zero exit, one message per problem, the offending VALUE never echoed):
  * any column whose name locates a person or a place (ZIP / postcode, county, FIPS, ZCTA, latitude / longitude,
    address, city, state), or identifies a person (names, date of birth, MRN, phone, e-mail, SSN, insurance ids),
    or carries an exact date or an exact age (supply `day_index` / `year` and `age_band` instead);
  * text values that look like e-mail addresses, phone numbers, SSNs, exact dates, ZIP codes, street addresses or
    coordinate pairs, in any text column of any file;
  * participant ids that look like MRNs (6+ digits, or an MRN / patient prefix), SSNs, e-mails, phone numbers,
    dates or personal names;
  * ages above 89 not top-coded as `90+`;
  * a manifest that does not confirm the data are de-identified and that this use is permitted;
  * fewer than 10 cases or 10 controls (configurable, never below 5) for the primary label;
  * undeclared columns in participants.csv, non-numeric feature columns, malformed LOINC / ICD-10-CM codes,
    participant ids in a data file that are not in participants.csv, and data files the contract does not name.

These checks are a safety net, not a de-identification certificate: the data owner remains responsible for
de-identification, IRB / data-use approval and the right to use the data here (docs/BRING_YOUR_OWN_DATA.md).
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ..config import load_config
from ..validate import GEO_COLUMN_RE
from . import common as K

# ------------------------------------------------------------------------------------------------------------------
# manifest schema (also shipped as templates/byod/manifest.schema.json; a test keeps the two identical)
# ------------------------------------------------------------------------------------------------------------------
_NAME = {"type": "string", "pattern": r"^[a-z][a-z0-9_]{0,40}$"}
MANIFEST_SCHEMA: dict = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://measure-it.local/byod/manifest.schema.json",
    "title": "measure-it bring-your-own-data manifest",
    "type": "object",
    "additionalProperties": False,
    "required": ["schema_version", "dataset_id", "title", "owner", "licence", "data_use", "labels", "comparator",
                 "measurement", "primary_analysis"],
    "properties": {
        "schema_version": {"const": K.SCHEMA_VERSION},
        "dataset_id": {"type": "string", "pattern": K.MANIFEST_ID_RE.pattern,
                       "description": "lower-case id; the engine dataset id becomes byod_<dataset_id>"},
        "title": {"type": "string", "minLength": 3},
        "owner": {"type": "string", "minLength": 2, "description": "organisation responsible for the data"},
        "description": {"type": "string"},
        "licence": {"type": "string", "minLength": 3, "description": "licence or data-use statement"},
        "synthetic": {"type": "boolean", "description": "true only for generated (not real) data"},
        "demo": {"type": "boolean", "description": "true: the performance record is excluded from rankings unless "
                                                   "deploy --demo"},
        "data_use": {
            "type": "object", "additionalProperties": False,
            "required": ["deidentified", "use_permitted", "statement"],
            "properties": {
                "deidentified": {"const": True},
                "use_permitted": {"const": True},
                "statement": {"type": "string", "minLength": 10},
                "irb_or_dua_reference": {"type": "string"}}},
        "labels": {
            "type": "array", "minItems": 1,
            "items": {"type": "object", "additionalProperties": False,
                      "required": ["column", "condition", "definition", "label_basis"],
                      "properties": {"column": _NAME, "condition": {"type": "string"},
                                     "definition": {"type": "string", "minLength": 10},
                                     "label_basis": {"enum": sorted(K.LABEL_BASES)}}}},
        "comparator": {
            "type": "object", "additionalProperties": False, "required": ["type", "description"],
            "properties": {"type": {"enum": sorted(K.COMPARATOR_TYPES)}, "description": {"type": "string"}}},
        "measurement": {
            "type": "object", "additionalProperties": False, "required": ["class"],
            "properties": {"class": {"type": "string"}, "bundle": {"type": "string"}, "device": {"type": "string"}}},
        "devices": {
            "type": "object",
            "additionalProperties": {
                "type": "object", "additionalProperties": False,
                "properties": {
                    "measurement_class": {"type": "string"}, "description": {"type": "string"},
                    "adapter": {"type": "object", "additionalProperties": False, "required": ["module", "class"],
                                "properties": {"module": {"type": "string"}, "class": {"type": "string"},
                                               "register_as": _NAME}}}}},
        "id_hashing": {"enum": ["none", "salted_sha256"]},
        "min_group_size": {"type": "object", "additionalProperties": False,
                           "properties": {"cases": {"type": "integer", "minimum": 5},
                                          "controls": {"type": "integer", "minimum": 5}}},
        "primary_analysis": {
            "type": "object", "additionalProperties": False, "required": ["label", "feature_blocks"],
            "properties": {
                "label": _NAME,
                "feature_blocks": {"type": "array", "minItems": 1, "items": {"type": "string"}},
                "combined_blocks": {"type": "array", "items": {"type": "string"}},
                "specificity": {"type": "number", "exclusiveMinimum": 0.5, "maximum": 0.99},
                "covariates": {"type": "array", "items": {"enum": ["age_band", "sex"]}},
                "n_splits": {"type": "integer", "minimum": 2, "maximum": 10},
                "n_repeats": {"type": "integer", "minimum": 1, "maximum": 50},
                "n_permutations": {"type": "integer", "minimum": 0, "maximum": 5000},
                "n_bootstrap": {"type": "integer", "minimum": 100, "maximum": 10000},
                "C": {"type": "number", "exclusiveMinimum": 0}}},
        "omics": {
            "type": "object",
            "description": "optional, per omics_<layer> file: the assay platform; platform nmr_nightingale with units "
                           "nightingale_standard lets clinical-chemistry-equivalent measures meet NHANES labs as LOINC "
                           "(configs/harmonize_nightingale_loinc.yaml)",
            "additionalProperties": {
                "type": "object", "additionalProperties": False,
                "properties": {"platform": _NAME, "units": {"type": "string"}, "description": {"type": "string"}}}},
        "subgroup_analysis": {
            "type": "object", "additionalProperties": False,
            "description": "PRE-SPECIFIED device-defined subgroup inside a label's cases (`measure-it byod subgroup`)",
            "required": ["name", "within_label", "device_positive"],
            "properties": {
                "name": _NAME,
                "description": {"type": "string"},
                "within_label": _NAME,
                "device_positive": {
                    "type": "object", "additionalProperties": False,
                    "properties": {
                        "column": _NAME,
                        "block": {"type": "string", "pattern": r"^device_[a-z][a-z0-9_]{0,40}$"},
                        "feature": {"type": "string"},
                        "threshold": {"type": "number"},
                        "direction": {"enum": [">=", ">", "<=", "<"]}},
                    "oneOf": [{"required": ["column"]}, {"required": ["block", "feature", "threshold"]}]},
                "ehr_domains": {"type": "array", "minItems": 1, "uniqueItems": True,
                                "items": {"enum": ["condition", "demographic", "drug", "measurement", "survey"]}},
                "max_features": {"type": "integer", "minimum": 1, "maximum": 100},
                "min_feature_participants": {"type": "integer", "minimum": 3},
                "exclude_feature_keys": {"type": "array", "items": {"type": "string"}},
                "n_splits": {"type": "integer", "minimum": 2, "maximum": 10},
                "n_repeats": {"type": "integer", "minimum": 1, "maximum": 50},
                "n_permutations": {"type": "integer", "minimum": 0, "maximum": 5000},
                "n_bootstrap": {"type": "integer", "minimum": 100, "maximum": 10000},
                "C": {"type": "number", "exclusiveMinimum": 0}}},
    },
}
DEFAULT_MIN_GROUP = 10
HARD_MIN_GROUP = 5

# ------------------------------------------------------------------------------------------------------------------
# identifying / geographic content
# ------------------------------------------------------------------------------------------------------------------
GEO_TOKENS = {"zip", "zipcode", "zip5", "zip3", "postcode", "postal", "county", "fips", "zcta", "latitude", "longitude",
              "lat", "lon", "lng", "gps", "geoid", "geocode", "address", "street", "city", "township", "tract",
              "census", "municipality", "neighborhood", "neighbourhood", "region", "country", "geolocation"}
GEO_EXACT = {"state", "state_code", "us_state", "state_abbr", "state_name", "province", "location", "residence",
             "home", "site_location", "hospital", "clinic", "facility", "site"}
ID_TOKENS = {"mrn", "ssn", "dob", "birthdate", "birthday", "phone", "telephone", "mobile", "fax", "email", "mail",
             "surname", "medicare", "medicaid", "passport", "license", "licence", "ip", "url", "npi", "vin", "insurance",
             "beneficiary", "account", "biometric", "photo", "face", "fingerprint", "initials"}
ID_SUBSTRINGS = ("birth", "social_security", "socialsecurity", "medical_record", "record_number", "health_plan",
                 "first_name", "last_name", "given_name", "family_name", "full_name", "patient_name", "middle_name",
                 "maiden", "firstname", "lastname", "fullname")
DATE_TOKENS = {"date", "datetime", "timestamp", "dod", "admission", "admitted", "discharge", "discharged", "death",
               "died", "deceased"}
AGE_EXACT = {"age", "age_years", "age_yrs", "age_at_visit", "age_at_enrollment", "age_at_enrolment", "agey", "age_y"}

VALUE_PATTERNS = {
    "e-mail address": re.compile(r"[^@\s]+@[^@\s]+\.[A-Za-z]{2,}"),
    "social security number": re.compile(r"\b\d{3}-\d{2}-\d{4}\b"),
    "phone number": re.compile(r"(?:\+?1[\s.-]?)?\(?\b\d{3}\)?[\s.-]\d{3}[\s.-]\d{4}\b"),
    "exact date": re.compile(r"\b(?:19|20)\d{2}[-/.](?:0?[1-9]|1[0-2])[-/.](?:0?[1-9]|[12]\d|3[01])\b"
                             r"|\b(?:0?[1-9]|1[0-2])[-/.](?:0?[1-9]|[12]\d|3[01])[-/.](?:19|20)?\d{2}\b"
                             r"|\b\d{1,2}[ -](?:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)[a-z]*[ -]\d{2,4}\b",
                             re.IGNORECASE),
    "ZIP code": re.compile(r"^\s*\d{5}(?:-\d{4})?\s*$"),
    "street address": re.compile(r"\b\d{1,6}\s+(?:[A-Za-z0-9.]+\s){0,4}(?:street|st|avenue|ave|road|rd|boulevard|"
                                 r"blvd|lane|ln|drive|dr|court|ct|way|place|pl|terrace|parkway|pkwy|highway|hwy)\b\.?",
                                 re.IGNORECASE),
    "coordinate pair": re.compile(r"-?\d{1,3}\.\d{3,}\s*[,;]\s*-?\d{1,3}\.\d{3,}"),
}
ID_PATTERNS = {
    "MRN-like number (6 or more digits)": re.compile(r"^\d{6,}$"),
    "MRN / patient-number prefix": re.compile(r"^(?:mrn|mr|pt|pat|patient|med|hosp|chart|acct)[-_ #:.]?\d+$",
                                             re.IGNORECASE),
    "personal name": re.compile(r"^[A-Z][a-z]+[ _,]+[A-Z][a-z]+$"),
}
ID_CHARSET = re.compile(r"^[A-Za-z0-9._-]{1,64}$")
LOINC_RE = re.compile(r"^\d{1,7}-\d$")
ICD10CM_RE = re.compile(r"^[A-Z]\d[0-9A-Z](?:\.?[0-9A-Z]{1,4})?$")
RXNORM_RE = re.compile(r"^\d{1,9}$")
AGE_BAND_RE = re.compile(r"^(\d{1,2})-(\d{1,2})$")
AGE_TOP_RE = re.compile(r"^(\d{2})\+$")
SEX_VALUES = {"female": "female", "f": "female", "male": "male", "m": "male", "other": "other",
              "unknown": "unknown", "": "unknown"}
FEATURE_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.\-]{0,80}$")
LONG_KEYS = ("day_index",)


def _tokens(col: str) -> list[str]:
    return [t for t in re.split(r"[^a-z0-9]+", str(col).lower()) if t]


def column_problem(col: str) -> str | None:
    """Why a column name is refused (None = allowed). Checked in every file."""
    c = str(col).strip()
    low = re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_")
    toks = set(_tokens(c))
    if GEO_COLUMN_RE.search(c) or toks & GEO_TOKENS or low in GEO_EXACT:
        return "locates a person or a place (geographic column); person-level data never carry geography"
    if low in ("name", "names", "patient", "person") or any(s in low for s in ID_SUBSTRINGS) or toks & ID_TOKENS:
        return "identifies a person (name, date of birth, MRN, phone, e-mail, SSN or similar)"
    if toks & DATE_TOKENS:
        return "carries a calendar date; supply a relative day number `day_index` (or `year`) instead"
    if low in AGE_EXACT or re.fullmatch(r"age(_?in)?_?(years?|yrs?|months?)", low):
        return "is an exact age; supply `age_band` (e.g. 40-49; 90+ top-coded) in participants.csv instead"
    return None


def _mask(v: str) -> str:
    return f"<{len(str(v))} characters, starting {str(v)[:1]!r}>"


def value_problems(s: pd.Series, col: str, file: str, max_examples: int = 1) -> list[str]:
    """Refusals for the text values of one column (numeric columns are not scanned). Values are never echoed."""
    vals = s.dropna().astype(str)
    vals = vals[vals.str.strip() != ""]
    if vals.empty:
        return []
    numeric = pd.to_numeric(vals, errors="coerce").notna()
    text = vals[~numeric]
    out = []
    for name, rx in VALUE_PATTERNS.items():
        pool = text
        if name == "ZIP code":
            # in a text column any 5-digit(-4) value is a ZIP candidate; in a numeric column only leading-zero ones
            # ('02139' is never a numeric feature value as supplied)
            pool = vals if not numeric.all() else vals[vals.str.fullmatch(r"\s*0\d{4}(?:-\d{4})?\s*")]
        if pool.empty:
            continue
        hit = pool[pool.map(lambda v: bool(rx.search(v)))]
        if len(hit):
            out.append(f"{file}: column {col!r}: {len(hit)} value(s) look like: {name} "
                       f"(e.g. {_mask(hit.iloc[0])}); remove or transform it before supplying the data")
    return out


def id_problems(ids: pd.Series, file: str) -> list[str]:
    out = []
    s = ids.astype("string")
    if s.isna().any() or (s.str.strip() == "").any():
        out.append(f"{file}: {int((s.isna() | (s.str.strip() == '')).sum())} empty participant_id value(s)")
    s = s.dropna()
    bad_chars = s[~s.map(lambda v: bool(ID_CHARSET.fullmatch(str(v))))]
    if len(bad_chars):
        out.append(f"{file}: {len(bad_chars)} participant_id value(s) use characters outside A-Z a-z 0-9 . _ - "
                   f"or exceed 64 characters (e.g. {_mask(bad_chars.iloc[0])})")
    for name, rx in {**ID_PATTERNS, **{k: v for k, v in VALUE_PATTERNS.items() if k != "ZIP code"}}.items():
        hit = s[s.map(lambda v: bool(rx.search(str(v))))]
        if len(hit):
            out.append(f"{file}: {len(hit)} participant_id value(s) look like: {name} (e.g. {_mask(hit.iloc[0])}): "
                       "use study codes, not medical-record or other identifying numbers")
    return out


def age_band_problem(v: str) -> str | None:
    v = str(v).strip()
    m = AGE_BAND_RE.fullmatch(v)
    if m:
        lo, hi = int(m.group(1)), int(m.group(2))
        if hi < lo or hi - lo < 4:
            return "age bands must span at least 5 years"
        if hi > 89:
            return "ages above 89 must be top-coded as 90+"
        return None
    m = AGE_TOP_RE.fullmatch(v)
    if m:
        return None if int(m.group(1)) <= 90 else "the top-coded band must start at 90 or below (90+)"
    return "not an age band like 40-49 or 90+"


# ------------------------------------------------------------------------------------------------------------------
# the report
# ------------------------------------------------------------------------------------------------------------------

@dataclass
class Report:
    root: Path
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    will_create: list[str] = field(default_factory=list)
    info: dict = field(default_factory=dict)
    manifest: dict = field(default_factory=dict)
    inputs: K.Inputs | None = None
    frames: dict = field(default_factory=dict)       # file key -> DataFrame as read (text)
    resolved: dict = field(default_factory=dict)     # label column -> canonical condition id, class, bundle, ...

    @property
    def ok(self) -> bool:
        return not self.errors

    def err(self, msg: str) -> None:
        self.errors.append(msg)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)

    def to_dict(self) -> dict:
        return {"ok": self.ok, "root": str(self.root.name), "errors": self.errors, "warnings": self.warnings,
                "will_create": self.will_create, "info": self.info, "resolved": self.resolved}

    def text(self) -> str:
        lines = [f"byod validate {self.root.name}: {'OK' if self.ok else 'REFUSED'} "
                 f"({len(self.errors)} error(s), {len(self.warnings)} warning(s))"]
        for e in self.errors:
            lines.append(f"  ERROR   {e}")
        for w in self.warnings:
            lines.append(f"  warning {w}")
        if self.info:
            lines.append("  summary:")
            for k, v in self.info.items():
                lines.append(f"    {k}: {v}")
        if self.ok and self.will_create:
            lines.append("  will create:")
            lines += [f"    {w}" for w in self.will_create]
        return "\n".join(lines)


def resolve_condition(text: str) -> tuple[str | None, str]:
    """Label condition -> canonical condition id (or a scored condition set) via measure_it.ontology.normalize."""
    from ..ontology.normalize import normalize_condition
    ids = {c["id"] for c in load_config("conditions")["conditions"]}
    try:
        from ..scoring.opportunity import condition_sets
        sets = condition_sets()
    except Exception:  # noqa: BLE001 - sets need configs only; never fatal for validation
        sets = {}
    t = str(text).strip()
    if t in ids:
        return t, "canonical condition id"
    if t in sets:
        return t, "scored condition set id"
    for sid, s in sets.items():
        if t.lower() in {a.lower() for a in s.get("aliases", [])} or t.lower() == str(s.get("label", "")).lower():
            return sid, f"condition set alias ({s.get('label')})"
    r = normalize_condition(t)
    cands = sorted({m.get("canonical_condition_id") for m in r.get("matches", []) or [] if m.get("canonical_condition_id")})
    if r.get("status") == "matched" and len(cands) == 1:
        return cands[0], f"measure_it.ontology.normalize: {r['matches'][0].get('match_reason', '')}"
    return None, (f"normalize_condition status {r.get('status')!r}"
                  + (f", candidates {cands}" if cands else "") + ": give a canonical condition id or an exact name")


def _measurement_classes() -> dict:
    return {m["id"]: m for m in load_config("measurements")["measurement_classes"]}


def check_manifest(rep: Report) -> None:
    m = rep.manifest
    try:
        import jsonschema
        v = jsonschema.Draft202012Validator(MANIFEST_SCHEMA)
        for e in sorted(v.iter_errors(m), key=lambda e: list(e.path)):
            where = "/".join(str(p) for p in e.path) or "(top level)"
            rep.err(f"manifest.yaml {where}: {e.message}")
    except ImportError:   # the explicit checks below still run
        for k in MANIFEST_SCHEMA["required"]:
            if k not in m:
                rep.err(f"manifest.yaml: missing required field {k!r}")
    du = m.get("data_use") or {}
    if du.get("deidentified") is not True:
        rep.err("manifest.yaml data_use.deidentified must be `true`: the engine accepts de-identified data only")
    if du.get("use_permitted") is not True:
        rep.err("manifest.yaml data_use.use_permitted must be `true`: confirm that your IRB / data-use agreement "
                "permits this use on this machine")
    mid = str(m.get("dataset_id", ""))
    if mid and not K.MANIFEST_ID_RE.fullmatch(mid):
        rep.err(f"manifest.yaml dataset_id {mid!r}: use 3-41 lower-case letters, digits or '_' (no '__')")
    if "__" in mid:
        rep.err("manifest.yaml dataset_id must not contain '__' (reserved for <table>__<source> partitions)")
    classes = _measurement_classes()
    meas = m.get("measurement") or {}
    cls = meas.get("class")
    if cls and cls not in classes:
        rep.err(f"manifest.yaml measurement.class {cls!r} is not a measurement class in configs/measurements.yaml "
                f"(one of: {', '.join(sorted(classes))})")
    bundle = meas.get("bundle")
    bundles = load_config("relevance").get("measurement_bundles", {}) or {}
    if bundle:
        if bundle not in bundles:
            rep.err(f"manifest.yaml measurement.bundle {bundle!r} is not a bundle in configs/relevance.yaml "
                    f"(one of: {', '.join(sorted(bundles))})")
        elif cls and cls not in (bundles[bundle] or {}).get("members", []):
            rep.err(f"manifest.yaml measurement.class {cls!r} is not a member of bundle {bundle!r}")
    rep.resolved["measurement_class"] = cls
    rep.resolved["bundles_containing_class"] = sorted(b for b, s in bundles.items() if cls in (s or {}).get("members", []))
    for name, spec in (m.get("devices") or {}).items():
        dc = (spec or {}).get("measurement_class")
        if dc and dc not in classes:
            rep.err(f"manifest.yaml devices.{name}.measurement_class {dc!r} is not in configs/measurements.yaml")
        ad = (spec or {}).get("adapter")
        if ad:
            p = (rep.root / str(ad.get("module", ""))).resolve()
            if not str(p).startswith(str(rep.root.resolve())) or not p.exists() or p.suffix != ".py":
                rep.err(f"manifest.yaml devices.{name}.adapter.module {ad.get('module')!r}: a .py file inside the "
                        "dataset folder is required")
            elif f"class {ad.get('class')}" not in p.read_text(errors="replace"):
                rep.err(f"manifest.yaml devices.{name}.adapter.class {ad.get('class')!r} not defined in {p.name}")
    labels = {}
    for lab in m.get("labels") or []:
        col = lab.get("column")
        cid, why = resolve_condition(lab.get("condition", ""))
        if cid is None:
            rep.err(f"manifest.yaml labels[{col}].condition {lab.get('condition')!r} does not map to one condition: "
                    f"{why}")
        labels[col] = {"condition_id": cid, "mapping": why, "label_basis": lab.get("label_basis"),
                       "definition": lab.get("definition")}
    rep.resolved["labels"] = labels
    if m.get("id_hashing") == "salted_sha256":
        salt = os.environ.get(K.SALT_ENV, "")
        if len(salt) < 16:
            rep.err(f"id_hashing: salted_sha256 needs a salt of >= 16 characters in the environment variable "
                    f"{K.SALT_ENV} (it is never written to disk; keep it with the data owner)")
    if m.get("synthetic") and not m.get("demo"):
        rep.err("a synthetic dataset must also set demo: true (its record never enters default rankings)")


def _read(rep: Report, f: K.InputFile) -> pd.DataFrame | None:
    try:
        df = K.read_frame(f.path)
    except Exception as exc:  # noqa: BLE001 - reported, not raised
        rep.err(f"{f.path.name}: cannot be read ({type(exc).__name__}: {str(exc)[:160]})")
        return None
    df.columns = [str(c).strip() for c in df.columns]
    rep.frames[f.path.name] = df
    return df


def _common_columns(rep: Report, df: pd.DataFrame, fname: str) -> None:
    for c in df.columns:
        why = column_problem(c)
        if why:
            rep.err(f"{fname}: column {c!r} refused: it {why}")
        if df[c].dtype == object or pd.api.types.is_string_dtype(df[c]):
            for p in value_problems(df[c], c, fname):
                rep.err(p)
    dup = df.columns[df.columns.duplicated()].tolist()
    if dup:
        rep.err(f"{fname}: duplicated column names {dup}")


def _numeric_columns(rep: Report, df: pd.DataFrame, fname: str, skip: set[str]) -> list[str]:
    feats = []
    for c in df.columns:
        if c in skip:
            continue
        if not FEATURE_NAME_RE.fullmatch(c):
            rep.err(f"{fname}: feature column name {c[:40]!r} must start with a letter and use letters, digits, "
                    "'_', '.' or '-' (<= 81 characters)")
            continue
        x = pd.to_numeric(df[c], errors="coerce")
        bad = df[c].notna() & x.isna()
        if bad.any():
            rep.err(f"{fname}: feature column {c!r} has {int(bad.sum())} non-numeric value(s); feature files hold "
                    "numbers only (text columns can carry identifying content)")
        feats.append(c)
    return feats


def check_dir(root: str | Path, min_cases: int | None = None, min_controls: int | None = None) -> Report:
    root = Path(root).resolve()
    rep = Report(root=root)
    if not root.is_dir():
        rep.err(f"{root}: not a directory")
        return rep
    inp = K.discover(root)
    rep.inputs = inp
    if not inp.manifest_path.exists():
        rep.err("manifest.yaml missing (copy templates/byod/manifest.yaml)")
        return rep
    try:
        rep.manifest = K.load_manifest(inp.manifest_path)
    except Exception as exc:  # noqa: BLE001
        rep.err(f"manifest.yaml does not parse: {exc}")
        return rep
    check_manifest(rep)
    for p in inp.unknown_data_files:
        rep.err(f"{p.name}: not a file name the contract knows (participants.csv, device_<name>.csv|parquet, "
                "ehr_labs.csv, ehr_diagnoses.csv, ehr_medications.csv, omics_<layer>.csv|parquet, "
                "survey_<instrument>.csv|parquet, self_report_history.csv); remove it from the folder")
    for f in inp.files:
        if f.name and not K.NAME_RE.fullmatch(f.name):
            rep.err(f"{f.path.name}: the part after device_/omics_/survey_ must be lower-case letters, digits or '_'")
    if inp.participants is None:
        rep.err("participants.csv missing")
        return rep
    m = rep.manifest
    label_cols = [lab.get("column") for lab in m.get("labels") or [] if lab.get("column")]
    # ---------------- participants.csv
    P = _read(rep, inp.participants)
    if P is None:
        return rep
    fname = "participants.csv"
    _common_columns(rep, P, fname)
    if "participant_id" not in P.columns:
        rep.err(f"{fname}: participant_id column missing")
        return rep
    sa = m.get("subgroup_analysis") or {}
    dp_col = (sa.get("device_positive") or {}).get("column")
    allowed = {"participant_id", "age_band", "sex", *label_cols, *([dp_col] if dp_col else [])}
    extra = [c for c in P.columns if c not in allowed]
    if extra:
        rep.err(f"{fname}: undeclared column(s) {extra}: participants.csv holds participant_id, the label columns "
                "declared in the manifest, optionally age_band and sex, and the subgroup_analysis device_positive "
                "column if one is declared (features go in device_/ehr_/omics_/survey_ files)")
    for p in id_problems(P["participant_id"], fname):
        rep.err(p)
    if P["participant_id"].duplicated().any():
        rep.err(f"{fname}: {int(P['participant_id'].duplicated().sum())} duplicated participant_id value(s)")
    for c in label_cols:
        if c not in P.columns:
            rep.err(f"{fname}: label column {c!r} (declared in the manifest) missing")
            continue
        v = P[c].dropna().astype(str).str.strip()
        bad = v[~v.isin(["0", "1", "0.0", "1.0"])]
        if len(bad):
            rep.err(f"{fname}: label column {c!r} must hold 1 (case), 0 (control) or blank; "
                    f"{len(bad)} other value(s)")
    if dp_col and dp_col in P.columns:
        v = P[dp_col].dropna().astype(str).str.strip()
        bad = v[~v.isin(["0", "1", "0.0", "1.0"])]
        if len(bad):
            rep.err(f"{fname}: subgroup device_positive column {dp_col!r} must hold 1, 0 or blank; "
                    f"{len(bad)} other value(s)")
    if "age_band" in P.columns:
        probs = P["age_band"].dropna().map(age_band_problem).dropna()
        if len(probs):
            rep.err(f"{fname}: {len(probs)} age_band value(s) refused: {probs.iloc[0]}")
    if "sex" in P.columns:
        bad = P["sex"].dropna().astype(str).str.strip().str.lower()
        bad = bad[~bad.isin(SEX_VALUES)]
        if len(bad):
            rep.err(f"{fname}: sex must be female / male / other / unknown (or f / m); {len(bad)} other value(s)")
    ids = set(P["participant_id"].dropna().astype(str))
    coverage: dict[str, set] = {}
    features: dict[str, list[str]] = {}
    # ---------------- data files
    for f in inp.files:
        if f.kind == "participants":
            continue
        df = _read(rep, f)
        if df is None:
            continue
        fn = f.path.name
        _common_columns(rep, df, fn)
        if "participant_id" not in df.columns:
            rep.err(f"{fn}: participant_id column missing")
            continue
        pid = df["participant_id"].astype("string")
        orphan = sorted(set(pid.dropna()) - ids)
        if orphan:
            rep.err(f"{fn}: {len(orphan)} participant_id value(s) are not in participants.csv")
        coverage[f.block] = set(pid.dropna()) & ids
        if "day_index" in df.columns:
            d = pd.to_numeric(df["day_index"], errors="coerce")
            if d.isna().any() or (d % 1 != 0).any() or (d < -3660).any() or (d > 3660).any():
                rep.err(f"{fn}: day_index must be an integer day number relative to the participant's own start "
                        "(-3660..3660), never a calendar date")
        if f.kind in ("device", "omics", "survey"):
            feats = _numeric_columns(rep, df, fn, {"participant_id", "day_index"})
            features[f.block] = feats
            if not feats:
                rep.err(f"{fn}: no feature columns")
            if "day_index" not in df.columns and pid.duplicated().any():
                rep.err(f"{fn}: {int(pid.duplicated().sum())} repeated participant_id row(s); a per-day file needs a "
                        "day_index column, a per-participant file one row per participant")
            if f.kind == "omics" and "day_index" in df.columns:
                rep.err(f"{fn}: omics files are one row per participant (no day_index)")
        elif f.kind == "ehr_labs":
            need = {"participant_id", "loinc", "value", "unit"}
            miss = need - set(df.columns)
            if miss:
                rep.err(f"{fn}: missing column(s) {sorted(miss)} (participant_id, loinc, value, unit)")
                continue
            extra = set(df.columns) - need - {"lab_name", "day_index"}
            if extra:
                rep.err(f"{fn}: undeclared column(s) {sorted(extra)} (allowed: participant_id, loinc, value, unit, "
                        "lab_name, day_index)")
            bad = df["loinc"].dropna().astype(str)
            bad = bad[~bad.str.strip().str.fullmatch(LOINC_RE)]
            if len(bad):
                rep.err(f"{fn}: {len(bad)} loinc value(s) are not LOINC codes like 1988-5")
            v = pd.to_numeric(df["value"], errors="coerce")
            if (df["value"].notna() & v.isna()).any():
                rep.err(f"{fn}: {int((df['value'].notna() & v.isna()).sum())} non-numeric lab value(s) (use numbers; "
                        "report below-detection values as the detection limit)")
            for c in ("unit", "lab_name"):
                if c in df.columns and (df[c].dropna().astype(str).str.len() > 80).any():
                    rep.err(f"{fn}: column {c!r} has values longer than 80 characters (free text is not accepted)")
            features[f.block] = sorted(df["loinc"].dropna().astype(str).str.strip().unique())
        elif f.kind == "ehr_diagnoses":
            need = {"participant_id", "icd10cm"}
            if need - set(df.columns):
                rep.err(f"{fn}: needs columns participant_id and icd10cm")
                continue
            extra = set(df.columns) - need - {"day_index"}
            if extra:
                rep.err(f"{fn}: undeclared column(s) {sorted(extra)} (allowed: participant_id, icd10cm, day_index)")
            codes = df["icd10cm"].dropna().astype(str).str.strip().str.upper()
            bad = codes[~codes.str.fullmatch(ICD10CM_RE)]
            if len(bad):
                rep.err(f"{fn}: {len(bad)} icd10cm value(s) are not ICD-10-CM codes like G93.32")
            features[f.block] = sorted(codes.unique())
        elif f.kind == "ehr_medications":
            need = {"participant_id", "rxnorm"}
            if need - set(df.columns):
                rep.err(f"{fn}: needs columns participant_id and rxnorm (RxNorm ingredient CUI)")
                continue
            extra = set(df.columns) - need - {"day_index"}
            if extra:
                rep.err(f"{fn}: undeclared column(s) {sorted(extra)} (allowed: participant_id, rxnorm, day_index)")
            codes = df["rxnorm"].dropna().astype(str).str.strip()
            bad = codes[~codes.str.fullmatch(RXNORM_RE)]
            if len(bad):
                rep.err(f"{fn}: {len(bad)} rxnorm value(s) are not RxNorm CUIs (digits, e.g. 6809)")
            features[f.block] = sorted(codes.unique())
        elif f.kind == "self_report_history":
            need = {"participant_id", "condition"}
            if need - set(df.columns):
                rep.err(f"{fn}: needs columns participant_id and condition (condition name as reported)")
                continue
            extra = set(df.columns) - need - {"icd10cm", "day_index"}
            if extra:
                rep.err(f"{fn}: undeclared column(s) {sorted(extra)} (allowed: participant_id, condition, icd10cm, "
                        "day_index)")
            if (df["condition"].dropna().astype(str).str.len() > 80).any():
                rep.err(f"{fn}: column 'condition' has values longer than 80 characters (free text is not accepted; "
                        "give the condition name only)")
            if "icd10cm" in df.columns:
                codes = df["icd10cm"].dropna().astype(str).str.strip().str.upper()
                codes = codes[codes != ""]
                bad = codes[~codes.str.fullmatch(ICD10CM_RE)]
                if len(bad):
                    rep.err(f"{fn}: {len(bad)} icd10cm value(s) are not ICD-10-CM codes like I73.0")
            from ..harmonize.self_report import map_history
            mapped = df.apply(lambda r: map_history(r.get("condition"), r.get("icd10cm"))["status"] == "mapped",
                              axis=1) if len(df) else pd.Series(dtype=bool)
            if len(df) and not mapped.all():
                rep.warn(f"{fn}: {int((~mapped).sum())} of {len(df)} row(s) do not map to ICD-10-CM "
                         "(configs/harmonize_self_report_icd10.yaml, exact aliases): kept, counted, not modelled")
            features[f.block] = sorted(df["condition"].dropna().astype(str).str.strip().str.lower().unique())
    check_extra_sections(rep, inp, P, label_cols)
    # ---------------- primary analysis + group sizes
    pa = K.analysis_defaults(m.get("primary_analysis"))
    blocks = inp.blocks()
    for key in ("feature_blocks", "combined_blocks"):
        for b in pa.get(key) or []:
            if b not in blocks:
                rep.err(f"manifest.yaml primary_analysis.{key}: block {b!r} has no file (blocks present: "
                        f"{sorted(blocks) or 'none'})")
    label = pa.get("label")
    if label and label not in label_cols:
        rep.err(f"manifest.yaml primary_analysis.label {label!r} is not one of the declared labels {label_cols}")
    for cov in pa.get("covariates") or []:
        if cov not in P.columns:
            rep.err(f"manifest.yaml primary_analysis.covariates: {cov!r} requested but participants.csv has no {cov} "
                    "column")
    mg = m.get("min_group_size") or {}
    n_min_cases = int(min_cases if min_cases is not None else mg.get("cases", DEFAULT_MIN_GROUP))
    n_min_ctrl = int(min_controls if min_controls is not None else mg.get("controls", DEFAULT_MIN_GROUP))
    if min(n_min_cases, n_min_ctrl) < HARD_MIN_GROUP:
        rep.err(f"minimum group size below the engine's floor of {HARD_MIN_GROUP}")
    if min(n_min_cases, n_min_ctrl) < DEFAULT_MIN_GROUP:
        rep.warn(f"minimum group size lowered below {DEFAULT_MIN_GROUP}: estimates will be very imprecise")
    if label in P.columns:
        y = pd.to_numeric(P.set_index("participant_id")[label], errors="coerce")
        with_data = set().union(*[coverage.get(b, set()) for b in pa.get("feature_blocks") or []]) \
            if pa.get("feature_blocks") else set()
        yy = y[y.notna() & y.index.isin(with_data)]
        n1, n0 = int((yy == 1).sum()), int((yy == 0).sum())
        rep.info["primary label"] = f"{label}: {n1} cases, {n0} controls with data in {pa.get('feature_blocks')}"
        if n1 < n_min_cases or n0 < n_min_ctrl:
            rep.err(f"primary label {label!r}: {n1} cases and {n0} controls with data in the primary blocks; at least "
                    f"{n_min_cases} cases and {n_min_ctrl} controls are required")
        if n1 < int(pa["n_splits"]) or n0 < int(pa["n_splits"]):
            rep.err(f"primary label {label!r}: fewer cases or controls than CV folds ({pa['n_splits']})")
        for c in label_cols:
            if c != label and c in P.columns:
                yc = pd.to_numeric(P[c], errors="coerce")
                if int((yc == 1).sum()) < n_min_cases or int((yc == 0).sum()) < n_min_ctrl:
                    rep.warn(f"label {c!r} has fewer than {n_min_cases} cases or {n_min_ctrl} controls: only described, "
                             "not analysed")
        rep.resolved["group_sizes"] = {"cases": n1, "controls": n0}
    # ---------------- summary + what will be created
    ds = K.engine_id(str(m.get("dataset_id", "unknown")))
    rep.info["dataset"] = f"{ds} ({m.get('title', '')}); owner {m.get('owner', '')}"
    rep.info["participants"] = int(P["participant_id"].nunique())
    rep.info["blocks"] = {b: len(features.get(b, [])) for b in blocks}
    rep.info["participants per block"] = {b: len(coverage.get(b, set())) for b in blocks}
    omics = [b for b in blocks if b.startswith("omics_")]
    if omics:
        others = set().union(*[coverage.get(b, set()) for b in blocks if not b.startswith("omics_")] or [set()])
        rep.info["omics overlap with other blocks (same participant_id)"] = {
            b: len(coverage.get(b, set()) & others) for b in omics}
    rep.info["comparator"] = (m.get("comparator") or {}).get("type")
    rep.info["demo / synthetic"] = f"{bool(m.get('demo'))} / {bool(m.get('synthetic'))}"
    if rep.ok:
        rep.will_create = will_create(rep, ds, blocks, pa)
    return rep


def check_extra_sections(rep: Report, inp: K.Inputs, P: pd.DataFrame, label_cols: list[str]) -> None:
    """manifest `omics` (platform declarations) and `subgroup_analysis` against the files present."""
    m = rep.manifest
    blocks = inp.blocks()
    for layer, spec in (m.get("omics") or {}).items():
        if f"omics_{layer}" not in blocks:
            rep.err(f"manifest.yaml omics.{layer}: no omics_{layer}.csv|parquet file")
        if (spec or {}).get("platform") == "nmr_nightingale" and (spec or {}).get("units") != "nightingale_standard":
            rep.warn(f"manifest.yaml omics.{layer}: platform nmr_nightingale without units: nightingale_standard; "
                     "no Nightingale measure will be mapped to LOINC (units are never guessed)")
    sa = m.get("subgroup_analysis")
    if not sa:
        return
    if sa.get("within_label") not in label_cols:
        rep.err(f"manifest.yaml subgroup_analysis.within_label {sa.get('within_label')!r} is not one of the declared "
                f"labels {label_cols}")
    dp = sa.get("device_positive") or {}
    if dp.get("column") and dp["column"] not in P.columns:
        rep.err(f"manifest.yaml subgroup_analysis.device_positive.column {dp['column']!r} is not a participants.csv "
                "column")
    if dp.get("block"):
        f = blocks.get(dp["block"])
        if f is None:
            rep.err(f"manifest.yaml subgroup_analysis.device_positive.block {dp['block']!r} has no file")
        else:
            cols = set(rep.frames.get(f.path.name, pd.DataFrame()).columns)
            has_adapter = bool(((m.get("devices") or {}).get(f.name) or {}).get("adapter"))
            feat = str(dp.get("feature", ""))
            base = feat.split("__")[0]
            if not has_adapter and feat not in cols and base not in cols:
                rep.err(f"manifest.yaml subgroup_analysis.device_positive.feature {feat!r} is not a column of "
                        f"{f.path.name} (per-day files are summarised as <column>__mean / __sd / __n_days)")
    if "device_positive" in sa and not dp:
        rep.err("manifest.yaml subgroup_analysis.device_positive: give `column`, or `block` + `feature` + `threshold`")


def will_create(rep: Report, ds: str, blocks: dict, pa: dict) -> list[str]:
    m = rep.manifest
    classes = _measurement_classes()
    out = [f"data/processed/participants__{ds}.parquet (person layer; ids {ds}:<"
           + ("salted hash" if m.get("id_hashing") == "salted_sha256" else "your participant_id") + ">)",
           f"data/processed/participant_conditions__{ds}.parquet (label rows"
           + (" + ICD-10-CM diagnosis rows)" if "ehr_diagnoses" in blocks else ")")]
    for b, f in blocks.items():
        if f.kind == "device":
            dc = ((m.get("devices") or {}).get(f.name) or {}).get("measurement_class") or m["measurement"]["class"]
            fam = classes.get(dc, {}).get("modality_family")
            t = "participant_wearable_features" if fam == "wearable" else "participant_device_features"
            out.append(f"data/processed/{t}__{ds}.parquet ({b}; class {dc})")
            if "day_index" in rep.frames[f.path.name].columns and fam == "wearable":
                out.append(f"data/processed/participant_wearable_daily__{ds}.parquet ({b}, per day)")
            out.append(f"data/processed/participant_adapter_scores__{ds}.parquet ({b} adapter scores)")
        elif f.kind == "ehr_labs":
            out.append(f"data/processed/participant_labs__{ds}.parquet (LOINC-coded labs)")
        elif f.kind == "omics":
            out.append(f"data/processed/participant_omics_linked__{ds}.parquet ({b}; participant-linked because the "
                       "owner supplies the same participant_id across files; overlap counted above)")
        elif f.kind == "ehr_medications":
            out.append(f"data/processed/participant_medications__{ds}.parquet (RxNorm ingredients)")
        elif f.kind == "survey":
            out.append(f"data/processed/participant_survey__{ds}.parquet ({b}; PRO items, long)")
        elif f.kind == "self_report_history":
            out.append(f"data/processed/participant_conditions__{ds}.parquet (self-reported history mapped to "
                       "ICD-10-CM, mapping_method self_report_map)")
    merged: dict[str, list[str]] = {}          # one line per table (several device files share a table)
    for line in out:
        path, _, why = line.partition(" (")
        merged.setdefault(path, []).append(why.rstrip(")"))
    out = [f"{p} ({'; '.join(w)})" for p, w in merged.items()]
    out += [f"data/processed/{K.DATASETS_TABLE}.parquet (one metadata row for {ds})",
            f"data/raw/{ds}/ (byod_registry_entry.yaml, DATA_AUDIT.md, input_manifest.json; local-only)",
            f"data/processed/person_concepts__{ds}.parquet (harmonized vocabulary: DEMOG, ICD10CM, LOINC, RXNORM, "
            "SURVEY, DEVICE, OMICS; docs/HARMONIZATION_CONTRACT.md)",
            "canonical unions and the Digital Phenotype Vector of this dataset only (never pooled)",
            f"after `evaluate`: results/byod/{K.manifest_id_of(ds)}/ (locked plan, evaluation), "
            f"data/processed/phenotype_signatures__{ds}.parquet, one record in {K.RECORDS_TABLE} "
            f"(tier user_supplied_own_computation; class {m['measurement']['class']}; primary blocks "
            f"{pa.get('feature_blocks')}; specificity {pa.get('specificity')})",
            *([f"after `subgroup`: results/byod/{K.manifest_id_of(ds)}/computable_phenotype.json + SUBGROUP.md, rows "
               f"in {K.PHENOTYPES_TABLE} and {K.STRATA_TABLE} (aggregate only)"] if m.get("subgroup_analysis") else []),
            "after `deploy`: measurement_performance and the affected deployment_opportunities rows for "
            f"{[v['condition_id'] for v in rep.resolved.get('labels', {}).values()]} x "
            f"{rep.resolved.get('bundles_containing_class')}"
            + (" (demo record: only with --demo)" if m.get("demo") else "")]
    return out
