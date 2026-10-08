"""Part B of the launch agent: apply an APPROVED mapping (docs/AGENT_CONTRACT.md §2-§3) to the loaded tables.

    apply_mapping(tables, mapping, out, approve=False) -> report dict
        <out>/.gitignore             `*` (everything under a run folder is local-only)
        <out>/byod/                  a `measure-it byod` folder (docs/BRING_YOUR_OWN_DATA.md), validated on write
        <out>/omop/                  All of Us-shaped OMOP CDM v5.4 parquet (measure_it.agent.omop)
        <out>/apply_report.json      what was written, dropped, refused, and the validator's full report

Privacy transforms (contract §0 rule 4), all deterministic and applied here, never by a model:
  * columns with role `drop` or `unmapped`, and columns absent from the mapping, are never read into an output;
    a column whose NAME the BYOD validator classes as identifying or geographic is dropped even when the mapping gives
    it a data role (`dropped_by_privacy_guard`);
  * participant ids stay study codes; ids that look identifying (MRN-like numbers, names, e-mails, ...: the patterns of
    `measure_it.byod.checks.id_problems`), an id column whose name looks identifying, or `id_hashing: salted_sha256` in
    the mapping are replaced by `h<16 hex of sha256(salt ":" id)>` with the salt from MEASURE_IT_BYOD_SALT (the same
    rule as `byod ingest`); without a salt such ids are REFUSED (nothing is written);
  * dates become `day_index` per participant (the participant's first observed date over every table = 0); the date
    itself is never written; a numeric "date" column is taken as an already-relative day number;
  * exact ages become 10-year `age_band` (90+ top-coded; supplied bands are widened to their decade);
  * text values written to an output (self-reported condition names, demographic categories) that match the
    validator's identifying-value patterns, or exceed 80 characters, are dropped and counted.

Data use (manifest `data_use`): a public dataset (`public=True`) gets `deidentified: true, use_permitted: true,
statement: "public dataset: <source>"`. For private data the engine never asserts permission: the mapping must carry a
user-confirmed `data_use` block (deidentified + use_permitted true, a statement), otherwise placeholders that FAIL
validation are written and `data_use_confirmation_required` is true in the report.
"""
from __future__ import annotations

import datetime as dt
import getpass
import hashlib
import json
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from ..byod import checks as C
from ..byod import common as K

WRITING_ROLES = {"label", "age", "age_band", "sex", "demographic", "self_report_history", "diagnosis_code", "lab",
                 "medication_code", "medication_flag", "survey", "device", "omics"}
KNOWN_ROLES = WRITING_ROLES | {"participant_id", "date", "drop", "unmapped"}
YES = {"1", "1.0", "yes", "y", "true", "t", "case", "positive", "pos"}
NO = {"0", "0.0", "no", "n", "false", "f", "control", "negative", "neg"}
SEX_DEFAULT = {"f": "female", "female": "female", "woman": "female", "w": "female", "m": "male", "male": "male",
               "man": "male", "other": "other", "intersex": "other", "nonbinary": "other", "non-binary": "other",
               "unknown": "unknown", "u": "unknown", "": "unknown"}
NAME_SAFE = re.compile(r"[^a-z0-9_]+")
FEATURE_SAFE = re.compile(r"[^A-Za-z0-9_.\-]+")
MAX_TEXT = 80


class ApplyRefused(PermissionError):
    """apply_mapping refuses to write anything; `.reasons` lists every reason."""

    def __init__(self, reasons: list[str]):
        self.reasons = list(reasons)
        super().__init__("; ".join(self.reasons))


# ------------------------------------------------------------------------------------------------------------------
# the harmonized bundle (shared by the BYOD writer and the OMOP writer)
# ------------------------------------------------------------------------------------------------------------------

@dataclass
class Harmonized:
    dataset_id: str
    participants: pd.DataFrame                                   # participant_id, <label cols>, age_band, sex
    labels: list[dict] = field(default_factory=list)            # column, condition, definition, label_basis
    demographics: pd.DataFrame | None = None                     # participant_id, variable, value (OMOP only)
    devices: dict[str, pd.DataFrame] = field(default_factory=dict)
    device_meta: dict[str, dict] = field(default_factory=dict)  # name -> measurement_class, description
    labs: pd.DataFrame | None = None                             # participant_id, loinc, value, unit, lab_name,
                                                                 # day_index, platform
    diagnoses: pd.DataFrame | None = None                        # participant_id, icd10cm, day_index
    medications: pd.DataFrame | None = None                      # participant_id, rxnorm, day_index
    history: pd.DataFrame | None = None                          # participant_id, condition, icd10cm, day_index
    surveys: dict[str, pd.DataFrame] = field(default_factory=dict)
    omics: dict[str, pd.DataFrame] = field(default_factory=dict)
    omics_meta: dict[str, dict] = field(default_factory=dict)   # layer -> platform, units, description
    warnings: list[str] = field(default_factory=list)


LAB_COLS = ["participant_id", "loinc", "value", "unit", "lab_name", "day_index", "platform"]
DX_COLS = ["participant_id", "icd10cm", "day_index"]
RX_COLS = ["participant_id", "rxnorm", "day_index"]
HX_COLS = ["participant_id", "condition", "icd10cm", "day_index"]


def _frame(rows: list[dict], cols: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)


# ------------------------------------------------------------------------------------------------------------------
# small deterministic transforms
# ------------------------------------------------------------------------------------------------------------------

def _key(v) -> str:
    s = str(v).strip()
    try:
        f = float(s)
        if np.isfinite(f) and f == int(f):
            return str(int(f))
    except ValueError:
        pass
    return s.lower()


def apply_value_map(s: pd.Series, value_map: dict | None) -> pd.Series:
    """Map values through a value_map (keys compared as trimmed, case-insensitive text; '1' == 1 == 1.0). Values not
    in the map become NaN."""
    if not value_map:
        return s
    vm = {_key(k): v for k, v in value_map.items()}
    return s.map(lambda v: vm.get(_key(v), np.nan) if pd.notna(v) else np.nan)


def to_binary(s: pd.Series, value_map: dict | None = None) -> pd.Series:
    """-> 1.0 / 0.0 / NaN. With a value_map its targets must be 1/0 (other targets -> NaN); without one, yes/no-like
    values are recognised."""
    if value_map:
        m = apply_value_map(s, value_map)
        return m.map(lambda v: np.nan if pd.isna(v) else 1.0 if _key(v) in YES else 0.0 if _key(v) in NO else np.nan)
    return s.map(lambda v: 1.0 if pd.notna(v) and _key(v) in YES else 0.0 if pd.notna(v) and _key(v) in NO else np.nan)


def normalize_sex(s: pd.Series, value_map: dict | None = None) -> pd.Series:
    m = apply_value_map(s, value_map) if value_map else s
    return m.map(lambda v: SEX_DEFAULT.get(str(v).strip().lower(), "unknown") if pd.notna(v) else "unknown")


def age_to_band(s: pd.Series) -> pd.Series:
    """Exact ages -> 10-year bands (90+); values already like '40-49' / '90+' are widened to their decade."""
    from ..harmonize.schema import age_band_from_supplied, age_band_from_years
    num = pd.to_numeric(s, errors="coerce")
    out = age_band_from_years(num)
    txt = s.where(num.isna())
    if txt.notna().any():
        sup = age_band_from_supplied(txt.astype(object))
        out = pd.Series(np.where(num.notna(), out, sup), index=s.index, dtype=object)
    out = out.map(lambda v: v if v is None or (isinstance(v, str) and C.age_band_problem(v) is None) else None)
    return out


def normalize_icd10(v) -> str | None:
    s = str(v or "").strip().upper().replace(" ", "")
    if not s or s in ("NAN", "NONE"):
        return None
    if "." not in s and len(s) > 3:
        s = s[:3] + "." + s[3:]
    return s if C.ICD10CM_RE.fullmatch(s) else None


def normalize_rxnorm(v) -> str | None:
    s = _key(v)
    return s if C.RXNORM_RE.fullmatch(s) else None


def normalize_loinc(v) -> str | None:
    s = str(v or "").strip()
    return s if C.LOINC_RE.fullmatch(s) else None


def safe_name(text: str, fallback: str = "x") -> str:
    s = NAME_SAFE.sub("_", str(text).lower()).strip("_")
    s = re.sub(r"_+", "_", s)
    if not s or not s[0].isalpha():
        s = f"{fallback}_{s}" if s else fallback
    return s[:41]


def safe_feature(text: str) -> str:
    s = FEATURE_SAFE.sub("_", str(text)).strip("_")
    if not s or not s[0].isalpha():
        s = "f_" + s
    return s[:81]


def _text_ok(v: str) -> bool:
    """A text value that may be written: short, and no identifying-value pattern (checks.VALUE_PATTERNS)."""
    s = str(v)
    if len(s) > MAX_TEXT:
        return False
    return not any(rx.search(s) for name, rx in C.VALUE_PATTERNS.items() if name != "ZIP code") \
        and not C.VALUE_PATTERNS["ZIP code"].search(s)


def _privacy_guard(col: str) -> str | None:
    why = C.column_problem(col)
    if why and ("geographic" in why or "identifies" in why):
        return why
    return None


def hash_id(salt: str, native: str) -> str:
    return "h" + hashlib.sha256(f"{salt}:{native}".encode()).hexdigest()[:16]


# ------------------------------------------------------------------------------------------------------------------
# reading the mapping
# ------------------------------------------------------------------------------------------------------------------

def is_approved(mapping: dict) -> bool:
    return bool(mapping.get("approved_by")) and bool(mapping.get("approved_at"))


def _table_specs(mapping: dict) -> dict:
    t = mapping.get("tables") or {}
    if not isinstance(t, dict) or not t:
        raise ApplyRefused(["mapping has no `tables` section (docs/AGENT_CONTRACT.md §2)"])
    return t


def _id_column(name: str, spec: dict) -> str | None:
    if spec.get("participant_id"):
        return spec["participant_id"]
    for c, cs in (spec.get("columns") or {}).items():
        if (cs or {}).get("role") == "participant_id":
            return c
    return None


def _day_column(spec: dict) -> str | None:
    if spec.get("day_index_from"):
        return spec["day_index_from"]
    for c, cs in (spec.get("columns") or {}).items():
        if (cs or {}).get("role") == "date":
            return c
    return None


def _long_lab_spec(spec: dict) -> dict | None:
    """Long lab tables: table-level `long_labs: {code, value, unit[, name]}` or column specs
    `{role: lab, part: code|value|unit|name}` (contract §2, long tables)."""
    ll = dict(spec.get("long_labs") or {})
    for c, cs in (spec.get("columns") or {}).items():
        cs = cs or {}
        if cs.get("role") == "lab" and cs.get("part") in ("code", "value", "unit", "name"):
            ll.setdefault(cs["part"], c)
            if cs["part"] == "code" and cs.get("code_map"):
                ll.setdefault("code_map", cs["code_map"])
    return ll if ll.get("code") and ll.get("value") else None


# ------------------------------------------------------------------------------------------------------------------
# harmonize
# ------------------------------------------------------------------------------------------------------------------

def harmonize(tables: dict[str, pd.DataFrame], mapping: dict) -> tuple[Harmonized, dict]:
    """Tables + mapping -> Harmonized bundle and an audit dict (dropped, unmapped, counts). Raises ApplyRefused."""
    specs = _table_specs(mapping)
    audit = {"dropped": [], "dropped_by_privacy_guard": [], "unmapped": [], "not_in_mapping": [],
             "values_dropped": Counter(), "tables_skipped": [], "participant_id_transform": "study codes kept"}
    warnings: list[str] = []
    refusals: list[str] = []
    ds = safe_name(mapping.get("dataset_id") or "agent_dataset", "ds")
    if len(ds) < 3:
        ds = (ds + "_ds")[:41]
    ds = re.sub(r"_{2,}", "_", ds)

    # ---- 1. participant ids (decide hashing over ALL tables first)
    id_cols: dict[str, str] = {}
    raw_ids: list[pd.Series] = []
    force_hash = mapping.get("id_hashing") == "salted_sha256"
    for tname, spec in specs.items():
        spec = spec or {}
        if tname not in tables:
            warnings.append(f"table {tname!r} is in the mapping but was not loaded; skipped")
            audit["tables_skipped"].append(tname)
            continue
        idc = _id_column(tname, spec)
        df = tables[tname]
        if not idc or idc not in df.columns:
            warnings.append(f"table {tname!r}: no participant_id column in the mapping (or it is absent); skipped")
            audit["tables_skipped"].append(tname)
            continue
        id_cols[tname] = idc
        raw_ids.append(df[idc].dropna().astype(str).str.strip())
        if _privacy_guard(idc) and idc.lower() not in ("participant_id", "record_id", "subject_id", "study_id"):
            force_hash = True
    if not id_cols:
        raise ApplyRefused(["no table has a mapped participant_id column"])
    all_ids = pd.Series(pd.concat(raw_ids).unique()) if raw_ids else pd.Series([], dtype=str)
    all_ids = all_ids[all_ids != ""]
    probs = C.id_problems(all_ids, "participant_id")
    if probs or force_hash:
        salt = os.environ.get(K.SALT_ENV, "")
        if len(salt) < 16:
            why = (["participant ids look identifying: " + "; ".join(p.split(": ", 1)[-1] for p in probs)]
                   if probs else ["participant ids must be hashed (identifying id column name or id_hashing)"])
            raise ApplyRefused(why + [f"set {K.SALT_ENV} (>= 16 characters, kept by the data owner, never written) "
                                      "to replace them with salted hashes, or map a study-code column instead"])
        id_map = {v: hash_id(salt, v) for v in all_ids}
        audit["participant_id_transform"] = "salted_sha256 (salt from MEASURE_IT_BYOD_SALT, not stored)"
    else:
        id_map = {v: v for v in all_ids}

    # ---- 2. dates -> day_index (first observed date per participant over every table = 0)
    parsed: dict[str, pd.Series] = {}
    relative: dict[str, pd.Series] = {}
    for tname, idc in id_cols.items():
        dcol = _day_column(specs[tname] or {})
        if not dcol:
            continue
        df = tables[tname]
        if dcol not in df.columns:
            warnings.append(f"table {tname!r}: date column {dcol!r} not found; no day_index")
            continue
        s = df[dcol]
        num = pd.to_numeric(s, errors="coerce")
        if s.notna().sum() and num.notna().sum() >= 0.9 * s.notna().sum():
            relative[tname] = num.round()
        else:
            parsed[tname] = pd.to_datetime(s.astype(str).where(s.notna()), errors="coerce", format="mixed")
            bad = int((s.notna() & parsed[tname].isna()).sum())
            if bad:
                warnings.append(f"table {tname!r}: {bad} value(s) of the date column could not be read as dates "
                                "(rows kept without day_index)")
    first: dict[str, pd.Timestamp] = {}
    for tname, d in parsed.items():
        pid = tables[tname][id_cols[tname]].astype(str).str.strip().map(id_map)
        g = pd.DataFrame({"p": pid, "d": d}).dropna().groupby("p")["d"].min()
        for p, v in g.items():
            first[p] = min(first.get(p, v), v)

    def day_index(tname: str, pid: pd.Series) -> pd.Series | None:
        """Aligned to pid's index (the kept rows of the table)."""
        if tname in relative:
            return relative[tname].loc[pid.index]
        if tname in parsed:
            return (parsed[tname].loc[pid.index] - pid.map(first)).dt.days.astype("float")
        return None

    # ---- 3. column roles
    labels: dict[str, dict] = {}
    label_vals: dict[str, list[pd.DataFrame]] = {}
    ages, sexes, demo = [], [], []
    devices: dict[str, list[pd.DataFrame]] = {}
    device_meta: dict[str, dict] = {}
    surveys: dict[str, list[pd.DataFrame]] = {}
    omics: dict[str, list[pd.DataFrame]] = {}
    omics_meta: dict[str, dict] = {}
    lab_rows, dx_rows, rx_rows, hx_rows = [], [], [], []
    nightingale = _nightingale_by_loinc()

    for tname, idc in id_cols.items():
        spec = specs[tname] or {}
        cols = spec.get("columns") or {}
        df = tables[tname]
        pid = df[idc].astype(str).str.strip().map(id_map)
        keep = pid.notna() & (pid != "")
        df, pid = df[keep], pid[keep]
        dix = day_index(tname, pid)
        repeated = bool(pid.duplicated().any())
        dcol = _day_column(spec)
        for c in df.columns:
            if c not in cols and c != idc and c != dcol:
                audit["not_in_mapping"].append(f"{tname}.{c}")
        ll = _long_lab_spec(spec)
        long_parts = set(ll.values()) if ll else set()
        if ll:
            code_map = ll.get("code_map") or {}
            codes = df[ll["code"]].map(lambda v: code_map.get(str(v).strip(), v) if pd.notna(v) else v)
            loinc = codes.map(normalize_loinc)
            val = pd.to_numeric(df[ll["value"]], errors="coerce")
            unit = df[ll["unit"]].astype("string") if ll.get("unit") in df.columns else pd.Series(pd.NA, index=df.index)
            name = df[ll["name"]].astype("string") if ll.get("name") in df.columns else pd.Series(pd.NA, index=df.index)
            ok = loinc.notna() & val.notna()
            audit["values_dropped"][f"{tname}: long lab rows without a LOINC code or numeric value"] += int((~ok).sum())
            for i in np.flatnonzero(ok.to_numpy()):
                n = name.iloc[i]
                lab_rows.append({"participant_id": pid.iloc[i], "loinc": loinc.iloc[i], "value": float(val.iloc[i]),
                                 "unit": None if pd.isna(unit.iloc[i]) else str(unit.iloc[i]).strip()[:MAX_TEXT],
                                 "lab_name": str(n)[:MAX_TEXT] if pd.notna(n) and _text_ok(n) else None,
                                 "day_index": None if dix is None or pd.isna(dix.iloc[i]) else int(dix.iloc[i]),
                                 "platform": None})
        for c, cs in cols.items():
            cs = dict(cs or {})
            role = cs.get("role", "unmapped")
            if c in long_parts or c == idc or role in ("participant_id", "date"):
                continue
            if role not in KNOWN_ROLES:
                warnings.append(f"{tname}.{c}: unknown role {role!r}; not written")
                audit["unmapped"].append({"column": f"{tname}.{c}", "reason": f"unknown role {role!r}"})
                continue
            if role == "drop":
                audit["dropped"].append({"column": f"{tname}.{c}", "reason": cs.get("reason", "drop")})
                continue
            if role == "unmapped":
                audit["unmapped"].append({"column": f"{tname}.{c}", "reason": cs.get("reason", "")})
                continue
            if c not in df.columns:
                warnings.append(f"{tname}.{c}: in the mapping but not in the table; skipped")
                continue
            guard = _privacy_guard(c)
            if guard and role in WRITING_ROLES - {"age", "age_band"}:
                audit["dropped_by_privacy_guard"].append({"column": f"{tname}.{c}", "role": role,
                                                          "reason": f"column name {guard}"})
                continue
            s = df[c]
            base = pd.DataFrame({"participant_id": pid.to_numpy()}, index=df.index)
            if role == "label":
                for lab in (cs.get("labels") or [cs]):
                    cond = lab.get("condition") or mapping.get("condition_hint")
                    if not cond:
                        warnings.append(f"{tname}.{c}: label without a condition; not written")
                        continue
                    col = safe_name(lab.get("column") or cond, "label")
                    vm = lab.get("value_map")
                    y = to_binary(s, vm)
                    labels.setdefault(col, {
                        "column": col, "condition": str(cond),
                        "label_basis": lab.get("label_basis") or cs.get("label_basis"),
                        "definition": lab.get("definition") or (
                            f"label {cond} from source column {c!r}"
                            + (f" (value_map {json.dumps({str(k): v for k, v in vm.items()})})" if vm else
                               " (yes/no-like values)"))})
                    label_vals.setdefault(col, []).append(base.assign(v=y.to_numpy()))
            elif role in ("age", "age_band"):
                ages.append(base.assign(v=age_to_band(s).to_numpy(), n=pd.to_numeric(s, errors="coerce").to_numpy()))
            elif role == "sex":
                sexes.append(base.assign(v=normalize_sex(s, cs.get("value_map")).to_numpy()))
            elif role == "demographic":
                v = apply_value_map(s, cs.get("value_map")) if cs.get("value_map") else s
                v = v.where(v.isna(), v.astype(str).str.strip())
                bad = v.notna() & ~v.map(lambda x: _text_ok(x) if pd.notna(x) else True)
                audit["values_dropped"][f"{tname}.{c}: demographic values with identifying patterns or > 80 chars"] \
                    += int(bad.sum())
                demo.append(base.assign(variable=safe_name(cs.get("variable") or c, "demo"),
                                        value=v.where(~bad).to_numpy()))
            elif role == "self_report_history":
                vm = cs.get("value_map")
                code = normalize_icd10(cs.get("icd10cm")) if cs.get("icd10cm") else None
                if vm or code:          # wide yes/no column for one condition
                    y = to_binary(s, vm)
                    name = str(cs.get("condition") or cs.get("label") or c)[:MAX_TEXT]
                    for i in np.flatnonzero((y == 1).to_numpy()):
                        hx_rows.append({"participant_id": pid.iloc[i], "condition": name, "icd10cm": code,
                                        "day_index": None if dix is None or pd.isna(dix.iloc[i]) else int(dix.iloc[i])})
                else:                   # long column of condition names
                    for i, v in enumerate(s.to_numpy()):
                        if pd.isna(v) or not str(v).strip():
                            continue
                        if not _text_ok(str(v).strip()):
                            audit["values_dropped"][f"{tname}.{c}: condition names > 80 chars or identifying"] += 1
                            continue
                        hx_rows.append({"participant_id": pid.iloc[i], "condition": str(v).strip(), "icd10cm": None,
                                        "day_index": None if dix is None or pd.isna(dix.iloc[i]) else int(dix.iloc[i])})
            elif role == "diagnosis_code":
                codes = s.map(normalize_icd10)
                audit["values_dropped"][f"{tname}.{c}: values that are not ICD-10-CM codes"] += \
                    int((s.notna() & codes.isna()).sum())
                for i in np.flatnonzero(codes.notna().to_numpy()):
                    dx_rows.append({"participant_id": pid.iloc[i], "icd10cm": codes.iloc[i],
                                    "day_index": None if dix is None or pd.isna(dix.iloc[i]) else int(dix.iloc[i])})
            elif role == "medication_code":
                codes = s.map(normalize_rxnorm)
                audit["values_dropped"][f"{tname}.{c}: values that are not RxNorm CUIs (digits)"] += \
                    int((s.notna() & codes.isna()).sum())
                for i in np.flatnonzero(codes.notna().to_numpy()):
                    rx_rows.append({"participant_id": pid.iloc[i], "rxnorm": codes.iloc[i],
                                    "day_index": None if dix is None or pd.isna(dix.iloc[i]) else int(dix.iloc[i])})
            elif role == "medication_flag":
                code = normalize_rxnorm(cs.get("rxnorm"))
                if not code:
                    audit["unmapped"].append({"column": f"{tname}.{c}", "reason": "medication_flag without a valid "
                                                                                 "rxnorm ingredient CUI"})
                    continue
                y = to_binary(s, cs.get("value_map"))
                for i in np.flatnonzero((y == 1).to_numpy()):
                    rx_rows.append({"participant_id": pid.iloc[i], "rxnorm": code,
                                    "day_index": None if dix is None or pd.isna(dix.iloc[i]) else int(dix.iloc[i])})
            elif role == "lab":
                loinc = normalize_loinc(cs.get("loinc"))
                if not loinc:
                    audit["unmapped"].append({"column": f"{tname}.{c}", "reason": f"lab without a valid LOINC code "
                                                                                 f"({cs.get('loinc')!r})"})
                    continue
                val = pd.to_numeric(apply_value_map(s, cs.get("value_map")) if cs.get("value_map") else s,
                                    errors="coerce")
                audit["values_dropped"][f"{tname}.{c}: non-numeric lab values"] += int((s.notna() & val.isna()).sum())
                platform = cs.get("platform")
                unit = cs.get("unit")
                ng = nightingale.get(loinc) if platform and ("nmr" in str(platform) or "nightingale" in str(platform)) \
                    else None
                if ng and unit and _unit_eq(unit, ng["source_unit"]):
                    # NMR measure in Nightingale's own units: the BYOD omics route keeps the platform flag
                    layer = "nightingale"
                    omics.setdefault(layer, []).append(
                        base.assign(**{ng["measure"]: val.to_numpy()}).groupby("participant_id")[[ng["measure"]]].mean())
                    omics_meta.setdefault(layer, {"platform": "nmr_nightingale", "units": "nightingale_standard",
                                                  "description": "NMR measures routed from lab-role columns"})
                    continue
                if platform:
                    warnings.append(f"{tname}.{c}: platform {platform!r} cannot be carried by ehr_labs.csv (the BYOD "
                                    "lab file has no platform column); it is kept in the OMOP measurement_source_value")
                for i in np.flatnonzero(val.notna().to_numpy()):
                    lab_rows.append({"participant_id": pid.iloc[i], "loinc": loinc, "value": float(val.iloc[i]),
                                     "unit": unit, "lab_name": cs.get("label") or c,
                                     "day_index": None if dix is None or pd.isna(dix.iloc[i]) else int(dix.iloc[i]),
                                     "platform": platform})
            elif role in ("survey", "device", "omics"):
                if role == "survey":
                    key = safe_name(cs.get("instrument") or "survey", "survey")
                    feat = safe_feature(cs.get("item") or c)
                elif role == "device":
                    key = safe_name(cs.get("block") or "device", "device")
                    feat = safe_feature(cs.get("feature") or c)
                    meta = device_meta.setdefault(key, {})
                    if cs.get("measurement_class"):
                        meta.setdefault("measurement_class", cs["measurement_class"])
                    if cs.get("description"):
                        meta.setdefault("description", cs["description"])
                else:
                    key = safe_name(cs.get("layer") or "omics", "omics")
                    feat = safe_feature(cs.get("feature") or c)
                    if cs.get("platform"):
                        omics_meta.setdefault(key, {})["platform"] = safe_name(cs["platform"], "p")
                    if cs.get("units"):
                        omics_meta.setdefault(key, {})["units"] = str(cs["units"])
                v = apply_value_map(s, cs.get("value_map")) if cs.get("value_map") else s
                num = pd.to_numeric(v, errors="coerce")
                audit["values_dropped"][f"{tname}.{c}: non-numeric {role} values"] += int((s.notna() & num.isna()).sum())
                frame = base.assign(**{feat: num.to_numpy()})
                if dix is not None and repeated and role != "omics":
                    frame["day_index"] = dix.to_numpy()
                    frame = frame.dropna(subset=["day_index"]).set_index(["participant_id", "day_index"])
                else:
                    if repeated:
                        warnings.append(f"{tname}.{c}: several rows per participant and no day_index; "
                                        f"{role} values averaged per participant")
                    frame = frame.groupby("participant_id")[[feat]].mean()
                if frame.index.duplicated().any():
                    frame = frame.groupby(level=list(range(frame.index.nlevels))).mean()
                {"survey": surveys, "device": devices, "omics": omics}[role].setdefault(key, []).append(frame)
    if refusals:
        raise ApplyRefused(refusals)

    # ---- 4. assemble
    ids = sorted(set(id_map.values()))
    P = pd.DataFrame({"participant_id": ids})
    for col, frames in label_vals.items():
        lv = pd.concat(frames)
        g = lv.dropna(subset=["v"]).groupby("participant_id")["v"]
        conflict = int((g.nunique() > 1).sum())
        if conflict:
            warnings.append(f"label {col!r}: {conflict} participant(s) with conflicting values; set blank")
        val = g.first().where(g.nunique() <= 1)
        P[col] = P["participant_id"].map(val)
    if ages:
        a = pd.concat(ages)
        a = a.sort_values("n", na_position="last").dropna(subset=["v"]).groupby("participant_id")["v"].first()
        P["age_band"] = P["participant_id"].map(a)
    if sexes:
        sx = pd.concat(sexes)
        sx = sx[sx["v"] != "unknown"].groupby("participant_id")["v"].first()
        P["sex"] = P["participant_id"].map(sx).fillna("unknown")
    for lab in labels.values():
        if not lab.get("label_basis"):
            lab["label_basis"] = "proxy"
            warnings.append(f"label {lab['column']!r}: label_basis not given in the mapping; recorded as `proxy` "
                            "(the conservative basis) - edit the mapping if a case definition applies")
        if len(lab["definition"]) < 10:
            lab["definition"] = f"label {lab['condition']}: {lab['definition']}"

    def _merge(frames: list[pd.DataFrame]) -> pd.DataFrame:
        if len({f.index.nlevels for f in frames}) > 1:     # per-day and per-participant columns in one block
            warnings.append("a block mixes per-day and per-participant columns; per-day values averaged per "
                            "participant")
            frames = [f.groupby(level=0).mean() if f.index.nlevels > 1 else f for f in frames]
        out = pd.concat(frames, axis=1)
        out = out.T.groupby(level=0).first().T if out.columns.duplicated().any() else out
        return out.reset_index()

    h = Harmonized(
        dataset_id=ds, participants=P, labels=list(labels.values()),
        demographics=pd.concat(demo, ignore_index=True).dropna(subset=["value"]) if demo else None,
        devices={k: _merge(v) for k, v in devices.items()}, device_meta=device_meta,
        labs=_frame(lab_rows, LAB_COLS), diagnoses=_frame(dx_rows, DX_COLS).drop_duplicates(),
        medications=_frame(rx_rows, RX_COLS).drop_duplicates(), history=_frame(hx_rows, HX_COLS).drop_duplicates(),
        surveys={k: _merge(v) for k, v in surveys.items()},
        omics={k: _merge(v) for k, v in omics.items()}, omics_meta=omics_meta, warnings=warnings)
    audit["values_dropped"] = {k: v for k, v in audit["values_dropped"].items() if v}
    return h, audit


def _unit_eq(a: str, b: str) -> bool:
    n = lambda u: str(u).strip().lower().replace("µ", "u").replace(" ", "")  # noqa: E731
    return n(a) == n(b)


def _nightingale_by_loinc() -> dict:
    try:
        from ..harmonize import schema as S
        return {r["loinc"]: r for r in S.config("harmonize_nightingale_loinc")["mappings"]}
    except Exception:  # noqa: BLE001 - config missing: no NMR routing
        return {}


# ------------------------------------------------------------------------------------------------------------------
# manifest
# ------------------------------------------------------------------------------------------------------------------

def _classes() -> dict:
    from ..config import load_config
    return {m["id"]: m for m in load_config("measurements")["measurement_classes"]}


def _bundles() -> dict:
    from ..config import load_config
    return load_config("relevance").get("measurement_bundles", {}) or {}


def resolve_measurement_class(hint: str | None, device_classes: list[str]) -> tuple[str | None, str | None, str]:
    """(class, bundle, how) from the mapping's measurement_hint and the device blocks' classes."""
    classes, bundles = _classes(), _bundles()
    if hint and hint in classes:
        return hint, None, "measurement_hint is a class"
    if hint:
        bundle = hint if hint in bundles else None
        members = list((bundles.get(hint) or {}).get("members", []))
        if not bundle:
            try:
                from ..measurements.resolve import resolve_measurement
                r = resolve_measurement(hint, single=True)
                if r.get("status") == "matched" and r.get("measurement_id"):
                    if r.get("kind") == "class" and r["measurement_id"] in classes:
                        return r["measurement_id"], None, f"measurement_hint resolved: {r.get('match_reason')}"
                    if r["measurement_id"] in bundles:
                        bundle, members = r["measurement_id"], list(r.get("member_ids") or [])
            except Exception:  # noqa: BLE001 - resolver needs configs only; fall through to device classes
                pass
        if bundle:
            pick = next((c for c in device_classes if c in members), None) or (members[0] if members else None)
            return pick, bundle, f"measurement_hint is bundle {bundle}; class {pick} (a member present in the devices)"
    if device_classes:
        top = Counter(device_classes).most_common(1)[0][0]
        return top, None, "most common device measurement_class (no usable measurement_hint)"
    return None, None, "no measurement_hint and no device measurement_class"


def build_manifest(h: Harmonized, mapping: dict, *, public: bool) -> tuple[dict, dict]:
    """(manifest, info) for the BYOD folder. info: data_use_source, data_use_confirmation_required, caveats."""
    caveats: list[str] = []
    classes = _classes()
    devices = {}
    for name in h.devices:
        meta = h.device_meta.get(name, {})
        mc = meta.get("measurement_class")
        if mc and mc not in classes:
            mc2, _, _ = resolve_measurement_class(mc, [])
            if mc2 in classes:
                mc = mc2
            else:
                h.warnings.append(f"device {name}: measurement_class {mc!r} is not in configs/measurements.yaml; "
                                  "omitted")
                mc = None
        d = {"description": meta.get("description") or f"device block {name} (from the approved mapping)"}
        if mc:
            d["measurement_class"] = mc
        devices[name] = d
    dev_classes = [d["measurement_class"] for d in devices.values() if d.get("measurement_class")]
    cls, bundle, how = resolve_measurement_class(mapping.get("measurement_hint"), dev_classes)
    source = str(mapping.get("source") or mapping.get("dataset_id") or h.dataset_id)
    du_in = mapping.get("data_use") or {}
    if public:
        data_use = {"deidentified": True, "use_permitted": True, "statement": f"public dataset: {source}"}
        du_source, confirm = "public", False
    elif du_in.get("deidentified") is True and du_in.get("use_permitted") is True \
            and len(str(du_in.get("statement") or "")) >= 10:
        data_use = {k: du_in[k] for k in ("deidentified", "use_permitted", "statement", "irb_or_dua_reference")
                    if k in du_in}
        du_source, confirm = "mapping (confirmed by the user)", False
    else:
        data_use = {"deidentified": False, "use_permitted": False,
                    "statement": "PLACEHOLDER: the data owner must confirm that these data are de-identified and that "
                                 "an IRB protocol / data-use agreement permits this analysis on this machine "
                                 "(add a data_use block to mapping.yaml and re-apply)"}
        du_source, confirm = "placeholder (not confirmed)", True
    comp = mapping.get("comparator")
    if isinstance(comp, str):
        comp = {"type": comp}
    if comp and comp.get("type") in K.COMPARATOR_TYPES:
        comparator = {"type": comp["type"], "description": comp.get("description") or f"{comp['type']} comparator "
                                                                                       "(from the mapping)"}
    else:
        comparator = {"type": "healthy", "description": "DEFAULT: controls assumed to be healthy (the mapping did not "
                                                        "state a comparator); edit if they are look-alike patients"}
        caveats.append("comparator defaulted to `healthy`: every performance number against healthy controls is "
                       "optimistic; state the comparator in the mapping if the controls are look-alike patients")
    labels = [{"column": lab["column"], "condition": lab["condition"], "definition": lab["definition"],
               "label_basis": lab["label_basis"] if lab["label_basis"] in K.LABEL_BASES else "proxy"}
              for lab in h.labels]
    label_cols = [lab["column"] for lab in labels]
    hint = safe_name(mapping.get("condition_hint") or "", "c") if mapping.get("condition_hint") else None
    primary_label = next((lab["column"] for lab in labels if lab["condition"] == mapping.get("condition_hint")
                          or lab["column"] == hint), label_cols[0] if label_cols else None)
    blocks = _blocks(h)
    dev_blocks = [f"device_{n}" for n, d in devices.items() if d.get("measurement_class") == cls] \
        or [f"device_{n}" for n in devices]
    fallback = [b for b in blocks if b.startswith(("ehr_labs", "omics_", "ehr_diagnoses"))][:1]
    pa = {"label": primary_label, "feature_blocks": dev_blocks or fallback or blocks[:1],
          "specificity": 0.90, "covariates": [c for c in ("age_band", "sex") if c in h.participants.columns],
          "n_splits": 5, "n_repeats": 10, "n_permutations": 200, "n_bootstrap": 1000}
    if mapping.get("primary_analysis"):
        pa.update(mapping["primary_analysis"])
    syn = bool(mapping.get("synthetic"))
    title = str(mapping.get("title") or f"{h.dataset_id}: harmonized by the measure-it launch agent")
    manifest = {
        "schema_version": K.SCHEMA_VERSION, "dataset_id": h.dataset_id, "title": title,
        "owner": str(mapping.get("owner") or (source if public else "data owner (not stated in the mapping)")),
        "description": str(mapping.get("description") or "Written by measure_it.agent.apply from an approved mapping "
                           f"(approved_by {mapping.get('approved_by')}, approved_at {mapping.get('approved_at')})."),
        "licence": str(mapping.get("licence") or (f"public dataset: {source} (see the source's licence)" if public
                                                   else "PLACEHOLDER: data-use statement to be supplied by the owner")),
        "synthetic": syn, "demo": syn or bool(mapping.get("demo")),
        "data_use": data_use, "labels": labels, "comparator": comparator,
        "measurement": {"class": cls} if cls else {},
        "id_hashing": "none",
        "primary_analysis": pa,
    }
    if bundle and cls and cls in (_bundles().get(bundle) or {}).get("members", []):
        manifest["measurement"]["bundle"] = bundle
    if mapping.get("measurement_hint"):
        manifest["measurement"]["device"] = str(mapping["measurement_hint"])
    if devices:
        manifest["devices"] = devices
    om = {}
    for layer in h.omics:
        meta = h.omics_meta.get(layer, {})
        e = {k: meta[k] for k in ("platform", "units", "description") if meta.get(k)}
        if e:
            om[layer] = e
    if om:
        manifest["omics"] = om
    sa = mapping.get("subgroup_analysis")
    if sa:
        sa = json.loads(json.dumps(sa))
        dp = sa.get("device_positive") or {}
        if dp.get("block") and not str(dp["block"]).startswith("device_"):
            dp["block"] = f"device_{dp['block']}"
        manifest["subgroup_analysis"] = sa
    info = {"data_use_source": du_source, "data_use_confirmation_required": confirm,
            "measurement_class_resolution": how, "caveats": caveats}
    return manifest, info


def _blocks(h: Harmonized) -> list[str]:
    out = [f"device_{n}" for n in h.devices]
    if h.labs is not None and len(h.labs):
        out.append("ehr_labs")
    if h.diagnoses is not None and len(h.diagnoses):
        out.append("ehr_diagnoses")
    if h.medications is not None and len(h.medications):
        out.append("ehr_medications")
    out += [f"survey_{n}" for n in h.surveys] + [f"omics_{n}" for n in h.omics]
    if h.history is not None and len(h.history):
        out.append("self_report_history")
    return out


# ------------------------------------------------------------------------------------------------------------------
# writing
# ------------------------------------------------------------------------------------------------------------------

def _int_day(df: pd.DataFrame) -> pd.DataFrame:
    if "day_index" in df.columns:
        df = df.copy()
        df["day_index"] = pd.to_numeric(df["day_index"], errors="coerce").astype("Int64")
        if df["day_index"].isna().all():
            df = df.drop(columns="day_index")
    return df


def write_byod(h: Harmonized, manifest: dict, byod_dir: Path) -> dict[str, int]:
    """Write the BYOD folder; returns {file name: rows}. Existing files of a previous apply are replaced."""
    byod_dir.mkdir(parents=True, exist_ok=True)
    for p in byod_dir.iterdir():           # a re-apply must not leave a stale block behind
        if p.is_file() and p.suffix in (".csv", ".parquet", ".yaml"):
            p.unlink()
    files: dict[str, int] = {}

    def put(name: str, df: pd.DataFrame) -> None:
        if df is None or df.empty:
            return
        df.to_csv(byod_dir / name, index=False)
        files[name] = int(len(df))

    P = h.participants.copy()
    for lab in h.labels:
        P[lab["column"]] = P[lab["column"]].map(lambda v: "" if pd.isna(v) else str(int(v)))
    put("participants.csv", P)
    for name, df in h.devices.items():
        put(f"device_{name}.csv", _int_day(df))
    if h.labs is not None and len(h.labs):
        labs = h.labs.drop(columns=["platform"])
        put("ehr_labs.csv", _int_day(labs))
    for name, df in (("ehr_diagnoses.csv", h.diagnoses), ("ehr_medications.csv", h.medications)):
        if df is not None and len(df):
            put(name, _int_day(df))
    if h.history is not None and len(h.history):
        hx = h.history.copy()
        if hx["icd10cm"].isna().all():
            hx = hx.drop(columns="icd10cm")
        put("self_report_history.csv", _int_day(hx))
    for name, df in h.surveys.items():
        put(f"survey_{name}.csv", _int_day(df))
    for name, df in h.omics.items():
        put(f"omics_{name}.csv", df)
    (byod_dir / "manifest.yaml").write_text(
        "# written by measure_it.agent.apply from an approved mapping (docs/AGENT_CONTRACT.md §3)\n"
        + yaml.safe_dump(manifest, sort_keys=False, width=110, allow_unicode=True))
    return files


def _gitignore(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    gi = out / ".gitignore"
    if not gi.exists():
        gi.write_text("# measure-it launch run: local-only (may hold private, user-supplied material); never commit\n*\n")


def approve_mapping(mapping: dict, approver: str | None = None) -> dict:
    """Stamp approved_by / approved_at (in place) and return the mapping."""
    if not is_approved(mapping):
        mapping["approved_by"] = approver or os.environ.get("MEASURE_IT_APPROVER") or getpass.getuser()
        mapping["approved_at"] = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    return mapping


def apply_mapping(tables: dict[str, pd.DataFrame], mapping: dict, out: Path, *, approve: bool, public: bool = False,
                  omop: bool = True) -> dict:
    """Apply an approved mapping: write <out>/byod/ (validated), <out>/omop/ and <out>/apply_report.json.

    Raises ApplyRefused when the mapping is not approved (and approve=False), has no usable tables, or the
    participant ids look identifying and no MEASURE_IT_BYOD_SALT is set; nothing is written then. A folder that fails
    the BYOD validator is still written and returned with `ok: False` and every validator reason.
    """
    out = Path(out)
    if not is_approved(mapping):
        if not approve:
            raise ApplyRefused(["mapping is not approved: set approved_by and approved_at in mapping.yaml, or pass "
                                "--approve (docs/AGENT_CONTRACT.md §0 rule 3)"])
        approve_mapping(mapping)
    h, audit = harmonize(tables, mapping)
    manifest, info = build_manifest(h, mapping, public=public)
    _gitignore(out)
    byod_dir = out / "byod"
    files = write_byod(h, manifest, byod_dir)
    rep = C.check_dir(byod_dir)
    val = rep.to_dict()
    du_errors = [e for e in rep.errors if "data_use" in e]
    report = {
        "ok": rep.ok, "dataset_id": h.dataset_id, "byod_dir": str(byod_dir), "files": files,
        "approved_by": mapping.get("approved_by"), "approved_at": mapping.get("approved_at"),
        "data_use_source": info["data_use_source"],
        "data_use_confirmation_required": info["data_use_confirmation_required"],
        "errors_other_than_data_use": [e for e in rep.errors if e not in du_errors],
        "measurement_class": manifest.get("measurement", {}).get("class"),
        "measurement_class_resolution": info["measurement_class_resolution"],
        "participant_id_transform": audit["participant_id_transform"],
        "n_participants": int(len(h.participants)),
        "labels": manifest["labels"], "comparator": manifest["comparator"],
        "primary_analysis": manifest["primary_analysis"],
        "dropped": audit["dropped"], "dropped_by_privacy_guard": audit["dropped_by_privacy_guard"],
        "unmapped": audit["unmapped"], "not_in_mapping": audit["not_in_mapping"],
        "values_dropped": audit["values_dropped"], "tables_skipped": audit["tables_skipped"],
        "warnings": h.warnings, "caveats": info["caveats"], "validation": val,
    }
    if omop:
        from .omop import write_omop
        try:
            report["omop"] = write_omop(h, out / "omop")
        except Exception as exc:  # noqa: BLE001 - the BYOD folder stands on its own; report the OMOP failure
            report["omop"] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
    (out / "mapping_applied.yaml").write_text(yaml.safe_dump(mapping, sort_keys=False, allow_unicode=True))
    (out / "apply_report.json").write_text(json.dumps(report, indent=1, default=str))
    return report


# ------------------------------------------------------------------------------------------------------------------
# a BYOD folder back into the bundle (for OMOP export of an existing folder)
# ------------------------------------------------------------------------------------------------------------------

def harmonized_from_byod(byod_dir: Path) -> Harmonized:
    """Read a BYOD folder (e.g. examples/byod_synthetic output) into a Harmonized bundle."""
    byod_dir = Path(byod_dir)
    m = K.load_manifest(byod_dir / "manifest.yaml")
    inp = K.discover(byod_dir)
    P = K.read_frame(inp.participants.path)
    for lab in m.get("labels") or []:
        if lab["column"] in P.columns:
            P[lab["column"]] = pd.to_numeric(P[lab["column"]], errors="coerce")
    h = Harmonized(dataset_id=m["dataset_id"], participants=P,
                   labels=[{k: lab.get(k) for k in ("column", "condition", "definition", "label_basis")}
                           for lab in m.get("labels") or []])
    h.device_meta = {k: dict(v or {}) for k, v in (m.get("devices") or {}).items()}
    h.omics_meta = {k: dict(v or {}) for k, v in (m.get("omics") or {}).items()}

    def num(df: pd.DataFrame) -> pd.DataFrame:
        for c in df.columns:
            if c != "participant_id":
                df[c] = pd.to_numeric(df[c], errors="coerce")
        return df

    for f in inp.files:
        if f.kind == "participants":
            continue
        df = K.read_frame(f.path)
        if f.kind == "device":
            h.devices[f.name] = num(df)
        elif f.kind == "survey":
            h.surveys[f.name] = num(df)
        elif f.kind == "omics":
            h.omics[f.name] = num(df)
        elif f.kind == "ehr_labs":
            df["value"] = pd.to_numeric(df["value"], errors="coerce")
            h.labs = df.reindex(columns=LAB_COLS)
        elif f.kind == "ehr_diagnoses":
            h.diagnoses = df.reindex(columns=DX_COLS)
        elif f.kind == "ehr_medications":
            h.medications = df.reindex(columns=RX_COLS)
        elif f.kind == "self_report_history":
            h.history = df.reindex(columns=HX_COLS)
    for attr in ("labs", "diagnoses", "medications", "history"):
        df = getattr(h, attr)
        if df is not None and "day_index" in df.columns:
            df["day_index"] = pd.to_numeric(df["day_index"], errors="coerce")
    return h
