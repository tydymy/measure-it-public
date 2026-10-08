"""Part B: the harmonized data as All of Us-shaped OMOP CDM v5.4 tables (docs/AGENT_CONTRACT.md §3).

    write_omop(harmonized, out_dir) -> summary dict          (called by measure_it.agent.apply.apply_mapping)
    export_omop_from_byod(byod_dir, out_dir) -> summary dict  (any existing `measure-it byod` folder)

Tables (parquet, one file each): person, observation_period, condition_occurrence, measurement, drug_exposure,
observation, concept; Fitbit-shaped activity_summary, heart_rate_summary, sleep_daily_summary when a wearable device
block has per-day rows. Column lists follow the All of Us curation schemas (github.com/all-of-us/curation,
data_steward/resource_files/schemas/cdm/clinical/*.json and wearables/fitbit/*.json, read 2026-10-07); omop/README.md
lists every extension column and every choice made here.

Concept ids: every (vocabulary_id, concept_code) used gets a LOCAL concept_id >= 2,000,000,000 (the OMOP range for
site-specific concepts), assigned in sorted (vocabulary_id, concept_code) order so a re-export is identical. Both
`*_concept_id` and `*_source_concept_id` hold these local ids. No standard OMOP concept id is written: inside All of Us
the same codes resolve through the CDR's own `concept` table (vocabulary_id + concept_code), which is how the stage-D
query packs (measure_it.similar.omop) look codes up.

Dates are synthetic: 2000-01-01 + day_index (day 0 when a row has no day_index). They carry order and spacing within a
participant, nothing else.
"""
from __future__ import annotations

import datetime as dt
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ANCHOR = dt.date(2000, 1, 1)
LOCAL_ID_START = 2_000_000_000
V_TYPE = "MEASURE_IT_TYPE"
V_SURVEY = "MEASURE_IT_SURVEY"
V_SELF = "MEASURE_IT_SELF_REPORT"
V_DEVICE = "MEASURE_IT_DEVICE"
V_DEMOG = "MEASURE_IT_DEMOGRAPHIC"
V_OMICS = "MEASURE_IT_OMICS"

PERSON_COLS = ["person_id", "gender_concept_id", "year_of_birth", "month_of_birth", "day_of_birth", "birth_datetime",
               "race_concept_id", "ethnicity_concept_id", "location_id", "provider_id", "care_site_id",
               "person_source_value", "gender_source_value", "gender_source_concept_id", "race_source_value",
               "race_source_concept_id", "ethnicity_source_value", "ethnicity_source_concept_id",
               # All of Us extension columns (person_ext in the curation pipeline; exposed on person in the Workbench)
               "sex_at_birth_concept_id", "sex_at_birth_source_concept_id", "sex_at_birth_source_value",
               # measure-it extension (not in OMOP or All of Us)
               "age_band_source_value"]
OBS_PERIOD_COLS = ["observation_period_id", "person_id", "observation_period_start_date",
                   "observation_period_end_date", "period_type_concept_id"]
CONDITION_COLS = ["condition_occurrence_id", "person_id", "condition_concept_id", "condition_start_date",
                  "condition_start_datetime", "condition_end_date", "condition_end_datetime",
                  "condition_type_concept_id", "condition_status_concept_id", "stop_reason", "provider_id",
                  "visit_occurrence_id", "visit_detail_id", "condition_source_value", "condition_source_concept_id",
                  "condition_status_source_value"]
MEASUREMENT_COLS = ["measurement_id", "person_id", "measurement_concept_id", "measurement_date",
                    "measurement_datetime", "measurement_time", "measurement_type_concept_id", "operator_concept_id",
                    "value_as_number", "value_as_concept_id", "unit_concept_id", "range_low", "range_high",
                    "provider_id", "visit_occurrence_id", "visit_detail_id", "measurement_source_value",
                    "measurement_source_concept_id", "unit_source_value", "value_source_value"]
DRUG_COLS = ["drug_exposure_id", "person_id", "drug_concept_id", "drug_exposure_start_date",
             "drug_exposure_start_datetime", "drug_exposure_end_date", "drug_exposure_end_datetime",
             "verbatim_end_date", "drug_type_concept_id", "stop_reason", "refills", "quantity", "days_supply", "sig",
             "route_concept_id", "lot_number", "provider_id", "visit_occurrence_id", "visit_detail_id",
             "drug_source_value", "drug_source_concept_id", "route_source_value", "dose_unit_source_value"]
OBSERVATION_COLS = ["observation_id", "person_id", "observation_concept_id", "observation_date",
                    "observation_datetime", "observation_type_concept_id", "value_as_number", "value_as_string",
                    "value_as_concept_id", "qualifier_concept_id", "unit_concept_id", "provider_id",
                    "visit_occurrence_id", "visit_detail_id", "observation_source_value",
                    "observation_source_concept_id", "unit_source_value", "qualifier_source_value",
                    "value_source_concept_id", "value_source_value", "questionnaire_response_id"]
CONCEPT_COLS = ["concept_id", "concept_name", "domain_id", "vocabulary_id", "concept_class_id", "standard_concept",
                "concept_code", "valid_start_date", "valid_end_date", "invalid_reason"]
# Fitbit tables (All of Us curation schemas wearables/fitbit/*.json; src_id omitted)
ACTIVITY_COLS = ["date", "activity_calories", "calories_bmr", "calories_out", "elevation", "fairly_active_minutes",
                 "floors", "lightly_active_minutes", "marginal_calories", "sedentary_minutes", "steps",
                 "very_active_minutes", "person_id"]
HR_SUMMARY_COLS = ["person_id", "date", "zone_name", "min_heart_rate", "max_heart_rate", "minute_in_zone",
                   "calorie_count", "mean_heart_rate_device"]          # last: measure-it extension
SLEEP_COLS = ["person_id", "sleep_date", "is_main_sleep", "minute_in_bed", "minute_to_fall_asleep", "minute_asleep",
              "minute_after_wakeup", "minute_awake", "minute_restless", "minute_deep", "minute_light", "minute_rem",
              "minute_wake"]
HR_ZONE = "device_daily_summary (not a Fitbit heart-rate zone)"

INT_COLS = re.compile(r"(_id|_concept_id|year_of_birth|month_of_birth|day_of_birth|refills|days_supply|steps|floors|"
                      r"min_heart_rate|max_heart_rate|minute_in_zone|^minute_.*)$")
DATE_COLS = re.compile(r"(_date|^date)$")

# per-day wearable columns -> Fitbit-shaped columns (patterns on the device feature name)
_NOT_HR = re.compile(r"hrv|rmssd|sdnn|variab", re.I)
FITBIT_RULES = [
    ("activity_summary", "steps", re.compile(r"(^|_)steps?($|_)|step_?count", re.I)),
    ("sleep_daily_summary", "minute_asleep",
     re.compile(r"sleep.*(min|dur|total|time)|total_?sleep|minutes?_asleep|(^|_)tst($|_)", re.I)),
    ("heart_rate_summary", "min_heart_rate",
     re.compile(r"(rest|resting|min|lowest).*(hr|heart)|(hr|heart_?rate).*(min|rest)|(^|_)rhr($|_)", re.I)),
    ("heart_rate_summary", "max_heart_rate", re.compile(r"(hr|heart_?rate).*max|max.*(hr|heart)", re.I)),
    ("heart_rate_summary", "mean_heart_rate_device",
     re.compile(r"(hr|heart_?rate).*(mean|avg|average)|(mean|avg|average).*(hr|heart)|^hr(_bpm)?$", re.I)),
]


class ConceptRegistry:
    """Collects every (vocabulary_id, concept_code) used; assigns local ids >= 2e9 at the end."""

    def __init__(self) -> None:
        self.rows: dict[tuple[str, str], dict] = {}

    def add(self, vocab: str, code: str, name: str | None = None, domain: str = "Observation",
            concept_class: str = "Local") -> tuple[str, str]:
        k = (str(vocab), str(code))
        r = self.rows.setdefault(k, {"concept_name": None, "domain_id": domain, "concept_class_id": concept_class})
        if name and not r["concept_name"]:
            r["concept_name"] = str(name)[:255]
        return k

    def ids(self) -> dict[tuple[str, str], int]:
        return {k: LOCAL_ID_START + i for i, k in enumerate(sorted(self.rows))}

    def frame(self) -> pd.DataFrame:
        ids = self.ids()
        rows = []
        for k, r in self.rows.items():
            rows.append({"concept_id": ids[k], "concept_name": r["concept_name"] or f"{k[0]} {k[1]} (name not looked up)",
                         "domain_id": r["domain_id"], "vocabulary_id": k[0], "concept_class_id": r["concept_class_id"],
                         "standard_concept": None, "concept_code": k[1], "valid_start_date": dt.date(1970, 1, 1),
                         "valid_end_date": dt.date(2099, 12, 31), "invalid_reason": None})
        return pd.DataFrame(rows, columns=CONCEPT_COLS).sort_values("concept_id").reset_index(drop=True)


def _date(day) -> dt.date:
    try:
        d = int(day)
    except (TypeError, ValueError):
        d = 0
    return ANCHOR + dt.timedelta(days=d)


def _finish(rows: list[dict], cols: list[str], reg_ids: dict, id_col: str | None) -> pd.DataFrame:
    """Resolve concept keys (tuples) to local ids, add the row id, order and type the columns."""
    df = pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)
    for c in cols:
        if c.endswith("concept_id") and len(df):
            df[c] = df[c].map(lambda v: reg_ids.get(v, 0) if isinstance(v, tuple) else v)
    if id_col and len(df):
        df[id_col] = np.arange(1, len(df) + 1)
    for c in cols:
        if DATE_COLS.search(c):
            df[c] = df[c].astype(object)
        elif c.endswith("datetime"):
            df[c] = pd.to_datetime(df[c], errors="coerce")
        elif INT_COLS.search(c):
            df[c] = pd.to_numeric(df[c], errors="coerce").round().astype("Int64")
        elif c in ("value_as_number", "quantity", "range_low", "range_high", "calorie_count", "activity_calories",
                   "calories_bmr", "calories_out", "elevation", "fairly_active_minutes", "lightly_active_minutes",
                   "marginal_calories", "sedentary_minutes", "very_active_minutes", "mean_heart_rate_device"):
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
        else:
            df[c] = df[c].astype("string")
    return df


def _ucum(unit):
    from ..similar.omop import _ucum as stage_d_ucum      # the same normalisation the query packs apply
    return stage_d_ucum(unit)


def _label_codes(condition: str) -> tuple[str | None, list[str], bool]:
    """Label condition -> (canonical id, ICD-10-CM codes of its base population, proxy?)."""
    from ..config import load_config
    from ..similar.phenotype_io import base_population
    cid = condition if condition in {c["id"] for c in load_config("conditions")["conditions"]} else None
    if cid is None:
        try:
            from ..byod.checks import resolve_condition
            cid, _ = resolve_condition(condition)
        except Exception:  # noqa: BLE001
            cid = None
    if not cid:
        return None, [], False
    bp = base_population({"condition_id": cid})
    return cid, bp["icd10cm"], bool(bp.get("proxy"))


def _wearable(meta: dict) -> bool:
    try:
        from ..config import load_config
        fam = {m["id"]: m.get("modality_family") for m in load_config("measurements")["measurement_classes"]}
        return fam.get((meta or {}).get("measurement_class")) == "wearable"
    except Exception:  # noqa: BLE001
        return False


def write_omop(h, out_dir: Path) -> dict:
    """Write the OMOP / All of Us-shaped export of a Harmonized bundle (measure_it.agent.apply.Harmonized)."""
    from ..harmonize.self_report import instrument_code, map_history
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    reg = ConceptRegistry()
    P = h.participants
    pids = sorted(P["participant_id"].astype(str))
    person_id = {p: i + 1 for i, p in enumerate(pids)}
    t_record = reg.add(V_TYPE, "study_record", "measure-it: value recorded by the study (local type concept)",
                       "Type Concept", "Type Concept")
    days: dict[int, list[int]] = {}

    def seen(pid: int, day) -> None:
        try:
            days.setdefault(pid, []).append(int(day))
        except (TypeError, ValueError):
            days.setdefault(pid, []).append(0)

    notes: list[str] = []
    # ---- person
    demo = h.demographics if h.demographics is not None else pd.DataFrame(columns=["participant_id", "variable", "value"])
    race = demo[demo["variable"].str.contains("race", na=False)].groupby("participant_id")["value"].first()
    eth = demo[demo["variable"].str.contains("ethnic|hispanic", na=False)].groupby("participant_id")["value"].first()
    persons = []
    for _, r in P.iterrows():
        pid = str(r["participant_id"])
        sex = str(r.get("sex") or "unknown") if "sex" in P.columns else "unknown"
        g = (reg.add("Gender", "F", "FEMALE", "Gender", "Gender") if sex == "female" else
             reg.add("Gender", "M", "MALE", "Gender", "Gender") if sex == "male" else
             reg.add(V_DEMOG, f"sex:{sex}", f"sex {sex} (measure-it local)", "Gender"))
        rv, ev = race.get(pid), eth.get(pid)
        rc = reg.add(V_DEMOG, f"race:{rv}", f"race: {rv} (as recorded by the study)", "Race") if rv else 0
        ec = reg.add(V_DEMOG, f"ethnicity:{ev}", f"ethnicity: {ev} (as recorded by the study)", "Ethnicity") if ev else 0
        persons.append({"person_id": person_id[pid], "gender_concept_id": g, "year_of_birth": None,
                        "race_concept_id": rc, "ethnicity_concept_id": ec, "person_source_value": pid,
                        "gender_source_value": sex, "gender_source_concept_id": g, "race_source_value": rv,
                        "race_source_concept_id": rc, "ethnicity_source_value": ev, "ethnicity_source_concept_id": ec,
                        "sex_at_birth_concept_id": g, "sex_at_birth_source_concept_id": g,
                        "sex_at_birth_source_value": sex,
                        "age_band_source_value": r.get("age_band") if "age_band" in P.columns else None})
    # ---- condition_occurrence: EHR diagnoses + study-label cases
    conds = []
    if h.diagnoses is not None:
        for _, r in h.diagnoses.iterrows():
            p = person_id.get(str(r["participant_id"]))
            if p is None:
                continue
            k = reg.add("ICD10CM", r["icd10cm"], None, "Condition", "ICD10 code")
            conds.append({"person_id": p, "condition_concept_id": k, "condition_source_concept_id": k,
                          "condition_start_date": _date(r.get("day_index")), "condition_type_concept_id": t_record,
                          "condition_source_value": r["icd10cm"]})
            seen(p, r.get("day_index"))
    t_label = reg.add(V_TYPE, "study_label", "measure-it: study case label recorded as the condition's ICD-10-CM "
                      "code (local type concept)", "Type Concept", "Type Concept")
    label_notes = {}
    for lab in h.labels:
        col = lab["column"]
        if col not in P.columns:
            continue
        cid, codes, proxy = _label_codes(lab["condition"])
        if not codes:
            label_notes[col] = f"condition {lab['condition']!r}: no ICD-10-CM code; cases not recorded as a condition"
            continue
        code = sorted(codes, key=lambda c: (-len(c), c))[0]
        label_notes[col] = (f"cases (label {col} = 1) recorded as condition_occurrence {code} with "
                            f"condition_type_concept_id = local 'study_label'" + (" (PROXY code)" if proxy else ""))
        k = reg.add("ICD10CM", code, None, "Condition", "ICD10 code")
        for _, r in P[pd.to_numeric(P[col], errors="coerce") == 1].iterrows():
            p = person_id[str(r["participant_id"])]
            conds.append({"person_id": p, "condition_concept_id": k, "condition_source_concept_id": k,
                          "condition_start_date": ANCHOR, "condition_type_concept_id": t_label,
                          "condition_source_value": code, "condition_status_source_value": f"study_label:{col}"})
            seen(p, 0)
    # ---- measurement: labs, Nightingale NMR (LOINC via the curated config), device features (local vocabulary)
    meas = []
    if h.labs is not None:
        for _, r in h.labs.iterrows():
            p = person_id.get(str(r["participant_id"]))
            if p is None or pd.isna(r.get("value")):
                continue
            k = reg.add("LOINC", r["loinc"], r.get("lab_name") if pd.notna(r.get("lab_name")) else None,
                        "Measurement", "Lab Test")
            unit = r.get("unit") if pd.notna(r.get("unit")) else None
            u = reg.add("UCUM", _ucum(unit), _ucum(unit), "Unit", "Unit") if unit else 0
            plat = r.get("platform") if "platform" in r and pd.notna(r.get("platform")) else None
            meas.append({"person_id": p, "measurement_concept_id": k, "measurement_source_concept_id": k,
                         "measurement_date": _date(r.get("day_index")), "measurement_type_concept_id": t_record,
                         "value_as_number": float(r["value"]), "unit_concept_id": u, "unit_source_value": unit,
                         "measurement_source_value": r["loinc"] + (f" | platform={plat}" if plat else "")})
            seen(p, r.get("day_index"))
    try:
        from ..harmonize.byod import nightingale_lookup, _norm
        ng = nightingale_lookup()
    except Exception:  # noqa: BLE001
        ng, _norm = {}, None
    omics_not_exported = {}
    for layer, df in h.omics.items():
        meta = h.omics_meta.get(layer, {})
        nmr = meta.get("platform") == "nmr_nightingale" and meta.get("units") == "nightingale_standard"
        mapped = []
        for c in df.columns:
            if c == "participant_id":
                continue
            row = ng.get(_norm(c)) if nmr and _norm else None
            if row is None:
                continue
            mapped.append(c)
            k = reg.add("LOINC", row["loinc"], row.get("loinc_name"), "Measurement", "Lab Test")
            u = reg.add("UCUM", _ucum(row["unit"]), row["unit"], "Unit", "Unit")
            note = (f"{c} | platform=nmr_nightingale; {row['source_unit']} x {row['factor']} -> {row['unit']} "
                    "(configs/harmonize_nightingale_loinc.yaml; NMR is not calibrated to routine clinical chemistry)")
            for pid, v in zip(df["participant_id"].astype(str), pd.to_numeric(df[c], errors="coerce")):
                p = person_id.get(pid)
                if p is None or pd.isna(v):
                    continue
                meas.append({"person_id": p, "measurement_concept_id": k, "measurement_source_concept_id": k,
                             "measurement_date": ANCHOR, "measurement_type_concept_id": t_record,
                             "value_as_number": float(v) * float(row["factor"]), "unit_concept_id": u,
                             "unit_source_value": row["unit"], "measurement_source_value": note,
                             "value_source_value": f"{v:g} {row['source_unit']}"})
                seen(p, 0)
        rest = [c for c in df.columns if c != "participant_id" and c not in mapped]
        if rest:
            omics_not_exported[layer] = len(rest)
    for name, df in h.devices.items():
        per_day = "day_index" in df.columns
        for c in df.columns:
            if c in ("participant_id", "day_index"):
                continue
            k = reg.add(V_DEVICE, f"{name}:{c}", f"device {name}: {c} (measure-it local)", "Measurement",
                        "Device feature")
            vals = pd.to_numeric(df[c], errors="coerce")
            dayv = df["day_index"] if per_day else pd.Series([0] * len(df), index=df.index)
            for pid, v, d in zip(df["participant_id"].astype(str), vals, dayv):
                p = person_id.get(pid)
                if p is None or pd.isna(v):
                    continue
                meas.append({"person_id": p, "measurement_concept_id": k, "measurement_source_concept_id": k,
                             "measurement_date": _date(d), "measurement_type_concept_id": t_record,
                             "value_as_number": float(v), "measurement_source_value": f"{name}:{c}"})
                seen(p, d)
    # ---- drug_exposure
    drugs = []
    if h.medications is not None:
        for _, r in h.medications.iterrows():
            p = person_id.get(str(r["participant_id"]))
            if p is None:
                continue
            k = reg.add("RxNorm", str(r["rxnorm"]), None, "Drug", "Ingredient")
            d = _date(r.get("day_index"))
            drugs.append({"person_id": p, "drug_concept_id": k, "drug_source_concept_id": k,
                          "drug_exposure_start_date": d, "drug_exposure_end_date": d,
                          "drug_type_concept_id": t_record, "drug_source_value": str(r["rxnorm"])})
            seen(p, r.get("day_index"))
    # ---- observation: survey answers, self-reported history, other demographics
    obs = []
    t_survey = reg.add(V_TYPE, "survey", "measure-it: participant-reported survey answer (local type concept)",
                       "Type Concept", "Type Concept")
    for name, df in h.surveys.items():
        instr = instrument_code(name)
        per_day = "day_index" in df.columns
        for c in df.columns:
            if c in ("participant_id", "day_index"):
                continue
            code = f"{instr}:{c}"
            k = reg.add(V_SURVEY, code, f"{instr} {c} (PPI-style question, measure-it local)", "Observation",
                        "Question")
            vals = pd.to_numeric(df[c], errors="coerce")
            dayv = df["day_index"] if per_day else pd.Series([0] * len(df), index=df.index)
            for pid, v, d in zip(df["participant_id"].astype(str), vals, dayv):
                p = person_id.get(pid)
                if p is None or pd.isna(v):
                    continue
                obs.append({"person_id": p, "observation_concept_id": k, "observation_source_concept_id": k,
                            "observation_date": _date(d), "observation_type_concept_id": t_survey,
                            "value_as_number": float(v), "observation_source_value": code})
                seen(p, d)
    if h.history is not None and len(h.history):
        q = reg.add(V_SELF, "personal_medical_history", "Personal medical history (self-reported; PPI-style, "
                    "measure-it local)", "Observation", "Question")
        for _, r in h.history.iterrows():
            p = person_id.get(str(r["participant_id"]))
            if p is None:
                continue
            icd = r.get("icd10cm") if pd.notna(r.get("icd10cm")) else None
            m = map_history(r.get("condition"), icd)
            a = reg.add("ICD10CM", m["icd10cm"], None, "Condition", "ICD10 code") if m["status"] == "mapped" else 0
            obs.append({"person_id": p, "observation_concept_id": q, "observation_source_concept_id": q,
                        "observation_date": _date(r.get("day_index")), "observation_type_concept_id": t_survey,
                        "value_as_string": str(r.get("condition"))[:80], "value_as_concept_id": a,
                        "value_source_concept_id": a, "value_source_value": m["icd10cm"],
                        "observation_source_value": "personal_medical_history",
                        "qualifier_source_value": f"self_report_map:{m['status']}"})
            seen(p, r.get("day_index"))
    other_demo = demo[~demo["variable"].str.contains("race|ethnic|hispanic", na=False)]
    for _, r in other_demo.iterrows():
        p = person_id.get(str(r["participant_id"]))
        if p is None:
            continue
        k = reg.add(V_DEMOG, str(r["variable"]), f"demographic {r['variable']} (measure-it local)", "Observation",
                    "Question")
        obs.append({"person_id": p, "observation_concept_id": k, "observation_source_concept_id": k,
                    "observation_date": ANCHOR, "observation_type_concept_id": t_survey,
                    "value_as_string": str(r["value"])[:80], "observation_source_value": str(r["variable"])})
    # ---- Fitbit-shaped tables from per-day wearable blocks
    act, hrs, slp = [], [], []
    fitbit_map: dict[str, dict] = {}
    fitbit_unmapped: dict[str, list[str]] = {}
    for name, df in h.devices.items():
        if "day_index" not in df.columns or not _wearable(h.device_meta.get(name, {})):
            continue
        cols = {}
        for c in df.columns:
            if c in ("participant_id", "day_index"):
                continue
            hit = None
            for table, target, rx in FITBIT_RULES:
                if rx.search(c) and not (table == "heart_rate_summary" and _NOT_HR.search(c)):
                    if (table, target) not in cols.values():
                        hit = (table, target)
                    break
            if hit:
                cols[c] = hit
            else:
                fitbit_unmapped.setdefault(name, []).append(c)
        if not cols:
            continue
        fitbit_map[name] = {c: f"{t}.{col}" for c, (t, col) in cols.items()}
        for _, r in df.iterrows():
            p = person_id.get(str(r["participant_id"]))
            if p is None or pd.isna(r["day_index"]):
                continue
            d = _date(r["day_index"])
            a, hr, s = {}, {}, {}
            for c, (t, col) in cols.items():
                v = pd.to_numeric(r[c], errors="coerce")
                if pd.isna(v):
                    continue
                {"activity_summary": a, "heart_rate_summary": hr, "sleep_daily_summary": s}[t][col] = float(v)
            if a:
                act.append({"person_id": p, "date": d, **a})
            if hr:
                hrs.append({"person_id": p, "date": d, "zone_name": HR_ZONE, **hr})
            if s:
                slp.append({"person_id": p, "sleep_date": d, "is_main_sleep": "true", **s})
    # ---- observation_period, then resolve ids and write
    periods = []
    for p in person_id.values():
        dd = days.get(p, [0])
        periods.append({"person_id": p, "observation_period_start_date": _date(min(dd)),
                        "observation_period_end_date": _date(max(dd)), "period_type_concept_id": t_record})
    ids = reg.ids()
    tables = {
        "person": _finish(persons, PERSON_COLS, ids, None),
        "observation_period": _finish(periods, OBS_PERIOD_COLS, ids, "observation_period_id"),
        "condition_occurrence": _finish(conds, CONDITION_COLS, ids, "condition_occurrence_id"),
        "measurement": _finish(meas, MEASUREMENT_COLS, ids, "measurement_id"),
        "drug_exposure": _finish(drugs, DRUG_COLS, ids, "drug_exposure_id"),
        "observation": _finish(obs, OBSERVATION_COLS, ids, "observation_id"),
    }
    if act:
        tables["activity_summary"] = _finish(act, ACTIVITY_COLS, ids, None)
    if hrs:
        tables["heart_rate_summary"] = _finish(hrs, HR_SUMMARY_COLS, ids, None)
    if slp:
        tables["sleep_daily_summary"] = _finish(slp, SLEEP_COLS, ids, None)
    tables["concept"] = reg.frame()
    for old in out_dir.glob("*.parquet"):
        old.unlink()
    for name, df in tables.items():
        df.to_parquet(out_dir / f"{name}.parquet", index=False)
    summary = {"ok": True, "dir": str(out_dir), "tables": {k: int(len(v)) for k, v in tables.items()},
               "anchor_date": ANCHOR.isoformat(), "local_concept_id_min": LOCAL_ID_START,
               "local_concept_id_max": int(tables["concept"]["concept_id"].max()) if len(tables["concept"]) else None,
               "vocabularies": tables["concept"]["vocabulary_id"].value_counts().to_dict(),
               "study_labels": label_notes, "fitbit_columns": fitbit_map, "fitbit_unmapped_columns": fitbit_unmapped,
               "omics_features_not_exported": omics_not_exported, "notes": notes}
    (out_dir / "README.md").write_text(readme(summary))
    (out_dir / "omop_export.json").write_text(json.dumps(summary, indent=1, default=str))
    return summary


def export_omop_from_byod(byod_dir: Path, out_dir: Path) -> dict:
    """OMOP export of an existing BYOD folder (no demographics beyond age band and sex)."""
    from .apply import harmonized_from_byod
    return write_omop(harmonized_from_byod(byod_dir), out_dir)


def readme(s: dict) -> str:
    t = s["tables"]
    rows = "\n".join(f"| `{k}` | {v} |" for k, v in t.items())
    fit = "\n".join(f"* device `{d}`: " + ", ".join(f"`{c}` -> `{v}`" for c, v in m.items())
                    for d, m in s["fitbit_columns"].items()) or "* none (no per-day wearable device block)"
    unf = "; ".join(f"`{d}`: {', '.join(cs)}" for d, cs in s["fitbit_unmapped_columns"].items()) or "none"
    labels = "\n".join(f"* `{k}`: {v}" for k, v in s["study_labels"].items()) or "* no labels"
    return f"""# OMOP CDM v5.4 export, shaped like the All of Us Curated Data Repository

Written by `measure_it.agent.omop` from the same harmonized data as the BYOD folder next to it
(docs/AGENT_CONTRACT.md §3). Local-only: this folder sits under a run folder whose `.gitignore` is `*`.

| table | rows |
|---|---|
{rows}

## Dates

Every date is synthetic: **{s['anchor_date']} + day_index** (day 0 when the source row had no day index). Dates keep
the order and spacing of one participant's records and nothing else; they are not calendar dates and must never be
read as such. `year_of_birth` is **null** for everyone (only a 10-year age band exists, by design); the band is in the
non-OMOP extension column `person.age_band_source_value` (e.g. `40-49`, `90+`).

## Concepts: what is local and what is standard

* Every `(vocabulary_id, concept_code)` used is a row of the local `concept` table with a **local concept_id >=
  {s['local_concept_id_min']:,}** (the OMOP range reserved for site-specific concepts), assigned in sorted order.
  Both `*_concept_id` and `*_source_concept_id` hold these local ids; `standard_concept` is null.
* **No standard OMOP concept id is invented.** The codes themselves are real where the vocabulary is real:
  `ICD10CM` (diagnoses), `LOINC` (labs), `RxNorm` (ingredient CUIs), `UCUM` (units), `Gender` (`F`/`M`).
  Inside All of Us (or any CDM with Athena vocabularies) the same codes resolve through the CDR's `concept` table by
  `vocabulary_id` + `concept_code`, which is exactly how the stage-D query packs look codes up.
* `MEASURE_IT_*` vocabularies are local and have no OMOP counterpart: `MEASURE_IT_TYPE` (type concepts:
  `study_record`, `study_label`, `survey`), `MEASURE_IT_SURVEY` (PRO items, code `<INSTRUMENT>:<item>`, e.g.
  `COMPASS31:vasomotor`), `MEASURE_IT_SELF_REPORT` (`personal_medical_history`), `MEASURE_IT_DEVICE`
  (`<device>:<feature>`), `MEASURE_IT_DEMOGRAPHIC` (race / ethnicity / other study categories as recorded; also sex
  values other than female / male).
* Vocabularies used: {json.dumps(s['vocabularies'])}.

## Table by table

* `person`: one row per participant; `person_source_value` = the study code written to the BYOD folder (a salted hash
  when the agent hashed identifying ids). `gender_concept_id`, `sex_at_birth_concept_id` (All of Us extension; in the
  curation pipeline it lives in `person_ext`, the Workbench exposes it on `person`) and their source columns all carry
  the study's single `sex` field (female / male / other / unknown): the study did not distinguish gender from sex at
  birth. Race / ethnicity only when the mapping had such a demographic column. No location, no year of birth.
* `observation_period`: first to last synthetic date of the person's rows.
* `condition_occurrence`: EHR diagnoses (ICD-10-CM source codes, type `study_record`) and **study-label cases**:
{labels}
  These rows make a query pack's base population (people with the condition's ICD-10-CM code) equal to the study's
  case group; exclude them with `condition_type_concept_id` = the local `study_label` concept if you need EHR-coded
  diagnoses only. Controls get no condition row.
* `measurement`: LOINC labs (`value_as_number`, `unit_source_value` as supplied, `unit_concept_id` = local UCUM concept
  of the unit normalised as the stage-D packs do); Nightingale NMR measures mapped to LOINC by
  `configs/harmonize_nightingale_loinc.yaml` with the explicit unit factor, the platform and the conversion stated in
  `measurement_source_value` and the original value + unit in `value_source_value` (**NMR is not calibrated to routine
  clinical chemistry**); device features as local `MEASURE_IT_DEVICE` concepts (not an All of Us domain; a CDM query
  that selects LOINC concepts never sees them). Other omics features are not exported (no OMOP domain):
  {json.dumps(s['omics_features_not_exported'])}.
* `drug_exposure`: RxNorm ingredient CUIs (start = end = the synthetic date).
* `observation`: survey answers as All of Us PPI-style question / answer rows (`observation_source_value` =
  `<INSTRUMENT>:<item>`, numeric answer in `value_as_number`); self-reported medical history as
  `personal_medical_history` with the reported condition in `value_as_string` and, when the curated alias table maps
  it, the ICD-10-CM code in `value_source_value` / `value_as_concept_id` (like the All of Us personal medical history
  survey, it is NOT a condition_occurrence); other demographic categories.
* Fitbit-shaped tables (`activity_summary`, `heart_rate_summary`, `sleep_daily_summary`), column names and types from
  the All of Us curation schemas (`data_steward/resource_files/schemas/wearables/fitbit/*.json` in
  github.com/all-of-us/curation, read 2026-10-07; `src_id` omitted). Built from per-day rows of wearable device
  blocks; **the study device is not a Fitbit** (other sensors and algorithms, values not interchangeable):
{fit}
  `heart_rate_summary` has one row per day with `zone_name` = "{HR_ZONE}"; the daily mean heart rate is in the
  extension column `mean_heart_rate_device` (Fitbit has no such column). `sleep_daily_summary.is_main_sleep` is
  `'true'` because the device reports one daily total. Per-day columns with no Fitbit counterpart (HRV, SpO2, ...): {unf}.
  The Workbench minute-level heart-rate table (`heart_rate_minute_level`; `heart_rate_intraday` in curation) is not
  produced: there are no minute-level data.

## Running the same analysis inside All of Us

1. Fit the subgroup / computable phenotype locally (`measure-it byod subgroup` on the BYOD folder).
2. `measure-it similar` writes the All of Us BigQuery query pack (`aou_bigquery__<id>.sql`); in the Researcher
   Workbench replace `` `{{CDR}}` `` with your CDR dataset (the `WORKSPACE_CDR` environment variable) and run it. It
   returns aggregate counts only, with All of Us small-cell suppression (1-20) built in.
3. `query_pack_check.json` in this folder (when present) is the result of running that same pack, transpiled to DuckDB,
   against THIS export: it shows the SQL executes and how many rows each step keeps here. Two differences from the
   real CDR: (a) local concept ids instead of standard ones (stubs for `concept_relationship` 'Maps to' and
   `concept_ancestor` map each local concept to itself); (b) `year_of_birth` is null here, so the check also runs a
   variant with a check-only year of birth derived from the age-band midpoint.
4. Self-reported history and PRO answers are observations here and in All of Us; the stage-D packs evaluate
   ICD-10-CM features only from `condition_occurrence`, so a phenotype built on self-reported history is evaluated in
   All of Us from coded diagnoses, which differ (under-coding of post-infectious conditions is common).
"""
