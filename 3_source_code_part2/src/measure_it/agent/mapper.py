"""Propose `mapping.yaml` from a profile (docs/AGENT_CONTRACT.md section 2), and the part-A command functions.

    mapping = propose_mapping(profile, dataset=ds, condition_hint="long_covid")
    write_mapping(mapping, out_dir)                      # <out>/mapping.yaml (human-editable; needs approval)
    profile_cmd(path, out);  map_cmd(path, out, llm=None|"claude", condition_hint=..., measurement_hint=...)

Three methods, in order, each recorded per column in ``proposed_by`` (method, confidence 0-1, one-line rationale):

1. ``rule`` - names, flags and layouts: participant-id candidates, age / sex / dates, identifier and geography columns
   (dropped), code columns (ICD-10-CM / LOINC / RxNorm values), REDCap structure (checkbox / yes-no history items,
   forms named after an instrument or a device, repeating instruments), long tables (code + value + unit columns),
   platform hints (Olink, SomaScan, LEGENDplex, Nightingale).
2. ``dictionary`` (exact, after normalisation) and ``fuzzy`` (rapidfuzz token_sort_ratio >= threshold, default 0.80)
   against vocabularies already in the repo: configs/harmonize_nightingale_loinc.yaml and harmonize_nhanes_loinc.yaml
   (LOINC), harmonize_self_report_icd10.yaml (ICD-10-CM), harmonize_survey_crosswalk.yaml (instruments),
   measurements.yaml / relevance.yaml (measurement classes and bundles), conditions.yaml + measure_it.ontology.normalize
   (cohort / label columns).
3. optional ``llm`` (measure_it.agent.llm) for columns still unmapped; every proposal is validated here
   (`validate_proposal`): code format and check digit, LOINC / ICD-10-CM present in a local list, RxNorm format,
   role in the contract's list, unit plausible. A failure leaves the column unmapped with the reason.

Nothing is invented: an unknown column stays ``unmapped`` with the reason.
"""
from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import pandas as pd
import yaml
from rapidfuzz import fuzz, process

from ..config import CONFIGS, PROCESSED, load_config
from . import llm as L
from .profile import icd10_ok, luhn_ok, profile_dataset, write_profile
from .read import Dataset, _code_str, read_input

DEFAULT_THRESHOLD = 0.80
LLM_MIN_CONFIDENCE = 0.6
ROLES = L.ROLES
MIN_FUZZY_LEN = 4

LABEL_NAMES = re.compile(r"^(cohort|cohort_name|cohort_group|group|arm|study_arm|study_group|case_control|casecontrol|"
                         r"case_status|diagnosis_group|dx_group|label|status|condition|disease|disease_group|phenotype|"
                         r"participant_group|participant_type|group_assignment|enrollment_group|enrolment_group|"
                         r"subject_group|population)$")
SEX_NAMES = re.compile(r"^(sex|gender|sex_at_birth|biological_sex|birth_sex|sex_assigned_at_birth|riagendr|"
                       r"gender_source_value|gender_concept_id|sex_cd|participant_sex)$")
AGE_BAND_NAMES = re.compile(r"^(age_band|ageband|age_group|agegroup|age_cat|age_category|age_range|age_bracket|"
                            r"age_decade|age_bin)$")
DEMOG_NAMES = re.compile(r"(^|_)(race|ethnicity|ethnic|hispanic|education|educ|marital|employment|income|"
                         r"race_concept_id|ethnicity_concept_id|race_source_value|ethnicity_source_value)(_|$)")
DISPLAY_NAMES = re.compile(r"(_text|_display|_description|_desc|_name)$")
LAB_CODE_NAMES = re.compile(r"^(loinc|loinc_code|loinc_num|test_code|lbtestcd|testcd|lbtest|analyte|test|test_name|"
                            r"lab_test|lab_name|measurement_source_value|parameter|component|assay|lab)$")
LAB_VALUE_NAMES = re.compile(r"^(value|result|value_as_number|lborres|lbstresn|aval|measurement_value|result_value|"
                             r"lab_value|numeric_value|value_numeric)$")
LAB_UNIT_NAMES = re.compile(r"^(unit|units|lborresu|lbstresu|unit_source_value|value_unit|result_unit|uom|unit_of_measure)$")
HISTORY_TEXT_NAMES = re.compile(r"(self_?reported_)?(condition|diagnosis|medical_history|history|reported_condition|"
                                r"condition_name)")
UNIT_STRIP = re.compile(r"\s*[\(\[][^\)\]]*[\)\]]\s*|_(mg_?dl|mmol_?l|umol_?l|nmol_?l|g_?l|g_?dl|mg_?l|ng_?ml|pg_?ml|"
                        r"pct|percent|bpm|ms|kg|cm|mmhg|npx)$", re.IGNORECASE)
DATE_PREFER = re.compile(r"visit|collect|effective|sample|assess|measure|exam|start|onset|enrol|baseline|record")
DEVICE_CONTEXT = re.compile(r"ring|wearable|oura|fitbit|evie|movano|watch|garmin|whoop|actigraph|accelerom|daily|"
                            r"capillar|nailfold|caprio|device|sensor|wear")
DEVICE_PREFIX = re.compile(r"^(ring|wear|wearable|oura|fitbit|evie|watch|garmin|whoop|cap|nvc|nailfold|capillary|"
                           r"capillaroscopy|caprio|stat|nasa_lean|lean)_")
CAP_TOKENS = re.compile(r"^(cap|nvc|nailfold|capillary|capillaroscopy|caprio)_|capillar|nailfold")
WEARABLE_TOKENS = [
    (re.compile(r"rmssd|sdnn|(^|_)hrv(_|$)|heart_rate_variability"), "hrv"),
    (re.compile(r"spo2|oxygen_sat|o2_?sat|(^|_)sat(_|$)"), "continuous_spo2"),
    (re.compile(r"(^|_)(rhr|hr|bpm|pulse|heart_rate|resting_hr|resting_heart_rate|hr_rest|hr_mean|hr_max|hr_min)(_|$)"),
     "wearable_heart_rate"),
    (re.compile(r"sleep|(^|_)(rem|deep|light|awake|wake|waso|tst)(_|$)"), "sleep_objective"),
    (re.compile(r"steps|(^|_)(activity|active|active_min|met|mets|calories|kcal|sedentary)(_|$)"), "accelerometry"),
    (re.compile(r"(^|_)(temp|temperature|skin_temp|temp_dev|temperature_deviation)(_|$)"), "continuous_temperature"),
    (re.compile(r"resp_?rate|(^|_)(rr|breathing_rate|respiratory_rate|resp)(_|$)"), "respiratory_rate"),
]
OLINK_RE = re.compile(r"olink|(^|_)npx(_|$)|^oid\d{5}$", re.IGNORECASE)
SOMA_RE = re.compile(r"somascan|somalogic|^seq\.\d+\.\d+$", re.IGNORECASE)
LEGEND_RE = re.compile(r"legendplex|legend_plex", re.IGNORECASE)
NMR_CONTEXT = re.compile(r"nightingale|(^|_)nmr(_|$)|metabolom", re.IGNORECASE)
OMOP_ROW_ID = re.compile(r"_(occurrence|exposure|era|period)_id$|^(measurement|observation|visit_occurrence|visit_detail|"
                         r"note|specimen|death)_id$")
CONTROL_RE = re.compile(r"^(hc|hcs|healthy|healthy controls?|healthy volunteers?|controls?|ctrl|ctl|hv|non ?cases?|"
                        r"unaffected|normal|reference|comparator|0)$")
COHORT_ABBREV = {"lc": "long_covid", "pasc": "long_covid", "longcovid": "long_covid", "long covid": "long_covid",
                 "post covid": "long_covid", "clyme": "ptlds", "ptlds": "ptlds", "chronic lyme": "ptlds",
                 "post treatment lyme": "ptlds", "post lyme": "ptlds", "plds": "ptlds", "alyme": "lyme_disease",
                 "acute lyme": "lyme_disease", "lyme": "lyme_disease", "mecfs": "me_cfs", "me cfs": "me_cfs",
                 "cfs": "me_cfs", "pots": "pots"}
YES = {"1", "yes", "y", "true", "t", "checked", "present", "positive"}
NO = {"0", "no", "n", "false", "f", "unchecked", "absent", "negative"}
UNKNOWN_ANSWERS = {"unknown", "dont know", "don t know", "do not know", "prefer not to say", "prefer not to answer",
                   "not sure", "unsure", "declined", "refused", "missing", "na", "n a", "not applicable"}
SEX_WORDS = {"female": "female", "f": "female", "woman": "female", "w": "female", "girl": "female",
             "male": "male", "m": "male", "man": "male", "boy": "male",
             "other": "other", "non binary": "other", "nonbinary": "other", "intersex": "other", "x": "other",
             "diverse": "other", "unknown": "unknown", "prefer not to say": "unknown", "prefer not to answer": "unknown",
             "not reported": "unknown", "u": "unknown", "declined": "unknown"}
OMOP_GENDER = {"8532": "female", "8507": "male"}          # OMOP standard concepts FEMALE / MALE (Gender vocabulary)


def nk(s) -> str:
    """Normalised key: lower case, apostrophes removed, other non-alphanumerics -> one space."""
    s = str(s or "").lower().replace("'", "").replace("’", "")
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def ck(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s or "").lower())


def slug(s) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(s or "").lower()).strip("_")


# ------------------------------------------------------------------------------------------------------------------
# vocabularies already in the repo
# ------------------------------------------------------------------------------------------------------------------
class Vocab:
    """Lazy view of the repo's vocabularies (configs + processed ontology tables when present)."""

    def __init__(self, processed: Path = PROCESSED, configs: Path = CONFIGS):
        self.processed = Path(processed)
        self.configs = Path(configs)

    @cached_property
    def nightingale(self) -> dict:
        cfg = load_config("harmonize_nightingale_loinc")
        out = {}
        for r in cfg["mappings"]:
            for a in [r["measure"], *r.get("aliases", [])]:
                out.setdefault(ck(a), r)
        return out

    @cached_property
    def nightingale_unmapped(self) -> set:
        return {ck(r["measure"]) for r in load_config("harmonize_nightingale_loinc").get("unmapped", [])
                if re.fullmatch(r"[A-Za-z0-9_]+", str(r["measure"]))}

    @cached_property
    def nhanes(self) -> list[dict]:
        return load_config("harmonize_nhanes_loinc")["mappings"]

    @cached_property
    def nhanes_by_variable(self) -> dict:
        out = {}
        for r in self.nhanes:
            out.setdefault(ck(r["variable"]), r)
        return out

    @cached_property
    def lab_lexicon(self) -> list[tuple[str, dict, str]]:
        """(normalised text, row, source) for fuzzy LOINC matching: Nightingale measures and aliases, NHANES LOINC
        component names (LONG_COMMON_NAME before the first '[' / ' in ')."""
        lex, seen = [], set()
        for r in load_config("harmonize_nightingale_loinc")["mappings"]:
            for a in [r["measure"], *r.get("aliases", [])]:
                k = nk(a.replace("_", " "))
                if k and (k, r["loinc"]) not in seen:
                    seen.add((k, r["loinc"]))
                    lex.append((k, r, "configs/harmonize_nightingale_loinc.yaml"))
        for r in self.nhanes:
            comp = re.split(r"\s\[|\sin\s(?:Serum|Plasma|Blood|Urine|Red)", r["loinc_name"])[0]
            k = nk(comp)
            if k and (k, r["loinc"]) not in seen:
                seen.add((k, r["loinc"]))
                lex.append((k, r, "configs/harmonize_nhanes_loinc.yaml"))
        return lex

    @cached_property
    def loinc_units(self) -> dict[str, set]:
        out: dict[str, set] = {}
        for r in self.nhanes:
            out.setdefault(r["loinc"], set()).update({r.get("unit"), *r.get("accept_units", [])} - {None})
        for r in load_config("harmonize_nightingale_loinc")["mappings"]:
            out.setdefault(r["loinc"], set()).update({r.get("unit"), r.get("source_unit")} - {None})
        return out

    @cached_property
    def loinc_codes(self) -> set:
        codes = set()
        for p in sorted(self.configs.glob("*.yaml")):
            codes.update(re.findall(r"\bloinc:\s*[\"']?(\d{1,7}-\d)\b", p.read_text()))
        pc = self.processed / "person_concepts.parquet"
        if pc.exists():
            try:
                d = pd.read_parquet(pc, columns=["vocabulary", "concept_code"])
                codes.update(d.loc[d["vocabulary"] == "LOINC", "concept_code"].astype(str).unique())
            except Exception:  # noqa: BLE001 - an unreadable optional table only narrows the list
                pass
        return {c for c in codes if luhn_ok(c)}

    @cached_property
    def icd_codes(self) -> tuple[set, str]:
        p = self.processed / "ontology_icd10cm_codes.parquet"
        if p.exists():
            try:
                d = pd.read_parquet(p, columns=["code_nodot"])
                return set(d["code_nodot"].astype(str).str.upper()), "data/processed/ontology_icd10cm_codes.parquet"
            except Exception:  # noqa: BLE001
                pass
        codes = set()
        for p in sorted(self.configs.glob("*.yaml")):
            codes.update(c.replace(".", "").upper() for c in re.findall(r"\bicd10cm:\s*[\"']?([A-Z]\d[0-9A-Z](?:\.[0-9A-Z]{1,4})?)",
                                                                      p.read_text()))
        return codes, "ICD-10-CM codes in configs/*.yaml (ontology_icd10cm_codes not built)"

    @cached_property
    def history_aliases(self) -> list[tuple[str, dict]]:
        out = []
        for r in load_config("harmonize_self_report_icd10")["mappings"]:
            for a in r["aliases"]:
                out.append((nk(a), r))
        return sorted(out, key=lambda x: -len(x[0]))

    @cached_property
    def instruments(self) -> list[dict]:
        out = []
        for r in load_config("harmonize_survey_crosswalk")["instruments"]:
            keys = sorted({slug(r["code"]), *[slug(f) for f in r["file_names"]]}, key=len, reverse=True)
            out.append({"code": r["code"], "label": r["label"], "keys": keys,
                        "label_key": ck(re.split(r"[(\[]", r["label"])[0])})
        return out

    @cached_property
    def instrument_codes(self) -> set:
        return {i["code"] for i in self.instruments}

    @cached_property
    def classes(self) -> dict:
        out = {}
        for m in load_config("measurements")["measurement_classes"]:
            pats = []
            for p in m.get("patterns", []) or []:
                try:
                    pats.append(re.compile(p, re.IGNORECASE))
                except re.error:
                    continue
            out[m["id"]] = {"name": m["name"], "aliases": m.get("aliases") or [], "patterns": pats,
                            "family": m.get("modality_family")}
        return out

    @cached_property
    def bundles(self) -> dict:
        return load_config("relevance").get("measurement_bundles", {})

    @cached_property
    def device_lexicon(self) -> list[tuple[str, str, str]]:
        """(normalised alias, class id, source) for device classes (not lab / omics families)."""
        lex = []
        for cid, c in self.classes.items():
            if c["family"] in ("lab", "omics"):
                continue
            for a in [cid.replace("_", " "), re.sub(r"\(.*?\)", "", c["name"]), *c["aliases"]]:
                if nk(a):
                    lex.append((nk(a), cid, "configs/measurements.yaml"))
        for bid, b in self.bundles.items():
            members = [m for m in b.get("members", []) if self.classes.get(m, {}).get("family") not in ("lab", "omics")]
            if not members:
                continue
            for a in [bid.replace("_", " "), *b.get("aliases", [])]:
                if nk(a):
                    lex.append((nk(a), members[0], "configs/relevance.yaml"))
        return lex

    @cached_property
    def conditions(self) -> dict[str, str]:
        """compact key -> canonical condition id (ids, preferred names, aliases, search terms)."""
        out = {}
        for c in load_config("conditions")["conditions"]:
            for a in [c["id"], c["preferred_name"], re.sub(r"\(.*?\)", "", c["preferred_name"]),
                      *(c.get("aliases") or []), *(c.get("search_terms") or [])]:
                k = ck(re.sub(r"\(.*?\)", "", str(a)))
                if k:
                    out.setdefault(k, c["id"])
        return out

    @cached_property
    def condition_ids(self) -> set:
        return {c["id"] for c in load_config("conditions")["conditions"]}


# ------------------------------------------------------------------------------------------------------------------
# small resolvers
# ------------------------------------------------------------------------------------------------------------------
def resolve_condition_text(text: str, vocab: Vocab) -> tuple[str | None, str]:
    """Cohort level / label text -> ('control' | condition id | None, how)."""
    t = nk(text)
    if not t:
        return None, "empty"
    if CONTROL_RE.match(t):
        return "control", "control wording"
    if t in COHORT_ABBREV:
        return COHORT_ABBREV[t], f"cohort abbreviation '{text}' (curated in measure_it.agent.mapper.COHORT_ABBREV)"
    if ck(text) in COHORT_ABBREV:
        return COHORT_ABBREV[ck(text)], f"cohort abbreviation '{text}'"
    cid = vocab.conditions.get(ck(re.sub(r"\(.*?\)", "", str(text))))
    if cid:
        return cid, "configs/conditions.yaml name or alias"
    if re.search(r"\bhealthy\b|\bcontrols?\b", t):
        return "control", "control wording"
    try:
        from ..ontology.normalize import normalize_condition
        r = normalize_condition(str(text))
        cands = sorted({m.get("canonical_condition_id") for m in r.get("matches", []) or [] if m.get("canonical_condition_id")})
        if r.get("status") == "matched" and len(cands) == 1:
            return cands[0], f"measure_it.ontology.normalize ({r['matches'][0].get('match_reason', '')})"
    except Exception:  # noqa: BLE001 - the ontology tables are optional here
        pass
    return None, "no condition match"


def resolve_condition_hint(hint: str | None, vocab: Vocab) -> str | None:
    if not hint:
        return None
    if hint in vocab.condition_ids:
        return hint
    cid, _ = resolve_condition_text(hint, vocab)
    if cid and cid != "control":
        return cid
    try:
        from ..byod.checks import resolve_condition
        cid, _ = resolve_condition(hint)
        return cid
    except Exception:  # noqa: BLE001
        return None


def yes_no_map(values: list[str], choices: dict | None) -> dict | None:
    """Codes / levels of a yes-no item -> {value: 1 | 0 | None}; None when the column is not yes-no."""
    src = {str(k): str(v) for k, v in (choices or {}).items()} or {v: v for v in values}
    if not src:
        return None
    out = {}
    for code, lab in src.items():
        k = nk(lab)
        if k in YES:
            out[code] = 1
        elif k in NO:
            out[code] = 0
        elif k in UNKNOWN_ANSWERS:
            out[code] = None
        else:
            return None
    return out if 1 in out.values() or 0 in out.values() else None


def sex_map(values: list[str], choices: dict | None, col: str) -> tuple[dict, list[str]]:
    if col.lower() == "gender_concept_id":
        vm = {v: OMOP_GENDER[v] for v in values if v in OMOP_GENDER}
        return vm, [v for v in values if v not in OMOP_GENDER]
    src = {str(k): str(v) for k, v in (choices or {}).items()} or {v: v for v in values}
    vm, unknown = {}, []
    for code, lab in src.items():
        s = SEX_WORDS.get(nk(lab))
        if s:
            vm[code] = s
        else:
            unknown.append(code)
    return vm, unknown


def looks_like_bands(values: list[str]) -> bool:
    vals = [v for v in values if v]
    return bool(vals) and all(re.fullmatch(r"\d{1,2}\s*[-–]\s*\d{1,2}|\d{2}\s*\+", v) for v in vals)


# ------------------------------------------------------------------------------------------------------------------
# the proposal engine
# ------------------------------------------------------------------------------------------------------------------
@dataclass
class Proposal:
    mapping: dict
    method: str
    confidence: float
    rationale: str
    no_llm: bool = False         # a rule decided it stays unmapped (do not ask a model)


def _p(mapping, method, confidence, rationale, no_llm=False) -> Proposal:
    return Proposal(mapping, method, round(float(confidence), 3), rationale, no_llm)


class Mapper:
    def __init__(self, profile: dict, dataset: Dataset | None = None, *, condition_hint: str | None = None,
                 measurement_hint: str | None = None, threshold: float = DEFAULT_THRESHOLD, vocab: Vocab | None = None):
        self.profile = profile
        self.dataset = dataset
        self.vocab = vocab or Vocab()
        self.threshold = threshold
        self.condition_hint = condition_hint
        self.condition_id = resolve_condition_hint(condition_hint, self.vocab)
        self.measurement_hint = measurement_hint

    # local (never sent anywhere) distinct values of a low-cardinality column, for value maps
    def local_levels(self, table: str, col: dict) -> list[str]:
        vals: list[str] = []
        if self.dataset is not None:
            try:
                s = self.dataset.table(table).df[col["name"]].dropna()
                u = s.map(_code_str).unique()
                if len(u) <= 60:
                    vals = sorted(map(str, u))
            except KeyError:
                pass
        if not vals:
            vals = [v for v, _ in (col.get("levels") or []) if not str(v).startswith("<")]
        for k in (col.get("choices") or {}):
            if k not in vals:
                vals.append(str(k))
        return vals

    # ---------------------------------------------------------------------------------------------- lab matching
    def lab_exact(self, text: str) -> tuple[dict, str] | None:
        k = ck(UNIT_STRIP.sub("", str(text or "")))
        if k in self.vocab.nightingale:
            return self.vocab.nightingale[k], "configs/harmonize_nightingale_loinc.yaml"
        if k in self.vocab.nhanes_by_variable:
            return self.vocab.nhanes_by_variable[k], "configs/harmonize_nhanes_loinc.yaml"
        return None

    def fuzzy(self, texts: list[str], lexicon: list[tuple]) -> tuple[float, tuple | None, str]:
        choices = [e[0] for e in lexicon]
        best, best_e, best_t = 0.0, None, ""
        for t in texts:
            q = nk(UNIT_STRIP.sub("", str(t)).replace("_", " "))
            if len(q) < MIN_FUZZY_LEN or not choices:
                continue
            r = process.extractOne(q, choices, scorer=fuzz.token_sort_ratio)
            if r and r[1] / 100 > best:
                best, best_e, best_t = r[1] / 100, lexicon[r[2]], t
        return best, best_e, best_t

    def lab_entry(self, row: dict, col: dict, source: str, nmr: bool) -> tuple[dict | None, str]:
        """A lab mapping with a plausible unit, or (None, why)."""
        unit = col.get("unit")
        e = {"role": "lab", "loinc": row["loinc"]}
        loinc_unit = row.get("unit")
        if "source_unit" in row:                              # Nightingale row
            if unit and unit not in (row["source_unit"], loinc_unit):
                return None, (f"unit {unit} matches neither the Nightingale unit {row['source_unit']} nor the LOINC "
                              f"unit {loinc_unit} in {source}")
            if nmr and (unit in (None, row["source_unit"])):
                e.update({"unit": row["source_unit"], "platform": "nmr_nightingale", "loinc_unit": loinc_unit,
                          "factor": row["factor"]})
            else:
                e["unit"] = unit or None
                if unit is None:
                    e["unit_note"] = "no unit in the source metadata; fill before apply"
        else:
            accept = {loinc_unit, *row.get("accept_units", [])}
            if unit and unit not in accept:
                return None, f"unit {unit} differs from the LOINC term's unit {loinc_unit} in {source} (no conversion)"
            e["unit"] = unit or None
            if not unit:
                e["unit_note"] = (f"no unit in the source metadata; the LOINC term's unit in {source} is {loinc_unit}: "
                                  "confirm before apply")
        e["loinc_name"] = row.get("loinc_name")
        return e, ""

    # ---------------------------------------------------------------------------------------------- per column
    def propose_column(self, t: dict, col: dict, ctx: dict) -> Proposal:
        name = col["name"]
        low = slug(name)
        flags = set(col.get("flags") or [])
        label = col.get("label") or ""
        dtype = col.get("dtype")
        choices = col.get("choices") or None

        # -- privacy first
        if "geography" in flags:
            return _p({"role": "drop", "reason": "geography"}, "rule", 0.99, "geography pattern (measure_it.byod.checks)")
        if "identifier" in flags:
            return _p({"role": "drop", "reason": "identifier"}, "rule", 0.99,
                      "identifier pattern (measure_it.byod.checks) or REDCap 'Identifier?' = y")
        if "exact_age" in flags:
            return _p({"role": "age", "transform": "age_band"}, "rule", 0.95, "exact age -> 10-year age_band")
        if "date" in flags:
            if dtype in ("date", "datetime"):
                return _p({"role": "date"}, "rule", 0.95, "calendar date -> per-participant day_index (never written)")
            return _p({"role": "drop", "reason": "date"}, "rule", 0.9, "date-like column that is not a date value")
        if "identifier?" in flags or "geography?" in flags:
            why = "identifier" if "identifier?" in flags else "geography"
            return _p({"role": "drop", "reason": why}, "rule", 0.6,
                      f"the label suggests {why}; dropped for safety, review and change if wrong")
        if "id_values_look_identifying" in flags and name == ctx.get("participant_id"):
            return _p({"role": "participant_id", "hash_ids": True}, "rule", 0.8,
                      "participant-id candidate whose values look like MRNs / identifying numbers: hash them "
                      "(id_hashing: salted_sha256)")

        # -- table structure
        if name == ctx.get("participant_id"):
            return _p({"role": "participant_id"}, "rule", 0.97, f"participant-id candidate ({t['kind'] or t['shape']} table)")
        if "structural" in flags or OMOP_ROW_ID.search(low):
            return _p({"role": "drop", "reason": "structural"}, "rule", 0.95,
                      "REDCap / OMOP structural column (event, repeat instance, form status, row key)")
        if low.endswith("_concept_id") and low not in ("gender_concept_id", "race_concept_id", "ethnicity_concept_id"):
            return _p({"role": "unmapped", "reason": "OMOP concept id without a local concept table; the "
                       "*_source_value column carries the source code"}, "rule", 0.0, "OMOP concept id", no_llm=True)
        if "free_text" in flags:
            if dtype == "category" and re.search(r"condition|diagnos|history", low):
                pass
            else:
                return _p({"role": "drop", "reason": "free_text"}, "rule", 0.95, "free text is never kept")

        # -- self-reported condition names in a long table
        if dtype == "category" and t["shape"] != "wide" and HISTORY_TEXT_NAMES.fullmatch(low):
            return _p({"role": "self_report_history", "from": "text"}, "rule", 0.85,
                      "condition names as reported; mapped per value through configs/harmonize_self_report_icd10.yaml "
                      "(exact aliases) on apply")

        # -- demographics
        if SEX_NAMES.match(low) or (re.search(r"\b(sex|gender)\b", label.lower()) and dtype in ("category", "bool", "int")
                                    and col.get("n_unique", 0) <= 6):
            vm, unk = sex_map(ctx["levels"], choices, name)
            m = {"role": "sex", "value_map": vm or None}
            why = "sex / gender column"
            if unk:
                m["unmapped_values"] = unk
                why += f"; values {unk} have no label: fill value_map by hand"
            return _p(m, "rule", 0.95 if vm and not unk else 0.7, why)
        if AGE_BAND_NAMES.match(low) or (dtype == "category" and looks_like_bands(ctx["levels"]) and "age" in low):
            return _p({"role": "age_band"}, "rule", 0.9, "age band column (bands are widened to decades on apply)")
        if dtype in ("category", "bool", "int") and (LABEL_NAMES.match(low) or re.search(r"\bcohort\b|study group|case.?control",
                                                                                         label.lower())):
            p = self.label_proposal(col, ctx)
            if p:
                return p
        if DEMOG_NAMES.search(low):
            return _p({"role": "demographic"}, "rule", 0.85, "demographic column (race / ethnicity / education ...)")

        # -- code columns and long tables
        cp = col.get("code_pattern")
        ll = ctx.get("long_lab") or {}
        if name in (ll.get("code_column"), ll.get("value_column"), ll.get("unit_column")):
            part = "code" if name == ll["code_column"] else "value" if name == ll["value_column"] else "unit"
            m = {"role": "lab", "part": part}
            if part == "code" and ll.get("code_map"):
                m["code_map"] = ll["code_map"]
                if ll.get("unmapped_levels"):
                    m["unmapped_levels"] = ll["unmapped_levels"]
            return _p(m, ll.get("method", "rule"), ll.get("confidence", 0.95), ll["rationale"])
        if cp == "icd10":
            return _p({"role": "diagnosis_code"}, "rule", 0.95, "values are ICD-10-CM codes")
        if cp == "rxnorm":
            return _p({"role": "medication_code"}, "rule", 0.9, "values are RxNorm codes (supply ingredient CUIs)")
        if cp == "snomed":
            return _p({"role": "unmapped", "reason": "SNOMED CT codes: no local SNOMED CT -> ICD-10-CM map"
                       + ("; the ICD-10-CM column of this table is used" if ctx.get("has_icd") else "")},
                      "rule", 0.0, "SNOMED CT code column", no_llm=True)
        if cp == "loinc":
            return _p({"role": "unmapped", "reason": "LOINC code column without a numeric value column in this table"},
                      "rule", 0.0, "LOINC code column", no_llm=True)
        if DISPLAY_NAMES.search(low) and (t.get("kind") or "").startswith(("fhir:", "omop:")):
            return _p({"role": "drop", "reason": "free_text"}, "rule", 0.9, "code display text (the code column is kept)")
        if (t.get("kind") or "").startswith("fhir:") and low in ("resource_type", "status", "clinical_status", "category",
                                                               "value_code"):
            return _p({"role": "drop", "reason": "structural"}, "rule", 0.8, "FHIR resource status / category")

        # -- surveys
        p = self.survey_proposal(col, ctx)
        if p:
            return p

        # -- self-reported history (yes / no item naming a condition)
        yn = yes_no_map(ctx["levels"], choices) if dtype in ("bool", "category", "int") else None
        if yn is not None:
            p = self.history_proposal(col, yn)
            if p:
                return p

        # -- omics platforms
        p = self.omics_proposal(col, ctx)
        if p:
            return p

        # -- labs: exact vocabulary names
        if dtype in ("int", "float"):
            for text in (name, label):
                hit = self.lab_exact(text) if text else None
                if hit:
                    e, why = self.lab_entry(hit[0], col, hit[1], ctx["nmr"])
                    if e:
                        conf = 0.9 if hit[0].get("confidence") != "low" else 0.8
                        return _p(e, "dictionary", conf, f"exact name match to {hit[1]} '{hit[0].get('measure') or hit[0].get('variable')}'")
                    return _p({"role": "unmapped", "reason": why}, "dictionary", 0.0, why, no_llm=True)

        # -- devices
        p = self.device_proposal(col, ctx)
        if p:
            return p

        # -- condition-named 0/1 column -> label
        if yn is not None:
            for text in (name, label):
                cid, how = resolve_condition_text(text, self.vocab) if text else (None, "")
                if cid and cid != "control" and nk(text) not in ("me",):
                    vm = {k: v for k, v in yn.items()}
                    return _p({"role": "label", "condition": cid, "value_map": vm,
                               "label_basis": None}, "dictionary", 0.8,
                              f"0/1 column named after a condition ({how}); set label_basis (the data owner knows it)")

        # -- fuzzy: labs, devices, instruments
        texts = [x for x in (name, label) if x]
        best = []
        if dtype in ("int", "float"):
            s, e, tx = self.fuzzy(texts, self.vocab.lab_lexicon)
            best.append((s, "lab", e, tx))
        s, e, tx = self.fuzzy(texts, self.vocab.device_lexicon)
        best.append((s, "device", e, tx))
        best.sort(key=lambda x: -x[0])
        s, kind, e, tx = best[0]
        if e is not None and s >= self.threshold:
            if kind == "lab":
                ent, why = self.lab_entry(e[1], col, e[2], ctx["nmr"])
                if ent:
                    return _p(ent, "fuzzy", min(s, 0.9), f"'{tx}' ~ '{e[0]}' in {e[2]} (score {s:.2f})")
                return _p({"role": "unmapped", "reason": why}, "fuzzy", 0.0, why)
            if dtype in ("int", "float", "bool"):
                return _p({"role": "device", "block": ctx["device_block"] or e[1], "measurement_class": e[1]}, "fuzzy",
                          min(s, 0.85), f"'{tx}' ~ '{e[0]}' in {e[2]} (score {s:.2f})")
        why = (f"no rule; best fuzzy score {s:.2f}" + (f" ('{tx}' ~ '{e[0]}')" if e is not None else "")
               + f" < {self.threshold:.2f}")
        return _p({"role": "unmapped", "reason": why}, "none", 0.0, why)

    # ---------------------------------------------------------------------------------------------- sub-proposals
    def label_proposal(self, col: dict, ctx: dict) -> Proposal | None:
        choices = col.get("choices") or {}
        levels = ctx["levels"]
        src = {str(k): str(v) for k, v in choices.items()} if choices else {v: v for v in levels}
        if not src or len(src) > 12:
            return None
        lc, hows = {}, []
        for code, lab in src.items():
            cid, how = resolve_condition_text(lab, self.vocab)
            lc[code] = cid
            if cid:
                hows.append(f"{lab}->{cid}")
        cases = [c for c in lc.values() if c and c != "control"]
        if not cases:
            return None
        counts = {}
        for v, n in (col.get("levels") or []):
            counts[str(v)] = n or 0
        if self.condition_id and self.condition_id in cases:
            primary = self.condition_id
            why = f"condition_hint {self.condition_hint!r}"
        else:
            freq: dict[str, int] = {}
            for code, cid in lc.items():
                if cid and cid != "control":
                    freq[cid] = freq.get(cid, 0) + int(counts.get(code, 0) or 0)
            primary = max(sorted(freq), key=lambda k: freq[k])
            why = "most frequent case level" + (f" (condition_hint {self.condition_hint!r} matches no level)"
                                                if self.condition_hint else "")
        vm = {code: (1 if cid == primary else 0 if cid == "control" else None) for code, cid in lc.items()}
        m = {"role": "label", "condition": primary, "value_map": vm, "label_basis": None,
             "level_conditions": {code: cid for code, cid in lc.items()}}
        unresolved = [code for code, cid in lc.items() if cid is None]
        conf = 0.9 if not unresolved else 0.7
        rationale = (f"cohort / group column; levels {', '.join(hows)}; primary condition from {why}; other conditions "
                     f"are null (not in this label); set label_basis (the data owner knows it)")
        if unresolved:
            rationale += f"; unresolved levels {unresolved}"
        return _p(m, "rule", conf, rationale)

    def history_proposal(self, col: dict, yn: dict) -> Proposal | None:
        label = col.get("label") or ""
        texts = [label, col["name"]]
        m = re.search(r"\(choice=([^)]*)\)", label)
        if m:
            texts.insert(0, m.group(1))
        from ..harmonize.self_report import map_history
        for i, text in enumerate(texts):
            if not text:
                continue
            r = map_history(text)
            if r["status"] == "mapped":
                return _p({"role": "self_report_history", "icd10cm": r["icd10cm"], "value_map": yn,
                           "label": r["label"]}, "dictionary", 0.9,
                          f"exact alias '{text}' in configs/harmonize_self_report_icd10.yaml (confidence {r['confidence']})")
        for text in texts:
            k = f" {nk(text)} "
            for alias, row in self.vocab.history_aliases:
                if len(alias) >= 5 and f" {alias} " in k:
                    return _p({"role": "self_report_history", "icd10cm": row["icd10cm"], "value_map": yn,
                               "label": row["label"]}, "dictionary", 0.85,
                              f"'{text}' contains the alias '{alias}' (configs/harmonize_self_report_icd10.yaml)")
        return None

    def survey_proposal(self, col: dict, ctx: dict) -> Proposal | None:
        low = slug(col["name"])
        form = slug(col.get("form") or "")
        tname = slug(ctx["table_name"])
        for inst in self.vocab.instruments:
            for key in inst["keys"]:
                if form and (form == key or form.startswith(key + "_") or form.endswith("_" + key)):
                    item = re.sub(rf"^{re.escape(key)}_?", "", low) or low
                    return _p({"role": "survey", "instrument": inst["code"], "item": item}, "rule", 0.9,
                              f"REDCap form '{col.get('form')}' is the instrument {inst['code']}")
                if low == key or re.match(rf"^{re.escape(key)}(_|$)", low):
                    item = low[len(key):].strip("_") or "total"
                    return _p({"role": "survey", "instrument": inst["code"], "item": item}, "dictionary", 0.85,
                              f"column prefix '{key}' = instrument {inst['code']} (configs/harmonize_survey_crosswalk.yaml)")
                if tname == key or tname.startswith(key + "_"):
                    return _p({"role": "survey", "instrument": inst["code"], "item": low}, "rule", 0.8,
                              f"table '{ctx['table_name']}' is the instrument {inst['code']}")
            lk = inst["label_key"]
            if col.get("label") and len(lk) >= 3 and ck(col["label"]).startswith(ck(inst["code"])):
                return _p({"role": "survey", "instrument": inst["code"], "item": low}, "dictionary", 0.8,
                          f"label starts with the instrument name {inst['code']}")
        return None

    def omics_proposal(self, col: dict, ctx: dict) -> Proposal | None:
        name, label = col["name"], col.get("label") or ""
        if col.get("dtype") not in ("int", "float"):
            return None
        hay = f"{name} {label} {ctx['context']}"
        if OLINK_RE.search(name) or OLINK_RE.search(label) or OLINK_RE.search(ctx["context"]) or col.get("unit") == "NPX":
            return _p({"role": "omics", "layer": "proteomics", "platform": "olink"}, "rule", 0.9,
                      "Olink platform hint (name / label / form / NPX unit)")
        if SOMA_RE.search(hay):
            return _p({"role": "omics", "layer": "proteomics", "platform": "somascan"}, "rule", 0.9, "SomaScan platform hint")
        if LEGEND_RE.search(hay):
            return _p({"role": "omics", "layer": "cytokines", "platform": "legendplex"}, "rule", 0.85, "LEGENDplex hint")
        if NMR_CONTEXT.search(ctx["context"]) and not self.lab_exact(name):
            return _p({"role": "omics", "layer": "metabolomics", "platform": "nmr_nightingale"}, "rule", 0.85,
                      "Nightingale / NMR table or form; not a clinical-chemistry-equivalent measure "
                      "(configs/harmonize_nightingale_loinc.yaml), so it stays OMICS")
        return None

    def device_proposal(self, col: dict, ctx: dict) -> Proposal | None:
        name, label = col["name"], col.get("label") or ""
        if col.get("dtype") not in ("int", "float", "bool"):
            return None
        low = slug(name)
        in_ctx = bool(DEVICE_CONTEXT.search(ctx["context"]) or DEVICE_PREFIX.match(low))
        block = ctx["device_block"]
        if CAP_TOKENS.search(low) or (in_ctx and re.search(r"capillar|nailfold|caprio", ctx["context"])):
            return _p({"role": "device", "block": block or "capillaroscopy", "measurement_class": "capillaroscopy"},
                      "rule", 0.9, "capillaroscopy feature (prefix / form)")
        if in_ctx:
            stem = DEVICE_PREFIX.sub("", low)
            for rx, cid in WEARABLE_TOKENS:
                if rx.search(stem):
                    return _p({"role": "device", "block": block or "wearable", "measurement_class": cid}, "rule", 0.85,
                              f"wearable feature token in '{name}' ({ctx['table_name']} is a device table)")
        hay = f"{name.replace('_', ' ')} {label}"
        for cid, c in self.vocab.classes.items():
            if c["family"] in ("lab", "omics"):
                continue
            for rx in c["patterns"]:
                if rx.search(hay):
                    return _p({"role": "device", "block": block or cid, "measurement_class": cid}, "dictionary", 0.8,
                              f"configs/measurements.yaml pattern for {cid} matches '{hay.strip()}'")
        if in_ctx and block:
            hint = self.measurement_hint if self.measurement_hint in self.vocab.classes else None
            if hint is None and self.measurement_hint in self.vocab.bundles:
                hint = (self.vocab.bundles[self.measurement_hint].get("members") or [None])[0]
            if hint:
                return _p({"role": "device", "block": block, "measurement_class": hint}, "rule", 0.6,
                          f"numeric feature of the device table '{ctx['table_name']}'; class from measurement_hint")
        return None

    # ---------------------------------------------------------------------------------------------- per table
    def long_lab(self, t: dict) -> dict | None:
        cols = {c["name"]: c for c in t["columns"]}
        code = next((c for c in t["columns"] if c["code_pattern"] == "loinc"), None)
        if code is None:
            code = next((c for c in t["columns"] if LAB_CODE_NAMES.match(slug(c["name"]))
                         and c["dtype"] in ("category", "text") and not (set(c["flags"]) & L.SENSITIVE_FLAGS)), None)
        val = next((c for c in t["columns"] if LAB_VALUE_NAMES.match(slug(c["name"])) and c["dtype"] in ("int", "float")),
                   None)
        if code is None or val is None or t["shape"] != "long":
            return None
        unit = next((c for c in t["columns"] if LAB_UNIT_NAMES.match(slug(c["name"]))), None)
        out = {"code_column": code["name"], "value_column": val["name"], "unit_column": unit["name"] if unit else None,
               "code_system": "loinc"}
        if code["code_pattern"] == "loinc":
            out["rationale"] = "long lab table: LOINC code + numeric value (+ unit) columns"
            return out
        levels = self.local_levels(t["name"], code)
        lmap, miss = {}, []
        for lv in levels:
            hit = self.lab_exact(lv)
            if hit:
                lmap[lv] = hit[0]["loinc"]
                continue
            s, e, _ = self.fuzzy([lv], self.vocab.lab_lexicon)
            if e is not None and s >= self.threshold:
                lmap[lv] = e[1]["loinc"]
            else:
                miss.append(lv)
        if not lmap:
            return None
        out.update({"code_map": lmap, "unmapped_levels": miss, "method": "fuzzy", "confidence": 0.8,
                    "rationale": f"long lab table: test names mapped to LOINC through configs (exact / fuzzy >= "
                                 f"{self.threshold:.2f}); {len(miss)} name(s) left unmapped"})
        _ = cols
        return out

    def table_context(self, t: dict) -> dict:
        names = [c["name"] for c in t["columns"]]
        forms = {c.get("form") for c in t["columns"] if c.get("form")}
        kind = t.get("kind") or ""
        tname = t["name"]
        n_ng = sum(1 for n in names if ck(n) in self.vocab.nightingale or ck(n) in self.vocab.nightingale_unmapped)
        nmr_table = bool(NMR_CONTEXT.search(f"{tname} {t.get('sheet') or ''} {kind}")) or n_ng >= 5
        device_block = None
        if DEVICE_CONTEXT.search(f"{tname} {kind} {t.get('sheet') or ''}"):
            device_block = slug(tname)
        ids = t.get("id_candidates") or []
        return {"participant_id": ids[0] if ids else None, "nmr_table": nmr_table, "device_block": device_block,
                "table_name": tname, "forms": forms, "has_icd": any(c["code_pattern"] == "icd10" for c in t["columns"]),
                "long_lab": self.long_lab(t)}

    def map_table(self, t: dict) -> tuple[dict, dict, list]:
        tctx = self.table_context(t)
        columns, proposed, pending = {}, {}, []
        for col in t["columns"]:
            form = col.get("form") or ""
            context = f"{t['name']} {t.get('kind') or ''} {t.get('sheet') or ''} {form}".strip()
            nmr = tctx["nmr_table"] or bool(NMR_CONTEXT.search(form))
            dev_block = tctx["device_block"]
            if not dev_block and form and DEVICE_CONTEXT.search(form):
                dev_block = slug(form)
            ctx = {**tctx, "context": context, "nmr": nmr, "device_block": dev_block,
                   "levels": self.local_levels(t["name"], col) if col.get("dtype") in ("bool", "category", "int") else []}
            p = self.propose_column(t, col, ctx)
            columns[col["name"]] = p.mapping
            proposed[col["name"]] = {"method": p.method, "confidence": p.confidence, "rationale": p.rationale}
            if p.mapping["role"] == "unmapped" and not p.no_llm:
                pending.append(col["name"])
        entry: dict = {"source_file": t["source_file"]}
        if t.get("sheet"):
            entry["sheet"] = t["sheet"]
        if t.get("kind"):
            entry["kind"] = t["kind"]
        entry["shape"] = t["shape"]
        pid = next((c for c, m in columns.items() if m["role"] == "participant_id"), None)
        entry["participant_id"] = pid
        dates = [c for c, m in columns.items() if m["role"] == "date"]
        if dates:
            entry["day_index_from"] = sorted(dates, key=lambda c: (0 if DATE_PREFER.search(c.lower()) else 1,
                                                                   dates.index(c)))[0]
        names = [c["name"] for c in t["columns"]]
        if "redcap_event_name" in names:
            entry["event_column"] = "redcap_event_name"
        if "redcap_repeat_instance" in names:
            entry["repeat_instance_column"] = "redcap_repeat_instance"
        if tctx["long_lab"]:
            entry["long_lab"] = {k: v for k, v in tctx["long_lab"].items()
                                 if k in ("code_column", "value_column", "unit_column", "code_system")}
        entry["columns"] = columns
        entry["proposed_by"] = proposed
        return entry, proposed, pending

    # ---------------------------------------------------------------------------------------------- LLM
    def validate_proposal(self, prop: dict, col: dict) -> tuple[dict | None, str]:
        """Deterministic validation of one model proposal -> (mapping, '') or (None, reason)."""
        role = prop.get("role")
        if role not in ROLES:
            return None, f"role {role!r} is not in the contract's list"
        if role not in L.LLM_ROLES:
            return None, f"role {role} is assigned only by rules"
        system = (prop.get("code_system") or "none").lower()
        code = str(prop.get("code") or "").strip()
        dtype = col.get("dtype")
        if role == "unmapped":
            return None, "model left it unmapped: " + str(prop.get("rationale", ""))[:200]
        if role == "drop":
            return {"role": "drop", "reason": "llm_proposed"}, ""
        if role == "demographic":
            return {"role": "demographic"}, ""
        if role == "lab":
            if dtype not in ("int", "float"):
                return None, "lab proposal for a non-numeric column"
            if system != "loinc" or not code:
                return None, "lab proposal without a LOINC code"
            if not luhn_ok(code):
                return None, f"{code} is not a well-formed LOINC code (format / check digit)"
            if code not in self.vocab.loinc_codes:
                return None, f"LOINC {code} is not in the local list (configs/*.yaml, person_concepts)"
            units = self.vocab.loinc_units.get(code)
            unit = col.get("unit")
            if unit and units and unit not in units:
                return None, f"unit {unit} is not a unit of LOINC {code} in configs ({sorted(units)}); no conversion"
            m = {"role": "lab", "loinc": code, "unit": unit or None}
            if not unit:
                m["unit_note"] = ("no unit in the source metadata" + (f"; units of this LOINC code in configs: "
                                                                     f"{sorted(units)}" if units else "")
                                  + ": confirm before apply")
            return m, ""
        if role == "self_report_history":
            yn = yes_no_map(self.local_levels_for(col), col.get("choices"))
            if yn is None:
                return None, "self_report_history proposal for a column that is not a yes/no item"
            if system != "icd10cm" or not code:
                return None, "self_report_history proposal without an ICD-10-CM code"
            c = code.upper()
            if not icd10_ok(c):
                return None, f"{code} is not a well-formed ICD-10-CM code"
            codes, src = self.vocab.icd_codes
            if c.replace(".", "") not in codes:
                return None, f"ICD-10-CM {code} is not in the local list ({src})"
            if "." not in c and len(c) > 3:
                c = c[:3] + "." + c[3:]
            return {"role": "self_report_history", "icd10cm": c, "value_map": yn}, ""
        if role == "diagnosis_code":
            if col.get("code_pattern") != "icd10":
                return None, "diagnosis_code proposal but the values are not ICD-10-CM codes"
            return {"role": "diagnosis_code"}, ""
        if role == "medication_code":
            if col.get("code_pattern") != "rxnorm":
                return None, "medication_code proposal but the values are not RxNorm codes"
            return {"role": "medication_code"}, ""
        if role == "medication_flag":
            yn = yes_no_map(self.local_levels_for(col), col.get("choices"))
            if yn is None:
                return None, "medication_flag proposal for a column that is not a yes/no item"
            if system != "rxnorm" or not re.fullmatch(r"\d{1,9}", code):
                return None, "medication_flag proposal without a well-formed RxNorm code"
            return {"role": "medication_flag", "rxnorm": code, "value_map": yn,
                    "note": "RxNorm code format checked only (no local RxNorm list)"}, ""
        if role == "survey":
            inst = L.sanitize_token(prop.get("instrument"))
            if not inst:
                return None, "survey proposal without an instrument"
            inst_u = inst.upper()
            known = next((i["code"] for i in self.vocab.instruments if ck(i["code"]) == ck(inst_u)
                          or ck(inst_u) in [ck(k) for k in i["keys"]]), None)
            item = L.sanitize_token(prop.get("item")) or slug(col["name"])
            if dtype not in ("int", "float", "bool"):
                return None, "survey items must be numeric (code text answers first)"
            m = {"role": "survey", "instrument": known or re.sub(r"[^A-Z0-9_]", "_", inst_u), "item": item}
            if not known:
                m["note"] = "instrument not in configs/harmonize_survey_crosswalk.yaml (no crosswalk)"
            return m, ""
        if role == "device":
            mc = prop.get("measurement_class")
            if mc not in self.vocab.classes:
                return None, f"measurement_class {mc!r} is not in configs/measurements.yaml"
            if self.vocab.classes[mc]["family"] in ("lab", "omics"):
                return None, f"{mc} is a lab / omics class, not a device"
            return {"role": "device", "block": slug(mc), "measurement_class": mc}, ""
        if role == "omics":
            layer = (prop.get("omics_layer") or "").lower()
            if layer not in L.OMICS_LAYERS:
                return None, f"omics_layer {layer!r} is not one of {L.OMICS_LAYERS}"
            if dtype not in ("int", "float"):
                return None, "omics features must be numeric"
            m = {"role": "omics", "layer": layer}
            plat = L.sanitize_token(prop.get("platform"), 40)
            if plat:
                m["platform"] = plat.lower()
            return m, ""
        return None, f"role {role} not handled"

    def local_levels_for(self, col: dict) -> list[str]:
        return self.local_levels(col.get("_table", ""), col)

    def run_llm(self, provider, tables: dict, out) -> dict:
        items, where = [], {}
        prof_tables = {t["name"]: t for t in self.profile["tables"]}
        for tname, entry in tables.items():
            t = prof_tables[tname]
            cols = {c["name"]: c for c in t["columns"]}
            for cname in entry.pop("_pending", []):
                key = f"c{len(items) + 1}"
                it = L.column_item(key, t, cols[cname])
                if it is None:
                    continue
                items.append(it)
                where[key] = (tname, cname)
        info = {"provider": getattr(provider, "name", None), "model": getattr(provider, "model", None),
                "n_columns_sent": len(items), "n_proposals": 0, "n_accepted": 0, "n_rejected": 0, "errors": [],
                "log": "llm_requests.jsonl"}
        if not items:
            return info
        context = {"format": self.profile.get("format"), "condition_hint": self.condition_hint,
                   "measurement_hint": self.measurement_hint, "instruments": sorted(self.vocab.instrument_codes),
                   "measurement_classes": sorted(k for k, v in self.vocab.classes.items()
                                                 if v["family"] not in ("lab", "omics")),
                   "omics_layers": L.OMICS_LAYERS}
        props, errors = L.run_proposals(provider, items, context, out)
        info["errors"] = errors
        info["n_proposals"] = len(props)
        for key, (tname, cname) in where.items():
            entry = tables[tname]
            col = {**{c["name"]: c for c in prof_tables[tname]["columns"]}[cname], "_table": tname}
            prop = props.get(key)
            if prop is None:
                old = entry["columns"][cname].get("reason", "")
                entry["columns"][cname]["reason"] = old + ("; llm: request failed" if errors else "; llm: no proposal")
                continue
            try:
                conf = max(0.0, min(1.0, float(prop.get("confidence", 0))))
            except (TypeError, ValueError):
                conf = 0.0
            m, why = self.validate_proposal(prop, col)
            if m is not None and conf < LLM_MIN_CONFIDENCE:
                m, why = None, f"model confidence {conf:.2f} < {LLM_MIN_CONFIDENCE}"
            rationale = str(prop.get("rationale", ""))[:300]
            if m is None:
                info["n_rejected"] += 1
                entry["columns"][cname] = {"role": "unmapped", "reason": f"llm proposal rejected: {why}"}
                entry["proposed_by"][cname] = {"method": "llm", "confidence": 0.0,
                                               "rationale": f"rejected ({why}); model said: {rationale}"}
            else:
                info["n_accepted"] += 1
                entry["columns"][cname] = m
                entry["proposed_by"][cname] = {"method": "llm", "confidence": round(conf, 3),
                                               "rationale": f"validated; model said: {rationale}"}
        return info

    # ---------------------------------------------------------------------------------------------- all
    def propose(self, *, dataset_id: str, llm=None, out=None, echo=print) -> dict:
        tables = {}
        for t in self.profile["tables"]:
            entry, _, pending = self.map_table(t)
            entry["_pending"] = pending
            tables[t["name"]] = entry
        llm_info = None
        provider = L.get_provider(llm, echo=echo)
        if provider is not None:
            llm_info = self.run_llm(provider, tables, out)
        for e in tables.values():
            e.pop("_pending", None)
        summary: dict = {"by_method": {}, "by_role": {}, "n_columns": 0}
        for e in tables.values():
            for c, m in e["columns"].items():
                meth = e["proposed_by"][c]["method"] if m["role"] != "unmapped" else "unmapped"
                summary["by_method"][meth] = summary["by_method"].get(meth, 0) + 1
                summary["by_role"][m["role"]] = summary["by_role"].get(m["role"], 0) + 1
                summary["n_columns"] += 1
        return {
            "dataset_id": dataset_id,
            "condition_hint": self.condition_hint,
            "condition_id": self.condition_id,
            "measurement_hint": self.measurement_hint,
            "approved_by": None,
            "approved_at": None,
            "input": self.profile.get("input"),
            "format": self.profile.get("format"),
            "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
            "generated_by": "measure_it.agent.mapper",
            "fuzzy_threshold": self.threshold,
            "llm": llm_info,
            "summary": summary,
            "tables": tables,
            "subgroup_analysis": None,
        }


# ------------------------------------------------------------------------------------------------------------------
# public API
# ------------------------------------------------------------------------------------------------------------------
def dataset_id_from_path(path: str | Path) -> str:
    p = Path(path)
    stem = p.name if p.is_dir() else p.stem
    s = slug(stem) or "dataset"
    if not s[0].isalpha():
        s = "ds_" + s
    s = s[:41]
    return s if len(s) >= 3 else (s + "_ds")[:41]


def propose_mapping(profile: dict, dataset: Dataset | None = None, *, dataset_id: str | None = None,
                    condition_hint: str | None = None, measurement_hint: str | None = None,
                    threshold: float = DEFAULT_THRESHOLD, llm=None, out: str | Path | None = None,
                    vocab: Vocab | None = None, echo=print) -> dict:
    """Profile (+ the local dataset, for complete value maps) -> mapping dict (contract section 2)."""
    m = Mapper(profile, dataset, condition_hint=condition_hint, measurement_hint=measurement_hint,
               threshold=threshold, vocab=vocab)
    return m.propose(dataset_id=dataset_id or dataset_id_from_path(profile.get("input", "dataset")), llm=llm, out=out,
                     echo=echo)


HEADER = """# mapping.yaml - PROPOSED by measure_it.agent.mapper (docs/AGENT_CONTRACT.md section 2).
# Review every column, edit freely, then approve (set approved_by / approved_at, or run apply with --approve).
# Roles: participant_id, label, age, age_band, sex, demographic, self_report_history, diagnosis_code, lab,
#        medication_code, medication_flag, survey, device, omics, date, drop, unmapped.
# value_map keys are source values as strings (integral numbers without '.0'); null = not in this label / unknown.
# proposed_by: method (rule | dictionary | fuzzy | llm | none), confidence 0-1, rationale. Nothing here was invented:
# codes come from configs/ vocabularies, the data, or a model proposal that passed local validation.
"""


def write_mapping(mapping: dict, out: str | Path) -> Path:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    p = out / "mapping.yaml"
    p.write_text(HEADER + yaml.safe_dump(mapping, sort_keys=False, allow_unicode=True, width=120,
                                         default_flow_style=None))
    return p


def _private_dir(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    gi = out / ".gitignore"
    if not gi.exists():
        gi.write_text("*\n")


def profile_cmd(path: str | Path, out: str | Path, echo=print) -> dict:
    """Read `path`, write `<out>/profile.json`; return the profile."""
    out = Path(out)
    _private_dir(out)
    ds = read_input(path)
    prof = profile_dataset(ds)
    p = write_profile(prof, out)
    echo(f"profile: {len(prof['tables'])} table(s), format {prof['format']} -> {p}")
    return prof


def map_cmd(path: str | Path, out: str | Path, llm=None, condition_hint: str | None = None,
            measurement_hint: str | None = None, *, dataset_id: str | None = None,
            threshold: float = DEFAULT_THRESHOLD, echo=print) -> dict:
    """Read + profile + propose; write `<out>/profile.json` and `<out>/mapping.yaml` (and `<out>/llm_requests.jsonl`
    when `llm="claude"` ran). Returns {"profile", "mapping", "profile_path", "mapping_path"}."""
    out = Path(out)
    _private_dir(out)
    ds = read_input(path)
    prof = profile_dataset(ds)
    pp = write_profile(prof, out)
    mapping = propose_mapping(prof, ds, dataset_id=dataset_id or dataset_id_from_path(path),
                              condition_hint=condition_hint, measurement_hint=measurement_hint, threshold=threshold,
                              llm=llm, out=out, echo=echo)
    mp = write_mapping(mapping, out)
    s = mapping["summary"]
    echo(f"mapping: {s['n_columns']} column(s); by method {s['by_method']} -> {mp}")
    return {"profile": prof, "mapping": mapping, "profile_path": pp, "mapping_path": mp}
