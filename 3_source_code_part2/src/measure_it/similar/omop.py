"""Query packs that apply a computable phenotype inside controlled-access OMOP CDM cohorts, run by the user, plus a
plain clinic export.

Nothing here touches person-level data: the generated SQL is run by an authorised user inside the All of Us
Researcher Workbench (BigQuery) or the N3C Enclave (Spark SQL), and returns AGGREGATE counts only, with the
small-cell rules of each programme enforced in the query (All of Us: no count of 1-20 may be published; N3C: no
count below 20 may be published; both: complementary suppression so a suppressed cell cannot be recovered by
subtraction).

Code mapping without an Athena download: source codes are resolved through the CDM's own `concept` table
(vocabulary_id + concept_code). ICD-10-CM codes match source concepts (condition_source_concept_id) and their standard
targets through `concept_relationship` 'Maps to' (condition_concept_id); LOINC codes are standard measurement
concepts; RxNorm ingredients expand to every descendant drug concept through `concept_ancestor`.

Features with no OMOP equivalent (PRO instruments, ring-wearable summaries, omics) are listed as NOT EVALUABLE and
held at their centre, folded into a reduced intercept (a reduced phenotype; the rule of
measure_it.harmonize.phenotype.score_phenotype); the coverage (share of |coefficient| mass evaluable) is written into
every file. All of Us Fitbit tables give an OPTIONAL sensitivity block for wearable features: Fitbit and the device
the phenotype was fitted on are different devices with different algorithms; values are not interchangeable.

Outputs (local-only, results/byod/<id>/similar/):
  aou_bigquery.sql, aou_fitbit_sensitivity_bigquery.sql (when wearable features exist), n3c_spark.sql,
  clinic_export/{icd10cm,loinc,rxnorm,demographics,not_evaluable_in_ehr}.csv + README.md, query_pack.json
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .phenotype_io import base_population, coverage, coverage_flag, is_nmr, reduced_intercept, threshold_lp

AOU_SUPPRESS_MAX = 20     # All of Us Data and Statistics Dissemination Policy: no participant counts of 1-20
N3C_SUPPRESS_MAX = 19     # N3C: counts below 20 must not be published
FEMALE_CONCEPT_N3C = 8532  # OMOP standard gender concept 'FEMALE'
AOU_STATE_OBS_CONCEPT = 1585249  # PPI 'StreetAddress_PIIState' (state of residence answer) in observation

UCUM_ALIASES = {"mg/dl": "mg/dL", "g/dl": "g/dL", "mmol/l": "mmol/L", "umol/l": "umol/L", "µmol/l": "umol/L",
                "nmol/l": "nmol/L", "u/l": "U/L", "iu/l": "[IU]/L", "ng/ml": "ng/mL", "pg/ml": "pg/mL",
                "ug/dl": "ug/dL", "mg/l": "mg/L", "%": "%", "fl": "fL", "1000 cells/ul": "10*3/uL",
                "10^3/ul": "10*3/uL", "10*3/ul": "10*3/uL", "mmhg": "mm[Hg]", "kg/m2": "kg/m2", "bpm": "/min",
                "beats/min": "/min", "ml/min/1.73m2": "mL/min/{1.73_m2}", "mol/l": "mol/L"}

# Wearable summaries All of Us Fitbit tables can approximate (column names per the AoU "Data Types and
# Organization" article and the CDR data dictionary; VERIFY against the dictionary of the CDR you query).
FITBIT_PROXIES = [
    (re.compile(r"step", re.I), "daily_steps",
     "SELECT person_id, AVG(steps) AS v FROM {p}activity_summary WHERE steps > 0 GROUP BY person_id",
     "mean daily Fitbit steps (activity_summary.steps)"),
    (re.compile(r"sleep.*(dur|total|min|time)|total_sleep|tst", re.I), "sleep_minutes",
     "SELECT person_id, AVG(minute_asleep) AS v FROM {p}sleep_daily_summary WHERE LOWER(CAST(is_main_sleep AS STRING)) IN ('true', '1') GROUP BY person_id",
     "mean main-sleep minutes asleep (sleep_daily_summary.minute_asleep)"),
    (re.compile(r"(resting|rest|night|sleep).*(hr|heart)|rhr|lowest_hr|min_hr", re.I), "resting_hr",
     "SELECT person_id, APPROX_QUANTILES(day_min, 2)[OFFSET(1)] AS v FROM (SELECT person_id, DATE(datetime) AS d, "
     "MIN(heart_rate_value) AS day_min FROM {p}heart_rate_minute_level GROUP BY person_id, DATE(datetime)) "
     "GROUP BY person_id",
     "median over days of the daily minimum minute heart rate (heart_rate_minute_level), a resting-HR proxy"),
]
FITBIT_NOT_AVAILABLE = re.compile(r"hrv|rmssd|sdnn|spo2|oxygen|temp|respir", re.I)


def _sql_str(s: str) -> str:
    return "'" + str(s).replace("'", "''") + "'"


def _ucum(unit: str | None) -> str | None:
    if not unit:
        return None
    return UCUM_ALIASES.get(unit.strip().lower(), unit.strip())


def _code_predicate(col: str, codes: list[str]) -> str:
    """Exact codes, plus prefix match for 3-character ICD-10-CM categories and codes ending in 'x'/'*'."""
    exact = [c for c in codes if not (len(c) == 3 or c.endswith(("x", "*")))]
    prefix = [c.rstrip("x*") for c in codes if c not in exact]
    parts = []
    if exact:
        parts.append(f"{col} IN ({', '.join(_sql_str(c) for c in exact)})")
    parts += [f"{col} LIKE {_sql_str(p + '%')}" for p in prefix]
    return "(" + " OR ".join(parts) + ")" if parts else "FALSE"


def icd_codes(f: dict) -> list[str]:
    """Codes to match for an ICD-10-CM feature: the 3-character category for code_match 'category_prefix' (the
    feature_matrix default), else the listed codes."""
    cat = f["feature_key"].split(":", 1)[1]
    if f.get("code_match", "category_prefix" if len(cat) == 3 else "exact") == "category_prefix" or not f.get("codes"):
        return [cat]
    return list(f["codes"])


class Dialect:
    def __init__(self, name: str, prefix: str):
        self.name, self.p = name, prefix

    def median(self, x: str) -> str:
        return f"APPROX_QUANTILES({x}, 2)[OFFSET(1)]" if self.name == "bigquery" else f"percentile_approx({x}, 0.5)"

    def current_year(self) -> str:
        return "EXTRACT(YEAR FROM CURRENT_DATE())" if self.name == "bigquery" else "year(current_date())"

    def ln(self, x: str) -> str:
        return f"LN({x})"


def classify_omop(f: dict) -> tuple[bool, str]:
    """(evaluable in an OMOP CDM, reason)."""
    v, t = f["vocabulary"], f["transform"]
    if v == "ICD10CM":
        return True, "condition_occurrence via ICD10CM source concepts + 'Maps to'"
    if v == "LOINC":
        note = "measurement (LOINC standard concept)"
        if is_nmr(f):
            note += "; NMR (Nightingale) values in the phenotype vs routine-lab values in the EHR: not calibrated"
        return t in ("standardized_value", "value", "log_standardized_value", "presence", "above_threshold",
                     "below_threshold"), note
    if v == "RXNORM":
        return True, "drug_exposure via RxNorm ingredient descendants (concept_ancestor)"
    if v == "DEMOG":
        key = f["feature_key"].split(":", 1)[1]
        return key in ("age_mid", "female", "male", "age_years"), "person table"
    if v == "SURVEY":
        return False, ("PRO instrument item/score: no OMOP concept; All of Us survey answers (observation, PPI "
                       "vocabulary) do not include this instrument")
    if v == "DEVICE":
        return False, "wearable/device summary: no OMOP concept (see the optional Fitbit sensitivity block for All of Us)"
    if v == "OMICS":
        return False, "omics feature: not in the OMOP clinical tables"
    return False, f"vocabulary {v} not mapped"


def _transform_sql(f: dict, x: str) -> str:
    """SQL expression of the transformed feature value given the raw per-person value column x (NULL = missing)."""
    t = f["transform"]
    if t == "presence":
        return f"CASE WHEN {x} IS NULL THEN 0 ELSE 1 END"
    if t == "standardized_value":
        return f"COALESCE(({x} - ({float(f['center'])})) / ({float(f['scale'])}), 0)"
    if t == "log_standardized_value":
        return f"COALESCE((LN(GREATEST({x}, 1e-9)) - ({float(f['center'])})) / ({float(f['scale'])}), 0)"
    if t == "value":
        return f"COALESCE({x}, 0)"
    if t == "above_threshold":
        return f"CASE WHEN {x} >= ({float(f['threshold'])}) THEN 1 ELSE 0 END"
    if t == "below_threshold":
        return f"CASE WHEN {x} <= ({float(f['threshold'])}) THEN 1 ELSE 0 END"
    raise ValueError(f"unsupported transform {t!r} for {f['feature_key']}")


def build_sql(ph: dict, *, dialect: str, prefix: str = "", fitbit: bool = False) -> tuple[str, dict]:
    """One query returning aggregate counts. Returns (sql, info)."""
    D = Dialect(dialect, prefix)
    p = prefix
    aou = dialect == "bigquery"
    bp = base_population(ph)
    lp_thr = threshold_lp(ph)
    evaluable, reasons, ctes, terms = {}, {}, [], []
    fitbit_used = []
    if not bp["icd10cm"]:
        raise ValueError(f"no base-population ICD-10-CM codes for condition {bp['condition_id']!r}")
    ctes.append(f"base_src AS (\n  SELECT concept_id FROM {p}concept\n  WHERE vocabulary_id = 'ICD10CM' AND "
                f"{_code_predicate('concept_code', bp['icd10cm'])}\n)")
    ctes.append(f"base_std AS (\n  SELECT cr.concept_id_2 AS concept_id FROM {p}concept_relationship cr\n"
                f"  JOIN base_src s ON cr.concept_id_1 = s.concept_id\n  WHERE cr.relationship_id = 'Maps to'\n)")
    ctes.append(f"base AS (\n  SELECT DISTINCT co.person_id FROM {p}condition_occurrence co\n"
                f"  WHERE co.condition_source_concept_id IN (SELECT concept_id FROM base_src)\n"
                f"     OR co.condition_concept_id IN (SELECT concept_id FROM base_std)\n)")
    sex_col = "sex_at_birth_concept_id" if aou else "gender_concept_id"
    female_expr = ("CASE WHEN LOWER(sx.concept_name) = 'female' THEN 1 ELSE 0 END" if aou
                   else f"CASE WHEN pe.gender_concept_id = {FEMALE_CONCEPT_N3C} THEN 1 ELSE 0 END")
    sex_join = f"\n  LEFT JOIN {p}concept sx ON sx.concept_id = pe.{sex_col}" if aou else ""
    ctes.append(f"cohort AS (\n  SELECT pe.person_id, {D.current_year()} - pe.year_of_birth AS age,\n"
                f"         {female_expr} AS female\n  FROM {p}person pe\n  JOIN base b ON b.person_id = pe.person_id"
                f"{sex_join}\n  WHERE {D.current_year()} - pe.year_of_birth >= 18\n)")
    joins = []
    for i, f in enumerate(ph["features"]):
        ok, why = classify_omop(f)
        reasons[f["feature_key"]] = why
        coef = float(f["coefficient"])
        name = f"f{i}"
        v = f["vocabulary"]
        if v == "DEVICE" and fitbit and aou and not FITBIT_NOT_AVAILABLE.search(f["feature_key"]):
            for rx, short, q, desc in FITBIT_PROXIES:
                if rx.search(f["feature_key"]) or rx.search(f.get("label", "")):
                    ctes.append(f"{name} AS (\n  {q.format(p=p)}\n)")
                    joins.append(f"LEFT JOIN {name} ON {name}.person_id = c.person_id")
                    terms.append(f"({coef}) * {_transform_sql(f, name + '.v')}")
                    fitbit_used.append({"feature_key": f["feature_key"], "fitbit_proxy": short, "definition": desc})
                    ok, reasons[f["feature_key"]] = True, f"Fitbit proxy ({desc}); different device, not calibrated"
                    break
            evaluable[f["feature_key"]] = ok
            continue
        evaluable[f["feature_key"]] = ok
        if not ok:
            continue
        if v == "ICD10CM":
            ctes.append(f"{name}_src AS (\n  SELECT concept_id FROM {p}concept WHERE vocabulary_id = 'ICD10CM' AND "
                        f"{_code_predicate('concept_code', icd_codes(f))}\n)")
            ctes.append(f"{name} AS (\n  SELECT DISTINCT co.person_id, 1 AS v FROM {p}condition_occurrence co\n"
                        f"  JOIN cohort c ON c.person_id = co.person_id\n"
                        f"  WHERE co.condition_source_concept_id IN (SELECT concept_id FROM {name}_src)\n"
                        f"     OR co.condition_concept_id IN (SELECT cr.concept_id_2 FROM {p}concept_relationship cr "
                        f"JOIN {name}_src s ON cr.concept_id_1 = s.concept_id WHERE cr.relationship_id = 'Maps to')\n)")
        elif v == "LOINC":
            unit = _ucum(f.get("unit"))
            unit_filter = (f"\n    AND m.unit_concept_id IN (SELECT concept_id FROM {p}concept WHERE vocabulary_id = "
                           f"'UCUM' AND concept_code = {_sql_str(unit)})") if unit and f["transform"] != "presence" else ""
            ctes.append(f"{name} AS (\n  SELECT m.person_id, {D.median('m.value_as_number')} AS v\n"
                        f"  FROM {p}measurement m\n  JOIN cohort c ON c.person_id = m.person_id\n"
                        f"  WHERE m.measurement_concept_id IN (SELECT concept_id FROM {p}concept WHERE vocabulary_id = "
                        f"'LOINC' AND concept_code IN ({', '.join(_sql_str(c) for c in f['codes'])}))\n"
                        f"    AND m.value_as_number IS NOT NULL{unit_filter}\n  GROUP BY m.person_id\n)")
        elif v == "RXNORM":
            ctes.append(f"{name} AS (\n  SELECT DISTINCT de.person_id, 1 AS v FROM {p}drug_exposure de\n"
                        f"  JOIN cohort c ON c.person_id = de.person_id\n"
                        f"  WHERE de.drug_concept_id IN (\n    SELECT ca.descendant_concept_id FROM {p}concept_ancestor ca\n"
                        f"    JOIN {p}concept i ON i.concept_id = ca.ancestor_concept_id\n"
                        f"    WHERE i.vocabulary_id = 'RxNorm' AND i.concept_code IN "
                        f"({', '.join(_sql_str(c) for c in f['codes'])}))\n)")
        elif v == "DEMOG":
            key = f["feature_key"].split(":", 1)[1]
            col = {"age_mid": "LEAST(c.age, 90)", "age_years": "LEAST(c.age, 90)", "female": "c.female",
                   "male": "(1 - c.female)"}[key]
            if key == "age_mid":   # age-band midpoint, 10-year bands, top band 90+ (harmonize.schema.age_mid)
                col = "CASE WHEN c.age >= 90 THEN 92.5 ELSE FLOOR(c.age / 10) * 10 + 4.5 END"
            expr = col if f["transform"] == "presence" else _transform_sql(f, col)
            terms.append(f"({coef}) * {expr}")
            continue
        joins.append(f"LEFT JOIN {name} ON {name}.person_id = c.person_id")
        terms.append(f"({coef}) * {_transform_sql(f, name + '.v')}")
    icpt = reduced_intercept(ph, evaluable)   # not-evaluable features held at their centre
    lp = " +\n         ".join([f"({icpt})"] + terms)
    ctes.append(f"scored AS (\n  SELECT c.person_id,\n         {lp} AS lp\n  FROM cohort c\n  "
                + "\n  ".join(joins) + "\n)")
    supp = AOU_SUPPRESS_MAX if aou else N3C_SUPPRESS_MAX
    if aou:
        ctes.append(f"state AS (\n  SELECT o.person_id, MAX(COALESCE(v.concept_name, o.value_source_value)) AS state\n"
                    f"  FROM {p}observation o\n  LEFT JOIN {p}concept v ON v.concept_id = o.value_source_concept_id\n"
                    f"  WHERE o.observation_source_concept_id = {AOU_STATE_OBS_CONCEPT}\n  GROUP BY o.person_id\n)")
        state_join = "LEFT JOIN state st ON st.person_id = s.person_id"
        state_col = "COALESCE(st.state, 'unknown')"
    else:
        state_join = (f"JOIN {p}person pe ON pe.person_id = s.person_id\n  LEFT JOIN {p}location l "
                      f"ON l.location_id = pe.location_id")
        state_col = "COALESCE(l.state, 'unknown')"
    ctes.append(f"by_level AS (\n  SELECT 'all' AS level, 'all' AS stratum, COUNT(*) AS n_base,\n"
                f"         COALESCE(SUM(CASE WHEN s.lp >= {lp_thr:.6f} THEN 1 ELSE 0 END), 0) AS n_positive\n  FROM scored s\n"
                f"  UNION ALL\n  SELECT 'state' AS level, {state_col} AS stratum, COUNT(*) AS n_base,\n"
                f"         COALESCE(SUM(CASE WHEN s.lp >= {lp_thr:.6f} THEN 1 ELSE 0 END), 0) AS n_positive\n  FROM scored s\n"
                f"  {state_join}\n  GROUP BY {state_col}\n)")
    # primary suppression (1..supp) on each count and on the phenotype-negative difference n_base - n_positive (it would
    # otherwise be recoverable by subtraction; found by the first local execution, agent.omop_run); complementary: if exactly one state cell of a column is
    # suppressed, also suppress the smallest unsuppressed state cell so it cannot be recovered from the total
    ctes.append(f"flagged AS (\n  SELECT level, stratum, n_base, n_positive,\n"
                f"    CASE WHEN n_base BETWEEN 1 AND {supp} THEN 1 ELSE 0 END AS sb,\n"
                f"    CASE WHEN n_positive BETWEEN 1 AND {supp} OR (n_base - n_positive) BETWEEN 1 AND {supp} "
                f"THEN 1 ELSE 0 END AS sp\n  FROM by_level\n)")
    ctes.append("ranked AS (\n  SELECT f.*,\n"
                "    SUM(sb) OVER (PARTITION BY level) AS nsb, SUM(sp) OVER (PARTITION BY level) AS nsp,\n"
                "    ROW_NUMBER() OVER (PARTITION BY level, sb ORDER BY n_base) AS rb,\n"
                "    ROW_NUMBER() OVER (PARTITION BY level, sp ORDER BY CASE WHEN n_positive = 0 THEN 1 ELSE 0 END, n_positive) AS rp\n  FROM flagged f\n)")
    sql = ("-- " + "\n-- ".join(_header(ph, bp, dialect, supp)) + "\nWITH\n" + ",\n".join(ctes) + "\n"
           "SELECT level, stratum,\n"
           "  CASE WHEN sb = 1 OR (level = 'state' AND nsb = 1 AND sb = 0 AND rb = 1) THEN NULL ELSE n_base END "
           "AS n_base,\n"
           "  CASE WHEN sp = 1 OR (level = 'state' AND nsp = 1 AND sp = 0 AND rp = 1) THEN NULL ELSE n_positive END "
           "AS n_phenotype_positive,\n"
           f"  CASE WHEN sb = 1 OR sp = 1 THEN 'suppressed: count 1-{supp}' "
           "WHEN level = 'state' AND ((nsb = 1 AND sb = 0 AND rb = 1) OR (nsp = 1 AND sp = 0 AND rp = 1)) "
           "THEN 'suppressed: complementary' ELSE '' END AS suppression\n"
           "FROM ranked\nORDER BY level, stratum")
    cov = coverage(ph, evaluable)
    flag = coverage_flag(cov["coverage"])
    if flag != "adequate":
        sql = (f"-- WARNING: feature coverage {cov['coverage']:.2f} ({flag}). Counts from this query describe a\n"
               "-- different, reduced rule; do not report them as the study phenotype.\n" + sql)
    cov["flag"] = flag
    return sql + "\n", {"evaluable": evaluable, "reasons": reasons, "coverage": cov, "base_population": bp,
                        "fitbit": fitbit_used, "threshold_lp": lp_thr}


def _header(ph: dict, bp: dict, dialect: str, supp: int) -> list[str]:
    target = ("All of Us Researcher Workbench (BigQuery; set the dataset prefix to your CDR, e.g. the "
              "WORKSPACE_CDR environment variable)" if dialect == "bigquery" else "N3C Enclave (Spark SQL)")
    return [f"measure-it similar: computable phenotype {ph['phenotype_id']} (version {ph.get('version', 1)})",
            f"Target: {target}.",
            f"Base population: {bp['label']}; adults (>= 18 by year of birth).",
            "Output: AGGREGATE counts only (overall and by state of residence).",
            f"Small cells: counts 1-{supp} are returned as NULL, with complementary suppression across states. "
            + ("All of Us policy forbids publishing participant counts of 1-20 (and derived values that reveal them)."
               if dialect == "bigquery" else "N3C forbids publishing counts below 20."),
            "Score: linear predictor = intercept + sum(coefficient x transformed feature); positive if lp >= "
            "logit(decision_threshold). Missing labs -> centre (contribution 0); absent codes / drugs -> 0.",
            "Features without an OMOP equivalent are held at their centre and folded into the intercept (reduced "
            "phenotype); see query_pack.json for coverage.",
            "Generated SQL was parsed (sqlglot) but has NOT been executed against a CDM; review before running."]


def clinic_export(ph: dict, out_dir: Path) -> dict:
    """CSV code lists + a one-page README a health-system report writer can implement in the clinic's EHR."""
    out_dir.mkdir(parents=True, exist_ok=True)
    bp = base_population(ph)
    rows = {"icd10cm": [], "loinc": [], "rxnorm": [], "demographics": [], "not_evaluable_in_ehr": []}
    evaluable = {}
    for f in ph["features"]:
        v, coef = f["vocabulary"], float(f["coefficient"])
        common = {"feature_key": f["feature_key"], "label": f.get("label", ""), "coefficient": coef,
                  "transform": f["transform"]}
        if v == "ICD10CM":
            for c in icd_codes(f):
                rows["icd10cm"].append({**common, "code": c, "match": "prefix" if len(c) == 3 else "exact",
                                        "codes_seen_in_study": " ".join(f.get("codes", [])),
                                        "note": f.get("ehr_note", "")})
            evaluable[f["feature_key"]] = True
        elif v == "LOINC":
            for c in f["codes"]:
                rows["loinc"].append({**common, "loinc": c, "unit": f.get("unit", ""), "center": f.get("center"),
                                      "scale": f.get("scale"), "threshold": f.get("threshold"),
                                      "aggregation": "median of the patient's results in a stated look-back window (e.g. 24 months)",
                                      "platform_caveat": ("phenotype fitted on NMR (Nightingale) values; routine "
                                                          "lab values are not calibrated to them") if is_nmr(f) else ""})
            evaluable[f["feature_key"]] = True
        elif v == "RXNORM":
            for c in f["codes"]:
                rows["rxnorm"].append({**common, "rxnorm_ingredient_cui": c,
                                       "match": "any product containing this ingredient"})
            evaluable[f["feature_key"]] = True
        elif v == "DEMOG":
            rows["demographics"].append({**common, "center": f.get("center"), "scale": f.get("scale"),
                                         "definition": {"DEMOG:age_mid": "midpoint of the 10-year age band "
                                                        "(e.g. 40-49 -> 44.5; 90+ -> 92.5)",
                                                        "DEMOG:female": "1 if sex is female else 0"}.get(
                                             f["feature_key"], "")})
            evaluable[f["feature_key"]] = True
        else:
            why = {"SURVEY": "patient-reported instrument; evaluable only if the clinic administers this exact "
                             "instrument and stores the score",
                   "DEVICE": "device/wearable summary; not in the EHR",
                   "OMICS": "research omics feature; not in the EHR"}.get(v, "not an EHR vocabulary")
            xw = "; ".join(f"{r.get('reference')} ({r.get('match_quality')})" for r in f.get("reference_crosswalk") or []
                           if r.get("reference"))
            rows["not_evaluable_in_ehr"].append({**common, "domain": f["domain"], "reason": why,
                                                 "center": f.get("center"), "scale": f.get("scale"),
                                                 "reference_crosswalk": xw})
            evaluable[f["feature_key"]] = False
    for k, r in rows.items():
        pd.DataFrame(r).to_csv(out_dir / f"{k}.csv", index=False)
    cov = coverage(ph, evaluable)
    lp_thr = threshold_lp(ph)
    readme = _clinic_readme(ph, bp, cov, lp_thr, rows, reduced_intercept(ph, evaluable))
    (out_dir / "README.md").write_text(readme)
    return {"coverage": cov, "files": sorted(p.name for p in out_dir.iterdir())}


def _fmt(x) -> str:
    return "" if x is None or (isinstance(x, float) and pd.isna(x)) else f"{x:g}" if isinstance(x, float) else str(x)


def _clinic_readme(ph, bp, cov, lp_thr, rows, icpt) -> str:
    dom = ", ".join(f"{d} {v['evaluable_share_of_domain']:.0%}" for d, v in sorted(cov["by_domain"].items()))
    L = [f"# Computable phenotype for a clinic EHR: `{ph['phenotype_id']}`", "",
         "A candidate profile of patients who resemble the device-positive subgroup of a research study. It is a "
         "research-readiness aid for finding people who might benefit from the measurement, not a diagnosis, and it "
         "does not identify anyone outside your own system. Run it inside your EHR reporting tool; publish only "
         "aggregate counts and follow your institution's small-cell rules.", "",
         "## 1. Base population", "",
         f"Adults (18+) with at least one diagnosis code in: {', '.join(bp['icd10cm'])} ({bp['label']}). Codes of 3 "
         "characters (or ending in x) match every code that starts with them.", "",
         "## 2. Features (code lists in this folder)", "",
         f"* `icd10cm.csv`: {len(rows['icd10cm'])} code rows. Feature = 1 if any listed code is on the record, else 0.",
         f"* `loinc.csv`: {len(rows['loinc'])} lab rows. Use results in the listed unit only (do not convert other "
         "units silently); take the median of the patient's results; transform as in the column `transform`.",
         f"* `rxnorm.csv`: {len(rows['rxnorm'])} ingredient rows. Feature = 1 if any order/administration contains "
         "the ingredient.",
         f"* `demographics.csv`: {len(rows['demographics'])} rows.",
         f"* `not_evaluable_in_ehr.csv`: {len(rows['not_evaluable_in_ehr'])} features your EHR cannot evaluate "
         "(PRO instruments, device summaries, omics). They are held at their study centre, which is already folded "
         "into the intercept below. If your clinic administers the same instrument, add it back with its "
         "centre/scale.", "",
         "## 3. Scoring rule", "",
         "```",
         f"lp = {icpt:g}   (model intercept {float(ph['intercept']):g} + not-evaluable features at their centre)",
         "     + sum over features of coefficient x value",
         "value = 1/0                              for transform 'presence'",
         "value = (median_result - center) / scale for transform 'standardized_value' (0 when no result)",
         f"candidate = lp >= {lp_thr:.4f}            (= probability >= {float(ph['decision_threshold']):g})",
         "```", "",
         "## 4. What your EHR can evaluate", "",
         f"Coverage = share of the model's |standardized coefficient| mass your EHR can evaluate: **{cov['coverage']:.0%}** "
         f"({cov['n_evaluable']} of {cov['n_features']} features). By domain: {dom}. Below 100% the rule is a "
         "reduced phenotype: expect fewer candidates and a different (unknown) accuracy than in the study.", "",
         "## 5. Caveats", "",
         f"* Study performance (in the research cohort, all features): {json.dumps(ph.get('performance', {}))}.",
         "* Lab values fitted on NMR metabolomics (rows with a `platform_caveat`) are not interchangeable with routine "
         "clinical chemistry without a calibration study.",
         "* Self-reported history in the study was mapped to ICD-10-CM; coded diagnoses in an EHR are recorded "
         "differently (under-coding of post-infectious conditions is common).",
         ] + [f"* {c}" for c in ph.get("caveats", [])] + [""]
    return "\n".join(L)


def write_query_pack(ph: dict, out_dir: Path) -> dict:
    """Write every query-pack file for one phenotype under out_dir; returns the manifest written as query_pack.json."""
    out_dir.mkdir(parents=True, exist_ok=True)
    sid = re.sub(r"[^A-Za-z0-9_.-]+", "_", ph["phenotype_id"].split(":", 1)[-1])
    manifest = {"phenotype_id": ph["phenotype_id"], "files": {}, "coverage": {}, "not_evaluable": {}}
    aou_sql, aou = build_sql(ph, dialect="bigquery", prefix="`{CDR}`.")
    (out_dir / f"aou_bigquery__{sid}.sql").write_text(aou_sql)
    manifest["files"]["all_of_us"] = f"aou_bigquery__{sid}.sql"
    manifest["coverage"]["all_of_us"] = aou["coverage"]
    manifest["not_evaluable"]["omop"] = {k: aou["reasons"][k] for k, ok in aou["evaluable"].items() if not ok}
    manifest["base_population"] = aou["base_population"]
    if any(f["vocabulary"] == "DEVICE" for f in ph["features"]):
        fb_sql, fb = build_sql(ph, dialect="bigquery", prefix="`{CDR}`.", fitbit=True)
        if fb["fitbit"]:
            (out_dir / f"aou_fitbit_sensitivity_bigquery__{sid}.sql").write_text(
                "-- SENSITIVITY ONLY: wearable features approximated from All of Us Fitbit tables. Fitbit and the\n"
                "-- study device are different devices (sensors, algorithms, wear); centre/scale were fitted on the\n"
                "-- study device, so the Fitbit-based contributions are uncalibrated. Report beside, never instead of,\n"
                "-- the main query.\n" + fb_sql)
            manifest["files"]["all_of_us_fitbit_sensitivity"] = f"aou_fitbit_sensitivity_bigquery__{sid}.sql"
            manifest["coverage"]["all_of_us_with_fitbit_proxies"] = fb["coverage"]
            manifest["fitbit_proxies"] = fb["fitbit"]
    n3c_sql, n3c = build_sql(ph, dialect="spark", prefix="")
    (out_dir / f"n3c_spark__{sid}.sql").write_text(n3c_sql)
    manifest["files"]["n3c"] = f"n3c_spark__{sid}.sql"
    manifest["coverage"]["n3c"] = n3c["coverage"]
    ce = clinic_export(ph, out_dir / f"clinic_export__{sid}")
    manifest["files"]["clinic_export"] = f"clinic_export__{sid}/"
    manifest["coverage"]["clinic_ehr"] = ce["coverage"]
    manifest["suppression"] = {"all_of_us": f"counts 1-{AOU_SUPPRESS_MAX} -> NULL (+ complementary)",
                               "n3c": f"counts 1-{N3C_SUPPRESS_MAX} -> NULL (+ complementary)"}
    manifest["validated"] = "parsed with sqlglot (bigquery / spark dialects); not executed against a CDM"
    (out_dir / f"query_pack__{sid}.json").write_text(json.dumps(manifest, indent=1, default=str))
    return manifest


def parse_check(sql: str, dialect: str) -> None:
    """Raise if sqlglot cannot parse the SQL in the given dialect (sqlglot is optional: uv run --with sqlglot)."""
    import sqlglot
    sqlglot.parse_one(sql.replace("`{CDR}`.", "`proj.cdr`."), read=dialect, error_level=sqlglot.ErrorLevel.RAISE)
