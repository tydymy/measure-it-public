"""All of Us public Data Browser: participant counts for the lab concepts behind the published lab biomarkers.

    uv run python -m measure_it.labs.aou_lab_concepts            # query (cached), save raw JSON, print the table
    uv run python -m measure_it.labs.aou_lab_concepts --append   # also append new rows to controlled_cohort_counts.csv

What it does
    POSTs the Data Browser's own search endpoint (the one its web page calls)
        https://public.api.researchallofus.org/v1/databrowser/searchConcepts
        body {"query": <term>, "domain": MEASUREMENT | PROCEDURE, "standardConceptFilter": "STANDARD_OR_CODE_ID_MATCH",
              "maxResults": 25, "minCount": 1}
    through measure_it.http (disk-cached), keeps the concepts listed in LAB_CONCEPTS (matched by vocabulary + code,
    never by name similarity), and records the terms that returned no concept at all (e.g. CGRP).

Reading the counts
    * `countValue` = participants with at least one record of the concept (standard concept, including records mapped
      to it from other source codes); `sourceCountValue` = participants whose record carried this code as the source
      code. Both are rounded to multiples of 20 by the Data Browser; the smallest value it shows is 20, which this
      module reports as "<=20" (a count of 1-20 is never shown, in line with the All of Us rule that no count of 1-20
      may be disseminated).
    * These are counts of people who HAD THE TEST in their EHR, not people with an abnormal result, and not people
      with a diagnosis. Who gets tested is decided by a clinician who already suspects the condition (ascertainment
      bias; see docs/DATA_ACCESS_PLAN.md "Lab biomarkers in All of Us").
    * The Data Browser does not publish intersections (test x diagnosis); those need the Registered Tier.

Outputs
    data/raw/controlled_cohorts/aou_databrowser/searchConcepts_labs.json   raw responses per query term
    results/tables/controlled_cohort_counts.csv                            appended rows (with --append; deduplicated)
"""
from __future__ import annotations

import argparse
import csv
import json
import re

from ..config import RAW, TABLES, UNKNOWN, utc_now_iso

SEARCH_URL = "https://public.api.researchallofus.org/v1/databrowser/searchConcepts"
CDR_URL = "https://public.api.researchallofus.org/v1/databrowser/cdrversion-used"
ANALYSIS_URL = "https://public.api.researchallofus.org/v1/databrowser/concept-analysis-results"
SEX_SPLIT_CONCEPTS = {"3019420": "serum tryptase (LOINC 21582-2)"}  # sex-at-birth split read from the analysis endpoint
RAW_PATH = RAW / "controlled_cohorts" / "aou_databrowser" / "searchConcepts_labs.json"
COUNTS_CSV = TABLES / "controlled_cohort_counts.csv"
MIN_SHOWN = 20  # the Data Browser's floor; 20 means 1-20

# (group, target biomarker, search term, domain, vocabulary, code). Codes were chosen by reading every concept each
# term returns (see the raw JSON) and keeping the ones that measure the target analyte.
LAB_CONCEPTS: list[tuple[str, str, str, str, str, str]] = [
    ("mast_cell", "serum tryptase", "tryptase", "MEASUREMENT", "LOINC", "21582-2"),
    ("mast_cell", "serum tryptase", "tryptase", "MEASUREMENT", "SNOMED", "121873004"),
    ("mast_cell", "serum tryptase", "tryptase", "MEASUREMENT", "LOINC", "51834-0"),
    ("mast_cell", "urinary leukotriene E4", "leukotriene", "MEASUREMENT", "LOINC", "33343-5"),
    ("mast_cell", "urinary leukotriene E4", "leukotriene", "MEASUREMENT", "LOINC", "33344-3"),
    ("mast_cell", "urinary leukotriene E4", "leukotriene", "MEASUREMENT", "LOINC", "101115-4"),
    ("mast_cell", "urinary 2,3-dinor-11beta-PGF2alpha", "prostaglandin", "MEASUREMENT", "LOINC", "94381-1"),
    ("mast_cell", "urinary 2,3-dinor-11beta-PGF2alpha", "prostaglandin", "MEASUREMENT", "LOINC", "97658-9"),
    ("mast_cell", "urinary prostaglandin D2", "prostaglandin", "MEASUREMENT", "LOINC", "14054-1"),
    ("mast_cell", "urinary N-methylhistamine", "methylhistamine", "MEASUREMENT", "LOINC", "44340-8"),
    ("mast_cell", "urinary N-methylhistamine", "methylhistamine", "MEASUREMENT", "LOINC", "13781-0"),
    ("mast_cell", "urinary N-methylhistamine", "methylhistamine", "MEASUREMENT", "LOINC", "26053-9"),
    ("mast_cell", "urinary N-methylhistamine", "methylhistamine", "MEASUREMENT", "LOINC", "12714-2"),
    ("mast_cell", "plasma/serum histamine", "histamine", "MEASUREMENT", "LOINC", "2416-6"),
    ("mast_cell", "KIT D816V (c.2447A>T)", "KIT gene", "MEASUREMENT", "LOINC", "88519-4"),
    ("mast_cell", "KIT targeted mutation analysis", "KIT gene", "MEASUREMENT", "LOINC", "55201-8"),
    ("mast_cell", "KIT gene analysis (CPT 81273, mastocytosis)", "KIT gene", "MEASUREMENT", "CPT4", "81273"),
    ("catecholamine", "plasma norepinephrine", "norepinephrine plasma", "MEASUREMENT", "LOINC", "2666-6"),
    ("catecholamine", "plasma norepinephrine", "norepinephrine plasma", "MEASUREMENT", "LOINC", "14852-8"),
    ("catecholamine", "plasma norepinephrine supine", "norepinephrine plasma", "MEASUREMENT", "LOINC", "1601-4"),
    ("catecholamine", "plasma norepinephrine standing", "norepinephrine plasma", "MEASUREMENT", "LOINC", "17368-2"),
    ("catecholamine", "plasma epinephrine standing", "standing", "MEASUREMENT", "LOINC", "95054-3"),
    ("catecholamine", "plasma catecholamines", "norepinephrine plasma", "MEASUREMENT", "LOINC", "2056-0"),
    ("catecholamine", "plasma catecholamines 3 panel", "norepinephrine plasma", "MEASUREMENT", "LOINC", "34551-2"),
    ("catecholamine", "plasma free fractionated catecholamines", "norepinephrine plasma", "MEASUREMENT", "LOINC",
     "42493-7"),
    ("catecholamine", "24-h urine norepinephrine", "norepinephrine", "MEASUREMENT", "LOINC", "2668-2"),
    ("cortisol", "serum/plasma cortisol", "cortisol", "MEASUREMENT", "LOINC", "2143-6"),
    ("cortisol", "serum/plasma cortisol AM", "cortisol", "MEASUREMENT", "LOINC", "9813-7"),
    ("cortisol", "24-h urine free cortisol", "cortisol", "MEASUREMENT", "LOINC", "2147-7"),
    ("cortisol", "salivary cortisol", "cortisol", "MEASUREMENT", "LOINC", "2142-8"),
    ("autonomic", "ganglionic AChR antibody", "ganglionic", "MEASUREMENT", "LOINC", "42233-7"),
    ("autonomic", "heart rate standing", "standing", "MEASUREMENT", "LOINC", "69001-6"),
    ("autonomic", "heart rate supine", "supine", "MEASUREMENT", "LOINC", "68999-2"),
    ("autonomic", "tilt table evaluation", "tilt table", "PROCEDURE", "CPT4", "93660"),
    ("autonomic", "autonomic function testing (combined)", "95924", "PROCEDURE", "CPT4", "95924"),
    ("blood_volume", "plasma volume, single sampling", "78110", "PROCEDURE", "CPT4", "78110"),
    ("blood_volume", "plasma volume, multiple sampling", "78111", "PROCEDURE", "CPT4", "78111"),
    ("blood_volume", "whole blood volume incl. plasma + red cell volume", "78122", "PROCEDURE", "CPT4", "78122"),
]

# Terms searched that should be recorded even when nothing matches (absence is a finding).
# A returned concept counts as relevant only if its name matches the regex (the search is fuzzy: 'adrenergic
# receptor' returns an ADRA1B germline-variant concept and a generic CPT molecular-pathology code, neither of which is
# an autoantibody assay).
ABSENCE_TERMS: list[tuple[str, str, str, str]] = [
    ("neuropeptide", "calcitonin gene-related peptide (CGRP)", "calcitonin gene", r"calcitonin gene|CGRP"),
    ("neuropeptide", "calcitonin gene-related peptide (CGRP)", "CGRP", r"calcitonin gene|CGRP"),
    ("autoantibody", "GPCR / adrenergic / muscarinic receptor autoantibodies", "muscarinic", r"\bAb\b|antibod"),
    ("autoantibody", "GPCR / adrenergic / muscarinic receptor autoantibodies", "adrenergic receptor",
     r"\bAb\b|antibod"),
    ("mast_cell", "TPSAB1 copy number (hereditary alpha-tryptasemia)", "TPSAB1", r"TPSAB1|alpha.?tryptase"),
]


def _body(term: str, domain: str) -> dict:
    return {"query": term, "domain": domain, "standardConceptFilter": "STANDARD_OR_CODE_ID_MATCH",
            "maxResults": 25, "minCount": 1}


def fetch_data() -> dict:
    """{'cdr': ..., 'retrieved_at': ..., 'responses': {'DOMAIN|term': response}} (cached HTTP)."""
    from ..http import get_json, post_json
    queries = sorted({(t, d) for _, _, t, d, _, _ in LAB_CONCEPTS} | {(t, "MEASUREMENT") for _, _, t, _ in ABSENCE_TERMS})
    out = {"endpoint": SEARCH_URL, "retrieved_at": utc_now_iso(), "cdr": get_json(CDR_URL), "responses": {}}
    for term, dom in queries:
        out["responses"][f"{dom}|{term}"] = {"body": _body(term, dom),
                                             "response": post_json(SEARCH_URL, json_body=_body(term, dom))}
    out["sex_analysis"] = {}
    for cid in SEX_SPLIT_CONCEPTS:
        d = get_json(ANALYSIS_URL, params={"concept-ids": cid, "domain-id": "Measurement"})
        item = (d.get("items") or [{}])[0]
        out["sex_analysis"][cid] = {"url": ANALYSIS_URL, "params": {"concept-ids": cid, "domain-id": "Measurement"},
                                    "genderAnalysis": item.get("genderAnalysis")}
    return out


def sex_rows(data: dict) -> list[dict]:
    rows = []
    for cid, label in SEX_SPLIT_CONCEPTS.items():
        ga = ((data.get("sex_analysis") or {}).get(cid) or {}).get("genderAnalysis") or {}
        for r in ga.get("results") or []:
            rows.append({"concept_id": cid, "label": label, "sex": r.get("analysisStratumName"),
                         "participants": r.get("countValue")})
    return rows


def save_raw(data: dict) -> None:
    RAW_PATH.parent.mkdir(parents=True, exist_ok=True)
    RAW_PATH.write_text(json.dumps(data, indent=1, sort_keys=True), encoding="utf-8")


def load_raw() -> dict | None:
    return json.loads(RAW_PATH.read_text(encoding="utf-8")) if RAW_PATH.exists() else None


def shown(v) -> str:
    """Data Browser value as displayed: 20 is the floor (1-20)."""
    if v is None:
        return UNKNOWN
    if int(v) == 0:
        return "0"
    return f"<={MIN_SHOWN}" if int(v) <= MIN_SHOWN else str(int(v))


def concept_rows(data: dict) -> list[dict]:
    rows = []
    for group, target, term, dom, vocab, code in LAB_CONCEPTS:
        items = (data["responses"].get(f"{dom}|{term}") or {}).get("response", {}).get("items", [])
        hit = next((it for it in items if it["vocabularyId"] == vocab and it["conceptCode"] == code), None)
        rows.append({
            "group": group, "target_biomarker": target, "search_term": term, "domain": dom, "vocabulary": vocab,
            "code": code, "concept_id": hit["conceptId"] if hit else UNKNOWN,
            "concept_name": hit["conceptName"] if hit else UNKNOWN,
            "participants": hit["countValue"] if hit else None,
            "participants_source_code": hit["sourceCountValue"] if hit else None,
            "participants_shown": shown(hit["countValue"]) if hit else UNKNOWN,
            "found": hit is not None,
        })
    return rows


def absence_rows(data: dict) -> list[dict]:
    rows = []
    for group, target, term, pattern in ABSENCE_TERMS:
        items = (data["responses"].get(f"MEASUREMENT|{term}") or {}).get("response", {}).get("items", [])
        relevant = [i for i in items if re.search(pattern, i["conceptName"], re.I)]
        rows.append({"group": group, "target_biomarker": target, "search_term": term, "n_concepts": len(items),
                     "n_relevant_concepts": len(relevant),
                     "concepts": [f"{i['vocabularyId']} {i['conceptCode']} {i['conceptName']}" for i in items]})
    return rows


def cohort_count_rows(data: dict, retrieved_date: str) -> list[dict]:
    """Rows in the controlled_cohort_counts.csv schema (cohort, quantity, value, source_url, retrieved_date, status)."""
    src = (f"{SEARCH_URL} (POST; raw: data/raw/controlled_cohorts/aou_databrowser/searchConcepts_labs.json)")
    out = []
    for r in concept_rows(data):
        if not r["found"]:
            continue
        v = r["participants"]
        value = str(v) if v > MIN_SHOWN else f"<={MIN_SHOWN}"
        out.append({"cohort": "All of Us",
                    "quantity": (f"participants with lab/measurement concept {r['vocabulary']} {r['code']} "
                                 f"({r['concept_name']}) [{r['target_biomarker']}] (rounded to 20; tested, not "
                                 f"abnormal)"),
                    "value": value, "source_url": src, "retrieved_date": retrieved_date, "status": "verified"})
    for r in sex_rows(data):
        out.append({"cohort": "All of Us",
                    "quantity": f"participants with {r['label']}, sex at birth = {r['sex']} (rounded to 20)",
                    "value": str(r["participants"]) if r["participants"] > MIN_SHOWN else f"<={MIN_SHOWN}",
                    "source_url": f"{ANALYSIS_URL}?concept-ids={r['concept_id']}&domain-id=Measurement",
                    "retrieved_date": retrieved_date, "status": "verified"})
    for a in absence_rows(data):
        out.append({"cohort": "All of Us",
                    "quantity": (f"Labs & Measurements concepts measuring {a['target_biomarker']} returned for "
                                 f"search '{a['search_term']}' (number of CONCEPTS, not participants; "
                                 f"{a['n_concepts']} concepts returned in total)"),
                    "value": str(a["n_relevant_concepts"]), "source_url": src, "retrieved_date": retrieved_date,
                    "status": "verified"})
    return out


def append_counts(rows: list[dict]) -> int:
    """Append rows whose (cohort, quantity) is not already in the CSV; keeps its columns. Returns rows added."""
    with open(COUNTS_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames
        have = {(r["cohort"], r["quantity"]) for r in reader}
    new = [r for r in rows if (r["cohort"], r["quantity"]) not in have]
    if new:
        with open(COUNTS_CSV, "a", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            for r in new:
                w.writerow({c: r.get(c, "") for c in cols})
    return len(new)


def run(fetch: bool = True, append: bool = True) -> dict:
    """Pipeline entry point (step `labs_aou_lab_concepts`). With fetch=False (what the pipeline passes offline) it reads
    the saved raw JSON only and never touches the network; with fetch=True it re-queries through the disk-cached HTTP
    layer and re-saves the raw JSON. New rows are appended to controlled_cohort_counts.csv only when their
    (cohort, quantity) is not there yet, so re-runs are idempotent."""
    data = None if fetch else load_raw()
    if data is None:
        if not fetch:
            raise FileNotFoundError(f"{RAW_PATH} missing: run `python -m measure_it.labs.aou_lab_concepts` online once")
        data = fetch_data()
        save_raw(data)
    rows = concept_rows(data)
    added = append_counts(cohort_count_rows(data, data["retrieved_at"][:10])) if append else 0
    out = {"cdr": data.get("cdr", {}), "concepts": len(rows), "concepts_found": sum(r["found"] for r in rows),
           "absence_checks": len(absence_rows(data)), "rows_appended": added, "raw": str(RAW_PATH)}
    print(json.dumps(out, default=str)[:800])
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--append", action="store_true", help="append rows to results/tables/controlled_cohort_counts.csv")
    ap.add_argument("--no-fetch", action="store_true", help="use the saved raw JSON only")
    a = ap.parse_args(argv)
    data = load_raw() if a.no_fetch else None
    if data is None:
        data = fetch_data()
        save_raw(data)
    cdr = data.get("cdr", {})
    print(f"CDR: {json.dumps(cdr)[:200]}")
    for r in concept_rows(data):
        print(f"{r['group']:13} {r['vocabulary']:6} {r['code']:>10} {r['participants_shown']:>7} "
              f"(source {shown(r['participants_source_code']) if r['found'] else '-'}) {r['concept_name']}")
    for r in sex_rows(data):
        print(f"SEX {r['label']}: {r['sex']} {shown(r['participants'])}")
    for x in absence_rows(data):
        print(f"ABSENCE CHECK {x['search_term']!r}: {x['n_relevant_concepts']} relevant of {x['n_concepts']} "
              f"concepts {x['concepts'][:3]}")
    if a.append:
        n = append_counts(cohort_count_rows(data, data["retrieved_at"][:10]))
        print(f"appended {n} rows to {COUNTS_CSV}")


if __name__ == "__main__":
    main()
