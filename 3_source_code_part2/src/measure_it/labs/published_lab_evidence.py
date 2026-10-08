"""Published lab / biomarker evidence for the target invisible illnesses (curated literature table).

Every row is a PUBLISHED CLAIM, NOT REPRODUCED BY THIS PROJECT. The table answers "what does the peer-reviewed
literature report about how well lab test X separates condition Y from controls, or how often it is abnormal?" and
keeps the design, the comparator (healthy vs symptomatic look-alike), the reference standard, whether the biomarker
is itself part of the diagnostic criteria (incorporation bias) and the risk-of-bias notes next to every number.
It is the lab-test companion of measure_it.measurements.published_evidence (device evidence) and reuses that
module's quote verifier and number tracer unchanged.

    uv run python -m measure_it.labs.published_lab_evidence              # verify quotes, write the table
    uv run python -m measure_it.labs.published_lab_evidence --no-fetch   # verify against stored texts only

Inputs
    data/raw/published_lab_evidence/published_lab_evidence.csv   the curated table (one row per claim)
    Europe PMC REST (search, resultType=core; <PMCID>/fullTextXML) abstracts, metadata, open-access full text

Anti-fabrication contract (enforced here and in tests/test_published_lab_evidence.py)
    * every numeric claim sits inside `source_quote`, a verbatim extract (segments joined by ' || ') of the text
      named by `quote_source` (europepmc_abstract / europepmc_fulltext);
    * `build()` re-fetches those texts through measure_it.http (disk-cached) and checks every quote segment as a
      substring after light normalisation; a row whose quote is not found is kept but flagged quote_verified=False;
    * AUROC / sensitivity / specificity (and counts marked "as stated in source_quote") must appear in the quote;
      counts the curator summed are explained in `count_derivation`;
    * the Europe PMC DOI must equal the curated DOI (doi_matches_source).

Output
    data/processed/published_lab_evidence.parquet   data_layer 'condition_molecular', evidence_type
                                                    'published_biomarker', object_id 'published_lab_evidence:<id>'
    data/raw/published_lab_evidence/source_texts.jsonl

Query
    get_published_lab_evidence(condition=None, biomarker=None) -> dict
        {'status': 'matched', 'claim_label', 'n_rows', 'rows': [...]} or {'status': UNKNOWN / NOT AVAILABLE, 'reason'}
"""
from __future__ import annotations

import argparse
import json
import re
from functools import lru_cache

import pandas as pd

from ..config import RAW, UNKNOWN, load_config, raw_dir, utc_now_iso
from ..download import load_manifest, record_file
from ..measurements import published_evidence as pe
from ..provenance import add_provenance
from ..store import read_table, table_exists, write_table

SOURCE_ID = "published_lab_evidence"
TABLE = "published_lab_evidence"
CSV_NAME = "published_lab_evidence.csv"
TEXTS_NAME = "source_texts.jsonl"
CLAIM_LABEL = pe.CLAIM_LABEL  # "published claim — not reproduced by this project"
CURATED_ON = "2026-09-24"
EPMC_SEARCH, EPMC_FULLTEXT = pe.EPMC_SEARCH, pe.EPMC_FULLTEXT

REQUIRED_COLUMNS = [
    "evidence_id", "condition_ids", "condition_definition_in_study", "biomarker", "specimen", "assay_method",
    "analyte_class", "measurement_classes", "measurement_class_gap", "claim_type", "finding_direction",
    "incorporation_bias", "comparator_type", "study_design", "evidence_tier", "n_cases", "n_controls",
    "count_derivation", "control_group", "reference_standard", "population_setting", "metric", "metric_value",
    "ci_or_spread", "auroc", "sensitivity", "specificity", "external_validation", "risk_of_bias_notes",
    "source_quote", "quote_source", "pmid", "pmcid", "doi", "first_author", "year", "source_url_override",
]
CLAIM_TYPES = pe.CLAIM_TYPES | {
    "no_validated_biomarker",        # the source states that no validated lab biomarker exists
    "treatment_response_prediction", # the lab value predicts response to a treatment
    "prognostic_association",        # the lab value predicts a later outcome (e.g. developing PTLDS)
    "passive_transfer_experiment",   # patient material transfers the phenotype to animals (mechanistic)
}
FINDING_DIRECTIONS = pe.FINDING_DIRECTIONS
QUOTE_SOURCES = {"europepmc_abstract", "europepmc_fulltext"}
EVIDENCE_TIERS = pe.EVIDENCE_TIERS
EXTERNAL_VALIDATION = pe.EXTERNAL_VALIDATION
INCORPORATION_BIAS = {"yes", "partial", "no"}
COMPARATOR_TYPES = {"healthy", "disease_or_symptomatic", "mixed", "none"}
ANALYTE_CLASSES = {
    "mast_cell_mediator", "genetic_variant", "catecholamine", "autoantibody", "blood_volume", "hormone",
    "cytokine_chemokine", "immune_cell_function", "immune_cell_phenotype", "metabolomic_panel", "proteomic_panel",
    "viral_antigen", "coagulation_microclot", "histology", "microrna", "epigenetic", "serology_antibody",
    "neuropeptide", "clinical_chemistry_panel", "tumour_marker", "biophysical_cell_assay", "other",
}
# measurement classes (configs/measurements.yaml) that a lab row may carry
LAB_MEASUREMENT_CLASSES = {"blood_biomarkers", "metabolomics", "proteomics", "transcriptomics", "immune_assays",
                           "small_fiber_testing", "cpet"}  # cpet: lab sampled around a 2-day CPET provocation
DISCRIMINATION_CLAIMS = pe.DISCRIMINATION_CLAIMS
_split = pe._split


# --------------------------------------------------------------------------- curated table
def csv_path():
    return RAW / SOURCE_ID / CSV_NAME


def load_curated() -> pd.DataFrame:
    return pd.read_csv(csv_path(), dtype=str, keep_default_na=False, encoding="utf-8")


def validate_curated(df: pd.DataFrame) -> list[str]:
    """Structural checks on the curated table; returns a list of problems (empty = valid)."""
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        return [f"missing columns {missing}"]
    problems: list[str] = []
    if df["evidence_id"].duplicated().any():
        problems.append(f"duplicate evidence_id {sorted(df.loc[df.evidence_id.duplicated(), 'evidence_id'])}")
    conds = {c["id"] for c in load_config("conditions")["conditions"]}
    classes = {m["id"] for m in load_config("measurements")["measurement_classes"]} & LAB_MEASUREMENT_CLASSES
    for r in df.itertuples(index=False):
        eid = r.evidence_id
        if not re.fullmatch(r"PLE-\d{3}", eid):
            problems.append(f"{eid}: evidence_id format")
        for c in _split(r.condition_ids) or ["<none>"]:
            if c not in conds:
                problems.append(f"{eid}: unknown condition id {c!r}")
        mcs = _split(r.measurement_classes)
        for m in mcs:
            if m not in classes:
                problems.append(f"{eid}: measurement class {m!r} not a lab class in configs/measurements.yaml")
        if not mcs and not r.measurement_class_gap.strip():
            problems.append(f"{eid}: no measurement class and no measurement_class_gap explanation")
        for col, allowed in (("claim_type", CLAIM_TYPES), ("finding_direction", FINDING_DIRECTIONS),
                             ("quote_source", QUOTE_SOURCES), ("evidence_tier", EVIDENCE_TIERS),
                             ("external_validation", EXTERNAL_VALIDATION), ("incorporation_bias", INCORPORATION_BIAS),
                             ("comparator_type", COMPARATOR_TYPES), ("analyte_class", ANALYTE_CLASSES)):
            if getattr(r, col) not in allowed:
                problems.append(f"{eid}: {col}={getattr(r, col)!r} not in {sorted(allowed)}")
        for col in ("auroc", "sensitivity", "specificity"):
            v = getattr(r, col)
            if v:
                try:
                    f = float(v)
                except ValueError:
                    problems.append(f"{eid}: {col}={v!r} is not a number")
                    continue
                if not 0.0 <= f <= 1.0:
                    problems.append(f"{eid}: {col}={v!r} outside [0, 1]")
        for col in ("n_cases", "n_controls"):
            v = getattr(r, col)
            if v and not v.isdigit():
                problems.append(f"{eid}: {col}={v!r} is not a non-negative integer")
        if not r.pmid.isdigit():
            problems.append(f"{eid}: pmid {r.pmid!r}")
        if not r.doi.strip():
            problems.append(f"{eid}: empty doi")
        if r.quote_source == "europepmc_fulltext" and not r.pmcid.startswith("PMC"):
            problems.append(f"{eid}: europepmc_fulltext needs a PMCID")
        if not r.source_quote.strip() or not r.risk_of_bias_notes.strip():
            problems.append(f"{eid}: empty source_quote or risk_of_bias_notes")
        if r.incorporation_bias in {"yes", "partial"} and "incorporation" not in r.risk_of_bias_notes.lower():
            problems.append(f"{eid}: incorporation_bias={r.incorporation_bias} but risk_of_bias_notes do not say so")
        if r.claim_type == "diagnostic_accuracy" and not (r.auroc or r.sensitivity or r.specificity
                                                          or "accuracy" in r.metric.lower()):
            problems.append(f"{eid}: diagnostic_accuracy row without AUROC, sensitivity, specificity or accuracy")
    return problems


# --------------------------------------------------------------------------- source texts
def write_source_texts(texts: dict[str, dict]) -> None:
    p = raw_dir(SOURCE_ID) / TEXTS_NAME
    content = "".join(json.dumps(texts[k], ensure_ascii=False, sort_keys=True) + "\n" for k in sorted(texts))
    if p.exists() and p.read_text(encoding="utf-8") == content and TEXTS_NAME in load_manifest(SOURCE_ID)["files"]:
        return  # unchanged: keep the recorded retrieved_at so offline rebuilds are deterministic
    p.write_text(content, encoding="utf-8")
    record_file(SOURCE_ID, TEXTS_NAME, url=EPMC_SEARCH,
                note="Europe PMC core records (abstract + metadata) and open-access fullTextXML where a quote needs "
                     "it; the texts every source_quote is checked against")


def read_source_texts() -> dict[str, dict]:
    p = RAW / SOURCE_ID / TEXTS_NAME
    if not p.exists():
        return {}
    return {d["pmid"]: d for d in map(json.loads, p.read_text(encoding="utf-8").splitlines()) if d}


# --------------------------------------------------------------------------- build
def build(fetch: bool = True) -> pd.DataFrame:
    df = load_curated()
    problems = validate_curated(df)
    if problems:
        raise ValueError("curated table invalid:\n  " + "\n  ".join(problems))
    texts = pe.fetch_source_texts(df) if fetch else read_source_texts()
    if fetch:
        write_source_texts(texts)
    out = pe.verify_quotes(df, texts)
    untraced = [pe.untraced_numbers(r) for _, r in df.iterrows()]
    out["untraced_numbers"] = ["; ".join(u) for u in untraced]
    out["numbers_traced_to_quote"] = [not u for u in untraced]
    for col in ("title", "journal", "pub_year", "pub_types"):
        out[f"source_{col}"] = [
            ("; ".join(texts.get(p, {}).get(col) or []) if col == "pub_types" else texts.get(p, {}).get(col, ""))
            for p in out.pmid]
    out["source_url"] = [u or f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for u, p in zip(out.source_url_override, out.pmid)]
    out["europepmc_url"] = [f"https://europepmc.org/article/MED/{p}" for p in out.pmid]
    out["claim_label"] = CLAIM_LABEL
    out["project_reproduction"] = "not attempted"
    out["reports_person_level_discrimination"] = out.claim_type.isin(DISCRIMINATION_CLAIMS)
    out["object_id"] = "published_lab_evidence:" + out.evidence_id
    for col in ("n_cases", "n_controls"):
        out[col] = pe._to_num(out[col], int)
    for col in ("auroc", "sensitivity", "specificity"):
        out[col] = pe._to_num(out[col])
    manifest = load_manifest(SOURCE_ID)["files"]
    retrieved = (manifest.get(TEXTS_NAME) or {}).get("retrieved_at") or utc_now_iso()
    out = add_provenance(
        out, data_layer="condition_molecular",
        source_name="Peer-reviewed literature via Europe PMC / PubMed (curated lab/biomarker-evidence extraction)",
        source_version=f"curated {CURATED_ON}; {len(out)} claims from {out.pmid.nunique()} publications",
        retrieved_at=retrieved, evidence_type="published_biomarker", source_record_id="evidence_id",
        evidence_level=out["evidence_tier"],
        provenance_notes=CLAIM_LABEL + ". Numbers are copied from source_quote (verbatim; checked in quote_verified); "
                                       "counts marked in count_derivation were summed by the curator. Condition-level "
                                       "published evidence: no row describes a participant of any dataset in this "
                                       "project.")
    write_table(out, TABLE, producer="measure_it.labs.published_lab_evidence",
                description="Published lab / biomarker evidence for target invisible illnesses; " + CLAIM_LABEL)
    return out


def summary(df: pd.DataFrame) -> dict:
    return {
        "rows": int(len(df)), "publications": int(df.pmid.nunique()),
        "quote_verified": int(df.quote_verified.sum()), "doi_matches_source": int(df.doi_matches_source.sum()),
        "numbers_traced_to_quote": int(df.numbers_traced_to_quote.sum()),
        "by_condition": df.condition_ids.str.split(";").explode().value_counts().to_dict(),
        "by_claim_type": df.claim_type.value_counts().to_dict(),
        "by_finding_direction": df.finding_direction.value_counts().to_dict(),
        "incorporation_bias": df.incorporation_bias.value_counts().to_dict(),
        "comparator_type": df.comparator_type.value_counts().to_dict(),
        "diagnostic_accuracy_rows": int((df.claim_type == "diagnostic_accuracy").sum()),
        "diagnostic_accuracy_rows_with_external_validation": int(
            ((df.claim_type == "diagnostic_accuracy") & (df.external_validation == "yes")).sum()),
    }


# --------------------------------------------------------------------------- query
PUBLIC_COLUMNS = [
    "object_id", "evidence_id", "claim_label", "condition_ids", "condition_definition_in_study", "biomarker",
    "specimen", "assay_method", "analyte_class", "measurement_classes", "measurement_class_gap", "claim_type",
    "finding_direction", "incorporation_bias", "comparator_type", "study_design", "evidence_tier", "n_cases",
    "n_controls", "count_derivation", "control_group", "reference_standard", "population_setting", "metric",
    "metric_value", "ci_or_spread", "auroc", "sensitivity", "specificity", "external_validation",
    "risk_of_bias_notes", "source_quote", "quote_source", "quote_verified", "pmid", "doi", "first_author", "year",
    "source_title", "source_url", "project_reproduction", "reports_person_level_discrimination",
]
_SEARCH_FIELDS = ("biomarker", "specimen", "assay_method", "analyte_class", "measurement_classes",
                  "measurement_class_gap")


@lru_cache(maxsize=1)
def _table(_stamp: float) -> pd.DataFrame:
    return read_table(TABLE)


def _load_table() -> pd.DataFrame | None:
    if not table_exists(TABLE):
        return None
    from ..store import processed_path
    return _table(processed_path(TABLE).stat().st_mtime)


def _clean(v):
    if v is None or v is pd.NA or (isinstance(v, float) and pd.isna(v)):
        return UNKNOWN
    if hasattr(v, "item"):
        v = v.item()
    return UNKNOWN if v == "" else v


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower())


def get_published_lab_evidence(condition: str | None = None, biomarker: str | None = None) -> dict:
    """Published lab/biomarker claims for a condition and/or a biomarker.

    `biomarker` matches an analyte_class or measurement class exactly (e.g. 'autoantibody', 'metabolomics'), otherwise
    every word (>= 3 characters, or a short token such as 'NK' / 'LTE4') must appear in the biomarker, specimen,
    assay, analyte class or measurement-class fields ('tryptase', 'urinary LTE4', 'CA-125', 'norepinephrine').
    Returns {'status': 'matched', 'claim_label', 'filters', 'n_rows', 'rows': [...]} with each row's design,
    comparator, incorporation-bias flag, risk-of-bias notes and verbatim quote, or {'status': UNKNOWN / NOT AVAILABLE,
    'reason'}. Missing values are UNKNOWN, never 0.
    """
    df = _load_table()
    if df is None:
        return {"status": UNKNOWN, "reason": f"table {TABLE} not built (run measure_it.labs.published_lab_evidence)"}
    if not condition and not biomarker:
        return {"status": UNKNOWN, "reason": "give a condition and/or a biomarker"}
    filters: dict = {}
    mask = pd.Series(True, index=df.index)
    if condition:
        from ..measurements.query import resolve_condition
        rc = resolve_condition(condition)
        if rc.get("status") != "matched":
            return {"status": UNKNOWN, "reason": rc.get("reason", "condition not resolved"), "condition": condition,
                    "candidates": rc.get("candidates", [])}
        cid = rc["condition_id"]
        filters["condition_id"] = cid
        mask &= df.condition_ids.map(lambda s: cid in _split(s))
    if biomarker:
        b = biomarker.strip().lower()
        if b in ANALYTE_CLASSES or b in LAB_MEASUREMENT_CLASSES:
            filters.update(biomarker=biomarker, match_method="class")
            mask &= (df.analyte_class == b) | df.measurement_classes.map(lambda s: b in _split(s))
        else:
            words = [w for w in _norm(biomarker).split() if len(w) > 2 or re.search(r"\d", w) or w in {"nk", "ca"}]
            if not words:
                return {"status": UNKNOWN, "reason": f"biomarker {biomarker!r} has no searchable words"}
            filters.update(biomarker=biomarker, match_method="text", words=words)
            hay = df[list(_SEARCH_FIELDS)].astype(str).agg(" ".join, axis=1).map(_norm)
            def _has(s: str, w: str) -> bool:  # short tokens ('nk', 'ca') must be whole words
                return (f" {w} " in f" {s} ") if len(w) <= 2 else (w in s)
            mask &= hay.map(lambda s: all(_has(s, w) for w in words))
    sub = df[mask]
    if sub.empty:
        return {"status": UNKNOWN, "reason": "no curated published lab-evidence row matches", "filters": filters,
                "claim_label": CLAIM_LABEL}
    sub = sub.sort_values(["condition_ids", "evidence_tier", "evidence_id"])
    rows = [{c: _clean(r[c]) for c in PUBLIC_COLUMNS if c in sub.columns} for _, r in sub.iterrows()]
    return {"status": "matched", "claim_label": CLAIM_LABEL, "filters": filters, "n_rows": len(rows),
            "n_null_or_mixed": int(sub.finding_direction.isin(["null", "mixed"]).sum()),
            "n_incorporation_bias": int(sub.incorporation_bias.isin(["yes", "partial"]).sum()),
            "n_symptomatic_comparator": int(sub.comparator_type.isin(["disease_or_symptomatic", "mixed"]).sum()),
            "n_with_external_validation": int((sub.external_validation == "yes").sum()), "rows": rows}


# --------------------------------------------------------------------------- registry / CLI
def registry_entry(df: pd.DataFrame) -> dict:
    return {
        "source_id": SOURCE_ID,
        "name": "Published lab and biomarker evidence for invisible illnesses (curated literature extraction)",
        "publisher": "Peer-reviewed journals, indexed by PubMed / Europe PMC (EMBL-EBI); extraction curated by this project",
        "landing_url": "https://europepmc.org/",
        "access_urls": [EPMC_SEARCH, EPMC_FULLTEXT],
        "license": "Abstracts and metadata via the Europe PMC REST API (open); full texts used are open-access "
                   "articles; quotes are short verbatim extracts for verification. Copyright stays with publishers.",
        "access_conditions": "Open; no registration or data-use agreement",
        "retrieved_at": df["retrieved_at"].iloc[0],
        "source_version": f"curated {CURATED_ON}",
        "update_date": CURATED_ON,
        "data_layer": "condition_molecular",
        "evidence_type": "published_biomarker",
        "unit_of_observation": "one published claim (publication x condition x biomarker x metric)",
        "sample_size": {"claims": int(len(df)), "publications": int(df.pmid.nunique()),
                        "quote_verified": int(df.quote_verified.sum())},
        "geographic_resolution": "none",
        "person_level": False,
        "geographic": False,
        "omics": False,  # published aggregate claims (some are metabolomic/proteomic), no omics data held
        "wearable": False,
        "participant_linkage": "none (aggregate published results)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested" if bool(df.quote_verified.all() and df.numbers_traced_to_quote.all()) else "partial",
        "processed_outputs": [TABLE],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.labs.published_lab_evidence",
        "limitations": "Curated, not a systematic review: searches were targeted (see DATA_AUDIT.md), so absence of a "
                       "row is not evidence of absence. Every number is a published claim, not reproduced by this "
                       "project. Most rows compare patients with healthy controls (spectrum bias), several biomarkers "
                       "are part of the diagnostic criteria they are scored against (incorporation bias, flagged per "
                       "row), and almost none report external validation.",
    }


def run(fetch: bool = True) -> pd.DataFrame:
    from ..registry import write_registry_entry
    df = build(fetch=fetch)
    write_registry_entry(registry_entry(df))
    print(json.dumps(summary(df), indent=2))
    bad = df.loc[~df.quote_verified, ["evidence_id", "quote_missing_segments"]]
    if len(bad):
        print("QUOTES NOT FOUND:\n" + bad.to_string(index=False))
    un = df.loc[~df.numbers_traced_to_quote, ["evidence_id", "untraced_numbers"]]
    if len(un):
        print("NUMBERS NOT TRACED TO THE QUOTE:\n" + un.to_string(index=False))
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-fetch", action="store_true", help="verify against data/raw/.../source_texts.jsonl only")
    a = ap.parse_args()
    run(fetch=not a.no_fetch)
