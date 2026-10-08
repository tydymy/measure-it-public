"""Tests for measure_it.labs.published_lab_evidence (curated published lab / biomarker evidence) and the All of Us
lab-concept counts (measure_it.labs.aou_lab_concepts).

Unit tests use synthetic rows (a fabricated number or a missing incorporation-bias note must fail). Data tests check
the curated CSV, the offline re-verification of every quote against the stored Europe PMC texts, the processed
table, the query helper, the doc and the All of Us rows.
"""
from __future__ import annotations

import pandas as pd
import pytest

from measure_it.config import DOCS, RAW, TABLES, UNKNOWN
from measure_it.labs import aou_lab_concepts as aou
from measure_it.labs import published_lab_evidence as ple
from measure_it.provenance import PROVENANCE_COLUMNS, check_provenance


def _row(**kw) -> dict:
    base = {c: "" for c in ple.REQUIRED_COLUMNS}
    base.update(evidence_id="PLE-900", condition_ids="mcas", condition_definition_in_study="x", biomarker="tryptase",
                specimen="serum", assay_method="immunoassay", analyte_class="mast_cell_mediator",
                measurement_classes="blood_biomarkers", claim_type="diagnostic_accuracy", finding_direction="positive",
                incorporation_bias="no", comparator_type="healthy", study_design="case-control",
                evidence_tier="small_primary", n_cases="40", n_controls="20",
                count_derivation="as stated in source_quote", control_group="healthy", reference_standard="x",
                population_setting="x", metric="sensitivity", metric_value="80%", sensitivity="0.80",
                external_validation="no", risk_of_bias_notes="healthy controls only (spectrum bias)",
                source_quote="40 patients and 20 controls; sensitivity 80%", quote_source="europepmc_abstract",
                pmid="123", doi="10.1/x", first_author="A", year="2020")
    base.update(kw)
    return base


# --------------------------------------------------------------------------- unit
def test_validate_accepts_a_clean_synthetic_row():
    assert ple.validate_curated(pd.DataFrame([_row()])) == []


def test_validate_rejects_bad_enums_and_missing_incorporation_note():
    df = pd.DataFrame([
        _row(evidence_id="PLE-901", analyte_class="vibes"),
        _row(evidence_id="PLE-902", comparator_type="friends"),
        _row(evidence_id="PLE-903", incorporation_bias="yes"),  # notes do not mention incorporation
        _row(evidence_id="PLE-904", measurement_classes="accelerometry"),  # a device class, not a lab class
        _row(evidence_id="PLE-905", quote_source="pdf_text"),
        _row(evidence_id="PDE-001"),
    ])
    problems = " | ".join(ple.validate_curated(df))
    for needle in ("vibes", "friends", "PLE-903: incorporation_bias=yes", "accelerometry", "pdf_text",
                   "PDE-001: evidence_id format"):
        assert needle in problems, needle


def test_number_tracer_rejects_a_fabricated_sensitivity():
    import measure_it.measurements.published_evidence as pe
    assert pe.untraced_numbers(_row()) == []
    assert pe.untraced_numbers(_row(sensitivity="0.85")) == ["sensitivity=0.85"]
    assert pe.untraced_numbers(_row(n_cases="41")) == ["n_cases=41"]


def test_aou_shown_uses_the_disclosure_floor():
    assert aou.shown(20) == "<=20" and aou.shown(0) == "0" and aou.shown(4420) == "4420"
    assert aou.shown(None) == UNKNOWN


def test_aou_concept_rows_match_by_code_not_name():
    data = {"responses": {"MEASUREMENT|tryptase": {"response": {"items": [
        {"vocabularyId": "LOINC", "conceptCode": "7748-7", "conceptId": 1, "conceptName": "Tryptase IgE Ab",
         "countValue": 40, "sourceCountValue": 20},
        {"vocabularyId": "LOINC", "conceptCode": "21582-2", "conceptId": 3019420,
         "conceptName": "Tryptase [Mass/volume] in Serum or Plasma", "countValue": 4420, "sourceCountValue": 3340},
    ]}}}}
    rows = {(r["vocabulary"], r["code"]): r for r in aou.concept_rows(data)}
    assert rows[("LOINC", "21582-2")]["participants"] == 4420
    assert rows[("LOINC", "51834-0")]["found"] is False  # absent from this fake response: UNKNOWN, not 0
    assert rows[("LOINC", "51834-0")]["participants_shown"] == UNKNOWN


# --------------------------------------------------------------------------- data: curated CSV and quotes
@pytest.mark.data
def test_curated_csv_is_valid():
    df = ple.load_curated()
    assert len(df) >= 50
    assert ple.validate_curated(df) == []


@pytest.mark.data
def test_every_target_condition_has_rows():
    df = ple.load_curated()
    conds = set(df.condition_ids.str.split(";").explode())
    for c in ("mcas", "pots", "long_covid", "me_cfs", "fibromyalgia", "eds_hsd", "lyme_disease", "ptlds",
              "migraine", "ibs", "gastroparesis", "endometriosis"):
        assert c in conds, c


@pytest.mark.data
def test_stored_source_texts_reverify_every_quote_offline():
    import measure_it.measurements.published_evidence as pe
    df = ple.load_curated()
    texts = ple.read_source_texts()
    assert set(df.pmid) <= set(texts), "source_texts.jsonl is missing PMIDs: re-run the module"
    v = pe.verify_quotes(df, texts)
    bad = v.loc[~v.quote_verified, ["evidence_id", "quote_missing_segments"]]
    assert bad.empty, bad.to_string()
    assert v.doi_matches_source.all()
    untraced = {r.evidence_id: pe.untraced_numbers(r._asdict()) for r in df.itertuples(index=False)}
    assert not {k: u for k, u in untraced.items() if u}


# --------------------------------------------------------------------------- data: processed table
@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    from measure_it.store import read_table
    return read_table(ple.TABLE)


@pytest.mark.data
def test_table_provenance_and_ids(table):
    check_provenance(table, ple.TABLE)
    assert set(PROVENANCE_COLUMNS) <= set(table.columns)
    assert set(table.data_layer) == {"condition_molecular"}
    assert set(table.evidence_type) == {"published_biomarker"}
    assert table.object_id.str.fullmatch(r"published_lab_evidence:PLE-\d{3}").all()
    assert table.object_id.is_unique


@pytest.mark.data
def test_table_every_claim_verified_and_labelled(table):
    assert table.quote_verified.all()
    assert table.numbers_traced_to_quote.all()
    assert table.doi_matches_source.all()
    assert (table.claim_label == ple.CLAIM_LABEL).all()
    assert (table.project_reproduction == "not attempted").all()
    assert ple.CLAIM_LABEL in table.provenance_notes.iloc[0]


@pytest.mark.data
def test_nulls_and_incorporation_bias_are_kept(table):
    exploded = table.assign(c=table.condition_ids.str.split(";")).explode("c")
    nulls = exploded[exploded.finding_direction.isin(["null", "mixed"])]
    for cond in ("me_cfs", "long_covid", "pots"):
        assert (nulls.c == cond).any(), f"no null/mixed row kept for {cond}"
    # the MCAS tryptase criterion is part of the MCAS definition: at least one MCAS row must be flagged
    mcas = exploded[exploded.c == "mcas"]
    assert mcas.incorporation_bias.isin(["yes", "partial"]).any()
    # hEDS: the table must say that no validated biomarker exists
    assert ((exploded.c == "eds_hsd") & (exploded.claim_type == "no_validated_biomarker")).any()


@pytest.mark.data
def test_diagnostic_accuracy_rows_carry_bias_notes(table):
    da = table[table.claim_type == "diagnostic_accuracy"]
    assert len(da) >= 5
    assert (da.risk_of_bias_notes.str.len() > 40).all()
    assert da.reports_person_level_discrimination.all()
    assert not table.loc[table.claim_type != "diagnostic_accuracy", "reports_person_level_discrimination"].any()


# --------------------------------------------------------------------------- data: query helper
@pytest.mark.data
def test_query_condition_and_biomarker_text():
    r = ple.get_published_lab_evidence("MCAS", "tryptase")
    assert r["status"] == "matched" and r["claim_label"] == ple.CLAIM_LABEL
    assert all("mcas" in row["condition_ids"].split(";") for row in r["rows"])
    assert all("tryptase" in (row["biomarker"] + row["assay_method"]).lower() for row in r["rows"])
    assert all(row["source_quote"] and row["risk_of_bias_notes"] for row in r["rows"])


@pytest.mark.data
def test_query_by_class_and_unknowns():
    r = ple.get_published_lab_evidence(biomarker="autoantibody")
    assert r["status"] == "matched" and r["filters"]["match_method"] == "class"
    assert {row["analyte_class"] for row in r["rows"]} == {"autoantibody"}
    assert ple.get_published_lab_evidence("definitely not a condition xyz")["status"] == UNKNOWN
    assert ple.get_published_lab_evidence("migraine", "unobtainium")["status"] == UNKNOWN
    assert ple.get_published_lab_evidence()["status"] == UNKNOWN


@pytest.mark.data
def test_query_missing_values_are_unknown_not_zero():
    r = ple.get_published_lab_evidence("EDS")
    assert r["status"] == "matched"
    row = next(x for x in r["rows"] if x["claim_type"] == "no_validated_biomarker")
    assert row["auroc"] == UNKNOWN and row["sensitivity"] == UNKNOWN


# --------------------------------------------------------------------------- data: doc, registry, All of Us
@pytest.mark.data
def test_doc_lists_every_row_and_the_label(table):
    doc = (DOCS / "PUBLISHED_LAB_EVIDENCE.md").read_text(encoding="utf-8")
    assert "published claim — not reproduced by this project" in doc
    missing = [e for e in table.evidence_id if e not in doc]
    assert not missing, missing


@pytest.mark.data
def test_registry_entry_written():
    import yaml
    e = yaml.safe_load((RAW / ple.SOURCE_ID / "registry_entry.yaml").read_text())
    assert e["evidence_type"] == "published_biomarker" and e["data_layer"] == "condition_molecular"
    assert e["person_level"] is False and e["true_participant_linkage_across_modalities"] is False
    assert (RAW / ple.SOURCE_ID / "DATA_AUDIT.md").exists()


@pytest.mark.data
def test_aou_lab_counts_in_cohort_table_match_the_raw_json():
    data = aou.load_raw()
    assert data is not None, "run python -m measure_it.labs.aou_lab_concepts"
    counts = pd.read_csv(TABLES / "controlled_cohort_counts.csv", dtype=str, keep_default_na=False)
    assert list(counts.columns) == ["cohort", "quantity", "value", "source_url", "retrieved_date", "status"]
    expected = aou.cohort_count_rows(data, data["retrieved_at"][:10])
    have = {(r.cohort, r.quantity): r.value for r in counts.itertuples(index=False)}
    for r in expected:
        assert have.get((r["cohort"], r["quantity"])) == r["value"], r["quantity"]
    tryptase = next(r for r in aou.concept_rows(data) if r["code"] == "21582-2")
    assert tryptase["participants"] >= 1000
    # CGRP has no Labs & Measurements concept: recorded as 0 concepts, never as 0 participants
    cgrp = [r for r in expected if "CGRP" in r["quantity"]]
    assert cgrp and all("number of CONCEPTS" in r["quantity"] for r in cgrp)


@pytest.mark.data
def test_data_access_plan_has_the_lab_section():
    doc = (DOCS / "DATA_ACCESS_PLAN.md").read_text(encoding="utf-8")
    assert "### Lab biomarkers in All of Us" in doc
    assert "21582-2" in doc and "ascertainment" in doc.lower()
