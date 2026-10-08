"""Tests for measure_it.measurements.published_evidence (curated published device evidence).

Unit tests check the quote verifier and the number tracer on synthetic text (a fabricated number must fail).
Data tests check the curated CSV, the processed table, the offline re-verification against the stored source texts,
the query helper and the doc.
"""
from __future__ import annotations

import pandas as pd
import pytest

from measure_it.config import DOCS, UNKNOWN
from measure_it.measurements import published_evidence as pe
from measure_it.provenance import PROVENANCE_COLUMNS, check_provenance


# --------------------------------------------------------------------------- unit: normalisation and quote checks
def test_normalize_strips_markup_but_keeps_inequalities():
    raw = "HR <i>P</i> &lt; .001; <h4>Results</h4>ME/CFS &#x2212;10.8 at Test 1, P < 0.05 and > 10 y"
    n = pe.normalize_for_match(raw)
    assert "p<.001" in n
    assert "results" in n and "<h4>" not in n
    assert "-10.8attest1" in n              # U+2212 minus and no-break space unified
    assert "p<0.05and>10y" in n            # literal < and > in a sentence are text, not a tag


def test_check_quote_finds_segments_and_rejects_fabrication():
    src = ("<h4>Results</h4>Patients with POTS (n=15) and healthy controls (n=15) underwent testing. "
           "Sp (specificity) was 67% for STAND10.")
    ok = "Patients with POTS (n=15) and healthy controls (n=15) || Sp (specificity) was 67%"
    assert pe.check_quote(ok, src) == (2, 2, [])
    fabricated = "Patients with POTS (n=15) || Sp (specificity) was 76%"
    n, found, missing = pe.check_quote(fabricated, src)
    assert (n, found) == (2, 1) and missing == ["Sp (specificity) was 76%"]


def test_untraced_numbers_flags_values_absent_from_the_quote():
    row = {"source_quote": "sensitivity of 0.950 ± 0.062 and an AUC of 0.900; Eighty-seven healthy subjects",
           "auroc": "0.900", "sensitivity": "0.950", "specificity": "", "n_cases": "", "n_controls": "87",
           "count_derivation": "as stated in source_quote"}
    assert pe.untraced_numbers(row) == []
    assert pe.untraced_numbers({**row, "specificity": "0.85"}) == ["specificity=0.85"]
    assert pe.untraced_numbers({**row, "n_controls": "88"}) == ["n_controls=88"]
    # a count the curator summed is not required to appear verbatim
    assert pe.untraced_numbers({**row, "n_controls": "88", "count_derivation": "31 + 57 (sum computed by curator)"}) == []


def test_rate_forms_cover_percent_and_decimal():
    forms = pe._rate_forms("0.9111")
    assert "91.11%" in forms and "0.9111" in forms
    assert "87%" in pe._rate_forms("0.87")


# --------------------------------------------------------------------------- data: curated CSV
@pytest.mark.data
def test_curated_csv_is_valid():
    df = pe.load_curated()
    assert len(df) >= 60
    assert pe.validate_curated(df) == []


@pytest.mark.data
def test_validate_curated_catches_bad_ids():
    df = pe.load_curated().head(3).copy()
    df.loc[0, "condition_ids"] = "not_a_condition"
    df.loc[1, "measurement_classes"] = "tricorder"
    df.loc[2, "claim_type"] = "proof"
    problems = " | ".join(pe.validate_curated(df))
    assert "not_a_condition" in problems and "tricorder" in problems and "proof" in problems


@pytest.mark.data
def test_stored_source_texts_reverify_every_quote_offline():
    """Guards against editing the CSV without re-running the build: every quote must still be in the stored text."""
    df = pe.load_curated()
    texts = pe.read_source_texts()
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
    return read_table(pe.TABLE)


@pytest.mark.data
def test_table_provenance_and_ids(table):
    check_provenance(table, pe.TABLE)
    assert set(PROVENANCE_COLUMNS) <= set(table.columns)
    assert set(table.data_layer) == {"measurement"}
    assert set(table.evidence_type) == {"published_biomarker"}
    assert table.object_id.str.fullmatch(r"published_evidence:PDE-\d{3}").all()
    assert table.object_id.is_unique


@pytest.mark.data
def test_table_every_claim_verified_and_labelled(table):
    assert table.quote_verified.all()
    assert table.numbers_traced_to_quote.all()
    assert table.doi_matches_source.all()
    assert (table.claim_label == pe.CLAIM_LABEL).all()
    assert pe.CLAIM_LABEL in table.provenance_notes.iloc[0]


@pytest.mark.data
def test_published_and_reproduced_kept_apart(table):
    rep = dict(zip(table.evidence_id, table.project_reproduction))
    assert rep["PDE-040"].startswith("attempted: NOT reproduced")
    assert {v for k, v in rep.items() if k != "PDE-040"} == {"not attempted"}


@pytest.mark.data
def test_nulls_are_reported_beside_positives(table):
    exploded = table.assign(c=table.condition_ids.str.split(";")).explode("c")
    nulls = exploded[exploded.finding_direction.isin(["null", "mixed"])]
    for cond in ("me_cfs", "long_covid", "pots"):
        assert (nulls.c == cond).any(), f"no null/mixed row kept for {cond}"
    # the independent 2-day CPET replication and the pooled Long COVID HRV result stay null
    fd = dict(zip(table.evidence_id, table.finding_direction))
    assert fd["PDE-015"] == "null" and fd["PDE-038"] == "null"


@pytest.mark.data
def test_diagnostic_accuracy_rows_carry_bias_notes(table):
    da = table[table.claim_type == "diagnostic_accuracy"]
    assert len(da) >= 10
    assert (da.risk_of_bias_notes.str.len() > 40).all()
    assert da.external_validation.isin(["yes", "no", "not_applicable"]).all()
    assert da.reports_person_level_discrimination.all()
    assert not table.loc[table.claim_type != "diagnostic_accuracy", "reports_person_level_discrimination"].any()


# --------------------------------------------------------------------------- data: query helper
@pytest.mark.data
def test_query_condition_and_measurement_class():
    r = pe.get_published_device_evidence("POTS", "wearable")
    assert r["status"] == "matched" and r["claim_label"] == pe.CLAIM_LABEL
    ids = {row["evidence_id"] for row in r["rows"]}
    assert "PDE-006" in ids
    assert all("pots" in row["condition_ids"].split(";") for row in r["rows"])
    assert all(row["claim_label"] == pe.CLAIM_LABEL and row["source_quote"] for row in r["rows"])


@pytest.mark.data
def test_query_device_text_fallback_and_unknowns():
    r = pe.get_published_device_evidence("Long COVID", "pupillometry")
    assert r["status"] == "matched" and r["filters"]["match_method"] == "device_text"
    assert [row["evidence_id"] for row in r["rows"]] == ["PDE-049"]
    assert pe.get_published_device_evidence("MCAS")["status"] == UNKNOWN
    assert pe.get_published_device_evidence("definitely not a condition xyz")["status"] == UNKNOWN
    cpet = pe.get_published_device_evidence("ME/CFS", "CPET")
    assert cpet["status"] == "matched" and cpet["n_null_or_mixed"] >= 1


@pytest.mark.data
def test_query_missing_values_are_unknown_not_zero():
    r = pe.get_published_device_evidence("fibromyalgia")
    row = next(x for x in r["rows"] if x["evidence_id"] == "PDE-031")
    assert row["auroc"] == UNKNOWN and row["n_cases"] == UNKNOWN


# --------------------------------------------------------------------------- data: doc and registry
@pytest.mark.data
def test_doc_lists_every_row_and_the_label(table):
    doc = (DOCS / "PUBLISHED_DEVICE_EVIDENCE.md").read_text(encoding="utf-8")
    assert "published claim — not reproduced by this project" in doc
    missing = [e for e in table.evidence_id if e not in doc]
    assert not missing, missing


@pytest.mark.data
def test_registry_entry_written():
    import yaml
    from measure_it.config import RAW
    e = yaml.safe_load((RAW / pe.SOURCE_ID / "registry_entry.yaml").read_text())
    assert e["evidence_type"] == "published_biomarker" and e["data_layer"] == "measurement"
    assert e["person_level"] is False and e["true_participant_linkage_across_modalities"] is False
    assert (RAW / pe.SOURCE_ID / "DATA_AUDIT.md").exists()
