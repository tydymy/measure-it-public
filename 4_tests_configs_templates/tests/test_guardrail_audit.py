"""Regression tests for the guardrail audit (docs/GUARDRAIL_AUDIT.md).

Each test pins one violation the audit found and fixed at its source:
* guesses instead of UNKNOWN (partial/ambiguous condition text in get_patient_phenotype_signature and the Condition
  Explorer; fuzzy word overlap for a measurement-bundle name in get_regulatory_context);
* 'nan' strings for level-D members of a condition set;
* agent answers that silently accepted an out-of-scope premise (individual diagnosis, 'best' clinic, participant
  location or omics, causes / treatment);
* a proxy series shown without its burden evidence level (geography Test 4 report);
* forbidden product terms in user-facing text that `measure-it validate` did not scan (MCP/API descriptions, README,
  captions, data audits, result JSON).
Pure tests run anywhere; tests marked `data` need data/processed.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from measure_it.config import PROCESSED, UNKNOWN
from measure_it.validate import forbidden_hits

HAS_DATA = (PROCESSED / "deployment_opportunities.parquet").exists()
needs_data = pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")


# ------------------------------------------------------------------------------------------------ pure

def test_agent_scope_notes_flag_out_of_scope_premises():
    from measure_it.agents.deployment_agent import (CAPILLAROSCOPY_QUESTION, PERSON_LOCATION_NOTE, SCOPE_RULES,
                                                    scope_notes)
    q = "Does the patient have POTS, and which clinic is best in Texas for diagnosing it with a smartwatch?"
    notes = scope_notes(q)
    assert notes[:2] == [SCOPE_RULES[0][1], SCOPE_RULES[1][1]]
    q2 = "What is the multi-omics profile of NHANES participants living in San Diego County, and which gene causes it?"
    notes2 = scope_notes(q2, {"geo_id": "06073"})
    assert SCOPE_RULES[2][1] in notes2 and SCOPE_RULES[3][1] in notes2 and PERSON_LOCATION_NOTE in notes2
    assert scope_notes(CAPILLAROSCOPY_QUESTION) == []
    assert scope_notes("Where are the diagnostic deserts for POTS?") == []       # 'diagnostic' is not 'diagnose'
    for rx, note in SCOPE_RULES:
        assert not forbidden_hits(note), note
    assert not forbidden_hits(PERSON_LOCATION_NOTE)


def test_test4_legend_names_every_series_level():
    from measure_it.geography.features import test4_series_legend
    corr = pd.DataFrame({"series_id": ["lyme_county", "long_covid_county_proxy", "long_covid_county_proxy",
                                       "desert_me_cfs_county"],
                         "label": ["Lyme incidence 2023 (A) / reported cases", "PLACES PHLTH crude (C proxy) / x",
                                   "PLACES PHLTH crude (C proxy) / x", "diagnostic desert, me_cfs"]})
    leg = test4_series_legend(corr)
    assert "`long_covid_county_proxy` = PLACES PHLTH crude (C proxy)" in leg and "(A)" in leg
    assert leg.count("long_covid_county_proxy") == 1 and "`desert_me_cfs_county`" not in leg
    assert test4_series_legend(pd.DataFrame()) == ""


def test_language_check_reads_json_strings_and_service_text():
    from measure_it import validate as V
    doc = {"a": ["never a \"best clinic\" here"], "b": {"c": "this facility is the best clinic"}}
    hits = [h[1] for s in V._json_strings(doc) for h in forbidden_hits(s)]
    assert hits == ["best clinic"]                                     # the quoted mention is excused, the claim is not
    svc = V.service_texts()
    assert "tool:find_candidate_clinics" in svc and "mcp.GUARDRAILS_MD" in svc
    bad = {k: forbidden_hits(v) for k, v in svc.items() if forbidden_hits(v)}
    assert not bad, bad


# ------------------------------------------------------------------------------------------------ data

@pytest.mark.data
@needs_data
@pytest.mark.parametrize("text", ["chronic", "syndrome", "headache", "HP:0012378"])
def test_phenotype_signature_does_not_guess_from_partial_matches(text):
    from measure_it import tools as T
    r = T.get_patient_phenotype_signature(text)
    assert r["status"] == UNKNOWN, (text, [d["phenotype_id"] for d in (r["data"] or {}).get("datasets") or []])


@pytest.mark.data
@needs_data
def test_normalize_condition_partial_or_ambiguous_is_unknown_with_candidates():
    from measure_it import tools as T
    for text, status in (("chronic", "partial"), ("HP:0012378", "ambiguous")):
        r = T.normalize_condition(text)
        assert r["status"] == UNKNOWN and r["data"]["match_status"] == status and r["reason"], text
        assert r["data"]["matches"] or r["data"]["other_candidates"]            # candidates kept, none chosen
    assert T.normalize_condition("G93.32")["status"] == T.OK


@pytest.mark.data
@needs_data
def test_phenotype_signature_still_answers_exact_inputs():
    from measure_it import tools as T
    for text, pid in (("ME/CFS", "mecfs_like_proxy"), ("PASC", "long_covid_self_report"),
                      ("depression", "depression_phq9"), ("fatigue", "fatigue")):
        r = T.get_patient_phenotype_signature(text)
        assert r["status"] == T.OK and pid in [d["phenotype_id"] for d in r["data"]["datasets"]], text


@pytest.mark.data
@needs_data
def test_regulatory_context_resolves_bundles_not_fuzzy_words():
    from measure_it import tools as T
    r = T.get_regulatory_context("wearable autonomic monitoring", max_records=2)
    assert r["status"] == T.OK
    got = [(m["measurement_id"], m["method"]) for m in r["data"]["matched_measurements"]]
    assert {m for m, _ in got} == {"accelerometry", "wearable_heart_rate", "hrv", "ecg_ambulatory", "ppg",
                                   "posture_detection"}
    assert all(meth == "exact" for _, meth in got) and "autonomic_testing" not in {m for m, _ in got}
    assert r["data"]["measurement_resolution"]["kind"] == "bundle"
    # an exact / pattern class match is kept as it was (the bundle alias 'tilt table' does not widen it)
    t = T.get_regulatory_context("tilt table", max_records=2)
    assert [m["measurement_id"] for m in t["data"]["matched_measurements"]] == ["autonomic_testing"]
    # fuzzy-only matches are labelled as candidates, not an identification
    b = T.get_regulatory_context("blood test", max_records=2)
    assert any("fuzzy vocabulary" in c for c in b["caveats"])
    assert T.get_regulatory_context("quantum banana")["status"] == UNKNOWN


@pytest.mark.data
@needs_data
def test_level_d_set_members_have_no_nan_strings():
    from measure_it import tools as T
    r = T.rank_deployment_opportunities("demo cluster", "wearable autonomic monitoring", "county", top_n=2)
    assert r["status"] == T.OK
    text = json.dumps(r)
    assert '"nan"' not in text
    for rec in r["data"]["recommendations"]:
        pots = rec["geography"]["burden"]["members"]["pots"]
        assert pots["burden_evidence_level"] == "D" and pots["burden_value"] is None
        assert pots["burden_source_resolution"] == UNKNOWN and pots["burden_measure_id"] is None


@pytest.mark.data
@needs_data
def test_agent_answer_states_scope_for_individual_questions():
    from measure_it.agents.deployment_agent import SCOPE_RULES, answer_question
    res = answer_question("Does the patient have POTS, and which clinic is best in Texas?", top_n=2)
    md = res["answer_markdown"]
    assert "## Scope (what this engine does not do)" in md
    assert SCOPE_RULES[0][1] in md and SCOPE_RULES[1][1] in md
    assert not forbidden_hits(md)


@pytest.mark.data
@needs_data
def test_condition_explorer_does_not_show_a_guessed_condition(monkeypatch):
    from pathlib import Path

    from streamlit.testing.v1 import AppTest
    monkeypatch.setenv("MEASURE_IT_OFFLINE", "1")
    page = Path(__file__).resolve().parents[1] / "dashboard" / "pages" / "2_Condition_Explorer.py"
    at = AppTest.from_file(str(page), default_timeout=300).run()
    at.text_input(key="cond_text").input("chronic").run()
    assert not at.exception
    warnings = [str(w.value) for w in at.warning]
    assert any(UNKNOWN in w and "partial" in w for w in warnings), warnings
    subs = [s.value for s in at.subheader]
    assert "1. Phenotype features" not in subs and "Limitations and provenance" in subs
