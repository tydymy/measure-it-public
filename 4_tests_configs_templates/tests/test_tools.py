"""Tests for the shared tool facade (measure_it.tools) and the deterministic deployment agent.

Pure tests run anywhere; tests marked `data` need data/processed. Data.gov calls run with MEASURE_IT_OFFLINE=1 so
they are answered from the HTTP cache (the demo query) or return UNKNOWN (never the network).
"""
from __future__ import annotations

import json
import re

import pytest

from measure_it import tools as T
from measure_it.config import UNKNOWN
from measure_it.validate import forbidden_hits

NONSENSE = "zzqxv blorf quux"

# tool -> (demo arguments, nonsense arguments)
CASES = {
    "search_condition": ({"condition": "long covid"}, {"condition": NONSENSE}),
    "normalize_condition": ({"condition_or_code": "G93.32"}, {"condition_or_code": NONSENSE}),
    "get_patient_phenotype_signature": ({"condition": "ME/CFS"}, {"condition": NONSENSE}),
    "get_molecular_context": ({"condition": "Long COVID"}, {"condition": NONSENSE}),
    "discover_candidate_measurements": ({"condition": "Long COVID", "phenotype": "orthostatic intolerance"},
                                        {"condition": NONSENSE}),
    "get_measurement_evidence": ({"measurement": "wearable autonomic monitoring", "condition": "Long COVID"},
                                 {"measurement": NONSENSE, "condition": "Long COVID"}),
    "get_regulatory_context": ({"technology": "wearable ECG patch"}, {"technology": NONSENSE}),
    "get_condition_burden": ({"condition": "ME/CFS", "geography_level": "state"}, {"condition": NONSENSE}),
    "get_geographic_context": ({"geography_id": "San Diego County, California"},
                               {"geography_id": "Zzqxv County, Nowhere"}),
    "find_candidate_clinics": ({"geography_id": "06073", "condition": "Long COVID",
                                "measurement": "wearable autonomic monitoring"},
                               {"geography_id": "99999", "condition": "Long COVID"}),
    "find_relevant_trials": ({"condition": "Long COVID", "measurement": "wearable autonomic monitoring"},
                             {"condition": NONSENSE}),
    "find_relevant_research_centers": ({"condition": "ME/CFS"}, {"condition": NONSENSE}),
    "rank_deployment_opportunities": ({"condition": "Long COVID or ME/CFS",
                                       "measurement": "wearable autonomic monitoring", "top_n": 3},
                                      {"condition": NONSENSE, "measurement": "wearable autonomic monitoring"}),
    "search_us_open_data": ({"query": "long COVID"}, {"query": NONSENSE}),
    "trace_evidence": ({"object_id": "condition:long_covid"}, {"object_id": "zzqxv:blorf"}),
    "get_phenotype_measurement_evidence": ({"phenotype": "orthostatic intolerance"}, {"phenotype": NONSENSE}),
    "get_measurable_biology": ({"condition": "ME/CFS"}, {"condition": NONSENSE}),
    "list_sources": ({}, {"data_layer": NONSENSE}),
}
ENVELOPE_KEYS = {"tool", "status", "query", "data", "caveats", "provenance", "truncation", "meta"}


@pytest.fixture
def offline(monkeypatch):
    monkeypatch.setenv("MEASURE_IT_OFFLINE", "1")


# ------------------------------------------------------------------------------------------------ pure

def test_every_spec_tool_is_exposed():
    spec = ["search_condition", "normalize_condition", "get_patient_phenotype_signature", "get_molecular_context",
            "discover_candidate_measurements", "get_measurement_evidence", "get_regulatory_context",
            "get_condition_burden", "get_geographic_context", "find_candidate_clinics", "find_relevant_trials",
            "find_relevant_research_centers", "rank_deployment_opportunities", "search_us_open_data", "trace_evidence"]
    assert T.SPEC_TOOLS == spec
    assert set(spec) | {"get_phenotype_measurement_evidence", "get_measurable_biology", "list_sources"} == set(T.TOOLS)
    assert set(CASES) == set(T.TOOLS)


def test_collect_object_ids_takes_whole_values_only():
    d = {"a": "condition:long_covid", "b": ["trial:NCT01", "see trial:NCT02 in text", "geo:06073"],
         "c": {"d": "fda:code:DQA", "e": "not:an_id", "f": "opportunity:x|y|06073"}}
    assert T.collect_object_ids(d) == ["condition:long_covid", "trial:NCT01", "geo:06073", "fda:code:DQA",
                                       "opportunity:x|y|06073"]


def test_cap_lists_records_truncation_and_protects_primary_lists():
    trunc = []
    out = T.cap_lists({"rows": list(range(10)), "x": {"y": list(range(4))}}, 3, "data", trunc,
                      protect=("data.rows",))
    assert out["rows"] == list(range(10)) and out["x"]["y"] == [0, 1, 2]
    assert trunc == [{"path": "data.x.y", "returned": 3, "total": 4}]
    merged = T._merge_trunc([{"path": "p", "returned": 3, "total": 5}, {"path": "p", "returned": 3, "total": 9}])
    assert merged == [{"path": "p", "returned": 3, "total": 9, "n_lists": 2}]


def test_fit_size_never_cuts_below_floor():
    trunc = []
    data = {"big": [{"s": "x" * 100} for _ in range(200)], "small": list(range(8))}
    out = T.fit_size(data, 2_000, trunc)
    assert len(out["big"]) == 10 and out["small"] == list(range(8))
    assert any(t["path"] == "data.big" and t["total"] == 200 for t in trunc)


# ------------------------------------------------------------------------------------------------ data: every tool

def _check_envelope(r: dict, tool: str) -> None:
    assert ENVELOPE_KEYS <= set(r), set(r) ^ ENVELOPE_KEYS
    assert r["tool"] == tool
    assert r["status"] in (T.OK, UNKNOWN, T.ERROR)
    json.dumps(r, allow_nan=False)                       # NaN / inf became null
    assert isinstance(r["caveats"], list) and T.GENERAL_CAVEAT in r["caveats"]
    assert isinstance(r["truncation"], list)
    assert {"object_ids", "tables", "sources", "producer"} <= set(r["provenance"])
    hits = forbidden_hits(json.dumps(r, indent=1))
    assert not hits, hits[:3]


@pytest.mark.data
@pytest.mark.parametrize("tool", list(CASES))
def test_tool_ok_on_demo_input(tool, offline):
    r = T.TOOLS[tool](**CASES[tool][0])
    _check_envelope(r, tool)
    assert r["status"] == T.OK, r.get("reason")
    assert "reason" not in r


@pytest.mark.data
@pytest.mark.parametrize("tool", list(CASES))
def test_tool_unknown_on_nonsense(tool, offline):
    r = T.TOOLS[tool](**CASES[tool][1])
    _check_envelope(r, tool)
    assert r["status"] == UNKNOWN, (tool, r["status"], r.get("reason"))
    assert r["reason"]


@pytest.mark.data
def test_level_d_burden_is_unknown_with_reason():
    r = T.get_condition_burden("POTS", "county")
    assert r["status"] == UNKNOWN
    assert r["data"]["burden_evidence_level"] == "D" and r["data"]["value"] == UNKNOWN
    assert "no" in r["reason"].lower()


@pytest.mark.data
def test_inherited_state_values_are_labelled_not_county_prevalence():
    r = T.get_condition_burden("Long COVID", "county", max_rows=5)
    assert r["status"] == T.OK
    summ = r["data"]["summary"]
    if summ["n_inherited_state_values"]:   # long_covid_burden: hps_inherited
        assert all(row["inherited"] is True and row["source_geographic_resolution"] == "state"
                   for row in r["data"]["rows"])
        assert any("never describe them as county prevalence" in c for c in r["caveats"])
    else:                                  # brfss_sae (default since 2026-10-07): modelled, labelled
        assert summ["n_modelled_small_area_values"] > 0
        assert all(row["derivation"] == "mrp_small_area_estimate" for row in r["data"]["rows"])
        assert any("Never describe them as observed county prevalence" in c for c in r["caveats"])
    assert any(t["path"] == "data.rows" and t["returned"] == 5 for t in r["truncation"])
    assert r["data"]["burden_evidence_level"] == "A"


@pytest.mark.data
def test_molecular_context_is_condition_level_enrichment():
    r = T.get_molecular_context("ME/CFS")
    assert r["status"] == T.OK
    assert r["data"]["label"].startswith("condition-level molecular enrichment")
    assert "coherence" in r["data"] and "post_hoc_flag_note" in r["data"]["coherence"]
    text = json.dumps(r).lower()
    for m in re.finditer(r"multi-omics", text):            # only ever mentioned as what it is NOT
        assert re.search(r"\b(not|never|no)\b", text[max(0, m.start() - 60):m.start()]), text[m.start() - 80:m.end()]


@pytest.mark.data
def test_ambiguous_condition_is_not_guessed():
    rc = T.resolve_condition_arg("fatigue")
    assert rc["status"] == UNKNOWN and rc["reason"]
    r = T.get_molecular_context("fatigue")
    assert r["status"] == UNKNOWN


@pytest.mark.data
def test_rank_compact_keeps_spec_schema_and_shared_context():
    r = T.rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=3)
    d = r["data"]
    spec = {"condition", "phenotype", "measurement", "technology", "geography", "candidate_sites",
            "research_evidence", "molecular_context", "uncertainties", "provenance", "recommended_next_step"}
    assert len(d["recommendations"]) == 3
    for rec in d["recommendations"]:
        assert spec <= set(rec)
        g = rec["geography"]
        assert g["burden_evidence_level"] in ("A", "B", "C", "D")
        assert {"name", "fips", "burden", "vulnerability", "diagnostic_desert"} <= set(g)
        assert rec["provenance"]["object_ids"][0].startswith("opportunity:")
        assert "candidate deployment opportunity" in rec["recommended_next_step"].lower()
    assert {"phenotype", "molecular_context", "technology", "condition", "measurement"} <= set(d["shared_context"])
    assert len(json.dumps(r, separators=(",", ":"))) < 120_000
    # equal weights: rank and the (equal-weight) Monte Carlo interval describe the same ranking; no basis block
    assert "monte_carlo_basis" not in d["recommendations"][0]["geography"]


@pytest.mark.data
def test_rank_compact_names_the_basis_of_a_non_equal_rank():
    """Under another weight set the compact geography block must not show that set's rank next to the equal-weight
    Monte Carlo interval without saying which is which."""
    r = T.rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=2,
                                        weight_set="burden_led")
    for rec in r["data"]["recommendations"]:
        g = rec["geography"]
        assert g["rank_basis"]["weight_set"] == "burden_led" and g["rank_basis"]["rank_column"] == "rank_burden_led"
        assert "EQUAL weights" in g["monte_carlo_basis"] and "not the rank under 'burden_led'" in g["monte_carlo_basis"]
        assert g["rank_equal_weights"] is not None and "rank_interval_5_95" in g


@pytest.mark.data
def test_find_relevant_trials_accepts_bundle_class_and_alias():
    b = T.find_relevant_trials("ME/CFS", "wearable autonomic monitoring")
    c = T.find_relevant_trials("ME/CFS", "accelerometry")
    assert b["status"] == c["status"] == T.OK
    assert b["data"]["summary"]["n_trials"] >= c["data"]["summary"]["n_trials"] > 0
    assert all(t["object_id"].startswith("trial:") for t in b["data"]["trials"])
    assert b["data"]["view"].startswith("literal")


@pytest.mark.data
def test_results_are_memoised():
    T.clear_cache()
    a = T.get_regulatory_context("pulse oximeter")
    b = T.get_regulatory_context("pulse oximeter")
    assert a["meta"]["cache_hit"] is False and b["meta"]["cache_hit"] is True
    assert b["meta"]["elapsed_s"] < 0.05
    assert {k: v for k, v in a.items() if k != "meta"} == {k: v for k, v in b.items() if k != "meta"}


# ------------------------------------------------------------------------------------------------ data: agent

def _cited_lines(res: dict) -> list[str]:
    out = []
    for key in ("what", "where", "who", "uncertain", "prov"):
        out += res["sections"][key]
    return out


@pytest.mark.data
def test_agent_capillaroscopy_question_every_sentence_cites_ids():
    from measure_it.agents.deployment_agent import CAPILLAROSCOPY_QUESTION, PERSON_NOTE, answer_question
    res = answer_question(CAPILLAROSCOPY_QUESTION, top_n=3)
    assert res["stopped"] is None
    assert res["parsed"]["measurement"]["measurement_id"] == "nailfold_capillaroscopy"
    assert res["parsed"]["phenotype"]["id"] == "orthostatic_intolerance"
    assert [c["condition_id"] for c in res["parsed"]["conditions"]] == ["dysautonomia"]
    for line in _cited_lines(res):
        if line == PERSON_NOTE:
            continue
        assert re.search(r"\[[a-z_]+:[^\]]+\]\s*$", line), line
    assert not forbidden_hits(res["answer_markdown"])
    # every cited id resolves
    for oid in res["cited_object_ids"]:
        t = T.trace_evidence(oid)
        assert t["status"] == T.OK, (oid, t.get("reason"))
    tools = [c["tool"] for c in res["tool_calls"]]
    assert tools.index("normalize_condition") < tools.index("rank_deployment_opportunities") < tools.index(
        "find_candidate_clinics") < tools.index("trace_evidence")
    assert any("burden" in u.lower() and "evidence level d" in u.lower() for u in res["unknowns"])


@pytest.mark.data
def test_agent_unknown_condition_stops_without_ranking():
    from measure_it.agents.deployment_agent import answer_question
    res = answer_question("Where should we deploy wearable monitoring for Zorblaxian drift syndrome?")
    assert res["stopped"]
    assert not res["parsed"]["conditions"]
    assert "rank_deployment_opportunities" not in [c["tool"] for c in res["tool_calls"]]
    assert any(UNKNOWN in u for u in res["unknowns"])


@pytest.mark.data
def test_agent_level_d_condition_keeps_burden_unknown():
    from measure_it.agents.deployment_agent import answer_question
    res = answer_question("Where should we deploy autonomic testing for POTS?", top_n=3)
    assert [c["condition_id"] for c in res["parsed"]["conditions"]] == ["pots"]
    assert res["parsed"]["measurement"]["measurement_id"] == "autonomic_function_testing"   # bundle alias preferred
    where = " ".join(res["sections"]["where"])
    assert f"{UNKNOWN} (evidence level D)" in where
    assert "burden percentile " + UNKNOWN in where
    for oid in res["cited_object_ids"]:
        assert T.trace_evidence(oid)["status"] == T.OK, oid


@pytest.mark.data
def test_agent_parses_named_geography_and_cites_its_stored_opportunity():
    from measure_it.agents.deployment_agent import answer_question, parse_geography
    assert parse_geography("Deploy CPET in St. Louis County, Missouri and nearby")["geo_id"] == "29189"
    assert parse_geography("Where in Michigan should we deploy it?")["geo_id"] == "26"
    assert parse_geography("Where in the United States should we deploy it?") is None
    res = answer_question("Could we deploy wearable monitoring for Long COVID in San Diego County, California?",
                          top_n=2)
    assert res["parsed"]["geography"]["geo_id"] == "06073"
    opp = "opportunity:long_covid|wearable_autonomic_activity_monitoring|06073"
    assert opp in res["cited_object_ids"]
    line = next(x for x in res["sections"]["where"] if opp in x)
    assert ("inherited by the county, not county prevalence" in line
            or "modelled small-area estimate, not observed county prevalence" in line)
    clinics = [c for c in res["tool_calls"] if c["tool"] == "find_candidate_clinics"]
    assert clinics and clinics[0]["args"]["geography_id"] == "06073"


@pytest.mark.data
def test_molecular_entities_carry_traceable_object_ids():
    """Genes, variants, GWAS studies and GEO series in the molecular context are cited with CONVENTIONS ids."""
    assert T.molecular_object_id("gene", "ENSG00000109670") == "gene:ENSG00000109670"
    assert T.molecular_object_id("variant", "rs141691232") == "variant:rs141691232"
    assert T.molecular_object_id("study", "GCST90480593") == "gwas_study:GCST90480593"
    assert T.molecular_object_id("study", "GSE227375") == "geo_series:GSE227375"
    assert T.molecular_object_id("drug", "CHEMBL1431") is None and T.molecular_object_id("gene", "CA14") is None
    r = T.get_molecular_context("ME/CFS", top_n=3)
    ids = [i for i in r["provenance"]["object_ids"] if not i.startswith("condition:")]
    assert {i.split(":", 1)[0] for i in ids} >= {"gene", "gwas_study", "geo_series", "variant"}
    for oid in ids:
        assert T.trace_evidence(oid)["status"] == T.OK, oid
    rank = T.rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", top_n=1)
    shared = T.collect_object_ids(rank["data"]["shared_context"]["molecular_context"])
    assert any(i.startswith("gene:") for i in shared)


@pytest.mark.data
def test_regulatory_context_cites_matched_measurement_classes():
    r = T.get_regulatory_context("nailfold capillaroscopy")
    assert r["status"] == T.OK
    assert "measurement:capillaroscopy" in r["provenance"]["object_ids"]
    assert T.trace_evidence("measurement:capillaroscopy")["status"] == T.OK


@pytest.mark.data
def test_agent_named_state_is_not_presented_as_a_within_state_ranking():
    """A named state: the county ranking is labelled national, the state's stored row is a rank among states, and the
    facility pool is the whole state (no radius)."""
    from measure_it.agents.deployment_agent import answer_question
    res = answer_question("Where should we deploy wearable autonomic monitoring for Long COVID or ME/CFS in Texas?",
                          top_n=3)
    assert res["parsed"]["geography"]["geo_id"] == "48"
    where = res["sections"]["where"]
    assert any("NATIONAL ranking" in x for x in where)
    assert any(x.startswith("Counties in Texas among the national top") or x.startswith("No county in Texas")
               for x in where)
    stored = next(x for x in where if "opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|48]" in x)
    assert "state-level ranking" in stored and "among states" in stored
    who = " ".join(res["sections"]["who"])
    assert "inside the state" in who and "within 50 km" not in who
    for oid in res["cited_object_ids"]:
        assert T.trace_evidence(oid)["status"] == T.OK, oid


@pytest.mark.data
def test_agent_spans_use_the_shared_vocabulary():
    """The agent finds measurement spans in the shared resolver's vocabulary (class aliases included) and phenotype
    spans in the curated axis aliases ('PEM' -> post_exertional_malaise); spans are confirmed by the same resolver."""
    from measure_it.agents.deployment_agent import parse_measurements, parse_phenotypes
    from measure_it.facilities.matching import resolve_measurement
    ph = parse_phenotypes("Which wearable could capture PEM in ME/CFS?")
    assert [(p["kind"], p["id"]) for p in ph] == [("axis", "post_exertional_malaise")]
    ms = parse_measurements("Could a Holter monitor or a pulse oximeter help?")
    assert [m["lexicon_target"] for m in ms] == ["ecg_ambulatory", "continuous_spo2"]
    for m in ms:
        assert resolve_measurement(m["lexicon_term"])["measurement_id"] == m["lexicon_target"]
    wb = parse_measurements("Where should we deploy wearable-based monitoring?")
    assert wb and resolve_measurement(wb[0]["lexicon_term"])["measurement_id"] == "wearable_autonomic_activity_monitoring"


# ------------------------------------------------------------------ review fixes (2026-09-24)
def test_agent_system_tiers_separate_flagged_and_literature_support():
    """The agent must not call drug-target-driven or literature-only systems 'supported' (results/MOLECULAR_COHERENCE.md)."""
    from measure_it.agents.deployment_agent import describe_system_tiers, system_support_tiers
    rows = [{"physiological_system": "immune", "source": "ot_literature", "post_hoc_flagged": False},
            {"physiological_system": "autonomic_cardiac", "source": "ot_other", "post_hoc_flagged": True},
            {"physiological_system": "neuronal", "source": "ot_literature", "post_hoc_flagged": False},
            {"physiological_system": "neuronal", "source": "ot_other", "post_hoc_flagged": True},
            {"physiological_system": "vascular_endothelial", "source": "ot_genetic", "post_hoc_flagged": False}]
    tiers = system_support_tiers(rows)
    assert tiers == {"non_literature_unflagged": ["vascular_endothelial"],
                     "drug_target_flagged": ["autonomic_cardiac", "neuronal"], "literature_only": ["immune"]}
    text = describe_system_tiers(tiers)
    assert "supported physiological systems" not in text and "post-hoc-flagged" in text
    assert system_support_tiers(UNKNOWN) == {"non_literature_unflagged": [], "drug_target_flagged": [],
                                            "literature_only": []}


@pytest.mark.data
def test_demo_conditions_have_no_unflagged_non_literature_system(offline):
    """What the report states in section 8, read through the tool the agent uses."""
    from measure_it.agents.deployment_agent import system_support_tiers
    for c in ("long_covid", "me_cfs", "pots"):
        r = T.get_molecular_context(c, top_n=3)
        tiers = system_support_tiers((r["data"].get("coherence") or {}).get("supported_systems"))
        assert tiers["non_literature_unflagged"] == [], c


@pytest.mark.data
def test_offline_search_matches_a_cached_query_ignoring_case(offline):
    exact = T.search_us_open_data("long COVID")
    variant = T.search_us_open_data("  Long   COVID ")
    assert exact["status"] == variant["status"] == T.OK
    assert variant["data"]["n_results"] == exact["data"]["n_results"]
    assert any("cached response to the query 'long COVID'" in c for c in variant["caveats"])
    miss = T.search_us_open_data("ME/CFS")
    assert miss["status"] == UNKNOWN and "cached queries" in miss["reason"]


@pytest.mark.data
def test_technology_evidence_carries_the_phenotype_signal_verdict():
    from measure_it.scoring import recommend as REC
    ctx = REC.technology_context("wearable_autonomic_activity_monitoring",
                                 ("accelerometry", "wearable_heart_rate", "hrv", "ecg_ambulatory", "ppg",
                                  "posture_detection"), ("long_covid", "me_cfs"))
    vals = [t.get("phenotype_signal_strength") for t in ctx["technology_evidence"] if "measurement_id" in t]
    assert vals and all(v is not None for v in vals)
    assert "null_result" in vals


def test_agent_never_calls_omics_signature_features_wearable():
    from measure_it.agents.deployment_agent import _feature_kind
    assert _feature_kind("nhanes") == "wearable features"
    assert _feature_kind("stanford_longcovid_uwakwe2025") == "wearable features"
    for ds in ("geo_GSE270045", "mapmecfs_nih_pi_mecfs", None):
        assert "not wearable" in _feature_kind(ds)
