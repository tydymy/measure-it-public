"""Tests for measure_it.measurements.{discovery,precision_audit,evidence,query}.

Unit tests need no data. Tests marked `data` read data/processed outputs (run
`uv run python -m measure_it.measurements.discovery` and `... .evidence` first).
"""
from __future__ import annotations

import json
import re

import pandas as pd
import pytest

from measure_it.config import TABLES, UNKNOWN, load_config
from measure_it.measurements import discovery as D
from measure_it.measurements import phenotype_evidence as PE
from measure_it.measurements import precision_audit as PA
from measure_it.measurements import query as Q
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table, table_exists

CLASS_IDS = [c["id"] for c in load_config("measurements")["measurement_classes"]]
BUNDLE_IDS = list(load_config("relevance")["measurement_bundles"])
CONDITION_IDS = [c["id"] for c in load_config("conditions")["conditions"]]


def _cc(mid: str, version: str = "final") -> D.CompiledClass:
    return next(c for c in D.compile_classes(D.measurement_classes(version)) if c.id == mid)


def _hit(text: str, mid: str, version: str = "final") -> dict | None:
    return D.scan_text(text, _cc(mid, version))


# ------------------------------------------------------------------------------------------ unit
@pytest.mark.parametrize("text", ["Fatigue Severity Scale", "PROMIS fatigue short form", "SF-36 physical function",
                                  "Pittsburgh Sleep Quality Index", "COMPASS-31 total score", "self-reported steps",
                                  "sleep diary", "Borg rating of perceived exertion", "PHQ-9"])
def test_pro_regex_flags_questionnaires(text):
    assert D.pro_hit(text) is not None


@pytest.mark.parametrize("text", ["large-scale cohort of actigraphy", "pro-inflammatory cytokines", "heart rate variability",
                                  "polysomnography", "vas deferens", "the time scale of recovery"])
def test_pro_regex_spares_objective_text(text):
    assert D.pro_hit(text) is None


def test_scan_text_objective_and_pro_flags():
    h = _hit("Daily step count measured by a wrist accelerometer for 7 days", "accelerometry")
    assert h["objective_flag"] and "[[" in h["snippet"]
    h = _hit("Physical activity measured with the IPAQ questionnaire and an accelerometer", "accelerometry")
    assert h is not None and not h["objective_flag"] and h["pro_hit"].lower() == "questionnaire"


def test_sleep_exclude_pattern_blocks_psqi():
    h = _hit("PSQI sleep efficiency component", "sleep_objective", "v1")
    assert h is not None and not h["objective_flag"] and h["exclude_hit"]


def test_nonhuman_context_flag():
    h = _hit("Cardiac function will be quantified by echocardiography in telemetry-instrumented rats.", "vascular_imaging")
    assert h is not None and h["nonhuman_context_hit"].lower() == "rats"


def test_local_window_ignores_hard_wraps_but_keeps_list_items():
    text = ("Background sentence one about something else entirely. " * 3
            + "In this aim telemetry instrumented\nrats will be studied with echocardiography of the heart\nand more text. "
            + "Another sentence that is unrelated to the measurement. " * 3)
    i = text.find("echocardiography")
    w = D.local_window(text, i, i + 16)
    assert "rats" in w and "Background" not in w
    elig = ("Inclusion Criteria:\n\n* Adults aged 18 or over with a long description that pads the field beyond the "
            "short-field limit so that windowing applies to this criterion text xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            "* No ECG abnormality at screening\n* Able to complete a questionnaire")
    i = elig.find("ECG")
    w = D.local_window(elig, i, i + 3)
    assert "ECG" in w and "questionnaire" not in w


@pytest.mark.parametrize("text,mid", [
    ("Mast Cell Activation Syndrome (MCAS)", "immune_assays"),           # T-cell pattern must need a word boundary
    ("Valutazione Degli Esiti a Medio-lungo Termine", "ecg_ambulatory"),  # 'zio' inside an Italian word
    ("16S rRNA sequencing of stool", "transcriptomics"),
    ("a digital biomarker of fatigue from wearables", "blood_biomarkers"),
    ("resting 12-lead ECG at screening", "ecg_ambulatory"),
    ("30 seconds sit-to-stand test", "digital_gait"),
    ("pulse wave velocity by applanation tonometry", "ppg"),
    ("vital signs including respiratory rate", "respiratory_rate"),
    ("DPA-714 PET/MRI of neuroinflammation", "vascular_imaging"),
])
def test_revised_patterns_reject_round1_false_positives(text, mid):
    h = _hit(text, mid)
    assert h is None or not h["objective_flag"]


@pytest.mark.parametrize("text,mid", [
    ("24-hour Holter ECG recording", "ecg_ambulatory"), ("Zio XT patch for 14 days", "ecg_ambulatory"),
    ("inertial measurement unit worn for one week", "digital_gait"), ("serum C-reactive protein", "blood_biomarkers"),
    ("lowest oxygen saturation during the 6-minute walk test", "continuous_spo2"),
    ("resting heart rate from a Garmin watch", "wearable_heart_rate"), ("PBMC RNA-seq", "transcriptomics"),
    ("intraepidermal nerve fiber density on skin biopsy", "small_fiber_testing"),
    ("sleep efficiency assessed by activity tracker", "sleep_objective"), ("flow cytometry of T cells", "immune_assays"),
])
def test_revised_patterns_keep_true_uses(text, mid):
    h = _hit(text, mid)
    assert h is not None and h["objective_flag"]


def test_v1_patterns_are_the_shipped_config():
    """The config is the single source of truth: live patterns = the audited revision, v1 = audit_revision.*_v1."""
    cfg = load_config("measurements")["measurement_classes"]
    live = {c["id"]: c["patterns"] for c in cfg}
    shipped = {c["id"]: (c.get("audit_revision") or {}).get("patterns_v1", c["patterns"]) for c in cfg}
    assert {c["id"]: c["patterns"] for c in D.measurement_classes("v1")} == shipped
    assert {c["id"]: c["patterns"] for c in D.measurement_classes("final")} == live
    revised = {c["id"]: c["patterns"] for c in D.measurement_classes("final")}
    assert revised["accelerometry"] == shipped["accelerometry"]          # passed round 1, untouched
    assert revised["ecg_ambulatory"] != shipped["ecg_ambulatory"]
    assert set(D.PATTERN_REVISIONS) == {c["id"] for c in cfg if c.get("audit_revision")}
    for mid, rev in D.PATTERN_REVISIONS.items():
        assert rev.get("changes"), f"{mid}: every revision must document its change and rationale"
        for p in rev.get("patterns", []) + rev.get("exclude_patterns", []):
            re.compile(p)
    # the revision hash (and so every mention's pattern_version) is the one the audit recorded
    assert D.pattern_version_label("final").endswith("sha1 12e9f79f)")


def test_regulatory_matching_uses_audited_patterns():
    from measure_it.measurements import regulatory as reg
    ms = D.measurement_classes("final")
    ids = lambda q: [m["measurement_id"] for m in reg.match_measurement_classes(q, ms)]  # noqa: E731
    assert "immune_assays" not in ids("Mast Cell Activation")          # shipped 'T[- ]cell activation' matched it
    assert ids("mast cell activation syndrome") == []
    assert ids("T-cell activation assay") == ["immune_assays"]
    for q in ("Pittsburgh Sleep Quality Index", "COMPASS-31 questionnaire", "sleep diary", "fatigue severity scale",
              "self-reported sleep quality"):
        assert ids(q) == [], q                                          # questionnaires are not measurement classes
        assert reg.get_regulatory_context(q)["status"] == reg.UNKNOWN


def test_wilson_and_seed():
    lo, hi = PA.wilson(30, 30)
    assert 0.88 < lo < 0.9 and hi == 1.0
    assert PA.wilson(0, 0)[0] != PA.wilson(0, 0)[0]  # nan
    assert PA.class_seed("hrv", 1) == PA.class_seed("hrv", 1) != PA.class_seed("hrv", 2)


def test_resolve_measurement_bundles_and_aliases():
    r = Q.resolve_measurement("wearable autonomic monitoring")
    assert r["status"] == "matched" and r["kind"] == "bundle" and "hrv" in r["member_ids"]
    # an exact class id wins over a bundle alias (reviewer fix): 'capillaroscopy' is the class, not the bundle that
    # also contains microvascular_function; 'CPET' is the cpet class; bundle-only aliases still resolve to bundles
    r = Q.resolve_measurement("CPET")
    assert r["status"] == "matched" and r["kind"] == "class" and r["measurement_ids"] == ["cpet"]
    r = Q.resolve_measurement("capillaroscopy")
    assert r["kind"] == "class" and r["measurement_ids"] == ["capillaroscopy"]
    r = Q.resolve_measurement("CAPRIO")
    assert r["kind"] == "bundle" and r["measurement_ids"] == ["nailfold_capillaroscopy"]
    assert Q.resolve_measurement("2-day CPET")["measurement_ids"] == ["exercise_capacity_testing"]
    r = Q.resolve_measurement("accelerometry")
    assert r["measurement_ids"] == ["accelerometry"]
    assert Q.resolve_measurement("Holter")["measurement_ids"] == ["ecg_ambulatory"]
    assert Q.resolve_measurement("tarot cards")["status"] == UNKNOWN
    assert Q.resolve_measurement("")["status"] == UNKNOWN
    amb = Q.resolve_measurement("blood test")      # only loose word matches to several classes -> no guess
    assert amb["status"] == UNKNOWN and len(amb["candidates"]) > 1


AGREE = ["wearable", "Wearable", "wearables", "wearable device", "wearable monitoring", "wearable autonomic monitoring",
         "CPET", "2-day CPET", "2 day cpet", "exercise testing", "cardiopulmonary exercise testing", "capillaroscopy",
         "CAPRIO", "tilt table", "autonomic testing", "autonomic", "heart rate", "Heart-rate variability", "HRV",
         "accelerometry", "accelerometer", "actigraphy", "Holter", "EKG", "ECG", "wearable ECG", "pulse oximeter",
         "SpO2", "polysomnography", "microvascular", "orthostatic", "blood test", "smartwatch", "tarot cards", "",
         "PEM", "not a measurement"]


@pytest.mark.parametrize("text", AGREE)
def test_every_consumer_resolves_measurements_the_same_way(text):
    """One resolver (measure_it.measurements.resolve): the scoring/facility path (single target) and the evidence/tool
    path agree on every phrase (before 2026-09-24 'wearable' was the wearable bundle for one and ambiguous for the
    other, 'accelerometer' / 'Holter' / 'pulse oximeter' resolved on one path only)."""
    from measure_it.facilities import matching as M
    from measure_it.measurements.resolve import resolve_measurement as R
    a, b, c = M.resolve_measurement(text), Q.resolve_measurement(text), R(text)
    assert a["status"] == b["status"] == c["status"], text
    if a["status"] == "matched":
        assert [a["measurement_id"]] == b["measurement_ids"] == c["measurement_ids"], text
        assert a["members"] == b["member_ids"]
    else:
        assert a.get("candidates") == b.get("candidates"), text


def test_obvious_aliases_resolve():
    from measure_it.measurements.resolve import resolve_measurement as R
    for text, target in [("wearable", "wearable_autonomic_activity_monitoring"),
                         ("wearable device", "wearable_autonomic_activity_monitoring"),
                         ("accelerometer", "accelerometry"), ("Holter", "ecg_ambulatory"), ("EKG", "ecg_ambulatory"),
                         ("pulse oximeter", "continuous_spo2"), ("SpO2", "continuous_spo2"),
                         ("polysomnography", "sleep_objective"), ("exercise testing", "exercise_capacity_testing")]:
        r = R(text)
        assert r["status"] == "matched" and r["measurement_id"] == target, (text, r)
    r = R("wearable")
    assert r["kind"] == "bundle" and r["match_step"] == "exact"
    for text in ("heart rate", "autonomic", "orthostatic"):          # ambiguous: never guessed
        r = R(text)
        assert r["status"] == UNKNOWN and len(r["candidates"]) > 1, text
    # the phenotype side: PEM / post-exertional malaise are the post_exertional_malaise axis
    from measure_it.store import read_table, table_exists
    if table_exists("phenotype_axes"):
        ax = read_table("phenotype_axes")
        for text in ("PEM", "post-exertional malaise", "PESE"):
            p = Q.resolve_phenotype(text, ax)
            assert p["status"] == "matched" and p["axis_id"] == "post_exertional_malaise", text
            assert p["match_reason"].startswith("exact")


def test_clustered_interval_and_adequacy():
    j = pd.Series(["yes"] * 10)
    one_record = PA.clustered_wilson(j, pd.Series(["trial:A"] * 10))
    spread = PA.clustered_wilson(j, pd.Series([f"trial:{i}" for i in range(10)]))
    assert one_record["n_records"] == 1 and one_record["n_eff_records"] == 1.0
    assert spread["n_eff_records"] == 10.0
    assert one_record["precision_strict_ci_low_clustered"] < spread["precision_strict_ci_low_clustered"]
    assert abs(spread["precision_strict_ci_low_clustered"] - PA.wilson(10, 10)[0]) < 1e-9
    assert PA.sample_adequacy(10, 5) == "adequate"
    assert PA.sample_adequacy(10, 4) == "insufficient" and PA.sample_adequacy(9, 9) == "insufficient"


def test_consumer_applicant_regex_excludes_hospital_ultrasound():
    from measure_it.measurements import evidence as E
    rx = re.compile(E.CONSUMER_APPLICANTS, re.IGNORECASE)
    assert rx.search("Samsung Electronics Co., Ltd.") and rx.search("Apple, Inc.") and rx.search("Withings")
    assert not rx.search("Samsung Medison Co., Ltd.")


def test_ratio_ci():
    from measure_it.measurements import evidence as E
    lo, hi = E.ratio_ci(74, 1174, 10, 5414)
    r = (74 / 1174) / (10 / 5414)
    assert lo < r < hi and lo > 1
    assert all(v != v for v in E.ratio_ci(0, 10, 3, 10))     # NaN when a count is 0


def test_resolve_condition_by_id_and_unknown():
    assert Q.resolve_condition("pots")["condition_id"] == "pots"
    assert Q.resolve_condition("")["status"] == UNKNOWN


# ------------------------------------------------------------------------------------------ data


def _need(*names):
    missing = [n for n in names if not table_exists(n)]
    if missing:
        pytest.skip(f"tables not built: {missing}")


@pytest.mark.data
def test_mentions_table_contract():
    _need("measurement_mentions")
    m = read_table("measurement_mentions")
    assert set(PROVENANCE_COLUMNS) <= set(m.columns)
    assert m["mention_id"].is_unique
    assert set(m["role"]) <= set(D.ROLES)
    assert set(m["measurement_id"]) <= set(CLASS_IDS)
    ok = m["object_id"].str.match(r"^(trial:NCT\d{8}|nih:\d+)$")
    assert ok.all(), m.loc[~ok, "object_id"].head()
    assert (m["evidence_flag"] <= m["objective_flag"]).all()
    assert (m["evidence_flag"] == (m["objective_flag"] & m["nonhuman_context_hit"].isna())).all()
    assert (m["data_layer"] == "measurement").all() and (m["evidence_type"] == "text_mined").all()


@pytest.mark.data
def test_registry_dimensions_and_contract():
    _need("measurement_registry")
    r = read_table("measurement_registry")
    assert set(PROVENANCE_COLUMNS) <= set(r.columns)
    assert sorted(r["measurement_id"]) == sorted(CLASS_IDS + BUNDLE_IDS)
    assert (r["object_id"] == "measurement:" + r["measurement_id"]).all()
    for dim in ("measurement_evidence_strength", "technology_maturity", "regulatory_visibility", "deployment_complexity",
                "phenotype_signal_strength"):
        assert dim in r.columns
        assert any(c.startswith(dim + "__") for c in r.columns), f"{dim} has no component columns"
    # filled by phenotype_evidence.py: a verdict, or UNKNOWN (never NaN / 0) with a reason
    assert r["phenotype_signal_strength"].notna().all()
    assert set(r["phenotype_signal_strength"]) <= {"supported", "mixed", "null_result", UNKNOWN}
    assert set(r["phenotype_signal_strength__status"]) <= {"supported", "mixed", "null_result", PE.NOT_OBS, PE.NOT_AN}
    assert set(r["measurement_evidence_strength"]) <= {"high", "moderate", "low", "none"}
    assert set(r["technology_maturity"]) <= {"established", "emerging", "category_only", "no_fda_category"}
    assert set(r["regulatory_visibility"]) <= {"visible", "partial", "not_visible"}
    no_code = r[r["kind"] == "class"].set_index("measurement_id").loc[["capillaroscopy", "posture_detection",
                                                                       "metabolomics", "proteomics"]]
    assert (no_code["regulatory_visibility"] == "not_visible").all()
    ev = "measurement_evidence_strength__"
    assert (r[ev + "n_trials_outcome_measure"] <= r[ev + "n_trials_objective"]).all()
    assert (r[ev + "n_trials_recruiting"] <= r[ev + "n_trials_objective"]).all()


@pytest.mark.data
def test_bundle_counts_are_unions_not_sums():
    _need("measurement_registry")
    r = read_table("measurement_registry").set_index("measurement_id")
    col = "measurement_evidence_strength__n_trials_objective"
    for bid, b in load_config("relevance")["measurement_bundles"].items():
        mem = r.loc[b["members"], col]
        assert mem.max() <= r.loc[bid, col] <= mem.sum()


@pytest.mark.data
def test_condition_evidence_grid_and_ids():
    _need("condition_measurement_evidence")
    e = read_table("condition_measurement_evidence")
    assert len(e) == len(CONDITION_IDS) * (len(CLASS_IDS) + len(BUNDLE_IDS))
    assert e["object_id"].is_unique
    assert (e["object_id"] == "measurement_evidence:" + e["condition_id"] + "|" + e["measurement_id"]).all()
    assert e["provenance_notes"].str.contains("not evidence that the measurement works").all()
    ev = "measurement_evidence_strength__"
    assert (e[ev + "n_trials_objective"] <= e["n_condition_trials_literal"]).all()
    assert (e[ev + "n_trials_objective"] <= e["n_trials_objective_expanded_view"]).all()
    for col in ("example_trial_object_ids", "example_grant_object_ids", "example_fda_object_ids", "molecular_object_ids"):
        for ids in e[col].map(json.loads):
            assert all(re.match(r"^(trial|nih|fda|geo_series):", x) for x in ids)
    lc = e[(e["condition_id"] == "long_covid") & (e["measurement_id"] == "wearable_autonomic_activity_monitoring")].iloc[0]
    assert lc[ev + "n_trials_objective"] > 0 and json.loads(lc["example_trial_object_ids"])


@pytest.mark.data
def test_trial_counts_match_mentions():
    _need("measurement_mentions", "condition_measurement_evidence")
    m = read_table("measurement_mentions")
    tc = read_table("trial_conditions")
    lit = set(tc.loc[(tc["condition_id"] == "me_cfs") & tc["condition_literal_match"], "nct_id"])
    d = m[(m["source"] == "ctgov") & m["evidence_flag"] & m["counted_role"] & (m["measurement_id"] == "hrv")]
    expected = len(set(d["record_id"]) & lit)
    e = read_table("condition_measurement_evidence")
    got = e.loc[(e["condition_id"] == "me_cfs") & (e["measurement_id"] == "hrv"),
                "measurement_evidence_strength__n_trials_objective"].iloc[0]
    assert got == expected


@pytest.mark.data
def test_precision_audit_files_complete():
    for r in (1, 2):
        s = pd.read_csv(PA.sample_path(r))
        j = pd.read_csv(PA.judgement_path(r))
        assert set(s["mention_id"].dropna()) == set(j["mention_id"])
        assert set(j["judgement"]) <= set(PA.JUDGEMENTS)
    s1 = pd.read_csv(PA.sample_path(1))
    assert (s1.groupby("measurement_id")["mention_id"].count() >= 30).sum() >= 20   # >= 30 per class where available
    s2 = pd.read_csv(PA.sample_path(2))
    assert not set(s2["mention_id"].dropna()) & set(s1["mention_id"].dropna())      # round 2 is a fresh sample
    summ = pd.read_csv(PA.SUMMARY_OUT)
    assert summ["reviewer"].str.contains("not a domain expert").all()
    assert set(summ.loc[summ["audit_round"] == 2, "measurement_id"]) == set(D.PATTERN_REVISIONS)


@pytest.mark.data
def test_discover_candidates_ordering_and_phenotype():
    _need("condition_measurement_evidence", "measurement_registry")
    r = Q.discover_candidate_measurements("Long COVID")
    assert r["status"] == "ok" and r["condition"]["condition_id"] == "long_covid"
    keys = [(-c["dimensions"]["measurement_evidence_strength"]["n_trials_objective"],
             -c["dimensions"]["measurement_evidence_strength"]["n_trials_outcome_measure"],
             -c["dimensions"]["measurement_evidence_strength"]["n_nih_core_projects"]) for c in r["candidates"]]
    assert keys == sorted(keys)
    assert [c["rank"] for c in r["candidates"]] == list(range(1, len(r["candidates"]) + 1))
    c0 = r["candidates"][0]
    assert c0["reasons"] and set(c0["dimensions"]) == {"measurement_evidence_strength", "technology_maturity",
                                                       "regulatory_visibility", "deployment_complexity",
                                                       "phenotype_signal_strength"}
    assert {b["bundle_id"] for b in r["bundles"]} == set(BUNDLE_IDS)
    p = Q.discover_candidate_measurements("POTS", phenotype="postural tachycardia")
    assert p["status"] == "ok" and p["phenotype"]["axis_id"] == "postural_tachycardia"
    allowed = {x["measurement_id"] for x in p["phenotype"]["observing_measurements"]}
    assert {c["measurement_id"] for c in p["candidates"]} <= allowed
    assert "blood_biomarkers" not in allowed
    assert Q.discover_candidate_measurements("not a real disease xyz")["status"] == UNKNOWN
    assert Q.discover_candidate_measurements("Long COVID", phenotype="qwerty zxcv")["status"] == UNKNOWN


@pytest.mark.data
def test_get_measurement_evidence_object_ids_and_unknown():
    _need("condition_measurement_evidence", "measurement_mentions")
    e = Q.get_measurement_evidence("wearable autonomic monitoring", "ME/CFS")
    assert e["status"] == "ok" and e["measurement"]["kind"] == "bundle"
    b = next(i for i in e["evidence"] if i["kind"] == "bundle")
    assert b["n_trial_object_ids"] == b["counts"]["n_trials_objective"]
    assert all(x.startswith("trial:") for x in b["trial_object_ids"])
    assert Q.get_measurement_evidence("tarot cards", "ME/CFS")["status"] == UNKNOWN
    none = Q.get_measurement_evidence("retinal imaging", "gastroparesis")
    assert none["status"] == UNKNOWN and "absence of a mention" in none["reason"]


@pytest.mark.data
def test_results_tables_written():
    for name in ("measurement_registry.csv", "measurement_condition_evidence.csv", "measurement_precision_audit.csv",
                 "measurement_precision_summary.csv", "measurement_pattern_changelog.csv",
                 "measurement_pro_filter_summary.csv", "measurement_person_level_datasets.csv"):
        assert (TABLES / name).exists(), name


@pytest.mark.data
def test_precision_status_is_not_claimed_on_tiny_samples():
    _need("measurement_registry")
    r = read_table("measurement_registry").set_index("measurement_id")
    cls = r[r["kind"] == "class"]
    assert set(cls["text_mining_precision_status"]) <= {"meets_0_8", "insufficient_sample", "below_0_8",
                                                         "design_data_only", "not_audited"}
    meets = cls[cls["text_mining_precision_status"] == "meets_0_8"]
    assert (meets["text_mining_precision_n_judged"] >= PA.MIN_JUDGED).all()
    assert (meets["text_mining_precision_n_records"] >= PA.MIN_RECORDS).all()
    assert (meets["text_mining_precision_strict"] >= 0.8).all()
    assert cls.loc["respiratory_rate", "text_mining_precision_status"] == "insufficient_sample"   # 1 fresh mention
    assert cls.loc["continuous_temperature", "text_mining_precision_status"] == "design_data_only"
    assert (cls["text_mining_precision_below_0_8"] == (cls["text_mining_precision_status"] == "below_0_8")).all()
    ok = cls["text_mining_precision_n_judged"] > 0
    assert (cls.loc[ok, "text_mining_precision_ci_low_clustered"] <= cls.loc[ok, "text_mining_precision_ci_low"] + 1e-9).all()


@pytest.mark.data
def test_vascular_imaging_has_no_consumer_clearances():
    _need("measurement_registry")
    r = read_table("measurement_registry").set_index("measurement_id")
    assert r.loc["vascular_imaging", "technology_maturity__n_consumer_device_clearances"] == 0


@pytest.mark.data
def test_fda_decision_totals_do_not_double_count_filtered_subsets():
    _need("measurement_registry", "measurement_regulatory_status")
    r = read_table("measurement_registry").set_index("measurement_id")
    st = read_table("measurement_regulatory_status")
    mem = load_config("relevance")["measurement_bundles"]["wearable_autonomic_activity_monitoring"]["members"]
    u = st[st["measurement_id"].isin(mem) & st["product_code"].notna()].copy()
    u["f"] = u["device_filter"].fillna("")
    u = u.drop_duplicates(["product_code", "f"])
    unfilt = set(u.loc[u["f"] == "", "product_code"])
    u = u[(u["f"] == "") | ~u["product_code"].isin(unfilt)]
    assert r.loc["wearable_autonomic_activity_monitoring", "fda_n_510k"] == int(pd.to_numeric(u["n_510k"]).sum())


@pytest.mark.data
def test_specific_term_and_hrv_sensitivity_columns():
    _need("condition_measurement_evidence")
    e = read_table("condition_measurement_evidence")
    ev = "measurement_evidence_strength__"
    assert (e["n_trials_objective_specific_terms"] <= e[ev + "n_trials_objective"]).all()
    assert (e["n_nih_core_projects_specific_terms"] <= e[ev + "n_nih_core_projects"]).all()
    plain = e[e["broad_literal_terms_excluded"] == "[]"]
    assert (plain["n_trials_objective_specific_terms"] == plain[ev + "n_trials_objective"]).all()
    assert (plain["n_nih_core_projects_specific_terms"] == plain[ev + "n_nih_core_projects"]).all()
    d = e[(e["condition_id"] == "dysautonomia") & (e["measurement_id"] == "hrv")].iloc[0]
    assert d["n_trials_objective_specific_terms"] < d[ev + "n_trials_objective"]
    b = e[e["kind"] == "bundle"]
    assert (b["n_trials_objective_excluding_hrv"] <= b[ev + "n_trials_objective"]).all()
    assert e.loc[e["kind"] == "class", "n_trials_objective_excluding_hrv"].isna().all()


@pytest.mark.data
def test_phenotype_resolution_refuses_ambiguous_words():
    _need("phenotype_axes")
    ax = read_table("phenotype_axes")
    for word in ("orthostatic", "intolerance", "sleep"):
        r = Q.resolve_phenotype(word, ax)
        assert r["status"] == UNKNOWN and len(r["candidates"]) > 1, word
    assert Q.resolve_phenotype("post-exertional malaise", ax)["axis_id"] == "post_exertional_malaise"
    assert Q.resolve_phenotype("", ax)["status"] == UNKNOWN


@pytest.mark.data
def test_capillaroscopy_query_is_the_class_not_the_bundle():
    _need("condition_measurement_evidence", "measurement_mentions")
    e = Q.get_measurement_evidence("capillaroscopy", "migraine")
    assert e["status"] == UNKNOWN          # the bundle's microvascular_function trials must not stand in for it
    lc = Q.get_measurement_evidence("capillaroscopy", "long_covid")
    assert lc["status"] == "ok" and [i["measurement_id"] for i in lc["evidence"]] == ["capillaroscopy"]


# ------------------------------------------------------------------------------------------ Phase 3: unit (no data)
@pytest.mark.parametrize("n,model,under,expected", [
    (3, "positive", False, "supported"), (3, "not_tested", False, "supported"), (3, "null_result", False, "mixed"),
    (3, "negative", False, "mixed"), (3, "positive", True, "mixed"), (3, "not_tested", True, "mixed"),
    (0, "positive", False, "mixed"), (0, "null_result", False, "null_result"), (0, "not_tested", True, "null_result"),
])
def test_verdict_rule(n, model, under, expected):
    assert PE.verdict(n, model, under) == expected


def test_model_status_rules():
    assert PE.model_status_delta(0.01, 0.03, 0.005) == "positive"
    assert PE.model_status_delta(0.01, 0.03, 0.2) == PE.NULL          # CI > 0 but the negative control fails
    assert PE.model_status_delta(0.01, 0.03, float("nan")) == PE.NULL  # no permutation p -> not positive
    assert PE.model_status_delta(-0.03, -0.01, 0.5) == "negative"
    assert PE.model_status_delta(-0.01, 0.03, 0.01) == PE.NULL
    assert PE.model_status_delta(float("nan"), float("nan"), 0.01) == "not_tested"
    assert PE.model_status_auroc(0.79, 0.92, 0.001) == "positive"
    assert PE.model_status_auroc(0.44, 0.63, 0.27) == PE.NULL
    assert PE.NULL != "null"      # the bare string 'null' is read back as missing by pandas' CSV reader


def test_summarise_features_effect_floor_and_like_units():
    rows = pd.DataFrame({"feature": ["a", "b", "c", "d"],
                         "effect_measure": ["hedges_g", "hedges_g", "paired_difference_fraction", "hedges_g"],
                         "effect_size": [0.15, -0.35, 0.9, 0.5], "ci_low": [0.1, -0.5, 0.1, -0.1],
                         "ci_high": [0.2, -0.2, 1.0, 1.1], "q": [0.01, 0.02, 0.2, 0.3], "object_id": list("abcd")})
    s = PE.summarise_features(rows)
    assert s["n_features_fdr"] == 2 and s["n_features_fdr_meaningful"] == 1       # 0.15 is below the 0.2 floor
    assert s["best_feature"] == "b" and s["best_effect_passes_fdr"]
    s2 = PE.summarise_features(rows.assign(q=[0.5, 0.5, 0.5, 0.5]))
    assert s2["n_features_fdr"] == 0 and s2["best_feature"] == "d"               # SD-unit effects compared first


def test_demo_phenotype_resolution_is_exact_only():
    assert PE.resolve_demo_phenotype("Orthostatic/autonomic dysfunction") == "orthostatic_autonomic"
    assert PE.resolve_demo_phenotype("sleep_circadian") == "sleep_circadian"
    assert PE.resolve_demo_phenotype("Post-infection physiological recovery") == "post_infection_recovery"
    assert PE.resolve_demo_phenotype("activity intolerance") == "activity_intolerance_pem"
    for word in ("sleep", "orthostatic", "", "tarot"):
        assert PE.resolve_demo_phenotype(word) is None
    assert set(PE.PHENOTYPES) == {"orthostatic_autonomic", "activity_intolerance_pem", "sleep_circadian",
                                  "post_infection_recovery"}


def test_nhanes_feature_families_partition_the_26_features():
    assert len(PE.NHANES_FAMILY) == 26
    sleep = [f for f, fam in PE.NHANES_FAMILY.items() if fam == "sleep_proxy"]
    assert len(sleep) == 6 and not {f for f, fam in PE.NHANES_FAMILY.items() if fam in PE.ACC_FAMILIES} & set(sleep)


# ------------------------------------------------------------------------------------------ Phase 3: data
@pytest.mark.data
def test_signal_results_contract_and_sanity_checks():
    _need("measurement_phenotype_signal")
    s = read_table("measurement_phenotype_signal")
    assert set(PROVENANCE_COLUMNS) <= set(s.columns) and s["object_id"].is_unique
    assert s["object_id"].str.match(r"^measurement_signal:[^|]+\|[^|]+\|[a-z_]+\|[a-z_]+$").all()
    assert set(s["verdict"]) <= set(PE.VERDICTS)
    assert (s["n_features_fdr_meaningful"] <= s["n_features_fdr"]).all()
    assert (s["n_features_fdr"] <= s["n_features_tested"]).all()
    for ids in s["signature_object_ids"].map(json.loads):
        assert ids and all(x.startswith("signature:") for x in ids)
    # rule sanity checks fixed in the plan
    chk = pd.read_csv(TABLES / "phase3_rule_sanity_checks.csv")
    assert chk["passed"].all(), chk
    # an acute COVID-19 cohort is never a condition label
    ac = s[s["dataset_id"] == PE.ACUTE]
    assert len(ac) and not ac["is_condition_label"].any() and ac["condition_id"].isna().all()
    # NHANES has no heart rate: no NHANES row for any HR-derived class
    assert set(s.loc[s["dataset_id"] == PE.NHANES, "measurement_id"]) <= {"accelerometry", "sleep_objective"}


@pytest.mark.data
def test_condition_fill_guardrails():
    _need("condition_measurement_evidence", "measurement_phenotype_signal")
    e = read_table("condition_measurement_evidence")
    ps = "phenotype_signal_strength"
    assert e[ps].notna().all()
    unk = e[e[ps] == UNKNOWN]
    assert (unk[ps + "__status"].isin([PE.NOT_OBS, PE.NOT_AN])).all()
    assert (unk[ps + "__reason"].str.len() > 10).all()                      # every UNKNOWN says why
    assert unk[ps + "__n_features_fdr"].isna().all()                         # UNKNOWN, not zero
    lc = e[e["condition_id"] == "long_covid"].set_index("measurement_id")
    assert lc.loc["wearable_heart_rate", ps + "__label_id"] == PE.LC_LABEL
    assert lc.loc["accelerometry", ps] == UNKNOWN                            # Uwakwe released no steps
    adj = json.loads(lc.loc["accelerometry", ps + "__adjacent_evidence"])
    assert adj and all("NOT a Long COVID label" in a["label"] for a in adj)  # acute cohort only as adjacent evidence
    assert not (e[ps + "__dataset_id"] == PE.ACUTE).any()                    # never the verdict of any condition
    me = e[e["condition_id"] == "me_cfs"].set_index("measurement_id")
    assert me.loc["accelerometry", ps + "__label_class"] == "constructed_proxy"
    assert json.loads(me.loc["accelerometry", ps + "__sensitivity_labels"])
    assert me.loc["proteomics", ps + "__status"] == PE.NOT_AN
    for cid in ("pots", "dysautonomia", "eds_hsd", "mcas"):                  # no public person-level label at all
        assert (e.loc[e["condition_id"] == cid, ps] == UNKNOWN).all()
    b = e[e["kind"] == "bundle"]
    for _, r in b.iterrows():                                                # bundle value = highest member verdict
        mem = [m["phenotype_signal_strength"] for m in json.loads(r[ps + "__by_member"])]
        best = PE.best_verdict(mem)
        assert r[ps] == (best if best else UNKNOWN)


@pytest.mark.data
def test_molecular_context_is_a_separate_group():
    _need("condition_measurement_evidence", "measurable_biology")
    e = read_table("condition_measurement_evidence")
    mc = "molecular_context__best_support_status"
    assert set(e[mc]) <= set(PE.MOL_RANK) | {"no_system_link", UNKNOWN}
    acc = e[(e["measurement_id"] == "accelerometry") & (e["kind"] == "class")]
    assert (acc[mc] == "no_system_link").all()
    mb = read_table("measurable_biology")
    row = e[(e["condition_id"] == "pots") & (e["measurement_id"] == "hrv")].iloc[0]
    expected = PE._best_mol(mb.loc[(mb["condition_id"] == "pots") & (mb["measurement_class"] == "hrv"), "support_status"])
    assert row[mc] == expected
    ids = json.loads(row["molecular_context__object_ids"])
    assert ids and all(x.startswith("molbio:") for x in ids)
    # never blended into another dimension: dimension values do not depend on it
    assert not any(c.startswith("molecular_context") for c in ("phenotype_signal_strength", "measurement_evidence_strength"))


@pytest.mark.data
def test_phase3_table_shape_and_dimensions():
    _need("phase3_measurement_evidence")
    p = read_table("phase3_measurement_evidence")
    assert set(PROVENANCE_COLUMNS) <= set(p.columns) and p["object_id"].is_unique
    assert set(p["phenotype_id"]) == set(PE.PHENOTYPES)
    for dim in ("phenotype_signal_strength", "measurement_evidence_strength", "technology_maturity",
                "regulatory_visibility", "deployment_complexity"):
        assert dim in p.columns
    assert not any("composite" in c or "overall_score" in c for c in p.columns)
    for pre in ("ctgov__", "nih__", "fda__", "molecular_context__"):                 # evidence per source
        assert any(c.startswith(pre) for c in p.columns), pre
    assert set(p["phenotype_signal_strength"]) <= set(PE.VERDICTS) | {UNKNOWN}
    o = p[p["phenotype_id"] == "orthostatic_autonomic"].set_index("measurement_id")
    assert o.loc["hrv", "phenotype_signal_strength"] == UNKNOWN                     # no public beat-to-beat HRV
    assert o.loc["posture_detection", "phenotype_signal_strength"] == UNKNOWN
    pi = p[p["phenotype_id"] == "post_infection_recovery"].set_index("measurement_id")
    assert pi.loc["accelerometry", "phenotype_signal_strength__label_class"] == "acute_infection_episode"
    reg = read_table("measurement_registry").set_index("measurement_id")
    for r in p.itertuples():                                                         # registry dimensions carried over
        assert r.technology_maturity == reg.loc[r.measurement_id, "technology_maturity"]
        assert r.regulatory_visibility == reg.loc[r.measurement_id, "regulatory_visibility"]
    for f in ("phase3_measurement_evidence.csv", "phase3_signal_results.csv", "phase3_condition_signal.csv",
              "phase3_phenotype_summary.csv", "phase3_verdict_sensitivity.csv"):
        assert (TABLES / f).exists(), f
    back = pd.read_csv(TABLES / "phase3_signal_results.csv")
    assert back["verdict"].notna().all()                    # no verdict string is read back as missing


@pytest.mark.data
def test_query_returns_phenotype_evidence():
    _need("condition_measurement_evidence", "phase3_measurement_evidence", "measurement_phenotype_signal")
    r = Q.discover_candidate_measurements("Long COVID")
    hr = next(c for c in r["candidates"] if c["measurement_id"] == "wearable_heart_rate")
    ps = hr["dimensions"]["phenotype_signal_strength"]
    assert ps["value"] == PE.NULL and ps["label_class"] == "self_reported_condition_case_definition"
    assert "molecular_context" in hr and "molecular_context" not in hr["dimensions"]
    acc = next(c for c in r["candidates"] if c["measurement_id"] == "accelerometry")
    assert acc["dimensions"]["phenotype_signal_strength"]["value"] == UNKNOWN
    assert acc["dimensions"]["phenotype_signal_strength"]["reason"]
    d = Q.discover_candidate_measurements("ME/CFS", phenotype="activity intolerance / post-exertional pattern")
    assert d["status"] == "ok" and d["phenotype"]["kind"] == "demo_phenotype"
    assert {c["measurement_id"] for c in d["candidates"]} <= {x["measurement_id"] for x in d["phenotype"]["observing_measurements"]}
    assert all("phenotype_evidence" in c for c in d["candidates"])
    # an exact axis still resolves to the axis, not to a demo phenotype
    assert Q.discover_candidate_measurements("POTS", phenotype="postural tachycardia")["phenotype"]["axis_id"] == "postural_tachycardia"
    e = Q.get_measurement_evidence("accelerometry", "long_covid")
    item = e["evidence"][0]
    roles = {x["role"] for x in item["phenotype_signal_results"]}
    assert roles and all("adjacent" in x for x in roles)                    # acute cohort only as adjacent evidence
    m = Q.get_measurement_evidence("accelerometry", "ME/CFS")["evidence"][0]
    assert any(x["role"] == "reference label" for x in m["phenotype_signal_results"])
    ph = Q.get_phenotype_measurement_evidence("sleep/circadian disruption")
    assert ph["status"] == "ok" and ph["phenotype"]["phenotype_id"] == "sleep_circadian"
    assert set(ph["technologies"][0]["dimensions"]) == {"phenotype_signal_strength", "measurement_evidence_strength",
                                                        "technology_maturity", "regulatory_visibility",
                                                        "deployment_complexity"}
    assert set(ph["technologies"][0]["evidence_by_source"]) == {"public_patient_datasets", "clinicaltrials_gov",
                                                                "nih_reporter", "openfda", "molecular_context"}
    assert Q.get_phenotype_measurement_evidence("postural tachycardia")["phenotype"]["phenotype_id"] == "orthostatic_autonomic"
    assert Q.get_phenotype_measurement_evidence("tarot")["status"] == UNKNOWN


# ------------------------------------------------------------------------------------------ Phase 3: reviewer additions
def test_scope_flags_are_read_from_scope_notes():
    assert "ACUTE" in PE._scope_flags(PE.PHENOTYPES["post_infection_recovery"]["scope_note"]["accelerometry"])
    assert "DEFINITIONAL" in PE._scope_flags(PE.PHENOTYPES["activity_intolerance_pem"]["scope_note"]["accelerometry"])
    assert "drug-defined" in PE._scope_flags(PE.PHENOTYPES["sleep_circadian"]["scope_note"]["sleep_objective"])
    assert PE._scope_flags(None) == "" and PE._scope_flags("") == ""


@pytest.mark.data
def test_phase3_caveats_travel_with_the_rows():
    _need("phase3_measurement_evidence", "measurement_phenotype_signal")
    p = read_table("phase3_measurement_evidence")
    assert not p.apply(lambda c: c.astype(str).str.contains("__none__")).any().any()   # no internal placeholder leaks
    pi = p[(p["phenotype_id"] == "post_infection_recovery") & (p["kind"] == "class")].set_index("measurement_id")
    for m in ("accelerometry", "wearable_heart_rate"):                                 # acute phase, not recovery
        assert "ACUTE PHASE ONLY" in pi.loc[m, "phenotype_signal_strength__scope_note"]
    s = read_table("measurement_phenotype_signal")
    nh = s[s["dataset_id"] == PE.NHANES]
    assert nh.loc[nh["label_class"] == "rx_reason_code", "notes"].str.contains("TAKING a prescription").all()
    assert nh.loc[nh["label_id"].isin(PE.NHANES_DEFINITIONAL), "notes"].str.contains("PARTLY DEFINITIONAL").all()
    assert not nh.loc[nh["label_id"] == "mortality", "notes"].str.contains("DEFINITIONAL").any()
    txt = PE.REPORT.read_text()
    assert "asymmetric" in txt and "winner's curse" in txt and "do not validate it independently" in txt


@pytest.mark.data
def test_demo_phenotype_evidence_is_not_passed_off_as_a_condition_result():
    _need("condition_measurement_evidence", "phase3_measurement_evidence", "measurement_phenotype_signal")
    d = Q.discover_candidate_measurements("ME/CFS", phenotype="post-infection physiological recovery")
    acc = next(c for c in d["candidates"] if c["measurement_id"] == "accelerometry")
    pe = acc["phenotype_evidence"]
    assert pe["dimensions"]["phenotype_signal_strength"]["value"] == "supported"      # the acute COVID-19 cohort
    assert pe["applies_to_queried_condition"] is False and "NOT a result for" in pe["scope"]
    assert acc["dimensions"]["phenotype_signal_strength"]["value"] == "mixed"          # the ME/CFS-level value
    lc = Q.discover_candidate_measurements("Long COVID", phenotype="sleep/circadian disruption")
    sl = next(c for c in lc["candidates"] if c["measurement_id"] == "sleep_objective")["phenotype_evidence"]
    assert sl["applies_to_queried_condition"] is False and sl["reference_label_condition_id"] is None  # rx_insomnia
    ai = Q.discover_candidate_measurements("ME/CFS", phenotype="activity intolerance")
    a2 = next(c for c in ai["candidates"] if c["measurement_id"] == "accelerometry")["phenotype_evidence"]
    assert a2["applies_to_queried_condition"] is True and a2["reference_label_condition_id"] == "me_cfs"
    unk = next(c for c in ai["candidates"] if c["measurement_id"] == "cpet")["phenotype_evidence"]
    assert unk["applies_to_queried_condition"] is None
    assert any("applies_to_queried_condition" in c for c in ai["caveats"])
