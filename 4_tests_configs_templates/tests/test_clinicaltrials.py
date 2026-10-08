"""Tests for measure_it.ingestion.clinicaltrials.

Unit tests exercise the pure helpers on synthetic records. Tests marked `data`
check the processed outputs written by `uv run python -m measure_it.ingestion.clinicaltrials`.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from measure_it.config import RAW, UNKNOWN, load_config
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import processed_path, read_table
from measure_it.ingestion import clinicaltrials as ct

TABLES = ["clinical_trials", "trial_conditions", "trial_interventions", "trial_outcomes", "trial_sites",
          "ctgov_facility_summary"]


# ------------------------------------------------------------------ unit tests

@pytest.mark.parametrize("term,expected", [
    ("PASC", True), ("ME/CFS", True), ("hEDS", True), ("POTS", True), ("MCAS", True),
    ("long COVID", False), ("Ehlers-Danlos", False), ("fibromyalgia", False), ("migraine", False),
    ("post-COVID-19 condition", False),
])
def test_is_bare_acronym(term, expected):
    assert ct.is_bare_acronym(term) is expected


def test_query_string_is_phrase_quoted():
    assert ct.query_string("long COVID") == '"long COVID"'
    assert ct.query_string('"migraine"') == '"migraine"'


def test_search_plan_never_queries_bare_acronyms():
    plan = ct.search_plan()
    conds = load_config("conditions")["conditions"]
    assert len(plan) == sum(len(c.get("search_terms", [])) for c in conds)
    assert set(plan.loc[plan["skipped"], "term"]) >= {"PASC", "ME/CFS"}
    assert plan.loc[~plan["skipped"], "query_string"].str.match(r'^".+"$').all()
    assert plan.loc[plan["skipped"], "query_string"].isna().all()


def test_literal_term_matching():
    assert ct.literal_term_in("post COVID syndrome", "Post-COVID Syndrome")
    assert ct.literal_term_in("migraine", "Chronic Migraines")
    assert ct.literal_term_in("Ehlers-Danlos", "Ehlers Danlos Syndrome, Hypermobility Type")
    assert not ct.literal_term_in("pots", "teapots")
    assert not ct.literal_term_in("long COVID", None)
    assert not ct.literal_term_in("long COVID", float("nan"))


def test_literal_matching_folds_covid19_and_sars_cov2_spellings():
    # regression: "Post-COVID-19 Syndrome" (33 long_covid trials) was flagged non-literal for "post COVID syndrome"
    assert ct.literal_term_in("post COVID syndrome", "Post-COVID-19 Syndrome")
    assert ct.literal_term_in("long COVID", "Long COVID-19")
    assert ct.literal_term_in("long COVID", "Long-COVID19 fatigue")
    assert ct.literal_term_in("post-COVID-19 condition", "Post COVID Condition")
    assert ct.literal_term_in("post-acute sequelae of SARS-CoV-2", "Post-Acute Sequelae of SARS-CoV2 Infection")
    # folding must not create matches across different words
    assert not ct.literal_term_in("long COVID", "long-term COVID-19 vaccine follow-up")
    assert not ct.literal_term_in("long COVID", "COVID-19 190 day follow-up")


def test_condition_literal_match_uses_all_condition_terms():
    # regression: a trial retrieved by "long COVID" (synonym expansion) whose conditions list
    # "Post-COVID Conditions" was counted as non-literal for long_covid
    trial = {"conditions": "Post-COVID Conditions | Obesity", "keywords": None, "brief_title": "Adiposity study",
             "official_title": None, "acronym": None, "condition_mesh_terms": None, "condition_mesh_ancestors": None}
    assert not (ct.LITERAL_FIELDS & set(ct.match_fields("long COVID", trial)))
    assert ct.condition_literal_match(["long COVID", "post-COVID condition"], trial)
    # only a MeSH-ancestor hit is still not literal
    t2 = dict(trial, conditions="Multiple System Atrophy", condition_mesh_ancestors="Autonomic Nervous System Diseases")
    assert not ct.condition_literal_match(["dysautonomia", "autonomic nervous system disorder"], t2)


def test_normalize_text_and_site_key_handle_nan():
    # regression: NaN facility normalized to "nan", merging unnamed sites into a pseudo-facility "nan|city|..."
    assert ct.normalize_text(float("nan")) == ""
    assert ct.normalize_text(pd.NA) == ""
    assert ct.site_key(float("nan"), "Phoenix", "Arizona", "United States") == "|phoenix|arizona|united states"
    assert not ct.is_generic_facility(float("nan"))


@pytest.mark.parametrize("name,generic", [
    ("Research Site", True), ("Pfizer Investigational Site", True), ("Local Institution - 0012", True),
    ("Investigational Site Number : 8400001", True), ("Site 101", True), ("Novartis Investigative Site", True),
    ("For additional information regarding investigative sites for this trial, contact 1-877-CTLILLY", True),
    ("Johns Hopkins University Clinical Research Site", False), ("Duke Clinical Research Site", False),
    ("Mayo Clinic in Rochester", False), ("Stanford University", False), (None, False),
])
def test_generic_facility(name, generic):
    assert ct.is_generic_facility(name) is generic


@pytest.mark.parametrize("name,generic", [
    # facilities-review placeholders (2026-09-23)
    ("Site Reference ID/Investigator# 56266", True), ("SIte reference ID 304", True),
    ("Contact Medtronic for Exact Locations", True), ("Contact Medtronic for exact location", True),
    ("IQVIA Virtual Site", True), ("Remote trial - anyone residing in the United States", True),
    ("Online Patient Enrollment System", True), ("Nationwide home based study program", True),
    ("Alethios Digital Research Platform (Decentralized)", True), ("Participants' homes", True),
    # real facilities that merely contain a keyword stay identifiable
    ("Remote Monitoring Clinic", False), ("Virginia Mason Medical Center", False), ("Contact Lens Institute", False),
])
def test_generic_facility_review_placeholders(name, generic):
    assert ct.is_generic_facility(name) is generic


def test_generic_facility_city_only_names():
    assert ct.is_generic_facility("Denver", "Denver", "Colorado")
    assert ct.is_generic_facility("Denver, Colorado", "Denver", "Colorado")
    assert not ct.is_generic_facility("Denver Health", "Denver", "Colorado")
    assert not ct.is_generic_facility("Denver", "Aurora", "Colorado")
    assert not ct.is_generic_facility("Denver")        # no city given: regex only


def test_date_year_and_lists():
    assert ct.date_year("2021-03") == 2021
    assert ct.date_year(None) is None
    assert ct.join_list(["a", None, " ", "b"]) == "a | b"
    assert ct.split_list("a | b") == ["a", "b"]
    assert ct.split_list(None) == []


SYNTHETIC = {
    "protocolSection": {
        "identificationModule": {"nctId": "NCT00000001", "briefTitle": "Wearable HRV in Long COVID",
                                 "officialTitle": "A Study of Heart Rate Variability in Post-COVID Syndrome"},
        "statusModule": {"overallStatus": "RECRUITING", "startDateStruct": {"date": "2023-05", "type": "ACTUAL"},
                         "lastUpdatePostDateStruct": {"date": "2025-01-02"}},
        "sponsorCollaboratorsModule": {"leadSponsor": {"name": "Uni A", "class": "OTHER"},
                                       "collaborators": [{"name": "NIH", "class": "NIH"}]},
        "descriptionModule": {"briefSummary": "We measure HRV."},
        "conditionsModule": {"conditions": ["Long COVID"], "keywords": ["dysautonomia"]},
        "designModule": {"studyType": "INTERVENTIONAL", "phases": ["NA"],
                         "designInfo": {"allocation": "RANDOMIZED"}, "enrollmentInfo": {"count": 40, "type": "ESTIMATED"}},
        "armsInterventionsModule": {
            "armGroups": [{"label": "A", "type": "EXPERIMENTAL", "interventionNames": ["Device: Watch"]}],
            "interventions": [{"type": "DEVICE", "name": "Watch", "description": "smartwatch PPG",
                               "otherNames": ["Apple Watch"], "armGroupLabels": ["A"]}]},
        "outcomesModule": {"primaryOutcomes": [{"measure": "RMSSD", "timeFrame": "12 weeks"}],
                           "secondaryOutcomes": [{"measure": "Fatigue Severity Scale", "description": "PRO",
                                                  "timeFrame": "12 weeks"}]},
        "eligibilityModule": {"eligibilityCriteria": "Inclusion: adults", "sex": "ALL", "minimumAge": "18 Years"},
        "contactsLocationsModule": {"locations": [
            {"facility": "Uni A Hospital", "city": "Boston", "state": "Massachusetts", "zip": "02114",
             "country": "United States", "status": "RECRUITING", "geoPoint": {"lat": 42.36, "lon": -71.06},
             "contacts": [{"name": "Dr X", "email": "x@example.org"}]},
            {"facility": "Research Site", "city": "Oslo", "country": "Norway"}]},
    },
    "derivedSection": {"conditionBrowseModule": {"meshes": [{"id": "D000094024", "term": "Post-Acute COVID-19 Syndrome"}],
                                                 "ancestors": [{"id": "D018352", "term": "Coronavirus Infections"}]},
                       "miscInfoModule": {"versionHolder": "2026-09-23"}},
    "hasResults": False,
}


def test_parse_study_flattens_all_requested_fields():
    p = ct.parse_study(SYNTHETIC)
    t = p["trial"]
    assert t["nct_id"] == "NCT00000001"
    assert t["overall_status"] == "RECRUITING" and t["start_year"] == 2023
    assert t["phases"] == "NA" and t["allocation"] == "RANDOMIZED" and t["enrollment_count"] == 40
    assert t["lead_sponsor_class"] == "OTHER" and t["collaborators"] == "NIH"
    assert t["condition_mesh_terms"] == "Post-Acute COVID-19 Syndrome"
    assert t["n_locations"] == 2 and t["n_us_locations"] == 1
    assert json.loads(t["arm_groups_json"])[0]["label"] == "A"
    assert p["interventions"][0]["other_names"] == "Apple Watch"
    assert {o["outcome_type"] for o in p["outcomes"]} == {"primary", "secondary"}
    assert p["sites"][0]["geopoint_lat"] == 42.36 and p["sites"][1]["geopoint_lat"] is None
    # contacts (named people, emails) must never be carried into the tables
    assert all("contacts" not in s and "email" not in json.dumps(s) for s in p["sites"])


def test_match_fields_literal_vs_expansion():
    t = ct.parse_study(SYNTHETIC)["trial"]
    assert ct.match_fields("long COVID", t)[:2] == ["condition", "brief_title"]
    assert ct.match_fields("post COVID syndrome", t) == ["official_title"]
    # only reachable through MeSH ancestor -> not a literal match
    hits = ct.match_fields("coronavirus infections", t)
    assert hits == ["condition_mesh_ancestor"] and not (ct.LITERAL_FIELDS & set(hits))


def test_measurement_regex_uses_class_excludes():
    inc, exc, kind = ct.measurement_regex("sleep_objective")
    assert kind.startswith("measurement class")
    text = "Pittsburgh Sleep Quality Index (PSQI) at baseline"
    assert inc.search(exc.sub(" ", text)) is None       # questionnaire is not objective sleep measurement
    assert inc.search("overnight polysomnography")
    inc2, exc2, kind2 = ct.measurement_regex(r"tilt[- ]table")
    assert kind2 == "regex" and exc2 is None and inc2.search("Tilt-Table test")


# ------------------------------------------------------------------ data tests

def _needs_outputs():
    if not all(processed_path(t).exists() for t in TABLES):
        pytest.skip("run `uv run python -m measure_it.ingestion.clinicaltrials` first")


@pytest.fixture(scope="module")
def tables():
    _needs_outputs()
    return {t: read_table(t) for t in TABLES}


@pytest.mark.data
def test_provenance_and_layer(tables):
    for name, df in tables.items():
        assert set(PROVENANCE_COLUMNS) <= set(df.columns), name
        assert (df["data_layer"] == "facility").all(), name
        assert (df["evidence_type"] == "trial_registry_record").all(), name
        assert df["source_version"].str.contains("dataTimestamp").all(), name
        assert df["provenance_notes"].str.contains("not evidence").all(), name


@pytest.mark.data
def test_clinical_trials_unique_and_linked(tables):
    trials, tc = tables["clinical_trials"], tables["trial_conditions"]
    assert trials["nct_id"].is_unique
    assert trials["nct_id"].str.match(r"^NCT\d{8}$").all()
    assert set(tc["nct_id"]) == set(trials["nct_id"])
    for child in ("trial_interventions", "trial_outcomes", "trial_sites"):
        assert set(tables[child]["nct_id"]) <= set(trials["nct_id"]), child
    assert not tc.duplicated(["nct_id", "condition_id", "matched_term"]).any()


@pytest.mark.data
def test_query_log_matches_trial_conditions(tables):
    qlog = json.loads((RAW / ct.SOURCE_ID / "query_log.json").read_text())
    tc = tables["trial_conditions"]
    ran = [q for q in qlog["queries"] if not q["skipped"]]
    skipped = [q for q in qlog["queries"] if q["skipped"]]
    assert {q["term"] for q in skipped} >= {"PASC", "ME/CFS"}
    assert not set(tc["matched_term"]) & {q["term"] for q in skipped}
    for q in ran:
        assert q["n_unique_ids"] == q["total_count"], q["query_string"]
        got = tc[(tc["condition_id"] == q["condition_id"]) & (tc["matched_term"] == q["term"])]
        assert got["nct_id"].nunique() == q["n_unique_ids"], q["query_string"]
        assert (got["query_string"] == q["query_string"]).all()
    cov = qlog["record_coverage"]
    assert cov["requested"] == cov["returned"] and not cov["missing"]


@pytest.mark.data
def test_every_condition_searched(tables):
    ids = {c["id"] for c in load_config("conditions")["conditions"]}
    assert set(tables["trial_conditions"]["condition_id"]) <= ids
    for cid in ("long_covid", "me_cfs", "pots", "dysautonomia", "fibromyalgia", "migraine"):
        assert (tables["trial_conditions"]["condition_id"] == cid).any(), cid


@pytest.mark.data
def test_trials_cover_all_statuses_and_years(tables):
    t = tables["clinical_trials"]
    assert {"COMPLETED", "RECRUITING"} <= set(t["overall_status"])
    assert t["overall_status"].nunique() >= 5
    assert t["start_year"].min() < 2005 and t["start_year"].max() >= 2025


@pytest.mark.data
def test_outcomes_and_interventions_shape(tables):
    oc, iv = tables["trial_outcomes"], tables["trial_interventions"]
    assert set(oc["outcome_type"]) <= {"primary", "secondary", "other"}
    assert oc["measure"].notna().mean() > 0.99
    assert iv["intervention_type"].notna().mean() > 0.99
    assert not oc.duplicated(["nct_id", "outcome_type", "outcome_idx"]).any()


@pytest.mark.data
def test_sites_geography(tables):
    s = tables["trial_sites"]
    assert not s.duplicated(["nct_id", "site_idx"]).any()
    assert set(s["source_geographic_resolution"]) <= {"city_centroid", "zcta_centroid", "state", "none"}
    assert set(s["country_group"]) <= {"US", "US_territory", "non_US", "unknown_country"}
    pt = s["geocode_method"] == ct.GEOPOINT_METHOD
    assert (s.loc[pt, "source_geographic_resolution"] == "city_centroid").all()  # geoPoint = city centroid
    assert s.loc[pt, ["lat", "lon"]].notna().all().all()
    us = s[s["country_group"] == "US"]
    assert us["county_fips"].notna().mean() > 0.95
    cf = us["county_fips"].dropna()
    assert cf.str.match(r"^\d{5}$").all()
    counties = set(read_table("geographies").query("geo_level == 'county'")["geo_id"])
    assert set(cf) <= counties
    assert set(s["county_assignment_method"]) <= {"zip5_as_zcta_internal_point", "geopoint_county_majority_of_zcta",
                                                   "geopoint_point_in_polygon_cb_2024_5m", "none"}
    maj = s["county_assignment_method"] == "geopoint_county_majority_of_zcta"
    assert (s.loc[maj, "pip_county_share_of_zcta"] >= 0.5).all()
    assert (s.loc[maj, "county_fips"] == s.loc[maj, "county_fips_pip"]).all()
    zip_first = s["county_assignment_method"] == "zip5_as_zcta_internal_point"
    assert (s.loc[zip_first, "county_fips"] == s.loc[zip_first, "county_fips_zip"]).all()
    pip = s["county_assignment_method"].isin(["geopoint_point_in_polygon_cb_2024_5m", "geopoint_county_majority_of_zcta"])
    assert (s.loc[pip, "county_fips"] == s.loc[pip, "county_fips_pip"]).all()
    # non-U.S. sites are kept but never given a U.S. county
    assert s.loc[s["country_group"] == "non_US", "county_fips"].isna().all()
    assert (s["country_group"] == "non_US").any()
    # U.S. coordinates fall inside a plausible U.S. bounding box
    ll = us.dropna(subset=["lat"])
    assert ll["lat"].between(17, 72).all() and ll["lon"].between(-180, -60).all()


@pytest.mark.data
def test_site_history_counts(tables):
    h, s = tables["ctgov_facility_summary"], tables["trial_sites"]
    assert h["site_key"].is_unique
    assert h["n_trials"].sum() >= s["nct_id"].nunique()
    assert (h["n_trials_site_recruiting"] <= h["n_trials"]).all()
    assert (h["n_device_or_diagnostic_trials"] <= h["n_trials"]).all()


@pytest.mark.data
def test_unnamed_sites_are_flagged_not_facilities(tables):
    h, s = tables["ctgov_facility_summary"], tables["trial_sites"]
    # no pseudo-facility keyed on the string "nan"
    assert not h["site_key"].str.startswith("nan|").any()
    assert (s["facility_missing"] == s["facility"].isna() | (s["facility"].fillna("").str.strip() == "")).all()
    assert (s["facility_identifiable"] == ~(s["facility_generic"] | s["facility_missing"])).all()
    unnamed = h[h["facility"].isna()]
    assert len(unnamed) and unnamed["facility_missing"].all() and not unnamed["facility_identifiable"].any()
    assert not h.loc[h["facility_generic"], "facility_identifiable"].any()
    sites = ct.trial_sites_for_condition("pots", exclude_generic=True)
    assert sites["facility_identifiable"].all() and sites["facility"].notna().all()


@pytest.mark.data
def test_condition_literal_match_column(tables):
    tc, t = tables["trial_conditions"], tables["clinical_trials"]
    # condition-level literal is implied by any per-term literal of the same (trial, condition)
    per = tc.groupby(["nct_id", "condition_id"]).agg(any_term=("literal_match", "any"),
                                                      cond=("condition_literal_match", "first"),
                                                      n=("condition_literal_match", "nunique"))
    assert (per["n"] == 1).all()
    assert (per["cond"] | ~per["any_term"]).all()
    anyl = tc.groupby("nct_id")["condition_literal_match"].any()
    assert (t.set_index("nct_id")["any_literal_match"] == anyl.reindex(t["nct_id"]).values).all()


@pytest.mark.data
def test_connecticut_uses_2024_planning_regions(tables):
    s = tables["trial_sites"]
    geo = read_table("geographies")
    legacy = set(geo.loc[geo["ct_legacy"].astype(bool), "geo_id"])
    ct_sites = s[(s["country_group"] == "US") & (s["reported_state_fips"] == "09")]
    assert len(ct_sites)
    assert not ct_sites["county_fips"].isin(legacy).any()
    assert ct_sites["county_fips"].dropna().str.match(r"^091[1-9]0$").all()


@pytest.mark.data
def test_quoting_check_reports_recall_cost():
    qlog = json.loads((RAW / ct.SOURCE_ID / "query_log.json").read_text())
    for r in qlog["quoting_check"]:
        assert r["n_quoted_not_in_unquoted"] == 0
        assert 0 <= r["n_of_those_naming_condition_term"] <= r["n_unquoted_only_outside_condition"]


@pytest.mark.data
def test_trials_for_condition_helper(tables):
    lc = ct.trials_for_condition("long_covid")
    assert len(lc) > 100 and lc["nct_id"].is_unique
    lit = ct.trials_for_condition("long_covid", literal_only=True)
    assert 0 < len(lit) <= len(lc)
    hrv = ct.trials_for_condition("long_covid", "hrv")
    assert 0 < len(hrv) < len(lc)
    assert set(hrv["nct_id"]) <= set(lc["nct_id"])
    assert hrv["measurement_match_fields"].notna().all()
    assert "not evidence" in hrv.attrs["guardrail"]
    rec = ct.trials_for_condition("me_cfs", statuses=["RECRUITING"])
    assert (rec["overall_status"] == "RECRUITING").all()
    with pytest.raises(ValueError):
        ct.trials_for_condition("not_a_condition")


@pytest.mark.data
def test_sites_and_get_trial_helpers(tables):
    sites = ct.trial_sites_for_condition("pots")
    assert len(sites) > 0 and (sites["country_group"] == "US").all()
    nct = tables["clinical_trials"]["nct_id"].iloc[0]
    rec = ct.get_trial(nct)
    assert rec["nct_id"] == nct and rec["trial_conditions"]
    assert ct.get_trial("NCT99999999")["status"] == UNKNOWN


@pytest.mark.data
def test_registry_and_audit_written():
    import yaml

    from measure_it.registry import REQUIRED_FIELDS
    entry = yaml.safe_load((RAW / ct.SOURCE_ID / "registry_entry.yaml").read_text())
    assert not [f for f in REQUIRED_FIELDS if f not in entry]
    assert entry["status"] in {"ingested", "partial", "blocked"}
    audit = (RAW / ct.SOURCE_ID / "DATA_AUDIT.md").read_text()
    assert "dataTimestamp" in audit and "not** evidence" in audit
    manifest = json.loads((RAW / ct.SOURCE_ID / "MANIFEST.json").read_text())
    assert "api_version.json" in manifest["files"]
    assert any(k.startswith("records/") for k in manifest["files"])


@pytest.mark.data
def test_geocode_sites_branches():
    if not processed_path("zcta_centroids").exists():
        pytest.skip("census_geography tables missing")
    sites = pd.DataFrame([
        # Emory: postal city Atlanta, geoPoint = Atlanta centroid (Fulton); ZIP 30322 -> DeKalb
        {"nct_id": "N1", "site_idx": 0, "facility": "Emory University", "city": "Atlanta", "state": "Georgia",
         "zip": "30322", "country": "United States", "site_status": None,
         "geopoint_lat": 33.749, "geopoint_lon": -84.38798},
        # no ZIP -> point-in-polygon fallback
        {"nct_id": "N1", "site_idx": 1, "facility": "X", "city": "Atlanta", "state": "Georgia",
         "zip": None, "country": "United States", "site_status": None,
         "geopoint_lat": 33.749, "geopoint_lon": -84.38798},
        # ZIP typo in another state -> geoPoint used, conflict flagged
        {"nct_id": "N1", "site_idx": 2, "facility": "Y", "city": "Atlanta", "state": "Georgia",
         "zip": "10065", "country": "United States", "site_status": None,
         "geopoint_lat": 33.749, "geopoint_lon": -84.38798},
        # non-U.S.
        {"nct_id": "N1", "site_idx": 3, "facility": "Z", "city": "Oslo", "state": None, "zip": "0372",
         "country": "Norway", "site_status": None, "geopoint_lat": 59.91, "geopoint_lon": 10.75},
        # territory without geoPoint, with ZIP
        {"nct_id": "N1", "site_idx": 4, "facility": "W", "city": "San Juan", "state": None, "zip": "00936",
         "country": "Puerto Rico", "site_status": None, "geopoint_lat": None, "geopoint_lon": None},
    ])
    g = ct.geocode_sites(sites)
    assert g.loc[0, "county_fips"] == "13089" and g.loc[0, "county_fips_pip"] == "13121"
    assert g.loc[0, "county_assignment_method"] == "zip5_as_zcta_internal_point"
    assert g.loc[0, "source_geographic_resolution"] == "city_centroid" and g.loc[0, "geocode_method"] == ct.GEOPOINT_METHOD
    assert g.loc[1, "county_fips"] == "13121"
    assert g.loc[1, "county_assignment_method"] == "geopoint_point_in_polygon_cb_2024_5m"
    assert bool(g.loc[2, "zip_state_conflict"]) and g.loc[2, "county_fips"] == "13121"
    assert g.loc[3, "country_group"] == "non_US" and g.loc[3, "county_fips"] is None
    assert g.loc[3, "source_geographic_resolution"] == "city_centroid"
    assert g.loc[4, "country_group"] == "US_territory" and g.loc[4, "state_fips"] == "72"


@pytest.mark.data
def test_geocode_majority_rule_uses_real_ormond_beach_rows(tables):
    s = tables["trial_sites"]
    ob = s[(s["zip"].fillna("").str.startswith("32174")) & (s["county_fips_pip"] == "12127")]
    if ob.empty:
        pytest.skip("no Ormond Beach 32174 sites in this snapshot")
    # ZCTA 32174 internal point lies in Flagler (12035) but 70.6% of its land is Volusia (12127)
    assert (ob["county_fips_zip"] == "12035").all()
    assert (ob["county_fips"] == "12127").all()
    assert (ob["county_assignment_method"] == "geopoint_county_majority_of_zcta").all()
