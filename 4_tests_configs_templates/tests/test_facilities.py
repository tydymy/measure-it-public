"""Tests for measure_it.facilities (registry, matching, query, shuffle_test)."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest

from measure_it.config import TABLES, UNKNOWN
from measure_it.facilities import matching as M
from measure_it.facilities import registry as R
from measure_it.provenance import check_provenance
from measure_it.store import read_table, table_exists

BANNED = re.compile(r"\bbest\b|\bproves?\b|diagnosed by ai|patient has|definitive biomarker|optimal treatment",
                    re.IGNORECASE)


# --------------------------------------------------------------------------- pure-function unit tests

def test_name_normalisation():
    assert R.norm_name("Headache Wellness Center, P.A.") == "headache wellness center"
    assert R.norm_name("St. Luke's Hosp.") == "saint lukes hospital"
    assert R.norm_name("The Mayo Clinic, Inc.") == "mayo clinic"
    assert R.norm_name("Univ of North Carolina Chapel Hill") == "university of north carolina chapel hill"
    assert R.norm_name(None) == "" and R.norm_name(float("nan")) == ""
    # legal suffixes are dropped only at the end
    assert R.norm_name("PA Medical Group LLC") == "pa medical group"


def test_informative_tokens_drop_generic_words():
    assert R.informative_tokens(R.name_tokens("University of Michigan Health System")) == {"michigan"}
    assert R.informative_tokens(R.name_tokens("Family Medical Clinic")) == set()


def _cp(a, b, level, state, city):
    na, nb = R.norm_name(a), R.norm_name(b)
    return R.classify_pair(na, nb, R.informative_tokens(na.split()), R.informative_tokens(nb.split()), level,
                           R._place_tokens(state, city), None, frozenset(city.split()))


@pytest.mark.data   # _place_tokens reads state names from data/processed/geographies.parquet
def test_classify_pair_rules():
    assert _cp("Mayo Clinic in Rochester", "MAYO CLINIC", "zip", "MN", "rochester") == "contained_name_same_zip"
    assert _cp("Gilbert Neurology", "GILBERT NEUROLOGY, LLC", "zip", "AZ", "gilbert") == "exact_name_same_zip"
    # word order matters for containment
    assert _cp("Department of Gastroenterology University of Utah Hospital", "UTAH GASTROENTEROLOGY LLC",
               "city", "UT", "salt lake city") is None
    # numbered sites are different sites
    assert _cp("OnSite Clinical Solutions - ClinEdge (Site 147)", "OnSite Clinical Solutions - ClinEdge (Site 146)",
               "city", "NC", "charlotte") is None
    # city-level links need >= 2 informative tokens
    assert _cp("Smith Clinic", "Smith Clinic", "city", "OH", "dayton") is None
    # generic-only names never link
    assert _cp("Family Medical Clinic", "Family Medical Clinic", "zip", "OH", "dayton") is None
    # a name that is only the city cannot anchor a containment link
    assert _cp("Akron Childrens Hospital", "AKRON GENERAL MEDICAL CENTER", "zip", "OH", "akron") is None


def test_haversine():
    assert abs(float(R.haversine_km(0, 0, 1, 0)) - 111.19) < 0.1
    assert float(R.haversine_km(40, -75, 40, -75)) == 0.0


def _nodes(rows):
    df = pd.DataFrame(rows, columns=["record_id", "source", "lat", "lon"])
    return df


def test_cluster_constraints():
    nodes = _nodes([("hrsa_site:A", "hrsa", 40.0, -75.0), ("hrsa_site:B", "hrsa", 40.0, -75.0),
                    ("npi:1", "nppes", 40.0, -75.0), ("nih_org:X", "nih", 40.0, -75.0),
                    ("nih_org:Y", "nih", 40.0, -75.0), ("ctgov_site:far", "ctgov", 41.0, -75.0)])
    edges = pd.DataFrame([("hrsa_site:A", "npi:1", "exact_name_same_zip", 0.0),
                          ("hrsa_site:B", "npi:1", "hrsa_org_name_same_zip", 0.0),
                          ("nih_org:X", "npi:1", "exact_name_same_zip", 0.0),
                          ("nih_org:Y", "npi:1", "fuzzy_name_same_zip", 0.0),
                          ("ctgov_site:far", "npi:1", "exact_name_same_city", 111.0)],
                         columns=["a", "b", "method", "dist_km"])
    edges["conf"] = edges["method"].map(R.CONFIDENCE)
    root, acc, rej = R.cluster(nodes, edges)
    assert root["hrsa_site:A"] == root["npi:1"] == root["nih_org:X"]
    assert root["hrsa_site:B"] != root["npi:1"]          # at most one HRSA site per facility
    assert root["nih_org:Y"] != root["npi:1"]            # at most one NIH organisation per facility
    assert root["ctgov_site:far"] != root["npi:1"]       # no facility wider than 25 km
    assert rej == {"two_hrsa_sites": 1, "two_nih_orgs": 1, "span_over_25km": 1}
    # never across states, whatever the link type
    n2 = _nodes([("npi:7", "nppes", 40.0, -75.0), ("npi:8", "nppes", 40.0, -75.0)]).assign(state=["NJ", "NY"])
    e2 = pd.DataFrame([("npi:7", "npi:8", "nppes_same_name_same_zip", 0.0, 0.95)], columns=["a", "b", "method", "dist_km", "conf"])
    root2, _, rej2 = R.cluster(n2, e2)
    assert root2["npi:7"] != root2["npi:8"] and rej2 == {"cross_state": 1}
    fac = R.assign_facilities(nodes.assign(native_id=nodes["record_id"].str.split(":").str[1],
                                           precision="point"), root)
    assert fac.loc["npi:1", "facility_id"] == "hrsa-A"   # HRSA anchor takes priority


def test_round_robin_and_orders():
    orders = {c: np.array([], dtype=int) for c in M.CHARACTERISTICS}
    orders["clinical_specialty_match"] = np.array([3, 1, 2])
    orders["relevant_trial_history"] = np.array([3, 4])
    orders["distance_to_target_population"] = np.array([1, 2, 3, 4, 5])
    sel, via = M.round_robin(orders, 4)
    assert sel == [3, 4, 1, 2]
    assert via[3] == "clinical_specialty_match" and via[4] == "relevant_trial_history"
    sel, _ = M.round_robin(orders, 50)
    assert sorted(sel) == [1, 2, 3, 4, 5]


def test_characteristic_orders_are_separate():
    n = 4
    z = np.zeros(n, dtype=int)
    F = M.Features(n_cond_groups=np.array([1, 0, 0, 0]), n_meas_groups=np.array([1, 0, 0, 0]),
                   n_union_groups=np.array([2, 0, 0, 1]), trials_cond=np.array([0, 5, 0, 0]),
                   trials_cond_open=z, trials_cond_latest=z, nih_core=np.array([0, 0, 3, 0]),
                   nih_active_core=z, nih_obligations=np.full(n, np.nan), tech=z, tech_cond=z,
                   intrinsic_access=z, eligible=np.array([True, True, True, True]), cond_groups=[], meas_groups=[])
    dist = np.array([10.0, 20.0, 30.0, 5.0])
    o = M.characteristic_orders(F, np.arange(n), dist, np.zeros(n, dtype=int))
    assert o["clinical_specialty_match"].tolist() == [0, 3]
    assert o["relevant_trial_history"].tolist() == [1]
    assert o["NIH_research_activity"].tolist() == [2]
    assert o["distance_to_target_population"].tolist() == [3, 0, 1, 2]
    assert len(o["community_access"]) == 0


def test_bh_qvalues_and_geo_bootstrap():
    from measure_it.facilities import shuffle_test as S
    q = S.bh_qvalues(np.array([0.01, 0.04, 0.03, 0.5, np.nan]))
    assert np.allclose(q[:4], [0.04, 0.16 / 3, 0.16 / 3, 0.5]) and np.isnan(q[4])
    assert (q[:4] >= np.array([0.01, 0.04, 0.03, 0.5])).all()
    # county bootstrap: resampling within stratum; a constant difference has a degenerate interval
    diff = np.zeros((6, 2))
    diff[:3, 0] = 1.0
    st = np.array(["metro"] * 3 + ["nonmetro"] * 3)
    gb = S.geo_bootstrap(diff, st, n_boot=200)
    est, lo, hi = gb["metro_minus_nonmetro"]
    assert est[0] == lo[0] == hi[0] == 1.0 and est[1] == 0.0
    assert gb["all"][0][0] == 0.5


def test_resolve_measurement_bundles_and_aliases():
    r = M.resolve_measurement("wearable autonomic monitoring")
    assert r["status"] == "matched" and r["measurement_id"] == "wearable_autonomic_activity_monitoring"
    assert "hrv" in r["members"] and "cardiology" in r["implementer_groups"]
    assert M.resolve_measurement("CPET")["measurement_id"] == "cpet"
    assert M.resolve_measurement("tilt table")["measurement_id"] == "autonomic_function_testing"
    assert M.resolve_measurement("Heart-rate variability")["measurement_id"] == "hrv"
    assert M.resolve_measurement("not a measurement")["status"] == UNKNOWN
    assert M.resolve_measurement("wearable")["measurement_id"] == "wearable_autonomic_activity_monitoring"
    assert M.resolve_measurement("Holter")["measurement_id"] == "ecg_ambulatory"      # shared resolver: class pattern
    amb = M.resolve_measurement("heart rate")
    assert amb["status"] == UNKNOWN and set(amb["candidates"]) == {"hrv", "wearable_heart_rate"}


# --------------------------------------------------------------------------- data tests

@pytest.fixture(scope="module")
def facilities():
    if not table_exists("facilities"):
        pytest.skip("facilities not built")
    return read_table("facilities")


@pytest.fixture(scope="module")
def links():
    if not table_exists("facility_source_links"):
        pytest.skip("facility_source_links not built")
    return read_table("facility_source_links")


@pytest.mark.data
def test_facility_tables_provenance_and_ids(facilities, links):
    for name in ("facilities", "facility_source_links", "facility_trials", "facility_nih_projects",
                 "clinic_registry", "research_site_registry"):
        check_provenance(read_table(name), name)
    assert facilities["facility_id"].is_unique
    assert (facilities["object_id"] == "facility:" + facilities["facility_id"]).all()
    assert set(facilities["geocode_precision"]) <= {"point", "zcta_centroid", "city_centroid", "none"}
    assert links["record_id"].is_unique
    assert links["facility_id"].isin(facilities["facility_id"]).all()
    ok = links["object_id"].fillna("").str.match(r"^(npi:\d{10}|hrsa_site:\S+|)$")
    assert ok.all()
    assert facilities["npi_object_ids"].fillna("").str.fullmatch(r"(npi:\d{10}(\|npi:\d{10})*)?").all()


@pytest.mark.data
def test_resolution_constraints(links):
    multi = links[links.groupby("facility_id")["record_id"].transform("size") > 1]
    per = multi.groupby("facility_id")
    assert (per["source"].apply(lambda s: (s == "hrsa").sum()) <= 1).all()
    assert (per["source"].apply(lambda s: (s == "nih").sum()) <= 1).all()
    # never across states
    st = multi.dropna(subset=["state"]).groupby("facility_id")["state"].nunique()
    assert (st <= 1).all()
    # no facility wider than 25 km
    g = multi.dropna(subset=["lat"])
    for fid, grp in list(g.groupby("facility_id"))[:3000]:
        la, lo = grp["lat"].to_numpy(), grp["lon"].to_numpy()
        d = R.haversine_km(la[:, None], lo[:, None], la[None, :], lo[None, :])
        assert d.max() <= R.MAX_SPAN_KM + 1e-6, fid
    # every link method is a known rule
    assert set(links["match_method"]) <= set(R.CONFIDENCE) | {"single_record", "anchor"}


@pytest.mark.data
def test_no_sponsor_placeholders_in_facilities(links):
    from measure_it.ingestion.clinicaltrials import GENERIC_FACILITY_RE
    ct = links[links["source"] == "ctgov"]
    assert len(ct) > 0
    assert not ct["name"].map(lambda x: bool(GENERIC_FACILITY_RE.match(str(x)))).any()
    assert not ct["name"].str.contains(R.NONFACILITY_RE).any()


@pytest.mark.data
def test_research_registry_counts():
    from measure_it.ingestion.nih_reporter import obligation_total
    rr = read_table("research_site_registry")
    ft = read_table("facility_trials")
    fn = read_table("facility_nih_projects")
    assert rr["facility_id"].is_unique
    # no burden variable in the research registry
    assert not any(c.startswith(("burden", "prevalence")) for c in rr.columns)
    # precision-view flag: literal-match trial or title/abstract NIH match
    assert (rr["has_precision_view_record"] ==
            ((rr["n_trials_literal"] > 0) | (rr["n_nih_projects_precision"] > 0))).all()
    n = ft.groupby("facility_id")["nct_id"].nunique()
    assert (rr.set_index("facility_id")["n_trials"].reindex(n.index) == n).all()
    # obligations counted once per appl_id via obligation_total
    top = rr.sort_values("n_nih_projects", ascending=False).iloc[0]
    assert np.isclose(top["nih_award_obligations_usd"], obligation_total(fn[fn["facility_id"] == top["facility_id"]]))
    # every U.S. identifiable NIH appl_id and trial is attached to a facility
    assert fn["facility_id"].notna().all() and ft["facility_id"].notna().all()


@pytest.mark.data
def test_find_candidate_clinics_demo(facilities):
    r = M.find_candidate_clinics("06073", "Long COVID", "wearable autonomic monitoring")
    assert r["status"] == "ok"
    assert 0 < len(r["candidates"]) <= r["max_sites"]
    assert r["geography"]["geo_id"] == "06073" and r["condition"]["condition_id"] == "long_covid"
    ids = [c["facility_id"] for c in r["candidates"]]
    assert len(ids) == len(set(ids))
    for c in r["candidates"]:
        assert c["framing"] == M.FRAMING
        assert set(c["characteristics"]) == set(M.CHARACTERISTICS)
        assert c["object_id"].startswith("facility:")
        assert c["reasons"] and not any(BANNED.search(x) for x in c["reasons"])
        assert c["distance_km"] <= r["radius_km"] or c["inside_geography"]
        via = c["characteristics"][c["selected_via"]]
        assert via["rank"] is not None and via["rank"] >= 1
    assert "score" not in {k for c in r["candidates"] for k in c}
    # aliases resolve identically
    r2 = M.find_candidate_clinics("San Diego County, California", "PASC", "wearables")
    assert [c["facility_id"] for c in r2["candidates"]] == ids
    df = M.candidates_frame(r)
    assert len(df) == len(ids) and "relevant_trial_history__rank" in df.columns


@pytest.mark.data
def test_find_candidate_clinics_unknowns(facilities):
    assert M.find_candidate_clinics("99999", "long_covid", "hrv")["status"] == UNKNOWN
    assert M.find_candidate_clinics("06073", "not a disease xyz", "hrv")["status"] == UNKNOWN
    assert M.find_candidate_clinics("06073", "long_covid", "not a measurement xyz")["status"] == UNKNOWN
    # a county with no geocoded candidate facility inside it and a tiny radius -> UNKNOWN
    cand = facilities[facilities["is_clinic_candidate"].astype(bool)]
    g = read_table("geographies")
    empty = g[(g["geo_level"] == "county") & ~g["ct_legacy"] & ~g["geo_id"].isin(set(cand["county_fips"].dropna()))]
    assert len(empty) > 0
    r = M.find_candidate_clinics(empty["geo_id"].iloc[0], "long_covid", "hrv", radius_km=0.01)
    assert r["status"] == UNKNOWN and "no geocoded facility" in r["reason"]


@pytest.mark.data
def test_find_relevant_research_centers():
    from measure_it.facilities.query import find_relevant_research_centers
    r = find_relevant_research_centers("ME/CFS")
    assert r["status"] == "ok" and r["n_returned"] > 0
    c = r["centers"][0]
    assert c["object_id"].startswith("facility:")
    assert all(x.startswith("trial:") for x in c["trial_object_ids"])
    assert all(x.startswith("nih:") for x in c["nih_object_ids"])
    assert find_relevant_research_centers("xyz not a condition")["status"] == UNKNOWN
    rm = find_relevant_research_centers("long_covid", measurement="wearables")
    assert rm["status"] == "ok" and "n_trials_mentioning_measurement" in rm["centers"][0]


@pytest.mark.data
def test_research_centers_read_tables_once_per_version(monkeypatch):
    """Table-level cache: repeated calls do not re-read facility_trials / facility_nih_projects / facilities; a new
    table version (mtime) is a new cache key; results are unchanged by the cache."""
    from measure_it.facilities import query as FQ
    FQ._tables.cache_clear()
    reads = []
    real = FQ.read_table
    monkeypatch.setattr(FQ, "read_table", lambda name, **kw: reads.append(name) or real(name, **kw))
    a = FQ.find_relevant_research_centers("ME/CFS", top_n=5)
    b = FQ.find_relevant_research_centers("long_covid", measurement="wearable", top_n=5)
    c = FQ.find_relevant_research_centers("ME/CFS", top_n=5)
    assert sorted(reads) == ["facilities", "facility_nih_projects", "facility_trials"]
    import json
    assert json.dumps(a, default=str, sort_keys=True) == json.dumps(c, default=str, sort_keys=True)   # NaN-safe
    assert b["status"] == "ok"
    monkeypatch.setattr(FQ, "_stamp", lambda: ("new version",))
    FQ.find_relevant_research_centers("ME/CFS", top_n=5)
    assert len(reads) == 6
    FQ._tables.cache_clear()


@pytest.mark.data
def test_test6_outputs():
    p = TABLES / "test6_clinic_shuffle.csv"
    if not p.exists():
        pytest.skip("Test 6 not run")
    s = pd.read_csv(p)
    assert {"observed", "null_mean", "null_p2_5", "null_p97_5", "p_two_sided_empirical", "q_bh_fdr_family",
            "geo_boot_ci_lo", "geo_boot_ci_hi", "robust_q05_and_geo_ci"} <= set(s.columns)
    t = s[s["metric"] != "jaccard_vs_observed"]
    assert (t["q_bh_fdr_family"] >= t["p_two_sided_empirical"] - 1e-12).all()
    assert (t.groupby(["condition_id", "measurement_id", "null"]).size() == 40).all()
    chk = pd.read_csv(TABLES / "test6_clinic_shuffle_checks.csv").set_index("check")["value"]
    assert float(chk["pool_size_invariance_violations"]) == 0
    assert float(chk["slot_path_vs_matching_select_mismatches_first25_geos"]) == 0
