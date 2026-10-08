"""Tests for measure_it.ingestion.nih_reporter.

Unit tests are pure functions (no network). @pytest.mark.data tests check the
processed outputs written by `uv run python -m measure_it.ingestion.nih_reporter`.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from measure_it.config import load_config, raw_dir
from measure_it.ingestion import nih_reporter as nr
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import table_exists


# --------------------------------------------------------------------------- unit tests
def test_phrase_and_or_group():
    assert nr.phrase("long COVID") == '"long COVID"'
    assert nr.phrase(' say "x" ') == '"say  x"'
    assert nr.or_group(["a"]) == '"a"'
    assert nr.or_group(["a b", "c"]) == '("a b" OR "c")'
    with pytest.raises(ValueError):
        nr.phrase('  ""  ')


def test_make_criteria_shape():
    c = nr.make_criteria('"x"', [2015, 2016])
    assert c == {"fiscal_years": [2015, 2016],
                 "advanced_text_search": {"operator": "advanced",
                                          "search_field": "projecttitle,abstracttext,terms",
                                          "search_text": '"x"'}}


def test_query_building_covers_config():
    cfg = load_config("conditions")
    qs = nr.all_queries()
    n_terms = sum(len(c["search_terms"]) for c in cfg["conditions"])
    cond = [q for q in qs if q.query_type == "condition_term"]
    aug = [q for q in qs if q.query_type == "augmented"]
    assert len(cond) == n_terms
    assert len(aug) == len(cfg["demo_cluster"]["members"]) * len(nr.AUGMENT_TERMS)
    assert len({q.query_id for q in qs}) == len(qs)
    for q in qs:
        assert q.fiscal_years == list(range(2015, 2027))
    a = next(q for q in aug if q.query_id == "aug:pots:wearable")
    assert a.search_text == ('("postural orthostatic tachycardia syndrome" OR "postural tachycardia syndrome") '
                             'AND ("wearable" OR "wearables")')


def test_page_offsets_and_split():
    assert nr.page_offsets(0) == [0]
    assert nr.page_offsets(500) == [0]
    assert nr.page_offsets(501) == [0, 500]
    assert nr.page_offsets(15000)[-1] == 14500
    assert max(nr.page_offsets(15000)) <= nr.MAX_OFFSET
    with pytest.raises(ValueError):
        nr.page_offsets(15001)
    assert not nr.needs_fiscal_year_split(15000)
    assert nr.needs_fiscal_year_split(15001)


def _fake_poster(per_fy: dict[int, int]):
    """A fake API: per_fy maps fiscal year -> number of records."""
    calls = []

    def poster(body):
        calls.append(body)
        fys = body["criteria"]["fiscal_years"]
        ids = [fy * 100000 + i for fy in fys for i in range(per_fy.get(fy, 0))]
        page = ids[body["offset"]: body["offset"] + body["limit"]]
        return {"meta": {"total": len(ids), "search_id": "x", "properties": {"URL": "u"}},
                "results": [{"appl_id": a, "fiscal_year": a // 100000} for a in page]}
    return poster, calls


def test_fetch_criteria_pages_without_split():
    poster, calls = _fake_poster({2020: 1200})
    recs, log = nr.fetch_criteria(nr.make_criteria('"x"', [2020]), poster=poster)
    assert len(recs) == 1200 and len({r["appl_id"] for r in recs}) == 1200
    assert [c["offset"] for c in calls] == [0, 500, 1000]
    assert all(c["limit"] == 500 and c["sort_field"] == "appl_id" for c in calls)


def test_fetch_criteria_splits_by_fiscal_year_over_offset_cap():
    poster, calls = _fake_poster({2019: 9000, 2020: 8000})
    recs, log = nr.fetch_criteria(nr.make_criteria('"x"', [2019, 2020]), poster=poster)
    assert len(recs) == 17000 and len({r["appl_id"] for r in recs}) == 17000
    assert log[0]["split_by_fiscal_year"] is True
    assert all(c["offset"] <= nr.MAX_OFFSET for c in calls)
    assert {tuple(c["criteria"]["fiscal_years"]) for c in calls[1:]} == {(2019,), (2020,)}


def test_fetch_criteria_single_year_too_large_raises():
    poster, _ = _fake_poster({2020: 16000})
    with pytest.raises(RuntimeError):
        nr.fetch_criteria(nr.make_criteria('"x"', [2020]), poster=poster)


def test_parse_terms():
    assert nr.parse_terms("<A><B C><>") == "A; B C"
    assert nr.parse_terms(None) is None
    assert nr.parse_terms("") is None


def test_phrase_regex_separators_and_boundaries():
    p = nr.phrase_regex("post-acute sequelae of SARS-CoV-2")
    assert p.search("the Post Acute Sequelae of SARS CoV 2 cohort")
    p2 = nr.phrase_regex("ME/CFS")
    assert p2.search("patients with ME/CFS and") and p2.search("ME-CFS")
    assert not p2.search("HOME/CFSX")
    p3 = nr.phrase_regex("PASC")
    assert p3.search("(PASC)") and not p3.search("PASCAL")


def test_phrase_regex_unicode_dash():
    assert nr.phrase_regex("Ehlers-Danlos").search("Classic Ehlers\u2013Danlos syndrome")


def test_acronym_context():
    assert nr.acronym_context("PASC", "PASC after SARS-CoV-2", None) is True
    assert nr.acronym_context("PASC", "pancreatic stellate cells (PASC)", float("nan")) is False
    assert nr.acronym_context("PASC", "Wrap PASC Cohorts", "OTA-21-015A") is True
    assert nr.acronym_context("long COVID", "anything") is None


def test_match_fields_tiers():
    m = nr.match_fields(["long COVID"], "Long COVID study", None, None, None)
    assert m["match_tier"] == "title_abstract" and m["matched_in_title"]
    m = nr.match_fields(["long COVID"], "x", "y", None, "COVID-19; Long COVID")
    assert m["match_tier"] == "terms_only"
    m = nr.match_fields(["long COVID"], "x", "y", "z", "w")
    assert m["match_tier"] == "not_found_locally"
    m = nr.match_fields(["wearable", "wearables"], None, None, "uses wearables daily", None)
    assert m["match_tier"] == "title_abstract" and m["matched_in_abstract"]


def _proj(rows):
    return pd.DataFrame(rows, columns=["appl_id", "project_num", "fiscal_year", "is_subproject", "award_amount"])


def test_obligation_total_parent_includes_subprojects():
    # Parent 100 = subprojects 60 + 40; counting all three would double count.
    df = _proj([(1, "P01X-01", 2020, False, 100.0), (2, "P01X-01", 2020, True, 60.0),
                (3, "P01X-01", 2020, True, 40.0), (4, "R01Y-02", 2020, False, 10.0)])
    assert nr.obligation_total(df) == 110.0
    assert nr.subproject_overlap_rows(df) == 2


def test_obligation_total_subprojects_without_parent_and_duplicates():
    df = _proj([(2, "P01X-01", 2020, True, 60.0), (3, "P01X-01", 2020, True, 40.0),
                (3, "P01X-01", 2020, True, 40.0),              # same appl_id matched by two queries
                (5, "P01X-02", 2021, False, np.nan), (6, "P01X-02", 2021, True, 7.0)])
    assert nr.obligation_total(df) == 107.0
    assert nr.obligation_total(df.iloc[:0]) == 0.0


def test_obligation_total_prefers_parent_when_subprojects_exceed_it():
    # RePORTER inconsistency seen in real data (e.g. 5U54NS105541-05 FY2021: parent 1,922,881 = IC funding,
    # subprojects sum to 3,992,540): the parent is still the counted obligation.
    df = _proj([(1, "U54X-05", 2021, False, 100.0), (2, "U54X-05", 2021, True, 90.0),
                (3, "U54X-05", 2021, True, 60.0)])
    assert nr.obligation_total(df) == 100.0
    df2 = pd.concat([df, _proj([(4, "P01Y-01", 2020, False, 50.0), (5, "P01Y-01", 2020, True, 50.0),
                                (6, "P50Z-02", 2020, False, 80.0), (7, "P50Z-02", 2020, True, 30.0),
                                (8, "R01Q-01", 2020, False, 9.0)])])
    c = nr.parent_subproject_consistency(df2)
    assert (c["groups"], c["equal"], c["parent_larger"], c["parent_smaller"]) == (3, 1, 1, 1)
    assert c["parent_smaller_excess"] == 50.0 and c["parent_smaller_groups"] == ["U54X-05 FY2021"]


def test_save_raw_is_byte_identical_across_runs(tmp_path, monkeypatch):
    import time
    paths = {"records": tmp_path / "records.jsonl.gz", "hits": tmp_path / "query_hits.csv",
             "queries": tmp_path / "queries.json"}
    monkeypatch.setattr(nr, "raw_paths", lambda: paths)
    monkeypatch.setattr(nr, "_manifest_record", lambda *a, **k: None)
    recs = {2: {"appl_id": 2, "x": "é"}, 1: {"appl_id": 1}}
    qlog = [{"query_id": "q", "requests": [{"fetched_at": "2026-09-23T00:00:00+00:00"}]}]
    nr.save_raw(recs, [("q", 2), ("q", 1)], qlog)
    first = nr.sha256_file(paths["records"])
    time.sleep(1.1)  # a gzip header timestamp would change here
    nr.save_raw(recs, [("q", 2), ("q", 1)], qlog)
    assert nr.sha256_file(paths["records"]) == first
    back, hits, _ = nr.load_raw()
    assert back == {1: {"appl_id": 1}, 2: {"appl_id": 2, "x": "é"}} and hits == [("q", 1), ("q", 2)]


def test_flatten_record_minimal():
    r = {"appl_id": 42, "fiscal_year": 2024, "subproject_id": None, "project_num": "5R01HL1-02",
         "core_project_num": "R01HL1", "project_title": "t", "terms": "<a><b>",
         "organization": {"org_name": "U", "org_city": "DALLAS", "org_state": "TX", "org_zipcode": "753909105",
                          "org_country": "UNITED STATES", "org_ipf_code": "578404"},
         "principal_investigators": [{"profile_id": 1, "full_name": "A B", "is_contact_pi": False},
                                     {"profile_id": 2, "full_name": "C D", "is_contact_pi": True}],
         "award_amount": 10, "is_active": True, "project_start_date": "2023-09-15T00:00:00",
         "geo_lat_lon": {"lat": 32.8, "lon": -96.8}}
    f = nr.flatten_record(r)
    assert f["appl_id"] == 42 and f["is_subproject"] is False
    assert f["terms"] == "a; b"
    assert f["contact_pi_profile_id"] == "2" and f["pi_profile_ids"] == "1; 2" and f["n_pis"] == 2
    assert f["org_zipcode"] == "753909105" and f["org_ipf_code"] == "578404"
    assert f["award_amount"] == 10.0 and np.isnan(f["direct_cost_amt"])
    assert str(f["project_start_date"]) == "2023-09-15"


# --------------------------------------------------------------------------- data tests
needs_data = pytest.mark.skipif(not table_exists("nih_projects"), reason="run the nih_reporter module first")


@pytest.fixture(scope="module")
def tables():
    from measure_it.store import read_table
    return read_table("nih_projects"), read_table("nih_project_conditions"), read_table("nih_project_query_hits")


@pytest.mark.data
@needs_data
def test_provenance_and_keys(tables):
    proj, pc, qh = tables
    for df in (proj, pc, qh):
        assert all(c in df.columns for c in PROVENANCE_COLUMNS)
        assert (df["data_layer"] == "facility").all()
    assert proj["appl_id"].is_unique
    assert set(pc["appl_id"]) <= set(proj["appl_id"])
    assert set(qh["appl_id"]) <= set(proj["appl_id"])
    assert not pc.duplicated(["appl_id", "condition_id", "query_id"]).any()
    assert not qh.duplicated(["appl_id", "query_id"]).any()
    assert proj["fiscal_year"].between(2015, 2026).all()
    assert (proj["evidence_type"] == "grant_record").all()
    assert set(pc["evidence_type"]) == {"text_mined"}


@pytest.mark.data
@needs_data
def test_every_project_came_from_a_query(tables):
    proj, pc, qh = tables
    assert set(proj["appl_id"]) == set(pc["appl_id"]) | set(qh["appl_id"])


@pytest.mark.data
@needs_data
def test_geocoding_resolution_consistent(tables):
    proj, _, _ = tables
    z = proj["geocode_method"] == nr.GEO_ZCTA
    pt = proj["geocode_method"] == nr.GEO_REPORTER
    none = proj["geocode_method"] == "not_geocoded"
    assert (z | pt | none).all()
    assert (proj.loc[z, "source_geographic_resolution"] == "zcta_centroid").all()
    assert (proj.loc[pt, "source_geographic_resolution"] == "point").all()
    assert (proj.loc[none, "source_geographic_resolution"] == "none").all()
    assert proj.loc[none, ["lat", "lon", "county_fips"]].isna().all().all()
    geo = z | pt
    assert proj.loc[geo, "lat"].between(-15, 72).all() and proj.loc[geo, "lon"].between(-180, 146).all()
    assert proj.loc[geo, "county_fips"].str.fullmatch(r"\d{5}").all()
    # the RePORTER-point fallback is only used when the ZIP failed, and must sit in the org's own state
    assert proj.loc[pt, "zip_geocode_failed"].all()
    assert proj.loc[pt, "geo_state_matches_org_state"].astype(bool).all()
    # foreign organizations are never geocoded (e.g. a Nigerian postal code read as a U.S. ZIP)
    us = proj["org_country"].eq("UNITED STATES")
    assert not geo[~us].any()
    assert geo[us & proj["org_zipcode"].notna()].mean() > 0.99  # measured, see DATA_AUDIT.md


@pytest.mark.data
@needs_data
def test_queries_json_records_exact_criteria(tables):
    q = json.loads((raw_dir(nr.SOURCE_ID) / "queries.json").read_text())
    ids = {x["query_id"] for x in q["queries"]}
    assert ids == {s.query_id for s in nr.all_queries()}
    for x in q["queries"]:
        crit = x["request_body_template"]["criteria"]
        assert crit["fiscal_years"] == list(range(2015, 2027))
        assert crit["advanced_text_search"]["search_text"] == x["search_text"]
        assert "warning" not in x, x.get("warning")
    manifest = json.loads((raw_dir(nr.SOURCE_ID) / "MANIFEST.json").read_text())
    assert {"records.jsonl.gz", "query_hits.csv", "queries.json"} <= set(manifest["files"])


@pytest.mark.data
@needs_data
def test_condition_links_match_query_hits(tables):
    _, pc, qh = tables
    import pandas as pd
    hits = pd.read_csv(raw_dir(nr.SOURCE_ID) / "query_hits.csv")
    cond_hits = hits[hits["query_id"].str.startswith("cond:")]
    aug_hits = hits[hits["query_id"].str.startswith("aug:")]
    assert len(cond_hits) == len(pc)
    assert len(aug_hits) == len(qh)
    assert set(pc["match_tier"]) <= {"title_abstract", "terms_only", "not_found_locally"}
    assert set(qh["condition_id"]) == set(load_config("conditions")["demo_cluster"]["members"])
    assert set(qh["augmented_term"]) <= set(nr.AUGMENT_TERMS)


@pytest.mark.data
@needs_data
def test_rcdc_categories_exist_in_data(tables):
    proj, pc, _ = tables
    cats = {c.strip() for s in proj["spending_categories_desc"].dropna() for c in s.split(";")}
    for cid, cat in nr.RCDC_CATEGORY_FOR_CONDITION.items():
        assert cat in cats, f"{cid}: RCDC category {cat!r} never observed"
        sub = pc[(pc["condition_id"] == cid)]
        assert sub["rcdc_category_assigned"].notna().all()


@pytest.mark.data
@needs_data
def test_helpers(tables):
    proj, pc, _ = tables
    p_all = nr.projects_for_condition("pots", include_terms_only=True, include_likely_false_positives=True)
    assert p_all["appl_id"].is_unique
    assert set(p_all["appl_id"]) == set(pc.loc[pc["condition_id"] == "pots", "appl_id"])
    assert set(p_all["best_match_tier"]) <= {"title_abstract", "terms_only", "not_found_locally"}
    p = nr.projects_for_condition("pots")  # precision view
    assert set(p["best_match_tier"]) == {"title_abstract"} and set(p["appl_id"]) <= set(p_all["appl_id"])
    lc = nr.projects_for_condition("long_covid", include_terms_only=True)
    fp_ids = set(pc.loc[(pc["condition_id"] == "long_covid") & pc["likely_false_positive"], "appl_id"])
    other = set(pc.loc[(pc["condition_id"] == "long_covid") & ~pc["likely_false_positive"], "appl_id"])
    assert fp_ids and not (fp_ids - other) & set(lc["appl_id"])
    rc = nr.research_centers_for_condition("pots")
    assert rc["n_projects"].sum() == len(p)
    assert np.isclose(rc["award_obligations_total_usd"].sum(), nr.obligation_total(p))
    me = nr.research_centers_for_condition("me_cfs")
    all_missing = me["n_award_amount_missing"] == me["n_projects"]
    assert me.loc[all_missing, "award_obligations_total_usd"].isna().all()  # VA: unknown, not $0
    assert (rc["award_obligations_total_usd"] <= rc["award_obligations_total_usd"].sum()).all()
    from measure_it.provenance import check_provenance
    check_provenance(rc, "research_centers_for_condition(pots)")  # full provenance contract, not a subset
    assert rc["source_record_id"].is_unique
    empty = nr.research_centers_for_condition("post_infectious_syndrome")
    assert empty.empty and "org_name" in empty.columns
    with pytest.raises(ValueError):
        nr.projects_for_condition("not_a_condition")


@pytest.mark.data
@needs_data
def test_condition_totals_not_double_counted(tables):
    proj, pc, _ = tables
    for cid in pc["condition_id"].unique():
        sub = proj[proj["appl_id"].isin(pc.loc[pc["condition_id"] == cid, "appl_id"])]
        dedup = nr.obligation_total(sub)
        assert dedup <= sub["award_amount"].sum() + 1e-6
    # union total is below the sum of overlapping per-condition totals
    union = nr.obligation_total(proj)
    per = sum(nr.obligation_total(proj[proj["appl_id"].isin(pc.loc[pc["condition_id"] == c, "appl_id"])])
              for c in pc["condition_id"].unique())
    assert union <= per + 1e-6


@pytest.mark.data
@needs_data
def test_connecticut_uses_2024_planning_regions(tables):
    proj, _, _ = tables
    ct = proj["county_fips"].fillna("").str.startswith("09")
    assert ct.any()
    assert proj.loc[ct, "county_fips"].str.fullmatch(r"091[1-9]0").all()      # 2024 planning regions only
    assert (proj.loc[ct, "county_assignment"] != "zcta_county_rel2020_max_land_share").all()
    audit = (raw_dir(nr.SOURCE_ID) / "DATA_AUDIT.md").read_text()
    assert f"{int(ct.sum())} carry a CT county FIPS" in audit and " and 0 are legacy counties" in audit


@pytest.mark.data
@needs_data
def test_audit_parent_subproject_numbers_match_table(tables):
    proj, _, _ = tables
    c = nr.parent_subproject_consistency(proj)
    audit = (raw_dir(nr.SOURCE_ID) / "DATA_AUDIT.md").read_text()
    assert (f"groups in the table that hold a parent and at least one subproject: parent = sum of the "
            f"retrieved subprojects in {c['equal']}; parent larger in {c['parent_larger']}") in audit
    assert f"parent **smaller** in {c['parent_smaller']}" in audit
    # the parent amount is the IC-funding total, which is why obligation_total prefers it
    ic = proj["ic_fundings"].dropna().map(lambda j: sum(float(f.get("total_cost") or 0) for f in json.loads(j)))
    amt = proj.loc[ic.index, "award_amount"]
    assert ((amt - ic).abs()[amt.notna()] <= 1).all()


@pytest.mark.data
@needs_data
def test_audit_and_registry_written():
    import yaml
    from measure_it.registry import REQUIRED_FIELDS
    d = raw_dir(nr.SOURCE_ID)
    audit = (d / "DATA_AUDIT.md").read_text()
    assert "award obligations" in audit.lower() and "not** a total" in audit
    reg = yaml.safe_load((d / "registry_entry.yaml").read_text())
    assert all(f in reg for f in REQUIRED_FIELDS)
    assert reg["data_layer"] == "facility" and reg["person_level"] is False
