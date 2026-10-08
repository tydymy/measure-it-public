"""Tests for measure_it.agents.datagov (Data.gov catalog discovery agent).

Unit tests use a fake HTTP layer (no network). @pytest.mark.data tests check the
processed outputs written by `uv run python -m measure_it.agents.datagov`.
"""
from __future__ import annotations

import json

import pytest

from measure_it import http
from measure_it.agents import datagov as dg
from measure_it.config import UNKNOWN
from measure_it.provenance import PROVENANCE_COLUMNS

# A trimmed record in the shape the catalog.data.gov/search API returns (DCAT-US).
REC = {
    "title": "Lyme disease public use aggregated data with geography, 2008-2021",
    "publisher": "Centers for Disease Control and Prevention",
    "description": "County-level Lyme disease case counts.",
    "identifier": "https://data.cdc.gov/api/views/qtbi-xd4i",
    "slug": "lyme-disease-public-use-aggregated-data-with-geography-2008-2021",
    "keyword": ["lyme", "tickborne"],
    "harvest_record": "https://catalog.data.gov/harvest_record/abc",
    "last_harvested_date": "2026-09-10T00:00:00",
    "organization": {"name": "U.S. Department of Health & Human Services", "slug": "hhs",
                     "organization_type": "Federal Government"},
    "dcat": {
        "accessLevel": "public",
        "identifier": "https://data.cdc.gov/api/views/qtbi-xd4i",
        "landingPage": "https://www.cdc.gov/lyme/data-research/facts-stats/surveillance-data-1.html",
        "modified": "2026-09-09",
        "issued": "2023-01-01",
        "license": "http://opendefinition.org/licenses/odc-odbl/",
        "publisher": {"@type": "org:Organization", "name": "Centers for Disease Control and Prevention",
                      "subOrganizationOf": {"name": "U.S. Department of Health & Human Services"}},
        "distribution": [
            {"downloadURL": "https://data.cdc.gov/api/views/qtbi-xd4i/rows.csv", "mediaType": "text/csv"},
            {"accessURL": "https://data.cdc.gov/resource/qtbi-xd4i.json", "format": "API", "title": "SODA API"},
            {"title": "no url here"},
        ],
    },
}
META = {"api_endpoint": "test", "request_url": "https://catalog.data.gov/search?q=x",
        "http_status": 200, "retrieved_at": "2026-09-23T00:00:00+00:00"}


def _fake_response(payload, status=200):
    return http.CachedResponse(url="https://catalog.data.gov/search?q=x", status=status,
                               headers={"Content-Type": "application/json"},
                               content=json.dumps(payload).encode(), fetched_at="2026-09-23T00:00:00+00:00",
                               from_cache=True)


# ----------------------------------------------------------------- unit tests
def test_normalise_record_fields_and_routing():
    r = dg.normalise_record(REC, rank=1, query="Lyme disease", meta=META, more_available=True)
    assert r["is_unknown"] is False and r["result_status"] == "ok"
    assert r["title"].startswith("Lyme disease public use")
    assert r["agency"] == r["publisher"] == "Centers for Disease Control and Prevention"
    assert r["parent_organization"] == "U.S. Department of Health & Human Services"
    assert "U.S. Department of Health & Human Services" in r["publisher_hierarchy"]
    assert r["access_level"] == "public"
    assert r["modified"] == "2026-09-09"
    assert r["identifier"] == "https://data.cdc.gov/api/views/qtbi-xd4i"
    assert r["landing_page"].startswith("https://www.cdc.gov/lyme")
    # routing goes to the agency landing page, never to the catalog
    assert r["retrieve_from"] == r["landing_page"]
    assert r["retrieve_from_basis"] == "dcat.landingPage"
    assert not dg.is_catalog_url(r["retrieve_from"])
    assert r["machine_readable_url"] == "https://data.cdc.gov/api/views/qtbi-xd4i/rows.csv"
    # distribution without any URL is dropped; accessURL is used when downloadURL is absent
    assert [x["url_kind"] for x in r["resources"]] == ["downloadURL", "accessURL"]
    assert all(set(x) >= {"url", "format", "name"} for x in r["resources"])
    assert r["resources"][1]["format"] == "API" and r["resources"][1]["name"] == "SODA API"
    assert r["resources"][0]["name"] == UNKNOWN
    assert r["datagov_url"].startswith("https://catalog.data.gov/dataset/")


def test_retrieve_from_skips_catalog_hosts():
    res = [{"url": "https://catalog.data.gov/dataset/x", "url_kind": "downloadURL"},
           {"url": "https://www.hrsa.gov/data.csv", "url_kind": "downloadURL"}]
    url, basis = dg.choose_retrieve_from("https://catalog.data.gov/dataset/x", res, "abc")
    assert url == "https://www.hrsa.gov/data.csv" and basis == "dcat.distribution.downloadURL"
    url, basis = dg.choose_retrieve_from(None, [], "https://data.cdc.gov/api/views/abcd-1234")
    assert url.startswith("https://data.cdc.gov") and basis == "dcat.identifier"
    url, _ = dg.choose_retrieve_from(None, [], "not-a-url")
    assert url == UNKNOWN


def test_gsa_records_may_route_to_data_gov_but_others_may_not():
    res = [{"url": "http://www.data.gov/developers/apis", "url_kind": "downloadURL"}]
    assert dg.choose_retrieve_from("https://www.data.gov/developers/harvesting", res, "x")[0] == UNKNOWN
    url, basis = dg.choose_retrieve_from("https://www.data.gov/developers/harvesting", res, "x", dg.GSA_OWN_HOSTS)
    assert url == "https://www.data.gov/developers/harvesting" and basis == "dcat.landingPage"
    rec = {"title": "Data.gov CKAN API", "publisher": "Data.gov",
           "organization": {"name": "General Services Administration"},
           "dcat": {"landingPage": "https://www.data.gov/developers/harvesting"}}
    r = dg.normalise_record(rec, rank=1, query="q", meta=META, more_available=False)
    assert r["retrieve_from"] == "https://www.data.gov/developers/harvesting"
    # catalog.data.gov pages are never an agency host, even for GSA
    assert dg.is_catalog_url("https://catalog.data.gov/dataset/x", dg.GSA_OWN_HOSTS)


def test_object_valued_landing_page_is_a_url():
    # data.va.gov records give landingPage as a DCAT Document object (seen in the raw responses)
    lp = {"@type": "Document", "accessURL": "https://www.data.va.gov/d/39pc-24dr",
          "title": "Emerging Pathogens Initiative (EPI) - Dataset Home"}
    rec = {"title": "Emerging Pathogens Initiative (EPI)", "publisher": "Department of Veterans Affairs",
           "dcat": {"identifier": "VA-VHA-PCS-007", "landingPage": lp, "accessLevel": "non-public"}}
    r = dg.normalise_record(rec, rank=1, query="q", meta=META, more_available=False)
    assert r["landing_page"] == "https://www.data.va.gov/d/39pc-24dr"  # not a JSON dump
    assert r["retrieve_from"] == r["landing_page"] and r["retrieve_from_basis"] == "dcat.landingPage"
    assert dg._url_of([None, {"downloadURL": "https://x.gov/a.csv"}]) == "https://x.gov/a.csv"
    assert dg._url_of({"title": "no url"}) is None


def test_single_distribution_object_is_not_dropped():
    rec = {"title": "T", "dcat": {"distribution": {"downloadURL": "https://www.ers.usda.gov/x.xlsx",
                                                   "mediaType": "application/vnd.ms-excel"}}}
    r = dg.normalise_record(rec, rank=1, query="q", meta=META, more_available=False)
    assert r["resource_urls"] == ["https://www.ers.usda.gov/x.xlsx"]
    assert r["retrieve_from"] == "https://www.ers.usda.gov/x.xlsx"


def test_healthdata_gov_stub_routes_to_agency_distribution():
    # FDA PMA record: landing page is an HHS healthdata.gov catalog page; the data are at accessdata.fda.gov
    res = [{"url": "http://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfPMA/pma.cfm", "url_kind": "downloadURL"}]
    url, basis = dg.choose_retrieve_from("https://healthdata.gov/d/798x-p6ne", res, "x")
    assert url.startswith("http://www.accessdata.fda.gov") and "healthdata.gov catalog page" in basis
    # HHS-hosted data (distributions on healthdata.gov itself) keep the healthdata.gov landing page
    hosted = [{"url": "https://healthdata.gov/api/v3/views/rxn6-qnx8/query.json", "url_kind": "downloadURL"}]
    assert dg.choose_retrieve_from("https://healthdata.gov/d/rxn6-qnx8", hosted, "x") == \
        ("https://healthdata.gov/d/rxn6-qnx8", "dcat.landingPage")
    # no distributions: the healthdata.gov page is the only URL there is
    assert dg.choose_retrieve_from("https://healthdata.gov/d/w6hy-npne", [], "x")[0] == "https://healthdata.gov/d/w6hy-npne"
    # other agency landing pages are untouched
    cdc = [{"url": "https://data.cdc.gov/x.csv", "url_kind": "downloadURL"}]
    assert dg.choose_retrieve_from("https://www.cdc.gov/lyme/", cdc, "x")[1] == "dcat.landingPage"


def test_machine_readable_excludes_documents_and_images():
    assert dg.is_machine_readable_media("text/csv")
    assert dg.is_machine_readable_media("application/json; charset=utf-8")
    for mt in ("application/pdf", "text/html", "application/msword", UNKNOWN, "", None, "image/jpeg",
               "application/vnd.openxmlformats-officedocument.wordprocessingml.document"):
        assert not dg.is_machine_readable_media(mt), mt
    rec = {"title": "T", "dcat": {"distribution": [
        {"downloadURL": "https://pasteur.epa.gov/guide.pdf", "mediaType": "application/pdf"},
        {"downloadURL": "https://pasteur.epa.gov/data.csv", "mediaType": "text/csv"}]}}
    r = dg.normalise_record(rec, rank=1, query="q", meta=META, more_available=False)
    assert r["machine_readable_url"] == "https://pasteur.epa.gov/data.csv"


def test_publisher_as_plain_string_and_missing_fields():
    rec = {"title": "T", "dcat": {"publisher": "Some Agency"}}
    r = dg.normalise_record(rec, rank=1, query="q", meta=META, more_available=False)
    assert r["publisher"] == "Some Agency"
    for f in ("access_level", "landing_page", "modified", "identifier", "retrieve_from", "datagov_url"):
        assert r[f] == UNKNOWN
    assert r["resources"] == []


def test_search_returns_unknown_on_empty(monkeypatch):
    monkeypatch.setattr(dg, "_polite_get", lambda *a, **k: _fake_response({"results": [], "sort": "relevance"}))
    out = dg.search_us_open_data("no such dataset", 5, endpoint="origin")
    assert len(out) == 1
    u = out[0]
    assert u["is_unknown"] is True and u["result_status"] == "no_results"
    assert u["title"] == UNKNOWN and u["agency"] == UNKNOWN and u["retrieve_from"] == UNKNOWN
    assert u["resources"] == []


def test_search_returns_unknown_when_unreachable(monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("connection refused")
    monkeypatch.setattr(dg, "_polite_get", boom)
    out = dg.search_us_open_data("long COVID", 5, endpoint="auto")
    assert len(out) == 1 and out[0]["is_unknown"] and out[0]["result_status"] == "unreachable"
    assert "connection refused" in out[0]["note"]


def test_every_query_failing_is_recorded_blocked_not_a_crash(monkeypatch, tmp_path):
    """A Data.gov outage or rate limit (every query fails) must be recorded as blocked (CONVENTIONS 3.4), not crash
    the step: on an empty result frame pandas keeps the str dtype and '&' with a bool Series used to raise."""
    def boom(*a, **k):
        raise RuntimeError("connection refused")
    captured = {}
    monkeypatch.setattr(dg, "_polite_get", boom)
    monkeypatch.setattr(dg, "raw_dir", lambda sid: tmp_path)
    monkeypatch.setattr(dg, "load_manifest", lambda sid: {"files": {}})
    monkeypatch.setattr(dg, "load_source_registry",
                        lambda: {"sources": [{"source_id": "cdc_lyme", "name": "Lyme disease cases"},
                                             {"source_id": "nppes", "name": "NPPES NPI file"}]})
    monkeypatch.setattr(dg, "write_registry_entry", lambda e: captured.setdefault("entry", e))
    log, res = dg.run_discovery(rows=5)
    assert (log["index_status"] == "query_failed").all()
    assert res["is_unknown"].astype(bool).all()
    stats = dg._write_registry_and_audit(log, res, {"probes": []}, "test")
    assert stats["status"] == "blocked" and stats["queries_failed"] == len(log)
    assert captured["entry"]["status"] == "blocked"
    assert any("queries failed" in x for x in captured["entry"]["limitations"])
    audit = (tmp_path / "DATA_AUDIT.md").read_text()
    assert "| Status | blocked:" in audit


def test_search_http_error_is_unknown_not_crash(monkeypatch):
    monkeypatch.setattr(dg, "_polite_get", lambda *a, **k: _fake_response({"detail": {}, "message": "Not Found"}, 404))
    out = dg.search_us_open_data("x", 5, endpoint="gateway")
    assert out[0]["result_status"] == "unreachable" and "HTTP 404" in out[0]["note"]


def test_search_ok_passes_rows_and_ranks(monkeypatch):
    seen = {}

    def fake(url, *, params, headers, min_interval, **kw):
        seen.update(url=url, params=params, headers=headers)
        return _fake_response({"results": [REC, REC], "after": "cursor", "sort": "relevance"})
    monkeypatch.setattr(dg, "_polite_get", fake)
    monkeypatch.delenv("DATA_GOV_API_KEY", raising=False)
    monkeypatch.delenv("API_DATA_GOV_KEY", raising=False)
    monkeypatch.delenv("DATAGOV_API_KEY", raising=False)
    out = dg.search_us_open_data("Lyme disease", rows=7)
    assert seen["url"] == dg.ORIGIN_SEARCH_URL and seen["params"] == {"q": "Lyme disease", "per_page": 7}
    assert [r["rank"] for r in out] == [1, 2] and all(r["more_available"] for r in out)


def test_gateway_key_goes_in_header_not_params(monkeypatch):
    seen = {}

    def fake(url, *, params, headers, min_interval, **kw):
        seen.update(url=url, params=params, headers=headers)
        return _fake_response({"results": [REC]})
    monkeypatch.setattr(dg, "_polite_get", fake)
    monkeypatch.setenv("DATA_GOV_API_KEY", "secret-test-key")
    dg.search_us_open_data("x", 3)
    assert seen["url"] == dg.GATEWAY_SEARCH_URL
    assert seen["headers"] == {"X-Api-Key": "secret-test-key"}
    assert "api_key" not in seen["params"]  # keeps the key out of the HTTP cache key/meta


def test_cached_error_is_refetched_not_served(monkeypatch):
    # a probe cached the no-key 403 under the same URL+params (the cache key ignores headers)
    calls = []

    def fake(url, *, params, headers, min_interval, refresh=False, **kw):
        calls.append(refresh)
        if not refresh:
            r = _fake_response({"error": {"code": "API_KEY_MISSING"}}, 403)
            return r
        return http.CachedResponse(url=url, status=200, headers={}, content=json.dumps({"results": [REC]}).encode(),
                                   fetched_at="2026-09-23T00:00:00+00:00", from_cache=False)
    monkeypatch.setattr(dg, "_polite_get", fake)
    monkeypatch.setenv("DATA_GOV_API_KEY", "k")
    out = dg.search_us_open_data("lyme", 4, endpoint="gateway")
    assert calls == [False, True]
    assert out[0]["is_unknown"] is False and out[0]["title"].startswith("Lyme")


def test_keyless_probe_cannot_share_a_cache_entry_with_a_search(monkeypatch):
    seen = []

    def fake(url, *, params, headers, **kw):
        seen.append((url, params, headers))
        return _fake_response({"results": []})
    monkeypatch.setattr(dg, "_polite_get", fake)
    dg.probe_endpoints()
    search_keys = {"q", "per_page", "org_type"}  # the only params search_raw sends
    for url, params, headers in seen:
        if url in (dg.GATEWAY_SEARCH_URL, dg.ORIGIN_SEARCH_URL) and not (headers or {}).get("X-Api-Key") \
                and url == dg.GATEWAY_SEARCH_URL:
            assert not set(params) <= search_keys, params  # the no-key 403 probe is keyed apart


def test_invalid_arguments():
    with pytest.raises(ValueError):
        dg.search_us_open_data("  ")
    with pytest.raises(ValueError):
        dg.search_us_open_data("x", rows=0)
    with pytest.raises(ValueError):
        dg.search_us_open_data("x", rows=5, endpoint="ckan")


def test_rules_match_expected_and_reject_wrong_publisher():
    r = dg.normalise_record(REC, rank=1, query="Lyme", meta=META, more_available=False)
    assert dg.matches_rule(r, dg.DISCOVERY_RULES["cdc_lyme"])
    assert dg.match_scope(r, dg.DISCOVERY_RULES["cdc_lyme"]) == "planned_dataset"
    assert not dg.matches_rule(r, dg.DISCOVERY_RULES["cdc_places"])
    fake = dict(r, title="Social Vulnerability Index 2022", publisher="FEMA", publisher_hierarchy="FEMA",
                parent_organization="Department of Homeland Security")
    assert dg.match_scope(fake, dg.DISCOVERY_RULES["cdc_svi"]) is None  # right title, wrong agency
    unk = dg.unknown_result("q", "no_results", "none")
    assert dg.match_scope(unk, dg.DISCOVERY_RULES["cdc_lyme"]) is None


def _as(r, **kw):
    base = dg.normalise_record(REC, rank=1, query="q", meta=META, more_available=False)
    base.update(kw)
    return base


def test_related_vs_planned_scope():
    cdc = dict(publisher="Centers for Disease Control and Prevention",
               publisher_hierarchy="Centers for Disease Control and Prevention",
               parent_organization="U.S. Department of Health & Human Services")
    svi_rule = dg.DISCOVERY_RULES["cdc_svi"]
    assert dg.match_scope(_as(REC, title="Provisional COVID-19 Deaths by Week and County Social Vulnerability Index",
                              **cdc), svi_rule) == "related_dataset"
    assert dg.match_scope(_as(REC, title="CDC/ATSDR Social Vulnerability Index 2022 - United States", **cdc),
                          svi_rule) == "planned_dataset"
    lyme = dg.DISCOVERY_RULES["cdc_lyme"]
    assert dg.match_scope(_as(REC, title="Lyme disease public use line-listed data without geography, 2022-2023",
                              **cdc), lyme) == "related_dataset"
    nh = dg.DISCOVERY_RULES["nhanes_2011_2014"]
    assert dg.match_scope(_as(REC, title="National Health and Nutrition Examination Survey (NHANES) - National "
                                         "Cardiovascular Disease Surveillance System", **cdc), nh) == "related_dataset"
    assert dg.match_scope(_as(REC, title="NHANES 2013-2014 Examination Data", **cdc), nh) == "planned_dataset"
    hrsa = dict(publisher="Health Resources and Services Administration", publisher_hierarchy="HRSA",
                parent_organization="U.S. Department of Health & Human Services")
    hpsa = dg.DISCOVERY_RULES["hrsa_hpsa"]
    assert dg.match_scope(_as(REC, title="Find Shortage Areas: HPSA & MUA/P by Address", **hrsa), hpsa) == "related_dataset"
    assert dg.match_scope(_as(REC, title="Health Professional Shortage Areas - Primary Care", **hrsa), hpsa) == "planned_dataset"


def test_blog_records_never_match_and_case_sensitive_tokens():
    fda = dict(publisher="Food and Drug Administration", publisher_hierarchy="Food and Drug Administration",
               parent_organization="U.S. Department of Health & Human Services")
    rule = dg.DISCOVERY_RULES["openfda_device"]
    assert dg.match_scope(_as(REC, title="Blog | OpenFDA Makes Medical Device-Related Data Easier", **fda), rule) is None
    assert dg.match_scope(_as(REC, title="Product Classification", **fda), rule) == "planned_dataset"
    nih = dict(publisher="National Institutes of Health", publisher_hierarchy="NIH",
               parent_organization="U.S. Department of Health & Human Services")
    rep = dg.DISCOVERY_RULES["nih_reporter"]
    assert dg.match_scope(_as(REC, title="Cre reporter strains produced by targeted insertion", **nih), rep) is None
    assert dg.match_scope(_as(REC, title="NIH Research Portfolio Online Reporting Tools (RePORTER)", **nih), rep) \
        == "planned_dataset"


def test_best_match_prefers_planned_over_better_ranked_related():
    cdc = dict(publisher="CDC", publisher_hierarchy="CDC", parent_organization="HHS")
    res = [_as(REC, rank=1, title="Lyme disease public use line-listed data without geography", **cdc),
           _as(REC, rank=2, title="Lyme disease public use aggregated data with geography, 2008-2021", **cdc)]
    m, scope = dg.best_match(res, dg.DISCOVERY_RULES["cdc_lyme"])
    assert scope == "planned_dataset" and m["rank"] == 2


def test_strip_html():
    assert dg.strip_html("<p>A &amp; B</p>\n<ul><li>C</li></ul>") == "A & B C"
    assert dg.strip_html(None) is None


def test_uncurated_rule_fallback():
    rule = dg._rule_for({"source_id": "new_src", "name": "Example Widget Survey (beta)"})
    assert rule["curated"] is False and rule["queries"] == ["Example Widget Survey"]
    ok = {"is_unknown": False, "title": "The Example Widget Survey 2025"}
    assert dg.matches_rule(ok, rule)
    assert not dg.matches_rule({"is_unknown": False, "title": "Example Survey"}, rule)


def test_every_registry_source_has_a_rule():
    from measure_it.registry import load_source_registry
    ids = {s["source_id"] for s in load_source_registry()["sources"]}
    assert ids <= set(dg.DISCOVERY_RULES), ids - set(dg.DISCOVERY_RULES)


# ----------------------------------------------------------------- data tests
@pytest.fixture(scope="module")
def log():
    from measure_it.store import read_table
    return read_table("data_gov_discovery_log")


@pytest.fixture(scope="module")
def results():
    from measure_it.store import read_table
    return read_table("data_gov_search_results")


@pytest.mark.data
def test_outputs_have_provenance(log, results):
    for df in (log, results):
        assert set(PROVENANCE_COLUMNS) <= set(df.columns)
        assert (df["evidence_type"] == "metadata_catalog").all()
        assert (df["source_geographic_resolution"] == "none").all()
        assert df["retrieved_at"].str.startswith("2026").all()


@pytest.mark.data
def test_log_covers_every_registry_source_and_demo_query(log):
    from measure_it.registry import load_source_registry
    ids = {s["source_id"] for s in load_source_registry()["sources"]}
    src = log[log["query_kind"] == "registry_source"]
    assert ids <= set(src["matched_source_id"])
    assert set(dg.DEMO_QUERIES) == set(log.loc[log["query_kind"] == "demo", "query"])
    required = {"query", "n_results", "top_title", "top_publisher", "top_resource_urls",
                "matched_source_id", "indexed", "retrieved_at"}
    assert required <= set(log.columns)
    assert set(log["indexed"]) <= {"yes", "no"}
    # exactly one final and one best attempt per registry source
    assert src.loc[src["is_final_attempt"], "matched_source_id"].is_unique
    assert src.loc[src["is_best_attempt"], "matched_source_id"].is_unique
    assert set(src.loc[src["is_best_attempt"], "matched_source_id"]) == set(src["matched_source_id"])
    assert set(log["index_status"]) <= set(dg.STATUS_ORDER)
    # indexed == yes exactly when the planned dataset was matched
    assert ((log["indexed"] == "yes") == (log["match_scope"] == "planned_dataset")).all()
    assert ((src["source_indexed"] == "yes") == (src["source_index_status"] == "dataset_indexed")).all()


@pytest.mark.data
def test_indexed_rows_route_to_agency_urls(log):
    yes = log[log["indexed"] == "yes"]
    assert len(yes) > 0
    assert (yes["match_title"] != UNKNOWN).all()
    assert not yes["match_title"].str.match(dg.NON_DATASET_TITLE, case=False).any()
    for sid, u in zip(yes["matched_source_id"], yes["match_retrieve_from"]):
        own = dg.GSA_OWN_HOSTS if sid == dg.SOURCE_ID else frozenset()  # GSA's own record may use data.gov
        assert u.startswith("http") and not dg.is_catalog_url(u, own), (sid, u)


@pytest.mark.data
def test_zero_result_queries_are_explicit_unknown(log, results):
    zero = log[log["n_results"] == 0]
    assert (zero["top_title"] == UNKNOWN).all() and (zero["indexed"] == "no").all()
    unk = results[results["is_unknown"].astype(bool)]
    assert (unk["title"] == UNKNOWN).all() and (unk["retrieve_from"] == UNKNOWN).all()
    # every UNKNOWN result row corresponds to a zero-result / failed query in the log
    zero_keys = set(zip(zero["query_kind"], zero["query"]))
    assert set(zip(unk["query_kind"], unk["query"])) <= zero_keys


@pytest.mark.data
def test_results_are_real_catalog_records(results, log):
    real = results[~results["is_unknown"].astype(bool)]
    assert len(real) > 0
    assert (real["title"] != UNKNOWN).all()
    assert real["request_url"].str.contains(r"catalog\.data\.gov/search|api\.gsa\.gov/technology/datagov/v4").all()
    gsa = real["publisher_hierarchy"].str.contains(dg.GSA_PUBLISHER, case=False, regex=True)
    on_catalog = real["retrieve_from"].map(lambda u: u.startswith("http") and dg.is_catalog_url(u))
    assert not (on_catalog & ~gsa).any()
    assert not real["retrieve_from"].map(lambda u: dg._host(u) == "catalog.data.gov").any()
    for rj in real["resources_json"].head(200):
        for res in json.loads(rj):
            assert res["url"].startswith("http") and {"url", "format", "name"} <= set(res)
    # counts in the log agree with the stored result rows
    counts = real.groupby(["query_kind", "query"]).size()
    for _, row in log.iterrows():
        assert counts.get((row["query_kind"], row["query"]), 0) == row["n_results"]


@pytest.mark.data
def test_raw_responses_match_log_and_results(log, results):
    from measure_it.config import raw_dir
    d = raw_dir(dg.SOURCE_ID)
    on_disk = {f"search_responses/{p.name}" for p in (d / "search_responses").glob("*.json")}
    assert on_disk == set(log["raw_response_file"])  # no stale or missing raw responses
    for _, row in log.iterrows():
        raw = json.loads((d / row["raw_response_file"]).read_text())
        assert raw["query"] == row["query"]
        n = len((raw["payload"] or {}).get("results") or [])
        assert n == row["n_results"]
        # the stored result rows are those records, in API order
        got = results[(results["query_kind"] == row["query_kind"]) & (results["query"] == row["query"])
                      & ~results["is_unknown"].astype(bool)].sort_values("rank")
        want = [str((r.get("dcat") or {}).get("identifier") or r.get("identifier")) for r in raw["payload"]["results"]] \
            if n else []
        assert got["identifier"].tolist() == want


@pytest.mark.data
def test_landing_pages_are_urls_and_healthdata_stubs_rerouted(results):
    real = results[~results["is_unknown"].astype(bool)]
    lp = real["landing_page"]
    # URLs as published (NASA uses ivo:// identifiers), never a serialised DCAT object
    assert not lp.str.match(r"\s*[\[{]").any()
    assert (lp.eq(UNKNOWN) | lp.str.match(r"https?://|ivo://")).all()
    on_hd = real[real["retrieve_from_host"].isin(dg.SECONDARY_CATALOG_HOSTS)]
    for rj in on_hd["resources_json"]:
        hosts = {dg._host(x["url"]) for x in json.loads(rj)}
        assert hosts <= dg.SECONDARY_CATALOG_HOSTS | dg.CATALOG_HOSTS, hosts


@pytest.mark.data
def test_raw_artifacts_present():
    from measure_it.config import raw_dir
    d = raw_dir(dg.SOURCE_ID)
    for f in ("MANIFEST.json", "DATA_AUDIT.md", "registry_entry.yaml", "endpoint_probes.json", "openapi.json"):
        assert (d / f).exists(), f
    probes = json.loads((d / "endpoint_probes.json").read_text())
    legacy = [p for p in probes["probes"] if p["label"].startswith("legacy CKAN package_search")]
    assert legacy and legacy[0].get("http_status") == 404
    audit = (d / "DATA_AUDIT.md").read_text()
    assert "TODO" not in audit and "| source_id | data_gov_catalog |" in audit


def test_census_county_file_not_confused_with_cd_within_county():
    census = dict(publisher="U.S. Department of Commerce, U.S. Census Bureau, Geography Division",
                  publisher_hierarchy="U.S. Census Bureau", parent_organization="U.S. Census Bureau")
    rule = dg.DISCOVERY_RULES["census_geography"]
    cd = ("2025 Cartographic Boundary File (SHP), 119th Congressional District within Current County and "
          "Equivalent Entities for United States, 1:500,000")
    assert dg.match_scope(_as(REC, title=cd, **census), rule) == "related_dataset"
    county = "2016 Cartographic Boundary File, Current County and Equivalent for United States, 1:5,000,000"
    assert dg.match_scope(_as(REC, title=county, **census), rule) == "planned_dataset"
    sub = "2023 Cartographic Boundary File (SHP), County Subdivision for Oregon, 1:500,000"
    assert dg.match_scope(_as(REC, title=sub, **census), rule) == "related_dataset"
