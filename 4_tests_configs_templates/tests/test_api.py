"""Tests for the FastAPI app (measure_it.api.app) and the display boundaries (measure_it.api.geo).

Every route is exercised with the FastAPI TestClient: one GET and one POST per tool (demo arguments and an input that
must come back UNKNOWN / NOT AVAILABLE), /health, /tools, /sources, /geojson/{level}. Tool routes must return exactly
the envelope of the same function in measure_it.tools (compared without `meta`, which carries timings). Tests marked
`data` need data/processed; Data.gov calls run with MEASURE_IT_OFFLINE=1 (HTTP cache only, never the network).
"""
from __future__ import annotations

import json
import re

import pytest
from fastapi.testclient import TestClient

from measure_it import tools as T
from measure_it.api import app as A
from measure_it.api import geo as G
from measure_it.config import PROCESSED, UNKNOWN

NONSENSE = "zzqxv blorf quux"
# tool -> (demo arguments, arguments that must give UNKNOWN / NOT AVAILABLE)
CASES = {
    "search_condition": ({"condition": "long covid"}, {"condition": NONSENSE}),
    "normalize_condition": ({"condition_or_code": "G93.32"}, {"condition_or_code": NONSENSE}),
    "get_patient_phenotype_signature": ({"condition": "ME/CFS", "top_n": 3}, {"condition": NONSENSE}),
    "get_molecular_context": ({"condition": "Long COVID", "top_n": 3}, {"condition": NONSENSE}),
    "discover_candidate_measurements": ({"condition": "Long COVID", "phenotype": "orthostatic intolerance"},
                                        {"condition": NONSENSE}),
    "get_measurement_evidence": ({"measurement": "wearable autonomic monitoring", "condition": "Long COVID"},
                                 {"measurement": NONSENSE, "condition": "Long COVID"}),
    "get_regulatory_context": ({"technology": "wearable ECG patch"}, {"technology": NONSENSE}),
    "get_condition_burden": ({"condition": "Long COVID", "geography_level": "county", "geo_id": "06073"},
                             {"condition": "POTS", "geography_level": "county"}),
    "get_geographic_context": ({"geography_id": "San Diego County, California"},
                               {"geography_id": "Zzqxv County, Nowhere"}),
    "find_candidate_clinics": ({"geography_id": "06073", "condition": "Long COVID",
                                "measurement": "wearable autonomic monitoring", "max_sites": 5},
                               {"geography_id": "99999", "condition": "Long COVID"}),
    "find_relevant_trials": ({"condition": "Long COVID", "measurement": "wearable autonomic monitoring",
                              "statuses": ["RECRUITING"], "max_trials": 5}, {"condition": NONSENSE}),
    "find_relevant_research_centers": ({"condition": "ME/CFS", "top_n": 5}, {"condition": NONSENSE}),
    "rank_deployment_opportunities": ({"condition": "Long COVID or ME/CFS",
                                       "measurement": "wearable autonomic monitoring", "top_n": 3},
                                      {"condition": NONSENSE, "measurement": "wearable autonomic monitoring"}),
    "search_us_open_data": ({"query": "long COVID"}, {"query": NONSENSE}),
    "trace_evidence": ({"object_id": "geo:06073"}, {"object_id": "zzqxv:blorf"}),
    "get_phenotype_measurement_evidence": ({"phenotype": "orthostatic intolerance"}, {"phenotype": NONSENSE}),
    "get_measurable_biology": ({"condition": "ME/CFS"}, {"condition": NONSENSE}),
    "list_sources": ({}, {"data_layer": NONSENSE}),
}
ENVELOPE_KEYS = {"tool", "status", "query", "data", "caveats", "provenance", "truncation", "meta"}
HAS_DATA = (PROCESSED / "condition_registry.parquet").exists()
HAS_BOUNDARIES = (PROCESSED / "county_boundaries_2024.geoparquet").exists()


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    monkeypatch.setenv("MEASURE_IT_OFFLINE", "1")


@pytest.fixture(scope="module")
def client():
    return TestClient(A.app)


def _no_meta(env: dict) -> dict:
    return {k: v for k, v in env.items() if k != "meta"}


def _params(args: dict) -> dict:
    return {k: v for k, v in args.items() if v is not None}


# ------------------------------------------------------------------------------------------------ pure

def test_every_tool_has_a_get_and_a_post_route():
    routes = {(r.path, m) for r in A.app.routes for m in getattr(r, "methods", set())}
    assert set(CASES) == set(T.TOOLS)
    for name in T.TOOLS:
        assert (f"/tools/{name}", "GET") in routes, name
        assert (f"/tools/{name}", "POST") in routes, name
    for path in ("/", "/health", "/tools", "/sources", "/sources/{source_id}", "/geojson/{level}"):
        assert (path, "GET") in routes, path


def test_tools_listing_describes_every_parameter(client):
    r = client.get("/tools")
    assert r.status_code == 200
    d = r.json()
    assert d["n_tools"] == len(T.TOOLS) == 18
    assert [t["name"] for t in d["tools"]] == T.SPEC_TOOLS + T.EXTRA_TOOLS
    by = {t["name"]: t for t in d["tools"]}
    rank = {p["name"]: p for p in by["rank_deployment_opportunities"]["parameters"]}
    assert rank["geography_level"]["enum"] == ["county", "state"] and rank["detail"]["enum"] == ["compact", "full"]
    assert rank["condition"]["required"] and not rank["top_n"]["required"] and rank["top_n"]["default"] == 10
    for t in d["tools"]:
        assert t["wraps"] and t["routes"] == {"get": f"/tools/{t['name']}", "post": f"/tools/{t['name']}"}
        assert all(p["description"] for p in t["parameters"])


def test_health(client):
    r = client.get("/health")
    assert r.status_code == 200
    d = r.json()
    assert d["status"] in ("ok", "degraded") and d["n_tools"] == 18 and d["unknown_sentinel"] == UNKNOWN
    assert d["offline"] is True


def test_validation_errors_are_422_before_the_tool_runs(client):
    assert client.get("/tools/search_condition").status_code == 422                       # missing required
    assert client.get("/tools/get_condition_burden", params={"condition": "x", "geography_level": "planet"}) \
        .status_code == 422                                                                # enum
    assert client.get("/tools/rank_deployment_opportunities", params={
        "condition": "x", "measurement": "y", "detail": "everything"}).status_code == 422
    assert client.get("/tools/search_condition", params={"condition": "x", "limit": "many"}).status_code == 422
    assert client.post("/tools/normalize_condition", json={}).status_code == 422            # missing in body
    assert client.post("/tools/normalize_condition").status_code == 422                    # no body at all
    assert client.post("/tools/list_sources", json={"bogus": 1}).status_code == 422        # unknown key
    assert client.get("/geojson/zcta").status_code == 422
    assert client.get("/geojson/county", params={"tolerance": 5}).status_code == 422


def test_error_envelope_maps_to_http_500_with_the_same_body():
    env = {"tool": "x", "status": "error", "reason": "boom", "data": None}
    r = A._respond(env)
    assert r.status_code == 500 and json.loads(r.body) == env and r.headers["x-tool-status"] == "error"
    assert A._respond({"status": UNKNOWN}).status_code == 200


def test_openapi_documents_the_tool_routes(client):
    spec = client.get("/openapi.json").json()
    assert "/tools/rank_deployment_opportunities" in spec["paths"]
    ops = spec["paths"]["/tools/get_condition_burden"]
    names = {p["name"] for p in ops["get"]["parameters"]}
    assert {"condition", "geography_level", "geo_id", "max_rows", "order"} <= names
    assert "requestBody" in ops["post"]


def test_geo_argument_checks():
    with pytest.raises(ValueError):
        G.boundaries_geojson("tract")
    with pytest.raises(ValueError):
        G.boundaries_geojson("county", tolerance=-1)
    with pytest.raises(ValueError):
        G.boundaries_geojson("county", state="CA")


# ------------------------------------------------------------------------------------------------ tools over HTTP

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
@pytest.mark.parametrize("name", list(CASES))
def test_tool_get_and_post_return_the_facade_envelope(client, name):
    demo, _ = CASES[name]
    direct = T.TOOLS[name](**demo)
    rg = client.get(f"/tools/{name}", params=_params(demo))
    rp = client.post(f"/tools/{name}", json=demo)
    for r in (rg, rp):
        assert r.status_code == 200, (name, r.text[:300])
        env = r.json()
        assert ENVELOPE_KEYS <= set(env), name
        assert env["tool"] == name and env["status"] == "ok", (name, env.get("reason"))
        assert r.headers["x-tool-status"] == "ok"
        assert _no_meta(env) == json.loads(json.dumps(_no_meta(direct))), name
    assert rg.json()["query"] == rp.json()["query"]


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
@pytest.mark.parametrize("name", list(CASES))
def test_tool_unknown_inputs_return_unknown_with_a_reason(client, name):
    _, bad = CASES[name]
    for r in (client.get(f"/tools/{name}", params=_params(bad)), client.post(f"/tools/{name}", json=bad)):
        assert r.status_code == 200, (name, r.text[:300])
        env = r.json()
        assert env["status"] == UNKNOWN, (name, env["status"])
        assert env.get("reason"), name
        assert r.headers["x-tool-status"] == UNKNOWN


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_rank_route_returns_the_spec_schema(client):
    r = client.get("/tools/rank_deployment_opportunities", params={
        "condition": "Long COVID or ME/CFS", "measurement": "wearable autonomic monitoring", "top_n": 2,
        "weight_set": "equal"})
    env = r.json()
    rec = env["data"]["recommendations"][0]
    for k in ("condition", "phenotype", "measurement", "technology", "geography", "candidate_sites",
              "research_evidence", "molecular_context", "uncertainties", "provenance", "recommended_next_step"):
        assert k in rec, k
    g = rec["geography"]
    assert g["burden_evidence_level"] in ("A", "B", "C") and g["rank_interval_5_95"][0] <= g["rank"]
    assert env["provenance"]["object_ids"][0].startswith(("opportunity:", "condition:", "geo:", "measurement:"))


# ------------------------------------------------------------------------------------------------ sources

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_sources(client):
    r = client.get("/sources")
    env = r.json()
    assert r.status_code == 200 and env["status"] == "ok" and env["tool"] == "list_sources"
    n = env["data"]["n_sources"]
    assert n == len(env["data"]["sources"]) >= 20
    full = client.get("/sources", params={"full": "true"}).json()
    assert len(full["data"]["sources"]) == n and full["data"]["detail"].startswith("full")
    assert len(full["data"]["sources"][0]) >= len(env["data"]["sources"][0])
    person = client.get("/sources", params={"data_layer": "person"}).json()
    assert person["status"] == "ok" and all(s["data_layer"] == "person" for s in person["data"]["sources"])
    bad = client.get("/sources", params={"data_layer": NONSENSE}).json()
    assert bad["status"] == UNKNOWN and bad["reason"]
    one = client.get("/sources/nppes")
    assert one.status_code == 200 and one.json()["source"]["source_id"] == "nppes"
    missing = client.get("/sources/zzqxv")
    assert missing.status_code == 404 and missing.json()["detail"]["status"] == UNKNOWN


# ------------------------------------------------------------------------------------------------ geojson

@pytest.mark.data
@pytest.mark.skipif(not HAS_BOUNDARIES, reason="needs the boundary geoparquets")
@pytest.mark.parametrize("level,width,n_min", [("county", 5, 3100), ("state", 2, 51)])
def test_geojson(client, level, width, n_min):
    r = client.get(f"/geojson/{level}")
    assert r.status_code == 200 and r.headers["content-type"].startswith("application/geo+json")
    fc = r.json()
    assert fc["type"] == "FeatureCollection" and len(fc["features"]) >= n_min
    assert fc["metadata"]["source"]["source_id"] == "census_geography" and "Not for spatial joins" in \
        fc["metadata"]["caveat"]
    ids = [f["id"] for f in fc["features"]]
    assert len(ids) == len(set(ids))
    for f in fc["features"][:200]:
        assert f["id"] == f["properties"]["geo_id"] and len(f["id"]) == width and f["id"].isdigit()
        assert f["properties"]["object_id"] == f"geo:{f['id']}"
        assert f["geometry"]["type"] in ("Polygon", "MultiPolygon")
    coords = re.findall(r"-?\d+\.\d+", json.dumps(fc["features"][0]["geometry"]))
    assert coords and all(len(c.split(".")[1]) <= 3 for c in coords)       # snapped to the 0.001-degree grid
    # cached: identical bytes, ETag -> 304
    r2 = client.get(f"/geojson/{level}", headers={"If-None-Match": r.headers["etag"]})
    assert r2.status_code == 304
    assert client.get(f"/geojson/{level}").content == r.content


@pytest.mark.data
@pytest.mark.skipif(not HAS_BOUNDARIES, reason="needs the boundary geoparquets")
def test_geojson_state_filter_and_simplification(client):
    ca = client.get("/geojson/county", params={"state": "06"}).json()
    assert len(ca["features"]) == 58 and all(f["id"].startswith("06") for f in ca["features"])
    coarse = client.get("/geojson/county", params={"state": "06", "tolerance": 0.05}).content
    fine = client.get("/geojson/county", params={"state": "06", "tolerance": 0}).content
    assert len(coarse) < len(fine)
    assert "06073" in {f["id"] for f in ca["features"]}
    universe = {f["id"]: f["properties"]["in_ranking_universe"] for f in G.boundaries_geojson("state")["features"]}
    assert universe["06"] is True and universe.get("72") is False
