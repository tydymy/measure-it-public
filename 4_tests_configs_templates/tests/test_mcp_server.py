"""Tests for the MCP server (measure_it.mcp.server) through the FastMCP in-memory client.

The client talks to the server object in-process (no subprocess, no network). Tests marked `data` need
data/processed. Data.gov calls run with MEASURE_IT_OFFLINE=1 (answered from the HTTP cache or UNKNOWN).
"""
from __future__ import annotations

import asyncio
import json

import pytest

from measure_it import tools as T
from measure_it.config import UNKNOWN
from measure_it.validate import forbidden_hits

from test_tools import CASES  # same demo / nonsense inputs as the facade tests

fastmcp = pytest.importorskip("fastmcp")


def _server():
    from measure_it.mcp.server import build_server
    return build_server()


def run(coro):
    return asyncio.run(coro)


async def _with_client(fn):
    from fastmcp import Client
    async with Client(_server()) as c:
        return await fn(c)


@pytest.fixture
def offline(monkeypatch):
    monkeypatch.setenv("MEASURE_IT_OFFLINE", "1")


# ------------------------------------------------------------------------------------------------ listing

def test_tool_list_contains_every_spec_tool_with_guardrail_docs():
    tools = run(_with_client(lambda c: c.list_tools()))
    names = {t.name for t in tools}
    assert set(T.SPEC_TOOLS) <= names
    assert {"get_phenotype_measurement_evidence", "get_measurable_biology", "list_sources"} <= names
    for t in tools:
        assert "never add facts beyond tool output" in t.description, t.name
        assert "UNKNOWN / NOT AVAILABLE" in t.description, t.name
        props = t.input_schema.get("properties", {})
        for p, spec in props.items():
            assert spec.get("description"), (t.name, p)
        assert t.annotations.read_only_hint is True
    rank = next(t for t in tools if t.name == "rank_deployment_opportunities")
    assert rank.input_schema["properties"]["geography_level"]["enum"] == ["county", "state"]
    assert rank.input_schema["required"] == ["condition", "measurement"]


def test_prompts_and_resources_describe_guardrails_and_layers():
    async def go(c):
        prompts = {p.name for p in await c.list_prompts()}
        res = {str(r.uri) for r in await c.list_resources()}
        g = await c.get_prompt("measure_it_guardrails")
        d = await c.get_prompt("deployment_question", {"question": "Where should we deploy CPET for ME/CFS?"})
        layers = await c.read_resource("measure-it://data-layers")
        ids = await c.read_resource("measure-it://object-ids")
        return prompts, res, g.messages[0].content.text, d.messages[0].content.text, layers[0].text, ids[0].text
    prompts, res, g, d, layers, ids = run(_with_client(go))
    assert {"measure_it_guardrails", "deployment_question"} <= prompts
    assert {"measure-it://guardrails", "measure-it://data-layers", "measure-it://object-ids",
            "measure-it://tools"} <= res
    assert "Never add facts beyond tool output" in g and "UNKNOWN stays UNKNOWN" in g
    assert "condition-level molecular enrichment" in g and "burden_evidence_level" in g
    assert "CPET for ME/CFS" in d and "rank_deployment_opportunities" in d
    assert "person" in layers and "ecological" in layers
    assert "participant" in ids and "opportunity" in ids
    assert not forbidden_hits(g) and not forbidden_hits(d)     # forbidden terms appear only quoted


# ------------------------------------------------------------------------------------------------ calls

async def _call(c, tool, args):
    r = await c.call_tool(tool, args, raise_on_error=False)
    assert not r.is_error, (tool, r.content[0].text[:300] if r.content else "")
    return r.structured_content


@pytest.mark.data
def test_every_tool_ok_on_demo_and_unknown_on_nonsense(offline):
    async def go(c):
        out = {}
        for tool, (demo, bad) in CASES.items():
            out[tool] = (await _call(c, tool, demo), await _call(c, tool, bad))
        return out
    res = run(_with_client(go))
    for tool, (ok, bad) in res.items():
        assert ok["status"] == "ok", (tool, ok.get("reason"))
        assert bad["status"] == UNKNOWN and bad["reason"], (tool, bad["status"])
        for r in (ok, bad):
            json.dumps(r, allow_nan=False)
            hits = forbidden_hits(json.dumps(r, indent=1))
            assert not hits, (tool, hits[:3])


@pytest.mark.data
def test_trace_round_trips_every_object_id_of_a_ranking():
    async def go(c):
        rank = await _call(c, "rank_deployment_opportunities",
                           {"condition": "Long COVID or ME/CFS", "measurement": "wearable autonomic monitoring"})
        ids = T.collect_object_ids(rank)
        traces = {}
        for oid in ids:
            traces[oid] = await _call(c, "trace_evidence", {"object_id": oid})
        return rank, ids, traces
    rank, ids, traces = run(_with_client(go))
    assert rank["status"] == "ok" and len(rank["data"]["recommendations"]) == 10
    namespaces = {i.split(":", 1)[0] for i in ids}
    assert {"opportunity", "geo", "condition", "measurement", "facility", "trial", "signature"} <= namespaces
    bad = {oid: t.get("reason") for oid, t in traces.items() if t["status"] != "ok"
           or t["data"]["object_id"] != oid or not t["data"]["rows"]}
    assert not bad, list(bad.items())[:5]
    # the opportunity trace walks back to burden rows and their raw source files
    opp = next(i for i in ids if i.startswith("opportunity:"))
    t = traces[opp]["data"]
    assert any(lk.get("table") == "geo_condition_burden" for lk in t["lineage"])
    assert any(s["raw"]["record_files"] for s in t["sources"])


@pytest.mark.data
def test_person_level_trace_refused_through_mcp():
    r = run(_with_client(lambda c: _call(c, "trace_evidence", {"object_id": "participant:nhanes:62161"})))
    assert r["status"] == UNKNOWN and r["data"]["refused"] is True


@pytest.mark.data
def test_invalid_enum_is_rejected_by_schema():
    async def go(c):
        return await c.call_tool("get_condition_burden", {"condition": "ME/CFS", "geography_level": "galaxy"},
                                 raise_on_error=False)
    r = run(_with_client(go))
    assert r.is_error
