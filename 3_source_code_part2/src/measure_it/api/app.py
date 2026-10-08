"""FastAPI app for the "Measure It to Cure It" public-data engine: every agentic-layer tool over HTTP.

    uv run measure-it serve-api                      # http://127.0.0.1:8000 (OpenAPI docs at /docs)
    uv run uvicorn measure_it.api.app:app --port 8000

Routes
    GET  /health                 liveness + data stamp, processed-table count, offline flag
    GET  /tools                  every tool: title, wrapped function, parameters (type, default, enum, description)
    GET  /tools/{tool}           one route per tool; arguments as query parameters
    POST /tools/{tool}           the same tool; arguments as a JSON object
    GET  /sources                SOURCE_REGISTRY.yaml summary (the list_sources tool); ?full=true for whole entries
    GET  /sources/{source_id}    one full registry entry
    GET  /geojson/{level}        simplified county / state boundaries (GeoJSON, cached, ETag)

Every tool route returns exactly the envelope of the same function in `measure_it.tools` (the MCP server and the
dashboard call the same functions): {tool, status, query, data, caveats, provenance, truncation, reason?, meta}.
HTTP status: 200 for status "ok" and for "UNKNOWN / NOT AVAILABLE" (a valid answer: the data are missing, with the
reason), 500 for status "error" (same body), 422 when an argument fails validation before the tool runs. The
`X-Tool-Status` header repeats the envelope status. Nothing here computes a result: the routes only validate
arguments and serialise the tool output.
"""
from __future__ import annotations

import inspect
import os
import time
import typing
from typing import Annotated, Any, Literal, Optional

from fastapi import Body, FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, Response
from pydantic import ConfigDict, Field, create_model

from .. import tools as T
from ..config import PROCESSED, SOURCE_REGISTRY_PATH, UNKNOWN
from ..mcp.server import PARAM_DOCS, PER_TOOL_DOCS, TITLES, TYPE_OVERRIDES
from ..mcp.trace import registry, to_jsonable
from . import geo as G

VERSION = "0.1.0"
TOOL_ORDER = T.SPEC_TOOLS + T.EXTRA_TOOLS
DESCRIPTION = """Public-data engine that answers: what objective measurement could make an invisible-illness phenotype
visible, where in the U.S. deploying it is a **candidate deployment opportunity**, and which facilities or research
sites could realistically implement or study it.

Every `/tools/{tool}` route returns the JSON envelope of the same function in `measure_it.tools` (shared with the MCP
server and the dashboard). `status` is `ok`, `UNKNOWN / NOT AVAILABLE` (data missing; `reason` says why; nothing is
filled in) or `error`. Use only facts in `data`; cite `provenance.object_ids`; `trace_evidence` resolves any of them.
Person-level data are group statistics only and never joined to places or facilities; molecular evidence is
condition-level molecular enrichment (never patient multi-omics); burden values carry `burden_evidence_level` and a
state value attached to a county is inherited context, not county prevalence."""

app = FastAPI(title="Measure It to Cure It: public-data API", version=VERSION, description=DESCRIPTION,
              openapi_tags=[{"name": "service"}, {"name": "spec tools", "description": "the SPEC agentic-layer tools"},
                            {"name": "extra tools"}, {"name": "sources"}, {"name": "geography"}])
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["GET", "POST"], allow_headers=["*"],
                   expose_headers=["X-Tool-Status", "ETag"])


# ------------------------------------------------------------------------------------------------ helpers

def _respond(env: dict) -> JSONResponse:
    status = env.get("status")
    return JSONResponse(content=env, status_code=500 if status == "error" else 200,
                        headers={"X-Tool-Status": str(status)})


def _param_types(name: str) -> dict[str, Any]:
    hints = typing.get_type_hints(T.TOOLS[name])
    return {p: TYPE_OVERRIDES.get((name, p), t) for p, t in hints.items() if p != "return"}


def _doc(name: str, param: str) -> str:
    return PER_TOOL_DOCS.get((name, param), PARAM_DOCS.get(param, param))


def _type_name(t: Any) -> str:
    if typing.get_origin(t) is Literal:
        return "enum"
    s = getattr(t, "__name__", None) or str(t)
    return s.replace("typing.", "")


def tool_spec(name: str) -> dict:
    """Machine-readable description of one tool (what GET /tools lists)."""
    fn = T.TOOLS[name]
    sig = inspect.signature(fn)
    types = _param_types(name)
    params = []
    for p in sig.parameters.values():
        t = types.get(p.name)
        item = {"name": p.name, "type": _type_name(t), "required": p.default is inspect.Parameter.empty,
                "default": None if p.default is inspect.Parameter.empty else p.default,
                "description": _doc(name, p.name)}
        if typing.get_origin(t) is Literal:
            item["enum"] = list(typing.get_args(t))
        params.append(item)
    return {"name": name, "title": TITLES.get(name), "kind": "spec" if name in T.SPEC_TOOLS else "extra",
            "description": inspect.getdoc(fn), "wraps": T.TOOL_PRODUCERS.get(name), "parameters": params,
            "routes": {"get": f"/tools/{name}", "post": f"/tools/{name}"}}


def _get_endpoint(name: str):
    """GET handler whose signature carries every tool argument as a typed, documented query parameter."""
    fn = T.TOOLS[name]
    types = _param_types(name)
    params = []
    for p in inspect.signature(fn).parameters.values():
        default = p.default   # inspect.Parameter.empty for a required argument
        params.append(inspect.Parameter(p.name, inspect.Parameter.KEYWORD_ONLY, default=default,
                                        annotation=Annotated[types[p.name], Query(description=_doc(name, p.name))]))

    def endpoint(**kwargs) -> JSONResponse:
        return _respond(fn(**kwargs))

    endpoint.__signature__ = inspect.Signature(params, return_annotation=JSONResponse)
    endpoint.__name__ = f"{name}_get"
    return endpoint


def _post_endpoint(name: str):
    """POST handler taking the tool arguments as one JSON object (unknown keys rejected)."""
    fn = T.TOOLS[name]
    types = _param_types(name)
    fields, required = {}, False
    for p in inspect.signature(fn).parameters.values():
        if p.default is inspect.Parameter.empty:
            required = True
            fields[p.name] = (types[p.name], Field(..., description=_doc(name, p.name)))
        else:
            fields[p.name] = (types[p.name], Field(p.default, description=_doc(name, p.name)))
    model = create_model(f"{''.join(w.title() for w in name.split('_'))}Request",
                         __config__=ConfigDict(extra="forbid"), **fields)

    def endpoint(body=None) -> JSONResponse:
        return _respond(fn(**(body.model_dump() if body is not None else {})))

    # real annotation objects (this module uses postponed annotations, which FastAPI cannot resolve for a local model)
    ann = Annotated[model, Body()] if required else Annotated[Optional[model], Body()]
    endpoint.__signature__ = inspect.Signature(
        [inspect.Parameter("body", inspect.Parameter.KEYWORD_ONLY, annotation=ann,
                           default=inspect.Parameter.empty if required else None)], return_annotation=JSONResponse)
    endpoint.__name__ = f"{name}_post"
    return endpoint, model


for _name in TOOL_ORDER:
    _tag = "spec tools" if _name in T.SPEC_TOOLS else "extra tools"
    _desc = (inspect.getdoc(T.TOOLS[_name]) or _name) + f"\n\nWraps `{T.TOOL_PRODUCERS.get(_name)}`. " + T.GENERAL_CAVEAT
    app.add_api_route(f"/tools/{_name}", _get_endpoint(_name), methods=["GET"], tags=[_tag], name=f"{_name}_get",
                      summary=TITLES.get(_name, _name), description=_desc, response_class=JSONResponse)
    _post, _ = _post_endpoint(_name)
    app.add_api_route(f"/tools/{_name}", _post, methods=["POST"], tags=[_tag], name=f"{_name}_post",
                      summary=f"{TITLES.get(_name, _name)} (JSON body)", description=_desc,
                      response_class=JSONResponse)


# ------------------------------------------------------------------------------------------------ service routes

@app.get("/", tags=["service"], summary="Index")
def index() -> dict:
    return {"service": app.title, "version": VERSION, "docs": "/docs", "openapi": "/openapi.json",
            "health": "/health", "tools": "/tools", "sources": "/sources", "geojson": ["/geojson/county",
                                                                                      "/geojson/state"]}


@app.get("/health", tags=["service"], summary="Liveness and data stamp")
def health() -> dict:
    n_tables = len(list(PROCESSED.glob("*.parquet"))) if PROCESSED.exists() else 0
    stamp = T.data_stamp()
    ok = n_tables > 0 and SOURCE_REGISTRY_PATH.exists()
    return {"status": "ok" if ok else "degraded",
            "reason": None if ok else "data/processed or SOURCE_REGISTRY.yaml missing: run the pipeline first",
            "service": app.title, "version": VERSION, "n_tools": len(TOOL_ORDER), "processed_tables": n_tables,
            "data_stamp": stamp,
            "data_stamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(stamp)) if stamp else UNKNOWN,
            "offline": os.environ.get("MEASURE_IT_OFFLINE") == "1", "unknown_sentinel": UNKNOWN}


@app.get("/tools", tags=["service"], summary="List every tool with its parameters")
def list_tools() -> dict:
    return {"n_tools": len(TOOL_ORDER), "spec_tools": T.SPEC_TOOLS, "extra_tools": T.EXTRA_TOOLS,
            "envelope": ["tool", "status", "query", "data", "caveats", "provenance", "truncation", "reason", "meta"],
            "status_values": ["ok", UNKNOWN, "error"], "tools": [tool_spec(n) for n in TOOL_ORDER]}


@app.get("/sources", tags=["sources"], summary="SOURCE_REGISTRY.yaml (the list_sources tool)")
def sources(data_layer: Annotated[Optional[str], Query(description=PARAM_DOCS["data_layer"])] = None,
            full: Annotated[bool, Query(description="Return whole registry entries instead of the summary.")] = False
            ) -> JSONResponse:
    env = T.list_sources(data_layer=data_layer)
    if full and env.get("status") == "ok":
        by_id = registry()["by_id"]
        env = {**env, "data": {**env["data"], "sources": [to_jsonable(by_id[s["source_id"]])
                                                          for s in env["data"]["sources"]],
                               "detail": "full SOURCE_REGISTRY.yaml entries"}}
    return _respond(env)


@app.get("/sources/{source_id}", tags=["sources"], summary="One SOURCE_REGISTRY.yaml entry")
def source(source_id: str) -> dict:
    s = registry()["by_id"].get(source_id)
    if s is None:
        raise HTTPException(status_code=404, detail={"status": UNKNOWN, "reason": f"no source {source_id!r} in "
                                                     "SOURCE_REGISTRY.yaml", "sources": sorted(registry()["by_id"])})
    return {"status": "ok", "source": to_jsonable(s), "data_audit": f"data/raw/{source_id}/DATA_AUDIT.md"}


@app.get("/geojson/{level}", tags=["geography"], summary="Simplified county or state boundaries (GeoJSON)",
         response_class=Response,
         responses={200: {"content": {"application/geo+json": {}}}, 304: {"description": "not modified (ETag)"}})
def geojson(request: Request, level: Literal["county", "state"],
            tolerance: Annotated[Optional[float], Query(ge=0, le=G.MAX_TOLERANCE_DEG, description=(
                "Simplification tolerance in degrees (default 0.01 county, 0.02 state; 0 = no simplification)"))]
            = None,
            state: Annotated[Optional[str], Query(description="Only features of this 2-digit state FIPS, e.g. '06'")]
            = None) -> Response:
    """Display geometry (Census 2024 cartographic boundaries, simplified; see measure_it.api.geo): feature id =
    FIPS; properties geo_id, name, state_abbr, object_id, in_ranking_universe. Not for spatial joins."""
    try:
        raw, etag = G.boundaries_geojson_bytes(level, tolerance, state)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e)) from e
    headers = {"ETag": f'"{etag}"', "Cache-Control": "public, max-age=86400"}
    if request.headers.get("if-none-match", "").strip('"') == etag:
        return Response(status_code=304, headers=headers)
    return Response(content=raw, media_type="application/geo+json", headers=headers)
