"""MCP server (FastMCP, stdio) exposing the "Measure It to Cure It" public-data tools.

    uv run measure-it serve-mcp                 # stdio server (what .mcp.json starts)
    uv run python -m measure_it.mcp.server      # same

Every tool is a thin typed wrapper around the function of the same name in ``measure_it.tools`` (one shared facade
for the MCP server, the deterministic agent and the tests). Each returns the JSON envelope
``{tool, status, query, data, caveats, provenance, truncation, reason?, meta}``; ``status`` is ``ok``,
``UNKNOWN / NOT AVAILABLE`` or ``error``.

Also exposed:
* prompts   ``measure_it_guardrails`` (the rules a calling model must follow) and ``deployment_question`` (the
            documented tool order for "what should we measure, where, and who could deploy it?");
* resources ``measure-it://guardrails``, ``measure-it://data-layers``, ``measure-it://object-ids`` (namespaces
            trace_evidence resolves) and ``measure-it://tools`` (tool -> wrapped function).
"""
from __future__ import annotations

import functools
import inspect
import json
from typing import Annotated, Literal

from pydantic import Field

from .. import tools as T
from ..config import UNKNOWN
from .trace import namespaces_doc

SERVER_NAME = "measure-it-public"

GUARDRAILS_MD = f"""# Measure It to Cure It (public-data alpha): rules for any model calling these tools

1. **Never add facts beyond tool output.** Every number, place, facility, trial, grant, device, gene or effect size
   you state must appear in a tool result. Do not estimate, extrapolate, round into new claims or fill gaps from
   general knowledge. No model-generated number is data.
2. **UNKNOWN stays UNKNOWN.** When a tool returns status `{UNKNOWN}` (or a field holds that value), say it is
   unknown / not available and give the tool's `reason`. Do not guess a value.
3. **Cite object ids.** Tool results carry object ids (`condition:`, `measurement:`, `opportunity:`, `geo:`,
   `facility:`, `trial:`, `nih:`, `signature:` ...). Cite them next to each fact; `trace_evidence(object_id)` resolves
   any of them to its rows, sources, DATA_AUDIT.md and raw files.
4. **Four data layers are never joined as if they were the same people.** Person-level wearable/clinical data
   (NHANES, Stanford wearable studies, MapMECFS) are summarised as group statistics; they are never located and never
   linked to a place, facility or trial. Geographic joins are ecological (places, not people).
5. **Condition-level molecular enrichment** (Open Targets, GWAS Catalog, GEO/SRA, mapMECFS supplements) is not measured
   on any participant: call it condition-level molecular enrichment, never patient or personal multi-omics, and never
   a mechanism.
6. **Burden evidence levels.** Always state `burden_evidence_level` with a burden value: A = direct condition measure,
   B = closely matching coded condition, C = symptom/comorbidity proxy, D = no usable estimate (burden UNKNOWN). A state
   value attached to a county (`inherited` = true, source resolution `state`) is inherited context, never county
   prevalence.
7. **Research and regulatory signals are not efficacy.** A registered trial or NIH grant shows research activity; an
   FDA record is a deployment-readiness signal. Neither shows that a measurement works or detects a condition.
8. **Product language.** Say candidate, deployment opportunity, measurement desert, research readiness, evidence-
   supported, measurable phenotype, implementation candidate. Never say "proves", "diagnosed by AI", "patient has",
   "definitive biomarker", "best clinic", "optimal treatment", "confirms mechanism" or "cures". A ranked region is a
   candidate deployment opportunity for pilot evaluation, not a validated diagnostic pathway; a facility "has
   characteristics suggesting it may be a viable implementation or study partner".
"""

DATA_LAYERS_MD = """# Data layers

| layer | what | examples | joins |
|---|---|---|---|
| person | person-level wearable, questionnaire, lab, diagnosis data | NHANES 2011-2014 accelerometry; Stanford smartwatch studies; MapMECFS participants | summarised only (signature:, measurement_signal:); never located, never joined to layers 3-4 |
| condition_molecular | condition-level molecular evidence | Open Targets, GWAS Catalog, GEO/SRA, mapMECFS supplements, Reactome | by condition concept only ("condition-level molecular enrichment") |
| geographic | population / place data | CDC HPS long COVID (state), CMS MMD (county, Medicare FFS), CDC PLACES, SVI, ACS, Lyme | ecological joins by FIPS; state values are never county prevalence |
| facility | providers, facilities, trial sites, NIH awardees, devices | NPPES, HRSA, ClinicalTrials.gov, NIH RePORTER, openFDA | by location (ZIP/city centroids) and registry ids |
| ontology / measurement / derived | condition registry, measurement registry, scores | Mondo/EFO/HPO/ICD-10-CM; measurement classes; deployment_opportunities | concept bridges and ecological scores; never person-level |

The ontology bridges layers by concept (condition -> phenotype -> measurement), never by person.
"""


def _tools_md() -> str:
    lines = ["# Tools", "", "| tool | wraps |", "|---|---|"]
    for n in T.SPEC_TOOLS + T.EXTRA_TOOLS:
        lines.append(f"| {n} | {T.TOOL_PRODUCERS.get(n)} |")
    return "\n".join(lines)


def _object_ids_md() -> str:
    lines = ["# Object-id namespaces resolved by trace_evidence", "", "| namespace | table | native id |",
             "|---|---|---|"]
    for r in namespaces_doc():
        lines.append(f"| {r['namespace']} | {r['table']} | {r['native_id']} |")
    return "\n".join(lines)


INSTRUCTIONS = (
    "Public-data engine that answers: what objective measurement could make an invisible-illness phenotype visible, "
    "where in the U.S. deploying it is a candidate opportunity, and which facilities or research sites could "
    "realistically implement or study it. Rules: use only facts present in tool output, never add facts beyond it; "
    "report 'UNKNOWN / NOT AVAILABLE' as unknown with the tool's reason; cite the object ids from tool output; state "
    "burden_evidence_level with every burden value; person-level data are group statistics only and never linked to "
    "places; molecular evidence is condition-level molecular enrichment, not patient multi-omics; say 'candidate "
    "deployment opportunity', never 'validated'. Suggested order for a deployment question: normalize_condition -> "
    "get_patient_phenotype_signature -> discover_candidate_measurements -> get_measurement_evidence -> "
    "get_regulatory_context -> get_condition_burden -> rank_deployment_opportunities -> find_candidate_clinics / "
    "find_relevant_research_centers -> trace_evidence. Read the measure_it_guardrails prompt or the "
    "measure-it://guardrails resource for the full rules.")

TOOL_PREAMBLE = ("Returns the JSON envelope {status: 'ok' | 'UNKNOWN / NOT AVAILABLE' | 'error', data, caveats, "
                 "provenance (object ids + sources), truncation}. Use ONLY facts in the returned data; never add "
                 "facts beyond tool output; if status is 'UNKNOWN / NOT AVAILABLE', report it as unknown with the "
                 "reason. Cite provenance.object_ids; trace_evidence resolves them.")

# ------------------------------------------------------------------------------------------------ parameter docs

COND = ("Condition name, alias, acronym, ICD-10-CM code or canonical id, e.g. 'Long COVID', 'ME/CFS', 'POTS', "
        "'dysautonomia', 'long_covid'. Resolved with the ontology normaliser; ambiguous text returns UNKNOWN with "
        "candidates.")
MEAS = ("Measurement class id/name, bundle id/label or alias, e.g. 'wearable autonomic monitoring', "
        "'nailfold capillaroscopy', 'autonomic testing', 'CPET', 'accelerometry', 'wearable_heart_rate'.")
PARAM_DOCS = {
    "condition": COND,
    "condition_or_code": "Free text, alias, ICD-10-CM code (G93.32), or a MONDO / EFO / MeSH / HPO id.",
    "limit": "Maximum number of candidates (1-25).",
    "top_n": "Maximum number of items returned.",
    "include_broad": "Include rows retrieved by broader ontology ids (Open Targets / GWAS Catalog).",
    "phenotype": ("Optional phenotype: a phenotype axis (e.g. 'orthostatic intolerance', 'post-exertional malaise', "
                  "HPO id) or a demo phenotype ('orthostatic/autonomic dysfunction', 'activity intolerance')."),
    "measurement": MEAS,
    "max_object_ids": "Maximum trial / grant object ids returned per measurement.",
    "technology": ("Technology or device name, e.g. 'wearable ECG patch', 'Zio', 'EndoPAT', 'nailfold "
                   "capillaroscopy', 'pulse oximeter'. Questionnaires are refused (not objective measurements)."),
    "max_records": "Maximum matching FDA device records returned.",
    "geography_level": "Geography level.",
    "geo_id": "Optional single geography (FIPS, 'geo:<fips>' or a name) to return only its row.",
    "include_alternates": "Also return the county proxy and pre-specified alternate measures (labelled by role).",
    "max_rows": "Maximum rows returned (a summary of all rows is always included).",
    "order": "Row order: 'geo_id' or 'value_desc' (presentation only).",
    "geography_id": ("State or county: FIPS ('06073', '06'), 'geo:<fips>', 'zcta:<ZCTA>', a state name/abbreviation "
                     "or a county name such as 'San Diego County, California'."),
    "radius_km": "Search radius around the geography's internal point (default 50 km).",
    "max_sites": "Maximum candidate facilities returned (default 15).",
    "literal_only": "Only trials whose record names the condition literally (precision view; default true).",
    "statuses": "Optional ClinicalTrials.gov overall_status filter, e.g. ['RECRUITING'].",
    "us_sites_only": "Only trials with at least one U.S. site.",
    "max_trials": "Maximum trials returned (a summary of all is included).",
    "include_expanded": "Include MeSH/terms-expanded matches (lower precision).",
    "weight_set": ("Named weight set from configs/scoring.yaml: equal (default), burden_led, equity_led, capacity_led, "
                   "study_partner, access_gap, burden_only, saturation_adjusted."),
    "include_small_population": "Also rank counties below the minimum population (unstable rates).",
    "max_sites_per_region": "Candidate facilities shown per region in compact mode.",
    "detail": "'compact' (shared context once, capped sites) or 'full' (everything, large).",
    "query": "Free-text Data.gov catalog query, e.g. 'long COVID', 'social vulnerability index'.",
    "rows": "Number of Data.gov results (default 20).",
    "object_id": ("Object id '<namespace>:<native id>', e.g. 'opportunity:long_covid_or_me_cfs|"
                  "wearable_autonomic_activity_monitoring|40109', 'facility:npi-1538213764', 'trial:NCT05172024', "
                  "'geo:06073', 'signature:nhanes|mecfs_like_proxy|interdaily_stability'."),
    "include_unsupported": "Also return tested-negative (not supported) rows.",
    "data_layer": "Optional layer filter: person, condition_molecular, geographic, facility, ontology, metadata.",
}
PER_TOOL_DOCS = {
    ("get_patient_phenotype_signature", "condition"): (
        "Condition or phenotype text, e.g. 'Long COVID', 'ME/CFS', 'depression', 'mortality', 'sleep disorder'."),
    ("rank_deployment_opportunities", "condition"): (
        COND + " Also condition sets: 'Long COVID or ME/CFS', 'demo cluster', or any 'A or B' of conditions."),
    ("find_relevant_trials", "condition"): COND + " Condition-set ids (e.g. long_covid_or_me_cfs) are expanded.",
    ("find_candidate_clinics", "measurement"): "Optional. " + MEAS,
    ("find_relevant_research_centers", "measurement"): "Optional. " + MEAS,
    ("find_relevant_trials", "measurement"): "Optional. " + MEAS + " Only trials whose registered text mentions it.",
    ("get_measurable_biology", "measurement"): "Optional measurement class or bundle id/alias filter.",
    ("rank_deployment_opportunities", "geography_level"): "Ranking level: 'county' (default) or 'state'.",
    ("get_condition_burden", "geography_level"): "'national', 'state' (default) or 'county'.",
}
TYPE_OVERRIDES = {
    ("get_condition_burden", "geography_level"): Literal["national", "state", "county"],
    ("get_condition_burden", "order"): Literal["geo_id", "value_desc"],
    ("rank_deployment_opportunities", "geography_level"): Literal["county", "state"],
    ("rank_deployment_opportunities", "detail"): Literal["compact", "full"],
}
TITLES = {
    "search_condition": "Search conditions", "normalize_condition": "Normalize a condition or code",
    "get_patient_phenotype_signature": "Wearable phenotype signature (public person-level data, group statistics)",
    "get_molecular_context": "Condition-level molecular enrichment",
    "discover_candidate_measurements": "Candidate objective measurements",
    "get_measurement_evidence": "Measurement evidence for a condition",
    "get_regulatory_context": "FDA regulatory context (deployment-readiness signal)",
    "get_condition_burden": "Condition burden by geography (with evidence level)",
    "get_geographic_context": "Geographic context", "find_candidate_clinics": "Candidate implementation partners",
    "find_relevant_trials": "Relevant ClinicalTrials.gov trials",
    "find_relevant_research_centers": "Research readiness by facility",
    "rank_deployment_opportunities": "Rank candidate deployment opportunities",
    "search_us_open_data": "Search Data.gov (metadata catalog)", "trace_evidence": "Trace an object id to sources",
    "get_phenotype_measurement_evidence": "Phase 3 phenotype measurement evidence",
    "get_measurable_biology": "Measurable biology (molecular support per measurement class)",
    "list_sources": "List public data sources",
}


def _typed_wrapper(name: str, fn):
    """Same signature as the facade function, with pydantic Field descriptions (and enums) for every parameter."""
    sig = inspect.signature(fn)
    params, ann = [], {}
    for p in sig.parameters.values():
        t = TYPE_OVERRIDES.get((name, p.name), p.annotation)
        doc = PER_TOOL_DOCS.get((name, p.name), PARAM_DOCS.get(p.name, p.name))
        a = Annotated[t, Field(description=doc)]
        params.append(p.replace(annotation=a))
        ann[p.name] = a
    ann["return"] = dict

    @functools.wraps(fn)
    def wrapper(**kwargs):
        return fn(**kwargs)

    wrapper.__signature__ = sig.replace(parameters=params, return_annotation=dict)
    wrapper.__annotations__ = ann
    del wrapper.__wrapped__
    return wrapper


def _description(name: str, fn) -> str:
    doc = inspect.getdoc(fn) or name
    return f"{doc}\n\nWraps {T.TOOL_PRODUCERS.get(name)}. {TOOL_PREAMBLE}"


BYOD_TOOLS = {"list_user_datasets": "List user-supplied (BYOD) datasets on this machine",
              "get_user_dataset_metric": "Locked plan, metrics and performance record of a user dataset"}


def _register_byod_tools(mcp) -> None:
    """User-supplied datasets (measure_it.byod; docs/BRING_YOUR_OWN_DATA.md): aggregate results only, read-only.
    Kept outside tools.TOOLS so the SPEC/extra tool list, the API routes and their counts are unchanged."""
    from ..byod import tools as BT
    PER_TOOL_DOCS.setdefault(("get_user_dataset_metric", "dataset_id"),
                             "user dataset id, with or without the byod_ prefix (see list_user_datasets)")
    for name, title in BYOD_TOOLS.items():
        fn = getattr(BT, name)
        mcp.tool(_typed_wrapper(name, fn), name=name, title=title,
                 description=f"{inspect.getdoc(fn) or name}\n\nWraps measure_it.byod.tools. {TOOL_PREAMBLE}",
                 annotations={"title": title, "readOnlyHint": True, "destructiveHint": False, "idempotentHint": True,
                              "openWorldHint": False}, tags={"byod"})


def build_server():
    from fastmcp import FastMCP

    mcp = FastMCP(name=SERVER_NAME, instructions=INSTRUCTIONS, version="0.1.0")
    for name in T.SPEC_TOOLS + T.EXTRA_TOOLS:
        fn = T.TOOLS[name]
        mcp.tool(_typed_wrapper(name, fn), name=name, title=TITLES.get(name), description=_description(name, fn),
                 annotations={"title": TITLES.get(name), "readOnlyHint": True, "destructiveHint": False,
                              "idempotentHint": True, "openWorldHint": name == "search_us_open_data"},
                 tags={"spec" if name in T.SPEC_TOOLS else "extra"})
    _register_byod_tools(mcp)

    @mcp.prompt(name="measure_it_guardrails", title="Guardrails for answering with measure-it tools")
    def measure_it_guardrails() -> str:
        """The rules any model must follow when it answers with these tools (no facts beyond tool output, UNKNOWN
        stays UNKNOWN, object-id citations, data layers, evidence levels, product language)."""
        return GUARDRAILS_MD + "\n" + DATA_LAYERS_MD

    @mcp.prompt(name="deployment_question", title="Answer a deployment question with the tools")
    def deployment_question(question: str) -> str:
        """Tool plan for 'what should we measure, where should we measure it, and who could realistically deploy
        it?' questions."""
        return (
            GUARDRAILS_MD + "\n\nAnswer this question using only tool output:\n\n" + str(question) + "\n\n"
            "Tool order (the deterministic agent `uv run measure-it ask` uses the same order):\n"
            "1. normalize_condition for every condition or phenotype phrase (search_condition when nothing matches).\n"
            "2. get_phenotype_measurement_evidence when a phenotype is named.\n"
            "3. get_patient_phenotype_signature per condition (group statistics from public person-level data).\n"
            "4. discover_candidate_measurements per condition (phenotype filter when given); "
            "get_measurement_evidence(measurement, condition); get_regulatory_context(measurement).\n"
            "5. get_molecular_context per condition (condition-level molecular enrichment only).\n"
            "6. get_condition_burden(condition, level) - state the burden_evidence_level; D means burden UNKNOWN.\n"
            "7. rank_deployment_opportunities(condition, measurement, 'county') - candidate deployment opportunities.\n"
            "8. find_candidate_clinics(top region, condition, measurement) and find_relevant_research_centers.\n"
            "9. trace_evidence on the top opportunity id to show the provenance chain.\n"
            "Write three sections (what to measure / where / who could deploy it) plus uncertainties; cite object "
            "ids after each factual sentence; keep UNKNOWN parts UNKNOWN.")

    @mcp.resource("measure-it://guardrails", name="guardrails", mime_type="text/markdown",
                  description="Rules for any model answering with these tools.")
    def guardrails_resource() -> str:
        return GUARDRAILS_MD

    @mcp.resource("measure-it://data-layers", name="data_layers", mime_type="text/markdown",
                  description="The four data layers (+ ontology/measurement/derived) and how they may be joined.")
    def data_layers_resource() -> str:
        return DATA_LAYERS_MD

    @mcp.resource("measure-it://object-ids", name="object_ids", mime_type="text/markdown",
                  description="Object-id namespaces that trace_evidence resolves (docs/CONVENTIONS.md section 7).")
    def object_ids_resource() -> str:
        return _object_ids_md()

    @mcp.resource("measure-it://tools", name="tools", mime_type="text/markdown",
                  description="Each MCP tool and the query function it wraps.")
    def tools_resource() -> str:
        return _tools_md()

    return mcp


@functools.lru_cache(maxsize=1)
def get_server():
    return build_server()


def serve(transport: str = "stdio", host: str | None = None, port: int | None = None) -> None:
    """Run the server (stdio by default). The MCP SDK points fd 1 at stderr while serving, so stray prints from
    wrapped modules never reach the JSON-RPC stream.

    transport "http" (streamable HTTP, endpoint /mcp) or "sse" serves over the network instead (the container
    deployment, docs/DEPLOYMENT.md) and adds GET /health (liveness + data stamp) for container health checks."""
    srv = get_server()
    if transport == "stdio":
        srv.run(transport=transport, show_banner=False)
        return
    _add_health_route(srv)
    srv.run(transport=transport, show_banner=False, host=host or "127.0.0.1", port=port or 8765)


def _add_health_route(srv) -> None:
    if getattr(srv, "_measure_it_health", False):
        return
    from starlette.responses import JSONResponse

    from ..config import PROCESSED, SOURCE_REGISTRY_PATH

    @srv.custom_route("/health", methods=["GET"])
    async def health(request):  # noqa: ARG001 - starlette handler signature
        n_tables = len(list(PROCESSED.glob("*.parquet"))) if PROCESSED.exists() else 0
        ok = n_tables > 0 and SOURCE_REGISTRY_PATH.exists()
        return JSONResponse({"status": "ok" if ok else "degraded", "service": SERVER_NAME,
                             "reason": None if ok else "data/processed or SOURCE_REGISTRY.yaml missing",
                             "n_tools": len(T.SPEC_TOOLS + T.EXTRA_TOOLS), "processed_tables": n_tables,
                             "data_stamp": T.data_stamp()})

    srv._measure_it_health = True


def tool_manifest() -> str:
    return json.dumps({n: T.TOOL_PRODUCERS.get(n) for n in T.SPEC_TOOLS + T.EXTRA_TOOLS}, indent=1)


if __name__ == "__main__":  # pragma: no cover
    serve()
