"""Shared tool facade: one plain-Python function per SPEC agentic-layer tool (plus a few extras).

Every function wraps an existing query function (it never re-derives a result) and returns one JSON-serialisable
envelope::

    {"tool": <name>, "status": "ok" | "UNKNOWN / NOT AVAILABLE" | "error", "query": {...arguments},
     "data": {...the wrapped function's result, cleaned},
     "caveats": [...], "provenance": {"object_ids": [...], "tables": [...], "sources": [...], "producer": ...},
     "truncation": [{"path", "returned", "total"}, ...],          # present whenever a list was capped
     "reason": "..." (when status is not ok), "meta": {"elapsed_s", "cache_hit", "data_stamp"}}

NaN / inf / NA become null. When the wrapped function has no data the status is ``UNKNOWN / NOT AVAILABLE`` and
``reason`` says why; nothing is filled in. The MCP server (measure_it.mcp.server), the deterministic agent
(measure_it.agents.deployment_agent) and the tests all call these functions.

Tool -> wrapped function
    search_condition / normalize_condition      measure_it.ontology.normalize
    get_patient_phenotype_signature             measure_it.wearables.signatures
    get_molecular_context                       measure_it.omics.query (+ omics.coherence summary)
    discover_candidate_measurements             measure_it.measurements.query
    get_measurement_evidence                    measure_it.measurements.query
    get_phenotype_measurement_evidence          measure_it.measurements.query (Phase 3)
    get_regulatory_context                      measure_it.measurements.regulatory
    get_condition_burden / get_geographic_context  measure_it.geography.query
    find_candidate_clinics                      measure_it.facilities.matching
    find_relevant_trials                        measure_it.ingestion.clinicaltrials.trials_for_condition
    find_relevant_research_centers              measure_it.facilities.query
    rank_deployment_opportunities               measure_it.scoring.recommend
    search_us_open_data                         measure_it.agents.datagov
    get_measurable_biology                      measure_it.omics.graph
    trace_evidence                              measure_it.mcp.trace
    list_sources                                SOURCE_REGISTRY.yaml

Caching: results are memoised per (tool, arguments, data stamp), where the data stamp is the newest mtime of
data/processed, configs, results/tables and SOURCE_REGISTRY.yaml (re-checked at most every 2 s), so a rebuilt table
invalidates every cached answer. The wrapped modules keep their own table caches.

    uv run python -m measure_it.tools --benchmark      # latency of every tool (cold and warm) on realistic inputs
"""
from __future__ import annotations

import json
import os
import re
import sys
import threading
import time
from collections import OrderedDict

from .config import CONFIGS, PROCESSED, SOURCE_REGISTRY_PATH, TABLES, UNKNOWN
from .mcp.trace import NAMESPACES, registry, sources_for_table, to_jsonable
from .mcp.trace import trace_evidence as _trace

OK, ERROR = "ok", "error"
SPEC_TOOLS = [
    "search_condition", "normalize_condition", "get_patient_phenotype_signature", "get_molecular_context",
    "discover_candidate_measurements", "get_measurement_evidence", "get_regulatory_context", "get_condition_burden",
    "get_geographic_context", "find_candidate_clinics", "find_relevant_trials", "find_relevant_research_centers",
    "rank_deployment_opportunities", "search_us_open_data", "trace_evidence",
]
EXTRA_TOOLS = ["get_phenotype_measurement_evidence", "get_measurable_biology", "list_sources"]
MAX_LIST = 100          # safety cap for any list inside a tool result
MAX_PROVENANCE_IDS = 150
MAX_BYTES = 60_000      # size guard: MCP clients truncate long tool output (Claude Code: ~25k tokens by default)
RANK_MAX_BYTES = 90_000  # rank_deployment_opportunities compact: 10 regions + shared context

MOLECULAR_LABEL = ("condition-level molecular enrichment (public databases; not participant-linked; not patient "
                   "multi-omics; not mechanism)")
GENERAL_CAVEAT = ("Tool output only: do not add facts, numbers or places that are not in it; report UNKNOWN / NOT "
                  "AVAILABLE as unknown.")

# tables each tool reads (for provenance -> SOURCE_REGISTRY sources)
TOOL_TABLES: dict[str, list[str]] = {
    "search_condition": ["condition_registry", "condition_ontology_mappings"],
    "normalize_condition": ["condition_registry", "condition_ontology_mappings", "condition_phenotypes",
                            "ontology_icd10cm_codes"],
    "get_patient_phenotype_signature": ["phenotype_signatures"],
    "get_molecular_context": ["condition_molecular_evidence", "condition_molecular_coverage", "test3_reactome_ora"],
    "discover_candidate_measurements": ["condition_measurement_evidence", "measurement_registry",
                                        "measurement_phenotype_signal", "measurable_biology"],
    "get_measurement_evidence": ["condition_measurement_evidence", "measurement_registry", "measurement_mentions",
                                 "measurement_phenotype_signal", "measurable_biology"],
    "get_phenotype_measurement_evidence": ["phase3_measurement_evidence", "measurement_phenotype_signal",
                                           "measurement_registry"],
    "get_regulatory_context": ["measurement_regulatory_status", "fda_510k", "fda_named_device_findings"],
    "get_condition_burden": ["geo_condition_burden", "geo_burden_definitions", "geo_condition_features"],
    "get_geographic_context": ["geographies", "geo_context", "geo_condition_features", "provider_density_county",
                               "facilities__hrsa", "trial_sites", "clinical_trials", "nih_projects"],
    "find_candidate_clinics": ["facilities", "facility_trials", "facility_nih_projects", "geographies"],
    "find_relevant_trials": ["clinical_trials", "trial_conditions", "trial_interventions", "trial_outcomes"],
    "find_relevant_research_centers": ["facility_trials", "facility_nih_projects", "facilities"],
    "rank_deployment_opportunities": ["deployment_opportunities", "geo_condition_features", "facilities",
                                      "facility_trials", "facility_nih_projects", "providers", "phenotype_signatures",
                                      "condition_molecular_evidence", "measurement_regulatory_status",
                                      "condition_measurement_evidence"],
    "search_us_open_data": ["data_gov_search_results"],
    "get_measurable_biology": ["measurable_biology"],
    "trace_evidence": [],
    "list_sources": [],
}
TOOL_PRODUCERS = {
    "search_condition": "measure_it.ontology.normalize.search_condition",
    "normalize_condition": "measure_it.ontology.normalize.normalize_condition",
    "get_patient_phenotype_signature": "measure_it.wearables.signatures.get_patient_phenotype_signature",
    "get_molecular_context": "measure_it.omics.query.get_molecular_context + measure_it.omics.coherence."
                             "get_molecular_coherence",
    "discover_candidate_measurements": "measure_it.measurements.query.discover_candidate_measurements",
    "get_measurement_evidence": "measure_it.measurements.query.get_measurement_evidence",
    "get_phenotype_measurement_evidence": "measure_it.measurements.query.get_phenotype_measurement_evidence",
    "get_regulatory_context": "measure_it.measurements.regulatory.get_regulatory_context",
    "get_condition_burden": "measure_it.geography.query.get_condition_burden",
    "get_geographic_context": "measure_it.geography.query.get_geographic_context",
    "find_candidate_clinics": "measure_it.facilities.matching.find_candidate_clinics",
    "find_relevant_trials": "measure_it.ingestion.clinicaltrials.trials_for_condition",
    "find_relevant_research_centers": "measure_it.facilities.query.find_relevant_research_centers",
    "rank_deployment_opportunities": "measure_it.scoring.recommend.rank_deployment_opportunities",
    "search_us_open_data": "measure_it.agents.datagov.search_us_open_data",
    "get_measurable_biology": "measure_it.omics.graph.get_measurable_biology",
    "trace_evidence": "measure_it.mcp.trace.trace_evidence",
    "list_sources": "SOURCE_REGISTRY.yaml (measure_it.registry)",
}


# ============================================================================================ envelope helpers

_OID_NS = sorted(NAMESPACES, key=len, reverse=True)
OBJECT_ID_FULL_RE = re.compile(r"^(?:" + "|".join(map(re.escape, _OID_NS)) + r"):\S+$")


def collect_object_ids(obj, out: list | None = None) -> list[str]:
    """Every string value that IS an object id (<namespace>:<id> with a CONVENTIONS namespace), in order, unique.
    A dict marked `"stored": False` (an on-the-fly opportunity row, provenance.opportunity_row) is skipped: its id has
    no stored row, so it is not a traceable object id."""
    out = [] if out is None else out
    seen = set(out)
    stack = [obj]
    while stack:
        v = stack.pop()
        if isinstance(v, dict):
            if v.get("stored") is False:
                continue
            stack.extend(reversed(list(v.values())))
        elif isinstance(v, (list, tuple)):
            stack.extend(reversed(list(v)))
        elif isinstance(v, str) and OBJECT_ID_FULL_RE.match(v) and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def cap_lists(obj, max_items: int, path: str, trunc: list, protect: tuple = ()) -> object:
    """Cap every list longer than max_items (recursively), recording each cut in `trunc`; lists at a path in
    `protect` keep their length (their items are still capped inside)."""
    if isinstance(obj, dict):
        return {k: cap_lists(v, max_items, f"{path}.{k}", trunc, protect) for k, v in obj.items()}
    if isinstance(obj, list):
        if len(obj) > max_items and path not in protect:
            trunc.append({"path": path, "returned": max_items, "total": len(obj)})
            obj = obj[:max_items]
        return [cap_lists(v, max_items, f"{path}[]", trunc, protect) for v in obj]
    return obj


def _cap(lst, n: int, path: str, trunc: list) -> list:
    lst = list(lst or [])
    if n is not None and len(lst) > n:
        trunc.append({"path": path, "returned": int(n), "total": len(lst)})
        return lst[:n]
    return lst


def _source_brief(source_id: str) -> dict:
    s = registry()["by_id"].get(source_id)
    if s is None:
        return {"source_id": source_id, "status": UNKNOWN}
    return {"source_id": source_id, "name": s.get("name"), "source_version": s.get("source_version"),
            "retrieved_at": s.get("retrieved_at"), "data_audit": f"data/raw/{source_id}/DATA_AUDIT.md"}


def _provenance(tool: str, data, extra_tables: list[str] | None = None) -> dict:
    tables = list(dict.fromkeys((TOOL_TABLES.get(tool) or []) + (extra_tables or [])))
    src = []
    for t in tables:
        for s in sources_for_table(t):
            if s not in src:
                src.append(s)
    ids = collect_object_ids(data)
    out = {"producer": TOOL_PRODUCERS.get(tool), "tables": tables, "sources": [_source_brief(s) for s in src],
           "object_ids": ids[:MAX_PROVENANCE_IDS], "n_object_ids": len(ids),
           "trace_with": "trace_evidence(object_id) resolves any of these ids to rows, sources and raw files"}
    return out


def _pop_caveats(d: dict) -> list[str]:
    out = []
    if not isinstance(d, dict):
        return out
    for k in ("caveats", "guardrails", "guardrail", "general_uncertainties", "regulatory_note", "interpretation",
              "framing"):
        v = d.get(k)
        if isinstance(v, str) and v and v != UNKNOWN:
            out.append(v)
        elif isinstance(v, list):
            out += [x for x in v if isinstance(x, str)]
    return out


def _dedupe(xs: list[str]) -> list[str]:
    return list(dict.fromkeys(x for x in xs if x))


# ============================================================================================ caching

_CACHE: "OrderedDict[str, str]" = OrderedDict()
_CACHE_MAX = 256
_CACHE_LOCK = threading.Lock()
_STAMP = {"value": None, "checked": 0.0}
STAMP_TTL_S = 2.0


def data_stamp() -> float:
    """Newest mtime of the inputs every tool reads (checked at most every STAMP_TTL_S seconds)."""
    now = time.monotonic()
    if _STAMP["value"] is not None and now - _STAMP["checked"] < STAMP_TTL_S:
        return _STAMP["value"]
    m = 0.0
    for d in (PROCESSED, CONFIGS, TABLES):
        try:
            with os.scandir(d) as it:
                for e in it:
                    try:
                        m = max(m, e.stat().st_mtime)
                    except OSError:
                        continue
        except OSError:
            continue
    try:
        m = max(m, SOURCE_REGISTRY_PATH.stat().st_mtime)
    except OSError:
        pass
    _STAMP.update(value=m, checked=now)
    return m


def clear_cache() -> None:
    with _CACHE_LOCK:
        _CACHE.clear()
    _STAMP.update(value=None, checked=0.0)


def _run(tool: str, query: dict, fn, *, cache: bool = True) -> dict:
    """Call fn() -> envelope (without meta); memoise the JSON by (tool, query, data stamp); never raise."""
    t0 = time.perf_counter()
    stamp = data_stamp()
    key = json.dumps([tool, query, stamp], sort_keys=True, default=str)
    if cache:
        with _CACHE_LOCK:
            hit = _CACHE.get(key)
            if hit is not None:
                _CACHE.move_to_end(key)
        if hit is not None:
            env = json.loads(hit)
            env["meta"] = {"elapsed_s": round(time.perf_counter() - t0, 4), "cache_hit": True, "data_stamp": stamp}
            return env
    try:
        env = fn()
    except Exception as exc:  # noqa: BLE001 - a tool never raises to the caller
        env = {"status": ERROR, "data": None, "reason": f"{type(exc).__name__}: {str(exc)[:500]}",
               "caveats": [GENERAL_CAVEAT], "provenance": {"producer": TOOL_PRODUCERS.get(tool)}, "truncation": []}
    env = {"tool": tool, "status": env.get("status", ERROR), "query": to_jsonable(query), **{
        k: v for k, v in env.items() if k not in ("tool", "status", "query")}}
    env.setdefault("caveats", [])
    env.setdefault("truncation", [])
    if env["status"] != OK and not env.get("reason"):
        env["reason"] = "no data returned by the wrapped function"
    env = to_jsonable(env)
    text = json.dumps(env, allow_nan=False)
    if cache and env["status"] != ERROR:
        with _CACHE_LOCK:
            _CACHE[key] = text
            _CACHE.move_to_end(key)
            while len(_CACHE) > _CACHE_MAX:
                _CACHE.popitem(last=False)
    env["meta"] = {"elapsed_s": round(time.perf_counter() - t0, 4), "cache_hit": False, "data_stamp": stamp}
    return env


def _merge_trunc(trunc: list) -> list:
    """One entry per path: the smallest cap applied and the largest original length (n_lists when several lists at
    the same path, e.g. data.rows[].reasons, were cut)."""
    out: dict[str, dict] = {}
    for t in trunc:
        p = t.get("path")
        if p not in out:
            out[p] = dict(t)
            continue
        o = out[p]
        o["returned"] = min(o.get("returned", 0), t.get("returned", 0))
        o["total"] = max(o.get("total", 0), t.get("total", 0))
        o["n_lists"] = o.get("n_lists", 1) + 1
    return list(out.values())


PROTECTED_LISTS = ("data.recommendations", "data.candidates", "data.centers", "data.trials", "data.results",
                   "data.rows", "data.matches", "data.technologies", "data.sources", "data.datasets", "data.evidence")


def fit_size(data, max_bytes: int, trunc: list, floor: int = 10):
    """Safety net against very long outputs: while the JSON is above max_bytes, cap lists at 50, 25 and then `floor`
    items (each cut is recorded in truncation). The primary result lists (PROTECTED_LISTS) are never cut here (their
    length is set by the tool's own top_n / max_rows argument), and no list is cut below `floor` items, so an answer
    can stay above max_bytes rather than lose its short lists (members, caveats)."""
    size = lambda d: len(json.dumps(d, allow_nan=False, separators=(",", ":")))  # noqa: E731
    for n in (50, 25, floor):
        if size(data) <= max_bytes:
            return data
        data = cap_lists(data, n, "data", trunc, PROTECTED_LISTS)
    return data


def _env(tool: str, status: str, data, *, caveats: list | None = None, reason: str | None = None,
         truncation: list | None = None, extra_tables: list[str] | None = None, max_list: int = MAX_LIST,
         max_bytes: int | None = None) -> dict:
    trunc = list(truncation or [])
    data = cap_lists(to_jsonable(data), max_list, "data", trunc)
    limit = MAX_BYTES if max_bytes is None else max_bytes
    if limit:
        data = fit_size(data, limit, trunc)
    env = {"status": status, "data": data, "caveats": _dedupe(list(caveats or []) + [GENERAL_CAVEAT]),
           "provenance": _provenance(tool, data, extra_tables), "truncation": _merge_trunc(trunc)}
    if reason:
        env["reason"] = reason
    return env


def _status(v) -> str:
    return OK if v not in (UNKNOWN, "unknown", "ambiguous", None, ERROR) else UNKNOWN


# ============================================================================================ condition / measurement args

def resolve_condition_arg(text: str, *, allow_sets: bool = False) -> dict:
    """Condition text -> canonical id through the ontology normaliser (the normalize_condition tool).

    Canonical ids (and, with allow_sets, condition-set ids) pass through. Only an unambiguous 'matched' result is
    used; ambiguous / partial / unknown input returns UNKNOWN with the candidates, never a guess."""
    from .config import load_config
    from .ontology.normalize import normalize_condition as _norm
    q = str(text or "").strip()
    if not q:
        return {"status": UNKNOWN, "query": q, "reason": "empty condition"}
    ids = {c["id"] for c in load_config("conditions")["conditions"]}
    if q in ids:
        return {"status": "matched", "query": q, "condition_id": q, "object_id": f"condition:{q}",
                "match_reason": "canonical condition id"}
    if allow_sets:
        from .scoring.opportunity import condition_sets
        if q in condition_sets():
            return {"status": "matched", "query": q, "condition_id": q, "object_id": f"condition:{q}",
                    "match_reason": "condition-set id", "members": condition_sets()[q]["members"]}
    r = _norm(q)
    if r.get("status") == "matched":
        m = r["matches"][0]
        return {"status": "matched", "query": q, "condition_id": m["canonical_condition_id"],
                "preferred_name": m.get("preferred_name"), "object_id": f"condition:{m['canonical_condition_id']}",
                "match_reason": m.get("match_reason"), "confidence": m.get("confidence"),
                "other_candidates": [c["canonical_condition_id"] for c in r.get("other_candidates", [])][:5]}
    cands = [c["canonical_condition_id"] for c in (r.get("matches") or []) + (r.get("other_candidates") or [])]
    reason = r.get("reason") or (f"normalize_condition status {r.get('status')!r}: not one unambiguous registry "
                                 f"condition (candidates {cands[:5]})")
    return {"status": UNKNOWN, "query": q, "reason": reason, "candidates": cands[:5],
            **({"phenotype_axis": r["phenotype_axis"]} if r.get("phenotype_axis") else {})}


def _cond_or_unknown(tool: str, condition: str, **kw):
    rc = resolve_condition_arg(condition, **kw)
    if rc["status"] != "matched":
        return None, _env(tool, UNKNOWN, {"condition_resolution": rc}, reason=f"condition: {rc['reason']}")
    return rc, None


# ============================================================================================ tools: ontology

def search_condition(condition: str, limit: int = 10) -> dict:
    """Ranked registry conditions for free text, each with ids (MONDO/EFO/MeSH/ICD-10-CM) and why it matched."""
    def go():
        from .ontology.normalize import search_condition as _search
        res = _search(condition, limit=max(1, min(int(limit), 25)))
        for r in res:
            r["object_id"] = f"condition:{r['canonical_condition_id']}"
        if not res:
            return _env("search_condition", UNKNOWN, {"query": condition, "candidates": []},
                        reason="no registry condition matches this text (exact, word-subset or close-spelling "
                               "match; acronyms match exactly only)")
        return _env("search_condition", OK, {"query": condition, "candidates": res,
                                             "scoring": "1.0 x kind weight exact; 0.8 word subset; 0.75 x similarity "
                                                        "close spelling"})
    return _run("search_condition", {"condition": condition, "limit": limit}, go)


def normalize_condition(condition_or_code: str) -> dict:
    """Free text, alias, ICD-10-CM code, or MONDO/EFO/MeSH/HPO id -> canonical registry condition(s)."""
    def go():
        from .ontology.normalize import normalize_condition as _norm
        r = _norm(condition_or_code)
        st = r.get("status")
        data = dict(r)
        data["match_status"] = data.pop("status", None)
        for m in data.get("matches", []) + data.get("other_candidates", []):
            m["object_id"] = f"condition:{m['canonical_condition_id']}"
        cav = []
        if st in ("ambiguous", "partial"):
            cav.append(f"match_status is {st!r}: more than one condition fits or the match is below the exact "
                       "threshold; do not treat the first match as the answer.")
        if data.get("input_type") == "hpo":
            cav.append("HPO ids return conditions annotated with the phenotype (a phenotype association, not an "
                       "equivalence).")
        if data.get("phenotype_axis"):
            cav.append("The text is also a phenotype axis (a phenotype, not a condition).")
        # only an unambiguous 'matched' result is an answer (as in resolve_condition_arg): a partial or ambiguous
        # match keeps its candidates in data but the envelope status is UNKNOWN, so no caller takes matches[0]
        if st in ("ambiguous", "partial"):
            cands = [m["canonical_condition_id"] for m in data.get("matches", []) + data.get("other_candidates", [])]
            return _env("normalize_condition", UNKNOWN, data, caveats=cav,
                        reason=(f"match_status {st!r}: not one unambiguous registry condition; candidates "
                                f"{cands[:5]} are listed in data.matches / data.other_candidates, none is chosen"))
        return _env("normalize_condition", UNKNOWN if st == UNKNOWN else OK, data, caveats=cav,
                    reason=r.get("reason") if st == UNKNOWN else None)
    return _run("normalize_condition", {"condition_or_code": condition_or_code}, go)


# ============================================================================================ tools: person-level

def get_patient_phenotype_signature(condition: str, top_n: int = 5) -> dict:
    """Wearable phenotype signatures (group differences in public person-level datasets) for a condition/phenotype."""
    def go():
        from .wearables.signatures import get_patient_phenotype_signature as _sig
        r = _sig(condition, top_n=max(1, min(int(top_n), 20)))
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        trunc = []
        for i, d in enumerate(data.get("datasets") or []):
            nr = d.get("null_results") or {}
            nr["features"] = _cap(nr.get("features"), 25, f"data.datasets[{i}].null_results.features", trunc)
            d["descriptive_rows_not_tested"] = _cap(d.get("descriptive_rows_not_tested"), 10,
                                                    f"data.datasets[{i}].descriptive_rows_not_tested", trunc)
            d["object_ids"] = _cap(d.get("object_ids"), 30, f"data.datasets[{i}].object_ids", trunc)
        cav = _pop_caveats(data) + [
            "Person-level datasets are summarised as group statistics only; these participants are not the people of "
            "any geography, facility or trial returned by other tools.",
            "is_proxy = true means the label is a proxy (e.g. an ME/CFS-like symptom proxy), not the named condition."]
        return _env("get_patient_phenotype_signature", OK if st == "found" else UNKNOWN, data, caveats=cav,
                    truncation=trunc, reason=r.get("reason") if st != "found" else None)
    return _run("get_patient_phenotype_signature", {"condition": condition, "top_n": top_n}, go)


# ============================================================================================ tools: molecular

_MOL_ID_RULES = (("gene", re.compile(r"ENSG\d{11}"), "gene"), ("variant", re.compile(r"rs\d+"), "variant"),
                 ("study", re.compile(r"GCST\d+"), "gwas_study"), ("study", re.compile(r"GSE\d+"), "geo_series"))


def molecular_object_id(entity_type, entity_id) -> str | None:
    """CONVENTIONS section 7 object id of a condition_molecular_evidence entity (Ensembl gene, rsID, GWAS Catalog study,
    GEO series), or None when the entity has no namespace (drugs, pathways, SRA studies, analytes, gene symbols)."""
    for et, rx, ns in _MOL_ID_RULES:
        if entity_type == et and isinstance(entity_id, str) and rx.fullmatch(entity_id):
            return f"{ns}:{entity_id}"
    return None


def _annotate_molecular_ids(evidence) -> None:
    """Add object_id to every top entity of get_molecular_context's evidence summary (in place) so the entities are
    cited and traceable."""
    if not isinstance(evidence, dict):
        return
    for block in evidence.values():
        te = block.get("top_entities") if isinstance(block, dict) else None
        groups = te.values() if isinstance(te, dict) else [te] if isinstance(te, list) else []
        for items in groups:
            for it in items or []:
                if isinstance(it, dict):
                    oid = molecular_object_id(it.get("entity_type"), it.get("entity_id"))
                    if oid:
                        it["object_id"] = oid


def get_molecular_context(condition: str, top_n: int = 5, include_broad: bool = True) -> dict:
    """Condition-level molecular enrichment (Open Targets, GWAS Catalog, GEO/SRA, mapMECFS) + Test 3 coherence."""
    def go():
        rc, unk = _cond_or_unknown("get_molecular_context", condition)
        if unk:
            return unk
        from .omics.coherence import get_molecular_coherence
        from .omics.query import get_molecular_context as _mol
        cid = rc["condition_id"]
        r = _mol(cid, top_n=max(1, min(int(top_n), 25)), include_broad=bool(include_broad))
        coh = get_molecular_coherence(cid)
        data = {"label": MOLECULAR_LABEL, "condition_resolution": rc,
                **{k: v for k, v in r.items() if k != "status"}, "molecular_status": r.get("status")}
        _annotate_molecular_ids(data.get("evidence"))
        cav = _pop_caveats(data)
        cstat = coh.get("status")
        data["coherence"] = {
            "status": cstat, "reason": coh.get("reason"), "source_counts": coh.get("source_counts"),
            "gene_agreement": coh.get("gene_agreement"), "supported_systems": coh.get("supported_systems"),
            "disagreements": coh.get("disagreements"),
            "post_hoc_flag_note": "supported_systems[].post_hoc_flagged = true: the system rests partly on drug-target "
                                  "gene families of trialled drugs or on one GWAS locus (omics.coherence.post_hoc_flagged)"
                                  "; read those rows with caution.",
            "method": "Test 3 (docs/ANALYSIS_PLAN_MOLECULAR.md; results/MOLECULAR_COHERENCE.md)"}
        cav += [c for c in (coh.get("caveats") or []) if isinstance(c, str)]
        cav.append("Always call this condition-level molecular enrichment; it is never patient or personal "
                   "multi-omics.")
        ok = r.get("status") == "available" or cstat == "available"
        return _env("get_molecular_context", OK if ok else UNKNOWN, data, caveats=cav,
                    reason=None if ok else (r.get("reason") or coh.get("reason") or "no molecular evidence rows"))
    return _run("get_molecular_context", {"condition": condition, "top_n": top_n, "include_broad": include_broad}, go)


def get_measurable_biology(condition: str, measurement: str | None = None, include_unsupported: bool = False) -> dict:
    """Measurement classes with condition-level molecular support (condition -> system data-derived, system -> class
    curated) for one condition, optionally one class or bundle."""
    def go():
        rc, unk = _cond_or_unknown("get_measurable_biology", condition)
        if unk:
            return unk
        from .omics.graph import get_measurable_biology as _mb
        r = _mb(rc["condition_id"], measurement, include_unsupported=bool(include_unsupported))
        st = r.get("status")
        data = {"label": MOLECULAR_LABEL, "condition_resolution": rc, **{k: v for k, v in r.items() if k != "status"},
                "result": st}
        cav = _pop_caveats(data)
        if st == "not_supported":
            cav.append("result = not_supported is a tested negative (molecular data exist, no mapped system is "
                       "supported), not missing data.")
        ok = st in ("available", "not_supported")
        return _env("get_measurable_biology", OK if ok else UNKNOWN, data, caveats=cav,
                    reason=None if ok else r.get("reason"))
    return _run("get_measurable_biology", {"condition": condition, "measurement": measurement,
                                           "include_unsupported": include_unsupported}, go)


# ============================================================================================ tools: measurements

def discover_candidate_measurements(condition: str, phenotype: str | None = None, top_n: int = 10) -> dict:
    """Candidate objective measurements for a condition (optionally one phenotype), in the documented research-activity
    order, with the five separate dimensions per candidate."""
    def go():
        rc, unk = _cond_or_unknown("discover_candidate_measurements", condition)
        if unk:
            return unk
        from .measurements.query import discover_candidate_measurements as _disc
        r = _disc(rc["condition_id"], phenotype or None, top_n=max(1, min(int(top_n), 30)))
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        data["condition_resolution"] = rc
        cav = _pop_caveats(data)
        return _env("discover_candidate_measurements", OK if st == "ok" else UNKNOWN, data, caveats=cav,
                    reason=r.get("reason") if st != "ok" else None)
    return _run("discover_candidate_measurements", {"condition": condition, "phenotype": phenotype, "top_n": top_n}, go)


def get_measurement_evidence(measurement: str, condition: str, max_object_ids: int = 25) -> dict:
    """Evidence (trial / NIH / FDA object ids, counts, snippets, phenotype-signal results) that research on a condition
    deploys a measurement class or bundle."""
    def go():
        rc, unk = _cond_or_unknown("get_measurement_evidence", condition)
        if unk:
            return unk
        from .measurements.query import get_measurement_evidence as _ev
        r = _ev(measurement, rc["condition_id"], max_object_ids=max(1, min(int(max_object_ids), 100)))
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        data["condition_resolution"] = rc
        trunc = []
        for i, e in enumerate(data.get("evidence") or []):
            for k in ("trial_object_ids", "grant_object_ids"):
                n = e.get(f"n_{k}")
                if n is not None and n > len(e.get(k) or []):
                    trunc.append({"path": f"data.evidence[{i}].{k}", "returned": len(e.get(k) or []), "total": n})
        cav = _pop_caveats(data)
        return _env("get_measurement_evidence", OK if st == "ok" else UNKNOWN, data, caveats=cav, truncation=trunc,
                    reason=r.get("reason") if st != "ok" else None)
    return _run("get_measurement_evidence", {"measurement": measurement, "condition": condition,
                                             "max_object_ids": max_object_ids}, go)


def get_phenotype_measurement_evidence(phenotype: str) -> dict:
    """Phase 3 table for a demo phenotype (e.g. 'orthostatic/autonomic dysfunction', 'orthostatic intolerance'):
    observable signals, candidate technologies, evidence per source and five separate dimensions."""
    def go():
        from .measurements.query import get_phenotype_measurement_evidence as _p3
        r = _p3(phenotype)
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        cav = _pop_caveats(data)
        return _env("get_phenotype_measurement_evidence", OK if st == "ok" else UNKNOWN, data, caveats=cav,
                    reason=r.get("reason") if st != "ok" else None)
    return _run("get_phenotype_measurement_evidence", {"phenotype": phenotype}, go)


def get_regulatory_context(technology: str, max_records: int = 10) -> dict:
    """FDA regulatory context (product codes, 510(k)/De Novo/PMA counts, example records) for a technology name."""
    def go():
        from .measurements.regulatory import get_regulatory_context as _reg
        r = _reg(technology)
        if not any(isinstance(m, dict) and m.get("method") != "class_vocabulary_fuzzy"
                   for m in r.get("matched_measurements") or []):
            # no exact / pattern class match: a bundle name resolves to its member classes rather than to fuzzy
            # word overlap ('wearable autonomic monitoring' otherwise matched autonomic_testing at score 0.4)
            r = _regulatory_for_bundle(technology, _reg) or r
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        trunc = []
        fuzzy_only = [m for m in data.get("matched_measurements") or [] if isinstance(m, dict)]
        fuzzy_only = bool(fuzzy_only) and all(m.get("method") == "class_vocabulary_fuzzy" for m in fuzzy_only)
        data["matching_fda_records"] = _cap(data.get("matching_fda_records"), max(0, int(max_records)),
                                            "data.matching_fda_records", trunc)
        for rec in data["matching_fda_records"]:
            if rec.get("clearance_id"):
                rec["object_id"] = f"fda:{rec['clearance_id']}"
        status_rows = data.get("regulatory_status")
        for row in status_rows if isinstance(status_rows, list) else []:
            if row.get("product_code") not in (None, UNKNOWN):
                row["product_code_object_id"] = f"fda:code:{row['product_code']}"
        # the matched measurement classes are cited too (a class without a product code still has a registry row)
        for key in ("matched_measurements", "summary", "regulatory_status"):
            rows = data.get(key)
            for row in rows if isinstance(rows, list) else []:
                if isinstance(row, dict) and isinstance(row.get("measurement_id"), str) and row["measurement_id"]:
                    row.setdefault("measurement_object_id", f"measurement:{row['measurement_id']}")
        cav = _pop_caveats(data)
        if fuzzy_only:
            cav.append("Every matched class was found only by fuzzy vocabulary overlap (method class_vocabulary_fuzzy, "
                       "score < 0.95): these are candidate classes for the text, not an identification of the "
                       "technology; name a measurement class or bundle for an exact match.")
        return _env("get_regulatory_context", OK if st == "matched" else UNKNOWN, data, caveats=cav, truncation=trunc,
                    reason=r.get("reason") if st != "matched" else None)
    return _run("get_regulatory_context", {"technology": technology, "max_records": max_records}, go)


def _regulatory_for_bundle(technology: str, reg) -> dict | None:
    """A measurement BUNDLE name ('wearable autonomic monitoring', 'nailfold capillaroscopy') -> the regulatory rows of
    its member classes (exact class matches), instead of fuzzy word overlap with unrelated classes. None when the text
    is not exactly a bundle id / label / alias (measurements.query.resolve_measurement)."""
    if not str(technology or "").strip():
        return None
    from .measurements.query import resolve_measurement
    rm = resolve_measurement(technology)
    if rm.get("status") != "matched" or rm.get("kind") != "bundle" or not rm.get("member_ids"):
        return None
    parts = {m: reg(m) for m in rm["member_ids"]}
    hit = {m: p for m, p in parts.items() if p.get("status") == "matched"}
    seen, records = set(), []
    for p in hit.values():
        for rec in p.get("matching_fda_records") or []:
            k = rec.get("clearance_id") or json.dumps(rec, sort_keys=True, default=str)
            if k not in seen:
                seen.add(k)
                records.append(rec)
    first = next(iter(parts.values()), {})
    out = {"query": technology, "regulatory_note": first.get("regulatory_note"),
           "status": "matched" if hit else UNKNOWN,
           "measurement_resolution": {k: rm.get(k) for k in ("measurement_id", "kind", "member_ids", "match_reason")},
           "matched_measurements": [x for p in hit.values() for x in p.get("matched_measurements") or []],
           "summary": [x for p in hit.values() for x in (p.get("summary") or [])],
           "regulatory_status": [x for p in hit.values() for x in (p.get("regulatory_status")
                                                                  if isinstance(p.get("regulatory_status"), list) else [])],
           "matching_fda_records": records,
           "members_without_regulatory_rows": {m: p.get("reason") for m, p in parts.items() if m not in hit},
           "provenance": next((p.get("provenance") for p in hit.values() if p.get("provenance")), None)}
    if not hit:
        out["reason"] = f"no member class of bundle {rm.get('measurement_id')} has FDA regulatory rows"
    return out


# ============================================================================================ tools: geography

def _burden_summary(rows: list[dict]) -> dict:
    vals = [r["value"] for r in rows if isinstance(r.get("value"), (int, float)) and not isinstance(r.get("value"), bool)]
    vals_sorted = sorted(vals)
    n = len(vals_sorted)
    med = (vals_sorted[n // 2] if n % 2 else (vals_sorted[n // 2 - 1] + vals_sorted[n // 2]) / 2) if n else None
    return {"n_rows": len(rows), "n_with_value": n, "n_value_unknown": len(rows) - n,
            "min": vals_sorted[0] if n else UNKNOWN, "median": med if n else UNKNOWN,
            "max": vals_sorted[-1] if n else UNKNOWN,
            "value_units": sorted({str(r.get("value_unit")) for r in rows if r.get("value_unit")}),
            "evidence_levels": sorted({str(r.get("burden_evidence_level")) for r in rows}),
            "source_resolutions": sorted({str(r.get("source_geographic_resolution")) for r in rows}),
            "n_inherited_state_values": sum(1 for r in rows if r.get("inherited") is True),
            "n_modelled_small_area_values": sum(1 for r in rows if r.get("derivation") == "mrp_small_area_estimate")}


def get_condition_burden(condition: str, geography_level: str = "state", geo_id: str | None = None,
                         include_alternates: bool = False, max_rows: int = 60, order: str = "geo_id") -> dict:
    """Burden estimates for a condition at 'national', 'state' or 'county' level, each row with its
    burden_evidence_level (A direct / B coded condition / C proxy / D none), source resolution and inherited flag."""
    def go():
        rc, unk = _cond_or_unknown("get_condition_burden", condition)
        if unk:
            return unk
        from .geography.query import get_condition_burden as _b
        r = _b(rc["condition_id"], geography_level, geo_id=geo_id, include_alternates=bool(include_alternates))
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        data["condition_resolution"] = rc
        rows = data.get("rows") or []
        trunc = []
        if rows:
            data["summary"] = _burden_summary(rows)
            if order == "value_desc":
                rows = sorted(rows, key=lambda x: (-(x["value"]) if isinstance(x.get("value"), (int, float)) else
                                                   float("inf"), str(x.get("geo_id"))))
            data["rows"] = _cap(rows, max(1, int(max_rows)), "data.rows", trunc)
            data["row_order"] = "value descending (presentation order only)" if order == "value_desc" else "geo_id"
        cav = _pop_caveats(data)
        lvl = data.get("burden_evidence_level")
        if lvl:
            cav.append(f"burden_evidence_level {lvl}: A = direct condition measure, B = closely matching coded "
                       "condition, C = symptom/comorbidity proxy, D = no usable burden estimate. Always state it.")
        if data.get("summary", {}).get("n_inherited_state_values"):
            cav.append("Rows with inherited = true are STATE estimates attached to counties (source resolution "
                       "'state'); never describe them as county prevalence.")
        if data.get("summary", {}).get("n_modelled_small_area_values") and r.get("geo_level", geography_level) == "county":
            cav.append("Rows with derivation = mrp_small_area_estimate are MODELLED county estimates (BRFSS 2023 "
                       "multilevel regression and post-stratification); their within-state differences come from "
                       "county composition and covariates only. Never describe them as observed county prevalence.")
        return _env("get_condition_burden", _status(st), data, caveats=cav, truncation=trunc,
                    reason=r.get("reason") if _status(st) != OK else None)
    return _run("get_condition_burden", {"condition": condition, "geography_level": geography_level, "geo_id": geo_id,
                                         "include_alternates": include_alternates, "max_rows": max_rows,
                                         "order": order}, go)


def get_geographic_context(geography_id: str) -> dict:
    """One state, county or ZCTA (FIPS, 'geo:<fips>', 'zcta:<zcta>' or a name such as 'San Diego County,
    California'): burden across conditions, vulnerability, ACS, providers, HRSA sites, trials, NIH projects."""
    def go():
        from .geography.query import get_geographic_context as _g
        r = _g(geography_id)
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        cav = _pop_caveats(data)
        return _env("get_geographic_context", _status(st), data, caveats=cav,
                    reason=r.get("reason") if _status(st) != OK else None)
    return _run("get_geographic_context", {"geography_id": geography_id}, go)


# ============================================================================================ tools: facilities / research

def find_candidate_clinics(geography_id: str, condition: str, measurement: str | None = None,
                           radius_km: float | None = None, max_sites: int | None = None) -> dict:
    """Candidate implementation / study-partner facilities around a geography, six characteristics ranked separately
    (never one score, never a quality ranking of facilities)."""
    def go():
        rc, unk = _cond_or_unknown("find_candidate_clinics", condition)
        if unk:
            return unk
        from .facilities.matching import find_candidate_clinics as _fcc
        r = _fcc(geography_id, rc["condition_id"], measurement or None, radius_km=radius_km,
                 max_sites=None if max_sites is None else max(1, min(int(max_sites), 50)))
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        data["result"] = st
        data["condition_resolution"] = rc
        cav = _pop_caveats(data)
        if data.get("absence_note"):
            cav.append(data["absence_note"])
        ok = st == "ok"
        reason = None if ok else (r.get("reason") or "no candidate facility")
        return _env("find_candidate_clinics", OK if ok else UNKNOWN, data, caveats=cav, reason=reason)
    return _run("find_candidate_clinics", {"geography_id": geography_id, "condition": condition,
                                           "measurement": measurement, "radius_km": radius_km,
                                           "max_sites": max_sites}, go)


TRIAL_COLS = ["nct_id", "brief_title", "overall_status", "study_type", "phases", "start_date",
              "primary_completion_date", "last_update_post_date", "enrollment_count", "lead_sponsor",
              "lead_sponsor_class", "conditions", "n_us_locations", "has_results", "study_url", "condition_terms",
              "condition_literal_match", "measurement_match_fields", "measurement_match_example", "source_version",
              "retrieved_at"]


def find_relevant_trials(condition: str, measurement: str | None = None, literal_only: bool = True,
                         statuses: list[str] | None = None, us_sites_only: bool = False, max_trials: int = 25) -> dict:
    """ClinicalTrials.gov trials registered for a condition (literal condition match by default), optionally only those
    whose registered text mentions a measurement class / bundle / alias."""
    def go():
        import pandas as pd
        rc, unk = _cond_or_unknown("find_relevant_trials", condition, allow_sets=True)
        if unk:
            return unk
        from .ingestion.clinicaltrials import GUARDRAIL_NOTE, trials_for_condition
        members = rc.get("members") or [rc["condition_id"]]
        rm = None
        classes: list[str | None] = [None]
        if measurement is not None and str(measurement).strip():
            from .measurements.query import resolve_measurement
            rm = resolve_measurement(measurement)
            if rm.get("status") != "matched":
                return _env("find_relevant_trials", UNKNOWN, {"condition_resolution": rc, "measurement_resolution": rm},
                            reason=f"measurement: {rm.get('reason')}")
            classes = list(dict.fromkeys(rm["member_ids"]))
        frames = []
        for cid in members:
            for cls in classes:
                df = trials_for_condition(cid, cls, literal_only=bool(literal_only), statuses=statuses or None,
                                          us_sites_only=bool(us_sites_only))
                if len(df):
                    frames.append(df.assign(condition_id=cid, measurement_class=cls))
        base = {"condition_resolution": rc, "measurement_resolution": rm,
                "view": "literal condition match (trial_conditions.condition_literal_match)" if literal_only
                        else "all matched terms incl. MeSH-expanded matches",
                "fields_searched_for_measurement": "titles, brief summary, keywords, intervention names/descriptions/"
                                                   "other names, outcome measures/descriptions (not eligibility)"}
        if not frames:
            return _env("find_relevant_trials", UNKNOWN, {**base, "n_trials": 0, "trials": []},
                        caveats=[GUARDRAIL_NOTE],
                        reason=(f"no registered trial for {members} "
                                + (f"whose registered text mentions {classes} " if rm else "")
                                + ("(literal condition view)" if literal_only else "")))
        allt = pd.concat(frames, ignore_index=True)
        agg = {"condition_id": lambda s: sorted(set(s)),
               "measurement_class": lambda s: sorted({x for x in s if isinstance(x, str)})}
        first = {c: "first" for c in TRIAL_COLS if c in allt.columns and c != "nct_id"}
        g = allt.groupby("nct_id", sort=False).agg({**first, **agg}).reset_index()
        g = g.sort_values(["last_update_post_date", "nct_id"], ascending=[False, True])
        g.insert(0, "object_id", "trial:" + g["nct_id"])
        g = g.rename(columns={"condition_id": "condition_ids", "measurement_class": "measurement_classes_matched"})
        trunc = []
        recs = to_jsonable(g.to_dict("records"), max_text=300)
        summary = {"n_trials": int(len(g)), "by_overall_status": to_jsonable(g["overall_status"].value_counts()
                                                                             .to_dict()),
                   "n_recruiting": int((g["overall_status"] == "RECRUITING").sum()),
                   "by_study_type": to_jsonable(g["study_type"].value_counts().to_dict())}
        data = {**base, "summary": summary, "order": "last_update_post_date descending",
                "trials": _cap(recs, max(1, int(max_trials)), "data.trials", trunc)}
        return _env("find_relevant_trials", OK, data, caveats=[GUARDRAIL_NOTE], truncation=trunc)
    return _run("find_relevant_trials", {"condition": condition, "measurement": measurement,
                                         "literal_only": literal_only, "statuses": statuses,
                                         "us_sites_only": us_sites_only, "max_trials": max_trials}, go)


def find_relevant_research_centers(condition: str, measurement: str | None = None, top_n: int = 15,
                                   include_expanded: bool = False) -> dict:
    """Facilities with ClinicalTrials.gov site history or NIH RePORTER awards for a condition (research readiness,
    kept separate from burden), counts and object ids per resolved facility."""
    def go():
        rc, unk = _cond_or_unknown("find_relevant_research_centers", condition)
        if unk:
            return unk
        from .facilities.query import find_relevant_research_centers as _frc
        r = _frc(rc["condition_id"], measurement or None, top_n=max(1, min(int(top_n), 100)),
                 include_expanded=bool(include_expanded))
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        data["condition_resolution"] = rc
        trunc = []
        for i, c in enumerate(data.get("centers") or []):
            for k in ("trial_object_ids", "nih_object_ids", "npi_object_ids", "nih_org_names",
                      "measurement_classes_registered"):
                if isinstance(c.get(k), list):
                    c[k] = _cap(c[k], 10, f"data.centers[{i}].{k}", trunc)
            c.pop("member_names", None)
        n_tot, n_ret = data.get("n_facilities"), data.get("n_returned")
        if isinstance(n_tot, int) and isinstance(n_ret, int) and n_tot > n_ret:
            trunc.insert(0, {"path": "data.centers", "returned": n_ret, "total": n_tot})
        cav = _pop_caveats(data) + ["The order is a display order (trials + NIH core projects), not a score or a "
                                    "quality ranking."]
        return _env("find_relevant_research_centers", OK if st == "ok" else UNKNOWN, data, caveats=cav,
                    truncation=trunc, reason=r.get("reason") if st != "ok" else None)
    return _run("find_relevant_research_centers", {"condition": condition, "measurement": measurement,
                                                   "top_n": top_n, "include_expanded": include_expanded}, go)


# ============================================================================================ tools: deployment ranking

SITE_KEEP = ("object_id", "facility_name", "primary_kind", "city", "state", "lat", "lon", "geocode_precision",
             "distance_km", "inside_geography", "selected_via")
GEO_KEEP = ("name", "fips", "level", "object_id", "state", "population", "burden_evidence_level",
            "burden_evidence_level_note", "burden_excluded_level_D", "composite", "rank", "n_regions_ranked",
            "rank_interval_5_95", "p_top10_monte_carlo")
PERF_KEEP = ("object_id", "performance_status", "quality_tier", "selected_record_id", "label_basis", "comparator",
             "auroc", "auroc_ci_low", "auroc_ci_high", "op_sensitivity", "op_specificity", "measurement_evidence",
             "evidence_conflict", "performance_note")
YIELD_KEEP = ("basis", "performance_status", "reach", "reached_adults", "expected_yield_component",
              "expected_detectable_cases", "expected_false_positives", "expected_ppv", "expected_yield_index",
              "expected_false_positives_upper", "language")
BURDEN_MEMBER_KEEP = ("burden_value", "burden_evidence_level", "burden_inherited", "burden_source_resolution",
                      "burden_measure_id", "contributes_to_set_burden")


def _sub(d, keys) -> dict:
    return {k: d.get(k) for k in keys if isinstance(d, dict) and k in d}


def _compact_ids(block, n: int = 5) -> dict:
    """{'n': .., 'object_ids': [...]} blocks of research_evidence, ids capped (the count stays)."""
    if not isinstance(block, dict):
        return block
    out = {k: v for k, v in block.items() if k not in ("object_ids", "examples")}
    ids = block.get("object_ids") or []
    out["object_ids"] = ids[:n]
    if len(ids) > n:
        out["object_ids_note"] = f"{n} of {len(ids)} shown"
    return out


def _compact_phenotype(ph) -> object:
    if not isinstance(ph, list):
        return ph
    out = []
    for p in ph:
        if not isinstance(p, dict):
            out.append(p)
            continue
        q = _sub(p, ("condition_id", "status", "reason", "note"))
        ds = []
        for d in p.get("datasets") or []:
            e = _sub(d, ("dataset_id", "phenotype_id", "is_proxy", "n_cases", "n_controls", "n_features_tested",
                         "n_features_fdr_significant"))
            e["top_features"] = [_sub(f, ("feature", "effect_size", "ci_low", "ci_high", "q_value_bh", "object_id"))
                                 for f in (d.get("top_features") or [])[:2]]
            ds.append(e)
        if ds:
            q["datasets"] = ds
        out.append(q)
    return out


def _compact_molecular(mol) -> object:
    if not isinstance(mol, list):
        return mol
    out = []
    for m in mol:
        if not isinstance(m, dict):
            out.append(m)
            continue
        q = _sub(m, ("condition_id", "label", "status", "reason"))
        ev = {}
        for et, v in (m.get("evidence_summary") or {}).items():
            top = (v.get("top_entities") or [])[:1]
            ev[et] = {**_sub(v, ("n_records", "n_distinct_entities")),
                      "sources": sorted((v.get("counts_by_source") or {}).keys()),
                      "top_entity": (f"{top[0].get('entity_id')} ({top[0].get('entity_label')})" if top else None)}
            oid = molecular_object_id(top[0].get("entity_type"), top[0].get("entity_id")) if top else None
            if oid:
                ev[et]["top_entity_object_id"] = oid
        if ev:
            q["evidence_summary"] = ev
        agr = m.get("cross_source_genetic_agreement")
        if isinstance(agr, dict):
            q["cross_source_genetic_agreement"] = _sub(agr, ("open_targets_genetic_genes", "gwas_catalog_mapped_genes",
                                                             "in_both"))
        out.append(q)
    return out


def _compact_technology(tech: dict) -> dict:
    reg = []
    for r in tech.get("regulatory_context") or []:
        if isinstance(r, dict):
            reg.append({**_sub(r, ("measurement_id", "status", "reason")),
                        **({"summary": _sub(r.get("summary"), ("n_mapped_product_codes", "confidence_of_best_mapping",
                                                              "regulatory_visibility", "n_510k_total", "n_denovo_total",
                                                              "n_pma_total"))} if isinstance(r.get("summary"), dict)
                           else {})})
    tev = []
    for e in tech.get("technology_evidence") or []:
        if isinstance(e, dict):
            mes = e.get("measurement_evidence_strength") if isinstance(e.get("measurement_evidence_strength"),
                                                                        dict) else {}
            tev.append({"evidence_object_id": e.get("evidence_object_id") or
                        f"{e.get('condition_id')}: {e.get('status')} ({e.get('reason')})",
                        "tier": mes.get("tier"), "n_trials_objective": mes.get("n_trials_objective"),
                        "n_nih_core_projects": mes.get("n_nih_core_projects"),
                        "technology_maturity": e.get("technology_maturity")})
    out = {"regulatory_context": reg, "technology_evidence": tev, "guardrail": tech.get("guardrail")}
    if isinstance(tech.get("measurement_performance"), dict):
        out["measurement_performance"] = _sub(tech["measurement_performance"], PERF_KEEP)
    if tech.get("member_condition_performance"):
        out["member_condition_performance"] = [_sub(x, PERF_KEEP) for x in tech["member_condition_performance"]
                                               if isinstance(x, dict)]
    return out


def _compact_burden(b) -> object:
    if not isinstance(b, dict):
        return b
    out = {k: v for k, v in b.items() if k in ("value", "unit", "measure_id", "ci_95", "percentile",
                                                "member_mean_percentile", "inherited", "source_resolution",
                                                "members_without_usable_burden")}
    if isinstance(b.get("members"), dict):
        out["members"] = {c: _sub(m, BURDEN_MEMBER_KEEP) for c, m in b["members"].items()}
    return out


def _compact_recommendation(rec: dict, max_sites: int, trunc: list, path: str) -> dict:
    out = {"condition": {"see": "data.shared_context.condition"},
           "phenotype": {"see": "data.shared_context.phenotype"},
           "measurement": {"see": "data.shared_context.measurement"},
           "technology": {"see": "data.shared_context.technology"}}
    geo = rec.get("geography") or {}
    g = _sub(geo, GEO_KEEP)
    # the Monte Carlo interval is always the equal-weight one: when `rank` is not the equal-weight rank among
    # population-eligible regions, the compact block says which is which (as recommended_next_step does)
    rb = geo.get("rank_basis") or {}
    if rb and (rb.get("weight_set") != "equal" or str(rb.get("rank_column", "")).endswith("_incl_small")):
        g["rank_basis"] = rb
        g["monte_carlo_basis"] = geo.get("monte_carlo_basis")
        g["rank_equal_weights"] = (geo.get("ranks_under_weight_sets") or {}).get("equal")
    g["burden"] = _compact_burden(geo.get("burden"))
    for k, keep in (("vulnerability", ("svi_overall", "percentile")), ("diagnostic_desert", ("index", "percentile")),
                    ("clinic_capacity", ("index", "percentile")),
                    ("research_readiness", ("index", "percentile", "condition_trials_in_pool",
                                            "condition_nih_core_projects_in_pool",
                                            "technology_experience_trials_in_pool")),
                    ("technology_saturation", ("percentile",))):
        if isinstance(geo.get(k), dict):
            g[k] = _sub(geo[k], keep)
    if isinstance(geo.get("measurement_evidence"), dict):
        g["measurement_evidence"] = _sub(geo["measurement_evidence"], ("value", "performance_status", "performance_tier",
                                                                       "object_id"))
    if isinstance(geo.get("expected_yield"), dict):
        g["expected_yield"] = _sub(geo["expected_yield"], YIELD_KEEP)
    if isinstance(geo.get("evidence_weighted"), dict):
        g["evidence_weighted"] = _sub(geo["evidence_weighted"], ("status", "rank", "rank_interval_5_95"))
    out["geography"] = g
    sites = rec.get("candidate_sites")
    if isinstance(sites, list):
        sites = _cap(sites, max(0, int(max_sites)), f"{path}.candidate_sites", trunc)
        out["candidate_sites"] = [_sub(s, SITE_KEEP) for s in sites]
    else:
        out["candidate_sites"] = sites
    out["candidate_site_query_status"] = [_sub(q, ("condition_id", "status", "n_eligible"))
                                          for q in rec.get("candidate_site_query_status") or [] if isinstance(q, dict)]
    rev = rec.get("research_evidence")
    if isinstance(rev, dict):
        out["research_evidence"] = {k: (_compact_ids(v, 2) if isinstance(v, dict) else v) for k, v in rev.items()
                                    if k not in ("note", "pool")}
    else:
        out["research_evidence"] = rev
    out["molecular_context"] = {"see": "data.shared_context.molecular_context"}
    out["uncertainties"] = rec.get("uncertainties")
    prov = rec.get("provenance") or {}
    ids = prov.get("object_ids") or []
    out["provenance"] = {"object_ids": _cap(ids, 8, f"{path}.provenance.object_ids", trunc), "n_object_ids": len(ids),
                         "opportunity_row": prov.get("opportunity_row")}
    out["recommended_next_step"] = rec.get("recommended_next_step")
    return out


def rank_deployment_opportunities(condition: str, measurement: str, geography_level: str = "county", top_n: int = 10,
                                  weight_set: str = "equal", include_small_population: bool = False,
                                  max_sites_per_region: int = 3, detail: str = "compact") -> dict:
    """Top-N candidate deployment opportunities (SPEC output schema per region): components, rank interval, candidate
    sites, research evidence, uncertainties, provenance and a recommended next step. detail='compact' shows the
    condition-level context once (data.shared_context) and a few candidate sites per region; detail='full' returns
    everything (large)."""
    def go():
        from .scoring.recommend import rank_deployment_opportunities as _rank
        n = max(1, min(int(top_n), 50))
        r = _rank(condition, measurement, geography_level, top_n=n, weight_set=weight_set,
                  include_small_population=bool(include_small_population), with_context=True)
        st = r.get("status")
        data = {k: v for k, v in r.items() if k != "status"}
        cav = _pop_caveats(data)
        if st != "ok":
            return _env("rank_deployment_opportunities", UNKNOWN, data, caveats=cav, reason=r.get("reason"))
        trunc = []
        recs = data.get("recommendations") or []
        full = detail == "full"
        if not full and recs:
            first = recs[0]
            data["shared_context"] = {
                "note": "condition- and measurement-level context; identical for every region, so shown once. "
                        "Person-level signatures are group statistics from public cohorts, not the people of any "
                        "region; molecular context is condition-level molecular enrichment.",
                "condition": _sub(first.get("condition") or {}, ("condition_id", "label", "kind", "members",
                                                                 "object_ids")),
                "measurement": _sub(first.get("measurement") or {}, ("measurement_id", "label", "kind", "members",
                                                                     "implementer_groups", "adapter", "adapter_status",
                                                                     "object_id")),
                "phenotype": _compact_phenotype(first.get("phenotype")),
                "molecular_context": _compact_molecular(first.get("molecular_context")),
                "technology": {"name": (first.get("technology") or {}).get("name"),
                               **_compact_technology(first.get("technology") or {})},
                "research_evidence_pool": (first.get("research_evidence") or {}).get("pool")
                if isinstance(first.get("research_evidence"), dict) else None}
            data["recommendations"] = [_compact_recommendation(rec, max_sites_per_region, trunc,
                                                               f"data.recommendations[{i}]")
                                       for i, rec in enumerate(recs)]
            if isinstance(data.get("excluded_incomplete_burden"), list):
                ex = _cap(data["excluded_incomplete_burden"], 5, "data.excluded_incomplete_burden", trunc)
                data["excluded_incomplete_burden"] = [_sub(x, ("name", "fips", "object_id", "members_without_value",
                                                               "components_kept", "burden", "rank"))
                                                      if isinstance(x, dict) else x for x in ex]
            data["detail"] = ("compact: shared context shown once, candidate sites capped per region, six-characteristic "
                              "ranks in find_candidate_clinics; detail='full' returns everything")
            data["candidate_sites_note"] = ("candidate_sites are facilities with characteristics suggesting they may be "
                                            "viable implementation or study partners (not a quality ranking); "
                                            "selected_via names the characteristic that selected each; "
                                            "find_candidate_clinics(fips, condition, measurement) gives the reasons "
                                            "and the six separate characteristic ranks")
            data["provenance_note"] = ("each recommendation's provenance.object_ids lists the first ids; "
                                       "trace_evidence(<opportunity id>) walks back to the component rows and sources")
        else:
            data["detail"] = "full"
        cav.append("Each region is a candidate deployment opportunity for pilot evaluation, not a validated "
                   "diagnostic pathway. State burden_evidence_level with every burden value; an inherited state "
                   "value is not county prevalence.")
        return _env("rank_deployment_opportunities", OK, data, caveats=cav, truncation=trunc,
                    max_list=MAX_LIST if not full else 10 ** 6, max_bytes=0 if full else RANK_MAX_BYTES)
    return _run("rank_deployment_opportunities", {"condition": condition, "measurement": measurement,
                                                  "geography_level": geography_level, "top_n": top_n,
                                                  "weight_set": weight_set,
                                                  "include_small_population": include_small_population,
                                                  "max_sites_per_region": max_sites_per_region, "detail": detail}, go)


# ============================================================================================ tools: open data, trace, sources

DATAGOV_KEEP = ("rank", "title", "agency", "publisher", "parent_organization", "description", "access_level",
                "resource_urls", "resource_formats", "modified", "identifier", "landing_page", "retrieve_from",
                "api_endpoint", "retrieved_at", "more_available", "note")


def _norm_query(q) -> str:
    return " ".join(str(q or "").casefold().split())


def _cached_datagov_queries() -> list[str]:
    """The catalog queries run at ingestion (data_gov_discovery_log): the ones an offline build can answer."""
    try:
        from .store import read_table
        return list(dict.fromkeys(read_table("data_gov_discovery_log", columns=["query"])["query"].astype(str)))
    except Exception:  # noqa: BLE001 - no processed log: nothing to offer
        return []


def search_us_open_data(query: str, rows: int = 20) -> dict:
    """Data.gov catalog search (metadata only): title, agency, description, access level, resource URLs, format,
    modified date, and the agency URL that holds the records."""
    def go():
        from .agents.datagov import MAX_ROWS, search_us_open_data as _s
        from .http import offline
        q = str(query or "").strip()
        if not q:
            return _env("search_us_open_data", UNKNOWN, {"results": []}, reason="empty query")
        n = max(1, min(int(rows), MAX_ROWS))
        res = _s(q, rows=n)
        cav = ["Data.gov is a metadata catalog: it does not host the records; retrieve_from names the agency URL that "
               "does. A catalog hit is not evidence about a condition.",
               "Results are in Data.gov's own relevance order; this engine does not evaluate catalog relevance (the "
               "top hit can be unrelated to the query, e.g. 'Long-term Care and COVID-19' for 'long COVID')."]
        unknown = len(res) == 1 and res[0].get("is_unknown")
        cached = []
        if unknown and offline() and res[0].get("result_status") == "unreachable":
            # offline, the HTTP cache is keyed on the exact query string: answer a case/whitespace variant of a
            # query that was run at ingestion from that cached response, and say so
            cached = _cached_datagov_queries()
            alt = next((c for c in cached if _norm_query(c) == _norm_query(q) and c != q), None)
            if alt is not None:
                res2 = _s(alt, rows=n)
                if not (len(res2) == 1 and res2[0].get("is_unknown")):
                    res, unknown = res2, False
                    cav.append(f"Offline mode: answered from the cached response to the query {alt!r} (same words, "
                               "different case or spacing).")
        trunc = []
        out = []
        for r in res:
            rec = {k: r.get(k) for k in DATAGOV_KEEP if k in r}
            rec["resource_urls"] = _cap(rec.get("resource_urls"), 5, f"data.results[{len(out)}].resource_urls", trunc)
            if isinstance(rec.get("description"), str) and len(rec["description"]) > 400:
                rec["description"] = rec["description"][:400] + " ...[cut]"
            out.append(rec)
        if unknown:
            reason = f"{res[0].get('result_status')}: {res[0].get('note')}"
            if offline() and res[0].get("result_status") == "unreachable":
                cached = cached or _cached_datagov_queries()
                reason += (" Offline mode answers only queries already in the HTTP cache (matched ignoring case and "
                           "spacing); cached queries: " + ("; ".join(repr(c) for c in cached[:60]) or UNKNOWN) + ".")
            return _env("search_us_open_data", UNKNOWN, {"results": []}, caveats=cav, reason=reason)
        return _env("search_us_open_data", OK, {"n_results": len(out), "results": out}, caveats=cav, truncation=trunc)
    # not memoised here: measure_it.http caches every catalog response on disk, and an 'unreachable' answer must not
    # outlive a network outage
    return _run("search_us_open_data", {"query": query, "rows": rows}, go, cache=False)


def trace_evidence(object_id: str, max_rows: int = 20) -> dict:
    """Resolve an object id (<namespace>:<id>, docs/CONVENTIONS.md section 7) to its rows, provenance, SOURCE_REGISTRY
    entries, DATA_AUDIT.md paths, raw files and one level of lineage. Person-level ids are refused."""
    def go():
        r = _trace(object_id, max_rows=max(1, min(int(max_rows), 100)))
        st = r.get("status")
        data = {k: v for k, v in r.items() if k not in ("status", "truncation", "caveats")}
        trunc = [{**t, "path": t["path"]} for t in (r.get("truncation") or [])]
        cav = list(r.get("caveats") or [])
        if r.get("refused"):
            cav.append("Person-level rows are summarised only; they are never joined to places or facilities and never "
                       "exposed individually.")
        env = _env("trace_evidence", OK if st == "ok" else UNKNOWN, data, caveats=cav, truncation=trunc,
                   reason=r.get("reason") if st != "ok" else None)
        # provenance of a trace is the trace itself: sources are listed in data.sources
        env["provenance"]["sources"] = [{"source_id": s.get("source_id"), "name": s.get("name"),
                                         "source_version": s.get("source_version"),
                                         "retrieved_at": s.get("retrieved_at")}
                                        for s in (data.get("sources") or []) if isinstance(s, dict)]
        return env
    return _run("trace_evidence", {"object_id": object_id, "max_rows": max_rows}, go)


SOURCE_KEEP = ("source_id", "name", "publisher", "data_layer", "status", "source_version", "retrieved_at",
               "update_date", "license", "access_conditions", "landing_url", "unit_of_observation",
               "geographic_resolution", "person_level", "geographic", "omics", "wearable",
               "true_participant_linkage_across_modalities", "processed_outputs", "audit", "ingestion_module",
               "limitations")


def list_sources(data_layer: str | None = None) -> dict:
    """SOURCE_REGISTRY.yaml summary: every public source with publisher, version, retrieval date, layer, status,
    licence, processed outputs and its DATA_AUDIT.md path."""
    def go():
        reg = registry()
        srcs = []
        for s in reg["by_id"].values():
            if data_layer and s.get("data_layer") != data_layer:
                continue
            srcs.append({k: s.get(k) for k in SOURCE_KEEP})
        layers: dict[str, int] = {}
        for s in reg["by_id"].values():
            layers[str(s.get("data_layer"))] = layers.get(str(s.get("data_layer")), 0) + 1
        data = {"n_sources": len(srcs), "by_data_layer": layers, "registry": str(SOURCE_REGISTRY_PATH.name),
                "registry_generated_at": reg.get("generated_at"), "sources": srcs}
        if not srcs:
            return _env("list_sources", UNKNOWN, data, reason=f"no source with data_layer {data_layer!r}; layers: "
                                                               f"{sorted(layers)}")
        return _env("list_sources", OK, data, caveats=[
            "Four data layers (person, condition_molecular, geographic, facility) describe different units and are "
            "never joined as if they were the same people."])
    return _run("list_sources", {"data_layer": data_layer}, go)


TOOLS = {name: globals()[name] for name in SPEC_TOOLS + EXTRA_TOOLS}


# ============================================================================================ benchmark

BENCHMARK_CALLS: list[tuple[str, dict]] = [
    ("search_condition", {"condition": "long covid"}),
    ("normalize_condition", {"condition_or_code": "G93.32"}),
    ("get_patient_phenotype_signature", {"condition": "ME/CFS"}),
    ("get_molecular_context", {"condition": "Long COVID"}),
    ("get_measurable_biology", {"condition": "ME/CFS"}),
    ("discover_candidate_measurements", {"condition": "Long COVID", "phenotype": "orthostatic intolerance"}),
    ("get_measurement_evidence", {"measurement": "wearable autonomic monitoring", "condition": "Long COVID"}),
    ("get_phenotype_measurement_evidence", {"phenotype": "orthostatic/autonomic dysfunction"}),
    ("get_regulatory_context", {"technology": "wearable ECG patch"}),
    ("get_condition_burden", {"condition": "Long COVID", "geography_level": "county"}),
    ("get_geographic_context", {"geography_id": "San Diego County, California"}),
    ("find_candidate_clinics", {"geography_id": "06073", "condition": "Long COVID",
                                "measurement": "wearable autonomic monitoring"}),
    ("find_relevant_trials", {"condition": "Long COVID", "measurement": "wearable autonomic monitoring"}),
    ("find_relevant_research_centers", {"condition": "ME/CFS"}),
    ("rank_deployment_opportunities", {"condition": "Long COVID or ME/CFS",
                                       "measurement": "wearable autonomic monitoring", "geography_level": "county"}),
    ("search_us_open_data", {"query": "long COVID"}),
    ("trace_evidence", {"object_id": "opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|40109"}),
    ("list_sources", {}),
]


def benchmark(calls: list[tuple[str, dict]] | None = None) -> list[dict]:
    """Cold (first call in this process, tool cache empty) and warm (memoised) latency of each tool."""
    out = []
    clear_cache()
    for name, args in calls or BENCHMARK_CALLS:
        t0 = time.perf_counter()
        r1 = TOOLS[name](**args)
        cold = time.perf_counter() - t0
        t0 = time.perf_counter()
        TOOLS[name](**args)
        warm = time.perf_counter() - t0
        out.append({"tool": name, "args": args, "status": r1["status"], "cold_s": round(cold, 3),
                    "warm_s": round(warm, 4), "bytes": len(json.dumps(r1))})
    return out


def _main(argv: list[str]) -> int:
    if "--benchmark" in argv:
        rows = benchmark()
        w = max(len(r["tool"]) for r in rows)
        print(f"{'tool':{w}s}  {'status':24s} {'cold_s':>8s} {'warm_s':>8s} {'kB':>7s}")
        for r in rows:
            print(f"{r['tool']:{w}s}  {r['status']:24s} {r['cold_s']:8.3f} {r['warm_s']:8.4f} {r['bytes'] / 1e3:7.1f}")
        return 0
    if len(argv) >= 1 and argv[0] in TOOLS:
        import contextlib
        args = json.loads(argv[1]) if len(argv) > 1 else {}
        with contextlib.redirect_stdout(sys.stderr):   # progress lines of the layers must not corrupt the JSON
            res = TOOLS[argv[0]](**args)
        print(json.dumps(res, indent=1))
        return 0
    print("usage: python -m measure_it.tools --benchmark | <tool> '<json args>'", file=sys.stderr)
    return 2


if __name__ == "__main__":  # pragma: no cover
    code = _main(sys.argv[1:])
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(code)   # skip native-library teardown (see measure_it.omics.query.exit_cleanly)
