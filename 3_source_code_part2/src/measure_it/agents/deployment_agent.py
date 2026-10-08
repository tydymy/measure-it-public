"""Deterministic, LLM-free deployment agent: answer_question(question) -> templated, cited answer.

    uv run measure-it ask "Given phenotype orthostatic intolerance and candidate measurement technology nailfold
                           capillaroscopy, where should we deploy it?"

It answers "what should we measure, where should we measure it, and who could realistically deploy it?" with the
tools of measure_it.tools (the same facade the MCP server exposes). No language model is involved: the question is
parsed with the project's normalisers, the tools are called in a fixed, documented order, and the answer is filled
from tool output into fixed sentence templates. Every factual sentence ends with the object ids it rests on
([condition:...], [opportunity:...], [facility:...] ...); anything a tool returns as UNKNOWN stays UNKNOWN.

Parsing (deterministic)
    conditions   every span of the question that matches a registry term (preferred names, aliases, search terms,
                 Mondo labels/synonyms of measure_it.ontology; acronyms such as POTS, ME/CFS case-sensitively), longest
                 non-overlapping spans first; each span is then passed to normalize_condition, and only its 'matched'
                 condition is used. No span -> search_condition on the phrase after 'for'/'in'/'with' (or on the whole
                 question); nothing found -> the condition is UNKNOWN and nothing is ranked.
    phenotype    spans matching a phenotype axis (phenotype_axes) or a Phase 3 demo phenotype label/alias.
    measurement  spans matching a measurement class id/name (parts split at '/' and brackets) or a bundle
                 id/label/alias of configs (measurements.yaml, relevance.yaml), confirmed with
                 facilities.matching.resolve_measurement; '<word>-based' is read as '<word>'. None named -> the top
                 candidate of discover_candidate_measurements for the first condition.
    geography    '<Name> County|Parish|Borough|..., <State>', 'in <State name>', or a FIPS / geo: id, resolved with
                 geography.query.resolve_geography; otherwise the whole U.S. The level is 'state' when the question
                 asks for states, else 'county'. 'N regions/counties/states' sets how many regions are listed (default
                 5).
    condition set  all parsed conditions together: a known condition set when the members match (e.g. the demo
                 cluster), otherwise 'A or B ...' (an ad hoc set scored by the same engine).

Tool order
    1 normalize_condition (each condition/phenotype span) | search_condition (fallback)
    2 get_phenotype_measurement_evidence (named phenotype)
    3 get_patient_phenotype_signature (each condition)
    4 discover_candidate_measurements (each condition; phenotype filter when named)
    5 get_measurement_evidence (measurement x each condition); get_regulatory_context (each member class)
    6 get_molecular_context (each condition)
    7 get_condition_burden (each condition, at the level)
    8 rank_deployment_opportunities (condition set x measurement x level). The ranking is national: for a named
      state at county level the national top 50 is requested and the state's counties among them are listed
      separately (a within-state ranking is UNKNOWN - the tool has no state filter)
    9 named geography: get_geographic_context + trace_evidence(opportunity id; its level is stated - a named state's
      stored row is a rank among states); find_candidate_clinics (named geography - a state's pool is the whole
      state - else the top-ranked region; each condition)
   10 find_relevant_research_centers and find_relevant_trials (each condition x measurement)
   11 trace_evidence (top opportunity -> component rows -> sources and raw files)
"""
from __future__ import annotations

import json
import re
import time
from functools import lru_cache

from .. import tools as T
from ..config import UNKNOWN, load_config

DEFAULT_TOP_N = 5
STATE_WINDOW = 50   # a named state at county level: its counties are looked up among the national top STATE_WINDOW
PERSON_NOTE = ("Person-level results are group statistics from public cohorts; those participants are not the people "
               "of any region, facility or trial named below.")


# ============================================================================================ parsing helpers

def _fold(text: str) -> str:
    """Lower-case ASCII with every other character replaced by a space: same length as the input (indices match)."""
    return "".join(c.lower() if c.isascii() else " " for c in str(text))


def _span_regex(tokens: list[str]) -> re.Pattern:
    return re.compile(r"(?<![a-z0-9])" + r"[^a-z0-9]+".join(map(re.escape, tokens)) + r"(?![a-z0-9])")


def _tokens(text: str) -> list[str]:
    return re.sub(r"[^a-z0-9]+", " ", _fold(text)).split()


def _prep(question: str) -> str:
    """'wearable-based' -> 'wearable      ' (same length)."""
    return re.sub(r"(?i)(\w)-based\b", lambda m: m.group(1) + " " * 6, question)


@lru_cache(maxsize=1)
def _condition_lexicon() -> tuple:
    from ..ontology.normalize import get_index, is_acronym
    out, seen = [], set()
    for e in get_index().lexicon:
        text = str(e.text).strip()
        key = (text, e.condition_id)
        if key in seen or not text:
            continue
        seen.add(key)
        if is_acronym(text):
            out.append((text, e.condition_id, re.compile(r"(?<![A-Za-z0-9])" + re.escape(text) + r"(?![A-Za-z0-9])"),
                        True))
        else:
            toks = _tokens(re.sub(r"\([^)]*\)", " ", text))
            if not toks or sum(len(t) for t in toks) < 4:
                continue
            out.append((text, e.condition_id, _span_regex(toks), False))
    return tuple(out)


@lru_cache(maxsize=1)
def _measurement_lexicon() -> tuple:
    out = []

    def add(text, target, prio):
        toks = _tokens(text)
        if toks and sum(len(t) for t in toks) >= 3:
            out.append((text, target, _span_regex(toks), prio))

    # the shared vocabulary of measure_it.measurements.resolve (class ids/names/aliases, bundle ids/labels/aliases);
    # the spans found here are confirmed with the shared resolver (facilities.matching.resolve_measurement)
    from ..measurements.resolve import vocabulary
    for text, target, prio in vocabulary():
        add(re.sub(r"\([^)]*\)", " ", text.replace("_", " ")), target, prio)
    for m in load_config("measurements")["measurement_classes"]:
        for part in re.split(r"[/()]", str(m.get("name") or "")):
            if len(part.strip()) >= 4:
                add(part, m["id"], 3)
    return tuple(out)


@lru_cache(maxsize=1)
def _phenotype_lexicon() -> tuple:
    out = []
    try:
        from ..store import read_table
        ax = read_table("phenotype_axes", columns=["axis_id", "axis_label", "hpo_label"])
        from ..measurements.query import phenotype_axis_aliases
        aliases = phenotype_axis_aliases()
        for r in ax.itertuples():
            for t in (r.axis_label, r.hpo_label, str(r.axis_id).replace("_", " ")):
                if isinstance(t, str) and len(t) >= 5:
                    out.append((t, "axis", r.axis_id, _span_regex(_tokens(t))))
            for t in aliases.get(r.axis_id, []):     # curated aliases (e.g. PEM), resolved as exact axis matches
                if _tokens(t):
                    out.append((t, "axis", r.axis_id, _span_regex(_tokens(t))))
    except FileNotFoundError:
        pass
    from ..measurements.phenotype_evidence import PHENOTYPES
    for pid, p in PHENOTYPES.items():
        for t in [p["label"]] + list(p.get("aliases", [])):
            out.append((t, "demo_phenotype", pid, _span_regex(_tokens(t))))
    return tuple(out)


def _longest_spans(matches: list[tuple[int, int, object]]) -> list[tuple[int, int, object]]:
    """Greedy longest-first selection of non-overlapping (start, end, payload) spans, returned in text order."""
    keep: list[tuple[int, int, object]] = []
    for s, e, p in sorted(matches, key=lambda x: (-(x[1] - x[0]), x[0])):
        if all(e <= ks or s >= ke for ks, ke, _ in keep):
            keep.append((s, e, p))
    return sorted(keep, key=lambda x: x[0])


def parse_conditions(question: str) -> list[dict]:
    q = _prep(question)
    f = _fold(q)
    hits = []
    for text, cid, rx, acr in _condition_lexicon():
        for m in rx.finditer(q if acr else f):
            hits.append((m.start(), m.end(), {"span": question[m.start():m.end()], "lexicon_term": text,
                                              "lexicon_condition": cid}))
    return [dict(p, start=s, end=e) for s, e, p in _longest_spans(hits)]


def parse_phenotypes(question: str) -> list[dict]:
    f = _fold(_prep(question))
    hits = []
    for term, kind, pid, rx in _phenotype_lexicon():
        for m in rx.finditer(f):
            hits.append((m.start(), m.end(), {"span": question[m.start():m.end()], "kind": kind, "id": pid,
                                              "lexicon_term": term}))
    return [dict(p, start=s, end=e) for s, e, p in _longest_spans(hits)]


def parse_measurements(question: str) -> list[dict]:
    f = _fold(_prep(question))
    hits = []
    for term, target, rx, prio in _measurement_lexicon():
        for m in rx.finditer(f):
            hits.append((m.start(), m.end(), {"span": question[m.start():m.end()], "lexicon_term": term,
                                              "lexicon_target": target, "priority": prio}))
    return [dict(p, start=s, end=e) for s, e, p in _longest_spans(hits)]


COUNTY_RE = re.compile(r"((?:[A-Z][A-Za-z.'\-]*\s+)+(?:County|Parish|Borough|Census Area|Planning Region|Municipality|"
                       r"City and Borough|city)),\s*([A-Z][A-Za-z]+(?:\s(?:of\s)?[A-Z][A-Za-z]+)*)")
FIPS_RE = re.compile(r"\bgeo:(\d{2}|\d{5})\b|\bFIPS\s*(\d{2}|\d{5})\b")


def parse_geography(question: str) -> dict | None:
    from ..geography.query import resolve_geography
    for m in COUNTY_RE.finditer(question):
        text = f"{m.group(1).strip()}, {m.group(2).strip()}"
        r = resolve_geography(text)
        if r.get("status") == "ok":
            return {"span": text, **r}
    m = FIPS_RE.search(question)
    if m:
        r = resolve_geography(m.group(1) or m.group(2))
        if r.get("status") == "ok":
            return {"span": m.group(0), **r}
    for name in sorted(_state_names(), key=len, reverse=True):
        mm = re.search(r"\b(?:in|within|across)\s+(?:the\s+state\s+of\s+)?" + re.escape(name) + r"\b", question)
        if mm:
            r = resolve_geography(name)
            if r.get("status") == "ok":
                return {"span": name, **r}
    return None


@lru_cache(maxsize=1)
def _state_names() -> tuple:
    from ..store import read_table
    g = read_table("geographies", columns=["geo_level", "name"])
    return tuple(sorted(set(g.loc[g["geo_level"] == "state", "name"].astype(str))))


def parse_level(question: str) -> str:
    return "state" if re.search(r"\b(which|what|top|\d+)\s+(U\.S\.\s+)?states\b|\bstate[- ]level\b", question, re.I) \
        else "county"


def parse_top_n(question: str) -> int:
    m = re.search(r"\b(\d{1,2})\s+(?:U\.S\.\s+|US\s+)?(?:regions|counties|states|places|locations|areas|sites)\b",
                  question, re.I)
    return max(1, min(int(m.group(1)), 25)) if m else DEFAULT_TOP_N


FALLBACK_PHRASE_RE = re.compile(r"\b(?:for|in|with|of)\s+(?:patients\s+with\s+)?([A-Za-z][A-Za-z0-9'/\- ]{2,60}?)"
                                r"(?=[?.,;:!]|\s+(?:and|or|where|which|in|at|to|using|with)\b|$)")

# Premises the engine cannot answer. A question that carries one is still answered as a deployment question, but the
# answer first says (fixed text, no model) what the engine does not do, so the premise is never silently accepted.
SCOPE_RULES: list[tuple[re.Pattern, str]] = [
    (re.compile(r"\b(?:does|do|did|is|has|have|can)\s+(?:the|this|that|my|our|a)\s+(?:patient|person|participant|"
                r"individual|child|son|daughter|wife|husband|partner)\b|\bpatient\s+(?:has|have)\b|"
                r"\b(?:i|we|he|she)\s+(?:have|has)\s+(?:it|this|that|long|me|pots|a)\b|\bdiagnos(?:e|es|ed|ing|is)\b",
                re.I),
     "The engine does not assess, diagnose or locate any individual: person-level public cohorts are summarised as "
     "group statistics, and places and facilities are ranked only as candidate deployment opportunities."),
    (re.compile(r"\b(?:best|optimal|top[- ]rated|leading)\b", re.I),
     "No clinic, doctor or site is ranked as best: facilities are listed as candidates with characteristics "
     "suggesting they may be viable implementation or study partners, with the reasons for each."),
    (re.compile(r"\b(?:multi-?omics|omics\s+profile|genome|genomes|genetic\s+profile)\b", re.I),
     "No public dataset used here links omics to wearable participants: molecular context is condition-level molecular "
     "enrichment from public databases, not patient multi-omics."),
    (re.compile(r"\b(?:causes?|caused|causal|mechanisms?|cure[sd]?|treat(?:ment|ments|ed|ing)?|therap(?:y|ies))\b",
                re.I),
     "The engine reports measurement and deployment opportunities; it does not establish causes or mechanisms and "
     "does not recommend treatments."),
]
PERSON_DATASET_RE = re.compile(r"\b(?:nhanes|stanford|mapmecfs|participants?|respondents?|subjects?)\b", re.I)
PERSON_LOCATION_NOTE = ("Person-level participants (NHANES, the Stanford wearable studies, MapMECFS) are never "
                        "located: their public files carry no geography below the nation, so no participant is placed "
                        "in any county, state or facility.")


def scope_notes(question: str, geography: dict | None = None) -> list[str]:
    """Fixed scope statements for premises the engine cannot answer (individual diagnosis, 'best' clinic, participant
    location or omics, causes / treatment), in rule order."""
    q = str(question or "")
    out = [note for rx, note in SCOPE_RULES if rx.search(q)]
    asks_where = re.search(r"\b(?:liv(?:e|es|ing)|located|reside|where)\b", q, re.I)
    if PERSON_DATASET_RE.search(q) and (geography or asks_where):
        out.append(PERSON_LOCATION_NOTE)
    return out


# ============================================================================================ agent

class _Log:
    def __init__(self):
        self.calls: list[dict] = []

    def call(self, step: str, tool: str, **args) -> dict:
        t0 = time.perf_counter()
        r = T.TOOLS[tool](**args)
        self.calls.append({"step": step, "tool": tool, "args": args, "status": r.get("status"),
                           "reason": (r.get("reason") or "")[:300] or None,
                           "elapsed_s": round(time.perf_counter() - t0, 3),
                           "cache_hit": (r.get("meta") or {}).get("cache_hit")})
        return r


def _ids(*xs) -> str:
    flat = []
    for x in xs:
        if isinstance(x, (list, tuple)):
            flat += [str(i) for i in x if i]
        elif x:
            flat.append(str(x))
    flat = list(dict.fromkeys(flat))
    return " [" + "; ".join(flat) + "]" if flat else ""


def _fmt(v, nd: int = 2) -> str:
    if v is None or v == UNKNOWN:
        return UNKNOWN
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int) or (isinstance(v, float) and float(v).is_integer() and abs(v) >= 1):
        return f"{int(v):,}"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def _condition_query(cids: list[str]) -> str:
    from ..scoring.opportunity import condition_sets
    for sid, s in condition_sets().items():
        if set(s["members"]) == set(cids) and len(cids) > 1:
            return sid
    return cids[0] if len(cids) == 1 else " or ".join(cids)


def answer_question(question: str, *, top_n: int | None = None, geography_level: str | None = None) -> dict:
    """Answer 'what should we measure, where, and who could deploy it?' from tool output only (no LLM)."""
    log = _Log()
    q = str(question or "").strip()
    unknowns: list[str] = []
    lines: dict[str, list[str]] = {"parse": [], "what": [], "where": [], "who": [], "uncertain": [], "prov": []}

    # ------------------------------------------------------------------ 1. parse + normalise
    level = geography_level or parse_level(q)
    n_regions = top_n or parse_top_n(q)
    phen_spans = parse_phenotypes(q)
    cond_spans = parse_conditions(q)
    meas_spans = parse_measurements(q)
    geo = parse_geography(q)
    lines["scope"] = scope_notes(q, geo)
    conditions: list[dict] = []
    for sp in cond_spans:
        r = log.call("1 normalise", "normalize_condition", condition_or_code=sp["span"])
        d = r.get("data") or {}
        if r["status"] == T.OK and d.get("match_status") == "matched":
            m = d["matches"][0]
            cid = m["canonical_condition_id"]
            if cid not in [c["condition_id"] for c in conditions]:
                conditions.append({"condition_id": cid, "preferred_name": m.get("preferred_name"), "span": sp["span"],
                                   "match_reason": m.get("match_reason"),
                                   "other_candidates": [c["canonical_condition_id"] for c in d.get("other_candidates", [])],
                                   "phenotype_axis": (d.get("phenotype_axis") or {}).get("axis_id")})
        else:
            unknowns.append(f"'{sp['span']}' did not normalise to one registry condition "
                            f"(normalize_condition: {d.get('match_status') or r['status']}).")
    if not conditions:
        phrase = None
        m = FALLBACK_PHRASE_RE.search(q)
        if m:
            phrase = m.group(1).strip()
        phrase = phrase or q
        r = log.call("1 normalise", "search_condition", condition=phrase)
        cands = (r.get("data") or {}).get("candidates") or []
        if r["status"] == T.OK and cands:
            lines["parse"].append(f"No condition in the question matched a registry term exactly; search_condition "
                                  f"on '{phrase}' returned candidates " + ", ".join(
                                      f"{c['canonical_condition_id']} ({c['match_reason']})" for c in cands[:3])
                                  + " - not used without an exact match." + _ids([c["object_id"] for c in cands[:3]]))
        unknowns.append(f"Condition: {UNKNOWN} - no registry condition matches '{phrase}' (normalize_condition / "
                        f"search_condition).")
    measurement = None
    for sp in meas_spans:
        from ..facilities.matching import resolve_measurement
        bundle = _bundle_for_phrase(sp["lexicon_term"])
        rm = resolve_measurement(sp["lexicon_term"])
        if bundle:
            mid, how = bundle, f"exact bundle id / label / alias '{sp['lexicon_term']}' (bundles are preferred for deployment)"
        elif rm.get("status") == "matched":
            mid, how = rm["measurement_id"], rm.get("match_reason")
        else:
            mid, how = sp["lexicon_target"], f"lexicon term '{sp['lexicon_term']}' (resolve_measurement: {rm.get('reason')})"
        measurement = {"measurement_id": mid, "span": sp["span"], "match": how}
        break
    phenotype = phen_spans[0] if phen_spans else None

    parsed = {"conditions": conditions, "phenotype": phenotype, "measurement": measurement,
              "geography": ({k: geo.get(k) for k in ("span", "geo_id", "geo_level", "name", "object_id")}
                            if geo else {"scope": "United States (no specific geography named)"}),
              "level": level, "n_regions": n_regions}
    if conditions:
        lines["parse"].append("Conditions: " + "; ".join(
            f"'{c['span']}' -> {c['preferred_name']} ({c['match_reason']}"
            + (f"; other candidates {c['other_candidates']}" if c["other_candidates"] else "") + ")"
            for c in conditions) + _ids([f"condition:{c['condition_id']}" for c in conditions]))
    if phenotype:
        lines["parse"].append(f"Phenotype: '{phenotype['span']}' ({phenotype['kind'].replace('_', ' ')} "
                              f"{phenotype['id']}).")
    lines["parse"].append(("Measurement: '" + measurement["span"] + "' -> " + measurement["measurement_id"]
                           + f" ({measurement['match']})" + _ids(f"measurement:{measurement['measurement_id']}"))
                          if measurement else "Measurement: none recognised in the question (no measurement class, "
                                              "bundle or alias matched; a device or brand the vocabulary does not "
                                              "know, e.g. 'smartwatch', is not guessed).")
    lines["parse"].append(f"Geography: {geo['name']} ({geo['geo_level']})" + _ids(geo["object_id"]) if geo
                          else f"Geography: the United States; regions ranked at {level} level.")

    if not conditions:
        return _finish(q, parsed, log, lines, unknowns, stopped="no condition could be resolved, so no measurement, "
                                                                "burden or deployment result can be attributed to it")

    cids = [c["condition_id"] for c in conditions]

    # ------------------------------------------------------------------ 2. phenotype (Phase 3)
    p3 = None
    if phenotype:
        p3 = log.call("2 phenotype", "get_phenotype_measurement_evidence", phenotype=phenotype["span"])
        if p3["status"] == T.OK:
            d = p3["data"]
            ph = d.get("phenotype") or {}
            techs = d.get("technologies") or []
            obs = d.get("observable_in_public_person_level_data") or []
            lines["what"].append(
                f"Phenotype '{phenotype['span']}' belongs to the demo phenotype '{ph.get('label')}' "
                f"({d.get('match_reason')}); its observable signals are "
                + ", ".join(map(str, ph.get("observable_signals") or [])) + "."
                + _ids([t.get("object_id") for t in techs[:3]]))
            lines["what"].append(
                f"Of its {len(techs)} candidate technologies, public person-level data give a signal result for "
                + (", ".join(obs) if obs else "none") + "; the others are "
                + UNKNOWN + " in public data." + _ids([t.get("object_id") for t in techs if t.get("measurement_id")
                                                       in obs] or [t.get("object_id") for t in techs[:1]]))
            if measurement:
                members = _members(measurement["measurement_id"])
                linked = [t for t in techs if t.get("measurement_id") in members + [measurement["measurement_id"]]]
                if linked:
                    lines["what"].append(f"{measurement['measurement_id']} is among the phenotype's technologies "
                                         f"({', '.join(t['measurement_id'] for t in linked)})."
                                         + _ids([t.get('object_id') for t in linked]))
                else:
                    lines["what"].append(
                        f"{measurement['measurement_id']} is NOT linked to this phenotype in the phenotype-axis map: "
                        f"the Phase 3 table has no row for it, so its phenotype-level evidence is {UNKNOWN}."
                        + _ids([t.get("object_id") for t in techs[:2]]))
                    unknowns.append(f"Phenotype-level evidence for {measurement['measurement_id']} x "
                                    f"'{phenotype['span']}': {UNKNOWN} (not linked in the phenotype-axis map).")
        else:
            unknowns.append(f"Phase 3 evidence for phenotype '{phenotype['span']}': {UNKNOWN} ({p3.get('reason')}).")

    # ------------------------------------------------------------------ 3. person-level signatures
    for c in conditions:
        r = log.call("3 signature", "get_patient_phenotype_signature", condition=c["condition_id"])
        if r["status"] != T.OK:
            unknowns.append(f"Wearable phenotype signature for {c['condition_id']}: {UNKNOWN} "
                            f"({(r.get('reason') or '')[:160]}...).")
            lines["what"].append(f"Public person-level wearable signature for {c['preferred_name']}: {UNKNOWN} - "
                                 f"{(r.get('reason') or '')[:200]}" + _ids(f"condition:{c['condition_id']}"))
            continue
        # wearable datasets first (the question is about measuring), then molecular person-level datasets
        dsets = sorted(r["data"].get("datasets") or [],
                       key=lambda d: not str(d.get("dataset_id") or "").startswith(WEARABLE_SIGNATURE_DATASETS))
        for d in dsets[:3]:
            top = (d.get("top_features") or [])[:2]
            proxy = " (a PROXY label, not the condition itself)" if d.get("is_proxy") is True else ""
            base = (f"{d.get('dataset_id')}: phenotype '{d.get('phenotype_label')}'{proxy}, "
                    f"{_fmt(d.get('n_cases'))} cases vs {_fmt(d.get('n_controls'))} non-cases; "
                    f"{_fmt(d.get('n_features_fdr_significant'))} of {_fmt(d.get('n_features_tested'))} "
                    f"{_feature_kind(d.get('dataset_id'))} passed BH-FDR")
            if "CONFOUND" in str(d.get("caveats") or "").upper():
                base += " (the analysis flags a possible confound; see the signature's caveats)"
            if top:
                base += "; e.g. " + "; ".join(
                    f"{f['feature']} {_fmt(f.get('effect_size'))} ({_fmt(f.get('effect_measure'))}, 95% CI "
                    f"{_fmt(f.get('ci_low'))} to {_fmt(f.get('ci_high'))})" for f in top)
            lines["what"].append(f"Measurable phenotype in public data for {c['preferred_name']} - " + base + "."
                                 + _ids([f.get("object_id") for f in top] or [f"condition:{c['condition_id']}"]))
    lines["what"].append(PERSON_NOTE)

    # ------------------------------------------------------------------ 4. candidate measurements
    first_candidate = None
    for c in conditions:
        r = log.call("4 candidates", "discover_candidate_measurements", condition=c["condition_id"],
                     phenotype=(phenotype["span"] if phenotype else None), top_n=5)
        if r["status"] != T.OK:
            unknowns.append(f"Candidate measurements for {c['condition_id']}"
                            + (f" x phenotype '{phenotype['span']}'" if phenotype else "") + f": {UNKNOWN} "
                            f"({(r.get('reason') or '')[:200]}).")
            continue
        cands = r["data"].get("candidates") or []
        if cands and first_candidate is None:
            first_candidate = cands[0]
        if cands:
            lines["what"].append(
                f"Candidate objective measurements for {c['preferred_name']} (ordered by registered research activity, "
                f"not by validity): " + "; ".join(
                    f"{x['measurement_id']} ({_mes_counts(x)})" for x in cands[:4]) + "."
                + _ids([x.get("evidence_object_id") for x in cands[:4]]))
    if measurement is None and first_candidate is not None:
        measurement = {"measurement_id": first_candidate["measurement_id"], "span": None,
                       "match": "top candidate of discover_candidate_measurements (none recognised in the question)"}
        parsed["measurement"] = measurement
        lines["what"].append("No measurement was recognised in the question, so the top candidate, "
                             f"{first_candidate['measurement_id']}, is used below."
                             + _ids(first_candidate.get("evidence_object_id")))
    if measurement is None:
        unknowns.append(f"Measurement: {UNKNOWN} - none named and no candidate found.")
        return _finish(q, parsed, log, lines, unknowns, stopped="no measurement to deploy")
    mid = measurement["measurement_id"]

    # ------------------------------------------------------------------ 5. measurement evidence + regulatory
    for c in conditions:
        r = log.call("5 evidence", "get_measurement_evidence", measurement=mid, condition=c["condition_id"],
                     max_object_ids=5)
        if r["status"] == T.OK:
            evs = r["data"].get("evidence") or []
            for e in ([x for x in evs if x.get("measurement_id") == mid] or evs)[:1]:
                cnt = e.get("counts") or {}
                lines["what"].append(
                    f"{e.get('name')} for {c['preferred_name']}: {_fmt(cnt.get('n_trials_objective'))} registered "
                    f"trial(s) describe objective use ({_fmt(cnt.get('n_trials_outcome_measure'))} as an outcome "
                    f"measure) and {_fmt(cnt.get('n_nih_core_projects'))} NIH core project(s) mention it; evidence "
                    f"tier {((e.get('dimensions') or {}).get('measurement_evidence_strength') or {}).get('tier')} "
                    "(research activity, not proof that it works)."
                    + _ids(e.get("evidence_object_id"), (e.get("trial_object_ids") or [])[:2]))
        else:
            lines["what"].append(f"Research evidence for {mid} in {c['preferred_name']}: {UNKNOWN} - "
                                 f"{(r.get('reason') or '')[:220]}" + _ids(f"condition:{c['condition_id']}",
                                                                          f"measurement:{mid}"))
            unknowns.append(f"Trial/NIH evidence for {mid} x {c['condition_id']}: {UNKNOWN}.")
    for m in _members(mid)[:4]:
        r = log.call("5 regulatory", "get_regulatory_context", technology=m, max_records=3)
        if r["status"] == T.OK:
            s = next((x for x in r["data"].get("summary") or [] if x.get("measurement_id") == m), None)
            if s:
                codes = [row.get("product_code_object_id") for row in r["data"].get("regulatory_status") or []
                         if isinstance(row, dict) and row.get("product_code_object_id")]
                lines["what"].append(
                    f"FDA context for {m}: {s.get('regulatory_visibility')}; {_fmt(s.get('n_mapped_product_codes'))} "
                    f"mapped product code(s), {_fmt(s.get('n_510k_total'))} 510(k), {_fmt(s.get('n_denovo_total'))} "
                    f"De Novo, {_fmt(s.get('n_pma_total'))} PMA decisions (a deployment-readiness signal, not "
                    "evidence that it detects the condition)." + _ids(f"measurement:{m}", codes[:2]))
        else:
            unknowns.append(f"FDA regulatory context for {m}: {UNKNOWN} ({r.get('reason')}).")

    # ------------------------------------------------------------------ 6. molecular (condition-level)
    for c in conditions:
        r = log.call("6 molecular", "get_molecular_context", condition=c["condition_id"], top_n=3)
        if r["status"] == T.OK:
            ev = r["data"].get("evidence") if isinstance(r["data"].get("evidence"), dict) else {}
            parts = [f"{k} {_fmt(v.get('n_records'))} records" for k, v in list(ev.items())[:4]]
            sysd = (r["data"].get("coherence") or {}).get("supported_systems")
            tiers = system_support_tiers(sysd)
            lines["what"].append(
                f"Condition-level molecular enrichment for {c['preferred_name']} (public databases; not patient "
                f"multi-omics, not mechanism): " + (", ".join(parts) or UNKNOWN) + "; " + describe_system_tiers(tiers)
                + _ids(f"condition:{c['condition_id']}"))
        else:
            unknowns.append(f"Condition-level molecular enrichment for {c['condition_id']}: {UNKNOWN}.")

    # ------------------------------------------------------------------ 7. burden
    for c in conditions:
        r = log.call("7 burden", "get_condition_burden", condition=c["condition_id"], geography_level=level,
                     max_rows=3)
        d = r.get("data") or {}
        lvl = d.get("burden_evidence_level", UNKNOWN)
        if r["status"] == T.OK:
            sm = d.get("summary") or {}
            inh = (" These are STATE estimates inherited by counties (source resolution state), not county "
                   "prevalence." if sm.get("n_inherited_state_values") else "")
            lines["where"].append(
                f"Burden of {c['preferred_name']} at {level} level: evidence level {lvl} "
                f"({_level_word(lvl)}), measure {d.get('primary_measure_id')} at {d.get('primary_source_resolution')} "
                f"resolution; {_fmt(sm.get('n_with_value'))} of {_fmt(sm.get('n_rows'))} geographies have a value "
                f"(range {_fmt(sm.get('min'))}-{_fmt(sm.get('max'))} {', '.join(sm.get('value_units') or [])})."
                + inh + _ids(f"condition:{c['condition_id']}"))
        else:
            lines["where"].append(
                f"Burden of {c['preferred_name']} at {level} level: {UNKNOWN} (evidence level {lvl}) - "
                f"{(r.get('reason') or '')[:260]}" + _ids(f"condition:{c['condition_id']}"))
            unknowns.append(f"Burden of {c['condition_id']} ({level}): {UNKNOWN}, evidence level {lvl}.")

    # ------------------------------------------------------------------ 8. ranking
    cq = _condition_query(cids)
    # a named state with a county-level question: the ranking tool ranks nationally (it has no within-state filter),
    # so ask for the national top STATE_WINDOW and report the named state's counties among them separately
    state_filter = geo.get("state_abbr") if geo and geo.get("geo_level") == "state" and level == "county" else None
    rk = log.call("8 rank", "rank_deployment_opportunities", condition=cq, measurement=mid, geography_level=level,
                  top_n=max(n_regions, STATE_WINDOW) if state_filter else n_regions, max_sites_per_region=3)
    recs = []
    recs_all: list[dict] = []
    stored = False
    if rk["status"] == T.OK:
        rd = rk["data"]
        recs_all = rd.get("recommendations") or []
        recs = recs_all[:n_regions]
        stored = str(rd.get("ranking_basis", "")).startswith("deployment_opportunities")
        if not stored:
            lines["where"].append(
                f"This combination is not precomputed, so the ranking was computed on the fly by the same engine "
                f"({rd.get('ranking_basis')}); its regions have no stored opportunity row, so they are cited by "
                "geography, condition and measurement ids." + _ids(f"measurement:{mid}", [f"condition:{c}" for c in cids]))
        sc = rd.get("shared_context") or {}
        cset = (sc.get("condition") or {}).get("label") or cq
        lines["where"].append(
            f"Candidate deployment opportunities for {cset} x {(sc.get('measurement') or {}).get('label') or mid} "
            f"({level} level, '{rd.get('weight_set')}' weights; {_fmt(rd.get('n_regions_ranked'))} regions ranked; "
            f"{rd.get('ranking_universe')}"
            + (f"; a NATIONAL ranking - the ranking tool has no within-state filter, {geo['name']} counties are "
               "listed after it" if state_filter else "") + "):"
            + _ids((sc.get("condition") or {}).get("object_ids"), f"measurement:{mid}"))
        for rec in recs:
            g = rec.get("geography") or {}
            ri = g.get("rank_interval_5_95") or [None, None]
            b = g.get("burden") or {}
            opp = _opp_id(rec) if stored else None
            lines["where"].append(
                f"  {_fmt(g.get('rank'))}. {g.get('name')}: composite {_fmt(g.get('composite'), 3)}, rank interval "
                f"{_fmt(ri[0])}-{_fmt(ri[1])} (P(top 10) {_fmt(g.get('p_top10_monte_carlo'))}); burden percentile "
                f"{_fmt(b.get('percentile'))} (evidence level {g.get('burden_evidence_level')}"
                + (", includes a state value inherited by the county" if b.get("inherited") else "") + "), "
                f"vulnerability pct {_fmt((g.get('vulnerability') or {}).get('percentile'))}, diagnostic desert pct "
                f"{_fmt((g.get('diagnostic_desert') or {}).get('percentile'))}, clinic capacity pct "
                f"{_fmt((g.get('clinic_capacity') or {}).get('percentile'))}, research readiness pct "
                f"{_fmt((g.get('research_readiness') or {}).get('percentile'))}." + _ids(opp, g.get("object_id")))
        if state_filter:
            ins = [r for r in recs_all if (r.get("geography") or {}).get("state") == state_filter]
            if ins:
                lines["where"].append(
                    f"Counties in {geo['name']} among the national top {len(recs_all)}: " + "; ".join(
                        f"{(r.get('geography') or {}).get('name')} (national rank "
                        f"{_fmt((r.get('geography') or {}).get('rank'))}, composite "
                        f"{_fmt((r.get('geography') or {}).get('composite'), 3)}, burden evidence level "
                        f"{(r.get('geography') or {}).get('burden_evidence_level')})" for r in ins) + "."
                    + _ids([_opp_id(r) if stored else None for r in ins],
                           [(r.get("geography") or {}).get("object_id") for r in ins]))
            else:
                lines["where"].append(
                    f"No county in {geo['name']} is among the national top {len(recs_all)} of this ranking; a ranking "
                    f"of {geo['name']} counties alone is {UNKNOWN} from these tools." + _ids(geo["object_id"]))
        for u in (rd.get("general_uncertainties") or [])[:3]:
            lines["uncertain"].append(u + _ids(f"measurement:{mid}", (sc.get("condition") or {}).get("object_ids")))
        if recs:
            u0 = recs[0].get("uncertainties") or []
            for u in u0[:4]:
                lines["uncertain"].append(f"{(recs[0].get('geography') or {}).get('name')}: {u}"
                                          + _ids((recs[0].get("geography") or {}).get("object_id")))
    else:
        lines["where"].append(f"Ranked deployment opportunities: {UNKNOWN} - {(rk.get('reason') or '')[:300]}"
                              + _ids([f"condition:{c}" for c in cids], f"measurement:{mid}"))
        unknowns.append(f"Deployment ranking for {cq} x {mid}: {UNKNOWN}.")

    # ------------------------------------------------------------------ 9. named geography / top region: clinics
    target_geo = None
    if geo:
        target_geo = geo["geo_id"]
        gc = log.call("9 geography", "get_geographic_context", geography_id=target_geo)
        if gc["status"] == T.OK:
            gd = gc["data"]
            pop = (gd.get("population") or {}).get("total_population")
            vul = (gd.get("vulnerability") or {}).get("svi_overall")
            lines["where"].append(f"Named geography {geo['name']}: population {_fmt(pop)} (ACS), SVI overall "
                                  f"{_fmt(vul)}." + _ids(geo["object_id"]))
        opp_id = f"opportunity:{_set_id(cq)}|{mid}|{target_geo}"
        tr = log.call("9 geography", "trace_evidence", object_id=opp_id)
        if tr["status"] == T.OK:
            comp = tr["data"].get("components") or {}
            bc = comp.get("burden") or {}
            row_level = ((tr["data"].get("rows") or [{}])[0] or {}).get("geo_level") or geo.get("geo_level")
            lines["where"].append(
                f"{geo['name']} in the stored {row_level}-level ranking: rank {_fmt(tr['data'].get('rank_equal'))} "
                f"among {'states' if row_level == 'state' else 'counties'} (equal weights), "
                f"composite {_fmt(tr['data'].get('composite_equal'), 3)}; burden percentile "
                f"{_fmt(bc.get('percentile'))} (evidence level {bc.get('evidence_level')}"
                + ("; includes a state value inherited by the county, not county prevalence" if bc.get("inherited")
                   else "")
                + ("; Long COVID burden is a modelled small-area estimate, not observed county prevalence"
                   if bc.get("modelled_small_area") else "") + ")."
                + (f" {tr['data']['composite_rank_note']}" if tr["data"].get("composite_rank_note") else "")
                + _ids(opp_id))
        else:
            lines["where"].append(f"Stored opportunity row for {geo['name']}: {UNKNOWN} - "
                                  f"{(tr.get('reason') or '')[:200]}" + _ids(geo["object_id"]))
    elif recs:
        target_geo = (recs[0].get("geography") or {}).get("fips")
    if target_geo:
        gname = geo["name"] if geo else (recs[0].get("geography") or {}).get("name")
        for c in conditions:
            fc = log.call("9 clinics", "find_candidate_clinics", geography_id=target_geo, condition=c["condition_id"],
                          measurement=mid, max_sites=6)
            fd = fc.get("data") or {}
            if fc["status"] == T.OK:
                # the facility pool of a state is the whole state (find_candidate_clinics ignores the radius there)
                state_pool = (fd.get("geography") or {}).get("geo_level") == "state"
                pool = (f"In {gname} for {c['preferred_name']}: {_fmt(fd.get('n_facilities_in_pool'))} geocoded "
                        f"facilities inside the state" if state_pool else
                        f"Around {gname} for {c['preferred_name']}: {_fmt(fd.get('n_facilities_in_pool'))} geocoded "
                        f"facilities within {_fmt(fd.get('radius_km'))} km")
                lines["who"].append(
                    f"{pool}, {_fmt(fd.get('n_eligible'))} with a relevance "
                    "signal. Facilities with characteristics suggesting they may be viable implementation or study "
                    "partners (one per characteristic, round-robin; not a quality ranking):"
                    + _ids(f"geo:{target_geo}", f"condition:{c['condition_id']}"))
                dist_note = " from the state's internal point" if state_pool else ""
                for s in (fd.get("candidates") or [])[:3]:
                    lines["who"].append(f"  - {s.get('facility_name')} ({s.get('city')}, {s.get('state')}; "
                                        f"{_fmt(s.get('distance_km'))} km{dist_note}; selected via "
                                        f"{s.get('selected_via')}): " + (s.get("reasons") or [""])[0]
                                        + _ids(s.get("object_id")))
                if fd.get("absence_note"):
                    lines["who"].append(fd["absence_note"] + _ids(f"geo:{target_geo}"))
            else:
                lines["who"].append(f"Candidate facilities around {gname} for {c['preferred_name']}: {UNKNOWN} - "
                                    f"{(fc.get('reason') or '')[:200]}" + _ids(f"geo:{target_geo}"))
                unknowns.append(f"Candidate facilities ({target_geo}, {c['condition_id']}): {UNKNOWN}.")

    # ------------------------------------------------------------------ 10. research centres + trials
    trial_count_explained = False
    for c in conditions:
        rc = log.call("10 research", "find_relevant_research_centers", condition=c["condition_id"], measurement=mid,
                      top_n=3)
        if rc["status"] == T.OK:
            cs = rc["data"].get("centers") or []
            lines["who"].append(
                f"Research readiness for {c['preferred_name']} nationally ({_fmt(rc['data'].get('n_facilities'))} "
                f"facilities with a trial site or NIH award; display order, not a ranking): " + "; ".join(
                    f"{x.get('facility_name')} ({x.get('city')}, {x.get('state')}: {_fmt(x.get('n_trials'))} trials, "
                    f"{_fmt(x.get('n_nih_core_projects'))} NIH core projects, "
                    f"{_fmt(x.get('n_trials_mentioning_measurement'))} trials mentioning the measurement)" for x in cs)
                + "." + _ids([x.get("object_id") for x in cs]))
        else:
            unknowns.append(f"Research centres for {c['condition_id']}: {UNKNOWN} ({(rc.get('reason') or '')[:160]}).")
        tr = log.call("10 trials", "find_relevant_trials", condition=c["condition_id"], measurement=mid, max_trials=3)
        if tr["status"] == T.OK:
            sm = tr["data"].get("summary") or {}
            ts = tr["data"].get("trials") or []
            how = ("a plain text search of titles, summaries, keywords, interventions and outcomes; it differs from "
                   "the 'describe objective use' count in section 1, which comes from the audited measurement mining "
                   "without questionnaire-only or non-human-context mentions" if not trial_count_explained
                   else "plain text search, as above")
            trial_count_explained = True
            lines["who"].append(
                f"{_fmt(sm.get('n_trials'))} registered {c['preferred_name']} trial(s) (literal condition match) whose "
                f"registered text mentions {mid} ({how}); {_fmt(sm.get('n_recruiting'))} recruiting. Registration is a "
                "research-activity signal, not evidence of effectiveness." + _ids([t.get("object_id") for t in ts]))
        else:
            lines["who"].append(f"Registered {c['preferred_name']} trials mentioning {mid}: {UNKNOWN} - "
                                f"{(tr.get('reason') or '')[:200]}" + _ids(f"condition:{c['condition_id']}"))

    # ------------------------------------------------------------------ 11. provenance chain of the top opportunity
    top_opp = None
    if recs and stored:
        top_opp = _opp_id(recs[0])
    if top_opp:
        tr = log.call("11 trace", "trace_evidence", object_id=top_opp, max_rows=1)
        if tr["status"] == T.OK:
            td = tr["data"]
            burden_links = [lk for lk in td.get("lineage") or [] if lk.get("table") == "geo_condition_burden"]
            for lk in burden_links[:3]:
                row = lk.get("row") or {}
                lines["prov"].append(
                    f"{top_opp} <- burden row {row.get('burden_row_id')} ({row.get('measure_label')}, value "
                    f"{_fmt(row.get('value'))} {row.get('value_unit') or ''}, evidence level "
                    f"{row.get('burden_evidence_level')}, source resolution {row.get('source_geographic_resolution')}) "
                    f"from {row.get('source_name')}." + _ids(top_opp))
            for s in (td.get("sources") or [])[:4]:
                raw = (s.get("raw") or {})
                files = [f.get("path") for f in raw.get("record_files") or []][:2]
                lines["prov"].append(
                    f"Source {s.get('source_id')}: {s.get('name')} ({s.get('source_version')}; retrieved "
                    f"{s.get('retrieved_at')}); audit {', '.join(a['path'] for a in s.get('data_audit') or [])}; raw "
                    + (", ".join(files) if files else f"see {raw.get('manifest')}") + "." + _ids(top_opp))
        else:
            lines["prov"].append(f"Trace of {top_opp}: {UNKNOWN} ({tr.get('reason')})")
    return _finish(q, parsed, log, lines, unknowns)


WEARABLE_SIGNATURE_DATASETS = ("nhanes", "stanford_")   # datasets whose signature features are wearable features


def _feature_kind(dataset_id) -> str:
    """Name the features of a signature by modality: omics signatures (GEO cohorts, mapMECFS) are never 'wearable'."""
    return ("wearable features" if str(dataset_id or "").startswith(WEARABLE_SIGNATURE_DATASETS)
            else "molecular (not wearable) features")


LITERATURE_SOURCES = frozenset({"ot_literature"})   # Europe PMC co-mention (omics.coherence.SOURCE_FAMILY)


def system_support_tiers(supported_systems) -> dict[str, list[str]]:
    """Group a condition's enriched physiological systems by what supports them (Test 3; results/MOLECULAR_COHERENCE.md):
    'non_literature_unflagged' (genetic, mapMECFS or Open Targets 'other' support without a post-hoc flag),
    'drug_target_flagged' (non-literature support only from post-hoc-flagged sources: mostly drug targets of
    trialled drugs, or one GWAS locus) and 'literature_only' (Europe PMC co-mention only)."""
    out = {"non_literature_unflagged": [], "drug_target_flagged": [], "literature_only": []}
    if not isinstance(supported_systems, list):
        return out
    by_sys: dict[str, list[dict]] = {}
    for x in supported_systems:
        if isinstance(x, dict) and x.get("physiological_system"):
            by_sys.setdefault(str(x["physiological_system"]), []).append(x)
    for sysid, rows in sorted(by_sys.items()):
        nonlit = [x for x in rows if x.get("source") not in LITERATURE_SOURCES]
        if any(not x.get("post_hoc_flagged") for x in nonlit):
            out["non_literature_unflagged"].append(sysid)
        elif nonlit:
            out["drug_target_flagged"].append(sysid)
        else:
            out["literature_only"].append(sysid)
    return out


def describe_system_tiers(tiers: dict[str, list[str]]) -> str:
    if not any(tiers.values()):
        return f"enriched physiological systems: {UNKNOWN}."
    parts = ["systems with non-literature support after removing post-hoc-flagged sources: "
             + (", ".join(tiers["non_literature_unflagged"]) or "none")]
    if tiers["drug_target_flagged"]:
        parts.append("carried only by post-hoc-flagged sources (mostly drug targets of trialled drugs, or one GWAS "
                     "locus; not disease biology): " + ", ".join(tiers["drug_target_flagged"]))
    if tiers["literature_only"]:
        parts.append("literature co-mention only (study attention): " + ", ".join(tiers["literature_only"]))
    return "; ".join(parts) + "."


def _opp_id(rec: dict) -> str | None:
    """The recommendation's own stored opportunity row id (None for an on-the-fly ranking: no stored row; a member
    condition's stored row in provenance.object_ids is a component id, not this region's row)."""
    row = (rec.get("provenance") or {}).get("opportunity_row") or {}
    return row.get("object_id") if row.get("stored") else None


def _set_id(cq: str) -> str:
    """The condition id used in opportunity object ids ('A or B' ad hoc sets are stored as 'a+b')."""
    if " or " in cq:
        return "+".join(sorted(p.strip() for p in cq.split(" or ")))
    return cq


def _bundle_for_phrase(text: str) -> str | None:
    """Bundle whose id, label or alias equals the phrase (token-wise), else None."""
    toks = _tokens(text)
    for bid, b in load_config("relevance").get("measurement_bundles", {}).items():
        for t in [bid.replace("_", " "), re.sub(r"\([^)]*\)", " ", b.get("label", ""))] + list(b.get("aliases", [])):
            if _tokens(t) == toks:
                return bid
    return None


def _members(mid: str) -> list[str]:
    b = load_config("relevance").get("measurement_bundles", {}).get(mid)
    return list(b["members"]) if b else [mid]


def _mes_counts(x: dict) -> str:
    mes = (x.get("dimensions") or {}).get("measurement_evidence_strength") or {}
    return (f"{_fmt(mes.get('n_trials_objective'))} trials, {_fmt(mes.get('n_nih_core_projects'))} NIH core "
            f"projects; tier {mes.get('tier')}")


def _level_word(lvl) -> str:
    return {"A": "direct condition measure", "B": "closely matching coded condition", "C": "symptom/comorbidity proxy",
            "D": "no usable burden estimate"}.get(str(lvl), "unknown")


def _finish(q: str, parsed: dict, log: _Log, lines: dict, unknowns: list[str], stopped: str | None = None) -> dict:
    md = [f"# {q}", "", "_Deterministic answer (no language model): every sentence is filled from tool output and ends "
          "with the object ids it rests on; `trace_evidence(<id>)` resolves any of them. UNKNOWN / NOT AVAILABLE means "
          "a tool had no data._", "", "## How the question was read"]
    md += [f"- {x}" for x in lines["parse"]]
    if lines.get("scope"):
        md += ["", "## Scope (what this engine does not do)"]
        md += [f"- {x}" for x in lines["scope"]]
    if stopped:
        md += ["", f"**Stopped: {stopped}.** Nothing is ranked or attributed without a resolved input."]
    sections = [("## 1. What to measure", "what"), ("## 2. Where to measure it", "where"),
                ("## 3. Who could realistically deploy it", "who"), ("## Uncertainties", "uncertain"),
                ("## Provenance (top opportunity traced to sources and raw files)", "prov")]
    for title, key in sections:
        if lines[key]:
            md += ["", title]
            md += [(x if x.startswith("  ") else f"- {x}") for x in lines[key]]
    if unknowns:
        md += ["", "## UNKNOWN / NOT AVAILABLE"]
        md += [f"- {u}" for u in dict.fromkeys(unknowns)]
    md += ["", "## Tool calls (in order)", "", "| # | step | tool | status | s |", "|---|---|---|---|---|"]
    for i, c in enumerate(log.calls, 1):
        md.append(f"| {i} | {c['step']} | {c['tool']}({json.dumps(c['args'])[1:-1][:90]}) | {c['status']} | "
                  f"{c['elapsed_s']} |")
    text = "\n".join(md)
    cited = []
    for group in re.findall(r"\[([^\[\]]+)\]\s*$", text, flags=re.M):     # the citation brackets at line ends
        for oid in group.split("; "):
            if T.OBJECT_ID_FULL_RE.match(oid.strip()) and oid.strip() not in cited:
                cited.append(oid.strip())
    return {"question": q, "parsed": parsed, "stopped": stopped, "sections": lines, "unknowns": list(dict.fromkeys(unknowns)),
            "tool_calls": log.calls, "cited_object_ids": cited, "answer_markdown": text}


SPEC_DEMO_QUESTION = ("Where in the United States would deploying wearable-based autonomic or activity monitoring be "
                      "most valuable for detecting objective abnormalities associated with Long COVID, ME/CFS, "
                      "dysautonomia/POTS-like phenotypes, or related invisible illnesses, and which clinics or research "
                      "sites are best positioned to implement the technology?")
CAPILLAROSCOPY_QUESTION = ("Given phenotype orthostatic intolerance and candidate measurement technology nailfold "
                           "capillaroscopy, where should we deploy it?")
