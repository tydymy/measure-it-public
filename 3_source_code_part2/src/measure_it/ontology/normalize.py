"""Condition search and normalization over the ontology layer.

    search_condition(text)            -> ranked candidate conditions with ids and match reasons
    normalize_condition(text_or_code) -> canonical condition(s) for free text, an alias,
                                         an ICD-10-CM code, or a MONDO / EFO / MeSH / HPO id

Both read the tables written by `measure_it.ontology.build_registry`
(condition_registry, condition_ontology_mappings, condition_phenotypes,
phenotype_axes, ontology_icd10cm_codes). Nothing here invents an identifier:
when no registry row supports a match the result is
``{"status": UNKNOWN, ...}`` with the reason.

The pure text helpers (norm_text, compact_key, singular_key, is_acronym) are
also used by build_registry so matching at build time and at query time is
identical.
"""
from __future__ import annotations

import re
import unicodedata
from dataclasses import dataclass, field
from functools import lru_cache

import pandas as pd
from rapidfuzz import fuzz

from ..config import UNKNOWN

# --------------------------------------------------------------------------- text helpers


def norm_text(text: str) -> str:
    """Lower-case ASCII, punctuation (hyphens, slashes, brackets) -> single spaces."""
    t = unicodedata.normalize("NFKD", str(text)).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]+", " ", t).strip()


def compact_key(text: str) -> str:
    """norm_text without spaces: 'Post-exertional malaise' == 'Postexertional malaise'."""
    return norm_text(text).replace(" ", "")


def _singular_word(w: str) -> str:
    return w[:-1] if len(w) > 3 and w.endswith("s") and not w.endswith("ss") else w


def singular_key(text: str) -> str:
    """Plural-insensitive key, applied identically to both sides of a comparison."""
    return " ".join(_singular_word(w) for w in norm_text(text).split())


def is_acronym(text: str) -> bool:
    """Short token with >=2 capitals and no spaces (PASC, CFS, hEDS, ME/CFS, POTS)."""
    s = str(text).strip()
    return " " not in s and 0 < len(s) <= 7 and sum(ch.isupper() for ch in s) >= 2


def strip_parenthetical(text: str) -> str:
    return re.sub(r"\s*\([^)]*\)\s*", " ", str(text)).strip()


# --------------------------------------------------------------------------- id helpers

ICD_RE = re.compile(r"^(?:ICD-?10-?CM[:\s]*)?([A-Z][0-9][0-9A-Z])(?:\.?([0-9A-Z]{1,4}))?$", re.I)
MESH_BARE_RE = re.compile(r"^[CD]\d{6}(?:\d{3})?$")
CURIE_RE = re.compile(r"^([A-Za-z][A-Za-z0-9.]*)[:_](\S+)$")

# Canonical CURIE prefixes used in every table this module writes.
PREFIX_CANON = {
    "MONDO": "MONDO", "EFO": "EFO", "HP": "HP", "HPO": "HP", "DOID": "DOID", "MESH": "MESH",
    "ORPHANET": "ORPHANET", "ORPHA": "ORPHANET", "ORDO": "ORPHANET", "UMLS": "UMLS", "NCIT": "NCIT",
    "SCTID": "SNOMEDCT", "SNOMEDCT": "SNOMEDCT", "SNOMED": "SNOMEDCT", "OMIM": "OMIM", "OMIMPS": "OMIMPS",
    "MEDGEN": "MEDGEN", "GARD": "GARD", "ICD10CM": "ICD10CM", "ICD9": "ICD9", "ICD9CM": "ICD9CM",
    "ICD10WHO": "ICD10WHO", "ICD11.FOUNDATION": "ICD11", "ICD11": "ICD11", "MEDDRA": "MEDDRA",
    "NORD": "NORD", "NANDO": "NANDO", "ICDO": "ICDO", "ONCOTREE": "ONCOTREE", "OMIA": "OMIA",
    "VENOM": "VENOM", "WIKIPEDIA": "WIKIPEDIA", "DECIPHER": "DECIPHER", "GTR": "GTR", "HGNC": "HGNC",
    "CSP": "CSP", "PMID": "PMID", "RO": "RO",
}


def canonical_curie(curie: str) -> str:
    """'mesh:D015673' / 'Orphanet:1983' / 'MONDO_0005404' / 'icd10cm:Q79.6' -> canonical CURIE."""
    s = str(curie).strip()
    m = CURIE_RE.match(s)
    if not m:
        return s
    prefix, local = m.group(1), m.group(2)
    canon = PREFIX_CANON.get(prefix.upper(), prefix.upper())
    if canon == "ICD10CM":
        local = icd_dotted(local) or local.upper()
    return f"{canon}:{local}"


def icd_nodot(code: str) -> str:
    return str(code).replace(".", "").strip().upper()


def icd_dotted(code: str) -> str | None:
    """'G9332' / 'g93.32' -> 'G93.32'; ranges such as 'B90-B94' are returned unchanged."""
    s = str(code).strip().upper()
    if re.fullmatch(r"[A-Z][0-9][0-9A-Z]-[A-Z][0-9][0-9A-Z]", s):
        return s
    m = ICD_RE.match(s)
    if not m:
        return None
    return m.group(1).upper() + (("." + m.group(2).upper()) if m.group(2) else "")


# --------------------------------------------------------------------------- index

PREDICATE_RANK = {"exact": 0, "close": 1, "narrow": 2, "broad": 3, "related": 4}
PRED_CONF = {"exact": 1.0, "close": 0.95, "narrow": 0.85, "broad": 0.7, "related": 0.6}
ROLE_CONF = {"primary": 1.0, "component": 0.95, "narrower": 0.85, "broader_context": 0.6, "related_distinct": 0.5}
REGISTRY_PREDICATES = {"exact", "close", "narrow"}

# How much a text match on each lexicon kind is worth (ties broken by this).
KIND_WEIGHT = {
    "preferred_name": 1.0, "alias": 0.98, "search_term": 0.97,
    "mondo_label": 0.96, "mondo_exact_synonym": 0.95, "efo_label": 0.94, "mesh_label": 0.94,
}


@dataclass
class LexEntry:
    condition_id: str
    text: str
    kind: str
    weight: float
    compact: str
    singular: str
    tokens: frozenset
    acronym: bool


@dataclass
class OntologyIndex:
    conditions: dict = field(default_factory=dict)          # cid -> summary dict
    lexicon: list = field(default_factory=list)             # LexEntry
    ids: dict = field(default_factory=dict)                 # canonical CURIE -> [mapping dict]
    icd: list = field(default_factory=list)                 # (code_nodot, dotted, cid, predicate, source)
    icd_blocks: list = field(default_factory=list)          # (start, end, cid, predicate, source)
    hpo: dict = field(default_factory=dict)                 # HP id -> [(cid, scope, label)]
    axes: dict = field(default_factory=dict)                # HP id -> [axis dict]
    axis_lexicon: dict = field(default_factory=dict)        # compact label -> axis dict
    cms_codes: dict = field(default_factory=dict)           # code_nodot -> (dotted, description, billable)

    @classmethod
    def from_tables(cls, registry: pd.DataFrame, mappings: pd.DataFrame,
                    phenotypes: pd.DataFrame | None = None, axes: pd.DataFrame | None = None,
                    cms_codes: pd.DataFrame | None = None) -> "OntologyIndex":
        idx = cls()

        def add_lex(cid, text, kind, factor=1.0, display=None):
            if not isinstance(text, str) or not text.strip():
                return
            nt = norm_text(text)
            if not nt:
                return
            idx.lexicon.append(LexEntry(
                cid, display or text, kind, KIND_WEIGHT.get(kind, 0.9) * factor, compact_key(text), singular_key(text),
                frozenset(singular_key(text).split()), is_acronym(text)))

        for r in registry.to_dict("records"):
            cid = r["canonical_condition_id"]
            idx.conditions[cid] = {
                "canonical_condition_id": cid,
                "preferred_name": r["preferred_name"],
                "primary_mondo_id": r.get("primary_mondo_id"),
                "mondo_ids": _as_list(r.get("mondo_ids")),
                "icd10cm_codes": _as_list(r.get("icd10cm_codes")),
                "efo_ids": _as_list(r.get("efo_ids")),
                "mesh_ids": _as_list(r.get("mesh_ids")),
                "grouping_only": bool(r.get("grouping_only", False)),
            }
            add_lex(cid, r["preferred_name"], "preferred_name")
            add_lex(cid, strip_parenthetical(r["preferred_name"]), "preferred_name")
            for inner in re.findall(r"\(([^)]*)\)", r["preferred_name"]):
                if len(inner.split()) >= 2:
                    add_lex(cid, inner, "preferred_name", display=r["preferred_name"])
            for a in _as_list(r.get("aliases")):
                add_lex(cid, a, "alias")
                if "(" in a:
                    add_lex(cid, strip_parenthetical(a), "alias", display=a)
            for s in _as_list(r.get("search_terms")):
                add_lex(cid, s, "search_term")
            # Mondo labels/synonyms of the classes in the registry list; down-weighted when the
            # primary class is curated as broader than the condition (e.g. POTS -> orthostatic intolerance).
            factor = 1.0 if r.get("primary_mondo_relation", "exact") in ("exact", "close", "narrow") else 0.9
            for s in _as_list(r.get("mondo_labels")):
                add_lex(cid, s, "mondo_label", factor)
            for s in _as_list(r.get("mondo_exact_synonyms")):
                add_lex(cid, s, "mondo_exact_synonym", factor)

        for m in mappings.to_dict("records"):
            tid = canonical_curie(m["target_id"])
            idx.ids.setdefault(tid, []).append(m)
            # Short-form alias (MONDO_0005404 / EFO_0000555) resolves to the same rows.
            pred = m.get("predicate_condition")
            cid = m["canonical_condition_id"]
            if m["target_ontology"] == "ICD10CM" and m.get("in_registry_list"):
                if m.get("code_kind") == "block":
                    a, b = str(m["target_id"]).split(":", 1)[1].split("-")
                    idx.icd_blocks.append((a, b, cid, pred, m.get("mapping_source")))
                else:
                    code = str(m["target_id"]).split(":", 1)[1]
                    idx.icd.append((icd_nodot(code), code, cid, pred, m.get("mapping_source")))
            if m["target_ontology"] in ("EFO", "MESH") and m.get("in_registry_list") and isinstance(m.get("target_label"), str):
                # a narrower id's label (e.g. MeSH 'Neurocirculatory Asthenia' for dysautonomia) is a subtype name:
                # down-weighted so an exact text hit on it is reported as partial, not matched.
                add_lex(cid, m["target_label"], "efo_label" if m["target_ontology"] == "EFO" else "mesh_label",
                        1.0 if pred in ("exact", "close") else 0.9)

        if phenotypes is not None and len(phenotypes):
            for p in phenotypes[["canonical_condition_id", "hpo_id", "hpo_label", "association_scope", "negated"]].to_dict("records"):
                if p.get("negated") is True:
                    continue
                idx.hpo.setdefault(p["hpo_id"], []).append((p["canonical_condition_id"], p["association_scope"], p["hpo_label"]))
        if axes is not None and len(axes):
            for a in axes.to_dict("records"):
                if isinstance(a.get("hpo_id"), str):
                    idx.axes.setdefault(a["hpo_id"], []).append(a)
                for t in (a.get("axis_label"), a.get("hpo_label")):
                    if isinstance(t, str):
                        idx.axis_lexicon[compact_key(t)] = a
        if cms_codes is not None and len(cms_codes):
            if "in_fy2026" in cms_codes.columns:
                cms_codes = cms_codes[cms_codes["in_fy2026"].astype(bool)]
            for c in cms_codes[["code_nodot", "code", "description", "billable"]].itertuples(index=False):
                idx.cms_codes[c.code_nodot] = (c.code, c.description, bool(c.billable))
        return idx


def _as_list(v) -> list:
    if v is None:
        return []
    if isinstance(v, float) and pd.isna(v):
        return []
    if isinstance(v, str):
        return [v]
    return [x for x in list(v) if x is not None]


@lru_cache(maxsize=1)
def get_index() -> OntologyIndex:
    """Load the processed ontology tables (cached for the process)."""
    from ..store import read_table, table_exists

    def opt(name):
        return read_table(name) if table_exists(name) else None

    return OntologyIndex.from_tables(
        read_table("condition_registry"), read_table("condition_ontology_mappings"),
        opt("condition_phenotypes"), opt("phenotype_axes"), opt("ontology_icd10cm_codes"))


# --------------------------------------------------------------------------- search


def _summary(idx: OntologyIndex, cid: str) -> dict:
    c = idx.conditions.get(cid, {})
    return {k: c.get(k) for k in ("canonical_condition_id", "preferred_name", "primary_mondo_id",
                                  "mondo_ids", "icd10cm_codes", "efo_ids", "mesh_ids")}


def search_condition(text: str, *, limit: int = 10, min_score: float = 0.6,
                     index: OntologyIndex | None = None) -> list[dict]:
    """Rank registry conditions against free text. Each candidate carries its ids and why it matched.

    Scores: 1.0 x kind weight for an exact (punctuation/plural-insensitive) match;
    0.8 when all query words appear in a registry term covering >=50% of it (or vice versa);
    0.75 x similarity for a close spelling (token_sort_ratio >= 88). Acronym queries
    (e.g. 'POTS', 'CFS') only match exactly, never fuzzily.
    """
    idx = index or get_index()
    q = str(text or "").strip()
    if not q:
        return []
    qc, qs = compact_key(q), singular_key(q)
    qn = norm_text(q)
    qtok = frozenset(qs.split())
    q_acr = is_acronym(q)
    best: dict[str, dict] = {}
    for e in idx.lexicon:
        score, reason = 0.0, ""
        if e.compact == qc or (not q_acr and not e.acronym and e.singular == qs):
            score, reason = 1.0 * e.weight, f"exact match to {e.kind} '{e.text}'"
        elif q_acr or e.acronym or not qtok:
            continue
        elif qtok <= e.tokens and len(qtok) / max(len(e.tokens), 1) >= 0.5:
            score, reason = 0.8 * e.weight, f"all query words appear in {e.kind} '{e.text}'"
        elif e.tokens <= qtok and len(e.tokens) / max(len(qtok), 1) >= 0.5 and max(len(w) for w in e.tokens) >= 4:
            score, reason = 0.8 * e.weight, f"query contains {e.kind} '{e.text}'"
        else:
            ratio = fuzz.token_sort_ratio(qn, norm_text(e.text))
            if ratio >= 88:
                score, reason = 0.75 * e.weight * ratio / 100, f"close spelling ({ratio:.0f}/100) of {e.kind} '{e.text}'"
        if score >= min_score and score > best.get(e.condition_id, {}).get("score", 0):
            best[e.condition_id] = {**_summary(idx, e.condition_id), "score": round(score, 4),
                                    "match_reason": reason, "matched_text": e.text, "matched_kind": e.kind}
    return sorted(best.values(), key=lambda d: (-d["score"], d["canonical_condition_id"]))[:limit]


# --------------------------------------------------------------------------- normalize


def classify_input(value: str, index: OntologyIndex | None = None) -> tuple[str, str]:
    """Return (input_type, canonical value). input_type in icd10cm/mondo/efo/mesh/hpo/curie/text."""
    s = str(value).strip()
    m = CURIE_RE.match(s)
    if m and not ICD_RE.match(s):
        prefix = PREFIX_CANON.get(m.group(1).upper())
        if prefix:
            cur = canonical_curie(s)
            kind = {"MONDO": "mondo", "EFO": "efo", "MESH": "mesh", "HP": "hpo", "ICD10CM": "icd10cm"}.get(prefix, "curie")
            return kind, (cur.split(":", 1)[1] if kind == "icd10cm" else cur)
    if s.upper().startswith("ICD10CM:") or s.upper().startswith("ICD-10-CM"):
        d = icd_dotted(s.split(":", 1)[-1] if ":" in s else s)
        if d:
            return "icd10cm", d
    if MESH_BARE_RE.match(s.upper()):
        idx = index
        if not (idx and icd_nodot(s) in idx.cms_codes):
            return "mesh", f"MESH:{s.upper()}"
    if ICD_RE.match(s):
        return "icd10cm", icd_dotted(s)
    return "text", s


def _unknown(query: str, input_type: str, reason: str, **extra) -> dict:
    return {"status": UNKNOWN, "query": query, "input_type": input_type, "matches": [], "reason": reason, **extra}


def _finish(query: str, input_type: str, matches: list[dict], exact_threshold: float = 0.9, **extra) -> dict:
    matches = sorted(matches, key=lambda d: (-d["confidence"], PREDICATE_RANK.get(d.get("predicate") or "related", 9),
                                             d["canonical_condition_id"]))
    # one row per condition, best first
    seen, uniq = set(), []
    for mt in matches:
        if mt["canonical_condition_id"] not in seen:
            seen.add(mt["canonical_condition_id"])
            uniq.append(mt)
    top = uniq[0]["confidence"]
    best = [m for m in uniq if m["confidence"] == top]
    if top >= exact_threshold:
        status = "matched" if len(best) == 1 else "ambiguous"
    else:
        status = "partial"
    return {"status": status, "query": query, "input_type": input_type, "matches": best,
            "other_candidates": [m for m in uniq if m not in best], **extra}


def normalize_condition(text_or_code: str, *, index: OntologyIndex | None = None) -> dict:
    """Map free text, an alias, an ICD-10-CM code, or a MONDO/EFO/MeSH/HPO id to registry condition(s).

    Returns {"status": "matched" | "ambiguous" | "partial" | UNKNOWN, "query", "input_type",
    "matches": [{canonical_condition_id, preferred_name, match_reason, matched_id, predicate,
    confidence, ...ids}], "other_candidates": [...], ...}. HPO ids return the conditions whose
    Monarch disease-phenotype annotations include the term (a phenotype link, not an equivalence).
    """
    idx = index or get_index()
    q = str(text_or_code or "").strip()
    if not q:
        return _unknown(q, "empty", "empty input")
    kind, val = classify_input(q, idx)

    if kind == "icd10cm":
        return _normalize_icd(idx, q, val)
    if kind == "hpo":
        return _normalize_hpo(idx, q, val)
    if kind in ("mondo", "efo", "mesh", "curie"):
        rows = idx.ids.get(val, [])
        if not rows:
            return _unknown(q, kind, f"{val} is not mapped to any registry condition")
        matches = []
        for r in rows:
            pred = r.get("predicate_condition")
            role = r.get("mondo_role")
            conf = PRED_CONF.get(pred, 0.5)
            if str(r.get("target_id", "")).startswith("MONDO:") and isinstance(role, str) and role in ROLE_CONF:
                # the Mondo class itself: its role caps the confidence, and so does the curated relation
                # (e.g. MONDO:0001315 is the primary class of POTS but curated as BROADER than POTS).
                conf = min(ROLE_CONF[role], conf)
            note = ""
            if r.get("target_status") and r["target_status"] != "current":
                conf = min(conf, 0.9)
                note = f"; target status: {r['target_status']}"
            matches.append({**_summary(idx, r["canonical_condition_id"]), "matched_id": val, "predicate": pred,
                            "confidence": conf,
                            "match_reason": f"{val} mapped via {r.get('mapping_source')} "
                                            f"(role {r.get('mondo_role') or 'n/a'}, predicate {pred}){note}"})
        return _finish(q, kind, matches)

    # free text / alias
    cands = search_condition(val, index=idx, min_score=0.6)
    axis = idx.axis_lexicon.get(compact_key(val))
    extra = {"phenotype_axis": {k: axis.get(k) for k in ("axis_id", "hpo_id", "hpo_label")}} if axis else {}
    good = [c for c in cands if c["score"] >= 0.7]
    if not good:
        reason = "no registry condition matches this text"
        if axis:
            reason += f"; it matches phenotype axis '{axis['axis_id']}' ({axis.get('hpo_id')}), which is a phenotype, not a condition"
        return _unknown(q, "text", reason, candidates_below_threshold=cands[:3], **extra)
    matches = [{**c, "matched_id": None, "predicate": None, "confidence": c["score"]} for c in good]
    return _finish(q, "text", matches, **extra)


def _normalize_icd(idx: OntologyIndex, q: str, code: str) -> dict:
    nd = icd_nodot(code)
    cms = idx.cms_codes.get(nd)
    extra = {"icd10cm_code": code,
             "cms_fy2026": ({"code": cms[0], "description": cms[1], "billable": cms[2]} if cms else
                            (UNKNOWN if not idx.cms_codes else "not a code in the CMS FY2026 code set"))}
    matches = []
    for mnd, mcode, cid, pred, src in idx.icd:
        if nd == mnd:
            conf, why = 1.0, f"exact ICD-10-CM code {mcode}"
        elif nd.startswith(mnd):
            conf, why = 0.92, f"{code} falls under mapped ICD-10-CM category {mcode}"
        elif mnd.startswith(nd):
            conf, why = 0.6, f"{code} is a parent category of mapped code {mcode}"
        else:
            continue
        matches.append({**_summary(idx, cid), "matched_id": f"ICD10CM:{mcode}", "predicate": pred,
                        "confidence": conf - 0.01 * PREDICATE_RANK.get(pred, 4),
                        "match_reason": f"{why} ({src}; predicate {pred})"})
    cat = nd[:3]
    for a, b, cid, pred, src in idx.icd_blocks:
        if a <= cat <= b:
            matches.append({**_summary(idx, cid), "matched_id": f"ICD10CM:{a}-{b}", "predicate": pred,
                            "confidence": 0.5, "match_reason": f"{code} is inside ICD-10-CM block {a}-{b} ({src}; predicate {pred})"})
    if not matches:
        if cms:
            reason = f"ICD-10-CM {cms[0]} ('{cms[1]}') is a valid FY2026 code but is not mapped to any registry condition"
        elif idx.cms_codes:
            reason = f"{code} is not a code in the CMS FY2026 ICD-10-CM code set"
        else:
            reason = f"{code} is not mapped to any registry condition"
        return _unknown(q, "icd10cm", reason, **extra)
    return _finish(q, "icd10cm", matches, exact_threshold=0.9, **extra)


def _normalize_hpo(idx: OntologyIndex, q: str, hp: str) -> dict:
    axes = idx.axes.get(hp, [])
    extra = {"phenotype_axes": [a.get("axis_id") for a in axes]} if axes else {}
    rows = idx.hpo.get(hp, [])
    if not rows:
        reason = f"{hp} has no Monarch disease-phenotype annotation on any registry condition"
        if axes:
            reason += f"; it is the HPO term of phenotype axis {[a.get('axis_id') for a in axes]}"
        return _unknown(q, "hpo", reason, **extra)
    matches = []
    for cid, scope, label in rows:
        conf = 0.8 if scope == "direct" else 0.5
        matches.append({**_summary(idx, cid), "matched_id": hp, "predicate": "has_phenotype", "confidence": conf,
                        "match_reason": f"Monarch annotates {'the condition class' if scope == 'direct' else 'a descendant class'} "
                                        f"with {hp} '{label}' (phenotype association, not equivalence)"})
    return _finish(q, "hpo", matches, exact_threshold=0.8, **extra)
