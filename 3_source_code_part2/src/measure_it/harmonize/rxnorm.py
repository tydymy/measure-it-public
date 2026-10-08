"""Drug names -> RxNorm ingredients (IN / MIN) through the NLM RxNav REST API (cached by measure_it.http).

    ingredient_for_name("METOPROLOL SUCCINATE") -> {"rxcui": "6918", "name": "metoprolol", "tty": "IN", ...}

Steps for one name (no fuzzy matching: an unmatched name stays unmapped and is counted):
  1. GET /REST/rxcui.json?name=<name>&search=2   (exact, then normalised string match; never approximate)
  2. GET /REST/rxcui/<cui>/properties.json       term type of the hit
  3. IN or MIN -> done; PIN / BN / SCD / SBD / ... -> GET /REST/rxcui/<cui>/related.json?tty=IN MIN and keep the hit
     only when exactly one MIN, or else exactly one IN, is related (a dose-form group such as "acyclovir Topical
     Product" resolves to its ingredient this way).
Combination products written "A; B" (NHANES generic names) are split and each component is mapped separately.
"""
from __future__ import annotations

import re
from functools import lru_cache

from .. import http

BASE = "https://rxnav.nlm.nih.gov/REST"
NOT_DRUG_CODES = {"55555", "77777", "99999", ""}     # NHANES: unknown / refused / don't know
# names that RxNav resolves by string match to one specific ingredient although they denote a class whose member is
# not recorded (NHANES "INSULIN" -> "insulin, regular, human"): left unmapped on purpose
CLASS_NAMES = {"INSULIN"}


def _get(path: str, params: dict | None = None) -> dict:
    return http.get_json(f"{BASE}/{path}", params=params, timeout=60)


@lru_cache(maxsize=None)
def ingredient_for_name(name: str) -> dict:
    """{'status': 'mapped'|'unmapped'|'unavailable', 'rxcui', 'name', 'tty', 'query', 'reason'}."""
    q = re.sub(r"\s+", " ", str(name or "")).strip()
    if q.upper() in NOT_DRUG_CODES or not q:
        return {"status": "unmapped", "query": q, "reason": "not a drug name (NHANES unknown / refused code)"}
    if q.upper() in CLASS_NAMES:
        return {"status": "unmapped", "query": q, "reason": "class name; the specific ingredient is not recorded"}
    try:
        r = _get("rxcui.json", {"name": q, "search": "2"})
        ids = ((r or {}).get("idGroup") or {}).get("rxnormId") or []
        if not ids:
            return {"status": "unmapped", "query": q, "reason": "no exact or normalised RxNorm match"}
        out = []
        for cui in ids:
            props = (_get(f"rxcui/{cui}/properties.json") or {}).get("properties") or {}
            tty = props.get("tty")
            if tty in ("IN", "MIN"):
                out.append({"rxcui": str(cui), "name": props.get("name"), "tty": tty})
                continue
            rel = _get(f"rxcui/{cui}/related.json", {"tty": "IN MIN"}) or {}
            groups = (rel.get("relatedGroup") or {}).get("conceptGroup") or []
            by_tty = {g.get("tty"): g.get("conceptProperties") or [] for g in groups}
            ins, mins = by_tty.get("IN", []), by_tty.get("MIN", [])
            pick = mins if len(mins) == 1 else (ins if len(ins) == 1 else [])
            if pick:
                out.append({"rxcui": str(pick[0]["rxcui"]), "name": pick[0]["name"], "tty": pick[0]["tty"],
                            "via": f"{cui} ({tty})"})
        uniq = {o["rxcui"]: o for o in out}
        if len(uniq) == 1:
            return {"status": "mapped", "query": q, **next(iter(uniq.values()))}
        if not uniq:
            return {"status": "unmapped", "query": q, "reason": "match is not linked to exactly one ingredient"}
        return {"status": "unmapped", "query": q, "reason": f"ambiguous: {len(uniq)} ingredients"}
    except http.OfflineCacheMiss:
        return {"status": "unavailable", "query": q, "reason": "offline and not in the HTTP cache"}
    except Exception as exc:  # noqa: BLE001 - a failed lookup is recorded, never invented
        return {"status": "unavailable", "query": q, "reason": f"RxNav error: {type(exc).__name__}"}


def split_components(generic_name: str) -> list[str]:
    """'HYDROCHLOROTHIAZIDE; LISINOPRIL' -> ['HYDROCHLOROTHIAZIDE', 'LISINOPRIL']."""
    return [p.strip() for p in str(generic_name or "").split(";") if p.strip()]


def map_names(names, echo=None) -> dict[str, dict]:
    """{component name: lookup result} for every distinct component of the given generic names."""
    comps = sorted({c for n in names for c in split_components(n)})
    out = {}
    for i, c in enumerate(comps):
        out[c] = ingredient_for_name(c)
        if echo and (i + 1) % 200 == 0:
            echo(f"  rxnav: {i + 1}/{len(comps)} names")
    return out
