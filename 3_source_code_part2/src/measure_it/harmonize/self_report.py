"""Self-reported medical history and PRO instruments: curated, exact-match mappings (no fuzzy matching).

    map_history("Raynaud's phenomenon") -> {"icd10cm": "I73.0", "confidence": "medium", "label": ..., "status": "mapped"}
    instrument_code("compass31")        -> "COMPASS31"
    crosswalk_for("FSS9:total")         -> [{"reference": "PHQ9:item4_tired", "match_quality": "related_construct", ...}]
"""
from __future__ import annotations

import fnmatch
import re
from functools import lru_cache

from . import schema as S

CONF_CAP = {"high": "medium", "medium": "medium", "low": "low"}


def norm_text(s) -> str:
    s = str(s or "").lower().replace("'", "").replace("’", "")
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


@lru_cache(maxsize=1)
def _history_index() -> dict:
    cfg = S.config("harmonize_self_report_icd10")
    idx = {}
    for r in cfg["mappings"]:
        for a in r["aliases"]:
            k = norm_text(a)
            if k in idx and idx[k]["icd10cm"] != r["icd10cm"]:
                raise ValueError(f"harmonize_self_report_icd10.yaml: alias {a!r} maps to two codes")
            idx[k] = r
    return idx


def map_history(text: str, icd10cm: str | None = None) -> dict:
    """One self-reported condition -> ICD-10-CM (mapping_method self_report_map; confidence capped at medium)."""
    code = str(icd10cm or "").strip().upper()
    if code:
        if "." not in code and len(code) > 3:
            code = code[:3] + "." + code[3:]
        if S.ICD10CM_RE.fullmatch(code):
            return {"status": "mapped", "icd10cm": code, "label": str(text or code), "confidence": "medium",
                    "how": "code supplied with the self-report"}
    r = _history_index().get(norm_text(text))
    if r is None:
        return {"status": "unmapped", "icd10cm": None, "label": str(text or ""), "confidence": None,
                "how": "no exact alias in configs/harmonize_self_report_icd10.yaml"}
    return {"status": "mapped", "icd10cm": r["icd10cm"], "label": r["label"], "confidence": CONF_CAP[r["confidence"]],
            "how": "curated alias"}


@lru_cache(maxsize=1)
def _instruments() -> dict:
    cfg = S.config("harmonize_survey_crosswalk")
    out = {}
    for r in cfg["instruments"]:
        for f in r["file_names"]:
            out[f.lower()] = r
    return out


def instrument_code(file_name: str) -> str:
    """survey_<name>.csv -> instrument code (configs/harmonize_survey_crosswalk.yaml) or NAME upper-cased."""
    r = _instruments().get(str(file_name).lower())
    return r["code"] if r else re.sub(r"[^A-Z0-9_]", "_", str(file_name).upper())


def instrument_label(code: str) -> str:
    for r in _instruments().values():
        if r["code"] == code:
            return r["label"]
    return code


def crosswalk_for(item_code: str) -> list[dict]:
    """Reference-cohort items for one study survey item ('<INSTRUMENT>:<item>'); [] when none is listed. Rows with
    match_quality 'none' are returned too (an explicit gap)."""
    cfg = S.config("harmonize_survey_crosswalk")
    out = []
    for r in cfg["crosswalk"]:
        if fnmatch.fnmatchcase(item_code, r["item"]):
            ref = r.get("reference")
            if ref and ref.endswith(":*"):
                ref = ref[:-1] + item_code.split(":", 1)[1]
            out.append({"item": item_code, "reference": ref, "match_quality": r["match_quality"],
                        "direction": r.get("direction"), "rationale": r["rationale"]})
    return out
