"""Optional language-model step of the mapper (docs/AGENT_CONTRACT.md section 0 rules 1-2).

A model may only PROPOSE how a still-unmapped column maps to a role and a code. `measure_it.agent.mapper` validates
every proposal deterministically before it is used; a failed proposal becomes `unmapped` with the reason.

What leaves the machine (`build_items`): table and column names, inferred type, data-dictionary label and choices,
unit, numeric min / median / max, categorical levels seen in >= 11 rows (smaller ones only as ``<11 other`` without a
count), share missing, REDCap field type / form, the detected code pattern. Never row-level values, participant ids,
free text, dates, or any column flagged by the `measure_it.byod.checks` identifier / geography patterns. Every request
payload is appended to ``<out>/llm_requests.jsonl`` BEFORE it is sent, so the user can audit it (no secrets are
logged). Without ``llm="claude"`` nothing here runs and no model is contacted.

Provider interface: any object with ``name``, ``model`` and ``request(payload: dict) -> dict`` returning
``{"proposals": [...]}`` (schema `PROPOSAL_SCHEMA`). `ClaudeProvider` implements it with the official ``anthropic``
SDK (optional extra ``agent``; imported lazily). Claude Fable 5.1 (``claude-fable-5-1``) runs under 30-day data
retention (not zero-retention unless Anthropic authorised it); choose another model with ``model=`` or
``MEASURE_IT_LLM_MODEL``.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Protocol

from ..byod import checks as C

DEFAULT_MODEL = "claude-fable-5-1"       # the most capable widely released Claude model (claude-api skill, 2026-09-25)
MODEL_ENV = "MEASURE_IT_LLM_MODEL"
FALLBACK_BETA = "server-side-fallback-2026-07-01"     # `fallbacks: "default"`: a policy decline is re-run server-side
BATCH_SIZE = 60
ROLES = ["participant_id", "label", "age", "age_band", "sex", "demographic", "self_report_history", "diagnosis_code",
         "lab", "medication_code", "medication_flag", "survey", "device", "omics", "date", "drop", "unmapped"]
LLM_ROLES = [r for r in ROLES if r not in ("participant_id", "label", "age", "age_band", "sex", "date")]
OMICS_LAYERS = ["proteomics", "metabolomics", "lipidomics", "transcriptomics", "epigenomics", "genomics", "glycomics",
                "cytokines", "antibodies", "microbiome", "steroidomics"]
SENSITIVE_FLAGS = {"identifier", "geography", "date", "exact_age", "free_text", "identifier?", "geography?",
                   "id_values_look_identifying"}

_NULLABLE = {"anyOf": [{"type": "string"}, {"type": "null"}]}
PROPOSAL_SCHEMA = {
    "type": "object",
    "properties": {"proposals": {"type": "array", "items": {
        "type": "object",
        "properties": {
            "key": {"type": "string"},
            "role": {"type": "string", "enum": LLM_ROLES},
            "code_system": {"type": "string", "enum": ["loinc", "icd10cm", "rxnorm", "none"]},
            "code": _NULLABLE,
            "instrument": _NULLABLE,
            "item": _NULLABLE,
            "measurement_class": _NULLABLE,
            "omics_layer": _NULLABLE,
            "platform": _NULLABLE,
            "confidence": {"type": "number"},
            "rationale": {"type": "string"},
        },
        "required": ["key", "role", "code_system", "code", "instrument", "item", "measurement_class", "omics_layer",
                     "platform", "confidence", "rationale"],
        "additionalProperties": False}}},
    "required": ["proposals"],
    "additionalProperties": False,
}

SYSTEM_PROMPT = """You help harmonize a research dataset to a fixed schema. You see only column METADATA (names, \
data-dictionary labels, units, value ranges, frequent levels); you never see participant records.

For every column in `columns`, return one proposal keyed by its `key`:
- role: one of the allowed roles. self_report_history = a yes/no item about a medical condition the participant \
reported (give its ICD-10-CM code); lab = a laboratory or clinical-chemistry measurement (give its LOINC code); \
medication_flag = a yes/no item about one drug (give its RxNorm ingredient code); diagnosis_code / medication_code = \
a column whose VALUES are ICD-10-CM / RxNorm codes; survey = a patient-reported questionnaire item or score (give \
`instrument` and `item`); device = a measurement-device feature (give `measurement_class` from the list); omics = an \
omics feature (give `omics_layer` and `platform`); demographic = other demographics; drop = should not be kept; \
unmapped = you are not confident.
- Give a code only if you are certain it is a real, current code of that system for exactly this concept. Never \
guess or construct a code: when unsure, use role "unmapped" (or no code) and say why. Every code is checked against \
local vocabularies and rejected if absent.
- confidence: your probability (0 to 1) that the proposal is correct. rationale: one short sentence.
Answer with JSON matching the schema; one proposal per column."""


class LLMUnavailable(RuntimeError):
    """The provider cannot run here (package missing, no credentials). The message says which."""


class LLMRefused(RuntimeError):
    """The model declined the request (stop_reason 'refusal')."""


class LLMProvider(Protocol):
    name: str
    model: str

    def request(self, payload: dict) -> dict: ...


# ------------------------------------------------------------------------------------------------------------------
# payload (contract rule 2)
# ------------------------------------------------------------------------------------------------------------------
def _clean_text(v) -> str | None:
    """Drop a string that looks identifying (e-mail, phone, date, ZIP, address, coordinates)."""
    if v is None:
        return None
    s = str(v)
    for name, rx in C.VALUE_PATTERNS.items():
        if rx.search(s):
            return None
    for rx in C.ID_PATTERNS.values():
        if rx.search(s):
            return None
    return s


def column_item(key: str, table: dict, col: dict, small_cell: int = 11) -> dict | None:
    """One column's metadata as allowed by contract rule 2, or None when the column must not be sent at all."""
    flags = set(col.get("flags") or [])
    if flags & SENSITIVE_FLAGS or col["name"] in (table.get("id_candidates") or [])[:1]:
        return None
    if C.column_problem(col["name"]):
        return None
    levels = None
    if col.get("levels") and col.get("dtype") in ("bool", "category"):
        levels = []
        for v, n in col["levels"]:
            v = str(v)
            if v == f"<{small_cell} other":
                levels.append([v, None])            # collapsed small levels: never a count (it could be < 11)
                continue
            if n is None or n < small_cell:
                continue
            if v.startswith("<"):
                levels.append([v, int(n)])
                continue
            cv = _clean_text(v)
            if cv is not None:
                levels.append([cv, int(n)])
    choices = None
    if col.get("choices"):
        choices = {str(k): _clean_text(v) for k, v in col["choices"].items() if _clean_text(v) is not None}
    item = {"key": key, "table": table["name"], "column": col["name"], "label": _clean_text(col.get("label")),
            "dtype": col.get("dtype"), "unit": col.get("unit"), "share_missing": col.get("share_missing"),
            "numeric": col.get("numeric") if col.get("dtype") in ("int", "float") else None,
            "levels": levels, "choices": choices, "code_pattern": col.get("code_pattern"),
            "field_type": col.get("field_type"), "form": col.get("form")}
    return {k: v for k, v in item.items() if v not in (None, [], {})}


def build_payload(items: list[dict], context: dict) -> dict:
    return {"dataset_context": context, "allowed_roles": LLM_ROLES, "columns": items}


def log_request(out: str | Path | None, record: dict) -> None:
    if out is None:
        return
    p = Path(out)
    p.mkdir(parents=True, exist_ok=True)
    with open(p / "llm_requests.jsonl", "a") as fh:
        fh.write(json.dumps(record, sort_keys=True, default=str) + "\n")


# ------------------------------------------------------------------------------------------------------------------
# Claude
# ------------------------------------------------------------------------------------------------------------------
def claude_available() -> tuple[bool, str]:
    try:
        import anthropic  # noqa: F401
    except ImportError:
        return False, ("the 'anthropic' package is not installed (optional extra: uv sync --extra agent); "
                       "skipping the LLM step")
    if not (os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN")):
        return False, "ANTHROPIC_API_KEY is not set; skipping the LLM step"
    return True, "ok"


class ClaudeProvider:
    """Claude through the official SDK, JSON via structured outputs (`output_config.format`, json_schema)."""

    name = "claude"

    def __init__(self, model: str | None = None, *, client=None, max_tokens: int = 16000, effort: str = "high",
                 fallbacks: bool = True):
        self.model = model or os.environ.get(MODEL_ENV) or DEFAULT_MODEL
        self.max_tokens = max_tokens
        self.effort = effort
        self.fallbacks = fallbacks
        self._client = client
        if client is None:
            ok, why = claude_available()
            if not ok:
                raise LLMUnavailable(why)

    @property
    def client(self):
        if self._client is None:
            import anthropic
            self._client = anthropic.Anthropic()
        return self._client

    def request_params(self, payload: dict) -> dict:
        params = {
            "model": self.model,
            "max_tokens": self.max_tokens,
            "system": SYSTEM_PROMPT,
            "messages": [{"role": "user", "content": json.dumps(payload, sort_keys=True)}],
            "output_config": {"effort": self.effort,
                              "format": {"type": "json_schema", "schema": PROPOSAL_SCHEMA}},
        }
        if self.fallbacks:
            params["betas"] = [FALLBACK_BETA]
            params["fallbacks"] = "default"
        return params

    def request(self, payload: dict) -> dict:
        resp = self.client.beta.messages.create(**self.request_params(payload))
        if getattr(resp, "stop_reason", None) == "refusal":
            det = getattr(resp, "stop_details", None)
            raise LLMRefused(f"model declined (category {getattr(det, 'category', None)})")
        if getattr(resp, "stop_reason", None) == "max_tokens":
            raise RuntimeError("model output hit max_tokens; proposals incomplete")
        text = next((b.text for b in resp.content if getattr(b, "type", None) == "text"), None)
        if text is None:
            raise RuntimeError("no text block in the model response")
        return json.loads(text)


def get_provider(spec, *, model: str | None = None, echo=print):
    """None / 'none' -> None (no model is contacted); 'claude' -> ClaudeProvider, or None with a message when it
    cannot run here; a provider object -> itself."""
    if spec is None or (isinstance(spec, str) and spec.lower() in ("", "none", "off", "false")):
        return None
    if isinstance(spec, str):
        if spec.lower() != "claude":
            raise ValueError(f"unknown LLM provider {spec!r}; supported: none, claude")
        try:
            return ClaudeProvider(model=model)
        except LLMUnavailable as e:
            echo(f"LLM step skipped: {e}")
            return None
    if hasattr(spec, "request"):
        return spec
    raise TypeError(f"not an LLM provider: {spec!r}")


def run_proposals(provider, items: list[dict], context: dict, out: str | Path | None,
                  batch_size: int = BATCH_SIZE) -> tuple[dict[str, dict], list[str]]:
    """Send `items` in batches; return ({key: proposal}, errors). Each payload is logged before it is sent."""
    got: dict[str, dict] = {}
    errors: list[str] = []
    for i in range(0, len(items), batch_size):
        batch = items[i:i + batch_size]
        payload = build_payload(batch, context)
        rec = {"time_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
               "provider": getattr(provider, "name", type(provider).__name__),
               "model": getattr(provider, "model", None), "batch": i // batch_size, "n_columns": len(batch),
               "system_prompt_sha256": hashlib.sha256(SYSTEM_PROMPT.encode()).hexdigest(),
               "response_schema_sha256": hashlib.sha256(json.dumps(PROPOSAL_SCHEMA, sort_keys=True).encode()).hexdigest(),
               "payload": payload}
        log_request(out, rec)
        try:
            resp = provider.request(payload)
        except Exception as e:  # noqa: BLE001 - any provider failure leaves the columns unmapped, with the reason
            errors.append(f"batch {i // batch_size}: {type(e).__name__}: {e}")
            continue
        for p in (resp or {}).get("proposals", []) or []:
            if isinstance(p, dict) and isinstance(p.get("key"), str):
                got.setdefault(p["key"], p)
    return got, errors


def sanitize_token(s: str | None, max_len: int = 60) -> str | None:
    if not s:
        return None
    t = re.sub(r"[^A-Za-z0-9_.:-]+", "_", str(s)).strip("_")[:max_len]
    return t or None
