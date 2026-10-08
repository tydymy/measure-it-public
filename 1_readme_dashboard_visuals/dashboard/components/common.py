"""Shared page scaffolding for the dashboard: page setup, tool calls with a per-page log, UNKNOWN handling, the
'Limitations and provenance' panel and the trace_evidence explorer.

Every answer on a page comes from `measure_it.tools` (imported, not over HTTP); the dense map layers read the same
processed tables through `components.data`. Nothing here computes a new number.
"""
from __future__ import annotations

import json
import time
from typing import Any

import pandas as pd
import streamlit as st

from measure_it import tools as T
from measure_it.config import UNKNOWN

EVIDENCE_LEVELS = {"A": "direct condition measure", "B": "closely matching coded condition",
                   "C": "symptom/comorbidity proxy", "D": "no usable burden estimate"}
FRAMING = ("Every ranked region is a **candidate deployment opportunity** for pilot evaluation, not a validated "
           "diagnostic pathway. Facilities are listed because they have characteristics suggesting they may be viable "
           "implementation or study partners, not as a quality ranking.")
LAYER_NOTE = ("Four data layers stay separate: person-level public cohorts are summarised as group statistics and are "
              "never located; molecular evidence is condition-level molecular enrichment; geographic joins are "
              "ecological (places, not people); facilities, trials and grants are registries.")
CSS = """
<style>
.mi-badge {display:inline-block; padding:0.05rem 0.45rem; border-radius:0.6rem; font-size:0.78rem; font-weight:600;
           border:1px solid rgba(128,128,128,0.45); margin-right:0.25rem; white-space:nowrap}
.mi-lvl-A {background:rgba(42,120,214,0.14)} .mi-lvl-B {background:rgba(42,120,214,0.09)}
.mi-lvl-C {background:rgba(235,104,52,0.14)} .mi-lvl-D {background:rgba(128,128,128,0.16)}
.mi-unknown {background:rgba(128,128,128,0.16)}
div[data-testid="stMetricValue"] {font-size:1.45rem}
</style>
"""


# ------------------------------------------------------------------------------------------------ page setup

def setup_page(title: str, intro: str, icon: str = ":material/monitoring:") -> None:
    try:
        st.set_page_config(page_title=f"{title} | Measure It to Cure It", page_icon=icon, layout="wide")
    except Exception:  # noqa: BLE001 - set_page_config may only run once per page in some hosts
        pass
    st.markdown(CSS, unsafe_allow_html=True)
    st.title(title)
    st.caption(intro)
    with st.sidebar:
        st.caption("Measure It to Cure It: public-data alpha. Every answer comes from the shared tool facade "
                   "(`measure_it.tools`); UNKNOWN / NOT AVAILABLE means the data are missing, never a guess.")
        st.caption(f"Data stamp: {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(T.data_stamp()))}")


# ------------------------------------------------------------------------------------------------ tool calls

class ToolLog:
    """Calls measure_it.tools functions and keeps every envelope, for the page's provenance panel."""

    def __init__(self) -> None:
        self.calls: list[dict] = []
        self.direct_tables: dict[str, str] = {}

    def call(self, name: str, **kwargs) -> dict:
        env = T.TOOLS[name](**kwargs)
        self.calls.append(env)
        return env

    def table(self, name: str, why: str) -> None:
        """Record a processed table the page reads directly (map layers)."""
        self.direct_tables[name] = why


def ok(env: dict | None) -> bool:
    return bool(env) and env.get("status") == "ok"


def unknown_box(env: dict, what: str) -> None:
    """Render a non-ok envelope: UNKNOWN / NOT AVAILABLE (or error) with the tool's own reason."""
    status = env.get("status", UNKNOWN)
    reason = env.get("reason") or "no reason returned"
    if status == "error":
        st.error(f"**{what}: tool error** ({env.get('tool')}): {reason}")
    else:
        st.warning(f"**{what}: {UNKNOWN}**. {reason}")
    data = env.get("data") or {}
    for key in ("condition_resolution", "measurement_resolution"):
        res = data.get(key) if isinstance(data, dict) else None
        if isinstance(res, dict) and res.get("candidates"):
            st.caption(f"Candidates from the normaliser: {', '.join(map(str, res['candidates']))}")


def is_unknown(v) -> bool:
    return v is None or v == UNKNOWN or (isinstance(v, float) and v != v)


def fmt(v, nd: int = 2, unit: str = "") -> str:
    """Number -> text; UNKNOWN / None stay 'UNKNOWN / NOT AVAILABLE'."""
    if is_unknown(v):
        return UNKNOWN
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, (int, float)):
        if float(v).is_integer() and abs(v) >= 1:
            return f"{int(v):,}{unit}"
        return f"{v:,.{nd}f}{unit}"
    return str(v)


def level_badge(level) -> str:
    lv = str(level) if level not in (None, "") else "D"
    lv = lv if lv in EVIDENCE_LEVELS else "D"
    return f'<span class="mi-badge mi-lvl-{lv}">level {lv}: {EVIDENCE_LEVELS[lv]}</span>'


def badge(text: str) -> str:
    return f'<span class="mi-badge">{text}</span>'


def records_df(rows: list[dict] | None, cols: list[str] | None = None) -> pd.DataFrame:
    df = pd.DataFrame(rows or [])
    if cols:
        df = df[[c for c in cols if c in df.columns]]
    for c in df.columns:  # lists/dicts -> readable text so st.dataframe can show them
        if df[c].map(lambda x: isinstance(x, (list, dict))).any():
            df[c] = df[c].map(lambda x: "; ".join(map(str, x)) if isinstance(x, list)
                              else (json.dumps(x) if isinstance(x, dict) else x))
    return df


def show_df(df: pd.DataFrame, height: int | str = "auto", key: str | None = None, **kw) -> None:
    if df is None or len(df) == 0:
        st.caption(f"{UNKNOWN}: no rows.")
        return
    df = df.copy()
    for c in df.columns:  # mixed columns (numbers + the UNKNOWN sentinel) are shown as text, not coerced
        if df[c].dtype == object:
            kinds = {type(x) for x in df[c] if x is not None and not (isinstance(x, float) and x != x)}
            if len(kinds) > 1 or kinds & {list, dict, tuple}:
                df[c] = df[c].map(lambda x: "" if x is None or (isinstance(x, float) and x != x) else str(x))
    if height == "auto":   # show every row of short tables; scroll longer ones
        height = 35 * (len(df) + 1) + 3 if len(df) <= 16 else 420
    st.dataframe(df, hide_index=True, height=height, key=key, **kw)


def caveat_list(caveats: list[str], limit: int | None = None) -> None:
    items = list(dict.fromkeys(c for c in caveats if isinstance(c, str) and c.strip()))
    if limit:
        items = items[:limit]
    if items:
        st.markdown("\n".join(f"- {c}" for c in items))


# ------------------------------------------------------------------------------------------------ provenance

def _sources(envs: list[dict]) -> pd.DataFrame:
    rows = {}
    for e in envs:
        for s in (e.get("provenance") or {}).get("sources") or []:
            if isinstance(s, dict) and s.get("source_id") and s["source_id"] not in rows:
                rows[s["source_id"]] = {k: s.get(k) for k in ("source_id", "name", "source_version", "retrieved_at",
                                                               "data_audit")}
    return pd.DataFrame(list(rows.values()))


def limitations_panel(log: ToolLog, notes: list[str]) -> None:
    """The 'Limitations and provenance' panel every page ends with."""
    st.divider()
    st.subheader("Limitations and provenance")
    st.markdown("\n".join(f"- {n}" for n in notes + [LAYER_NOTE]))
    envs = log.calls
    unknowns = [e for e in envs if e.get("status") != "ok"]
    if unknowns:
        st.markdown("**Missing data on this page (returned as UNKNOWN / NOT AVAILABLE, never filled in):**")
        st.markdown("\n".join(f"- `{e.get('tool')}`: {e.get('reason')}" for e in unknowns))
    cav = [c for e in envs for c in e.get("caveats") or []]
    with st.expander(f"Caveats returned by the tools ({len(set(cav))})"):
        caveat_list(cav)
    with st.expander(f"Tool calls ({len(envs)}), sources and tables"):
        calls = pd.DataFrame([{
            "tool": e.get("tool"), "status": e.get("status"),
            "query": json.dumps(e.get("query"), default=str)[:160],
            "producer": (e.get("provenance") or {}).get("producer"),
            "tables": ", ".join((e.get("provenance") or {}).get("tables") or []),
            "n_object_ids": (e.get("provenance") or {}).get("n_object_ids"),
            "lists_capped": len(e.get("truncation") or []),
            "elapsed_s": (e.get("meta") or {}).get("elapsed_s")} for e in envs])
        show_df(calls, key=f"calls_{id(log)}")
        src = _sources(envs)
        if len(src):
            st.markdown("**Sources (SOURCE_REGISTRY.yaml; each has a DATA_AUDIT.md):**")
            show_df(src, key=f"src_{id(log)}")
        if log.direct_tables:
            st.markdown("**Processed tables read directly for map layers** (same tables the tools read; display "
                        "only):")
            st.markdown("\n".join(f"- `data/processed/{t}`: {why}" for t, why in log.direct_tables.items()))
        trunc = [dict(t, tool=e.get("tool")) for e in envs for t in e.get("truncation") or []]
        if trunc:
            st.markdown("**Lists capped by the tools** (the full count is kept):")
            show_df(pd.DataFrame(trunc)[["tool", "path", "returned", "total"]].head(40), key=f"trunc_{id(log)}")
        st.caption(f"Data stamp {T.data_stamp():.0f} (newest mtime of data/processed, configs, results/tables, "
                   "SOURCE_REGISTRY.yaml). Every object id resolves with trace_evidence.")


def trace_explorer(object_ids: list[str], key: str, log: ToolLog | None = None, label: str | None = None) -> None:
    """Expander: pick an object id and resolve it with trace_evidence (rows, sources, DATA_AUDIT.md, lineage)."""
    ids = [i for i in dict.fromkeys(object_ids) if isinstance(i, str) and i]
    with st.expander(label or f"Trace evidence: resolve any of {len(ids)} cited object ids"):
        if not ids:
            st.caption(f"{UNKNOWN}: no object ids to trace.")
            return
        oid = st.selectbox("Object id", ids, key=key)
        env = (log.call if log else lambda n, **k: T.TOOLS[n](**k))("trace_evidence", object_id=oid, max_rows=5)
        if not ok(env):
            unknown_box(env, f"trace_evidence({oid})")
            return
        d = env["data"]
        st.markdown(f"`{d.get('object_id')}`: table **{d.get('primary_table')}**, {d.get('n_rows')} row(s), "
                    f"lookup via {d.get('lookup_backend')}")
        rows = d.get("rows") or []
        if rows:
            r0 = {k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in rows[0].items()}
            show_df(pd.DataFrame({"field": list(r0.keys()), "value": [str(v)[:300] for v in r0.values()]}),
                    height=260, key=f"{key}_row")
        src = d.get("sources") or []
        if src:
            show_df(records_df(src, ["source_id", "name", "source_version", "retrieved_at", "data_audit",
                                     "manifest"]), key=f"{key}_src")
        lin = d.get("lineage") or []
        lin_ids = [x.get("object_id") for x in lin if isinstance(x, dict) and x.get("object_id")]
        if lin_ids:
            st.caption("Lineage (one level): " + ", ".join(f"`{i}`" for i in lin_ids[:20])
                       + (" ..." if len(lin_ids) > 20 else ""))
        for n in (d.get("notes") or [])[:4] if isinstance(d.get("notes"), list) else []:
            st.caption(str(n))


def burden_caveat_box(log: ToolLog, members: list[str], level: str, labels: dict | None = None) -> list[dict]:
    """Always-visible burden caveats: per member condition the evidence level, the source resolution, whether a state
    value is inherited by counties, and whether the measure is a proxy (from get_condition_burden)."""
    labels = labels or {}
    lines, out = [], []
    for cid in members:
        env = log.call("get_condition_burden", condition=cid, geography_level=level, max_rows=1)
        d = env.get("data") or {}
        lvl = d.get("burden_evidence_level") or "D"
        row = (d.get("rows") or [{}])[0] if d.get("rows") else {}
        name = labels.get(cid, cid)
        if lvl == "D" or env.get("status") != "ok":
            txt = (f"**{name}: level D, no usable burden estimate** ({UNKNOWN}). {d.get('rationale') or env.get('reason')}"
                   " The burden layer is empty, the diagnostic desert is access-only and the opportunity score drops "
                   "the burden weight (the other weights are renormalised).")
        else:
            src = row.get("source_name") or ""
            meas = row.get("measure_label") or d.get("primary_measure_id")
            txt = f"**{name}: level {lvl} ({EVIDENCE_LEVELS.get(lvl, '')})**, {meas}"
            txt += f" [{src}]" if src else ""
            txt += f"; source resolution: {d.get('primary_source_resolution')}."
            if d.get("inherited") and level == "county":
                txt += (" At county level every county carries its **state** estimate: inherited context with no "
                        "within-state variation, **not county prevalence**.")
            if lvl == "C":
                txt += " This is a **symptom/comorbidity proxy**, not the condition's prevalence."
            elif lvl == "B":
                txt += " A closely matching coded condition (claims), not a direct measure."
            if d.get("rationale"):
                txt += f" {d['rationale']}"
        lines.append(txt)
        out.append({"condition_id": cid, "level": lvl, "inherited": bool(d.get("inherited")),
                    "status": env.get("status")})
    st.warning("**Burden evidence (always shown).** " + "  \n".join(lines), icon=":material/info:")
    return out


def sections(names: list[str], numbered: bool = True) -> list:
    """Stacked, numbered page sections (containers with a subheader) instead of tabs, so every part of a page is
    visible at once (and in a full-page screenshot). A short index line links the reader through them."""
    st.markdown(" | ".join(f"**{i + 1}.** {n}" if numbered else n for i, n in enumerate(names)))
    out = []
    for i, n in enumerate(names):
        c = st.container()
        c.subheader(f"{i + 1}. {n}" if numbered else n, divider="gray")
        out.append(c)
    return out


def json_bytes(obj: Any) -> bytes:
    return json.dumps(obj, indent=1, allow_nan=False, default=str).encode()
