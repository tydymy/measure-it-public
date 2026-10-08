"""Static, read-only reviewer snapshot: one self-contained HTML page with the engine's key views.

    uv run measure-it export-static                    # dist/static_snapshot/index.html

The page reproduces, read-only: the overview (what the engine is and is not, the four data layers and the guardrails,
from README.md), the primary deployment query (Long COVID or ME/CFS x wearable autonomic monitoring: ranked table with
every component, Monte Carlo rank intervals, the evidence_weighted vs equal comparison and the ranks under every weight
set), a county choropleth of the equal-weight composite with the top regions and their candidate sites, evidence
cards for the top 10 (candidate sites with reasons, research evidence ids, uncertainties), the measurement performance
table of the metric link (UNKNOWN shown as UNKNOWN), results/FINAL_REPORT.md §14 (what separates cases from controls
and what does not), the validation tests of §13 and the limitations of §19.

Every number is read at build time: the ranking from ``data/processed/deployment_opportunities.parquet`` and the tool
facade (``measure_it.tools.rank_deployment_opportunities``, which reads the same table), performance from
``data/processed/measurement_performance.parquet``, the report sections from results/FINAL_REPORT.md, boundaries from
the Census geoparquet. Nothing is typed in by hand. The only external request is one pinned script (d3, for the map)
from cdn.jsdelivr.net; without it the page still shows every table, and the map says it could not load.
"""
from __future__ import annotations

import html
import json
import math
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import PROCESSED, PROJECT_ROOT, RESULTS, SOURCE_REGISTRY_PATH, UNKNOWN, load_config

OUT = PROJECT_ROOT / "dist" / "static_snapshot" / "index.html"
TITLE = "Measure It Snapshot"
D3_URL = "https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"
REASONS = {"P": "not ranked: fewer than the minimum population", "I": "not ranked: incomplete burden",
           "N": "not ranked"}
ALLOWED_SCRIPT_HOSTS = ("cdn.jsdelivr.net", "cdnjs.cloudflare.com")
MAX_BYTES = 12 * 1024 * 1024
PRIMARY = {"condition": "Long COVID or ME/CFS", "measurement": "wearable autonomic monitoring",
           "geography_level": "county", "weight_set": "equal", "top_n": 10}
PRIMARY_IDS = ("long_covid_or_me_cfs", "wearable_autonomic_activity_monitoring")
COMPONENTS = [("burden_pct", "burden"), ("vulnerability_pct", "vulnerability"), ("diagnostic_desert_pct", "desert"),
              ("clinic_capacity_pct", "capacity"), ("research_readiness_pct", "readiness")]
SITES_PER_CARD = 5
SITES_ON_MAP = 3
NON_STATE_FIPS = {"60", "66", "69", "72", "78"}
REPORT = RESULTS / "FINAL_REPORT.md"
README = PROJECT_ROOT / "README.md"


# ------------------------------------------------------------------------------------------------ small helpers

def esc(x) -> str:
    return html.escape("" if x is None else str(x), quote=True)


def _isnum(x) -> bool:
    return isinstance(x, (int, float, np.integer, np.floating)) and not (isinstance(x, float) and math.isnan(x)) \
        and not (isinstance(x, np.floating) and np.isnan(x))


def f3(x) -> str:
    return f"{float(x):.3f}" if _isnum(x) else "UNKNOWN"


def f2(x) -> str:
    return f"{float(x):.2f}" if _isnum(x) else "UNKNOWN"


def fint(x) -> str:
    return f"{int(round(float(x))):,}" if _isnum(x) else "UNKNOWN"


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if math.isnan(float(o)) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def _script_json(o) -> str:
    """JSON safe inside <script type=application/json> (no '</' sequences)."""
    return json.dumps(_jsonable(o), separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")


def _git() -> dict:
    def run(*a):
        try:
            r = subprocess.run(["git", *a], cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=20)
            return r.stdout.strip() if r.returncode == 0 else None
        except (OSError, subprocess.SubprocessError):
            return None
    return {"commit": run("rev-parse", "--short", "HEAD") or UNKNOWN, "describe": run("describe", "--always", "--dirty")}


# ------------------------------------------------------------------------------------------------ markdown

def _md():
    from markdown_it import MarkdownIt
    return MarkdownIt("commonmark", {"html": False, "linkify": False, "typographer": False}).enable("table")


LINK_RE = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")


def md_to_html(text: str) -> str:
    """Markdown -> HTML with links flattened to their text (the page makes no requests and links nowhere)."""
    text = LINK_RE.sub(r"\1", text)
    # a wrapped line that starts with "+ 0.05" or "+ BMI" is prose, not a list item
    text = re.sub(r"(?m)^(\s+)\+ ", r"\1\\+ ", text)
    out = _md().render(text)
    # tables scroll inside their own container
    return out.replace("<table>", '<div class="scroll"><table class="md">').replace("</table>", "</table></div>")


def md_section(text: str, heading_re: str, stop_re: str = r"^#{1,3} ") -> str:
    """The markdown under the first heading matching heading_re, up to the next heading matching stop_re."""
    lines = text.splitlines()
    start = next((i for i, l in enumerate(lines) if re.match(heading_re, l)), None)
    if start is None:
        raise ValueError(f"heading {heading_re!r} not found")
    end = next((j for j in range(start + 1, len(lines)) if re.match(stop_re, lines[j])), len(lines))
    body = "\n".join(lines[start + 1:end]).strip()
    return re.sub(r"\n---\s*$", "", body).strip()


def readme_parts() -> dict:
    t = README.read_text()
    lines = t.splitlines()
    # the introduction: from "The engine asks" up to the paragraph that starts with **It is not** (inclusive)
    i0 = next(i for i, l in enumerate(lines) if l.startswith("The engine asks"))
    i1 = next(i for i, l in enumerate(lines) if l.startswith("**It is not**"))
    i2 = next(j for j in range(i1, len(lines)) if not lines[j].strip())
    intro = "\n".join(lines[i0:i2])
    layers = md_section(t, r"^## The four data layers", r"^(## |---)")
    return {"intro": intro, "layers": layers}


def _blockquote(text: str, start_re: str) -> list[str]:
    """The block-quoted lines starting at the first line matching start_re (the report's dated update note)."""
    lines = text.splitlines()
    i = next((k for k, l in enumerate(lines) if re.match(start_re, l)), None)
    if i is None:
        return []
    out = []
    for l in lines[i:]:
        if not l.startswith(">"):
            break
        out.append(l)
    return out


def changes_html() -> str:
    """Section 'What changed': the capillaroscopy top 10 under the current defaults and the S7-S9 agreement table."""
    d = pd.read_parquet(PROCESSED / "deployment_opportunities.parquet",
                        columns=["condition_id", "measurement_id", "geo_level", "geo_name", "rank_equal",
                                 "composite_equal", "burden_pct", "clinic_capacity_pct", "clinic_capacity_basis",
                                 "research_readiness_pct", "rank_mc_p05", "rank_mc_p95", "p_top10_mc"])
    d = d[(d["condition_id"] == "long_covid_or_me_cfs") & (d["measurement_id"] == "nailfold_capillaroscopy")
          & (d["geo_level"] == "county") & (d["rank_equal"] <= 10)].sort_values("rank_equal")
    rows = "".join(
        f"<tr><td class='num'>{int(r.rank_equal)}</td><td class='l'><b>{esc(r.geo_name)}</b></td>"
        f"<td class='num strong'>{r.composite_equal:.3f}</td><td class='num'>{r.burden_pct:.2f}</td>"
        f"<td class='num'>{r.clinic_capacity_pct:.2f}</td><td class='num'>{r.research_readiness_pct:.2f}</td>"
        f"<td class='num'>{int(round(r.rank_mc_p05))}-{int(round(r.rank_mc_p95))}</td>"
        f"<td class='num'>{r.p_top10_mc:.2f}</td></tr>" for r in d.itertuples())
    cap = "<table class='data'><thead><tr><th>rank</th><th class='l'>county</th><th>composite</th><th>burden pct</th>" \
          "<th>capacity pct</th><th>readiness pct</th><th>rank interval 5-95</th><th>P(top 10)</th></tr></thead>" \
          f"<tbody>{rows}</tbody></table>"
    sl_path = RESULTS / "tables" / "stage_layers_sensitivity.csv"
    sl = ""
    if sl_path.exists():
        t = pd.read_csv(sl_path)
        lab = {"S7_legacy_hps_inherited": "S7: old Long COVID burden (inherited state value)",
               "S8_legacy_density": "S8: old clinic capacity (provider density)",
               "S9_legacy_both": "S9: both old rules (the ranking before 2026-10-07)"}
        body = "".join(
            f"<tr><td class='l'>{esc(lab.get(r.analysis, r.analysis))}</td><td class='l'>{esc(r.condition_id)} x "
            f"{esc(r.measurement_id)}</td><td class='num'>{int(r.top10_overlap)} of 10</td>"
            f"<td class='num'>{int(r.top25_overlap)} of 25</td><td class='num'>{r.kendall_tau_b:.2f}</td></tr>"
            for r in t.itertuples())
        sl = "<table class='data compact'><thead><tr><th class='l'>variant vs current ranking</th>" \
             "<th class='l'>query</th><th>top-10 kept</th><th>top-25 kept</th><th>Kendall tau-b</th></tr></thead>" \
             f"<tbody>{body}</tbody></table>"
    basis = (d["clinic_capacity_basis"].iloc[0] if len(d) else UNKNOWN)
    return (f"<h3>Long COVID or ME/CFS &times; nailfold capillaroscopy: top 10 counties</h3>"
            f"<p class='muted'>Equal weights. Clinic capacity basis: <code>{esc(basis)}</code> (capillaroscopy has no "
            "billing code, so it keeps provider density). Long COVID burden is the modelled BRFSS 2023 small-area "
            "estimate, not observed county prevalence.</p>"
            f"<div class='scroll'>{cap}</div>"
            "<h3>The old rules against the current ranking</h3>"
            "<p class='muted'>results/STAGE_LAYERS_SCORING.md: each variant puts a rule from before 2026-10-07 back. "
            "For capillaroscopy S8 changes nothing because its capacity is still density.</p>"
            f"<div class='scroll'>{sl}</div>")


def report_parts() -> dict:
    t = REPORT.read_text()
    return {
        "s13": md_section(t, r"^## 13\. ", r"^## "),
        "s14_intro": md_section(t, r"^## 14\. ", r"^#{2,3} "),
        "s14_2": md_section(t, r"^### 14\.2 ", r"^#{2,3} "),
        "s14_3": md_section(t, r"^### 14\.3 ", r"^#{2,3} "),
        "s14_4": md_section(t, r"^### 14\.4 ", r"^#{2,3} "),
        "s14_5": md_section(t, r"^### 14\.5 ", r"^#{2,3} "),
        "s19": md_section(t, r"^## 19\. ", r"^## "),
        "update": "\n".join(l[2:] if l.startswith("> ") else l[1:] for l in
                            _blockquote(t, r"^> \*\*Update ")),
        "headings": {k: next((l.lstrip("# ").strip() for l in t.splitlines() if re.match(p, l)), k)
                     for k, p in (("s13", r"^## 13\. "), ("s14", r"^## 14\. "), ("s14_2", r"^### 14\.2 "),
                                  ("s14_3", r"^### 14\.3 "), ("s14_4", r"^### 14\.4 "), ("s14_5", r"^### 14\.5 "),
                                  ("s19", r"^## 19\. "))},
    }


# ------------------------------------------------------------------------------------------------ data

def opportunity_rows() -> pd.DataFrame:
    cid, mid = PRIMARY_IDS
    return pd.read_parquet(PROCESSED / "deployment_opportunities.parquet",
                           filters=[("condition_id", "==", cid), ("measurement_id", "==", mid),
                                    ("geo_level", "==", "county")])


def primary_envelope() -> dict:
    from .. import tools as T
    env = T.rank_deployment_opportunities(condition=PRIMARY["condition"], measurement=PRIMARY["measurement"],
                                          geography_level=PRIMARY["geography_level"], top_n=PRIMARY["top_n"],
                                          weight_set=PRIMARY["weight_set"], detail="full")
    if env.get("status") != "ok":
        raise RuntimeError(f"rank_deployment_opportunities returned {env.get('status')}: {env.get('reason')}")
    return env


def performance_table() -> pd.DataFrame:
    p = PROCESSED / "measurement_performance.parquet"
    df = pd.read_parquet(p) if p.exists() else pd.read_csv(RESULTS / "tables" / "measurement_performance.csv")
    return df


def labels() -> tuple[dict, dict]:
    cond = {c["id"]: c["preferred_name"] for c in load_config("conditions")["conditions"]}
    try:
        from ..scoring.opportunity import condition_sets
        cond.update({k: v["label"] for k, v in condition_sets().items()})
    except Exception:  # noqa: BLE001
        pass
    meas = {}
    p = PROCESSED / "measurement_registry.parquet"
    if p.exists():
        m = pd.read_parquet(p, columns=["measurement_id", "name"])
        meas = dict(zip(m["measurement_id"], m["name"]))
    return cond, meas


def map_geometry(state_tol: float = 0.15, nd: int = 1) -> dict:
    """Embedded map geometry: simplified 2024 state outlines plus one dot per 2024 county at its Census internal point.

    County polygons are not embedded (they made the page 1.7 MB); counties are drawn as dots at their internal points.
    """
    states = boundaries(state_tol, nd)["state"]
    g = pd.read_parquet(PROCESSED / "geographies.parquet", columns=["geo_id", "geo_level", "name", "state_abbr",
                                                                    "lat", "lon", "vintage"])
    g = g[(g["geo_level"] == "county") & (g["vintage"] == "2024") & ~g["geo_id"].str[:2].isin(NON_STATE_FIPS)]
    dots = [[gid, round(float(lo), nd), round(float(la), nd)] for gid, lo, la in zip(g["geo_id"], g["lon"], g["lat"])]
    return {"state": states, "dots": dots}


def boundaries(tol: float = 0.015, nd: int = 2, levels: tuple = ("state",)) -> dict:
    """Simplified county + state polygons (clockwise exterior rings, as d3-geo expects), 2-decimal coordinates."""
    import geopandas as gpd
    import shapely
    from shapely.geometry import mapping
    from shapely.geometry.polygon import orient

    def load(level: str, t: float):
        g = gpd.read_parquet(PROCESSED / f"{level}_boundaries_2024.geoparquet")
        if g.crs is not None and g.crs.to_epsg() != 4326:
            g = g.to_crs(4326)
        g = g[~g["geo_id"].str[:2].isin(NON_STATE_FIPS)].sort_values("geo_id").reset_index(drop=True)
        geoms = shapely.simplify(g.geometry.values, t, preserve_topology=True)
        out = []
        for gid, name, geom in zip(g["geo_id"], g["name"], geoms):
            if geom is None or geom.is_empty:
                continue
            polys = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
            coords = []
            for p in polys:
                p = orient(p, sign=-1.0)   # clockwise exterior: d3's spherical winding
                m = mapping(p)["coordinates"]
                coords.append([[[round(x, nd), round(y, nd)] for x, y in ring] for ring in m])
            out.append({"type": "Feature", "id": gid, "properties": {"n": name},
                        "geometry": {"type": "MultiPolygon", "coordinates": coords}})
        return {"type": "FeatureCollection", "features": out}

    return {lv: load(lv, tol) for lv in levels}


def collect() -> dict:
    """Everything the page shows, read from the processed tables and results files."""
    rows = opportunity_rows()
    env = primary_envelope()
    d = env["data"]
    ranked = rows[rows["rank_equal"].notna()].copy()
    top = ranked.sort_values("rank_equal").head(PRIMARY["top_n"])
    top_ew = rows[rows["rank_evidence_weighted"].notna()].sort_values("rank_evidence_weighted").head(PRIMARY["top_n"])
    weight_sets = list(load_config("scoring")["weight_sets"])
    cond_labels, meas_labels = labels()
    perf = performance_table()
    try:
        import yaml
        n_sources = len((yaml.safe_load(SOURCE_REGISTRY_PATH.read_text()) or {}).get("sources", []))
    except Exception:  # noqa: BLE001
        n_sources = None
    n_tables = len(list(PROCESSED.glob("*.parquet")))
    return {"rows": rows, "ranked": ranked, "top": top, "top_ew": top_ew, "env": env, "data": d,
            "weight_sets": weight_sets, "cond_labels": cond_labels, "meas_labels": meas_labels, "perf": perf,
            "readme": readme_parts(), "report": report_parts(), "n_sources": n_sources, "n_tables": n_tables,
            "git": _git(), "built_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat()}


# ------------------------------------------------------------------------------------------------ rendering pieces

def interval_svg(p05, p95, rank, n: int, width: int = 150) -> str:
    """Rank interval on a log scale from 1 to n: a bar from p05 to p95 and a dot at the rank."""
    if not (_isnum(p05) and _isnum(p95) and _isnum(rank)) or n <= 1:
        return "UNKNOWN"
    lo, hi = math.log10(1), math.log10(n)

    def x(v):
        return 4 + (width - 8) * (math.log10(max(1.0, min(float(v), n))) - lo) / (hi - lo)
    ticks = [t for t in (1, 10, 100, 1000) if t <= n]
    tk = "".join(f'<line class="tick" x1="{x(t):.1f}" x2="{x(t):.1f}" y1="12" y2="16"/>' for t in ticks)
    return (f'<svg class="ivl" viewBox="0 0 {width} 18" width="{width}" height="18" role="img" '
            f'aria-label="rank interval {fint(p05)} to {fint(p95)} of {n}">'
            f'<line class="axis" x1="4" x2="{width - 4}" y1="14" y2="14"/>{tk}'
            f'<line class="band" x1="{x(p05):.1f}" x2="{x(p95):.1f}" y1="8" y2="8"/>'
            f'<circle class="dot" cx="{x(rank):.1f}" cy="8" r="3.2"/></svg>')


def pct_cell(v) -> str:
    if not _isnum(v):
        return '<td class="num unk">UNKNOWN</td>'
    w = max(0.0, min(1.0, float(v))) * 100
    return f'<td class="num pct"><span class="bar" style="--w:{w:.1f}%"></span><span>{float(v):.2f}</span></td>'


def burden_chip(r) -> str:
    lvl = r.get("burden_evidence_level") or "D"
    inh = bool(r.get("burden_inherited")) if r.get("burden_inherited") is not None else False
    mc = r.get("member_components")
    members = []
    try:
        for cid, m in (json.loads(mc) if isinstance(mc, str) else {}).items():
            members.append(f"{cid} {m.get('burden_evidence_level') or 'D'}"
                           + (" (inherited state value)" if m.get("burden_inherited") else ""))
    except ValueError:
        pass
    title = "; ".join(members)
    return (f'<span class="chip lvl-{esc(lvl)}" title="{esc(title)}">level {esc(lvl)}</span>'
            + (' <span class="chip warn" title="a member value is the state estimate carried by the county">'
               'inherited</span>' if inh or "inherited" in title else ""))


def ranked_table(c: dict) -> str:
    top, n = c["top"], int(c["data"]["n_regions_ranked"])
    head = ("<tr><th>rank</th><th class='l'>region</th><th>composite</th>"
            + "".join(f"<th>{esc(lbl)} pct</th>" for _, lbl in COMPONENTS)
            + "<th class='l'>rank interval 5-95 (1,000 draws)</th><th>P(top 10)</th><th class='l'>burden evidence</th>"
            + "<th>evidence_weighted rank</th><th>evidence_weighted composite</th></tr>")
    body = []
    for r in top.to_dict("records"):
        body.append(
            f'<tr data-geo="{esc(r["geo_id"])}" data-rank="{fint(r["rank_equal"])}" '
            f'data-composite="{f3(r["composite_equal"])}">'
            f'<td class="num">{fint(r["rank_equal"])}</td>'
            f'<td class="l region" title="{esc(r["object_id"])}"><b>{esc(r["geo_name"])}</b><br>'
            f'<code>{esc(r["geo_object_id"])}</code></td>'
            f'<td class="num strong">{f3(r["composite_equal"])}</td>'
            + "".join(pct_cell(r.get(col)) for col, _ in COMPONENTS)
            + f'<td class="l">{interval_svg(r.get("rank_mc_p05"), r.get("rank_mc_p95"), r.get("rank_equal"), n)}'
              f'<span class="ivl-t">{fint(r.get("rank_mc_p05"))}-{fint(r.get("rank_mc_p95"))}</span></td>'
            f'<td class="num">{f3(r.get("p_top10_mc"))}</td>'
            f'<td class="l">{burden_chip(r)}</td>'
            f'<td class="num">{fint(r.get("rank_evidence_weighted"))}</td>'
            f'<td class="num">{f3(r.get("composite_evidence_weighted"))}</td></tr>')
    return f'<div class="scroll"><table class="data" id="ranked"><thead>{head}</thead><tbody>{"".join(body)}</tbody></table></div>'


def comparison(c: dict) -> str:
    top, top_ew = c["top"], c["top_ew"]
    a, b = list(top["geo_id"]), list(top_ew["geo_id"])
    both = [g for g in a if g in b]
    n = len(a)

    def tbl(df, rank_col, other_col, comp_col, caption):
        rows = "".join(
            f'<tr data-geo="{esc(r["geo_id"])}"><td class="num">{fint(r[rank_col])}</td>'
            f'<td class="l">{esc(r["geo_name"])}{" <span class=\'chip ok\'>in both</span>" if r["geo_id"] in both else ""}</td>'
            f'<td class="num">{f3(r[comp_col])}</td><td class="num">{fint(r[other_col])}</td></tr>'
            for r in df.to_dict("records"))
        other = "equal rank" if rank_col != "rank_equal" else "evidence_weighted rank"
        return (f'<figure class="cmp"><figcaption>{esc(caption)}</figcaption><div class="scroll"><table class="data">'
                f'<thead><tr><th>rank</th><th class="l">region</th><th>composite</th><th>{other}</th></tr></thead>'
                f'<tbody>{rows}</tbody></table></div></figure>')
    ws = load_config("scoring")["weight_sets"]
    ew = ws.get("evidence_weighted", {})
    eq = ws.get("equal", {})
    fmt = lambda w: ", ".join(f"{k} {v:g}" for k, v in w.items())  # noqa: E731
    return (f'<p class="lede">{len(both)} of {n} regions are in both top-{n} lists. '
            f'<b>equal</b>: {esc(fmt(eq))}. <b>evidence_weighted</b>: {esc(fmt(ew))}.</p>'
            f'<div class="two">{tbl(top, "rank_equal", "rank_evidence_weighted", "composite_equal", "Top " + str(n) + " under equal weights")}'
            f'{tbl(top_ew, "rank_evidence_weighted", "rank_equal", "composite_evidence_weighted", "Top " + str(n) + " under evidence_weighted")}</div>')


def weightset_table(c: dict) -> str:
    recs = c["data"]["recommendations"]
    ws = c["weight_sets"]
    head = "<tr><th class='l'>region</th>" + "".join(f"<th>{esc(w)}</th>" for w in ws) + "</tr>"
    body = ""
    for r in recs:
        g = r["geography"]
        ranks = g.get("ranks_under_weight_sets") or {}
        body += (f'<tr><td class="l">{esc(g["name"])}</td>'
                 + "".join(f'<td class="num">{fint(ranks.get(w))}</td>' for w in ws) + "</tr>")
    return f'<div class="scroll"><table class="data compact"><thead>{head}</thead><tbody>{body}</tbody></table></div>'


def perf_callout(c: dict) -> str:
    mp = c["data"].get("measurement_performance") or {}
    status = mp.get("performance_status") or UNKNOWN
    if status != "known":
        return f'<div class="callout warn"><b>Measurement performance: {esc(status)}</b></div>'
    ci = f'{f3(mp.get("auroc"))} ({f3(mp.get("auroc_ci_low"))}-{f3(mp.get("auroc_ci_high"))})'
    op = (f'{f2(mp.get("op_sensitivity"))} ({f2(mp.get("op_sensitivity_ci_low"))}-'
          f'{f2(mp.get("op_sensitivity_ci_high"))}) at specificity {f2(mp.get("op_specificity"))}')
    return (f'<div class="callout"><div class="kv"><span>performance record</span><code>{esc(mp.get("selected_record_id"))}'
            f'</code></div><div class="kv"><span>AUROC (95% CI)</span><b>{ci}</b></div>'
            f'<div class="kv"><span>sensitivity</span><b>{op}</b></div>'
            f'<div class="kv"><span>cases / controls</span><b>{fint(mp.get("n_cases"))} / {fint(mp.get("n_controls"))}</b></div>'
            f'<div class="kv"><span>measurement_evidence</span><b>{f3(mp.get("measurement_evidence"))}</b></div>'
            f'<p class="note">{esc(mp.get("performance_note"))}</p></div>')


def cards(c: dict) -> str:
    out = []
    for i, r in enumerate(c["data"]["recommendations"]):
        g = r["geography"]
        ivl = g.get("rank_interval_5_95") or [None, None]
        sites = r.get("candidate_sites") or []
        site_html = []
        for s in sites[:SITES_PER_CARD]:
            reasons = "".join(f"<li>{esc(x)}</li>" for x in s.get("reasons") or [])
            site_html.append(
                f'<li class="site"><div class="site-h"><b>{esc(s.get("facility_name"))}</b> '
                f'<span class="muted">{esc(s.get("city"))}, {esc(s.get("state"))} &middot; '
                f'{f2(s.get("distance_km"))} km &middot; {esc(s.get("geocode_precision"))}</span></div>'
                f'<code>{esc(s.get("object_id"))}</code><ul class="reasons">{reasons}</ul></li>')
        more = len(sites) - SITES_PER_CARD
        re_ = r.get("research_evidence") or {}
        ids = []
        for key, lbl in (("condition_trials", "condition trials"), ("technology_experience_trials",
                                                                     "technology-experience trials"),
                         ("nih_projects", "NIH projects")):
            v = re_.get(key) or {}
            oids = v.get("object_ids") or []
            n = v.get("n", v.get("n_appl_ids", len(oids)))
            chips = "".join(f"<code>{esc(o)}</code>" for o in oids[:12])
            extra = f' <span class="muted">+{len(oids) - 12} more</span>' if len(oids) > 12 else ""
            ids.append(f'<div class="ids"><span class="lbl">{esc(lbl)} ({fint(n)})</span>'
                       f'{chips or "<span class=muted>none in the pool</span>"}{extra}</div>')
        unc = "".join(f"<li>{esc(u)}</li>" for u in r.get("uncertainties") or [])
        out.append(
            f'<details class="card" {"open" if i == 0 else ""} data-geo="{esc(g.get("fips"))}">'
            f'<summary><span class="rk">{fint(g.get("rank"))}</span><span class="nm">{esc(g.get("name"))}</span>'
            f'<span class="meta">composite {f3(g.get("composite"))} &middot; interval {fint(ivl[0])}-{fint(ivl[1])} '
            f'&middot; burden level {esc(g.get("burden_evidence_level"))} &middot; {len(sites)} candidate sites</span>'
            f'</summary><div class="card-b"><p class="step">{esc(r.get("recommended_next_step"))}</p>'
            f'<h4>Candidate sites (first {min(SITES_PER_CARD, len(sites))} of {len(sites)})</h4>'
            f'<ol class="sites">{"".join(site_html)}</ol>'
            + (f'<p class="muted">{more} further candidate sites in the tool output.</p>' if more > 0 else "")
            + f'<h4>Research evidence</h4>{"".join(ids)}<p class="muted">{esc(re_.get("note"))}</p>'
            f'<h4>Uncertainties</h4><ul class="unc">{unc}</ul></div></details>')
    return "".join(out)


def perf_table(c: dict) -> str:
    df = c["perf"].copy()
    cl, ml = c["cond_labels"], c["meas_labels"]
    order = list(dict.fromkeys(df["condition_id"]))
    first = PRIMARY_IDS[0]
    if first in order:
        order = [first] + [o for o in order if o != first]
    chips = "".join(f'<button type="button" class="fchip" data-cond="{esc(o)}" id="perf-{esc(o)}">'
                    f'{esc(cl.get(o, o))}</button>' for o in order)
    body = []
    for cid in order:
        g = df[df["condition_id"] == cid]
        n_known = int((g["performance_status"] == "known").sum())
        body.append(f'<tr class="grp" data-cond="{esc(cid)}"><th colspan="9" class="l">{esc(cl.get(cid, cid))} '
                    f'<span class="muted">({n_known} of {len(g)} measurements with known performance)</span></th></tr>')
        for r in g.to_dict("records"):
            st = r.get("performance_status") or UNKNOWN
            unk = st == UNKNOWN
            auroc = (f'{f3(r.get("auroc"))} ({f3(r.get("auroc_ci_low"))}-{f3(r.get("auroc_ci_high"))})'
                     if _isnum(r.get("auroc")) else "UNKNOWN")
            sens = (f'{f2(r.get("op_sensitivity"))} @ {f2(r.get("op_specificity"))}'
                    if _isnum(r.get("op_sensitivity")) else "UNKNOWN")
            ncc = (f'{fint(r.get("n_cases"))} / {fint(r.get("n_controls"))}' if _isnum(r.get("n_cases"))
                   else "UNKNOWN")
            body.append(
                f'<tr data-cond="{esc(cid)}" data-meas="{esc(r["measurement_id"])}" class="{"unkrow" if unk else ""}">'
                f'<td class="l">{esc(ml.get(r["measurement_id"], r["measurement_id"]))}</td>'
                f'<td class="l"><span class="chip {"warn" if unk else ("ok" if st == "known" else "mid")}">{esc(st)}</span></td>'
                f'<td class="l">{esc(r.get("quality_tier_label") or "UNKNOWN")}</td>'
                f'<td class="l"><code>{esc(r.get("selected_record_id") if isinstance(r.get("selected_record_id"), str) else "UNKNOWN")}</code></td>'
                f'<td class="num">{auroc}</td><td class="num">{sens}</td><td class="num">{ncc}</td>'
                f'<td class="num">{f3(r.get("measurement_evidence"))}</td>'
                f'<td class="l">{esc(r.get("comparator_kind") if isinstance(r.get("comparator_kind"), str) else "UNKNOWN")}</td></tr>')
    head = ("<tr><th class='l'>measurement</th><th class='l'>performance</th><th class='l'>quality tier</th>"
            "<th class='l'>record</th><th>AUROC (95% CI)</th><th>sensitivity @ specificity</th>"
            "<th>cases / controls</th><th>measurement_evidence</th><th class='l'>comparator</th></tr>")
    n_known = int((df["performance_status"] == "known").sum())
    n_partial = int((df["performance_status"] == "partial").sum())
    n_unk = int((df["performance_status"] == UNKNOWN).sum())
    return (f'<p class="lede">{len(df)} condition x measurement pairs: {n_known} known, {n_partial} partial, '
            f'{n_unk} {esc(UNKNOWN)}. A pair whose performance is UNKNOWN or partial gets no evidence_weighted '
            f'rank; it is never scored as 0 or 1.</p>'
            f'<div class="filters" role="group" aria-label="Filter by condition">'
            f'<button type="button" class="fchip on" data-cond="" id="perf-all">All</button>{chips}</div>'
            f'<div class="scroll tall"><table class="data" id="perf"><thead>{head}</thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>')


def map_payload(c: dict) -> dict:
    rows, n = c["rows"], int(c["data"]["n_regions_ranked"])
    vals = {}
    for r in rows.to_dict("records"):
        if _isnum(r.get("rank_equal")):
            vals[r["geo_id"]] = [round(float(r["composite_equal"]), 3), int(r["rank_equal"])]
        else:
            if bool(r.get("burden_incomplete")):
                why = "I"
            elif not bool(r.get("rank_eligible", True)):
                why = "P"
            else:
                why = "N"
            vals[r["geo_id"]] = [None, None, why]
    top_ids = [str(g) for g in c["top"]["geo_id"]]
    sites = []
    for rec in c["data"]["recommendations"]:
        for s in (rec.get("candidate_sites") or [])[:SITES_ON_MAP]:
            if _isnum(s.get("lat")) and _isnum(s.get("lon")):
                sites.append({"lat": round(float(s["lat"]), 4), "lon": round(float(s["lon"]), 4),
                              "name": s.get("facility_name"), "region": rec["geography"]["name"],
                              "precision": s.get("geocode_precision")})
    comp = [v[0] for v in vals.values() if v[0] is not None]
    names = {str(g): n for g, n in zip(c["top"]["geo_id"], c["top"]["geo_name"])}
    return {"values": vals, "top": top_ids, "sites": sites, "n_ranked": n, "names": names,
            "min": min(comp) if comp else None, "max": max(comp) if comp else None,
            "reasons": REASONS}


# ------------------------------------------------------------------------------------------------ page

CSS = r"""
:root{
  --bg:#f3f5f2; --surface:#ffffff; --surface-2:#eef1ec; --ink:#16212b; --muted:#56636e; --line:#d6dcd3;
  --accent:#0d6a70; --accent-ink:#ffffff; --accent-soft:#d5ebe9; --warn:#8a5300; --warn-bg:#f8ecd6;
  --ok:#1d6b3a; --ok-bg:#dcefe1; --mid:#5b4a8a; --mid-bg:#e7e2f4; --null:#8e3b46;
  --ramp-lo:#e4eff0; --ramp-hi:#0a4a66; --unranked:#d3d7cf; --stateline:#ffffff; --countyline:#ffffff;
  --outline:#c2410c; --site:#7c2d12; --focus:#c2410c;
  --serif: Charter, "Bitstream Charter", "Sitka Text", Cambria, Georgia, serif;
  --sans: system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  --mono: ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
  color-scheme: light;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --bg:#0f1519; --surface:#162027; --surface-2:#1b262e; --ink:#e3e9ed; --muted:#9aa8b3; --line:#2a3842;
    --accent:#5cc3c1; --accent-ink:#06282a; --accent-soft:#163a3c; --warn:#e7a948; --warn-bg:#35270f;
    --ok:#7fd49b; --ok-bg:#15301f; --mid:#b9a8ee; --mid-bg:#2a2340; --null:#e58a96;
    --ramp-lo:#1b2d36; --ramp-hi:#8ee0ea; --unranked:#2c353b; --stateline:#0f1519; --countyline:#0f1519;
    --outline:#fb923c; --site:#fdba74; --focus:#fb923c; color-scheme: dark;
  }
}
:root[data-theme="dark"]{
  --bg:#0f1519; --surface:#162027; --surface-2:#1b262e; --ink:#e3e9ed; --muted:#9aa8b3; --line:#2a3842;
  --accent:#5cc3c1; --accent-ink:#06282a; --accent-soft:#163a3c; --warn:#e7a948; --warn-bg:#35270f;
  --ok:#7fd49b; --ok-bg:#15301f; --mid:#b9a8ee; --mid-bg:#2a2340; --null:#e58a96;
  --ramp-lo:#1b2d36; --ramp-hi:#8ee0ea; --unranked:#2c353b; --stateline:#0f1519; --countyline:#0f1519;
  --outline:#fb923c; --site:#fdba74; --focus:#fb923c; color-scheme: dark;
}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 var(--sans); overflow-x:hidden}
.wrap{max-width:1180px; margin:0 auto; padding-inline:16px; padding-block:0 64px}
header.top{padding-block:40px 20px; display:grid; gap:10px}
.eyebrow{font:600 12px/1.2 var(--sans); letter-spacing:.08em; text-transform:uppercase; color:var(--accent)}
h1{font:600 clamp(30px,5vw,46px)/1.08 var(--serif); margin:0; text-wrap:balance; letter-spacing:-.01em}
h2{font:600 clamp(22px,3vw,28px)/1.2 var(--serif); margin:0 0 6px; text-wrap:balance}
h3{font:600 18px/1.3 var(--serif); margin:22px 0 8px}
h4{font:600 12px/1.3 var(--sans); letter-spacing:.06em; text-transform:uppercase; color:var(--muted); margin:16px 0 6px}
p{margin:0 0 10px; max-width:75ch}
.lede{color:var(--muted); max-width:80ch}
.muted{color:var(--muted)}
code{font:12.5px/1.4 var(--mono); background:var(--surface-2); border-radius:4px; padding:1px 5px; overflow-wrap:anywhere}
.stamp{font:12.5px/1.5 var(--mono); color:var(--muted)}
.facts{display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:0; background:var(--surface);
  border:1px solid var(--line); border-radius:8px; overflow:hidden; margin-top:6px}
.fact{background:var(--surface); padding:12px 14px; display:grid; gap:2px; align-content:start}
.fact b{font:600 22px/1.2 var(--serif); font-variant-numeric:tabular-nums}
.fact span{font-size:12.5px; color:var(--muted)}
nav.toc{position:sticky; top:env(safe-area-inset-top,0px); z-index:5; background:var(--bg);
  border-bottom:1px solid var(--line); margin-inline:-16px; padding:8px 16px}
nav.toc ul{list-style:none; margin:0; padding:0; display:flex; gap:6px; overflow-x:auto; scrollbar-width:thin}
nav.toc a{display:block; white-space:nowrap; padding:5px 10px; border-radius:999px; color:var(--ink);
  text-decoration:none; font-size:13px; border:1px solid var(--line); background:var(--surface)}
nav.toc a:hover{border-color:var(--accent)}
a:focus-visible, button:focus-visible, summary:focus-visible{outline:2px solid var(--focus); outline-offset:2px}
section{padding-block:36px 8px; border-bottom:1px solid var(--line); scroll-margin-top:56px}
section:last-of-type{border-bottom:0}
.sec-h{display:grid; gap:4px; margin-bottom:14px}
.scroll{overflow-x:auto; max-width:100%; -webkit-overflow-scrolling:touch; border:1px solid var(--line);
  border-radius:8px; background:var(--surface); margin:10px 0 14px}
.scroll.tall{max-height:640px; overflow-y:auto}
table{border-collapse:collapse; width:100%; font-size:13.5px}
th,td{padding:7px 10px; border-bottom:1px solid var(--line); vertical-align:top; text-align:right}
th{font:600 12px/1.3 var(--sans); color:var(--muted); background:var(--surface-2); position:sticky; top:0; z-index:1;
  white-space:nowrap}
th.l,td.l{text-align:left}
td.num{font-variant-numeric:tabular-nums; white-space:nowrap}
td.strong{font-weight:700}
td.unk{color:var(--warn)}
tbody tr:last-child td{border-bottom:0}
table.data td.l b{font-weight:600}
td.region{min-width:12rem}
table.data td.l code{display:inline-block; margin-top:2px; font-size:11px}
tr.grp th{background:var(--accent-soft); color:var(--ink); font-size:13px; position:static}
tr.unkrow td{color:var(--muted)}
td.pct{position:relative; min-width:74px}
td.pct .bar{position:absolute; left:6px; right:6px; bottom:4px; height:3px; background:var(--line); border-radius:2px}
td.pct .bar::after{content:""; position:absolute; left:0; top:0; bottom:0; width:var(--w); background:var(--accent);
  border-radius:2px}
table.md td, table.md th{text-align:left; min-width:9rem; white-space:normal}
table.md td:first-child{min-width:11rem}
.compact td,.compact th{padding:5px 8px}
svg.ivl{vertical-align:middle; overflow:visible}
svg.ivl .axis{stroke:var(--line); stroke-width:1}
svg.ivl .tick{stroke:var(--muted); stroke-width:1}
svg.ivl .band{stroke:var(--accent); stroke-width:5; stroke-linecap:round; opacity:.45}
svg.ivl .dot{fill:var(--accent); stroke:var(--surface); stroke-width:1.2}
.ivl-t{display:inline-block; margin-left:8px; font:12px var(--mono); color:var(--muted); font-variant-numeric:tabular-nums}
.chip{display:inline-block; padding:1px 8px; border-radius:999px; font:600 11.5px/1.6 var(--sans);
  background:var(--surface-2); color:var(--ink); white-space:nowrap}
.chip.warn{background:var(--warn-bg); color:var(--warn)}
.chip.ok{background:var(--ok-bg); color:var(--ok)}
.chip.mid{background:var(--mid-bg); color:var(--mid)}
.chip.lvl-A{background:var(--ok-bg); color:var(--ok)}
.chip.lvl-B{background:var(--accent-soft); color:var(--accent)}
.chip.lvl-C,.chip.lvl-D{background:var(--warn-bg); color:var(--warn)}
.callout{background:var(--surface); border:1px solid var(--line); border-radius:8px; padding:14px 16px;
  display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:10px 18px; margin:12px 0}
.callout.warn{border-color:var(--warn)}
.callout .note{grid-column:1/-1; color:var(--muted); font-size:13px; margin:0}
.kv{display:grid; gap:2px}
.kv span{font-size:12px; color:var(--muted)}
.kv b{font-variant-numeric:tabular-nums}
.two{display:grid; grid-template-columns:repeat(auto-fit,minmax(min(100%,420px),1fr)); gap:16px}
figure{margin:0}
figcaption{font:600 13px/1.3 var(--sans); margin-top:8px}
.map-wrap{background:var(--surface); border:1px solid var(--line); border-radius:8px; padding:10px; position:relative}
#county-map{width:100%; aspect-ratio:16/10; max-width:100%}
#county-map svg{width:100%; height:100%; display:block}
#county-map .county{stroke:none}
#county-map .statefill{fill:var(--surface-2); stroke:var(--line); stroke-width:.8}
#county-map .state{fill:none; stroke:var(--stateline); stroke-width:1}
#county-map .top{fill:none; stroke:var(--outline); stroke-width:2}
#county-map .site{fill:var(--site); stroke:var(--surface); stroke-width:.8}
#map-msg{padding:40px 10px; text-align:center; color:var(--muted)}
.legend{display:flex; flex-wrap:wrap; align-items:center; gap:8px 18px; font-size:12.5px; color:var(--muted); margin-top:8px}
.ramp{display:inline-block; width:150px; height:10px; border-radius:2px;
  background:linear-gradient(90deg,var(--ramp-lo),var(--ramp-hi)); vertical-align:middle; margin-inline:6px}
.sw{display:inline-block; width:12px; height:12px; border-radius:2px; vertical-align:middle; margin-right:5px}
.sw.un{background:var(--unranked)} .sw.top{border:2px solid var(--outline)} .sw.site{background:var(--site); border-radius:50%}
#tip{position:absolute; pointer-events:none; background:var(--surface); color:var(--ink); border:1px solid var(--line);
  border-radius:6px; padding:6px 9px; font-size:12.5px; box-shadow:0 4px 14px rgba(0,0,0,.12); max-width:260px}
details.card{background:var(--surface); border:1px solid var(--line); border-radius:8px; margin:10px 0}
details.card summary{cursor:pointer; list-style:none; display:grid; grid-template-columns:auto 1fr; gap:2px 12px;
  padding:12px 14px; align-items:baseline}
details.card summary::-webkit-details-marker{display:none}
details.card summary .rk{grid-row:span 2; font:600 26px/1 var(--serif); color:var(--accent); min-width:2ch;
  font-variant-numeric:tabular-nums}
details.card summary .nm{font:600 17px/1.3 var(--serif)}
details.card summary .meta{font-size:12.5px; color:var(--muted)}
details.card[open] summary{border-bottom:1px solid var(--line)}
.card-b{padding:6px 14px 14px}
.step{max-width:90ch}
ol.sites{margin:0; padding-left:20px; display:grid; gap:10px}
.site-h{display:flex; flex-wrap:wrap; gap:4px 10px; align-items:baseline}
ul.reasons{margin:4px 0 0; padding-left:18px; color:var(--muted); font-size:13px}
ul.unc{margin:0; padding-left:18px; font-size:13.5px; display:grid; gap:4px}
.ids{display:flex; flex-wrap:wrap; gap:4px; align-items:center; margin:4px 0}
.ids .lbl{font-size:12.5px; color:var(--muted); margin-right:4px}
.filters{display:flex; flex-wrap:wrap; gap:6px; margin-top:8px}
.fchip{font:13px var(--sans); padding:4px 11px; border-radius:999px; border:1px solid var(--line);
  background:var(--surface); color:var(--ink); cursor:pointer}
.fchip.on{background:var(--accent); color:var(--accent-ink); border-color:var(--accent)}
.md{max-width:100%}
.prose ul{padding-left:20px; display:grid; gap:6px; max-width:95ch}
.prose li p{margin:0}
.prose{font-size:14.5px}
.layers table.md td{min-width:10rem}
footer{padding-block:28px; color:var(--muted); font-size:12.5px}
@media (max-width:640px){
  body{font-size:14.5px}
  header.top{padding-block:28px 14px}
  .fact b{font-size:19px}
  th,td{padding:6px 8px}
}
@media (prefers-reduced-motion: reduce){*{scroll-behavior:auto!important}}
"""

JS = r"""
(function(){
  // performance table filter
  var chips = document.querySelectorAll('.fchip');
  chips.forEach(function(b){ b.addEventListener('click', function(){
    var c = b.getAttribute('data-cond');
    chips.forEach(function(x){ x.classList.toggle('on', x === b); });
    document.querySelectorAll('#perf tbody tr').forEach(function(tr){
      tr.hidden = !!c && tr.getAttribute('data-cond') !== c;
    });
  }); });

  // county map (d3 from cdn.jsdelivr.net; the page works without it)
  var msg = document.getElementById('map-msg');
  var box = document.getElementById('county-map');
  if (typeof d3 === 'undefined') { msg.textContent = 'The map needs one script from cdn.jsdelivr.net (d3), which did not load. Every number is in the tables below.'; return; }
  var M = JSON.parse(document.getElementById('map-data').textContent);
  var geo = JSON.parse(document.getElementById('geo-data').textContent);
  var topSet = new Set(M.top);
  var tip = document.getElementById('tip');
  function tok(n){ return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
  function fmt(x){ return x == null ? 'UNKNOWN' : x.toFixed(3); }
  function draw(){
    box.innerHTML = '';
    var W = box.clientWidth || 900, H = Math.round(W * 10 / 16);
    var proj = d3.geoAlbersUsa().fitExtent([[6, 6], [W - 6, H - 6]], geo.state);
    var path = d3.geoPath(proj);
    var color = d3.scaleLinear().domain([M.min, M.max]).range([tok('--ramp-lo'), tok('--ramp-hi')])
                  .interpolate(d3.interpolateLab).clamp(true);
    var svg = d3.select(box).append('svg').attr('viewBox', '0 0 ' + W + ' ' + H)
      .attr('role', 'img').attr('aria-label', 'County dot map of the equal-weight composite; top regions ringed');
    svg.append('g').selectAll('path').data(geo.state.features).join('path').attr('class', 'statefill').attr('d', path);
    var r = Math.max(1.6, W / 520);
    var pts = geo.dots.map(function(d){ var p = proj([d[1], d[2]]); return p ? {id: d[0], name: M.names[d[0]] || ('County FIPS ' + d[0]), x: p[0], y: p[1]} : null; })
                      .filter(Boolean);
    svg.append('g').selectAll('circle').data(pts).join('circle').attr('class', 'county')
      .attr('cx', function(d){ return d.x; }).attr('cy', function(d){ return d.y; }).attr('r', r)
      .attr('fill', function(d){ var v = M.values[d.id]; return v && v[0] != null ? color(v[0]) : tok('--unranked'); })
      .on('pointermove', function(ev, d){
        var v = M.values[d.id] || [null, null, 'not in the ranking universe'];
        var why = M.reasons[v[2]] || v[2];
        tip.innerHTML = '<b>' + d.name + '</b><br>' + (v[0] != null ? ('composite ' + fmt(v[0]) + ' &middot; rank ' + v[1] + ' of ' + M.n_ranked) : why);
        var b = box.parentNode.getBoundingClientRect();
        tip.style.left = Math.min(ev.clientX - b.left + 12, b.width - 270) + 'px';
        tip.style.top = (ev.clientY - b.top + 12) + 'px'; tip.hidden = false;
      })
      .on('pointerleave', function(){ tip.hidden = true; });
    svg.append('g').selectAll('circle').data(pts.filter(function(d){ return topSet.has(d.id); })).join('circle')
      .attr('class', 'top').attr('cx', function(d){ return d.x; }).attr('cy', function(d){ return d.y; }).attr('r', r * 3.2)
      .append('title').text(function(d){ return d.name + ' (top ' + M.top.length + ')'; });
    svg.append('g').selectAll('circle').data(M.sites.filter(function(s){ return proj([s.lon, s.lat]); }))
      .join('circle').attr('class', 'site').attr('r', Math.max(2.2, W / 420))
      .attr('cx', function(s){ return proj([s.lon, s.lat])[0]; }).attr('cy', function(s){ return proj([s.lon, s.lat])[1]; })
      .append('title').text(function(s){ return s.name + ' (' + s.region + '; location: ' + s.precision + ')'; });
  }
  draw();
  var t; window.addEventListener('resize', function(){ clearTimeout(t); t = setTimeout(draw, 150); });
  try { window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', draw); } catch (e) {}
  try { new MutationObserver(draw).observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']}); } catch (e) {}
})();
"""


def render(c: dict) -> str:
    d, top, rep, rd = c["data"], c["top"], c["report"], c["readme"]
    n_ranked = int(d["n_regions_ranked"])
    mp = map_payload(c)
    geo = map_geometry()
    top1 = top.iloc[0]
    excl = d.get("excluded_incomplete_burden") or []
    snapshot = {
        "primary_query": PRIMARY, "n_regions_ranked": n_ranked,
        "deployment_opportunities_version": d.get("deployment_opportunities_version"),
        "top10_equal": [{"geo_id": r["geo_id"], "geo_name": r["geo_name"], "rank_equal": r["rank_equal"],
                         "composite_equal": r["composite_equal"], "rank_mc_p05": r["rank_mc_p05"],
                         "rank_mc_p95": r["rank_mc_p95"], "p_top10_mc": r["p_top10_mc"],
                         "rank_evidence_weighted": r["rank_evidence_weighted"],
                         **{col: r[col] for col, _ in COMPONENTS}} for r in top.to_dict("records")],
        "top10_evidence_weighted": [{"geo_id": r["geo_id"], "rank_evidence_weighted": r["rank_evidence_weighted"],
                                     "composite_evidence_weighted": r["composite_evidence_weighted"],
                                     "rank_equal": r["rank_equal"]} for r in c["top_ew"].to_dict("records")],
        "git": c["git"], "built_at": c["built_at"],
    }
    facts = [
        (fint(n_ranked), "counties ranked (population at least the minimum, complete burden)"),
        (fint(len(excl)), "counties not ranked: incomplete burden"),
        (esc(top1["geo_name"]), f"rank 1 under equal weights (interval {fint(top1['rank_mc_p05'])}-{fint(top1['rank_mc_p95'])})"),
        (fint(c["n_sources"]), "sources in SOURCE_REGISTRY.yaml"),
        (fint(c["n_tables"]), "processed tables"),
    ]
    facts_html = "".join(f'<div class="fact"><b>{v}</b><span>{esc(k)}</span></div>' for v, k in facts)
    lo, hi = mp["min"], mp["max"]
    framing = esc(d.get("framing"))
    meth = esc(d.get("method"))
    unc = "".join(f"<li>{esc(u)}</li>" for u in d.get("general_uncertainties") or [])
    toc = [("overview", "Overview"), ("changes", "What changed"), ("query", "Primary query"), ("map", "Map"), ("cards", "Top 10 evidence"),
           ("performance", "Measurement performance"), ("works", "What works"), ("validation", "Validation"),
           ("limits", "Limitations")]
    toc_html = "".join(f'<li><a href="#{a}">{esc(b)}</a></li>' for a, b in toc)
    h = rep["headings"]
    return f"""<!doctype html>
<meta charset="utf-8">
<title>{esc(TITLE)}</title>
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="description" content="Read-only snapshot of the Measure It to Cure It public-data engine: primary deployment query, county map, evidence cards, measurement performance and validation.">
<style>{CSS}</style>
<div class="wrap">
<header class="top">
  <div class="eyebrow">Measure It to Cure It &middot; public-data alpha &middot; read-only snapshot</div>
  <h1>{esc(TITLE)}</h1>
  <p class="lede">Where could an objective measurement for Long COVID or ME/CFS be piloted, what supports it, and what
  is still unknown. Every figure below is read from the engine's processed tables and results files when this page
  was built; nothing on it is live, and nothing describes a person.</p>
  <div class="stamp">built {esc(c["built_at"])} &middot; code {esc(c["git"].get("describe") or c["git"].get("commit"))}
  &middot; deployment_opportunities version {esc(d.get("deployment_opportunities_version"))}</div>
  <div class="facts">{facts_html}</div>
</header>
<nav class="toc" aria-label="Sections"><ul>{toc_html}</ul></nav>

<section id="overview">
  <div class="sec-h"><div class="eyebrow">Overview</div><h2>What the engine is, and is not</h2></div>
  <div class="prose">{md_to_html(rd["intro"])}</div>
  <h3>The four data layers and the guardrails</h3>
  <div class="prose layers">{md_to_html(rd["layers"])}</div>
</section>

<section id="changes">
  <div class="sec-h"><div class="eyebrow">Update 2026-10-07</div><h2>What changed since the 2026-09-24 build</h2></div>
  <div class="prose">{md_to_html(re.sub(r"Sections below describe the build of 2026-09-24 and are not\s+rewritten\.\s*",
                                         "The ranking, map, cards and performance sections on this page are rebuilt "
                                         "from the current tables; the report sections 13, 14 and 19 further down "
                                         "are from the 2026-09-24 build. ", rep.get("update") or ""))}</div>
  {changes_html()}
</section>

<section id="query">
  <div class="sec-h"><div class="eyebrow">Primary deployment query</div>
  <h2>{esc(d["condition_resolution"].get("query") or PRIMARY["condition"])} &times; {esc(d["measurement_resolution"].get("label"))}</h2>
  <p class="lede">{framing}</p></div>
  <p class="muted">Universe: {esc(d.get("ranking_universe"))}. Basis: {esc(d.get("ranking_basis"))}. Weight set
  <b>{esc(d.get("weight_set"))}</b> ({esc(", ".join(f"{k} {v:g}" for k, v in (d.get("weights") or {}).items()))}).
  Method: {meth}.</p>
  <h3>The measurement's own evidence</h3>
  {perf_callout(c)}
  <h3>Top {len(top)} regions under equal weights</h3>
  <p class="muted">Component columns are percentiles (0-1, bar = value). The rank interval is the 5th-95th percentile
  of the rank over Monte Carlo draws of burden and weights, drawn on a log scale from 1 to {fint(n_ranked)}.</p>
  {ranked_table(c)}
  <h3>evidence_weighted compared with equal</h3>
  {comparison(c)}
  <h3>Rank of the top {len(top)} under every weight set</h3>
  <p class="muted">The same regions under the {len(c["weight_sets"])} weight sets of configs/scoring.yaml.</p>
  {weightset_table(c)}
  <h3>Uncertainties that apply to every region</h3>
  <ul class="unc">{unc}</ul>
</section>

<section id="map">
  <div class="sec-h"><div class="eyebrow">County map</div><h2>Equal-weight composite by county</h2>
  <p class="lede">The colour runs from low to high composite (legend under the map). The top {len(top)} are ringed; dots are their first
  {SITES_ON_MAP} candidate sites (ZIP/ZCTA or city centroids). Grey counties are not ranked; hover or tap for the reason.</p></div>
  <div class="map-wrap"><div id="county-map"><div id="map-msg">Loading the map&hellip;</div></div><div id="tip" hidden></div>
  <div class="legend"><span>composite {f3(lo)}<span class="ramp"></span>{f3(hi)}</span>
  <span><span class="sw un"></span>not ranked</span><span><span class="sw top"></span>top {len(top)}</span>
  <span><span class="sw site"></span>candidate site</span></div></div>
  <p class="muted">Each dot is one 2024 county at its Census internal point; state outlines are simplified Census 2024
  cartographic boundaries, for display only. Alaska and Hawaii are drawn as insets.</p>
</section>

<section id="cards">
  <div class="sec-h"><div class="eyebrow">Evidence cards</div><h2>The top {len(top)} regions, one card each</h2>
  <p class="lede">Candidate sites are facilities with characteristics suggesting they may be viable implementation or
  study partners, not a quality ranking. Trials and grants are research-activity signals, not evidence that the
  measurement works.</p></div>
  {cards(c)}
</section>

<section id="performance">
  <div class="sec-h"><div class="eyebrow">Metric link</div><h2>Measurement performance</h2>
  <p class="muted">data/processed/measurement_performance: the scored performance record per condition and
  measurement, with its quality tier and comparator.</p></div>
  {perf_table(c)}
</section>

<section id="works">
  <div class="sec-h"><div class="eyebrow">results/FINAL_REPORT.md &sect;14</div><h2>{esc(h["s14"])}</h2></div>
  <div class="prose">{md_to_html(rep["s14_intro"])}</div>
  <h3>{esc(h["s14_2"])}</h3><div class="prose">{md_to_html(rep["s14_2"])}</div>
  <h3>{esc(h["s14_3"])}</h3><div class="prose">{md_to_html(rep["s14_3"])}</div>
  <h3>{esc(h["s14_4"])}</h3><div class="prose">{md_to_html(rep["s14_4"])}</div>
  <h3>{esc(h["s14_5"])}</h3><div class="prose">{md_to_html(rep["s14_5"])}</div>
</section>

<section id="validation">
  <div class="sec-h"><div class="eyebrow">results/FINAL_REPORT.md &sect;13</div><h2>{esc(h["s13"])}</h2></div>
  <p class="muted"><b>Describes the 2026-09-24 build.</b> Current ranking stability and the clinic-shuffle result are
  under What changed (the report's section 13 has not been rewritten).</p>
  <div class="prose">{md_to_html(rep["s13"])}</div>
</section>

<section id="limits">
  <div class="sec-h"><div class="eyebrow">results/FINAL_REPORT.md &sect;19</div><h2>{esc(h["s19"])}</h2></div>
  <p class="muted"><b>Written for the 2026-09-24 build.</b> Two limitations have since changed: Long COVID county
  burden is now a modelled BRFSS 2023 small-area estimate (its within-state ordering is modelled, not observed),
  not the inherited state value, and clinic capacity is measurement-specific Medicare billing where a billing code
  exists. See What changed.</p>
  <div class="prose">{md_to_html(rep["s19"])}</div>
</section>

<footer>
  Built by measure_it.export.static_site from data/processed/deployment_opportunities.parquet,
  data/processed/measurement_performance.parquet, the tool rank_deployment_opportunities, results/FINAL_REPORT.md,
  README.md and Census 2024 geography. Every ranked region is a candidate deployment opportunity for pilot
  evaluation, not a validated diagnostic pathway.
</footer>
</div>
<script type="application/json" id="snapshot-data">{_script_json(snapshot)}</script>
<script type="application/json" id="map-data">{_script_json(mp)}</script>
<script type="application/json" id="geo-data">{_script_json(geo)}</script>
<script src="{D3_URL}"></script>
<script>{JS}</script>
"""


def build_static_site(out: Path | None = None) -> dict:
    out = Path(out or OUT)
    out.parent.mkdir(parents=True, exist_ok=True)
    page = render(collect())
    raw = page.encode("utf-8")
    if len(raw) > MAX_BYTES:
        raise RuntimeError(f"static snapshot is {len(raw) / 1e6:.1f} MB, above the {MAX_BYTES / 1e6:.0f} MB limit")
    out.write_bytes(raw)
    # Publishable variant: the artifact host wraps pages in its own doctype/head, so drop those lines.
    art = out.with_name("artifact.html")
    body = "\n".join(l for l in page.splitlines()
                     if not l.startswith(("<!doctype", '<meta charset', '<meta name="viewport"')))
    art.write_text(body, encoding="utf-8")
    return {"path": str(out), "bytes": len(raw), "artifact_path": str(art), "artifact_bytes": len(body.encode())}
