"""Plotly map builders (geo subplot, Albers USA; no tile server). Colours: one-hue blue ramp for every choropleth
(magnitude), point overlays in validated categorical slots with distinct symbols (identity is never colour alone),
specialist density as neutral hollow circles sized by value."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from measure_it.config import UNKNOWN

from .data import subset_geojson

BLUES = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf",
         "#1c5cab", "#184f95", "#104281", "#0d366b"]
COLORSCALE = [[i / (len(BLUES) - 1), c] for i, c in enumerate(BLUES)]
NODATA = "#dcdad3"
EDGE = "rgba(255,255,255,0.55)"
STATE_EDGE = "rgba(60,60,60,0.55)"
# categorical slots validated all-pairs (orange, aqua, violet; dataviz validator) + neutral for magnitude bubbles
OVERLAY = {"hrsa": {"color": "#eb6834", "symbol": "circle", "name": "HRSA health-center sites"},
           "trials": {"color": "#1baf7a", "symbol": "diamond", "name": "Trial sites (city centroids)"},
           "nih": {"color": "#4a3aa7", "symbol": "square", "name": "NIH-funded organisations"},
           "specialists": {"color": "#52514e", "symbol": "circle-open", "name": "Specialist density"},
           "sites": {"color": "#eb6834", "symbol": "diamond", "name": "Candidate sites"},
           "clinics": {"color": "#4a3aa7", "symbol": "square", "name": "Candidate facilities"}}


def _layout(fig: go.Figure, height: int, legend: bool = True) -> go.Figure:
    fig.update_geos(scope="usa", projection_type="albers usa", showland=True, landcolor="rgba(0,0,0,0)",
                    showlakes=False, showcoastlines=False, showframe=False, showsubunits=False, showcountries=False,
                    bgcolor="rgba(0,0,0,0)")
    fig.update_layout(height=height, margin=dict(l=0, r=0, t=30, b=0), paper_bgcolor="rgba(0,0,0,0)",
                      showlegend=legend, legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="left", x=0,
                                                     bgcolor="rgba(0,0,0,0)"),
                      hoverlabel=dict(align="left"), dragmode=False)
    return fig


def choropleth_map(fc: dict, df: pd.DataFrame, value_col: str, hover: list[str], colorbar_title: str,
                   nodata_hover: list[str] | None = None, state_fc: dict | None = None,  # columns for grey hover
                   zmin: float | None = None, zmax: float | None = None, height: int = 600,
                   highlight: list[str] | None = None, reversescale: bool = False) -> go.Figure:
    """One data trace for geographies with a value; one grey trace (with its reason on hover) for those without."""
    fig = go.Figure()
    vals = pd.to_numeric(df[value_col], errors="coerce") if len(df) else pd.Series(dtype=float)
    has = vals.notna() & np.isfinite(vals)
    d_ok, d_no = df[has], df[~has]
    known_ids = set(df["geo_id"]) if len(df) else set()
    outside = [f["id"] for f in fc["features"] if f["id"] not in known_ids]
    no_ids = list(d_no["geo_id"]) + outside
    if no_ids:
        text = hover_text(d_no, nodata_hover) if (nodata_hover and len(d_no)) else [UNKNOWN] * len(d_no)
        text = text + ["not in this layer"] * len(outside)
        fig.add_trace(go.Choropleth(
            geojson=subset_geojson(fc, no_ids), locations=no_ids, z=[0] * len(no_ids),
            colorscale=[[0, NODATA], [1, NODATA]], showscale=False, marker_line_color=EDGE, marker_line_width=0.3,
            hovertext=text, hoverinfo="text", name="no value (UNKNOWN / NOT AVAILABLE)", showlegend=True))
    if len(d_ok):
        fig.add_trace(go.Choropleth(
            geojson=subset_geojson(fc, d_ok["geo_id"]), locations=list(d_ok["geo_id"]), z=vals[has].tolist(),
            colorscale=COLORSCALE, reversescale=reversescale, zmin=zmin, zmax=zmax, marker_line_color=EDGE,
            marker_line_width=0.3, hovertext=hover_text(d_ok, hover), hoverinfo="text", name=colorbar_title,
            showlegend=False,
            colorbar=dict(title=dict(text=colorbar_title, side="right"), thickness=12, len=0.62, x=1.0, y=0.45)))
    if state_fc is not None:
        fig.add_trace(go.Choropleth(
            geojson=state_fc, locations=[f["id"] for f in state_fc["features"]], z=[0] * len(state_fc["features"]),
            colorscale=[[0, "rgba(0,0,0,0)"], [1, "rgba(0,0,0,0)"]], showscale=False, marker_line_color=STATE_EDGE,
            marker_line_width=0.7, hoverinfo="skip", showlegend=False, name="state borders"))
    if highlight:
        hl = [h for h in highlight if h in set(d_ok["geo_id"]) | set(no_ids)]
        if hl:
            fig.add_trace(go.Choropleth(
                geojson=subset_geojson(fc, hl), locations=hl, z=[0] * len(hl),
                colorscale=[[0, "rgba(0,0,0,0)"], [1, "rgba(0,0,0,0)"]], showscale=False,
                marker_line_color="#e34948", marker_line_width=2.2, hoverinfo="skip", name="top ranked (red outline)",
                showlegend=True))
    return _layout(fig, height)


def hover_text(df: pd.DataFrame, cols: list[str]) -> list[str]:
    rows = []
    for r in df[cols].itertuples(index=False):
        parts = []
        for c, v in zip(cols, r):
            if c == cols[0]:
                parts.append(f"<b>{v}</b>")
            elif isinstance(v, str) and v.startswith("<br>"):
                parts.append(v[4:])
            else:
                parts.append(f"{c}: {_short(v)}")
        rows.append("<br>".join(p for p in parts if p))
    return rows


def _short(v) -> str:
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return UNKNOWN
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (float, np.floating)):
        return f"{v:.3g}" if abs(v) < 1000 else f"{v:,.0f}"
    return str(v)


def add_points(fig: go.Figure, df: pd.DataFrame, kind: str, hover: list[str], size_col: str | None = None,
               max_size: float = 22.0, base_size: float = 5.0, name: str | None = None) -> go.Figure:
    if df is None or not len(df):
        return fig
    sty = OVERLAY[kind]
    if size_col:
        v = pd.to_numeric(df[size_col], errors="coerce").fillna(0).to_numpy(float)
        top = np.nanpercentile(v[v > 0], 98) if (v > 0).any() else 1.0
        size = base_size + (max_size - base_size) * np.sqrt(np.clip(v, 0, top) / max(top, 1e-9))
    else:
        size = np.full(len(df), base_size)
    line = dict(width=1.6, color=sty["color"]) if sty["symbol"].endswith("open") else dict(width=0.8, color="white")
    fig.add_trace(go.Scattergeo(
        lat=df["lat"], lon=df["lon"], mode="markers", name=name or sty["name"],
        marker=dict(size=size, color=sty["color"], symbol=sty["symbol"], line=line,
                    opacity=0.85 if not sty["symbol"].endswith("open") else 0.9),
        hovertext=hover_text(df, hover), hoverinfo="text"))
    return fig


_BBOX_CACHE: dict[int, tuple[dict, dict]] = {}


def _coords_bbox(coords) -> tuple[float, float, float, float]:
    a = np.asarray([c for c in _flatten(coords)], dtype=float)
    return float(a[:, 0].min()), float(a[:, 1].min()), float(a[:, 0].max()), float(a[:, 1].max())


def _flatten(coords):
    if coords and isinstance(coords[0], (int, float)):
        yield coords
        return
    for c in coords:
        yield from _flatten(c)


def _feature_bboxes(fc: dict) -> dict[str, tuple[float, float, float, float]]:
    """id -> (minlon, minlat, maxlon, maxlat), memoised per FeatureCollection object (the cached, never-mutated
    boundaries of components.data.geojson; the collection is kept referenced so its id cannot be reused)."""
    hit = _BBOX_CACHE.get(id(fc))
    if hit is not None and hit[0] is fc:
        return hit[1]
    out = {f["id"]: _coords_bbox(f["geometry"]["coordinates"]) for f in fc["features"] if f.get("geometry")}
    _BBOX_CACHE[id(fc)] = (fc, out)
    return out


def local_map(fc_area: dict, points: list[tuple[pd.DataFrame, str, list[str]]], height: int = 480,
              area_hover: str = "", context_fc: dict | None = None) -> go.Figure:
    """Zoomed map of one geography (its polygon) with point layers. The view is the bounding box of the polygon and
    the points (padded); neighbouring units come from `context_fc` (the same simplified Census boundaries, drawn in
    grey), so no basemap coastline at another resolution is drawn over them. Falls back to fitting the polygon when it
    crosses the antimeridian (Alaska)."""
    fig = go.Figure()
    ids = [f["id"] for f in fc_area["features"]]
    boxes = [_coords_bbox(f["geometry"]["coordinates"]) for f in fc_area["features"] if f.get("geometry")]
    lons = [b[0] for b in boxes] + [b[2] for b in boxes]
    lats = [b[1] for b in boxes] + [b[3] for b in boxes]
    for df, _, _ in points:
        if df is not None and len(df) and {"lat", "lon"} <= set(df.columns):
            lons += pd.to_numeric(df["lon"], errors="coerce").dropna().tolist()
            lats += pd.to_numeric(df["lat"], errors="coerce").dropna().tolist()
    view = None
    if lons and lats and max(lons) - min(lons) < 60:
        pad_x = max((max(lons) - min(lons)) * 0.06, 0.05)
        pad_y = max((max(lats) - min(lats)) * 0.06, 0.05)
        view = (min(lons) - pad_x, min(lats) - pad_y, max(lons) + pad_x, max(lats) + pad_y)
    if context_fc is not None and view is not None:
        bb, own = _feature_bboxes(context_fc), set(ids)
        near = [i for i, (x0, y0, x1, y1) in bb.items()
                if i not in own and x1 >= view[0] and x0 <= view[2] and y1 >= view[1] and y0 <= view[3]]
        if near:
            fig.add_trace(go.Choropleth(
                geojson=subset_geojson(context_fc, near), locations=near, z=[0] * len(near),
                colorscale=[[0, "rgba(128,128,128,0.07)"], [1, "rgba(128,128,128,0.07)"]], showscale=False,
                marker_line_color="#b4b2a9", marker_line_width=0.8, hoverinfo="skip", name="neighbouring units",
                showlegend=False))
    fig.add_trace(go.Choropleth(geojson=fc_area, locations=ids, z=[0] * len(ids),
                                colorscale=[[0, "rgba(42,120,214,0.10)"], [1, "rgba(42,120,214,0.10)"]],
                                showscale=False, marker_line_color="#2a78d6", marker_line_width=1.6,
                                hovertext=[area_hover] * len(ids), hoverinfo="text", name="geography",
                                showlegend=False))
    for df, kind, hover in points:
        add_points(fig, df, kind, hover, base_size=9)
    common = dict(projection_type="mercator", showland=False, showlakes=False, showocean=False, showcoastlines=False,
                  showsubunits=False, showcountries=False, showframe=False, bgcolor="rgba(0,0,0,0)")
    if view is not None:
        fig.update_geos(fitbounds=False, lonaxis_range=[view[0], view[2]], lataxis_range=[view[1], view[3]],
                        **common)
    else:
        fig.update_geos(fitbounds="locations", **common)
    fig.update_layout(height=height, margin=dict(l=0, r=0, t=30, b=0), paper_bgcolor="rgba(0,0,0,0)",
                      legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="left", x=0),
                      hoverlabel=dict(align="left"))
    return fig


def ranking_map(rank_data: dict, weight_set: str, level: str, cat: dict, stamp: float, log=None,
                height: int = 540) -> tuple[go.Figure, str]:
    """Choropleth for a rank_deployment_opportunities answer: every ranked region's composite (precomputed
    combinations, from deployment_opportunities) with the returned top regions outlined; for combinations scored on the
    fly only the returned regions are drawn (coloured by rank). Returns (figure, caption)."""
    from . import data as D
    recs = rank_data.get("recommendations") or []
    top_ids = [r["geography"]["fips"] for r in recs]
    cid = (rank_data.get("condition_resolution") or {}).get("condition_id")
    mid = (rank_data.get("measurement_resolution") or {}).get("measurement_id")
    rank_col = rank_data.get("rank_column") or f"rank_{weight_set}"
    fc = D.geojson(level, 0.02 if level == "county" else None)
    state_fc = D.geojson("state") if level == "county" else None
    if cid in cat["scored_conditions"] and mid in cat["scored_measurements"] and rank_col == f"rank_{weight_set}":
        u = D.opportunity_rows(cid, mid, level, stamp)
        if log is not None:
            log.table("deployment_opportunities", f"composite_{weight_set} of every ranked region (map colour)")
        u = u[u[rank_col].notna()]
        n = len(u)
        uni = pd.DataFrame({"geo_id": u["geo_id"], "region": u["geo_name"], "value": u[f"composite_{weight_set}"],
                            "rank": [f"{x:.0f} of {n}" for x in u[rank_col]], "composite": u[f"composite_{weight_set}"],
                            "burden level": u["burden_evidence_level"].fillna("D") + np.where(
                                u["burden_inherited"].fillna(False).astype(bool), " (inherited state value)", "")})
        fig = choropleth_map(fc, uni, "value", ["region", "rank", "composite", "burden level"],
                             f"composite ({weight_set})", state_fc=state_fc, highlight=top_ids, height=height)
        note = (f"Blue = composite of every ranked region under the '{weight_set}' weights (darker = higher); red "
                f"outline = the {len(top_ids)} returned regions; grey = not ranked (small population or incomplete "
                "burden)")
        return fig, note
    rows = [{"geo_id": r["geography"]["fips"], "region": r["geography"]["name"],
             "value": r["geography"].get("rank"), "rank": r["geography"].get("rank"),
             "composite": r["geography"].get("composite"),
             "burden level": r["geography"].get("burden_evidence_level")} for r in recs]
    fig = choropleth_map(fc, pd.DataFrame(rows, columns=["geo_id", "region", "value", "rank", "composite",
                                                         "burden level"]),
                         "value", ["region", "rank", "composite", "burden level"], "rank (1 = top)",
                         state_fc=state_fc, reversescale=True, highlight=top_ids, height=height)
    return fig, (f"Blue = the {len(top_ids)} returned regions (darker = higher rank); this combination is scored on the "
                 "fly, so only the returned regions are drawn")
