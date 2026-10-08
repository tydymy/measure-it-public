"""Small Plotly charts used across pages (single axis, thin marks, categorical slots in fixed order)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SEQ = "#2a78d6"
GRID = "rgba(128,128,128,0.18)"


def _base(fig: go.Figure, height: int, xtitle: str = "", ytitle: str = "") -> go.Figure:
    fig.update_layout(height=height, margin=dict(l=10, r=10, t=36, b=10), paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="rgba(0,0,0,0)", legend=dict(orientation="h", yanchor="bottom", y=1.0, x=0),
                      hoverlabel=dict(align="left"), bargap=0.25)
    fig.update_xaxes(title=xtitle, gridcolor=GRID, zerolinecolor="rgba(128,128,128,0.5)")
    fig.update_yaxes(title=ytitle, gridcolor=GRID)
    return fig


def forest(df: pd.DataFrame, group_col: str, label_col: str, effect: str = "effect", lo: str = "ci_low",
           hi: str = "ci_high", xtitle: str = "effect (group difference)", hover_cols: list[str] | None = None,
           height: int | None = None) -> go.Figure:
    """Effect sizes with 95% CIs, one colour per dataset/phenotype (fixed slot order), zero line."""
    fig = go.Figure()
    groups = list(dict.fromkeys(df[group_col]))
    for i, g in enumerate(groups):
        d = df[df[group_col] == g]
        hover = ["<br>".join(f"{c}: {r[c]}" for c in (hover_cols or []) if c in r) for _, r in d.iterrows()]
        fig.add_trace(go.Scatter(
            x=d[effect], y=d[label_col], mode="markers", name=str(g), marker=dict(size=9, color=SLOTS[i % len(SLOTS)],
                                                                                  line=dict(width=1, color="white")),
            error_x=dict(type="data", symmetric=False, array=(d[hi] - d[effect]).tolist(),
                         arrayminus=(d[effect] - d[lo]).tolist(), thickness=1.6, width=0,
                         color=SLOTS[i % len(SLOTS)]),
            hovertext=hover, hoverinfo="text+x"))
    fig.add_vline(x=0, line_width=1, line_color="rgba(128,128,128,0.8)")
    fig.update_yaxes(autorange="reversed", type="category")
    _base(fig, height or max(240, 28 * len(df) + 60 + 22 * len(groups)), xtitle)
    fig.update_layout(showlegend=True, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))
    return fig


def hbar(df: pd.DataFrame, x: str, y: str, xtitle: str = "", hover_cols: list[str] | None = None,
         height: int | None = None, color: str = SEQ, text: str | None = None) -> go.Figure:
    hover = ["<br>".join(f"{c}: {r[c]}" for c in (hover_cols or []) if c in r) for _, r in df.iterrows()]
    fig = go.Figure(go.Bar(x=df[x], y=df[y], orientation="h", marker=dict(color=color, line=dict(width=0)),
                           hovertext=hover or None, hoverinfo="text+x" if hover_cols else "x+y",
                           text=df[text] if text else None, textposition="outside" if text else None,
                           cliponaxis=False))
    fig.update_yaxes(autorange="reversed", type="category")
    return _base(fig, height or max(220, 24 * len(df) + 80), xtitle)


def error_bars(df: pd.DataFrame, x: str, y: str, lo: str, hi: str, xtitle: str = "", height: int | None = None,
               hover_cols: list[str] | None = None) -> go.Figure:
    hover = ["<br>".join(f"{c}: {r[c]}" for c in (hover_cols or []) if c in r) for _, r in df.iterrows()]
    fig = go.Figure(go.Scatter(
        x=df[x], y=df[y], mode="markers", marker=dict(size=9, color=SEQ, line=dict(width=1, color="white")),
        error_x=dict(type="data", symmetric=False, array=(df[hi] - df[x]).tolist(),
                     arrayminus=(df[x] - df[lo]).tolist(), thickness=1.6, width=0, color=SEQ),
        hovertext=hover or None, hoverinfo="text+x" if hover_cols else "x+y", showlegend=False))
    fig.update_yaxes(autorange="reversed", type="category")
    return _base(fig, height or max(220, 26 * len(df) + 80), xtitle)


def stacked_contributions(df: pd.DataFrame, label_col: str, parts: list[tuple[str, str]], xtitle: str,
                          height: int | None = None) -> go.Figure:
    """Horizontal stacked bars: one segment per (column, legend name), fixed slot order, 2px surface gaps."""
    fig = go.Figure()
    for i, (col, name) in enumerate(parts):
        fig.add_trace(go.Bar(x=df[col], y=df[label_col], orientation="h", name=name,
                             marker=dict(color=SLOTS[i % len(SLOTS)], line=dict(width=1.5, color="white")),
                             hovertemplate=f"{name}: %{{x:.3f}}<extra>%{{y}}</extra>"))
    fig.update_layout(barmode="stack")
    fig.update_yaxes(autorange="reversed", type="category")
    _base(fig, height or max(260, 30 * len(df) + 100), xtitle)
    fig.update_layout(legend=dict(orientation="h", traceorder="normal", yanchor="bottom", y=1.0, xanchor="left", x=0,
                                  entrywidth=0, font=dict(size=11)))
    return fig


def finite(v) -> float:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return np.nan
    return f if np.isfinite(f) else np.nan
