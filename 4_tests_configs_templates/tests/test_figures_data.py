"""Tests for the SPEC data Figures 3-9 (measure_it.figures.data_figures and its fig_* / data_* helpers).

Unit tests need nothing but the code. Data tests (marked `data`) render every figure from the processed and result
tables and check (1) layout: no text outside the canvas or overlapping, fonts >= 5.5 pt, double-column width,
height <= 247 mm; (2) product language and guardrail labels (proxy, Rx-defined, NOT Long COVID, inherited
burden, city centroids, UNKNOWN / NOT AVAILABLE); (3) that what the figures show equals what the tables and the
measure_it query functions return; (4) the published files (PNG 300 dpi + SVG, interactive maps, CAPTIONS.md block).
"""
from __future__ import annotations

import hashlib
import re
import struct
import warnings

import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from measure_it import pipeline as P
from measure_it.config import FIGURES, MAPS, RESULTS, TABLES, UNKNOWN
from measure_it.figures import data_checks as DC
from measure_it.figures import data_common as DCm
from measure_it.figures import style as S

STEMS = {
    "fig03": "fig03_wearable_phenotype_differences", "fig04": "fig04_clinical_vs_wearable_models",
    "fig05": "fig05_molecular_measurable_physiology", "fig06": "fig06_long_covid_burden_unmet_need",
    "fig07": "fig07_infrastructure_overlay", "fig08": "fig08_deployment_opportunities",
    "fig09": "fig09_san_diego_drilldown",
}
MAP_STEMS = ["fig06_long_covid_burden_unmet_need", "fig07_infrastructure_overlay", "fig08_deployment_opportunities",
             "fig09_san_diego_drilldown"]


# ================================================================================================= unit tests
def _fig(w=3.0, h=2.0):
    f = Figure(figsize=(w, h), dpi=100)
    FigureCanvasAgg(f)
    return f


def test_layout_checker_finds_overlap_and_outside_text():
    f = _fig()
    f.text(0.1, 0.5, "overlapping label one", fontsize=8)
    f.text(0.12, 0.5, "overlapping label two", fontsize=8)
    f.text(0.95, 0.5, "runs off the right edge", fontsize=8)
    probs = DC.layout_problems(f)
    assert any(p.startswith("text overlap") for p in probs)
    assert any(p.startswith("outside figure") for p in probs)


def test_layout_checker_ignores_undrawn_ticks_and_axis_off():
    f = _fig()
    ax = f.add_axes((0.2, 0.2, 0.6, 0.6))
    ax.set_xticks([-5, 0, 0.5, 1, 7], ["-5", "0", "0.5", "1", "7"], fontsize=6)
    ax.set_xlim(0, 1)                            # -5 and 7 are now outside the view: never drawn
    ax.set_yticks([])
    off = f.add_axes((0, 0, 1, 1))
    off.set_axis_off()                           # full-figure canvas axes: its ticks are not drawn either
    assert DC.layout_problems(f) == []


def test_formatting_helpers():
    assert DCm.num(-0.5, 2) == "−0.50"
    assert DCm.ci_text(0.1, -0.2, 0.45, 2, True) == "+0.10 (−0.20 to +0.45)"
    assert DCm.p_text(1 / 201, floor=1 / 201) == "p ≤ 0.005"
    assert DCm.p_text(0.0001) == "p < 0.001"
    assert DCm.p_text(0.53) == "p = 0.53"
    assert DCm.intc(2340.0) == "2,340"


def test_quantile_classes_and_labels():
    from measure_it.figures.fig_geography import _class_labels, _quantile_bins
    v = pd.Series([0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 50.0])
    e = _quantile_bins(v, [0.25, 0.5, 0.75], zero_class=True)
    assert e[0] == 0 and e[1] < 1e-6 and e[-1] > 50 and e == sorted(e)
    labs = _class_labels(e, True)
    assert labs[0] == "0" and labs[1].startswith(">0") and labs[-1].startswith("≥")
    e2 = _quantile_bins(pd.Series(np.arange(100.0)), [0.0, 0.2, 0.4, 0.6, 0.8], zero_class=False)
    assert len(e2) == 6 and e2[-1] > 99                   # five classes, the last ends at the maximum


def test_leader_labels_do_not_collide():
    from measure_it.figures.fig_geography import _place_labels
    pts = [(1.0, 1.0), (1.02, 1.01), (1.05, 0.98), (3.0, 2.0)]
    lab = _place_labels(pts, pts, (0, 5, 0, 4), r_label=0.06)
    for i in range(len(lab)):
        for j in range(i + 1, len(lab)):
            assert np.hypot(lab[i][0] - lab[j][0], lab[i][1] - lab[j][1]) >= 2.3 * 0.06 - 1e-9


def test_leader_lines_do_not_cross_in_a_dense_cluster():
    from measure_it.figures.fig_geography import _place_labels, _seg_point_dist, _segments_cross
    rng = np.random.default_rng(1)
    pts = [tuple(x) for x in rng.uniform([1.0, 1.0], [1.6, 1.5], size=(12, 2))]
    lab = _place_labels(pts, pts, (0, 4, 0, 3), r_label=0.06)
    segs = list(zip(pts, lab))
    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            assert not _segments_cross(*segs[i], *segs[j]), (i, j)
            assert np.hypot(lab[i][0] - lab[j][0], lab[i][1] - lab[j][1]) >= 2.3 * 0.06 - 1e-9
    assert _seg_point_dist((0, 0), (2, 0), (1, 1)) == pytest.approx(1.0)
    assert _segments_cross((0, 0), (1, 1), (0, 1), (1, 0))
    assert not _segments_cross((0, 0), (1, 1), (1, 1), (2, 0))               # a shared end point is not a crossing


def test_plotly_rings_are_clockwise():
    import shapely

    from measure_it.figures.data_html import d3_orient
    ccw = shapely.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])          # RFC 7946 winding
    out = d3_orient([ccw])[0]
    assert not out.exterior.is_ccw


def test_pipeline_registers_figure_steps_after_scoring():
    names = [s.name for s in P.STEPS]
    for step in ("figures_schematics", "figures_data"):
        assert step in names
        assert names.index(step) > names.index("scoring_report")
        assert names.index(step) < names.index("validate")
    assert P.STEP_BY_NAME["figures_data"].target == "measure_it.figures.data_figures:run"
    assert P.STEP_BY_NAME["figures_schematics"].target == "measure_it.figures.schematics:run"
    assert {"scoring_report", "omics_graph", "stanford_models", "wearable_models"} <= set(
        P.STEP_BY_NAME["figures_data"].deps)
    assert "registry_merge" in P.STEP_BY_NAME["figures_schematics"].deps
    assert names[-1] == "validate"


def test_caption_block_is_inserted_once_before_figure_10(tmp_path, monkeypatch):
    from measure_it.figures import data_figures as D
    p = tmp_path / "CAPTIONS.md"
    p.write_text("# Figure captions\n\n## Figure 1. A\n\ntext\n\n---\n\n## Figure 10. B\n\ntext\n")
    monkeypatch.setattr(D, "caption_text", lambda res: f"{D.MARK_START}\nblock v1\n{D.MARK_END}\n")
    D.write_captions({}, p)
    t = p.read_text()
    assert t.index("block v1") < t.index("## Figure 10.") and t.index("## Figure 1.") < t.index("block v1")
    monkeypatch.setattr(D, "caption_text", lambda res: f"{D.MARK_START}\nblock v2\n{D.MARK_END}\n")
    D.write_captions({}, p)
    t = p.read_text()
    assert "block v1" not in t and t.count(D.MARK_START) == 1 and "block v2" in t and "## Figure 10. B" in t


# ================================================================================================= data tests
@pytest.fixture(scope="module")
def figs():
    from measure_it.figures import data_figures as D
    out = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        for key in D.BUILDERS:
            out[key] = D.render(key, check=False)
    glyph = [str(x.message) for x in w if "glyph" in str(x.message).lower() or "missing from" in str(x.message).lower()]
    out["_glyph_warnings"] = glyph
    return out


def _text(res) -> str:
    return " ".join(DC.figure_text(res.fig))


@pytest.mark.data
@pytest.mark.parametrize("key", list(STEMS))
def test_layout_is_clean(figs, key):
    res = figs[key]
    assert res.stem == STEMS[key]
    assert DC.layout_problems(res.fig) == []
    assert DC.min_font_size(res.fig) >= S.FONT_SIZES["min"] - 1e-9
    w, h = res.fig.get_size_inches()
    assert abs(w - S.DOUBLE_COL) < 1e-9 and h <= S.MAX_HEIGHT + 1e-9


@pytest.mark.data
@pytest.mark.parametrize("key", list(STEMS))
def test_text_keeps_a_margin_from_the_canvas_edge(figs, key):
    f = figs[key].fig
    f.canvas.draw()
    r = f.canvas.get_renderer()
    for t in DC._visible_texts(f):
        bb = t.get_window_extent(r)
        if bb.width <= 0:
            continue
        m = min(bb.x0 - f.bbox.x0, f.bbox.x1 - bb.x1, bb.y0 - f.bbox.y0, f.bbox.y1 - bb.y1) / f.dpi
        assert m >= 0.015, (t.get_text(), m)


@pytest.mark.data
def test_no_missing_glyphs(figs):
    assert figs["_glyph_warnings"] == []


@pytest.mark.data
@pytest.mark.parametrize("key", list(STEMS))
def test_product_language(figs, key):
    assert DC.banned_phrases_in(_text(figs[key])) == []


@pytest.mark.data
def test_fig3_labels_proxy_rx_and_acute_covid(figs):
    from measure_it.figures.fig_wearable import FEATURE_GROUPS
    res = figs["fig03"]
    t = _text(res)
    assert "ME/CFS-like PROXY" in t and "not ME/CFS" in t
    assert "Fibromyalgia (Rx-defined)" in t and f"{res.facts['phenotypes']['rx_fibromyalgia']['n_cases']} cases" in t
    assert "NOT Long COVID" in t and "Long COVID (Uwakwe 2025): null" in t
    grouped = sorted(f for _, feats in FEATURE_GROUPS for f, _ in feats)
    assert grouped == res.facts["features"]                      # every table feature drawn exactly once


@pytest.mark.data
def test_fig3_numbers_match_tables(figs):
    f = figs["fig03"].facts
    ep = pd.read_csv(TABLES / "stanford_acute_endpoints.csv")
    r = ep[(ep["endpoint"] == "mean step z, days 0..14") & (ep["group"] == "pooled_placebo_eligible")
           & (ep["onset_type"] == "paired difference true minus placebo")].iloc[0]
    assert f["steps_0_14_diff"]["estimate"] == pytest.approx(r["estimate"])
    summ = pd.read_csv(TABLES / "nhanes_signature_summary.csv").set_index("phenotype_id")
    for p, v in f["phenotypes"].items():
        assert v["n_cases"] == summ.loc[p, "n_cases"] and v["n_fdr"] == summ.loc[p, "n_fdr_significant"]
    assert f["uwakwe_wearable_auroc"]["auroc"] == pytest.approx(0.5, abs=0.1)   # the documented null


@pytest.mark.data
def test_fig4_readings_match_the_nhanes_writeup(figs):
    from measure_it.wearables.nhanes_cohort import TARGETS
    md = (RESULTS / "WEARABLE_NHANES_RESULTS.md").read_text().split("## 2.")[0]
    verd = figs["fig04"].facts["verdicts"]
    for t, v in verd.items():
        row = next(line for line in md.splitlines() if line.startswith(f"| {TARGETS[t].label} |"))
        assert ("Test 1 SUPPORTED" in row) == (v["test1"]["verdict"] == "supported"), t
        exp2 = ("supported, robust" if "Test 2 SUPPORTED" in row else "not robust" if "NOT ROBUST" in row
                else "negative" if "NEGATIVE" in row else "null")
        assert v["test2"]["verdict"] == exp2, t
    t = _text(figs["fig04"])
    assert "not estimable:" in t and "block perm." in t          # the Uwakwe Long COVID null is shown


@pytest.mark.data
def test_fig5_draws_the_table_edges_and_flags(figs):
    from measure_it.omics.graph import DEMO
    e = pd.read_csv(TABLES / "fig5_edges.csv")
    cp = e[(e["edge_type"] == "condition_pathway") & e["condition_id"].isin(DEMO)]
    d = figs["fig05"].facts["demo"]
    assert d["n_condition_pathway_edges"] == len(cp) and d["n_pathways"] == cp["target_node"].nunique()
    assert d["n_flagged_condition_pathway_edges"] == int(cp["post_hoc_flag"].fillna("").ne("").sum())
    t = _text(figs["fig05"])
    assert "not patient multi-omics" in t and "curated" in t and "post-hoc flag" in t


@pytest.mark.data
def test_fig6_state_level_values_and_inherited_caveat(figs):
    from measure_it.geography.query import get_condition_burden
    f = figs["fig06"].facts
    sb = get_condition_burden("long_covid", "state")
    vals = [r["value"] for r in sb["rows"] if r["value"] != UNKNOWN]
    assert f["n_states"] == len(vals) == 51 and f["state_level"] == "A" and f["state_resolution"] == "state"
    assert f["max"][1] == pytest.approx(max(vals)) and f["min"][1] == pytest.approx(min(vals))
    t = _text(figs["fig06"])
    if f["county_inherited"]:   # long_covid_burden: hps_inherited
        assert "Inherited burden" in t and "not county prevalence" in t and "STATE" in t
    else:                       # brfss_sae (default since 2026-10-07): the caveat box says it is modelled
        assert "Modelled burden" in t and "not observed county prevalence" in t
    for lv in "ABCD":
        assert re.search(rf"\b{lv}\b", t)                        # evidence-level key lists A-D


@pytest.mark.data
def test_fig7_trial_sites_are_city_centroids(figs):
    f = figs["fig07"].facts
    assert set(f["trial_site_resolution"]) == {"city_centroid"}
    assert f["n_trial_locations_drawn"] == f["n_trial_locations"]
    assert f["n_hrsa_drawn"] <= f["n_hrsa_sites"] and f["n_nih_orgs_drawn"] <= f["n_nih_orgs"]
    assert 0 < f["specialist_share"] < 0.5
    t = _text(figs["fig07"])
    assert "CITY CENTROIDS" in t and "Specialist-only variant" in t


@pytest.mark.data
@pytest.mark.parametrize("key", ["fig08", "fig09"])
def test_map_leader_lines_are_untangled(figs, key):
    from measure_it.figures.fig_geography import _seg_point_dist, _segments_cross
    segs = figs[key].facts["leaders"]
    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            assert not _segments_cross(*segs[i], *segs[j]), (key, i + 1, j + 1)
            assert _seg_point_dist(*segs[i], segs[j][1]) >= 0.06, (key, i + 1, j + 1)   # no leader under a label


@pytest.mark.data
def test_fig8_top10_incomplete_and_weight_sets_match_scoring(figs):
    f = figs["fig08"].facts
    demo = pd.read_csv(TABLES / "demo_top10_regions.csv", dtype={"geo_id": str})
    assert [t["geo_id"] for t in f["top10"]] == list(demo["geo_id"])
    assert [t["p95"] for t in f["top10"]] == pytest.approx(list(demo["rank_interval_p95"]))
    inc = pd.read_csv(TABLES / "scoring_incomplete_burden.csv", dtype={"geo_id": str})
    inc = inc[(inc["condition_id"] == "long_covid_or_me_cfs")
              & (inc["measurement_id"] == "wearable_autonomic_activity_monitoring") & (inc["geo_level"] == "county")]
    assert f["n_incomplete"] == len(inc)
    st = pd.read_csv(TABLES / "test7_region_stability_primary.csv", dtype={"geo_id": str})
    st = st[(st["condition_id"] == "long_covid_or_me_cfs") & (st["geo_level"] == "county")].set_index("geo_id")
    for gid, ranks in f["weight_set_ranks"].items():
        for ws, rk in ranks.items():
            assert rk == st.loc[gid, f"rank_{ws}"], (gid, ws)
    t = _text(figs["fig08"])
    assert "incomplete set burden" in t and "not a validated diagnostic pathway" in t


@pytest.mark.data
def test_fig9_candidates_are_find_candidate_clinics_output(figs):
    from measure_it.facilities.matching import find_candidate_clinics
    from measure_it.geography.query import get_geographic_context
    f = figs["fig09"].facts
    ids = []
    for c in ("long_covid", "me_cfs"):
        r = find_candidate_clinics("06073", c, "wearable_autonomic_activity_monitoring")
        assert r["radius_km"] == f["radius_km"]
        ids += [x["facility_id"] for x in r["candidates"]]
    assert f["site_ids"] == list(dict.fromkeys(ids))
    ctx = get_geographic_context("06073")
    lv = {b["condition_id"]: b["burden_evidence_level"] for b in ctx["burden_by_condition"]}
    for cid, _, level, _ in f["burden_rows"]:
        assert lv[cid] == level
    assert len(f["level_d"]) == sum(v == "D" for v in lv.values())
    t = _text(figs["fig09"])
    assert UNKNOWN in t and ("state, inherited" in t or "county, modelled" in t) and "never one combined score" in t
    # NIH counts are scoped (all registry conditions vs the two query conditions) and match the context
    assert "any registry condition" in t
    per = {b["condition_id"]: b["nih_projects_n"] for b in ctx["burden_by_condition"]}
    assert f["nih_member_projects"] == {c: int(per[c]) for c in ("long_covid", "me_cfs")}
    # a characteristic that find_candidate_clinics reports absent from the pool is stated on the figure
    for cid, c in f["characteristics_absent"]:
        assert c in find_candidate_clinics("06073", cid, "wearable_autonomic_activity_monitoring")[
            "characteristics_absent_in_pool"]
    if f["characteristics_absent"]:
        assert "absence in public registries" in t


# ------------------------------------------------------------------------------------------ published files
def _png_dpi(path) -> float:
    data = path.read_bytes()
    i = data.find(b"pHYs")
    assert i > 0, "no pHYs chunk"
    ppu_x, _, unit = struct.unpack(">IIB", data[i + 4:i + 13])
    assert unit == 1
    return ppu_x * 0.0254


@pytest.mark.data
@pytest.mark.parametrize("stem", list(STEMS.values()))
def test_published_png_and_svg(stem):
    png, svg = FIGURES / f"{stem}.png", FIGURES / f"{stem}.svg"
    assert png.exists() and svg.exists()
    assert round(_png_dpi(png)) == 300
    w = struct.unpack(">I", png.read_bytes()[16:20])[0]
    assert abs(w / 300 - S.DOUBLE_COL) < 0.01
    assert "<svg" in svg.read_text()[:2000]          # SVG text is converted to paths (style.save_figure)


@pytest.mark.data
@pytest.mark.parametrize("stem", MAP_STEMS)
def test_interactive_maps_open_offline(stem):
    p = MAPS / f"{stem}.html"
    assert p.exists()
    html = p.read_text()
    assert 'src="plotly.min.js"' in html and "cdn.plot.ly" not in html
    assert (MAPS / "plotly.min.js").exists()
    assert f'id="{stem}"' in html


@pytest.mark.data
def test_captions_block_is_current(figs):
    from measure_it.figures import data_figures as D
    text = (FIGURES / "CAPTIONS.md").read_text()
    assert text.count(D.MARK_START) == 1 and text.count(D.MARK_END) == 1
    block = text[text.index(D.MARK_START):text.index(D.MARK_END) + len(D.MARK_END)]
    assert block.strip() == D.caption_text({k: figs[k] for k in D.BUILDERS}).strip()
    assert DC.banned_phrases_in(block) == []
    for n in range(3, 10):
        assert f"## Figure {n}." in block


@pytest.mark.data
def test_rebuild_is_byte_identical(tmp_path, figs):
    from measure_it.figures import data_figures as D
    hashes = []
    for run in ("a", "b"):
        res = D.render("fig04")
        with S.style_context():
            paths = S.save_figure(res.fig, res.stem, tmp_path / run)
        hashes.append([hashlib.sha256(p.read_bytes()).hexdigest() for p in paths])
    assert hashes[0] == hashes[1]
