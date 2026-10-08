"""Tests for the shared figure style and the SPEC Figure 1 / 2 / 10 schematics.

The schematics read SOURCE_REGISTRY.yaml, docs/SPEC.md and the adapter code (all in the repository), not
data/processed, so these tests need no pipeline run. `Canvas.layout_problems()` measures every text
artist: text outside the figure, outside the box it belongs to, or overlapping other text fails.
"""
from __future__ import annotations

import copy
import hashlib
import warnings
from pathlib import Path

import matplotlib
import pytest

from measure_it.figures import schematics as SC
from measure_it.figures import style as S
from measure_it.measurements import adapters as A
from measure_it.measurements.adapter import MeasurementAdapter, MeasurementDescription
from measure_it.provenance import DATA_LAYERS


# ------------------------------------------------------------------------------------------ style
def test_categorical_palette_is_fixed_and_never_cycles():
    assert len(S.CATEGORICAL) == 6 and len(set(S.CATEGORICAL)) == 6
    assert S.categorical(3) == list(S.CATEGORICAL[:3])
    with pytest.raises(ValueError):
        S.categorical(7)
    assert S.OKABE_ITO["yellow"] not in S.CATEGORICAL          # fails lightness band / contrast


def test_every_data_layer_has_a_colour_and_label():
    assert set(DATA_LAYERS) <= set(S.LAYER_COLORS)
    assert set(DATA_LAYERS) <= set(S.LAYER_LABELS)
    four = [S.LAYER_COLORS[k] for k in ("person", "condition_molecular", "geographic", "facility")]
    assert len(set(four)) == 4


def test_sequential_and_diverging_ramps_are_monotone_in_lightness():
    L = [S.oklab(c)[0] for c in S.SEQUENTIAL_BLUE]
    assert all(a > b for a, b in zip(L, L[1:]))
    lo = [S.oklab(c)[0] for c in (*S.DIVERGING_LOW, S.DIVERGING_MID)]
    hi = [S.oklab(c)[0] for c in (S.DIVERGING_MID, *S.DIVERGING_HIGH)]
    assert all(a < b for a, b in zip(lo, lo[1:])) and all(a > b for a, b in zip(hi, hi[1:]))
    # arms matched step for step (within 0.01 OKLab L)
    for a, b in zip(S.DIVERGING_LOW, reversed(S.DIVERGING_HIGH)):
        assert abs(S.oklab(a)[0] - S.oklab(b)[0]) < 0.01


def test_text_ink_contrast():
    for k in ("primary", "secondary", "muted"):
        assert S.contrast_ratio(S.INK[k], S.INK["surface"]) >= 4.5, k
    assert S.contrast_ratio(S.INK["muted"], S.INK["panel"]) >= 4.5


def test_fonts_are_bundled_with_matplotlib():
    from matplotlib import font_manager as fm
    mpl_data = Path(matplotlib.get_data_path()).resolve()
    for fam in (S.FONT_FAMILY, S.MONO_FAMILY):
        path = Path(fm.findfont(fm.FontProperties(family=fam), fallback_to_default=False)).resolve()
        assert mpl_data in path.parents, path


def test_save_figure_is_deterministic(tmp_path):
    from matplotlib.figure import Figure
    hashes = []
    for run in ("a", "b"):
        with S.style_context():
            fig = Figure(figsize=(S.SINGLE_COL, 1.5))
            ax = fig.add_subplot()
            ax.plot([0, 1, 2], [1, 0, 2], color=S.CATEGORICAL[0])
            paths = S.save_figure(fig, "tiny", tmp_path / run)
        hashes.append([hashlib.sha256(p.read_bytes()).hexdigest() for p in paths])
    assert [p.suffix for p in paths] == [".png", ".svg"]
    assert hashes[0] == hashes[1]


# ------------------------------------------------------------------------------------------ parsing helpers
@pytest.mark.parametrize("res, expected", [
    ("none", []),
    ("none (catalog metadata; spatial fields not used)", []),
    ("national, state (no county)", ["state", "national"]),
    ("county (CT legacy counties 2001-2022, CT planning regions 2023)", ["county"]),
    ("county (2024 vintage, CT planning regions); zcta; national reference row", ["zcta", "county", "national"]),
    ("practice ZIP5 -> ZCTA internal point (zcta_centroid); 2024 county aggregates", ["zcta", "county"]),
    ("site point (geoPoint, city-level) or ZCTA centroid; county_fips (2024) from the site ZIP's ZCTA first, "
     "else geoPoint point-in-polygon", ["point", "zcta", "county"]),
    ("census tract / county subdivision / county / facility point; county summary (2024)",
     ["point", "zcta", "county"]),
])
def test_resolution_levels(res, expected):
    assert SC.resolution_levels(res) == expected


def test_canonical_table_folds_partitions_and_rejects_sentences():
    assert SC.canonical_table("condition_molecular_evidence__gwas_catalog") == "condition_molecular_evidence"
    assert SC.canonical_table("geo_context") == "geo_context"
    assert SC.canonical_table("county_boundaries_2024.geoparquet") == "county_boundaries_2024.geoparquet"
    assert SC.canonical_table("results/tables/x.csv") == "results/tables/x.csv"
    assert SC.canonical_table("none written to data/processed; used by measure_it.omics.graph") is None


def test_access_flags():
    assert SC.access_flags({"access_conditions": "open download, no registration"}) == []
    gated = SC.access_flags({"access_conditions": "mapMECFS data files: registered account required "
                                                  "(anonymous requests return HTTP 403 'Login Required')"})
    assert gated == ["login-gated files: not automated"]
    keyed = SC.access_flags({"access_conditions": "api.census.gov data requests required an API key"})
    assert keyed == ["key-gated API not used; open route used"]


def test_spec_mcp_tools_parsed_from_spec():
    tools = SC.spec_mcp_tools()
    assert len(tools) == len(set(tools)) >= 15
    assert tools[0] == "search_condition" and "trace_evidence" in tools
    assert "rank_deployment_opportunities" in tools


# ------------------------------------------------------------------------------------------ figures
def _build(builder, *args):
    with S.style_context(), warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        cv = builder(*args)
        probs = cv.layout_problems()
    glyph = [str(x.message) for x in w if "missing from" in str(x.message).lower() or "glyph" in str(x.message).lower()]
    return cv, probs, glyph


@pytest.fixture(scope="module")
def registry():
    """One registry snapshot per module: other contributors may re-merge SOURCE_REGISTRY.yaml mid-run."""
    return SC.load_registry()


@pytest.fixture(scope="module")
def fig1(registry):
    return _build(SC.figure1, copy.deepcopy(registry))


@pytest.fixture(scope="module")
def fig2(registry):
    return _build(SC.figure2, copy.deepcopy(registry))


@pytest.fixture(scope="module")
def fig10():
    return _build(SC.figure10)


@pytest.mark.parametrize("name", ["fig1", "fig2", "fig10"])
def test_layout_is_clean(name, request):
    cv, probs, glyph = request.getfixturevalue(name)
    assert probs == [], probs[:10]
    assert glyph == []
    assert abs(cv.W - S.DOUBLE_COL) < 1e-9
    assert cv.H <= S.MAX_HEIGHT + 1e-9, f"{name} is {cv.H:.2f} in tall (> 247 mm)"


@pytest.mark.parametrize("name", ["fig1", "fig2", "fig10"])
def test_product_language(name, request):
    cv, _, _ = request.getfixturevalue(name)
    assert SC.banned_phrases_in(" ".join(cv.all_text())) == []


def test_banned_phrase_matching_uses_whole_words():
    assert SC.banned_phrases_in("wearable features improves AUROC; secures funding") == []
    assert SC.banned_phrases_in("This PROVES it; the best clinic") == ["best clinic", "proves"]


def test_fig1_shows_chain_layers_joins_seam_and_tools(fig1, registry):
    cv, _, _ = fig1
    text = " ".join(cv.all_text())
    for i, (step, _) in enumerate(SC.CHAIN, 1):
        assert f"{i}  {step}" in text
    for join in ("PERSON-LEVEL", "CONCEPT-LEVEL", "ECOLOGICAL", "concept-level join", "ecological join"):
        assert join in text
    assert "MeasurementAdapter" in text and "MCP / agent layer" in text
    for tool in SC.spec_mcp_tools():
        assert tool in text.replace("\n", ""), tool
    reg = registry
    n = {}
    for s in reg["sources"]:
        n[s["data_layer"]] = n.get(s["data_layer"], 0) + 1
    for i, layer in enumerate(("person", "condition_molecular", "geographic", "facility"), 1):
        k = n.get(layer, 0)
        assert f"LAYER {i} · {k} SOURCE{'S' if k != 1 else ''}" in text
    assert f"indexes {len(reg['sources'])} public sources" in text
    other = [s for s in reg["sources"] if s["data_layer"] not in ("person", "condition_molecular", "geographic",
                                                                  "facility", "ontology")]
    if other:
        assert "(including " in text
    for a in A.list_adapters():
        assert a["class"] in text


def test_fig2_is_generated_from_the_registry(fig2, registry):
    cv, _, _ = fig2
    reg = registry
    text = " ".join(cv.all_text())
    flat = text.replace(" ", " ")
    for s in reg["sources"]:
        assert SC.display_name(s) in text, s["source_id"]
    for st in SC.STATUS_ORDER:
        k = sum(1 for s in reg["sources"] if s.get("status") == st)
        assert f"{st} ({k})" in text
    for tabs in SC.source_table_map(reg["sources"]).values():
        for t in tabs:
            assert SC.table_label(t) in flat, t
    gated = [s for s in reg["sources"] if SC.access_flags(s)]
    assert text.count("login-gated files: not automated") == sum(
        "login-gated files: not automated" in SC.access_flags(s) for s in gated)
    partial = [s for s in reg["sources"] if s.get("status") == "partial"]
    assert all(SC.display_name(s) in text for s in partial)
    assert f"({len(reg['sources'])} sources" in text


def test_fig2_generalises_to_other_registries():
    """A synthetic registry (no real source) with a blocked source, a sentence output and a new layer."""
    base = {"publisher": "Example agency", "unit_of_observation": "county x year", "access_conditions": "open",
            "geographic_resolution": "county", "status": "ingested"}
    reg = {"generated_at": "2026-09-23T00:00:00+00:00", "sources": [
        {**base, "source_id": "a", "name": "Source A", "data_layer": "geographic", "processed_outputs": ["t1__a", "t2"]},
        {**base, "source_id": "b", "name": "Source B", "data_layer": "geographic", "processed_outputs": ["t1__b"],
         "status": "blocked", "access_conditions": "registered account required"},
        {**base, "source_id": "c", "name": "Source C", "data_layer": "person", "processed_outputs": ["t1__c"],
         "geographic_resolution": "none"},
        {**base, "source_id": "d", "name": "Source D", "data_layer": "metadata",
         "processed_outputs": ["none written to data/processed; used by an analysis"]},
    ]}
    cv, probs, glyph = _build(SC.figure2, copy.deepcopy(reg))
    assert probs == [] and glyph == []
    text = " ".join(cv.all_text())
    assert "blocked (1)" in text and "ingested (3)" in text and "partial (0)" in text
    assert "login-gated files: not automated" in text
    assert "no processed table written" in text
    chips = SC.table_chips(reg["sources"])
    t1 = next(c for c in chips if "t1" in c["tables"])
    assert t1["sources"] == ["a", "b", "c"] and t1["band"] == "geographic" and t1["cross_layer"] == ["person"]


def test_fig10_matches_the_adapter_code(fig10):
    cv, _, _ = fig10
    text = " ".join(cv.all_text())
    for m in MeasurementAdapter.__abstractmethods__:
        assert f"{m}(" in text
    for f in (f.name for f in __import__("dataclasses").fields(MeasurementDescription)):
        assert f in text.replace(" ", " "), f
    for info in A.list_adapters():
        d = A.get_adapter(info["name"]).describe_measurement()
        assert info["class"] in text
        for ph in d.phenotypes:
            assert ph in text, ph
        for mid in d.measurement_ids:
            assert mid in text
    stub = [i for i in A.list_adapters() if i["status"] == "stub"]
    assert stub, "the capillaroscopy stub should be registered"
    assert SC._probe_stub_methods(A.get_adapter(stub[0]["name"])) == {
        m: "raises NotImplementedError" for m in ("preprocess", "embed", "phenotype_score")}
    assert text.count("raises NotImplementedError") == 3
    assert "CAPRIO DATA NOT USED" in text


def test_build_all_writes_png_and_svg_deterministically(tmp_path):
    first = SC.build_all(tmp_path / "a")
    second = SC.build_all(tmp_path / "b")
    for key, stem in SC.FIG_STEMS.items():
        assert [p.name for p in first[key]] == [f"{stem}.png", f"{stem}.svg"]
        for p, q in zip(first[key], second[key]):
            assert p.read_bytes() == q.read_bytes(), p.name
    from PIL import Image
    with Image.open(first["fig01"][0]) as im:
        dpi = im.info.get("dpi", (0, 0))
    assert round(dpi[0]) == S.DPI


# ------------------------------------------------------------------------------------------ review regressions
def test_stacked_lines_keep_1_2_leading_so_underscores_clear_the_next_line():
    """kvlist items and stacked key/value rows are set at >= 1.2 x size baseline-to-baseline. At the old ~1.05 x
    step, an underscore on one line touched the ascenders of the line below (read as an overbar)."""
    cv = SC.Canvas(3.0, 3.0)
    card = SC.Card("T", S.LAYER_COLORS["person"], rows=[
        ("kvlist", ("phenotypes (score columns):", ["low_activity", "circadian_disruption", "activity_fragmentation"])),
        ("kv", ("a_rather_long_key_name_here:", "value_that_is_long_enough_to_force_stacking_below")),
    ])
    with S.style_context():
        SC.render_card(cv, SC.Box(0.1, 0.1, 2.0, 2.5), card)
    ys = {a.get_text(): a.get_position()[1] for a, _ in cv.texts}
    step = SC.STACK_LEADING * 5.6 * SC.PT - 1e-6
    items = ["phenotypes (score columns):", "low_activity", "circadian_disruption", "activity_fragmentation"]
    for upper, lower in zip(items, items[1:]):
        assert ys[upper] - ys[lower] >= step, (upper, lower, ys[upper] - ys[lower])
    key = "a_rather_long_key_name_here:"
    value_lines = [t for t in ys if t.startswith("value_")]
    assert value_lines and ys[key] - ys[value_lines[0]] >= step


def test_card_text_keeps_clear_of_the_border():
    """An unwrapped overline that runs into the card padding is a layout problem, not just one outside the box."""
    cv = SC.Canvas(3.0, 1.5)
    card = SC.Card("Title", S.LAYER_COLORS["geographic"], badge="5", overline="LAYER 3 · 9 SOURCES · EXTRA")
    with S.style_context():
        SC.render_card(cv, SC.Box(0.1, 0.1, 1.25, 1.0), card)
        probs = cv.layout_problems()
    assert any("LAYER 3" in p for p in probs)


def test_fig10_draws_one_implements_arrow_per_adapter(fig10):
    cv, _, _ = fig10
    text = " ".join(cv.all_text())
    infos = A.list_adapters()
    n_stub = sum(i["status"] == "stub" for i in infos)
    assert text.count("implements (stub)") == min(n_stub, 1)       # the figure draws the first stub only
    assert text.count("implements") - text.count("implements (stub)") - text.count("implements these") == sum(
        i["status"] == "implemented" for i in infos)
