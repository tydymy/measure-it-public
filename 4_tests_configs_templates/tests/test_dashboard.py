"""Headless tests for the Streamlit dashboard (dashboard/app.py + dashboard/pages/*).

streamlit.testing.v1.AppTest runs the home page and each page script in-process with the default inputs and with
alternative inputs (other levels, layers, weight sets, overlays, conditions, and inputs that must come back UNKNOWN /
NOT AVAILABLE). Each run must finish without an exception and render the page's key elements, including the
'Limitations and provenance' panel; the rendered text is scanned with the project's product-language check.
Tests marked `data` need data/processed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

from measure_it.config import PROCESSED, UNKNOWN
from measure_it.validate import forbidden_hits

ROOT = Path(__file__).resolve().parents[1]
DASH = ROOT / "dashboard"
PAGES = {
    "map": DASH / "pages" / "1_National_Opportunity_Map.py",
    "condition": DASH / "pages" / "2_Condition_Explorer.py",
    "measurement": DASH / "pages" / "3_Measurement_Explorer.py",
    "geography": DASH / "pages" / "4_Geography_Explorer.py",
    "recommendation": DASH / "pages" / "5_Deployment_Recommendation.py",
}
HAS_DATA = (PROCESSED / "deployment_opportunities.parquet").exists()
TIMEOUT = 300
needs_data = [pytest.mark.data, pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")]


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    monkeypatch.setenv("MEASURE_IT_OFFLINE", "1")


def _app(path: Path):
    from streamlit.testing.v1 import AppTest
    return AppTest.from_file(str(path), default_timeout=TIMEOUT)


def _texts(at) -> list[str]:
    out = []
    for kind in ("title", "header", "subheader", "markdown", "caption", "info", "warning", "error", "success"):
        out += [str(e.value) for e in getattr(at, kind)]
    return out


def _check(at, title: str) -> list[str]:
    assert not at.exception, [e.message for e in at.exception]
    assert [t.value for t in at.title][:1] == [title]
    subs = [s.value for s in at.subheader]
    assert "Limitations and provenance" in subs, subs
    texts = _texts(at)
    bad = [(t[:80], h) for t in texts for h in forbidden_hits(t)]
    assert not bad, bad
    return texts


def _has(texts: list[str], needle: str) -> bool:
    return any(needle.lower() in t.lower() for t in texts)


def _plotly(at) -> int:
    return len(at.get("plotly_chart"))


# ------------------------------------------------------------------------------------------------ pure

def test_every_page_file_exists_and_has_a_limitations_panel():
    assert (DASH / "app.py").exists()
    for p in PAGES.values():
        assert p.exists(), p
        assert "limitations_panel(" in p.read_text(), p
    assert "limitations_panel(" in (DASH / "app.py").read_text()


def test_dashboard_source_text_uses_product_language():
    for p in [DASH / "app.py", *PAGES.values(), *(DASH / "components").glob("*.py")]:
        assert not forbidden_hits(p.read_text()), p


def test_choropleth_builder_on_a_tiny_geojson():
    sys.path.insert(0, str(DASH))
    from components import maps as M
    sq = lambda x: {"type": "Polygon", "coordinates": [[[x, 30], [x + 1, 30], [x + 1, 31], [x, 31], [x, 30]]]}  # noqa
    fc = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "id": g, "properties": {"geo_id": g}, "geometry": sq(-100 + i)}
        for i, g in enumerate(["01001", "01003", "01005"])]}
    df = pd.DataFrame({"geo_id": ["01001", "01003"], "Region": ["A", "B"], "value": [0.2, float("nan")],
                       "Why no value": ["", "level D"]})
    fig = M.choropleth_map(fc, df, "value", ["Region"], "test", nodata_hover=["Region", "Why no value"],
                           highlight=["01001"])
    names = [t.name for t in fig.data]
    assert names[0].startswith("no value") and "top ranked (red outline)" in names
    grey = fig.data[0]
    assert list(grey.locations) == ["01003", "01005"]          # NaN row + a feature outside the layer
    assert "level D" in grey.hovertext[0] and grey.hovertext[1] == "not in this layer"
    data_trace = fig.data[1]
    assert list(data_trace.locations) == ["01001"] and [f["id"] for f in data_trace.geojson["features"]] == ["01001"]


# ------------------------------------------------------------------------------------------------ home

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_home():
    at = _app(DASH / "app.py").run()
    texts = _check(at, "Measure It to Cure It: public-data alpha")
    assert _has(texts, "candidate deployment opportunity")
    assert any(m.label == "Sources" and int(m.value) >= 20 for m in at.metric)
    assert "Four data layers, never joined as the same people" in [s.value for s in at.subheader]


# ------------------------------------------------------------------------------------------------ 1 map

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_map_default():
    at = _app(PAGES["map"]).run()
    texts = _check(at, "National Opportunity Map")
    assert _plotly(at) == 1
    assert _has(texts, "Burden evidence (always shown)") and _has(texts, "not county prevalence")
    assert _has(texts, "symptom/comorbidity proxy")
    assert len(at.dataframe) >= 1 and at.dataframe[0].value.shape[0] >= 10
    fig = at.get("plotly_chart")[0]
    assert fig is not None


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_map_alternatives():
    at = _app(PAGES["map"]).run()
    # state level, burden layer, every overlay on, another weight set
    at.radio(key="map_level").set_value("state").run()
    at.radio(key="map_layer").set_value("Burden (with evidence level)").run()
    for k in ("ov_hrsa", "ov_trials", "ov_nih", "ov_spec"):
        at.checkbox(key=k).check()
    at.selectbox(key="map_weights").set_value("burden_led").run()
    _check(at, "National Opportunity Map")
    assert _plotly(at) == 1
    # long COVID burden on counties: the colour bar and a warning say it is the inherited STATE value
    at.radio(key="map_level").set_value("county").run()
    at.selectbox(key="map_condition").set_value("long_covid").run()
    texts = _check(at, "National Opportunity Map")
    fig = at.get("plotly_chart")[0].proto.spec
    if _has(texts, "of counties carry their **state** estimate"):   # long_covid_burden: hps_inherited
        assert _has(texts, "not county prevalence") and "inherited STATE value on counties" in fig
    else:                                                          # brfss_sae (default since 2026-10-07)
        assert _has(texts, "**modelled** small-area estimate") and _has(texts, "not observed county prevalence")
        assert "MODELLED small-area estimate on counties" in fig
    # a level-D condition: burden UNKNOWN, desert access-only, caveat says so
    at.selectbox(key="map_condition").set_value("pots").run()
    texts = _check(at, "National Opportunity Map")
    assert _has(texts, "level D, no usable burden estimate")
    at.radio(key="map_layer").set_value("Diagnostic desert").run()
    texts = _check(at, "National Opportunity Map")
    assert _has(texts, "access-only")
    # a condition outside the precomputed ranking: opportunity layer is UNKNOWN, burden layer still works
    at.selectbox(key="map_condition").set_value("fibromyalgia").run()
    at.radio(key="map_layer").set_value("Deployment opportunity").run()
    texts = _check(at, "National Opportunity Map")
    assert _has(texts, f"Deployment opportunity: {UNKNOWN}")


# ------------------------------------------------------------------------------------------------ 2 condition

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_condition_default_and_alternatives():
    at = _app(PAGES["condition"]).run()
    texts = _check(at, "Condition Explorer")
    subs = [s.value for s in at.subheader]
    for s in ("1. Phenotype features", "2. Wearable abnormalities", "3. Molecular enrichment",
              "4. Candidate measurements", "5. Technologies (FDA)", "6. Burden sources"):
        assert s in subs, s
    assert _has(texts, "Long COVID (post-COVID-19 condition)")
    assert _has(texts, "not patient multi-omics")
    assert _has(texts, "null (not significant) features")
    assert len(at.dataframe) >= 8 and _plotly(at) >= 1
    dims = ["phenotype signal strength", "measurement evidence strength", "technology maturity",
            "regulatory visibility", "deployment complexity"]
    cand = next(df.value for df in at.dataframe if "technology maturity" in df.value.columns)
    assert list(cand.columns[2:7]) == dims
    assert "null result" in set(cand["phenotype signal strength"])          # Long COVID wearable HR: null, shown
    fda = next(df.value for df in at.dataframe if "device classes" in df.value.columns)
    assert fda["measurement"].nunique() >= 3
    # ME/CFS (proxy label with FDR-significant NHANES features) + a phenotype focus
    at.text_input(key="cond_text").input("ME/CFS").run()
    at.selectbox(key="cond_phenotype").set_value("orthostatic_intolerance").run()
    texts = _check(at, "Condition Explorer")
    assert _has(texts, "Myalgic encephalomyelitis") and _plotly(at) >= 2
    # an ICD-10-CM code resolves
    at.text_input(key="cond_text").input("G90.A").run()
    texts = _check(at, "Condition Explorer")
    assert _has(texts, "Postural orthostatic tachycardia syndrome")
    # nonsense -> UNKNOWN, page stops cleanly
    at.text_input(key="cond_text").input("zzqxv blorf").run()
    texts = _check(at, "Condition Explorer")
    assert _has(texts, UNKNOWN)


# ------------------------------------------------------------------------------------------------ 3 measurement

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_measurement_default_and_alternatives():
    at = _app(PAGES["measurement"]).run()
    texts = _check(at, "Measurement Explorer")
    subs = [s.value for s in at.subheader]
    for s in ("1. Performance evidence (metric link)", "2. Conditions where it is studied",
              "3. Signals & public datasets", "4. Relevant trials", "5. FDA device classes",
              "6. Regions with unmet need"):
        assert s in subs, s
    assert _has(texts, "Wearable autonomic / activity monitoring") and _has(texts, "candidate deployment opportunity")
    assert _plotly(at) == 3 and len(at.dataframe) >= 7     # + the metric-link performance chart and table
    perf = next(df.value for df in at.dataframe if "sensitivity @ spec" in df.value.columns)
    assert "OWN-MM-STEPS-POOLED" in set(perf["record"]) and UNKNOWN in set(perf["status"])
    unmet = next(df.value for df in at.dataframe if "readiness pct" in df.value.columns)
    assert "rank interval 5-95 (equal-weight MC)" in unmet.columns and len(unmet) == 10
    assert _has(texts, "Post-infectious syndromes (grouping)**: UNKNOWN / NOT AVAILABLE")
    # another class, state level, another weight set (scored on the fly)
    at.text_input(key="meas_text").input("CPET").run()
    at.radio(key="meas_level").set_value("state").run()
    at.selectbox(key="meas_weights").set_value("equal").run()
    _check(at, "Measurement Explorer")
    assert _plotly(at) == 3
    at.text_input(key="meas_text").input("zzqxv blorf").run()
    texts = _check(at, "Measurement Explorer")
    assert _has(texts, UNKNOWN)


# ------------------------------------------------------------------------------------------------ 4 geography

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_geography_default_and_alternatives():
    at = _app(PAGES["geography"]).run()
    texts = _check(at, "Geography Explorer")
    subs = [s.value for s in at.subheader]
    for s in ("1. Burden by condition", "2. Vulnerability & context", "3. Providers & specialists",
              "4. Trials & NIH research", "5. Candidate technologies & facilities", "6. Nearby health centers (map)"):
        assert s in subs, s
    assert _has(texts, "San Diego County, California")
    assert _has(texts, "inherited by the county") or _has(texts, "**modelled** small-area estimate")
    assert _plotly(at) >= 3 and any(m.label.startswith("Population") for m in at.metric)
    burden = next(df.value for df in at.dataframe if "evidence level" in df.value.columns)
    assert set(burden["evidence level"]) <= {"A", "B", "C", "D"} and "unit" in burden.columns
    inh = burden[burden["inherited state value"].fillna(False).astype(bool)]
    mod = burden[burden["modelled estimate"].fillna(False).astype(bool)]
    assert len(inh) + len(mod) >= 1
    assert inh["burden"].str.contains("(CA STATE value, inherited)", regex=False).all()
    assert mod["burden"].str.contains("(modelled small-area estimate)", regex=False).all()
    spec = next(df.value for df in at.dataframe if "specialists (excl. primary care)" in df.value.columns)
    assert list(spec.columns[:3]) == ["condition", "specialists (excl. primary care)", "specialists per 100k"]
    # a state, another condition, a smaller radius
    at.text_input(key="geo_text").input("Mississippi").run()
    at.selectbox(key="geo_condition").set_value("me_cfs").run()
    at.slider(key="geo_radius").set_value(30).run()
    texts = _check(at, "Geography Explorer")
    assert _has(texts, "Mississippi")
    at.text_input(key="geo_text").input("Zzqxv County, Nowhere").run()
    texts = _check(at, "Geography Explorer")
    assert _has(texts, UNKNOWN)


# ------------------------------------------------------------------------------------------------ 5 recommendation

@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_recommendation_default():
    at = _app(PAGES["recommendation"]).run()
    texts = _check(at, "Deployment Recommendation")
    assert _has(texts, "Find 10 U.S. regions where wearable monitoring")
    table = at.dataframe[0].value
    assert len(table) == 10
    for col in ("rank", "composite", "rank interval 5-95", "P(top 10)", "burden level", "burden pct",
                "vulnerability pct", "desert pct", "capacity pct", "readiness pct", "rank range across weight sets",
                "object_id"):
        assert col in table.columns, col
    assert table["object_id"].str.startswith("opportunity:").all()
    assert list(table.columns[2:8]) == ["composite", "burden pct", "vulnerability pct", "desert pct", "capacity pct",
                                        "readiness pct"]
    phen = next(df.value for df in at.dataframe if "FDR-significant" in df.value.columns)
    assert phen["null results"].str.contains("did not pass").any()           # null results are visible
    assert _plotly(at) == 3                               # contributions, rank intervals, map
    assert len(at.get("download_button")) == 2
    assert _has(texts, "Candidate deployment opportunity (rank 1 of")
    assert any("Trace evidence" in e.label for e in at.expander)


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_recommendation_alternatives():
    at = _app(PAGES["recommendation"]).run()
    at.radio(key="rec_level").set_value("state").run()
    at.selectbox(key="rec_weights").set_value("burden_led").run()
    at.number_input(key="rec_topn").set_value(5).run()
    _check(at, "Deployment Recommendation")
    assert len(at.dataframe[0].value) == 5
    assert {"rank interval 5-95 (equal-weight MC)", "P(top 10) (equal-weight MC)"} <= set(at.dataframe[0].value.columns)
    at.selectbox(key="rec_region").set_value(2).run()
    _check(at, "Deployment Recommendation")
    at.text_input(key="rec_condition").input("zzqxv blorf").run()
    texts = _check(at, "Deployment Recommendation")
    assert _has(texts, UNKNOWN)
