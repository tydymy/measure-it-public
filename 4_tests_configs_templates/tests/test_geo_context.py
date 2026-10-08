"""Tests for the geographic context sources: CDC/ATSDR SVI, ACS 5-year, USDA RUCC (+ census_geography audit).

Unit tests use synthetic inputs. @pytest.mark.data tests read data/processed and data/raw and are skipped
when a table has not been built.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
import yaml

from measure_it.config import PROCESSED, RAW
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.registry import REQUIRED_FIELDS
from measure_it.ingestion import census_acs, rucc, svi


# ============================================================================ unit tests
def test_moe_sum_uses_largest_zero_estimate_moe_only():
    est = pd.DataFrame({"a": [10, 5], "b": [0, 0], "c": [0, 7]})
    moe = pd.DataFrame({"a": [3, 1], "b": [4, 2], "c": [5, 2]})
    out = census_acs.moe_sum(est, moe)
    # row 0: nonzero a (3) + largest zero MOE (5) -> sqrt(9 + 25)
    assert out[0] == pytest.approx(math.sqrt(9 + 25))
    # row 1: nonzero a (1) and c (2) + zero b (2) -> sqrt(1 + 4 + 4)
    assert out[1] == pytest.approx(3.0)


def test_moe_sum_propagates_missing_moe():
    est = pd.DataFrame({"a": [10.0], "b": [3.0]})
    moe = pd.DataFrame({"a": [3.0], "b": [np.nan]})
    assert np.isnan(census_acs.moe_sum(est, moe)[0])


def test_proportion_formula_and_ratio_fallback():
    p, m = census_acs.proportion_with_moe(pd.Series([20.0]), pd.Series([5.0]), pd.Series([100.0]), pd.Series([10.0]))
    assert p[0] == pytest.approx(0.2)
    assert m[0] == pytest.approx(math.sqrt(25 - 0.04 * 100) / 100)
    # negative radicand -> ratio formula (plus sign)
    p, m = census_acs.proportion_with_moe(pd.Series([90.0]), pd.Series([2.0]), pd.Series([100.0]), pd.Series([10.0]))
    assert m[0] == pytest.approx(math.sqrt(4 + 0.81 * 100) / 100)
    # zero denominator -> null, not inf
    p, m = census_acs.proportion_with_moe(pd.Series([0.0]), pd.Series([1.0]), pd.Series([0.0]), pd.Series([1.0]))
    assert np.isnan(p[0]) and np.isnan(m[0])


def test_clean_sentinels_estimates_null_controlled_moe_zero():
    df = pd.DataFrame({"X_001E": [5.0, -666666666.0, 7.0], "X_001M": [-555555555.0, -222222222.0, 2.0]})
    out, counts = census_acs.clean_sentinels(df, ["X_001E", "X_001M"])
    assert np.isnan(out.loc[1, "X_001E"]) and out.loc[0, "X_001E"] == 5
    assert out.loc[0, "X_001M"] == 0.0            # controlled estimate: no sampling error
    assert np.isnan(out.loc[1, "X_001M"])         # MOE not computable
    assert counts == {"X_001E": {"-666666666": 1}, "X_001M": {"-555555555": 1, "-222222222": 1}}


def _shells() -> pd.DataFrame:
    ages = ["Under 5 years", "5 to 9 years", "15 to 17 years", "18 and 19 years", "62 to 64 years",
            "65 and 66 years", "85 years and over"]
    rows = [("B01001", "0", "B01001_001", "Total:", "Total population")]
    rows += [("B01001", "2", f"B01001_{i + 3:03d}", a, "Total population") for i, a in enumerate(ages)]
    return pd.DataFrame(rows, columns=["Table ID", "Indent", "Unique ID", "Label", "Universe"])


def test_b01001_age_bands_from_labels():
    sh = _shells()
    assert census_acs._b01001_age_cells(sh, 18, 64) == ["B01001_006", "B01001_007"]
    assert census_acs._b01001_age_cells(sh, 65, None) == ["B01001_008", "B01001_009"]
    assert census_acs._b01001_age_cells(sh, 0, 17) == ["B01001_003", "B01001_004", "B01001_005"]


def _age_measures():
    """65+ from one published cell and 18-64 as a difference, each cross-checked against B01001 band sums."""
    m65 = census_acs.Measure("age_65_plus", "B09020", "share", "65+", ["B09020_001"], "B01001_001",
                             label_rule="Total:", check_cells=["B01001_020", "B01001_044"],
                             check_label_rule=r"65 and 66 years")
    m1864 = census_acs.Measure("age_18_64", "B01001", "difference", "18-64", ["B01001_001"], "B01001_001",
                               label_rule="Total:", subtract=["B09001_001", "B09020_001"],
                               check_cells=["B01001_007", "B01001_031"], check_label_rule=r"18 and 19 years")
    return [m65, m1864]


def _age_frames(b20_e=30.0):
    geo = ["0500000US01001", "0500000US01003"]
    b01001 = pd.DataFrame({"GEO_ID": geo, "B01001_001E": [100.0, 0.0], "B01001_001M": [-555555555.0, 11.0],
                           "B01001_020E": [10.0, 0.0], "B01001_020M": [4.0, 11.0],
                           "B01001_044E": [20.0, 0.0], "B01001_044M": [5.0, 11.0],
                           "B01001_007E": [25.0, 0.0], "B01001_007M": [6.0, 11.0],
                           "B01001_031E": [25.0, 0.0], "B01001_031M": [6.0, 11.0]})
    b09001 = pd.DataFrame({"GEO_ID": geo, "B09001_001E": [20.0, 0.0], "B09001_001M": [3.0, 11.0]})
    b09020 = pd.DataFrame({"GEO_ID": geo, "B09020_001E": [b20_e, 0.0], "B09020_001M": [2.0, 11.0]})
    return {"B01001": b01001, "B09001": b09001, "B09020": b09020}


def test_age_numerators_use_single_cells_and_small_moe():
    out, _, checks = census_acs.assemble(_age_frames(), _age_measures())
    r = out.iloc[0]
    assert r["n_age_65_plus"] == 30 and r["n_age_65_plus_moe"] == 2.0        # published single-cell MOE
    assert r["n_age_18_64"] == 50                                            # 100 - 20 - 30
    assert r["n_age_18_64_moe"] == pytest.approx(math.sqrt(0 + 9 + 4))       # controlled total MOE = 0
    assert r["pct_age_65_plus"] == pytest.approx(30.0) and r["pct_age_18_64"] == pytest.approx(50.0)
    # denominator 0 -> no share, not 0 %
    assert np.isnan(out.iloc[1]["pct_age_65_plus"]) and np.isnan(out.iloc[1]["pct_age_18_64"])
    assert checks["age_65_plus"]["rows_differing"] == 0 and checks["age_18_64"]["rows_compared"] == 2


def test_age_numerator_disagreeing_with_band_sum_raises():
    with pytest.raises(ValueError, match="differs from the sum"):
        census_acs.assemble(_age_frames(b20_e=31.0), _age_measures())


def test_verify_labels_checks_cell_universe():
    sh = pd.DataFrame([("B09020", "0", "B09020_001", "Total:", "T", "Population 65 years and over"),
                       ("B01001", "0", "B01001_001", "Total:", "Sex by Age", "Total population")],
                      columns=["Table ID", "Indent", "Unique ID", "Label", "Title", "Universe"])
    ok = census_acs.Measure("age_65_plus", "B09020", "share", "65+", ["B09020_001"], "B01001_001",
                            label_rule="Total:", cell_universes={"B09020_001": "Population 65 years and over"})
    assert census_acs.verify_labels([ok], sh)["age_65_plus"]["denominator"] == "B01001_001"
    bad = census_acs.Measure("age_65_plus", "B09020", "share", "65+", ["B09020_001"], "B01001_001",
                             label_rule="Total:", cell_universes={"B09020_001": "Population under 18 years"})
    with pytest.raises(ValueError, match="universe"):
        census_acs.verify_labels([bad], sh)


def test_rucc_adjacency_parsing():
    assert rucc.adjacency_from_definition("Urban population of 20,000 or more, adjacent to a metro area",
                                          "nonmetro") is True
    assert rucc.adjacency_from_definition("Urban population of fewer than 5,000, not adjacent to a metro area",
                                          "nonmetro") is False
    assert rucc.adjacency_from_definition("Counties in metro areas of 1 million population or more", "metro") is None
    with pytest.raises(ValueError):
        rucc.adjacency_from_definition("something else", "nonmetro")


def _codes() -> pd.DataFrame:
    return pd.DataFrame({"rucc_2023": [1, 6], "metro_status": ["metro", "nonmetro"],
                         "code_definition": ["Counties in metro areas of 1 million population or more",
                                             "Urban population of 5,000 to 20,000, adjacent to a metro area"]})


def test_rucc_transform_flags_and_missing_code():
    wide = pd.DataFrame({
        "FIPS": ["01001", "01005", "60030"], "State": ["AL", "AL", "AS"],
        "County_Name": ["A", "B", "Rose Island"], "Population_2020": ["100", "50", "0"],
        "RUCC_2023": ["1", "6", ""],
        "Description": ["Metro - Counties in metro areas of 1 million population or more",
                        "Nonmetro - Urban population of 5,000 to 20,000, adjacent to a metro area",
                        "Not Applicable"]})
    out = rucc.transform(wide, _codes())
    assert out["is_metro"].tolist()[:2] == [True, False]
    assert pd.isna(out.loc[2, "is_metro"]) and pd.isna(out.loc[2, "rucc_2023"])
    assert out.loc[1, "nonmetro_adjacent_to_metro"] == True  # noqa: E712 (pandas boolean)
    assert pd.isna(out.loc[0, "nonmetro_adjacent_to_metro"])


def test_rucc_transform_rejects_flag_description_conflict():
    wide = pd.DataFrame({"FIPS": ["01001"], "State": ["AL"], "County_Name": ["A"], "Population_2020": ["1"],
                         "RUCC_2023": ["1"], "Description": ["Nonmetro - mislabelled"]})
    with pytest.raises(ValueError):
        rucc.transform(wide, _codes())


def test_svi_clean_converts_minus999_and_counts():
    cols = svi.keep_columns()
    row = {c: 0.5 for c in cols}
    row.update({"ST": "01", "STATE": "Alabama", "ST_ABBR": "AL", "STCNTY": "01001", "COUNTY": "A",
                "FIPS": "1001", "LOCATION": "A, Alabama"})
    df = pd.DataFrame([row, dict(row, FIPS="01003", RPL_THEMES=-999.0, EP_NOVEH=-999.0)])
    out, counts = svi.clean(df)
    assert out["FIPS"].tolist() == ["01001", "01003"]     # zero-padded
    assert np.isnan(out.loc[1, "RPL_THEMES"]) and np.isnan(out.loc[1, "EP_NOVEH"])
    assert counts["__n_cells_equal_minus999"] == 2
    assert counts["RPL_THEMES"] == {"-999": 1}


def test_svi_keep_columns_cover_all_themes_and_components():
    cols = set(svi.keep_columns())
    for t in (1, 2, 3, 4):
        assert {f"RPL_THEME{t}", f"SPL_THEME{t}", f"F_THEME{t}"} <= cols
    assert {"RPL_THEMES", "SPL_THEMES", "F_TOTAL"} <= cols
    n_vars = sum(len(v) for v in svi.VARIABLES.values())
    assert n_vars == 16
    for vs in svi.VARIABLES.values():
        for v in vs:
            assert {f"E_{v}", f"EP_{v}", f"MP_{v}", f"EPL_{v}"} <= cols


# ============================================================================ data tests
def _table(name: str) -> pd.DataFrame:
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p)


def _canonical_counties() -> set[str]:
    g = _table("geographies")
    return set(g.loc[(g["geo_level"] == "county") & (~g["ct_legacy"]), "geo_id"])


CT_REGIONS = {f"09{c}" for c in ("110", "120", "130", "140", "150", "160", "170", "180", "190")}


def _check_provenance(df: pd.DataFrame, evidence_type: str) -> None:
    for c in PROVENANCE_COLUMNS:
        assert c in df.columns, c
    assert (df["data_layer"] == "geographic").all()
    assert (df["evidence_type"] == evidence_type).all()
    assert df["retrieved_at"].str.match(r"\d{4}-\d{2}-\d{2}T").all()


@pytest.mark.data
def test_geo_vulnerability_table():
    d = _table("geo_vulnerability")
    _check_provenance(d, "composite_index")
    assert len(d) == 3144
    assert d["geo_id"].is_unique and d["geo_id"].str.fullmatch(r"\d{5}").all()
    assert not d["geo_id"].str.startswith("72").any()          # PR ranked separately
    for c in ["rpl_themes", "rpl_theme1", "rpl_theme2", "rpl_theme3", "rpl_theme4"]:
        assert d[c].notna().all() and d[c].between(0, 1).all()
    assert set(d["geo_id"]) <= _canonical_counties()
    assert CT_REGIONS <= set(d["geo_id"])
    assert (d["source_geographic_resolution"] == "county").all()
    assert d["provenance_notes"].str.contains("NOT disease prevalence").all()
    # composite components are retained (guardrail: every score keeps its components)
    for c in ["spl_theme1", "epl_pov150", "ep_uninsur", "mp_uninsur", "f_total"]:
        assert c in d.columns


@pytest.mark.data
def test_geo_vulnerability_puerto_rico_is_separate():
    d = _table("geo_vulnerability_puerto_rico")
    assert len(d) == 78 and d["geo_id"].str.startswith("72").all()
    assert (d["svi_ranking_universe"] == "Puerto Rico municipios only").all()
    assert d["provenance_notes"].str.contains("NOT comparable").all()


@pytest.mark.data
def test_geo_context_acs_county_state():
    d = _table("geo_context__acs")
    _check_provenance(d, "survey_estimate")
    cty, st = d[d["geo_level"] == "county"], d[d["geo_level"] == "state"]
    assert len(cty) == 3222 and len(st) == 52
    assert not d.duplicated(["geo_id", "geo_level"]).any()
    assert (cty["source_geographic_resolution"] == "county").all()
    assert (st["source_geographic_resolution"] == "state").all()
    assert set(cty["geo_id"]) <= _canonical_counties()
    assert CT_REGIONS <= set(cty["geo_id"])
    assert not cty["geo_id"].isin([f"090{i:02d}" for i in range(1, 16, 2)]).any()
    assert cty["total_population"].sum() == st["total_population"].sum()
    pct = [c for c in d.columns if c.startswith("pct_") and not c.endswith("_moe")]
    assert len(pct) == 7
    for c in pct:
        assert d[c].dropna().between(0, 100).all(), c
        assert (d[f"{c}_moe"].dropna() >= 0).all(), c
        assert cty[c].isna().sum() == 0, c
    assert cty["median_household_income"].isna().sum() <= 5
    assert (cty["pop_density_per_km2"] > 0).all()
    # shared universes are identical estimates
    assert (cty["with_disability_universe"] == cty["uninsured_universe"]).all()
    assert (cty["households_no_vehicle_universe"] == cty["households_broadband_universe"]).all()
    assert (d["acs_vintage"] == d["acs_vintage"].iloc[0]).all()


@pytest.mark.data
def test_acs_consistent_with_svi_inputs():
    """Derived ACS shares track SVI 2022's independently published ACS inputs (earlier period)."""
    a = _table("geo_context__acs")
    s = _table("geo_vulnerability")
    m = a[a["geo_level"] == "county"].merge(s, on="geo_id")
    assert len(m) == 3144
    rho = lambda x, y: m[x].corr(m[y], method="spearman")  # noqa: E731
    assert rho("pct_uninsured", "ep_uninsur") > 0.85
    assert rho("pct_with_disability", "ep_disabl") > 0.85
    assert rho("pct_age_65_plus", "ep_age65") > 0.9
    assert rho("pct_households_no_vehicle", "ep_noveh") > 0.8
    assert rho("pct_households_broadband", "ep_noint") < -0.85


@pytest.mark.data
def test_acs_moes_close_to_published_subject_table_moes():
    """Approximated share MOEs stay near the replicate-weight MOEs SVI 2022 took from S/DP tables.

    Summing 12 B01001 age bands gave a 65+ MOE ~5x the published S0101_C02_030M; the single-cell B09020_001
    numerator brings it back near 1. Different ACS periods, so the band is wide.
    """
    a = _table("geo_context__acs")
    s = _table("geo_vulnerability")
    m = a[a["geo_level"] == "county"].merge(s, on="geo_id")
    for ours, published in [("pct_age_65_plus_moe", "mp_age65"), ("pct_uninsured_moe", "mp_uninsur"),
                            ("pct_with_disability_moe", "mp_disabl"), ("pct_households_no_vehicle_moe", "mp_noveh"),
                            ("pct_households_broadband_moe", "mp_noint")]:
        ratio = (m[ours] / m[published].where(m[published] > 0)).median()
        assert 0.6 < ratio < 1.6, (ours, ratio)


@pytest.mark.data
def test_acs_age_numerators_match_b01001_band_sums():
    import json
    p = RAW / "census_acs" / "build_stats.json"
    if not p.exists():
        pytest.skip("census_acs not built")
    checks = json.loads(p.read_text())["numerator_checks"]
    assert set(checks) == {"age_18_64", "age_65_plus"}
    for v in checks.values():
        assert v["rows_compared"] > 30000 and v["rows_differing"] == 0


@pytest.mark.data
def test_geo_context_zcta():
    z = _table("geo_context_zcta__acs")
    _check_provenance(z, "survey_estimate")
    assert len(z) > 33000 and z["zcta"].is_unique and z["zcta"].str.fullmatch(r"\d{5}").all()
    cz = _table("zcta_centroids")
    assert set(z["zcta"]) <= set(cz["zcta"])
    assert (z["source_geographic_resolution"] == "zcta").all()
    # zero-population ZCTAs have no shares, never a fabricated 0
    zero = z["total_population"] == 0
    assert z.loc[zero, "pct_uninsured"].isna().all()
    assert z["pct_uninsured"].dropna().between(0, 100).all()


@pytest.mark.data
def test_geo_context_rucc():
    r = _table("geo_context__rucc")
    _check_provenance(r, "composite_index")
    assert len(r) == 3235 and r["geo_id"].is_unique
    assert set(r["geo_id"]) == _canonical_counties()
    coded = r["rucc_2023"].notna()
    assert r.loc[coded, "rucc_2023"].between(1, 9).all()
    assert (~coded).sum() == 2 and r.loc[~coded, "geo_id"].str.startswith("60").all()
    assert (r.loc[coded, "is_metro"] == (r.loc[coded, "rucc_2023"] <= 3)).all()
    assert r.loc[coded & r["is_metro"].fillna(False), "nonmetro_adjacent_to_metro"].isna().all()
    assert CT_REGIONS <= set(r["geo_id"])


@pytest.mark.data
@pytest.mark.parametrize("source_id", ["cdc_svi", "census_acs", "usda_rucc"])
def test_registry_entries_complete(source_id):
    p = RAW / source_id / "registry_entry.yaml"
    if not p.exists():
        pytest.skip(f"{source_id} not built")
    e = yaml.safe_load(p.read_text())
    assert not [f for f in REQUIRED_FIELDS if f not in e]
    assert e["person_level"] is False and e["geographic"] is True
    assert e["true_participant_linkage_across_modalities"] is False
    assert (RAW / source_id / "MANIFEST.json").exists()


@pytest.mark.data
@pytest.mark.parametrize("source_id", ["cdc_svi", "census_acs", "usda_rucc", "census_geography"])
def test_data_audits_follow_template(source_id):
    p = RAW / source_id / "DATA_AUDIT.md"
    if not p.exists():
        pytest.skip(f"{source_id} audit not written")
    text = p.read_text()
    for heading in ["## Files / endpoints retrieved", "## Key variables", "## Missingness", "## Linkage strategy",
                    "## Limitations and caveats", "## Processed outputs", "## Reproduce"]:
        assert heading in text, heading
    for field in ["source_id", "Retrieval date", "License", "Sample size", "True participant linkage"]:
        assert field in text
    assert "{" + "s[" not in text        # no unrendered template fields
