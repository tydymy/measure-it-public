"""Tests for the CDC geographic sources: PLACES, Long COVID (HPS), Lyme surveillance.

Unit tests exercise pure functions on synthetic inputs. `@pytest.mark.data` tests check the
processed outputs written by
  uv run python -m measure_it.ingestion.cdc_places
  uv run python -m measure_it.ingestion.cdc_long_covid
  uv run python -m measure_it.ingestion.cdc_lyme
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
import yaml

from measure_it.config import RAW
from measure_it.ingestion import cdc_long_covid as lc
from measure_it.ingestion import cdc_lyme as ly
from measure_it.ingestion import cdc_places as pl
from measure_it.provenance import PROVENANCE_COLUMNS, check_provenance
from measure_it.registry import REQUIRED_FIELDS
from measure_it.store import read_table, table_exists

# --------------------------------------------------------------------------- PLACES units


def test_parse_release_name():
    assert pl.parse_release_name("PLACES: Local Data for Better Health, County Data, 2025 release") == ("county", 2025)
    assert pl.parse_release_name("PLACES: Local Data for Better Health, ZCTA Data 2023 release") == ("zcta", 2023)
    assert pl.parse_release_name("PLACES: Local Data for Better Health, Census Tract Data 2024 release") == (
        "census_tract", 2024)
    assert pl.parse_release_name("PLACES: County Data (GIS Friendly Format), 2025 release") is None
    assert pl.parse_release_name("500 Cities: Local Data for Better Health, 2019 release") is None


def test_latest_release_requires_county_and_zcta():
    rows = [("a", "PLACES: Local Data for Better Health, County Data, 2026 release"),  # no ZCTA for 2026
            ("b", "PLACES: Local Data for Better Health, County Data, 2025 release"),
            ("c", "PLACES: Local Data for Better Health, ZCTA Data, 2025 release"),
            ("d", "PLACES: Local Data for Better Health, County Data 2024 release"),
            ("e", "PLACES: Local Data for Better Health, ZCTA Data 2024 release")]
    assert pl.latest_release(rows) == {"year": 2025, "county": "b", "zcta": "c"}
    with pytest.raises(RuntimeError):
        pl.latest_release(rows[:1])


def test_relevance_flags_measures_and_whole_categories():
    assert pl.relevance("PHLTH", "HLTHSTAT") == (True, "health_status")
    assert pl.relevance("LPA", "RISKBEH") == (True, "risk_behavior")
    assert pl.relevance("NEWNEED", "SOCLNEED") == (True, "social_needs")   # future social-needs measure
    assert pl.relevance("NEWDIS", "DISABLT") == (True, "disability")
    assert pl.relevance("MAMMOUSE", "PREVENT") == (False, "other")


def test_long_covid_detector():
    m = pd.DataFrame({"measure_id": ["PHLTH", "LONGCOV", "X"],
                      "measure_label": ["Frequent physical distress", "Long COVID among adults", "Post-acute sequelae"],
                      "short_question_text": ["", "", ""]})
    assert pl.long_covid_measures(m) == ["LONGCOV", "X"]
    assert pl.long_covid_measures(m.iloc[:1]) == []


def test_ct_vintage():
    assert pl.ct_vintage(pd.Series(["09110", "06073"])) == "ct_planning_regions_09110_09190"
    assert pl.ct_vintage(pd.Series(["09001"])) == "ct_legacy_counties_09001_09015"
    assert pl.ct_vintage(pd.Series(["09001", "09190"])) == "mixed"
    assert pl.ct_vintage(pd.Series(["06073"])) == "no_ct_rows"


def test_brfss_years_from_description():
    desc = ("The 2025 release uses 2023 BRFSS data for 35 measures and 2022 BRFSS data for 5 measures "
            "(all teeth lost, ...)")
    assert pl.brfss_years_from_description(desc) == "BRFSS 2023 (35 measures); BRFSS 2022 (5 measures)"


def test_zcta_query_is_deterministic():
    q1, q2 = pl.zcta_query(["PHLTH", "GHLTH"]), pl.zcta_query(["GHLTH", "PHLTH"])
    assert q1 == q2 and q1["$where"] == "measureid in('GHLTH','PHLTH')"

# --------------------------------------------------------------------------- Long COVID units


def test_indicator_map():
    assert lc.indicator_map("Currently experiencing long COVID, as a percentage of all adults")[0] == \
        "lc_current_pct_all_adults"
    assert lc.indicator_map("Ever had COVID") is None
    with pytest.raises(ValueError):
        lc.indicator_map("Something new")
    assert len({v[0] for v in lc.INDICATORS.values()}) == 8


def test_hps_phase_label():
    assert lc.hps_phase_label("3.1", "08/23/2023") == "3.10"
    assert lc.hps_phase_label("3.1", "05/01/2021") == "3.1"
    assert lc.hps_phase_label("4.2", "08/20/2024") == "4.2"


def _raw_lc(**kw):
    base = {"Indicator": "Currently experiencing long COVID, as a percentage of all adults", "Group": "By State",
            "State": "California", "Subgroup": "California", "Phase": "4.2", "Time Period": "72",
            "Time Period Label": "Aug 20 - Sep 16, 2024", "Time Period Start Date": "08/20/2024",
            "Time Period End Date": "09/16/2024", "Value": "3.8", "LowCI": "2.8", "HighCI": "5.0",
            "Quartile range": "", "Quartile number": "1", "Suppression Flag": ""}
    base.update(kw)
    return base


def test_lc_tidy_drops_gap_and_context_rows_and_maps_states():
    raw = pd.DataFrame([
        _raw_lc(),
        _raw_lc(Group="National Estimate", State="United States", Subgroup="United States", Value="5.3"),
        _raw_lc(Phase="-1", Group="By Age", State="United States", Subgroup="18 - 29 years", Value=""),
        _raw_lc(Indicator="Ever had COVID"),
        _raw_lc(State="Texas", Subgroup="Texas", Value="", LowCI="", HighCI="", **{"Suppression Flag": "1"}),
    ]).astype(str).replace({"": np.nan})
    out, dropped = lc.tidy(raw, {"California": "06", "Texas": "48"})
    assert dropped == {"collection_gap_rows": 1, "ever_had_covid_rows": 1}
    assert len(out) == 3
    assert set(out["geo_level"]) == {"state", "national"}
    assert out.set_index("geo_name").loc["California", "geo_id"] == "06"
    assert out.set_index("geo_name").loc["United States", "geo_id"] == "US"
    tx = out[out["geo_name"] == "Texas"].iloc[0]
    assert tx["suppressed"] and math.isnan(tx["value"])
    assert (out["burden_evidence_level"] == "A").all() and (out["condition_id"] == "long_covid").all()
    assert out["county_fips"].isna().all()
    # 'US' is not a row in `geographies`; states are
    canon = out.set_index("geo_name")["in_canonical_geographies"]
    assert not canon["United States"] and canon["California"]
    with pytest.raises(ValueError):
        lc.tidy(raw, {"California": "06"})  # Texas has no FIPS -> fail loudly


def test_scan_text_for_long_covid():
    assert lc.scan_text_for_long_covid("FALLVAC ... COVID (1) Flu (2)")["long_covid_hits"] == 0
    txt = "Did you have any symptoms lasting 3 months or longer that you did not have prior to COVID?"
    assert lc.scan_text_for_long_covid(txt)["long_covid_hits"] >= 1

# --------------------------------------------------------------------------- Lyme units


def test_case_definition_period():
    assert ly.case_definition_period(2007).startswith("2001-2007")
    assert ly.case_definition_period(2008).startswith("2008-2021")
    assert ly.case_definition_period(2021).startswith("2008-2021")
    assert ly.case_definition_period(2022).startswith("2022+")


def test_normalize_status():
    assert ly.normalize_status("High Incidenc") == "High Incidence"
    assert ly.normalize_status("Low Incidence") == "Low Incidence"


def test_ct_placeholders_dropped_not_reallocated():
    raw = pd.DataFrame({
        "Ctyname": ["Fairfield County", "Capitol", "Autauga County"],
        "stname": ["Connecticut", "Connecticut", "Alabama"], "ststatus": ["High Incidenc"] * 2 + ["Low Incidence"],
        "stcode": ["9", "9", "1"], "ctycode": ["1", "110", "1"],
        "Cases2022": ["503", None, "0"], "cases2023": ["0", "507", "0"]})
    long, drops = ly.drop_placeholders(ly.parse_county_file(raw))
    got = {(r.fips, r.year): r.cases for r in long.itertuples()}
    assert got == {("09001", 2022): 503, ("01001", 2022): 0, ("09110", 2023): 507, ("01001", 2023): 0}
    assert drops["ct_legacy_placeholder_cells"] == 1 and drops["blank_cells"] == 1


def test_county_file_coverage():
    pu = pd.DataFrame({
        "Year": ["2012"] * 7,
        "State": ["CT", "CT", "CT", "NY", "NY", "DC", "Suppressed"],
        "FIPS": ["09001", "Unknown", "Suppressed", "36027", "Suppressed", "Unknown", "Suppressed"],
        "Case_status": ["Confirmed", "Confirmed", "Probable", "Confirmed", "Probable", "Confirmed", "Confirmed"],
        "Frequency": ["60", "50", "10", "30", "10", "8", "5"]})
    county = pd.DataFrame({"state_abbr": ["CT", "CT", "NY", "DC", "WY"], "year": [2012] * 5,
                           "cases": [55, 15, 43, 8, 2]})
    s = ly.county_file_coverage(pu, county).set_index("state_abbr")
    # CT: 120 public-use cases (50 with unknown county), 70 in the county file -> 50 missing
    assert s.loc["CT", "pu_state_cases"] == 120 and s.loc["CT", "pu_unknown_county_cases"] == 50
    assert s.loc["CT", "county_file_state_cases"] == 70
    assert s.loc["CT", "state_cases_missing_from_county_file_pct"] == pytest.approx(100 * 50 / 120)
    # county file above the public-use total (state-suppressed public-use rows) -> 0, never negative
    assert s.loc["NY", "state_cases_missing_from_county_file_pct"] == 0
    # unknown-county cases that ARE in the county file (DC-like) -> 0, not the unknown share
    assert s.loc["DC", "pu_unknown_county_cases"] == 8 and s.loc["DC", "state_cases_missing_from_county_file_pct"] == 0
    # no public-use row for the state -> blank, not 0
    assert math.isnan(s.loc["WY", "state_cases_missing_from_county_file_pct"])
    assert "Suppressed" not in s.index


def test_suspected_state_nonreporting():
    county = pd.DataFrame({
        "state_abbr": ["MN"] * 6 + ["WY"] * 3 + ["NY"] * 2 + ["VA"] * 3,
        "year": [2019, 2019, 2020, 2020, 2021, 2021, 2019, 2020, 2021, 2022, 2023, 2010, 2011, 2012],
        "cases": [1000, 528, 0, 0, 1900, 2, 1, 0, 3, 0, 500, 900, 0, 19]})
    f = ly.suspected_state_nonreporting(county).set_index(["state_abbr", "year"])["suspected_state_nonreporting"]
    assert f[("MN", 2020)] and not f[("MN", 2019)] and not f[("MN", 2021)]
    assert not f[("WY", 2020)]          # low-incidence state: a zero between small counts is plausible
    assert f[("NY", 2022)]              # edge year: only the adjacent years present are checked (2023 = 500)
    assert not f[("VA", 2011)]          # one adjacent year below 20 -> not flagged
    assert f.dtype == bool


def test_incidence_per_100k():
    r = ly.incidence_per_100k(pd.Series([5, 0, 3]), pd.Series([50_000, 1_000, 0]))
    assert r.iloc[0] == pytest.approx(10.0) and r.iloc[1] == 0 and math.isnan(r.iloc[2])

# --------------------------------------------------------------------------- data tests


def _need(name):
    if not table_exists(name):
        pytest.skip(f"{name} not built")
    return read_table(name)


REQUIRED_CONTEXT = {"GHLTH", "PHLTH", "MHLTH", "DEPRESSION", "LPA", "OBESITY", "SLEEP", "ACCESS2", "CHECKUP",
                    "DISABILITY", "COGNITION", "HEARING", "INDEPLIVE", "MOBILITY", "SELFCARE", "VISION",
                    "ARTHRITIS", "CASTHMA", "COPD", "DIABETES", "CHD", "STROKE", "BPHIGH", "CSMOKING"}


@pytest.mark.data
def test_places_county_output():
    d = _need("geo_context__cdc_places")
    check_provenance(d, "geo_context__cdc_places")
    assert (d["evidence_type"] == "modeled_small_area_estimate").all()
    assert set(d["source_geographic_resolution"]) == {"county", "national"}
    assert d["measure_id"].nunique() >= 40
    assert not d.duplicated(["geo_id", "measure_id", "data_value_type_id"]).any()
    assert set(d["data_value_type_id"]) == {"CrdPrv", "AgeAdjPrv"}
    assert pl.long_covid_measures(d[["measure_id", "measure_label", "short_question_text"]].drop_duplicates()) == []
    flagged = set(d.loc[d["relevant_measure"], "measure_id"])
    assert REQUIRED_CONTEXT <= flagged
    assert set(d.loc[d["category_id"] == "SOCLNEED", "measure_id"]) <= flagged
    cty = d[d["geo_level"] == "county"]
    assert cty["geo_id"].str.fullmatch(r"\d{5}").all()
    assert cty["in_canonical_geographies"].all()
    assert cty["geo_id"].nunique() > 3100
    ct = set(cty.loc[cty["geo_id"].str.startswith("09"), "geo_id"])
    assert ct == pl.CT_PLANNING_REGIONS
    v = d.dropna(subset=["value"])
    assert ((v["ci_low"] <= v["value"]) & (v["value"] <= v["ci_high"])).all()
    assert v["value"].between(0, 100).all()
    assert d.loc[d["suppressed"], "value"].isna().all()
    nat = d[d["geo_level"] == "national"]
    assert (nat["geo_id"] == "US").all() and (nat["geo_name"] == "United States").all()
    assert not nat["in_canonical_geographies"].any()


@pytest.mark.data
def test_places_zcta_output():
    z = _need("geo_context_zcta__cdc_places")
    check_provenance(z, "geo_context_zcta__cdc_places")
    assert (z["source_geographic_resolution"] == "zcta").all()
    assert z["geo_id"].str.fullmatch(r"\d{5}").all()
    assert (z["data_value_type_id"] == "CrdPrv").all()
    assert z["relevant_measure"].all()
    assert not z.duplicated(["geo_id", "measure_id"]).any()
    assert z["geo_id"].nunique() > 30000
    assert z["in_canonical_zctas"].all()
    assert z["measure_label"].notna().all()


@pytest.mark.data
def test_long_covid_output():
    d = _need("geo_condition_burden__cdc_long_covid")
    check_provenance(d, "geo_condition_burden__cdc_long_covid")
    assert set(d["geo_level"]) == {"national", "state"}
    assert "county" not in set(d["source_geographic_resolution"])
    assert (d["geo_level"] == d["source_geographic_resolution"]).all()
    assert d["county_fips"].isna().all()
    assert (d["condition_id"] == "long_covid").all() and (d["burden_evidence_level"] == "A").all()
    assert (d["evidence_type"] == "survey_estimate").all()
    assert set(d["measure_id"]) == {v[0] for v in lc.INDICATORS.values()}
    assert not d["measure_label"].eq("Ever had COVID").any()
    st = d[d["geo_level"] == "state"]
    assert st["geo_id"].str.fullmatch(r"\d{2}").all() and st["geo_id"].nunique() == 51
    assert (st["subgroup_type"] == "By State").all()           # subgroups are national only
    assert d.loc[d["suppressed"], "value"].isna().all()
    assert d.loc[~d["suppressed"], "value"].notna().all()
    assert not d.duplicated(["measure_id", "geo_id", "subgroup_type", "subgroup", "time_period_id"]).any()
    prim = d[d["primary_burden_measure"]]
    assert set(prim["geo_id"]) == set(st["geo_id"]) | {"US"}
    v = d.dropna(subset=["value", "ci_low", "ci_high"])
    assert ((v["ci_low"] <= v["value"]) & (v["value"] <= v["ci_high"])).all()
    assert d["period_end"].max() <= "2024-12-31"
    assert d.loc[d["geo_level"] == "state", "in_canonical_geographies"].all()
    assert not d.loc[d["geo_level"] == "national", "in_canonical_geographies"].any()


@pytest.mark.data
def test_lyme_output():
    d = _need("geo_condition_burden__cdc_lyme")
    check_provenance(d, "geo_condition_burden__cdc_lyme")
    assert (d["source_geographic_resolution"] == "county").all()
    assert (d["condition_id"] == "lyme_disease").all() and (d["burden_evidence_level"] == "A").all()
    assert (d["evidence_type"] == "surveillance_case_count").all()
    assert not d["condition_id"].eq("ptlds").any()
    assert set(d["measure_id"]) == {"lyme_reported_cases", "lyme_incidence_per_100k"}
    assert not d.duplicated(["geo_id", "year", "measure_id"]).any()
    assert d["year"].min() == 2001 and d["year"].max() == 2023
    counts = d[d["measure_id"] == "lyme_reported_cases"]
    rates = d[d["measure_id"] == "lyme_incidence_per_100k"]
    assert (counts["value"] >= 0).all() and (counts["value"] == counts["numerator"]).all()
    r = rates.dropna(subset=["value"])
    assert np.allclose(r["value"], r["numerator"] / r["denominator"] * 1e5)
    assert rates["value"].isna().mean() < 0.01
    # Connecticut: legacy counties through 2022, planning regions only in 2023, never both in a year
    ct = counts[counts["state_fips"] == "09"]
    assert set(ct.loc[ct["year"] < 2023, "geo_id"]) == ly.CT_LEGACY
    assert set(ct.loc[ct["year"] == 2023, "geo_id"]) == ly.CT_PLANNING
    assert (~counts.loc[counts["year"] >= 2022, "comparable_with_pre_2022"]).all()
    assert counts.loc[counts["year"] < 2022, "comparable_with_pre_2022"].all()
    assert counts.groupby("year")["geo_id"].nunique().min() > 3000
    # county-file coverage vs the public-use files: bounded, one value per state-year, Connecticut 2012 the known
    # heavy case (2,657 public-use cases, 1,255 in the county file)
    col = "state_cases_missing_from_county_file_pct"
    assert d[col].dropna().between(0, 100).all()
    ct12 = counts[(counts["state_fips"] == "09") & (counts["year"] == 2012)][col]
    assert ct12.nunique() == 1 and ct12.iloc[0] == pytest.approx(100 * (2657 - 1255) / 2657)
    assert not d.groupby(["state_fips", "year"])[col].nunique().gt(1).any()
    # high-volume states always have public-use rows
    assert counts.loc[counts["state_abbr"].isin(["NY", "PA", "WI"]), col].notna().all()
    # Minnesota reported 0 cases in 2020 (CDC: 2019-2020 data incomplete); flagged, and only a handful of state-years
    flagged = counts.loc[counts["suspected_state_nonreporting"], ["state_abbr", "year"]].drop_duplicates()
    assert ("MN", 2020) in set(map(tuple, flagged.to_numpy().tolist())) and len(flagged) <= 5
    assert (counts.loc[counts["suspected_state_nonreporting"], "value"] == 0).all()


@pytest.mark.data
@pytest.mark.parametrize("source_id", ["cdc_places", "cdc_long_covid", "cdc_lyme"])
def test_registry_and_audit(source_id):
    p = RAW / source_id / "registry_entry.yaml"
    if not p.exists():
        pytest.skip(f"{source_id} not ingested")
    entry = yaml.safe_load(p.read_text())
    assert not [f for f in REQUIRED_FIELDS if f not in entry]
    assert entry["person_level"] is False and entry["geographic"] is True
    assert entry["true_participant_linkage_across_modalities"] is False
    audit = (RAW / source_id / "DATA_AUDIT.md").read_text()
    assert "## Missingness" in audit and "## Linkage strategy" in audit and "{" not in audit.split("\n")[0]
    assert (RAW / source_id / "MANIFEST.json").exists()


def test_provenance_columns_constant():
    # guard against an upstream rename that would silently drop provenance from these outputs
    assert {"data_layer", "source_geographic_resolution", "evidence_type"} <= set(PROVENANCE_COLUMNS)
