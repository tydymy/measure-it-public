"""Tests for measure_it.ingestion.hrsa (HRSA health-center sites, HPSA and MUA/P designations)."""
from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import yaml
from shapely.geometry import box

from measure_it.config import RAW
from measure_it.ingestion import hrsa
from measure_it.provenance import check_provenance
from measure_it.registry import REQUIRED_FIELDS
from measure_it.store import read_table, table_exists


# --------------------------------------------------------------------------- pure-function unit tests

def test_clean_str_maps_placeholders_to_na():
    s = pd.Series([" Active ", "N/A", "XX", "XXXXX", "", "Not Determined", "06037"])
    out = hrsa.clean_str(s)
    assert out.tolist()[0] == "Active"
    assert out.iloc[1:6].isna().all()
    assert out.iloc[6] == "06037"


def test_zip5_extracts_five_digits():
    out = hrsa.zip5(pd.Series(["94952-3388", "02115", "", "ABCDE", " 99501"]))
    assert out.iloc[0] == "94952" and out.iloc[1] == "02115" and out.iloc[4] == "99501"
    assert out.iloc[2:4].isna().all()


def test_npi_check_digit():
    assert hrsa.npi_is_valid("1234567893")      # CMS worked example
    assert not hrsa.npi_is_valid("1234567890")  # wrong check digit
    assert not hrsa.npi_is_valid("123456789")
    assert not hrsa.npi_is_valid(None)
    out = hrsa.clean_npi(pd.Series(["1234567893", "1234567890", "", "abc"]))
    assert out.iloc[0] == "1234567893" and out.iloc[1:].isna().all()


def test_component_county_rules():
    geo = pd.Series(["09003406100", "06037101110", "POINT", "46113", "POINT"])
    pub = pd.Series(["XXXXX", "06037", "XXXXX", "XXXXX", "12086"])
    typ = pd.Series(["CT", "CT", "POINT", "SCTY", "POINT"])
    out = hrsa.component_county(geo, pub, typ)
    assert out.iloc[0] == "09003"          # legacy CT tract -> GEOID prefix (reported, not remapped)
    assert out.iloc[1] == "06037"
    assert pd.isna(out.iloc[2])            # facility with undetermined county stays NA
    assert out.iloc[3] == "46113"
    assert out.iloc[4] == "12086"


def test_hpsa_category():
    dt = pd.Series(["Geographic HPSA", "High Needs Geographic HPSA", "HPSA Population", "Rural Health Clinic",
                    "Federally Qualified Health Center"])
    ct = pd.Series(["SCTY", "CT", "CSD", "POINT", "POINT"])
    assert hrsa.hpsa_category(dt, ct).tolist() == [
        "geographic", "geographic_high_needs", "population", "facility", "facility"]


def test_component_geography_uses_type_sets():
    # 'SCTY' contains the substring 'CT'; a mixed whole-county + tract designation is not "census_tract"
    s = pd.Series(["POINT", "SCTY", "CT", "CSD", "CSD;CT", "CT;SCTY", "CSD;SCTY", "CT;POINT"])
    assert hrsa.component_geography(s).tolist() == [
        "facility_point", "whole_county", "census_tract", "county_subdivision", "mixed_subcounty",
        "mixed_county_and_subcounty", "mixed_county_and_subcounty", "mixed"]


def test_bridge_1to1_renames_only_renames():
    """Shannon SD / Wade Hampton AK are 2015 renames (same boundary): placed on the 2024 code. The 02261 split and the
    51515 merger are boundary changes and stay as published (unplaced candidates)."""
    s = pd.Series(["02270", "46113", "02261", "51515", "06073", None], dtype="string")
    out = hrsa.bridge_1to1_renames(s)
    assert out.tolist()[:5] == ["02158", "46102", "02261", "51515", "06073"] and pd.isna(out.iloc[5])
    assert set(hrsa.FIPS_RENAMES_1TO1) <= set(hrsa.RETIRED_COUNTY_SUCCESSORS)
    assert all(hrsa.RETIRED_COUNTY_SUCCESSORS[k] == (v,) for k, v in hrsa.FIPS_RENAMES_1TO1.items())


def test_legacy_county_candidates_overlap_and_retired():
    # two 2024 "regions" side by side; legacy county 09001 straddles both, 09003 lies in one,
    # 09005 touches region B only by a sliver (< 1% of its area)
    regions = gpd.GeoDataFrame({"geo_id": ["09110", "09120"]},
                               geometry=[box(-73.0, 41.0, -72.5, 41.5), box(-72.5, 41.0, -72.0, 41.5)], crs=4326)
    legacy = gpd.GeoDataFrame({"geo_id": ["09001", "09003", "09005"]},
                              geometry=[box(-72.8, 41.1, -72.2, 41.4), box(-72.9, 41.1, -72.6, 41.4),
                                        box(-72.999, 41.1, -72.499, 41.4)], crs=4326)
    c2024 = {"09110", "09120", "02063", "02066", "02158"}
    cand = hrsa.legacy_county_candidates(c2024, boundaries=regions, legacy_ct=legacy)
    assert cand["09001"] == ("09110", "09120")
    assert cand["09003"] == ("09110",)
    assert cand["09005"] == ("09110",)          # sliver below LEGACY_CT_MIN_SHARE ignored
    assert cand["02261"] == ("02063", "02066")  # Valdez-Cordova split
    assert cand["02270"] == ("02158",)          # Wade Hampton -> Kusilvak
    assert "51515" not in cand                  # successor 51019 not in this c2024 -> dropped


def test_unplaced_legacy_counts_and_undetermined_coverage():
    comp = pd.DataFrame([
        ("M1", "Designated", "09001"), ("M1", "Designated", "09001"),   # one designation, two components
        ("M2", "Designated", "09003"),
        ("M3", "Proposed For Withdrawal", "09003"),
        ("M4", "Designated", "06037"),                                  # placed normally: not counted
    ], columns=["mua_id", "status", "county_fips"])
    cand = {"09001": ("09110", "09120"), "09003": ("09110",)}
    out = hrsa.unplaced_legacy_counts(comp, "mua_id", "status", "mua", cand).set_index("county_fips")
    assert out.loc["09110", "mua_n_unplaced_legacy_candidates"] == 2
    assert out.loc["09120", "mua_n_unplaced_legacy_candidates"] == 1
    assert out.loc["09110", "mua_n_unplaced_legacy_candidates_pfw"] == 1
    assert "06037" not in out.index
    cov = hrsa.mark_undetermined(pd.Series(["none", "none", "partial_county"]), pd.Series([2, 0, 3]))
    assert cov.tolist() == [hrsa.UNDETERMINED, "none", "partial_county"]


def test_coverage_category_precedence():
    assert hrsa.coverage_category(True, True, True, True) == "whole_county_geographic"
    assert hrsa.coverage_category(False, True, True, True) == "whole_county_population"
    assert hrsa.coverage_category(False, False, True, True) == "partial_county"
    assert hrsa.coverage_category(False, False, False, True) == "facility_only"
    assert hrsa.coverage_category(False, False, False, False) == "none"


def _comp(rows):
    return pd.DataFrame(rows, columns=["hpsa_id", "hpsa_status", "designation_type", "designation_category",
                                       "hpsa_score", "component_type", "county_fips"])


def test_summarize_hpsa_by_county_synthetic():
    comp = _comp([
        # county A: whole-county geographic (score 12) + an FQHC facility (score 20)
        ("G1", "Designated", "Geographic HPSA", "geographic", 12, "SCTY", "00001"),
        ("F1", "Designated", "Federally Qualified Health Center", "facility", 20, "POINT", "00001"),
        # county B: one population HPSA with three tract components (score counted once) + a PFW
        ("P1", "Designated", "HPSA Population", "population", 16, "CT", "00002"),
        ("P1", "Designated", "HPSA Population", "population", 16, "CT", "00002"),
        ("P1", "Designated", "HPSA Population", "population", 16, "CSD", "00002"),
        ("W1", "Proposed For Withdrawal", "Geographic HPSA", "geographic", 25, "SCTY", "00002"),
        # county C: facility only; county D: only proposed-for-withdrawal
        ("R1", "Designated", "Rural Health Clinic", "facility", 7, "POINT", "00003"),
        ("W2", "Proposed For Withdrawal", "HPSA Population", "population", 9, "CT", "00004"),
        # a component with no county is ignored
        ("X1", "Designated", "Correctional Facility", "facility", 3, "POINT", None),
    ])
    out = hrsa.summarize_hpsa_by_county(comp, "pc").set_index("county_fips")
    a, b, c, d = (out.loc[k] for k in ["00001", "00002", "00003", "00004"])
    assert a["pc_whole_county_geographic"] and a["pc_coverage"] == "whole_county_geographic"
    assert a["pc_n_designations"] == 2 and a["pc_max_score"] == 20 and a["pc_mean_score"] == 16.0
    assert a["pc_max_score_geographic"] == 12 and a["pc_max_score_facility"] == 20
    assert b["pc_partial_county"] and not b["pc_whole_county"] and b["pc_coverage"] == "partial_county"
    assert b["pc_n_designations"] == 1 and b["pc_mean_score"] == 16.0   # PFW score 25 not used
    assert b["pc_n_proposed_for_withdrawal"] == 1
    assert c["pc_coverage"] == "facility_only" and c["pc_n_facility"] == 1
    assert d["pc_coverage"] == "none" and d["pc_n_designations"] == 0 and d["pc_n_proposed_for_withdrawal"] == 1
    assert pd.isna(d["pc_max_score"])
    assert set(out.index) == {"00001", "00002", "00003", "00004"}
    assert out.loc["00001", "pc_hpsa_ids"] == "F1;G1"


def test_summarize_mua_by_county_ignores_missing_and_zero_imu():
    comp = pd.DataFrame([
        ("M1", "Designated", "Medically Underserved Area", 55.0, "SCTY", "00001"),
        ("M2", "Designated", "Medically Underserved Area", 0.0, "CT", "00001"),
        ("M3", "Designated", "Medically Underserved Population-Governor's Exception", np.nan, "CT", "00002"),
        ("M4", "Proposed For Withdrawal", "Medically Underserved Area", 40.0, "CT", "00003"),
    ], columns=["mua_id", "status", "designation_type", "imu_score", "component_type", "county_fips"])
    out = hrsa.summarize_mua_by_county(comp).set_index("county_fips")
    assert out.loc["00001", "mua_min_imu_score"] == 55.0
    assert out.loc["00001", "mua_coverage"] == "whole_county" and out.loc["00001", "mua_n_designations"] == 2
    assert pd.isna(out.loc["00002", "mua_min_imu_score"]) and out.loc["00002", "mua_coverage"] == "partial_county"
    assert out.loc["00003", "mua_coverage"] == "none" and out.loc["00003", "mua_n_proposed_for_withdrawal"] == 1


def test_point_in_county_within_nearest_outside():
    b = gpd.GeoDataFrame({"geo_id": ["00001", "00002"]},
                         geometry=[box(-100.0, 40.0, -99.0, 41.0), box(-99.0, 40.0, -98.0, 41.0)], crs=4326)
    lat = pd.Series([40.5, 40.5, 40.5, np.nan], index=[10, 11, 12, 13])
    lon = pd.Series([-99.5, -97.99, -90.0, np.nan], index=[10, 11, 12, 13])  # -97.99 is ~0.85 km east of the box
    out = hrsa.point_in_county(lat, lon, b)
    assert out.loc[10, "county_fips_pip"] == "00001" and out.loc[10, "pip_method"] == "within_polygon"
    assert out.loc[11, "county_fips_pip"] == "00002" and out.loc[11, "pip_method"] == "nearest_polygon_within_5km"
    assert pd.isna(out.loc[12, "county_fips_pip"]) and out.loc[12, "pip_method"] == "outside_all_2024_counties"
    assert out.loc[13, "pip_method"] == "no_coordinates"


# --------------------------------------------------------------------------- data tests on processed outputs

def _need(name):
    if not table_exists(name):
        pytest.skip(f"{name} not built; run `uv run python -m measure_it.ingestion.hrsa`")
    return read_table(name)


def _counties_2024():
    g = read_table("geographies")
    return set(g.loc[(g["geo_level"] == "county") & ~g["ct_legacy"], "geo_id"])


@pytest.mark.data
def test_facilities_table():
    f = _need("facilities__hrsa")
    check_provenance(f, "facilities__hrsa")
    raw = pd.read_csv(RAW / hrsa.HC_SOURCE / hrsa.HC_FILES["sites"], dtype=str, keep_default_na=False)
    assert len(f) == len(raw) > 10_000
    assert f["facility_id"].is_unique
    assert (f["data_layer"] == "facility").all()
    assert f["is_look_alike"].sum() == (raw["Health Center Type"].str.contains("Look-Alike")).sum()
    # coordinates: plausible and resolution matches geocode method
    src = f["geocode_method"] == "source_address_geocode"
    assert src.mean() > 0.99
    assert f.loc[src, "lat"].between(-20, 72).all() and f.loc[src, "lon"].between(-180, 180).all()
    assert (f.loc[src, "source_geographic_resolution"] == "point").all()
    assert (f.loc[f["geocode_method"] == "zip5_as_zcta_internal_point", "source_geographic_resolution"]
            == "zcta_centroid").all()
    assert f.loc[f["geocode_method"] == "not_geocoded", "lat"].isna().all()
    # county: every assigned county is a 2024 county; published FIPS agrees with PIP almost always
    c2024 = _counties_2024()
    assigned = f["county_fips_method"] != "not_assigned"
    assert f.loc[assigned, "county_fips"].isin(c2024).all()
    assert f.loc[~assigned, "county_fips"].isna().all()
    agree = f["county_fips_pip_agrees"].dropna()
    assert agree.mean() > 0.99
    # linked facility HPSA scores are valid HPSA scores and point at in-effect designations
    des = _need("hpsa_designations")
    ids = set(des.loc[des["discipline_code"] == "pc", "hpsa_id"])
    linked = f["pc_facility_hpsa_id"].dropna()
    assert len(linked) > 0.5 * len(f) and linked.isin(ids).all()
    assert f["pc_facility_hpsa_score"].dropna().between(0, 26).all()


@pytest.mark.data
def test_hpsa_designations_table():
    d = _need("hpsa_designations")
    check_provenance(d, "hpsa_designations")
    assert d["hpsa_record_id"].is_unique
    assert set(d["hpsa_status"]) <= set(hrsa.IN_EFFECT_STATUSES)
    assert set(d["discipline_code"]) == {"pc", "mh", "dh"}
    assert d["hpsa_score"].between(0, 26).all()
    fac = d["designation_category"] == "facility"
    assert (d.loc[fac, "component_types"] == "POINT").all()
    assert not d.loc[~fac, "component_types"].str.contains("POINT").any()
    # independent recount from the raw CSV: Designated PC IDs
    raw = pd.read_csv(RAW / hrsa.HPSA_SOURCE / hrsa.HPSA_FILES["hpsa_pc"], dtype=str, keep_default_na=False,
                      usecols=["HPSA ID", "HPSA Status"])
    n_raw = raw.loc[raw["HPSA Status"] == "Designated", "HPSA ID"].nunique()
    assert n_raw == int((d["is_designated"] & (d["discipline_code"] == "pc")).sum())


@pytest.mark.data
def test_hpsa_components_reference_designations():
    c = _need("hpsa_components")
    d = _need("hpsa_designations")
    check_provenance(c, "hpsa_components")
    assert c["component_record_id"].is_unique
    keys = set(zip(d["discipline_code"], d["hpsa_id"]))
    assert set(zip(c["discipline_code"], c["hpsa_id"])) <= keys
    tract = c["component_type"] == "CT"
    assert c.loc[tract, "component_geo_id"].str.fullmatch(r"\d{11}").all()


@pytest.mark.data
def test_geo_context_hpsa_table():
    g = _need("geo_context__hpsa")
    check_provenance(g, "geo_context__hpsa")
    assert g["county_fips"].is_unique
    assert set(g["county_fips"]) == _counties_2024()
    assert (g["source_geographic_resolution"] == "county").all()
    for p in ["pc", "mh", "dh"]:
        n_ids = g[f"{p}_hpsa_ids"].fillna("").map(lambda s: len(s.split(";")) if s else 0)
        assert (n_ids == g[f"{p}_n_designations"]).all()
        assert (g[f"{p}_any_designation"] == (g[f"{p}_n_designations"] > 0)).all()
        assert (g.loc[g[f"{p}_whole_county_geographic"], f"{p}_coverage"] == "whole_county_geographic").all()
        assert (g.loc[g[f"{p}_coverage"] == "none", f"{p}_n_designations"] == 0).all()
        assert g[f"{p}_max_score"].dropna().between(0, 26).all()
    # a county never reads "none" while a Designated designation on legacy/retired county codes may cover it
    for p in ["pc", "mh", "dh", "mua"]:
        und = g[f"{p}_coverage"] == hrsa.UNDETERMINED
        assert (g.loc[und, f"{p}_n_designations"] == 0).all()
        assert (g.loc[und, f"{p}_n_unplaced_legacy_candidates"] > 0).all()
        assert not ((g[f"{p}_coverage"] == "none") & (g[f"{p}_n_unplaced_legacy_candidates"] > 0)).any()
    # Connecticut: every Designated MUA/P in the raw file is on legacy county/tract codes, so no planning
    # region may claim "no MUA/P"
    mr = pd.read_csv(RAW / hrsa.HPSA_SOURCE / hrsa.HPSA_FILES["mua"], dtype=str, keep_default_na=False,
                     usecols=["MUA/P Status Description", "MUA/P Area Code"])
    ct_legacy = mr[(mr["MUA/P Status Description"] == "Designated")
                   & mr["MUA/P Area Code"].str.match(r"^090(0[1-9]|1[0-5])")]
    if len(ct_legacy):
        ct = g[g["county_fips"].str.startswith("09")]
        assert (ct["mua_coverage"] != "none").all()
    # cross-check county flags against components (Designated PC, 2024 counties)
    c = _need("hpsa_components")
    pc = c[(c["discipline_code"] == "pc") & (c["hpsa_status"] == "Designated") & c["county_in_2024_set"]]
    assert set(pc["county_fips"]) == set(g.loc[g["pc_any_designation"], "county_fips"])
    whole = set(pc.loc[(pc["component_type"] == "SCTY"), "county_fips"])
    assert whole == set(g.loc[g["pc_whole_county"], "county_fips"])


@pytest.mark.data
def test_mua_designations_table():
    m = _need("mua_designations")
    check_provenance(m, "mua_designations")
    assert m["mua_id"].is_unique
    assert set(m["status"]) <= set(hrsa.IN_EFFECT_STATUSES)
    ge = m["designation_type_code"].str.endswith("-GE")
    assert m.loc[ge, "imu_score"].isna().all()
    assert m["imu_score"].dropna().between(0, 100).all()
    assert (m.loc[~ge, "imu_score"].fillna(-1) == m.loc[~ge, "imu_score_published"].fillna(-1)).all()
    assert (m.loc[m["component_types"] == "CT", "source_geographic_resolution"] == "census_tract").all()


@pytest.mark.data
def test_fips_renames_placed_on_2024_code():
    """MUA/P 00104 (Wade Hampton Census Area, whole county) is published on the retired code 02270: it is placed on
    Kusilvak 02158 (county_fips_bridge names the rename), so 02158's MUA coverage is not 'undetermined'."""
    m = _need("mua_designations")
    assert not m["county_fips_list"].str.contains(r"\b(?:02270|46113)\b").any()
    br = m[m["county_fips_bridge"] != ""]
    assert set(br["county_fips_bridge"]) <= {"02270->02158", "46113->46102"}
    for old, new in hrsa.FIPS_RENAMES_1TO1.items():
        assert (br.loc[br["county_fips_bridge"].str.contains(old), "county_fips_list"].str.contains(new)).all()
    c = _need("hpsa_components")
    assert not c["county_fips"].isin(list(hrsa.FIPS_RENAMES_1TO1)).any()
    assert (c.loc[c["county_fips_method"] == "fips_rename_1to1", "county_fips_published"]
            .isin(list(hrsa.FIPS_RENAMES_1TO1))).all()
    g = _need("geo_context__hpsa").set_index("county_fips")
    if (br["county_fips_bridge"] == "02270->02158").any():
        assert g.loc["02158", "mua_coverage"] != hrsa.UNDETERMINED and g.loc["02158", "mua_n_designations"] >= 1


@pytest.mark.data
def test_facility_points_without_coordinates_not_labelled_point():
    d = _need("hpsa_designations")
    pt = d["component_geography"] == "facility_point"
    assert (d.loc[pt & d["facility_lat"].isna(), "source_geographic_resolution"] != "point").all()
    assert (d.loc[pt & d["facility_lat"].notna(), "source_geographic_resolution"] == "point").all()
    c = _need("hpsa_components")
    n_pt = int((c["source_geographic_resolution"] == "point").sum())
    assert n_pt == int((pt & d["facility_lat"].notna()).sum())


@pytest.mark.data
def test_registry_entries_and_audits():
    for sid in (hrsa.HC_SOURCE, hrsa.HPSA_SOURCE):
        p = RAW / sid / "registry_entry.yaml"
        if not p.exists():
            pytest.skip("registry entries not written yet")
        e = yaml.safe_load(p.read_text())
        assert all(k in e for k in REQUIRED_FIELDS)
        assert e["status"] == "ingested" and e["person_level"] is False
        audit = (RAW / sid / "DATA_AUDIT.md").read_text()
        assert "Retrieval date" in audit and "Missingness" in audit and "Linkage strategy" in audit
        assert (RAW / sid / "MANIFEST.json").exists()
    f = read_table("facilities__hrsa", columns=["facility_id"])
    assert f"{len(f):,} sites" in (RAW / hrsa.HC_SOURCE / "DATA_AUDIT.md").read_text()
