"""HRSA health-center sites (facility layer) and shortage designations (geographic layer).

Sources (data.hrsa.gov -> Data Downloads; HRSA Data Warehouse, refreshed daily):

hrsa_health_centers
    Health_Center_Service_Delivery_and_LookAlike_Sites.csv (+ data dictionary XLSX).
    One row per Health Center Program (FQHC) or FQHC Look-Alike site in scope.
hrsa_hpsa
    BCD_HPSA_FCT_DET_{PC,MH,DH}.csv: Health Professional Shortage Areas, one row per
    designation *component* (census tract, county subdivision, whole county, or a
    facility point); SDMS_AUTO_HPSA_SITE_{PC,MH,DH}.csv: auto-HPSA site lists that
    link each health-center site (BPHC assigned number) to its facility HPSA;
    MUA_DET.csv: Medically Underserved Areas/Populations (MUA/P), one row per component.
    Data dictionaries: HPSA_DATAMART_METADATA.XLSX, DD_AUTO_HPSA_SITE_METADATA.XLSX,
    MUA_DATAMART_METADATA.XLSX.

Outputs (data/processed):
    facilities__hrsa   one row per HRSA health-center site (facility layer)
    hpsa_designations  one row per (discipline, HPSA ID) still in effect
                       (status Designated or Proposed For Withdrawal; Withdrawn dropped)
    hpsa_components    one row per in-effect HPSA component with its county
    mua_designations   one row per in-effect MUA/P designation
    geo_context__hpsa  one row per 2024 county: HPSA (PC/MH/DH) + MUA/P summary

Guardrails
    * A shortage designation is an administrative federal designation (provider-to-
      population ratio + need indicators), not a disease-burden measure.
    * County flags record the *presence* of a designation touching the county. A
      "partial" designation covers only some tracts/subdivisions (or a low-income
      sub-population) of the county; it is not a county-wide statement.
    * Health-center sites are an administrative listing of in-scope sites; the file
      carries no service-line or clinical-capability detail.
    * Legacy Connecticut geographies are reported as unmatched, never forced onto the
      2024 planning regions.

Reproduce: uv run python -m measure_it.ingestion.hrsa
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from ..config import PROCESSED, RAW, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..geography.crosswalk import zip_to_geo
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

HC_SOURCE = "hrsa_health_centers"
HPSA_SOURCE = "hrsa_hpsa"
BASE_URL = "https://data.hrsa.gov/DataDownload/DD_Files/"
LANDING_URL = "https://data.hrsa.gov/data/download"
PRODUCER = "ingestion.hrsa"

HC_FILES = {
    "sites": "Health_Center_Service_Delivery_and_LookAlike_Sites.csv",
    "sites_dictionary": "Health_Center_Service_Delivery_and_LookAlike_Sites_Data_Download_Metadata.xlsx",
}
HPSA_FILES = {
    "hpsa_pc": "BCD_HPSA_FCT_DET_PC.csv",
    "hpsa_mh": "BCD_HPSA_FCT_DET_MH.csv",
    "hpsa_dh": "BCD_HPSA_FCT_DET_DH.csv",
    "hpsa_dictionary": "HPSA_DATAMART_METADATA.XLSX",
    "auto_pc": "SDMS_AUTO_HPSA_SITE_PC.csv",
    "auto_mh": "SDMS_AUTO_HPSA_SITE_MH.csv",
    "auto_dh": "SDMS_AUTO_HPSA_SITE_DH.csv",
    "auto_dictionary": "DD_AUTO_HPSA_SITE_METADATA.XLSX",
    "mua": "MUA_DET.csv",
    "mua_dictionary": "MUA_DATAMART_METADATA.XLSX",
}
DISCIPLINES = {"pc": "Primary Care", "mh": "Mental Health", "dh": "Dental Health"}

IN_EFFECT_STATUSES = ("Designated", "Proposed For Withdrawal")
COUNTY_BOUNDARIES = PROCESSED / "county_boundaries_2024.geoparquet"
PIP_SNAP_METERS = 5000.0

# Values the HRSA extracts use for "no value" in identifier fields.
PLACEHOLDERS = {"", "N/A", "NA", "n/a", "na", "NONE", "None", "none", "0", "00000", "000000", "0000000000",
                "Not Applicable", "Not Determined", "XX", "XXX", "XXXXX", "Unknown"}

HPSA_CATEGORY = {
    "Geographic HPSA": "geographic",
    "High Needs Geographic HPSA": "geographic_high_needs",
    "HPSA Population": "population",
}
UNDETERMINED = "undetermined_legacy_geography"
COVERAGE_ORDER = ["whole_county_geographic", "whole_county_population", "partial_county", "facility_only",
                  UNDETERMINED, "none"]

# Census 2020 cartographic counties (read-only; shipped by the census_geography source) give the
# legacy Connecticut county polygons used to find which 2024 planning regions a legacy code overlaps.
LEGACY_CT_COUNTIES_ZIP = RAW / "census_geography" / "cb_2020_us_county_5m.zip"
LEGACY_CT_MIN_SHARE = 0.01  # overlap must be >= 1% of the legacy county's area (drops 1:5M slivers)
# Retired county codes still used by HRSA components -> the 2024 county/ies covering the same area
# (Census "Substantial Changes to Counties and County Equivalent Entities").
RETIRED_COUNTY_SUCCESSORS = {"02261": ("02063", "02066"),  # Valdez-Cordova split 2019 -> Chugach, Copper River
                             "02270": ("02158",),           # Wade Hampton renamed Kusilvak 2015
                             "46113": ("46102",),           # Shannon renamed Oglala Lakota 2015
                             "51515": ("51019",)}           # Bedford city merged into Bedford County 2013
# The 1:1 renames among them (new code, unchanged boundary) are identities, not overlaps: a component on the old code
# IS the 2024 county, so it is placed there (2026-09-24; the CMS MMD and CDC Lyme ingestions bridge the same two
# renames). `county_fips_published` keeps the code HRSA published; HPSA components get county_fips_method
# 'fips_rename_1to1' and MUA/P designations name the rename in `county_fips_bridge`. The boundary changes (02261 split,
# 51515 merger) stay unplaced candidates.
FIPS_RENAMES_1TO1 = {"02270": "02158", "46113": "46102"}


def bridge_1to1_renames(county_fips: pd.Series) -> pd.Series:
    """County codes with the two 1:1 renames (FIPS_RENAMES_1TO1) replaced by their 2024 code; others unchanged."""
    return county_fips.where(~county_fips.isin(list(FIPS_RENAMES_1TO1)), county_fips.map(FIPS_RENAMES_1TO1))


# --------------------------------------------------------------------------- helpers (pure)

def clean_str(s: pd.Series) -> pd.Series:
    """Strip whitespace; map HRSA placeholder strings to <NA>."""
    out = s.astype("string").str.strip()
    return out.mask(out.isin(PLACEHOLDERS))


def zip5(s: pd.Series) -> pd.Series:
    """First five digits of a ZIP / ZIP+4 string, else <NA>."""
    return s.astype("string").str.extract(r"^\s*(\d{5})")[0].astype("string")


def npi_is_valid(npi: str | None) -> bool:
    """NPI check digit (Luhn over '80840' + first 9 digits), per CMS NPI standard."""
    if npi is None or not isinstance(npi, str) or not re.fullmatch(r"\d{10}", npi):
        return False
    digits = [int(c) for c in "80840" + npi[:9]]
    total = 0
    for i, d in enumerate(reversed(digits)):
        if i % 2 == 0:  # rightmost payload digit is doubled (the check digit is excluded)
            d *= 2
            if d > 9:
                d -= 9
        total += d
    return (10 - total % 10) % 10 == int(npi[9])


def clean_npi(s: pd.Series) -> pd.Series:
    v = s.astype("string").str.strip()
    ok = v.map(lambda x: npi_is_valid(x) if isinstance(x, str) else False).astype(bool)
    return v.where(ok & (v != "1234567890"))


def parse_date(s: pd.Series, fmt: str) -> pd.Series:
    return pd.to_datetime(s.astype("string").str.strip().str[:10], format=fmt, errors="coerce")


def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype("string").str.replace(",", "", regex=False).str.strip(), errors="coerce")


def valid_county_fips(s: pd.Series) -> pd.Series:
    v = s.astype("string").str.strip()
    return v.where(v.str.fullmatch(r"\d{5}").fillna(False))


def hpsa_category(designation_type: pd.Series, component_type: pd.Series) -> pd.Series:
    """geographic / geographic_high_needs / population / facility."""
    cat = designation_type.map(HPSA_CATEGORY)
    cat = cat.where(component_type != "POINT", "facility")
    return cat.fillna("facility").astype("string")


def component_county(geo_id: pd.Series, county_col: pd.Series, component_type: pd.Series) -> pd.Series:
    """County FIPS of an HPSA/MUA component.

    Uses the published state+county FIPS when it is a 5-digit code; for area
    components whose county field is undetermined ('XXXXX', e.g. legacy
    Connecticut tracts) falls back to the first 5 characters of the tract / county-
    subdivision / county GEOID. Facility points with an undetermined county stay <NA>.
    """
    pub = valid_county_fips(county_col)
    prefix = geo_id.astype("string").str.strip().str[:5]
    prefix = prefix.where(prefix.str.fullmatch(r"\d{5}").fillna(False) & (component_type != "POINT"))
    return pub.fillna(prefix)


def component_geography(component_types: pd.Series) -> pd.Series:
    """Geography of a designation from its ';'-joined component types (CT/CSD/SCTY/POINT).

    Compares type *sets* (a substring test would let 'SCTY' match 'CT')."""
    sets = component_types.fillna("").astype(str).str.split(";").map(lambda t: frozenset(x for x in t if x))
    single = {frozenset({"POINT"}): "facility_point", frozenset({"SCTY"}): "whole_county",
              frozenset({"CT"}): "census_tract", frozenset({"CSD"}): "county_subdivision",
              frozenset({"CT", "CSD"}): "mixed_subcounty"}
    return sets.map(lambda t: single.get(t, "mixed_county_and_subcounty" if "SCTY" in t else "mixed")).astype("string")


def legacy_county_candidates(c2024: set[str], boundaries: gpd.GeoDataFrame | None = None,
                             legacy_ct: gpd.GeoDataFrame | None = None) -> dict[str, tuple[str, ...]]:
    """Non-2024 county code -> the 2024 counties whose area it overlaps (candidates, not a mapping).

    Retired codes use RETIRED_COUNTY_SUCCESSORS. Legacy Connecticut counties (09001-09015) use the
    geometric overlap of the Census 2020 county polygon with the 2024 planning regions (share of the
    legacy county's area >= LEGACY_CT_MIN_SHARE); if those polygons are unavailable every CT planning
    region is a candidate. Used only to flag counties whose value is undetermined, never to assign.
    """
    cand = {k: tuple(v for v in vs if v in c2024) for k, vs in RETIRED_COUNTY_SUCCESSORS.items()}
    ct_regions = tuple(sorted(c for c in c2024 if c.startswith("09")))
    if legacy_ct is None and LEGACY_CT_COUNTIES_ZIP.exists():
        lg = gpd.read_file(f"zip://{LEGACY_CT_COUNTIES_ZIP}")
        legacy_ct = lg.loc[lg["STATEFP"] == "09", ["GEOID", "geometry"]].rename(columns={"GEOID": "geo_id"})
    if legacy_ct is None or legacy_ct.empty:
        cand.update({f"09{c:03d}": ct_regions for c in range(1, 16, 2)})
        return {k: v for k, v in cand.items() if v}
    if boundaries is None:
        boundaries = gpd.read_parquet(COUNTY_BOUNDARIES)
    reg = boundaries.loc[boundaries["geo_id"].isin(ct_regions), ["geo_id", "geometry"]] \
        .rename(columns={"geo_id": "region"}).to_crs(5070)
    lg = legacy_ct[["geo_id", "geometry"]].to_crs(5070)
    lg = lg.assign(legacy_area=lg.area)
    ov = gpd.overlay(lg, reg, how="intersection", keep_geom_type=True)
    ov = ov[ov.area / ov["legacy_area"] >= LEGACY_CT_MIN_SHARE]
    for code, sub in ov.groupby("geo_id"):
        cand[str(code)] = tuple(sorted(sub["region"].astype(str).unique()))
    return {k: v for k, v in cand.items() if v}


def unplaced_legacy_counts(comp: pd.DataFrame, id_col: str, status_col: str, prefix: str,
                           candidates: dict[str, tuple[str, ...]]) -> pd.DataFrame:
    """Per candidate 2024 county: distinct designations (Designated / Proposed For Withdrawal) that have
    a component on a legacy or retired county code overlapping that county. comp.county_fips holds
    the component's county code as published (not remapped), except the two 1:1 renames (FIPS_RENAMES_1TO1), which
    are already on their 2024 code and so are never candidates."""
    x = comp[comp["county_fips"].isin(list(candidates))]
    x = x.assign(candidate_county=x["county_fips"].map(candidates)).explode("candidate_county")
    out = {}
    for suffix, status in [("", "Designated"), ("_pfw", "Proposed For Withdrawal")]:
        s = x[x[status_col] == status].drop_duplicates(["candidate_county", id_col])
        out[f"{prefix}_n_unplaced_legacy_candidates{suffix}"] = s.groupby("candidate_county").size()
    res = pd.DataFrame(out).fillna(0).astype("int64")
    res.index.name = "county_fips"
    return res.reset_index()


def coverage_category(whole_geo: bool, whole_pop: bool, partial: bool, facility: bool) -> str:
    if whole_geo:
        return "whole_county_geographic"
    if whole_pop:
        return "whole_county_population"
    if partial:
        return "partial_county"
    if facility:
        return "facility_only"
    return "none"


def _type_counts(values: pd.Series) -> str:
    vc = values.value_counts()
    return "; ".join(f"{k}={int(v)}" for k, v in sorted(vc.items()))


def summarize_hpsa_by_county(comp: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """County summary of HPSA components for one discipline.

    comp columns: hpsa_id, hpsa_status, designation_type, designation_category,
    hpsa_score, component_type (CT/CSD/SCTY/POINT), county_fips.
    Flags and scores use status == 'Designated'; 'Proposed For Withdrawal' rows are
    counted separately. Scores are taken once per designation (not per component).
    """
    comp = comp.dropna(subset=["county_fips"])
    d = comp[comp["hpsa_status"] == "Designated"]
    pfw = comp[comp["hpsa_status"] == "Proposed For Withdrawal"]
    area = d["designation_category"] != "facility"
    whole = d["component_type"] == "SCTY"
    geo_cat = d["designation_category"].isin(["geographic", "geographic_high_needs"])

    out = pd.DataFrame(index=pd.Index(sorted(set(comp["county_fips"])), name="county_fips"))
    out[f"{prefix}_whole_county"] = d[area & whole].groupby("county_fips").size().reindex(out.index).notna()
    out[f"{prefix}_whole_county_geographic"] = (
        d[area & whole & geo_cat].groupby("county_fips").size().reindex(out.index).notna())
    out[f"{prefix}_whole_county_population"] = (
        d[area & whole & (d["designation_category"] == "population")].groupby("county_fips").size()
        .reindex(out.index).notna())
    out[f"{prefix}_partial_county"] = (
        d[area & d["component_type"].isin(["CT", "CSD"])].groupby("county_fips").size().reindex(out.index).notna())
    out[f"{prefix}_facility_designation"] = (
        d[~area].groupby("county_fips").size().reindex(out.index).notna())

    per_des = d.drop_duplicates(["county_fips", "hpsa_id"])
    gd = per_des.groupby("county_fips")
    out[f"{prefix}_n_designations"] = gd.size().reindex(out.index).fillna(0).astype("int64")
    for cat_name, cats in [("geographic", ["geographic", "geographic_high_needs"]), ("population", ["population"]),
                           ("facility", ["facility"])]:
        sub = per_des[per_des["designation_category"].isin(cats)]
        out[f"{prefix}_n_{cat_name}"] = sub.groupby("county_fips").size().reindex(out.index).fillna(0).astype("int64")
        out[f"{prefix}_max_score_{cat_name}"] = sub.groupby("county_fips")["hpsa_score"].max().reindex(out.index)
    out[f"{prefix}_max_score"] = gd["hpsa_score"].max().reindex(out.index)
    out[f"{prefix}_mean_score"] = gd["hpsa_score"].mean().reindex(out.index).round(3)
    out[f"{prefix}_designation_types"] = gd["designation_type"].agg(_type_counts).reindex(out.index)
    out[f"{prefix}_hpsa_ids"] = gd["hpsa_id"].agg(lambda s: ";".join(sorted(s))).reindex(out.index)
    out[f"{prefix}_n_proposed_for_withdrawal"] = (
        pfw.drop_duplicates(["county_fips", "hpsa_id"]).groupby("county_fips").size()
        .reindex(out.index).fillna(0).astype("int64"))
    out[f"{prefix}_coverage"] = [
        coverage_category(a, b, c, e) for a, b, c, e in zip(
            out[f"{prefix}_whole_county_geographic"], out[f"{prefix}_whole_county_population"],
            out[f"{prefix}_partial_county"], out[f"{prefix}_facility_designation"])]
    return out.reset_index()


def summarize_mua_by_county(comp: pd.DataFrame) -> pd.DataFrame:
    """County summary of MUA/P components (status Designated; PFW counted separately).

    comp columns: mua_id, status, designation_type, imu_score, component_type, county_fips.
    IMU score: lower = more underserved (dictionary: must be <= 62.0 to qualify; Governor's
    exceptions carry no IMU). mua_min_imu_score ignores missing and exactly-zero IMU values.
    """
    comp = comp.dropna(subset=["county_fips"])
    d = comp[comp["status"] == "Designated"]
    pfw = comp[comp["status"] == "Proposed For Withdrawal"]
    out = pd.DataFrame(index=pd.Index(sorted(set(comp["county_fips"])), name="county_fips"))
    out["mua_whole_county"] = d[d["component_type"] == "SCTY"].groupby("county_fips").size().reindex(out.index).notna()
    out["mua_partial_county"] = (
        d[d["component_type"].isin(["CT", "CSD"])].groupby("county_fips").size().reindex(out.index).notna())
    per = d.drop_duplicates(["county_fips", "mua_id"])
    g = per.groupby("county_fips")
    out["mua_n_designations"] = g.size().reindex(out.index).fillna(0).astype("int64")
    out["mua_min_imu_score"] = per[per["imu_score"] > 0].groupby("county_fips")["imu_score"].min().reindex(out.index)
    out["mua_designation_types"] = g["designation_type"].agg(_type_counts).reindex(out.index)
    out["mua_ids"] = g["mua_id"].agg(lambda s: ";".join(sorted(s))).reindex(out.index)
    out["mua_n_proposed_for_withdrawal"] = (
        pfw.drop_duplicates(["county_fips", "mua_id"]).groupby("county_fips").size()
        .reindex(out.index).fillna(0).astype("int64"))
    out["mua_coverage"] = np.select([out["mua_whole_county"], out["mua_partial_county"]],
                                    ["whole_county", "partial_county"], "none")
    return out.reset_index()


def point_in_county(lat: pd.Series, lon: pd.Series, boundaries: gpd.GeoDataFrame | None = None) -> pd.DataFrame:
    """2024 county containing each point; unmatched points snap to the nearest county
    polygon within PIP_SNAP_METERS (the 1:5M cartographic boundaries are generalised and
    clipped to the shoreline). Returns county_fips_pip + pip_method, index-aligned."""
    if boundaries is None:
        boundaries = gpd.read_parquet(COUNTY_BOUNDARIES)
    b = boundaries[["geo_id", "geometry"]].to_crs(4326)
    out = pd.DataFrame({"county_fips_pip": pd.Series(pd.NA, index=lat.index, dtype="string"),
                        "pip_method": pd.Series("no_coordinates", index=lat.index, dtype="string")})
    ok = lat.notna() & lon.notna()
    if not ok.any():
        return out
    pts = gpd.GeoDataFrame(index=lat.index[ok], geometry=gpd.points_from_xy(lon[ok], lat[ok]), crs=4326)
    j = gpd.sjoin(pts, b, how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    out.loc[j.index, "county_fips_pip"] = j["geo_id"].astype("string")
    out.loc[j.index[j["geo_id"].notna()], "pip_method"] = "within_polygon"
    miss = j.index[j["geo_id"].isna()]
    if len(miss):
        mp = pts.loc[miss].to_crs(3857)
        nb = gpd.sjoin_nearest(mp, b.to_crs(3857), how="left", max_distance=PIP_SNAP_METERS, distance_col="d_m")
        nb = nb[~nb.index.duplicated(keep="first")]
        hit = nb.index[nb["geo_id"].notna()]
        out.loc[hit, "county_fips_pip"] = nb.loc[hit, "geo_id"].astype("string")
        out.loc[hit, "pip_method"] = "nearest_polygon_within_5km"
        out.loc[nb.index[nb["geo_id"].isna()], "pip_method"] = "outside_all_2024_counties"
    return out


# --------------------------------------------------------------------------- fetch

def fetch() -> dict[str, Path]:
    paths = {}
    for key, fn in HC_FILES.items():
        paths[key] = download_file(BASE_URL + fn, HC_SOURCE, fn, min_bytes=1000)
    for key, fn in HPSA_FILES.items():
        paths[key] = download_file(BASE_URL + fn, HPSA_SOURCE, fn, min_bytes=1000)
    return paths


def _manifest_entry(source_id: str, filename: str) -> dict:
    return load_manifest(source_id)["files"][filename]


def _read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    return df.loc[:, ~df.columns.str.startswith("Unnamed")]


def dictionary_undefined_columns(dictionary_xlsx: Path, csv_path: Path) -> list[str]:
    """CSV columns not named in the HRSA dictionary (Attribute Name or any User Friendly Name)."""
    m = pd.read_excel(dictionary_xlsx)
    names = set(m.iloc[:, 0].astype(str).str.strip())
    for v in m.iloc[:, 2].dropna().astype(str):
        names |= {x.strip() for x in v.split(";")}
    cols = [c for c in pd.read_csv(csv_path, nrows=1).columns if not c.startswith("Unnamed")]
    return [c for c in cols if c not in names]


def _counties_2024() -> pd.DataFrame:
    g = read_table("geographies")
    return g[(g["geo_level"] == "county") & (~g["ct_legacy"])][
        ["geo_id", "name", "state_fips", "state_abbr", "state_name"]].reset_index(drop=True)


def _pairs(s: pd.Series) -> dict:
    """Two-level MultiIndex Series -> {"a|b": int}."""
    return {f"{a}|{b}": int(v) for (a, b), v in s.items()}


def _py(o):
    """Recursively convert numpy/pandas scalars to plain Python (for YAML/JSON)."""
    if isinstance(o, dict):
        return {str(k): _py(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_py(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, float) and np.isnan(o):
        return None
    if o is pd.NA:
        return None
    return o


def _latest_date(s: pd.Series) -> str:
    """Latest HRSA Data Warehouse record-create date in a column, as ISO yyyy-mm-dd."""
    d = pd.to_datetime(s.astype("string").str.strip().str[:10], errors="coerce", format="mixed")
    return d.max().date().isoformat()


def _pct(n: int, d: int) -> float:
    return round(100.0 * n / d, 2) if d else float("nan")


# --------------------------------------------------------------------------- health-center sites

def build_facilities(paths: dict[str, Path], hpsa_des: pd.DataFrame, auto_links: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    raw = _read_csv(paths["sites"])
    stats: dict = {"rows_raw": int(len(raw)), "n_columns": int(raw.shape[1]),
                   "dictionary_undefined_columns": dictionary_undefined_columns(paths["sites_dictionary"], paths["sites"])}
    c = {k: raw[k] for k in raw.columns}

    df = pd.DataFrame(index=raw.index)
    df["site_bphc_number"] = clean_str(c["BPHC Assigned Number"])
    df["facility_id"] = "hrsa_site:" + df["site_bphc_number"]
    df["facility_name"] = clean_str(c["Site Name"])
    df["facility_source"] = "HRSA Health Center Program site"
    df["health_center_type"] = clean_str(c["Health Center Type"])
    df["is_look_alike"] = df["health_center_type"].str.contains("Look-Alike", regex=False).fillna(False).astype(bool)
    df["health_center_number"] = clean_str(c["Health Center Number"])
    df["bhcmis_org_id"] = clean_str(c["BHCMIS Organization Identification Number"])
    df["health_center_name"] = clean_str(c["Health Center Name"])
    df["health_center_org_address"] = clean_str(c["Health Center Organization Street Address"])
    df["health_center_org_city"] = clean_str(c["Health Center Organization City"])
    df["health_center_org_state"] = clean_str(c["Health Center Organization State"])
    df["health_center_org_zip5"] = zip5(c["Health Center Organization ZIP Code"])
    df["grantee_org_type"] = clean_str(c["Grantee Organization Type Description"])
    df["site_type"] = clean_str(c["Health Center Type Description"])
    df["is_service_delivery_site"] = df["site_type"].isin(
        ["Service Delivery Site", "Administrative/Service Delivery Site"]).astype(bool)
    df["location_setting"] = c["Health Center Service Delivery Site Location Setting Description"].str.strip()
    df["location_type"] = c["Health Center Location Type Description"].str.strip()
    df["operator_type"] = c["Health Center Operator Description"].str.strip()
    df["operating_schedule"] = c["Health Center Operational Schedule Description"].str.strip()
    df["operating_calendar"] = c["Health Center Operating Calendar"].str.strip()
    df["operating_hours_per_week"] = to_num(c["Operating Hours per Week"])
    df["site_status"] = c["Site Status Description"].str.strip()
    df["is_active"] = (df["site_status"] == "Active").astype(bool)
    df["site_added_to_scope_date"] = parse_date(c["Site Added to Scope this Date"], "%m/%d/%Y")
    df["address"] = clean_str(c["Site Address"])
    df["city"] = clean_str(c["Site City"])
    df["state_abbr"] = clean_str(c["Site State Abbreviation"])
    df["state_fips"] = clean_str(c["State FIPS Code"])
    df["zip5"] = zip5(c["Site Postal Code"])
    df["site_phone"] = clean_str(c["Site Telephone Number"])
    df["site_url"] = clean_str(c["Site Web Address"])
    npi_raw = clean_str(c["FQHC Site NPI Number"])
    df["site_npi"] = clean_npi(npi_raw.fillna(""))
    df["site_medicare_billing_number"] = clean_str(c["FQHC Site Medicare Billing Number"])
    df["hhs_region"] = clean_str(c["HHS Region Code"])
    df["congressional_district"] = clean_str(c["Congressional District Code"])
    df["us_mexico_border_100km"] = c["U.S. - Mexico Border 100 Kilometer Indicator"].str.strip().map({"Y": True, "N": False}).astype("boolean")
    df["us_mexico_border_county"] = c["U.S. - Mexico Border County Indicator"].str.strip().map({"Y": True, "N": False}).astype("boolean")
    df["county_name_source"] = clean_str(c["Complete County Name"])
    df["county_fips_source"] = valid_county_fips(c["State and County Federal Information Processing Standard Code"])
    stats["npi_nonempty"] = int(npi_raw.notna().sum())
    stats["npi_valid"] = int(df["site_npi"].notna().sum())

    # --- coordinates: the source's own address geocode (X = lon, Y = lat)
    lon = to_num(c["Geocoding Artifact Address Primary X Coordinate"])
    lat = to_num(c["Geocoding Artifact Address Primary Y Coordinate"])
    bad = (lat.abs() > 90) | (lon.abs() > 180) | ((lat == 0) & (lon == 0))
    stats["source_coords_invalid"] = int(bad.sum())
    lat, lon = lat.mask(bad), lon.mask(bad)
    has_pt = lat.notna() & lon.notna()
    df["lat"], df["lon"] = lat, lon
    df["geocode_method"] = pd.Series("source_address_geocode", index=df.index, dtype="string").where(has_pt)

    # fallback: ZIP -> ZCTA internal point
    zg = zip_to_geo(df["zip5"].fillna(""))
    use_zip = ~has_pt & zg["lat"].notna()
    df.loc[use_zip, "lat"] = zg.loc[use_zip, "lat"]
    df.loc[use_zip, "lon"] = zg.loc[use_zip, "lon"]
    df.loc[use_zip, "geocode_method"] = "zip5_as_zcta_internal_point"
    df["geocode_method"] = df["geocode_method"].fillna("not_geocoded")

    # --- county: source FIPS, checked against point-in-polygon on the 2024 boundaries
    counties = _counties_2024()
    c2024 = set(counties["geo_id"])
    pip = point_in_county(lat.where(has_pt), lon.where(has_pt))
    df["county_fips_pip"] = pip["county_fips_pip"]
    df["pip_method"] = pip["pip_method"]
    src_ok = df["county_fips_source"].isin(c2024)
    df["county_fips_source_in_2024_set"] = src_ok.astype(bool)
    agree = pd.Series(pd.NA, index=df.index, dtype="boolean")
    both = df["county_fips_source"].notna() & df["county_fips_pip"].notna()
    agree[both] = (df.loc[both, "county_fips_source"] == df.loc[both, "county_fips_pip"]).to_numpy()
    df["county_fips_pip_agrees"] = agree
    dis = (agree == False).fillna(False).astype(bool)  # noqa: E712
    df["county_fips"] = df["county_fips_source"].where(src_ok)
    df["county_fips_method"] = pd.Series("source_published", index=df.index, dtype="string").where(src_ok)
    need = df["county_fips"].isna() & df["county_fips_pip"].notna()
    df.loc[need, "county_fips"] = df.loc[need, "county_fips_pip"]
    df.loc[need, "county_fips_method"] = "point_in_polygon_2024"
    need = df["county_fips"].isna() & zg["county_fips"].notna()
    df.loc[need, "county_fips"] = zg.loc[need, "county_fips"]
    df.loc[need, "county_fips_method"] = "zip5_zcta_county"
    df["county_fips_method"] = df["county_fips_method"].fillna("not_assigned")
    df["county_fips"] = df["county_fips"].astype("string")

    # --- auto-HPSA facility designation linked to this site (per discipline)
    for disc in DISCIPLINES:
        lk = auto_links[auto_links["discipline_code"] == disc]
        des = hpsa_des[hpsa_des["discipline_code"] == disc][["hpsa_id", "hpsa_status", "hpsa_score", "designation_type"]]
        lk = lk.merge(des, on="hpsa_id", how="inner")  # in-effect designations only
        lk = lk.sort_values(["site_bphc_number", "hpsa_status", "hpsa_score"], ascending=[True, True, False])
        n_links = lk.groupby("site_bphc_number").size()
        lk1 = lk.drop_duplicates("site_bphc_number").set_index("site_bphc_number")
        df[f"{disc}_facility_hpsa_id"] = df["site_bphc_number"].map(lk1["hpsa_id"]).astype("string")
        df[f"{disc}_facility_hpsa_status"] = df["site_bphc_number"].map(lk1["hpsa_status"]).astype("string")
        df[f"{disc}_facility_hpsa_score"] = df["site_bphc_number"].map(lk1["hpsa_score"]).astype("Int64")
        df[f"{disc}_facility_hpsa_n_linked"] = df["site_bphc_number"].map(n_links).fillna(0).astype("int64")
        stats[f"sites_with_{disc}_facility_hpsa"] = int(df[f"{disc}_facility_hpsa_id"].notna().sum())

    # --- stats
    n = len(df)
    stats.update({
        "sites_total": n,
        "sites_unique_ids": int(df["site_bphc_number"].nunique()),
        "sites_active": int(df["is_active"].sum()),
        "site_status_counts": df["site_status"].value_counts().to_dict(),
        "health_center_type_counts": df["health_center_type"].value_counts().to_dict(),
        "site_type_counts": df["site_type"].value_counts().to_dict(),
        "location_setting_counts": df["location_setting"].value_counts().to_dict(),
        "location_type_counts": df["location_type"].value_counts().to_dict(),
        "operator_type_counts": df["operator_type"].value_counts().to_dict(),
        "service_delivery_sites": int(df["is_service_delivery_site"].sum()),
        "look_alike_sites": int(df["is_look_alike"].sum()),
        "health_center_organizations": int(df["health_center_number"].nunique()),
        "geocode_method_counts": df["geocode_method"].value_counts().to_dict(),
        "county_fips_method_counts": df["county_fips_method"].value_counts().to_dict(),
        "pip_method_counts": df["pip_method"].value_counts().to_dict(),
        "county_fips_source_missing": int(df["county_fips_source"].isna().sum()),
        "county_fips_source_not_in_2024": int((df["county_fips_source"].notna() & ~src_ok).sum()),
        "county_fips_source_not_in_2024_states": df.loc[df["county_fips_source"].notna() & ~src_ok, "state_abbr"]
            .value_counts().to_dict(),
        "pip_compared": int(both.sum()),
        "pip_agree": int(agree.sum()),
        "pip_disagree": int(dis.sum()),
        "pip_disagree_examples": df.loc[dis, ["site_bphc_number", "state_abbr", "county_fips_source",
                                              "county_fips_pip", "pip_method"]].head(12).astype(str).to_dict("records"),
        "pip_disagree_by_method": df.loc[dis, "pip_method"].value_counts().to_dict(),
        "pip_disagree_same_state": int((dis & (df["county_fips_source"].str[:2] == df["county_fips_pip"].str[:2])
                                        .fillna(False)).sum()),
        "pip_disagree_by_state": df.loc[dis, "state_abbr"].value_counts().head(15).to_dict(),
        "sites_with_county": int(df["county_fips"].notna().sum()),
        "counties_with_site": int(df.loc[df["county_fips"].isin(c2024), "county_fips"].nunique()),
        "counties_with_service_delivery_site": int(
            df.loc[df["county_fips"].isin(c2024) & df["is_service_delivery_site"], "county_fips"].nunique()),
        "counties_2024_total": len(c2024),
        "states_territories": int(df["state_abbr"].nunique()),
        "sites_outside_2024_county_set_by_state": df.loc[~df["county_fips"].isin(c2024), "state_abbr"]
            .fillna("no site location (XX)").value_counts().to_dict(),
        "missingness": {col: {"n_missing": int(df[col].isna().sum()), "pct": _pct(int(df[col].isna().sum()), n)}
                        for col in ["facility_name", "address", "city", "zip5", "lat", "county_fips_source",
                                    "county_fips", "operating_hours_per_week", "site_npi", "site_url", "site_phone",
                                    "site_medicare_billing_number"]},
        "dw_record_create_date": _latest_date(c["Data Warehouse Record Create Date"]),
    })

    # per-row geographic resolution
    res = np.select([df["geocode_method"] == "source_address_geocode",
                     df["geocode_method"] == "zip5_as_zcta_internal_point", df["state_abbr"].notna().to_numpy()],
                    ["point", "zcta_centroid", "state"], "none")
    no_addr = df["address"].isna()
    stats["sites_without_address"] = int(no_addr.sum())
    stats["sites_without_address_by_setting"] = df.loc[no_addr, "location_setting"].value_counts().to_dict()
    stats["sites_without_site_state"] = int(df["state_abbr"].isna().sum())
    stats["source_geographic_resolution_counts"] = pd.Series(res).value_counts().to_dict()
    # how far do disagreeing points sit from the county HRSA published?
    if dis.any():
        b = gpd.read_parquet(COUNTY_BOUNDARIES).set_index("geo_id").to_crs(5070)
        sub = df.loc[dis]
        pts = gpd.GeoSeries(gpd.points_from_xy(sub["lon"], sub["lat"]), index=sub.index, crs=4326).to_crs(5070)
        stats["pip_disagree_distance_m"] = [int(round(pts[i].distance(b.loc[sub.at[i, "county_fips_source"], "geometry"])))
                                            for i in sub.index if sub.at[i, "county_fips_source"] in b.index]
    else:
        stats["pip_disagree_distance_m"] = []
    m = _manifest_entry(HC_SOURCE, HC_FILES["sites"])
    df = add_provenance(
        df, data_layer="facility",
        source_name="HRSA Health Center Service Delivery and Look-Alike Sites (data.hrsa.gov)",
        source_version=f"HRSA Data Warehouse extract {stats['dw_record_create_date']} (daily refresh)",
        retrieved_at=m["retrieved_at"], evidence_type="facility_registry", source_record_id="site_bphc_number",
        source_geographic_resolution="point",
        provenance_notes=("Administrative listing of in-scope HRSA health-center sites; no service-line or "
                          "clinical-capability detail. lat/lon = HRSA's address geocode (X/Y) unless "
                          "geocode_method says otherwise. county_fips = HRSA's published FIPS when it is a 2024 "
                          "county, else 2024 point-in-polygon, else ZIP->ZCTA county."),
    )
    df["source_geographic_resolution"] = res
    return df.reset_index(drop=True), stats


# --------------------------------------------------------------------------- HPSA

def _load_hpsa_components(paths: dict[str, Path]) -> tuple[pd.DataFrame, dict]:
    frames, stats = [], {"component_rows_by_discipline_status": {}, "designations_by_discipline_status": {}}
    for disc, label in DISCIPLINES.items():
        raw = _read_csv(paths[f"hpsa_{disc}"])
        comp_type = raw["HPSA Component Type Code"].str.strip().replace({"UNK": "POINT"})
        comp_type = comp_type.where(raw["HPSA Geography Identification Number"].str.strip() != "POINT", "POINT")
        f = pd.DataFrame({
            "discipline_code": disc,
            "discipline": raw["HPSA Discipline Class"].str.strip(),
            "hpsa_id": raw["HPSA ID"].str.strip(),
            "hpsa_name": raw["HPSA Name"].str.strip(),
            "designation_type": raw["Designation Type"].str.strip(),
            "hpsa_type_code": raw["HPSA Type Code"].str.strip(),
            "hpsa_status": raw["HPSA Status"].str.strip(),
            "hpsa_score": to_num(raw["HPSA Score"]).astype("Int64"),
            "mcta_score": (to_num(raw["PC MCTA Score"]).astype("Int64") if "PC MCTA Score" in raw
                           else pd.Series(pd.NA, index=raw.index, dtype="Int64")),
            "designation_date": parse_date(raw["HPSA Designation Date"], "%m/%d/%Y"),
            "last_update_date": parse_date(raw["HPSA Designation Last Update Date"], "%m/%d/%Y"),
            "withdrawn_date": parse_date(raw["Withdrawn Date"], "%m/%d/%Y"),
            "metropolitan_indicator": raw["Metropolitan Indicator"].str.strip(),
            "degree_of_shortage": clean_str(raw["HPSA Degree of Shortage"]).replace({"Not applicable": pd.NA}),
            "hpsa_fte": to_num(raw["HPSA FTE"]),
            "designation_population": to_num(raw["HPSA Designation Population"]),
            "pct_population_below_poverty": to_num(raw["% of Population Below 100% Poverty"]),
            "formal_ratio": clean_str(raw["HPSA Formal Ratio"]),
            "provider_ratio_goal": clean_str(raw["HPSA Provider Ratio Goal"]),
            "population_type": clean_str(raw["HPSA Population Type"]),
            "estimated_served_population": to_num(raw["HPSA Estimated Served Population"]),
            "estimated_underserved_population": to_num(raw["HPSA Estimated Underserved Population"]),
            "hpsa_shortage_fte": to_num(raw["HPSA Shortage"]),
            "primary_state_abbr": clean_str(raw["Primary State Abbreviation"]),
            "component_type": comp_type,
            "component_type_desc": raw["HPSA Component Type Description"].str.strip(),
            "component_geo_id": raw["HPSA Geography Identification Number"].str.strip(),
            "component_name": raw["HPSA Component Name"].str.strip(),
            "component_rural_status": clean_str(raw["Rural Status"]),
            "county_fips_published": raw["State and County Federal Information Processing Standard Code"].str.strip(),
            "facility_address": clean_str(raw["HPSA Address"]),
            "facility_city": clean_str(raw["HPSA City"]),
            "facility_zip5": zip5(raw["HPSA Postal Code"]),
            "facility_lat": to_num(raw["Latitude"]),
            "facility_lon": to_num(raw["Longitude"]),
            "facility_bhcmis_org_id": clean_str(raw["BHCMIS Organization Identification Number"]),
            "dw_record_create_date": raw["Data Warehouse Record Create Date"].str.strip(),
        })
        # For facility rows the published county is 'Common State County FIPS Code'
        f["county_fips_published"] = f["county_fips_published"].where(
            f["component_type"] != "POINT", raw["Common State County FIPS Code"].str.strip())
        f["designation_category"] = hpsa_category(f["designation_type"], f["component_type"])
        stats["component_rows_by_discipline_status"][label] = f["hpsa_status"].value_counts().to_dict()
        stats["designations_by_discipline_status"][label] = (
            f.drop_duplicates(["hpsa_id", "hpsa_status"])["hpsa_status"].value_counts().to_dict())
        stats.setdefault("raw_rows", {})[label] = int(len(f))
        reused = f.groupby("hpsa_id")[["hpsa_status", "designation_date", "hpsa_name"]].nunique().max(axis=1) > 1
        stats.setdefault("raw_ids_reused_by_multiple_designations", {})[label] = {
            "n_ids": int(reused.sum()),
            "statuses": f[f["hpsa_id"].isin(reused[reused].index)]["hpsa_status"].value_counts().to_dict()}
        stats.setdefault("raw_unique_ids", {})[label] = int(f["hpsa_id"].nunique())
        frames.append(f)
    comp = pd.concat(frames, ignore_index=True)
    stats["ids_in_more_than_one_discipline"] = int((comp.groupby("hpsa_id")["discipline_code"].nunique() > 1).sum())
    comp["county_fips"] = component_county(comp["component_geo_id"], comp["county_fips_published"], comp["component_type"])
    return comp, stats


def build_hpsa(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    comp_all, stats = _load_hpsa_components(paths)
    stats["dictionary_undefined_columns"] = dictionary_undefined_columns(paths["hpsa_dictionary"], paths["hpsa_pc"])
    stats["n_columns_pc"] = int(len(pd.read_csv(paths["hpsa_pc"], nrows=1).columns.drop(
        [c for c in pd.read_csv(paths["hpsa_pc"], nrows=1).columns if c.startswith("Unnamed")])))
    counties = _counties_2024()
    c2024 = set(counties["geo_id"])

    comp = comp_all[comp_all["hpsa_status"].isin(IN_EFFECT_STATUSES)].copy()
    stats["withdrawn_component_rows_dropped"] = int((comp_all["hpsa_status"] == "Withdrawn").sum())

    # facility points whose published county is missing or not a 2024 county (e.g. legacy CT
    # counties): 2024 point-in-polygon on the facility's own coordinates
    comp["county_fips_method"] = pd.Series("published_or_geoid_prefix", index=comp.index, dtype="string").where(
        comp["county_fips"].notna())
    renamed = comp["county_fips"].isin(list(FIPS_RENAMES_1TO1))
    comp["county_fips"] = bridge_1to1_renames(comp["county_fips"])
    comp.loc[renamed, "county_fips_method"] = "fips_rename_1to1"
    stats["in_effect_components_fips_rename_1to1"] = int(renamed.sum())
    need = ((comp["component_type"] == "POINT") & ~comp["county_fips"].isin(c2024)
            & comp["facility_lat"].notna() & comp["facility_lon"].notna())
    if need.any():
        pip = point_in_county(comp.loc[need, "facility_lat"], comp.loc[need, "facility_lon"])
        hit = pip.index[pip["county_fips_pip"].notna()]
        comp.loc[hit, "county_fips"] = pip.loc[hit, "county_fips_pip"]
        comp.loc[hit, "county_fips_method"] = "point_in_polygon_2024"
    stats["facility_points_reassigned_by_pip"] = int((comp["county_fips_method"] == "point_in_polygon_2024").sum())
    comp["county_fips_method"] = comp["county_fips_method"].fillna("not_assigned")
    comp["county_in_2024_set"] = comp["county_fips"].isin(c2024)

    # designation-level attributes must be constant within (discipline, hpsa_id)
    key = ["discipline_code", "hpsa_id"]
    des_cols = ["discipline", "hpsa_name", "designation_type", "hpsa_type_code", "designation_category",
                "hpsa_status", "hpsa_score", "mcta_score", "designation_date", "last_update_date",
                "metropolitan_indicator", "degree_of_shortage", "hpsa_fte", "designation_population",
                "pct_population_below_poverty", "formal_ratio", "provider_ratio_goal", "population_type",
                "estimated_served_population", "estimated_underserved_population", "hpsa_shortage_fte",
                "primary_state_abbr"]
    nun = comp.groupby(key)[des_cols].nunique(dropna=True)
    stats["in_effect_ids_with_inconsistent_attributes"] = {k: int(v) for k, v in (nun > 1).sum().items() if v}

    g = comp.groupby(key, sort=True)
    des = g[des_cols].first()
    des["n_components"] = g.size()
    des["component_types"] = g["component_type"].agg(lambda s: ";".join(sorted(set(s))))
    des["county_fips_list"] = g["county_fips"].agg(lambda s: ";".join(sorted(set(s.dropna()))))
    des["n_counties"] = g["county_fips"].nunique()
    whole = comp[comp["component_type"] == "SCTY"].groupby(key)["county_fips"].agg(lambda s: ";".join(sorted(set(s.dropna()))))
    des["whole_county_fips_list"] = whole.reindex(des.index).fillna("")
    des["n_components_county_not_in_2024_set"] = comp[~comp["county_in_2024_set"]].groupby(key).size() \
        .reindex(des.index).fillna(0).astype("int64")
    des["rural_status_values"] = g["component_rural_status"].agg(lambda s: ";".join(sorted(set(s.dropna()))))
    fac = comp[comp["component_type"] == "POINT"].drop_duplicates(key).set_index(key)
    for col in ["facility_address", "facility_city", "facility_zip5", "facility_lat", "facility_lon",
                "facility_bhcmis_org_id"]:
        des[col] = fac[col].reindex(des.index)
    des = des.reset_index()
    des["hpsa_record_id"] = des["discipline_code"] + ":" + des["hpsa_id"]
    des["is_designated"] = (des["hpsa_status"] == "Designated").astype(bool)
    des["is_proposed_for_withdrawal"] = (des["hpsa_status"] == "Proposed For Withdrawal").astype(bool)
    des["component_geography"] = component_geography(des["component_types"])
    is_pt = (des["component_geography"] == "facility_point").to_numpy()
    res = np.select([is_pt & des["facility_lat"].notna().to_numpy(),
                     is_pt & (des["county_fips_list"] != "").to_numpy(),
                     is_pt & des["primary_state_abbr"].notna().to_numpy(), is_pt,
                     (des["component_geography"] == "census_tract").to_numpy()],
                    ["point", "county", "state", "none", "census_tract"], "county")

    dw = _latest_date(comp_all["dw_record_create_date"])
    m = _manifest_entry(HPSA_SOURCE, HPSA_FILES["hpsa_pc"])
    version = f"HRSA Data Warehouse extract {dw} (daily refresh)"
    des = add_provenance(
        des, data_layer="geographic", source_name="HRSA Health Professional Shortage Areas (data.hrsa.gov)",
        source_version=version, retrieved_at=m["retrieved_at"], evidence_type="composite_index",
        source_record_id="hpsa_record_id", source_geographic_resolution="county",
        provenance_notes=("Federal shortage designation (administrative). hpsa_score is HRSA's composite HPSA "
                          "score (PC/MH 0-25, DH 0-26; higher = greater need). Not a disease-burden measure. "
                          "Status Designated or Proposed For Withdrawal; Withdrawn dropped. "
                          "source_geographic_resolution 'county' for county-subdivision/mixed components means "
                          "'located within the listed counties' (no county_subdivision resolution in the vocabulary)."),
    )
    des["source_geographic_resolution"] = res

    comp_out = comp[["discipline_code", "hpsa_id", "hpsa_status", "designation_type", "designation_category",
                     "hpsa_score", "component_type", "component_type_desc", "component_geo_id", "component_name",
                     "component_rural_status", "county_fips_published", "county_fips", "county_fips_method",
                     "county_in_2024_set"]].copy()
    # facility points without coordinates are only known to their county (or state)
    pt = comp_out["component_type"] == "POINT"
    has_xy = comp["facility_lat"].notna() & comp["facility_lon"].notna()
    pt_res = np.select([has_xy.to_numpy(), comp_out["county_fips"].notna().to_numpy(),
                        comp["primary_state_abbr"].notna().to_numpy()], ["point", "county", "state"], "none")
    comp_out["component_record_id"] = (comp_out["discipline_code"] + ":" + comp_out["hpsa_id"] + ":"
                                       + comp_out["component_type"] + ":" + comp_out["component_geo_id"])
    keep = ~comp_out["component_record_id"].duplicated()
    comp_out, pt, pt_res = comp_out[keep].reset_index(drop=True), pt[keep].to_numpy(), pt_res[keep.to_numpy()]
    cres = comp_out["component_type"].map({"CT": "census_tract", "SCTY": "county", "CSD": "county"})
    cres = cres.where(~pt, pd.Series(pt_res, index=comp_out.index))
    comp_out = add_provenance(
        comp_out, data_layer="geographic", source_name="HRSA Health Professional Shortage Areas (data.hrsa.gov)",
        source_version=version, retrieved_at=m["retrieved_at"], evidence_type="composite_index",
        source_record_id="component_record_id", source_geographic_resolution="county",
        provenance_notes=("One row per HPSA component; CSD components recorded at 'county' resolution (containing "
                          "county). Facility points are 'point' only when HRSA gives coordinates, else 'county' / "
                          "'state' / 'none' by what is known."),
    )
    comp_out["source_geographic_resolution"] = cres.fillna("county").to_numpy()

    # --- stats
    stats["in_effect_designations"] = int(len(des))
    stats["in_effect_by_discipline_status"] = des.groupby(["discipline", "hpsa_status"]).size() \
        .unstack(fill_value=0).to_dict("index")
    stats["in_effect_by_discipline_type"] = des.groupby(["discipline", "designation_type"]).size() \
        .unstack(fill_value=0).to_dict("index")
    stats["designated_by_discipline_category"] = des[des["is_designated"]].groupby(
        ["discipline", "designation_category"]).size().unstack(fill_value=0).to_dict("index")
    stats["in_effect_component_rows"] = int(len(comp))
    stats["in_effect_components_county_not_in_2024"] = int((~comp["county_in_2024_set"]).sum())
    stats["in_effect_components_county_not_in_2024_by_status_state"] = _pairs(
        comp.loc[~comp["county_in_2024_set"]].fillna({"primary_state_abbr": "NA"})
        .groupby(["hpsa_status", "primary_state_abbr"]).size())
    stats["area_components_with_unmapped_type"] = int(
        ((comp["component_type"] != "POINT") & ~comp["designation_type"].isin(HPSA_CATEGORY)).sum())
    stats["in_effect_components_county_not_in_2024_by_fips"] = comp.loc[~comp["county_in_2024_set"], "county_fips"] \
        .fillna("NA").value_counts().head(25).to_dict()
    stats["county_fips_method_counts"] = comp["county_fips_method"].value_counts().to_dict()
    stats["hpsa_score_designated_describe"] = {
        lab: des.loc[des["is_designated"] & (des["discipline_code"] == d), "hpsa_score"].astype(float)
        .describe().round(2).to_dict() for d, lab in DISCIPLINES.items()}
    stats["dw_record_create_date"] = dw
    fac_des = des[des["designation_category"] == "facility"]
    stats["facility_designations"] = int(len(fac_des))
    stats["facility_designations_missing_coords"] = int(fac_des["facility_lat"].isna().sum())
    stats["missingness_designations"] = {col: {"n_missing": int(des[col].isna().sum()), "pct": _pct(int(des[col].isna().sum()), len(des))}
                                         for col in ["hpsa_score", "designation_population", "hpsa_fte",
                                                     "pct_population_below_poverty", "formal_ratio",
                                                     "estimated_underserved_population", "facility_lat"]}
    return des, comp_out, comp, stats


def load_auto_links(paths: dict[str, Path]) -> tuple[pd.DataFrame, dict]:
    frames, stats = [], {}
    for disc in DISCIPLINES:
        a = _read_csv(paths[f"auto_{disc}"])
        stats[f"auto_{disc}_rows"] = int(len(a))
        stats[f"auto_{disc}_id_types"] = a["Site Source Identification Number Type"].value_counts().to_dict()
        a = a[a["Site Source Identification Number Type"].str.strip() == "BPHC Assigned Number"]
        frames.append(pd.DataFrame({"discipline_code": disc,
                                    "hpsa_id": a["HPSA Source Identification Number"].str.strip(),
                                    "site_bphc_number": a["Auto-HPSA Facility Site Source Identification Number"].str.strip()}))
    return pd.concat(frames, ignore_index=True).drop_duplicates(), stats


# --------------------------------------------------------------------------- MUA/P

def build_mua(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    raw = _read_csv(paths["mua"])
    ctype = raw["Medically Underserved Area/Population (MUA/P) Component Geographic Type Code"].str.strip()
    f = pd.DataFrame({
        "mua_id": raw["MUA/P ID"].str.strip(),
        "service_area_name": raw["MUA/P Service Area Name"].str.strip(),
        "designation_type_code": raw["Designation Type Code"].str.strip(),
        "designation_type": raw["Designation Type"].str.strip(),
        "status": raw["MUA/P Status Description"].str.strip(),
        "designation_date": parse_date(raw["Designation Date"], "%Y-%m-%d"),
        "update_date": parse_date(raw["MUA/P Update Date"], "%Y-%m-%d"),
        "imu_score": to_num(raw["IMU Score"]),
        "population_type": raw["Population Type"].str.strip(),
        "designation_population": to_num(raw["Designation Population in a Medically Underserved Area/Population (MUA/P)"]),
        "pct_below_poverty": to_num(raw["Percent of Population with Incomes at or Below 100 Percent of the U.S. Federal Poverty Level"]),
        "pct_age_65_plus": to_num(raw["Percentage of Population Age 65 and Over"]),
        "infant_mortality_rate": to_num(raw["Infant Mortality Rate"]),
        "providers_per_1000": to_num(raw["Providers per 1000 Population"]),
        "primary_state_abbr": raw["Primary State Abbreviation"].str.strip(),
        "component_type": ctype,
        "component_geo_id": raw["MUA/P Area Code"].str.strip(),
        "component_rural_status": clean_str(raw["Rural Status Description"]),
        "county_fips_published": raw["State and County Federal Information Processing Standard Code"].str.strip(),
        "dw": raw["Data Warehouse Record Create Date"].str.strip(),
    })
    stats = {"raw_rows": int(len(f)), "raw_unique_ids": int(f["mua_id"].nunique()),
             "dictionary_undefined_columns": dictionary_undefined_columns(paths["mua_dictionary"], paths["mua"]),
             "component_rows_by_status": f["status"].value_counts().to_dict(),
             "designations_by_status": f.drop_duplicates(["mua_id", "status", "designation_date"])["status"]
                 .value_counts().to_dict()}
    # Dictionary: "a Governor designation ... does not receive an IMU score" -> NA for *-GE rows.
    ge = f["designation_type_code"].str.endswith("-GE")
    stats["ge_rows_imu_set_to_na"] = int(ge.sum())
    f["imu_score_published"] = f["imu_score"]  # value as published (GE rows too), for traceability
    f["imu_score"] = f["imu_score"].mask(ge)
    stats["raw_ids_reused_across_statuses"] = int((f.groupby("mua_id")["status"].nunique() > 1).sum())
    f["county_fips"] = component_county(f["component_geo_id"], f["county_fips_published"], f["component_type"])
    f["county_fips_as_published"] = f["county_fips"]
    f["county_fips"] = bridge_1to1_renames(f["county_fips"])      # the two 1:1 renames (FIPS_RENAMES_1TO1)
    c2024 = set(_counties_2024()["geo_id"])
    comp = f[f["status"].isin(IN_EFFECT_STATUSES)].copy()
    comp["county_in_2024_set"] = comp["county_fips"].isin(c2024)
    renamed = comp["county_fips_as_published"].isin(list(FIPS_RENAMES_1TO1))
    stats["in_effect_components_fips_rename_1to1"] = int(renamed.sum())
    key = ["mua_id"]
    cols = ["service_area_name", "designation_type_code", "designation_type", "status", "designation_date",
            "update_date", "imu_score", "imu_score_published", "population_type", "designation_population",
            "pct_below_poverty",
            "pct_age_65_plus", "infant_mortality_rate", "providers_per_1000", "primary_state_abbr"]
    nun = comp.groupby(key)[cols].nunique(dropna=True)
    stats["in_effect_ids_with_inconsistent_attributes"] = {k: int(v) for k, v in (nun > 1).sum().items() if v}
    g = comp.groupby(key)
    des = g[cols].first()
    des["n_components"] = g.size()
    des["component_types"] = g["component_type"].agg(lambda s: ";".join(sorted(set(s))))
    des["county_fips_list"] = g["county_fips"].agg(lambda s: ";".join(sorted(set(s.dropna()))))
    des["n_counties"] = g["county_fips"].nunique()
    des["n_components_county_not_in_2024_set"] = comp[~comp["county_in_2024_set"]].groupby(key).size() \
        .reindex(des.index).fillna(0).astype("int64")
    br = comp[renamed]
    des["county_fips_bridge"] = (br["county_fips_as_published"] + "->" + br["county_fips"]).groupby(
        br["mua_id"]).agg(lambda s: ";".join(sorted(set(s)))).reindex(des.index).fillna("").astype(str)
    des = des.reset_index()
    des["mua_record_id"] = "mua:" + des["mua_id"]
    des["component_geography"] = component_geography(des["component_types"])
    stats["in_effect_non_ge_imu_exactly_zero"] = int((des["imu_score"] == 0).sum())
    stats["in_effect_imu_na"] = int(des["imu_score"].isna().sum())
    stats["in_effect_imu_describe"] = des.loc[des["imu_score"] > 0, "imu_score"].describe().round(2).to_dict()
    stats["missingness_designations"] = {
        col: {"n_missing": int(des[col].isna().sum()), "pct": _pct(int(des[col].isna().sum()), len(des))}
        for col in ["imu_score", "designation_population", "pct_below_poverty", "pct_age_65_plus",
                    "infant_mortality_rate", "providers_per_1000"]}
    stats.update({
        "in_effect_designations": int(len(des)),
        "in_effect_by_status_type": _pairs(des.groupby(["status", "designation_type"]).size()),
        "in_effect_component_rows": int(len(comp)),
        "in_effect_components_county_not_in_2024": int((~comp["county_in_2024_set"]).sum()),
        "in_effect_components_county_not_in_2024_by_fips": comp.loc[~comp["county_in_2024_set"], "county_fips"]
            .fillna("NA").value_counts().head(25).to_dict(),
        "in_effect_components_county_not_in_2024_by_state": comp.loc[~comp["county_in_2024_set"], "primary_state_abbr"]
            .value_counts().to_dict(),
        "dw_record_create_date": _latest_date(f["dw"]),
    })
    m = _manifest_entry(HPSA_SOURCE, HPSA_FILES["mua"])
    des = add_provenance(
        des, data_layer="geographic", source_name="HRSA Medically Underserved Areas/Populations (data.hrsa.gov)",
        source_version=f"HRSA Data Warehouse extract {stats['dw_record_create_date']}", retrieved_at=m["retrieved_at"],
        evidence_type="composite_index", source_record_id="mua_record_id", source_geographic_resolution="county",
        provenance_notes=("Federal MUA/P designation (administrative). IMU = Index of Medical Underservice "
                          "(0-100; <= 62 qualifies; lower = more underserved). Not a disease-burden measure. "
                          "imu_score is NA for Governor's exceptions (dictionary: no IMU); imu_score_published "
                          "keeps the value in the file. source_geographic_resolution 'county' for county-"
                          "subdivision/mixed designations means 'located within the listed counties'."),
    )
    des["source_geographic_resolution"] = np.where(des["component_geography"] == "census_tract",
                                                   "census_tract", "county")
    return des, comp, stats


# --------------------------------------------------------------------------- county context

def mark_undetermined(coverage: pd.Series, n_unplaced: pd.Series) -> pd.Series:
    """'none' becomes UNDETERMINED when a designation on legacy/retired county geography may cover the county."""
    return coverage.mask((coverage == "none") & (n_unplaced > 0), UNDETERMINED)


def build_geo_context(hpsa_comp: pd.DataFrame, mua_comp: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    counties = _counties_2024()
    c2024 = set(counties["geo_id"])
    candidates = legacy_county_candidates(c2024)
    out = counties.rename(columns={"geo_id": "county_fips", "name": "county_name"}).copy()
    for disc in DISCIPLINES:
        sub = hpsa_comp[hpsa_comp["discipline_code"] == disc]
        out = out.merge(summarize_hpsa_by_county(sub, disc), on="county_fips", how="left")
        out = out.merge(unplaced_legacy_counts(sub, "hpsa_id", "hpsa_status", disc, candidates),
                        on="county_fips", how="left")
    out = out.merge(summarize_mua_by_county(mua_comp), on="county_fips", how="left")
    out = out.merge(unplaced_legacy_counts(mua_comp, "mua_id", "status", "mua", candidates),
                    on="county_fips", how="left")
    for col in out.columns:
        if col.endswith(("_whole_county", "_whole_county_geographic", "_whole_county_population",
                         "_partial_county", "_facility_designation")):
            out[col] = out[col].astype("boolean").fillna(False).astype(bool)
        elif re.search(r"_n_(designations|geographic|population|facility|proposed_for_withdrawal"
                       r"|unplaced_legacy_candidates(_pfw)?)$", col):
            out[col] = out[col].fillna(0).astype("int64")
        elif col.endswith("_coverage"):
            out[col] = out[col].fillna("none")
    for p in [*DISCIPLINES, "mua"]:
        out[f"{p}_coverage"] = mark_undetermined(out[f"{p}_coverage"], out[f"{p}_n_unplaced_legacy_candidates"])
    for disc in DISCIPLINES:
        out[f"{disc}_any_designation"] = out[f"{disc}_n_designations"] > 0
        out[f"{disc}_any_incl_proposed_for_withdrawal"] = (out[f"{disc}_n_designations"]
                                                          + out[f"{disc}_n_proposed_for_withdrawal"]) > 0
    out["geo_id"] = out["county_fips"]
    out["geo_level"] = "county"
    stats = {"counties": int(len(out))}
    for disc, lab in DISCIPLINES.items():
        stats[lab] = {
            "counties_any_designated": int(out[f"{disc}_any_designation"].sum()),
            "counties_whole_county": int(out[f"{disc}_whole_county"].sum()),
            "counties_whole_county_geographic": int(out[f"{disc}_whole_county_geographic"].sum()),
            "counties_partial_county": int(out[f"{disc}_partial_county"].sum()),
            "coverage_counts": out[f"{disc}_coverage"].value_counts().reindex(COVERAGE_ORDER, fill_value=0).to_dict(),
            "counties_any_incl_pfw": int(out[f"{disc}_any_incl_proposed_for_withdrawal"].sum()),
            "counties_with_unplaced_legacy_candidates": int((out[f"{disc}_n_unplaced_legacy_candidates"] > 0).sum()),
            "counties_with_unplaced_legacy_candidates_pfw": int(
                (out[f"{disc}_n_unplaced_legacy_candidates_pfw"] > 0).sum()),
        }
    stats["MUA/P"] = {"counties_any_designated": int((out["mua_n_designations"] > 0).sum()),
                      "coverage_counts": out["mua_coverage"].value_counts().to_dict(),
                      "counties_with_unplaced_legacy_candidates": int(
                          (out["mua_n_unplaced_legacy_candidates"] > 0).sum()),
                      "undetermined_counties": sorted(out.loc[out["mua_coverage"] == UNDETERMINED, "county_fips"])}
    stats["legacy_county_candidates"] = {k: list(v) for k, v in sorted(candidates.items())}
    return out, stats


# --------------------------------------------------------------------------- audit + registry

def _fmt_counts(d: dict) -> str:
    return ", ".join(f"{k}: {v:,}" if isinstance(v, (int, np.integer)) else f"{k}: {v}" for k, v in d.items())


def _write_hc_audit(stats: dict, meta: dict) -> None:
    man = load_manifest(HC_SOURCE)["files"]
    files = "\n".join(f"| `{k}` | {v['url']} | {v['bytes']:,} | `{v['sha256'][:16]}…` | {v['retrieved_at']} | {v.get('last_modified')} |"
                      for k, v in sorted(man.items()))
    miss = "\n".join(f"| `{k}` | {v['n_missing']:,} | {v['pct']}% |" for k, v in stats["missingness"].items())
    txt = f"""# DATA AUDIT — HRSA Health Center Service Delivery and Look-Alike Sites

*Generated by `measure_it.ingestion.hrsa` at {meta['generated_at']}; every number below is computed from the retrieved file.*

| Field | Value |
|---|---|
| source_id | `{HC_SOURCE}` |
| Source (dataset/API name, exact files/endpoints) | HRSA Data Downloads → "Health Center Program Service Delivery and Look–Alike Sites": `{HC_FILES['sites']}` + data dictionary `{HC_FILES['sites_dictionary']}` |
| Publishing organization | Health Resources and Services Administration (HRSA), Bureau of Primary Health Care (BPHC); HRSA Data Warehouse |
| Retrieval date (UTC) | {man[HC_FILES['sites']]['retrieved_at']} |
| Source version / release | HRSA Data Warehouse record create date {stats['dw_record_create_date']} (file Last-Modified {man[HC_FILES['sites']].get('last_modified')}) |
| Source update date / cadence | Landing page as read on 2026-09-23: "Updated: 9/23/2026, Refresh Cycle: Daily" |
| License / access conditions | U.S. federal government data; landing page lists "Usage limitations: None"; open download, no registration |
| Unit of observation | One Health Center Program (FQHC) or FQHC Look-Alike site in scope (BPHC assigned site number) |
| Sample size (actual, as ingested) | {stats['sites_total']:,} sites ({stats['sites_unique_ids']:,} unique site numbers) from {stats['health_center_organizations']:,} health-center grants/look-alike designations |
| Geography (resolution, vintage) | Point (HRSA address geocode) + HRSA-published state/county FIPS (2024 vintage incl. CT planning regions); 118th Congressional Districts |
| Person-level? | no |
| Geographic? | yes (facility locations) |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — facility registry; no people |

## Files / endpoints retrieved
Landing page: {LANDING_URL}

| File | URL | Bytes | sha256 | Retrieved | Last-Modified |
|---|---|---|---|---|---|
{files}

## Key variables
| Output column | Source column | Meaning / use |
|---|---|---|
| `facility_id`, `site_bphc_number` | BPHC Assigned Number | Unique site id (`hrsa_site:BPS-H80-…` / `BPS-LAL-…`) |
| `facility_name` | Site Name | Site name |
| `health_center_number`, `health_center_name`, `bhcmis_org_id` | Health Center Number, Health Center Name, BHCMIS ID | Grantee (H80CS…) or look-alike (LALCS…) organisation |
| `is_look_alike` | Health Center Type | FQHC Look-Alike (no Section 330 grant) vs funded FQHC |
| `site_type`, `is_service_delivery_site` | Health Center Type Description | Service Delivery / Administrative / both; administrative-only sites do not see patients |
| `location_setting` | …Location Setting Description | All Other Clinic Types / School / Hospital / … |
| `location_type` | Health Center Location Type Description | Permanent / Seasonal / Mobile Van |
| `operating_hours_per_week` | Operating Hours per Week | hours |
| `site_status`, `is_active` | Site Status Description | status |
| `lat`, `lon`, `geocode_method` | Geocoding Artifact Address Primary Y / X | HRSA's own address geocode (`source_address_geocode`); fallback ZIP→ZCTA internal point |
| `county_fips`, `county_fips_source`, `county_fips_pip`, `county_fips_pip_agrees`, `county_fips_method` | State and County FIPS Code | Published county FIPS, verified by point-in-polygon against 2024 cartographic county boundaries |
| `site_npi` | FQHC Site NPI Number | 10-digit NPI kept only when the NPI check digit validates |
| `pc/mh/dh_facility_hpsa_id/_status/_score` | SDMS_AUTO_HPSA_SITE_*.csv → HPSA files | The organisation's automatic facility HPSA and its HPSA score, linked on BPHC site number |

The file carries **no service-line, specialty, patient-volume or population-served fields**. Grantee-level patient counts exist
in the separate UDS release (not ingested here).

## Counts (measured)
* Sites: {stats['sites_total']:,}; active: {stats['sites_active']:,} (status values in file: {_fmt_counts(stats['site_status_counts'])}). The download contains only currently active in-scope sites, so terminated sites cannot be counted from it.
* Health-center type: {_fmt_counts(stats['health_center_type_counts'])}.
* Site type: {_fmt_counts(stats['site_type_counts'])}; service-delivery sites (incl. admin/service): {stats['service_delivery_sites']:,}.
* Location setting: {_fmt_counts(stats['location_setting_counts'])}.
* Location type: {_fmt_counts(stats['location_type_counts'])}.
* Operator: {_fmt_counts(stats['operator_type_counts'])}.
* States/territories represented: {stats['states_territories']}.
* Geocoding: {_fmt_counts(stats['geocode_method_counts'])}. Invalid source coordinates nulled: {stats['source_coords_invalid']}.
* County assignment: {_fmt_counts(stats['county_fips_method_counts'])}.
* Point-in-polygon check (2024 1:5M boundaries): method {_fmt_counts(stats['pip_method_counts'])}.
  Of {stats['pip_compared']:,} sites with both a published FIPS and a PIP county, {stats['pip_agree']:,} agree and **{stats['pip_disagree']:,} disagree** ({_pct(stats['pip_disagree'], stats['pip_compared'])}%);
  {stats['pip_disagree_same_state']} of the disagreements are within the same state. Distance from the HRSA point to the polygon of the HRSA-published county (EPSG:5070): {stats['pip_disagree_distance_m']} m (per disagreeing site): every disagreeing point lies within {max(stats['pip_disagree_distance_m'] or [0]):,} m of the published county's generalised polygon. Disagreements by PIP method: {_fmt_counts(stats['pip_disagree_by_method'])}; by state: {_fmt_counts(stats['pip_disagree_by_state'])}. The published FIPS is kept as `county_fips`; the PIP value is retained in `county_fips_pip` and `county_fips_pip_agrees` flags each row.
  Examples: {stats['pip_disagree_examples'][:6]}
* Published county FIPS missing: {stats['county_fips_source_missing']}; published but not in the 2024 county set: {stats['county_fips_source_not_in_2024']} (by state: {stats['county_fips_source_not_in_2024_states']}).
* Sites whose final county is not a 2024 Census county (freely associated states FM/MH/PW have no Census counties; unassigned sites): {stats['sites_outside_2024_county_set_by_state']}.
* County coverage vs `data/processed/geographies`: {stats['counties_with_site']:,} of {stats['counties_2024_total']:,} 2024 counties have ≥1 site; {stats['counties_with_service_delivery_site']:,} have ≥1 service-delivery site.
* NPI: {stats['npi_nonempty']:,} non-empty, {stats['npi_valid']:,} pass the check digit and are not the placeholder 1234567890.
* Facility-HPSA link (auto-HPSA site lists): sites linked to an in-effect PC facility HPSA {stats['sites_with_pc_facility_hpsa']:,}; MH {stats['sites_with_mh_facility_hpsa']:,}; DH {stats['sites_with_dh_facility_hpsa']:,}.

## Missingness
| Column | Missing | % |
|---|---|---|
{miss}

{stats['sites_without_address']} sites have no street address ({_fmt_counts(stats['sites_without_address_by_setting'])}); {stats['sites_without_site_state']} of them also carry state code `XX` and no city, ZIP, coordinates or county, so their location is unknown (`source_geographic_resolution='none'`; the parent organisation's address is in `health_center_org_*`). HRSA's landing-page note says address components are suppressed for domestic-violence shelters; the other no-address sites are not explained by the source.
Resolution of `lat`/`lon`: {_fmt_counts(stats['source_geographic_resolution_counts'])}.

## Linkage strategy
* `county_fips` (2024 vintage) → `geographies`, `geo_context__hpsa`, and other county tables. Joins are ecological: a site in a county says nothing about any individual.
* `site_bphc_number` → HRSA auto-HPSA site lists → `hpsa_designations` (`hpsa_id`).
* `health_center_number` / `bhcmis_org_id` → UDS grantee data (not ingested).
* `site_npi` → NPPES (organisation NPI) when present; {stats['npi_valid']:,} sites only. Name/address matching to NPPES is not done here.
* Not joinable to any person-level dataset.

## Limitations and caveats
* Administrative listing of sites in a grant's scope; says nothing about which services, specialties or equipment a site offers.
* Administrative-only sites ({stats['site_type_counts'].get('Administrative', 0):,}) are included and flagged (`is_service_delivery_site=False`).
* Mobile-van sites ({stats['location_type_counts'].get('Mobile Van', 0):,}) carry a registered site address; the dictionary says X/Y is geocoded "based on its address", so the point is that address, not where the van operates.
* HRSA's geocode quality is not documented per row; PIP disagreement is reported above.
* Daily-refreshed: counts change day to day; the extract date is in `source_version`.

## Processed outputs
* `data/processed/facilities__hrsa.parquet` — {stats['sites_total']:,} rows.

## Reproduce
`uv run python -m measure_it.ingestion.hrsa`
"""
    (raw_dir(HC_SOURCE) / "DATA_AUDIT.md").write_text(txt)


RETIRED_COUNTY_CODES = {"02261": "Valdez-Cordova Census Area AK (split 2019)",
                        "02270": "Wade Hampton Census Area AK (renamed Kusilvak 2015, now 02158)",
                        "46113": "Shannon County SD (renamed Oglala Lakota 2015, now 46102)",
                        "51515": "Bedford city VA (merged into Bedford County 2013)"}


def _retired(fips_counts: dict) -> str:
    hits = [f"{k} {v}" for k, v in RETIRED_COUNTY_CODES.items() if k in fips_counts]
    return (" and retired county codes: " + "; ".join(hits)) if hits else ""


def _write_hpsa_audit(hs: dict, ms: dict, gs: dict, meta: dict) -> None:
    man = load_manifest(HPSA_SOURCE)["files"]
    files = "\n".join(f"| `{k}` | {v['url']} | {v['bytes']:,} | `{v['sha256'][:16]}…` | {v['retrieved_at']} | {v.get('last_modified')} |"
                      for k, v in sorted(man.items()))
    rows = []
    for lab in DISCIPLINES.values():
        rows.append(f"| {lab} | {hs['raw_rows'][lab]:,} | {hs['raw_unique_ids'][lab]:,} | "
                    f"{_fmt_counts(hs['designations_by_discipline_status'][lab])} | "
                    f"{_fmt_counts(hs['component_rows_by_discipline_status'][lab])} |")
    status_tbl = "\n".join(rows)
    cov_rows = []
    for lab in DISCIPLINES.values():
        s = gs[lab]
        cc = s["coverage_counts"]
        cov_rows.append(f"| {lab} | {s['counties_any_designated']:,} | {s['counties_whole_county']:,} | "
                        f"{s['counties_whole_county_geographic']:,} | {s['counties_partial_county']:,} | "
                        f"{cc['whole_county_geographic']:,} / {cc['whole_county_population']:,} / {cc['partial_county']:,} / "
                        f"{cc['facility_only']:,} / {cc[UNDETERMINED]:,} / {cc['none']:,} | {s['counties_any_incl_pfw']:,} | "
                        f"{s['counties_with_unplaced_legacy_candidates']:,} / "
                        f"{s['counties_with_unplaced_legacy_candidates_pfw']:,} |")
    cov_tbl = "\n".join(cov_rows)
    miss = "\n".join(f"| `{k}` | {v['n_missing']:,} | {v['pct']}% |" for k, v in hs["missingness_designations"].items())
    type_tbl = "\n".join(f"| {d} | {_fmt_counts({k: v for k, v in t.items() if v})} |"
                         for d, t in hs["in_effect_by_discipline_type"].items())
    score_tbl = "\n".join(f"| {d} | {s.get('count', 0):.0f} | {s.get('mean', float('nan'))} | {s.get('min', float('nan'))} | "
                          f"{s.get('50%', float('nan'))} | {s.get('max', float('nan'))} |" for d, s in hs["hpsa_score_designated_describe"].items())
    txt = f"""# DATA AUDIT — HRSA shortage designations (HPSA primary care / mental health / dental; MUA/P)

*Generated by `measure_it.ingestion.hrsa` at {meta['generated_at']}; every number below is computed from the retrieved files.*

| Field | Value |
|---|---|
| source_id | `{HPSA_SOURCE}` |
| Source (dataset/API name, exact files/endpoints) | HRSA Data Downloads → "Shortage Areas": HPSA "All HPSAs" CSVs for Primary Care, Mental Health, Dental Health; "Auto-HPSA scoring" site lists; MUA/P CSV; data dictionaries |
| Publishing organization | HRSA, Bureau of Health Workforce (BHW), Division of Policy and Shortage Designation (DPSD); HRSA Data Warehouse |
| Retrieval date (UTC) | {man[HPSA_FILES['hpsa_pc']]['retrieved_at']} |
| Source version / release | HRSA Data Warehouse record create date {hs['dw_record_create_date']} (HPSA), {ms['dw_record_create_date']} (MUA/P) |
| Source update date / cadence | Landing page as read on 2026-09-23: "Updated: 9/23/2026, Refresh Cycle: Daily" (Shortage Areas block) |
| License / access conditions | U.S. federal government data; landing page lists "Usage limitations: None"; open download, no registration |
| Unit of observation | Raw: one row per designation component (census tract, county subdivision, whole county, or facility point). Processed: designation (`hpsa_designations`, `mua_designations`), component (`hpsa_components`), 2024 county (`geo_context__hpsa`) |
| Sample size (actual, as ingested) | HPSA raw component rows: {_fmt_counts(hs['raw_rows'])}; in-effect designations kept: {hs['in_effect_designations']:,}; MUA/P raw rows {ms['raw_rows']:,}, in-effect designations {ms['in_effect_designations']:,}; counties summarised: {gs['counties']:,} |
| Geography (resolution, vintage) | Census tract / county subdivision / county components (GEOIDs as published; tract vintage not stated in the file; current Designated CT components use planning-region county codes 091xx) and facility points; summarised to 2024 counties |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — area/facility designations |

## Files / endpoints retrieved
Landing page: {LANDING_URL}

| File | URL | Bytes | sha256 | Retrieved | Last-Modified |
|---|---|---|---|---|---|
{files}

## Status breakdown (measured, before filtering)
Designation counts are unique (HPSA ID, status); component rows are raw CSV rows.

| Discipline | Raw component rows | Unique HPSA IDs | Designations by status | Component rows by status |
|---|---|---|---|---|
{status_tbl}

Kept: status **Designated** and **Proposed For Withdrawal** (not withdrawn); dropped {hs['withdrawn_component_rows_dropped']:,} Withdrawn component rows.
County flags/scores in `geo_context__hpsa` use **Designated only**; Proposed-For-Withdrawal designations are counted separately (`*_n_proposed_for_withdrawal`, `*_any_incl_proposed_for_withdrawal`).
In-effect designations by discipline and status: {hs['in_effect_by_discipline_status']}.
In-effect IDs with inconsistent designation-level attributes across their component rows: {hs['in_effect_ids_with_inconsistent_attributes'] or 'none'}.
Raw HPSA IDs reused by more than one designation (differing status, date or name): {hs['raw_ids_reused_by_multiple_designations']}. HPSA IDs occurring in more than one discipline file: {hs['ids_in_more_than_one_discipline']} (hence the key is discipline + ID).

### In-effect designations by type
| Discipline | Designation type counts |
|---|---|
{type_tbl}

Designated only, by category: {hs['designated_by_discipline_category']}.

### HPSA score (status Designated)
| Discipline | n | mean | min | median | max |
|---|---|---|---|---|---|
{score_tbl}

## County coverage (2024 counties = {gs['counties']:,} rows of `geographies`)
| Discipline | Counties with ≥1 Designated HPSA | Whole-county component (any type) | Whole-county geographic | Partial (tract/CSD) | coverage: whole-geo / whole-pop / partial / facility-only / undetermined-legacy / none | Any incl. Proposed-For-Withdrawal | Counties that may hold an unplaced legacy-coded designation: Designated / PFW |
|---|---|---|---|---|---|---|---|
{cov_tbl}

MUA/P: counties with ≥1 Designated MUA/P {gs['MUA/P']['counties_any_designated']:,}; coverage {gs['MUA/P']['coverage_counts']}; counties that may hold an unplaced legacy-coded Designated MUA/P: {gs['MUA/P']['counties_with_unplaced_legacy_candidates']} (coverage `{UNDETERMINED}`: {', '.join(gs['MUA/P']['undetermined_counties']) or 'none'}).

The two 1:1 county renames (same boundary, new code: {', '.join(f'{k} -> {v}' for k, v in FIPS_RENAMES_1TO1.items())}) are identities, so their components are placed on the 2024 code (since 2026-09-24; `county_fips_published` keeps the published code, HPSA `county_fips_method = fips_rename_1to1`, MUA/P `county_fips_bridge`; in-effect components bridged: HPSA {hs.get('in_effect_components_fips_rename_1to1', 0)}, MUA/P {ms.get('in_effect_components_fips_rename_1to1', 0)}). Every other legacy/retired county code (legacy Connecticut; boundary changes such as 02261 and 51515) is never remapped. To avoid a false "none", each 2024 county whose area overlaps a legacy
code that carries in-effect components gets `*_n_unplaced_legacy_candidates` (Designated) / `*_n_unplaced_legacy_candidates_pfw`
(distinct designations), and a coverage that would otherwise be `none` becomes `{UNDETERMINED}`. Candidate 2024 counties per code
(CT legacy counties: overlap ≥ {LEGACY_CT_MIN_SHARE:.0%} of the legacy county's area between the Census 2020 and 2024 1:5M polygons;
retired codes: Census successor counties): {gs['legacy_county_candidates']}.

Components whose county is not a 2024 county (not forced onto 2024 geography):
* HPSA in-effect component rows: {hs['in_effect_components_county_not_in_2024']:,}; by status|state: {hs['in_effect_components_county_not_in_2024_by_status_state']}; by county code: {hs['in_effect_components_county_not_in_2024_by_fips']}.
  These are Proposed-For-Withdrawal Connecticut components still on legacy CT county/tract codes (09001–09015), facility points with an undetermined county, and freely-associated-state / outlying-area codes.
* MUA/P in-effect component rows: {ms['in_effect_components_county_not_in_2024']:,}; by state: {ms['in_effect_components_county_not_in_2024_by_state']}; by county code: {ms['in_effect_components_county_not_in_2024_by_fips']}.
  Includes Connecticut MUA/Ps on legacy county/tract codes (09001–09015){_retired(ms['in_effect_components_county_not_in_2024_by_fips'])}; 64xxx/68xxx/70xxx are FM / MH / PW codes with no 2024 Census county.
* HPSA county-assignment method: {hs['county_fips_method_counts']}.

MUA/P status (raw designations): {ms['designations_by_status']}; in-effect by status|type: {ms['in_effect_by_status_type']}.
MUA/P in-effect IDs with inconsistent attributes: {ms['in_effect_ids_with_inconsistent_attributes'] or 'none'}; raw MUA/P IDs appearing under more than one status: {ms['raw_ids_reused_across_statuses']} (only in-effect rows are kept). Governor's-exception rows whose IMU was set to NA: {ms['ge_rows_imu_set_to_na']:,}.

## Key variables
| Column | Meaning |
|---|---|
| `hpsa_id`, `discipline_code` | HPSA source ID; key is (`discipline_code`, `hpsa_id`) → `hpsa_record_id` |
| `designation_type`, `designation_category` | Geographic / High Needs Geographic / Population (low-income, Medicaid, migrant farmworker, …) / facility (FQHC, look-alike, RHC, ITU, correctional, state mental hospital, other) |
| `hpsa_status` | Designated / Proposed For Withdrawal |
| `hpsa_score` | HRSA HPSA score; dictionary: integers 0–26, higher = greater priority (observed maxima in the table above); `mcta_score` = PC Maternity Care Target Area supplementary score (0–25) |
| `component_type` | CT tract / CSD county subdivision / SCTY whole county / POINT facility |
| `*_whole_county` | a Designated area or population designation has a whole-county (SCTY) component in the county |
| `*_partial_county` | a Designated designation covers ≥1 tract or county subdivision of the county |
| `*_coverage` | whole_county_geographic > whole_county_population > partial_county > facility_only > `{UNDETERMINED}` (no placed designation, but a legacy-coded one may cover the county) > none |
| `*_n_unplaced_legacy_candidates`, `*_n_unplaced_legacy_candidates_pfw` | distinct Designated / Proposed-For-Withdrawal designations with a component on a legacy CT or retired county code whose area overlaps this county (not counted in `*_n_designations`) |
| `imu_score_published` (`mua_designations`) | IMU as published, including Governor's-exception rows whose `imu_score` is set to NA |
| `*_max_score`, `*_mean_score` | max / mean HPSA score over distinct Designated designations touching the county (also split by geographic / population / facility) |
| `mua_min_imu_score` | lowest Index of Medical Underservice among Designated MUA/Ps in the county; dictionary: 0 = highest need, 100 = lowest, must be ≤ 62.0 to qualify, Governor's exceptions get no IMU (set to NA). Exactly-zero IMUs are excluded from the minimum ({ms['in_effect_non_ge_imu_exactly_zero']} in-effect non-exception MUA/Ps have IMU = 0.00, which may encode "not scored"). In-effect IMU > 0: {ms['in_effect_imu_describe']} |

{len(hs['dictionary_undefined_columns'])} of the {hs['n_columns_pc']} PC CSV columns are not named in `HPSA_DATAMART_METADATA.XLSX` (neither attribute nor user-friendly name): {', '.join(hs['dictionary_undefined_columns'])}. Their meaning is taken from the column names and value sets and was checked for internal consistency (above); the component-type and status value sets match the MUA/P dictionary's definitions of the same fields. The MUA/P CSV has {len(ms['dictionary_undefined_columns'])} undefined columns.

## Missingness (in-effect `hpsa_designations`)
| Column | Missing | % |
|---|---|---|
{miss}

MUA/P (in-effect `mua_designations`):

| Column | Missing | % |
|---|---|---|
{chr(10).join(f"| `{k}` | {v['n_missing']:,} | {v['pct']}% |" for k, v in ms['missingness_designations'].items())}

## Linkage strategy
* `geo_context__hpsa.county_fips` (2024) → `geographies` and other county tables (ecological join).
* `hpsa_components.component_geo_id` gives tract/CSD GEOIDs for sub-county work.
* `hpsa_designations.hpsa_id` ← auto-HPSA site list ← `facilities__hrsa.site_bphc_number` (facility HPSA of the health-center organisation).
* Nothing here is person-level; no join to participants.

## Limitations and caveats
* A shortage designation is an administrative federal determination (provider-to-population ratio plus need indicators); it is **not** a disease-prevalence or burden estimate and must not be used as one.
* A "partial" county flag means some tracts/subdivisions or a sub-population are designated; it is not a county-wide statement. Tract-level coverage fraction is not computed.
* Population HPSAs designate a sub-population (e.g. low-income) of the area, not every resident.
* Facility HPSAs (FQHC, look-alike, RHC, ITU, correctional, state mental hospital, other) are points: {hs['facility_designations']:,} in-effect facility designations, {hs['facility_designations_missing_coords']} without coordinates (so `facility_lat` missingness above is almost entirely area designations, which have no point). The auto-HPSA facility designation of a health center is named after the organisation and shared by its sites.
* Proposed-For-Withdrawal status is volatile; daily refresh means counts drift.
* Legacy-CT codes and retired codes with a boundary change are reported, not remapped; the two 1:1 renames (02270 -> 02158, 46113 -> 46102) are placed on the 2024 code.

## Processed outputs
* `data/processed/hpsa_designations.parquet` — {meta['n_hpsa_designations']:,} rows
* `data/processed/hpsa_components.parquet` — {meta['n_hpsa_components']:,} rows
* `data/processed/mua_designations.parquet` — {meta['n_mua_designations']:,} rows
* `data/processed/geo_context__hpsa.parquet` — {meta['n_geo_context']:,} rows

## Reproduce
`uv run python -m measure_it.ingestion.hrsa`
"""
    (raw_dir(HPSA_SOURCE) / "DATA_AUDIT.md").write_text(txt)


def _write_registry(hcs: dict, hs: dict, ms: dict, gs: dict) -> None:
    hc_man = load_manifest(HC_SOURCE)["files"]
    hp_man = load_manifest(HPSA_SOURCE)["files"]
    write_registry_entry({
        "source_id": HC_SOURCE,
        "name": "HRSA Health Center Service Delivery and Look-Alike Sites",
        "publisher": "HRSA Bureau of Primary Health Care (HRSA Data Warehouse)",
        "landing_url": LANDING_URL,
        "access_urls": [BASE_URL + f for f in HC_FILES.values()],
        "license": "U.S. federal government data (public domain); landing page: 'Usage limitations: None'",
        "access_conditions": "open download, no registration",
        "retrieved_at": hc_man[HC_FILES["sites"]]["retrieved_at"],
        "source_version": f"HRSA Data Warehouse extract {hcs['dw_record_create_date']}",
        "update_date": "daily refresh (landing page 'Updated: 9/23/2026')",
        "data_layer": "facility",
        "unit_of_observation": "health-center service-delivery / look-alike site",
        "sample_size": {"sites": hcs["sites_total"], "active_sites": hcs["sites_active"],
                        "service_delivery_sites": hcs["service_delivery_sites"],
                        "look_alike_sites": hcs["look_alike_sites"],
                        "health_center_organizations": hcs["health_center_organizations"],
                        "geocoded_source_point": hcs["geocode_method_counts"].get("source_address_geocode", 0),
                        "counties_with_site": hcs["counties_with_site"]},
        "geographic_resolution": "point (HRSA address geocode); county FIPS (2024)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (facility registry)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": ["facilities__hrsa"],
        "audit": f"data/raw/{HC_SOURCE}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.hrsa",
        "limitations": [
            "Administrative listing of in-scope sites; no service lines, specialties or patient volumes",
            "Download contains only active sites",
            f"{hcs['sites_without_address']} sites lack a street address; {hcs['sites_without_site_state']} have no "
            "site location at all (state 'XX'; incl. domestic-violence shelters, whose address HRSA suppresses)",
            f"{hcs['pip_disagree']} sites where HRSA's county FIPS disagrees with 2024 point-in-polygon (reported, source kept)",
        ],
    })
    write_registry_entry({
        "source_id": HPSA_SOURCE,
        "name": "HRSA shortage designations: HPSA (primary care, mental health, dental) and MUA/P",
        "publisher": "HRSA Bureau of Health Workforce, Division of Policy and Shortage Designation (HRSA Data Warehouse)",
        "landing_url": LANDING_URL,
        "access_urls": [BASE_URL + f for f in HPSA_FILES.values()],
        "license": "U.S. federal government data (public domain); landing page: 'Usage limitations: None'",
        "access_conditions": "open download, no registration",
        "retrieved_at": hp_man[HPSA_FILES["hpsa_pc"]]["retrieved_at"],
        "source_version": f"HRSA Data Warehouse extract {hs['dw_record_create_date']}",
        "update_date": "daily refresh (landing page 'Updated: 9/23/2026')",
        "data_layer": "geographic",
        "unit_of_observation": "shortage designation (components: tract / county subdivision / county / facility); summarised to 2024 county",
        "sample_size": {"hpsa_raw_component_rows": hs["raw_rows"],
                        "hpsa_in_effect_designations": hs["in_effect_designations"],
                        "hpsa_in_effect_by_discipline_status": hs["in_effect_by_discipline_status"],
                        "mua_in_effect_designations": ms["in_effect_designations"],
                        "counties": gs["counties"],
                        "counties_with_designated_pc_hpsa": gs["Primary Care"]["counties_any_designated"],
                        "counties_with_designated_mh_hpsa": gs["Mental Health"]["counties_any_designated"]},
        "geographic_resolution": "census tract / county subdivision / county / facility point; county summary (2024)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (area/facility designations)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": ["hpsa_designations", "hpsa_components", "mua_designations", "geo_context__hpsa"],
        "audit": f"data/raw/{HPSA_SOURCE}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.hrsa",
        "limitations": [
            "Administrative shortage designation, not disease burden",
            "Partial-county flag = some tracts/subdivisions or a sub-population designated",
            "Withdrawn designations dropped; Proposed For Withdrawal kept but excluded from county flags",
            "Legacy Connecticut and retired county codes with a boundary change reported as unmatched, not "
            "remapped (the 1:1 renames 02270 -> 02158 and 46113 -> 46102 are placed on the 2024 code); overlapping 2024 "
            f"counties that would otherwise read 'none' are marked '{UNDETERMINED}' in geo_context__hpsa",
            "evidence_type recorded as composite_index pending an 'administrative_designation' vocabulary entry",
        ],
    })


# --------------------------------------------------------------------------- entry point

def run() -> dict:
    paths = fetch()
    auto_links, auto_stats = load_auto_links(paths)
    hpsa_des, hpsa_comp_out, hpsa_comp, hs = build_hpsa(paths)
    hs.update(auto_stats)
    mua_des, mua_comp, ms = build_mua(paths)
    fac, hcs = build_facilities(paths, hpsa_des, auto_links)
    hs, ms, hcs = _py(hs), _py(ms), _py(hcs)

    geo, gs = build_geo_context(hpsa_comp[["discipline_code", "hpsa_id", "hpsa_status", "designation_type",
                                           "designation_category", "hpsa_score", "component_type", "county_fips"]],
                                mua_comp[["mua_id", "status", "designation_type", "imu_score", "component_type",
                                          "county_fips"]])
    gs = _py(gs)
    m = _manifest_entry(HPSA_SOURCE, HPSA_FILES["hpsa_pc"])
    geo = add_provenance(
        geo, data_layer="geographic",
        source_name="HRSA HPSA (PC/MH/DH) + MUA/P designations summarised to 2024 counties",
        source_version=f"HRSA Data Warehouse extract {hs['dw_record_create_date']}", retrieved_at=m["retrieved_at"],
        evidence_type="composite_index", source_record_id="county_fips", source_geographic_resolution="county",
        provenance_notes=("Presence/score of federal shortage designations touching the county (status Designated; "
                          "Proposed For Withdrawal counted separately). Partial = some tracts/subdivisions or a "
                          "sub-population. Administrative designation, not disease burden. IDs listed in *_hpsa_ids / "
                          "mua_ids for tracing. Designations published on legacy Connecticut or retired county "
                          "codes with a boundary change are not remapped (the 1:1 renames 02270 -> 02158 and "
                          "46113 -> 46102 are placed on the 2024 code): counts/flags exclude them, "
                          "*_n_unplaced_legacy_candidates "
                          "counts those whose legacy area overlaps this county, and a coverage of 'none' is "
                          f"replaced by '{UNDETERMINED}' there."),
    )

    write_table(fac, "facilities__hrsa", producer=PRODUCER,
                description="HRSA health-center service-delivery and look-alike sites (one row per site)")
    write_table(hpsa_des, "hpsa_designations", producer=PRODUCER,
                description="In-effect HRSA HPSA designations (PC/MH/DH), one row per (discipline, HPSA ID)")
    write_table(hpsa_comp_out, "hpsa_components", producer=PRODUCER,
                description="Components (tract/CSD/county/facility) of in-effect HPSA designations with county FIPS")
    write_table(mua_des, "mua_designations", producer=PRODUCER,
                description="In-effect HRSA MUA/P designations, one row per MUA/P ID")
    write_table(geo, "geo_context__hpsa", producer=PRODUCER,
                description="2024-county summary of HPSA (PC/MH/DH) and MUA/P designations")

    meta = {"generated_at": utc_now_iso(), "n_hpsa_designations": len(hpsa_des),
            "n_hpsa_components": len(hpsa_comp_out), "n_mua_designations": len(mua_des), "n_geo_context": len(geo)}
    _write_hc_audit(hcs, meta)
    _write_hpsa_audit(hs, ms, gs, meta)
    _write_registry(hcs, hs, ms, gs)
    stats = {"health_centers": hcs, "hpsa": hs, "mua": ms, "geo_context": gs, "meta": meta}
    (raw_dir(HC_SOURCE) / "audit_stats.json").write_text(json.dumps(hcs, indent=2, default=str))
    (raw_dir(HPSA_SOURCE) / "audit_stats.json").write_text(json.dumps({"hpsa": hs, "mua": ms, "geo_context": gs},
                                                                      indent=2, default=str))
    return stats


if __name__ == "__main__":
    s = run()
    print(json.dumps({"sites": s["health_centers"]["sites_total"],
                      "hpsa_in_effect": s["hpsa"]["in_effect_designations"],
                      "mua_in_effect": s["mua"]["in_effect_designations"],
                      "counties": s["geo_context"]["counties"]}, indent=2))
