"""Canonical U.S. geographies and the ZIP/ZCTA -> county crosswalk.

Canonical county set = Census 2024 vintage (Connecticut uses the 9 planning
regions, FIPS 09110-09190). The 8 legacy Connecticut counties (09001-09015)
are kept as a flagged alternate set because several sources (CMS, older CDC
releases) still publish on them. Nothing here reassigns legacy-CT values to
planning regions; mismatches must be reported, not forced.

ZIP geocoding: a USPS ZIP is treated as the ZCTA with the same code and placed
at the ZCTA's Census internal point. PO-box and unique ZIPs have no ZCTA and
stay ungeocoded. This is an approximation of a facility's location
(geocode_method = 'zip5_as_zcta_internal_point'); it is not an address geocode.

Outputs (data/processed):
  geographies            one row per state and county (2024) + legacy CT counties
  zcta_centroids         ZCTA internal points with the 2024 county containing them
  zcta_county_crosswalk  Census 2020 ZCTA-county relationship with land-area shares
"""
from __future__ import annotations

import io
import zipfile

import geopandas as gpd
import pandas as pd

from ..config import PROCESSED, raw_dir
from ..download import download_file, load_manifest
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, table_exists, write_table

SOURCE_ID = "census_geography"
URLS = {
    "gaz_counties": "https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2024_Gazetteer/2024_Gaz_counties_national.zip",
    "gaz_zcta": "https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2024_Gazetteer/2024_Gaz_zcta_national.zip",
    "rel_zcta_county": "https://www2.census.gov/geo/docs/maps-data/data/rel2020/zcta520/tab20_zcta520_county20_natl.txt",
    "cb_county_2024": "https://www2.census.gov/geo/tiger/GENZ2024/shp/cb_2024_us_county_5m.zip",
    "cb_state_2024": "https://www2.census.gov/geo/tiger/GENZ2024/shp/cb_2024_us_state_5m.zip",
    "cb_county_2020": "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_5m.zip",
}

STATE_FIPS = {
    "01": "AL", "02": "AK", "04": "AZ", "05": "AR", "06": "CA", "08": "CO", "09": "CT", "10": "DE", "11": "DC",
    "12": "FL", "13": "GA", "15": "HI", "16": "ID", "17": "IL", "18": "IN", "19": "IA", "20": "KS", "21": "KY",
    "22": "LA", "23": "ME", "24": "MD", "25": "MA", "26": "MI", "27": "MN", "28": "MS", "29": "MO", "30": "MT",
    "31": "NE", "32": "NV", "33": "NH", "34": "NJ", "35": "NM", "36": "NY", "37": "NC", "38": "ND", "39": "OH",
    "40": "OK", "41": "OR", "42": "PA", "44": "RI", "45": "SC", "46": "SD", "47": "TN", "48": "TX", "49": "UT",
    "50": "VT", "51": "VA", "53": "WA", "54": "WV", "55": "WI", "56": "WY", "60": "AS", "66": "GU", "69": "MP",
    "72": "PR", "78": "VI",
}
STATE_ABBR_TO_FIPS = {v: k for k, v in STATE_FIPS.items()}


def _fetch_all() -> dict:
    return {k: download_file(u, SOURCE_ID) for k, u in URLS.items()}


def _read_gazetteer(path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as zf:
        name = [n for n in zf.namelist() if n.endswith(".txt")][0]
        df = pd.read_csv(io.BytesIO(zf.read(name)), sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    return df


def _retrieved(filename: str) -> str:
    return load_manifest(SOURCE_ID)["files"][filename]["retrieved_at"]


def build() -> None:
    files = _fetch_all()

    # --- counties (2024) and states
    counties = gpd.read_file(f"zip://{files['cb_county_2024']}")
    counties = counties.rename(columns={"GEOID": "geo_id", "NAME": "name", "STATEFP": "state_fips",
                                        "STUSPS": "state_abbr", "STATE_NAME": "state_name", "ALAND": "aland_m2"})
    rp = counties.to_crs(5070).representative_point().to_crs(4326)
    counties["lat"], counties["lon"] = rp.y, rp.x
    gaz = _read_gazetteer(files["gaz_counties"])[["GEOID", "INTPTLAT", "INTPTLONG"]]
    counties = counties.merge(gaz.rename(columns={"GEOID": "geo_id"}), on="geo_id", how="left")
    counties["lat"] = pd.to_numeric(counties["INTPTLAT"]).fillna(counties["lat"])
    counties["lon"] = pd.to_numeric(counties["INTPTLONG"]).fillna(counties["lon"])
    cty = pd.DataFrame({
        "geo_id": counties["geo_id"], "geo_level": "county", "name": counties["NAMELSAD"],
        "state_fips": counties["state_fips"], "state_abbr": counties["state_abbr"],
        "state_name": counties["state_name"], "county_fips": counties["geo_id"],
        "lat": counties["lat"], "lon": counties["lon"], "aland_m2": counties["aland_m2"].astype("int64"),
        "vintage": "2024", "ct_legacy": False,
    })

    states = gpd.read_file(f"zip://{files['cb_state_2024']}")
    srp = states.to_crs(5070).representative_point().to_crs(4326)
    st = pd.DataFrame({
        "geo_id": states["GEOID"], "geo_level": "state", "name": states["NAME"],
        "state_fips": states["STATEFP"], "state_abbr": states["STUSPS"], "state_name": states["NAME"],
        "county_fips": None, "lat": srp.y.values, "lon": srp.x.values,
        "aland_m2": states["ALAND"].astype("int64"), "vintage": "2024", "ct_legacy": False,
    })

    legacy = gpd.read_file(f"zip://{files['cb_county_2020']}")
    legacy = legacy[legacy["STATEFP"] == "09"]
    lrp = legacy.to_crs(5070).representative_point().to_crs(4326)
    ct = pd.DataFrame({
        "geo_id": legacy["GEOID"], "geo_level": "county", "name": legacy["NAMELSAD"],
        "state_fips": "09", "state_abbr": "CT", "state_name": "Connecticut", "county_fips": legacy["GEOID"],
        "lat": lrp.y.values, "lon": lrp.x.values, "aland_m2": legacy["ALAND"].astype("int64"),
        "vintage": "2020_ct_legacy", "ct_legacy": True,
    })
    geos = pd.concat([st, cty, ct], ignore_index=True)
    geos = add_provenance(
        geos, data_layer="geographic", source_name="U.S. Census Bureau cartographic boundary files + Gazetteer",
        source_version="cb_2024_5m; 2024 Gazetteer; cb_2020_5m (CT legacy)", retrieved_at=_retrieved("cb_2024_us_county_5m.zip"),
        evidence_type="census_estimate", source_record_id="geo_id",
        source_geographic_resolution="county",
        provenance_notes="Internal point from Gazetteer where available; CT legacy counties flagged ct_legacy=True",
    )
    geos.loc[geos["geo_level"] == "state", "source_geographic_resolution"] = "state"
    write_table(geos, "geographies", producer="geography.crosswalk.build",
                description="States + 2024 counties + legacy CT counties with internal points")

    # --- county boundaries for maps (GeoParquet)
    gcty = counties[["geo_id", "NAMELSAD", "state_abbr", "geometry"]].rename(columns={"NAMELSAD": "name"})
    gcty.to_parquet(PROCESSED / "county_boundaries_2024.geoparquet")
    gst = states[["GEOID", "NAME", "STUSPS", "geometry"]].rename(columns={"GEOID": "geo_id", "NAME": "name", "STUSPS": "state_abbr"})
    gst.to_parquet(PROCESSED / "state_boundaries_2024.geoparquet")

    # --- ZCTA internal points -> containing 2024 county
    z = _read_gazetteer(files["gaz_zcta"])
    z = z.rename(columns={"GEOID": "zcta"})
    z["lat"] = pd.to_numeric(z["INTPTLAT"])
    z["lon"] = pd.to_numeric(z["INTPTLONG"])
    zg = gpd.GeoDataFrame(z[["zcta", "lat", "lon", "ALAND"]], geometry=gpd.points_from_xy(z["lon"], z["lat"]), crs=4326)
    joined = gpd.sjoin(zg, counties[["geo_id", "state_fips", "state_abbr", "geometry"]].to_crs(4326),
                       how="left", predicate="within")
    joined = joined.drop_duplicates("zcta")
    joined["county_assignment"] = joined["geo_id"].notna().map({True: "internal_point_in_polygon", False: None})
    # ~30 coastal ZCTAs have internal points in water outside the generalised 5m polygons
    # (e.g. 60611 Chicago). Assign the nearest county polygon and flag it.
    miss = joined["geo_id"].isna()
    if miss.any():
        pts = joined.loc[miss, ["zcta", "geometry"]].to_crs(5070)
        near = gpd.sjoin_nearest(pts, counties[["geo_id", "state_fips", "state_abbr", "geometry"]].to_crs(5070),
                                 how="left", distance_col="nearest_county_m").drop_duplicates("zcta").set_index("zcta")
        idx = joined.loc[miss, "zcta"]
        for col in ("geo_id", "state_fips", "state_abbr"):
            joined.loc[miss, col] = idx.map(near[col]).values
        joined.loc[miss, "county_assignment"] = "nearest_polygon (" + idx.map(near["nearest_county_m"]).round(0).astype(int).astype(str).values + " m)"
    zc = pd.DataFrame({
        "zcta": joined["zcta"], "lat": joined["lat"], "lon": joined["lon"],
        "aland_m2": pd.to_numeric(joined["ALAND"]).astype("int64"),
        "county_fips": joined["geo_id"], "state_fips": joined["state_fips"], "state_abbr": joined["state_abbr"],
        "county_assignment": joined["county_assignment"],
    })
    zc = add_provenance(
        zc, data_layer="geographic", source_name="U.S. Census Bureau 2024 ZCTA Gazetteer + cb_2024 county boundaries",
        source_version="2024", retrieved_at=_retrieved("2024_Gaz_zcta_national.zip"),
        evidence_type="census_estimate", source_record_id="zcta", source_geographic_resolution="zcta_centroid",
        provenance_notes="county = 2024 county polygon containing the ZCTA internal point; see county_assignment for nearest-polygon fallbacks",
    )
    write_table(zc, "zcta_centroids", producer="geography.crosswalk.build",
                description="ZCTA internal points and containing 2024 county")

    # --- 2020 ZCTA-county relationship with land-area shares (for apportionment)
    rel = pd.read_csv(files["rel_zcta_county"], sep="|", dtype=str)
    rel = rel.dropna(subset=["GEOID_ZCTA5_20", "GEOID_COUNTY_20"])
    rel["arealand_part"] = pd.to_numeric(rel["AREALAND_PART"])
    rel["arealand_zcta"] = pd.to_numeric(rel["AREALAND_ZCTA5_20"])
    xw = pd.DataFrame({
        "zcta": rel["GEOID_ZCTA5_20"], "county_fips_2020": rel["GEOID_COUNTY_20"],
        "land_share_of_zcta": (rel["arealand_part"] / rel["arealand_zcta"]).where(rel["arealand_zcta"] > 0),
    })
    xw = add_provenance(
        xw, data_layer="geographic", source_name="U.S. Census Bureau 2020 ZCTA5-county relationship file",
        source_version="rel2020", retrieved_at=_retrieved("tab20_zcta520_county20_natl.txt"),
        evidence_type="census_estimate", source_record_id=lambda d: d["zcta"] + ":" + d["county_fips_2020"],
        source_geographic_resolution="zcta",
        provenance_notes="2020 county vintage (legacy CT counties); land-area share, not population share",
    )
    write_table(xw, "zcta_county_crosswalk", producer="geography.crosswalk.build",
                description="2020 ZCTA-county relationship with land-area shares")

    write_registry_entry({
        "source_id": SOURCE_ID, "name": "Census cartographic boundaries, Gazetteer and ZCTA-county relationship files",
        "publisher": "U.S. Census Bureau", "landing_url": "https://www.census.gov/geographies/reference-files.html",
        "access_urls": list(URLS.values()), "license": "U.S. federal government work (public domain)",
        "access_conditions": "open download, no registration", "retrieved_at": _retrieved("cb_2024_us_county_5m.zip"),
        "source_version": "2024 cartographic boundary + Gazetteer; 2020 relationship file", "update_date": "annual",
        "data_layer": "geographic", "unit_of_observation": "state / county / ZCTA",
        "sample_size": {"counties_2024": int((cty["geo_level"] == "county").sum()), "states": int(len(st)),
                        "zctas": int(len(zc)), "zctas_without_county": int(zc["county_fips"].isna().sum())},
        "geographic_resolution": "state, county, zcta", "person_level": False, "geographic": True, "omics": False,
        "wearable": False, "participant_linkage": "not applicable",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": ["geographies", "zcta_centroids", "zcta_county_crosswalk",
                              "county_boundaries_2024.geoparquet", "state_boundaries_2024.geoparquet"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": "measure_it.geography.crosswalk",
        "limitations": ["ZIP != ZCTA; PO-box ZIPs are not geocoded",
                        "Connecticut county vintages differ (planning regions vs legacy counties)"],
    })


def zip_to_geo(zips: pd.Series) -> pd.DataFrame:
    """Map U.S. ZIPs (5 digits or ZIP+4, with or without hyphen) to ZCTA internal point + 2024 county.

    Anything else (foreign postal codes such as 6-digit codes, 4-digit ZIPs that lost
    their leading zero) is left ungeocoded rather than guessed. Callers should still
    restrict to U.S. addresses.
    """
    zc = read_table("zcta_centroids", columns=["zcta", "lat", "lon", "county_fips", "state_fips"])
    z5 = zips.astype(str).str.strip().str.extract(r"^(\d{5})(?:-?\d{4})?$")[0]
    out = pd.DataFrame({"zip5": z5}).merge(zc, left_on="zip5", right_on="zcta", how="left")
    out["geocode_method"] = out["lat"].notna().map({True: "zip5_as_zcta_internal_point", False: "not_geocoded"})
    out.index = zips.index
    return out[["zip5", "lat", "lon", "county_fips", "state_fips", "geocode_method"]]


def county_lookup() -> pd.DataFrame:
    g = read_table("geographies")
    return g[g["geo_level"] == "county"]


if __name__ == "__main__":
    build()
