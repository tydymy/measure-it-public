"""CDC PLACES (latest release on data.cdc.gov): county long-format context + ZCTA drill-down.

PLACES values are MODELED small-area estimates (multilevel regression and
post-stratification of BRFSS responses onto Census/ACS populations). They are
population context and, at most, labelled proxies. They are never observed
county prevalence of any invisible illness, and PLACES has no long-COVID,
ME/CFS, POTS or Lyme measure (the module checks for long COVID on every run).

Release discovery: the Socrata discovery API
(https://api.us.socrata.com/api/catalog/v1?domains=data.cdc.gov&q=PLACES) is queried and the
newest "PLACES: Local Data for Better Health, {County|ZCTA} Data, <year> release" pair is used.
Each dataset's /api/views/<id>.json metadata is saved next to the data (release year, BRFSS
data years, row-update time).

Raw files (data/raw/cdc_places, recorded in MANIFEST.json):
  places_county_<release>_<id>.csv                    full county file, all measures, crude + age-adjusted
  places_zcta_<release>_<id>_selected_measures.csv    SoQL extract, flagged measures only (crude; PLACES
                                                      publishes crude prevalence only for ZCTAs)
  places_data_dictionary_m35w-spkz.csv, metadata_<id>.json

Outputs (data/processed):
  geo_context__cdc_places        county x measure x value type (+ the US national row), long format
  geo_context_zcta__cdc_places   ZCTA x flagged measure, crude prevalence

Run: uv run python -m measure_it.ingestion.cdc_places
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd

from ..config import raw_dir
from ..download import download_file, load_manifest
from ..http import get_json
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

SOURCE_ID = "cdc_places"
SOURCE_NAME = "CDC PLACES: Local Data for Better Health"
PRODUCER = "measure_it.ingestion.cdc_places"
CATALOG_URL = "https://api.us.socrata.com/api/catalog/v1"
DOMAIN = "data.cdc.gov"
DICTIONARY_ID = "m35w-spkz"
LANDING_URL = "https://www.cdc.gov/places/"

# Measures kept flagged for the invisible-illness use case (SPEC section E), with a context group.
# Whole categories are flagged too, so newly added disability / health-status / social-needs
# measures are picked up without editing this list.
RELEVANT_MEASURES: dict[str, str] = {
    "GHLTH": "health_status", "PHLTH": "health_status", "MHLTH": "health_status",
    "DEPRESSION": "mental_health",
    "LPA": "risk_behavior", "SLEEP": "risk_behavior", "CSMOKING": "risk_behavior",
    "OBESITY": "chronic_disease",
    "ACCESS2": "access_to_care", "CHECKUP": "access_to_care",
    "ARTHRITIS": "chronic_disease", "CASTHMA": "chronic_disease", "COPD": "chronic_disease",
    "DIABETES": "chronic_disease", "CHD": "chronic_disease", "STROKE": "chronic_disease",
    "BPHIGH": "chronic_disease",
}
RELEVANT_CATEGORIES: dict[str, str] = {
    "DISABLT": "disability", "HLTHSTAT": "health_status", "SOCLNEED": "social_needs",
}
LONG_COVID_PATTERN = re.compile(r"covid|sars-cov|coronavirus|post-acute|\bpasc\b|long[- ]haul", re.IGNORECASE)
RELEASE_PATTERN = re.compile(
    r"^PLACES: Local Data for Better Health, (County|ZCTA|Census Tract|Place) Data,? (\d{4}) release$")
CT_PLANNING_REGIONS = {f"09{c}" for c in range(110, 200, 10)}
CT_LEGACY = {f"09{c:03d}" for c in range(1, 16, 2)}
ZCTA_COLUMNS = ("year,locationid,locationname,datasource,categoryid,measureid,datavaluetypeid,"
                "data_value_unit,data_value,data_value_footnote_symbol,data_value_footnote,"
                "low_confidence_limit,high_confidence_limit,totalpopulation,totalpop18plus")


# --------------------------------------------------------------------------- pure helpers
def parse_release_name(name: str) -> tuple[str, int] | None:
    """'PLACES: Local Data for Better Health, County Data, 2025 release' -> ('county', 2025)."""
    m = RELEASE_PATTERN.match(name.strip())
    if not m:
        return None
    return m.group(1).lower().replace(" ", "_"), int(m.group(2))


def latest_release(catalog_rows: list[tuple[str, str]]) -> dict:
    """Pick the newest release year that has both a County and a ZCTA dataset.

    catalog_rows: (dataset_id, dataset_name). Returns {'year', 'county', 'zcta'}.
    """
    by_year: dict[int, dict[str, str]] = {}
    for dsid, name in catalog_rows:
        parsed = parse_release_name(name)
        if parsed:
            level, year = parsed
            by_year.setdefault(year, {})[level] = dsid
    complete = [y for y, d in by_year.items() if {"county", "zcta"} <= set(d)]
    if not complete:
        raise RuntimeError("no PLACES release with both County and ZCTA datasets found in the catalog")
    year = max(complete)
    return {"year": year, "county": by_year[year]["county"], "zcta": by_year[year]["zcta"]}


def relevance(measure_id: str, category_id: str) -> tuple[bool, str]:
    if measure_id in RELEVANT_MEASURES:
        return True, RELEVANT_MEASURES[measure_id]
    if category_id in RELEVANT_CATEGORIES:
        return True, RELEVANT_CATEGORIES[category_id]
    return False, "other"


def long_covid_measures(measures: pd.DataFrame) -> list[str]:
    """Measure ids whose id, label or question text mentions COVID / post-acute sequelae."""
    text = measures.astype(str).agg(" ".join, axis=1)
    hits = measures.loc[text.str.contains(LONG_COVID_PATTERN), "measure_id"]
    return sorted(hits.unique().tolist())


def ct_vintage(fips: pd.Series) -> str:
    """Which Connecticut county-equivalent set a FIPS column uses."""
    s = set(fips.astype(str))
    has_new, has_old = bool(s & CT_PLANNING_REGIONS), bool(s & CT_LEGACY)
    if has_new and has_old:
        return "mixed"
    if has_new:
        return "ct_planning_regions_09110_09190"
    if has_old:
        return "ct_legacy_counties_09001_09015"
    return "no_ct_rows"


def brfss_years_from_description(desc: str) -> str:
    m = re.search(r"uses (\d{4}) BRFSS data for (\d+) measures and (\d{4}) BRFSS data for (\d+) measures", desc)
    if m:
        return f"BRFSS {m.group(1)} ({m.group(2)} measures); BRFSS {m.group(3)} ({m.group(4)} measures)"
    years = sorted(set(re.findall(r"BRFSS\)? (\d{4})", desc)))
    return "BRFSS " + ", ".join(years) if years else "UNKNOWN"


def _num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


# --------------------------------------------------------------------------- retrieval
def discover() -> dict:
    rows: list[tuple[str, str]] = []
    for offset in (0, 100, 200):
        d = get_json(CATALOG_URL, params={"domains": DOMAIN, "q": "PLACES", "limit": 100, "offset": offset})
        rows += [(r["resource"]["id"], r["resource"]["name"]) for r in d["results"]]
        if len(d["results"]) < 100:
            break
    rel = latest_release(rows)
    rel["catalog_hits"] = len(rows)
    years = sorted({p[1] for p in (parse_release_name(n) for _, n in rows) if p})
    rel["release_years"] = f"{years[0]}-{years[-1]}" if years else "UNKNOWN"
    return rel


def fetch_metadata(dsid: str) -> dict:
    p = download_file(f"https://{DOMAIN}/api/views/{dsid}.json", SOURCE_ID, f"metadata_{dsid}.json")
    return json.loads(p.read_text())


def fetch_county(dsid: str, year: int) -> Path:
    return download_file(f"https://{DOMAIN}/api/views/{dsid}/rows.csv?accessType=DOWNLOAD", SOURCE_ID,
                         f"places_county_{year}_{dsid}.csv", min_bytes=1_000_000)


def zcta_query(measure_ids: list[str]) -> dict:
    return {"$select": ZCTA_COLUMNS,
            "$where": "measureid in(" + ",".join(repr(m) for m in sorted(measure_ids)) + ")",
            "$order": "locationid,measureid,datavaluetypeid", "$limit": "5000000"}


def fetch_zcta(dsid: str, year: int, measure_ids: list[str]) -> tuple[Path, int]:
    q = zcta_query(measure_ids)
    expected = int(get_json(f"https://{DOMAIN}/resource/{dsid}.json",
                            params={"$select": "count(*)", "$where": q["$where"]})[0]["count"])
    url = f"https://{DOMAIN}/resource/{dsid}.csv?" + urlencode(q)
    p = download_file(url, SOURCE_ID, f"places_zcta_{year}_{dsid}_selected_measures.csv",
                      min_bytes=1_000_000, timeout=1800)
    return p, expected


# --------------------------------------------------------------------------- tidy
def tidy_county(raw: pd.DataFrame, geos: pd.DataFrame, release_year: int) -> pd.DataFrame:
    cty = geos[geos["geo_level"] == "county"][["geo_id", "vintage", "state_fips"]]
    st = geos[geos["geo_level"] == "state"][["state_abbr", "state_fips"]].drop_duplicates()
    out = pd.DataFrame({
        "geo_id": raw["LocationID"].astype(str),
        "geo_name": raw["LocationName"],
        "state_abbr": raw["StateAbbr"],
        "state_name": raw["StateDesc"],
        "measure_id": raw["MeasureId"],
        "measure_label": raw["Measure"],
        "short_question_text": raw["Short_Question_Text"],
        "category_id": raw["CategoryID"],
        "category": raw["Category"],
        "data_value_type_id": raw["DataValueTypeID"],
        "data_value_type": raw["Data_Value_Type"],
        "value": _num(raw["Data_Value"]),
        "value_unit": raw["Data_Value_Unit"],
        "ci_low": _num(raw["Low_Confidence_Limit"]),
        "ci_high": _num(raw["High_Confidence_Limit"]),
        "total_population": _num(raw["TotalPopulation"]).astype("Int64"),
        "total_pop_18plus": _num(raw["TotalPop18plus"]).astype("Int64"),
        "brfss_year": _num(raw["Year"]).astype("Int64"),
        "footnote_symbol": raw["Data_Value_Footnote_Symbol"],
        "footnote": raw["Data_Value_Footnote"],
    })
    national = out["state_abbr"].eq("US")
    out["geo_level"] = national.map({True: "national", False: "county"})
    out.loc[national, "geo_id"] = "US"
    out.loc[national, "geo_name"] = out.loc[national, "state_name"]  # raw LocationName is blank on the US row
    out["county_fips"] = out["geo_id"].where(~national)
    out = out.merge(st, on="state_abbr", how="left")
    out = out.merge(cty.rename(columns={"vintage": "geo_vintage", "state_fips": "_sf"}), on="geo_id", how="left")
    out["in_canonical_geographies"] = out["geo_vintage"].eq("2024")
    out["geo_vintage"] = out["geo_vintage"].fillna("not_in_geographies").where(~national, "national")
    out = out.drop(columns="_sf")
    out["places_release"] = str(release_year)
    out["suppressed"] = out["value"].isna()
    rel = [relevance(m, c) for m, c in zip(out["measure_id"], out["category_id"])]
    out["relevant_measure"] = [r[0] for r in rel]
    out["relevance_group"] = [r[1] for r in rel]
    out["estimate_kind"] = "modeled_small_area_estimate"
    return out


def tidy_zcta(raw: pd.DataFrame, labels: pd.DataFrame, zc: pd.DataFrame, release_year: int) -> pd.DataFrame:
    out = pd.DataFrame({
        "geo_id": raw["locationid"].astype(str).str.zfill(5),
        "measure_id": raw["measureid"],
        "category_id": raw["categoryid"],
        "data_value_type_id": raw["datavaluetypeid"],
        "value": _num(raw["data_value"]),
        "value_unit": raw["data_value_unit"],
        "ci_low": _num(raw["low_confidence_limit"]),
        "ci_high": _num(raw["high_confidence_limit"]),
        "total_population": _num(raw["totalpopulation"]).astype("Int64"),
        "total_pop_18plus": _num(raw["totalpop18plus"]).astype("Int64"),
        "brfss_year": _num(raw["year"]).astype("Int64"),
        "footnote_symbol": raw["data_value_footnote_symbol"],
    })
    out["geo_level"] = "zcta"
    out["zcta"] = out["geo_id"]
    out = out.merge(labels[["measure_id", "measure_label", "short_question_text", "category"]],
                    on="measure_id", how="left")
    z = zc[["zcta", "state_fips", "state_abbr", "county_fips"]].rename(columns={"zcta": "geo_id"})
    out = out.merge(z, on="geo_id", how="left")
    out["in_canonical_zctas"] = out["geo_id"].isin(zc["zcta"])
    out["county_assignment"] = "2024 county containing the ZCTA internal point (zcta_centroids); ecological, not a split"
    out["places_release"] = str(release_year)
    out["suppressed"] = out["value"].isna()
    rel = [relevance(m, c) for m, c in zip(out["measure_id"], out["category_id"])]
    out["relevant_measure"] = [r[0] for r in rel]
    out["relevance_group"] = [r[1] for r in rel]
    out["estimate_kind"] = "modeled_small_area_estimate"
    return out


# --------------------------------------------------------------------------- audit
def _pct(n: int, d: int) -> str:
    return f"{n:,} ({100 * n / d:.2f}%)" if d else f"{n:,}"


def missingness_table(df: pd.DataFrame, cols: list[str]) -> str:
    lines = ["| column | missing | of rows |", "|---|---|---|"]
    for c in cols:
        lines.append(f"| {c} | {_pct(int(df[c].isna().sum()), len(df))} | {len(df):,} |")
    return "\n".join(lines)


def write_audit(stats: dict) -> Path:
    s = stats
    measure_rows = "\n".join(
        f"| {r.measure_id} | {r.category_id} | {r.measure_label} | {r.brfss_year} | {r.n_counties:,} | "
        f"{'yes' if r.relevant_measure else 'no'} | {r.relevance_group} |"
        for r in s["measure_table"].itertuples())
    files = "\n".join(f"| {k} | {v['url'][:160]}{'...' if len(v['url']) > 160 else ''} | {v['bytes']:,} | "
                      f"`{v['sha256'][:16]}...` | {v['retrieved_at']} |" for k, v in s["manifest"].items())
    text = f"""# DATA AUDIT — CDC PLACES (county + ZCTA), {s['release_year']} release

| Field | Value |
|---|---|
| source_id | cdc_places |
| Source (dataset/API name, exact files/endpoints) | "{s['county_name']}" (data.cdc.gov `{s['county_id']}`, full CSV export); "{s['zcta_name']}" (`{s['zcta_id']}`, SoQL extract of {s['n_zcta_measures']} flagged measures); PLACES data dictionary (`{DICTIONARY_ID}`) |
| Publishing organization | CDC, Division of Population Health (NCCDPHP); funded with the Robert Wood Johnson Foundation / CDC Foundation |
| Retrieval date (UTC) | {s['retrieved_at']} |
| Source version / release | PLACES {s['release_year']} release. {s['brfss_years']} (parsed from the dataset description). County rows updated {s['county_rows_updated']}; ZCTA rows updated {s['zcta_rows_updated']} (Socrata `rowsUpdatedAt`). |
| Source update date / cadence | Annual releases (the discovery API returned {s['catalog_hits']} PLACES/500 Cities catalog entries; PLACES county/ZCTA/tract/place releases found: {s['release_years']}; the newest with both County and ZCTA data is used) |
| License / access conditions | {s['license']}; open download, no registration |
| Unit of observation | county (or county-equivalent) x measure x value type (crude / age-adjusted prevalence); ZCTA x measure (crude only) |
| Sample size (actual, as ingested) | county file: {s['county_raw_rows']:,} rows, {s['n_counties']:,} counties + 1 national row, {s['n_measures']} measures; ZCTA extract: {s['zcta_raw_rows']:,} rows (server count {s['zcta_expected_rows']:,}), {s['n_zctas']:,} ZCTAs, {s['n_zcta_measures']} measures |
| Geography (resolution, vintage) | County: {s['county_ct_vintage']}; {s['county_match_2024']:,}/{s['n_counties']:,} county ids match `geographies` 2024 counties, {s['county_match_legacy']} match legacy CT, {s['county_unmatched']} unmatched. ZCTA: {s['zcta_match']:,}/{s['n_zctas']:,} match `zcta_centroids` (2020-vintage ZCTAs in the 2024 Gazetteer) |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — ecological model-based estimates; no person identifiers |

## What these numbers are
**Modeled small-area estimates** (`evidence_type = modeled_small_area_estimate`). CDC fits multilevel
models to BRFSS responses and post-stratifies them onto Census population and ACS covariates. The
dataset description states: "Because the small area model cannot detect effects due to local
interventions, users are cautioned against using these estimates for program or policy evaluations."
Values are **population context** (and at most labelled proxies of symptom/functional burden).
They are not county prevalence of Long COVID, ME/CFS, POTS, Lyme disease or any other target
condition. Model-based county estimates borrow strength from state-level BRFSS and county
covariates, so between-county differences partly reflect demographic composition.

## Does the release contain a long-COVID measure?
**{s['long_covid_answer']}.** Evidence (measured on the files retrieved): regex
`{LONG_COVID_PATTERN.pattern}` matched {s['long_covid_hits']} of the {s['n_measures']} measure ids/labels/question texts
in the county file; the string "COVID" occurs in {s['covid_string_hits']} of {s['county_raw_rows']:,} county rows;
the PLACES data dictionary (`{DICTIONARY_ID}`, {s['dict_rows']} measures across releases 2016-2024) matched
{s['dict_hits']}. The measure list is below.

## Files / endpoints retrieved
Discovery: `{CATALOG_URL}?domains={DOMAIN}&q=PLACES` (limit 100, paged). Metadata: `https://{DOMAIN}/api/views/<id>.json`.

| file | url | bytes | sha256 | retrieved_at |
|---|---|---|---|---|
{files}

ZCTA SoQL: `$select={ZCTA_COLUMNS}`; `$where=measureid in(<{s['n_zcta_measures']} flagged ids>)`; `$order=locationid,measureid,datavaluetypeid`; `$limit=5000000`. The row count of the
extract ({s['zcta_raw_rows']:,}) was checked against `count(*)` with the same filter ({s['zcta_expected_rows']:,}).

## Key variables
| processed column | raw | meaning / units | use downstream |
|---|---|---|---|
| geo_id / county_fips | LocationID | 5-char county FIPS (US national row: geo_id `US`, raw LocationID `59`) | join to `geographies` (ecological) |
| measure_id, measure_label, short_question_text | MeasureId, Measure, Short_Question_Text | PLACES measure | filter |
| category_id | CategoryID | HLTHOUT, HLTHSTAT, DISABLT, PREVENT, RISKBEH, SOCLNEED | grouping |
| data_value_type_id | DataValueTypeID | CrdPrv = crude prevalence, AgeAdjPrv = age-adjusted prevalence (county only) | use AgeAdjPrv for between-county comparison, CrdPrv for counts |
| value, ci_low, ci_high | Data_Value, Low/High_Confidence_Limit | percent of adults, 95% CI of the model estimate | context variable |
| total_population, total_pop_18plus | TotalPopulation, TotalPop18plus | Census population used by PLACES | denominators |
| brfss_year | Year | BRFSS data year behind the estimate (2023 or 2022, see table) | vintage |
| suppressed, footnote_symbol, footnote | Data_Value_Footnote* | `*` = suppressed (population < 50); `#` = social-needs estimate based on 39 states + DC | missingness |
| relevant_measure, relevance_group | derived | flag for the invisible-illness context set (SPEC section E); whole DISABLT/HLTHSTAT/SOCLNEED categories flagged | selection |

## Measures ({s['n_measures']}) and county coverage
| measure_id | category | label | BRFSS year | counties with rows | flagged | group |
|---|---|---|---|---|---|---|
{measure_rows}

## Missingness
Coverage gaps are rows that are **absent**, not blank values:
* Kentucky and Pennsylvania: no rows for any measure built on 2023 BRFSS ({s['n_counties_2023']:,} of {s['n_counties']:,} counties have them;
  {s['ky_pa_counties']} KY+PA counties lack them). CDC PLACES release notes (www.cdc.gov/places/current-release-notes):
  "Estimates for Kentucky and Pennsylvania were not available for measures based on the 2023 BRFSS. However, the
  five measures based on the 2022 BRFSS were carried over for Kentucky and Pennsylvania."
* Health-related social needs (SOCLNEED): {s['n_counties_soclneed']:,} counties in {s['n_states_soclneed']} states + DC (footnote `#`:
  "The estimates are based on 39 states and DC."). Counties in the other states have no social-needs rows.
* Territories: {s['territory_counties_missing']} county-equivalents in `geographies` (PR, GU, VI, AS, MP) are not in PLACES.
* Suppressed values (county table): {_pct(s['county_suppressed'], s['county_rows'])} (`*`, population < 50).
* ZCTA: {s['zcta_with_2023']:,} of {s['n_zctas']:,} ZCTAs have 2023-BRFSS measures (KY/PA ZCTAs lack them);
  {s['zcta_with_soclneed']:,} have social-needs measures; {s['zcta_centroids_missing']:,} ZCTAs in `zcta_centroids` are absent from PLACES;
  {s['zcta_no_county']} PLACES ZCTAs have no `county_fips`/`state_fips` because their Census internal point lies outside every
  2024 county polygon in `zcta_centroids` (the PLACES value itself is present).

County table (all rows):

{s['county_missing_md']}

ZCTA table:

{s['zcta_missing_md']}

## Linkage strategy
* County: `geo_id` = 5-char FIPS joins `geographies.geo_id` (2024 vintage). Connecticut: PLACES {s['release_year']} uses the
  **9 planning regions (09110-09190)**; {s['ct_rows_match']} CT ids match `geographies` 2024 rows; 0 legacy-county rows.
  Row level: {s['rows_match_2024']:,} of {s['rows_county']:,} county rows match a 2024 county in `geographies`
  ({s['rows_ct']:,} Connecticut rows, all on planning regions); {s['rows_legacy']} rows match legacy CT counties.
* ZCTA: `geo_id` = ZCTA5 joins `zcta_centroids.zcta`; `county_fips` is the 2024 county containing the ZCTA internal point
  (a drill-down convenience; ZCTAs can straddle counties and are not apportioned here).
* Joins to any other layer are **ecological**. Nothing here describes individual people.

## Limitations and caveats
* Modeled, not observed: estimates are synthetic predictions with model CIs; they smooth toward covariate-predicted values.
* Two BRFSS years are mixed within one release (2023 for most measures; 2022 for SLEEP, DENTAL, MAMMOUSE, COLON_SCREEN, TEETHLOST).
* Measure definitions change between releases (e.g. LONELINESS replaced ISOLATION in the 2025 release per CDC release notes;
  the data dictionary `{DICTIONARY_ID}` predates the 2025 release and lists ISOLATION, not LONELINESS). Do not stack releases.
* Crude vs age-adjusted: age-adjusted values are for comparison between areas; crude values x population give counts.
* No invisible-illness condition measure: PHLTH, GHLTH, DISABILITY etc. are symptom/function context. If they are ever used as
  a burden proxy, the row must carry `burden_evidence_level = C` and be labelled as a proxy.

## Processed outputs
| table | rows |
|---|---|
| geo_context__cdc_places | {s['county_rows']:,} |
| geo_context_zcta__cdc_places | {s['zcta_rows']:,} |

## Reproduce
`uv run python -m measure_it.ingestion.cdc_places`
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


# --------------------------------------------------------------------------- main
def run() -> dict:
    rel = discover()
    year, cid, zid = rel["year"], rel["county"], rel["zcta"]
    cmeta, zmeta = fetch_metadata(cid), fetch_metadata(zid)
    dict_path = download_file(f"https://{DOMAIN}/api/views/{DICTIONARY_ID}/rows.csv?accessType=DOWNLOAD",
                              SOURCE_ID, f"places_data_dictionary_{DICTIONARY_ID}.csv")
    county_path = fetch_county(cid, year)
    manifest = load_manifest(SOURCE_ID)["files"]
    retrieved = manifest[county_path.name]["retrieved_at"]
    brfss_years = brfss_years_from_description(cmeta.get("description", ""))
    version = f"PLACES {year} release (data.cdc.gov {cid} / {zid}); {brfss_years}"

    raw = pd.read_csv(county_path, dtype=str)
    geos = read_table("geographies")
    county = tidy_county(raw, geos, year)

    measures = (county[["measure_id", "measure_label", "short_question_text", "category_id", "category"]]
                .drop_duplicates("measure_id"))
    lc_hits = long_covid_measures(measures)
    dictionary = pd.read_csv(dict_path, dtype=str)
    dict_hits = int(dictionary.apply(lambda c: c.str.contains(LONG_COVID_PATTERN, na=False)).any(axis=1).sum())
    covid_rows = int(raw.apply(lambda c: c.str.contains("COVID", case=False, na=False)).any(axis=1).sum())

    zcta_ids = sorted(county.loc[county["relevant_measure"], "measure_id"].unique())
    zcta_path, expected = fetch_zcta(zid, year, zcta_ids)
    manifest = load_manifest(SOURCE_ID)["files"]
    zraw = pd.read_csv(zcta_path, dtype=str)
    if len(zraw) != expected:
        raise RuntimeError(f"ZCTA extract has {len(zraw)} rows, server count is {expected}; delete the file and rerun")
    zc = read_table("zcta_centroids")
    zcta = tidy_zcta(zraw, measures, zc, year)

    notes = ("Modeled small-area estimate (BRFSS MRP); population context, not observed condition prevalence. "
             "PLACES has no long-COVID/ME-CFS/POTS/Lyme measure.")
    county_out = add_provenance(
        county, data_layer="geographic", source_name=f"{SOURCE_NAME}, County Data {year} release",
        source_version=version, retrieved_at=retrieved, evidence_type="modeled_small_area_estimate",
        source_record_id=lambda d: f"{cid}:" + d["geo_id"] + ":" + d["measure_id"] + ":" + d["data_value_type_id"],
        source_geographic_resolution="county", evidence_level="population_context", provenance_notes=notes)
    county_out.loc[county_out["geo_level"] == "national", "source_geographic_resolution"] = "national"
    write_table(county_out, "geo_context__cdc_places", producer=PRODUCER,
                description=f"PLACES {year} county long format: all {len(measures)} measures, crude + age-adjusted, modeled")

    zcta_out = add_provenance(
        zcta, data_layer="geographic", source_name=f"{SOURCE_NAME}, ZCTA Data {year} release",
        source_version=version, retrieved_at=manifest[zcta_path.name]["retrieved_at"],
        evidence_type="modeled_small_area_estimate",
        source_record_id=lambda d: f"{zid}:" + d["geo_id"] + ":" + d["measure_id"] + ":" + d["data_value_type_id"],
        source_geographic_resolution="zcta", evidence_level="population_context", provenance_notes=notes)
    write_table(zcta_out, "geo_context_zcta__cdc_places", producer=PRODUCER,
                description=f"PLACES {year} ZCTA crude prevalence for {len(zcta_ids)} flagged measures, modeled")

    # ---- measured stats for the audit and registry
    cty = county[county["geo_level"] == "county"]
    gc = geos[geos["geo_level"] == "county"]
    ids = pd.Series(cty["geo_id"].unique())
    mt = (cty.groupby(["measure_id", "category_id", "measure_label", "relevant_measure", "relevance_group"])
          .agg(brfss_year=("brfss_year", "max"), n_counties=("geo_id", "nunique")).reset_index()
          .sort_values(["category_id", "measure_id"]))
    n23 = int(cty[cty["brfss_year"] == 2023]["geo_id"].nunique())
    sn = cty[cty["category_id"] == "SOCLNEED"]
    ts = lambda m: pd.Timestamp(m.get("rowsUpdatedAt", 0), unit="s").isoformat() if m.get("rowsUpdatedAt") else "UNKNOWN"
    stats = {
        "release_year": year, "county_id": cid, "zcta_id": zid, "county_name": cmeta["name"], "zcta_name": zmeta["name"],
        "catalog_hits": rel["catalog_hits"], "release_years": rel["release_years"],
        "zcta_no_county": int(zcta.loc[zcta["county_fips"].isna(), "geo_id"].nunique()), "brfss_years": brfss_years, "retrieved_at": retrieved,
        "county_rows_updated": ts(cmeta), "zcta_rows_updated": ts(zmeta),
        "license": (cmeta.get("license") or {}).get("name", "UNKNOWN"),
        "county_raw_rows": len(raw), "n_counties": len(ids), "n_measures": len(measures),
        "zcta_raw_rows": len(zraw), "zcta_expected_rows": expected, "n_zctas": zcta["geo_id"].nunique(),
        "n_zcta_measures": len(zcta_ids),
        "county_ct_vintage": ct_vintage(ids),
        "county_match_2024": int(ids.isin(gc[~gc["ct_legacy"]]["geo_id"]).sum()),
        "county_match_legacy": int(ids.isin(gc[gc["ct_legacy"]]["geo_id"]).sum()),
        "county_unmatched": int((~ids.isin(gc["geo_id"])).sum()),
        "ct_rows_match": int(ids[ids.str.startswith("09")].isin(gc[~gc["ct_legacy"]]["geo_id"]).sum()),
        "rows_county": len(cty), "rows_match_2024": int(cty["geo_id"].isin(gc[~gc["ct_legacy"]]["geo_id"]).sum()),
        "rows_legacy": int(cty["geo_id"].isin(gc[gc["ct_legacy"]]["geo_id"]).sum()),
        "rows_ct": int(cty["geo_id"].str.startswith("09").sum()),
        "zcta_match": int(pd.Series(zcta["geo_id"].unique()).isin(zc["zcta"]).sum()),
        "long_covid_hits": len(lc_hits), "long_covid_answer": "Yes: " + ", ".join(lc_hits) if lc_hits else "No",
        "covid_string_hits": covid_rows, "dict_rows": len(dictionary), "dict_hits": dict_hits,
        "measure_table": mt, "n_counties_2023": n23,
        "ky_pa_counties": int(gc[~gc["ct_legacy"] & gc["state_abbr"].isin(["KY", "PA"])].shape[0]),
        "n_counties_soclneed": int(sn["geo_id"].nunique()), "n_states_soclneed": int(sn.loc[sn["state_abbr"] != "DC", "state_abbr"].nunique()),
        "territory_counties_missing": int((~gc[~gc["ct_legacy"]]["geo_id"].isin(ids)).sum()),
        "county_suppressed": int(county["suppressed"].sum()), "county_rows": len(county),
        "zcta_rows": len(zcta),
        "zcta_with_2023": int(zcta[zcta["brfss_year"] == 2023]["geo_id"].nunique()),
        "zcta_with_soclneed": int(zcta[zcta["category_id"] == "SOCLNEED"]["geo_id"].nunique()),
        "zcta_centroids_missing": int((~zc["zcta"].isin(zcta["geo_id"])).sum()),
        "manifest": manifest,
        "county_missing_md": missingness_table(county, ["value", "ci_low", "ci_high", "total_population",
                                                        "total_pop_18plus", "brfss_year", "state_fips", "county_fips"]),
        "zcta_missing_md": missingness_table(zcta, ["value", "ci_low", "ci_high", "total_population",
                                                    "state_fips", "county_fips", "measure_label"]),
    }
    write_audit(stats)
    write_registry_entry({
        "source_id": SOURCE_ID, "name": f"{SOURCE_NAME} ({year} release: county + ZCTA)",
        "publisher": "CDC Division of Population Health (NCCDPHP)", "landing_url": LANDING_URL,
        "access_urls": [f"https://{DOMAIN}/d/{cid}", f"https://{DOMAIN}/d/{zid}", f"https://{DOMAIN}/d/{DICTIONARY_ID}",
                        CATALOG_URL + f"?domains={DOMAIN}&q=PLACES"],
        "license": stats["license"], "access_conditions": "open download, no registration",
        "retrieved_at": retrieved, "source_version": version,
        "update_date": f"county rows {stats['county_rows_updated']}; zcta rows {stats['zcta_rows_updated']}; annual releases",
        "data_layer": "geographic",
        "unit_of_observation": "county x measure x value type (crude/age-adjusted); ZCTA x measure (crude)",
        "sample_size": {"county_rows": len(county), "counties": len(ids), "measures": len(measures),
                        "zcta_rows": len(zcta), "zctas": stats["n_zctas"], "zcta_measures": len(zcta_ids),
                        "counties_with_2023_brfss_measures": n23,
                        "counties_with_social_needs": stats["n_counties_soclneed"]},
        "geographic_resolution": "county (2024 vintage, CT planning regions); zcta; national reference row",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (ecological modeled estimates)",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": ["geo_context__cdc_places", "geo_context_zcta__cdc_places"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": PRODUCER,
        "evidence_type": "modeled_small_area_estimate",
        "contains_long_covid_measure": bool(lc_hits),
        "ct_county_vintage": stats["county_ct_vintage"],
        "limitations": [
            "Modeled small-area estimates (BRFSS MRP), not observed prevalence; not for program evaluation",
            "No long-COVID, ME/CFS, POTS or Lyme measure; symptom/function measures are context or labelled level-C proxies",
            "Kentucky and Pennsylvania lack all 2023-BRFSS measures in the 2025 release",
            "Social-needs measures cover 39 states + DC only",
            "Two BRFSS years mixed within the release (2023 for 35 measures, 2022 for 5)",
        ],
    })
    summary = {k: stats[k] for k in ("release_year", "county_id", "zcta_id", "county_rows", "zcta_rows", "n_counties",
                                     "n_zctas", "n_measures", "long_covid_answer", "county_ct_vintage")}
    print(json.dumps(summary, indent=2, default=str))
    return summary


if __name__ == "__main__":
    run()
