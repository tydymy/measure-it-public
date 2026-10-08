"""CDC/ATSDR Social Vulnerability Index (county) -> geo_vulnerability.

Release detection: the newest SVI county CSV served by the official download
service (svi2.cdc.gov/webapi, the back end of the ATSDR "SVI Data &
Documentation Download" page). On 2026-09-23 the 2024 request returned HTTP 404
and 2022 (ACS 2018-2022) was the newest release.

Outputs
  geo_vulnerability              U.S. county database (50 states + DC), counties
                                 ranked against all U.S. counties.
  geo_vulnerability_puerto_rico  Puerto Rico county database. CDC ranks Puerto
                                 Rico municipios only against each other, so
                                 its percentiles are NOT comparable with the
                                 U.S. table and live in a separate table.

State-level SVI: CDC does not publish SVI with the state as the unit. The
per-state "county" downloads rank a state's counties against each other; they
are not state-level indices and are not ingested.

Columns keep SVI's own names in lower case so they trace to the SVI 2022 data
dictionary: rpl_themes (overall percentile), rpl_theme1 socioeconomic status,
rpl_theme2 household characteristics, rpl_theme3 racial & ethnic minority
status, rpl_theme4 housing type & transportation; their component spl_/epl_
columns, the e_/ep_/mp_ inputs, and the f_ flags. -999 (and ACS special
values <= -999) are converted to null and counted.

SVI measures relative social vulnerability for emergency planning. It is not
disease prevalence or burden and must not be used as a proxy for either.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import requests

from ..config import raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..http import HTMLInsteadOfData
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from .census_acs import county_join_coverage

SOURCE_ID = "cdc_svi"
LANDING_URL = "https://www.atsdr.cdc.gov/place-health/php/svi/svi-data-documentation-download.html"
DOWNLOAD = "https://svi2.cdc.gov/webapi/Documents/download?year={year}&type=csv&category=states_counties&name={name}"
DOC = "https://svi2.cdc.gov/webapi/documents/publication?filename=SVI{year}Documentation_ZCTA"
YEAR_CANDIDATES = (2024, 2022)
ACS_PERIOD = {2022: "2018-2022", 2024: "2020-2024"}

THEMES = {
    "rpl_theme1": "Socioeconomic Status",
    "rpl_theme2": "Household Characteristics",
    "rpl_theme3": "Racial & Ethnic Minority Status",
    "rpl_theme4": "Housing Type & Transportation",
    "rpl_themes": "Overall",
}
# The 16 SVI variables by theme (SVI 2022 data dictionary) + the internet adjunct.
VARIABLES = {
    1: ["POV150", "UNEMP", "HBURD", "NOHSDP", "UNINSUR"],
    2: ["AGE65", "AGE17", "DISABL", "SNGPNT", "LIMENG"],
    3: ["MINRTY"],
    4: ["MUNIT", "MOBILE", "CROWD", "NOVEH", "GROUPQ"],
}
ADJUNCT = ["NOINT"]
ID_COLS = ["ST", "STATE", "ST_ABBR", "STCNTY", "COUNTY", "FIPS", "LOCATION", "AREA_SQMI"]


def keep_columns() -> list[str]:
    cols = ID_COLS + ["E_TOTPOP", "M_TOTPOP", "E_HU", "E_HH", "E_DAYPOP"]
    for t, vs in VARIABLES.items():
        for v in vs:
            cols += [f"E_{v}", f"EP_{v}", f"MP_{v}", f"EPL_{v}", f"F_{v}"]
        cols += [f"SPL_THEME{t}", f"RPL_THEME{t}", f"F_THEME{t}"]
    for v in ADJUNCT:
        cols += [f"E_{v}", f"EP_{v}", f"MP_{v}"]
    cols += ["SPL_THEMES", "RPL_THEMES", "F_TOTAL"]
    return cols


def fetch(years=YEAR_CANDIDATES) -> tuple[int, dict[str, Path], dict]:
    """Download the newest available release: U.S. county CSV, PR county CSV, documentation PDF."""
    probes = {}
    for year in years:
        us_name = f"SVI_{year}_US_COUNTY"
        try:
            us = download_file(DOWNLOAD.format(year=year, name=us_name), SOURCE_ID, f"{us_name}.csv",
                               min_bytes=100_000)
        except (requests.RequestException, HTMLInsteadOfData, RuntimeError) as exc:
            # A newer release that is absent (404), unreachable (offline rebuild) or served as an HTML page
            # is recorded and skipped, so a cached older release still rebuilds.
            resp = getattr(exc, "response", None)
            probes[year] = {"url": DOWNLOAD.format(year=year, name=us_name),
                            "status": resp.status_code if resp is not None else None,
                            "error": None if resp is not None else f"{type(exc).__name__}: {exc}"[:300],
                            "checked_at": utc_now_iso()}
            continue
        probes[year] = {"url": DOWNLOAD.format(year=year, name=us_name), "status": 200}
        pr = download_file(DOWNLOAD.format(year=year, name="PUERTORICO_COUNTY"), SOURCE_ID,
                           f"SVI_{year}_PUERTORICO_COUNTY.csv", min_bytes=10_000)
        doc = download_file(DOC.format(year=year), SOURCE_ID, f"SVI{year}Documentation_ZCTA.pdf",
                            min_bytes=10_000, allow_html=False)
        return year, {"us": us, "pr": pr, "doc": doc}, probes
    raise RuntimeError(f"no SVI county release found: {probes}")


def read_svi(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"ST": str, "STCNTY": str, "FIPS": str})


def clean(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Keep the analysis columns; convert -999 and ACS special values (<= -999) to null; count them."""
    missing_cols = [c for c in keep_columns() if c not in df.columns]
    if missing_cols:
        raise ValueError(f"SVI file lacks expected columns {missing_cols}")
    counts_all = {"__n_cells_equal_minus999": int((df.select_dtypes("number") == -999).sum().sum())}
    for c in df.select_dtypes("number").columns:
        v = df[c]
        bad = v <= -999
        if bad.any():
            counts_all[c] = {str(int(k)): int(n) for k, n in v[bad].value_counts().items()}
    out = df[keep_columns()].copy()
    for c in out.select_dtypes("number").columns:
        out[c] = out[c].mask(out[c] <= -999)
    out["FIPS"] = out["FIPS"].str.zfill(5)
    out["ST"] = out["ST"].str.zfill(2)
    return out, counts_all


def transform(df: pd.DataFrame, year: int, universe: str) -> pd.DataFrame:
    out = df.rename(columns=str.lower).rename(columns={
        "fips": "geo_id", "st": "state_fips", "st_abbr": "state_abbr", "state": "state_name",
        "county": "county_name", "location": "location_name", "area_sqmi": "svi_area_sqmi"})
    out = out.drop(columns=["stcnty"]).copy()
    out.insert(1, "geo_level", "county")
    extra = pd.DataFrame({"svi_release": year, "acs_period": ACS_PERIOD.get(year, ""),
                          "svi_ranking_universe": universe}, index=out.index)
    out = pd.concat([out, extra], axis=1)
    lead = ["geo_id", "geo_level", "state_fips", "state_abbr", "state_name", "county_name", "location_name",
            "svi_release", "acs_period", "svi_ranking_universe", "rpl_themes", "rpl_theme1", "rpl_theme2",
            "rpl_theme3", "rpl_theme4"]
    return out[lead + [c for c in out.columns if c not in lead]]


def _provenance(df: pd.DataFrame, year: int, manifest_entry: dict, universe_note: str) -> pd.DataFrame:
    return add_provenance(
        df, data_layer="geographic",
        source_name=f"CDC/ATSDR Social Vulnerability Index {year}",
        source_version=f"SVI {year} (ACS {ACS_PERIOD.get(year, '?')} 5-year); county database",
        retrieved_at=manifest_entry["retrieved_at"], evidence_type="composite_index",
        source_record_id="geo_id", source_geographic_resolution="county",
        provenance_notes=("Percentile ranks (0-1, higher = more socially vulnerable) of 16 ACS variables in 4 themes. "
                          + universe_note + " Social vulnerability is context for emergency/health planning; it is "
                          "NOT disease prevalence or burden."),
    )


def _missing(df: pd.DataFrame, cols: list[str]) -> dict:
    return {c: {"n_missing": int(df[c].isna().sum()), "pct_missing": round(100 * df[c].isna().mean(), 3)}
            for c in cols}


def build() -> dict:
    year, files, probes = fetch()
    manifest = load_manifest(SOURCE_ID)["files"]
    raw_us, raw_pr = read_svi(files["us"]), read_svi(files["pr"])
    us, us_counts = clean(raw_us)
    pr, pr_counts = clean(raw_pr)
    us_t = transform(us, year, "US counties (50 states + DC)")
    pr_t = transform(pr, year, "Puerto Rico municipios only")

    us_out = _provenance(us_t, year, manifest[files["us"].name],
                         "Ranked against all U.S. counties (50 states + DC); Puerto Rico excluded from this ranking.")
    write_table(us_out, "geo_vulnerability", producer="ingestion.svi.build",
                description=f"CDC/ATSDR SVI {year} U.S. county database (overall + 4 theme percentiles, components)")
    pr_out = _provenance(pr_t, year, manifest[files["pr"].name],
                         "Ranked only among Puerto Rico municipios; NOT comparable with U.S.-ranked percentiles.")
    write_table(pr_out, "geo_vulnerability_puerto_rico", producer="ingestion.svi.build",
                description=f"CDC/ATSDR SVI {year} Puerto Rico county database (PR-only ranking)")

    rpl = list(THEMES)
    ep = [f"ep_{v.lower()}" for vs in VARIABLES.values() for v in vs] + ["ep_noint"]
    stats = {
        "release": year, "acs_period": ACS_PERIOD.get(year), "probes": probes,
        "files": {p.name: {k: manifest[p.name].get(k) for k in ("url", "bytes", "sha256", "retrieved_at",
                                                                "last_modified")} for p in files.values()},
        "raw_shape": {"us": list(raw_us.shape), "pr": list(raw_pr.shape)},
        "rows": {"us": int(len(us_out)), "pr": int(len(pr_out))},
        "columns_kept": int(len(keep_columns())),
        "states_in_us_file": int(us_t["state_abbr"].nunique()),
        "sentinels_all_columns": {"us": us_counts, "pr": pr_counts},
        "missing_us": _missing(us_t, rpl + ep + ["e_totpop"]),
        "missing_pr": _missing(pr_t, rpl),
        "rpl_summary_us": us_t[rpl].describe().round(4).to_dict(),
        "population_sum_us": int(us_t["e_totpop"].sum()),
        "population_sum_pr": int(pr_t["e_totpop"].sum()),
        "zero_population_counties_us": int((us_t["e_totpop"] == 0).sum()),
        "theme_correlations_us": us_t[rpl].corr(method="spearman").round(3).to_dict(),
        "coverage_us": county_join_coverage(us_t["geo_id"], f"SVI {year} US county"),
        "coverage_us_plus_pr": county_join_coverage(pd.concat([us_t["geo_id"], pr_t["geo_id"]]),
                                                    f"SVI {year} US + PR county"),
    }
    (raw_dir(SOURCE_ID) / "build_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    _registry(stats, year, manifest, files)
    write_audit(stats)
    return stats


def _registry(stats: dict, year: int, manifest: dict, files: dict) -> None:
    write_registry_entry({
        "source_id": SOURCE_ID, "name": f"CDC/ATSDR Social Vulnerability Index {year} (county)",
        "publisher": "CDC / ATSDR Geospatial Research, Analysis, and Services Program (GRASP)",
        "landing_url": LANDING_URL, "access_urls": [manifest[p.name]["url"] for p in files.values()],
        "license": "U.S. federal government work (public domain); CDC suggests citing the SVI database",
        "access_conditions": "open download, no registration",
        "retrieved_at": manifest[files["us"].name]["retrieved_at"],
        "source_version": f"SVI {year} (ACS {ACS_PERIOD.get(year)} 5-year); newer release probed: "
                          + ", ".join(f"{y}: HTTP {p['status']}" for y, p in stats["probes"].items()),
        "update_date": "biennial releases; ATSDR download page dated 2026-05-11 lists 2022 as the newest; "
                       "2022 MP_CROWD correction 2024-12-11",
        "data_layer": "geographic", "unit_of_observation": "county / county-equivalent",
        "sample_size": {"us_counties": stats["rows"]["us"], "pr_municipios": stats["rows"]["pr"]},
        "geographic_resolution": "county (2022 county set; CT planning regions)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable", "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": ["geo_vulnerability", "geo_vulnerability_puerto_rico"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": "measure_it.ingestion.svi",
        "limitations": ["Relative percentile ranks, not absolute levels; not comparable across SVI releases",
                        "Puerto Rico ranked separately and not comparable with U.S. ranks",
                        "No state-level SVI is published",
                        "Vulnerability is not disease prevalence or burden"],
    })


def _fmt_missing(d: dict) -> str:
    rows = ["| Variable | Missing (n) | Missing (%) |", "|---|---|---|"]
    return "\n".join(rows + [f"| {k} | {v['n_missing']} | {v['pct_missing']} |" for k, v in d.items()])


def write_audit(s: dict) -> Path:
    files = "\n".join(f"| {f} | {v['url']} | {v['bytes']:,} | `{v['sha256'][:16]}…` | {v['last_modified'] or 'not sent'} |"
                      for f, v in s["files"].items())
    probes = "; ".join(f"{y}: HTTP {p['status']} ({p['url']})" for y, p in s["probes"].items())
    sent_us = s["sentinels_all_columns"]["us"]
    sent_pr = s["sentinels_all_columns"]["pr"]
    cov, cov2 = s["coverage_us"], s["coverage_us_plus_pr"]
    rs = s["rpl_summary_us"]
    rpl_tab = "\n".join(f"| {c} | {THEMES[c]} | {rs[c]['count']:.0f} | {rs[c]['min']} | {rs[c]['50%']} | {rs[c]['max']} |"
                        for c in THEMES)
    var_rows = []
    for t, vs in VARIABLES.items():
        for v in vs:
            var_rows.append(f"| {THEMES[f'rpl_theme{t}']} | e_{v.lower()}, ep_{v.lower()}, mp_{v.lower()}, "
                            f"epl_{v.lower()}, f_{v.lower()} |")
    text = f"""# DATA AUDIT — CDC/ATSDR Social Vulnerability Index {s['release']} (county)

| Field | Value |
|---|---|
| source_id | cdc_svi |
| Source (dataset/API name, exact files/endpoints) | SVI {s['release']} U.S. county CSV and Puerto Rico county CSV from the ATSDR SVI download service (svi2.cdc.gov/webapi), plus the SVI {s['release']} documentation PDF |
| Publishing organization | CDC / ATSDR Geospatial Research, Analysis, and Services Program (GRASP) |
| Retrieval date (UTC) | {min(v['retrieved_at'] for v in s['files'].values())} |
| Source version / release | SVI {s['release']}, built on ACS {s['acs_period']} 5-year estimates |
| Source update date / cadence | Download page dated May 11, 2026 lists 2022 as the newest release. 2022 MP_CROWD values corrected 2024-12-11 (the file here was retrieved after that correction date). Newer release probe: {probes} |
| License / access conditions | Public domain (U.S. federal work); open download, no registration |
| Unit of observation | County / county-equivalent |
| Sample size (actual, as ingested) | U.S. file {s['rows']['us']} counties ({s['states_in_us_file']} state codes: 50 states + DC); Puerto Rico file {s['rows']['pr']} municipios |
| Geography (resolution, vintage) | County; ACS 2022 geography (Connecticut = 9 planning regions) |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — area-level index |

## Files / endpoints retrieved

| File | URL | Bytes | sha256 | Last-Modified |
|---|---|---|---|---|
{files}

The download page itself (`{LANDING_URL}`) builds these URLs in JavaScript
(`https://svi2.cdc.gov/data-downloads/assets/main-DVz7aho5.js`); the year menu lists 2022, 2020, 2018, 2016,
2014, 2010, 2000. The legacy path `https://svi.cdc.gov/Documents/Data/2022/csv/states_counties/SVI_2022_US_county.csv`
served a byte-identical file (same sha256) on the retrieval date.

**State-level SVI:** not published. The per-state downloads rank one state's counties against each other;
they are county rows, not a state index, and were not ingested. **ZCTA-level SVI 2022** exists on the same
service (not requested; not ingested).

## Key variables

| Column | Meaning | Count / min / median / max (U.S. file) |
|---|---|---|
{rpl_tab}

(The theme columns are rpl_theme1 socioeconomic status, rpl_theme2 household characteristics,
rpl_theme3 racial & ethnic minority status, rpl_theme4 housing type & transportation; rpl_themes is the overall rank.)

Percentile ranks run 0–1 (higher = more vulnerable). Each theme rank is the percentile of `spl_theme<k>`, the sum of
the variable percentiles `epl_<var>`; `rpl_themes` ranks `spl_themes`. All component columns are kept.

| Theme | Kept inputs per variable |
|---|---|
{chr(10).join(var_rows)}
| Adjunct | e_noint, ep_noint, mp_noint (households without an internet subscription) |

Also kept: e_totpop, m_totpop, e_hu, e_hh, e_daypop, svi_area_sqmi, f_theme1–4, f_total. `e_` = ACS count
estimate, `ep_` = percent, `mp_` = 90% MOE of the percent, `epl_` = percentile rank, `f_<var>` = flag (1 when
`epl_<var>` ≥ 0.90), `f_theme<k>` / `f_total` = number of flagged variables in the theme / overall (counts, not 0/1).
Adjunct race/ethnicity breakdowns (e_afam, e_hisp, …) were dropped; theme 3 keeps its single input ep_minrty.

## Missingness

`-999` values (SVI's "no data") and ACS special values ≤ −999 in any column of the raw files, before cleaning
(`__n_cells_equal_minus999` = number of cells exactly −999 across all numeric columns):

* U.S. file: {sent_us}
* Puerto Rico file: {sent_pr}

(`-555555555` in an MOE column is the ACS code for a controlled estimate, MOE not appropriate. M_GROUPQ is not in the
kept column set.) After cleaning, missing values in the U.S. table:

{_fmt_missing(s['missing_us'])}

Puerto Rico table percentiles:

{_fmt_missing(s['missing_pr'])}

Counties with e_totpop = 0 in the U.S. file: {s['zero_population_counties_us']}. Sum of e_totpop: U.S. {s['population_sum_us']:,};
Puerto Rico {s['population_sum_pr']:,}.

## Linkage strategy

* Join to `geographies` on `geo_id` (5-digit county FIPS). U.S. file: {cov['source_counties']} counties, matched
  {cov['matched']} of {cov['canonical_counties_2024']} canonical 2024 counties; SVI rows not in geographies =
  {cov['source_unmatched']} {cov['source_unmatched_fips']}; canonical counties without a U.S. SVI row =
  {cov['canonical_without_source_row']} (by state/territory: {cov['canonical_without_source_row_by_state']}).
  Adding the Puerto Rico table: matched {cov2['matched']}; still without SVI =
  {cov2['canonical_without_source_row']} ({cov2['canonical_without_source_row_by_state']}; SVI does not cover
  American Samoa, Guam, Northern Mariana Islands or U.S. Virgin Islands).
* Connecticut: SVI {s['release']} is on **{cov['connecticut_vintage']}** ({', '.join(cov['connecticut_source_fips'])});
  {len(cov['connecticut_matched'])} match the canonical planning regions. No legacy CT county rows exist; legacy-CT
  sources cannot be joined to SVI without a crosswalk (none applied).
* Never pool `geo_vulnerability` and `geo_vulnerability_puerto_rico` percentiles into one ranking.
* Ecological join only: an SVI value describes a county, not any person in a person-level dataset.

## Limitations and caveats

* SVI is a relative ranking (percentiles within the release), not an absolute measure; ranks are not comparable
  across SVI releases or between the U.S. and Puerto Rico databases.
* SVI {s['release']} is built on ACS {s['acs_period']}; the ACS context table in this project is a later period, so SVI
  inputs and ACS columns differ in vintage.
* Inputs are ACS survey estimates with MOEs (mp_ columns); small counties have wide MOEs, so small rank
  differences are not meaningful.
* Theme 3 (racial & ethnic minority status) is a single-variable theme; it is included because it is part of the
  published SVI construct, not as a proxy for health status.
* Social vulnerability is not disease prevalence, not a burden measure, and not a proxy for either.

## Processed outputs

| Table | Rows |
|---|---|
| geo_vulnerability | {s['rows']['us']} |
| geo_vulnerability_puerto_rico | {s['rows']['pr']} |

`build_stats.json` in this directory holds every number above in machine-readable form.

## Reproduce

`uv run python -m measure_it.ingestion.svi`
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


def run() -> dict:
    return build()


if __name__ == "__main__":
    st = run()
    print(json.dumps({"release": st["release"], "rows": st["rows"], "probes": st["probes"],
                      "coverage_us": {k: st["coverage_us"][k] for k in ("matched", "source_unmatched",
                                                                         "canonical_without_source_row",
                                                                         "connecticut_vintage")}}, indent=2))
