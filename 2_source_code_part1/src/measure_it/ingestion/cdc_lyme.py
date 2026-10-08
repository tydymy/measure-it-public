"""CDC Lyme disease surveillance: county x year reported cases (2001-2023) + incidence per 100,000.

Primary file: CDC "Lyme disease case map" data set, county of residence, 2001-2023
  https://www.cdc.gov/lyme/media/files/2025/02/LD_Case_Counts_by_County_2023_updated.csv
Cross-check: data.cdc.gov "Lyme disease public use aggregated data with geography" (84rx-ksgd
1992-2007, qtbi-xd4i 2008-2021, x5j9-wybp 2022-2023; privacy-suppressed cells). The run reconciles
yearly totals of the two and writes the comparison into DATA_AUDIT.md.

Denominator for incidence (documented, not CDC's): U.S. Census Bureau Population Estimates Program,
July 1 resident population by county:
  2001-2009  2000-2010 intercensal (co-est00int-tot.csv)
  2010-2019  2010-2020 intercensal (cc-est2020int-agesex-all.csv, YEAR codes 2..11 = July 1 2010..2019)
  2020-2023  Vintage 2025 (co-est2025-alldata.csv; Connecticut = planning regions)
  2020-2021  Vintage 2021 (co-est2021-alldata.csv) for the 8 legacy Connecticut counties only
Renamed-only FIPS are carried to the current code for 2001-2009 (02270->02158 Kusilvak,
46113->46102 Oglala Lakota). No other county is re-mapped; missing denominators stay missing.

Connecticut: the CDC file reports legacy counties (09001-09015) for 2001-2022 and planning regions
(09110-09190) for 2023; the other set is blank or 0 in each year. Those placeholder cells are dropped
(counted in the audit), never re-allocated between vintages.

Guardrails: condition lyme_disease, burden_evidence_level 'A' (surveillance case counts), with heavy
under-ascertainment and a 2022 surveillance change (lab-only reporting in high-incidence states) that
breaks comparability with earlier years. No PTLDS burden is produced (PTLDS = level D).

Output (data/processed): geo_condition_burden__cdc_lyme

Run: uv run python -m measure_it.ingestion.cdc_lyme
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from ..config import raw_dir
from ..download import download_file, load_manifest
from ..geography.crosswalk import STATE_FIPS
from ..http import get_json
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

SOURCE_ID = "cdc_lyme"
SOURCE_NAME = "CDC Lyme disease surveillance, reported cases by county of residence"
PRODUCER = "measure_it.ingestion.cdc_lyme"
LANDING_URL = "https://www.cdc.gov/lyme/data-research/facts-stats/lyme-disease-case-map.html"
COUNTY_URL = "https://www.cdc.gov/lyme/media/files/2025/02/LD_Case_Counts_by_County_2023_updated.csv"
# www.cdc.gov (Akamai) answers the project User-Agent with HTTP 403; the plain requests UA is accepted.
CDC_WWW_HEADERS = {"User-Agent": f"python-requests/{requests.__version__}"}
PUBLIC_USE = {"84rx-ksgd": "1992_2007", "qtbi-xd4i": "2008_2021", "x5j9-wybp": "2022_2023"}
PEP = {
    "co-est00int-tot.csv":
        "https://www2.census.gov/programs-surveys/popest/datasets/2000-2010/intercensal/county/co-est00int-tot.csv",
    "cc-est2020int-agesex-all.csv":
        "https://www2.census.gov/programs-surveys/popest/datasets/2010-2020/intercensal/county/asrh/cc-est2020int-agesex-all.csv",
    "co-est2025-alldata.csv":
        "https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/counties/totals/co-est2025-alldata.csv",
    "co-est2021-alldata.csv":
        "https://www2.census.gov/programs-surveys/popest/datasets/2020-2021/counties/totals/co-est2021-alldata.csv",
}
FIPS_RENAMES_2000S = {"02270": "02158", "46113": "46102"}
CT_LEGACY = {f"09{c:03d}" for c in range(1, 16, 2)}
CT_PLANNING = {f"09{c}" for c in range(110, 200, 10)}


# --------------------------------------------------------------------------- pure helpers
def case_definition_period(year: int) -> str:
    if year <= 2007:
        return "2001-2007: confirmed cases (no probable category before 2008)"
    if year <= 2021:
        return "2008-2021: confirmed + probable"
    return "2022+: 2022 surveillance case definition; high-incidence states report on laboratory evidence alone"


def normalize_status(s: str) -> str:
    s = str(s).strip()
    if s.lower().startswith("high"):
        return "High Incidence"
    if s.lower().startswith("low"):
        return "Low Incidence"
    return s


def parse_county_file(raw: pd.DataFrame) -> pd.DataFrame:
    """Wide CDC county file -> long (fips, year, cases). Blank cells stay NaN."""
    cols = [c for c in raw.columns if c.lower().startswith("cases")]
    d = raw.copy()
    d["fips"] = d["stcode"].astype(str).str.zfill(2) + d["ctycode"].astype(str).str.zfill(3)
    long = d.melt(id_vars=["fips", "Ctyname", "stname", "ststatus"], value_vars=cols,
                  var_name="col", value_name="cases")
    long["year"] = long["col"].str[-4:].astype(int)
    long["cases"] = pd.to_numeric(long["cases"], errors="coerce")
    long["state_incidence_category"] = long["ststatus"].map(normalize_status)
    return long.drop(columns=["col", "ststatus"]).rename(columns={"Ctyname": "geo_name", "stname": "state_name"})


def drop_placeholders(long: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Drop blank cells and the Connecticut vintage that is not in use in a given year.

    In a year where any CT planning region has a value, legacy-county cells are placeholders
    (the file carries 0 there), and vice versa.
    """
    n0 = len(long)
    blank = long["cases"].isna()
    ct_new_years = set(long.loc[long["fips"].isin(CT_PLANNING) & long["cases"].notna(), "year"])
    ct_old_years = set(long.loc[long["fips"].isin(CT_LEGACY) & long["cases"].notna() & (long["cases"] > 0), "year"])
    ct_legacy_placeholder = long["fips"].isin(CT_LEGACY) & long["year"].isin(ct_new_years - ct_old_years)
    keep = ~blank & ~ct_legacy_placeholder
    out = long[keep].copy()
    return out, {"cells": n0, "blank_cells": int(blank.sum()),
                 "blank_cells_ct_planning_regions": int((blank & long["fips"].isin(CT_PLANNING)).sum()),
                 "ct_legacy_placeholder_cells": int((ct_legacy_placeholder & ~blank).sum()),
                 "ct_planning_region_years": sorted(ct_new_years), "kept_cells": int(keep.sum())}


def incidence_per_100k(cases: pd.Series, population: pd.Series) -> pd.Series:
    pop = population.astype("float64")
    return (cases.astype("float64") / pop * 1e5).where(pop > 0)


def _std_fips(df: pd.DataFrame) -> pd.Series:
    return (pd.to_numeric(df["STATE"]).astype(int).astype(str).str.zfill(2)
            + pd.to_numeric(df["COUNTY"]).astype(int).astype(str).str.zfill(3))


def population_table(paths: dict[str, Path]) -> pd.DataFrame:
    """County x year July-1 population with its source file (see module docstring)."""
    frames = []
    a = pd.read_csv(paths["co-est00int-tot.csv"], dtype=str, encoding="latin-1")
    a = a[pd.to_numeric(a["SUMLEV"]) == 50]
    a["fips"] = _std_fips(a).replace(FIPS_RENAMES_2000S)
    for y in range(2001, 2010):
        frames.append(pd.DataFrame({"fips": a["fips"], "year": y, "population": pd.to_numeric(a[f"POPESTIMATE{y}"]),
                                    "denominator_source": "Census PEP 2000-2010 intercensal (co-est00int-tot), July 1"}))
    b = pd.read_csv(paths["cc-est2020int-agesex-all.csv"], dtype=str, encoding="latin-1",
                    usecols=["SUMLEV", "STATE", "COUNTY", "YEAR", "POPESTIMATE"])
    b = b[(pd.to_numeric(b["SUMLEV"]) == 50) & pd.to_numeric(b["YEAR"]).between(2, 11)]
    frames.append(pd.DataFrame({"fips": _std_fips(b), "year": pd.to_numeric(b["YEAR"]) + 2008,
                                "population": pd.to_numeric(b["POPESTIMATE"]),
                                "denominator_source": "Census PEP 2010-2020 intercensal (cc-est2020int), July 1"}))
    c = pd.read_csv(paths["co-est2025-alldata.csv"], dtype=str, encoding="latin-1")
    c = c[pd.to_numeric(c["SUMLEV"]) == 50]
    c["fips"] = _std_fips(c)
    for y in range(2020, 2024):
        frames.append(pd.DataFrame({"fips": c["fips"], "year": y, "population": pd.to_numeric(c[f"POPESTIMATE{y}"]),
                                    "denominator_source": "Census PEP Vintage 2025 (co-est2025-alldata), July 1"}))
    d = pd.read_csv(paths["co-est2021-alldata.csv"], dtype=str, encoding="latin-1")
    d = d[pd.to_numeric(d["SUMLEV"]) == 50]
    d["fips"] = _std_fips(d)
    d = d[d["fips"].isin(CT_LEGACY)]
    for y in (2020, 2021):
        frames.append(pd.DataFrame({"fips": d["fips"], "year": y, "population": pd.to_numeric(d[f"POPESTIMATE{y}"]),
                                    "denominator_source": "Census PEP Vintage 2021 (co-est2021-alldata), July 1; "
                                                          "last vintage on legacy CT counties"}))
    pop = pd.concat(frames, ignore_index=True)
    dup = pop.duplicated(["fips", "year"], keep=False)
    if dup.any():
        raise ValueError(f"overlapping denominators for {pop[dup].head().to_dict('records')}")
    return pop


def suspected_state_nonreporting(county_long: pd.DataFrame, min_adjacent: int = 20) -> pd.DataFrame:
    """State-years whose county-file total is 0 while every adjacent year present has >= min_adjacent cases.

    Such a zero is almost certainly a state that did not report that year (CDC: "2019 and 2020 data
    from some jurisdictions are incomplete"), not an absence of disease. Returns state_abbr, year, flag.
    """
    tot = county_long.groupby(["state_abbr", "year"])["cases"].sum()
    rows = []
    for (st, yr), v in tot.items():
        adj = [tot.get((st, yr + k)) for k in (-1, 1)]
        adj = [a for a in adj if a is not None and pd.notna(a)]
        rows.append({"state_abbr": st, "year": yr,
                     "suspected_state_nonreporting": bool(v == 0 and adj and min(adj) >= min_adjacent)})
    return pd.DataFrame(rows)


def county_file_coverage(public_use: pd.DataFrame, county_long: pd.DataFrame) -> pd.DataFrame:
    """State x year: how many of the state's reported cases the county file leaves out.

    public_use: the data.cdc.gov aggregated rows (Year, State, FIPS, Frequency).
    county_long: county-file rows with state_abbr, year, cases.
    Cases whose county of residence is unknown are generally absent from the county file (see
    reconcile), so county counts and rates in a state-year with a large gap are understated.
    state_cases_missing_from_county_file_pct = max(public-use state total - county-file state total, 0)
    / public-use state total. It is a lower bound: public-use rows whose State is itself suppressed
    or unknown are not attributed to any state. NaN when the public-use files have no row for the
    state-year.
    """
    pu = public_use.assign(f=pd.to_numeric(public_use["Frequency"]), year=pd.to_numeric(public_use["Year"]).astype(int),
                           unk=public_use["FIPS"].astype(str).str.strip().str.lower().eq("unknown"))
    tot = pu.groupby(["State", "year"])["f"].sum().rename("pu_state_cases")
    unk = pu[pu["unk"]].groupby(["State", "year"])["f"].sum().rename("pu_unknown_county_cases")
    out = pd.concat([tot, unk], axis=1).fillna({"pu_unknown_county_cases": 0}).reset_index()
    out = out.rename(columns={"State": "state_abbr"})
    cf = county_long.groupby(["state_abbr", "year"])["cases"].sum().rename("county_file_state_cases").reset_index()
    out = cf.merge(out, on=["state_abbr", "year"], how="left")
    gap = (out["pu_state_cases"] - out["county_file_state_cases"]).clip(lower=0)
    out["state_cases_missing_from_county_file_pct"] = (100 * gap / out["pu_state_cases"]).where(
        out["pu_state_cases"] > 0)
    return out


def reconcile(public_use: pd.DataFrame, county_long: pd.DataFrame) -> pd.DataFrame:
    """Yearly totals: CDC county file vs public-use aggregated datasets (by FIPS category)."""
    pu = public_use.copy()
    pu["Frequency"] = pd.to_numeric(pu["Frequency"])
    pu["year"] = pd.to_numeric(pu["Year"]).astype(int)
    fips = pu["FIPS"].astype(str).str.strip()
    pu["fips_kind"] = np.where(fips.str.fullmatch(r"\d{5}"), "county",
                               np.where(fips.str.lower().eq("unknown"), "unknown", "suppressed_or_other"))
    t = pu.pivot_table(index="year", columns="fips_kind", values="Frequency", aggfunc="sum").fillna(0)
    cs = pu.pivot_table(index="year", columns="Case_status", values="Frequency", aggfunc="sum").fillna(0)
    t = t.join(cs.add_prefix("case_status_"))
    t["public_use_total"] = t[[c for c in ("county", "unknown", "suppressed_or_other") if c in t]].sum(axis=1)
    t["county_file_total"] = county_long.groupby("year")["cases"].sum()
    t["public_use_minus_unknown"] = t["public_use_total"] - t.get("unknown", 0)
    t["county_file_minus_pu_known"] = t["county_file_total"] - t["public_use_minus_unknown"]
    return t.reset_index()


# --------------------------------------------------------------------------- audit
def _pct(n: int, d: int) -> str:
    return f"{n:,} ({100 * n / d:.2f}%)" if d else f"{n:,}"


def write_audit(s: dict) -> Path:
    t, rec, py = s["table"], s["reconcile"], s["per_year"]
    rec_md = "\n".join(
        f"| {int(r.year)} | {int(r.public_use_total):,} | {int(r.case_status_Confirmed):,} | {int(r.case_status_Probable):,} | "
        f"{int(r.unknown):,} | {int(r.suppressed_or_other):,} | {int(r.county_file_total) if pd.notna(r.county_file_total) else 'n/a'} | "
        f"{r.county_file_minus_pu_known:+.0f} |" for r in rec[rec["year"] >= 2001].itertuples())
    py_md = "\n".join(
        f"| {int(r.year)} | {int(r.counties):,} | {int(r.cases):,} | {_pct(int(r.zero), int(r.counties))} | "
        f"{_pct(int(r.no_denominator), int(r.counties))} | {r.crude_rate:.1f} | {r.definition} |" for r in py.itertuples())
    files = "\n".join(f"| {k} | {v['url'][:150]} | {v['bytes']:,} | `{v['sha256'][:16]}...` | {v['retrieved_at']} |"
                      for k, v in s["manifest"].items())
    dr = s["drops"]
    nodenom = s["no_denominator"]
    nd_md = "\n".join(f"| {r.fips} | {r.geo_name}, {r.state_name} | {r.years} |" for r in nodenom.itertuples())
    text = f"""# DATA AUDIT — CDC Lyme disease county case counts (2001-2023)

| Field | Value |
|---|---|
| source_id | cdc_lyme |
| Source (dataset/API name, exact files/endpoints) | CDC Lyme disease case-map data set `LD_Case_Counts_by_County_2023_updated.csv` (linked from {LANDING_URL}); cross-check: data.cdc.gov public-use aggregated datasets {', '.join(PUBLIC_USE)}; denominators: Census PEP county files ({', '.join(PEP)}) |
| Publishing organization | CDC NCEZID Division of Vector-Borne Diseases (NNDSS data reported by states); Census Bureau (denominators) |
| Retrieval date (UTC) | {s['retrieved_at']} |
| Source version / release | county file Last-Modified {s['last_modified']}; public-use datasets rows updated {s['pu_updated']} |
| Source update date / cadence | annual, after final verification of state surveillance data |
| License / access conditions | U.S. Government work (public domain); open download. www.cdc.gov returns HTTP 403 to the project User-Agent; the file is fetched with the default python-requests User-Agent |
| Unit of observation | county of residence x year: reported Lyme disease cases (confirmed; + probable from 2008) |
| Sample size (actual, as ingested) | raw {s['raw_rows']:,} county rows x 23 years = {dr['cells']:,} cells; {dr['kept_cells']:,} county-years kept; {s['n_counties']:,} distinct county FIPS; {int(t.loc[t.measure_id == 'lyme_reported_cases', 'value'].sum()):,} cases total 2001-2023 |
| Geography (resolution, vintage) | county of residence. Connecticut: legacy counties 09001-09015 for 2001-2022, planning regions 09110-09190 for {', '.join(map(str, dr['ct_planning_region_years']))}. Matches to `geographies`: {s['match_2024']:,} FIPS in 2024 counties, {s['match_legacy']} legacy CT, {s['unmatched']} unmatched ({s['unmatched_list']}) |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — aggregated surveillance counts |

## Burden labelling
`condition_id = lyme_disease`, `burden_evidence_level = A` (surveillance case counts of the condition itself),
`evidence_type = surveillance_case_count`, `source_geographic_resolution = county`. Two measures per county-year:
`lyme_reported_cases` (metric_type case_count) and `lyme_incidence_per_100k` (metric_type incidence_per_100k,
`primary_burden_measure = True`). Both carry `numerator` (cases) and `denominator` (population) so analyses can use
rates or absolute counts.

**Level A here means "direct measure of reported disease", not "complete".** Reported cases undercount true infections
substantially (CDC's own data set description: "Not every case of Lyme disease is reported to CDC"), reporting practice
varies between states and years, and cases are assigned to county of residence, not place of exposure.

**2022 break:** "In 2022, states with a high incidence of Lyme disease started reporting cases based on laboratory
evidence alone without requirement for a clinical investigation, precluding comparison with historical data"
(data.cdc.gov dataset description). Rows carry `case_definition_period` and `comparable_with_pre_2022`; do not trend
across 2021/2022. `state_incidence_category` (High/Low Incidence) is CDC's state classification from the file.

**PTLDS (post-treatment Lyme disease syndrome): burden level D.** No public surveillance or survey source estimates
PTLDS by geography; this module writes no PTLDS rows. Lyme incidence is not a PTLDS measure; if a later step uses it
as a PTLDS proxy it must be labelled level C proxy explicitly.

## Denominator (documented choice, not CDC's)
Incidence per 100,000 = cases / July-1 resident population x 100,000. Population: Census PEP 2000-2010 intercensal
(2001-2009), 2010-2020 intercensal (2010-2019; YEAR codes 2-11 are July 1 2010-2019, code 1 is the April 2010 base,
code 12 April 2020), Vintage 2025 (2020-2023), and Vintage 2021 for legacy Connecticut counties in 2020-2021.
Renamed FIPS carried forward for 2001-2009: {FIPS_RENAMES_2000S}. County-years with no matching denominator keep the
count and get a blank rate:

| fips | county | years without denominator |
|---|---|---|
{nd_md}

## Per-year summary (kept county-years)
| year | counties | cases | counties with 0 cases | no denominator | crude rate per 100k (sum cases / sum pop, counties with pop) | case definition |
|---|---|---|---|---|---|---|
{py_md}

## Reconciliation with the data.cdc.gov public-use aggregated datasets
Public-use rows are cases by year x state x county FIPS x case status x sex x age group; FIPS is `Unknown`
when county of residence was not known and `Suppressed` when CDC's privacy protection (the Lee et al. 2021
method cited in the dataset description) removed the county.

| year | public-use total | confirmed | probable | FIPS unknown | FIPS suppressed/other | county file total | county file - (public-use - unknown) |
|---|---|---|---|---|---|---|---|
{rec_md}

Reading: in every year the county-file total equals the public-use total (confirmed + probable) minus cases with
unknown county to within {s['second_gap']:.0f} cases, except {s['max_gap_year']} ({s['max_gap']:+.0f}). So the county file counts
confirmed + probable cases with a known county of residence; cases with unknown county are generally not in it
(not always: in some state-years, e.g. DC, the county file matches the full state total). The state-year table
below measures the gap directly.
County-years where the public-use county sum exceeds the county file or is absent from it: {s['cell_diffs']}.

**Cases missing from the county file are concentrated in a few state-years**, where every county count and rate
is understated. Each row carries `state_cases_missing_from_county_file_pct` = max(public-use state total - county-file
state total, 0) / public-use state total x 100. It is a lower bound: public-use rows whose State is itself suppressed or
unknown are attributed to no state. State-years with a gap above 10% and at least 20 cases:

| state | year | public-use state cases | of which FIPS unknown | county-file state cases | % missing |
|---|---|---|---|---|---|
{s['unknown_share_md']}

{s['no_pu_md']}; for those the column is blank.

## Files / endpoints retrieved
| file | url | bytes | sha256 | retrieved_at |
|---|---|---|---|---|
{files}

## Key variables
| column | raw | meaning |
|---|---|---|
| geo_id / county_fips | stcode + ctycode | 5-char FIPS of county of residence |
| year | CasesYYYY / casesYYYY column | report year |
| value (lyme_reported_cases) / numerator | CasesYYYY | reported cases |
| value (lyme_incidence_per_100k) | derived | cases / population x 1e5 |
| denominator, denominator_source | Census PEP | July-1 population and the file it came from |
| state_incidence_category | ststatus | CDC High/Low Incidence state class (raw text truncated as "High Incidenc") |
| state_cases_missing_from_county_file_pct | derived: public-use vs county file | % of the state's public-use cases that year absent from the county file (lower bound) |
| case_definition_period, comparable_with_pre_2022 | derived from year | surveillance era |

## Missingness
* Raw blank cells: {dr['blank_cells']:,} of {dr['cells']:,} ({dr['blank_cells_ct_planning_regions']:,} are CT planning regions before 2023;
  other FIPS with blank cells: {s['other_blank']}).
* CT legacy-county placeholder cells dropped (0 in years reported on planning regions): {dr['ct_legacy_placeholder_cells']}.
* Incidence without denominator: {_pct(s['incidence_na'], s['incidence_rows'])} incidence rows.
* No suppression flags in the county file: every kept cell is a numeric count. A 0 is what the file reports; the file
  does not distinguish "no case reported" from "jurisdiction did not report". `suspected_state_nonreporting` = True
  where the state's county-file total is 0 in a year while every adjacent year present has >= 20 cases: {s['nonreport_md']}.
  Those zeros (and 0 incidence rates) are not evidence of absent disease.
* 2024 counties (50 states + DC) absent from the file: {s['absent_counties']}. Territories are absent.
  The file keeps Valdez-Cordova (02261, split in 2019) through 2023 alongside Chugach (02063, 2023 only), and keeps
  Bedford city (51515, merged into Bedford County 51019 in 2013) as 0 through 2023; those later zeros are rows for an
  entity that no longer exists (`in_canonical_geographies=False`), not observed absence of disease.
* `state_cases_missing_from_county_file_pct` missing: {_pct(s['unknown_pct_na'], len(t))} rows (state-years with no
  public-use row for the state: no cases, or the State field suppressed/unknown in the public-use files).

## Linkage strategy
`geo_id` joins `geographies.geo_id` (2024 vintage) for {s['match_2024']:,} FIPS. Legacy CT counties (2001-2022) match the
`ct_legacy=True` rows only; they must not be mapped onto planning regions. `02261` (Valdez-Cordova, dissolved 2019) and `51515`
(Bedford city, merged 2013) have no 2024 geography (`in_canonical_geographies=False`). Joins to other layers are ecological.
Row level ({len(t):,} rows): {s['rows_2024']:,} match a 2024 county in `geographies`; {s['rows_legacy']:,} match legacy CT counties
(`ct_legacy=True`, years 2001-2022); {s['rows_unmatched']:,} match nothing. Connecticut rows: {s['rows_ct']:,}, of which
{s['rows_ct_2024']:,} are planning regions (2023) and {s['rows_legacy']:,} legacy counties.

## Limitations and caveats
* Under-ascertainment and state-to-state differences in surveillance intensity; year-to-year changes can be artefacts.
* 2022 surveillance change makes 2022-2023 counts not comparable with 2001-2021 (national count {s['cases_2021']:,} in 2021 vs {s['cases_2022']:,} in 2022).
* "Due to the coronavirus disease 2019 (COVID-19) pandemic, 2019 and 2020 data from some jurisdictions are incomplete"
  (CDC Lyme surveillance data page); 2020 counts are low, and see `suspected_state_nonreporting`.
* County of residence is not county of exposure.
* Cases with unknown county of residence are mostly not in the county file: up to {s['unknown_max']} of a state's
  cases in a state-year with at least 100 cases (see the table under the reconciliation); use
  `state_cases_missing_from_county_file_pct`.
* Small denominators make rates unstable; use `denominator` to filter (scoring config min_population = 10,000).
* Denominator vintages differ across eras (intercensal vs postcensal); rates are this project's computation, not CDC-published county rates.

## Processed outputs
| table | rows |
|---|---|
| geo_condition_burden__cdc_lyme | {len(t):,} ({len(t) // 2:,} county-years x 2 measures) |

## Reproduce
`uv run python -m measure_it.ingestion.cdc_lyme`
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


# --------------------------------------------------------------------------- main
def run() -> dict:
    county_path = download_file(COUNTY_URL, SOURCE_ID, headers=CDC_WWW_HEADERS, min_bytes=100_000)
    (raw_dir(SOURCE_ID) / "denominators").mkdir(exist_ok=True)
    (raw_dir(SOURCE_ID) / "public_use").mkdir(exist_ok=True)
    pep_paths = {k: download_file(u, SOURCE_ID, f"denominators/{k}") for k, u in PEP.items()}
    pu_frames, pu_meta = [], {}
    for dsid, span in PUBLIC_USE.items():
        p = download_file(f"https://data.cdc.gov/api/views/{dsid}/rows.csv?accessType=DOWNLOAD", SOURCE_ID,
                          f"public_use/lyme_aggregated_geo_{span}_{dsid}.csv")
        pu_frames.append(pd.read_csv(p, dtype=str))
        m = get_json(f"https://data.cdc.gov/api/views/{dsid}.json")
        pu_meta[dsid] = pd.Timestamp(m["rowsUpdatedAt"], unit="s").isoformat() if m.get("rowsUpdatedAt") else "UNKNOWN"

    raw = pd.read_csv(county_path, dtype=str, encoding="latin-1")
    long, drops = drop_placeholders(parse_county_file(raw))
    pop = population_table(pep_paths)
    long = long.merge(pop, on=["fips", "year"], how="left")
    long["incidence"] = incidence_per_100k(long["cases"], long["population"])
    pu_all = pd.concat(pu_frames, ignore_index=True)
    long["state_abbr"] = long["fips"].str[:2].map(STATE_FIPS)
    coverage = county_file_coverage(pu_all, long)
    nonreport = suspected_state_nonreporting(long)
    n_before = len(long)
    long = long.merge(coverage[["state_abbr", "year", "state_cases_missing_from_county_file_pct"]],
                      on=["state_abbr", "year"], how="left").merge(nonreport, on=["state_abbr", "year"], how="left")
    if len(long) != n_before or long["suspected_state_nonreporting"].isna().any():
        raise ValueError("state-year flag merge changed the row set")

    geos = read_table("geographies")
    gc = geos[geos["geo_level"] == "county"][["geo_id", "vintage"]].rename(columns={"geo_id": "fips"})
    long = long.merge(gc, on="fips", how="left")
    long["in_canonical_geographies"] = long["vintage"].eq("2024")
    long["geo_vintage"] = long["vintage"].fillna("not_in_geographies")

    base = pd.DataFrame({
        "geo_id": long["fips"], "geo_level": "county", "geo_name": long["geo_name"],
        "state_fips": long["fips"].str[:2], "state_abbr": long["state_abbr"],
        "state_name": long["state_name"], "county_fips": long["fips"],
        "geo_vintage": long["geo_vintage"], "in_canonical_geographies": long["in_canonical_geographies"],
        "condition_id": "lyme_disease", "year": long["year"].astype("Int64"),
        "period_start": long["year"].astype(str) + "-01-01", "period_end": long["year"].astype(str) + "-12-31",
        "period_label": long["year"].astype(str),
        "numerator": long["cases"].astype("Int64"), "denominator": long["population"].astype("Int64"),
        "denominator_source": long["denominator_source"].fillna("UNKNOWN / NOT AVAILABLE"),
        "ci_low": np.nan, "ci_high": np.nan, "subgroup_type": "All", "subgroup": "All", "suppressed": False,
        "state_incidence_category": long["state_incidence_category"],
        "state_cases_missing_from_county_file_pct": long["state_cases_missing_from_county_file_pct"],
        "suspected_state_nonreporting": long["suspected_state_nonreporting"].astype(bool),
        "case_definition_period": long["year"].map(case_definition_period),
        "comparable_with_pre_2022": long["year"] < 2022, "burden_evidence_level": "A",
    })
    counts = base.assign(measure_id="lyme_reported_cases", measure_label="Reported Lyme disease cases (county of residence)",
                         metric_type="case_count", value=long["cases"].astype("float64"), value_unit="cases",
                         primary_burden_measure=False)
    rates = base.assign(measure_id="lyme_incidence_per_100k",
                        measure_label="Reported Lyme disease incidence per 100,000 population (Census PEP denominator)",
                        metric_type="incidence_per_100k", value=long["incidence"], value_unit="per 100,000 population",
                        primary_burden_measure=True)
    table = pd.concat([counts, rates], ignore_index=True)

    manifest = load_manifest(SOURCE_ID)["files"]
    cm = manifest[county_path.name]
    version = f"LD_Case_Counts_by_County_2023_updated.csv (Last-Modified {cm.get('last_modified')}); cases 2001-2023"
    out = add_provenance(
        table, data_layer="geographic", source_name=SOURCE_NAME, source_version=version, retrieved_at=cm["retrieved_at"],
        evidence_type="surveillance_case_count",
        source_record_id=lambda d: "LD_Case_Counts_by_County:" + d["geo_id"] + ":" + d["year"].astype(str) + ":" + d["measure_id"],
        source_geographic_resolution="county", evidence_level="A",
        provenance_notes="Reported surveillance cases (under-ascertained; 2022 case-definition break). Incidence uses "
                         "Census PEP July-1 population (see denominator_source). No PTLDS burden (level D).")
    write_table(out, "geo_condition_burden__cdc_lyme", producer=PRODUCER,
                description="CDC Lyme reported cases + incidence per 100k by county of residence, 2001-2023, level A")

    # ---- measured stats
    rec = reconcile(pd.concat(pu_frames, ignore_index=True), long)
    r01 = rec[rec["year"].between(2001, 2023)]
    gap = r01.loc[r01["county_file_minus_pu_known"].abs().idxmax()]
    pu = pd.concat(pu_frames, ignore_index=True)
    pu = pu[pu["FIPS"].astype(str).str.fullmatch(r"\d{5}")]
    puc = pu.assign(Frequency=pd.to_numeric(pu["Frequency"]), year=pd.to_numeric(pu["Year"]).astype(int)) \
        .groupby(["FIPS", "year"])["Frequency"].sum().rename("pu").reset_index().rename(columns={"FIPS": "fips"})
    cmp_ = long[["fips", "year", "cases"]].merge(puc, on=["fips", "year"], how="right")
    cmp_ = cmp_[cmp_["year"].between(2001, 2023)]
    diffs = cmp_[(cmp_["pu"] > cmp_["cases"]) | cmp_["cases"].isna()]
    cell_diffs = "; ".join(f"{r.fips} {r.year}: public-use {int(r.pu)} vs county file "
                           f"{'absent' if pd.isna(r.cases) else int(r.cases)}" for r in diffs.itertuples()) or "none"
    per_year = long.groupby("year").apply(lambda g: pd.Series({
        "counties": len(g), "cases": g["cases"].sum(), "zero": int((g["cases"] == 0).sum()),
        "no_denominator": int(g["population"].isna().sum()),
        "crude_rate": g.loc[g["population"].notna(), "cases"].sum() / g["population"].sum() * 1e5,
    }), include_groups=False).reset_index()
    per_year["definition"] = per_year["year"].map(lambda y: case_definition_period(int(y)).split(":")[0])
    nd = long[long["population"].isna()].groupby(["fips", "geo_name", "state_name"])["year"] \
        .apply(lambda s: ", ".join(map(str, sorted(s)))).rename("years").reset_index()
    fips_all = pd.Series(long["fips"].unique())
    rl = parse_county_file(raw)
    ob = rl[rl["cases"].isna() & ~rl["fips"].isin(CT_PLANNING)].groupby(["fips", "geo_name"])["year"].agg(["min", "max", "size"])
    s_other_blank = "; ".join(f"{f} {n} ({r['size']} cells, {r['min']}-{r['max']})" for (f, n), r in ob.iterrows()) or "none"
    g24 = geos[(geos["geo_level"] == "county") & ~geos["ct_legacy"] & geos["state_abbr"].isin(
        [a for f, a in STATE_FIPS.items() if int(f) <= 56])]
    ab = g24[~g24["geo_id"].isin(fips_all)]
    s_absent = "; ".join(f"{r.geo_id} {r.name}, {r.state_abbr}" for r in ab.itertuples()) or "none"
    g_all = geos[geos["geo_level"] == "county"]
    unmatched = long.loc[~long["fips"].isin(g_all["geo_id"]), ["fips", "geo_name"]].drop_duplicates()
    cov = coverage.assign(gap=coverage["pu_state_cases"] - coverage["county_file_state_cases"])
    hi = cov[(cov["state_cases_missing_from_county_file_pct"] > 10) & (cov["gap"] >= 20)].sort_values(
        ["state_abbr", "year"])
    unknown_share_md = "\n".join(
        f"| {r.state_abbr} | {r.year} | {int(r.pu_state_cases):,} | {int(r.pu_unknown_county_cases):,} | "
        f"{int(r.county_file_state_cases):,} | {r.state_cases_missing_from_county_file_pct:.1f} |"
        for r in hi.itertuples()) or "| none | | | | | |"
    big = cov[cov["pu_state_cases"] >= 100]
    top = big.loc[big["state_cases_missing_from_county_file_pct"].idxmax()]
    unknown_max = f"{top['state_cases_missing_from_county_file_pct']:.0f}% ({top['state_abbr']} {int(top['year'])})"
    no_pu = cov[cov["pu_state_cases"].isna() & (cov["county_file_state_cases"] > 0)]
    no_pu_md = (f"{len(no_pu)} state-years with county-file cases have no public-use row for the state (State suppressed "
                f"or unknown in the public-use files), largest: " + ", ".join(
                    f"{r.state_abbr} {r.year} ({int(r.county_file_state_cases)} cases)"
                    for r in no_pu.nlargest(3, "county_file_state_cases").itertuples()))
    stats = {
        "table": table, "reconcile": rec, "per_year": per_year, "drops": drops, "no_denominator": nd,
        "raw_rows": len(raw), "n_counties": len(fips_all), "retrieved_at": cm["retrieved_at"],
        "last_modified": cm.get("last_modified"), "pu_updated": ", ".join(f"{k} {v}" for k, v in pu_meta.items()),
        "match_2024": int(fips_all.isin(g_all[~g_all["ct_legacy"]]["geo_id"]).sum()),
        "match_legacy": int(fips_all.isin(g_all[g_all["ct_legacy"]]["geo_id"]).sum()),
        "unmatched": len(unmatched), "unmatched_list": ", ".join(f"{r.fips} {r.geo_name}" for r in unmatched.itertuples()),
        "incidence_na": int(rates["value"].isna().sum()), "incidence_rows": len(rates),
        "rows_2024": int(table["geo_id"].isin(g_all[~g_all["ct_legacy"]]["geo_id"]).sum()),
        "rows_legacy": int(table["geo_id"].isin(g_all[g_all["ct_legacy"]]["geo_id"]).sum()),
        "rows_unmatched": int((~table["geo_id"].isin(g_all["geo_id"])).sum()),
        "rows_ct": int((table["state_fips"] == "09").sum()),
        "rows_ct_2024": int(((table["state_fips"] == "09") & table["geo_id"].isin(CT_PLANNING)).sum()),
        "max_gap_year": int(gap["year"]), "max_gap": float(gap["county_file_minus_pu_known"]), "cell_diffs": cell_diffs,
        "second_gap": float(r01.loc[r01["year"] != gap["year"], "county_file_minus_pu_known"].abs().max()),
        "other_blank": s_other_blank, "absent_counties": s_absent,
        "unknown_share_md": unknown_share_md, "unknown_max": unknown_max, "no_pu_md": no_pu_md,
        "unknown_pct_na": int(table["state_cases_missing_from_county_file_pct"].isna().sum()),
        "nonreport_md": "; ".join(f"{r.state_abbr} {r.year} ({int(((long['state_abbr'] == r.state_abbr) & (long['year'] == r.year)).sum())} counties)"
                                  for r in nonreport[nonreport["suspected_state_nonreporting"]].itertuples()) or "none",
        "cases_2021": int(long.loc[long.year == 2021, "cases"].sum()),
        "cases_2022": int(long.loc[long.year == 2022, "cases"].sum()),
        "manifest": manifest,
    }
    write_audit(stats)
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "CDC Lyme disease surveillance: reported cases by county of residence, 2001-2023",
        "publisher": "CDC NCEZID Division of Vector-Borne Diseases (NNDSS)", "landing_url": LANDING_URL,
        "access_urls": [COUNTY_URL, *(f"https://data.cdc.gov/d/{k}" for k in PUBLIC_USE), *PEP.values()],
        "license": "U.S. Government work (public domain)",
        "access_conditions": "open download; www.cdc.gov returns HTTP 403 to non-default User-Agents, fetched with python-requests UA",
        "retrieved_at": cm["retrieved_at"], "source_version": version,
        "update_date": f"county file Last-Modified {cm.get('last_modified')}; annual",
        "data_layer": "geographic", "unit_of_observation": "county of residence x year (reported cases)",
        "sample_size": {"county_years": drops["kept_cells"], "counties": len(fips_all), "years": "2001-2023",
                        "rows": len(table), "total_cases_2001_2023": int(long["cases"].sum()),
                        "incidence_rows_without_denominator": stats["incidence_na"]},
        "geographic_resolution": "county (CT legacy counties 2001-2022, CT planning regions 2023)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (aggregate surveillance counts)",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": ["geo_condition_burden__cdc_lyme"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": PRODUCER,
        "burden_evidence_level": "A", "condition_ids": ["lyme_disease"],
        "ptlds_burden_evidence_level": "D (no burden rows written)",
        "denominator": "Census PEP July-1 county population (2000-2010 intercensal, 2010-2020 intercensal, V2025; V2021 for CT legacy 2020-2021)",
        "limitations": [
            "Reported surveillance cases: substantial under-ascertainment; varies by state and year",
            "2022 surveillance change (lab-only reporting in high-incidence states): not comparable with 2001-2021",
            "County of residence, not county of exposure",
            f"Zero state totals between reporting years are suspected non-reporting, not absence: {stats['nonreport_md']} "
            "(suspected_state_nonreporting)",
            f"Cases with unknown county of residence are mostly not in the county file (up to {unknown_max} of a state's "
            "cases); see state_cases_missing_from_county_file_pct",
            "Incidence denominators are this project's choice (Census PEP), not CDC-published county rates",
            "CT legacy counties 2001-2022 vs planning regions 2023; not re-mapped",
            "No PTLDS burden (level D)",
        ],
    })
    summary = {"rows": len(table), "county_years": drops["kept_cells"], "counties": len(fips_all),
               "incidence_na": stats["incidence_na"], "max_gap": [stats["max_gap_year"], stats["max_gap"]]}
    print(json.dumps(summary, indent=2, default=str))
    return summary


if __name__ == "__main__":
    run()
