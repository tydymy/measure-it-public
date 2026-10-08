"""American Community Survey 5-year geographic context (state, county, ZCTA).

Source route
------------
The brief asks for api.census.gov. As of 2026-09 every *data* request to
api.census.gov without a key is answered with HTTP 302 -> /data/missing_key.html
(response header ``X-DataWebAPI-KeyError: 1``; the page says "A valid key must
be included with each data API request"). Requesting a key is an account
sign-up, which this pipeline does not do on anyone's behalf. The metadata
endpoints (dataset list, variables, groups) still answer without a key, and
are used to detect the newest 5-year vintage.

The same estimates are published keylessly in the Census Bureau's
*table-based Summary File* (one pipe-delimited file per detailed table, every
summary level). We download the detailed (B) tables that underlie the requested
subject/profile tables and derive the percentages ourselves:

  requested (subject/profile)     derived here from detailed-table cells
  S1810_C03_001E  disability      B18101 "With a disability" cells / B18101_001
  S2701_C05_001E  uninsured       B27020 "No health insurance coverage" cells / B27020_001
  S1701_C03_001E  poverty         B17001_002 / B17001_001
  S2801_C02_014E  broadband       B28002_004 / B28002_001
  DP04_0058PE     no vehicle      B08201_002 / B08201_001
  S0101_C02_030E  age 65+         B09020_001 (population 65+, one published cell) / B01001_001
  S0101 age 18-64                 (B01001_001 - B09001_001 under 18 - B09020_001 65+) / B01001_001

The age numerators use single published cells instead of summing the 12 / 26
B01001 sex-by-age bands. The point estimates are identical (checked in code
against the B01001 band sums for every row), but the handbook approximation
for a sum of many cells overstated the 65+ share MOE about six-fold against the
published S0101_C02_030M (measured against SVI 2022's MP_AGE65, 2026-09-23).

Margins of error are 90% MOEs. Published MOEs are used for single cells
(total population, median age, median household income, one-cell numerators
and denominators). Sums and differences use the ACS handbook
root-sum-of-squares approximation (only the largest MOE among zero-estimate
cells), proportions use the handbook proportion formula (ratio formula when
the radicand is negative). These are approximations; the published
subject-table MOEs, which are computed with replicate weights, can differ.

Outputs
-------
geo_context__acs       one row per state (52) and county (3222), wide
geo_context_zcta__acs  one row per ZCTA, smaller variable set

Population density uses land area from the ``geographies`` / ``zcta_centroids``
tables (Census 2024 cartographic boundary + Gazetteer).
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from ..config import PROCESSED, RAW, raw_dir, utc_now_iso
from ..download import download_file, load_manifest
from ..http import get
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import read_table, write_table

SOURCE_ID = "census_acs"
SF_ROOT = "https://www2.census.gov/programs-surveys/acs/summary_file/{year}/table-based-SF"
API_DATASET = "https://api.census.gov/data/{year}/acs/acs5"
API_PROBE = "https://api.census.gov/data/{year}/acs/acs5?get=NAME,B01003_001E&for=state:01"
LANDING_URL = "https://www.census.gov/programs-surveys/acs/data/summary-file.html"
VINTAGE_CANDIDATES = (2025, 2024, 2023)

SUMLEVEL_PREFIX = {"state": "0400000US", "county": "0500000US", "zcta": "860Z200US"}

# ACS special values (census.gov "Notes on ACS estimate and annotation values").
ACS_SENTINELS = {
    -666666666: "estimate not computable (insufficient sample / open-ended median)",
    -999999999: "estimate or MOE not displayable (insufficient sample cases)",
    -888888888: "not applicable / not available",
    -555555555: "MOE not appropriate: estimate controlled to an independent population/housing estimate",
    -333333333: "MOE not computable: median in open-ended interval",
    -222222222: "MOE not computable: insufficient sample observations",
}
CONTROLLED_MOE = -555555555
# Open-ended median jam values in the Summary File.
INCOME_TOPCODE = 250001
INCOME_BOTTOMCODE = 2499

SQM_PER_SQMI = 2_589_988.110336


@dataclass
class Measure:
    """One output concept built from detailed-table cells.

    kind="estimate":   a single published cell (value + published MOE).
    kind="share":      numerator cells / denominator cell -> count, universe and percent.
    kind="difference": (cells[0] - sum(subtract)) / denominator -> count, universe and percent.

    cell_universes pins the table-shell Universe of a cell (used where the label alone,
    e.g. "Total:", does not identify the concept). check_cells is an independent cell
    list whose sum must equal the numerator in every row (verified in assemble()).
    """
    name: str
    table: str
    kind: str
    description: str
    cells: list[str] = field(default_factory=list)          # estimate cell(s) or numerator cells
    denominator: str | None = None
    label_rule: str = ""                                     # regex the numerator cell labels must match
    equivalent: str = ""                                     # subject/profile variable the brief named
    zcta: bool = True
    subtract: list[str] = field(default_factory=list)       # difference kind: cells subtracted from cells[0]
    cell_universes: dict[str, str] = field(default_factory=dict)
    check_cells: list[str] = field(default_factory=list)
    check_label_rule: str = ""

    def all_cells(self) -> list[str]:
        return [*self.cells, *self.subtract, *([self.denominator] if self.denominator else []), *self.check_cells]


def _b01001_age_cells(shells: pd.DataFrame, lo: int, hi: int | None) -> list[str]:
    """B01001 sex-by-age cells whose lower age bound is in [lo, hi] (hi=None: open)."""
    rows = shells[(shells["Table ID"] == "B01001") & (shells["Indent"] == "2")]
    out = []
    for uid, label in zip(rows["Unique ID"], rows["Label"]):
        m = re.match(r"(Under )?(\d+)", label)
        if not m:
            raise ValueError(f"unparsed B01001 label {label!r}")
        low = 0 if m.group(1) else int(m.group(2))
        if low >= lo and (hi is None or low <= hi):
            out.append(uid)
    return out


def _cells_matching(shells: pd.DataFrame, table: str, pattern: str) -> list[str]:
    rows = shells[(shells["Table ID"] == table) & shells["Label"].str.fullmatch(pattern)]
    return rows["Unique ID"].tolist()


def build_measures(shells: pd.DataFrame) -> list[Measure]:
    """Resolve cell lists from the table shells so labels, not memory, define each numerator."""
    return [
        Measure("total_population", "B01003", "estimate", "Total population", ["B01003_001"],
                label_rule=r"Total", equivalent="B01003_001E"),
        Measure("median_age", "B01002", "estimate", "Median age (years), total population", ["B01002_001"],
                label_rule=r"Total:", equivalent="B01002_001E"),
        Measure("age_18_64", "B01001", "difference",
                "Population aged 18-64 (total population minus under 18 minus 65 and over; share of total)",
                ["B01001_001"], "B01001_001", label_rule=r"Total:",
                equivalent="S0101 (18 years and over minus 65+)", zcta=False,
                subtract=["B09001_001", "B09020_001"],
                cell_universes={"B01001_001": "Total population", "B09001_001": "Population under 18 years",
                                "B09020_001": "Population 65 years and over"},
                check_cells=_b01001_age_cells(shells, 18, 64),
                check_label_rule=r"(18 and 19|20|21|22 to 24|25 to 29|30 to 34|35 to 39|40 to 44|45 to 49|"
                                 r"50 to 54|55 to 59|60 and 61|62 to 64) years"),
        Measure("age_65_plus", "B09020", "share", "Population aged 65 and over (share of total population)",
                ["B09020_001"], "B01001_001", label_rule=r"Total:", equivalent="S0101_C02_030E",
                cell_universes={"B09020_001": "Population 65 years and over"},
                check_cells=_b01001_age_cells(shells, 65, None),
                check_label_rule=r"(65 and 66|67 to 69|70 to 74|75 to 79|80 to 84) years|85 years and over"),
        Measure("with_disability", "B18101", "share",
                "Civilian noninstitutionalized population with a disability",
                _cells_matching(shells, "B18101", "With a disability"), "B18101_001",
                label_rule=r"With a disability", equivalent="S1810_C03_001E"),
        Measure("uninsured", "B27020", "share",
                "Civilian noninstitutionalized population with no health insurance coverage",
                _cells_matching(shells, "B27020", "No health insurance coverage"), "B27020_001",
                label_rule=r"No health insurance coverage", equivalent="S2701_C05_001E"),
        Measure("below_poverty", "B17001", "share",
                "Population for whom poverty status is determined, income below poverty level",
                ["B17001_002"], "B17001_001",
                label_rule=r"Income in the past 12 months below poverty level:", equivalent="S1701_C03_001E"),
        Measure("median_household_income", "B19013", "estimate",
                "Median household income in the past 12 months (inflation-adjusted dollars of the final year)",
                ["B19013_001"], label_rule=r"Median household income in the past 12 months.*",
                equivalent="B19013_001E"),
        Measure("households_no_vehicle", "B08201", "share", "Households with no vehicle available",
                ["B08201_002"], "B08201_001", label_rule=r"No vehicle available", equivalent="DP04_0058PE"),
        Measure("households_broadband", "B28002", "share",
                "Households with a broadband Internet subscription of any type",
                ["B28002_004"], "B28002_001", label_rule=r"Broadband of any type", equivalent="S2801_C02_014E"),
    ]


EXPECTED_DENOMINATOR_UNIVERSE = {
    "B01001": "Total population",
    "B18101": "Civilian noninstitutionalized population",
    "B27020": "Civilian noninstitutionalized population",
    "B17001": "Population for whom poverty status is determined",
    "B08201": "Households",
    "B28002": "Households",
}


def verify_labels(measures: list[Measure], shells: pd.DataFrame) -> dict:
    """Check every cell against the published table shells; raise on any mismatch."""
    by_uid = shells.set_index("Unique ID")
    report = {}
    for m in measures:
        if not m.cells:
            raise ValueError(f"{m.name}: no cells resolved from table shells")
        for uid in [*m.cells, *m.subtract]:
            if uid not in by_uid.index:
                raise ValueError(f"{m.name}: {uid} not in the table shells")
            label = by_uid.loc[uid, "Label"]
            if not re.fullmatch(m.label_rule, label):
                raise ValueError(f"{m.name}: {uid} label {label!r} !~ {m.label_rule!r}")
        for uid, universe in m.cell_universes.items():
            if by_uid.loc[uid, "Universe"] != universe:
                raise ValueError(f"{m.name}: {uid} universe {by_uid.loc[uid, 'Universe']!r} != {universe!r}")
        for uid in m.check_cells:
            label = by_uid.loc[uid, "Label"]
            if not re.fullmatch(m.check_label_rule, label):
                raise ValueError(f"{m.name}: check cell {uid} label {label!r} !~ {m.check_label_rule!r}")
        if m.denominator:
            d = by_uid.loc[m.denominator]
            den_table = m.denominator.split("_")[0]
            if d["Label"] != "Total:" or d["Universe"] != EXPECTED_DENOMINATOR_UNIVERSE[den_table]:
                raise ValueError(f"{m.name}: denominator {m.denominator} is {d['Label']!r}/{d['Universe']!r}")
        report[m.name] = {
            "table": m.table,
            "table_title": str(by_uid.loc[m.cells[0], "Title"]),
            "universe": str(by_uid.loc[m.cells[0], "Universe"]),
            "cells": {uid: f"{by_uid.loc[uid, 'Label']} [{by_uid.loc[uid, 'Title']}; {by_uid.loc[uid, 'Universe']}]"
                      if m.cell_universes else str(by_uid.loc[uid, "Label"]) for uid in m.cells},
            "subtract": {uid: f"{by_uid.loc[uid, 'Label']} [{by_uid.loc[uid, 'Title']}; {by_uid.loc[uid, 'Universe']}]"
                         for uid in m.subtract},
            "check_cells": list(m.check_cells),
            "denominator": m.denominator,
            "equivalent_requested_variable": m.equivalent,
        }
    return report


# ----------------------------------------------------------------------------- MOE arithmetic
def moe_sum(est: pd.DataFrame, moe: pd.DataFrame) -> pd.Series:
    """ACS handbook approximation for the MOE of a sum of cells.

    Root-sum-of-squares, except that among cells whose estimate is 0 only the
    largest MOE is included (Census ACS General Handbook 2020, ch. 8).
    """
    e = est.to_numpy(dtype=float)
    m = moe.to_numpy(dtype=float)
    zero = e == 0
    nonzero_sq = np.where(zero, 0.0, m ** 2).sum(axis=1)
    zero_max = np.where(zero, m, -np.inf).max(axis=1)
    zero_max = np.where(np.isfinite(zero_max), zero_max, 0.0)
    out = np.sqrt(nonzero_sq + zero_max ** 2)
    out[np.isnan(m).any(axis=1)] = np.nan
    return pd.Series(out, index=est.index)


def proportion_with_moe(num: pd.Series, num_moe: pd.Series, den: pd.Series, den_moe: pd.Series
                        ) -> tuple[pd.Series, pd.Series]:
    """Proportion (0-1) and its approximate MOE (ACS handbook proportion/ratio formulas)."""
    num, num_moe, den, den_moe = (s.astype(float) for s in (num, num_moe, den, den_moe))
    den_safe = den.where(den > 0)
    p = num / den_safe
    rad = num_moe ** 2 - p ** 2 * den_moe ** 2
    rad = rad.where(rad >= 0, num_moe ** 2 + p ** 2 * den_moe ** 2)
    return p, np.sqrt(rad) / den_safe


# ----------------------------------------------------------------------------- retrieval
def select_vintage(candidates=VINTAGE_CANDIDATES) -> tuple[int, dict]:
    """Newest ACS 5-year vintage whose (keyless) API metadata endpoint exists."""
    tried = {}
    for year in candidates:
        url = API_DATASET.format(year=year)
        try:
            r = get(url, reject_html=False)
        except Exception as exc:  # offline / unreachable: record and fall back to an older (cached) vintage
            tried[year] = f"error: {type(exc).__name__}"
            continue
        tried[year] = r.status
        if r.status == 200:
            ds = r.json()["dataset"][0]
            return year, {"api_dataset_url": url, "title": ds.get("title"), "c_vintage": ds.get("c_vintage"),
                          "modified": ds.get("modified"), "probe_status": tried}
    raise RuntimeError(f"no ACS 5-year vintage found among {candidates}: {tried}")


def probe_api_key_requirement(year: int) -> dict:
    """Record whether api.census.gov data requests work without a key (they did not on 2026-09-23)."""
    url = API_PROBE.format(year=year)
    try:
        r = get(url, reject_html=False, refresh=True)
    except Exception as exc:  # network failure: record, do not fail the build
        return {"url": url, "error": repr(exc), "checked_at": utc_now_iso()}
    needs_key = r.url.rstrip("/").endswith("missing_key.html") or r.headers.get("X-DataWebAPI-KeyError") == "1"
    return {"url": url, "final_url": r.url, "final_status": r.status, "needs_key": bool(needs_key),
            "checked_at": utc_now_iso(),
            "page_text": "A valid key must be included with each data API request." if needs_key else ""}


def sf_urls(year: int, tables: list[str]) -> dict[str, str]:
    root = SF_ROOT.format(year=year)
    urls = {f"ACS{year}5YR_Table_Shells.txt": f"{root}/documentation/ACS{year}5YR_Table_Shells.txt"}
    for t in tables:
        urls[f"acsdt5y{year}-{t.lower()}.dat"] = f"{root}/data/5YRData/acsdt5y{year}-{t.lower()}.dat"
    return urls


TABLES = ["B01001", "B01002", "B01003", "B08201", "B09001", "B09020", "B17001", "B18101", "B19013", "B27020",
          "B28002"]


def fetch(year: int | None = None) -> dict[str, Path]:
    """Download table shells + detailed-table Summary File files (skips files already in MANIFEST)."""
    if year is None:
        year, _ = select_vintage()
    out = {}
    for fname, url in sf_urls(year, TABLES).items():
        out[fname] = download_file(url, SOURCE_ID, fname, min_bytes=10_000, timeout=900)
    return out


def read_shells(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="|", dtype=str, keep_default_na=False)


def read_sf_table(path: Path, table: str, levels=("state", "county", "zcta")) -> pd.DataFrame:
    """Read one Summary File table, keeping only the requested summary levels, numeric, sentinels intact."""
    prefixes = [SUMLEVEL_PREFIX[lv] for lv in levels]
    con = duckdb.connect()
    df = con.execute(
        "SELECT * FROM read_csv(?, delim='|', header=true, all_varchar=true) WHERE left(GEO_ID, 9) IN ("
        + ",".join("?" * len(prefixes)) + ")", [str(path), *prefixes]).df()
    con.close()
    rename = {}
    for c in df.columns:
        m = re.fullmatch(rf"({table})_([EM])(\d{{3}})", c)
        if m:
            rename[c] = f"{m.group(1)}_{m.group(3)}{m.group(2)}"   # B01003_E001 -> B01003_001E
    df = df.rename(columns=rename)
    for c in rename.values():
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def clean_sentinels(df: pd.DataFrame, cols: list[str]) -> tuple[pd.DataFrame, dict]:
    """Estimates: every ACS special value -> NaN. MOEs: controlled (-555555555) -> 0, others -> NaN.

    Returns the cleaned frame and a count of every special value per column.
    """
    out = df.copy()
    counts: dict[str, dict[str, int]] = {}
    for c in cols:
        v = out[c]
        special = v.isin(list(ACS_SENTINELS)) | (v <= -100_000_000)
        if special.any():
            counts[c] = {str(int(k)): int(n) for k, n in v[special].value_counts().items()}
        if c.endswith("M"):
            out[c] = v.mask(v == CONTROLLED_MOE, 0.0).mask(special & (v != CONTROLLED_MOE))
        else:
            out[c] = v.mask(special)
    return out, counts


# ----------------------------------------------------------------------------- build
def assemble(frames: dict[str, pd.DataFrame], measures: list[Measure]) -> tuple[pd.DataFrame, dict, dict]:
    """Merge per-table frames on GEO_ID and compute the output columns for every measure.

    Returns (wide frame, special-value counts, numerator cross-checks). A cross-check
    compares a measure's numerator with the sum of its independent check_cells; any
    row where both are present and differ raises, because it would mean the chosen
    cells do not measure the concept the check cells define.
    """
    base = None
    sentinel_counts = {}
    for table, df in frames.items():
        needed = sorted({c for m in measures for c in m.all_cells() if c.split("_")[0] == table})
        if not needed:
            continue
        cols = [f"{c}{s}" for c in needed for s in "EM"]
        sub, counts = clean_sentinels(df[["GEO_ID", *cols]], cols)
        sentinel_counts.update(counts)
        base = sub if base is None else base.merge(sub, on="GEO_ID", how="outer", validate="one_to_one")
    out = pd.DataFrame({"GEO_ID": base["GEO_ID"]})
    checks = {}
    for m in measures:
        if m.kind == "estimate":
            c = m.cells[0]
            out[m.name] = base[f"{c}E"]
            out[f"{m.name}_moe"] = base[f"{c}M"]
            continue
        if m.kind == "difference":
            parts = [m.cells[0], *m.subtract]
            est = base[[f"{c}E" for c in parts]]
            moe = base[[f"{c}M" for c in parts]]
            n = (est.iloc[:, 0] - est.iloc[:, 1:].sum(axis=1)).where(est.notna().all(axis=1))
            n_moe = moe_sum(est, moe)       # handbook: MOE of a difference = MOE of the sum
        elif m.kind == "share":
            est = base[[f"{c}E" for c in m.cells]]
            moe = base[[f"{c}M" for c in m.cells]]
            n = est.sum(axis=1, min_count=len(m.cells)).where(est.notna().all(axis=1))
            n_moe = base[f"{m.cells[0]}M"] if len(m.cells) == 1 else moe_sum(est, moe)
        else:
            raise ValueError(f"{m.name}: unknown kind {m.kind!r}")
        if m.check_cells:
            chk_est = base[[f"{c}E" for c in m.check_cells]]
            chk = chk_est.sum(axis=1, min_count=len(m.check_cells)).where(chk_est.notna().all(axis=1))
            both = n.notna() & chk.notna()
            differ = both & (n != chk)
            checks[m.name] = {"check_cells": m.check_cells, "rows_compared": int(both.sum()),
                              "rows_differing": int(differ.sum())}
            if differ.any():
                raise ValueError(f"{m.name}: numerator differs from the sum of {m.check_cells} in "
                                 f"{int(differ.sum())} rows, e.g. {base.loc[differ, 'GEO_ID'].head(5).tolist()}")
        den, den_moe = base[f"{m.denominator}E"], base[f"{m.denominator}M"]
        p, p_moe = proportion_with_moe(n, n_moe, den, den_moe)
        out[f"n_{m.name}"] = n
        out[f"n_{m.name}_moe"] = n_moe
        out[f"{m.name}_universe"] = den
        out[f"{m.name}_universe_moe"] = den_moe
        out[f"pct_{m.name}"] = 100 * p
        out[f"pct_{m.name}_moe"] = 100 * p_moe
    if "median_household_income" in out:
        inc = out["median_household_income"]
        out["median_household_income_open_ended"] = np.select(
            [inc == INCOME_TOPCODE, inc == INCOME_BOTTOMCODE], ["top_coded_250000_plus", "bottom_coded_under_2500"],
            default="")
    return out, sentinel_counts, checks


def _geo_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    g = read_table("geographies")
    z = read_table("zcta_centroids", columns=["zcta", "aland_m2", "county_fips", "state_fips", "state_abbr"])
    return g, z


def county_join_coverage(fips: pd.Series, source_label: str) -> dict:
    """Coverage of a set of 5-digit county FIPS against the canonical 2024 counties (+ legacy CT)."""
    g = read_table("geographies")
    cty = g[(g["geo_level"] == "county") & (~g["ct_legacy"])]
    legacy = set(g.loc[g["ct_legacy"], "geo_id"])
    canon = set(cty["geo_id"])
    src = set(fips.dropna().astype(str))
    matched = src & canon
    src_unmatched = sorted(src - canon)
    canon_unmatched = sorted(canon - src)
    abbr = dict(zip(cty["geo_id"], cty["state_abbr"]))
    canon_unmatched_by_state = pd.Series([abbr[f] for f in canon_unmatched]).value_counts().to_dict() \
        if canon_unmatched else {}
    ct_src = sorted(f for f in src if f.startswith("09"))
    return {
        "source": source_label,
        "source_counties": len(src),
        "canonical_counties_2024": len(canon),
        "matched": len(matched),
        "source_unmatched": len(src_unmatched),
        "source_unmatched_fips": src_unmatched[:50],
        "source_unmatched_legacy_ct": sorted(set(src_unmatched) & legacy),
        "canonical_without_source_row": len(canon_unmatched),
        "canonical_without_source_row_by_state": {k: int(v) for k, v in canon_unmatched_by_state.items()},
        "canonical_without_source_row_fips": canon_unmatched[:60],
        "connecticut_source_fips": ct_src,
        "connecticut_vintage": ("planning regions (09110-09190), the canonical 2024 CT set" if ct_src and all(f >= "09100" for f in ct_src)
                                else "legacy counties (09001-09015)" if ct_src and all(f < "09100" for f in ct_src)
                                else "mixed" if ct_src else "absent"),
        "connecticut_matched": sorted(set(ct_src) & canon),
    }


def build(year: int | None = None) -> dict:
    vintage_year, ds_meta = select_vintage() if year is None else (year, {"api_dataset_url": API_DATASET.format(year=year)})
    api_probe = probe_api_key_requirement(vintage_year)
    files = fetch(vintage_year)
    manifest = load_manifest(SOURCE_ID)["files"]
    shells = read_shells(files[f"ACS{vintage_year}5YR_Table_Shells.txt"])
    measures = build_measures(shells)
    label_report = verify_labels(measures, shells)

    frames = {t: read_sf_table(files[f"acsdt5y{vintage_year}-{t.lower()}.dat"], t) for t in TABLES}
    wide, sentinel_counts, numerator_checks = assemble(frames, measures)
    wide["sumlevel_prefix"] = wide["GEO_ID"].str[:9]
    wide["geo_code"] = wide["GEO_ID"].str[9:]
    acs_label = f"{vintage_year - 4}-{vintage_year}"
    retrieved = min(manifest[f]["retrieved_at"] for f in files)
    sf_last_modified = sorted({manifest[f].get("last_modified") or "" for f in files if f.endswith(".dat")})

    g, z = _geo_frames()
    # ---------------- state + county
    sc = wide[wide["sumlevel_prefix"].isin([SUMLEVEL_PREFIX["state"], SUMLEVEL_PREFIX["county"]])].copy()
    sc["geo_level"] = np.where(sc["sumlevel_prefix"] == SUMLEVEL_PREFIX["state"], "state", "county")
    sc["geo_id"] = sc["geo_code"]
    canon = g[~g["ct_legacy"]][["geo_id", "geo_level", "name", "state_fips", "state_abbr", "aland_m2"]]
    sc = sc.merge(canon, on=["geo_id", "geo_level"], how="left", validate="one_to_one")
    sc["state_fips"] = sc["state_fips"].fillna(sc["geo_id"].str[:2])
    sc["county_fips"] = np.where(sc["geo_level"] == "county", sc["geo_id"], None)
    _add_density(sc)
    sc["acs_vintage"] = acs_label
    sc["acs_release_year"] = vintage_year
    measure_cols = [c for c in wide.columns if c not in ("GEO_ID", "sumlevel_prefix", "geo_code")]
    sc = sc[["geo_id", "geo_level", "name", "state_fips", "state_abbr", "county_fips", "acs_vintage",
             "acs_release_year", *measure_cols, "aland_m2", "aland_km2", "pop_density_per_km2",
             "pop_density_per_sqmi", "GEO_ID"]].sort_values(["geo_level", "geo_id"], ascending=[False, True])
    notes = (f"ACS {acs_label} 5-year estimates from the Census table-based Summary File (detailed tables). "
             "pct_* derived from detailed-table cells (see DATA_AUDIT.md for cell lists); *_moe are 90% MOEs, "
             "approximated with ACS handbook formulas for sums and proportions. Population density = "
             "total_population / 2024 Census land area. Ecological (area-level) context, not person-level data.")
    sc = add_provenance(
        sc.reset_index(drop=True), data_layer="geographic",
        source_name="U.S. Census Bureau, American Community Survey 5-year table-based Summary File",
        source_version=f"ACS {acs_label} 5-year (release {vintage_year}); SF files last-modified {sf_last_modified[-1]}",
        retrieved_at=retrieved, evidence_type="survey_estimate", source_record_id="GEO_ID",
        source_geographic_resolution="county", provenance_notes=notes,
    )
    sc.loc[sc["geo_level"] == "state", "source_geographic_resolution"] = "state"
    sc = sc.drop(columns=["GEO_ID"])
    write_table(sc, "geo_context__acs", producer="ingestion.census_acs.build",
                description=f"ACS {acs_label} 5-year context, states + counties, wide with MOE columns")

    # ---------------- ZCTA (smaller set)
    zc = wide[wide["sumlevel_prefix"] == SUMLEVEL_PREFIX["zcta"]].copy()
    zc["zcta"] = zc["geo_code"]
    z_measures = [m for m in measures if m.zcta]
    zcols = []
    for m in z_measures:
        zcols += [m.name, f"{m.name}_moe"] if m.kind == "estimate" else [f"pct_{m.name}", f"pct_{m.name}_moe"]
    zcols.append("median_household_income_open_ended")
    zc = zc.merge(z.rename(columns={"county_fips": "county_fips_2024_of_internal_point"}),
                  on="zcta", how="left", validate="one_to_one")
    _add_density(zc)
    zc["acs_vintage"] = acs_label
    zc["acs_release_year"] = vintage_year
    zc = zc[["zcta", "state_fips", "state_abbr", "county_fips_2024_of_internal_point", "acs_vintage",
             "acs_release_year", *zcols, "aland_m2", "aland_km2", "pop_density_per_km2", "pop_density_per_sqmi",
             "GEO_ID"]].sort_values("zcta")
    zc = add_provenance(
        zc.reset_index(drop=True), data_layer="geographic",
        source_name="U.S. Census Bureau, American Community Survey 5-year table-based Summary File",
        source_version=f"ACS {acs_label} 5-year (release {vintage_year}); 2020 ZCTA definitions",
        retrieved_at=retrieved, evidence_type="survey_estimate", source_record_id="GEO_ID",
        source_geographic_resolution="zcta",
        provenance_notes=notes + " ZCTA estimates have large MOEs for small populations; state/county "
                                 "columns locate the ZCTA internal point only (ZCTAs can span counties/states).",
    ).drop(columns=["GEO_ID"])
    write_table(zc, "geo_context_zcta__acs", producer="ingestion.census_acs.build",
                description=f"ACS {acs_label} 5-year context by ZCTA (smaller variable set)")

    stats = _stats(sc, zc, z, measures, label_report, sentinel_counts, files, manifest, ds_meta, api_probe,
                   vintage_year, acs_label)
    stats["numerator_checks"] = numerator_checks
    (raw_dir(SOURCE_ID) / "build_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    _write_registry(stats, files, manifest, vintage_year, acs_label, retrieved)
    write_audit(stats)
    return stats


def _add_density(df: pd.DataFrame) -> None:
    df["aland_km2"] = df["aland_m2"] / 1e6
    area = df["aland_m2"].where(df["aland_m2"] > 0)
    df["pop_density_per_km2"] = df["total_population"] / (area / 1e6)
    df["pop_density_per_sqmi"] = df["total_population"] / (area / SQM_PER_SQMI)


def _missing(df: pd.DataFrame, cols: list[str]) -> dict:
    return {c: {"n_missing": int(df[c].isna().sum()), "pct_missing": round(100 * df[c].isna().mean(), 3)}
            for c in cols}


def _stats(sc, zc, z, measures, label_report, sentinel_counts, files, manifest, ds_meta, api_probe,
           year, acs_label) -> dict:
    cty = sc[sc["geo_level"] == "county"]
    st = sc[sc["geo_level"] == "state"]
    key = ["total_population", "median_age", "pct_age_18_64", "pct_age_65_plus", "pct_with_disability",
           "pct_uninsured", "pct_below_poverty", "median_household_income", "pct_households_no_vehicle",
           "pct_households_broadband", "pop_density_per_km2"]
    zkey = [c for c in key if c in zc.columns]
    zset = set(zc["zcta"])
    zgaz = set(z["zcta"])
    return {
        "acs_vintage": acs_label, "release_year": year, "api_metadata": ds_meta, "api_probe": api_probe,
        "files": {f: {k: manifest[f].get(k) for k in ("url", "bytes", "sha256", "retrieved_at", "last_modified")}
                  for f in files},
        "rows": {"county": int(len(cty)), "state": int(len(st)), "zcta": int(len(zc))},
        "state_ids": sorted(st["geo_id"]),
        "county_by_state_count": int(cty["state_fips"].nunique()),
        "missing_county": _missing(cty, key + [f"{k}_moe" for k in key if f"{k}_moe" in cty.columns]),
        "missing_state": _missing(st, key),
        "missing_zcta": _missing(zc, zkey),
        "sentinel_counts_all_levels": sentinel_counts,
        "income_open_ended": {
            "county": cty["median_household_income_open_ended"].replace("", np.nan).value_counts().to_dict(),
            "zcta": zc["median_household_income_open_ended"].replace("", np.nan).value_counts().to_dict()},
        "zcta_zero_population": int((zc["total_population"] == 0).sum()),
        "zcta_population_lt_100": int((zc["total_population"] < 100).sum()),
        "zcta_median_rel_moe_pct_uninsured": float((zc["pct_uninsured_moe"] / zc["pct_uninsured"]).replace(
            [np.inf], np.nan).median()),
        "county_median_rel_moe_pct_uninsured": float((cty["pct_uninsured_moe"] / cty["pct_uninsured"]).replace(
            [np.inf], np.nan).median()),
        "county_summary": cty[key].describe().round(3).to_dict(),
        "national_population_sum_counties": int(cty["total_population"].sum()),
        "national_population_sum_states": int(st["total_population"].sum()),
        "county_coverage": county_join_coverage(cty["geo_id"], "ACS counties"),
        "state_coverage": {
            "acs_states": len(st), "matched_geographies": int(st["name"].notna().sum()),
            "geographies_states_without_acs": sorted(
                set(read_table("geographies").query("geo_level == 'state'")["geo_id"]) - set(st["geo_id"])),
        },
        "zcta_coverage": {"acs_zctas": len(zset), "gazetteer_2024_zctas": len(zgaz),
                          "matched": len(zset & zgaz), "acs_not_in_gazetteer": sorted(zset - zgaz)[:50],
                          "gazetteer_not_in_acs": len(zgaz - zset),
                          "gazetteer_not_in_acs_sample": sorted(zgaz - zset)[:30]},
        "label_report": label_report,
        "county_missing_median_income": cty.loc[cty["median_household_income"].isna(),
                                                ["geo_id", "name", "state_abbr", "total_population"]].to_dict("records"),
        "cross_check_svi": _cross_check_svi(cty),
        "zcta_without_internal_point_state": sorted(zc.loc[zc["state_fips"].isna(), "zcta"]),
    }


# SVI 2022 inputs (ACS 2018-2022) that measure the same concepts as columns here (ACS later period).
SVI_PAIRS = [("total_population", "e_totpop"), ("pct_uninsured", "ep_uninsur"),
             ("pct_with_disability", "ep_disabl"), ("pct_age_65_plus", "ep_age65"),
             ("pct_households_no_vehicle", "ep_noveh"), ("pct_below_poverty", "ep_pov150"),
             ("pct_households_broadband", "ep_noint")]


# SVI 2022 publishes the subject/profile-table MOE of each percent (mp_*); comparing them with the
# approximated pct_*_moe here measures how far the handbook approximation drifts (different period, so ~1 is
# the expectation, not an identity).
SVI_MOE_PAIRS = [("pct_uninsured_moe", "mp_uninsur"), ("pct_with_disability_moe", "mp_disabl"),
                 ("pct_age_65_plus_moe", "mp_age65"), ("pct_households_no_vehicle_moe", "mp_noveh"),
                 ("pct_households_broadband_moe", "mp_noint")]


def _cross_check_svi(cty: pd.DataFrame) -> dict:
    """Spearman correlation with SVI 2022 inputs, when geo_vulnerability has been built (else empty)."""
    try:
        svi = read_table("geo_vulnerability",
                         columns=["geo_id"] + [b for _, b in SVI_PAIRS] + [b for _, b in SVI_MOE_PAIRS])
    except FileNotFoundError:
        return {}
    m = cty.merge(svi, on="geo_id")
    moe = {}
    for a, b in SVI_MOE_PAIRS:
        r = (m[a] / m[b].where(m[b] > 0)).dropna()
        moe[f"{a} / svi.{b}"] = {"median_ratio": round(float(r.median()), 3),
                                 "p25": round(float(r.quantile(0.25)), 3), "p75": round(float(r.quantile(0.75)), 3),
                                 "n": int(len(r))}
    return {"n_counties": int(len(m)),
            "spearman": {f"{a} vs svi.{b}": round(float(m[a].corr(m[b], method="spearman")), 3)
                         for a, b in SVI_PAIRS},
            "moe_ratio_vs_published": moe}


def _write_registry(stats, files, manifest, year, acs_label, retrieved) -> None:
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": f"American Community Survey 5-year estimates ({acs_label}), table-based Summary File",
        "publisher": "U.S. Census Bureau",
        "landing_url": LANDING_URL,
        "access_urls": [manifest[f]["url"] for f in files],
        "license": "U.S. federal government work (public domain)",
        "access_conditions": ("Summary File: open download, no registration. api.census.gov data requests "
                              "required an API key at retrieval (HTTP 302 -> missing_key.html); no key was "
                              "requested, so the API was not used for data."),
        "retrieved_at": retrieved,
        "source_version": f"ACS {acs_label} 5-year (release {year})",
        "update_date": "annual (5-year release each December)",
        "data_layer": "geographic",
        "unit_of_observation": "state / county / ZCTA (area-level survey estimates)",
        "sample_size": {"states": stats["rows"]["state"], "counties": stats["rows"]["county"],
                        "zctas": stats["rows"]["zcta"]},
        "geographic_resolution": "state, county (2024 vintage; CT planning regions), zcta (2020 ZCTAs)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (area-level estimates)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": ["geo_context__acs", "geo_context_zcta__acs"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.census_acs",
        "limitations": [
            "Percentages derived from detailed-table cells, not downloaded subject/profile tables; MOEs approximated",
            "5-year pooled estimates (period estimate, not a single year)",
            "ZCTA estimates carry large MOEs for small populations",
            "Area-level context only; never person-level",
        ],
    })


# ----------------------------------------------------------------------------- audit
def _fmt_missing(d: dict) -> str:
    rows = ["| Variable | Missing (n) | Missing (%) |", "|---|---|---|"]
    rows += [f"| {k} | {v['n_missing']} | {v['pct_missing']} |" for k, v in d.items()]
    return "\n".join(rows)


def write_audit(stats: dict) -> Path:
    s = stats
    files = "\n".join(f"| {f} | {v['url']} | {v['bytes']:,} | `{v['sha256'][:16]}…` | {v['last_modified']} |"
                      for f, v in s["files"].items())
    lab = []
    for name, r in s["label_report"].items():
        cells = ", ".join(r["cells"])
        if r.get("subtract"):
            cells = f"{cells} minus " + " minus ".join(f"{k} ({v})" for k, v in r["subtract"].items())
        elif r["cells"] and any("[" in v for v in r["cells"].values()):
            cells = ", ".join(f"{k} ({v})" for k, v in r["cells"].items())
        lab.append(f"| {name} | {r['table']} ({r['table_title']}) | {r['universe']} | {cells} | "
                   f"{r['denominator'] or '—'} | {r['equivalent_requested_variable']} |")
    chk = "\n".join(f"| {k} | {', '.join(v['check_cells'])} | {v['rows_compared']} | {v['rows_differing']} |"
                     for k, v in s.get("numerator_checks", {}).items()) or "| none | | | |"
    moe_cmp = s["cross_check_svi"].get("moe_ratio_vs_published", {})
    moe_tab = "\n".join(f"| {k} | {v['median_ratio']} | {v['p25']}–{v['p75']} | {v['n']} |"
                         for k, v in moe_cmp.items()) or "| not computed (geo_vulnerability not built yet) | | | |"
    cov = s["county_coverage"]
    zc = s["zcta_coverage"]
    probe = s["api_probe"]
    sent = "\n".join(f"| {c} | {v} |" for c, v in sorted(s["sentinel_counts_all_levels"].items())) or "| none | |"
    text = f"""# DATA AUDIT — American Community Survey 5-year ({s['acs_vintage']})

| Field | Value |
|---|---|
| source_id | census_acs |
| Source (dataset/API name, exact files/endpoints) | ACS 5-year table-based Summary File, detailed tables {', '.join(TABLES)} (files listed below). api.census.gov used only for keyless metadata (vintage detection). |
| Publishing organization | U.S. Census Bureau |
| Retrieval date (UTC) | {min(v['retrieved_at'] for v in s['files'].values())} |
| Source version / release | ACS {s['acs_vintage']} 5-year (release year {s['release_year']}); API dataset metadata `modified` = {s['api_metadata'].get('modified')} |
| Source update date / cadence | Summary File .dat Last-Modified {sorted({v['last_modified'] for v in s['files'].values() if v['last_modified']})[-1]}; new 5-year release each December |
| License / access conditions | Public domain (U.S. federal work). Summary File: open download. api.census.gov data endpoints required a key at retrieval (see below). |
| Unit of observation | Geographic area (state, county, ZCTA) — survey estimates |
| Sample size (actual, as ingested) | {s['rows']['state']} states (50 + DC + PR), {s['rows']['county']} counties/county-equivalents, {s['rows']['zcta']} ZCTAs |
| Geography (resolution, vintage) | State; county on the 2024 vintage (Connecticut = 9 planning regions); ZCTA on 2020 ZCTA definitions |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — area-level estimates; no persons |

## Why the Summary File and not api.census.gov

The brief asked for `api.census.gov` without a key. On {probe.get('checked_at')} the request
`{probe.get('url')}` ended at `{probe.get('final_url')}` (HTTP {probe.get('final_status')}; the first hop is a
302 with header `X-DataWebAPI-KeyError: 1`). The page says: "{probe.get('page_text')}". needs_key = {probe.get('needs_key')}.
Requesting a key requires registering an e-mail address, which this pipeline does not do on the user's behalf.
The Census Bureau publishes the same ACS estimates keylessly as the table-based Summary File, which is used here.
Vintage detection: the keyless metadata endpoint `https://api.census.gov/data/<year>/acs/acs5` returned
{s['api_metadata'].get('probe_status')} (HTTP status by year); the newest available vintage was used.

## Files / endpoints retrieved

| File | URL | Bytes | sha256 | Last-Modified |
|---|---|---|---|---|
{files}

Full sha256 values are in `MANIFEST.json`. Each .dat file contains every summary level; only
`0400000US` (state), `0500000US` (county) and `860Z200US` (ZCTA) rows are read.

## Key variables

Cell lists were resolved from the published table shells (`ACS{s['release_year']}5YR_Table_Shells.txt`) and
every label was checked in code (`verify_labels`) before use.

| Output concept | Table | Universe | Numerator / estimate cells | Denominator | Variable named in the brief |
|---|---|---|---|---|---|
{chr(10).join(lab)}

Output columns (geo_context__acs): for each share concept `n_<x>`, `n_<x>_moe`, `<x>_universe`,
`<x>_universe_moe`, `pct_<x>` (0–100), `pct_<x>_moe`; for single estimates `<x>` and `<x>_moe`.
Also `median_household_income_open_ended` (Summary File jam values 250001 = "250,000+", 2499 = "2,500-"),
`aland_m2`, `aland_km2`, `pop_density_per_km2`, `pop_density_per_sqmi` (land area from the `geographies`
table, 2024 cartographic boundaries). ZCTA table keeps: total_population, median_age, pct_age_65_plus,
pct_with_disability, pct_uninsured, pct_below_poverty, median_household_income, pct_households_no_vehicle,
pct_households_broadband (+ MOEs) and density (land area from the 2024 ZCTA Gazetteer).

MOEs are 90% margins. Single cells use the published MOE. Sums and differences of cells use root-sum-of-squares
with only the largest MOE among zero-estimate cells; proportions use the ACS handbook proportion formula (ratio
formula when the radicand is negative). These are approximations of the replicate-weight MOEs Census publishes in
S0101, S1810, S2701, S1701, S2801 and DP04, and can differ from them. `-555555555` MOEs (controlled estimates, "MOE
not appropriate") are stored as 0; every other special value is set to null.

**Age numerators.** The 65+ count is the single published cell B09020_001 (population 65 years and over) and the
18–64 count is B01001_001 − B09001_001 (population under 18) − B09020_001, so each MOE combines at most three
published MOEs. Summing the 12 (65+) or 26 (18–64) B01001 sex-by-age bands gives the same point estimates but a
root-sum-of-squares MOE that ignores the population controls on those bands; in the first build of this table
that overstated the 65+ share MOE about six-fold against the published S0101_C02_030M (median ratio 5.3 versus
SVI 2022's MP_AGE65). The B01001 band sums are kept as an in-code cross-check of the point estimates (the build
fails if any row differs):

| Measure | Cross-check cells (B01001 bands) | Rows compared | Rows differing |
|---|---|---|---|
{chk}

Approximated MOEs versus the MOEs CDC/ATSDR SVI 2022 took from the published subject/profile tables (ACS
2018–2022, so a ratio near 1 is expected, not an identity; county rows with a published MOE > 0):

| Ratio (this table / SVI published) | Median | IQR | Counties |
|---|---|---|---|
{moe_tab}

## Missingness

Special values found in the Summary File cells used (all three summary levels, before cleaning):

| Cell | Counts by special value |
|---|---|
{sent}

Counties (n = {s['rows']['county']}):

{_fmt_missing(s['missing_county'])}

States (n = {s['rows']['state']}):

{_fmt_missing(s['missing_state'])}

ZCTAs (n = {s['rows']['zcta']}):

{_fmt_missing(s['missing_zcta'])}

County with no median household income (Census special value -666666666): {s['county_missing_median_income']}.

Median household income open-ended values: counties {s['income_open_ended']['county'] or 'none'};
ZCTAs {s['income_open_ended']['zcta'] or 'none'}. ZCTAs with total population 0: {s['zcta_zero_population']};
below 100: {s['zcta_population_lt_100']}. Median relative MOE of pct_uninsured: counties
{s['county_median_rel_moe_pct_uninsured']:.3f}, ZCTAs {s['zcta_median_rel_moe_pct_uninsured']:.3f}.
Sum of county total_population = {s['national_population_sum_counties']:,}; sum of the 52 state rows =
{s['national_population_sum_states']:,}.

Consistency check against the CDC/ATSDR SVI 2022 inputs (ACS 2018-2022, an earlier overlapping period;
Spearman rho across {s['cross_check_svi'].get('n_counties', 'n/a')} counties; `ep_noint` is households *without*
internet, so a negative rho is expected; `ep_pov150` is below 150% of poverty, not 100%):
{s['cross_check_svi'].get('spearman', 'not computed (geo_vulnerability not built yet)')}.

## Linkage strategy

* County rows join to `geographies` on `geo_id` (5-digit FIPS, 2024 vintage). Coverage:
  ACS counties = {cov['source_counties']}, canonical 2024 counties = {cov['canonical_counties_2024']},
  matched = {cov['matched']}, ACS counties not in geographies = {cov['source_unmatched']} {cov['source_unmatched_fips']},
  canonical counties without an ACS row = {cov['canonical_without_source_row']}
  (by state: {cov['canonical_without_source_row_by_state']} — the ACS does not cover American Samoa, Guam,
  Northern Mariana Islands or U.S. Virgin Islands).
* Connecticut: ACS rows are on **{cov['connecticut_vintage']}** ({', '.join(cov['connecticut_source_fips'])});
  all {len(cov['connecticut_matched'])} match the canonical 2024 planning regions. No legacy CT county
  (09001–09015) rows exist in this vintage, so legacy-CT sources cannot be joined to these rows without an
  explicit crosswalk (none is applied).
* State rows join on 2-digit state FIPS: {s['state_coverage']['matched_geographies']} of {s['state_coverage']['acs_states']}
  matched; geographies states without ACS rows: {s['state_coverage']['geographies_states_without_acs']}.
* ZCTA rows join to `zcta_centroids` on the 5-digit ZCTA: {zc['matched']} of {zc['acs_zctas']} ACS ZCTAs are in the
  2024 Gazetteer; {zc['gazetteer_not_in_acs']} Gazetteer ZCTAs have no ACS row. `state_fips`/`state_abbr`/
  `county_fips_2024_of_internal_point` locate the ZCTA's internal point only (taken from `zcta_centroids`).
  {len(s['zcta_without_internal_point_state'])} ACS ZCTAs have null state and county there because the internal
  point falls outside the generalized county polygons (a `census_geography` limitation, not an ACS gap):
  {', '.join(s['zcta_without_internal_point_state'])}.
* All joins are ecological (area-level). ACS rows describe places, not the individuals in any person-level
  dataset in this project.

## Limitations and caveats

* ACS 5-year estimates are period estimates pooled over {s['acs_vintage']}; they are survey estimates with
  sampling error (MOE columns), not counts.
* Percentages are derived from detailed-table cells rather than read from subject/profile tables. The point
  estimates follow the same definitions (same universes and cells) but the MOEs are approximations.
* Disability (B18101) uses the ACS six-question disability definition (hearing, vision, cognitive, ambulatory,
  self-care, independent living). It is not a measure of any specific condition and is not disease prevalence.
* Uninsured and disability universes are the civilian noninstitutionalized population; poverty universe excludes
  institutionalized people, military group quarters and unrelated individuals under 15; vehicle and broadband
  shares are per household (occupied housing unit).
* Median household income is in {s['release_year']} inflation-adjusted dollars; medians cannot be aggregated
  across areas.
* Small counties and most ZCTAs have large relative MOEs; use the MOE columns before ranking areas.
* Population density uses 2024 land area; water area is excluded.

## Processed outputs

| Table | Rows |
|---|---|
| geo_context__acs | {s['rows']['state'] + s['rows']['county']} ({s['rows']['state']} state + {s['rows']['county']} county) |
| geo_context_zcta__acs | {s['rows']['zcta']} |

`build_stats.json` in this directory holds every number above in machine-readable form.

## Reproduce

`uv run python -m measure_it.ingestion.census_acs`
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


# ----------------------------------------------------------------------------- census_geography audit
def census_geography_stats() -> dict:
    """Measure the outputs of measure_it.geography.crosswalk (that module is not modified here)."""
    import io
    import zipfile

    import geopandas as gpd

    geo_src = "census_geography"
    manifest = load_manifest(geo_src)["files"]
    g = read_table("geographies")
    z = read_table("zcta_centroids")
    x = read_table("zcta_county_crosswalk")
    cty = g[(g["geo_level"] == "county") & (~g["ct_legacy"])]
    legacy = g[g["ct_legacy"]]
    st = g[g["geo_level"] == "state"]
    with zipfile.ZipFile(RAW / geo_src / "2024_Gaz_counties_national.zip") as zf:
        name = [n for n in zf.namelist() if n.endswith(".txt")][0]
        gaz = pd.read_csv(io.BytesIO(zf.read(name)), sep="\t", dtype=str)
    gaz.columns = [c.strip() for c in gaz.columns]
    no_gaz_point = sorted(set(cty["geo_id"]) - set(gaz["GEOID"]))
    top = x.sort_values("land_share_of_zcta", ascending=False).drop_duplicates("zcta")
    m = z.merge(top[["zcta", "county_fips_2020"]], on="zcta", how="left")
    non_ct = m[(m["state_fips"] != "09") & m["county_fips"].notna()]
    xw_counties = set(x["county_fips_2020"])
    cb = gpd.read_parquet(PROCESSED / "county_boundaries_2024.geoparquet")
    sb = gpd.read_parquet(PROCESSED / "state_boundaries_2024.geoparquet")
    share_sum = x.groupby("zcta")["land_share_of_zcta"].sum()
    return {
        "files": {f: {k: v.get(k) for k in ("url", "bytes", "sha256", "retrieved_at", "last_modified")}
                  for f, v in manifest.items()},
        "geographies_rows": int(len(g)),
        "states": int(len(st)), "state_abbrs": sorted(st["state_abbr"]),
        "counties_2024": int(len(cty)),
        "counties_by_group": {
            "50 states + DC": int((~cty["state_abbr"].isin(["PR", "AS", "GU", "MP", "VI"])).sum()),
            "PR": int((cty["state_abbr"] == "PR").sum()),
            "island areas (AS, GU, MP, VI)": int(cty["state_abbr"].isin(["AS", "GU", "MP", "VI"]).sum())},
        "ct_planning_regions": sorted(cty.loc[cty["state_fips"] == "09", "geo_id"]),
        "ct_legacy_counties": sorted(legacy["geo_id"]),
        "duplicate_geo_id": int(g["geo_id"].duplicated().sum()),
        "missing_geographies": {c: int(g[c].isna().sum()) for c in g.columns if g[c].isna().any()},
        "aland_nonpositive": int((g["aland_m2"] <= 0).sum()),
        "counties_without_gazetteer_point": no_gaz_point,
        "zcta_rows": int(len(z)), "zcta_duplicates": int(z["zcta"].duplicated().sum()),
        "zcta_without_county": int(z["county_fips"].isna().sum()),
        "zcta_without_state": int(z["state_fips"].isna().sum()),
        "zcta_without_state_but_with_county": int((z["state_fips"].isna() & z["county_fips"].notna()).sum()),
        "zcta_without_county_list": sorted(z.loc[z["county_fips"].isna(), "zcta"]),
        "zcta_by_territory": z["state_abbr"].value_counts().reindex(["PR", "VI", "GU", "MP", "AS"]).fillna(0)
                              .astype(int).to_dict(),
        "zcta_internal_point_vs_max_land_county_agreement_non_ct": round(float(
            (non_ct["county_fips"] == non_ct["county_fips_2020"]).mean()), 4),
        "zcta_internal_point_vs_max_land_county_n": int(len(non_ct)),
        "crosswalk_rows": int(len(x)), "crosswalk_zctas": int(x["zcta"].nunique()),
        "crosswalk_counties_2020": len(xw_counties),
        "crosswalk_zctas_multi_county": int((x.groupby("zcta").size() > 1).sum()),
        "crosswalk_land_share_sum_min_max": [round(float(share_sum.min()), 6), round(float(share_sum.max()), 6)],
        "crosswalk_land_share_missing": int(x["land_share_of_zcta"].isna().sum()),
        "crosswalk_counties_not_in_2024": sorted(xw_counties - set(cty["geo_id"])),
        "counties_2024_not_in_crosswalk": sorted(set(cty["geo_id"]) - xw_counties),
        "county_boundaries_rows": int(len(cb)), "state_boundaries_rows": int(len(sb)),
        "county_boundaries_valid_share": round(float(cb.geometry.is_valid.mean()), 4),
        "boundaries_crs": str(cb.crs.to_epsg()),
    }


def write_census_geography_audit() -> Path:
    s = census_geography_stats()
    files = "\n".join(f"| {f} | {v['url']} | {v['bytes']:,} | `{v['sha256'][:16]}…` | {v['last_modified']} |"
                      for f, v in s["files"].items())
    retrieved = min(v["retrieved_at"] for v in s["files"].values())
    text = f"""# DATA AUDIT — Census geography (boundaries, Gazetteer, ZCTA–county relationship)

Audit of the outputs built by `measure_it.geography.crosswalk` (measured by
`measure_it.ingestion.census_acs.census_geography_stats`; the crosswalk module itself was not modified).

| Field | Value |
|---|---|
| source_id | census_geography |
| Source (dataset/API name, exact files/endpoints) | 2024 cartographic boundary files (county, state, 1:5,000,000), 2024 Gazetteer (counties, ZCTAs), 2020 cartographic county boundaries (Connecticut legacy counties only), 2020 ZCTA5–county relationship file |
| Publishing organization | U.S. Census Bureau (Geography Division) |
| Retrieval date (UTC) | {retrieved} |
| Source version / release | {'; '.join(f"{f} (Last-Modified {v['last_modified']})" for f, v in s['files'].items())} |
| Source update date / cadence | Boundaries and Gazetteer annually; relationship files per decennial census |
| License / access conditions | Public domain (U.S. federal work); open download, no registration |
| Unit of observation | State, county / county-equivalent, ZCTA, ZCTA×county part |
| Sample size (actual, as ingested) | {s['states']} states/state-equivalents; {s['counties_2024']} counties (2024) + {len(s['ct_legacy_counties'])} legacy CT counties; {s['zcta_rows']:,} ZCTAs; {s['crosswalk_rows']:,} ZCTA×county rows |
| Geography (resolution, vintage) | County 2024 (Connecticut = 9 planning regions); legacy CT counties from 2020; ZCTA 2020 definitions (2024 Gazetteer); relationship file on 2020 counties |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — reference geography |

## Files / endpoints retrieved

| File | URL | Bytes | sha256 | Last-Modified |
|---|---|---|---|---|
{files}

## Key variables

* `geographies`: geo_id (2-digit state / 5-digit county FIPS, strings), geo_level, name, state_fips, state_abbr,
  state_name, county_fips, lat/lon (Gazetteer internal point; representative point of the 5m polygon when the
  Gazetteer lacks the county), aland_m2 (land area, m²), vintage, ct_legacy.
* `zcta_centroids`: zcta, lat/lon (Gazetteer internal point), aland_m2, county_fips/state_fips/state_abbr of the
  2024 county polygon that contains the internal point.
* `zcta_county_crosswalk`: zcta, county_fips_2020, land_share_of_zcta (part land area / ZCTA land area — a land
  share, not a population share).
* `county_boundaries_2024.geoparquet`, `state_boundaries_2024.geoparquet`: 5m generalized polygons for maps
  (CRS EPSG:{s['boundaries_crs']}).

## Missingness

* `geographies` ({s['geographies_rows']} rows): duplicate geo_id = {s['duplicate_geo_id']}; null counts
  {s['missing_geographies']} (county_fips is null by design on state rows); land area ≤ 0 = {s['aland_nonpositive']}.
* Counties without a 2024 Gazetteer internal point ({len(s['counties_without_gazetteer_point'])}; lat/lon from the
  polygon representative point instead): {', '.join(s['counties_without_gazetteer_point'])} (the island areas
  AS, GU, MP, VI).
* `zcta_centroids`: {s['zcta_rows']:,} ZCTAs, duplicates {s['zcta_duplicates']}; {s['zcta_without_county']} ZCTAs have
  no county because their internal point falls outside every generalized (1:5m) county polygon — coastal or
  island ZCTAs: {', '.join(s['zcta_without_county_list'])}. `zip_to_geo` returns no county for these.
  state_fips/state_abbr are null for {s['zcta_without_state']} ZCTAs ({s['zcta_without_state_but_with_county']} of them
  with a county), i.e. the same point-in-polygon miss also drops the state (e.g. 60611 Chicago, 33109 Miami Beach,
  04108 Portland ME), although the relationship file places every one of them in a county.
* `zcta_county_crosswalk`: land_share missing = {s['crosswalk_land_share_missing']}; per-ZCTA land shares sum to
  {s['crosswalk_land_share_sum_min_max'][0]}–{s['crosswalk_land_share_sum_min_max'][1]}.

## Linkage strategy

* County sources join on the 5-digit FIPS string. Canonical set = {s['counties_2024']} counties of the 2024 vintage:
  {s['counties_by_group']}.
* Connecticut: canonical rows are the 9 planning regions {', '.join(s['ct_planning_regions'])}; the 8 legacy counties
  {', '.join(s['ct_legacy_counties'])} are kept with `ct_legacy=True` and `vintage='2020_ct_legacy'`. No values are
  reassigned between the two sets.
* ZIP → ZCTA → county: `zip_to_geo` treats a ZIP as the ZCTA of the same code and uses the 2024 county containing
  the ZCTA internal point. For non-CT ZCTAs with a county, that county equals the county holding the largest land
  share in the 2020 relationship file for {100 * s['zcta_internal_point_vs_max_land_county_agreement_non_ct']:.2f}% of
  {s['zcta_internal_point_vs_max_land_county_n']:,} ZCTAs; {s['crosswalk_zctas_multi_county']:,} of
  {s['crosswalk_zctas']:,} ZCTAs span more than one county, so a single county per ZCTA is an approximation.
* The relationship file is on 2020 counties: its {s['crosswalk_counties_2020']} counties include
  {len(s['crosswalk_counties_not_in_2024'])} not in the 2024 set ({', '.join(s['crosswalk_counties_not_in_2024'])}, legacy CT) and
  lack {len(s['counties_2024_not_in_crosswalk'])} 2024 counties ({', '.join(s['counties_2024_not_in_crosswalk'])}: the CT planning
  regions plus county-equivalents with no ZCTA).

## Limitations and caveats

* ZIP codes are USPS delivery routes, not areas; PO-box and unique ZIPs have no ZCTA and stay ungeocoded. A ZCTA
  internal point approximates a location; it is not an address geocode.
* 1:5m generalized polygons move coastlines; this is why {s['zcta_without_county']} coastal ZCTA points fall outside
  every county polygon. A fallback to the relationship file (largest land share) would resolve them but is on 2020
  counties and is not applied by the crosswalk module.
* Land-area shares are not population shares; apportioning counts from ZCTA to county by land share is a rough
  approximation, especially in rural ZCTAs.
* County vintages differ across federal sources (CT planning regions from 2022 onward; legacy CT counties in older
  products); unmatched rows must be reported, not forced.

## Processed outputs

| Table | Rows |
|---|---|
| geographies | {s['geographies_rows']} ({s['states']} state + {s['counties_2024']} county 2024 + {len(s['ct_legacy_counties'])} legacy CT) |
| zcta_centroids | {s['zcta_rows']} |
| zcta_county_crosswalk | {s['crosswalk_rows']} |
| county_boundaries_2024.geoparquet | {s['county_boundaries_rows']} (valid geometries {100 * s['county_boundaries_valid_share']:.1f}%) |
| state_boundaries_2024.geoparquet | {s['state_boundaries_rows']} |

## Reproduce

Build: `uv run python -m measure_it.geography.crosswalk`.
Re-measure this audit: `uv run python -m measure_it.ingestion.census_acs geography-audit`.
"""
    p = RAW / "census_geography" / "DATA_AUDIT.md"
    p.write_text(text)
    (raw_dir(SOURCE_ID) / "census_geography_audit_stats.json").write_text(json.dumps(s, indent=2, default=str))
    return p


def run() -> dict:
    stats = build()
    write_census_geography_audit()
    return stats


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "fetch":
        print(fetch())
    elif len(sys.argv) > 1 and sys.argv[1] == "geography-audit":
        print(write_census_geography_audit())
    else:
        st = run()
        print(json.dumps({"rows": st["rows"], "county_coverage": {k: st["county_coverage"][k] for k in
                          ("matched", "source_unmatched", "canonical_without_source_row", "connecticut_vintage")}},
                         indent=2))
