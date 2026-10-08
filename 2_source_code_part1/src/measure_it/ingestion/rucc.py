"""USDA ERS Rural-Urban Continuum Codes (2023) -> geo_context__rucc.

One row per county / county-equivalent: the 2023 RUCC (1-9), ERS's own
description, and flags derived from the ERS documentation sheet:

  metro_status                "metro" (codes 1-3) or "nonmetro" (codes 4-9),
                              read from the documentation sheet's
                              "Metropolitan counties" / "Nonmetropolitan
                              counties" sections, and cross-checked against
                              the "Metro - " / "Nonmetro - " prefix of every
                              row's Description.
  nonmetro_adjacent_to_metro  True / False for nonmetro codes (from the code
                              description text), null for metro codes.

Metro areas follow the OMB July 2023 delineation (Bulletin 23-01); populations
are 2020 Census DHC counts. Connecticut is on the 9 planning regions.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from ..config import raw_dir
from ..download import load_manifest, download_file
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from .census_acs import county_join_coverage

SOURCE_ID = "usda_rucc"
LANDING_URL = "https://www.ers.usda.gov/data-products/rural-urban-continuum-codes"
DOC_URL = "https://www.ers.usda.gov/data-products/rural-urban-continuum-codes/documentation"
URLS = {
    "2023-rural-urban-continuum-codes.csv": "https://www.ers.usda.gov/media/5768/2023-rural-urban-continuum-codes.csv",
    "2023-rural-urban-continuum-codes.xlsx": "https://www.ers.usda.gov/media/5767/2023-rural-urban-continuum-codes.xlsx",
}


def fetch() -> dict[str, Path]:
    return {f: download_file(u, SOURCE_ID, f, min_bytes=10_000) for f, u in URLS.items()}


def parse_code_definitions(xlsx: Path) -> pd.DataFrame:
    """Code -> (metro_status, definition) from the ERS 'Documentation' sheet."""
    doc = pd.read_excel(xlsx, sheet_name="Documentation", header=None, dtype=str)
    section, rows = None, []
    for a, b in zip(doc[0].fillna(""), doc[1].fillna("")):
        a, b = a.strip(), b.strip()
        if a.lower().startswith("metropolitan counties"):
            section = "metro"
        elif a.lower().startswith("nonmetropolitan counties"):
            section = "nonmetro"
        elif section and re.fullmatch(r"\d", a) and b:
            rows.append({"rucc_2023": int(a), "metro_status": section, "code_definition": b})
        elif a.lower().startswith("note"):
            section = None
    out = pd.DataFrame(rows)
    if sorted(out["rucc_2023"]) != list(range(1, 10)):
        raise ValueError(f"documentation sheet parsed to codes {sorted(out['rucc_2023'])}, expected 1-9")
    return out


def adjacency_from_definition(defn: str, metro_status: str) -> bool | None:
    if metro_status == "metro":
        return None
    d = defn.lower()
    if "not adjacent to a metro area" in d:
        return False
    if "adjacent to a metro area" in d:
        return True
    raise ValueError(f"cannot read adjacency from {defn!r}")


def load_wide(csv: Path) -> pd.DataFrame:
    long = pd.read_csv(csv, dtype=str, encoding="latin-1", keep_default_na=False)
    expected = {"FIPS", "State", "County_Name", "Attribute", "Value"}
    if set(long.columns) != expected:
        raise ValueError(f"unexpected RUCC columns {list(long.columns)}")
    wide = long.pivot(index=["FIPS", "State", "County_Name"], columns="Attribute", values="Value").reset_index()
    wide.columns.name = None
    wide["FIPS"] = wide["FIPS"].str.zfill(5)
    return wide


def transform(wide: pd.DataFrame, codes: pd.DataFrame) -> pd.DataFrame:
    df = pd.DataFrame({
        "geo_id": wide["FIPS"],
        "geo_level": "county",
        "state_fips": wide["FIPS"].str[:2],
        "state_abbr": wide["State"],
        "county_name": wide["County_Name"],
        "population_2020": pd.to_numeric(wide["Population_2020"].replace("", None)).astype("Int64"),
        "rucc_2023": pd.to_numeric(wide.get("RUCC_2023", pd.Series(index=wide.index, dtype=str)).replace("", None),
                                   errors="coerce").astype("Int64"),
        "rucc_description": wide["Description"],
    })
    df = df.merge(codes, on="rucc_2023", how="left")
    # Cross-check: ERS row descriptions carry "Metro - " / "Nonmetro - " prefixes.
    prefix = df["rucc_description"].str.extract(r"^(Metro|Nonmetro) - ")[0].str.lower()
    coded = df["rucc_2023"].notna()
    bad = coded & (prefix != df["metro_status"])
    if bad.any():
        raise ValueError(f"metro flag disagrees with description for {df.loc[bad, 'geo_id'].tolist()}")
    df["is_metro"] = df["metro_status"].map({"metro": True, "nonmetro": False}).astype("boolean")
    df["nonmetro_adjacent_to_metro"] = pd.array(
        [adjacency_from_definition(d, m) if isinstance(m, str) else None
         for d, m in zip(df["code_definition"], df["metro_status"])], dtype="boolean")
    return df.drop(columns=["code_definition"]).sort_values("geo_id").reset_index(drop=True)


def build() -> dict:
    files = fetch()
    manifest = load_manifest(SOURCE_ID)["files"]
    codes = parse_code_definitions(files["2023-rural-urban-continuum-codes.xlsx"])
    df = transform(load_wide(files["2023-rural-urban-continuum-codes.csv"]), codes)
    retrieved = manifest["2023-rural-urban-continuum-codes.csv"]["retrieved_at"]
    last_mod = manifest["2023-rural-urban-continuum-codes.csv"].get("last_modified") or "not sent by server"
    out = add_provenance(
        df, data_layer="geographic", source_name="USDA ERS Rural-Urban Continuum Codes 2023",
        source_version=f"RUCC 2023 (released Jan 2024; OMB July 2023 delineation, 2020 Census); CSV Last-Modified {last_mod}",
        retrieved_at=retrieved, evidence_type="composite_index", source_record_id="geo_id",
        source_geographic_resolution="county",
        provenance_notes=("County classification (metro size / urban population / metro adjacency). "
                          "metro_status and nonmetro_adjacent_to_metro derived from the ERS documentation sheet. "
                          "A rurality classification, not a health or burden measure."),
    )
    write_table(out, "geo_context__rucc", producer="ingestion.rucc.build",
                description="USDA ERS 2023 Rural-Urban Continuum Codes by county with metro/nonmetro flag")
    stats = _stats(df, codes, manifest)
    (raw_dir(SOURCE_ID) / "build_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    _registry(stats, manifest, retrieved)
    write_audit(stats)
    return stats


def _stats(df: pd.DataFrame, codes: pd.DataFrame, manifest: dict) -> dict:
    missing = df[df["rucc_2023"].isna()]
    by_code = (df.groupby("rucc_2023", dropna=False)
                 .agg(n=("geo_id", "size"), population_2020=("population_2020", "sum"))
                 .reset_index())
    return {
        "files": {f: {k: manifest[f].get(k) for k in ("url", "bytes", "sha256", "retrieved_at", "last_modified")}
                  for f in URLS},
        "rows": int(len(df)),
        "states_and_territories": int(df["state_abbr"].nunique()),
        "by_state_territory_rows": df["state_abbr"].value_counts().reindex(["PR", "AS", "GU", "MP", "VI"]).fillna(0)
                                                   .astype(int).to_dict(),
        "rucc_missing": int(df["rucc_2023"].isna().sum()),
        "rucc_missing_rows": missing[["geo_id", "state_abbr", "county_name", "population_2020",
                                      "rucc_description"]].to_dict("records"),
        "population_missing": int(df["population_2020"].isna().sum()),
        "by_code": [{"rucc_2023": (None if pd.isna(r.rucc_2023) else int(r.rucc_2023)), "n": int(r.n),
                     "population_2020": int(r.population_2020)} for r in by_code.itertuples()],
        "metro_counts": df["metro_status"].fillna("unclassified").value_counts().to_dict(),
        "metro_population_2020": df.groupby("metro_status")["population_2020"].sum().astype(int).to_dict(),
        "adjacent_counts": df["nonmetro_adjacent_to_metro"].astype(object)
                              .map({True: "adjacent", False: "not adjacent"})
                              .fillna("not applicable (metro or unclassified)").value_counts().to_dict(),
        "code_definitions": codes.to_dict("records"),
        "coverage": county_join_coverage(df["geo_id"], "RUCC 2023"),
    }


def _registry(stats: dict, manifest: dict, retrieved: str) -> None:
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "USDA ERS Rural-Urban Continuum Codes 2023",
        "publisher": "USDA Economic Research Service", "landing_url": LANDING_URL,
        "access_urls": list(URLS.values()) + [DOC_URL],
        "license": "U.S. federal government work (public domain); cite USDA ERS",
        "access_conditions": "open download, no registration",
        "retrieved_at": retrieved,
        "source_version": "RUCC 2023 (released January 2024; OMB Bulletin 23-01 delineation; 2020 Census)",
        "update_date": "2023 file last updated 2024-01-22 (ERS page); data product page updated 2025-12-30; "
                       "codes revised after each decennial census",
        "data_layer": "geographic", "unit_of_observation": "county / county-equivalent",
        "sample_size": {"counties": stats["rows"], "with_rucc_code": stats["rows"] - stats["rucc_missing"]},
        "geographic_resolution": "county (2023 county set; CT planning regions)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable", "true_participant_linkage_across_modalities": False,
        "status": "ingested", "processed_outputs": ["geo_context__rucc"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md", "ingestion_module": "measure_it.ingestion.rucc",
        "limitations": ["Whole-county classification hides within-county rural/urban variation",
                        "Nonmetro Virginia independent cities were combined with their counties of origin by ERS",
                        "Rurality is context, not disease burden or access to care"],
    })


def write_audit(s: dict) -> Path:
    files = "\n".join(f"| {f} | {v['url']} | {v['bytes']:,} | `{v['sha256'][:16]}…` | {v['last_modified'] or 'not sent'} |"
                      for f, v in s["files"].items())
    defs = "\n".join(f"| {d['rucc_2023']} | {d['metro_status']} | {d['code_definition']} |"
                     for d in s["code_definitions"])
    by_code = "\n".join(f"| {r['rucc_2023'] if r['rucc_2023'] is not None else 'missing'} | {r['n']} | "
                        f"{r['population_2020']:,} |" for r in s["by_code"])
    miss = "; ".join(f"{r['geo_id']} {r['county_name']} ({r['state_abbr']}, population_2020 = "
                     f"{r['population_2020']}, description '{r['rucc_description']}')" for r in s["rucc_missing_rows"])
    cov = s["coverage"]
    text = f"""# DATA AUDIT — USDA ERS Rural-Urban Continuum Codes 2023

| Field | Value |
|---|---|
| source_id | usda_rucc |
| Source (dataset/API name, exact files/endpoints) | 2023 Rural-Urban Continuum Codes, CSV + XLSX from the ERS data product page |
| Publishing organization | USDA Economic Research Service (ERS) |
| Retrieval date (UTC) | {s['files']['2023-rural-urban-continuum-codes.csv']['retrieved_at']} |
| Source version / release | RUCC 2023 ("updated January 2024" per the XLSX Documentation sheet); OMB metro delineation of July 2023 (Bulletin 23-01); 2020 Census DHC populations |
| Source update date / cadence | HTTP Last-Modified: {s['files']['2023-rural-urban-continuum-codes.csv']['last_modified'] or 'not sent by the /media/ endpoint'}; ERS page shows the data product updated 12/30/2025 and the 2023 file "Last Updated 1/22/2024"; codes are revised after each decennial census |
| License / access conditions | Public domain (U.S. federal work); open download, no registration |
| Unit of observation | County / county-equivalent |
| Sample size (actual, as ingested) | {s['rows']} rows ({s['states_and_territories']} states/territories incl. DC); {s['rows'] - s['rucc_missing']} with a code |
| Geography (resolution, vintage) | County; Connecticut on the 9 planning regions; includes PR ({s['by_state_territory_rows']['PR']}), AS ({s['by_state_territory_rows']['AS']}), GU ({s['by_state_territory_rows']['GU']}), MP ({s['by_state_territory_rows']['MP']}), VI ({s['by_state_territory_rows']['VI']}) |
| Person-level? | no |
| Geographic? | yes |
| Omics? | no |
| Wearable? | no |
| True participant linkage across modalities? | no — county classification |

## Files / endpoints retrieved

| File | URL | Bytes | sha256 | Last-Modified |
|---|---|---|---|---|
{files}

The CSV is long format (FIPS, State, County_Name, Attribute, Value) with attributes Population_2020,
RUCC_2023, Description; it is pivoted to one row per county. The XLSX is used for its Documentation sheet
(code definitions). The CSV served at the older ERS path
`https://ers.usda.gov/sites/default/files/_laserfiche/DataFiles/53251/Ruralurbancontinuumcodes2023.csv`
had the identical sha256 (ec455ee2…) when checked on 2026-09-23 and reported Last-Modified
Mon, 03 Mar 2025 17:25:44 GMT.

## Key variables

| Column | Meaning | Source |
|---|---|---|
| geo_id | 5-digit county FIPS (string) | CSV FIPS |
| population_2020 | 2020 Census population | CSV Population_2020 |
| rucc_2023 | Rural-Urban Continuum Code 1–9 | CSV RUCC_2023 |
| rucc_description | ERS description of the code | CSV Description |
| metro_status / is_metro | metro (1–3) / nonmetro (4–9) | derived: XLSX Documentation sheet sections, cross-checked against each row's "Metro - "/"Nonmetro - " prefix (0 disagreements, enforced in code) |
| nonmetro_adjacent_to_metro | nonmetro county adjacent to a metro area | derived from the code definition text; null for metro |

Code definitions parsed from the Documentation sheet:

| Code | Group | Definition |
|---|---|---|
{defs}

## Missingness

* rucc_2023 missing: {s['rucc_missing']} row(s): {miss or 'none'}. ERS documentation: two American Samoa
  county-equivalents were not classified because they reported zero population.
* population_2020 missing: {s['population_missing']}.

Rows and 2020 population by code:

| RUCC 2023 | Counties | Population 2020 |
|---|---|---|
{by_code}

Metro/nonmetro rows: {s['metro_counts']}; population: {s['metro_population_2020']}.
Nonmetro adjacency (nonmetro_adjacent_to_metro): {s['adjacent_counts']}.

## Linkage strategy

* Joins to `geographies` on `geo_id` (5-digit FIPS). RUCC counties = {cov['source_counties']}; canonical 2024
  counties = {cov['canonical_counties_2024']}; matched = {cov['matched']}; RUCC rows not in geographies =
  {cov['source_unmatched']} {cov['source_unmatched_fips']}; canonical counties without a RUCC row =
  {cov['canonical_without_source_row']} {cov['canonical_without_source_row_fips']}.
* Connecticut: RUCC rows are on **{cov['connecticut_vintage']}** ({', '.join(cov['connecticut_source_fips'])});
  {len(cov['connecticut_matched'])} match the canonical planning regions. Legacy CT counties (09001–09015) have
  no RUCC 2023 code; legacy-CT sources cannot be given a RUCC without a crosswalk (none applied).
* Ecological join only; a county's code describes the county, not residents in any person-level dataset.

## Limitations and caveats

* Whole-county classification: large metro counties contain rural areas and vice versa.
* In Virginia, nonmetro independent cities were combined with their counties of origin when computing codes
  (ERS Documentation sheet lists the nine pairs); each keeps its own row.
* Metro status follows OMB 2023 core-based statistical area delineations, which are commuting-based, not a
  measure of health-care access.
* Rurality is geographic context, not disease prevalence or burden.

## Processed outputs

| Table | Rows |
|---|---|
| geo_context__rucc | {s['rows']} |

## Reproduce

`uv run python -m measure_it.ingestion.rucc`
"""
    p = raw_dir(SOURCE_ID) / "DATA_AUDIT.md"
    p.write_text(text)
    return p


def run() -> dict:
    return build()


if __name__ == "__main__":
    st = run()
    print(json.dumps({k: st[k] for k in ("rows", "rucc_missing", "metro_counts")}
                     | {"coverage": {k: st["coverage"][k] for k in ("matched", "source_unmatched",
                                                                     "canonical_without_source_row",
                                                                     "connecticut_vintage")}}, indent=2))
