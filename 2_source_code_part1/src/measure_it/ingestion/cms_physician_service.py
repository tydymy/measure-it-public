"""Medicare Physician & Other Practitioners - by Provider and Service (data.cms.gov) -> provider_measurement_activity.

source_id: ``cms_physician_service`` (stage F, docs/HARMONIZATION_CONTRACT.md section 5).

Which clinicians actually PERFORM a measurement? NPPES taxonomies say what a clinician calls themself; Medicare
Part B claims say what they billed. This module pulls, for every HCPCS/CPT code in ``configs/measurement_hcpcs.yaml``,
the rows of the latest "by Provider and Service" file through the data.cms.gov data API with a server-side
``HCPCS_Cd`` filter (data-API helpers below; the ~3 GB national CSV is not downloaded), and writes
one row per rendering NPI x measurement class x HCPCS code (places of service summed, listed in
``place_of_service``). A code mapped to several classes (e.g. QSART 95923 -> autonomic_testing and
small_fiber_testing) gives one row per class. Measurement ids are those of configs/measurements.yaml, so a new
measurement is a config entry, not code; bundles are resolved downstream (facilities.activity) through their members.

Guardrails (in provenance_notes and DATA_AUDIT.md):
* Medicare fee-for-service Part B only (no Medicare Advantage, Medicaid, commercial, VA); rendering NPI, one year.
* CMS suppression: a provider x code x place-of-service row is published only for >= 11 beneficiaries, so
  low-volume performers are absent (absence is not evidence of absence).
* Billing a code does not mean the clinician sees Long COVID, ME/CFS or any target condition.
* ``code_role`` = dedicated / related / proxy (config); proxy rows are labelled, never presented as the measurement.

Run: ``uv run python -m measure_it.ingestion.cms_physician_service``
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .. import http
from ..config import load_config, raw_dir, utc_now_iso
from ..download import load_manifest, record_file
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table

# ---------------------------------------------------------------- data.cms.gov data API (shared with cms_partd_prescriber)
# The full national CSVs (~3-4 GB) are not downloaded: the data API filters server side. Endpoints:
#   {API}/slug?path=<dataset path>                                   -> the current version's dataset uuid
#   {API}/dataset/<uuid>/data/stats?filter[<col>]=<value>            -> found_rows
#   {API}/dataset/<uuid>/data?filter[<col>]=<value>&size=5000&offset=k
# Every response goes through measure_it.http (disk cache; offline mode answers from the cache only). The slug
# response is cached too: the first run pins the dataset version (data/raw/<source_id>/release.json);
# refresh_release=True re-reads the current version.

API = "https://data.cms.gov/data-api/v1"
PAGE_SIZE = 5000


def discover_release(source_id: str, slug_path: str, refresh: bool = False) -> dict:
    """Current dataset version for a data.cms.gov dataset path (uuid, version label, data year, modified date)."""
    js = http.get_json(f"{API}/slug", params={"path": slug_path}, refresh=refresh)
    cur = js["data"]["current_dataset"]
    version = str(cur.get("version") or "")
    rel = {
        "dataset_name": js["data"]["name"],
        "dataset_path": slug_path,
        "landing_url": f"https://data.cms.gov{slug_path}",
        "version_uuid": cur["uuid"],
        "version_name": cur["name"],
        "version": version,
        "data_year": int(version[:4]) if version[:4].isdigit() else None,
        "last_modified_date": cur.get("last_modified_date"),
        "api_url": f"{API}/dataset/{cur['uuid']}/data",
    }
    (raw_dir(source_id) / "release.json").write_text(json.dumps(rel, indent=2))
    return rel


def filtered_count(uuid: str, column: str, value: str) -> int:
    js = http.get_json(f"{API}/dataset/{uuid}/data/stats", params={f"filter[{column}]": value})
    return int(js["found_rows"])


def fetch_filtered(uuid: str, column: str, value: str, page_size: int = PAGE_SIZE) -> tuple[pd.DataFrame, dict]:
    """All rows with ``column == value`` (exact match), paged; returns (frame, fetch log)."""
    expected = filtered_count(uuid, column, value)
    rows: list[dict] = []
    offset = 0
    pages = 0
    while offset < expected:
        page = http.get_json(f"{API}/dataset/{uuid}/data",
                             params={f"filter[{column}]": value, "size": page_size, "offset": offset})
        pages += 1
        if not page:
            break
        rows.extend(page)
        offset += page_size
    df = pd.DataFrame(rows, dtype=str)
    log = {"filter_column": column, "filter_value": value, "found_rows_reported": expected,
           "rows_fetched": int(len(df)), "pages": pages, "complete": int(len(df)) == expected}
    return df, log


def write_raw_csv(df: pd.DataFrame, source_id: str, name: str) -> Path:
    p = raw_dir(source_id) / name
    df.to_csv(p, index=False)
    return p


SOURCE_ID = "cms_physician_service"
PRODUCER = "ingestion.cms_physician_service"
TABLE = "provider_measurement_activity"
METHODOLOGY_URL = "https://data.cms.gov/sites/default/files/2025-04/MUP_PHY_RY25_202350312_Methodology_508.pdf"
DICTIONARY_URL = ("https://data.cms.gov/resources/"
                  "medicare-physician-other-practitioners-by-provider-and-service-data-dictionary")

CAVEAT = (
    "Medicare fee-for-service Part B claims by rendering NPI, one calendar year (no Medicare Advantage, Medicaid, "
    "commercial or VA care). CMS publishes a provider x HCPCS x place-of-service row only for >= 11 beneficiaries, "
    "so low-volume performers are absent. Billing a code shows the clinician performs/interprets the procedure for "
    "Medicare patients; it does NOT mean the clinician sees Long COVID, ME/CFS, POTS or any target condition. "
    "code_role 'proxy' rows are indirect signals for a measurement without a dedicated code (configs/"
    "measurement_hcpcs.yaml). tot_benes_sum_over_pos can double-count a beneficiary seen in both settings."
)

ROLE_ORDER = {"dedicated": 0, "related": 1, "proxy": 2}


def load_codes_config() -> dict:
    return load_config("measurement_hcpcs")


def code_map(cfg: dict) -> pd.DataFrame:
    """One row per (measurement_class, code) with role, confidence and rationale from the config."""
    rows = []
    for mid, spec in (cfg.get("measurements") or {}).items():
        for key, role_default in (("codes", "dedicated"), ("related_codes", "related")):
            for c in spec.get(key) or []:
                role = "proxy" if c.get("proxy") else role_default
                rows.append({"measurement_class": mid, "hcpcs_cd": str(c["code"]), "code_role": role,
                             "proxy": bool(c.get("proxy", False)), "confidence": c.get("confidence", ""),
                             "config_descriptor": c.get("descriptor", ""), "rationale": c.get("rationale", "")})
    m = pd.DataFrame(rows)
    if m.empty:
        return m
    m["_o"] = m["code_role"].map(ROLE_ORDER)
    # a code listed twice for one class keeps its strongest role
    m = m.sort_values(["measurement_class", "hcpcs_cd", "_o"]).drop_duplicates(["measurement_class", "hcpcs_cd"])
    return m.drop(columns="_o").reset_index(drop=True)


def raw_name(rel: dict) -> str:
    return f"prov_svc_{rel['data_year']}_configured_codes.csv"


def fetch(rel: dict, codes: list[str]) -> tuple[pd.DataFrame, list[dict]]:
    frames, logs = [], []
    for c in sorted(set(codes)):
        df, log = fetch_filtered(rel["version_uuid"], "HCPCS_Cd", c)
        frames.append(df)
        logs.append(log)
    frames = [f for f in frames if len(f)]
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), logs


def shape(raw: pd.DataFrame, cmap: pd.DataFrame, rel: dict) -> pd.DataFrame:
    df = raw.copy()
    for c in ("Tot_Benes", "Tot_Srvcs", "Tot_Bene_Day_Srvcs", "Avg_Mdcr_Alowd_Amt"):
        df[c] = pd.to_numeric(df[c].replace("", np.nan), errors="coerce")
    df["allowed_usd"] = df["Tot_Srvcs"] * df["Avg_Mdcr_Alowd_Amt"]
    keys = ["Rndrng_NPI", "HCPCS_Cd"]
    first = ["Rndrng_Prvdr_Ent_Cd", "Rndrng_Prvdr_Type", "Rndrng_Prvdr_Crdntls", "Rndrng_Prvdr_Zip5",
             "Rndrng_Prvdr_State_Abrvtn", "Rndrng_Prvdr_State_FIPS", "Rndrng_Prvdr_Cntry", "Rndrng_Prvdr_RUCA",
             "HCPCS_Desc", "Rndrng_Prvdr_Mdcr_Prtcptg_Ind"]
    g = df.groupby(keys, sort=True)
    agg = g.agg(tot_services=("Tot_Srvcs", "sum"), tot_benes_sum_over_pos=("Tot_Benes", "sum"),
                tot_benes_max_pos=("Tot_Benes", "max"), tot_bene_day_services=("Tot_Bene_Day_Srvcs", "sum"),
                medicare_allowed_usd=("allowed_usd", "sum"), n_pos_rows=("Tot_Srvcs", "size"),
                **{c: (c, "first") for c in first}).reset_index()
    agg["place_of_service"] = g["Place_Of_Srvc"].agg(lambda s: "|".join(sorted(set(s)))).values
    m = agg.merge(cmap, left_on="HCPCS_Cd", right_on="hcpcs_cd", how="inner")
    out = pd.DataFrame({
        "npi": m["Rndrng_NPI"].astype(str),
        "measurement_class": m["measurement_class"],
        "hcpcs_cd": m["HCPCS_Cd"], "hcpcs_desc": m["HCPCS_Desc"],
        "code_role": m["code_role"], "proxy": m["proxy"], "code_confidence": m["confidence"],
        "tot_services": m["tot_services"], "tot_benes_sum_over_pos": m["tot_benes_sum_over_pos"].astype("Int64"),
        "tot_benes_max_pos": m["tot_benes_max_pos"].astype("Int64"),
        "tot_bene_day_services": m["tot_bene_day_services"], "medicare_allowed_usd": m["medicare_allowed_usd"].round(2),
        "place_of_service": m["place_of_service"], "n_pos_rows": m["n_pos_rows"].astype(int),
        "entity_code": m["Rndrng_Prvdr_Ent_Cd"], "cms_provider_type": m["Rndrng_Prvdr_Type"],
        "cms_credentials": m["Rndrng_Prvdr_Crdntls"], "medicare_participating": m["Rndrng_Prvdr_Mdcr_Prtcptg_Ind"],
        "cms_zip5": m["Rndrng_Prvdr_Zip5"], "cms_state": m["Rndrng_Prvdr_State_Abrvtn"],
        "cms_state_fips": m["Rndrng_Prvdr_State_FIPS"], "cms_country": m["Rndrng_Prvdr_Cntry"],
        "cms_ruca": m["Rndrng_Prvdr_RUCA"],
        "data_year": rel["data_year"], "payer_scope": "Medicare FFS Part B",
        "suppression_rule": "rows with < 11 beneficiaries per NPI x HCPCS x place of service are not published",
    })
    out["object_id"] = "npi:" + out["npi"]
    out = out.sort_values(["measurement_class", "hcpcs_cd", "npi"]).reset_index(drop=True)
    return add_provenance(
        out, data_layer="facility",
        source_name="CMS Medicare Physician & Other Practitioners - by Provider and Service",
        source_version=f"{rel['version_name']} (modified {rel['last_modified_date']})",
        retrieved_at=load_manifest(SOURCE_ID)["files"][raw_name(rel)]["retrieved_at"],
        evidence_type="administrative_claims",
        source_record_id=lambda d: d["npi"] + "|" + d["hcpcs_cd"] + "|" + d["measurement_class"],
        source_geographic_resolution="zcta_centroid", evidence_level="administrative_claims_ffs",
        provenance_notes=CAVEAT)


def write_audit(rel: dict, logs: list[dict], raw: pd.DataFrame, tab: pd.DataFrame, cmap: pd.DataFrame) -> None:
    m = load_manifest(SOURCE_ID)["files"][raw_name(rel)]
    benes = pd.to_numeric(raw["Tot_Benes"], errors="coerce") if len(raw) else pd.Series(dtype=float)
    per_code = (tab.drop_duplicates(["npi", "hcpcs_cd"]).groupby("hcpcs_cd")
                .agg(npis=("npi", "nunique"), services=("tot_services", "sum"),
                     individuals=("entity_code", lambda s: int((s == "I").sum()))))
    desc = tab.drop_duplicates("hcpcs_cd").set_index("hcpcs_cd")["hcpcs_desc"]
    per_class = (tab.groupby(["measurement_class", "code_role"])
                 .agg(npis=("npi", "nunique"), services=("tot_services", "sum")).reset_index())
    found = {lg["filter_value"]: lg["found_rows_reported"] for lg in logs}
    lines = [
        "# DATA AUDIT — Medicare Physician & Other Practitioners - by Provider and Service (configured codes only)", "",
        "| Field | Value |", "|---|---|",
        f"| source_id | {SOURCE_ID} |",
        f"| Source | {rel['dataset_name']}, version `{rel['version_name']}`; data API {rel['api_url']} filtered by "
        f"`HCPCS_Cd` (exact) for the {cmap['hcpcs_cd'].nunique()} codes in configs/measurement_hcpcs.yaml |",
        "| Publishing organization | CMS Office of Enterprise Data and Analytics |",
        f"| Retrieval date (UTC) | {m['retrieved_at']} |",
        f"| Source version / release | data year {rel['data_year']} (calendar year); dataset modified {rel['last_modified_date']} |",
        "| Source update date / cadence | annual |",
        "| License / access conditions | public, no registration (CMS public use file) |",
        "| Unit of observation | source: rendering NPI x HCPCS x place of service (F facility / O office); table: NPI x "
        "measurement class x HCPCS |",
        f"| Sample size (actual, as ingested) | {len(raw):,} source rows; {len(tab):,} table rows; "
        f"{tab['npi'].nunique():,} distinct NPIs |",
        "| Geography | provider practice ZIP5 and state as reported to CMS; county assigned downstream from NPPES |",
        "| Person-level? | no (provider aggregates) |", "| Geographic? | yes |", "| Omics? | no |",
        "| Wearable? | no |", "| True participant linkage across modalities? | no (not participants) |", "",
        "## Files / endpoints retrieved", "",
        f"* `{raw_name(rel)}`: {m['bytes']:,} bytes, sha256 `{m['sha256']}` (MANIFEST.json), assembled from paged "
        "API responses (cached under data/_http_cache).",
        f"* landing page {rel['landing_url']}; data dictionary {DICTIONARY_URL}; methodology {METHODOLOGY_URL}.", "",
        "## Suppression rule (Medicare FFS only)", "",
        "CMS publishes a provider x HCPCS x place-of-service row only when it covers 11 or more beneficiaries. Measured "
        f"here: minimum `Tot_Benes` over the {len(raw):,} source rows = "
        f"{int(benes.min()) if len(benes) else 'n/a'}; rows with `Tot_Benes` < 11: {int((benes < 11).sum())}. A "
        "clinician who performed a test for 1-10 Medicare FFS beneficiaries is therefore absent, as is all care paid "
        "by Medicare Advantage (about half of Medicare enrollees), Medicaid, commercial insurers and the VA. Facility "
        "(technical-component) billing by hospital outpatient departments is in institutional claims and not in this "
        "file; the professional component billed by the interpreting clinician is (place of service F).", "",
        "## Measured counts per code", "",
        "| HCPCS | CMS descriptor (2024 file) | API found_rows | NPIs | of which individual (entity I) | services |",
        "|---|---|---|---|---|---|",
    ]
    for code in sorted(cmap["hcpcs_cd"].unique()):
        if code in per_code.index:
            r = per_code.loc[code]
            lines.append(f"| {code} | {desc.get(code, '')} | {found.get(code, 0):,} | {int(r.npis):,} | "
                         f"{int(r.individuals):,} | {int(r.services):,} |")
        else:
            lines.append(f"| {code} | (no rows in this data year) | {found.get(code, 0):,} | 0 | 0 | 0 |")
    lines += ["", "Configured codes with 0 rows are kept in the config (they may be new, deleted, rarely billed or "
              "always below the suppression threshold); see the config's notes.", "",
              "## Measured counts per measurement class and code role", "",
              "| measurement_class | code_role | NPIs | services |", "|---|---|---|---|"]
    lines += [f"| {r.measurement_class} | {r.code_role} | {r.npis:,} | {int(r.services):,} |"
              for r in per_class.itertuples()]
    no_rows = sorted(set((load_codes_config().get("measurements") or {})) - set(tab["measurement_class"]))
    lines += ["", f"Measurement classes with no configured code or no rows: {', '.join(no_rows) or 'none'}.", "",
              "## Key variables", "",
              "* `Rndrng_NPI` -> `npi`; `HCPCS_Cd`/`HCPCS_Desc`; `Tot_Srvcs` -> `tot_services` (summed over places of "
              "service); `Tot_Benes` -> `tot_benes_sum_over_pos` (may double-count across F/O) and `tot_benes_max_pos` "
              "(a lower bound); `Rndrng_Prvdr_Ent_Cd` I individual / O organisation; `Rndrng_Prvdr_Type` CMS "
              "specialty from claims; `Rndrng_Prvdr_Zip5`, state; `Avg_Mdcr_Alowd_Amt` x services -> "
              "`medicare_allowed_usd`.",
              "* `measurement_class`, `code_role` (dedicated/related/proxy), `code_confidence` come from the config.", "",
              "## Missingness", "",
              f"* `cms_zip5` empty: {int((tab['cms_zip5'].fillna('') == '').sum()):,} rows; non-US country: "
              f"{int((tab['cms_country'].fillna('US') != 'US').sum()):,} rows.",
              f"* `tot_benes_sum_over_pos` null: {int(tab['tot_benes_sum_over_pos'].isna().sum()):,} rows.", "",
              "## Linkage strategy", "",
              "`npi` joins NPPES `providers.npi` (facilities.activity) for county, specialty group and contact fields; "
              "NPIs not in `providers` fall back to `cms_zip5` -> ZCTA -> 2024 county (labelled). Nothing links to "
              "people.", "",
              "## Limitations and caveats", "",
              "* Medicare FFS only, one year, >= 11 beneficiaries per row: a lower bound on who performs a test.",
              "* Billing is not diagnosis: a cardiologist billing tilt-table tests is not thereby a Long COVID clinician.",
              "* No CPT code exists for nailfold capillaroscopy (config); its rows are labelled proxies.",
              "* Several codes are shared across classes (RPM device supply/management is not specific to heart rate, "
              "activity or any sensor).", "",
              "## Processed outputs", "", f"* `{TABLE}`: {len(tab):,} rows.", "", "## Reproduce", "",
              "`uv run python -m measure_it.ingestion.cms_physician_service` (add `--refresh-release` to re-read the "
              "current CMS version)."]
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text("\n".join(lines) + "\n")


def write_registry(rel: dict, raw: pd.DataFrame, tab: pd.DataFrame) -> None:
    m = load_manifest(SOURCE_ID)["files"][raw_name(rel)]
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "Medicare Physician & Other Practitioners - by Provider and Service (configured HCPCS codes)",
        "publisher": "CMS (Office of Enterprise Data and Analytics)", "landing_url": rel["landing_url"],
        "access_urls": [rel["api_url"]], "license": "U.S. government public use file",
        "access_conditions": "open API, no registration", "retrieved_at": m["retrieved_at"],
        "source_version": rel["version_name"], "update_date": rel["last_modified_date"],
        "data_layer": "facility", "unit_of_observation": "rendering NPI x HCPCS x place of service (calendar year)",
        "sample_size": {"source_rows": int(len(raw)), "table_rows": int(len(tab)),
                        "npis": int(tab["npi"].nunique()),
                        "npis_per_measurement_class": {k: int(v) for k, v in
                                                       tab.groupby("measurement_class")["npi"].nunique().items()}},
        "geographic_resolution": "provider ZIP5 (county via NPPES join downstream)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (providers, not participants)",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": [TABLE], "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.cms_physician_service",
        "limitations": ["Medicare fee-for-service Part B only", "rows with < 11 beneficiaries suppressed by CMS",
                        "billing a code does not mean the clinician sees the target condition",
                        "no CPT code exists for nailfold capillaroscopy; proxies are labelled"],
        "data_year": rel["data_year"],
    })


def run(refresh_release: bool = False) -> dict:
    cfg = load_codes_config()
    rel = discover_release(SOURCE_ID, cfg["source"]["slug_path"], refresh=refresh_release)
    cmap = code_map(cfg)
    raw, logs = fetch(rel, cmap["hcpcs_cd"].tolist())
    incomplete = [lg for lg in logs if not lg["complete"]]
    if incomplete:
        raise RuntimeError(f"incomplete API pulls: {incomplete}")
    name = raw_name(rel)
    write_raw_csv(raw, SOURCE_ID, name)
    record_file(SOURCE_ID, name, url=rel["api_url"],
                note="rows of the configured HCPCS_Cd values, assembled from paged data-API responses",
                filters=[lg["filter_value"] for lg in logs], data_year=rel["data_year"])
    (raw_dir(SOURCE_ID) / "fetch_log.json").write_text(json.dumps({"fetched_at": utc_now_iso(), "pulls": logs},
                                                                   indent=2))
    tab = shape(raw, cmap, rel)
    write_table(tab, TABLE, producer=PRODUCER,
                description=f"Medicare FFS rendering NPIs billing configured measurement codes, {rel['data_year']}")
    write_audit(rel, logs, raw, tab, cmap)
    write_registry(rel, raw, tab)
    return {"data_year": rel["data_year"], "source_rows": int(len(raw)), "rows": int(len(tab)),
            "npis": int(tab["npi"].nunique()),
            "npis_per_class": tab.groupby("measurement_class")["npi"].nunique().to_dict()}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Medicare FFS providers billing configured measurement codes")
    ap.add_argument("--refresh-release", action="store_true", help="re-read the current CMS dataset version")
    a = ap.parse_args(argv)
    print(json.dumps(run(refresh_release=a.refresh_release), indent=2, default=str))


if __name__ == "__main__":
    main()
