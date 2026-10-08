"""Medicare Part D Prescribers - by Provider and Drug (data.cms.gov) -> provider_condition_rx_signals.

source_id: ``cms_partd_prescriber`` (stage F, docs/HARMONIZATION_CONTRACT.md section 5).

Only the generic names configured in ``configs/condition_rx_signals.yaml`` are pulled, through the data.cms.gov data
API with a server-side ``Gnrc_Name`` filter (helpers in ``measure_it.ingestion.cms_physician_service``); the ~4 GB national CSV is not
downloaded. Output: one row per prescriber NPI x signal x drug (generic x brand as CMS publishes it), with claims,
30-day fills, beneficiaries (blank in the source when 1-10) and the prescriber type CMS derived from claims.

Guardrails (also in the table's provenance_notes and DATA_AUDIT.md):
* A prescription is not a diagnosis; no configured drug is specific to Long COVID, ME/CFS or POTS (see the config).
* Medicare Part D enrollees only (PDP and MA-PD plans); rows with fewer than 11 claims are not published by CMS.

Run: ``uv run python -m measure_it.ingestion.cms_partd_prescriber``
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from ..config import load_config, raw_dir, utc_now_iso
from ..download import load_manifest, record_file
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from . import cms_physician_service as api

SOURCE_ID = "cms_partd_prescriber"
PRODUCER = "ingestion.cms_partd_prescriber"
TABLE = "provider_condition_rx_signals"
METHODOLOGY_URL = "https://data.cms.gov/sites/default/files/2026-05/MUP_DPR_RY26_20260421_Methodology_508.pdf"
DICTIONARY_URL = "https://data.cms.gov/resources/medicare-part-d-prescribers-by-provider-and-drug-data-dictionary"

CAVEAT = (
    "Medicare Part D prescribing by the prescriber NPI, one calendar year. A prescription is NOT a diagnosis: none of "
    "the configured drugs is specific to POTS, Long COVID or ME/CFS (configs/condition_rx_signals.yaml gives each "
    "drug's specificity). Part D enrollees only (stand-alone PDP and Medicare Advantage drug plans; mostly >= 65 or "
    "disabled). CMS publishes a prescriber x drug row "
    "only with >= 11 claims; Tot_Benes is blank when 1-10 beneficiaries. Absence of a row is not absence of "
    "prescribing."
)

NUM = ["Tot_Clms", "Tot_30day_Fills", "Tot_Day_Suply", "Tot_Drug_Cst", "Tot_Benes"]


def load_signals() -> dict:
    return load_config("condition_rx_signals")


def drug_rows(cfg: dict) -> pd.DataFrame:
    """One row per (signal, drug) from the config, with each condition's relation."""
    out = []
    for sid, s in cfg["signals"].items():
        for d in s["drugs"]:
            out.append({"signal_id": sid, "signal_label": s["label"], "gnrc_name": d["gnrc_name"],
                        "ingredient": d["ingredient"], "specificity": d["specificity"],
                        "conditions": "|".join(sorted(s["relation"])),
                        "condition_relation": json.dumps(s["relation"], sort_keys=True)})
    return pd.DataFrame(out)


def shape(raw: pd.DataFrame, drugs: pd.DataFrame, rel: dict) -> pd.DataFrame:
    """Raw API rows -> provider_condition_rx_signals rows (one per NPI x signal x generic x brand)."""
    df = raw.copy()
    for c in NUM:
        df[c] = pd.to_numeric(df[c].replace("", np.nan), errors="coerce")
    df = df.merge(drugs, left_on="Gnrc_Name", right_on="gnrc_name", how="inner")
    out = pd.DataFrame({
        "npi": df["Prscrbr_NPI"].astype(str),
        "signal_id": df["signal_id"], "signal_label": df["signal_label"],
        "conditions": df["conditions"], "condition_relation": df["condition_relation"],
        "drug_generic": df["Gnrc_Name"], "drug_brand": df["Brnd_Name"], "ingredient": df["ingredient"],
        "specificity": df["specificity"],
        "tot_claims": df["Tot_Clms"].astype("Int64"), "tot_30day_fills": df["Tot_30day_Fills"],
        "tot_day_supply": df["Tot_Day_Suply"], "tot_drug_cost_usd": df["Tot_Drug_Cst"],
        "tot_benes": df["Tot_Benes"].astype("Int64"),
        "benes_suppressed": df["Tot_Benes"].isna(),
        "prescriber_type": df["Prscrbr_Type"], "prescriber_type_source": df["Prscrbr_Type_Src"],
        "prescriber_state": df["Prscrbr_State_Abrvtn"], "prescriber_state_fips": df["Prscrbr_State_FIPS"],
        "data_year": rel["data_year"], "payer_scope": "Medicare Part D",
    })
    out["object_id"] = "npi:" + out["npi"]
    out = out.sort_values(["signal_id", "drug_generic", "npi", "drug_brand"]).reset_index(drop=True)
    return add_provenance(
        out, data_layer="facility", source_name="CMS Medicare Part D Prescribers - by Provider and Drug",
        source_version=f"{rel['version_name']} (modified {rel['last_modified_date']})",
        retrieved_at=load_manifest(SOURCE_ID)["files"][raw_name(rel)]["retrieved_at"],
        evidence_type="administrative_claims",
        source_record_id=lambda d: d["npi"] + "|" + d["drug_generic"] + "|" + d["drug_brand"],
        source_geographic_resolution="state", evidence_level="administrative_claims_partd",
        provenance_notes=CAVEAT)


def raw_name(rel: dict) -> str:
    return f"partd_prescriber_drug_{rel['data_year']}_configured_drugs.csv"


def fetch(rel: dict, drugs: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    frames, logs = [], []
    for g in sorted(drugs["gnrc_name"].unique()):
        df, log = api.fetch_filtered(rel["version_uuid"], "Gnrc_Name", g)
        frames.append(df)
        logs.append(log)
    raw = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return raw, logs


def write_audit(rel: dict, logs: list[dict], tab: pd.DataFrame, drugs: pd.DataFrame) -> None:
    by_drug = tab.groupby(["signal_id", "drug_generic"]).agg(
        rows=("npi", "size"), prescribers=("npi", "nunique"), claims=("tot_claims", "sum"),
        min_claims=("tot_claims", "min"), benes_suppressed=("benes_suppressed", "sum")).reset_index()
    types = (tab.groupby(["drug_generic", "prescriber_type"]).size().rename("n").reset_index()
             .sort_values(["drug_generic", "n"], ascending=[True, False]).groupby("drug_generic").head(5))
    m = load_manifest(SOURCE_ID)["files"][raw_name(rel)]
    lines = [
        "# DATA AUDIT — Medicare Part D Prescribers - by Provider and Drug (configured drugs only)", "",
        "| Field | Value |", "|---|---|",
        f"| source_id | {SOURCE_ID} |",
        f"| Source | {rel['dataset_name']}, version `{rel['version_name']}`; data API {rel['api_url']} filtered by "
        f"`Gnrc_Name` (exact) for the {drugs['gnrc_name'].nunique()} generic names in configs/condition_rx_signals.yaml |",
        "| Publishing organization | CMS Office of Enterprise Data and Analytics |",
        f"| Retrieval date (UTC) | {m['retrieved_at']} |",
        f"| Source version / release | data year {rel['data_year']} (calendar year); dataset modified {rel['last_modified_date']} |",
        "| Source update date / cadence | annual |",
        "| License / access conditions | public, no registration (CMS public use file) |",
        "| Unit of observation | prescriber NPI x brand name x generic name (one calendar year) |",
        f"| Sample size (actual, as ingested) | {len(tab):,} rows, {tab['npi'].nunique():,} distinct prescriber NPIs |",
        "| Geography | prescriber state only in this table (county comes from the NPPES join in facilities.activity) |",
        "| Person-level? | no (prescriber aggregates) |", "| Geographic? | yes (prescriber state) |",
        "| Omics? | no |", "| Wearable? | no |",
        "| True participant linkage across modalities? | no (not participants) |", "",
        "## Files / endpoints retrieved", "",
        f"* `{raw_name(rel)}`: {m['bytes']:,} bytes, sha256 `{m['sha256']}` (MANIFEST.json), assembled from paged "
        "API responses (cached under data/_http_cache).",
        f"* landing page {rel['landing_url']}; data dictionary {DICTIONARY_URL}; methodology {METHODOLOGY_URL}.", "",
        "Per filter value (API `found_rows` vs rows fetched):", "",
        "| Gnrc_Name | found_rows | fetched | pages | complete |", "|---|---|---|---|---|",
    ]
    lines += [f"| {lg['filter_value']} | {lg['found_rows_reported']:,} | {lg['rows_fetched']:,} | {lg['pages']} | "
              f"{lg['complete']} |" for lg in logs]
    lines += ["", "## Key variables", "",
              "* `Prscrbr_NPI` -> `npi`; `Gnrc_Name`/`Brnd_Name` -> `drug_generic`/`drug_brand`; `Tot_Clms` -> "
              "`tot_claims` (claims incl. refills); `Tot_30day_Fills`; `Tot_Benes` -> `tot_benes` (blank in source when "
              "1-10 -> null, `benes_suppressed`); `Prscrbr_Type` (CMS-derived specialty: from claims or NPPES, "
              "`Prscrbr_Type_Src`).",
              "* `signal_id`, `conditions`, `condition_relation`, `specificity` come from the config, not from CMS.",
              "", "## Measured counts per drug", "",
              "| signal | Gnrc_Name | rows | prescriber NPIs | claims | min claims in a row | rows with beneficiaries suppressed |",
              "|---|---|---|---|---|---|---|"]
    lines += [f"| {r.signal_id} | {r.drug_generic} | {r.rows:,} | {r.prescribers:,} | {int(r.claims):,} | "
              f"{int(r.min_claims)} | {int(r.benes_suppressed):,} |" for r in by_drug.itertuples()]
    lines += ["", f"Minimum `Tot_Clms` across all rows: {int(tab['tot_claims'].min())} (the CMS rule: a prescriber x "
              "drug row is published only with 11 or more claims).", "",
              "Top prescriber types per drug (CMS `Prscrbr_Type`):", "", "| Gnrc_Name | Prscrbr_Type | rows |",
              "|---|---|---|"]
    lines += [f"| {r.drug_generic} | {r.prescriber_type} | {r.n:,} |" for r in types.itertuples()]
    lines += ["", "## Missingness", "",
              f"* `tot_benes` null (suppressed 1-10): {int(tab['benes_suppressed'].sum()):,} of {len(tab):,} rows "
              f"({tab['benes_suppressed'].mean():.1%}).",
              f"* `prescriber_type` empty: {int((tab['prescriber_type'].fillna('') == '').sum()):,} rows.", "",
              "## Linkage strategy", "",
              "`npi` joins NPPES `providers.npi` (facilities.activity), which supplies county, specialty and contact "
              "fields. Nothing here links to people.", "",
              "## Limitations and caveats", "",
              "* A prescription is not a diagnosis; no configured drug is specific to POTS, Long COVID or ME/CFS. "
              "Midodrine is mostly dialysis/orthostatic hypotension in older adults, pyridostigmine myasthenia gravis, "
              "ivabradine heart failure, fludrocortisone adrenal insufficiency; `Cromolyn Sodium` mixes ophthalmic and "
              "oral forms (brand `Gastrocrom`, the unambiguous oral product, has 5 rows in 2024).",
              "* Part D enrollees only (stand-alone PDP and Medicare Advantage drug plans, per the CMS methodology; mostly "
              ">= 65 or disabled); POTS is mainly a condition of younger adults.",
              "* CMS suppression: rows with < 11 claims are absent; beneficiary counts 1-10 are blank. Low-volume "
              "prescribers are invisible, so absence is not evidence of absence.",
              "* Long COVID, ME/CFS, Lyme disease and PTLDS have no specific drug signal (the config says why).", "",
              "## Processed outputs", "",
              f"* `{TABLE}`: {len(tab):,} rows.", "", "## Reproduce", "",
              "`uv run python -m measure_it.ingestion.cms_partd_prescriber` (add `--refresh-release` to re-read the "
              "current CMS version; otherwise the cached version is reused)."]
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text("\n".join(lines) + "\n")


def write_registry(rel: dict, tab: pd.DataFrame) -> None:
    m = load_manifest(SOURCE_ID)["files"][raw_name(rel)]
    write_registry_entry({
        "source_id": SOURCE_ID, "name": "Medicare Part D Prescribers - by Provider and Drug (configured drugs)",
        "publisher": "CMS (Office of Enterprise Data and Analytics)", "landing_url": rel["landing_url"],
        "access_urls": [rel["api_url"]], "license": "U.S. government public use file",
        "access_conditions": "open API, no registration", "retrieved_at": m["retrieved_at"],
        "source_version": rel["version_name"], "update_date": rel["last_modified_date"],
        "data_layer": "facility", "unit_of_observation": "prescriber NPI x drug (calendar year)",
        "sample_size": {"rows": int(len(tab)), "prescriber_npis": int(tab["npi"].nunique()),
                        "per_drug_rows": {k: int(v) for k, v in tab["drug_generic"].value_counts().items()}},
        "geographic_resolution": "prescriber state (county via NPPES join downstream)",
        "person_level": False, "geographic": True, "omics": False, "wearable": False,
        "participant_linkage": "not applicable (prescribers, not participants)",
        "true_participant_linkage_across_modalities": False, "status": "ingested",
        "processed_outputs": [TABLE], "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": "measure_it.ingestion.cms_partd_prescriber",
        "limitations": ["a prescription is not a diagnosis; configured drugs are not specific to the target conditions",
                        "Medicare Part D enrollees only (PDP and MA-PD)",
                        "rows with < 11 claims are suppressed by CMS; beneficiary counts 1-10 blank",
                        "generic name does not distinguish dosage forms (cromolyn ophthalmic vs oral)"],
        "data_year": rel["data_year"],
    })


def run(refresh_release: bool = False) -> dict:
    cfg = load_signals()
    rel = api.discover_release(SOURCE_ID, cfg["source"]["slug_path"], refresh=refresh_release)
    drugs = drug_rows(cfg)
    raw, logs = fetch(rel, drugs)
    incomplete = [lg for lg in logs if not lg["complete"]]
    if incomplete:
        raise RuntimeError(f"incomplete API pulls: {incomplete}")
    name = raw_name(rel)
    api.write_raw_csv(raw, SOURCE_ID, name)
    record_file(SOURCE_ID, name, url=rel["api_url"],
                note="rows of the configured Gnrc_Name values, assembled from paged data-API responses",
                filters=[lg["filter_value"] for lg in logs], data_year=rel["data_year"])
    (raw_dir(SOURCE_ID) / "fetch_log.json").write_text(json.dumps({"fetched_at": utc_now_iso(), "pulls": logs},
                                                                   indent=2))
    tab = shape(raw, drugs, rel)
    write_table(tab, TABLE, producer=PRODUCER,
                description=f"Medicare Part D prescribers of configured condition-signal drugs, {rel['data_year']}")
    write_audit(rel, logs, tab, drugs)
    write_registry(rel, tab)
    return {"data_year": rel["data_year"], "rows": int(len(tab)), "npis": int(tab["npi"].nunique()),
            "per_drug": tab["drug_generic"].value_counts().to_dict()}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Medicare Part D prescribers of configured signal drugs")
    ap.add_argument("--refresh-release", action="store_true", help="re-read the current CMS dataset version")
    a = ap.parse_args(argv)
    print(json.dumps(run(refresh_release=a.refresh_release), indent=2, default=str))


if __name__ == "__main__":
    main()
