"""Reference data for the molecular-coherence analysis: Reactome (pinned release) and HGNC.

Two lightweight sources, each with raw files under data/raw/<source_id>/ (MANIFEST.json via download_file),
a DATA_AUDIT.md and a registry_entry.yaml:

* ``reactome`` -- Reactome pathway knowledgebase, release REACTOME_RELEASE, from the versioned archive
  https://download.reactome.org/<release>/ (never ``/download/current/``, which moves with every release).
  Files: ReactomePathways.gmt.zip (pathway -> gene symbols), ReactomePathways.txt (names, species),
  ReactomePathwaysRelation.txt (hierarchy). The release-97 archive files are byte-identical (sha256) to the
  ``current`` files the first coherence build fetched on 2026-09-23, so pinning changed no result.
  The optional Reactome AnalysisService cross-check (coherence.reactome_service_crosscheck) calls a live service
  that always runs the current release; its responses are cached in data/_http_cache.
* ``hgnc`` -- HGNC complete set (approved symbols, previous symbols, aliases, Ensembl ids, locus group). HGNC
  publishes this file without release numbers and updates it continuously; the copy used is the file in
  data/raw/hgnc/ (sha256 and Last-Modified in MANIFEST.json). It was seeded from the HTTP-cache response the first
  coherence build used (fetched 2026-09-23T23:08:00Z, Last-Modified 18 Sep 2026), so results are unchanged. A cold
  checkout without that file downloads the then-current HGNC set (a known non-reproducible input).

Both are reference vocabularies used for condition-level molecular enrichment only; neither describes a person.

    uv run python -m measure_it.omics.reference_data     # ensure files, write registry entries + audits
"""
from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd

from ..config import raw_dir, utc_now_iso
from ..download import download_file, load_manifest, record_file
from ..registry import write_registry_entry

REACTOME_SOURCE = "reactome"
HGNC_SOURCE = "hgnc"
REACTOME_RELEASE = 97
REACTOME_DL = f"https://download.reactome.org/{REACTOME_RELEASE}/"
REACTOME_FILES = {
    "gmt": "ReactomePathways.gmt.zip",
    "pathways": "ReactomePathways.txt",
    "relation": "ReactomePathwaysRelation.txt",
}
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_FILE = "hgnc_complete_set.txt"
CONSUMERS = ("measure_it.omics.coherence (results/tables/test3_*.csv, Test 3 molecular coherence) and "
             "measure_it.omics.graph (measurable_biology, fig5_*.csv)")


def reactome_url(key: str) -> str:
    return REACTOME_DL + REACTOME_FILES[key]


def _seed_hgnc_from_http_cache() -> Path | None:
    """Copy the HGNC response the first coherence build used from the HTTP cache into data/raw/hgnc (once)."""
    from .. import http

    key = http._cache_key("GET", HGNC_URL, None, None, None)
    meta_path, body_path = http._cache_paths(key)
    if not (meta_path.exists() and body_path.exists()):
        return None
    meta = json.loads(meta_path.read_text())
    if int(meta.get("status", 0)) != 200:
        return None
    dest = raw_dir(HGNC_SOURCE) / HGNC_FILE
    dest.write_bytes(body_path.read_bytes())
    hdr = {k.lower(): v for k, v in (meta.get("headers") or {}).items()}
    entry = record_file(HGNC_SOURCE, HGNC_FILE, url=HGNC_URL,
                        note="seeded from the data/_http_cache response used by the first molecular-coherence build",
                        last_modified=hdr.get("last-modified"), etag=hdr.get("etag"),
                        content_type=hdr.get("content-type"))
    # keep the original fetch time as the retrieval time of these bytes
    from ..download import _record
    entry["retrieved_at"] = meta.get("fetched_at", entry["retrieved_at"])
    entry["seeded_from_http_cache_at"] = utc_now_iso()
    _record(HGNC_SOURCE, HGNC_FILE, entry)
    return dest


def ensure_files() -> dict[str, Path]:
    """Local paths of every reference file (downloaded once; MEASURE_IT_OFFLINE=1 requires them to exist)."""
    out = {k: download_file(reactome_url(k), REACTOME_SOURCE, fn, min_bytes=10_000) for k, fn in REACTOME_FILES.items()}
    dest = raw_dir(HGNC_SOURCE) / HGNC_FILE
    if not (dest.exists() and HGNC_FILE in load_manifest(HGNC_SOURCE)["files"]):
        from ..http import offline
        if not offline():
            _seed_hgnc_from_http_cache()
    out["hgnc"] = download_file(HGNC_URL, HGNC_SOURCE, HGNC_FILE, min_bytes=1_000_000)
    return out


def file_meta(source_id: str, filename: str) -> dict:
    m = load_manifest(source_id)["files"][filename]
    return {"url": m["url"], "fetched_at": m["retrieved_at"], "last_modified": m.get("last_modified") or "",
            "sha256": m.get("sha256", ""), "bytes": m.get("bytes")}


# --------------------------------------------------------------------------- registry + audit
def _stats(paths: dict[str, Path]) -> dict:
    with zipfile.ZipFile(paths["gmt"]) as z:
        gmt = z.read(z.namelist()[0]).decode()
    gmt_rows = [ln.split("\t") for ln in gmt.splitlines() if ln.strip()]
    pw = pd.read_csv(paths["pathways"], sep="\t", header=None, names=["stid", "name", "species"], dtype=str)
    rel = pd.read_csv(paths["relation"], sep="\t", header=None, names=["parent", "child"], dtype=str)
    hs = pw[pw["species"] == "Homo sapiens"]
    hrel = rel[rel["parent"].str.startswith("R-HSA") & rel["child"].str.startswith("R-HSA")]
    hg = pd.read_csv(paths["hgnc"], sep="\t", dtype=str, usecols=["hgnc_id", "symbol", "locus_group", "status",
                                                                     "ensembl_gene_id"])
    appr = hg[hg["status"].fillna("Approved") == "Approved"]
    return {
        "reactome": {"gmt_pathways": len(gmt_rows), "gmt_pathways_hsa": sum(1 for r in gmt_rows if len(r) > 1 and
                                                                            r[1].startswith("R-HSA")),
                     "gmt_distinct_genes": len({g for r in gmt_rows for g in r[2:]}),
                     "pathways_all_species": len(pw), "pathways_homo_sapiens": len(hs),
                     "relations_all": len(rel), "relations_homo_sapiens": len(hrel),
                     "top_level_homo_sapiens": len(set(hs["stid"]) - set(hrel["child"]))},
        "hgnc": {"rows": len(hg), "approved": len(appr),
                 "approved_protein_coding": int((appr["locus_group"] == "protein-coding gene").sum()),
                 "approved_missing_ensembl_id": int(appr["ensembl_gene_id"].isna().sum()),
                 "approved_missing_ensembl_id_pct": round(100 * float(appr["ensembl_gene_id"].isna().mean()), 2)},
    }


def _files_md(source_id: str) -> str:
    rows = ["| file | url | bytes | sha256 | retrieved_at (UTC) | Last-Modified |", "|---|---|---|---|---|---|"]
    for fn, m in sorted(load_manifest(source_id)["files"].items()):
        rows.append(f"| {fn} | {m['url']} | {m.get('bytes', ''):,} | `{m.get('sha256', '')[:16]}...` | "
                    f"{m.get('retrieved_at', '')} | {m.get('last_modified') or ''} |")
    return "\n".join(rows)


def write_audits_and_registry(paths: dict[str, Path] | None = None) -> dict:
    paths = paths or ensure_files()
    st = _stats(paths)
    rm = {k: file_meta(REACTOME_SOURCE, fn) for k, fn in REACTOME_FILES.items()}
    hm = file_meta(HGNC_SOURCE, HGNC_FILE)
    r_ret = min(m["fetched_at"] for m in rm.values())
    r = st["reactome"]
    h = st["hgnc"]
    common = {"geographic_resolution": "none", "person_level": False, "geographic": False, "wearable": False,
              "participant_linkage": "not applicable (reference vocabulary)",
              "true_participant_linkage_across_modalities": False, "status": "ingested",
              "processed_outputs": [], "used_by": CONSUMERS}
    write_registry_entry({
        "source_id": REACTOME_SOURCE, "name": f"Reactome pathway knowledgebase, release {REACTOME_RELEASE} (GMT, pathway list, hierarchy)",
        "publisher": "Reactome (OICR, EMBL-EBI, NYU Langone, OHSU)", "landing_url": "https://reactome.org/download-data",
        "access_urls": [reactome_url(k) for k in REACTOME_FILES] + [
            "https://reactome.org/AnalysisService/identifiers/projection (live cross-check only; current release)"],
        "license": "CC BY 4.0 (Reactome data)", "access_conditions": "open download, no registration",
        "retrieved_at": r_ret, "source_version": f"Reactome release {REACTOME_RELEASE} (versioned archive {REACTOME_DL})",
        "update_date": "quarterly releases; release 97 archive files Last-Modified " + ", ".join(
            sorted({m["last_modified"] for m in rm.values() if m["last_modified"]})),
        "data_layer": "condition_molecular", "unit_of_observation": "pathway (gene set, name, parent-child relation)",
        "sample_size": {k: int(v) for k, v in r.items()}, "omics": True, **common,
        "audit": f"data/raw/{REACTOME_SOURCE}/DATA_AUDIT.md", "ingestion_module": "measure_it.omics.reference_data",
        "limitations": [
            "Pathway annotation, not measurement: over-representation shows which curated pathways a gene list touches",
            "Nested, overlapping pathways; hierarchy-based system classification uses a curated anchor table",
            "The AnalysisService cross-check runs the live current release (cached), not the pinned archive",
        ],
    })
    write_registry_entry({
        "source_id": HGNC_SOURCE, "name": "HGNC complete set (HUGO Gene Nomenclature Committee)",
        "publisher": "HGNC (EMBL-EBI)", "landing_url": "https://www.genenames.org/download/",
        "access_urls": [HGNC_URL], "license": "CC0 (HGNC data are freely available without restriction)",
        "access_conditions": "open download, no registration", "retrieved_at": hm["fetched_at"],
        "source_version": f"HGNC complete set, Last-Modified {hm['last_modified'] or 'UNKNOWN'} (sha256 {hm['sha256'][:12]})",
        "update_date": "continuous (file regenerated as nomenclature changes); no release numbers",
        "data_layer": "ontology", "unit_of_observation": "gene (approved symbol and identifiers)",
        "sample_size": {k: (float(v) if isinstance(v, float) else int(v)) for k, v in h.items()}, "omics": False, **common,
        "audit": f"data/raw/{HGNC_SOURCE}/DATA_AUDIT.md", "ingestion_module": "measure_it.omics.reference_data",
        "limitations": [
            "Unversioned file: a cold checkout without data/raw/hgnc gets the then-current set",
            "Symbol harmonisation maps previous/alias symbols only when unambiguous (coherence.SymbolMapper)",
        ],
    })
    (raw_dir(REACTOME_SOURCE) / "DATA_AUDIT.md").write_text(f"""# DATA AUDIT — Reactome pathway knowledgebase (release {REACTOME_RELEASE})

| Field | Value |
|---|---|
| source_id | {REACTOME_SOURCE} |
| Source (dataset/API name, exact files/endpoints) | {', '.join(REACTOME_FILES.values())} from the versioned archive {REACTOME_DL}; AnalysisService projection endpoint for the optional cross-check only |
| Publishing organization | Reactome (OICR, EMBL-EBI, NYU Langone, OHSU) |
| Retrieval date (UTC) | {r_ret} |
| Source version / release | Reactome release {REACTOME_RELEASE}, pinned by URL (not `/download/current/`) |
| Source update date / cadence | quarterly releases; archive files Last-Modified {', '.join(sorted({m['last_modified'] for m in rm.values() if m['last_modified']}))} |
| License / access conditions | CC BY 4.0; open download |
| Unit of observation | pathway (gene set; name/species; parent-child relation) |
| Sample size (actual, as ingested) | GMT: {r['gmt_pathways']:,} pathways ({r['gmt_pathways_hsa']:,} R-HSA), {r['gmt_distinct_genes']:,} distinct gene symbols; pathway list: {r['pathways_all_species']:,} (all species), {r['pathways_homo_sapiens']:,} Homo sapiens; relations: {r['relations_all']:,} ({r['relations_homo_sapiens']:,} R-HSA to R-HSA); {r['top_level_homo_sapiens']:,} human top-level pathways |
| Geography (resolution, vintage) | none |
| Person-level? | no |
| Geographic? | no |
| Omics? | yes (pathway annotation reference) |
| Wearable? | no |
| True participant linkage across modalities? | no (reference vocabulary; no participants) |

## Files / endpoints retrieved
{_files_md(REACTOME_SOURCE)}

The three release-97 archive files are byte-identical (sha256) to the `https://reactome.org/download/current/` files
that the first molecular-coherence build fetched on 2026-09-23 (current was release 97 then), so pinning the URL
changed no Test 3 number.

## Key variables
GMT: pathway name, stable id (R-HSA-...), member gene symbols (-> over-representation analysis). Pathway list:
stable id -> name, species (Homo sapiens kept). Relation file: parent -> child stable ids (-> top-level ancestors ->
physiological-system classification with the curated anchor table in docs/ANALYSIS_PLAN_MOLECULAR.md).

## Missingness
No missing identifiers in the three files (every row has a stable id; GMT rows have >= 1 gene).

## Linkage strategy
Gene symbols join to HGNC approved symbols (source `hgnc`); pathways join to Open Targets `reactome_top_level_term`
by stable id for a hierarchy cross-check. Condition-level only: nothing joins to a person.

## Limitations and caveats
* Pathway annotation, not measurement; enrichment of a public gene list is condition-level molecular context.
* Pathways are nested and overlap; pathway-level agreement p-values are anti-conservative (noted in Test 3).
* The live AnalysisService cross-check always runs the current release; its cached responses are what the tables use.

## Processed outputs
None in data/processed. Used by {CONSUMERS}.

## Reproduce
`uv run python -m measure_it.omics.reference_data` (files are downloaded once and recorded in MANIFEST.json), then
`uv run python -m measure_it.omics.coherence`.
""")
    (raw_dir(HGNC_SOURCE) / "DATA_AUDIT.md").write_text(f"""# DATA AUDIT — HGNC complete set

| Field | Value |
|---|---|
| source_id | {HGNC_SOURCE} |
| Source (dataset/API name, exact files/endpoints) | {HGNC_FILE} ({HGNC_URL}) |
| Publishing organization | HUGO Gene Nomenclature Committee (HGNC), EMBL-EBI |
| Retrieval date (UTC) | {hm['fetched_at']} |
| Source version / release | unversioned file; Last-Modified {hm['last_modified'] or 'UNKNOWN'}; sha256 {hm['sha256']} |
| Source update date / cadence | continuous (regenerated as nomenclature changes) |
| License / access conditions | CC0; open download |
| Unit of observation | gene (HGNC id, approved symbol, previous and alias symbols, locus group, Ensembl gene id) |
| Sample size (actual, as ingested) | {h['rows']:,} rows; {h['approved']:,} approved; {h['approved_protein_coding']:,} approved protein-coding genes (the Test 3 background) |
| Geography (resolution, vintage) | none |
| Person-level? | no |
| Geographic? | no |
| Omics? | no (nomenclature reference) |
| Wearable? | no |
| True participant linkage across modalities? | no (reference vocabulary; no participants) |

## Files / endpoints retrieved
{_files_md(HGNC_SOURCE)}

The file was seeded from the HTTP-cache response that the first molecular-coherence build used (same bytes), so
moving it under data/raw changed no result.

## Key variables
`symbol` (approved symbol), `prev_symbol` and `alias_symbol` (harmonisation of source gene names), `ensembl_gene_id`
(Open Targets / MapMECFS Ensembl ids -> symbols), `locus_group` ("protein-coding gene" defines the background),
`status` (Approved kept).

## Missingness
Approved genes without an Ensembl gene id: {h['approved_missing_ensembl_id']:,} ({h['approved_missing_ensembl_id_pct']}%).

## Linkage strategy
Symbols and Ensembl ids from Open Targets, GWAS Catalog and MapMECFS are mapped to approved symbols; previous/alias
symbols only when they map to one approved symbol. Condition-level only: nothing joins to a person.

## Limitations and caveats
* No release numbers: a checkout without data/raw/hgnc/{HGNC_FILE} downloads the then-current file (non-reproducible input).
* Ambiguous previous/alias symbols are left unmapped (counted in results/tables/test3_gene_mapping_qc.csv).

## Processed outputs
None in data/processed. Used by {CONSUMERS}.

## Reproduce
`uv run python -m measure_it.omics.reference_data`, then `uv run python -m measure_it.omics.coherence`.
""")
    return st


def run() -> dict:
    return write_audits_and_registry(ensure_files())


if __name__ == "__main__":  # pragma: no cover
    print(json.dumps(run(), indent=1))
