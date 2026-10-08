# Build conventions (read before writing any module)

This repository is built by several contributors working in parallel. These
rules keep their work composable and keep the scientific guardrails enforceable
in code.

## 1. Shared infrastructure — use it, do not re-implement or edit it

| Need | Use | Notes |
|---|---|---|
| Paths, seed, UNKNOWN sentinel | `measure_it.config` | `raw_dir(source_id)`, `PROCESSED`, `SEED=20260923`, `UNKNOWN="UNKNOWN / NOT AVAILABLE"` |
| API calls | `measure_it.http.get_json / post_json / get / post` | Disk-cached under `data/_http_cache`; throttled per host; raises `HTMLInsteadOfData` when an endpoint serves an HTML page (CDC does this with HTTP 200). Never call `requests` directly for API traffic. |
| File downloads | `measure_it.download.download_file(url, source_id, filename)` | Writes `data/raw/<source_id>/MANIFEST.json` (url, bytes, sha256, retrieved_at). Re-runs skip existing files. |
| Provenance | `measure_it.provenance.add_provenance(...)` | Every processed table must carry `PROVENANCE_COLUMNS`. `write_table` enforces it. Extend `EVIDENCE_TYPES` only by reporting it back, not by editing the file. |
| Writing tables | `measure_it.store.write_table(df, name, producer=...)` | Parquet in `data/processed/<name>.parquet` + `<name>.meta.json`. |
| Source registry | `measure_it.registry.write_registry_entry({...})` | Writes `data/raw/<source_id>/registry_entry.yaml`. Never edit `SOURCE_REGISTRY.yaml` directly. |
| ZIP -> ZCTA point + county | `measure_it.geography.crosswalk.zip_to_geo(series)` | Tables `geographies`, `zcta_centroids`, `zcta_county_crosswalk` already exist. |
| State FIPS <-> abbreviation | `measure_it.geography.crosswalk.STATE_FIPS`, `STATE_ABBR_TO_FIPS` | |

Shared files you must NOT modify: `pyproject.toml`, `uv.lock`, `SOURCE_REGISTRY.yaml`,
`src/measure_it/{config,http,download,provenance,store,registry}.py`,
`src/measure_it/measurements/adapter.py`, `src/measure_it/geography/crosswalk.py`,
anything in `configs/` that you did not create, and other contributors' modules.
If a shared file needs a change, report it in your final summary.
Do not run `git commit`, `git add`, or `uv add`. If you need a package that is
not installed, report it; for exploration you may use `uv run --with <pkg>`.
To stop a background process, kill it by PID; `pkill -f`/`pgrep -f` can match
and kill the tool's own shell.

## 2. Table naming

* Canonical tables: `data/processed/<table>.parquet`.
* Per-source partitions: `data/processed/<table>__<source>.parquet`
  (e.g. `participant_wearable_features__nhanes`). `store.union_partitions(table)`
  builds the canonical table later.
* A table that summarises one source but is not a partition of a canonical table must not use the
  `<table>__<source>` form (e.g. the ClinicalTrials.gov per-site history is `ctgov_facility_summary`, not
  `research_site_registry__clinicaltrials_gov`). The pipeline's `canonical_unions` step unions every
  `<table>__<source>` family except those whose canonical table a module builds (`pipeline.OWNED_CANONICALS`).
* Person-level ids are namespaced: `participant_id = "<dataset_id>:<native id>"`
  (e.g. `nhanes:62161`). There is no cross-dataset participant key, by design.
* User-supplied datasets (`measure-it byod`, docs/BRING_YOUR_OWN_DATA.md) use the dataset id `byod_<id>` and the
  partitions `<table>__byod_<id>`; their metadata and aggregate records live in the non-partition tables
  `byod_datasets` and `byod_performance_records` (never `measurement_performance_records__...`, which the union step
  would turn into the canonical table).

## 3. Scientific guardrails that code must enforce

1. Four layers — person, condition-molecular, geographic, facility — plus
   ontology/measurement/derived. The `data_layer` provenance column records
   which. Never join a person-level row to a geographic or facility row as if
   they were the same people. Geographic joins are ecological.
2. Never infer an individual's location from a public clinical dataset.
3. Omics are "participant-linked" ONLY when the source publishes shared
   participant identifiers across modalities, and the overlap was counted.
   Otherwise the output is "condition-level molecular enrichment".
4. Do not fabricate data. If a source is unreachable, gated, or lacks a
   variable, record that (status `blocked`/`partial` in the registry entry and
   in DATA_AUDIT.md) and move on. Do not invent HR/HRV for NHANES; NHANES
   2011-2014 has no heart-rate stream.
5. Do not downscale state estimates to counties. A state value attached to a
   county must keep `source_geographic_resolution="state"` and be labelled as
   inherited context, never as county prevalence.
6. Proxies are labelled. Burden rows carry `burden_evidence_level`
   A direct / B closely matching coded condition / C symptom-comorbidity proxy /
   D none.
7. An FDA record is a deployment-readiness signal, not evidence the device
   diagnoses the target illness. A trial using a measurement is not evidence the
   measurement works.
8. No LLM-generated numbers as data. Tools return `UNKNOWN / NOT AVAILABLE`
   when data are missing.
9. Every score keeps its component columns.
10. Language: "candidate", "evidence-supported", "measurable phenotype",
    "deployment opportunity", "measurement desert", "research readiness".
    Avoid "diagnosed by AI", "proves", "confirms mechanism", "patient has",
    "definitive biomarker", "best clinic", "optimal treatment".

## 4. Geography

* Canonical county vintage is 2024 (Connecticut = 9 planning regions,
  09110-09190). Legacy CT counties (09001-09015) exist in `geographies` with
  `ct_legacy=True`. Sources on legacy CT FIPS will not match 2024 CT regions;
  report the unmatched rows rather than forcing a mapping.
* FIPS codes are strings with leading zeros (state 2 chars, county 5 chars).
* Facility coordinates: use the source's own lat/lon when provided
  (`source_geographic_resolution="point"`); otherwise `zip_to_geo`
  (`"zcta_centroid"`, `geocode_method` column).

## 5. Each data source delivers

1. `src/measure_it/<package>/<module>.py` with a `run()` (or `build()`)
   entry point and `if __name__ == "__main__": run()`; idempotent; uses the cache.
2. Raw files under `data/raw/<source_id>/` (+ MANIFEST.json via download_file).
3. `data/raw/<source_id>/DATA_AUDIT.md` following
   `docs/templates/DATA_AUDIT_TEMPLATE.md`, with measured (not guessed)
   sample sizes and missingness.
4. `data/raw/<source_id>/registry_entry.yaml` via `write_registry_entry`.
5. Processed tables via `write_table`.
6. `tests/test_<module>.py` — fast tests that run on the processed outputs
   (mark with `@pytest.mark.data`) and pure-function unit tests.

## 6. Reproduce

One step: `uv run python -m measure_it.<package>.<module>` (each module's own CLI), or through the pipeline
`uv run measure-it pipeline --only <step>`. Everything (docs/REPRODUCIBILITY.md has timings and the determinism
table):

```bash
uv sync --frozen                                     # the environment pinned in uv.lock
uv run measure-it pipeline --list                    # the ordered DAG: 75 steps, entry points, dependencies
uv run measure-it pipeline --offline --workers 16    # all steps from data/raw + data/_http_cache only
uv run measure-it pipeline --offline --from scoring_report        # re-run a step and everything after it
uv run measure-it pipeline --only measurement_report --only build_db
uv run measure-it build-db                           # DuckDB: every table + SPEC views (spec.*) + _sources
uv run measure-it validate                           # provenance, required datasets, person-layer geography,
                                                     # object ids, product language, data audits (exit 1 on failure)
uv run measure-it merge-registry                     # data/raw/*/registry_entry.yaml -> SOURCE_REGISTRY.yaml
uv run pytest -m "not network"                       # tests (data tests need the processed tables)
uv run python -m measure_it.reproduce compare --baseline <copy of data/processed + results/tables> --since <run start>
```

`--offline` sets `MEASURE_IT_OFFLINE=1`: `measure_it.http.request` answers only from `data/_http_cache`,
`measure_it.download.download_file` only returns files already recorded in `MANIFEST.json`, both raise
`http.OfflineCacheMiss` otherwise, and every process installs a socket guard that refuses outbound connections.
Modules with their own offline switch get it automatically (`nih_reporter.run(offline=True)`,
`mapmecfs.run(fetch=False)`). Without `--offline` the same steps fill a cold cache from the network.
Steps run as subprocesses with logs in `results/pipeline_runs/<run id>/<step>.log` and one JSON line per
start/end in `results/pipeline_runs/<run id>.jsonl`. New modules add a `Step` to `measure_it.pipeline.STEPS`
(entry point, dependencies, worker parameter) rather than a separate script.

## 7. Object ids for traceability (`trace_evidence`)

Any row a recommendation cites must be addressable as `<namespace>:<native id>`:

| namespace | native id | table |
|---|---|---|
| `condition` | canonical_condition_id (e.g. `long_covid`) | condition_registry |
| `measurement` | measurement class id or bundle id | measurement_registry |
| `trial` | NCT id | clinical_trials / trial_sites |
| `nih` | appl_id | nih_projects |
| `npi` | NPI | providers |
| `hrsa_site` | BPHC site number | facilities__hrsa |
| `facility` | facility_id from clinic_registry / research_site_registry | clinic_registry |
| `fda` | K/P/DEN number or product code (`fda:code:DQA`) | fda_510k / measurement_regulatory_status |
| `gene` | Ensembl gene id | condition_molecular_evidence |
| `variant` | rsID | condition_molecular_evidence |
| `gwas_study` | GCST accession | condition_molecular_evidence |
| `geo_series` | GSE accession | geo_study_catalog |
| `geo` | FIPS (`geo:06073`, `geo:06`) | geographies / geo_condition_features |
| `opportunity` | deterministic id `<condition>|<measurement>|<geo_id>` | deployment_opportunities |
| `signature` | `<dataset>|<phenotype>|<feature>` | phenotype_signatures |
| `measurement_evidence` | `<condition_id>|<measurement_id>` | condition_measurement_evidence |
| `measurement_signal` | `<dataset>|<label>|<class>|<scope>` | measurement_phenotype_signal |
| `phenotype_evidence` | `<phenotype>|<measurement>` | phase3_measurement_evidence |
| `molbio` | `<condition>|<system>|<measurement_class>[|<catalogue>]` | measurable_biology |
| `reactome` | Reactome stable id (R-HSA-...) | results/tables/test3_reactome_ora.csv |
| `measurement_performance` | `<condition_id>|<measurement_id>` | measurement_performance (scoring.metric_link) |

Tables that recommendations cite should carry an `object_id` column built this way.

## 8. Curated relevance (configs/relevance.yaml)

`condition_specialties`, `measurement_implementers` (with `deployment_complexity`)
and `measurement_bundles` are curated assumptions (evidence_type `curated_config`).
Use them; do not re-invent them per module. Bundles (e.g.
`wearable_autonomic_activity_monitoring`) resolve to member measurement classes.
