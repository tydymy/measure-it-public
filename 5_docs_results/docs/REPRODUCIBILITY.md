# Reproducibility

This page says how to rebuild every processed table, result table, results/*.md write-up and the DuckDB of this
repository from the raw files and the HTTP cache, without touching the network, and what the reproduction run of
2026-09-24 found when it compared the rebuild with the build that was committed before it.

## 1. The command

```bash
uv sync --frozen                                    # the exact environment in uv.lock (section 2)
uv run measure-it pipeline --offline --workers 16   # 85 steps (the recorded runs below had 75), raw files + data/_http_cache only
uv run pytest -m "not network"                      # unit + data tests on the rebuilt outputs
```

`uv run measure-it pipeline --list` prints the DAG: every step, the entry point it calls (the same function the
module's own `python -m` command calls), its dependencies and its worker parameter. The steps, in order:
geography crosswalk; ontology; the 16 ingestion steps (SVI, ACS, RUCC, HRSA, PLACES, HPS long COVID, Lyme, CMS MMD,
NPPES, ClinicalTrials.gov, NIH RePORTER, openFDA, NHANES incl. accelerometry features, Stanford wearables, MapMECFS,
Data.gov); Open Targets, GWAS Catalog, GEO/SRA, the molecular-evidence union, the Reactome/HGNC reference files,
molecular coherence and the measurable-biology graph; NHANES models, signatures, unsupervised structure and report;
Stanford models; digital person; wearable adapter; measurement discovery (as-shipped and audited patterns), the
precision-audit summaries, Phase 3 measurement evidence, the measurement report and the measurement-adapter registry;
geography features; facility registry, its location-shuffle control and its example outputs; the scoring report
(deployment opportunities, Tests 5-7, sensitivity analyses, demo recommendation); the data figures (SPEC Figures 3-9,
results/maps and their CAPTIONS.md block); canonical unions of partitioned tables; registry merge; the schematic
figures (Figures 1, 2 and 10; Figure 2 is drawn from the SOURCE_REGISTRY.yaml that the registry merge writes); DuckDB
build; validation. The two figure steps joined the DAG after the reproduction runs of section 3 (47 steps then). On
2026-09-24 the DAG grew to 75 steps: five NHANES steps (`ingest_nhanes_labs_extended`, `ingest_nhanes_2003_2006`,
`wearable_nhanes0306_features`, `wearable_nhanes_layers`, `wearable_nhanes0306_models`) and 20 steps for the labs, omics
and device analyses (section 5.8): `ingest_muscle_me_charlton`, `ingest_fm_thermography`, `ingest_appelman_lc_pem`,
`omics_geo_cohorts`, `omics_person_linked`, `wearable_muscle_me_steps`, `wearable_fm_thermography`, `labs_plan`,
`labs_vs_diagnosis`, `labs_appelman_metabolomics`, `labs_report`, `labs_klein2023_mylc`, `labs_endo_arg1_repod`,
`labs_heds_hsd_olink_cinquina2026`, `labs_published_evidence`, `labs_dataset_discovery`, `labs_aou_lab_concepts`,
`measurement_published_evidence`, `measurement_device_dataset_discovery`, all ordered before `canonical_unions`.

Useful variants: `--only STEP` (repeatable; inputs must exist), `--from STEP` (that step and every later one),
no `--offline` (the same steps fill a cold cache from the network), `--workers N` (worker budget).

Each step runs in its own subprocess. `results/pipeline_runs/<run id>.jsonl` gets one JSON line per step start and
end (start/end time, duration, worker count, return code, status); `results/pipeline_runs/<run id>/<step>.log` holds
the step's stdout/stderr (git-ignored). A failed step stops its dependants and the run exits non-zero; the final
`validate` step exits non-zero when any check fails (`uv run measure-it validate` runs it alone).

To compare a rebuild with an earlier build, copy `data/processed` and `results/tables` to a directory first, run the
pipeline, then:

```bash
uv run python -m measure_it.reproduce compare --baseline <copy> --since <run start, UTC ISO> --out <csv>
```

## 2. Environment

| item | value |
|---|---|
| host | Linux 6.17 aarch64 (NVIDIA Grace, 20 cores, 121 GB RAM), glibc 2.39 |
| Python | 3.12.13 (`requires-python >= 3.12`) |
| uv | 0.12.2; environment = `uv.lock` (lock revision 3, 140 packages at the reproduction runs; `xlrd` 2.0.2 and `mat-io` 1.0.0 (+ `h5py`) added 2026-09-24 for the `.xls` / MATLAB supplements), `uv sync --frozen` |
| key packages (uv.lock) | numpy 2.5.3, pandas 3.0.6, pyarrow 25.0.1, duckdb 1.5.5, polars 1.44.2, scipy 1.18.1, scikit-learn 1.9.1, statsmodels 0.15.0, joblib 1.6.0, threadpoolctl 3.7.0, geopandas 1.1.4, shapely 2.1.2, pyproj 3.8.0, umap-learn 0.5.12, numba 0.67.0, networkx 3.7, rapidfuzz 3.14.6, pyreadstat 1.3.6, requests 2.34.2, pyyaml 6.0.3, typer 0.27.2, matplotlib 3.11.2 |
| seeds | `config.SEED = 20260923` for every stochastic step; parallel tasks derive their own seeds from it |
| threads | `wearables.models` pins OMP/MKL/OpenBLAS threads to 1 per worker, `wearables.unsupervised` to 4 |
| inputs | `data/raw` (37 GB and 27 sources at the reproduction runs; 40 GB and 37 sources on 2026-09-24, each with MANIFEST.json: url, bytes, sha256, retrieved_at), `data/_http_cache` (503 MB, 4,644 cached API responses; 550 MB on 2026-09-24), and three derived caches in `data/interim` (section 4) |

## 3. Step timings (reproduction run)

Run `20260924T062847Z` (`--offline --workers 16`): wall time 91.4 min, status `ok`; sum of step durations 197.8 min. The first run `20260924T045326Z` took 92.0 min (its only failure was the final `validate` step, section 5.5). Durations in seconds; `workers` is the step's own worker count within the budget of 16. The last column is the 2026-09-24 fixups (section 5.6): run `20260924T114249Z` (`--offline --workers 16`, 12 steps, wall time 290 s: `facilities_shuffle_test` ran beside `scoring_report` and the steps after it) plus the two single-step runs named in the cell; `figures_data` also took 18 s in its own runs `20260924T103316Z` and `20260924T110035Z`. The two figure steps were not in the DAG of the two reproduction runs. The review re-run `20260924T124442Z` (`--offline --workers 16`, 13 steps, wall time 489 s, ok) re-ran every step downstream of `ingest_hrsa` and `ingest_cms_mmd` (section 5.6).

| # | step | workers | 20260924T045326Z (s) | 20260924T062847Z (s) | 2026-09-24 fixups (s) | review re-run `20260924T124442Z` (s) |
|---|---|---|---|---|---|---|
| 1 | `geography_crosswalk` | 1 | 1.0 | 1.0 | — | — |
| 2 | `ontology` | 1 | 7.0 | 5.0 | — | — |
| 3 | `ingest_svi` | 1 | 1.0 | 1.0 | — | — |
| 4 | `ingest_census_acs` | 1 | 10 | 10 | — | — |
| 5 | `ingest_rucc` | 1 | 1.0 | 1.0 | — | — |
| 6 | `ingest_hrsa` | 1 | 7.0 | 7.0 | — | 7.0 |
| 7 | `ingest_cdc_places` | 1 | 9.0 | 9.0 | — | — |
| 8 | `ingest_cdc_long_covid` | 1 | 1.0 | 1.0 | — | — |
| 9 | `ingest_cdc_lyme` | 1 | 2.0 | 1.0 | — | — |
| 10 | `ingest_cms_mmd` | 1 | 19 | 14 | 14 (`20260924T113143Z`) | 14 |
| 11 | `ingest_nppes` | 4 | 43 | 39 | — | — |
| 12 | `ingest_clinicaltrials` | 1 | 57 | 55 | — | — |
| 13 | `ingest_nih_reporter` | 1 | 30 | 29 | — | — |
| 14 | `ingest_openfda` | 1 | 2.0 | 2.0 | — | — |
| 15 | `ingest_nhanes` | 8 | 14 | 14 | — | — |
| 16 | `ingest_stanford_wearables` | 8 | 479 | 484 | — | — |
| 17 | `ingest_mapmecfs` | 1 | 23 | 22 | — | — |
| 18 | `datagov` | 1 | 1.0 | 1.0 | — | — |
| 19 | `omics_open_targets` | 1 | 5.0 | 5.0 | — | — |
| 20 | `omics_gwas_catalog` | 1 | 85 | 82 | — | — |
| 21 | `omics_geo_sra` | 1 | 2.0 | 2.0 | — | — |
| 22 | `omics_union` | 1 | 1.0 | 1.0 | — | — |
| 23 | `omics_reference_data` | 1 | 1.0 | 1.0 | — | — |
| 24 | `omics_coherence` | 1 | 115 | 115 | — | — |
| 25 | `omics_graph` | 1 | 3.0 | 3.0 | — | — |
| 26 | `wearable_models` | 6 | 5,097 | 5,018 | — | — |
| 27 | `wearable_signatures` | 1 | 2.0 | 2.0 | — | — |
| 28 | `wearable_unsupervised` | 1 | 74 | 65 | — | — |
| 29 | `wearable_nhanes_report` | 1 | 1.0 | 1.0 | — | — |
| 30 | `stanford_models` | 6 | 4,857 | 4,817 | — | — |
| 31 | `digital_person` | 1 | 21 | 17 | — | — |
| 32 | `wearable_adapter` | 4 | 45 | 46 | — | — |
| 33 | `measurement_discovery_v1` | 1 | 61 | 62 | — | — |
| 34 | `measurement_discovery` | 1 | 118 | 117 | — | — |
| 35 | `measurement_precision_audit` | 1 | 1.0 | 1.0 | — | — |
| 36 | `measurement_phenotype_evidence` | 1 | 8.0 | 7.0 | 7 | — |
| 37 | `measurement_report` | 1 | 1.0 | 1.0 | 1.0 | — |
| 38 | `measurement_adapters` | 1 | 1.0 | 1.0 | 1.0 | — |
| 39 | `geography_features` | 1 | 110 | 113 | 74 (`20260924T113255Z`) | 77 |
| 40 | `facilities_registry` | 1 | 216 | 206 | — | 190 |
| 41 | `facilities_shuffle_test` | 6 | 308 | 310 | 290 | 292 |
| 42 | `facilities_examples` | 1 | 7.0 | 7.0 | 6.0 | 6.0 |
| 43 | `scoring_report` | 6 | 122 | 119 | 146 | 147 |
| 44 | `figures_data` | 1 | — | — | 21 | 21 |
| 45 | `canonical_unions` | 1 | 2.0 | 2.0 | 2.0 | 2.0 |
| 46 | `registry_merge` | 1 | 1.0 | 1.0 | 1.0 | 1.0 |
| 47 | `figures_schematics` | 1 | — | — | 4.0 | 4.0 |
| 48 | `build_db` | 1 | 49 | 48 | 53 | 55 |
| 49 | `validate` | 1 | 3.0 | 4.0 | 4.0 | 4.0 |

The two long poles are the NHANES models (5 x 5 CV, 1,200 label permutations, 120 wearable-row shuffles, 1,000 bootstrap resamples; 6 workers) and the Stanford models (6 workers); they run in parallel, and everything else fits around them.

## 4. Offline guarantee

`--offline` sets `MEASURE_IT_OFFLINE=1` for the pipeline and every step subprocess. With it:

* `measure_it.http.request` answers only from `data/_http_cache` and raises `http.OfflineCacheMiss` on a miss; a
  `refresh=True` call is served from the cached response when one exists (a refresh cannot happen offline) and raises
  otherwise.
* `measure_it.download.download_file` returns a raw file only when it exists and is recorded in the source's
  `MANIFEST.json`; anything else raises `OfflineCacheMiss` instead of downloading.
* The two places that stream bulk files with `requests` directly (the NHANES segmented PAXMIN download and the NPPES
  V.1 byte-range reader) call `http.ensure_online()` first.
* `measure_it/__init__.py` installs a socket guard in every process that imports the package, joblib/loky and
  ProcessPool workers included: `socket.getaddrinfo` and `connect` to any non-loopback address raise
  `http.OfflineNetworkBlocked`. A step that completes under the guard did no network I/O; both runs completed every
  step under it, so the reproduction used `data/raw`, `data/_http_cache` and `data/interim` only.
* Modules with their own offline switch get it automatically: `nih_reporter.run(offline=True)` (rebuild from the saved
  RePORTER pages) and `mapmecfs.run(fetch=False)` (rebuild from data/raw and the saved access log).
* Release / vintage probes that the online build answered from the network are recorded as not reachable in an
  offline run's audits (`data/raw/census_acs/DATA_AUDIT.md`: ACS 2025 `error: OfflineCacheMiss` where the online build
  recorded HTTP 404; `data/raw/cdc_svi/registry_entry.yaml`: SVI 2024 status `None` instead of 404). The vintage used
  is the same; only the probe record differs.

Derived caches reused from `data/interim` (not network data; each is rebuilt from `data/raw` when absent):

| cache | built by | holds | invalidation |
|---|---|---|---|
| `nhanes_2011_2014/paxmin_{G,H}.parquet` | `wearables.nhanes_features.convert_paxmin` | lossless conversion of the 8-9 GB PAXMIN XPT files (MIMS as float32, restored to the published 3 decimals) | delete the file |
| `stanford_longcovid_wearables/2026-09-23.3/` | `ingestion.stanford_wearables` | label-free per-participant daily summaries (reference-dependent features are recomputed every run) | `PIPELINE_VERSION` in the module |
| `measurement_discovery/measurement_mentions_v1.parquet` | step `measurement_discovery_v1` | as-shipped pattern mentions for the precision audit | rewritten every run |

The first two were rebuilt from the raw files in a separate run and reproduced the outputs (section 5.3).

`data/interim/nhanes_wearable_models/*.npz` are intermediate model outputs rewritten by every run of
`wearables.models`. `data/interim/repro_backup/` holds the comparison baselines of this page (git-ignored).

## 5. Determinism

### 5.1 Run-to-run (the determinism check)

The pipeline was run twice with identical code and inputs (`20260924T045326Z` and `20260924T062847Z`, both `--offline --workers 16`, fresh processes, randomised PYTHONHASHSEED per process). Every processed table and result table of the second run was compared with the first (`measure_it.reproduce compare`: row count, column set and dtypes, order-independent hash of all rows; per-column hashes when they differ):

| table class | identical | identical (run-stamped columns normalised) | identical except runtime fields | all |
|---|---|---|---|---|
| processed canonical / derived table | 66 | 15 | 0 | 81 |
| processed partition (`<table>__<source>`) | 36 | 5 | 0 | 41 |
| result table (results/tables *.csv, *.json) | 129 | 26 | 3 | 158 |
| all | 231 | 46 | 3 | 280 |

Tables that are not identical between the two runs:

| table | status | columns | note |
|---|---|---|---|
| `nhanes_run_metadata.csv` | identical except runtime fields | value | runtime fields only (seconds_* per stage, last_run_at) |
| `nhanes_unsupervised_summary.csv` | identical except runtime fields | value | runtime field only (seconds) |
| `test6_clinic_shuffle_checks.csv` | identical except runtime fields | value | runtime field only (runtime_s) |

Run-stamped columns (48 tables) are compared after replacing every ISO timestamp by a placeholder; apart from the timestamp their content is identical. They are `retrieved_at` of derived tables (= build time), `source_version` strings that embed a build time (`scoring build <time>`, `audit built <time>`), the `generated_at` inside `deployment_candidates.recommendation_json`, and the run-time fields of three result tables (`nhanes_run_metadata.csv` seconds per stage and `last_run_at`, `nhanes_unsupervised_summary.csv` seconds, `test6_clinic_shuffle_checks.csv` runtime_s). The `.meta.json` sidecars (`created_at`), the step logs and the results/*.md generation lines also carry run times.

### 5.2 Rebuild vs the committed build (before the integration fixes)

`data/processed` and `results/tables` were copied before any integration change (commit ec3e5b2 state) and compared with the first full rebuild `20260924T045326Z`. Every difference was traced to a cause; none is unexplained:

| table class | changed as expected | identical | identical (run-stamped columns normalised) | identical except runtime fields | all |
|---|---|---|---|---|---|
| processed canonical / derived table | 23 | 52 | 6 | 0 | 81 |
| processed partition (`<table>__<source>`) | 6 | 33 | 3 | 0 | 42 |
| result table (results/tables *.csv, *.json) | 40 | 93 | 22 | 3 | 158 |
| all | 69 | 178 | 31 | 3 | 281 |

Cause codes: **A1** audited measurement patterns in configs/measurements.yaml used by every consumer (the facility registry now tags trials with the audited patterns, so technology experience, research readiness and saturation change for the two bundles that contain a revised class, `wearable_autonomic_activity_monitoring` and `autonomic_function_testing`; nothing changes for `nailfold_capillaroscopy` or `exercise_capacity_testing`); **A2** bundle alias + several adapters per bundle; **A3** inherited-burden uncertainty multiplier (x2) in the scoring Monte Carlo (only `burden_mc_sd` and Monte Carlo rank intervals of combinations with long COVID as a member change; observed composites and ranks do not); **A4** ClinicalTrials.gov `city_centroid` resolution, placeholder facility names, `ctgov_facility_summary` rename (the facility registry keeps the identical 17,685 U.S. site rows); **A6** Reactome release-97 archive (byte-identical to the files used before) and HGNC from data/raw/hgnc: version strings only; **A7** `data_layer = metadata`; **A8** NHANES SLQ050/SLQ060 sleep-label signatures (new rows only; every earlier row identical); **B** new canonical tables from the union step; **S** stale committed table: the PLACES, ACS, NPPES, NIH RePORTER and ClinicalTrials.gov partitions had been built at 20:53-21:06 UTC on 2026-09-23, before the geography crosswalk was rebuilt at 22:54 UTC with the fix committed in 2478d0a (the 30 ZCTAs whose Census internal point lies outside every 2024 county polygon get their nearest county), and were never rebuilt; the rebuild fills state/county for those ZCTAs (757 PLACES ZCTA rows, 30 ACS ZCTAs), relabels the county method of 3,272 providers and 96 NIH awards (same county), adds a ZIP county for 60 trial sites (mostly ZIP 60611) and so moves trial counts in 44 county rows of `geo_condition_features` and a few Test 4 correlations in the 6th decimal. The scoring components were not affected by S (burden, vulnerability, desert and clinic-capacity percentiles are identical in all 48 combinations).

| table | rows before -> after | cause |
|---|---|---|
| `condition_burden_coverage.parquet` | new (14) | B |
| `ctgov_facility_summary.parquet` | new (28,874; = the old partition's rows) | A4 rename of `research_site_registry__clinicaltrials_gov` (content differs only in facility_generic / facility_identifiable / source_geographic_resolution) |
| `data_gov_discovery_log.parquet` | 51 -> 54 | A7 (all rows), A6 (3 new query rows for reactome / hgnc) |
| `data_gov_search_results.parquet` | 396 -> 399 | A7 (all rows), A6 (3 new rows) |
| `deployment_candidates.parquet` | 10 -> 10 | A1 (5 of the 10 demo regions differ), A3 |
| `deployment_opportunities.parquet` | 76,680 -> 76,680 | A1 (research readiness / saturation of the wearable and autonomic-testing bundles), A3 (Monte Carlo of long-COVID combinations; 2 new columns) |
| `facility_trials.parquet` | 17,598 -> 17,598 | A1 (`measurement_classes` of 3,513 rows) |
| `geo_condition_features.parquet` | 46,074 -> 46,074 | S (trial counts of 44 county rows) |
| `geo_context_zcta.parquet` | new (925,539) | B |
| `geo_context_zcta__acs.parquet` | 33,772 -> 33,772 | S (state/county of 30 ZCTAs) |
| `geo_context_zcta__cdc_places.parquet` | 891,767 -> 891,767 | S (state/county of 757 ZCTA rows) |
| `measurement_adapter_registry.parquet` | 10 -> 10 | A2 (`bundles_naming_this_adapter`) |
| `measurement_phenotype_signal.parquet` | 60 -> 68 | A8 |
| `measurement_registry.parquet` | 30 -> 30 | A2 (alias; new column `adapters_declared`), A8 (non-condition labels) |
| `modality_inventory.parquet` | new (109) | B |
| `modality_inventory__mapmecfs.parquet` | 109 -> 109 | A7 (was derived) |
| `modality_overlap.parquet` | new (5,886) | B |
| `modality_overlap__mapmecfs.parquet` | 5,886 -> 5,886 | A7 (was derived) |
| `nih_projects.parquet` | 9,342 -> 9,342 | S (county method label of 96 awards; county unchanged) |
| `participant_cluster_assignments.parquet` | new (17,732) | B |
| `participant_mortality.parquet` | new (19,931) | B |
| `phase3_measurement_evidence.parquet` | 32 -> 32 | A8 |
| `phenotype_signatures.parquet` | new (441) | B |
| `phenotype_signatures__nhanes.parquet` | 364 -> 416 | A8 |
| `providers.parquet` | 1,456,346 -> 1,456,346 | S (county method label of 3,272 rows; county unchanged) |
| `research_site_registry.parquet` | 8,722 -> 8,722 | A1 (measurement classes studied) |
| `research_site_registry__clinicaltrials_gov.parquet` | 28,874 -> removed | A4 (renamed; old file removed) |
| `trial_sites.parquet` | 47,482 -> 47,482 | A4 (resolution of 46,490 geoPoint rows; 857 placeholder names), S (ZIP county of 60 sites) |
| `wearable_reference.parquet` | new (399) | B |
| `demo_top10_candidate_sites.csv` | 163 -> 171 | A1, A3 |
| `demo_top10_regions.csv` | 10 -> 10 | A1, A3 |
| `demo_top10_research_evidence.csv` | 105 -> 79 | A1, A3 |
| `facility_candidate_examples.csv` | 60 -> 60 | A1 (matcher technology experience) |
| `facility_resolution_summary.csv` | 75 -> 75 | A4 (532 excluded U.S. site rows move from the non-facility / place-name categories to `us_rows_sponsor_placeholder`; the 17,685 kept rows are identical) |
| `geo_desert_sensitivity.csv` | 49 -> 49 | S |
| `geo_trial_radius_qa.csv` | 14 -> 14 | S |
| `measurement_registry.csv` | 30 -> 30 | A2, A8 |
| `nhanes_phenotype_signatures.csv` | 364 -> 416 | A8 |
| `nhanes_run_metadata.csv` | 14 -> 14 | runtime fields only |
| `nhanes_signature_summary.csv` | 14 -> 16 | A8 |
| `nhanes_targets.csv` | 14 -> 16 | A8 |
| `nhanes_unsupervised_summary.csv` | 10 -> 10 | runtime field only |
| `phase3_measurement_evidence.csv` | 32 -> 32 | A8 |
| `phase3_phenotype_summary.csv` | 4 -> 4 | A8 |
| `phase3_signal_results.csv` | 60 -> 68 | A8 |
| `phase3_verdict_sensitivity.csv` | 60 -> 68 | A8 |
| `test4_partial_correlations.csv` | 22 -> 22 | S |
| `test4_population_correlations.csv` | 22 -> 22 | S |
| `test4_verdicts.csv` | 22 -> 22 | S (6th decimal of one correlation) |
| `test5_component_correlations.csv` | 720 -> 720 | A1, A3 |
| `test5_contrasts.csv` | 288 -> 288 | A1, A3 |
| `test5_movements.csv` | 5,430 -> 5,348 | A1, A3 |
| `test5_readiness_zero_counts.csv` | 144 -> 144 | A1, A3 |
| `test6_clinic_shuffle.csv` | 258 -> 258 | A1 |
| `test6_clinic_shuffle_checks.csv` | 9 -> 9 | runtime field only |
| `test6_clinic_shuffle_draws.csv` | 36,900 -> 36,900 | A1 |
| `test6_clinic_shuffle_geographies.csv` | 603 -> 603 | A1 |
| `test6_scoring_clinic_shuffle.csv` | 16 -> 16 | A1, A3 |
| `test6_scoring_clinic_shuffle_draws.csv` | 1,200 -> 1,200 | A1, A3 |
| `test6_scoring_label_shuffle.csv` | 26 -> 26 | A1, A3 |
| `test6_scoring_label_shuffle_draws.csv` | 4,400 -> 4,400 | A1, A3 |
| `test7_monte_carlo_decomposition.csv` | 48 -> 48 | A1, A3 |
| `test7_pairwise_tau.csv` | 1,344 -> 1,344 | A1, A3 |
| `test7_region_stability_primary.csv` | 2,406 -> 2,406 | A1, A3 |
| `test7_region_stability_top25_any.csv` | 2,951 -> 2,908 | A1, A3 |
| `test7_sensitivity_comparisons.csv` | 10 -> 10 | A1, A3 |
| `test7_sensitivity_s1_tau.csv` | 25 -> 25 | A1, A3 |
| `test7_sensitivity_s1b_inherited_multiplier.csv` | new (25) | A3 (new sensitivity S1b) |
| `test7_summary.csv` | 48 -> 48 | A1, A3 |
| `test7_unstable_regions.csv` | 2,199 -> 2,174 | A1, A3 |
| `facility_build_stats.json` | (JSON) | A4 (ctgov_filter counts move between categories; 17,685 rows kept either way) |
| `test3_run_metadata.json` | (JSON) | A6 (version strings; built_at) |

### 5.2.1 Published numbers that changed

All results/*.md files are regenerated by their modules. The numbers that moved (every one follows from A1, A3 or
A8; results/SCORING_RESULTS.md has the full, regenerated text). **These are the values of the first rebuild
(`20260924T045326Z`); the scoring numbers below were superseded later the same day (plan deviations 7 and 8 and the
fixups of section 5.6); the current values are in section 5.6 and results/SCORING_RESULTS.md.**

* **Demo recommendation** (Long COVID or ME/CFS x wearable autonomic/activity monitoring, county, equal weights):
  5 of the 10 regions changed (A1: research readiness of the wearable bundle now counts only trials whose registered
  text matches the audited patterns). Now: Oklahoma County OK, Wayne County MI, Franklin County AL, George County MS,
  Crittenden County AR, Bossier Parish LA, Webster Parish LA, Scotland County NC, Greater Bridgeport Planning Region CT,
  Forrest County MS (before: Scotland, Robeson, Forrest, Oklahoma, Yazoo, Washington Parish, Hoke, Wayne, George,
  Perry). Greater Bridgeport is one of the 12 counties whose set burden rests on the inherited long-COVID value alone
  (no CMS value for Connecticut planning regions), which the results already flag as an artefact.
* **Monte Carlo** (A1 + A3): median 5th-95th percentile rank-interval width of the primary top-25 228 -> 327 ranks;
  top-10 counties that stay in the top-10 in at least half of the draws 4 -> 0. The new sensitivity S1b separates the
  two: with the old multiplier (1.0) the width is 255 and 1 top-10 county is stable, so A1 accounts for 228 -> 255 and
  the inherited-burden multiplier for 255 -> 327.
* **Test 5**: burden-only vs full top-25 overlap 3 -> 4 (Jaccard 0.064 -> 0.087), Kendall tau-b 0.459 -> 0.467;
  need + desert vs full keeps 3 -> 8 of the top-25 (Jaccard 0.064 -> 0.190).
* **Test 7**: Kendall's W over the 8 named weight sets 0.793 -> 0.805 (without burden_only 0.833 -> 0.846), over 1,000
  random vectors 0.704 -> 0.714; default top-25 regions meeting the instability rule 18 -> 23.
* **Test 6**: (b) research-readiness Spearman of shuffled vs observed 0.582 -> 0.471; (a) composite vs real burden
  0.650 -> 0.658 observed (null mean 0.102 -> 0.100).
* **Facility Test 6, clinic-location shuffle** (A1; hand-written in `docs/FACILITY_MATCHING.md` §4-§5, updated at
  review): technology experience of matched Long COVID sites 0.089 -> 0.022, so eligibility, the selected sites and
  every Test 6 metric moved. The supported nonmetro technology-experience deficit (0.024 vs 0.051) and overall deficit
  (0.089 vs 0.106) are no longer supported (0.003 vs 0.005, p 0.26; 0.022 vs 0.018, p 0.13); the overall FQHC
  difference became supported and the overall HPSA difference sample-specific; ME/CFS x wearable supported
  metro-nonmetro contrasts 7 -> 8 of 10. Example pools: San Diego facilities with technology experience 41 -> 3,
  Cook County 20 -> 0.
* **Phase 3** (A8): 68 signal results instead of 60 (the two SLQ labels x 4 scopes); no verdict of an earlier row
  changed; the sleep/circadian gap note now says the signatures exist.
* Unchanged: every NHANES and Stanford model, permutation and signature (earlier rows) table, every
  molecular-coherence table, the geography tables except the S rows, and `facilities`, `clinic_registry`,
  `facility_source_links`, `facility_nih_projects` (A1 changes only `facility_trials` and `research_site_registry`).

### 5.3 Cold derived caches

Run `20260924T080151Z` moved the NHANES PAXMIN parquet cache and the Stanford per-participant cache out of
`data/interim` and re-ran `ingest_nhanes` and `ingest_stanford_wearables` offline (56 s and 588 s), so both caches were
rebuilt from the raw XPT and ZIP files. The rebuilt PAXMIN parquet files equal the old ones row for row
(78,126,856 and 88,223,479 minute rows, order-independent hash), and every NHANES output is identical. The Stanford
values are identical too, but `participant_wearable_daily__stanford_covid` had the column `step_stream_present` in a
different position: the old cache predates the stepless-day rule, so `derive_features` appended the column at the end,
while a fresh cache has it before `rhr_computable`: the output schema depended on the cache state. Fixed in `ingestion.stanford_wearables.derive_features` (the column is placed before `rhr_computable` in
both cases); re-running the step from the OLD cache with the fix (`20260924T081401Z`) reproduced the cold output
exactly, and `digital_person`, `wearable_adapter`, the unions, the registry merge, the DuckDB and validation were re-run
(`20260924T082119Z`). `stanford_models` reads the table by column name, so its outputs are unaffected.

### 5.4 Final state vs the reproduction run

After the two fixes below (Stanford column order; desert top-25 tie order, `geography_features` re-run as
`20260924T083324Z`), the committed tables differ from run 2 in exactly three files
(`results/pipeline_runs/final_vs_20260924T062847Z.csv`): `participant_wearable_daily__stanford_covid` and the canonical
`participant_wearable_daily` (column order only: same columns, values and row order) and
`results/tables/geo_desert_top25.csv` (row order among exact ties, and one tie at the top-25 boundary: PTLDS state
list now ends with Michigan instead of Mississippi, both at desert 0.5123; nothing reads this file). Every other table
is identical to run 2.

### 5.4.1 Independent review re-run

The review fixed three things and re-ran `--offline --workers 6 --from scoring_report` (run `20260924T085911Z`,
175 s; scoring_report, canonical_unions, registry_merge, build_db, validate):
* the S1b sentence of results/SCORING_RESULTS.md section 8 named the multipliers in the opposite order to its numbers
  (now "multiplier 2 (configured) vs 1: width 327 vs 255");
* `canonical_unions` cast any column whose dtype differed across partitions to string, so `phenotype_signatures.n_cases`
  / `n_controls` (int64 vs Int64) and `geo_context_zcta.total_population` were strings in the canonical tables and the
  DuckDB; integer kinds now become Int64 and numeric kinds Float64 (`pipeline._common_dtype`);
* `validate`'s product-language check excused every forbidden term on a line that contained a rule word anywhere
  ("never", "avoid"), so a sentence that used a forbidden term and also contained "never" passed; a term is now excused only when quoted or when it is a bare
  item of a list introduced by the rule word ("Avoid (unless directly supported): diagnosed by AI; proves; ...").

Compared with run 2 (`reproduce compare`, window from 2026-09-24T08:00:30Z;
`results/pipeline_runs/20260924T085911Z_review_vs_20260924T062847Z.csv`): 275 identical, the 3 tables of 5.4, and
`phenotype_signatures` / `geo_context_zcta` numerically equal (dtype only). The scoring re-run reproduced every scoring
table and the demo recommendation exactly.

### 5.5 Nondeterminism found and fixed

* The cache-state-dependent column order of `participant_wearable_daily__stanford_covid` (5.3; fixed in
  `ingestion.stanford_wearables.derive_features`).
* `geography.features` built `geo_desert_top25.csv` with an unstable single-key sort (`sort_values("diagnostic_desert")`,
  quicksort), so the order of exact ties, and which tied region makes a top-25 list, depended on the sort algorithm and
  the input row order. The committed file differed from the rebuild in tie order although every value was the same.
  Now a stable sort with explicit tie-breakers (condition, level, geo id).
* No run-to-run nondeterminism: every table of the second run equals the first (section 5.1), including the Monte Carlo,
  permutation, bootstrap and shuffle outputs, the DuckDB-based NPPES aggregation and the parallel steps (NHANES models
  and Stanford models with 6 loky workers, the Stanford ingestion with 8 processes instead of the module default 12,
  the facility and scoring shuffles with 6 processes). Python's per-process string-hash randomisation (PYTHONHASHSEED)
  was left on; no output depends on set iteration order.
* The first full run exposed two false failures in the new `validate` step, fixed before the second run:
  `facility_source_links` rows from ClinicalTrials.gov and NIH carry an empty `object_id` by design (no CONVENTIONS
  namespace) and are now counted as "not addressable" rather than malformed; the MapMECFS registry `audit` field lists
  two files separated by `;` and each is now checked.
* The comparison itself first treated future-dated study dates in the Stanford tables (dates run to 2030-09-03 in the
  release) as run timestamps; `reproduce.compare` now counts a timestamp as run-stamped only inside the run window,
  and those tables compare identical without any normalisation.
* Build-order staleness is the third reproducibility defect of the committed state (cause S in 5.2): tables built
  before an upstream shared table was fixed. The pipeline's explicit DAG removes this class of problem, because a
  rebuild always runs every consumer after its inputs.

### 5.6 Changes after the reproduction run (2026-09-24)

Three changes after the reproduction run moved scoring numbers: the set-burden completeness rule (plan deviation 7:
regions where a member with a defined burden has no value are not ranked), the S6 `specialist_only` sensitivity (plan
deviation 8, new tables only) and the fixups below. The two figure steps joined the DAG (49 steps). New result tables
since the reproduction run: `scoring_incomplete_burden.csv`, `scoring_set_burden_rule_before_after.csv`,
`scoring_set_burden_rule_top10.csv` (deviation 7) and `test7_sensitivity_specialist_only.csv`, `_qa.csv`, `_top25.csv`
(S6); the figure steps write results/figures/fig01-fig10 (PNG + SVG), results/figures/CAPTIONS.md and results/maps.
No processed table was added; `geo_condition_burden__cms_mmd` and `geo_context__cms_mmd` gained the columns
`fips_as_published` and `fips_bridge`, and `mua_designations` gained `county_fips_bridge` (review follow-up below).

**Fixups (2026-09-24, runs `20260924T113143Z` ingest_cms_mmd, `20260924T113255Z` geography_features,
`20260924T114249Z` the 12 steps after them, all `--offline`; after the documentation edits `build_db` and `validate`
ran again as `20260924T122439Z`, 48 s + 3 s, ok).** (1) The CMS MMD ingestion bridges the two 1:1 FIPS
renames CMS still publishes under their pre-2015 codes (46113 -> 46102 Oglala Lakota SD, 02270 -> 02158 Kusilvak AK;
docs/BURDEN_DEFINITIONS.md); (2) `geo_condition_features.burden_unknown_reason` states the measured source reason
instead of "no published value ... (suppressed, cell omitted by the source, or entity absent from the source)"; (3)
one shared measurement resolver (`measure_it.measurements.resolve`) with curated aliases in configs/measurements.yaml
and configs/relevance.yaml (the `aliases` column of `measurement_registry` now lists class aliases too); (4)
provenance, next-step wording and candidate-site coordinates of `rank_deployment_opportunities`, the incomplete-burden
fields of `trace_evidence`, table-level caching of two query functions (no table effect). Compared with the state
before them (`reproduce compare`, window from 11:30 UTC): 240 tables identical, 46 changed (7 processed, 39 result).
Processed: the two CMS partitions (bridge columns and 2 x 12 years of rows moved to the 2024 codes; row counts equal),
`geo_condition_burden` and `geo_condition_features` (the two counties now carry CMS values; reason texts),
`deployment_opportunities` and `deployment_candidates` (one more ranked county, percentiles re-normalised, next-step
text), `measurement_registry` (`aliases`). Result tables: the geography Test 4 / Test 6 tables (n + 2 in the CMS
series; the permutation nulls of the other series moved within Monte Carlo error because Test 6 draws every series
from one random stream, which now consumes two more values), the scoring Test 5-7 tables and the demo tables;
`test6_clinic_shuffle_checks.csv` differs in its runtime field only, and every other facility table is identical (the
resolver change does not alter any facility or shuffle output).

**Review follow-up (2026-09-24, run `20260924T124442Z`, the 13 steps downstream of `ingest_hrsa` and `ingest_cms_mmd`, `--offline --workers 16`, 489 s, ok).** The HRSA ingestion now also places components published on the two 1:1 renames on the 2024 code (`hrsa.FIPS_RENAMES_1TO1`; before, it treated them as legacy codes that only "may" cover the successor county). One in-effect designation is affected: MUA/P 00104 (Wade Hampton Census Area, whole county, published on 02270), so Kusilvak AK's MUA/P coverage is `whole_county` instead of `undetermined_legacy_geography` (counties with an undetermined MUA/P coverage 12 -> 11; counties with >= 1 Designated MUA/P 2,830 -> 2,831; no HPSA component uses either renamed code). MUA/P coverage is context only and enters no burden, desert or opportunity score. Compared with the state after the fixups (`reproduce compare`, window from 12:44 UTC): 279 tables identical, 7 changed, all through that one county: `mua_designations` (new column `county_fips_bridge`; `county_fips_list`), `geo_context__hpsa`, `geo_context` and `geo_condition_features` (MUA/P columns), `facilities` and `clinic_registry` (`county_mua_coverage`), and `test6_clinic_shuffle_checks.csv` (runtime field only). Every scoring table, result table, results/*.md file and CAPTIONS.md reproduced unchanged.

Scoring numbers, first rebuild (section 5.2.1) -> current (primary query, county, equal weights):

| quantity | first rebuild (`20260924T045326Z`) | current |
|---|---|---|
| counties ranked for the primary query | 2,406 | 2,395 (2,394 after deviation 7; +1 with the FIPS bridge) |
| incomplete-burden counties (population-eligible) | — (rule not in force) | 24 (11); 26 (12) before the FIPS bridge |
| demo top-10 | Oklahoma OK, Wayne MI, Franklin AL, George MS, Crittenden AR, Bossier LA, Webster LA, Scotland NC, Greater Bridgeport CT, Forrest MS | Oklahoma OK, Wayne MI, Franklin AL, George MS, Crittenden AR, Bossier LA, Webster LA, Scotland NC, Forrest MS, Lauderdale MS (Greater Bridgeport CT not ranked: incomplete burden) |
| median 5th-95th rank-interval width, top-25 | 327 | 319 |
| top-10 counties in the top-10 in >= half of the draws | 0 | 0 |
| S1b width, multiplier 2 vs 1 | 327 vs 255 | 319 vs 249 (248 before the FIPS bridge) |
| Test 5 burden-only vs full: top-25 overlap / Jaccard / Kendall tau-b | 4 / 0.087 / 0.467 | 4 / 0.087 / 0.467 (0.466 before the FIPS bridge) |
| Test 5 need + desert vs full: top-25 kept / Jaccard | 8 / 0.190 | 8 / 0.190 |
| Test 7 Kendall's W named sets (without burden_only) / random | 0.805 (0.846) / 0.714 | 0.806 (0.847) / 0.715 |
| default top-25 regions meeting the instability rule | 23 | 22 |
| Test 6(b) research-readiness Spearman, shuffled vs observed | 0.471 | 0.465 |
| Test 6(a) composite vs real burden, observed (null mean) | 0.658 (0.100) | 0.657 (0.098; 95% band 0.072-0.122, was 0.076-0.124 before the FIPS bridge) |
| S6 specialist_only: top-10 / top-25 kept, tau-b | — (not declared) | 5 / 17, 0.881 over 2,395 counties |

The FIPS bridge itself left every equal-weight top-10 membership unchanged (two adjacent swaps inside a top-10, in the
demo cluster x autonomic testing and me_cfs x wearable) and changed top-25 membership in 2 of 48 combinations (demo
cluster x capillaroscopy and x wearable); Oglala Lakota County SD is now ranked 443 (burden-only 307) for the primary
query and Kusilvak AK is complete but below the population minimum (results/SCORING_RESULTS.md section 5).

### 5.7 Final full run on the final code (2026-09-24, after the report review)

Run `20260924T150209Z` (`--offline --workers 16`, all 49 steps, 104.4 min wall time on a host shared with other jobs)
re-ran the whole DAG on the final code, after the review fixes (scoring progress on stderr, the matcher's absence note
naming its condition, the phenotype-signal verdict in `technology_evidence`, Figure 10's S5 sentence, run logs with
relative paths). `data/processed` and `results/tables` were copied to `data/interim/repro_backup/pre_final_20260924`
first; `reproduce compare` (window from 15:02:09 UTC; `results/pipeline_runs/20260924T150209Z_vs_pre_final.csv`):

| outcome | tables |
|---|---|
| identical (run timestamps normalised) | 279 |
| changed by the review fixes only | 4: `deployment_opportunities` (`member_components`: the 12,876 pre-fix `"nan"` strings cleared, 0 left), `deployment_candidates` (`recommendation_json`), `demo_top10_candidate_sites.csv` and `facility_candidate_examples.csv` (`reasons`, `absence_note` wording) |
| runtime fields only | 3: `nhanes_run_metadata.csv`, `nhanes_unsupervised_summary.csv`, `test6_clinic_shuffle_checks.csv` |
| changed by concurrent work | 1: the canonical `phenotype_signatures` (441 -> 28,516 rows), because other contributors wrote `phenotype_signatures__geo_cohorts` and `phenotype_signatures__mapmecfs_omics` into the same tree during the run and `canonical_unions` unions every partition |
| new, not pipeline outputs | 29 tables of that concurrent work (GEO cohorts, person-linked mapMECFS omics, published device evidence) |

Every results/*.md write-up the pipeline writes (including `GEOGRAPHY_RESULTS.md`, previously synced by hand to the
generator text) regenerated identical apart from generation times. 48 steps ended ok; the final `validate` step failed
on the concurrent work's registry entry `published_device_evidence` (its `audit` field names `DATA_AUDIT.md` without
the `data/raw/<source_id>/` prefix), not on a pipeline output. The registry merge of this run also picked up the two
new sources of that work, so `SOURCE_REGISTRY.yaml` and Figure 2 now list 29 sources. From the next run on, each
`step_end` line records `max_rss_mb` (peak RSS of the step's largest process).

### 5.8 Integration of the labs, omics and device analyses (2026-09-24, after the report review)

The modules behind `results/LINKED_OMICS_RESULTS.md`, the GEO cohort tables, `results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md`,
`results/FM_THERMOGRAPHY_RESULTS.md`, `results/LABS_VS_DIAGNOSIS_RESULTS.md`, the three lab-dataset results files, the
two published-evidence tables, the two discovery tables and the All of Us lab-concept counts had been run as standalone
commands. They are now the 20 steps listed in section 1. Offline switches: `ingest_muscle_me_charlton` and
`ingest_fm_thermography` use `skip_download=True`, the two published-evidence steps `fetch=False` (verify against the
saved source texts), `labs_dataset_discovery` `download=False` (inspect the files on disk), and `labs_aou_lab_concepts`
`fetch=False` (read the saved Data Browser JSON; new rows of `controlled_cohort_counts.csv` are appended only when
absent). `data/processed` and `results/tables` were copied to `data/interim/repro_backup/pre_integration_0924b` first.

Run `20260924T215649Z` (`--offline --workers 16`, the 20 new steps plus `digital_person`): 24.2 min wall time, all ok;
`omics_geo_cohorts` and `omics_person_linked` (16 workers each) ran alone first, so the wall time is close to the sum of
step durations. Then run `20260924T222227Z` (`datagov`, `scoring_report`, `figures_data`, `canonical_unions`,
`registry_merge`, `figures_schematics`, `build_db`, `validate`; 3.3 min, ok) and, after two online `datagov` runs for the
new sources' catalogue queries, `20260924T223309Z` (`datagov`, `registry_merge`, `figures_schematics`, `build_db`,
`validate`; 68 s, ok).

| step | s | workers | peak RSS (MB) |
|---|---|---|---|
| `omics_geo_cohorts` | 149 | 16 | 5,973 |
| `omics_person_linked` | 353 | 16 | 1,039 |
| `labs_heds_hsd_olink_cinquina2026` | 485 | 8 | 463 |
| `labs_vs_diagnosis` | 357 | 16 | 4,386 |
| `labs_klein2023_mylc` | 119 | 8 | 371 |
| `wearable_fm_thermography` | 56 | 8 | 370 |
| `labs_dataset_discovery` | 18 | 1 | 239 |
| `labs_appelman_metabolomics` | 15 | 16 | 290 |
| `wearable_muscle_me_steps` | 13 | 1 | 364 |
| `labs_endo_arg1_repod` | 6 | 1 | 275 |
| the other 10 new steps | 1 each | 1 | 30-541 |
| `digital_person` (now 11 datasets) | 16 | 1 | 10,083 |
| `scoring_report` (re-run) | 125 | 6 | 3,114 |
| `build_db` | 60 | 1 | 4,809 |

`reproduce compare` against the copy (`results/pipeline_runs/20260924_integration_vs_pre_integration.csv`, window from
21:56 UTC): 422 tables identical and 3 numerically equal (`phenotype_signatures__nhanes_labs` and two
`labs_dx_nhanes_*.csv`, float noise within 1e-9): every output the new steps write reproduced what the standalone runs
had written. 12 changed, all expected: the canonical `participants`, `participant_labs`, `participant_mortality` and
`phenotype_signatures` unions and `digital_person_union_checks.csv` (the lab-dataset and NHANES 2003-2006 partitions
now enter them); `deployment_candidates.recommendation_json` (its phenotype block lists the new signature datasets; no
rank, score or site changed); `data_gov_discovery_log` / `data_gov_search_results` (Data.gov rules for the 10 new
sources); `device_dataset_candidates.csv` (the `access_detail` / `caveats` text no longer says to run with
`uv run --with`); three run-metadata files (runtime fields). `controlled_cohort_counts.csv` is byte-identical (the
comparison tool could not sniff its CSV dialect and reports an error). No full 75-step run has been made since.

## 6. What needs the network on a cold cache

Without `data/raw` and `data/_http_cache`, `uv run measure-it pipeline` (no `--offline`) downloads and queries the
services below. **This cold path has not been run end to end from an empty checkout**: the raw files and cached
responses were fetched step by step while the build was developed, and every run recorded on this page is offline. Its
duration and memory peak are therefore unknown; keep `--workers 16` (the default is 1). A Data.gov outage or rate
limit no longer stops it: since 2026-09-24 the `datagov` step records the source as `blocked` (every query failed) or
`partial` and exits 0, so `registry_merge`, `build_db` and `validate` still run.

| step | network use (cold) |
|---|---|
| `geography_crosswalk` | Census cartographic boundaries 2024, 2020 ZCTA-county relationship file, Gazetteer |
| `ontology` | Mondo release (GitHub: OBO, OWL header, SSSOM; latest-release API), EBI OLS4 API (3,239 cached calls), Monarch API, CMS ICD-10-CM files |
| `ingest_nhanes` | CDC NHANES XPT files incl. PAXMIN 8-9 GB per cycle (segmented range download), documentation pages, NCHS linked mortality |
| `ingest_stanford_wearables` | Stanford wearable releases, about 6.0 GB of ZIP files |
| `ingest_nppes` | CMS NPPES monthly V.2 ZIP (~1 GB, 11.7 GB CSV), deactivation report, NUCC taxonomy CSV, listing pages; V.1 layout via HTTP range reads |
| `ingest_clinicaltrials` | ClinicalTrials.gov API v2 (version, search areas, one query per condition term, record batches) |
| `ingest_nih_reporter` | NIH RePORTER API v2 (175 cached POST pages), data dictionary |
| `ingest_openfda` | openFDA device APIs (588 cached calls; within the anonymous 1,000/day limit) + classification bulk file |
| `ingest_cdc_places`, `ingest_cdc_long_covid`, `ingest_cdc_lyme` | data.cdc.gov Socrata downloads and metadata, www.cdc.gov Lyme county file, Census PEP denominators |
| `ingest_census_acs`, `ingest_svi`, `ingest_rucc`, `ingest_hrsa`, `ingest_cms_mmd` | Census ACS summary files and API probe, CDC/ATSDR SVI CSVs, USDA ERS RUCC, HRSA Data Warehouse files, data.cms.gov MMD API |
| `ingest_mapmecfs` | Nature supplementary data, GEO series matrices, mapmecfs.org probes, E-utilities, Pennsieve API |
| `datagov` | catalog.data.gov search (10 s crawl delay), api.gsa.gov gateway probes |
| `omics_open_targets`, `omics_gwas_catalog`, `omics_geo_sra` | Open Targets GraphQL (295 calls), GWAS Catalog REST v2, NCBI E-utilities |
| `omics_reference_data`, `omics_coherence` | Reactome release-97 archive files, HGNC complete set, Reactome AnalysisService (cross-check, 3 calls) |
| `ingest_nhanes_labs_extended`, `ingest_nhanes_2003_2006` | CDC NHANES XPT files (2011-2014 laboratory catalogue; 2003-2006 cycles incl. PAXRAW), NCHS linked mortality |
| `omics_geo_cohorts` | GEO FTP series matrices / supplementary files, GPL annotation, NCBI gene_info |
| `ingest_muscle_me_charlton`, `ingest_fm_thermography`, `ingest_appelman_lc_pem`, `labs_klein2023_mylc`, `labs_endo_arg1_repod`, `labs_heds_hsd_olink_cinquina2026` | journal / repository supplements (Nature Communications, PLOS ONE figshare, Nature, RepOD Dataverse, Europe PMC) |
| `measurement_published_evidence`, `labs_published_evidence` | Europe PMC abstracts / open full texts, PDF supplements (online runs only; offline they verify the saved texts) |
| `measurement_device_dataset_discovery`, `labs_dataset_discovery` | candidate files and landing pages (figshare, Zenodo, Dryad, Pennsieve, OSF, Dataverse, ...); the lab module skips the two records that appear to hold patient identifiers (`PRIVACY_EXCLUDED`) |
| `labs_aou_lab_concepts` | All of Us public Data Browser API (online runs only) |

All other steps read only data/processed, configs and results/tables. The HTTP layer throttles per host and caches every
response, so a second online run does not query the services again.

## 7. API keys

None are required. Optional: `DATA_GOV_API_KEY` (or `API_DATA_GOV_KEY`, `DATAGOV_API_KEY`) makes the Data.gov step try the
api.gsa.gov gateway first. Do not set it for an offline reproduction: the key is part of the HTTP cache key, so the
cached keyless searches would not be found.

## 8. Inputs that do not reproduce on a cold cache

The raw files and cached responses in this checkout are the reproducible inputs (MANIFEST.json records url, bytes,
sha256 and retrieval time of every raw file; each cached response records its fetch time). Fetching again later gets
different data from sources that change:

| source | how it changes | what is pinned here |
|---|---|---|
| ClinicalTrials.gov API v2 | registry updated daily by sponsors | saved API pages (dataTimestamp in `clinical_trials.source_version`) |
| NIH RePORTER | weekly refresh; FY2026 still open at retrieval | saved POST pages |
| openFDA device | new 510(k)/PMA decisions weekly | cached responses + classification bulk file |
| Open Targets Platform | quarterly releases (26.06 used) | cached GraphQL responses |
| GWAS Catalog | releases every few weeks (2026-09-13 used) | cached REST responses |
| NCBI GEO / SRA | continuous deposits | cached E-utilities responses |
| Mondo, OLS4, Monarch | monthly Mondo releases; live APIs | versioned Mondo files in data/raw, cached API responses; the "latest release" check is informational |
| HRSA, NPPES, CDC PLACES / SVI / HPS / Lyme, CMS MMD, ACS | daily (HRSA), monthly (NPPES), annual releases | raw files in data/raw with their MANIFEST.json |
| Data.gov catalog | harvested continuously; relevance ranking changes | cached search responses |
| Reactome | quarterly releases | pinned: release-97 archive URLs (`https://download.reactome.org/97/`), sha256 in MANIFEST.json; the AnalysisService cross-check is live (current release; cached) |
| HGNC complete set | continuous, no release numbers | the file in `data/raw/hgnc` (Last-Modified 18 Sep 2026, sha256 in MANIFEST.json); a cold checkout gets the then-current file |
| Stanford, NHANES, MapMECFS publication deposits | static releases | raw files with sha256 |

## 9. Re-running parts

* One step: `uv run measure-it pipeline --offline --only <step>` (its inputs must exist); a step and everything after
  it: `--from <step>`. Steps do not check whether their inputs are newer than their outputs; after changing a
  shared input (a config, the crosswalk, an ingestion module), re-run from the first consumer (`--list` shows them),
  or the whole pipeline, as the S cases in 5.2 show.
* `datagov` reads SOURCE_REGISTRY.yaml, which `registry_merge` writes at the end of the run (and `registry_merge`
  depends on `datagov`, which writes its own registry fragment), so a newly added source reaches the Data.gov
  discovery tables only on the next run. With an unchanged source list this is a fixed point (both runs of section 5.1
  are identical).
* Scratch runs of single modules (`uv run python -m measure_it.<package>.<module>`) are equivalent to the step with
  the module's default arguments; the pipeline passes only the worker counts and the offline switches listed by
  `--list`.
* The figures are two DAG steps: `figures_data` (SPEC Figures 3-9, results/maps, the Figure 3-9 block of
  results/figures/CAPTIONS.md; after `scoring_report`) and `figures_schematics` (Figures 1, 2 and 10; after
  `registry_merge`, because Figure 2 is drawn from the merged SOURCE_REGISTRY.yaml). The draft figures under
  results/figures/drafts and results/maps/drafts are written by the analysis steps themselves.
