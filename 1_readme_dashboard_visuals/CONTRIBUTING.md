# Contributing

Thank you for helping. This project is an alpha research prototype: its value rests on every number being traceable to
a public source and on a small set of scientific guardrails being enforced in code. Most of this page is about keeping
those two properties when you add something.

Before a large change, open an issue (or a draft pull request) so the approach can be agreed first. Participation is
governed by the [code of conduct](CODE_OF_CONDUCT.md); privacy and security problems go through
[SECURITY.md](SECURITY.md), never a public issue.

## Setting up

```bash
git clone <this repository> && cd measure-it-public
uv sync --frozen                                   # Python >= 3.12, the environment pinned in uv.lock
uv run pytest -m "not network and not data"        # the unit tests: need no data, run in CI (about 20 s)
uv run measure-it pipeline --list                  # the 75-step DAG
```

A clone holds code, configs, audits, manifests and result write-ups, but **no data**. To work on anything that reads
`data/processed`, build it first with `uv run measure-it pipeline --workers 16` (network, about 40 GB of downloads; see
README "Quickstart" and [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) section 6). After that,
`--offline` rebuilds from the local copies without the network. To run your own measurements through the engine, see
[docs/BRING_YOUR_OWN_DATA.md](docs/BRING_YOUR_OWN_DATA.md) and [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md).

Tests are marked: `@pytest.mark.data` needs the processed tables, `@pytest.mark.network` calls a live API. CI runs
`-m "not network and not data"`, so **any test that reads `data/processed` must carry `@pytest.mark.data`**.

## Guardrails every contribution keeps

These come from [docs/CONVENTIONS.md](docs/CONVENTIONS.md) section 3 and are checked by `uv run measure-it validate`
and the tests. A pull request that breaks one will not be merged, however useful the feature.

1. **Four data layers stay separate.** Person-level rows (layer 1) are never joined to places (layer 3) or facilities
   (layer 4) as if they were the same people; geographic joins are ecological. Participant ids are namespaced
   `<dataset_id>:<native id>` and datasets are never pooled into one person.
2. **No inferred locations.** Never infer where a participant lives from a clinical or research dataset.
3. **"Participant-linked" omics only when the source publishes shared participant ids** and the overlap was counted;
   otherwise it is condition-level molecular enrichment.
4. **No fabricated data.** An unreachable, gated or incomplete source is recorded as `blocked` / `partial` in its
   registry entry and DATA_AUDIT.md. Missing values stay missing: tools return `UNKNOWN / NOT AVAILABLE` with a reason.
   No number produced by a language model is data.
5. **No downscaling.** A state value attached to a county keeps `source_geographic_resolution="state"` and is labelled
   inherited context. Proxies carry a burden evidence level (A direct, B closely matching coded condition, C
   symptom/comorbidity proxy, D none).
6. **Signals are not evidence of validity.** An FDA record is a deployment-readiness signal; a trial or grant is a
   research-activity signal; neither shows that a measurement works or diagnoses anything.
7. **Scores keep their components**, and rankings are reported with their sensitivity analyses.
8. **Product language.** Use "candidate", "evidence-supported", "measurable phenotype", "deployment opportunity",
   "measurement desert", "research readiness". Do not write "diagnosed by AI", "proves", "confirms mechanism",
   "patient has", "definitive biomarker", "best clinic" or "optimal treatment" (`validate` scans for these).
9. **Analysis plans come first.** A new group comparison or model gets a `docs/ANALYSIS_PLAN_<NAME>.md` written and
   committed *before* the comparison is run; later changes are recorded as dated deviations in that file.
10. **Privacy.** Only public, de-identified or aggregate data enter this repository, and never raw data or processed
    tables (they are rebuilt by the pipeline and git-ignored). Private data stay outside the repository; a private
    adapter declares `requires_private_data = True` and only group-level results may reach public outputs. Nothing
    that contains direct identifiers is downloaded, even from a public repository.

## Adding a data source

A source module delivers the six things listed in [docs/CONVENTIONS.md](docs/CONVENTIONS.md) section 5:

1. **Module** `src/measure_it/<package>/<module>.py` (`ingestion/` for most sources) with an idempotent `run()` (or
   `build()`) entry point and `if __name__ == "__main__": run()`. Use the shared infrastructure, do not re-implement
   it: `measure_it.http.get_json/post_json/get/post` for API calls (disk-cached, throttled, offline-aware),
   `measure_it.download.download_file` for files (writes `MANIFEST.json` with url, bytes, sha256 and retrieval time),
   `measure_it.provenance.add_provenance` and `measure_it.store.write_table` for tables (the provenance columns are
   enforced), `measure_it.geography.crosswalk.zip_to_geo` for ZIP codes. Never call `requests` directly for API traffic:
   the offline mode and the cache depend on it.
2. **Raw files** under `data/raw/<source_id>/` (git-ignored) with their `MANIFEST.json` (tracked).
3. **`data/raw/<source_id>/DATA_AUDIT.md`** from [docs/templates/DATA_AUDIT_TEMPLATE.md](docs/templates/DATA_AUDIT_TEMPLATE.md),
   with *measured* sample sizes and missingness, the licence quoted from the provider with its URL, and the access
   conditions.
4. **`data/raw/<source_id>/registry_entry.yaml`** written with `measure_it.registry.write_registry_entry`; then
   `uv run measure-it merge-registry` regenerates `SOURCE_REGISTRY.yaml` (never edit that file by hand).
5. **Processed tables** via `write_table` into `data/processed/<table>.parquet` (canonical) or
   `<table>__<source>.parquet` (a partition that `canonical_unions` merges); naming rules in CONVENTIONS section 2.
   Rows that recommendations cite need an `object_id` (CONVENTIONS section 7).
6. **Tests** `tests/test_<module>.py`: pure-function unit tests (unmarked, run in CI) and checks on the processed
   outputs marked `@pytest.mark.data`.

Then add a `Step` to `measure_it.pipeline.STEPS` (entry point, dependencies, worker parameter) rather than a separate
script, run it (`uv run measure-it pipeline --only <step>`, then `--from` the first downstream step that reads it), and
run `uv run measure-it validate`. Add the source to [docs/DATA_LICENSES.md](docs/DATA_LICENSES.md).

**Before you write the module, check the licence.** Sources must be downloadable without an agreement that forbids
redistributing derived tables, and must contain no direct identifiers. When a provider's download redirects to a
pre-signed cloud URL, strip its query string before committing `MANIFEST.json` (the CI hygiene job rejects signed
URLs and absolute local paths).

## Adding a measurement adapter

Only the first link of the chain (person-level phenotype) depends on the measurement modality, and it sits behind one
interface, `measure_it.measurements.adapter.MeasurementAdapter` (`preprocess`, `embed`, `phenotype_score`,
`describe_measurement`). The full contract, the wearable adapters and the capillaroscopy stub are in
[docs/ADAPTERS.md](docs/ADAPTERS.md); README "Adding a measurement modality" has the five steps. In short: implement
the four methods; use `measurement_ids` that exist in `configs/measurements.yaml`; give each phenotype an ontology
hook (an HPO term) or its condition links resolve to UNKNOWN; register it with
`measure_it.measurements.adapters.register_adapter`; add a pipeline step that writes
`participant_adapter_scores__<dataset>`; and make `uv run pytest tests/test_adapters.py` pass (the contract test runs a
synthetic adapter through every downstream function). For private data, follow
[docs/BRING_YOUR_OWN_DATA.md](docs/BRING_YOUR_OWN_DATA.md) and keep the data and person-level outputs outside the
repository.

## Adding a condition

1. Add an entry to `configs/conditions.yaml`: `id` (snake_case slug), `preferred_name`, `aliases`, `search_terms`,
   `parent`, `flags`. Search terms must be specific: never a bare acronym that is ambiguous in free text (the file's
   header explains why "POTS" and "EDS" are not searched). Do **not** add ontology identifiers there: the `ontology`
   step resolves them from Mondo, OLS4, Monarch and ICD-10-CM and records where each came from.
2. Add its specialties to `condition_specialties` in `configs/relevance.yaml`, with a rationale.
3. Re-run from the ontology step **with the network** (the new terms are not in the HTTP cache, so an `--offline` run
   stops with `OfflineCacheMiss`): `uv run measure-it pipeline --from ontology --workers 16`. Search-driven sources
   (ClinicalTrials.gov, RePORTER, Open Targets, GWAS Catalog, GEO) then query the new terms.
4. Check the condition's burden rows: most conditions have no direct county measure; say which evidence level (A-D)
   each burden value has, and leave it D (null) rather than inventing a proxy.
5. Update the tests that enumerate conditions and run `uv run measure-it validate`.

## Pull requests

* Keep a pull request to one purpose. Describe what changed and, for data or scoring changes, which pipeline steps you
  re-ran and which reported numbers moved (`uv run python -m measure_it.reproduce compare` produces the table).
* `uv run pytest -m "not network and not data"` must pass (CI runs it together with a hygiene check and a lint for
  syntax errors and undefined names). On a built checkout, `uv run pytest -m "not network"` and
  `uv run measure-it validate` should pass too.
* Changes to shared infrastructure (`src/measure_it/{config,http,download,provenance,store,registry}.py`,
  `measurements/adapter.py`, `geography/crosswalk.py`, `pyproject.toml`, `uv.lock`) need a sentence explaining why;
  add dependencies with `uv add`, which updates `uv.lock`.
* Do not commit data payloads, notebooks with outputs containing per-person rows, credentials, absolute local paths or
  generated caches.
* By submitting a contribution you agree that it is licensed under the repository's [licence](LICENSE) (MIT). Data
  you add remain under their providers' terms ([docs/DATA_LICENSES.md](docs/DATA_LICENSES.md)).
