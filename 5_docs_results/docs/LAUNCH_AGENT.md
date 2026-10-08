# Launch agent: any dataset -> harmonized layers -> candidate launch sites

`measure-it launch` takes a dataset in (almost) any format, public or private (a device study such as MAESTRO),
harmonizes it to this engine's schema (the `measure-it byod` folder, [BRING_YOUR_OWN_DATA.md](BRING_YOUR_OWN_DATA.md))
and to the shape of the All of Us Curated Data Repository (OMOP CDM v5.4), runs the engine's layers on it, and returns
a ranked, traceable list of candidate launch counties for the measurement device, with the clinicians and sites to
contact. The binding interface is [AGENT_CONTRACT.md](AGENT_CONTRACT.md); the manual version of the same run is
[MAESTRO_RUNBOOK.md](MAESTRO_RUNBOOK.md). Code: `src/measure_it/agent/`.

Every ranked county is a candidate deployment opportunity for pilot evaluation, not a validated diagnostic pathway.

## 1. Use

```bash
# 1. propose a mapping (nothing is written to the engine; no model is contacted)
uv run measure-it launch map path/to/export_folder --condition long_covid --measurement nailfold_capillaroscopy \
    --out launch_runs/maestro_v1
# 2. review and edit launch_runs/maestro_v1/mapping.yaml (roles, codes, labels, comparator, subgroup_analysis,
#    and, for private data, the data_use block), then approve and apply it
uv run measure-it launch apply launch_runs/maestro_v1 --approve
# 3. run everything (re-uses the approved mapping of the run folder)
uv run measure-it launch run path/to/export_folder --condition long_covid --measurement nailfold_capillaroscopy \
    --out launch_runs/maestro_v1
# re-render the report from launch.json
uv run measure-it launch report launch_runs/maestro_v1
```

`launch run ... --approve` does steps 1-3 in one go (the approval is recorded as `approved_by` / `approved_at`).
`launch profile <path>` writes only `profile.json`. Options of `launch run`:

| option | meaning |
|---|---|
| `--out DIR` | run folder (default `launch_runs/<stamp>_<condition>_<measurement>`); re-use it to resume after editing `mapping.yaml` |
| `--llm none\|claude` | `none` (default): no model is contacted. `claude`: metadata-only mapping proposals for columns the rules could not map (needs `ANTHROPIC_API_KEY` and the optional `anthropic` package) |
| `--approve` | approve `mapping.yaml` |
| `--device-positive COLUMN[>=X]` | pre-specify the device-positive subgroup from a column mapped to role `device`: a 1/0 flag (`cap_positive`) or a threshold (`'cap_abnormality_score>=2'`); written into `mapping.yaml` as `subgroup_analysis` before approval and locked by `byod subgroup` before any outcome. Without it (or a hand-written `subgroup_analysis`) the subgroup, similar-people and subgroup-burden steps are skipped and the launch ranking falls back to the engine ranking |
| `--public` | the input is a public, already de-identified dataset (data use is stated as public) |
| `--synthetic` | the input is generated test data: the dataset is `synthetic: true, demo: true` |
| `--confirm-deidentified --confirm-use-permitted --data-use-statement TEXT [--irb-or-dua REF]` | the data owner's confirmation for private data (written into `mapping.yaml`'s `data_use` block) |
| `--keep` | leave the dataset ingested at the end (default: `byod remove`) |
| `--demo-rankings` | let a demo record enter the rankings for this run only (`byod deploy --demo`; the teardown restores them; not with `--keep`) |
| `--top N`, `--per-county N` | counties ranked, clinicians per county in the outreach workbook |
| `--fast` | few CV repeats / permutations / bootstraps (smoke runs only; recorded in the locked plan) |
| `--no-outreach`, `--remap`, `--verbose` | skip the workbook; discard the run folder's mapping and propose again; echo the engine log |

Exit codes: 0 completed, 1 completed with a failed step, 2 refused (bad arguments), 3 stopped (approval or data use).

## 2. The step plan

`<out>/plan.json` records every step as `done` (with its outputs), `skipped` + why, `failed` + error, or `stopped` +
message, with timings. What the data hold decides the steps:

| step | what | skipped when |
|---|---|---|
| profile | read any format (CSV/TSV/Excel/Parquet/JSON/SAS/SPSS/Stata/REDCap/OMOP/FHIR/ZIP) -> `profile.json` | the input already is a BYOD folder |
| map | rule / dictionary / fuzzy (optionally Claude) proposals -> `mapping.yaml`; an existing `mapping.yaml` in the run folder is never overwritten (`--remap` discards it) | BYOD folder input |
| approve | the approval gate | **stops** the run when the mapping is not approved |
| apply | `<out>/byod/` (identifier and geography columns dropped, dates -> `day_index`, ages -> 10-year bands), `<out>/omop/` | BYOD folder input |
| omop_export | All of Us-shaped OMOP tables (from apply, or from a BYOD folder input) | |
| byod_validate | `measure-it byod validate` on the written folder -> `byod_validation.json` | |
| data_use_gate | private data must confirm `data_use.deidentified`, `data_use.use_permitted`, `data_use.statement` | **stops** before any ingest and names the fields; also stops when the dataset id is already ingested (a run never touches data it did not create) |
| byod_ingest | person-layer partitions (shared processed tables) | validation failed |
| byod_evaluate | locked primary analysis -> performance record | no label |
| byod_subgroup | locked subgroup plan -> computable phenotype + strata fractions | no `subgroup_analysis.device_positive` |
| byod_similar | NHANES "similar people" estimates + All of Us / N3C / clinic query packs | no phenotype |
| query_pack_check | the stage-D All of Us pack transpiled to DuckDB and run against the local export | no phenotype or export; `sqlglot` not installed (run with `uv run --with sqlglot ...`) |
| byod_deploy | the performance record enters the metric link and the evidence-weighted ranking | no record; a demo record (unless `--demo-rankings`) |
| geo_subgroup_burden | device-subgroup people per county (in memory; written only to `<out>/geo_subgroup_burden.csv`) | no phenotype |
| launch_ranking | section 3 | |
| outreach | outreach workbook under `<out>/outreach/` | `--no-outreach` |
| teardown | engine results copied to `<out>/engine_results/`, then `byod remove <id>` and a check that nothing of the dataset is left | nothing ingested; `--keep` |
| report | `LAUNCH_REPORT.md` + `launch.json` | never (also after a stop) |

The teardown runs even when an earlier step failed, so a launch run never silently leaves user data in the engine.

## 3. "Greatest impact"

For the condition x measurement:

* **(a) engine ranking**: `deployment_opportunities` at county level under `evidence_weighted` when the measurement has
  a known performance record for the condition in the rankings (after `byod deploy`, or a public record), else under
  `equal`. Components: burden, vulnerability, diagnostic desert, measurement-specific clinic capacity (activity where
  the measurement has a dedicated billing code, implementer density otherwise; nailfold capillaroscopy has no code),
  research readiness, and the device evidence / expected yield when known; the engine's Monte Carlo rank interval.
* **(b) device-subgroup re-score**, when a device subgroup exists: the same engine combination
  (`opportunity.make_combo` + `opportunity.evaluate`) with the subgroup condition's burden replaced by the estimated
  device-subgroup rate per 100k adults per county (`geography.subgroup`: base burden x the cohort's device-positive
  fraction by age band x sex), under the same weight set; the estimated number of adults in the subgroup is reported
  next to it; an equal-weight Monte Carlo interval draws weights and subgroup rates.

The launch ranking is (b) where a subgroup exists, else (a); the report names which and lists the counties entering
and leaving the top N. Conditions outside the scored grid (e.g. `lyme_disease`, `ptlds`) use the BURDEN-ONLY fallback
of `measure_it.outreach` (and with a subgroup, a subgroup-burden-only ranking), labelled as such.

For the top counties: billing clinicians (Medicare FFS, `geo_measurement_capacity`), measurement experience sites,
candidate facilities (clinic candidates with a measurement implementer specialty). Names, addresses and NPIs appear
only in the outreach workbook.

## 4. Privacy and language-model boundaries

* A model may only **propose** how a column maps to a role and a code, and only with `--llm claude`. It sees the
  profile's metadata (names, types, dictionary labels and choices, units, numeric min/median/max, categorical levels
  with counts >= 11, share missing), never row-level records, participant ids, free text, dates or any column flagged
  as identifying or geographic. Every request is logged to `<out>/llm_requests.jsonl` before it is sent. Every
  proposal is validated locally (code format, code in a local vocabulary or config, unit plausible) and needs the
  human approval of `mapping.yaml`.
* **No language-model numbers.** Every value, score and ranking comes from the data and the engine; the report is a
  deterministic template.
* The data owner confirms de-identification and permission (`data_use`); the engine never asserts it for private data.
  `--public` and `--synthetic` state it for public and generated data.
* Identifier and geography columns are never written; exact dates become per-participant `day_index`; exact ages become
  10-year bands; the BYOD validator must pass before anything is ingested.
* Everything under the run folder is local-only (`.gitignore` of `*`): the BYOD folder, the OMOP export, query-pack
  results, the outreach workbook. Do not commit or share it without the data owner's agreement.
* Synthetic / demo datasets keep `demo: true`; their record never enters the default rankings
  (`MEASURE_IT_BYOD_DEMO` semantics, docs/BRING_YOUR_OWN_DATA.md section 7).

## 5. Run folder

```text
launch_runs/<run>/
  .gitignore                 *
  plan.json                  steps, statuses, reasons, timings (earlier plans of the folder under logs/)
  profile.json, mapping.yaml, mapping_applied.yaml, apply_report.json
  byod/                      the BYOD folder (manifest.yaml, participants.csv, device_*.csv, survey_*.csv, ...)
  byod_validation.json
  omop/                      person, observation_period, condition_occurrence, measurement, drug_exposure,
                             observation, concept (.parquet), README.md, query_pack_check.json
  engine_results/            copy of results/byod/<id>/ (plan locks, EVALUATION.md, SUBGROUP.md, phenotype, packs)
  geo_subgroup_burden.csv    device-subgroup people per county / state / nation (aggregate)
  launch_ranking.json        the two rankings with every component
  outreach/                  outreach workbook (named clinicians)
  launch.json, LAUNCH_REPORT.md
  logs/engine.log
```

## 6. Example (synthetic)

`tests/test_agent_launch.py` builds a SYNTHETIC, MAESTRO-like REDCap export with a data dictionary (cohorts Long COVID /
chronic Lyme / acute Lyme / healthy, capillaroscopy features with a reader's device-positive call, COMPASS-31 / FSS-9 /
DSQ-PEM, Nightingale-like NMR columns, a ZIP code, initials, dates and free-text notes that must be dropped) and runs:

```python
from measure_it.agent.orchestrate import run_launch
run_launch("maestro_redcap/", "long_covid", "nailfold_capillaroscopy", out="launch_runs/synthetic",
           stop_after="map")                     # profile + mapping.yaml; the test then edits it like a reviewer
run_launch("maestro_redcap/", "long_covid", "nailfold_capillaroscopy", out="launch_runs/synthetic",
           approve=True, synthetic=True, fast=True)
```

The run validates the BYOD folder, exports OMOP, ingests, evaluates the capillaroscopy record, finds the planted
device subgroup (COMPASS-31 vasomotor, GlycA and Raynaud history), estimates the subgroup per county, ranks counties
under the engine (equal weights: the demo record stays out of the rankings) and under the subgroup re-score, writes the
outreach workbook and the report, and removes the dataset again. Its numbers describe generated data only.
