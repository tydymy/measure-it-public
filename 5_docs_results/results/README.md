# results/

Git keeps the write-ups and the few files no pipeline step can rebuild. The figures the write-ups show are published
copies in [`../visuals/`](../visuals/) (refresh them from `figures/` after a rebuild). Everything else in this folder is written by
`uv run measure-it pipeline` (README "Building the data") and is not in git.

| path | in git | what it is |
|---|---|---|
| `FINAL_REPORT.md`, `*_RESULTS.md` and the other `*.md` | yes | the write-ups of the reference build (2026-09-25 run `20260925T194852Z` and later single steps) |
| `tables/labs_dx_plan_lock.json` | yes | the lab analysis plan's sha256 lock, written before any outcome was computed (`measure_it.labs.plan`) |
| `tables/measurement_precision_*_round{1,2}.csv`, `tables/measurement_pro_filter_{sample,judgements}.csv` | yes | manual precision judgements of the measurement text patterns |
| `tables/controlled_cohort_counts.csv` | yes | All of Us Data Browser counts, appended by hand-run queries (`measure_it.labs.aou_lab_concepts --append`) |
| `tables/*_run1_before_rule_fixes.*`, `tables/_invalid_runs/`, `tables/number_audit.csv` | yes | archived runs and the number audit the write-ups cite |
| `tables/` (everything else), `figures/`, `maps/`, `pipeline_runs/`, `*.json` | no | pipeline outputs: the CSV tables the write-ups cite by path, interactive maps, run logs, example recommendations |

The write-ups cite tables as `results/tables/<name>.csv`. On a fresh clone those paths exist after the pipeline has
run; the reference build's copies are in the `reference-results` archive attached to the GitHub release, which unpacks
at the project root.
