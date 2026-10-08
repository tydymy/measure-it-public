# Walkthrough: from install to your own data on the map (and back)

Every command below was run on 2026-09-28 on the build host (Linux aarch64, a checkout holding `data/processed`), in
this order; the outputs are real, abbreviated with `...`. The dataset added in steps 5-9 is **SYNTHETIC** (generated,
no real person) and tagged `demo`, so it only reaches a ranking when asked to (`--demo`), and step 10 removes it. The
contract behind the `byod` commands: [BRING_YOUR_OWN_DATA.md](BRING_YOUR_OWN_DATA.md).

## 0. Install, with the data

**uv on a checkout.** Python >= 3.12 and [uv](https://docs.astral.sh/uv/). The code is in git; the data are not
(README "Quickstart"): build them with the pipeline, or unpack a data bundle at the project root
([DEPLOYMENT.md](DEPLOYMENT.md) section 3).

```bash
uv sync --frozen
uv run measure-it bundle --verify dist/bundle      # only if you use an unpacked bundle: re-hash every file
```

```text
Checked 138 packages in 1ms
{
 "ok": true,
 "n_files": 1045,
 "missing": [], "size_mismatch": [], "sha256_mismatch": [],
 ...
}
```

**Docker + data bundle** (read-only containers for the API, the dashboard and the MCP server; `byod` itself runs with
uv on a checkout, because the containers mount the bundle read-only):

```bash
docker compose up -d --build
docker compose ps
curl -s http://127.0.0.1:8000/health
docker compose down
```

```text
 Container measure-it-public-api-1 Started
 Container measure-it-public-dashboard-1 Started
 Container measure-it-public-mcp-1 Started
NAME                            ...   STATUS                        PORTS
measure-it-public-api-1         ...   Up About a minute (healthy)   ..., 127.0.0.1:8000->8000/tcp, ...
measure-it-public-dashboard-1   ...   Up About a minute (healthy)   ..., 127.0.0.1:8501->8501/tcp, ...
measure-it-public-mcp-1         ...   Up About a minute (healthy)   ..., 127.0.0.1:8765->8765/tcp
{"status":"ok","reason":null,"service":"Measure It to Cure It: public-data API","version":"0.1.0","n_tools":18,
 "processed_tables":159,...,"offline":true,"unknown_sentinel":"UNKNOWN / NOT AVAILABLE"}
```

Check the build before changing anything:

```bash
uv run measure-it validate
```

```text
measure-it validate (2026-09-28T19:10:09+00:00): PASS
  [PASS] provenance: 160/160 processed tables pass check_provenance; ...
  [PASS] person_layer_no_geography: 65/65 person-layer tables carry no geographic identifier column
  ...
```

## 1. Dashboard

```bash
uv run measure-it dashboard --offline --headless --port 8599      # stop with Ctrl-C or kill <pid>
```

```text
measure-it dashboard: http://localhost:8599 (pid 2124425; stop with Ctrl-C or kill 2124425)
```

Open http://localhost:8599: home, the five public-data pages, and page 6 **Your data** (empty until step 6; it says
`UNKNOWN / NOT AVAILABLE: no user dataset has been ingested`).

## 2. Ask

```bash
uv run measure-it ask "Where should we deploy autonomic testing for ME/CFS?" --top-n 3
```

```text
# Where should we deploy autonomic testing for ME/CFS?
...
- Candidate deployment opportunities for Myalgic encephalomyelitis/chronic fatigue syndrome x Clinic autonomic function
  testing (county level, 'equal' weights; 2,395 regions ranked; ...):
  1. Okmulgee County, Oklahoma: composite 0.902, rank interval 1-14 (P(top 10) 0.92); burden percentile 0.93
     (evidence level C), ... [opportunity:me_cfs|autonomic_function_testing|40111; geo:40111]
...
## Uncertainties
- Okmulgee County, Oklahoma: Measurement performance UNKNOWN / NOT AVAILABLE for this condition and measurement (no
  record in this project's results or the curated published evidence): not ranked under evidence_weighted, never
  scored 0 or 1. The equal-weight rank carries no measurement-performance information. [geo:40111]
```

No public performance record exists for autonomic testing in ME/CFS, so this combination is not ranked under the
`evidence_weighted` weights.

## 3. Trace

```bash
uv run python -m measure_it.tools trace_evidence '{"object_id": "measurement_performance:me_cfs|autonomic_function_testing"}'
```

```text
{ "tool": "trace_evidence", "status": "ok", ...
  "rows": [{"condition_id": "me_cfs", "measurement_id": "autonomic_function_testing", ...
            "records": [{"record_id": "PDE-022", "quality_tier": 3, ..., "performance_status": "UNKNOWN / NOT AVAILABLE"}],
            "quality_tier": 4, "quality_tier_label": "none", "performance_status": "UNKNOWN / NOT AVAILABLE", ...
            "performance_note": "UNKNOWN / NOT AVAILABLE: no performance record of this project or of the curated
             published evidence fo..." }] ... }
```

## 4. Make the SYNTHETIC dataset

```bash
uv run python examples/byod_synthetic/make_example.py
head -3 examples/byod_synthetic/synthetic_demo/participants.csv
head -3 examples/byod_synthetic/synthetic_demo/device_synthetic_hrv_patch.csv
```

```text
SYNTHETIC tutorial dataset written to .../examples/byod_synthetic/synthetic_demo (90 cases, 110 controls)
participant_id,me_cfs,age_band,sex
SYNTHETIC-0001,0,30-39,female
SYNTHETIC-0002,0,20-29,male
participant_id,day_index,rmssd_ms,sdnn_ms,resting_hr_bpm,upright_minutes
SYNTHETIC-0001,0,66.0,73.6,63.8,556.0
SYNTHETIC-0001,1,45.1,60.7,61.2,512.0
```

The folder holds `manifest.yaml` (`synthetic: true`, `demo: true`, label `me_cfs` -> ME/CFS, comparator healthy,
measurement class `hrv`, primary analysis = the ECG-patch block), `participants.csv`, two device files (a 14-day ECG
patch through a partner-written adapter in `adapters/hrv_patch_adapter.py`, and per-participant wrist-actigraph
summaries through the engine's `TabularFeatureAdapter`), `ehr_labs.csv` (LOINC), `ehr_diagnoses.csv` (ICD-10-CM),
`omics_synthetic_proteomics.csv` and `README_SYNTHETIC.md`.

## 5. Validate (and see a refusal)

```bash
uv run measure-it byod validate examples/byod_synthetic/synthetic_demo
```

```text
byod validate synthetic_demo: OK (0 error(s), 0 warning(s))
  summary:
    primary label: me_cfs: 90 cases, 110 controls with data in ['device_synthetic_hrv_patch']
    ...
    participants per block: {'device_synthetic_hrv_patch': 200, 'device_synthetic_wrist_actigraphy': 184,
                             'ehr_diagnoses': 139, 'ehr_labs': 200, 'omics_synthetic_proteomics': 140}
    omics overlap with other blocks (same participant_id): {'omics_synthetic_proteomics': 140}
    comparator: healthy
    demo / synthetic: True / True
  will create:
    data/processed/participants__byod_synthetic_demo.parquet (person layer; ids byod_synthetic_demo:<your participant_id>)
    data/processed/participant_conditions__byod_synthetic_demo.parquet (label rows + ICD-10-CM diagnosis rows)
    data/processed/participant_wearable_features__byod_synthetic_demo.parquet (device_synthetic_hrv_patch; class hrv; ...)
    ...
    after `deploy`: measurement_performance and the affected deployment_opportunities rows for ['me_cfs'] x
    ['autonomic_function_testing', 'wearable_autonomic_activity_monitoring'] (demo record: only with --demo)
```

A copy with a ZIP code, a date of birth and an MRN-like id is refused (exit 1), and no offending value is printed:

```bash
cp -r examples/byod_synthetic/synthetic_demo examples/byod_synthetic/bad_demo
uv run python -c "
import pandas as pd
d = 'examples/byod_synthetic/bad_demo/'
p = pd.read_csv(d + 'participants.csv', dtype=str)
p['zip'] = '02139'; p['date_of_birth'] = '1984-05-02'; p.loc[0, 'participant_id'] = 'MRN0012345'
p.to_csv(d + 'participants.csv', index=False)
"
uv run measure-it byod validate examples/byod_synthetic/bad_demo; echo "exit $?"
rm -rf examples/byod_synthetic/bad_demo
```

```text
byod validate bad_demo: REFUSED (11 error(s), 0 warning(s))
  ERROR   participants.csv: column 'zip' refused: it locates a person or a place (geographic column); person-level data never carry geography
  ERROR   participants.csv: column 'zip': 200 value(s) look like: ZIP code (e.g. <5 characters, starting '0'>); remove or transform it before supplying the data
  ERROR   participants.csv: column 'date_of_birth' refused: it identifies a person (name, date of birth, MRN, phone, e-mail, SSN or similar)
  ...
  ERROR   participants.csv: 1 participant_id value(s) look like: MRN / patient-number prefix (e.g. <10 characters, starting 'M'>): use study codes, not medical-record or other identifying numbers
  ERROR   device_synthetic_hrv_patch.csv: 1 participant_id value(s) are not in participants.csv
  ...
exit 1
```

## 6. Ingest

```bash
uv run measure-it byod ingest examples/byod_synthetic/synthetic_demo
```

```text
byod ingest byod_synthetic_demo: 6 file(s), manifest sha256 42d02d225547
  device synthetic_hrv_patch: adapter synthetic_hrv_patch_v1 (SyntheticHrvPatchAdapter), class hrv, 200 participants x 6 features -> participant_wearable_features
  device synthetic_wrist_actigraphy: adapter byod_tabular_byod_synthetic_demo_synthetic_wrist_actigraphy (TabularFeatureAdapter), class accelerometry, 184 participants x 3 features -> participant_wearable_features
  wrote 7 partition(s): participants__byod_synthetic_demo, participant_conditions__byod_synthetic_demo, participant_labs__byod_synthetic_demo, participant_omics_linked__byod_synthetic_demo, participant_wearable_features__byod_synthetic_demo, participant_wearable_daily__byod_synthetic_demo, participant_adapter_scores__byod_synthetic_demo
  rebuilt canonical unions (phenotype_signatures, ..., participants, ..., participant_omics_linked, digital_person) and the Digital Phenotype Vectors
ingested byod_synthetic_demo: next `measure-it byod evaluate examples/byod_synthetic/synthetic_demo`
```

`data/raw/byod_synthetic_demo/` now holds `byod_registry_entry.yaml`, a generated `DATA_AUDIT.md` (measured sizes and
missingness, linkage, limitations) and `input_manifest.json`, with a `.gitignore` of `*`. The dataset has its own
Digital Phenotype Vector (`multimodal_person_representation`: demographics, labs, symptoms_conditions, wearable and a
precomputed proteomics PCA block), never pooled with NHANES, Stanford or mapMECFS.

## 7. Evaluate (the plan is locked first)

```bash
uv run measure-it byod evaluate examples/byod_synthetic/synthetic_demo
```

```text
byod evaluate byod_synthetic_demo: plan locked (sha256 80e73f2648e5, locked 2026-09-28T19:11:53+00:00)
  population: 90 cases, 110 controls with data in ['device_synthetic_hrv_patch']
  primary                      AUROC 0.701 (0.620-0.771), sens 0.30 at spec 0.90, perm p 0.0050 [supported]
  combined                     AUROC 0.714 (0.642-0.783), sens 0.34 at spec 0.90, perm p 0.0099 [supported]
  block:device_synthetic_wrist_actigraphy AUROC 0.673 (0.596-0.742), sens 0.25 at spec 0.90, perm p 0.0099 [supported]
  block:ehr_diagnoses          AUROC 0.745 (0.677-0.810), sens 0.63 at spec 0.90, perm p 0.0099 [supported]
  block:ehr_labs               AUROC 0.556 (0.476-0.633), sens 0.11 at spec 0.90, perm p 0.1980 [null]
  block:omics_synthetic_proteomics AUROC 0.492 (0.409-0.578), sens 0.10 at spec 0.90, perm p 0.5842 [null]
  covariates                   AUROC 0.575 (0.507-0.642), sens 0.10 at spec 0.90, perm p 0.0594 [null]
  primary+covariates           AUROC 0.726 (0.651-0.794), sens 0.32 at spec 0.90, perm p 0.0099 [supported]
  Delta AUROC (primary+covariates - covariates) +0.152 (+0.079 to +0.225)
  record BYOD-SYNTHETIC_DEMO-PRIMARY: tier user_supplied_own_computation, class hrv, measurement-only True [demo: excluded from default rankings]
evaluated byod_synthetic_demo: next `measure-it byod deploy examples/byod_synthetic/synthetic_demo --demo`
```

Only `primary` becomes the performance record. The diagnosis-code block looks strong because R53.83 (fatigue) is
simulated as label-adjacent; the label-defining code G93.32 (ME/CFS) was removed from every model. The proteomics layer
is null. All numbers are against **healthy** controls and therefore optimistic.

**The lock holds.** Change the plan after the lock (here the specificity), re-ingest, and evaluate again:

```bash
sed -i 's/  specificity: 0.9$/  specificity: 0.95/' examples/byod_synthetic/synthetic_demo/manifest.yaml   # macOS: sed -i ''
uv run measure-it byod evaluate examples/byod_synthetic/synthetic_demo
uv run measure-it byod ingest examples/byod_synthetic/synthetic_demo
uv run measure-it byod evaluate examples/byod_synthetic/synthetic_demo
```

```text
byod evaluate refused: manifest.yaml changed since byod_synthetic_demo was ingested: re-run `measure-it byod ingest` first
...
  an earlier evaluation of this dataset was discarded: run `byod evaluate` (and `byod deploy`) again
ingested byod_synthetic_demo: next `measure-it byod evaluate examples/byod_synthetic/synthetic_demo`
byod evaluate refused: the analysis plan of byod_synthetic_demo was locked on 2026-09-28T19:11:53+00:00 (sha256 80e73f2648e5) and the manifest now describes a different plan. Evaluate under the locked plan, or re-run with --amend to record an amendment (the performance record will say the plan was amended).
```

Put the original manifest back (the generator is seeded, so it rewrites the same files) and evaluate under the locked
plan; the record now also says the dataset was evaluated once before (`results/byod/LEDGER.jsonl`):

```bash
uv run python examples/byod_synthetic/make_example.py
uv run measure-it byod ingest examples/byod_synthetic/synthetic_demo
uv run measure-it byod evaluate examples/byod_synthetic/synthetic_demo
```

```text
byod evaluate byod_synthetic_demo: plan verified (sha256 80e73f2648e5, locked 2026-09-28T19:11:53+00:00)
  population: 90 cases, 110 controls with data in ['device_synthetic_hrv_patch']
  primary                      AUROC 0.701 (0.620-0.771), sens 0.30 at spec 0.90, perm p 0.0050 [supported]
  ...
```

## 8. Deploy

Without `--demo` the SYNTHETIC record stays out of every ranking:

```bash
uv run measure-it byod deploy examples/byod_synthetic/synthetic_demo
```

```text
byod deploy byod_synthetic_demo: this is a demo record; it stays EXCLUDED from the rankings (pass --demo to include it). Rebuilding the metric link without it.
[metric_link] 21 records, 102 condition x measurement rows (16 known, 3 partial, 83 UNKNOWN) (6s)
  metric link rebuilt (102 rows); scored combinations whose performance changed: none
  guard: 73/73 person-layer tables carry no geographic identifier column; participant ids in derived tables: none
```

With `--demo`:

```bash
uv run measure-it byod deploy examples/byod_synthetic/synthetic_demo --demo
```

```text
byod deploy byod_synthetic_demo: record BYOD-SYNTHETIC_DEMO-PRIMARY enters the metric link (demo mode)
[metric_link] 22 records, 102 condition x measurement rows (18 known, 3 partial, 81 UNKNOWN) (6s)
  metric link rebuilt (102 rows); scored combinations whose performance changed: [('me_cfs', 'autonomic_function_testing')]
  re-scored me_cfs x autonomic_function_testing x county: 3144 rows
  re-scored me_cfs x autonomic_function_testing x state: 51 rows
  me_cfs x autonomic_function_testing (county): ranked under evidence_weighted 0 -> 2385 regions; best joint rank None -> 328; equal-weight ranks unchanged True
  me_cfs x autonomic_function_testing (state): ranked under evidence_weighted 0 -> 51 regions; best joint rank None -> 28; equal-weight ranks unchanged True
  guard: 73/73 person-layer tables carry no geographic identifier column; participant ids in derived tables: none
  report: .../results/byod/synthetic_demo/DEPLOY_REPORT.md
```

From `results/byod/synthetic_demo/DEPLOY_REPORT.md`:

```text
| condition | measurement | before | after | why |
| me_cfs | wearable_autonomic_activity_monitoring | known / OWN-MM-STEPS-ME / evidence 0.451 | known / OWN-MM-STEPS-ME / evidence 0.451 | the selected record stays OWN-MM-STEPS-ME (tier order: a better tier is present) |
| me_cfs | autonomic_function_testing | UNKNOWN / NOT AVAILABLE / - | known / BYOD-SYNTHETIC_DEMO-PRIMARY / evidence 0.145 | the user record is the selected record |
| me_cfs | hrv | UNKNOWN / NOT AVAILABLE / - | known / BYOD-SYNTHETIC_DEMO-PRIMARY / evidence 0.145 | the user record is the selected record |

## Recommendation under evidence_weighted: me_cfs x autonomic_function_testing (county, top 3)
1. Okmulgee County, Oklahoma: candidate partners CARDIOLOGY CLINIC OF MUSKOGEE, INC., Movement Disorder Clinic of Oklahoma PLLC, FAMILY CLINIC OF WELEETKA
2. Laurens County, South Carolina: candidate partners SELF REGIONAL HEALTHCARE, Absher Neurology, PA, The Pendergrass Family Health Center
3. Bradford County, Florida: candidate partners SIMEDHEALTH, L.L.C., University of Florida, UNIVERSITY OF FLORIDA

### me_cfs x autonomic_function_testing (county)
* regions ranked under evidence_weighted: 0 -> 2385; equal-weight ranks unchanged: True
* evidence_weighted top-10 after: Okmulgee County, Oklahoma | Laurens County, South Carolina | Bradford County, Florida | ...
* joint (region x measurement) top-25 under evidence_weighted, by measurement: before {'exercise_capacity_testing': 17,
  'wearable_autonomic_activity_monitoring': 8}, after {'exercise_capacity_testing': 17, 'wearable_autonomic_activity_monitoring': 8};
  best joint rank of this measurement None -> 328
```

How to read it: a tier-2.5 record (AUROC lower bound 0.62, evidence 0.6 x 0.24 = 0.145) makes autonomic testing for
ME/CFS rankable under `evidence_weighted`, but it is weak evidence next to the public tier-1 MUSCLE-ME records, so it
enters the joint region x measurement list at rank 328 (county) and does not displace CPET or wearable pairs from the
joint top 25. It does not replace the tier-1 steps record for the wearable bundle. Equal-weight ranks do not move (by
design they carry no performance information).

## 9. See it: dashboard, heatmap, embedding, ask, trace

* Dashboard page 6 **Your data** (http://localhost:8599/Your_Data): the dataset row, SYNTHETIC / demo / comparator
  badges, the healthy-control warning, and tabs for the locked plan, the metrics, the ranking change and
  `heatmap_measurement_evidence.png` (condition x measurement evidence, the user cells outlined) and
  `embedding_dpv_pca.png` (the dataset's own Digital Phenotype Vector PCA, labels colouring points only). The CSVs
  behind both are in `results/byod/synthetic_demo/`.

```bash
uv run measure-it byod list
uv run measure-it ask "Where should we deploy autonomic testing for ME/CFS?" --top-n 3
uv run python -m measure_it.tools trace_evidence '{"object_id": "measurement_performance:me_cfs|autonomic_function_testing"}'
uv run measure-it byod show synthetic_demo          # get_user_dataset_metric: plan, models, record, deploy report (JSON)
uv run measure-it validate
```

```text
byod_synthetic_demo: SYNTHETIC tutorial dataset: ME/CFS vs healthy, ECG patch + actigraph + EHR + proteomics [deployed] owner measure-it examples (SYNTHETIC); me_cfs -> me_cfs; class hrv; comparator healthy; demo True; plan 80e73f2648e5
...
- Okmulgee County, Oklahoma: Measurement performance (tier 2.5: computed by this engine under a locked plan on
  user-supplied, local-only data (not public; not re-analysable by others)): BYOD-SYNTHETIC_DEMO-PRIMARY AUROC 0.701
  (0.620-0.771), clinical_case_definition: SYNTHETIC label: ... vs SYNTHETIC healthy controls; sensitivity 0.30
  (0.18-0.40) at specificity 0.90. Performance measured against healthy (or recovered) controls overstates real-world
  performance: ... [geo:40111]
...
{"performance_status": "known", "quality_tier": 2.5, "quality_tier_label": "user_supplied_own_computation",
 "selected_record_id": "BYOD-SYNTHETIC_DEMO-PRIMARY", "auroc": 0.701..., "measurement_evidence": 0.1445..., ...
 "lineage": [..., {"object_id": "signature:byod_synthetic_demo|me_cfs|model_auroc:primary",
                   "relation": "result behind the scored performance record"}]}
{
 "tool": "get_user_dataset_metric",
 "status": "ok",
 "query": {"dataset_id": "byod_synthetic_demo"},
 "data": {
  "dataset": {"dataset_id": "byod_synthetic_demo", "title": "SYNTHETIC tutorial dataset: ...", "synthetic": true,
              "demo": true, "comparator_type": "healthy", "primary_condition_id": "me_cfs", "measurement_class": "hrv",
              "n_cases_primary": 90, "n_controls_primary": 110, "status": "deployed",
              "plan_sha256": "80e73f2648e583c9...", "deployed_with_demo": true, ...},
  "locked_plan": {...}, "models": [...],
  "performance_record": {"record_id": "BYOD-SYNTHETIC_DEMO-PRIMARY", "quality_tier": 2.5, ...,
    "caveats": ["SYNTHETIC tutorial data (generated, not real people): demonstrates the mechanics only.",
                "User-supplied, local-only data: ...", "Comparator: healthy controls. Performance against healthy
                controls is optimistic: ...", "This dataset id was evaluated 1 time(s) before
                (results/byod/LEDGER.jsonl).", ...]}, ...
measure-it validate (2026-09-28T19:15:35+00:00): PASS
  [PASS] provenance: 170/170 processed tables pass check_provenance; ...
  [PASS] person_layer_no_geography: 73/73 person-layer tables carry no geographic identifier column
  ...
```

`trace_evidence` on `signature:byod_synthetic_demo|me_cfs|model_auroc:primary` resolves the model row (AUROC 0.701,
permutation p 0.005, 90 cases / 110 controls, `evidence_level supported_single_dataset`) with an empty SOURCE_REGISTRY
source list: user data are not a registered public source.

## 10. Remove

```bash
uv run measure-it byod remove synthetic_demo
uv run measure-it byod list
uv run measure-it validate
```

```text
byod remove byod_synthetic_demo: deleted 18 item(s)
  rebuilt canonical unions (...) and the Digital Phenotype Vectors
[metric_link] 21 records, 102 condition x measurement rows (16 known, 3 partial, 83 UNKNOWN) (6s)
  metric link rebuilt (102 rows); scored combinations whose performance changed: [('me_cfs', 'autonomic_function_testing')]
  re-scored me_cfs x autonomic_function_testing x county: 3144 rows
  re-scored me_cfs x autonomic_function_testing x state: 51 rows
removed byod_synthetic_demo; re-scored 2 combination(s)
no user datasets (measure-it byod ingest <dir>)
measure-it validate (2026-09-28T19:16:48+00:00): PASS
  [PASS] provenance: 160/160 processed tables pass check_provenance; ...
  [PASS] person_layer_no_geography: 65/65 person-layer tables carry no geographic identifier column
```

After removal, `measurement_performance`, `measurement_performance_records`, `deployment_opportunities`,
`deployment_candidates`, `participants`, `phenotype_signatures`, `participant_phenotype_embeddings` and
`digital_phenotype_datasets` hold the same values and dtypes as before step 6 (checked column by column; only the
build-time provenance of the 3,195 re-scored rows differs). `results/byod/LEDGER.jsonl` keeps the history
(`--purge-ledger` deletes it); the input folder `examples/byod_synthetic/synthetic_demo/` is yours and is left alone.

## Your own data next

Copy `templates/byod/`, replace the example files and every value in `manifest.yaml` (pre-specify the primary analysis
before you look at any result), and run steps 5-8 on your folder. For a real dataset leave `demo: false`; its record
then enters the rankings at `deploy` and in every later pipeline run until `byod remove`. Read
[BRING_YOUR_OWN_DATA.md](BRING_YOUR_OWN_DATA.md) section 1 first: de-identification, permission and any sharing of
results remain the data owner's responsibility.
