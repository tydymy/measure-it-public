# Bring your own data (BYOD)

An outside group can add its own person-level data (device or tool measurements, EHR labs and diagnoses, omics layers)
to this engine, evaluate a measurement under a plan locked before any outcome is computed, and see how the resulting
performance record moves the deployment ranking, without editing core code. Everything happens on the group's own
machine; nothing is uploaded, published or redistributed.

```bash
uv run measure-it byod validate <folder>     # schema + privacy + geography checks; refuses with every reason
uv run measure-it byod ingest   <folder>     # person-layer partitions, provenance, registry entry, DATA_AUDIT.md
uv run measure-it byod evaluate <folder>     # locks the primary analysis (sha256), then runs it
uv run measure-it byod subgroup <folder>     # locks the subgroup plan; device-defined subgroup -> computable phenotype
uv run measure-it byod deploy   <folder>     # metric link + evidence_weighted scoring with the new record
uv run measure-it byod list                  # user datasets on this machine
uv run measure-it byod show     <dataset_id> # locked plan, models, performance record, deploy report (JSON)
uv run measure-it byod remove   <dataset_id> # removes the dataset and every record; restores the rankings
```

A worked, fully synthetic example that runs end to end: [WALKTHROUGH.md](WALKTHROUGH.md). Template folder:
`templates/byod/` (manifest, JSON schema, one example file of each kind). Code: `src/measure_it/byod/`.

## 1. Responsibilities of the data owner

The engine checks for obvious identifying and geographic content, but it cannot certify de-identification. Before a
folder is given to `byod`, the data owner is responsible for:

* **De-identification** under the rules that apply to the data (for US health data, e.g. HIPAA Safe Harbor or an
  expert determination). The manifest must state `data_use.deidentified: true`.
* **Permission**: the IRB protocol, data-use agreement or consent that allows this secondary analysis on this machine
  (`data_use.use_permitted: true`, a `statement`, and ideally an `irb_or_dua_reference`). Nothing here replaces that
  review.
* **Labels**: the case definition, who applied it, and the comparator. The engine records them as supplied.
* **Where results go.** Results stay in local-only folders (below). Sharing any of them (including aggregate numbers)
  is the data owner's decision under the same agreement; small cells, rare categories and per-feature statistics can
  still be disclosive.

## 2. The input folder (contract)

```text
my_dataset/
  manifest.yaml              required: who, what, permission, labels, comparator, measurement, primary analysis
  participants.csv           required: participant_id, one column per label, optional age_band, sex
  device_<name>.csv|parquet  any number: per-participant features, or per-day rows with day_index
  ehr_labs.csv               optional: participant_id, loinc, value, unit[, lab_name, day_index]
  ehr_diagnoses.csv          optional: participant_id, icd10cm[, day_index]
  ehr_medications.csv        optional: participant_id, rxnorm[, day_index]
  survey_<instrument>.csv|parquet  any number: participant_id[, day_index] + numeric item / score columns
  self_report_history.csv    optional: participant_id, condition[, icd10cm][, day_index]
  omics_<layer>.csv|parquet  any number: one row per participant, one numeric column per feature
  adapters/*.py              optional: your own MeasurementAdapter for a device file
```

Any other data file (`.xlsx`, `.json`, `.txt`, ...) is refused, so nothing is silently skipped. `<name>` and `<layer>`
are lower-case letters, digits and `_`.

**participants.csv.** `participant_id` (study codes: letters, digits, `.`, `_`, `-`; at most 64 characters), one
column per label declared in the manifest holding `1` (case), `0` (control) or blank (not in that analysis),
optionally `age_band` (`40-49`, bands of at least 5 years; everything above 89 as `90+`) and `sex` (`female`, `male`,
`other`, `unknown`). No other column: features belong in the data files.

**device files.** Either one row per participant (`participant_id` + numeric features) or one row per participant and
day (`participant_id`, `day_index`, numeric features). `day_index` is a day number relative to the participant's own
start, never a calendar date. Without an adapter, the engine's `TabularFeatureAdapter` uses per-participant rows as
they are and summarises per-day rows per participant (mean, SD, number of days). Nothing is fitted.

**ehr_labs.csv.** One row per result: `loinc` (e.g. `1988-5`), numeric `value`, `unit`. The engine keeps one unit per
LOINC code (rows in a minority unit are dropped, never converted, and counted in DATA_AUDIT.md) and uses the median per
participant and code.

**ehr_diagnoses.csv.** One row per code: `icd10cm` (e.g. `G93.32`). Models use 3-character categories present in at
least 5 participants. A code that the ontology normaliser maps to the condition being analysed (e.g. G93.32 for ME/CFS)
is **label-defining**: it is excluded from every model and flagged (and left out of the PCA fit) in the Digital
Phenotype Vector.

**ehr_medications.csv.** One row per medication: `rxnorm`, the RxNorm **ingredient** CUI (digits, e.g. `6918`
metoprolol). The codes are not sent to any service to check them (user data stay on the machine), so supply
ingredient-level CUIs (IN / MIN), not clinical drugs. Models use CUIs present in at least 5 participants.

**survey_<instrument>.csv.** Patient-reported outcome instruments: one row per participant (or per participant and
`day_index` when repeated), one numeric column per item or score. The instrument code comes from the file name
(`survey_compass31.csv` -> `COMPASS31`, `survey_fss9.csv` -> `FSS9`; known names in
`configs/harmonize_survey_crosswalk.yaml`, otherwise the upper-cased name) and every column becomes the harmonized
concept `SURVEY:<INSTRUMENT>:<column>` (e.g. `SURVEY:COMPASS31:vasomotor`). Text answers are refused: code them.

**self_report_history.csv.** Self-reported medical history: one row per reported condition, `condition` as the
participant or the study form wrote it (at most 80 characters; the same identifying-content checks as every text
column), optionally `icd10cm` when the study already coded it. Names are mapped to ICD-10-CM through the curated,
exact-alias table `configs/harmonize_self_report_icd10.yaml` (`mapping_method = self_report_map`, confidence at most
`medium`); an unlisted wording stays unmapped, is counted in DATA_AUDIT.md and never enters a model. A reported
condition that maps to the label's condition (e.g. "Long COVID" for a Long COVID label) is label-defining and
excluded from models, like an EHR code.

**omics files.** One row per participant; numeric features. They count as participant-linked because the owner gives
the same `participant_id` in every file; ingestion counts the overlap and writes it into DATA_AUDIT.md.

**manifest.yaml** (schema `templates/byod/manifest.schema.json`, template `templates/byod/manifest.yaml`):

| field | meaning |
|---|---|
| `schema_version` | `1` |
| `dataset_id` | lower-case id (3-41 characters); the engine id becomes `byod_<dataset_id>` |
| `title`, `owner`, `description` | shown in the registry entry, the audit and the dashboard |
| `licence` | licence or data-use statement for this local analysis |
| `data_use` | `deidentified: true`, `use_permitted: true`, `statement`, optional `irb_or_dua_reference` |
| `synthetic`, `demo` | `true` only for generated data; `demo: true` keeps the record out of rankings unless `deploy --demo` |
| `labels[]` | `column` (in participants.csv), `condition` (name, alias, ICD-10-CM code or id, mapped by `measure_it.ontology.normalize`; a scored condition set id such as `long_covid_or_me_cfs` is also accepted), `definition`, `label_basis` (`clinical_case_definition`, `clinical_diagnosis`, `ehr_codes`, `self_report`, `proxy`) |
| `comparator` | `type`: `healthy`, `look_alike` (other causes of the same symptoms) or `general_population`; `description` |
| `measurement` | `class` (a class in `configs/measurements.yaml`, e.g. `hrv`), optional `bundle` (from `configs/relevance.yaml`, must contain the class), `device` text |
| `devices.<name>` | optional per device file: `measurement_class`, `description`, `adapter: {module, class, register_as}` |
| `id_hashing` | `none` or `salted_sha256` (below) |
| `min_group_size` | `{cases, controls}`, default 10 each, never below 5 |
| `omics.<layer>` | optional: `platform` (e.g. `nmr_nightingale`, `olink`), `units`, `description`. With `platform: nmr_nightingale` and `units: nightingale_standard` the clinical-chemistry-equivalent Nightingale measures (total-C, HDL-C, clinical LDL-C, TG, glucose, creatinine, albumin, apoB, apoA1, GlycA) become LOINC `measurement` concepts with explicit unit factors (`configs/harmonize_nightingale_loinc.yaml`); everything else stays OMICS |
| `subgroup_analysis` | optional, **pre-specified** (section 5b): `name`, `within_label`, `device_positive` (`column`, or `block` + `feature` + `threshold` + `direction`), `ehr_domains` (`condition`, `measurement`, `drug`, `demographic`, `survey`), `max_features` (10), `min_feature_participants` (5), `exclude_feature_keys`, `n_splits`, `n_repeats`, `n_permutations`, `n_bootstrap`, `C` |
| `primary_analysis` | **pre-specified**: `label`, `feature_blocks` (`device_<name>`, `ehr_labs`, `ehr_diagnoses`, `omics_<layer>`), optional `combined_blocks`, `specificity` (operating point, default 0.90), `covariates` (`age_band`, `sex`), `n_splits` (5), `n_repeats` (10), `n_permutations` (200), `n_bootstrap` (1000), `C` (1.0) |

**Your own adapter.** A device file can go through a partner-written class implementing
`measure_it.measurements.adapter.MeasurementAdapter` (`preprocess`, `embed`, `phenotype_score`,
`describe_measurement`), in a `.py` file inside the folder. Ingestion imports it, registers it with
`measurements.adapters.register_adapter` and checks the contract (`validate_description`; one numeric row per
participant from `embed`). Phenotype names that are `phenotype_axes` ids (e.g. `decreased_hrv`) resolve to HPO terms
and conditions; others resolve to UNKNOWN. The synthetic example ships one (`adapters/hrv_patch_adapter.py`).
Ingestion runs that code on your machine: only use adapters you trust.

## 3. What `validate` refuses

`byod validate` exits 1 and lists every reason; it never prints an offending value (only its length and first
character). It refuses:

* a column whose name locates a place or a person: ZIP / postcode, county, FIPS, ZCTA, census tract, latitude /
  longitude, GPS, address, street, city, state, region, country, site / clinic / hospital;
* a column that identifies a person: names, date of birth, MRN / medical record, phone, fax, e-mail, SSN, insurance or
  beneficiary numbers, licence / passport / account numbers, URLs, IP addresses, photos, initials;
* a column carrying a calendar date (`*_date`, admission, discharge, death, timestamp) or an exact age (`age`,
  `age_years`): use `day_index` / `year` and `age_band`;
* text values that look like e-mail addresses, phone numbers, SSNs, exact dates, ZIP codes, street addresses or
  coordinate pairs, in any text column of any file (numbers are not scanned, except 5-digit values with a leading zero);
* participant ids that look like MRNs (6 or more digits, or an MRN / patient / chart prefix), SSNs, e-mails, phone
  numbers, dates or personal names;
* ages above 89 not top-coded as `90+`; bands narrower than 5 years;
* a manifest that does not confirm de-identification and permission, fails the schema, names an unknown measurement
  class or bundle, or a label condition the ontology cannot map to exactly one condition;
* fewer than 10 cases or 10 controls (configurable, never below 5) with data in the primary blocks, or fewer than the
  number of CV folds;
* undeclared columns in participants.csv, non-numeric feature columns, free text longer than 80 characters in the lab
  file, malformed LOINC, ICD-10-CM or RxNorm codes, self-reported condition names over 80 characters, a
  `subgroup_analysis` that names an unknown label, device block or feature, participant ids missing from participants.csv, unrecognised data files,
  `synthetic: true` without `demo: true`, and `id_hashing: salted_sha256` without a salt.

These checks are a safety net: a clean report does not mean the data are de-identified.

## 4. `ingest`: what is created

| created | content |
|---|---|
| `data/processed/participants__byod_<id>` | participant_id `byod_<id>:<native id>`, labels, age band (`age_range`), sex |
| `participant_conditions__byod_<id>` | one row per label value (study_label), per ICD-10-CM code (ehr_diagnosis) and per self-reported condition (self_report_history, with its mapping status), with the normalised condition where the code maps exactly |
| `participant_labs__byod_<id>` | LOINC-coded labs (`lab_variable = loinc:<code>`) |
| `participant_wearable_features__byod_<id>` or `participant_device_features__byod_<id>` | adapter embeddings (wearable classes vs clinic devices), columns `<device>__<feature>` |
| `participant_wearable_daily__byod_<id>` | per-day rows of wearable devices |
| `participant_adapter_scores__byod_<id>` | adapter phenotype scores (`measurements.adapters.scores_to_long`): descriptive, not diagnoses |
| `participant_omics_linked__byod_<id>` | omics layers (long) |
| `participant_medications__byod_<id>` | EHR medications (RxNorm CUIs) |
| `participant_survey__byod_<id>` | PRO instruments (long: instrument, item, value, `SURVEY` concept code) |
| `person_concepts__byod_<id>` | the harmonized vocabulary shared with NHANES (DEMOG, ICD10CM, LOINC, RXNORM, SURVEY, DEVICE, OMICS; docs/HARMONIZATION_CONTRACT.md); study labels are outcomes and are not in it |
| `byod_datasets` | one metadata row per user dataset (no person rows) |
| `data/raw/byod_<id>/` | `byod_registry_entry.yaml` (all `SOURCE_REGISTRY` fields, `user_supplied: true`, `local_only: true`), generated `DATA_AUDIT.md` (measured sizes and missingness), `input_manifest.json` (file names, sizes, sha256; the data are not copied) |

Every person-layer row carries the provenance columns with `data_layer = person`, `evidence_type` `person_*`,
`source_geographic_resolution = none` and `provenance_notes` beginning `user-supplied; not redistributed`.
Ingestion then rebuilds the canonical unions (`participants`, `participant_labs`, ..., `phenotype_signatures`) and the
Digital Phenotype Vectors; the user dataset gets its own representation (within-dataset z-scores and per-block PCA;
omics as a precomputed PCA block), never pooled with another dataset. Re-ingesting replaces the dataset's partitions.

**Namespacing and hashing.** Participant ids become `byod_<id>:<your id>`. With `id_hashing: salted_sha256` they become
`byod_<id>:h<first 16 hex digits of sha256(salt ":" your id)>`; the salt is read from the environment variable
`MEASURE_IT_BYOD_SALT` (at least 16 characters) and is never written to disk. Keep the salt with the data owner: the
same salt gives the same tokens (so a re-ingest lines up), a different salt gives unrelated tokens, and without it
nobody can map a token back to your id. There is no cross-dataset participant key.

**Why the registry entry is not in SOURCE_REGISTRY.yaml.** `SOURCE_REGISTRY.yaml` is committed to git, and the Data.gov
discovery step sends each registered source's name to an external catalog. The BYOD entry is therefore written as
`byod_registry_entry.yaml` (not `registry_entry.yaml`), which `measure-it merge-registry` never merges.

## 5. `evaluate`: the locked analysis

1. **Lock.** The plan (label and condition, label basis, comparator, measurement class, feature blocks, combined
   blocks, covariates, specificity, CV / bootstrap / permutation settings, model, decision rule, engine version, seed)
   is hashed (sha256) and written to `results/byod/<id>/PLAN.md` and `plan_lock.json` before any outcome is computed,
   as `measure_it.labs.plan` does for the lab analyses. A later run whose manifest describes a different plan is
   refused; `--amend` accepts it as a recorded amendment, and the performance record then says the plan was amended.
   Every lock, evaluation and removal is appended to `results/byod/LEDGER.jsonl`, which `remove` keeps, so a dataset
   that is removed and evaluated again under another plan shows its history (the record counts earlier evaluations).
2. **Models** (the conventions of `measure_it.labs.lab_cc_engine`): L2 logistic regression; median imputation and
   scaling fitted inside each training fold; participant-level stratified k-fold CV repeated `n_repeats` times.
   `primary` = the pre-specified feature blocks (the only model that becomes a performance record); `combined` = the
   pre-specified combined blocks; each supplied block alone; `covariates` (age band midpoint + female indicator) and
   `primary+covariates` with a paired Delta AUROC.
3. **Metrics.** AUROC = mean over repeats of the out-of-fold AUROC, 95% CI from a participant bootstrap (people
   resampled within class, out-of-fold predictions held fixed); AUPRC of the repeat-averaged out-of-fold probabilities
   (same bootstrap, prevalence shown as its baseline); sensitivity at the manifest's specificity (threshold on the
   controls' out-of-fold probabilities, re-chosen in every bootstrap resample; `metric_link.sens_at_spec_boot`);
   label-permutation null (labels shuffled, CV re-run with one repeat per shuffle, one-sided p = (1 + #null >=
   observed) / (1 + n)). Decision rule: *supported* = CI lower bound > 0.5 and permutation p <= 0.05, else *null*.
4. **Outputs.** `phenotype_signatures__byod_<id>` (model rows and descriptive per-feature Hedges g with BH q),
   `byod_performance_records` (one aggregate record), `results/byod/<id>/EVALUATION.md`, `evaluation_models.csv`,
   `evaluation_features.csv`, `performance_record.json`.

If the manifest's specificity is below 0.90, the record's operating point is recomputed at 0.90 (the metric link's
fixed operating point, docs/ANALYSIS_PLAN_METRIC_LINK.md section 2.4) and both are reported.

## 5b. `subgroup`: a device-defined subgroup and its computable phenotype

The question: inside a condition's cases, the device finds a subgroup (e.g. Long COVID cases with an abnormal
nailfold capillaroscopy pattern). Can the subgroup be recognised **without the device**, from features that public
cohorts, OMOP cohorts or a clinic's EHR hold? `byod subgroup` answers it under a locked plan:

1. **Lock.** The `subgroup_analysis` section (with engine defaults) is hashed and written to
   `results/byod/<id>/SUBGROUP_PLAN.md` and `subgroup_plan_lock.json` before any outcome is computed; a changed plan
   is refused unless `--amend` (recorded; the outputs say so). LEDGER.jsonl records `subgroup_plan_locked` and
   `subgroup_evaluated`.
2. **Population and outcome.** The cases (1) of `within_label` with a defined device finding; device-positive is the
   owner's 1/0 column or `feature direction threshold` on the device block's features.
3. **Features.** `measure_it.harmonize.features.feature_matrix` on `person_concepts__byod_<id>`, restricted to
   `ehr_domains`; device and omics concepts never enter. Removed (label-free): label-defining ICD-10-CM categories
   (codes that map to the label's condition, and the condition's own ICD-10-CM codes), `exclude_feature_keys`,
   binary features with fewer than `min_feature_participants` on either side, numeric features observed in under half
   of the cases or constant (`subgroup_feature_audit.csv` lists every candidate and why).
4. **Model.** L1 logistic regression with at most `max_features` non-zero coefficients (largest C on the grid
   C x 2^-k), median imputation and scaling inside each training fold; repeated stratified k-fold CV, participant
   bootstrap CI, label-permutation null and the `evaluate` decision rule (functions shared with `evaluate`).
5. **Computable phenotype.** Refit on all cases -> `results/byod/<id>/computable_phenotype.json` (schema
   `templates/byod/computable_phenotype.schema.json`, contract section 2): intercept, features with codes,
   transform, center / scale, coefficients, and per feature `ehr_evaluable` (ICD-10-CM, LOINC, RxNorm and
   demographics yes; survey instruments no), the decision threshold (Youden J, in-sample), performance, and strata
   fractions (age band x sex; cells under 10 cases suppressed; pooled `all` row with a Wilson CI). Aggregate tables
   `computable_phenotypes` and `subgroup_strata_fractions` (no participant id, no geography) and
   `results/byod/<id>/SUBGROUP.md`. `measure_it.harmonize.phenotype.score_phenotype(phenotype, feature_matrix)` scores
   the phenotype on any other dataset and reports which features it could not observe (never imputing silently).

Nothing from `subgroup` enters the deployment ranking; `byod remove` deletes its rows and files.

## 6. The new tier: `user_supplied_own_computation` (2.5)

| tier | label | evidence factor | why it sits there |
|---|---|---|---|
| 1 | own_computed_case_definition | 1.0 | this project's own computation on public data with a clinical / study case definition |
| 2 | own_computed_proxy_label | 0.75 | this project's own computation on public data with a proxy or self-reported label |
| **2.5** | **user_supplied_own_computation** | **0.6** | **computed by this engine under a locked plan (CV, bootstrap, permutation null), but on data that are not public: nobody else can re-analyse them, and de-identification and labels are attested by the owner** |
| 3 | published_claim | 0.5 | a published result this project has not re-analysed |
| 4 | none | UNKNOWN | no record |

`measurement_evidence = 0.6 x clip((AUROC lower bound - 0.5) / 0.5, 0, 1)`. Selection per (condition, measurement)
is unchanged (docs/ANALYSIS_PLAN_METRIC_LINK.md section 2.5): the best tier first, so a user record replaces a
published claim or fills an UNKNOWN cell, but never replaces tier 1-2 evidence (it is still listed in the row's
records, and an own-computed null next to a positive raises `evidence_conflict`). A record is `measurement_only` only
when every primary block is a measurement of the declared class (device blocks for a device class, labs for a lab
class, omics for an omics class); a primary model that includes diagnosis codes is kept but not scored, as PDE-040 is.

**Comparator.** The comparator type is carried into the record (`comparator_kind`) and its caveats. Against healthy
controls every number is flagged as optimistic; a look-alike comparator (people with other causes of the same
symptoms) is the clinically relevant contrast. The tier factor does not change with the comparator; the caveat travels
with the record into every recommendation's uncertainties.

With no user dataset evaluated, the metric link, the scores and every committed output are bit-identical to the build
without BYOD (the hook adds nothing).

## 7. `deploy`: what crosses to the map

`byod deploy` rebuilds the metric link with the user record(s), re-scores every scored (condition, measurement)
combination whose performance row changed (county and state; the unchanged engine, same seeds and Monte Carlo), splices
the rows into `deployment_opportunities`, recomputes the joint (region x measurement) ranks, refreshes
`deployment_candidates` if the primary demo combination changed, and writes `results/byod/<id>/DEPLOY_REPORT.md`:

* the relevant cells (the record's condition x every scored bundle containing its class, and the class itself), before
  and after, and why (user record selected, a better tier kept, or demo record excluded);
* per re-scored combination: regions ranked under `evidence_weighted` before -> after, top-10 before and after,
  whether the equal-weight ranks are unchanged (they must be: `equal` carries no performance information), the joint
  top-25 bundle composition and the measurement's best joint rank;
* the top-5 recommendation under `evidence_weighted` for each changed county combination;
* heatmap data (`heatmap_measurement_evidence.csv` / `.png`: condition x measurement evidence, user cells outlined) and
  the dataset's own embedding (`embedding_dpv_pca.csv` / `.png`: its Digital Phenotype Vector PCA, within-dataset).

**What crosses into geography: only the metric.** The only object that reaches the scoring is the aggregate record
(AUROC and CI, sensitivity / specificity at the operating point, n, label basis, comparator) through
`measurement_performance`. Enforced by: validation refusing geographic columns; ingestion re-checking every person
table against `validate.GEO_COLUMN_RE`; `deploy` running `validate.check_person_geography()` and scanning the derived
tables (`deployment_opportunities`, `measurement_performance`, `measurement_performance_records`,
`deployment_candidates`, `byod_performance_records`) for any `byod_<id>:` participant id (it fails if one is found);
and tests that the scoring, facility and geography code never reads a person-layer table
(`tests/test_byod.py::test_scoring_code_reads_no_person_table`). A performance record is never carried to another
condition, condition set or class.

**Demo records.** A record from a manifest with `demo: true` (the synthetic tutorial) enters the metric link only when
`MEASURE_IT_BYOD_DEMO=1`, which `byod deploy --demo` sets for its run and remembers in `results/byod/state.json` (so
`remove` restores the same mode). Any other rebuild of the metric link (the pipeline step `scoring_metric_link`, or
`byod deploy` without `--demo`) leaves it out again.

`deploy` does not rewrite the committed `results/tables` CSVs or `results/SCORING_RESULTS.md`
(`metric_link.build(write_csv=False)`); `--export-static` writes a static snapshot to `results/byod/<id>/static_snapshot/`
and `--build-db` rebuilds the DuckDB. A later full pipeline run includes every non-demo user dataset that is still
ingested.

## 8. `remove`

`byod remove <id>` deletes the dataset's partitions and sidecars, its rows in `byod_datasets`,
`byod_performance_records`, `computable_phenotypes` and `subgroup_strata_fractions`, `data/raw/byod_<id>/` and `results/byod/<id>/`; rebuilds the canonical unions and the
Digital Phenotype Vectors; rebuilds the metric link and re-scores every combination whose performance row changed. The
processed tables then hold the same values as before the dataset was added (checked in
`tests/test_byod.py::test_end_to_end_synthetic`; only the build-time provenance of re-scored rows differs). The ledger
keeps a `removed` line (`--purge-ledger` deletes the dataset's lines).

## 9. Where things live, and git

Local-only folders carry their own `.gitignore` (`*`), so nothing in them can be committed: `data/raw/byod_<id>/`,
`results/byod/`. Processed partitions follow the repository rule (`data/processed/*.parquet` is ignored), but the
`.meta.json` sidecars are tracked by default; the repository `.gitignore` should also ignore
`data/processed/*byod_*.meta.json`. Rebuilt canonical tables (e.g. `participants.meta.json`) change their tracked
sidecars while a user dataset is ingested (row counts and the list of partitions); `remove` restores the content (a
new `created_at`). The synthetic example folder generated by `examples/byod_synthetic/make_example.py` is
`examples/byod_synthetic/synthetic_demo/` (synthetic, safe to commit, but regenerated on demand); `--kind maestro`
writes the MAESTRO-like example (capillaroscopy, ring wearable, PROs, Nightingale-like NMR, Olink-like panel,
self-reported history, no EHR files, a planted subgroup and a `subgroup_analysis` section) to
`examples/byod_synthetic/synthetic_maestro/`.

A **full pipeline run** while a user dataset is ingested includes it wherever the hooks reach (canonical unions, the
Digital Phenotype Vectors, the metric link, the scores), and the pipeline's own steps also write it into committed
summary files (e.g. `results/tables/digital_person_summary.csv`, `results/tables/measurement_performance*.csv`).
Remove user datasets (`byod remove`) before a pipeline run whose outputs will be committed. Recommended root
`.gitignore` lines:

```gitignore
# bring-your-own-data (measure-it byod): user-supplied, local-only material (docs/BRING_YOUR_OWN_DATA.md)
data/raw/byod_*/
data/processed/*byod_*.meta.json
results/byod/
examples/byod_synthetic/synthetic_demo/
examples/byod_synthetic/synthetic_maestro/
```

## 10. Access from tools, the dashboard and the MCP server

* Dashboard page 6 **Your data** (`dashboard/pages/6_Your_Data.py`): user datasets, the locked plan, metrics, the
  performance record and its caveats, the ranking change, heatmap and embedding.
* `measure_it.byod.tools.list_user_datasets()` and `get_user_dataset_metric(dataset_id)`: tool-facade envelopes
  (aggregate only), also registered on the MCP server with the tag `byod`. They are not in `measure_it.tools.TOOLS`,
  so the HTTP API's tool list and its count are unchanged.
* The `measurement_performance` row of the user cell is what `rank_deployment_opportunities`, `ask` and the
  recommendation JSON cite (`measurement_performance:<condition>|<measurement>`, tier 2.5 text in the uncertainties).

## 11. Limits

* One performance record per dataset (the primary analysis). Secondary models are reported, not scored.
* Conditions outside the scored grid (the four demo conditions and two condition sets) get a `measurement_performance`
  row but no ranking; the deployment grid itself is unchanged.
* The validation is pattern-based and conservative; it can refuse legitimate columns (rename them) and cannot detect
  every identifier (free-text is refused for that reason).
* Docker images mount the data bundle read-only: run `byod` with `uv` on a checkout, not inside the containers.
