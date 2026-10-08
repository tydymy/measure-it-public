# Analysis plan: digital person representation and measurement adapters

Scope: SPEC Phase 2 ("Represent the digital person") and the CAPRIO / capillaroscopy adapter section.
Written 2026-09-23 before any analysis in this plan was run. Deviations are listed at the end,
dated, as they happen. Every number quoted in `docs/ADAPTERS.md` or in results tables comes from
a file written by `measure_it.wearables.digital_person`, `measure_it.wearables.adapter` or
`measure_it.measurements.adapters`.

Code: `src/measure_it/wearables/digital_person.py`, `src/measure_it/wearables/adapter.py`,
`src/measure_it/measurements/capillaroscopy_adapter.py`, `src/measure_it/measurements/adapters.py`.
Interface (shared, unchanged): `src/measure_it/measurements/adapter.py`.

## 0. Inputs (all already processed; audits read)

| input | rows (meta.json) | audit caveats that bind this plan |
|---|---:|---|
| participants__nhanes | 19,931 | no geography; age top-coded at 80; children 6-17 included |
| participant_clinical_features__nhanes | 19,931 | skip patterns (PFQ061, DIQ070), cycle-only items (DLQ H only, HSQ G only), fasting/thyroid subsamples; mortality columns are outcomes |
| participant_conditions__nhanes | 43,810 | self-report diagnoses; rx reason codes are 2013-2014 only; ME/CFS, Long COVID, POTS, EDS essentially absent |
| participant_labs__nhanes | 735,587 | subsample labs; below-LOD fill values |
| participant_medications__nhanes | 27,057 | |
| participant_wearable_features__nhanes / _daily | 12,955 / 130,186 | no HR stream; single MIMS cut-point; sleep = classifier proxy |
| data/interim/nhanes_2011_2014/paxmin_{G,H}.parquet | 78.1 M / 88.2 M minutes | MIMS stored float32, restore to 3 decimals |
| participants__stanford_covid (+ conditions, wearable features, daily) | 2,367 | three sub-datasets with different id schemes and cohorts; only Uwakwe 2025 has a Long COVID label; Alavi non-positives are not verified negatives; dates shifted per person |
| participants__mapmecfs, participant_omics_linked__mapmecfs | 47 / 902,186 | omics-omics linkage only; no wearable/HRV/tilt/CPET data with identifiers in open files |
| phenotype_axes, condition_phenotypes, condition_registry | 13 / 2,709 / 14 | axis-condition annotations are sparse (Monarch) |
| configs/measurements.yaml, configs/relevance.yaml | | curated assumptions (evidence_type curated_config) |

## 1. Canonical person-level unions

Tables: participants, participant_conditions, participant_labs, participant_medications,
participant_wearable_features, participant_wearable_daily, participant_clinical_features,
participant_omics_linked.

* Method: stack every `<table>__<source>` partition; add `dataset_id` where a partition lacks it
  (NHANES -> `nhanes`; mapMECFS -> `mapmecfs`; Stanford keeps its three sub-dataset ids) and
  `source_partition`. Columns absent in a partition are null. Where one column has incompatible
  dtypes across partitions it is cast to string and the cast is recorded in the table description.
* Harmonised keys added to `participants` only: `sex_harmonized` (female / male / UNKNOWN),
  `age_years_numeric` (NaN where only an age band is released), `age_band_harmonized`.
* No cross-dataset participant key is created. `participant_id` keeps its `<dataset_id>:` prefix.
* Checks (written to `results/tables/digital_person_union_checks.csv`): rows = sum of partition
  rows; `participant_id` unique in `participants`; every child-table participant appears in
  `participants` with the same dataset_id (orphans counted, not dropped); prefix of every
  participant_id == dataset_id.

## 2. Digital Phenotype Vector (DPV), one per dataset, never pooled

Datasets: `nhanes`, `stanford_covid_mishra2020`, `stanford_covid_alavi2022`,
`stanford_longcovid_uwakwe2025` (DPV), `mapmecfs` (omics-only representation, section 2.4).

### 2.1 Populations
* nhanes: participants aged >= 18 at screening (adult-targeted questionnaires, labs and the MIMS
  cut-point derived in adults). Children are excluded from the DPV (not from the unions).
* Stanford sub-datasets: all rostered participants (Alavi 2,123 incl. 7 without processed
  wearable data; their wearable block is simply unavailable).
* mapmecfs: all 47 study numbers.

### 2.2 Blocks and features (curated lists in code, all recorded in `digital_phenotype_blocks`)
Blocks: `demographics`, `symptoms_conditions`, `labs`, `medications`, `wearable`,
`functional_status`, plus `clinical_exam` (anthropometry, blood pressure, exam pulse), which the
SPEC list does not name but which is neither a lab nor a symptom.
* Outcome/label columns are NEVER dimensions: Stanford `long_covid_label`, `cohort_group`,
  COVID-positive status; mapMECFS `group_label`; NHANES mortality. They are carried as metadata
  only, so the vector can be used as model input without trivial label leakage.
* Label-defining items are dimensions but flagged `label_defining=True`: Uwakwe post-acute symptom
  reports (the Long COVID label is defined by symptoms >= 12 weeks after infection).
* Stanford wearable block uses only label-free, record-relative features. Features aligned to an
  infection reference date (baseline_*, rhr_acute_*, *_pre14/_post*, recovery, cusum_*_pre/_post,
  infection alarm day) exist only for COVID-positive participants, so even their missingness
  indicator would encode the label; they are excluded.
* NHANES survey design variables, weights and ids are not dimensions.

### 2.3 Transformations (deterministic, no labels used)
1. binary -> 0/1; ordinal -> its numeric code; categorical -> one-hot; clock-time hours -> sin/cos.
2. Non-negative numeric features with sample skewness > 2 in the dataset population -> log1p.
3. Exclude a feature if missing in > 95 % of the population or zero variance (reason recorded).
4. z-score within the dataset population (mean / SD of observed values; unweighted: the DPV
   describes sample participants, it is not a population estimate).
5. Missingness indicator for every retained feature with 0 < missing fraction < 1.
6. Block available = >= 1 observed feature in the block.
7. Per-block PCA (scikit-learn, `random_state=SEED`) fitted on participants with the block
   available, missing z imputed as 0 (the dataset mean); K = min(5, n_features, n_fit - 1).
   Explained-variance ratios and the top loadings go to `digital_phenotype_block_pca`.
   PCA summaries are descriptive coordinates, not phenotypes or subtypes.

### 2.4 mapMECFS omics-only representation
* Blocks: `demographics` (sex, age), `omics_availability` (has_* flags), and one block per open
  omics modality: pbmc_rnaseq, muscle_rnaseq (read counts -> log2(CPM + 1), genes with CPM >= 1 in
  >= 50 % of that modality's samples, top 2,000 by variance), csf_somalogic, serum_somalogic
  (log2 RFU), csf_metabolomics (log2 of the batch-normalised values). Feature-wise z within
  dataset, PCA K = min(5, n - 1). Only PCA coordinates are stored for omics blocks (not 18k z-scores).
* Required label, stored verbatim in `digital_phenotype_datasets.representation_label` and in the
  omics rows: "participant-linked omics exist for 47 participants but NO wearable/physiology data
  is linkable in open data, so this is not a wearable+omics multimodal person" (47 and the
  wearable flag are computed from participants__mapmecfs, not typed).
* NHANES and Stanford datasets carry the statement that molecular context is condition-level
  molecular enrichment reached through the ontology, not individual multi-omics.

### 2.5 Outputs
`participant_phenotype_embeddings` (participant_id, dataset_id, block, dimension, dimension_type in
{zscore, missing_indicator, pca, block_available}, value), `digital_phenotype_blocks`,
`digital_phenotype_block_pca`, `digital_phenotype_datasets`.

## 3. WearableAdapter (wrist accelerometry)

### 3.1 Contract
raw minute frame: participant_id, a minute index (`abs_min` from 00:00 of recording day 1, or
`timestamp`), `mims`, `pred` (1 wake wear, 2 sleep wear, 3 non-wear, 4 unknown), optional `qf`
(QC flag count, default 0) and `dow` (1 = Sunday). Day summaries are rebuilt from minutes with the
PAXDAY definitions (verified on 60 participants before this plan: wake/sleep/non-wear/unknown
minutes = classifier minutes with qf == 0; day MIMS = sum over minutes with qf == 0 and MIMS >= 0).
`embed` calls `nhanes_features.day_level_features` and `nhanes_features.participant_minute_features`
(imported, not re-implemented); same valid-day rule (full 1,440-minute day with >= 600 wake-wear
minutes; >= 4 valid days).

### 3.2 Reproduction validation (pre-specified pass criterion)
Sample with `numpy.random.default_rng(SEED)`: 150 participants per cycle from
participant_wearable_features__nhanes plus 25 per cycle that have PAM day data but fail the rule.
Load their minutes from data/interim, run `embed`, compare every feature column.
Pass = identical NaN pattern and max |adapter - stored| <= 1e-9 * max(1, |stored|) for every
feature, and every failing participant has passes_valid_wear_rule = False. Output:
`results/tables/wearable_adapter_validation.csv`. A test re-runs this on a smaller sample.

### 3.3 NHANES age/sex reference (`wearable_reference__nhanes`)
* Population: participants in participant_wearable_features__nhanes (>= 4 valid days), age >= 6.
* Strata: sex x age band {6-11, 12-17, 18-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+ (top-coded)}.
  A stratum with n < 100 falls back to the sex-pooled age band (flagged).
* Statistics per stratum x feature: unweighted n, mean, SD; survey-weighted (wtmec4yr_pooled)
  p5, p25, median, p75, p95; robust scale = (p75 - p25) / 1.349. Caveat carried in the table: the
  valid-wear subset is not re-weighted for non-wear, so weighted values carry non-wear selection.
* Axis composites get their own reference rows (feature = `axis:<name>`).

### 3.4 Phenotype axes (fixed before looking at any label)
| axis | components (sign) | phenotype_axes link |
|---|---|---|
| low_activity | sedentary_fraction (+), mean_daily_mims (-), active_min_per_day (-) | reduced_physical_activity |
| circadian_disruption | interdaily_stability (-), intradaily_variability (+), relative_amplitude (-) | abnormal_circadian_rhythm |
| activity_fragmentation | astp (+) | fatigue (indirect behavioural correlate) |
| sleep_proxy_irregularity | sleep_proxy_midpoint_sd_h (+) | sleep_disturbance (algorithm-estimated sleep proxy) |

component z = sign * (x - stratum median) / stratum robust scale; composite = mean of available
component z (needs >= ceil(k/2) components); axis score = (composite - stratum composite median) /
stratum composite robust scale. Higher = further from the reference in the named direction.
Scores describe physiology relative to a population reference; they are not diagnoses.
Participants failing the valid-wear rule get NaN with reason `insufficient_valid_wear`; no age/sex
-> the pooled all-ages reference, flagged `reference_stratum = pooled_fallback`.

### 3.5 Secondary, descriptive construct check (NHANES adults >= 18 with scores)
Pre-specified axis -> self-reported contrast: low_activity -> any_functional_limitation_pfq;
low_activity -> PHQ-9 item 4 (tired / little energy) >= 2 vs 0; activity_fragmentation -> PHQ-9
item 4 >= 2 vs 0; circadian_disruption -> any_functional_limitation_pfq;
sleep_proxy_irregularity -> told_doctor_trouble_sleeping; sleep_proxy_irregularity -> PHQ-9 item 3
(sleep) >= 2 vs 0. Metric: AUROC of the axis score as a ranker (no fitting, no thresholds),
bootstrap 95 % CI (1,000 resamples), label-permutation null (1,000 permutations, two-sided p).
Comparators: AUROC of age alone and of the unadjusted primary raw feature (to show how much of
a raw-feature association is age). Unweighted, six contrasts, no multiplicity correction (stated).
These are construct checks, not diagnostic performance; null results are reported as such.
Output: `results/tables/wearable_adapter_construct_check.csv`.

## 4. WearableHeartRateAdapter (optional; Stanford daily RHR)
* Input: participant x day frame as in participant_wearable_daily__stanford_covid (or minute HR,
  which is summarised with `stanford_features.daily_summary`). `embed` = `stanford_features.
  longitudinal_summary` (imported) plus two label-free personal-baseline deviation features.
* Personal baseline (label-free, because Uwakwe has no infection date): participant's own median
  daily RHR and 1.4826 x MAD over the whole record; daily robust z; 7-day centred rolling mean
  (>= 4 observed days). Axes: `rhr_elevation_vs_personal_baseline` = 95th percentile of the
  7-day rolling z (a percentile, not the maximum, so the score does not grow mechanically with
  record length); `activity_reduction_vs_personal_baseline` = -(5th percentile of the 7-day
  rolling steps z) (NaN without steps). Requires >= 28 RHR (or step) days. No population reference exists for wearable RHR in this project
  (NHANES has none), so these are within-person scores only.
* Validation: `embed` on the stored daily table reproduces the longitudinal_summary columns of
  participant_wearable_features__stanford_covid (same pass criterion as 3.2).
* Descriptive check: Uwakwe Long COVID label vs each axis, AUROC + bootstrap CI + permutation p
  (exploratory; the audit already reports near-null record-relative HR features).

## 5. Capillaroscopy stub, registry, contract test
* CapillaroscopyAdapter: no data, no numbers; preprocess/embed/phenotype_score raise
  NotImplementedError naming what CAPRIO must supply; describe_measurement status 'stub'.
* Registry `measure_it.measurements.adapters`: get_adapter(name), list_adapters(),
  adapter_for_bundle(bundle_or_alias) from relevance.yaml measurement_bundles[*].adapter; downstream
  resolution functions consume only MeasurementDescription + participant scores.
* Contract test: an in-test `SyntheticFakeAdapter` (synthetic data only inside the test) passes
  through the same registry/resolution/score-standardisation functions unchanged.

## 6. Adapter outputs on public data
`participant_adapter_scores__nhanes` (WearableAdapter), `participant_adapter_scores__stanford_covid`
(HR adapter) and the union `participant_adapter_scores`; `measurement_adapter_registry`
(one row per adapter x phenotype, object_id `measurement:<bundle id>`).

## 7. What would count as failure / null
* Reproduction mismatch in 3.2 or 4 -> the adapter is not a faithful wrapper; reported, not tuned.
* Construct-check AUROCs whose CI includes 0.5 are reported as null.
* Any PCA block whose first component explains < 10 % of variance is reported as weakly structured.

## Deviations (all dated 2026-09-23, made after the first run unless stated)

1. **Result classification made stricter (after seeing the first run).** Section 7 said a CI that
   includes 0.5 is null. The first run's Uwakwe check had a bootstrap CI that just excluded 0.5 while
   the pre-specified permutation test gave p = 0.07. An AUROC is now called an "association" only
   when the bootstrap CI excludes 0.5 AND the permutation p < 0.05; if only the CI excludes 0.5 it
   is "inconclusive". This can only weaken claims (`adapter.classify_result`).
2. **Ages below the reference range are not scored.** NHANES 2013-2014 recorded accelerometry from
   age 3: 410 children aged 3-5 pass the valid-wear rule. Section 3.4 only defined the pooled
   fallback for missing age/sex. Scoring 3-5-year-olds against a pooled 6+ reference would compare
   them with older children and adults, so a known age below 6 now gets no score, with status
   `outside_reference_age_range`. The reference population (age >= 6) is unchanged: 12,545 people.
3. **Age-band boundaries** use [lo, hi + 1), so fractional ages such as 17.9 fall in a band
   (bug fix; NHANES ages are integers, so no NHANES number changes).
4. **Stanford wearable block includes `device_class`** (one-hot; Fitbit minute vs sparse HealthKit),
   which section 2.2 did not list. HR sampling density depends on the device, so the device is kept
   as a measurement-context dimension. It has zero variance in Uwakwe (one class) and is excluded
   there automatically.
5. Not deviations, recorded for clarity: the WearableAdapter claims only the `accelerometry`
   measurement class. The sleep proxy is not polysomnography, so `sleep_objective` is not claimed.
   NHANES medication flags are the 15 most common generics in the DPV population; rx-reason groups
   need >= 10 positives; Uwakwe symptom items need >= 5 positives.
6. Process note: in this build environment the file-writing tool refused `.md` files, so this plan
   and docs/ADAPTERS.md were written with a shell heredoc. The content is unchanged by that.

## Review corrections (independent review, 2026-09-23, after all results above existed)

These were made by a reviewer after the first complete run. None of them changes a pre-specified
axis, contrast, sample or pass criterion. Items 7-10 fix bugs; items 11-13 add sensitivity checks or
edge-case checks.

7. **Pooled-fallback composite reference.** The pooled `axis:*` reference rows were computed from
   composites standardised against each participant's own stratum. A pooled-fallback participant's
   composite, though, is standardised against the pooled component reference. The pooled rows now
   use pooled-standardised composites. Before the fix, 3-component axes of participants without
   age/sex came out on the wrong scale (robust SD 0.70 instead of 1 on a synthetic age-dependent
   population). No NHANES participant uses the pooled row, so NHANES scores are unchanged.
8. **`conditions_for_description` research-activity route.** The code expected
   `canonical_condition_id`, but `condition_measurement_evidence` uses `condition_id`, so the route
   always returned UNKNOWN. It now reads either column. It skips grid cells with
   `measurement_evidence_strength = none` and carries each cell's `object_id`.
9. **`molecular_context_route`.** The code filtered on a `match_type` key that does not exist, so
   partial text matches were routed. For example, 32 NHANES "Headache" (R51) prescription-reason
   rows went to migraine. Now the ICD-10-CM code decides alone when present; otherwise only exact
   matches count. Partial matches are listed as unconfirmed candidates and not used.
10. **DPV blocks outside the SPEC list** (mapMECFS `omics_availability`) now get `available:` rows
    and a PCA like every other block. `participant_phenotype_embeddings` now carries a
    `label_defining` column.
11. **Construct-check design sensitivity.** Added a Rao-Wu PSU bootstrap CI (n_h - 1 PSUs drawn per
    masked-variance stratum) and a survey-weighted AUROC next to the pre-specified unweighted
    participant-bootstrap results. The classification rule still uses the pre-specified columns.
12. **NHANES reproduction edge case.** Every rule-passing participant with fewer than 9 recorded days
    (235) was added. For these records the adapter's minute grid length differs from the pipeline's.
    The random sample and its seed stream are unchanged.
13. **HR adapter raw-minute reproduction.** The pre-specified Stanford check re-ran the same function
    on the same stored daily table (pass-through only). A second check now runs the adapter's minute
    path on the released Uwakwe 2025 minute files (126 participants). The Uwakwe label check also
    reports RHR-day counts by label.
