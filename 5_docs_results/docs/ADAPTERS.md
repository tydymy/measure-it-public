# Measurement adapters: wearables today, CAPRIO capillaroscopy tomorrow

The engine asks: which measurable phenotype, which measurement, where to deploy it, and with whom.
Only the first step depends on the measurement modality. That step sits behind one interface,
`MeasurementAdapter`. Every later layer (measurement registry, ontology, molecular context,
geography, clinics / research sites, scoring, MCP tools) reads only what an adapter declares
(`MeasurementDescription`) and the per-participant scores it emits. Replacing wrist accelerometry
with nailfold capillaroscopy therefore means writing and registering a new adapter. No downstream
module is rewritten.

Every number below comes from a file written by this code; the file is named next to it.
Analysis plan and deviations: `docs/ANALYSIS_PLAN_DIGITAL_PERSON.md`.

## 1. The interface (shared, unchanged: `src/measure_it/measurements/adapter.py`)

```python
class MeasurementAdapter(ABC):
    def preprocess(self, raw_data) -> pd.DataFrame        # raw person-level data -> cleaned frame keyed by participant_id
    def embed(self, raw_data) -> pd.DataFrame             # -> one row per participant_id, numeric embedding columns
    def phenotype_score(self, embedding) -> pd.DataFrame  # -> one row per participant_id, one column per phenotype axis
    def describe_measurement(self) -> MeasurementDescription

@dataclass
class MeasurementDescription:
    adapter_id, modality, measurement_ids,   # measurement_ids = ids in configs/measurements.yaml
    signals, input_unit_of_observation, embedding_features, phenotypes,
    public_training_data, requires_private_data, status ("implemented" | "stub"), limitations
```

Scores describe measured physiology relative to a reference. They are not diagnoses.

## 2. The adapters

| registry name | class (module) | status | measurement_ids | phenotypes (score columns) | reference |
|---|---|---|---|---|---|
| `wearable` | `WearableAdapter` (`measure_it.wearables.adapter`) | implemented | accelerometry | low_activity, circadian_disruption, activity_fragmentation, sleep_proxy_irregularity | NHANES 2011-2014 age/sex reference (`wearable_reference__nhanes`) |
| `wearable_heart_rate` | `WearableHeartRateAdapter` (same module) | implemented | wearable_heart_rate, accelerometry | rhr_elevation_vs_personal_baseline, activity_reduction_vs_personal_baseline | the person's own whole-record baseline (label-free) |
| `capillaroscopy` | `CapillaroscopyAdapter` (`measure_it.measurements.capillaroscopy_adapter`) | **stub**, requires_private_data | capillaroscopy, microvascular_function | planned: capillary_rarefaction, abnormal_capillary_morphology, microhemorrhage_burden, reduced_capillary_flow | to be supplied by CAPRIO |

`configs/relevance.yaml` already routes the bundles: `wearable_autonomic_activity_monitoring.adapters = [wearable,
wearable_heart_rate]` (a bundle may declare several adapters; the first is its primary adapter, returned by
`adapter_for_bundle`; `adapters_for_bundle` returns all of them; the older single-name form `adapter: <name>` is still
read), `nailfold_capillaroscopy.adapter = capillaroscopy`; `autonomic_function_testing` and
`exercise_capacity_testing` have `adapter: null` (no adapter, reported as UNKNOWN with the reason).

### 2.1 WearableAdapter (wrist accelerometry)
* **Input contract** (raw minute frame): `participant_id` (namespaced `<dataset_id>:<native id>`),
  `abs_min` (minutes since 00:00 of recording day 1) or `timestamp`, `mims`, `pred` (1 wake wear,
  2 sleep wear, 3 non-wear, 4 unknown), optional `qf` (QC flags) and `dow` (1 = Sunday).
  `load_nhanes_minutes(seqns, cycle)` produces exactly this frame from the NHANES interim PAXMIN parquet.
* **preprocess**: validates, restores MIMS to its published 3-decimal precision (float32 storage
  would turn a minute at exactly the 10.558 MIMS cut-point inactive), derives day index and flags.
  `day_table` rebuilds the PAXDAY day summaries from minutes. Classifier minute counts use QC-valid
  minutes. Day MIMS is the sum over minutes with qf = 0 and MIMS >= 0. These definitions were
  checked against PAXDAY before the plan was written.
* **embed**: calls `nhanes_features.day_level_features` and `nhanes_features.participant_minute_features`
  (which call `circadian.py`). Nothing is re-implemented, so the 38 features are by construction
  the definitions in `participant_wearable_features__nhanes`. The valid-wear rule is identical.
* **Validation** (pre-specified; `results/tables/wearable_adapter_validation.csv`, `_summary.json`):
  `embed` was run on real NHANES minutes for 300 rule-passing participants (150 per cycle, seed
  20260923). It reproduced all 38 features: identical NaN pattern, maximum absolute difference
  2.9e-11. All 50 sampled rule-failing participants were flagged `passes_valid_wear_rule = False`.
  Added in review: the adapter sizes its minute grid from the last recorded day, the pipeline from
  max(last day, 9). The two differ only for records shorter than 9 days, so every rule-passing
  participant with such a record (235; 12 were already in the random sample, 223 added, rows
  `sample = all_rule_passing_with_fewer_than_9_recorded_days`) was also checked. All 38 features were
  reproduced for them too.
  This is a real raw-to-feature test: the stored day-level features come from the PAXDAY files, and
  the adapter rebuilds its day table from the minute file alone.
* **phenotype_score**: each component feature becomes a robust z against its sex x age-band stratum
  (survey-weighted median, scale (p75 - p25)/1.349). The signed components are averaged into a
  composite, which is standardised against the same stratum's composite distribution. The axis
  definitions were fixed in the plan before any label was examined:

  | axis | components (sign) |
  |---|---|
  | low_activity | sedentary_fraction (+), mean_daily_mims (-), active_min_per_day (-) |
  | circadian_disruption | interdaily_stability (-), intradaily_variability (+), relative_amplitude (-) |
  | activity_fragmentation | astp (+) |
  | sleep_proxy_irregularity | sleep_proxy_midpoint_sd_h (+) |

  Each component z is returned next to the axis score (`<axis>__z__<feature>`).
* **Reference** `wearable_reference__nhanes` (399 rows): 12,545 rule-passing participants aged >= 6,
  in 18 sex x age-band strata plus a pooled row. Every stratum had >= 100 participants, so no
  sex-pooled fallback was needed. Weights: WTMEC2YR/2. Caveat carried in the table: the valid-wear
  subset is not re-weighted for non-wear. Participants without age or sex use the pooled row
  (flagged). Corrected in review: the pooled `axis:*` rows are now built from composites whose
  components were standardised against the pooled component reference. That is how a
  pooled-fallback participant's composite is built. Before the fix they were built from
  stratum-standardised composites, so pooled-fallback scores of 3-component axes came out on the
  wrong scale: on a synthetic age-dependent population the robust SD was 0.70 instead of 1
  (`test_pooled_fallback_reference_is_consistent_with_pooled_scoring`). No NHANES participant uses
  the pooled row, so no NHANES score changed. A known age below 6 is not scored. NHANES 2013-2014 recorded 410 rule-passing
  children aged 3-5, who get status `outside_reference_age_range`.
* Scores on NHANES (`participant_adapter_scores__nhanes`): 12,545 participants scored on three axes
  and 12,228 on sleep_proxy_irregularity (317 lack >= 2 sleep-proxy nights). 410 are outside the
  reference age range.

### 2.2 WearableHeartRateAdapter (Stanford consumer-wearable resting HR)
* Input: participant x day table (as `participant_wearable_daily__stanford_covid`), or minute HR
  (+ steps). Minute HR is summarised with `stanford_features.daily_summary`.
* embed = `stanford_features.longitudinal_summary`, imported, plus label-free personal-baseline
  features. The baseline is the person's own median and 1.4826 x MAD. The daily z is smoothed as a
  7-day centred rolling mean. The score is the 95th percentile of the rolling z (RHR), or minus the
  5th percentile (steps). It needs >= 28 observed days.
* Validation (`results/tables/wearable_hr_adapter_validation.csv`, `_summary.json`). There are two checks,
  and only the second tests anything beyond pass-through:
  1. `check = stored_daily_wrapper`: `embed` on the stored daily table reproduces all 15 longitudinal
     columns of `participant_wearable_features__stanford_covid` exactly, for all 2,360 participants.
     The stored columns were computed by the same `longitudinal_summary` function from the same
     daily rows. This check therefore only shows that the adapter's daily-input path passes the data
     through unchanged. It is not an independent reproduction.
  2. `check = uwakwe_raw_minutes` (added in review): `embed` on the released Uwakwe 2025 minute files
     (`longcovid_sdr_RHR_Data`, 126 participants) goes through the adapter's own minute path (minute
     median HR -> `daily_summary` -> `longitudinal_summary`). It reproduces all 15 columns for all 126
     participants, with maximum absolute difference 0.
  The step-masked Fitbit / HealthKit minute path (Mishra 2020, Alavi 2022, with step data) was not
  re-run from raw files. It remains unvalidated end to end.
* There is no population reference for wearable RHR in this project (NHANES has no HR stream).
  These are within-person scores. Because the score is an upper-tail statistic, typical values are
  about 1.1-1.3, not 0. Medians: Alavi 1.13, Mishra 1.18, Uwakwe 1.33
  (`participant_adapter_scores__stanford_covid`).

### 2.3 CapillaroscopyAdapter (stub)
No data, no numbers. `preprocess`, `embed` and `phenotype_score` raise `NotImplementedError` with
a message that says what CAPRIO must supply (section 5). `describe_measurement()` works today:
status `stub`, `requires_private_data = True`, measurement_ids `[capillaroscopy, microvascular_function]`,
planned signals (capillary density, morphology, microhemorrhages, flow) and the four planned
phenotypes. Because of that, `measurement_context()` can already resolve its implementers
(rheumatology, vascular, primary care, FQHC, cardiology, physiological laboratory), curated
deployment complexity (2-3) and FDA context. The FDA context: `microvascular_function` has 3
product codes with decisions on record; `capillaroscopy` has no product code found
(`measurement_regulatory_status`). Its phenotypes return UNKNOWN until the ontology layer has
axes for them. Through `condition_measurement_evidence` it also resolves the conditions in which
these measurement classes appear in registered trials and grants. That is research activity, not
evidence that the measurement works. `capillaroscopy` appears for long_covid (moderate: 3 trials)
and dysautonomia (low: 1 trial). `microvascular_function` appears for 9 conditions.

## 3. Data flow

```
                raw person-level data                 describe_measurement()
   (minutes of MIMS | daily RHR | nailfold images)           |
                      |                                       v
   preprocess -> embed -> phenotype_score        measurement_context(desc)          [measure_it.measurements.adapters]
                      |                            |-- resolve_measurement_ids: configs/measurements.yaml, curated
                      v                            |     implementers + deployment_complexity (relevance.yaml),
   scores_to_long(desc, scores, dataset_id)        |     bundles, FDA product codes (measurement_regulatory_status)
                      |                            |-- resolve_phenotypes: phenotype -> phenotype_axes axis -> HPO
                      v                            |     -> Monarch-annotated conditions (condition_phenotypes)
   participant_adapter_scores                      |-- conditions_for_description: + condition_measurement_evidence
   (participant_id, dataset_id, adapter_id,        |     (trials / grants), + implementer overlap with
    phenotype, score, components_json, ...)        |     condition_specialties
                      |                            v
                      +------> phenotype signal / condition -> molecular context (condition-level)
                               -> geography (burden, vulnerability, deserts; ecological)
                               -> clinics / research sites (implementer groups, trials, NIH)
                               -> deployment_opportunities -> MCP tools / dashboard
```

CAPRIO: nailfold image -> (preprocess) QC'd capillary fields -> (embed) microvascular embedding ->
(phenotype_score) capillary phenotype scores -> `scores_to_long` -> the same pipeline. Its
`describe_measurement()` routes the geography, clinic and scoring layers to `capillaroscopy` /
`microvascular_function` and the `nailfold_capillaroscopy` bundle, instead of `accelerometry` and
the wearable bundle.

Contract rules enforced in code (`validate_description`, `scores_to_long`):
* measurement_ids must be ids in `configs/measurements.yaml`, and phenotypes must be non-empty and unique;
* scores have one row per `participant_id`, namespaced with the dataset id. Scores are never pooled
  across datasets, and there is no cross-dataset participant key;
* every described phenotype is a column. Components are kept (`<phenotype>__*` columns become
  `components_json`);
* a stub has no participant scores;
* anything that cannot be resolved returns `UNKNOWN / NOT AVAILABLE` with a reason;
* `condition_measurement_evidence` is a full condition x measurement grid. A cell counts as a
  condition link only if trials or grants were actually counted
  (`measurement_evidence_strength != "none"`). Each link carries the cell's `object_id`
  (`measurement_evidence:<condition>|<measurement>`). Fixed in review: before, the route looked for a
  `canonical_condition_id` column, but the table has `condition_id`, so the route always returned
  UNKNOWN.

The contract test (`tests/test_adapters.py`) runs an in-test `SyntheticFakeAdapter` through
`validate_description`, `resolve_measurement_ids`, `resolve_phenotypes`, `measurement_context`,
`scores_to_long` and `summarize_scores` without modification. Its data are synthetic and exist only
in the test. The fake declares `retinal_imaging`, which no real adapter uses. The test also asserts
that the source of these downstream functions names no modality (no "mims", "nhanes", "stanford",
"capillar", adapter class names), and that they return the same structure for the fake, wearable,
heart-rate and stub adapters.

## 4. How a new adapter plugs in

```python
from measure_it.measurements.adapters import register_adapter, adapter_for_bundle, measurement_context, scores_to_long
register_adapter("capillaroscopy", CaprioCapillaroscopyAdapter, overwrite=True)  # or edit _BUILTIN
adapter = adapter_for_bundle("nailfold capillaroscopy")      # relevance.yaml: nailfold_capillaroscopy.adapter
ctx = measurement_context(adapter.describe_measurement())    # implementers, complexity, FDA, conditions
long = scores_to_long(adapter.describe_measurement(), adapter.phenotype_score(adapter.embed(images)), "caprio_study")
```

**From a partner's own folder, without editing code (bring your own data).** A dataset folder can name its own adapter
class in its manifest (`devices.<name>.adapter: {module, class, register_as}`, a `.py` file inside the folder);
`measure-it byod ingest` imports it, registers it with `register_adapter`, checks `validate_description` and the
one-row-per-participant `embed` contract, and writes `participant_adapter_scores__byod_<id>` through `scores_to_long`.
Device files without an adapter go through `measure_it.byod.adapters.TabularFeatureAdapter` (pre-computed features;
per-day rows summarised per participant). Contract and example: docs/BRING_YOUR_OWN_DATA.md, docs/WALKTHROUGH.md
(`examples/byod_synthetic/`, whose ECG-patch adapter emits the `decreased_hrv` phenotype axis).

## 5. Exactly what CAPRIO must provide

1. **Images and metadata** for `preprocess`: per-participant nailfold capillaroscopy images, with
   participant_id namespaced `<dataset_id>:<native id>`, finger, field, magnification, device,
   acquisition date, and the image-quality rules that decide which fields are usable.
2. **The microvascular embedding** for `embed`: versioned model weights and preprocessing that map
   usable fields to one row per participant (e.g. capillary density per mm, morphology class
   fractions, microhemorrhage counts, flow metrics). The dimension names go in `embedding_features`.
3. **Capillary phenotype definitions and a reference** for `phenotype_score`: one score per planned
   phenotype, computed against a healthy reference cohort (with age/sex strata where they matter),
   using a documented scoring rule. Components are kept as `<phenotype>__*` columns.
4. **A non-stub `describe_measurement()`**: status `implemented`, `public_training_data` or
   `requires_private_data = True`, limitations.
5. **Ontology hooks**: each phenotype must be either (a) a new `phenotype_axes` row with an HPO term,
   added by the ontology layer, or (b) an entry in `PHENOTYPE_AXIS_LINKS` pointing to an existing axis
   with a stated link type. Until then the phenotypes resolve to UNKNOWN, and no condition link is
   invented.
6. **Registration**: `register_adapter("capillaroscopy", ...)`. The bundle already names the adapter.
7. **Its own validation**: reproducibility of the embedding and the reference statistics, as done
   here for the wearable adapters.

## 6. What would NOT change

`configs/measurements.yaml` and the `nailfold_capillaroscopy` bundle (both already exist), the
measurement registry and discovery layer, the ontology tables and normalisation, the molecular
evidence layer, the geography layer (burden, vulnerability, deserts), the clinic and research-site
registries and matching, deployment scoring and its weight sensitivity, the MCP tools, the dashboard,
`participant_adapter_scores` (schema), and every function in `measure_it.measurements.adapters`
downstream of the registry. Deployment opportunities for capillaroscopy would change because the
measurement ids change (implementers, complexity, FDA visibility, conditions). That is the point of
the design; the code stays the same.

## 7. Results on public data (positive, null and negative)

**Construct check** (descriptive, pre-specified; `results/tables/wearable_adapter_construct_check.csv`).
NHANES adults. AUROC of each age/sex-referenced axis score used as a ranker; 1,000 bootstrap
resamples and 1,000 label permutations; no fitting and no thresholds.

| axis | self-reported contrast | n (pos) | axis AUROC (95 % CI) | age alone | unadjusted raw feature |
|---|---|---|---|---|---|
| low_activity | any functional limitation | 8,511 (2,352) | 0.603 (0.589-0.619) | 0.697 | mean_daily_mims 0.676 |
| low_activity | PHQ-9 tired >= 2 vs 0 | 5,565 (1,329) | 0.571 (0.551-0.588) | 0.516 | mean_daily_mims 0.565 |
| activity_fragmentation | PHQ-9 tired >= 2 vs 0 | 5,565 (1,329) | 0.574 (0.556-0.591) | 0.516 | astp 0.557 |
| circadian_disruption | any functional limitation | 8,511 (2,352) | 0.590 (0.577-0.604) | 0.697 | relative_amplitude 0.609 |
| sleep_proxy_irregularity | told doctor of trouble sleeping | 8,714 (2,199) | 0.537 (0.523-0.552) | 0.587 | sleep_proxy_midpoint_sd_h 0.522 |
| sleep_proxy_irregularity | PHQ-9 sleep item >= 2 vs 0 | 6,387 (1,242) | 0.589 (0.571-0.605) | 0.509 | sleep_proxy_midpoint_sd_h 0.584 |

All six CIs exclude 0.5 and all permutation p <= 0.001. The associations are small, though.
The primary CIs come from a participant-level bootstrap and ignore NHANES clustering. A design
sensitivity was added in review (columns `axis_ci_*_design_bootstrap` and
`axis_auroc_survey_weighted`). It uses a Rao-Wu bootstrap: in each masked-variance stratum, draw
n_h - 1 PSUs with replacement (61 PSUs). The design CIs are within 0.004 of the participant CIs; the
weakest is sleep_proxy_irregularity vs trouble sleeping, 0.522-0.551. The survey-weighted
(wtmec4yr_pooled) AUROCs are within 0.01 of the unweighted ones: 0.600, 0.572, 0.574, 0.581, 0.538
and 0.591. No classification changes.
* For functional limitation, **age alone ranks better (0.697) than any axis score**. Most of the
  unadjusted MIMS association (0.676) is age: after age/sex referencing, 0.603 remains.
* For fatigue and sleep items, which barely depend on age, the axis scores slightly exceed the raw
  features.
* The sleep-irregularity association with a doctor-told trouble-sleeping item is weak (0.537).

These are construct checks of a population-referenced score, not diagnostic performance. They
are unweighted and have no multiplicity correction across the 6 contrasts. Draft figure:
`results/figures/drafts/wearable_adapter_construct_check.png`.

**Heart-rate adapter vs the Long COVID label** (Uwakwe 2025, n = 126, 31 positive; exploratory;
`results/tables/wearable_hr_adapter_label_check.csv`). rhr_elevation_vs_personal_baseline has
AUROC 0.387 (bootstrap 95 % CI 0.284-0.494; permutation p = 0.069). The result is
**inconclusive**: the CI and the permutation test disagree. Its direction (lower within-person RHR
elevation among Long COVID participants) has no obvious explanation. The score is not aligned to
infection, because the release has no infection date. Observation density also differs by label:
Long COVID participants contribute fewer RHR days (median 355 vs 491; rhr_days alone gives AUROC
0.423). The score correlates only weakly with rhr_days (Spearman 0.08) and with record span (0.17).
So density does not obviously explain the direction, but it is an uncontrolled difference.
activity_reduction_vs_personal_baseline cannot be computed for Uwakwe, since no step data were
released (UNKNOWN).

**Negative / limiting findings to keep in view**
* The stub has no data, so nothing about capillaroscopy performance is known or implied.
* The two wrist-accelerometry axes with the most direct ontology links (low_activity ->
  reduced_physical_activity, circadian_disruption -> abnormal_circadian_rhythm) have **no Monarch
  condition annotations**, so their condition link is UNKNOWN (`measurement_adapter_registry`).
  Only the indirect / partial links (fatigue, sleep disturbance, tachycardia) reach conditions
  (dysautonomia, eds_hsd, lyme_disease, post_infectious_syndrome) through HPO. Neither Long COVID
  nor ME/CFS is reached through the HPO route. The research-activity route
  (`condition_measurement_evidence`, live since the review fix) links accelerometry to 12 of the 14
  registry conditions, all but post_infectious_syndrome and gastroparesis. That records trial and
  grant activity only; it does not validate the measurement.
* NHANES 2011-2014 has no heart-rate stream and no ME/CFS, Long COVID or POTS items. The wearable
  adapter's population reference is therefore not an invisible-illness cohort.

## 8. Digital person representation (the adapter's person-level context)

`measure_it.wearables.digital_person` builds canonical person-level unions and one Digital
Phenotype Vector per dataset. Datasets are never pooled.

* Unions (`results/tables/digital_person_union_checks.csv`): participants 22,345 rows across 5
  datasets, participant_id unique; participant_wearable_features 15,315; participant_wearable_daily
  457,120; participant_labs 735,587; participant_omics_linked 902,186. Row counts equal the partition
  sums. There are 0 prefix mismatches and 0 participants missing from `participants`. No
  cross-dataset key exists.
* DPV (`results/tables/digital_person_summary.csv`, `digital_phenotype_block_pca.csv`):

  | dataset | participants | blocks | z-score dims (excl. indicators, PCA) | PC1 explained variance across blocks |
  |---|---:|---|---:|---|
  | nhanes (adults >= 18) | 11,977 | 7 (demographics, symptoms_conditions, labs, medications, wearable, functional_status, clinical_exam) | 172 | 0.116 (symptoms_conditions) - 0.365 (clinical_exam) |
  | stanford_covid_alavi2022 | 2,123 | wearable | 15 | 0.339 |
  | stanford_covid_mishra2020 | 118 | wearable | 13 | 0.396 |
  | stanford_longcovid_uwakwe2025 | 126 | demographics, symptoms_conditions, wearable, clinical_exam | 166 | 0.243 - 0.645 |
  | mapmecfs (omics-only) | 47 | demographics, omics availability, 5 omics PCA blocks | 8 (+ omics PCA only) | omics blocks 0.180 (serum SomaLogic) - 0.321 (CSF SomaLogic); omics_availability 0.502 |

  No block has PC1 < 10 %. PCA coordinates are descriptive, not phenotypes or subtypes. The
  PCA and z-scores are fitted on each whole dataset population without labels. Anyone who uses
  them as model input must refit them inside training folds.
  Corrected in review: the mapMECFS `omics_availability` block had z dimensions but no
  `available:` rows and no PCA, because the block loop only covered the SPEC block names. As a
  result, `get_digital_phenotype` reported the block as unavailable. It now has both
  (`participant_phenotype_embeddings` 4,152,767 rows; `digital_phenotype_block_pca` 95 rows).
  Of the 47 mapMECFS participants, 43 have >= 2 open omics modalities and so are actually linked
  across modalities; 4 have one modality only.
  Outcome labels are never dimensions. The 41 Uwakwe symptom items reported >= 12 weeks after
  infection define the Long COVID label. They are kept as dimensions, flagged `label_defining`,
  and excluded from the PCA fit. Since the review, the flag is also a column of
  `participant_phenotype_embeddings` itself (5,166 rows = 41 items x 126 participants), not only of
  `digital_phenotype_blocks`. A consumer reading the embeddings alone can therefore drop them. Infection-aligned Stanford features are excluded, because even
  their missingness encodes the label.
* Labels (`digital_phenotype_datasets`):
  * mapMECFS: "participant-linked omics exist for 47 participants but NO wearable/physiology data
    is linkable in open data, so this is not a wearable+omics multimodal person".
  * NHANES and Stanford: molecular context is condition-level molecular enrichment reached through
    the ontology, not individual multi-omics.
* Query functions (UNKNOWN instead of guessing): `get_digital_phenotype(participant_id)`,
  `describe_dataset_representation(dataset_id)`, `molecular_context_route(participant_id)`.
  Corrected in review: `molecular_context_route` filtered on a `match_type` key that
  `normalize_condition` never returns, so partial text matches were routed. 32 NHANES
  prescription-reason rows for "Headache" / "Prevent headache" (ICD-10-CM R51) went to migraine
  molecular evidence that way. Now an ICD-10-CM code decides alone when the source has one (R51 has
  no registry condition, so it is UNKNOWN; G43 is migraine). Otherwise only exact normalisations
  (`status = matched`) are routes. Partial matches are returned as `unconfirmed_candidates` with the
  reason, and they are not used.

## 9. Reproduce

```
uv run python -m measure_it.wearables.digital_person     # unions + DPV tables
uv run python -m measure_it.wearables.adapter            # reference, validation, scores, checks, draft figures
uv run python -m measure_it.measurements.adapters        # measurement_adapter_registry
uv run pytest tests/test_adapters.py tests/test_digital_person.py
```
