# Harmonization contract: from a device or omics dataset to people, places and clinicians

Status: draft 2026-10-07. Binding for the four build stages below; change it only by
editing this file and saying so in the stage report.

Goal: the engine stays universal. Any new measurement or omics layer (the first one is MAESTRO/CAPRIO nailfold
capillaroscopy, private, local-only through `measure-it byod`) should be able to answer:

1. **Which people resemble the device-positive patients?** In public or controlled data, never by naming anyone:
   aggregate prevalence in reference cohorts (NHANES), query packs that controlled-access OMOP cohorts (All of Us, N3C)
   run themselves, and a computable phenotype a partner clinic runs inside its own EHR.
2. **Where are many such people, and which clinicians could use the measurement there?**

| stage | what | package (owner) |
|---|---|---|
| B + C | shared person-level vocabulary; device-defined subgroup inside a condition; computable phenotype | `measure_it.harmonize`, `measure_it.byod` (stage BC) |
| D | similar people: NHANES survey-weighted matching, OMOP query packs, clinic export | `measure_it.similar` (stage D) |
| E | modelled county burden (BRFSS small-area estimates) and subgroup burden | `measure_it.geography.sae`, `measure_it.geography.subgroup` (stage E) |
| F | clinicians and sites with measurement-specific activity; generic outreach list | `measure_it.facilities.activity`, `measure_it.outreach` (stage F) |

`pipeline.py`, `configs/scoring.yaml` and `src/measure_it/scoring/` are integrated by the coordinator after the stages
report; stages expose `run()`/`build()` entry points and say which pipeline steps and dependencies they need.

## 0. What the first real dataset holds (MAESTRO, public description, researched 2026-10-07)

MIT Tal lab, "Mucosal And systEmic Signatures Triggered by Responses to infectious Organisms"; single site (MIT CCTR);
cohorts **acute Lyme, chronic symptoms after Lyme, Long COVID, healthy controls** (no ME/CFS or POTS cohort; healthy
comparator only); target 240-300. Layers per participant under one study code: optional nailfold capillaroscopy;
Movano Evie ring (HR, HRV, SpO2, sleep, activity); STAT in-ear cerebral blood flow; NASA Lean Test; RightEye, WAVi EEG,
BrainCheck; Nevisense, Beighton; PROs (PROMIS Cognitive, COMPASS-31, FSS-9, DSQ-PEM, FUNCAP27, BIVSS, VAFS, THRIVE,
self-reported medical history); FLIP antibody profiling; Olink + LEGENDplex; **Nightingale blood NMR metabolomics**,
urine metabolomics; glycomics, steroidomics, epigenetics, 30x WGS; pathogen metagenomics/ddPCR.
**No EHR extract: no ICD-10, LOINC or RxNorm codes.** Consequences for this contract:

* The bridge to public data is: demographics; **self-reported history mapped to ICD-10-CM** (`mapping_method="self_report_map"`,
  confidence at most `medium`); **Nightingale clinical-chemistry-equivalent measures mapped to LOINC** (lipids, glucose,
  creatinine, albumin, GlycA ...; `mapping_method="curated_map"`, with a `platform="nmr_nightingale"` note because NMR and
  routine-lab values are not interchangeable without calibration); **PRO instruments** as `domain="survey"`,
  `vocabulary="SURVEY"`, `concept_code="<instrument>:<item or score>"` (e.g. `COMPASS31:total`), mapped where possible
  to items that exist in reference cohorts (NHANES/BRFSS symptom and limitation items) via a curated crosswalk;
  wearable summaries as `domain="device"`.
* Subgroup analysis `ehr_domains` therefore includes `survey`, and "EHR-type" means "any domain observable outside the
  device": a phenotype built only from survey + NMR features is valid, and its exported code lists say which features
  a clinic EHR can and cannot evaluate.
* Conditions in scope add `lyme_disease` (A69.2) and `ptlds` (no ICD-10-CM code; see ontology mappings) to Long COVID.

## 1. `person_concepts` (stage BC produces; D consumes)

Long person-level table, one row per participant x concept (x day when known). Partitions
`person_concepts__<dataset_id>`; canonical `person_concepts` via `store.union_partitions`.

| column | type | meaning |
|---|---|---|
| participant_id | str | namespaced `<dataset_id>:<native id>` (unchanged from the source tables) |
| dataset_id | str | e.g. `nhanes`, `nhanes0306`, `byod_maestro` |
| domain | str | `demographic`, `condition`, `measurement` (labs/vitals), `drug`, `device`, `omics`, `survey` |
| vocabulary | str | `ICD10CM`, `LOINC`, `RXNORM`, `DEMOG`, `DEVICE`, `OMICS`, `SURVEY` |
| concept_code | str | the code in that vocabulary: `G93.32`, `2160-0`, RxNorm ingredient CUI `6809`, `age_band`, `sex`, ... |
| concept_label | str | human label |
| value_as_number | float/null | numeric value (lab value; 1.0 for presence of a condition or drug) |
| unit | str/null | UCUM-like unit for labs; null otherwise |
| value_as_string | str/null | categorical value (`40-49`, `female`) |
| day_index | int/null | day relative to the participant's own start; never a calendar date |
| mapping_method | str | `native` (code supplied by the source), `curated_map` (a config row), `self_report_map` (self-reported history -> ICD-10-CM), `rxnav_api`, `rule` |
| mapping_confidence | str | `high`, `medium`, `low` |
| source_variable | str | the source column/variable the row came from |
| platform | str/null | assay platform note when the value is platform-specific, e.g. `nmr_nightingale; Total_C mmol/L x 38.67 -> mg/dL` (stage BC addition) |
| + `PROVENANCE_COLUMNS` | | `data_layer="person"`, `evidence_type` one of the existing `person_*` types |

Rules: conditions keep the full ICD-10-CM code; consumers roll up to the 3-character category with
`measure_it.harmonize.features.icd10_category`. Drugs are RxNorm **ingredients** (IN/MIN). Labs keep one unit per LOINC
code per dataset (minority units dropped and counted, never converted silently). Demographics are `DEMOG:age_band`
(bands of 10 years, `90+` top-coded; values like `40-49`) and `DEMOG:sex` (`female`, `male`, `other`, `unknown`).
Person rows never carry geography (guardrail 2, `validate.GEO_COLUMN_RE`).

`measure_it.harmonize.features.feature_matrix(dataset_ids, feature_keys=None, *, condition_level="category")` returns a
wide participant x feature frame. Feature keys are `<VOCAB>:<code>`: `ICD10CM:G93` (presence 0/1),
`LOINC:2160-0` (median value), `RXNORM:6809` (presence), `DEMOG:age_mid`, `DEMOG:female`.

## 2. Computable phenotype (stage BC produces; D and E consume)

JSON file `results/byod/<id>/computable_phenotype.json` (local-only) and one aggregate row per phenotype in the
non-partition table `computable_phenotypes` (no person rows). Schema `templates/byod/computable_phenotype.schema.json`.

```json
{
  "phenotype_id": "byod_<id>:<subgroup name>",
  "version": 1,
  "condition_id": "long_covid",
  "base_population": {"description": "...", "icd10cm": ["U09.9"], "label_column": "lc"},
  "subgroup_definition": "device-positive = capillaroscopy score >= threshold (locked plan)",
  "model": "l1_logistic",
  "intercept": -1.23,
  "features": [
    {"feature_key": "ICD10CM:I73", "domain": "condition", "vocabulary": "ICD10CM", "codes": ["I73.0", "I73.00"],
     "transform": "presence", "coefficient": 0.81, "label": "Raynaud syndrome"},
    {"feature_key": "LOINC:2160-0", "domain": "measurement", "vocabulary": "LOINC", "codes": ["2160-0"],
     "transform": "standardized_value", "center": 0.9, "scale": 0.2, "unit": "mg/dL", "coefficient": 0.12}
  ],
  "decision_threshold": 0.5,
  "performance": {"auroc": 0.71, "auroc_ci_low": 0.62, "auroc_ci_high": 0.80, "n_positive": 40, "n_negative": 60,
                   "comparator": "condition cases without the device finding"},
  "strata_fractions": [{"age_band": "40-49", "sex": "female", "n_cases": 31, "n_positive": 14,
                        "fraction": 0.45, "ci_low": 0.29, "ci_high": 0.62, "suppressed": false}],
  "plan_sha256": "...",
  "caveats": ["..."]
}
```

`strata_fractions` is also written as the non-partition table `subgroup_strata_fractions` (columns as in the JSON plus
`phenotype_id`, `condition_id`). Cells with fewer than 10 cases are `suppressed=true` with fraction null; consumers fall
back to the pooled fraction (row with `age_band="all"`, `sex="all"`, always present).

Stage BC additions (2026-10-07, stage BC report; implemented in `measure_it.harmonize`, `measure_it.byod.subgroup`):

* person_concepts: survey-style sources emit explicit `value_as_number = 0.0` condition rows for an answered "no"
  (NHANES self-report), so "answered no" differs from "not asked". `feature_matrix(dataset_ids, feature_keys=None, *,
  condition_level="category", domains=None, absent="auto", concepts=None)`: presence keys (ICD10CM, RXNORM) are 0 for
  a participant without a row only in a dataset that records the key and has no explicit 0 rows for it; a key a
  dataset never records is NaN for all its participants (`absent="zero"|"nan"` override). Extra key families:
  `SURVEY:<instrument>:<item>`, `DEVICE:<device>:<feature>`, `OMICS:<layer>:<feature>` (medians).
  `feature_catalog(concepts)` gives per key: codes, unit, platform, mapping methods/confidences, `ehr_evaluable`.
* Ages: a source top-coded below 90 gets the top band `<decade>+` (NHANES: `80+`); user-supplied bands of 5 years are
  widened to the containing decade (`40-44` -> `40-49`).
* Computable phenotype features carry, besides the fields above: `code_match` (`category_prefix` for ICD-10-CM keys:
  any code in the category; `exact`; `derived` for DEMOG), `center` and `scale` for every feature (presence features
  too: the case-population prevalence and SD), `coefficient_standardized` (the fitted L1 coefficient on the
  standardized scale; |value| is the importance weight used for coverage), `ehr_evaluable` (bool) + `ehr_note`,
  `mapping_methods`, `mapping_confidences`, `platform`, and for SURVEY features `reference_crosswalk`
  (configs/harmonize_survey_crosswalk.yaml). Linear predictor = `intercept + sum(coefficient * t(x))` with
  `t(x) = x` for `presence` (0/1; the centering is folded into the intercept) and `(x - center) / scale` for
  `standardized_value`. Top-level extras: `dataset_id`, `decision_threshold_rule` (Youden J on out-of-fold
  probabilities), `ehr_evaluable_summary`, `synthetic`, `demo`; `performance` adds `perm_p`, `decision`, `auprc`,
  threshold sensitivity/specificity. `phenotype_id` is `byod_<id>:<subgroup name>`.
* Suppressed strata rows have `n_cases`, `n_positive`, `fraction` and CIs null (small cells are disclosive).
  `subgroup_strata_fractions` also carries `dataset_id` and `stratum_id` (`<phenotype_id>|<age_band>|<sex>`).
  `computable_phenotypes` columns: `phenotype_id, dataset_id, condition_id, subgroup_name, version, model,
  subgroup_definition, n_features, n_features_ehr_evaluable, ehr_evaluable_coefficient_share, feature_keys` (JSON),
  `auroc, auroc_ci_low, auroc_ci_high, perm_p, decision, n_positive, n_negative, positive_fraction,
  decision_threshold, plan_sha256, plan_amended, synthetic, demo, phenotype_json` (the full JSON), `caveats` (JSON).
  Both tables: `data_layer="measurement"`, `evidence_type="derived_evidence_summary"`.
* `measure_it.harmonize.phenotype.score_phenotype(phenotype, feature_matrix, *, missing="center",
  output="probability") -> pd.Series`: a missing value is set to the feature's `center` (neutral relative to the
  training cases) and reported in `attrs` (`coverage` per feature, `missing_features`, `coefficient_mass_available`,
  `n_rows_with_any_missing`); `missing="nan"` returns NaN for any row missing a feature. `phenotype_coverage()` gives
  the per-row share of |coefficient_standardized| mass observed; `validate_phenotype()` checks the JSON schema.

## 3. Similar-people outputs (stage D)

Aggregate only: `similar_reference_prevalence` (phenotype_id, reference dataset, stratum, survey-weighted prevalence
and 95% CI, unweighted n, feature coverage = share of |coefficient| mass available in the reference),
`similar_reference_profile` (feature means device-positive-like vs others). Query packs and clinic exports are files
under `results/byod/<id>/similar/` (local-only).

Stage D additions (2026-10-07, stage D report):

* `similar_reference_prevalence` columns: `phenotype_id, reference_dataset, variant, age_band, sex, prevalence, se,
  ci_low, ci_high` (Taylor, Wald-t on the design df), `kg_ci_low, kg_ci_high` (Korn-Graubard), `design_df,
  n_unweighted, n_positive_unweighted, weighted_population, nchs_reliability, feature_coverage, n_features_used,
  n_features, coverage_by_domain` (JSON: domain -> evaluable share of that domain's |coefficient| mass),
  `weighted_share_complete_features, weight_variable, reference_population, condition_id, notes, object_id` +
  provenance (`data_layer="derived"`, `evidence_type="survey_estimate"`, resolution `national`). Strata: `all/all`,
  age band x sex, sex, age band (NHANES bands; adults 18-19 are `18-19`; the top-code decade is `80+`).
  `variant` is `all_available_features`, plus `excluding_nmr_mapped` when a LOINC feature has
  `platform` containing `nmr`/`nightingale` (NMR and routine-lab values are not calibrated).
* `similar_reference_profile`: one row per phenotype feature x reference: `available_in_reference,
  nmr_platform_mismatch, mean_/se_/n_phenotype_positive, mean_/se_/n_others` (survey-weighted).
* New table `similar_reference_distance`: Gower distance of reference adults to the device-positive centroid over the
  shared features, as weighted quantiles (`q05 ... q95`) and weighted mean, for `all_adults`, `phenotype_positive`,
  `others` (cells with n < 10 are not reported).
* Optional phenotype JSON key (stage BC may add it): `"device_positive_centroid": {"<feature_key>": mean, ...}`
  (aggregate device-positive feature means, raw units; presence features as proportions). Without it the distance is
  not computed and the table says so. Also read: `device_positive_profile`, or `feature_profile` rows with
  `feature_key` + `mean_positive`.
* Reduced phenotype rule (all stage-D targets, = `harmonize.phenotype.score_phenotype(missing="center")`): a feature
  the target cannot evaluate, or an item-missing value, is held at its `center` (standardized: contribution 0;
  presence: coefficient x training prevalence). The OMOP / clinic packs fold the not-evaluable features into a reduced
  intercept; inside an EHR an evaluable code or drug that is not recorded is 0. Coverage weights are
  |`coefficient_standardized`| (|`coefficient`| when absent), as in `harmonize.phenotype`. ICD-10-CM features with
  `code_match="category_prefix"` match every code of the 3-character category. SURVEY features are evaluated in a
  reference only through a `reference_crosswalk` row with `match_quality="same_item"`; same/related-construct rows
  are reported, never substituted. The SQL / clinic packs also understand the transforms `log_standardized_value`,
  `value`, `above_threshold`, `below_threshold` (with `threshold`), not used by the current schema.

## 4. Geography outputs (stage E)

`geo_condition_burden__brfss_long_covid_sae` (model-based county estimates, `evidence_type="modeled_small_area_estimate"`,
`source_geographic_resolution="county"`, method and validation in DATA_AUDIT.md) and `geo_subgroup_burden`
(phenotype_id, geo_id, geo_level, estimated adults in the subgroup and rate, 95% interval, inputs). A state value is
never relabelled as county prevalence (guardrail 5): the SAE is a model with its own validation.

## 5. Clinician activity outputs (stage F)

`provider_measurement_activity` (npi, measurement_class, hcpcs codes, services, beneficiaries, year, Medicare FFS only,
suppression noted), `provider_condition_rx_signals` (npi, signal id, drug, claims), `geo_measurement_capacity` (geo_id,
measurement_id, active clinicians, services, per 100k), all keyed by the measurement ids of `configs/measurements.yaml`
so any new measurement is a config entry, not code. Outreach lists (named clinicians) are written only under
`outreach/` (untracked) and never committed.
