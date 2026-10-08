# MUSCLE-ME: does a wearable step count separate Long COVID / ME/CFS from healthy controls?

**Dataset:** `charlton_lc_mecfs_cpet_source`: Charlton et al., Nat Commun 2026 (doi:10.1038/s41467-026-75725-y), Source Data sheet `Main`, MUSCLE-ME rows. **Plan (locked before any comparison):** `docs/ANALYSIS_PLAN_CHARLTON_LC_MECFS_CPET_SOURCE.md` (sha256 `dd66de0263eec659...`). **Modules:** `measure_it.ingestion.muscle_me_charlton`, `measure_it.wearables.muscle_me_steps_analysis`. **Computed:** 2026-09-25T21:24:55+00:00. **Seed:** 20260923.

## Answer

* **Primary (locked):** on the 73 complete cases (51 patients, 22 controls), daily steps separate patients from healthy controls with AUROC **0.83 (95% CI 0.73-0.92)**, permutation p 0.0001 (10,000 permutations). The CI lower bound is at or above the locked 0.70 threshold, so the rule returns **"wearable captures the illness-associated functional limitation"**.
* **That verdict does not survive the missing controls.** 8 of 30 controls have no step data. If they all had the lowest observed step count (the locked worst case), AUROC falls to **0.61 (0.46-0.76)**, and the CI includes 0.5 (null); the best case is 0.88 (0.79-0.94). The complete-case result assumes those 8 controls are missing at random. A post-hoc check (not in the plan) shows they are not missing completely at random: 7 of the 8 are women (vs 8 of the 22 controls with steps; Fisher p 0.035). Their peak power and VO2_rel are also lower (median 225 vs 266 W, Mann-Whitney p 0.10; 33.9 vs 39.5 (VO2_rel, units as released), p 0.28), but that is largely what their sex predicts, so the data do not show which way the complete-case AUROC is biased. The worst case is deliberately extreme (all 8 set to the lowest value in the whole sample, below every observed control), so the truth probably lies between the bounds, but the locked rule cannot be called robust.
* **The wearable does no better than the lab test it would replace.** In the same people, steps AUROC 0.82 vs relative VO2peak (CPET) 0.87: difference -0.04 (-0.19 to 0.09), i.e. **no detectable difference either way**. A step counter is cheaper than a maximal CPET and does not require a maximal exercise bout, which can itself provoke PEM. On these data it is about equally informative; it is not more informative.
* **Steps do not tell Long COVID from ME/CFS:** AUROC 0.58 (0.41-0.73). Whatever steps capture is shared by both illnesses, as expected for a measure of activity limitation.
* **Long COVID alone vs controls** is **weak** under the locked rule: AUROC 0.81 (0.67-0.92). ME/CFS vs controls: 0.85 (0.73-0.94), "wearable captures the illness-associated functional limitation".
* **What this shows:** a hip-worn research accelerometer records lower daily activity in people *already diagnosed* with Long COVID (Canadian Consensus Criteria) or ME/CFS than in *healthy* people (enrolled as age- and sex-matched controls, although the 22 controls with steps are mostly men; see (e)). That limitation is part of the case definition (incorporation bias). **What it does not show:** that steps detect undiagnosed illness, or that they separate these illnesses from deconditioning, depression or other chronic disease. There is no such comparator in the data. It is a measurable phenotype, not a diagnostic test.

## Nulls and weak results (reported first)

| analysis | result | verdict |
|---|---|---|
| primary, worst-case bound (8 missing controls = lowest steps) | AUROC 0.61 (0.46-0.76) | null (CI includes 0.5) |
| (c) Long COVID vs ME/CFS by steps | AUROC 0.58 (0.41-0.73), Mann-Whitney p 0.35 | includes 0.5 (steps do not separate LC from ME/CFS) |
| (d) steps minus VO2peak/kg, same people | -0.04 (-0.19 to 0.09) | no difference detected |
| (a) Long COVID vs controls | AUROC 0.81 (0.67-0.92) | weak (CI lower bound below 0.70) |
| (e) sex-only baseline | CV AUROC 0.51 (repeats 0.47-0.55) | null (age not released; all 30 controls were sex-matched, but 7 of the 8 controls missing steps are women, so the complete-case controls are 8 / 14 F/M vs patients 26 / 25) |
| (g) blood lactate (rest, baseline cycling, task failure), exploratory | Hedges g +0.25, -0.17, +0.04; min q 0.11 | no difference |

## Published claim vs our own computation

| | published (paper Table 1, median (IQR) steps/day) | this module (same Source Data) |
|---|---|---|
| healthy control | 7153 (5063-8405); vs ME/CFS p < 0.001 | 7153.6 (5063-8405), n = 22 |
| Long COVID | 4718 (2794-5832); vs control p < 0.001 | 4718.8 (2794-5832), n = 25 |
| ME/CFS | 3704 (2286-5052); vs control p < 0.001 | 3704.8 (2286-5052), n = 26 |

The group medians match Table 1 to within one step/day (the paper appears to truncate the decimals), so the ingested rows are the paper's rows. The paper reports group differences (p-values) in step counts as cohort description. In the text we read, it reports no discrimination statistic (AUROC), no missing-data bound and no claim that steps detect illness. The AUROCs, bounds and head-to-head comparison here are this module's own computation.

## Primary analysis (locked)

| population | n patients / controls | AUROC (95% CI) | AUPRC (95% CI; no-skill) | Mann-Whitney p | decision rule |
|---|---|---|---|---|---|
| complete case | 51 / 22 | 0.831 (0.725-0.918) | 0.926 (0.880-0.965; 0.70) | 8.4e-06 | wearable captures the illness-associated functional limitation |
| worst case (missing controls = lowest observed steps); fill value 733 | 51 / 30 | 0.612 (0.459-0.763) | 0.642 (0.568-0.765; 0.63) | 0.095 | null |
| best case (missing controls = highest observed steps); fill value 16688 | 51 / 30 | 0.876 (0.795-0.942) | 0.926 (0.879-0.966; 0.63) | 1.9e-08 | wearable captures the illness-associated functional limitation |

Label-permutation null (complete case, 10,000 permutations): mean 0.501, 95th percentile 0.623; p = 0.0001. Direction fixed in advance (lower steps = patient); nothing fitted. Stratified person bootstrap, 2,000 resamples.

![primary](../visuals/figures/charlton_lc_mecfs_cpet_source_steps_primary.png)

## Pre-specified secondaries

| | contrast | feature | n cases / controls | AUROC (95% CI) | Mann-Whitney p | locked rule applied |
|---|---|---|---|---|---|---|
| (a) | Long COVID vs healthy controls | Steps | 25 / 22 | 0.813 (0.675-0.918) | 0.00026 | weak |
| (a) | ME/CFS vs healthy controls | Steps | 26 / 22 | 0.848 (0.726-0.942) | 4e-05 | wearable captures the illness-associated functional limitation |
| (b) | patients (LC + ME/CFS) vs healthy controls | VO2_rel | 49 / 30 | 0.848 (0.758-0.926) | 2.5e-07 | lab CPET meets the same 0.70 rule |
| (c) | Long COVID vs ME/CFS (specificity check); AUROC = P(LC steps > ME/CFS steps) | Steps | 25 / 26 | 0.577 (0.414-0.726) | 0.35 | includes 0.5 (steps do not separate LC from ME/CFS) |
| (d) | patients vs healthy controls, people with both Steps and VO2_rel | Steps minus VO2_rel | 49 / 22 | 0.824 vs 0.866; difference -0.043 (-0.192 to 0.092) | bootstrap p 0.53 | no difference detected |

### (e) Cross-validated models and the demographics-only baseline

Stratified group 5-fold CV (one row per person, so person-level), 20 repeats, L2 logistic regression (C = 1). Age is not released, so the demographic baseline is **sex only**. All 30 controls were sex-matched to patients, but 7 of the 8 controls without steps are women, so in the 73 complete cases controls are 8 / 14 (female/male) and patients 26 / 25. Sex is therefore slightly imbalanced in the primary population (patients more often female); the pre-specified sex + steps model shows that adjusting for it does not remove the steps signal.

| model | CV AUROC mean (2.5-97.5% over repeats) | AUROC of repeat-averaged OOF | CV AUPRC mean (no-skill) |
|---|---|---|---|
| (1) sex | 0.509 (0.471-0.545) | 0.344 | 0.700 (0.70) |
| (2) log10(Steps) | 0.819 (0.810-0.828) | 0.824 | 0.921 (0.70) |
| (3) sex + log10(Steps) | 0.839 (0.814-0.853) | 0.840 | 0.933 (0.70) |

| Delta AUROC | role | value | interval |
|---|---|---|---|
| (3) sex + log10(Steps) minus (1) sex | pre-specified Delta (model 3 - model 1) | +0.496 | 0.335-0.660 |
| (2) log10(Steps) minus (1) sex | descriptive (model 2 - model 1) | +0.480 | 0.280-0.667 |
| (3) sex + log10(Steps) minus (1) sex | POST HOC: mean of per-repeat paired AUROC differences (range = 2.5/97.5 percentiles over 20 repeats, not a CI) | +0.330 | 0.301-0.360 |

**Read the pre-specified Delta with care.** It is computed on repeat-averaged out-of-fold predictions. For the sex-only model those predictions have AUROC 0.34, below chance. That is a known cross-validation artefact for an uninformative feature, and it inflates the Delta. The post-hoc per-repeat difference (+0.33) is the honest size. Either way, the increment over sex is almost entirely the steps signal: sex alone carries little here (CV AUROC about 0.5). **Do not compare this Delta with the NHANES Delta AUROC (+0.02-0.03):** there the baseline was a 9-variable demographic model that was itself informative; here the baseline is sex only and near chance.

* Steps-only model, label-permutation null (1000 permutations, one CV repeat each): observed 0.812, null mean 0.478, 95th percentile 0.622, p = 0.0010.
* Sensitivity / specificity of the steps-only model, with a Youden threshold chosen inside each training fold: sensitivity 0.55 (0.51-0.60 over repeats), specificity 0.91 (0.82-1.00). At this cut-off about half the patients are missed and about 1 in 11 healthy controls is flagged; these are against healthy people only.

### (f) Calibration of sex + log10(steps) (n = 73, crude)

Brier 0.148 (no-skill 0.211); calibration intercept +0.00; slope 1.30 (95% CI 0.59-2.02). Too few people to judge the slope. The probabilities also reflect the study's 70% patient share, not any real-world prevalence.

| quartile of predicted risk | n | mean predicted | observed patient rate |
|---|---|---|---|
| 1 | 19 | 0.39 | 0.32 |
| 2 | 18 | 0.64 | 0.67 |
| 3 | 18 | 0.81 | 0.83 |
| 4 | 18 | 0.97 | 1.00 |

## (g) Exploratory: which person-level layers in this file differ (no device claim)

Every numeric column, patients vs healthy controls. Hedges g with a bootstrap CI, AUROC = P(patient > control) with no direction chosen (0.5 = no separation; far from 0.5 either way = separation), Mann-Whitney p, and BH q across 56 variables. This is descriptive. It must not be used to pick a "best" measurement: that would be the many-variants-report-the-best pattern the plan rules out.

| layer | variables | BH q < 0.05 | largest separation, abs(AUROC - 0.5) | variable |
|---|---|---|---|---|
| blood_lactate | 3 | 0 | 0.12 | Rest_lactate (AUROC 0.62) |
| cpet | 20 | 15 | 0.39 | VO2_perc_pred_wasserman (AUROC 0.11) |
| muscle_capillarisation | 5 | 4 | 0.32 | TypeI_CF (AUROC 0.18) |
| muscle_histology | 19 | 9 | 0.26 | Total_I_Positive (AUROC 0.24) |
| muscle_mitochondria | 8 | 6 | 0.41 | S_linked (AUROC 0.09) |
| wearable_accelerometry | 1 | 1 | 0.33 | Steps (AUROC 0.17) |

Sanity check: `VO2_pred`, the predicted VO2peak from demographics and body size, does not differ (g -0.19, q 0.68), which fits the matching. The largest exploratory separations are in muscle mitochondrial respiration and in CPET variables (top AUROC about 0.9 in the reversed direction for each), with fibre-type variables weaker (top about 0.76). These are the maxima of many descriptive comparisons and are therefore optimistic. The muscle variables need a vastus lateralis **biopsy**. They are person-linked tissue physiology, not omics, and not a deployable device. Full table: `results/tables/charlton_lc_mecfs_cpet_source_exploratory_feature_map.csv`.

![forest](../visuals/figures/charlton_lc_mecfs_cpet_source_exploratory_forest.png)

## Deviations from the plan

* **Post hoc (labelled):** the per-repeat Delta AUROC row in (e), added after the pre-specified Delta turned out to be inflated by the below-chance sex baseline. The pre-specified number is still reported, unchanged.
* **Post hoc (labelled), added by the adversarial review:** a missing-data diagnostic comparing the 8 controls without steps with the 22 controls with steps on sex, peak power and VO2_rel (no step values involved; table `results/tables/charlton_lc_mecfs_cpet_source_posthoc_missing_controls.csv`). It changes no locked number; it shows the missing controls are not missing completely at random (7 of 8 are women); it does not show the direction of any bias.
* No other deviations. No feature, threshold, model or subgroup was added. Published Table 1 medians were read after the plan was locked.

## Limitations

* Single-site case-control study (Amsterdam UMC, NCT05225688); healthy controls only, so accuracy is against healthy people, not against other fatigued or deconditioned patients; activity limitation/PEM is part of both case definitions (incorporation bias); patients able to complete maximal CPET and a muscle biopsy (milder end); probably the same people as Appelman et al. 2024; age not released; 8 of 30 controls lack steps.
* The steps device is a research-grade hip ActiGraph (wGT3X-BT; method from Appelman 2024, same trial), not a consumer wrist wearable. The averaging window and wear-time rule are not in the release, and neither is whether the window overlapped the post-CPET PEM period.
* Controls had recovered from SARS-CoV-2 without residual symptoms (per the paper). That is a fair comparator for Long COVID, but it is still a healthy group.
* n = 73 with steps. One site, the Netherlands. Not independent of Appelman 2024.

## Measurement-class and ontology mapping

* `Steps` -> `accelerometry` (configs/measurements.yaml). CPET columns -> `cpet`. Lactate -> `blood_biomarkers`. Muscle biopsy histology and respirometry -> no class exists (reported, not added).
* Session LC -> `long_covid`, ME -> `me_cfs` via `measure_it.ontology.normalize.normalize_condition` (status matched).

## Outputs

* `results/tables/charlton_lc_mecfs_cpet_source_primary.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_secondary.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_head_to_head.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_cv_performance.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_cv_delta_auroc.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_cv_youden_sens_spec.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_cv_permutation_null.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_calibration.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_calibration_bins.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_exploratory_feature_map.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_group_medians.csv`
* `results/tables/charlton_lc_mecfs_cpet_source_posthoc_missing_controls.csv`
* `visuals/figures/charlton_lc_mecfs_cpet_source_steps_primary.png` (pipeline output: `results/figures/drafts/`)
* `visuals/figures/charlton_lc_mecfs_cpet_source_exploratory_forest.png`
* `data/processed/phenotype_signatures__charlton_lc_mecfs_cpet_source.parquet` (63 rows; object ids `signature:charlton_lc_mecfs_cpet_source|<phenotype>|<feature>`)

## Reproduce

```bash
uv run python -m measure_it.ingestion.muscle_me_charlton
uv run python -m measure_it.wearables.muscle_me_steps_analysis
```
Pipeline: `uv run measure-it pipeline --offline --only ingest_muscle_me_charlton --only wearable_muscle_me_steps`.
