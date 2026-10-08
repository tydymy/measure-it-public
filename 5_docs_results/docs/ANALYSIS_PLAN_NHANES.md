# Analysis plan — NHANES 2011-2014 wearable phenotype modelling (SPEC Phase 1; Tests 1, 2, 6)

_Written 2026-09-23, before any model in this plan was run. Owner modules:
`src/measure_it/wearables/{nhanes_cohort,models,unsupervised,signatures}.py`. Results:
`results/WEARABLE_NHANES_RESULTS.md`, `results/tables/nhanes_*.csv`, `results/figures/drafts/nhanes_*`.
Deviations made after this plan was written are listed in §10 and are never edited out._

**Disclosure.** Before this plan was finalised, one timing benchmark was run (one 80/20 split, fatigue
and ME/CFS-like proxy, wearable and clinical+wearable sets, the four model specifications below). It
printed single-split AUROCs. No specification, grid, feature or target was changed after it. The
HGB specification below is the one written before the benchmark, even though HGB looked weaker than
logistic regression on that split.

## 1. Data and population

Inputs (read-only): `participants__nhanes`, `participant_clinical_features__nhanes`,
`participant_conditions__nhanes`, `participant_medications__nhanes`, `participant_mortality__nhanes`,
`participant_wearable_features__nhanes`. All joins use `participant_id = nhanes:<SEQN>` only.

Population: participants passing the valid-wear rule (>= 4 valid days: full 1440-min recorded day
with >= 600 valid wake-wear minutes) **and** aged >= 18 at screening, **excluding** those pregnant at
the MEC exam (activity and fatigue differ systematically in pregnancy). Measured flow: 19,931 -> 14,693
with PAM data -> 12,955 valid wear -> 8,954 aged >= 18 -> 88 pregnant excluded -> **8,866**.
Each target then restricts to participants with that target answered (§2).

No heart rate or HRV exists in NHANES 2011-2014. The 60-second exam pulse (BPXPLS) is a clinical
feature and is never a wearable feature.

## 2. Pre-specified targets (binary)

| id | definition | population | role |
|---|---|---|---|
| `functional_limitation` | `any_functional_limitation_pfq == 1` (yes to PFQ049, PFQ051, PFQ054, PFQ057 or PFQ059) | aged >= 20 (PFQ asked at 20+), composite answered | primary |
| `fatigue` | DPQ040 (tired / little energy, past 2 weeks) >= 2 (more than half the days) | DPQ040 answered | primary |
| `fair_poor_health` | HSD010 (MEC private self-rated general health) in {4 fair, 5 poor} | HSD010 answered | primary |
| `depression_phq9` | PHQ-9 total >= 10 (all nine items answered) | PHQ-9 complete | secondary |
| `mecfs_like_proxy` | fatigue AND functional limitation AND **no** 'yes' to any major explanatory diagnosis: congestive heart failure, coronary heart disease, angina, heart attack, stroke, emphysema, chronic bronchitis, COPD (2013-14 only), cancer ever, liver condition, diabetes, thyroid problem, anemia treatment in past 3 months, sleep disorder (told by a doctor), rheumatoid or psoriatic arthritis. Missing answers do not exclude. Non-cases = everyone else in the population. | aged >= 20, DPQ040 and PFQ composite answered | primary. **This is an ME/CFS-LIKE PROXY, never ME/CFS**: there is no PEM item, no 6-month duration and no clinical assessment. |
| `mortality` | death by 2019-12-31 in the NCHS public-use Linked Mortality File | LMF-eligible | **positive control** |

Signature-only variants (no supervised model): `mecfs_like_proxy_excl_arthritis` (also excludes any
arthritis), `mecfs_like_proxy_excl_depression` (also excludes PHQ-8 without the fatigue item >= 10),
`mecfs_like_proxy_vs_healthy_reference` (proxy cases vs participants with neither fatigue nor
limitation and no exclusion diagnosis). Prescription reason-code phenotypes (2013-2014 only; cases =
>= 1 prescription whose NCHS-coded reason maps to the group; non-cases = all other 2013-2014
participants in the population): fibromyalgia (M79.7, 34 cases), migraine (G43, 69), IBS (K58, 12),
plus insomnia (G47.0/F51.0) and myalgia (M79.1) as within-data comparators. They are too small for
supervised models and are used only for signatures with CIs.

n, cases and prevalence (unweighted; and survey-weighted with `wtmec4yr_pooled` as a descriptive
figure only) go to `nhanes_targets.csv`.

## 3. Feature sets

* **demographics**: age (top-coded 80), female, race/ethnicity (5 indicators, ref. non-Hispanic
  White), income-poverty ratio, education (ordinal 1-5; missing for ages 18-19 by design).
* **clinical**: demographics + BMI, waist, mean systolic/diastolic BP, exam pulse (60 s), 32 labs
  (CBC, standard biochemistry, HbA1c, total/HDL cholesterol, 25-OH vitamin D, B12; nine skewed labs
  log1p), comorbidity count (18 'ever told' diagnoses), medication count (distinct generic
  prescriptions), current and former smoking. Fasting-subsample labs and the 2011-12 thyroid
  subsample are excluded because they are missing by subsample design (>50%).
* **wearable**: 26 accelerometry features from `participant_wearable_features__nhanes` (activity
  volume, day-to-day variability, weekend-weekday difference, wake/sleep-wear minutes, IS, IV, M10,
  L5, RA, onset times, sedentary fraction and minutes, active minutes, ASTP/SATP, bout lengths,
  active-period duration, sleep-proxy timing/regularity/duration). Clock times are expressed as hours
  after a fixed cut at the population trough (12:00; 18:00 for M10 onset). Five skewed features are
  log1p. Wear/QC bookkeeping (day counts, non-wear minutes and similar) is excluded.
* **demographics_wearable**, **clinical_wearable**: unions.

**Predictor exclusions (overlap with a target).** Never in any predictor set: every PFQ item and
the PFQ composite, every PHQ-9 item/total and DPQ100, both self-rated-health items (HSD010, HUQ010),
DLQ disability items, HSQ unhealthy/inactive-day items, self-reported physical activity (PAQ) and
self-reported sleep (SLQ). Per target: `mecfs_like_proxy` also drops `comorbidity_count`, because it
counts the exclusion diagnoses that define the proxy. Written to `nhanes_predictor_exclusions.csv`.

## 4. Models (inside one sklearn Pipeline per fit)

Preprocessing is fitted on the training fold only: median imputation with missingness indicators,
standard scaling, then:
1. `lr`: unpenalised logistic regression (C = inf, lbfgs).
2. `pen_lr` (**primary model**): penalised logistic regression tuned inside the training fold by
   stratified 3-fold CV on log-loss over C in {0.001, 0.01, 0.1, 1, 10} x l1_ratio in {0 (L2),
   0.5 (elastic net), 1 (L1)} (saga). Chosen hyperparameters are logged per fold.
3. `hgb`: histogram gradient boosting (max_iter 300, learning rate 0.05, 15 leaves, min 40 samples
   per leaf, L2 1.0, early stopping on a 15% validation split of the training fold). Not tuned further.

## 5. Splits

* Primary: participant-level repeated stratified 5-fold CV, 5 repeats (seeds from `config.SEED`).
  Each participant has one out-of-fold (OOF) prediction per repeat.
* Secondary (temporal): train on 2011-2012 (G), test on 2013-2014 (H). For mortality, H has ~2 years
  less follow-up, so prevalence and calibration-in-the-large are expected to shift.

## 6. Metrics and inference

For each target x feature set x model, on OOF predictions:
AUROC, AUPRC (reported next to prevalence), Brier score, calibration intercept (calibration-in-the-large:
logistic model of y with offset logit(p)) and calibration slope (logistic of y on logit(p)),
sensitivity and specificity at a threshold chosen inside the training fold (**rule: predicted
probability >= training-fold prevalence**), and reliability bins (deciles of predicted risk, pooled
over repeats).
Point estimate = mean over the 5 repeats. 95% CI = participant-level percentile bootstrap (B = 1000
multinomial resamples of participants, the metric recomputed on each repeat and averaged). The same
resamples are used for every model/feature set of a target, so differences are paired.

* **Test 1 (wearable signal):** AUROC of `wearable` and of `demographics_wearable` vs `demographics`;
  paired bootstrap Delta-AUROC with 95% CI.
* **Test 2 (multimodal increment):** Delta-AUROC = `clinical_wearable` - `clinical`, paired bootstrap 95%
  CI, per target and model (primary: `pen_lr`). Delta-AUPRC and Delta-Brier reported alongside.
* Temporal split: same metrics on H with bootstrap CIs over H participants.

A result counts as "evidence-supported" only if its 95% CI excludes 0 (Delta) or 0.5 (AUROC) **and**
it survives the negative controls in §7. Positive and null results are reported with equal prominence.

## 7. Negative and positive controls

* **Label permutation (Test 6):** per target, 200 permutations of the label over the population;
  primary model (`pen_lr`) on the `wearable` set, one 5-fold CV per permutation; statistic = OOF
  AUROC. Report null mean, null 95th percentile and empirical p = (1 + #{null >= observed}) / 201
  with observed = repeat-1 OOF AUROC (the same single-repeat design).
* **Wearable row shuffle (increment control):** 20 shuffles per target in which the wearable block is
  permuted across participants (marginals and within-block correlations kept, person linkage broken);
  `pen_lr` on clinical + shuffled wearable with the repeat-1 folds; Delta-AUROC vs `clinical` on the
  same folds. Expected ~0. Report the shuffle distribution vs the observed repeat-1 Delta.
* **Positive control:** mortality. Accelerometry is known to predict all-cause mortality in NHANES;
  if wearable-only AUROC is not clearly > 0.5 and the increment over demographics is not positive,
  the pipeline is treated as broken.

## 8. Unsupervised analysis (wearable features only)

* Input: the 26 wearable features (transformed as §3), median-imputed, z-scored over the analysis
  population (n = 8,866).
* PCA: variance explained, loadings (PC1-PC5).
* Gaussian mixture (full covariance, 5 initialisations) on the PC scores retaining >= 90% variance,
  k = 1..8, BIC. Selected k = argmin BIC; if the argmin is at the k = 8 boundary this is reported as
  "BIC keeps favouring more components" (a continuum, not discrete groups) and the k = 8 solution is
  characterised anyway.
* HDBSCAN (sklearn) on the same PC scores, min_cluster_size = 1% of n; noise fraction reported.
* Structure null: the same GMM BIC scan on data whose columns were permuted independently (keeps
  marginals, destroys joint structure).
* Stability: 50 bootstrap resamples. GMM refit at the selected k and used to label everyone, then
  ARI vs the full-data labels. HDBSCAN refit and ARI computed on the resampled participants. Also
  ARI across 10 GMM random initialisations on the full data.
* Characterisation: per cluster n, age, % female, activity; prevalence of every target with Wilson
  95% CI; age- and sex-adjusted odds ratio vs the largest cluster (logistic regression, 95% CI);
  chi-square test. Clusters are **never called subtypes**.
* UMAP (n_neighbors 30, min_dist 0.3, seed) for visualisation only.
* Output: `participant_cluster_assignments__nhanes` (processed; long, one row per participant per
  method).

## 9. Phenotype signatures

For every target, proxy variant and rx-code phenotype, and every wearable feature: z-score the
feature within the phenotype's population and fit OLS `z ~ case + age + age^2 + female`, HC3 robust SE.
The case coefficient is the age- and sex-adjusted standardized difference (SD units) with 95% CI and
p. BH-FDR across the 26 features within each phenotype. Record n cases / controls, raw group means,
label basis, and whether n_cases < 50 (underpowered: minimal detectable |d| ~ 0.5 at n = 34).
Written to processed `phenotype_signatures__nhanes` with `object_id =
signature:nhanes|<phenotype>|<feature>`. `condition_id`: fibromyalgia / migraine / ibs for the rx
groups; `me_cfs` ONLY for the proxy rows, whose `phenotype_label` says "ME/CFS-like proxy"; null for
everything else. `evidence_level` describes the label basis.
`signatures.get_patient_phenotype_signature(condition)` reads every `phenotype_signatures__*`
partition and returns per dataset: definition, label basis, n, top features with effect sizes and
CIs, null results and caveats, or `UNKNOWN / NOT AVAILABLE` with a reason.

## Known limitations (stated in advance)

Cross-sectional (except mortality); self-reported labels; rx reason codes identify treated people in
one cycle, not prevalence; no HR/HRV; wrist accelerometry with a single MIMS cut-point; algorithm
sleep proxy; survey design (weights, strata, PSU) not used for prediction metrics, so metrics
describe this unweighted wearable-valid sample, not the U.S. population; the valid-wear subset is not
re-weighted for non-wear selection.

## 10. Deviations from this plan

_(appended during/after the analysis; nothing above this section was edited after the analysis started)_

1. **Metrics stage re-run (no refit).** The first full run finished the cv, temporal, perm and shuffle stages
   and then crashed at the start of the metrics stage: joblib could not unpickle scipy's `logit` ufunc
   sent from `__main__`. `logit`/`expit` were replaced by mathematically identical numpy functions (unit tests
   compare the calibration fit with statsmodels to 1e-6), the entry point now runs through the importable
   module, and only the metrics stage was re-run on the cached out-of-fold predictions. No model was refitted
   and no specification changed.
2. **HDBSCAN stability reporting.** HDBSCAN labelled every participant as noise (0 clusters). The ARI between
   two all-noise labelings is trivially 1.0, so the stability code was changed to record ARI as not
   estimable (NaN) in that case and to report the number of clusters per bootstrap replicate. The
   unsupervised run was stopped and repeated with this change; PCA/GMM/UMAP use fixed seeds and are
   unaffected. The HDBSCAN bootstrap refits on the unique resampled participants (duplicated points
   distort density estimates), consistent with §8.
3. **Descriptive additions not in the plan** (none changes a pre-specified estimate or decision rule):
   a robustness qualifier in the results headline (agreement of the Test-2 CI across the three model
   families and the temporal split); a `context_clinical_vs_demographics` Delta-AUROC; a table describing
   the composition of the ME/CFS-like proxy and what each exclusion diagnosis removed; coefficients of the
   primary model refit on the whole population (descriptive only); an explicit caveat that Rx reason-code
   signatures may reflect the medication itself.
4. **Resolution of the shuffle control.** As planned, 20 shuffles were run; the smallest attainable
   empirical p is 1/21 = 0.048, so a "passing" shuffle control means that the observed increment exceeded
   every shuffled one, not that p is small. Stated in the results.
5. **Independent review (2026-09-23; no estimate changed).** The full pipeline was re-run end to end: every
   out-of-fold prediction, permutation and shuffle draw, and every `nhanes_*.csv` result table was
   reproduced byte for byte (only `nhanes_run_metadata.csv` timestamps differ). Changes made by the reviewer:
   (a) *Query fixes* in `signatures.get_patient_phenotype_signature`: text matching now uses whole tokens of
   the phenotype id/label with sensitivity/contrast qualifiers removed (previously "arthritis" returned only
   `mecfs_like_proxy_excl_arthritis`, a phenotype that EXCLUDES arthritis, and "health" matched the healthy-
   reference contrast); model-level AUROC rows stored by other producers are reported under
   `model_level_results` instead of as top wearable features (previously the only "top feature" for Long
   COVID was a demographics-only model AUROC); rows without a q value are listed as descriptive, not as null
   results; a missing `null_value` is reported as UNKNOWN instead of being assumed to be 0.
   (b) *Provenance*: `source_version`/`retrieved_at` are now per cycle (previously every row, including all
   2013-2014 participants in `participant_cluster_assignments__nhanes` and the 2013-2014-only Rx signatures,
   cited the 2011-2012 release). `phenotype_signatures__nhanes` gains a `null_value` column (0).
   (c) *Caveats*: the ME/CFS-like proxy rows now carry data-derived caveats (share of cases with PHQ-9 >= 10
   and arthritis; the null predictive increment); the healthy-reference contrast is flagged as an
   extreme-group contrast with inflated effects.
   (d) *Reporting, post hoc and descriptive*: Holm adjustment across the 6 targets of the paired-bootstrap
   p for Test 1 and Test 2. Test 2 for fatigue and depression meets the pre-specified rule (CI > 0) but does
   not survive Holm, is CI > 0 in 1/3 model families and not in the temporal split, so it is reported as
   NOT ROBUST / unconfirmed rather than "supported". The wearable-only-vs-demographics statement now separates
   point estimates (6/6 lower) from CIs (5/6 exclude 0). The positive-control text notes that most of the
   wearable-only mortality AUROC is age. A note that shuffled increments are negative on average (so beating
   every shuffle is a weak test). An HDBSCAN parameter-sensitivity table (`nhanes_hdbscan_sensitivity.csv`).
   Legend placement in two draft figures.
