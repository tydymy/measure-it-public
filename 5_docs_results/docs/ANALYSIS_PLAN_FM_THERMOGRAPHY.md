# Analysis plan: infrared thermography in fibromyalgia (`fm_thermography`)

**Written and saved 2026-09-24, before any case-control comparison was computed.** The primary analysis is the
one locked in `docs/DEVICE_DATASET_DISCOVERY.md` (section "Pre-specified primary analyses", item 2) and stored in
`results/tables/device_dataset_candidates.csv` (`recommended_primary_analysis`). It is copied here unchanged. The
secondary analyses below are added now and are fixed before any group statistic is seen. Anything not listed here is
labelled **exploratory** wherever it is reported.

Question: *can a portable infrared camera help uncover fibromyalgia (FM), a clinically diagnosed invisible illness?*
Operationally: does resting regional skin temperature separate women with FM from women without FM, and does it add
anything to age and BMI?

## What had been seen before this plan was written

* **Schema only.** 178 rows x 27 columns, one sheet (`Hoja5`). Columns: `Participant`, `Group`, `Age`, `BMI`,
  `{Neck, Upper_back, Lower_back, Chest, Knee, Elbow}_{min, max, ave}`, `Civil_status`, `Employment_status`,
  `Studies`, `Tobacco`, `Alcohol` (integer codes 1-4/1-5, labels not in the file).
* Label counts: Group 1 = 86 (all IDs `FM###`), Group 2 = 92 (all IDs `CG###`). One FM row lacks the knee and elbow
  temperatures, so complete cases on the six `*_ave` columns are 85 FM vs 92 controls. No other missing values.
* Pooled (not by group) min/max of each column, to check units: temperatures 25.9-36.3 degC, age 38-70 y,
  BMI 19.5-39.9 kg/m2.
* The paper's abstract and Methods (article XML retrieved to `data/raw/fm_thermography/`). The abstract reports
  Mann-Whitney differences at the lower back and knees (p < 0.05) "that did not reach a minimum of clinically
  detectable change", no difference elsewhere, and concludes thermography "is not an effective supplementary
  assessment tool". The Results tables (group means) were **not** read before this plan; they are read afterwards
  only for the concordance check in S1.
* From Methods: all participants are **post-menopausal women** (so sex is constant and there is no sex covariate);
  FM = ACR 2010 criteria assessed by rheumatologists, recruited from FM associations and specialised units;
  controls = age-matched women without FM symptoms recruited by advertisement and snowball sampling; measurements
  15 min after undressing in an 8 m2 room kept at 24 degC / 44% RH, 3-6 pm; FLIR E60bx, 1 m, emissivity 0.98; knee
  and elbow values are the mean of both sides.

## Primary analysis (locked 2026-09-24; unchanged)

* **Sample:** complete cases on the six `*_ave` temperatures: 85 FM (y = 1) vs 92 controls (y = 0). `*_min` and
  `*_max` are not used.
* **Models** (all L2 logistic regression, C = 1.0, features standardised inside each training fold):
  1. `demo`: age + BMI;
  2. `device`: `Neck_ave`, `Upper_back_ave`, `Lower_back_ave`, `Chest_ave`, `Knee_ave`, `Elbow_ave`;
  3. `demo_device`: both.
* **Cross-validation:** stratified 5-fold, 20 repeats; repeat r uses `StratifiedKFold(shuffle=True,
  random_state=20260923 + r)`. Folds depend only on y, so the three models share identical folds (paired). One row
  per participant, so participant-level CV needs no grouping (no repeated measures).
* **Primary statistic:** Delta AUROC = AUROC(model 3) - AUROC(model 1). AUROC is computed within each repeat on that
  repeat's pooled out-of-fold predictions (all 177 people) and averaged over the 20 repeats.
* **Uncertainty:** 95% percentile CI from a paired person bootstrap, 2,000 resamples (seed 20260923), resampling
  people with replacement within each group (stratified, so both classes are always present); the same resampled
  people are used for both models, the OOF predictions are held fixed, and each resample's statistic is the
  repeat-averaged AUROC difference.
* **Permutation null:** device-only AUROC (model 2) against 1,000 label shuffles (seed 20260923), each shuffle
  re-running the full 20 x 5 CV with the same fold scheme applied to the shuffled labels; one-sided
  p = (1 + #{null >= observed}) / 1,001.
* **Decision rule:** "thermography adds information beyond age and BMI" **only if** the Delta AUROC CI lower bound
  is > 0 **and** the device-only permutation p < 0.05. Otherwise the result is reported as a replication of the
  authors' published null. Whatever the outcome, the result is reported with the caveat below.
* **Caveat carried with every result:** the room was held at 24 degC by protocol, but per-session room temperature,
  date and season are not in the file; recruitment channels differed by group (associations/clinics vs
  advertisements), so any difference may partly reflect selection; this is one research group's single-site sample
  of post-menopausal women and generalises to nothing else without replication.

## Secondary analyses (fixed now; none can replace the primary)

* **S1. Per-feature standardised differences.** For each of the 18 temperature columns: Hedges' g (FM minus
  control) with a 2,000-resample stratified bootstrap 95% CI, two-sided Mann-Whitney U p, and Benjamini-Hochberg q
  across the 18 (the FDR family). Age and BMI get the same statistics as descriptive context, outside the family.
  Uses all available rows per feature (86 vs 92, or 85 vs 92 for knee/elbow), which matches the paper. Concordance
  with the paper's reported per-region results (direction and p < 0.05 or not) is tabulated after the fact.
* **S2. Age- and BMI-adjusted differences.** For each of the six `*_ave` temperatures, OLS
  `temp ~ FM + age + BMI`; the FM coefficient divided by the pooled within-group SD (an adjusted standardised
  difference), its 95% CI and p, BH q across the six. Reason: adiposity lowers skin surface temperature, so a
  raw difference may be a BMI difference.
* **S3. Discrimination and calibration of all three models:** AUROC and AUPRC (with the bootstrap CI above; AUPRC
  chance = 85/177 = 0.48), Brier score vs the no-skill Brier, calibration-in-the-large and calibration slope
  (averaged over repeats), and a 10-bin reliability plot for model 3 using the per-person mean of the 20 OOF
  predictions.
* **S4. Operating point for a screening use:** sensitivity at 90% specificity for models 2 and 3 (threshold set
  within each repeat on that repeat's OOF control scores; averaged over repeats; bootstrap CI). Reported because
  "uncovering" an illness is a screening question; it is not a decision statistic.
* **S5. Sensitivity of the primary to the feature set:** the same Delta AUROC with all 18 temperatures (`*_min`,
  `*_max`, `*_ave`) in model 3 instead of the six averages. This is the only feature-set variant that will be run.
  It is reported beside the primary whatever it shows and cannot replace it.

Not done, by design: no other classifiers, no tuning of C, no feature selection, no threshold search, no
sub-grouping, no use of the sociodemographic codes (their labels are undocumented in the file).

## Outputs

* `docs/ANALYSIS_PLAN_FM_THERMOGRAPHY.md` (this file), `results/FM_THERMOGRAPHY_RESULTS.md`
* `results/tables/fm_thermography_*.csv`, `results/figures/drafts/fm_thermography_*.png`
* `data/processed/participants__fm_thermography`, `participant_device_features__fm_thermography`,
  `phenotype_signatures__fm_thermography`
* Condition mapping: `measure_it.ontology.normalize.normalize_condition("fibromyalgia")` ->
  `fibromyalgia` (ICD-10-CM M79.7). Measurement class: none of `configs/measurements.yaml` is an exact fit; the
  nearest are `continuous_temperature` (skin temperature signal, but that class is a wearable) and
  `microvascular_function` (skin blood flow). The outputs record `continuous_temperature` as the signal match with
  the modality mismatch stated, and a new class `infrared_thermography` is proposed rather than editing the shared
  config.
