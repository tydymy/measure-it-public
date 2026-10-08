# Analysis plan: serum cortisol and plasma immune mediators in Long COVID (MY-LC, Klein et al., Nature 2023)

**Dataset id:** `klein2023_mylc_ml_table`. It is rank 3 in `docs/LAB_DATASET_DISCOVERY.md`, tied with rank 2 within
0.01, so both tied datasets are analysed. **Written 2026-09-24, before any analyte value was compared between
groups.** Module: `src/measure_it/labs/klein2023_mylc_ml_table.py`. Nothing here may change after results are seen,
and anything else is exploratory.

## What was seen before this plan

* **The file, at schema level only.** Supplementary Table 3 (MOESM4, `Sheet1`) has 185 rows x 7,431 columns, one row
  per person, keyed by `x0_LC_ID`. What was looked at:
  * the `x0_*` column names;
  * the counts of `x0_Censor_Cohort_ID` (3 = 101, 1 = 42, 2 = 42);
  * cross-tabs of cohort ID against **non-biomarker** columns: `x0_Censor_Complete`, the infection-test fields, the
    steroid and pituitary/adrenal flags, sex and hospitalisation;
  * how many people have non-missing cortisol, draw time, age and BMI in each cohort.
* **Cohort codes.** 1 = healthy uninfected (no infection-test fields). 2 = convalescent (35 of 42 have a positive
  test). 3 = Long COVID. Keeping `x0_Censor_Complete == 0` leaves 40 / 39 / 99, which matches the paper's analysed
  groups.
* **The paper's claims**, read from its text and kept separate from our computation:
  * The MY-LC groups differ in median cortisol (Kruskal-Wallis P < 0.0001).
  * Cortisol was lower in LC and similar in HC and CC.
  * In the Gale-Shapley-matched set, "cortisol alone achieved an AUC of 0.96 (95% CI 0.93-0.99)".
  * The cortisol difference "remained significant after accounting for variations in demographics and
    sample-collection times".
  * C4b, CCL19, CCL20, galectin-1, CCL4, APRIL and LH were higher in LC, and IL-5 was lower.

## Primary analysis (one)

* **Question:** does serum cortisol add information about Long COVID status beyond time of blood draw, age, sex and
  BMI, when evaluated on held-out people?
* **Sample:**
  * people with `x0_Censor_Complete == 0`;
  * cases: LC (cohort 3) = 1; controls: HC and CC pooled (cohorts 1 and 2) = 0;
  * complete cases on cortisol (`x1_Cytokines_Cortisol_Obs_Conc_ng/mL_GG21.1_ML`), draw time (`x0_Sample_Time_Min`),
    age, sex and BMI.
* **Models:** all are L2 logistic regressions with C = 1.0, standardised inside each fold.
  1. Covariates: draw time (linear), age, sex, BMI.
  2. log10(cortisol) alone.
  3. Covariates + log10(cortisol).
* **Evaluation:**
  * Stratified 5-fold CV x 20 repeats (seed 20260923 + repeat).
  * AUROC per repeat on the pooled out-of-fold predictions, averaged over repeats.
  * **Primary statistic: Delta AUROC = model 3 minus model 1**, with a 95% CI from a 2,000-resample paired,
    group-stratified person bootstrap.
  * Model 2 is tested against a label-permutation null of 1,000 shuffles, with the full CV re-run for each shuffle.
* **Decision rule:** "cortisol adds information beyond draw time and demographics" only if the Delta CI lower bound
  is above 0 **and** the model-2 permutation p is below 0.05. Otherwise "not shown".
* **Interpretation, fixed now:** a positive result replicates a published group difference in the same released
  rows. It is not an independent replication. It also does not make cortisol a diagnostic test:
  * the controls are volunteers;
  * cortisol is diurnal and assay-dependent;
  * there is one draw per person.

## Pre-specified secondary analyses (all reported whatever they show)

* **S1 per-analyte effects.** Each of the 144 `x1_Cytokines_*` analytes, LC vs pooled controls, on the primary
  analysis set with that analyte non-missing:
  * Hedges' g on log10(value + half the smallest positive value);
  * a stratified bootstrap CI (2,000 resamples);
  * a two-sided Mann-Whitney p;
  * BH q across the 144.
* **S2 cortisol by control type.** Fit-free AUROC = P(control cortisol > LC cortisol), computed for LC vs HC and for
  LC vs CC separately, each with a stratified bootstrap CI and a 1,000-shuffle permutation p.
* **S3 multi-analyte panel.** All 144 analytes (log10, with training-fold median imputation) in an L2 logistic model
  with C = 0.1, using the same CV, a bootstrap CI and a 1,000-shuffle permutation null.
* **S4 sensitivity.** The primary repeated after also excluding anyone with `x0_Censor_Oral_Steroid == 1` or
  `x0_Censor_Pit_Adre_Dysfunction == 1`.
* **S5 sensitivity at 90% specificity** for models 2 and 3.
* **S6 descriptive.** The single-feature, fit-free cortisol AUROC on the primary sample, placed next to the paper's
  0.96 (which was computed on a matched subset).

## Outputs

* Processed tables:
  * `participants__klein2023_mylc_ml_table`;
  * `participant_labs__klein2023_mylc_ml_table` (long format; 144 plasma analytes);
  * `phenotype_signatures__klein2023_mylc_ml_table`.
* Result tables: `results/tables/klein2023_mylc_ml_table_*`.
* Report: `results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md`.

## Deviation 1 (2026-09-24, after the first run; forced by the released value scale, not by the result)

**What went wrong.** The plan specified log10(cortisol) and log10(analyte + half the smallest positive value). The
column names say concentrations (for example `..._Obs_Conc_ng/mL_...`). In fact the released `x1_Cytokines_*`
values are **batch-integrated z-scores**: every analyte has mean 0 and SD 1 over the analysis set, and all 144 have
negative values. This matches the paper's Methods: "After batch integration, each feature was z-scored".

**Effect on the first run.** The log transform turned every negative value into NaN, and complete-case filtering
then silently dropped those people. For cortisol, 100 of 177 values are negative. The first run's primary sample was
therefore a biased subset:
* **7 LC vs 68 controls** instead of about 98 vs 79. Low cortisol is the Long COVID direction, so the LC cases were
  removed preferentially;
* results: Delta 0.050 (-0.055 to 0.172), cortisol-only CV AUROC 0.692, permutation p 0.034.

That run is **invalid**. Its outputs are kept for transparency, renamed `*_INVALID_log10_run_*`, in
`results/tables/_invalid_runs/`.

**Change.** Every analyte enters as the released z-score, untransformed. Everything else in the plan is unchanged:
models, CV, bootstrap, permutations, decision rule, secondaries and seeds.
* The fit-free AUROCs (S2, S6) are invariant to any monotone transform, so the only change for them is that nobody
  is dropped.
* Hedges' g in S1 is computed on the released z-scores.
