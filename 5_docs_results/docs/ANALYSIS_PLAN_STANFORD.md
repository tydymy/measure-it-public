# Analysis plan — Stanford Snyder Lab wearable analyses (SPEC Phase 1; Tests 1, 2 and 6)

Written 2026-09-23, **before** any of the analyses below were run. Code:
`src/measure_it/wearables/stanford_models.py` (`uv run python -m measure_it.wearables.stanford_models`).
Results: `results/WEARABLE_STANFORD_RESULTS.md`, `results/tables/stanford_*.csv`,
`results/figures/drafts/stanford_*`, processed table `phenotype_signatures__stanford`.
Seed: `config.SEED` (20260923) for every stochastic step. At most 6 worker processes.

Inputs (all produced by `measure_it.ingestion.stanford_wearables`; caveats read from
`data/raw/stanford_longcovid_wearables/DATA_AUDIT.md`):
`participants__stanford_covid`, `participant_conditions__stanford_covid`,
`participant_wearable_features__stanford_covid`, `participant_wearable_daily__stanford_covid`.

What was looked at before writing this plan: table schemas, row counts, missingness
and eligibility counts (who has a reference date / baseline / pre-infection data),
and the already-published exploratory rows of `stanford_reproduction_benchmark.csv`
(univariate whole-record RHR AUROCs 0.45-0.58, none with permutation p < 0.2; a
demographic-only RF AUROC 0.60). No model below was fitted before this plan.

---------------------------------------------------------------------------

## Part A — Uwakwe 2025 Long COVID discrimination (Test 1 on an invisible-illness label; Test 2 increment)

**Population.** All 126 participants of `stanford_longcovid_uwakwe2025` (31 Long COVID,
95 COVID-19 recovered without Long COVID). Unit = participant.

**Target.** `long_covid_label` (self-reported symptoms >= 12 weeks after a confirmed
COVID-19 diagnosis; CDC definition; end-of-study survey). This is a self-report label,
not a clinical diagnosis.

**Feature sets (fixed in advance).**

| set id | features | role |
|---|---|---|
| `demographics` | age bin (ordinal 0-5 from 10-year bins), female, BMI, White (vs any other/unreported), vaccinated before infection | primary baseline ("demographic/clinical": vaccination is a clinical-history item, not a symptom) |
| `wearable` | rhr_mean_record, rhr_night_mean_record, rhr_sd_daily, rhr_iqr_daily, rhr_rmssd_daily, rhr_trend_bpm_per_30d, cusum_alarms_per_30d | primary wearable physiology (whole-record resting-HR features; NOT aligned to infection: the release has no diagnosis date) |
| `demographics_plus_wearable` | union of the two | primary comparison |
| `demographics_age_sex` | age bin, female | sensitivity (the benchmark's exploratory demographic set) |
| `wear_pattern` | hr_minutes_per_day_mean, wear_day_fraction, record_days_span | sensitivity / artefact check: device sampling density and engagement, not physiology |
| `symptoms_reference_LEAKAGE` | the 7 symptom indicators of the paper's RFSM model (`RFSM_FEATURES` in the ingestion module) | **leakage-labelled reference only.** Symptoms define the label; never part of the primary comparison |

Comorbidity indicators are excluded (missing for 84/126). Missing values (BMI 9,
rhr_night_mean_record 4, cusum_alarms_per_30d 2) are median-imputed inside each training
fold. No feature selection.

**Models (fixed hyper-parameters; no tuning on test folds).** Every model is a
pipeline `median imputer -> standard scaler -> classifier` fitted on the training fold only.
1. `logistic`: unpenalised logistic regression (C = inf, lbfgs).
2. `l2_logistic`: L2 logistic; C chosen by an inner 5-fold stratified CV on the training
   fold only (11 values, 1e-3..1e2, log-loss). **Primary model.**
3. `random_forest`: 300 trees, min_samples_leaf = 3, max_features = sqrt, fixed seed.

**Validation.** 10 x 5 repeated stratified K-fold (participant = unit; each participant is
in exactly one test fold per repeat). Per repeat, metrics are computed on the pooled
out-of-fold (OOF) predictions of all 126 participants; the reported point estimate is the
mean over the 10 repeats.

**Metrics.** AUROC (primary), AUPRC (chance = 31/126 = 0.246), Brier score (and the
no-skill Brier of predicting the prevalence), calibration intercept (calibration-in-the-
large, logistic recalibration with logit(p) as offset) and slope, 5-bin calibration table,
sensitivity and specificity at the pre-specified threshold = training-fold prevalence.

**Uncertainty.** Participant-level bootstrap (2,000 resamples of participants; for each
resample the metric is recomputed for every repeat's OOF predictions and averaged) ->
percentile 95% CI. This CI covers evaluation-sample variability, not model-refit
variability; the across-repeat range is reported next to it.

**Primary comparisons (pre-specified).**
* Test 1: `wearable` (l2_logistic) AUROC vs chance (label-permutation p) and vs
  `demographics` (paired bootstrap ΔAUROC).
* Test 2 (multimodal increment): ΔAUROC = `demographics_plus_wearable` − `demographics`
  (l2_logistic), paired participant bootstrap CI, plus a **wearable-block permutation
  null** (rows of the 7 wearable columns permuted across participants, labels and
  demographics fixed; 1,000 permutations; full 10x5 procedure each) -> one-sided p.
* All three models are reported for all primary sets; the RF and unpenalised logistic
  results are secondary.

**Negative control (Test 6).** Label-permutation null: 1,000 permutations of the Long
COVID label, each re-running the full 10x5 CV for the same model/feature set;
statistic = mean-over-repeats AUROC; one-sided p = (1 + #null >= observed)/(1 + 1000).
Run for {demographics, wearable, demographics_plus_wearable} x {3 models}. Expected: null
AUROC centred at 0.5.

**Univariate feature effects (for the signature table).** For each of the 7 wearable
features: Hedges' g (Long COVID minus no Long COVID) with 5,000-resample stratified
bootstrap CI; two-sided label-permutation p (10,000) of the mean difference; BH q across
the 7 features; per-feature AUROC. Secondary: age-bin- and sex-adjusted difference
(OLS feature ~ LC + age bin + female, HC3 SE), expressed in SD units of the feature.

**Interpretation rule (fixed now).** A wearable signal is called "evidence-supported in
this dataset" only if the l2_logistic `wearable` AUROC has permutation p < 0.05 AND its
95% CI excludes 0.5. A multimodal increment is claimed only if the ΔAUROC CI excludes 0
AND the wearable-block permutation p < 0.05. Otherwise the result is reported as a null
result. n = 126 with 31 cases: CIs are expected to be wide (about ±0.1 AUROC).

---------------------------------------------------------------------------

## Part B — Acute post-infection trajectories (Mishra 2020 + Alavi 2022; NOT Long COVID)

**Population.** COVID-19-positive participants of `stanford_covid_mishra2020` (32) and
`stanford_covid_alavi2022` (84) with a reference day (symptom onset, else test/diagnosis
date; `reference_basis` kept) and a valid pre-infection baseline (>= 7 RHR days in days
−42..−15; the ingestion rule). Pooled analysis is primary; per-dataset and per-device-class
(Fitbit minute vs HealthKit sparse) results are reported as sensitivity. Mishra
"other illness" and Alavi unlabelled participants are not used (the latter are not
verified negatives).

**Daily deviations.** Within-person RHR z = (daily RHR − baseline mean)/baseline SD and step
z likewise, baseline = days −42..−15 relative to the reference. I recompute these from
`rhr_mean`/`steps` with the ingestion's own `baseline_stats`/`add_deviations` so the same
code produces the true-onset and placebo deviations; a test asserts equality with the
table's `rhr_z` for true onsets.

**Trajectories.** Mean RHR z and step z by day −28..+60 relative to onset, participant-level
bootstrap (2,000) 95% CI, n contributing per day.

**Acute elevation.** Peak daily RHR z in days −7..+14 > 2 (ingestion definition). Fraction
with elevation, Wilson CI.

**Recovery time.** Among participants with acute elevation: days from the reference day
to the first of 3 consecutive observed days with z <= 1 at or after max(peak day, 0);
right-censored at the last observed day (<= +90). Kaplan-Meier (statsmodels
`SurvfuncRight`), median with CI, events/censored counts. Recomputed from my z and
checked against the table's `days_to_rhr_recovery`.

**Negative control — placebo onsets (pre-specified).** For each participant, a placebo day
P is eligible iff (all in days relative to the true reference T):
1. P <= T − 29, so the placebo analysis window [P − 42, P + 14] (its baseline + its acute
   window) ends at or before T − 15, i.e. it never overlaps the true window
   [T − 14, T + 90] (the 14-day pre-symptomatic window onward);
2. its own baseline [P − 42, P − 15] has >= 7 RHR days (same rule as the true onset);
3. its acute window [P − 7, P + 14] has >= 1 observed RHR day (same rule as the true onset);
4. Mishra only: no recorded symptom-onset / diagnosis / recovery date of any episode of
   that participant falls in [P − 7, P + 14] (other illnesses).
Deviations for a placebo are computed exactly as for the true onset, relative to P.
* Primary placebo summary: per participant, the average over ALL eligible placebo days
  (each participant weight 1), compared with the same participant's true onset (paired, in
  the placebo-eligible subset). Paired difference with participant bootstrap CI; two-sided
  sign-flip permutation p (10,000).
* One random eligible placebo per participant (seeded) is used for the detection model and
  as a sensitivity for the elevation fraction (exact McNemar test).
* Endpoints compared: mean RHR z days 0..14, mean step z days 0..14, acute-elevation
  fraction, trajectory curves (placebo days −28..+14 only, because placebo windows must end
  before T − 15), and recovery time with both true and placebo follow-up truncated at +14.
* Expected: elevation fraction and mean z near their chance levels for placebo onsets.
  If placebo elevation is common (a single z > 2 in a 22-day window occurs by chance), the
  true-onset elevation fraction must be read against it, not against 0.

**Detection analysis (retrospective window discrimination, not prospective early warning).**
Units = windows: one true-onset window + one random eligible placebo window per
placebo-eligible participant (1:1). Features, computed from each window's own baseline:
* window A (acute, days 0..+7): mean RHR z, max RHR z, mean step z;
* window B (pre-symptomatic, days −7..−1): the same three features (separate model).
Windows with < 3 observed RHR days in the feature window are excluded (same rule for true
and placebo). Model: L2 logistic, fixed C = 1.0 (3 features; no tuning), imputer + scaler
fitted in the training fold. CV: 10 x 5 StratifiedGroupKFold grouped by participant (no
participant in both train and test). Also the model-free AUROC of mean RHR z alone.
Participant-cluster bootstrap CI (2,000). Null: 1,000 within-participant label swaps
(each participant's true/placebo labels randomly exchanged), full CV each. Sensitivity:
all eligible placebo windows instead of one.

**Temporal-leakage checks (asserted in code, written to
`results/tables/stanford_temporal_leakage_checks.csv`).**
1. Every baseline day (true and placebo) is <= reference − 15, strictly before any day used
   for features or trajectories of that window (feature windows start at reference − 7).
2. Every placebo analysis window ends before the true pre-symptomatic window (<= T − 15).
3. Grouped CV: no participant contributes windows to both train and test folds.
4. Imputation/scaling/tuning are fitted inside training folds only (sklearn pipelines).
5. The reference day itself (symptom onset / test date) is known only retrospectively;
   window A uses days after onset, so this is window discrimination, not a claim that the
   infection could be flagged at onset.

---------------------------------------------------------------------------

## Part C — `phenotype_signatures__stanford`

One row per (dataset, phenotype, feature); columns: dataset_id, phenotype_id,
phenotype_label, condition_id, feature, effect_size, effect_unit, ci_low, ci_high,
p_value, q_value, n_cases, n_controls, label_basis, object_id
(`signature:<dataset_id>|<phenotype_id>|<feature>`), evidence_level + provenance.
* Uwakwe rows -> condition_id `long_covid`: 7 univariate Hedges' g rows (+ 7 age/sex-adjusted
  rows) and model rows (l2_logistic AUROC for demographics, wearable,
  demographics_plus_wearable, and the ΔAUROC increment). The leakage-labelled symptom model
  is NOT written to the signature table.
* Acute rows -> condition_id null, phenotype_label 'post-infection recovery trajectory
  (acute COVID-19 cohort; not a Long COVID label)': paired true-vs-placebo differences,
  elevation fractions, KM median recovery, detection AUROCs.
* q_value: BH within each (dataset_id, phenotype_id) family over rows with a p-value.
* evidence_level (this table): `supported_single_dataset` (q < 0.05 and CI excludes the
  null value), `nominal_single_dataset` (p < 0.05 but not the former), `null_single_dataset`
  (p >= 0.05), `descriptive` (no test).

## Deviations from this plan
Recorded 2026-09-23 after the first full run of Part B (none of these changes the primary endpoints, their windows, the models or the tests):

1. Detection table: added descriptive univariate AUROCs of each of the three pre-specified window features (no fitting, no test), to show which feature carries the window-A model. Added after the first run showed the model and the model-free RHR score disagree; used only to describe which feature carries the model, never as a test.
2. Endpoints table: added DESCRIPTIVE mean RHR z and mean step z over days 15..90 for true onsets only (no placebo comparator is possible because placebo windows end at +14). Added after the trajectory plot showed RHR z above the placebo band in weeks 3-5; not tested, not written to the signature table.
3. Detection windows: the one random placebo window per participant is drawn among that participant's eligible placebo onsets that also have >= 3 observed RHR-z days in the feature window, and only complete true/placebo pairs are analysed (so the within-participant label-swap null is well defined). The plan said 'one random eligible placebo window' and '>= 3 observed days' without fixing the order of the two rules.

## Post-review changes (independent review, 2026-09-23; after all results above were seen)

None of these changes a pre-specified endpoint, window, model or test, and every pre-specified number is unchanged
(the full pipeline was re-run and diffed against the pre-review outputs). The additions below are sensitivity or
descriptive analyses, labelled as such in the tables and in `results/WEARABLE_STANFORD_RESULTS.md`.

1. **Signature table (Part C).** The demographic-only and demographics+wearable model AUROC rows are no longer written.
   `measure_it.wearables.signatures.get_patient_phenotype_signature("long_covid")` lists every row with q < 0.05 as a top
   feature of the wearable signature, so it returned `model_auroc:demographics:l2_logistic` (a demographics model) as the
   only "feature" of the Long COVID wearable signature. Those AUROCs are quoted as context in each Uwakwe row's `caveats`.
   The Uwakwe BH family therefore has 16 rows (7 unadjusted + 7 adjusted features, wearable-only AUROC, wearable increment).
   Added the columns the query reads (`feature_label`, `effect_measure`, `phenotype_definition`, `is_proxy`, `caveats`),
   which it previously reported as UNKNOWN.
2. **Recovery time (Part B).** Kaplan-Meier rows timed from the acute peak (true vs placebo, follow-up capped at +14),
   plus a participant-weighted KM over every elevated placebo onset (participant-bootstrap CI), and the peak-day
   distribution. Timed from the reference day, the true-vs-placebo contrast mostly reflects where the peak falls in
   days -7..+14 (placebo peaks spread across the window), not recovery speed.
3. **Univariate Uwakwe features (Part A).** A studentized (Welch-t) permutation p over the same 10,000 permutations. The
   pre-specified raw mean-difference permutation test is anti-conservative when the smaller group has the larger
   variance. The pre-specified p stays the primary column.
4. **Negative-control sensitivity (Part B).** Placebo onsets are enumerated from the start of each record, so more of
   them have short baselines. The elevation fraction and the 0..14-day mean RHR z are repeated with true AND placebo
   baselines of >= 20 RHR days.
5. **Wear-time check (Part B, descriptive).** True-minus-placebo change, days 0..14 vs own baseline, in (a) the share of
   calendar days inside the record span that are valid wear days and (b) HR minutes per valid wear day. Tests whether the
   post-onset step drop is reduced walking or reduced wearing.
6. **Robustness.** `km_summary` drops and counts rows with a missing recovery time/event instead of casting NaN to an
   integer.
