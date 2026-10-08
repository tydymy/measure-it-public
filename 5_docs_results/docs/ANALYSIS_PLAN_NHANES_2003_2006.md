# Analysis plan — NHANES 2003-2006 hip accelerometry: replication of the NHANES 2011-2014 Tests 1/2/6

_Written 2026-09-24 while the hip-accelerometry features were being derived and BEFORE any model on 2003-2006 was
fitted or any target was tabulated against a hip feature. Owner module: `src/measure_it/wearables/nhanes0306_models.py`.
Results: `results/NHANES_2003_2006_RESULTS.md`, `results/tables/nhanes0306_*.csv`. Deviations are appended in §8._

The 2011-2014 plan (`docs/ANALYSIS_PLAN_NHANES.md`) is re-run as closely as the data allow. What differs, and why:

## 1. Data and population

`participants__nhanes0306`, `participant_clinical_features__nhanes0306`, `participant_mortality__nhanes0306`,
`participant_wearable_features__nhanes0306` (joins on `participant_id = nhanes0306:<SEQN>` only). Population:
reliable monitor, passes the Troiano valid-wear rule (>= 4 days with >= 600 wear minutes), aged >= 18, not pregnant
at the MEC exam. The device is a hip-worn ActiGraph 7164 (counts/min, waking hours) — never pooled with, and never
put on the same scale as, the 2011-2014 wrist MIMS; the two cohorts are compared only at the level of conclusions.

## 2. Targets (same definitions; availability differs)

| id | 2003-2004 (C) | 2005-2006 (D) | note |
|---|---|---|---|
| functional_limitation | yes | yes | PFQ composite, ages >= 20 |
| fair_poor_health | yes | yes | HSD010 in {4, 5} |
| fatigue (DPQ040 >= 2) | — | yes | PHQ-9 starts 2005-2006 |
| depression_phq9 | — | yes | |
| mecfs_like_proxy | — | yes | fatigue + limitation + no exclusion diagnosis; COPD is not asked in 2003-2006 (the item is missing, and missing never excludes); arthritis type from MCQ190 (rheumatoid counts; no psoriatic category). ME/CFS-LIKE PROXY, never ME/CFS |
| mortality (positive control) | yes | yes | public-use LMF, follow-up to 2019-12-31 (longer than for 2011-2014) |

## 3. Feature sets

* **demographics**: identical construction (age — top-coded at 85 in these cycles —, female, 5 race indicators with
  the 2011-2014 levels; RIDRETH1 has no Non-Hispanic Asian category, so that indicator is constant 0 and Asian
  participants fall in "Other race including multi-racial"; income-poverty ratio; education).
* **clinical**: the 2011-2014 `clinical` list with the same names and transforms. Creatine kinase is not measured in
  2003-2006 (all-missing; dropped by the imputer). Comorbidity count uses the same 18 items (gout, celiac and
  psoriasis are not asked in these cycles and count as 0 — stated, not imputed).
* **clinical_crp** (secondary): clinical + log1p(CRP, mg/dL) (CRP exists in 2003-2006 only).
* **wearable (hip)**, 23 features fixed now: mean_daily_counts (log1p), sd_daily_counts (log1p), cv_daily_counts,
  weekend_minus_weekday_counts, mean_wear_min, mean_cpm_wear (log1p), sedentary_min_per_day, sedentary_fraction,
  light_min_per_day, mvpa_min_per_day (log1p), vigorous_min_per_day (log1p), mvpa_bout10_min_per_day (log1p), astp,
  satp, mean_active_bout_min (log1p), mean_sedentary_bout_min (log1p), interdaily_stability,
  intradaily_variability, m10_counts (log1p), l5_counts (log1p), relative_amplitude, m10_onset_h and l5_onset_h
  (hours after 00:00 — the waking-wear profile has no nocturnal trough to cut at). `steps_per_day` exists only in
  2005-2006 and is excluded from the modelled set (descriptive only). No heart rate/HRV exists.
* Sets: demographics, clinical, wearable, demographics_wearable, clinical_wearable (+ clinical_crp,
  clinical_crp_wearable). Predictor exclusions identical to the 2011-2014 plan (no PFQ, PHQ, self-rated health,
  self-reported activity or sleep item is ever a predictor; `comorbidity_count` dropped for the proxy).

## 4. Models, splits, metrics, controls — unchanged

`lr`, `pen_lr` (primary), `hgb` exactly as in `measure_it.wearables.models`; repeated stratified 5-fold CV x 5
repeats; secondary temporal split train 2003-2004 -> test 2005-2006 (targets present in both cycles); AUROC, AUPRC,
Brier, calibration, participant bootstrap (B = 1000) with paired Deltas; label permutation (200, `pen_lr`,
wearable-only); wearable-row shuffle (20) for the clinical -> clinical+wearable increment; mortality as the
positive control.

## 5. Test definitions and decision rule (as 2011-2014)

Test 1: demographics_wearable vs demographics (and wearable-only vs demographics). Test 2: clinical_wearable vs
clinical (primary) and clinical_crp_wearable vs clinical_crp (secondary). SUPPORTED = primary-model 95% CI of the
Delta excludes 0 AND the matching control p <= 0.05 (permutation for Test 1's wearable-only signal, shuffle for
Test 2). Holm across targets is reported descriptively, as in the 2011-2014 review.

## 6. Replication criterion (pre-specified)

For each target present in both cohorts and each test, the 2011-2014 verdict (from `results/tables/nhanes_*.csv`,
unchanged) is compared with the 2003-2006 verdict under the same rule: **replicated** = same verdict (SUPPORTED in
both, or NULL/not supported in both); **not replicated** = different verdicts. Also reported: whether each cohort's
point estimate lies inside the other's 95% CI. Because devices differ, a replicated verdict means "the same
conclusion holds with a different accelerometer", not that effect sizes are comparable.

## 7. Known limitations

Hip vs wrist, counts vs MIMS, waking-wear vs 24-h wear; PHQ-9-based targets only in one cycle (smaller n; the
ME/CFS-like proxy will have few cases); different era (2003-2006) and longer mortality follow-up; no survey weights.

## 8. Deviations from this plan

_(appended during/after the analysis)_

1. **Metrics stage re-run (no refit).** The first metrics stage crashed because, for the ME/CFS-like proxy (36
   cases), some penalised fits returned constant predictions, leaving the calibration slope undefined in every
   bootstrap resample; `numpy.quantile` of an empty array raised. CIs of a metric with no finite bootstrap value
   are now reported as missing; only the metrics stage was re-run on the cached out-of-fold predictions. No
   model, feature or rule changed.
2. **Replication summary split (post hoc, descriptive).** §6 compares verdicts for every test. In the report the
   headline count uses the two comparisons that carry the SUPPORTED/NULL rule (demographics+wearable vs
   demographics; clinical+wearable vs clinical); the wearable-only vs demographics comparison (part of "Test 1" in
   §5) is reported in full as "Test 1b" with a direction-only verdict (lower / higher / no difference), because the
   decision rule of §5 does not define SUPPORTED for a comparison that is expected to be negative.
3. **Descriptive additions not in the plan**: number of model families (lr / pen_lr / hgb) with CI > 0 per test and
   cohort; a "Reading" section; a note that the 2011-2014 Test 2 results for fatigue and depression were already
   labelled NOT ROBUST in `results/WEARABLE_NHANES_RESULTS.md`. None changes an estimate or a verdict.
