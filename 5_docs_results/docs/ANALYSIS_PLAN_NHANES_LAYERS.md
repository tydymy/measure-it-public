# Analysis plan — NHANES 2011-2014 laboratory layers: layered increment analysis

_Written 2026-09-24, after the extended laboratory catalogue was ingested and the (label-free) overlap map was
computed, and BEFORE any model in this plan was fitted. Owner module: `src/measure_it/wearables/nhanes_layers.py`.
Results: `results/NHANES_LAYERS_RESULTS.md`, `results/tables/nhanes_layers_*.csv`. Deviations made after this
plan was written are appended in §9 and never edited out._

**What was seen before writing.** Layer sizes, cycle coverage, subsample weights, detection rates and the overlap
counts (`nhanes_lab_layers.csv`, `nhanes_layer_overlap.csv`); no outcome was tabulated against any lab layer and no
model with a lab layer was fitted. The existing NHANES analysis (docs/ANALYSIS_PLAN_NHANES.md) is known: wearable
adds over demographics for 5/6 targets, and over the clinical set robustly only for functional limitation.

## 1. Population, targets, joins

Population and targets are exactly those of `docs/ANALYSIS_PLAN_NHANES.md` §1-2 (via
`measure_it.wearables.nhanes_cohort.load_cohort`): wearable-valid adults aged >= 18, not pregnant (n = 8,866);
targets `functional_limitation`, `fatigue` (DPQ040 >= 2), `fair_poor_health`, `depression_phq9` (PHQ-9 >= 10),
`mecfs_like_proxy` (ME/CFS-LIKE PROXY, never ME/CFS), `mortality` (positive control). Lab layers are joined on
`participant_id = nhanes:<SEQN>` only, from `participant_labs_extended__nhanes`, `participant_labs__nhanes` and
`participant_layer_membership__nhanes`.

**Restriction rule (no whole-layer imputation).** Every comparison is restricted to participants in the target's
population who HAVE every layer involved (>= 1 measured value in the layer; membership table). A subsample layer is
therefore only ever analysed among people sampled for it. Missing individual analytes WITHIN a layer a person has
are imputed inside the training fold (median + missingness indicator), as in the existing pipeline. n and cases are
reported for every comparison. A comparison with < 30 cases in its restricted population is reported as "not
estimable" and not fitted.

## 2. Layer features (label-free rules, computed once on each layer's members in the analysis population)

* Analytes: `analyte_role == "analyte"` (plus urinary creatinine, `URXUCR`, as one extra feature in layers built from
  urine specimens, to adjust for dilution). When a participant has the same variable in a main-sample and a
  special-sample file, the main-sample value is used.
* Eligible if measured in >= 70% of the layer's members (removes cycle-only or sex-only variables inside a layer,
  e.g. vaginal-swab HPV types, so no structurally unsampled sub-block is imputed).
* Continuous analytes with a detection-limit comment code: eligible only if detected (not below LOD) in >= 50% of
  measured values. Released values (NCHS fill LLOD/sqrt(2) below the limit) are used; no re-imputation.
* Continuous analytes need >= 10 distinct values. Transform: if all values >= 0 and sample skewness > 1,
  log(x + c) with c = half the smallest positive value; otherwise untransformed.
* Categorical results: the 0/1 `result_positive` indicator (Positive/Weakly positive/Reactive = 1, Negative = 0,
  others missing); eligible if >= 20 positives and >= 20 negatives among members. Categorical variables without a
  positive/negative coding are not used.
* Oral microbiome: the four alpha-diversity indices, plus genus relative abundance for genera present (> 0) in
  >= 20% of layer members, transformed log10(x + 1e-4).
* The same analyte released twice in different units (e.g. urine albumin in ug/mL and mg/L; serum folate in ng/mL and
  nmol/L) is kept once: a feature whose values are an exact linear copy (|r| > 0.999) of an earlier feature of the
  layer is dropped (added before any model was fitted, after inspecting the feature list).
* A layer with no eligible feature is reported and skipped. The feature list per layer is written to
  `nhanes_layers_features.csv`.

Base set **clinical** = the pre-specified `clinical` set of the existing plan (demographics + BMI, waist, BP, exam
pulse, 32 core labs, comorbidity count, medication count, smoking; `comorbidity_count` dropped for the proxy).
Wearable = the pre-specified 26 wearable features. The core lab layers already inside `clinical`
(`core_cbc`, `core_biochemistry_profile`, `lipids_glycaemia_full`, `nutrition_vitamin_d_b12`) are not tested as
increments.

**Full-sample lab layers** (defined before modelling, label-free): layers with >= 85% coverage of the analysis
population in each cycle. Besides the four core layers these are `hormones_sex_steroids`, `immune_celiac_serology`,
`infectious_hepatitis`, `nutrition_folate_b12_status`, `tobacco_cotinine`, `urine_albumin_creatinine`,
`urine_flow_rate` = **ALL_FULL** (7 layers).

## 3. Comparisons

Model: `pen_lr` only (primary model of the existing plan; `models.make_model("pen_lr")`, preprocessing inside each
fit).

**A. Per-layer screen (does layer L add over demographics + clinical?).** For every lab layer not in `clinical`
with >= 1,000 members in the analysis population and >= 1 eligible feature, and every target: clinical+L vs clinical
in (target population ∩ members of L). Repeated stratified 5-fold CV with **3 repeats** (reduced from 5 for compute;
~55 layers x 6 targets), seeds `SEED + repeat`.

**B. Headline comparisons (5-fold x 5 repeats):**
* B1 `labs_full`: clinical + ALL_FULL vs clinical, among people with all ALL_FULL layers.
* B2 `wear_over_full`: clinical + ALL_FULL + wearable vs clinical + ALL_FULL (same population). **This is the
  "does accelerometry still add on top of all full-sample lab layers" test.**
* For each subsample group g, population = people with ALL_FULL and every layer of g:
  * B3_g `labs_<g>`: clinical + ALL_FULL + g vs clinical + ALL_FULL;
  * B4_g `wear_over_<g>`: clinical + ALL_FULL + g + wearable vs clinical + ALL_FULL + g
    (**accelerometry on top of all lab layers available to that subsample**).
  Groups (fixed now): `envA` = env_metals_urine, env_arsenic_urine, env_perchlorate_nitrate_thiocyanate_urine,
  env_voc_metabolites_urine, env_pahs_urine, nutrition_iodine_urine, nutrition_trace_elements_serum (one-third
  subsample A, both cycles); `fasting` = metabolic_fasting_subsample; `fatty_acids` = nutrition_fatty_acids_serum;
  `blood_vocs` = env_vocs_blood; `oral_microbiome` = oral_microbiome_16s (2011-2012 only); `klotho` =
  aging_klotho_surplus (ages 40-79); `maximal` = envA + fasting + fatty_acids + blood_vocs (the largest set of
  subsample layers that co-occur; ~1,000 people, so several targets will be not estimable).

## 4. Metrics and inference (identical machinery to the existing plan)

OOF predictions per repeat; AUROC, AUPRC, Brier (mean over repeats); participant-level bootstrap (B = 1000
multinomial resamples; the same resamples for both sets of a comparison) for 95% CIs and the paired Delta
(AUROC primary; AUPRC and Brier alongside); two-sided bootstrap p = 2 x min(share <= 0, share >= 0).

**Multiplicity.** Per-layer screen (A): Benjamini-Hochberg FDR across all estimable layer x target Delta-AUROC tests.
Headline (B): Holm across the 6 targets within each comparison.

## 5. Controls

* **Block-shuffle control** (increment null): the added block (layer L, group g, or the wearable block) is permuted
  across participants (rows kept together, person linkage broken), repeat-1 folds, 20 shuffles; statistic = Delta
  AUROC vs the base on the same folds; empirical p = (1 + #{shuffled >= observed}) / 21. Run for **every B
  comparison**, and for every A comparison whose primary 95% CI excludes 0 (the control gates positives).
* **Label-permutation null** (Test-6 style): the added block alone (layer-only or wearable-only model, `pen_lr`, one
  5-fold CV) under 50 label permutations; p = (1 + #{null >= observed}) / 51, observed = repeat-1 AUROC of the same
  block-only model. Run for the same set as the shuffle control.
* **Positive control:** mortality, as before.

## 6. Decision rules

* A layer increment (A) is **evidence-supported** only if: 95% CI of Delta AUROC excludes 0, BH q < 0.05, and the
  block-shuffle p <= 0.05 (observed beats all 20 shuffles). CI > 0 but failing q or shuffle = "not robust".
  Otherwise null.
* "Accelerometry still adds on top of all lab layers" (B2/B4) is **supported** for a target and population if the
  95% CI of Delta AUROC excludes 0 and the wearable-block shuffle p <= 0.05; Holm-adjusted p reported alongside.
* Positive and null results are reported with equal prominence; n and cases shown for every row.

## 7. Outputs

`nhanes_layers_features.csv`, `nhanes_layers_populations.csv`, `nhanes_layers_screen.csv` (A),
`nhanes_layers_headline.csv` (B), `nhanes_layers_performance.csv`, `nhanes_layers_controls.csv`, a draft forest
figure `results/figures/drafts/nhanes_layers_increments.png`, and `results/NHANES_LAYERS_RESULTS.md`.

## 8. Known limitations (stated in advance)

Cross-sectional (except mortality); self-reported targets; subsample populations differ from the full population
(and from each other), so increments are not comparable across groups; no survey weights in prediction metrics;
many layers x targets inflate chance positives (hence FDR + shuffle gates); below-LOD substitution values; the oral
microbiome and several surplus layers exist in one cycle only; 3-repeat CV for the per-layer screen is less stable
than 5 repeats.

## 9. Deviations from this plan

_(appended during/after the analysis; nothing above this section was edited after the analysis started)_

1. **Pipeline fixes during the metrics/controls stages (no estimate changed).** (a) The controls stage first
   crashed merging an empty control table (no controls had run yet); the empty table now has its columns. (b) The
   verdict string "NULL" was read back by pandas as a missing value; it is now written "NULL (CI includes 0)".
   (c) Comparisons that were not fitted (< 30 cases) were dropped from the headline table because their metric rows
   are empty; they are now taken from the population table and listed as not estimable. (d) Bootstrap results are
   cached (`data/interim/nhanes_layers/eval_rows.csv`) so the metrics stage is not recomputed; seeds unchanged.
2. **Layers screened.** 50 layers entered the per-layer screen (the plan said "~55"). Not screened, by the rules
   of §2-3: `immune_autoantibodies_surplus` (every indirect-immunofluorescence intensity score has < 10 distinct
   values and the 0/1 pattern indicators are coded as ranges, so no feature met the continuous or binary rule),
   `infectious_sti_urine` (chlamydia/trichomonas positives < 20 or coverage < 70%), `env_volatile_nitrosamines_urine`
   (every analyte detected in < 50%), `env_fluoride_plasma` (212 members < 1,000).
3. **Descriptive additions to the report** (not in the plan; no decision rule changed): verdict counts; a side-by-side
   of the 2011-2014 Test 2 increments (accelerometry over clinical alone) and the `wear_over_full` increments; the
   list of layers not screened.
