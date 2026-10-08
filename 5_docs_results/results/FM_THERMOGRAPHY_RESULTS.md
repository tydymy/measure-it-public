# Infrared thermography in fibromyalgia: pre-specified case-control analysis

Dataset `fm_thermography`: PLOS ONE 2021 0253281 S1 Data (Sempere-Rubio et al., University of Valencia), 86 women with
fibromyalgia (ACR 2010, rheumatologist-assessed) and 92 age-matched women without FM symptoms, all post-menopausal.
One resting infrared image session per person (FLIR E60bx, 24 degC room, 15 min acclimatisation); region-of-interest
min/max/mean skin temperature at six sites. Plan: `docs/ANALYSIS_PLAN_FM_THERMOGRAPHY.md` (saved 2026-09-24 17:12 UTC, before any group statistic
was computed; the primary analysis was locked earlier in `docs/DEVICE_DATASET_DISCOVERY.md`). The ordering is
verified from file timestamps only: the plan predates every result file, and the discovery doc predates this
download; what was viewed interactively before the plan cannot be verified from files. Analysis
`fm_thermography_analysis 2026-09-24.1`, computed 2026-09-25T22:05:42+00:00.

## Answer first

**The pre-specified primary test is not met.** Adding the six regional average skin temperatures to age + BMI
changes cross-validated AUROC by **0.053 (95% CI -0.028 to 0.129)**. The
interval includes 0 (9.2% of bootstrap resamples are at or below 0). By the
locked rule this is reported as **a replication of the authors' published null** ("thermal imaging ruled out as a
supplementary assessment").

The two criteria of the rule split, and both results are stated:

| criterion (pre-specified) | result | met? |
|---|---|---|
| Delta AUROC (age + BMI + 6 temperatures minus age + BMI), CI lower bound > 0 | 0.053 (-0.028 to 0.129) | **no** |
| device-only AUROC beats 1,000 label shuffles, p < 0.05 | AUROC 0.630 (0.553 to 0.708); null mean 0.496, null 95th percentile 0.582; p = 0.0050 (4 of 1000 shuffles reached the observed AUROC) | yes |

What that means in plain terms:
* There **is** a modest group-level difference: the camera alone ranks a random FM woman above a random control
  about 63% of the time, more than label shuffling produces. The difference runs the
  **opposite way to the authors' hypothesis** (cooler skin from autonomic impairment): women with FM have *warmer*
  lower back and knees; the other four regions do not differ unadjusted. The paper itself reports this direction
  (Discussion: lower back and knee "higher" in FM), so this is a reproduction of their finding, not a new one.
* It is **not shown to add** to what age and BMI already give, and it is **far from useful for finding an individual
  case**: at 90% specificity the device-only model flags **25%
  (16-36%)** of women with FM.
* The primary result is **inconclusive rather than a demonstrated absence**: the Delta interval runs from a small
  loss to a gain of 0.13, so a useful increment is not excluded; 177 women cannot resolve it. The
  locked rule calls this a null, and it is reported as one.
* The design compares FM with healthy women only. It cannot say whether warm knees mark FM or any chronic pain,
  lower activity, or the recruitment difference between the groups.

**Verdict for the lead's question ("can a device other than capillaroscopy be shown, with public data, to uncover an
invisible illness?"): this dataset does not show it for thermography.** It shows a small group difference
that we reproduce from the same released rows (18/18 per-region significance calls match the paper; this is a
reproduction of the authors' analysis, not an independent replication) that falls short of the pre-registered bar
for added value, with low screening sensitivity.

## Published claim vs our own demonstration

| | authors (published, PLOS ONE 2021) | this analysis (our own, from the released rows) |
|---|---|---|
| Question | group difference in skin temperature (Mann-Whitney per region) | discrimination beyond age + BMI (cross-validated), plus per-region effects |
| Direction | FM warmer at lower back and knee (Discussion); hypothesis had been cooler | same direction (FM warmer) |
| Per-region differences | lower back and knee min/max/mean p < 0.05; rest n.s. | same 6/18 regions at p < 0.05 (**18/18 concordant**); after BH over 18, only knee mean survives (q = 0.010) |
| Size | effect r 0.15-0.26; "did not exceed 0.5 degC", called not clinically meaningful | knee mean +0.59 degC, g = 0.58; lower back mean +0.48 degC, g = 0.37 |
| Conclusion | "not an effective supplementary assessment tool" | incremental value not demonstrated (pre-specified rule unmet); device-only signal above chance but weak |

## Primary and secondary model results (S3, S4, S5)

177 women (85 FM, 92 controls), L2 logistic (C = 1), stratified 5-fold x 20 repeats; 95% CIs from a
group-stratified person bootstrap (2,000). AUPRC chance = 0.480.

| feature set | role | AUROC (95% CI) | AUPRC (95% CI) | Brier | calibration slope | sensitivity at 90% specificity |
|---|---|---|---|---|---|---|
| demo | primary (model 1, baseline) | 0.582 (0.497 to 0.663) | 0.581 (0.509 to 0.671) | 0.245 (no-skill 0.250) | 0.72 | 0.191 (0.103 to 0.356) |
| device | primary (model 2, device only) | 0.630 (0.553 to 0.708) | 0.640 (0.572 to 0.713) | 0.237 (no-skill 0.250) | 0.69 | 0.248 (0.158 to 0.359) |
| demo_device | primary (model 3) | 0.635 (0.557 to 0.709) | 0.646 (0.573 to 0.724) | 0.239 (no-skill 0.250) | 0.63 | 0.286 (0.197 to 0.388) |
| demo_device18 | secondary S5 (all 18 temperatures) | 0.683 (0.605 to 0.756) | 0.677 (0.604 to 0.757) | 0.229 (no-skill 0.250) | 0.64 | 0.296 (0.187 to 0.426) |

* **Primary Delta** (model 3 minus model 1): 0.053 (-0.028 to 0.129); across
  the 20 repeats it ranged 0.005 to 0.092.
* **S5, all 18 temperatures** (min, max and mean) instead of the six means: Delta 0.100
  (0.011 to 0.185). This interval excludes 0. It is the one pre-specified feature-set
  sensitivity analysis; it cannot replace the primary, and a secondary that clears the bar where the primary does
  not is weak evidence (two feature sets were examined, and the 18 columns are three summaries of the same six
  images). It is reported because the plan said it would be, whatever it showed.
* **Calibration (S3):** calibration-in-the-large is ~0 for every model, but slopes of
  0.63-0.72 (< 1) mean the predicted probabilities are too
  spread out for the weak signal: the models over-state how sure they are. The decile reliability plot
  (`results/figures/drafts/fm_thermography_calibration.png`) is noisy at n = 177. No predicted probability from these
  models should be read as an individual's risk.
* Age + BMI alone: AUROC 0.582 (0.497 to 0.663). That comes from BMI (FM mean 27.8 vs control
  26.0 kg/m2, g = 0.42); the groups were age-matched (g = 0.00).

## Per-region effects (S1): unadjusted, all available rows

FM minus control. BH q across the 18 temperature columns.

| region_stat | n FM / control | mean FM / control (degC) | difference (degC) | Hedges g (95% CI) | Mann-Whitney p | BH q |
|---|---|---|---|---|---|---|
| `neck_min` | 86 / 92 | 31.58 / 31.89 | -0.31 | -0.26 (-0.56 to 0.03) | 0.1152 | 0.259 |
| `neck_max` | 86 / 92 | 34.54 / 34.43 | 0.10 | 0.11 (-0.18 to 0.40) | 0.4132 | 0.465 |
| `neck_ave` | 86 / 92 | 33.56 / 33.57 | -0.01 | -0.01 (-0.32 to 0.27) | 0.8808 | 0.881 |
| `upper_back_min` | 86 / 92 | 30.54 / 30.40 | 0.14 | 0.11 (-0.19 to 0.40) | 0.3969 | 0.465 |
| `upper_back_max` | 86 / 92 | 33.75 / 33.63 | 0.12 | 0.11 (-0.18 to 0.40) | 0.3793 | 0.465 |
| `upper_back_ave` | 86 / 92 | 32.16 / 32.12 | 0.04 | 0.04 (-0.26 to 0.33) | 0.7754 | 0.821 |
| `lower_back_min` | 86 / 92 | 29.75 / 29.25 | 0.50 | 0.35 (0.05 to 0.65) | 0.0300 | 0.108 |
| `lower_back_max` | 86 / 92 | 33.99 / 33.65 | 0.35 | 0.31 (0.01 to 0.62) | 0.0460 | 0.138 |
| `lower_back_ave` | 86 / 92 | 31.69 / 31.21 | 0.48 | 0.37 (0.08 to 0.67) | 0.0207 | 0.095 |
| `chest_min` | 86 / 92 | 30.83 / 30.64 | 0.18 | 0.15 (-0.16 to 0.44) | 0.2820 | 0.423 |
| `chest_max` | 86 / 92 | 34.28 / 34.10 | 0.18 | 0.21 (-0.08 to 0.50) | 0.2389 | 0.423 |
| `chest_ave` | 86 / 92 | 32.34 / 32.19 | 0.15 | 0.14 (-0.17 to 0.41) | 0.3418 | 0.465 |
| `knee_min` | 85 / 92 | 29.70 / 29.24 | 0.46 | 0.43 (0.14 to 0.72) | 0.0105 | 0.095 |
| `knee_max` | 85 / 92 | 31.86 / 31.49 | 0.36 | 0.37 (0.10 to 0.67) | 0.0211 | 0.095 |
| `knee_ave` | 85 / 92 | 30.80 / 30.21 | 0.59 | 0.58 (0.31 to 0.88) | 0.0006 | 0.010 |
| `elbow_min` | 85 / 92 | 30.28 / 30.05 | 0.23 | 0.23 (-0.07 to 0.53) | 0.1047 | 0.259 |
| `elbow_max` | 85 / 92 | 32.79 / 32.61 | 0.18 | 0.22 (-0.06 to 0.51) | 0.2688 | 0.423 |
| `elbow_ave` | 85 / 92 | 31.52 / 31.34 | 0.18 | 0.21 (-0.08 to 0.51) | 0.1968 | 0.394 |

## Age- and BMI-adjusted differences (S2), six means, primary sample

OLS `temperature ~ FM + age + BMI`; standardised by the pooled within-group SD; BH q across the six.

| region | adjusted difference, degC (95% CI) | adjusted standardised difference (95% CI) | p | BH q | BMI coefficient (degC per kg/m2) |
|---|---|---|---|---|---|
| `neck_ave` | 0.14 (-0.15 to 0.43) | 0.14 (-0.15 to 0.42) | 0.3439 | 0.344 | -0.084 |
| `upper_back_ave` | 0.26 (-0.05 to 0.57) | 0.23 (-0.04 to 0.50) | 0.0943 | 0.113 | -0.120 |
| `lower_back_ave` | 0.54 (0.16 to 0.93) | 0.43 (0.12 to 0.73) | 0.0060 | 0.018 | -0.035 |
| `chest_ave` | 0.30 (0.00 to 0.61) | 0.29 (0.00 to 0.57) | 0.0483 | 0.072 | -0.083 |
| `knee_ave` | 0.64 (0.33 to 0.94) | 0.64 (0.33 to 0.94) | 0.0001 | 0.000 | -0.027 |
| `elbow_ave` | 0.31 (0.07 to 0.54) | 0.36 (0.08 to 0.65) | 0.0127 | 0.025 | -0.069 |

Adjustment **enlarges** the FM differences rather than removing them. Higher BMI goes with cooler skin at most
sites, and the FM group has the higher BMI, so the unadjusted comparison understates how much warmer the FM group
is. After adjustment 3 of 6 regions have q < 0.05 (lower_back, knee, elbow). This is
a group-level description. It does not change the primary result, which already included BMI.

## Checks against the paper

* **Numbers reproduce.** All 18 group means match paper Table 2 within 0.016
  degC. The p < 0.05 calls agree for 18/18 regions (`results/tables/fm_thermography_paper_concordance.csv`).
* **The paper's "under 0.5 degC" statement does not hold for the knee mean.** The Discussion says the knee and
  lower-back differences "did not exceed 0.5°, either in the average or maximum temperature". Paper Table 2's
  own knee means (30.80 vs 30.21) differ by 0.59 degC, and so do the released rows (+0.59 degC).
  The lower-back mean difference (+0.48 degC) and all maxima are under 0.5 degC.
* **Group coding verified two ways.** All 23 sociodemographic category counts in paper Table 1 are reproduced
  exactly by the file's codes under Group 1 = FM, and so are the Table 2 FM-column means. The codes'
  labels were recovered this way; they were **not** used in any analysis (not in the plan).
* **Paper text error found.** The Results text says BMI was "27.83 (4.75)" for controls and "25.97 (4.00)" for FM.
  In the file, **FM has the higher BMI (27.83)**, and the file's group coding is confirmed by both tables. So the
  text has the two BMI values swapped. That matters here because BMI works against the temperature difference
  (see S2).
* **Education differs sharply by group:** university education 12/86 FM vs 36/92 controls (paper Table 1,
  chi-square p < 0.05). This is consistent with the different recruitment channels. It was not adjusted for, by
  design; it is one reason the group difference may not be specific to FM.

## Independent review checks (2026-09-24)

* Recomputed from the raw xlsx with separate code (sklearn `roc_auc_score`, statsmodels BH): device-only AUROC
  0.630, Delta 0.053, knee-mean difference +0.589 degC (g 0.584, BH q 0.0105), adjusted knee 0.636 SD, sensitivity
  at 90% specificity 0.248, BMI 27.83 (FM) vs 25.97 (controls), university 12/86 vs 36/92: all match the tables.
* **Cross-validation partition sensitivity.** Folds are drawn on rows sorted by participant ID. Run once on the
  file's original row order (same seeds), the numbers move by Monte Carlo CV noise: device-only AUROC 0.636, Delta
  0.060 (bootstrap CI with an independent seed -0.026 to 0.142). The decision is unchanged. The reported values
  are the pre-specified run; this check is reported only to show the size of partition noise (about 0.005-0.007).
* **The bootstrap CIs hold the out-of-fold predictions fixed** (as pre-specified), so they include person sampling
  but not model re-fitting variability. They are, if anything, too narrow. That cannot turn the primary null into
  a positive, but it does weaken S5: its lower bound (0.011) is close enough to 0 that a refit-aware interval
  could include 0.

## Caveats carried with every number

Single-site case-control sample of post-menopausal women (Valencia); FM recruited from associations/specialised units, controls by advertisement, so selection may contribute; no non-FM pain comparator; room held at 24 degC by protocol but per-session room temperature, date and season are not recorded. The authors' own published conclusion is negative. All participants are post-menopausal women aged 38-70 from one region of Spain. Only ROI summaries were
released, not images, so no other thermal feature (asymmetry, hot spots, texture) can be tested from this source.

## What this adds to the project's device evidence

| device, illness | label | our result | status |
|---|---|---|---|
| NHANES wrist accelerometry, general-population fatigue/limitation | self-report proxy, not a diagnosis | modest (Delta AUROC +0.02-0.03) | earlier round |
| Stanford smartwatch resting HR, Long COVID (Uwakwe 2025) | self-reported Long COVID | null (AUROC ~0.50) | earlier round |
| **Infrared thermography, fibromyalgia (this report)** | **ACR 2010, rheumatologist-assessed** | **primary rule unmet (Delta 0.053, CI -0.028 to 0.129); device-only AUROC 0.630, permutation p 0.005; sensitivity 25% at 90% specificity** | **replication of the published null by the locked rule** |

Measurement mapping: no class in `configs/measurements.yaml` is an exact fit. The rows carry
`measurement_class_id = continuous_temperature` (a signal match: skin temperature; that class is a wearable,
this is a clinic camera). A new class `infrared_thermography` (modality_family `imaging`) is proposed to the config
owner, not added. Condition: `normalize_condition("fibromyalgia")` gives `fibromyalgia` (ICD-10-CM M79.7).

## Outputs

* Processed: `participants__fm_thermography` (178), `participant_device_features__fm_thermography` (178),
  `phenotype_signatures__fm_thermography` (26 rows; object ids
  `signature:fm_thermography|fibromyalgia|<feature>`).
* Tables: `results/tables/fm_thermography_{model_performance, delta_auroc, primary_decision, permutation_null,
  feature_effects, adjusted_effects, full_data_coefficients, calibration_bins, oof_predictions, paper_concordance}.csv`,
  `fm_thermography_run_metadata.json`.
* Figures (drafts): `results/figures/drafts/fm_thermography_auroc_and_null.png`, `results/figures/drafts/fm_thermography_feature_effects.png`, `results/figures/drafts/fm_thermography_calibration.png`.

## Reproduce

```bash
uv run python -m measure_it.ingestion.fm_thermography
uv run python -m measure_it.wearables.fm_thermography_analysis --jobs 16   # ~1 min; deterministic (seed 20260923)
uv run pytest tests/test_fm_thermography.py
```
