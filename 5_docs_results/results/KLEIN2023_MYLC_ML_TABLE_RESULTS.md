# Serum cortisol and plasma mediators in Long COVID (MY-LC): pre-specified re-analysis

Dataset `klein2023_mylc_ml_table`: Klein et al., Nature 2023, Supplementary Table 3, one row per participant. Analysis set
(`x0_Censor_Complete == 0`): 99 LC, 40 HC, 39 CC. Plan: `docs/ANALYSIS_PLAN_KLEIN2023_MYLC_ML_TABLE.md` (written before any analyte was compared). Analysis
`klein2023_mylc_ml_table 2026-09-24.1`, computed 2026-09-25T22:07:05+00:00.

## Answer first

**Primary: adding serum cortisol to draw time, age, sex and BMI changes cross-validated AUROC by
0.253 (95% CI 0.177 to 0.328); cortisol alone: AUROC
0.963 (0.933 to 0.987), label-permutation
p <= 0.0010 (no shuffle reached the observed value) (1000 shuffles, null 95th percentile 0.556).** Sample:
98 LC vs 77 controls with complete cortisol and covariates.
Decision by the locked rule: **cortisol adds information beyond draw time and demographics (pre-specified rule met)**.

* Fit-free cortisol AUROC (lower = LC) on the analysis set (S6): 0.966 (0.937 to
  0.988). The paper reports 0.96 (0.93-0.99) on its Gale-Shapley-matched subset.
* Excluding steroid / pituitary-adrenal flags (S4): Delta 0.266 (0.191 to 0.344).
* 61 of 144 plasma analytes differ at BH q < 0.05 (S1).

## How to read this

* The released rows **reproduce the published cortisol result**: fit-free AUROC 0.966 and CV AUROC
  0.963 against the paper's 0.96. The primary increment over draw time and demographics is
  large, because cortisol alone already carries almost all the separation (model 3 0.958 vs model 2
  0.963).
* **The size of the effect is itself a caveat.** Cortisol in LC vs controls has Hedges g -2.78 (on the
  released z-scores). That is a very large separation for a single draw of a diurnal, pulsatile hormone. The table
  releases **batch-integrated z-scores**, not concentrations, and no assay-batch or plate column. So a group-by-batch
  or recruitment-channel difference can be **neither checked nor excluded** from the public file.
* The covariates alone separate the groups (AUROC 0.704), so draw time and demographics differ by
  group. **Post hoc, exploratory (E1, not part of the decision):** within draw-time tertiles the fit-free cortisol
  AUROC is 0.969 (0.942 to 0.989). Group counts per tertile are in
  `results/tables/klein2023_mylc_ml_table_drawtime_tertile_composition.csv`.
* This is the same released data re-analysed, not an independent replication. Whether hypocortisolaemia marks Long
  COVID outside this recruitment design (clinic patients vs advertised volunteers) needs other cohorts. The open
  candidates for that check are listed in `docs/LAB_DATASET_DISCOVERY.md`: Appelman 2024, with cortisol and minutes
  since waking; and a PLOS ONE cohort with 10 vs 7 people.

## Published claim vs our computation

| | authors (published) | this analysis (same released rows) |
|---|---|---|
| Cortisol | lower in LC, similar in HC and CC; "cortisol alone achieved an AUC of 0.96 (95% CI 0.93-0.99)" in a matched subset; significant after adjusting for demographics and sample-collection time | fit-free AUROC 0.966 on the full analysis set; CV cortisol-only 0.963; increment over draw time + demographics 0.253 (0.177 to 0.328) |
| Other mediators | C4b, CCL19, CCL20, galectin-1, CCL4, APRIL, LH higher; IL-5 lower (Kruskal-Wallis across 3 groups) | 61 of 144 at q < 0.05, LC vs pooled controls (table below) |

## Models (CV, 5-fold x 20)

| model | cases / controls | AUROC (95% CI) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
| 1 covariates (draw time, age, sex, BMI) | 98 / 77 | 0.704 (0.626 to 0.780) | not run | 0.26 |
| 2 cortisol alone (released z-score) | 98 / 77 | 0.963 (0.933 to 0.987) | <= 0.0010 (no shuffle reached the observed value) | 0.91 |
| 3 covariates + cortisol | 98 / 77 | 0.958 (0.927 to 0.983) | not run | 0.90 |
| S3 144-analyte panel (released z-scores, L2 C=0.1) | 99 / 79 | 0.947 (0.903 to 0.980) | <= 0.0010 (no shuffle reached the observed value) | 0.91 |
| S4 covariates + cortisol, steroid / pituitary-adrenal flags excluded | 93 / 76 | 0.955 (0.922 to 0.982) | not run | 0.89 |

## Fit-free cortisol AUROC by control type (S2, S6)

| contrast | LC / controls | AUROC (lower cortisol = LC) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
| S6 LC vs HC+CC (primary sample) | 98 / 79 | 0.966 (0.937 to 0.988) | <= 0.0010 (no shuffle reached the observed value) | 0.91 |
| S2 LC vs HC | 98 / 40 | 0.973 (0.947 to 0.993) | <= 0.0010 (no shuffle reached the observed value) | 0.93 |
| S2 LC vs CC | 98 / 39 | 0.958 (0.921 to 0.987) | <= 0.0010 (no shuffle reached the observed value) | 0.91 |
| E1 POST HOC exploratory: LC vs HC+CC within draw-time tertiles (pair-weighted) | 98 / 77 | 0.969 (0.942 to 0.989) | <= 0.0010 (no shuffle reached the observed value) | 0.91 |

## Top 15 plasma analytes (S1, sorted by p)

| analyte | unit named in column | n LC / controls | median z LC / controls | Hedges g on z (95% CI) | Mann-Whitney p | BH q |
|---|---|---|---|---|---|---|
| Cortisol | ng/mL | 98 / 79 | -0.80 / 0.87 | -2.78 (-3.42 to -2.32) | 1.89e-26 | 0.0000 |
| IL-12 | pg/mL | 98 / 79 | -0.58 / 0.43 | -0.81 (-1.50 to -0.41) | 6.96e-11 | 0.0000 |
| LIF | pg/mL | 98 / 79 | -0.56 / 0.25 | -0.55 (-0.92 to -0.24) | 5.50e-08 | 0.0000 |
| IL-25 | pg/mL | 98 / 79 | -0.41 / 0.02 | -0.46 (-0.82 to -0.23) | 8.70e-08 | 0.0000 |
| IFNγ | pg/mL | 98 / 79 | -0.66 / 0.12 | -0.46 (-0.83 to -0.16) | 1.01e-07 | 0.0000 |
| IL-13 | pg/mL | 98 / 79 | -0.52 / 0.16 | -0.60 (-1.05 to -0.29) | 2.51e-07 | 0.0000 |
| IL-4 | pg/mL | 98 / 79 | -0.51 / -0.11 | -0.51 (-0.78 to -0.24) | 4.32e-07 | 0.0000 |
| Complement C4b | ng/mL | 98 / 79 | 0.32 / -0.51 | 0.84 (0.55 to 1.14) | 4.80e-07 | 0.0000 |
| GCSF | pg/mL | 97 / 79 | -0.48 / 0.00 | -0.47 (-0.83 to -0.15) | 7.82e-07 | 0.0000 |
| IL-21 | pg/mL | 98 / 79 | -0.40 / -0.01 | -0.36 (-0.70 to -0.06) | 8.79e-07 | 0.0000 |
| GDF15 | pg/mL | 98 / 79 | -0.03 / -0.50 | 0.55 (0.29 to 0.85) | 1.06e-06 | 0.0000 |
| IL-2 | pg/mL | 98 / 79 | -0.37 / -0.18 | -0.42 (-0.80 to -0.10) | 2.43e-06 | 0.0000 |
| FGF2 | pg/mL | 98 / 79 | -0.39 / -0.12 | -0.59 (-0.82 to -0.37) | 3.48e-06 | 0.0000 |
| CCL23 | pg/mL | 98 / 79 | 0.18 / -0.47 | 0.70 (0.42 to 1.00) | 5.20e-06 | 0.0001 |
| XCL1 | pg/mL | 94 / 76 | -0.23 / -0.38 | 0.15 (-0.17 to 0.50) | 5.72e-06 | 0.0001 |

## Caveats carried with every number

Case-control by recruitment channel (LC clinics vs advertisement); one blood draw per person; cortisol is diurnal (draw time adjusted linearly) and measured on a multiplex immunoassay; values and labels are the authors' released rows, so this re-analyses the same data, it is not an independent replication. The cohort-code mapping (1 = HC, 2 = CC, 3 = LC) and draw-time unit are inferred from the file (no codebook);
the analysed n reproduces the paper's. Bootstrap CIs hold out-of-fold predictions fixed.

## Outputs

`participants__klein2023_mylc_ml_table`, `participant_labs__klein2023_mylc_ml_table`, `phenotype_signatures__klein2023_mylc_ml_table` (146 rows);
`results/tables/klein2023_mylc_ml_table_{model_performance, delta_auroc, cortisol_fitfree, primary_decision, analyte_effects,
permutation_null}.csv`; figure results/figures/drafts/klein2023_mylc_ml_table_auroc_and_null.png.

## Reproduce

```bash
uv run python -m measure_it.labs.klein2023_mylc_ml_table --jobs 16     # deterministic, seed 20260923
uv run pytest tests/test_klein2023_mylc_ml_table.py
```
