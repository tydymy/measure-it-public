# Serum arginase-1 in endometriosis: pre-specified re-analysis of the released rows

Dataset `endo_arg1_repod` (RepOD doi:10.18150/PY1P9X; Pliszkiewicz et al., J Clin Med 2024). Plan: `docs/ANALYSIS_PLAN_ENDO_ARG1_REPOD.md` (written
before any arginase value was compared). Analysis `endo_arg1_repod 2026-09-24.1`, computed 2026-09-25T21:24:47+00:00.

## Answer first

**Primary: pre-operative serum ARG1, endometriosis (120) vs women operated on for other benign
gynaecological conditions (32): fit-free AUROC 0.846 (95% CI 0.762 to
0.918), label-permutation p <= 0.0001 (no shuffle reached the observed value) (10000 shuffles).** Decision by the locked
rule: **ARG1 separates endometriosis from surgical controls (pre-specified rule met)**. Published (paper, GB vs K1, 105 vs 22): AUC 0.848 (0.769-0.926).

* **Batch check (S1):** within ELISA plates that hold both groups, AUROC 0.854 (0.775 to
  0.928), permutation p <= 0.0001 (no shuffle reached the observed value). The pooled result is therefore
  not explained by plate: surgical controls and patients were assayed on shared plates.
* **Healthy controls (S2, batch unknown):** AUROC 0.918 (0.863 to 0.964); the
  paper reports 0.912. No ELISA plate is recorded for any healthy control, so this contrast cannot be checked for batch.
* **Age (S3, CV):** age alone AUROC 0.478; age + log ARG1 0.857; Delta
  0.378 (0.252 to 0.503). The Delta is large because age alone carries almost no
  information here.
* **This reproduces the published estimate in a larger released sample.** The file has 120 vs 32 people with the
  primary ARG1 value; the paper had 105 vs 22. We get 0.846 against the paper's 0.848, and 0.918 against 0.912 vs
  healthy.
* **Stage:** ARG1 does not track rASRM stage, which is noted here as a descriptive finding.
* **Stage (S5):** Spearman rho -0.00 (p = 0.991, n = 120).

## Published claim vs our computation

| | authors (published) | this analysis (released rows) |
|---|---|---|
| Groups | 105 GB / 22 K1 / 53 K2 | 120 endometriosis / 32 surgical controls (Ctrl + Ctrl2) with ARG1 / 53 healthy |
| ARG1 vs surgical controls | AUC 0.848 (0.769-0.926) | 0.846 (0.762 to 0.918); plate-stratified 0.854 |
| ARG1 vs healthy | AUC 0.912 (0.859-0.965) | 0.918 (0.863 to 0.964), batch unknown |
| Arginase activity | AUC < 0.7, "not clinically useful" | see panel |

## Fit-free ARG1 AUROC

| analysis | cases / controls | AUROC (95% CI) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
| PRIMARY endometriosis vs surgical controls (Ctrl + Ctrl2), pooled | 120 / 32 | 0.846 (0.762 to 0.918) | <= 0.0001 (no shuffle reached the observed value) | 0.66 |
| S1 same, ELISA-plate-stratified (plates holding both groups) | 120 / 30 | 0.854 (0.775 to 0.928) | <= 0.0001 (no shuffle reached the observed value) | 0.66 |
| S2 endometriosis vs healthy (Ctrl3; batch unknown) | 120 / 53 | 0.918 (0.863 to 0.964) | <= 0.0001 (no shuffle reached the observed value) | 0.79 |

ELISA plate composition of the primary sample (plates with both groups enter S1):

| plate | endometriosis | surgical controls |
|---|---|---|
| 0 | 4 | 4 |
| 1 | 12 | 5 |
| 2 | 18 | 2 |
| 3 | 18 | 8 |
| 4 | 20 | 0 |
| 5 | 23 | 0 |
| 6 | 11 | 11 |
| 7 | 14 | 0 |

## Panel (S4): endometriosis vs surgical controls

| analyte | n | median case / control | AUROC (higher = case) | Hedges g on log10 (95% CI) | Mann-Whitney p | BH q (7) |
|---|---|---|---|---|---|---|
| ARG1, pre-operative | 120 / 32 | 89.94 / 31.18 | 0.846 | 1.50 (1.02 to 2.10) | 1.87e-09 | 0.0000 |
| arginase activity, pre-operative | 118 / 22 | 0.80 / 0.72 | 0.602 | 0.56 (-0.02 to 1.09) | 1.30e-01 | 0.3031 |
| ARG2, pre-operative (authors: not reproducible) | 106 / 22 | 110.19 / 120.08 | 0.423 | -0.21 (-0.82 to 0.18) | 2.58e-01 | 0.4519 |
| glycaemia | 127 / 25 | 87.60 / 90.20 | 0.381 | -0.42 (-0.90 to 0.05) | 6.09e-02 | 0.2132 |
| creatinine | 127 / 25 | 0.80 / 0.80 | 0.469 | -0.10 (-0.48 to 0.27) | 6.23e-01 | 0.6282 |
| AST | 26 / 3 | 19.80 / 17.60 | 0.596 | 0.24 (-0.45 to 1.12) | 6.16e-01 | 0.6282 |
| ALT | 27 / 3 | 13.00 / 12.00 | 0.593 | 0.33 (-0.08 to 0.75) | 6.28e-01 | 0.6282 |

## Caveats carried with every number

Single-centre (Warsaw) case-control; the healthy controls (Ctrl3) have no ELISA plate recorded, so a batch difference cannot be excluded for patient-vs-healthy contrasts; the released file has more people than the paper (127/34/54 vs 105/22/53); surgical controls are few. Controls operated on for myomas, cysts or pain are the clinically closer comparator, but there are only
32. Nothing here evaluates ARG1 in women with pelvic pain referred for suspected endometriosis,
which is the intended-use population.

## Outputs

`participants__endo_arg1_repod`, `participant_labs__endo_arg1_repod`, `phenotype_signatures__endo_arg1_repod` (10
rows); `results/tables/endo_arg1_repod_*.csv`.

## Reproduce

```bash
uv run python -m measure_it.labs.endo_arg1_repod      # deterministic, seed 20260923
uv run pytest tests/test_endo_arg1_repod.py
```
