# Analysis plan: serum arginase-1 in endometriosis vs surgical and healthy controls (RepOD PY1P9X)

**Dataset id:** `endo_arg1_repod`. It is rank 2 in `docs/LAB_DATASET_DISCOVERY.md`, tied within 0.01 with rank 3, so
both tied datasets are analysed.

**Data:** RepOD doi:10.18150/PY1P9X, the per-person file behind Pliszkiewicz et al., J Clin Med 2024;13:1489
(PMC10933979).

**Written 2026-09-24, before any arginase value was compared between groups.** Module:
`src/measure_it/labs/endo_arg1_repod.py`. Nothing here may change after results are seen.

## What was seen before this plan

* **Schema only:**
  * 217 rows x 50 columns.
  * `Group` counts: Patient 127, Ctrl 25, Ctrl2 9, Ctrl3 54, blank 2.
  * Cross-tabs of `Group` against the **non-biomarker** columns `diagnosis`, `Grading`, `platenr` and `ELISAkitnr`.
  * Non-missing counts by group for `Arg1PREng/ml`, `Arg1PREng/mlold`, age, glycaemia and creatinine.
* **What those show:**
  * Ctrl = 25 women operated on for other benign gynaecological conditions (myoma, cysts, pain, diagnostic
    laparoscopy and others).
  * Ctrl2 = 9 women with uterine myomas.
  * Ctrl3 = 51 women recorded as 'healthy', plus 3 others.
  * **None of the 54 Ctrl3 samples has an ELISA plate or kit number.** Patients and surgical controls share plates
    0-3 and 6.
  * Ctrl2 has no glycaemia or creatinine values.
* **Published claims** (paper text, kept separate from our computation):
  * Groups: 105 endometriosis (GB), 22 surgical controls (K1), 53 healthy (K2).
  * Preoperative Arg-1 was higher in GB than in either control group.
  * ROC AUC 0.848 (95% CI 0.769-0.926) for GB vs K1, and 0.912 (0.859-0.965) for GB vs K2.
  * The Arg-2 assay was not reproducible and was excluded.
  * Arginase activity had AUC < 0.7.
* The released file has more people than the paper (127 / 34 / 54 against 105 / 22 / 53), so our numbers are not
  expected to match exactly.

## Primary analysis (one)

* **Question:** does preoperative serum arginase-1 separate endometriosis from women operated on for other benign
  gynaecological conditions? This is the clinically relevant comparison.
* **Sample:**
  * Cases (label 1): `Group == "Patient"`, surgically confirmed endometriosis.
  * Controls (label 0): `Group in {"Ctrl", "Ctrl2"}`, surgical non-endometriosis controls.
  * Complete cases on `Arg1PREng/ml` (the column without the `old` suffix), with `n/a` treated as missing.
* **Statistic:** there is one pre-specified analyte, and the direction is fixed in advance from the paper (higher =
  endometriosis). Nothing is fitted, so no cross-validation is needed.
  * AUROC = P(ARG1 in a case > ARG1 in a control), with ties counted as 0.5.
  * 95% CI from 2,000 group-stratified person bootstraps (seed 20260923).
  * One-sided label-permutation p from 10,000 shuffles.
* **Decision rule:** "ARG1 separates endometriosis from surgical controls" if the CI lower bound is above 0.5
  **and** the permutation p is below 0.05. The published AUC (0.848) is placed next to ours, not merged with it.

## Pre-specified secondary analyses (all reported whatever they show)

* **S1 batch check (key):** a plate-stratified AUROC for the primary contrast.
  * The AUROC is computed within each ELISA plate that holds both groups.
  * The plate AUROCs are combined as a mean weighted by the number of case-control pairs (van Elteren-style).
  * Its CI comes from a stratified bootstrap that resamples within plate x group.
  * If the plate-stratified CI lower bound is at or below 0.5 while the pooled one is above 0.5, the pooled result
    is reported as batch-sensitive.
* **S2 healthy controls:** Patient vs Ctrl3, with the same fit-free AUROC, CI and permutation p. Labelled
  **batch-unknown** whatever it shows, because no plate is recorded for any Ctrl3 sample.
* **S3 adjustment for age:** CV Delta AUROC (age + log ARG1, minus age alone) for the primary contrast.
  * L2 logistic, C = 1.
  * Stratified 5-fold x 20 repeats.
  * Paired person bootstrap CI (2,000 resamples).
* **S4 per-analyte panel:** Patient vs surgical controls for seven analytes:
  * ARG1 (pre); arginase activity (pre); ARG2 (pre; flagged unreliable by the authors); glycaemia; creatinine; AST;
    ALT.
  * For each: fit-free AUROC, Hedges' g on log values, two-sided Mann-Whitney p, and BH q across the 7.
* **S5 stage:** Spearman correlation of ARG1 (pre) with rASRM stage (I-IV) among patients.
* **S6 sensitivity at 90% specificity**, with the threshold at the 90th percentile of control ARG1.

## Outputs

* Processed tables: `participants__endo_arg1_repod`, `participant_labs__endo_arg1_repod`,
  `phenotype_signatures__endo_arg1_repod`.
* Result tables: `results/tables/endo_arg1_repod_*`.
* Report: `results/ENDO_ARG1_REPOD_RESULTS.md`.
