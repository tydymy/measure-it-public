# Analysis plan — MUSCLE-ME daily steps vs Long COVID / ME/CFS (charlton_lc_mecfs_cpet_source)

**Locked:** 2026-09-24, before any case-control comparison was computed in this module. The file is written by
`measure_it.wearables.muscle_me_steps_analysis.lock_plan()`; the analysis refuses to run if this file is missing or
has been edited (sha256 check against the plan text in the module).

**Question:** can a wearable (accelerometer daily step count) separate people with Long COVID or ME/CFS from matched
healthy controls, and how does it compare with the laboratory comparator (maximal cycle CPET, `VO2_rel`)?

## What had been seen before this plan was written

* The primary analysis, the decision rule and secondaries (a)-(c) were locked on 2026-09-24 in
  `docs/DEVICE_DATASET_DISCOVERY.md` ("Pre-specified primary analyses ... locked 2026-09-24, before any outcome was
  examined"). They are restated here unchanged in substance.
* In this module, before this file was written, only the **structure** of the Source Data workbook was read: sheet
  names; the 61 column names of sheet `Main` (`Subject, Session, Group, Sex, Sx_duration, Steps`, 20 CPET columns,
  muscle fibre-type / fibre-size / capillary / mitochondrial-respiration columns, 3 lactate columns); the
  `Group x Session` table (AGBRESA BDC 24 / HDT55 24; MUSCLE-ME CON 30, LC 25, ME 26); the `Session x Sex` table
  (CON 15F/15M, LC 13F/12M, ME 13F/13M); and non-missing counts per column per session (Steps: CON 22/30,
  LC 25/25, ME 26/26; VO2_rel: CON 30/30, LC 23/25, ME 26/26). Each `Subject` appears in one `Session` only, so
  there are no repeated measures among the MUSCLE-ME rows.
* **No step, VO2 or other outcome value was printed or summarised by group.** The only outcome number already known
  is the patient step range quoted in the paper's Results (733-8,609 steps/day), which was in the discovery table.
* **Age is not in the workbook.** The demographics-only baseline can therefore use sex only.

## Data

* File `41467_2026_75725_MOESM4_ESM.xlsx` (Charlton et al., Nat Commun 2026, doi:10.1038/s41467-026-75725-y),
  sheet `Main`, rows `Group == "MUSCLE-ME"` (81 people). The AGBRESA bed-rest rows (24 people x 2 sessions) are
  excluded from every analysis here (separate cohort, no step data, not an invisible illness).
* Labels (`Session`): `LC` = Long COVID meeting the Canadian Consensus Criteria for ME/CFS; `ME` = ME/CFS diagnosed
  before 2019; `CON` = age- and sex-matched healthy controls (per the paper). Patients = LC + ME.
* Participant id: `charlton_lc_mecfs_cpet_source:<Subject>`.

## Primary analysis (one test)

| item | specification |
|---|---|
| population | MUSCLE-ME rows with non-missing `Steps` (complete case): patients (LC or ME) vs CON |
| label | 1 = patient (LC or ME), 0 = CON |
| feature | `Steps` (daily step count as released), one feature |
| model | none. Direction fixed in advance: lower steps = patient. Nothing is fitted. |
| metric | AUROC = P(control Steps > patient Steps), ties counted 0.5 |
| uncertainty | 95% percentile CI from 2,000 stratified person bootstraps (patients and controls resampled separately, with replacement), seed 20260923 |
| decision | **"wearable captures the illness-associated functional limitation"** if the CI lower bound is >= 0.70; **"null"** if the CI includes 0.50; **"weak"** otherwise |
| missing controls | the 8 controls without `Steps` are also analysed as bounds: worst case = assign them the lowest observed `Steps` value among the MUSCLE-ME rows; best case = the highest observed value; same AUROC and bootstrap. The decision is taken on the complete-case result; the bounds are reported next to it. |

Also reported for the primary (descriptions of the same test, not additional tests): AUPRC with the same bootstrap
(no-skill baseline = patient prevalence); the two-sided Mann-Whitney U p; and a label-permutation null (10,000
permutations of the patient/control labels among the complete cases, seed 20260923,
p = (1 + #{perm AUROC >= observed}) / (1 + 10,000)).

## Pre-specified secondary analyses (all reported, whatever they show)

* **(a)** LC vs CON and ME/CFS vs CON separately, `Steps`, same method (fixed direction, bootstrap CI; the decision
  rule is applied descriptively).
* **(b)** Laboratory comparator: `VO2_rel` (relative VO2peak), patients vs CON, same method, direction fixed
  (lower = patient).
* **(c)** Specificity check: LC vs ME/CFS by `Steps`, reported as P(LC Steps > ME/CFS Steps) with a bootstrap CI.
  Expected near 0.5; a result near 0.5 means steps do not tell the two illnesses apart.
* **(d)** Head-to-head device vs lab: among people with both `Steps` and `VO2_rel`, AUROC(Steps) minus
  AUROC(VO2_rel), both with the fixed direction, paired stratified person bootstrap (2,000, seed 20260923).
* **(e)** Demographics-only baseline, cross-validated. Features: sex only (age is not released). Models: L2 logistic
  regression (C = 1.0, standardised inside each fold) on (1) sex, (2) log10(Steps), (3) sex + log10(Steps), in the
  73 complete cases. Stratified group 5-fold CV (groups = participant; each person has one row, so this is
  person-level CV), repeated 20 times (seed 20260923 + repeat). AUROC and AUPRC on out-of-fold predictions (mean and
  2.5/97.5 percentiles across repeats), and **Delta AUROC = model 3 minus model 1** with a paired stratified person
  bootstrap CI (2,000) on the repeat-averaged out-of-fold predictions. Sensitivity and specificity of model 2 at a
  Youden threshold chosen **inside each training fold** and applied to the held-out fold (pooled over repeats).
  Label-permutation null for model 2: 1,000 permutations, one CV repeat each.
* **(f)** Calibration of model 3 (n = 73, crude): Brier score, calibration intercept and slope (logistic
  recalibration of repeat-averaged out-of-fold logits), and a 4-bin reliability table. Descriptive only.
* **(g)** Exploratory and descriptive only: every numeric `Main` column (CPET, muscle histology, mitochondrial
  respiration, capillarisation, lactate; `Sx_duration` excluded because controls have none) compared patients vs CON:
  Hedges g (patient minus control, pooled SD) with a 2,000-bootstrap CI, AUROC = P(patient > control) with **no
  direction chosen**, Mann-Whitney p, Benjamini-Hochberg q across all columns. **No claim about any device or
  measurement is made from (g).** It maps which person-level layers in this file differ, for the project lead's
  question about layers beyond condition-level enrichment. Muscle-biopsy variables are physiology/histology, not
  omics.

No other feature, threshold, model, subgroup or transformation will be tried. Anything added after this date is
labelled post hoc in the results.

## Interpretation fixed in advance

Activity limitation is part of the case definition for both groups (incorporation bias). A high AUROC would show
that an accelerometer measures the **defining** limitation in **diagnosed** patients versus **healthy** people. It
would **not** show detection of undiagnosed illness, and it would not show specificity against deconditioning,
depression or other chronic disease (the comparison group is healthy). Patients able to complete maximal CPET are the
milder end of the spectrum. The cohort probably overlaps Appelman et al. 2024 (same trial NCT05225688); the two are
not counted as independent evidence. The paper's own statistics (published claims) are kept separate from this
module's computation.

## Outputs

* `data/processed/participants__charlton_lc_mecfs_cpet_source`,
  `participant_wearable_features__charlton_lc_mecfs_cpet_source` (steps),
  `participant_exam_features__charlton_lc_mecfs_cpet_source` (CPET, muscle, lactate),
  `phenotype_signatures__charlton_lc_mecfs_cpet_source`.
* `results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md`, `results/tables/charlton_lc_mecfs_cpet_source_*.csv`,
  `results/figures/drafts/charlton_lc_mecfs_cpet_source_*`.
* Modules: `measure_it.ingestion.muscle_me_charlton` (ingestion), `measure_it.wearables.muscle_me_steps_analysis`
  (analysis).
