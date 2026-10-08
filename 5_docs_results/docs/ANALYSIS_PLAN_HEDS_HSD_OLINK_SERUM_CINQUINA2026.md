# Analysis plan: serum Olink proteome in hEDS/HSD vs controls (Cinquina et al., Clin Proteomics 2026)

**Dataset id:** `heds_hsd_olink_serum_cinquina2026` (rank 1 of the verified-open lab datasets in
`docs/LAB_DATASET_DISCOVERY.md`). **Written 2026-09-24, before any protein value was compared between groups.**
Module: `src/measure_it/labs/heds_hsd_olink_serum_cinquina2026.py` (ingest + analysis). Nothing in this plan may be
changed after results are seen; anything not listed here is labelled exploratory.

## What was seen before this plan

* Schema-level only:
  * sheet names and the NPX matrix layout: 460 assays x 352 samples, plus 10 Olink plate-control rows and the LOD,
    missing-frequency and normalisation rows;
  * sample-ID prefix counts: P = 88, H = 88, C = 176;
  * the column headers of the patient clinical sheets, and the sex and relationship counts in them.
* The mapping P = hEDS and H = HSD comes from the `ID NPx OLINK` column of MOESM6: sheet `hEDS` lists P1 onwards
  and sheet `HSD` lists H1 onwards.
* The paper text (Methods, Results, Abstract). Its published claims, kept separate from our computation:
  * 54 differentially expressed proteins (DEPs) in hEDS vs controls, 49 in HSD vs controls, none in hEDS vs HSD, and
    69 in the pooled patients vs controls. All used Wilcoxon + BH on the whole cohort.
  * No Italian-vs-US differences within hEDS or within HSD.
  * LASSO, then a random forest and a classification tree on the 69 DEPs. These were fitted on the full data, with
    no held-out evaluation reported. MYOC, COMP, NUDT5 and CLEC10A are among the top-ranked proteins.
  * The first rows of the paper's DEP tables (MOESM8) were visible when the sheets were listed, so the names of the
    top-ranked proteins are known. No values were seen.

## The design problem the plan is built around

The paper states that **all 176 controls are Italian** (Brescia). The patients are **44 + 44 Italian (P1-P44,
H1-H44)** and **44 + 44 American (P45-P88, H45-H88)**. The American samples were drawn at home and shipped at 4 degC
for up to 24 h before processing in Baltimore. Any patient-vs-control contrast that includes the American patients
is therefore confounded by country, recruitment channel and pre-analytical handling.

**The primary analysis uses the Italian patients only**, whose site and sample handling match the controls'. The
country mapping comes from the ID ranges in the paper text, not from a column in the file.

## Primary analysis (one)

* **Question:** does the serum Olink proteome separate Italian hEDS/HSD patients from Italian healthy controls, when
  evaluated on held-out people?
* **Sample:** Italian patients P1-P44 + H1-H44 (n = 88, label 1) vs controls C1-C176 (n = 176, label 0). The Olink
  plate-control, LOD, missing-frequency and normalisation rows are dropped.
* **Features:** all NPX assays, prepared in three steps:
  1. The two proteins duplicated across panels are reduced to one (the first occurrence is kept).
  2. Assays with more than 25% missing values across the whole NPX matrix (all 352 samples, labels ignored) are
     dropped.
  3. Remaining missing values are imputed with the training-fold median inside each fold.
* **Model:** L2 logistic regression, C = 0.1, features standardised inside each fold. No tuning and no feature
  selection.
* **Evaluation:**
  * Stratified 5-fold CV repeated 20 times (fold seed 20260923 + repeat). AUROC is computed on the pooled
    out-of-fold predictions of each repeat and averaged over repeats.
  * 95% CI from 2,000 group-stratified person bootstraps of the out-of-fold predictions (seed 20260923).
  * **Label-permutation null:** 1,000 label shuffles, with the whole CV re-run for each shuffle. One-sided
    p = (1 + number of null AUROCs >= observed) / 1,001.
* **Decision rule:** "the serum proteome discriminates hEDS/HSD from controls (Italian, same site)" if the AUROC CI
  lower bound is above 0.5 **and** the permutation p is below 0.05. Otherwise "not shown".
* **What a positive result would not mean:** it is **not** a diagnostic claim. The controls are healthy volunteers,
  not people with other causes of joint pain or hypermobility.

## Pre-specified secondary analyses (all reported whatever they show)

* **S1 per-protein effects (Italian primary sample):**
  * For every retained assay: Hedges' g (patients minus controls; NPX is already log2), a 2,000-resample stratified
    bootstrap CI, a two-sided Mann-Whitney p, and BH q across all retained assays.
  * Report the number of assays with q < 0.05, and their overlap with the paper's 69 DEP names from the pooled
    cohort. The paper's list is a published claim and is read only after S1 has been computed.
* **S2 site negative control:** American vs Italian **patients** (88 vs 88), with the same model and CV. The paper
  reports no per-protein differences here. If this AUROC's CI lower bound is above 0.5, country or sample handling is
  detectable in the proteome. Every analysis that pools American patients against Italian controls is then labelled
  site-confounded.
* **S3 whole-cohort replication of the paper's contrast:** all 176 patients vs 176 controls, same model and CV.
  Labelled "site-confounded by design" whatever it shows.
* **S4 hEDS vs HSD (Italian only, 44 vs 44):** same model and CV. The paper reports no difference, so an AUROC near
  0.5 is expected.
* **S5 sensitivity:** the primary repeated after excluding patients whose `RELATIONSHIP WITH PROBAND` is not `P`
  (relatives of a proband), to remove within-family dependence.
* **S6 sensitivity at 90% specificity** for the primary model. The threshold is the 90th percentile of control scores
  in each repeat.

## Covariates

Patient age and sex are in the clinical sheets, but **control age and sex are not released**. No age or sex
adjustment of the patient-vs-control contrast is therefore possible.

The paper reports that controls are older on average (mean 42 vs 36-38 years) and that about 91% of every group is
female. This limitation is stated with every result.

## Outputs

* Processed tables:
  * `participants__heds_hsd_olink_serum_cinquina2026`: 352 people, with group, site (from the ID range), patient
    age, sex and relationship, and comorbidity flags;
  * `participant_labs__heds_hsd_olink_serum_cinquina2026`: long format, one row per participant x assay NPX;
  * `phenotype_signatures__heds_hsd_olink_serum_cinquina2026`.
* Result tables: `results/tables/heds_hsd_olink_serum_cinquina2026_*`.
* Report: `results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md`.

## Amendment 1 (2026-09-24 21:14 UTC, before any protein value was compared between groups; the run started ~21:16)

The NPX sheet has two column blocks that were not visible when the plan above was written:
* per-panel `Plate ID` (four plates, SP230158-SP230161, 88 samples each);
* per-panel `QC Warning` (Pass / Warning).

Two decisions follow from them:
* **QC warnings.** As in the paper, an NPX value is set to missing when its panel's `QC Warning` for that sample is
  `Warning` (1-8 samples per panel). This happens before the > 25% missingness filter, and the imputation is
  unchanged.
* **Plate balance.** A plate x group cross-tab (labels only) shows every group on all four plates:

  | plate | controls | Italian hEDS | Italian HSD | US hEDS | US HSD |
  |---|---|---|---|---|---|
  | SP230158 | 42 | 14 | 9 | 11 | 12 |
  | SP230159 | 45 | 16 | 8 | 12 | 7 |
  | SP230160 | 40 | 5 | 14 | 15 | 14 |
  | SP230161 | 49 | 9 | 13 | 6 | 11 |

  Plate is therefore not confounded with diagnosis. As a pre-specified extra check (**S7**), the primary model is
  also run with plate indicators added as covariates. The primary statistic is unchanged.
