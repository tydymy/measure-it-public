# Analysis plan: person-level modelling of the mapMECFS participant-linked omics

Written 2026-09-24, **before** any analysis below was run. Code: `src/measure_it/omics/person_linked.py`
(`uv run python -m measure_it.omics.person_linked`). Results: `results/LINKED_OMICS_RESULTS.md`,
`results/tables/linked_omics_*.csv`, `results/figures/drafts/linked_omics_*`, processed partition
`phenotype_signatures__mapmecfs_omics`. Seed: `config.SEED` (20260923) for every stochastic step.

## What was looked at before writing this plan

Only structure, never outcomes: `docs/MAP_MECFS_LINKAGE_AUDIT.md`, `data/raw/mapmecfs_nih_pi_mecfs/DATA_AUDIT.md`,
the schemas and row counts of `participants__mapmecfs`, `participant_omics_linked__mapmecfs` and
`condition_molecular_evidence__mapmecfs`, the number of people per layer and group (from the audit), sex by group
(HV 11 F / 12 M / 1 unknown; PI-ME/CFS 12 F / 10 M / 1 discordant), matrix density (no missing cells in any layer;
PBMC has 25 pseudo-autosomal genes listed twice with identical values; CSF metabolomics has no zeros), value scales
(RNA-seq raw read counts, SomaScan RFU, metabolomics batch-normalised imputed relative abundance) and the
feature-name format. No classifier, test statistic or group comparison was computed before this file was written.

## Why this analysis, and what it can and cannot say

This is the one dataset in the project in which several omics layers belong to the **same people** (NIH study
number; `docs/MAP_MECFS_LINKAGE_AUDIT.md`). The question is whether any layer carries person-level signal that
separates PI-ME/CFS from healthy volunteers under honest cross-validation, not whether the published group-level
differences exist.

It cannot say anything about wearables: **no wearable, actigraphy, HRV, orthostatic or CPET measurement is linked to
these people in open data** (0 wearable x omics identifier pairs). It is not a device evaluation. It is also a
small, single-site cohort (NIH Clinical Center), and adjudication status is not in the open files: the deposits hold
21 PI-ME/CFS + 21 HV SomaScan samples while the paper's analytic cohort is 17 PI-ME/CFS + 21 HV, so up to 4 deposited
"PI-ME/CFS" participants may not be in the adjudicated cohort, and no open file says which. The published
per-analyte results (`condition_molecular_evidence__mapmecfs`) were computed on overlapping people, so agreement with
them is a **reproduction check, not independent replication**. With n = 21-42, most layers are expected to be null.

## Population and labels

Unit = participant (`mapmecfs:MECFS_<study number>`). Label: `group_label` (PI-ME/CFS = 1, HV = 0) as deposited;
concordant across all deposits for all 47 people.

| layer | people | HV | PI-ME/CFS | features | value |
|---|---|---|---|---|---|
| csf_somalogic | 42 | 21 | 21 | 1,317 aptamers | RFU (hybNorm, plateScale, medNorm) |
| serum_somalogic | 42 | 21 | 21 | 1,317 aptamers | RFU |
| pbmc_rnaseq | 27 | 15 | 12 | 18,369 genes (18,344 after de-duplicating PAR genes) | read counts |
| muscle_rnaseq | 25 | 12 | 13 | 11,423 genes | read counts |
| csf_metabolomics | 21 | 11 | 10 | 445 metabolites | batch-normalised imputed relative abundance; only ID-bearing rows (all Birth Sex = 1) |

Covariates for the baseline: **sex** as recorded in that layer's own deposit (`sex_by_source[layer]`; MECFS_311 is
Female in SomaScan and Male in RNA-seq and is used as recorded per layer); CSF metabolomics has only Birth Sex = 1
rows, so sex is constant there and the baseline is age only. **Age**: `age_years_somalogic` (the only open age; missing
for 5 people outside the SomaScan runs), median-imputed inside each training fold.

## Preprocessing (all inside the training fold; applied unchanged to the test fold)

* RNA-seq: counts -> log2(CPM + 1) with each sample's own library size (per-sample, no cross-sample information);
  keep genes with CPM >= 1 in >= 50% of the **training** samples; drop PAR duplicates (keep the first).
* SomaScan: log2(RFU).
* Metabolomics: log2(value) (no zeros in the data).
* Then: univariate Welch t-test on the training fold -> keep the **top k = 50** features by |t| -> standardise (mean
  and SD from the training fold) -> classifier.

## Models (fixed hyper-parameters, no tuning)

1. **Primary: L2 logistic regression**, C = 1.0 (sklearn default), lbfgs, class_weight = balanced.
2. Secondary: elastic-net logistic, C = 1.0, l1_ratio = 0.5, saga, max_iter 10,000, class_weight = balanced.
3. Baseline: L2 logistic (C = 1.0, balanced) on [age, female] (age only for CSF metabolomics).

Hyper-parameters are fixed because an inner tuning loop on 17-34 training people is itself noise, and AUROC of an
L2 logistic model depends only on the direction of the weight vector, which is insensitive to C at this scale.

## Validation

* 10 x 5-fold repeated stratified CV (each person in exactly one test fold per repeat). Per repeat, AUROC on the
  pooled out-of-fold (OOF) predictions of all people; point estimate = mean over the 10 repeats.
* 95% CI: stratified participant bootstrap of the OOF predictions (2,000 resamples; cases and controls resampled
  separately; for each resample the AUROC of every repeat is recomputed and averaged) -> percentile interval. This
  covers evaluation-sample variability, not refit variability; the across-repeat range is reported too.
* **Label-permutation null: 1,000 permutations** of the label, each re-running the full 10 x 5 procedure (including
  preprocessing and feature selection) -> one-sided p = (1 + #null >= observed) / 1,001. Run for the primary model,
  the elastic net, the baseline and the fusion models.
* A layer **passes** if the primary L2 model's permutation p < 0.01 (Bonferroni over the 5 layers). Raw p is always
  reported. **Top features are reported only for a layer that passes.**
* Omics vs baseline: paired bootstrap of ΔAUROC (primary L2 minus age/sex) on the same people.

## Pre-specified sensitivity analyses (reported in full whatever they show)

* S1 **sex/age-residualised features**: each log-scale feature is regressed on [female, age] in the training fold
  and replaced by its residual (the training fit is applied to the test fold) before selection; primary L2 model;
  1,000 permutations. Guards against the classifier reading sex-linked genes or age.
* S2 **k = 10 and k = 200**, primary L2 model, AUROC + bootstrap CI only (no permutation). Shows whether the
  primary conclusion depends on k; the primary result remains k = 50 whatever S2 shows.

## Late fusion (people who have every fused layer)

Within each outer training fold, each layer's primary pipeline is fitted on the training people; test-fold
predictions are fused by the **unweighted mean of the per-layer predicted log-odds**. No stacking (too few people).
10 x 5 CV, bootstrap CI, 1,000 permutations on the fusion cohort. Pre-specified fusion sets and expected n:

| fusion | people (HV / PI) from the audit overlaps | status |
|---|---|---|
| F1 csf_somalogic + serum_somalogic | 42 (21 / 21) | primary fusion |
| F2 csf_somalogic + serum_somalogic + pbmc_rnaseq | 27 | secondary |
| F3 csf_somalogic + serum_somalogic + muscle_rnaseq | 21 | secondary |
| F4 pbmc_rnaseq + muscle_rnaseq | 17 | secondary |
| F5 all five layers | 9 | reported as n only: not analysed (fewer than 6 people in a class is below any usable 5-fold CV) |

Each fusion is compared with the best single layer **on the same people** (its layers' primary models refitted on
the fusion cohort), reported as ΔAUROC with paired bootstrap CI; no permutation for Δ.

## Univariate differential analysis per layer (all people in the layer)

* Primary test: Welch t-test on the log-scale values (RNA-seq: log2(CPM + 1) of genes with CPM >= 1 in >= 50% of
  the layer's samples). Effect: difference in means of log2 values (PI-ME/CFS minus HV) and Hedges' g with a normal
  approximation 95% CI. BH-FDR within layer; FDR-significant = q < 0.05.
* Secondary: OLS of each feature on [group, female, age (median-imputed)] -> group coefficient p, BH-FDR within
  layer (CSF metabolomics: [group, age]).
* Count of FDR-significant features per layer is reported regardless of the classifier result.

## Agreement with the paper's published direction (reproduction check on overlapping people)

Join our per-feature effect to `condition_molecular_evidence__mapmecfs`: PBMC -> Supplementary Data 16A (Ensembl
id), muscle -> 19A (gene symbol), serum SomaScan -> 17A and CSF SomaScan -> 17B (SomaLogic SeqId; all tested
aptamers published), CSF metabolomics -> 14A (metabolite name; significant results only). Statistics:

* sign concordance of our mean difference with the published direction among the paper's nominally significant
  analytes (p < 0.05), with an exact binomial test against 0.5;
* for 17A/17B (all analytes published): Spearman correlation of our effect with the published effect across all
  shared aptamers.

The paper used 15 PI / 18 HV for 17A/17B and 14 / 15 for 16A, a subset of our people plus possibly others; this is
labelled a reproduction check.

## Cross-layer consistency: CSF vs serum SomaScan (same 42 people, same 1,317 aptamers)

* Spearman correlation of per-aptamer Welch t statistics (CSF vs serum). Null: **1,000 joint label permutations**
  (the same permuted labels applied to both layers, preserving each person's CSF-serum coupling) -> one-sided p.
* Overlap of the top-50 aptamers by |t| in each fluid, with the same joint-permutation null for the overlap count.
* Among aptamers FDR-significant in either fluid: count significant in both and sign agreement.
* Descriptive: per-aptamer across-person Spearman correlation of CSF with serum log2 RFU (how much CSF level tracks
  serum level at all), summarised as a distribution.

## Outputs

* `results/tables/linked_omics_model_performance.csv` (layer x model: n, AUROC, CI, repeat range, permutation p)
* `linked_omics_permutation_null.csv`, `linked_omics_fusion.csv`, `linked_omics_sensitivity.csv`,
  `linked_omics_univariate.csv` (all features), `linked_omics_univariate_summary.csv`,
  `linked_omics_paper_concordance.csv`, `linked_omics_csf_serum_consistency.csv`, `linked_omics_top_features.csv`
  (only passing layers), `linked_omics_cohort.csv`.
* Figures: `linked_omics_auroc_forest.png` (AUROC with CI and permutation null band per layer/model),
  `linked_omics_csf_vs_serum.png` (t-statistic scatter).
* `phenotype_signatures__mapmecfs_omics`: one row per layer x feature (Hedges' g, CI, Welch p, BH q) plus one model
  row per layer x model (feature `model_auroc:<layer>:<model>`, effect_measure "cross-validated AUROC"), columns as
  the canonical `phenotype_signatures`; condition_id `me_cfs`; label_basis "NIH intramural adjudicated PI-ME/CFS vs
  healthy volunteer (per paper)"; data_layer person (aggregate statistics of person-level data).

## Language

Results are "person-level omics signal in one small single-site cohort", never "biomarker", "diagnostic" or
"patient multi-omics with wearables". A passing layer is a candidate for replication, not a validated test.
