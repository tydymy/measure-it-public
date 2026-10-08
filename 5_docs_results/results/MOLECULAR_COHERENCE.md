# Molecular coherence (Test 3) and the phenotype -> measurable-biology graph

Generated 2026-09-25T22:08:18+00:00 by `measure_it.omics.graph.write_report` from files written by `measure_it.omics.coherence` (built 2026-09-25T22:05:21+00:00) and `measure_it.omics.graph`. Plan (written before the analysis): `docs/ANALYSIS_PLAN_MOLECULAR.md`. Every number below is read from `results/tables/test3_*.csv`, `results/tables/fig5_*.csv` or `data/processed/measurable_biology.parquet`.

**Scope.** Condition-level molecular enrichment only: gene lists that public databases attach to a disease concept. Nothing here is measured on, or linked to, any participant, and no result describes a person. No combined omics score is computed.

Versions: open_targets: Open Targets Platform data release 26.06; API 26.6.3; gwas_catalog: GWAS Catalog data release 2026-09-13; REST API 2.0 (api release 2025-08-01); EFO v3.93.0; gene build GRCh38.p14; mapmecfs: Walitt B, et al. Deep phenotyping of post-infectious myalgic encephalomyelitis/chronic fatigue syndrome. Nat Commun 15, 907 (2024); Supplementary Data 1-24 (MOESM4); reactome: Reactome release 97 (https://download.reactome.org/97/ReactomePathways.gmt.zip, fetched 2026-09-24T04:34:10+00:00); hgnc: HGNC complete set (Last-Modified Fri, 18 Sep 2026 13:36:42 GMT, fetched 2026-09-23T23:08:00+00:00).

**Independent review (2026-09-23).** The gene-level and pathway-level statistics were recomputed independently and match. Corrections made in review: (1) the per-condition disagreement list used a fixed empirical-p threshold of 0.05 that cannot be reached with 16-18 cross-condition pairs, so pairs whose same-condition overlap is above every cross pair were listed as disagreements; the criterion is now p <= max(0.05, 1/(1 + n_cross)). (2) The flag-concordance result is now reported with a Holm adjustment over the two pre-specified definitions, the sources behind each concordant pair, and exploratory definitions without literature and post-hoc-flagged evidence. (3) The drug records behind Open Targets 'other' support are named. (4) Post-hoc flags are drawn in the Figure 5 drafts. (5) The query tool separates a tested negative from missing data. Details are in the plan's 'Deviations' section.

## Headline results

Positive and negative results are listed together. Permutation p-values test whether sources agree more for the **same** condition than for **different** conditions (condition-label shuffle).

* Gene level, ot_genetic x gwas_catalog (7 conditions): same-condition mean log2 fold 2.849 vs cross-condition 0.213; difference 2.604 (95% CI 1.498 to 3.560); permutation p = 2.0e-04 (BH across pairs 3.0e-04) -> agreement **above the cross-condition baseline**.
* Gene level, ot_genetic x ot_literature (10 conditions): same-condition mean log2 fold 1.729 vs cross-condition 0.680; difference 1.039 (95% CI 0.428 to 1.615); permutation p = 1.0e-04 (BH across pairs 2.0e-04) -> agreement **above the cross-condition baseline**.
* Gene level, ot_genetic x ot_other (10 conditions): same-condition mean log2 fold 1.510 vs cross-condition 0.281; difference 1.228 (95% CI 0.436 to 2.072); permutation p = 1.0e-04 (BH across pairs 2.0e-04) -> agreement **above the cross-condition baseline**.
* Gene level, gwas_catalog x ot_literature (8 conditions): same-condition mean log2 fold 1.197 vs cross-condition 0.443; difference 0.700 (95% CI -0.072 to 1.327); permutation p = 0.0063 (BH across pairs 0.0076) -> agreement **above the cross-condition baseline by the permutation test, but the bootstrap CI of the difference includes 0**.
* Gene level, gwas_catalog x ot_other (8 conditions): same-condition mean log2 fold 0.039 vs cross-condition 0.149; difference -0.101 (95% CI -0.675 to 0.379); permutation p = 0.6272 (BH across pairs 0.6272) -> agreement **NOT above the cross-condition baseline**.
* Gene level, ot_literature x ot_other (12 conditions): same-condition mean log2 fold 2.594 vs cross-condition 1.586; difference 1.008 (95% CI 0.571 to 1.473); permutation p = 5.0e-05 (BH across pairs 2.0e-04) -> agreement **above the cross-condition baseline**.
* Pathway level, ot_genetic x gwas_catalog (5 conditions): same-condition Spearman rho 0.254 vs 0.077 (difference 95% CI 0.018 to 0.322); permutation p = 0.0167 (BH across pairs 0.0250) -> **above the cross-condition baseline**.
* Pathway level, ot_genetic x ot_literature (8 conditions): same-condition Spearman rho 0.177 vs 0.120 (difference 95% CI 0.007 to 0.112); permutation p = 4.7e-04 (BH across pairs 9.4e-04) -> **above the cross-condition baseline**.
* Pathway level, ot_genetic x ot_other (8 conditions): same-condition Spearman rho 0.194 vs 0.092 (difference 95% CI 0.051 to 0.151); permutation p = 5.0e-05 (BH across pairs 1.5e-04) -> **above the cross-condition baseline**.
* Pathway level, gwas_catalog x ot_literature (6 conditions): same-condition Spearman rho 0.108 vs 0.103 (difference 95% CI -0.027 to 0.047); permutation p = 0.3958 (BH across pairs 0.4750) -> **NOT above the cross-condition baseline**.
* Pathway level, gwas_catalog x ot_other (6 conditions): same-condition Spearman rho 0.036 vs 0.067 (difference 95% CI -0.058 to -0.001); permutation p = 0.8361 (BH across pairs 0.8361) -> **NOT above the cross-condition baseline**.
* Pathway level, ot_literature x ot_other (11 conditions): same-condition Spearman rho 0.286 vs 0.204 (difference 95% CI 0.016 to 0.151); permutation p = 5.0e-05 (BH across pairs 1.5e-04) -> **above the cross-condition baseline**.
* mapMECFS (me_cfs) x ot_genetic: overlap 3 genes, log2 fold -0.660; me_cfs ranks 9 of 10 conditions (share of conditions matching at least as well = 0.900; best match: dysautonomia) -> me_cfs is **not** the best-matching condition.
* mapMECFS (me_cfs) x gwas_catalog: overlap 0 genes, log2 fold -1.198; me_cfs ranks 7 of 8 conditions (share of conditions matching at least as well = 1; best match: lyme_disease) -> me_cfs is **not** the best-matching condition.
* mapMECFS (me_cfs) x ot_literature: overlap 138 genes, log2 fold 0.274; me_cfs ranks 3 of 12 conditions (share of conditions matching at least as well = 0.250; best match: ptlds) -> me_cfs is **not** the best-matching condition.
* mapMECFS (me_cfs) x ot_other: overlap 7 genes, log2 fold -0.152; me_cfs ranks 4 of 12 conditions (share of conditions matching at least as well = 0.333; best match: eds_hsd) -> me_cfs is **not** the best-matching condition.
* Systems vs curated clinical flags (any source; pre-specified): 20 of 28 expected condition-system pairs observed (null mean 18.9); permutation p = 0.3391; Holm-adjusted over the two pre-specified definitions 0.3391.
* Systems vs curated clinical flags (condition specific; pre-specified): 18 of 28 expected condition-system pairs observed (null mean 14.0); permutation p = 0.0247; Holm-adjusted over the two pre-specified definitions 0.0495.
* Systems vs curated clinical flags (non literature sources; exploratory (analyst, post hoc)): 14 of 28 expected condition-system pairs observed (null mean 12.6); permutation p = 0.2917.
* Systems vs curated clinical flags (condition specific non literature; exploratory (added in review, post hoc)): 13 of 28 expected condition-system pairs observed (null mean 10.9); permutation p = 0.1593.
* Systems vs curated clinical flags (condition specific non literature unflagged; exploratory (added in review, post hoc)): 3 of 28 expected condition-system pairs observed (null mean 2.4); permutation p = 0.4500.
* **Reading of the flag-concordance result (review).** The only nominally significant pre-specified definition is the 'condition-specific' one, and 'condition-specific' here means only that the system's fold is above the median of the other conditions (about half of all tests pass by construction). Of its 18 concordant pairs, 5 rest on Europe PMC co-mention alone, 10 include Open Targets 'other' support that is mostly drug targets of trialled drugs (post-hoc flag), and 0 involve genetic evidence. Literature co-mention and trialled-drug targets both reflect what is already believed clinically about a condition, which is what the curated flags also encode, so this concordance is partly circular. Restricted to condition-specific, non-literature, unflagged support: 3 of 28 expected pairs observed (null mean 2.4; p = 0.4500; exploratory, added in review). The public molecular resources therefore do **not** independently recover the curated clinical profile. Across the sensitivity analyses the condition-specific concordance p ranges 0.007-0.126; it is >= 0.05 in S1_ot_indirect (p = 0.126), S8_ot_other_without_drug_targets (p = 0.107) (S8 is post hoc).
* Enriched Reactome pathways (FDR < 0.05, overlap >= 2) summed over conditions, by source: ot_genetic: 59 across 8 condition sets; gwas_catalog: 13 across 6 condition sets; ot_literature: 2568 across 11 condition sets; ot_other: 682 across 11 condition sets; mapmecfs: 0 across 1 condition sets.
* Demo cluster, long_covid: literature_only: immune, vascular/endothelial, connective tissue. Systems with non-literature support after removing post-hoc-flagged evidence: none.
* Demo cluster, me_cfs: multi_source: neuronal [ot_other_mostly_drug_targets]; single_source: autonomic/cardiac [ot_other_mostly_drug_targets], mitochondrial/metabolic [ot_other_mostly_drug_targets]; literature_only: immune, vascular/endothelial, connective tissue. Systems with non-literature support after removing post-hoc-flagged evidence: none.
* Demo cluster, pots: single_source: autonomic/cardiac [ot_other_mostly_drug_targets], neuronal [ot_other_mostly_drug_targets]; literature_only: immune, vascular/endothelial. Systems with non-literature support after removing post-hoc-flagged evidence: none.

## 1. Source counts per condition

Gene counts are HGNC protein-coding genes in the primary scope (plan section 3). GEO/SRA counts are study metadata from text queries and are never used as gene evidence.

| condition_id | ot_unique_targets_all_ids | ot_unique_targets_direct | gwas_associations | gwas_studies | geo_series_in_scope | sra_studies | mapmecfs_rows | genes_ot_genetic | genes_gwas_catalog | genes_ot_literature | genes_ot_other | genes_mapmecfs | n_gene_sources_nonempty |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| long_covid | 347 | 347 | 0 | 4 | 52 | 45 | 0 | 2 | 0 | 300 | 46 | 0 | 3 |
| me_cfs | 1775 | 1775 | 9 | 19 | 28 | 8 | 14565 | 70 | 9 | 1587 | 109 | 1387 | 5 |
| pots | 231 | 231 | 0 | 0 | 0 | 3 | 0 | 1 | 0 | 135 | 101 | 0 | 3 |
| dysautonomia | 3246 | 351 | 3 | 8 | 49 | 10 | 0 | 28 | 1 | 313 | 42 | 0 | 4 |
| fibromyalgia | 595 | 595 | 50 | 41 | 13 | 9 | 0 | 0 | 13 | 474 | 147 | 0 | 3 |
| eds_hsd | 3273 | 923 | 0 | 3 | 18 | 11 | 0 | 41 | 0 | 595 | 537 | 0 | 3 |
| mcas | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| lyme_disease | 1297 | 1297 | 3 | 3 | 27 | 12 | 0 | 11 | 4 | 436 | 635 | 0 | 4 |
| ptlds | 4 | 4 | 0 | 0 | 4 | 3 | 0 | 0 | 0 | 3 | 1 | 0 | 2 |
| post_infectious_syndrome | 8240 | 1 | 0 | 0 | 12 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| migraine | 2930 | 2218 | 854 | 92 | 20 | 13 | 0 | 526 | 380 | 1770 | 374 | 0 | 4 |
| ibs | 2050 | 2050 | 483 | 85 | 21 | 15 | 0 | 175 | 188 | 1833 | 120 | 0 | 4 |
| gastroparesis | 236 | 236 | 10 | 5 | 5 | 4 | 0 | 47 | 9 | 168 | 33 | 0 | 4 |
| endometriosis | 3561 | 3511 | 1029 | 71 | 218 | 123 | 0 | 268 | 248 | 3016 | 419 | 0 | 4 |

Identifier harmonisation (distinct raw identifiers per source and their outcome):

| source | n_distinct_raw_identifiers | n_ensembl | n_unmapped | n_non_protein_coding | n_approved | n_ambiguous_alias | n_previous_symbol | n_ambiguous_previous_symbol |
|---|---|---|---|---|---|---|---|---|
| ot_genetic | 1106 | 1100 | 5 | 1 | 0 | 0 | 0 | 0 |
| gwas_catalog | 1299 | 0 | 2 | 533 | 757 | 3 | 3 | 1 |
| ot_literature | 5606 | 5417 | 8 | 175 | 6 | 0 | 0 | 0 |
| ot_other | 2046 | 1765 | 98 | 183 | 0 | 0 | 0 | 0 |
| mapmecfs | 1987 | 536 | 219 | 355 | 863 | 0 | 14 | 0 |

Conditions with no protein-coding gene from any source: mcas, post_infectious_syndrome (reported as UNKNOWN / NOT AVAILABLE downstream, never imputed).

## 2. Gene-level agreement (Q1)

Background N = 19297 HGNC approved protein-coding genes. log2 fold = log2((k + 0.5)/(E + 0.5)).

| source_a | source_b | n_conditions | same_mean | cross_mean_all_pairs | mean_diff_same_minus_cross | diff_ci95_low | diff_ci95_high | perm_p | n_perm | perm_exact | perm_p_bh_across_pairs | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ot_genetic | gwas_catalog | 7 | 2.849 | 0.213 | 2.604 | 1.498 | 3.560 | 2.0e-04 | 5040 | True | 3.0e-04 | ok |
| ot_genetic | ot_literature | 10 | 1.729 | 0.680 | 1.039 | 0.428 | 1.615 | 1.0e-04 | 20000 | False | 2.0e-04 | ok |
| ot_genetic | ot_other | 10 | 1.510 | 0.281 | 1.228 | 0.436 | 2.072 | 1.0e-04 | 20000 | False | 2.0e-04 | ok |
| ot_genetic | mapmecfs | 1 | -0.660 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| gwas_catalog | ot_literature | 8 | 1.197 | 0.443 | 0.700 | -0.072 | 1.327 | 0.006 | 40320 | True | 0.008 | ok |
| gwas_catalog | ot_other | 8 | 0.039 | 0.149 | -0.101 | -0.675 | 0.379 | 0.627 | 40320 | True | 0.627 | ok |
| gwas_catalog | mapmecfs | 1 | -1.198 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| ot_literature | ot_other | 12 | 2.594 | 1.586 | 1.008 | 0.571 | 1.473 | 5.0e-05 | 20000 | False | 2.0e-04 | ok |
| ot_literature | mapmecfs | 1 | 0.274 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| ot_other | mapmecfs | 1 | -0.152 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |

Jaccard (size-dependent; for reading only):

| source_a | source_b | n_conditions | same_mean | cross_mean_all_pairs | perm_p |
|---|---|---|---|---|---|
| ot_genetic | gwas_catalog | 7 | 0.1277 | 0.0033 | 2.0e-04 |
| ot_genetic | ot_literature | 10 | 0.0239 | 0.0073 | 5.0e-05 |
| ot_genetic | ot_other | 10 | 0.0258 | 0.0044 | 5.0e-05 |
| ot_genetic | mapmecfs | 1 | 0.0021 | n/a | n/a |
| gwas_catalog | ot_literature | 8 | 0.0199 | 0.0062 | 5.0e-05 |
| gwas_catalog | ot_other | 8 | 0.0066 | 0.0034 | 0.0375 |
| gwas_catalog | mapmecfs | 1 | 0.0000 | n/a | n/a |
| ot_literature | ot_other | 12 | 0.0695 | 0.0254 | 5.0e-05 |
| ot_literature | mapmecfs | 1 | 0.0487 | n/a | n/a |
| ot_other | mapmecfs | 1 | 0.0047 | n/a | n/a |

Same-condition overlaps (every condition with both sets non-empty):

| source_a | source_b | condition_a | n_a | n_b | overlap | jaccard | log2_fold | hypergeom_p | hypergeom_fdr_same_condition |
|---|---|---|---|---|---|---|---|---|---|
| gwas_catalog | mapmecfs | me_cfs | 9 | 1387 | 0 | 0 | -1.198 | 1 | 1 |
| gwas_catalog | ot_literature | dysautonomia | 1 | 313 | 0 | 0 | -0.046 | 1 | 1 |
| gwas_catalog | ot_literature | endometriosis | 248 | 3016 | 115 | 0.037 | 1.557 | 1.5e-30 | 6.6e-30 |
| gwas_catalog | ot_literature | fibromyalgia | 13 | 474 | 3 | 0.006 | 2.095 | 0.004 | 0.006 |
| gwas_catalog | ot_literature | gastroparesis | 9 | 168 | 1 | 0.006 | 1.375 | 0.076 | 0.114 |
| gwas_catalog | ot_literature | ibs | 188 | 1833 | 60 | 0.031 | 1.721 | 8.1e-18 | 2.7e-17 |
| gwas_catalog | ot_literature | lyme_disease | 4 | 436 | 2 | 0.005 | 2.082 | 0.003 | 0.006 |
| gwas_catalog | ot_literature | me_cfs | 9 | 1587 | 0 | 0 | -1.311 | 1 | 1 |
| gwas_catalog | ot_literature | migraine | 380 | 1770 | 151 | 0.076 | 2.099 | 2.9e-59 | 2.9e-58 |
| gwas_catalog | ot_other | dysautonomia | 1 | 42 | 0 | 0 | -0.006 | 1 | 1 |
| gwas_catalog | ot_other | endometriosis | 248 | 419 | 13 | 0.020 | 1.198 | 0.003 | 0.006 |
| gwas_catalog | ot_other | fibromyalgia | 13 | 147 | 0 | 0 | -0.261 | 1 | 1 |
| gwas_catalog | ot_other | gastroparesis | 9 | 33 | 0 | 0 | -0.044 | 1 | 1 |
| gwas_catalog | ot_other | ibs | 188 | 120 | 0 | 0 | -1.739 | 1 | 1 |
| gwas_catalog | ot_other | lyme_disease | 4 | 635 | 0 | 0 | -0.337 | 1 | 1 |
| gwas_catalog | ot_other | me_cfs | 9 | 109 | 0 | 0 | -0.140 | 1 | 1 |
| gwas_catalog | ot_other | migraine | 380 | 374 | 24 | 0.033 | 1.639 | 4.5e-07 | 1.1e-06 |
| ot_genetic | gwas_catalog | dysautonomia | 28 | 1 | 1 | 0.036 | 1.581 | 0.001 | 0.003 |
| ot_genetic | gwas_catalog | endometriosis | 268 | 248 | 83 | 0.192 | 4.404 | 2.7e-94 | 3.1e-93 |
| ot_genetic | gwas_catalog | gastroparesis | 47 | 9 | 0 | 0 | -0.062 | 1 | 1 |
| ot_genetic | gwas_catalog | ibs | 175 | 188 | 70 | 0.239 | 4.999 | 3.1e-98 | 4.6e-97 |
| ot_genetic | gwas_catalog | lyme_disease | 11 | 4 | 2 | 0.154 | 2.315 | 1.8e-06 | 4.2e-06 |
| ot_genetic | gwas_catalog | me_cfs | 70 | 9 | 3 | 0.039 | 2.716 | 3.8e-06 | 8.6e-06 |
| ot_genetic | gwas_catalog | migraine | 526 | 380 | 172 | 0.234 | 3.990 | 1.2e-172 | 2.4e-171 |
| ot_genetic | mapmecfs | me_cfs | 70 | 1387 | 3 | 0.002 | -0.660 | 0.888 | 1 |
| ot_genetic | ot_literature | dysautonomia | 28 | 313 | 4 | 0.012 | 2.238 | 0.001 | 0.002 |
| ot_genetic | ot_literature | eds_hsd | 41 | 595 | 34 | 0.056 | 4.290 | 3.1e-45 | 2.0e-44 |
| ot_genetic | ot_literature | endometriosis | 268 | 3016 | 119 | 0.038 | 1.495 | 2.1e-29 | 8.7e-29 |
| ot_genetic | ot_literature | gastroparesis | 47 | 168 | 1 | 0.005 | 0.722 | 0.337 | 0.463 |
| ot_genetic | ot_literature | ibs | 175 | 1833 | 52 | 0.027 | 1.616 | 3.3e-14 | 9.9e-14 |
| ot_genetic | ot_literature | long_covid | 2 | 300 | 1 | 0.003 | 1.498 | 0.031 | 0.048 |
| ot_genetic | ot_literature | lyme_disease | 11 | 436 | 3 | 0.007 | 2.225 | 0.002 | 0.003 |
| ot_genetic | ot_literature | me_cfs | 70 | 1587 | 5 | 0.003 | -0.186 | 0.692 | 0.908 |
| ot_genetic | ot_literature | migraine | 526 | 1770 | 173 | 0.081 | 1.832 | 8.4e-54 | 7.0e-53 |
| ot_genetic | ot_literature | pots | 1 | 135 | 1 | 0.007 | 1.565 | 0.007 | 0.011 |
| ot_genetic | ot_other | dysautonomia | 28 | 42 | 5 | 0.077 | 3.294 | 3.6e-09 | 9.7e-09 |
| ot_genetic | ot_other | eds_hsd | 41 | 537 | 34 | 0.062 | 4.394 | 8.6e-47 | 6.4e-46 |
| ot_genetic | ot_other | endometriosis | 268 | 419 | 13 | 0.019 | 1.095 | 0.006 | 0.010 |
| ot_genetic | ot_other | gastroparesis | 47 | 33 | 1 | 0.013 | 1.370 | 0.077 | 0.114 |
| ot_genetic | ot_other | ibs | 175 | 120 | 3 | 0.010 | 1.140 | 0.096 | 0.138 |
| ot_genetic | ot_other | long_covid | 2 | 46 | 0 | 0 | -0.014 | 1 | 1 |
| ot_genetic | ot_other | lyme_disease | 11 | 635 | 0 | 0 | -0.786 | 1 | 1 |
| ot_genetic | ot_other | me_cfs | 70 | 109 | 1 | 0.006 | 0.744 | 0.328 | 0.461 |
| ot_genetic | ot_other | migraine | 526 | 374 | 52 | 0.061 | 2.295 | 2.3e-22 | 9.1e-22 |
| ot_genetic | ot_other | pots | 1 | 101 | 1 | 0.010 | 1.570 | 0.005 | 0.009 |
| ot_literature | mapmecfs | me_cfs | 1587 | 1387 | 138 | 0.049 | 0.274 | 0.010 | 0.016 |
| ot_literature | ot_other | dysautonomia | 313 | 42 | 29 | 0.089 | 4.642 | 7.0e-43 | 4.1e-42 |
| ot_literature | ot_other | eds_hsd | 595 | 537 | 224 | 0.247 | 3.718 | 2.2e-205 | 6.3e-204 |
| ot_literature | ot_other | endometriosis | 3016 | 419 | 171 | 0.052 | 1.378 | 4.5e-36 | 2.2e-35 |
| ot_literature | ot_other | fibromyalgia | 474 | 147 | 31 | 0.053 | 2.938 | 2.1e-20 | 7.4e-20 |
| ot_literature | ot_other | gastroparesis | 168 | 33 | 12 | 0.063 | 3.989 | 3.9e-17 | 1.2e-16 |
| ot_literature | ot_other | ibs | 1833 | 120 | 66 | 0.035 | 2.483 | 3.6e-36 | 1.9e-35 |
| ot_literature | ot_other | long_covid | 300 | 46 | 5 | 0.015 | 2.178 | 7.1e-04 | 0.002 |
| ot_literature | ot_other | lyme_disease | 436 | 635 | 59 | 0.058 | 2.003 | 1.7e-20 | 6.1e-20 |
| ot_literature | ot_other | me_cfs | 1587 | 109 | 26 | 0.016 | 1.485 | 5.6e-07 | 1.4e-06 |
| ot_literature | ot_other | migraine | 1770 | 374 | 291 | 0.157 | 3.066 | 6.6e-231 | 3.9e-229 |
| ot_literature | ot_other | pots | 135 | 101 | 11 | 0.049 | 3.253 | 1.2e-10 | 3.4e-10 |
| ot_literature | ot_other | ptlds | 3 | 1 | 0 | 0 | -4.5e-04 | 1 | 1 |
| ot_other | mapmecfs | me_cfs | 109 | 1387 | 7 | 0.005 | -0.152 | 0.675 | 0.906 |

mapMECFS (me_cfs only) against each condition's set of another source:

| other_source | status | n_conditions_compared | me_cfs_overlap | me_cfs_log2_fold | me_cfs_hypergeom_p | me_cfs_rank_of_n | rank_p | other_conditions_log2_fold_median | best_condition |
|---|---|---|---|---|---|---|---|---|---|
| ot_genetic | ok | 10 | 3 | -0.660 | 0.888 | 9 | 0.900 | -0.151 | dysautonomia |
| gwas_catalog | ok | 8 | 0 | -1.198 | 1 | 7 | 1 | -0.194 | lyme_disease |
| ot_literature | ok | 12 | 138 | 0.274 | 0.010 | 3 | 0.250 | 0.116 | ptlds |
| ot_other | ok | 12 | 7 | -0.152 | 0.675 | 4 | 0.333 | -0.283 | eds_hsd |

## 3. Pathway-level agreement (Q2, Reactome)

Local over-representation analysis over 1736 Reactome pathways (10-500 genes), background 11303 Reactome-annotated protein-coding genes, BH-FDR per gene set.

| condition_id | source | genes_in_background | enriched_pathways |
|---|---|---|---|
| dysautonomia | ot_genetic | 16 | 0 |
| dysautonomia | ot_literature | 266 | 35 |
| dysautonomia | ot_other | 40 | 10 |
| eds_hsd | ot_genetic | 35 | 48 |
| eds_hsd | ot_literature | 512 | 156 |
| eds_hsd | ot_other | 441 | 205 |
| endometriosis | gwas_catalog | 158 | 0 |
| endometriosis | ot_genetic | 161 | 1 |
| endometriosis | ot_literature | 2524 | 627 |
| endometriosis | ot_other | 324 | 117 |
| fibromyalgia | gwas_catalog | 10 | 0 |
| fibromyalgia | ot_literature | 418 | 148 |
| fibromyalgia | ot_other | 141 | 49 |
| gastroparesis | gwas_catalog | 5 | 13 |
| gastroparesis | ot_genetic | 27 | 0 |
| gastroparesis | ot_literature | 154 | 29 |
| gastroparesis | ot_other | 33 | 14 |
| ibs | gwas_catalog | 118 | 0 |
| ibs | ot_genetic | 94 | 0 |
| ibs | ot_literature | 1587 | 453 |
| ibs | ot_other | 114 | 49 |
| long_covid | ot_literature | 266 | 143 |
| long_covid | ot_other | 44 | 17 |
| lyme_disease | ot_genetic | 6 | 6 |
| lyme_disease | ot_literature | 391 | 193 |
| lyme_disease | ot_other | 487 | 60 |
| me_cfs | gwas_catalog | 5 | 0 |
| me_cfs | mapmecfs | 838 | 0 |
| me_cfs | ot_genetic | 35 | 0 |
| me_cfs | ot_literature | 1356 | 424 |
| me_cfs | ot_other | 103 | 49 |
| migraine | gwas_catalog | 232 | 0 |
| migraine | ot_genetic | 322 | 4 |
| migraine | ot_literature | 1456 | 327 |
| migraine | ot_other | 334 | 66 |
| pots | ot_literature | 127 | 33 |
| pots | ot_other | 95 | 46 |

| source_a | source_b | statistic | n_conditions | same_mean | cross_mean_all_pairs | mean_diff_same_minus_cross | diff_ci95_low | diff_ci95_high | perm_p | n_perm | status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ot_genetic | gwas_catalog | spearman_rho_neglog10p | 5 | 0.254 | 0.077 | 0.170 | 0.018 | 0.322 | 0.017 | 120 | ok |
| ot_genetic | gwas_catalog | jaccard_enriched | 5 | 0 | 0 | 0 | 0 | 0 | 1 | 120 | ok |
| ot_genetic | ot_literature | spearman_rho_neglog10p | 8 | 0.177 | 0.120 | 0.056 | 0.007 | 0.112 | 4.7e-04 | 40320 | ok |
| ot_genetic | ot_literature | jaccard_enriched | 8 | 0.037 | 0.006 | 0.031 | -0.002 | 0.092 | 0.013 | 40320 | ok |
| ot_genetic | ot_other | spearman_rho_neglog10p | 8 | 0.194 | 0.092 | 0.103 | 0.051 | 0.151 | 5.0e-05 | 40320 | ok |
| ot_genetic | ot_other | jaccard_enriched | 8 | 0.041 | 0.004 | 0.037 | -2.8e-04 | 0.087 | 0.009 | 40320 | ok |
| ot_genetic | mapmecfs | spearman_rho_neglog10p | 1 | 0.116 | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| ot_genetic | mapmecfs | jaccard_enriched | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 0 condition(s) with both sets; a cross-condition control needs >= 2 |
| gwas_catalog | ot_literature | spearman_rho_neglog10p | 6 | 0.108 | 0.103 | 0.006 | -0.027 | 0.047 | 0.396 | 720 | ok |
| gwas_catalog | ot_literature | jaccard_enriched | 6 | 0 | 0.003 | -0.003 | -0.006 | -7.6e-04 | 1 | 720 | ok |
| gwas_catalog | ot_other | spearman_rho_neglog10p | 6 | 0.036 | 0.067 | -0.029 | -0.058 | -0.001 | 0.836 | 720 | ok |
| gwas_catalog | ot_other | jaccard_enriched | 6 | 0 | 1.3e-04 | -1.7e-04 | -3.4e-04 | 0 | 1 | 720 | ok |
| gwas_catalog | mapmecfs | spearman_rho_neglog10p | 1 | -0.051 | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| gwas_catalog | mapmecfs | jaccard_enriched | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 0 condition(s) with both sets; a cross-condition control needs >= 2 |
| ot_literature | ot_other | spearman_rho_neglog10p | 11 | 0.286 | 0.204 | 0.082 | 0.016 | 0.151 | 5.0e-05 | 20000 | ok |
| ot_literature | ot_other | jaccard_enriched | 11 | 0.124 | 0.067 | 0.056 | 0.017 | 0.104 | 1.5e-04 | 20000 | ok |
| ot_literature | mapmecfs | spearman_rho_neglog10p | 1 | -0.059 | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| ot_literature | mapmecfs | jaccard_enriched | 1 | 0 | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| ot_other | mapmecfs | spearman_rho_neglog10p | 1 | -0.073 | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |
| ot_other | mapmecfs | jaccard_enriched | 1 | 0 | n/a | n/a | n/a | n/a | n/a | n/a | UNKNOWN / NOT AVAILABLE: 1 condition(s) with both sets; a cross-condition control needs >= 2 |

Cross-checks: hierarchy-derived top level agrees with the Open Targets `reactome_top_level_term` for 2246 of 2246 pathway ids (1). Reactome AnalysisService vs local ORA:

| condition_id | source | status | n_genes_submitted | n_pathways_compared | spearman_rho_p_values | n_enriched_local | n_fdr05_service | n_in_both | jaccard |
|---|---|---|---|---|---|---|---|---|---|
| me_cfs | ot_genetic | ok | 70 | 276 | 0.591 | 0 | 0 | 0 | n/a |
| migraine | gwas_catalog | ok | 380 | 948 | 0.867 | 0 | 1 | 0 | 0 |
| long_covid | ot_literature | ok | 300 | 975 | 0.913 | 143 | 158 | 128 | 0.740 |

## 4. Physiological systems (Q3)

A system is supported by a source when its gene set is enriched for the system's Reactome gene universe (BH FDR < 0.05) and at least one pathway of the system is enriched. `*` = condition-specific (above the median of the other conditions for the same source).

| condition_id | source | immune | neuronal | vascular_endothelial | autonomic_cardiac | mitochondrial_metabolic | connective_tissue | other |
|---|---|---|---|---|---|---|---|---|
| dysautonomia | ot_genetic | . | . | . | . | .* | . | . |
| dysautonomia | ot_literature | S | S* | S | S* | .* | .* | .* |
| dysautonomia | ot_other | .* | . | .* | S* | .* | .* | .* |
| eds_hsd | ot_genetic | .* | .* | .* | . | . | S* | . |
| eds_hsd | ot_literature | S | . | S | . | .* | S* | .* |
| eds_hsd | ot_other | .* | S | S* | S | . | S* | S* |
| endometriosis | gwas_catalog | .* | .* | .* | .* | . | .* | . |
| endometriosis | ot_genetic | . | .* | .* | .* | . | . | S* |
| endometriosis | ot_literature | S | . | S | . | . | S* | . |
| endometriosis | ot_other | .* | . | . | S | S* | .* | .* |
| fibromyalgia | gwas_catalog | . | . | . | . | .* | . | .* |
| fibromyalgia | ot_literature | S* | S* | S | .* | .* | S | S* |
| fibromyalgia | ot_other | . | S* | . | S* | S* | . | . |
| gastroparesis | gwas_catalog | S* | . | .* | . | . | . | . |
| gastroparesis | ot_genetic | .* | . | . | .* | .* | . | .* |
| gastroparesis | ot_literature | S | .* | . | .* | .* | . | .* |
| gastroparesis | ot_other | .* | S* | . | . | . | . | . |
| ibs | gwas_catalog | .* | .* | . | . | . | .* | . |
| ibs | ot_genetic | .* | .* | . | . | .* | .* | . |
| ibs | ot_literature | S* | S* | S | S | . | S | . |
| ibs | ot_other | . | S* | . | S* | . | . | . |
| long_covid | ot_literature | S* | . | S* | . | . | S* | . |
| long_covid | ot_other | . | . | .* | . | . | . | .* |
| lyme_disease | ot_genetic | .* | . | . | . | . | . | . |
| lyme_disease | ot_literature | S* | . | S* | . | . | S* | . |
| lyme_disease | ot_other | S* | . | S* | . | . | S* | . |
| me_cfs | gwas_catalog | . | .* | . | . | .* | . | .* |
| me_cfs | mapmecfs | . | . | . | . | . | . | . |
| me_cfs | ot_genetic | . | . | .* | .* | . | .* | .* |
| me_cfs | ot_literature | S* | S | S* | . | . | S* | . |
| me_cfs | ot_other | . | S* | . | S* | S* | . | . |
| migraine | gwas_catalog | . | . | .* | .* | .* | .* | .* |
| migraine | ot_genetic | . | .* | .* | S* | .* | .* | .* |
| migraine | ot_literature | S | S* | S* | S* | .* | S | .* |
| migraine | ot_other | .* | S* | S* | S* | . | .* | .* |
| pots | ot_literature | S | . | S* | .* | .* | . | S* |
| pots | ot_other | . | S* | . | S* | .* | . | . |

Flag concordance (expected systems from curated `condition_registry` flags vs data-supported systems):

| observed_definition | definition_status | n_conditions | n_expected_condition_system_pairs | n_observed_condition_system_pairs | n_concordant | perm_null_mean_sum | perm_p | perm_p_holm_prespecified | per_system |
|---|---|---|---|---|---|---|---|---|---|
| any_source | pre-specified | 11 | 28 | 40 | 20 | 18.906 | 0.339 | 0.339 | {"immune": {"expected": 3, "observed": 11, "both": 3}, "autonomic_cardiac": {"expected": 6, "observed": 8, "both": 4}, "vascular_endothelial": {"expected": 5, "observed": 10, "both": 5}, "neuronal": {"expected": 7, "observed": 8, "both": 5}, "mitochondrial_metabolic": {"expected": 7, "observed": 3, "both": 3}} |
| condition_specific | pre-specified | 11 | 28 | 28 | 18 | 14.000 | 0.025 | 0.049 | {"immune": {"expected": 3, "observed": 6, "both": 3}, "autonomic_cardiac": {"expected": 6, "observed": 6, "both": 3}, "vascular_endothelial": {"expected": 5, "observed": 6, "both": 4}, "neuronal": {"expected": 7, "observed": 7, "both": 5}, "mitochondrial_metabolic": {"expected": 7, "observed": 3, "both": 3}} |
| non_literature_sources | exploratory (analyst, post hoc) | 11 | 28 | 23 | 14 | 12.642 | 0.292 | n/a | {"immune": {"expected": 3, "observed": 2, "both": 1}, "autonomic_cardiac": {"expected": 6, "observed": 8, "both": 4}, "vascular_endothelial": {"expected": 5, "observed": 3, "both": 2}, "neuronal": {"expected": 7, "observed": 7, "both": 4}, "mitochondrial_metabolic": {"expected": 7, "observed": 3, "both": 3}} |
| condition_specific_non_literature | exploratory (added in review, post hoc) | 11 | 28 | 20 | 13 | 10.912 | 0.159 | n/a | {"immune": {"expected": 3, "observed": 2, "both": 1}, "autonomic_cardiac": {"expected": 6, "observed": 6, "both": 3}, "vascular_endothelial": {"expected": 5, "observed": 3, "both": 2}, "neuronal": {"expected": 7, "observed": 6, "both": 4}, "mitochondrial_metabolic": {"expected": 7, "observed": 3, "both": 3}} |
| condition_specific_non_literature_unflagged | exploratory (added in review, post hoc) | 11 | 28 | 5 | 3 | 2.370 | 0.450 | n/a | {"immune": {"expected": 3, "observed": 1, "both": 1}, "autonomic_cardiac": {"expected": 6, "observed": 1, "both": 0}, "vascular_endothelial": {"expected": 5, "observed": 2, "both": 1}, "neuronal": {"expected": 7, "observed": 1, "both": 1}, "mitochondrial_metabolic": {"expected": 7, "observed": 0, "both": 0}} |

Concordant pairs and the sources behind them (`[post-hoc flag]` = drug-target gene family or single GWAS locus):

| observed_definition | concordant_pairs |
|---|---|
| any_source | dysautonomia:autonomic_cardiac(ot_literature+ot_other[post-hoc flag]); dysautonomia:vascular_endothelial(ot_literature); dysautonomia:neuronal(ot_literature); eds_hsd:autonomic_cardiac(ot_other); eds_hsd:vascular_endothelial(ot_literature+ot_other); endometriosis:mitochondrial_metabolic(ot_other[post-hoc flag]); fibromyalgia:neuronal(ot_literature+ot_other[post-hoc flag]); fibromyalgia:mitochondrial_metabolic(ot_other[post-hoc flag]); long_covid:immune(ot_literature); long_covid:vascular_endothelial(ot_literature); lyme_disease:immune(ot_literature+ot_other); me_cfs:immune(ot_literature); me_cfs:autonomic_cardiac(ot_other[post-hoc flag]); me_cfs:neuronal(ot_literature+ot_other[post-hoc flag]); me_cfs:mitochondrial_metabolic(ot_other[post-hoc flag]); migraine:vascular_endothelial(ot_literature+ot_other[post-hoc flag]); migraine:neuronal(ot_literature+ot_other); pots:autonomic_cardiac(ot_other[post-hoc flag]); pots:vascular_endothelial(ot_literature); pots:neuronal(ot_other[post-hoc flag]) |
| condition_specific | dysautonomia:autonomic_cardiac(ot_literature+ot_other[post-hoc flag]); dysautonomia:neuronal(ot_literature); eds_hsd:vascular_endothelial(ot_other); endometriosis:mitochondrial_metabolic(ot_other[post-hoc flag]); fibromyalgia:neuronal(ot_literature+ot_other[post-hoc flag]); fibromyalgia:mitochondrial_metabolic(ot_other[post-hoc flag]); long_covid:immune(ot_literature); long_covid:vascular_endothelial(ot_literature); lyme_disease:immune(ot_literature+ot_other); me_cfs:immune(ot_literature); me_cfs:autonomic_cardiac(ot_other[post-hoc flag]); me_cfs:neuronal(ot_other[post-hoc flag]); me_cfs:mitochondrial_metabolic(ot_other[post-hoc flag]); migraine:vascular_endothelial(ot_literature+ot_other[post-hoc flag]); migraine:neuronal(ot_literature+ot_other); pots:autonomic_cardiac(ot_other[post-hoc flag]); pots:vascular_endothelial(ot_literature); pots:neuronal(ot_other[post-hoc flag]) |
| non_literature_sources | dysautonomia:autonomic_cardiac(ot_other[post-hoc flag]); eds_hsd:autonomic_cardiac(ot_other); eds_hsd:vascular_endothelial(ot_other); endometriosis:mitochondrial_metabolic(ot_other[post-hoc flag]); fibromyalgia:neuronal(ot_other[post-hoc flag]); fibromyalgia:mitochondrial_metabolic(ot_other[post-hoc flag]); lyme_disease:immune(ot_other); me_cfs:autonomic_cardiac(ot_other[post-hoc flag]); me_cfs:neuronal(ot_other[post-hoc flag]); me_cfs:mitochondrial_metabolic(ot_other[post-hoc flag]); migraine:vascular_endothelial(ot_other[post-hoc flag]); migraine:neuronal(ot_other); pots:autonomic_cardiac(ot_other[post-hoc flag]); pots:neuronal(ot_other[post-hoc flag]) |
| condition_specific_non_literature | dysautonomia:autonomic_cardiac(ot_other[post-hoc flag]); eds_hsd:vascular_endothelial(ot_other); endometriosis:mitochondrial_metabolic(ot_other[post-hoc flag]); fibromyalgia:neuronal(ot_other[post-hoc flag]); fibromyalgia:mitochondrial_metabolic(ot_other[post-hoc flag]); lyme_disease:immune(ot_other); me_cfs:autonomic_cardiac(ot_other[post-hoc flag]); me_cfs:neuronal(ot_other[post-hoc flag]); me_cfs:mitochondrial_metabolic(ot_other[post-hoc flag]); migraine:vascular_endothelial(ot_other[post-hoc flag]); migraine:neuronal(ot_other); pots:autonomic_cardiac(ot_other[post-hoc flag]); pots:neuronal(ot_other[post-hoc flag]) |
| condition_specific_non_literature_unflagged | eds_hsd:vascular_endothelial(ot_other); lyme_disease:immune(ot_other); migraine:neuronal(ot_other) |

Source-level bias check: Open Targets literature sets support the immune system for 11 of 11 conditions tested.

## 5. Disagreement between sources (explicit)

Counts by kind: system_supported_by_some_sources_only: 47; overlap_not_condition_specific: 27; agreement: 19; presence_absence: 16; no_shared_genes: 13.

Presence/absence disagreements (one source names genes, the other none):

| condition_id | source_a | source_b | n_a | n_b | detail |
|---|---|---|---|---|---|
| long_covid | ot_genetic | gwas_catalog | 2 | 0 | Open Targets genetic_association targets name 2 protein-coding genes; GWAS Catalog mapped genes name none |
| long_covid | gwas_catalog | ot_literature | 0 | 300 | Open Targets literature (Europe PMC co-mention) targets name 300 protein-coding genes; GWAS Catalog mapped genes name none |
| long_covid | gwas_catalog | ot_other | 0 | 46 | Open Targets other datatypes (drug targets, animal models, expression, genetic literature, somatic, pathway) name 46 protein-coding genes; GWAS Catalog mapped genes name none |
| pots | ot_genetic | gwas_catalog | 1 | 0 | Open Targets genetic_association targets name 1 protein-coding genes; GWAS Catalog mapped genes name none |
| pots | gwas_catalog | ot_literature | 0 | 135 | Open Targets literature (Europe PMC co-mention) targets name 135 protein-coding genes; GWAS Catalog mapped genes name none |
| pots | gwas_catalog | ot_other | 0 | 101 | Open Targets other datatypes (drug targets, animal models, expression, genetic literature, somatic, pathway) name 101 protein-coding genes; GWAS Catalog mapped genes name none |
| fibromyalgia | ot_genetic | gwas_catalog | 0 | 13 | GWAS Catalog mapped genes name 13 protein-coding genes; Open Targets genetic_association targets name none |
| fibromyalgia | ot_genetic | ot_literature | 0 | 474 | Open Targets literature (Europe PMC co-mention) targets name 474 protein-coding genes; Open Targets genetic_association targets name none |
| fibromyalgia | ot_genetic | ot_other | 0 | 147 | Open Targets other datatypes (drug targets, animal models, expression, genetic literature, somatic, pathway) name 147 protein-coding genes; Open Targets genetic_association targets name none |
| eds_hsd | ot_genetic | gwas_catalog | 41 | 0 | Open Targets genetic_association targets name 41 protein-coding genes; GWAS Catalog mapped genes name none |
| eds_hsd | gwas_catalog | ot_literature | 0 | 595 | Open Targets literature (Europe PMC co-mention) targets name 595 protein-coding genes; GWAS Catalog mapped genes name none |
| eds_hsd | gwas_catalog | ot_other | 0 | 537 | Open Targets other datatypes (drug targets, animal models, expression, genetic literature, somatic, pathway) name 537 protein-coding genes; GWAS Catalog mapped genes name none |
| ptlds | ot_genetic | ot_literature | 0 | 3 | Open Targets literature (Europe PMC co-mention) targets name 3 protein-coding genes; Open Targets genetic_association targets name none |
| ptlds | ot_genetic | ot_other | 0 | 1 | Open Targets other datatypes (drug targets, animal models, expression, genetic literature, somatic, pathway) name 1 protein-coding genes; Open Targets genetic_association targets name none |
| ptlds | gwas_catalog | ot_literature | 0 | 3 | Open Targets literature (Europe PMC co-mention) targets name 3 protein-coding genes; GWAS Catalog mapped genes name none |
| ptlds | gwas_catalog | ot_other | 0 | 1 | Open Targets other datatypes (drug targets, animal models, expression, genetic literature, somatic, pathway) name 1 protein-coding genes; GWAS Catalog mapped genes name none |

Fibromyalgia: GWAS Catalog mapped genes name 13 protein-coding genes; Open Targets genetic_association targets name none. The GWAS Catalog audit lists possible reasons (release timing, Open Targets credible-set/L2G assignment vs the Catalog's positional mapping, sub-threshold Catalog associations); which one applies was not verified here.

Per-condition specificity: the same-condition fold is compared with the cross-condition pairs that share one side with it. With 7-22 cross pairs per condition the smallest attainable empirical p is 0.043-0.125, so a pair counts as condition-specific when its p is <= max(0.05, 1/(1 + n_cross)), i.e. when the same-condition fold is above every cross pair whenever 0.05 cannot be reached (review correction: the pre-review build used a fixed 0.05 threshold, which listed pairs such as migraine and IBS Open Targets genetic x GWAS Catalog as disagreements although their same-condition fold was the highest of all pairs).

Gene-level pairs with shared genes that are not condition-specific, or with no shared genes:

| condition_id | source_a | source_b | n_a | n_b | overlap | empirical_p_vs_cross | min_attainable_empirical_p | kind |
|---|---|---|---|---|---|---|---|---|
| long_covid | ot_genetic | ot_literature | 2 | 300 | 1 | 0.190 | 0.048 | overlap_not_condition_specific |
| long_covid | ot_genetic | ot_other | 2 | 46 | 0 | 0.333 | 0.048 | no_shared_genes |
| long_covid | ot_literature | ot_other | 300 | 46 | 5 | 0.261 | 0.043 | overlap_not_condition_specific |
| me_cfs | ot_genetic | ot_literature | 70 | 1587 | 5 | 1 | 0.048 | overlap_not_condition_specific |
| me_cfs | ot_genetic | ot_other | 70 | 109 | 1 | 0.381 | 0.048 | overlap_not_condition_specific |
| me_cfs | ot_genetic | mapmecfs | 70 | 1387 | 3 | 0.900 | 0.100 | overlap_not_condition_specific |
| me_cfs | gwas_catalog | ot_literature | 9 | 1587 | 0 | 0.947 | 0.053 | no_shared_genes |
| me_cfs | gwas_catalog | ot_other | 9 | 109 | 0 | 0.684 | 0.053 | no_shared_genes |
| me_cfs | gwas_catalog | mapmecfs | 9 | 1387 | 0 | 1 | 0.125 | no_shared_genes |
| me_cfs | ot_literature | ot_other | 1587 | 109 | 26 | 0.522 | 0.043 | overlap_not_condition_specific |
| me_cfs | ot_literature | mapmecfs | 1587 | 1387 | 138 | 0.250 | 0.083 | overlap_not_condition_specific |
| me_cfs | ot_other | mapmecfs | 109 | 1387 | 7 | 0.333 | 0.083 | overlap_not_condition_specific |
| pots | ot_genetic | ot_literature | 1 | 135 | 1 | 0.095 | 0.048 | overlap_not_condition_specific |
| pots | ot_genetic | ot_other | 1 | 101 | 1 | 0.095 | 0.048 | overlap_not_condition_specific |
| pots | ot_literature | ot_other | 135 | 101 | 11 | 0.130 | 0.043 | overlap_not_condition_specific |
| dysautonomia | ot_genetic | ot_literature | 28 | 313 | 4 | 0.095 | 0.048 | overlap_not_condition_specific |
| dysautonomia | gwas_catalog | ot_literature | 1 | 313 | 0 | 0.632 | 0.053 | no_shared_genes |
| dysautonomia | gwas_catalog | ot_other | 1 | 42 | 0 | 0.263 | 0.053 | no_shared_genes |
| fibromyalgia | gwas_catalog | ot_other | 13 | 147 | 0 | 0.842 | 0.053 | no_shared_genes |
| fibromyalgia | ot_literature | ot_other | 474 | 147 | 31 | 0.130 | 0.043 | overlap_not_condition_specific |
| lyme_disease | ot_genetic | ot_other | 11 | 635 | 0 | 1 | 0.048 | no_shared_genes |
| lyme_disease | gwas_catalog | ot_other | 4 | 635 | 0 | 0.895 | 0.053 | no_shared_genes |
| lyme_disease | ot_literature | ot_other | 436 | 635 | 59 | 0.087 | 0.043 | overlap_not_condition_specific |
| ptlds | ot_literature | ot_other | 3 | 1 | 0 | 0.261 | 0.043 | no_shared_genes |
| migraine | ot_genetic | ot_literature | 526 | 1770 | 173 | 0.095 | 0.048 | overlap_not_condition_specific |
| migraine | gwas_catalog | ot_other | 380 | 374 | 24 | 0.105 | 0.053 | overlap_not_condition_specific |
| ibs | ot_genetic | ot_literature | 175 | 1833 | 52 | 0.095 | 0.048 | overlap_not_condition_specific |
| ibs | ot_genetic | ot_other | 175 | 120 | 3 | 0.238 | 0.048 | overlap_not_condition_specific |
| ibs | gwas_catalog | ot_other | 188 | 120 | 0 | 0.947 | 0.053 | no_shared_genes |
| ibs | ot_literature | ot_other | 1833 | 120 | 66 | 0.304 | 0.043 | overlap_not_condition_specific |
| gastroparesis | ot_genetic | gwas_catalog | 47 | 9 | 0 | 0.706 | 0.059 | no_shared_genes |
| gastroparesis | ot_genetic | ot_literature | 47 | 168 | 1 | 0.333 | 0.048 | overlap_not_condition_specific |
| gastroparesis | ot_genetic | ot_other | 47 | 33 | 1 | 0.190 | 0.048 | overlap_not_condition_specific |
| gastroparesis | gwas_catalog | ot_literature | 9 | 168 | 1 | 0.263 | 0.053 | overlap_not_condition_specific |
| gastroparesis | gwas_catalog | ot_other | 9 | 33 | 0 | 0.474 | 0.053 | no_shared_genes |
| endometriosis | ot_genetic | ot_literature | 268 | 3016 | 119 | 0.143 | 0.048 | overlap_not_condition_specific |
| endometriosis | ot_genetic | ot_other | 268 | 419 | 13 | 0.143 | 0.048 | overlap_not_condition_specific |
| endometriosis | gwas_catalog | ot_literature | 248 | 3016 | 115 | 0.105 | 0.053 | overlap_not_condition_specific |
| endometriosis | gwas_catalog | ot_other | 248 | 419 | 13 | 0.368 | 0.053 | overlap_not_condition_specific |
| endometriosis | ot_literature | ot_other | 3016 | 419 | 171 | 0.348 | 0.043 | overlap_not_condition_specific |

Gene-level pairs whose same-condition overlap is above every cross-condition pair:

| condition_id | source_a | source_b | n_a | n_b | overlap | empirical_p_vs_cross | min_attainable_empirical_p |
|---|---|---|---|---|---|---|---|
| me_cfs | ot_genetic | gwas_catalog | 70 | 9 | 3 | 0.059 | 0.059 |
| dysautonomia | ot_genetic | gwas_catalog | 28 | 1 | 1 | 0.059 | 0.059 |
| dysautonomia | ot_genetic | ot_other | 28 | 42 | 5 | 0.048 | 0.048 |
| dysautonomia | ot_literature | ot_other | 313 | 42 | 29 | 0.043 | 0.043 |
| fibromyalgia | gwas_catalog | ot_literature | 13 | 474 | 3 | 0.053 | 0.053 |
| eds_hsd | ot_genetic | ot_literature | 41 | 595 | 34 | 0.048 | 0.048 |
| eds_hsd | ot_genetic | ot_other | 41 | 537 | 34 | 0.048 | 0.048 |
| eds_hsd | ot_literature | ot_other | 595 | 537 | 224 | 0.043 | 0.043 |
| lyme_disease | ot_genetic | gwas_catalog | 11 | 4 | 2 | 0.059 | 0.059 |
| lyme_disease | ot_genetic | ot_literature | 11 | 436 | 3 | 0.048 | 0.048 |
| lyme_disease | gwas_catalog | ot_literature | 4 | 436 | 2 | 0.053 | 0.053 |
| migraine | ot_genetic | gwas_catalog | 526 | 380 | 172 | 0.059 | 0.059 |
| migraine | ot_genetic | ot_other | 526 | 374 | 52 | 0.048 | 0.048 |
| migraine | gwas_catalog | ot_literature | 380 | 1770 | 151 | 0.053 | 0.053 |
| migraine | ot_literature | ot_other | 1770 | 374 | 291 | 0.043 | 0.043 |
| ibs | ot_genetic | gwas_catalog | 175 | 188 | 70 | 0.059 | 0.059 |
| ibs | gwas_catalog | ot_literature | 188 | 1833 | 60 | 0.053 | 0.053 |
| gastroparesis | ot_literature | ot_other | 168 | 33 | 12 | 0.043 | 0.043 |
| endometriosis | ot_genetic | gwas_catalog | 268 | 248 | 83 | 0.059 | 0.059 |

## 5b. Post-hoc diagnostics (NOT pre-specified; added after inspecting the results)

Two features of the evidence were noticed after the primary analysis ran: (i) positional GWAS mapping can place several genes of one locus into the same pathway; (ii) the Open Targets clinical (drug-target) datatype lists every target of a trialled drug mechanism, so one mechanism can contribute a whole gene family. The pre-specified support rule is unchanged; these diagnostics are reported beside it.

GWAS-supported systems and the number of independent GWAS loci (positions > 1 Mb apart) behind the genes of the system's enriched pathways:

| condition_id | physiological_system | overlap | n_independent_gwas_loci | genes_in_enriched_pathways | top_enriched_pathways |
|---|---|---|---|---|---|
| gastroparesis | immune | 4 | 1 | HLA-DPB1;HLA-DQA1;HLA-DQB1 | Phosphorylation of CD3 and TCR zeta chains; Translocation of ZAP-70 to Immunological synapse; Generation of second messenger molecules; Downstream TCR signaling; Interferon gamma signaling |

Open Targets 'other'-supported condition-system links: 25; in 15 of them at least half of the genes of the system's enriched pathways carry only drug-target (clinical) evidence among the non-literature, non-genetic datatypes:

| condition_id | physiological_system | overlap | share_overlap_drug_target_only | top_drug_mechanisms_behind_overlap | top_enriched_pathways |
|---|---|---|---|---|---|
| dysautonomia | autonomic_cardiac | 14 | 1 | n/a | Phase 0 - rapid depolarisation |
| fibromyalgia | autonomic_cardiac | 28 | 1 | ESLICARBAZEPINE ACETATE: Sodium channel alpha subunit blocker (10 of 18 genes); LACOSAMIDE: Sodium channel alpha subunit blocker (10 of 18 genes); LIDOCAINE: Sodium channel alpha subunit blocker (10 of 18 genes) | Phase 0 - rapid depolarisation; Cardiac conduction; Phase 2 - plateau phase |
| pots | neuronal | 35 | 1 | GABAPENTIN: Voltage-gated calcium channel modulator (20 of 34 genes); MEMANTINE: Glutamate [NMDA] receptor negative allosteric modulator (7 of 34 genes); MECAMYLAMINE: Neuronal acetylcholine receptor (2 of 34 genes) | Presynaptic depolarization and calcium channel opening; Transmission across Chemical Synapses; Neuronal System; NCAM1 interactions; NCAM signaling for neurite out-growth |
| me_cfs | mitochondrial_metabolic | 61 | 1 | METFORMIN: Mitochondrial complex I (NADH dehydrogenase) inhibitor (49 of 58 genes); PREGABALIN: Voltage-gated calcium channel modulator (7 of 58 genes); CLONIDINE: Adrenergic receptor alpha-2 agonist (2 of 58 genes) | Complex I biogenesis; Respiratory electron transport; Aerobic respiration and respiratory electron transport; Adrenaline,noradrenaline inhibits insulin secretion; Regulation of insulin secretion |
| me_cfs | autonomic_cardiac | 18 | 1 | PREGABALIN: Voltage-gated calcium channel modulator (8 of 9 genes); TILARGININE: Nitric oxide synthase inhibitor (1 of 9 genes) | Phase 2 - plateau phase; Phase 0 - rapid depolarisation; Cardiac conduction |
| me_cfs | neuronal | 32 | 1 | PREGABALIN: Voltage-gated calcium channel modulator (20 of 32 genes); KETAMINE: Glutamate [NMDA] receptor negative allosteric modulator (7 of 32 genes); TRIMETHAPHAN: Neuronal acetylcholine receptor (2 of 32 genes) | Presynaptic depolarization and calcium channel opening; Transmission across Chemical Synapses; Neuronal System; NCAM1 interactions; NCAM signaling for neurite out-growth |
| ibs | autonomic_cardiac | 21 | 1 | LIDOCAINE: Sodium channel alpha subunit blocker (10 of 18 genes); MEXILETINE: Sodium channel alpha subunit blocker (10 of 18 genes); PREGABALIN: Voltage-gated calcium channel modulator (8 of 18 genes) | Phase 0 - rapid depolarisation; Cardiac conduction; Phase 2 - plateau phase |
| gastroparesis | neuronal | 12 | 1 | MOSAPRIDE: Serotonin 3 (5-HT3) receptor antagonist (5 of 12 genes); SIMPINICLINE: Neuronal acetylcholine receptor (4 of 12 genes); BOTULINUM TOXIN TYPE A: Synaptosomal nerve-associated protein 25 (SNAP-25) hydrolytic enzyme (1 of 12 genes) | Transmission across Chemical Synapses; Neuronal System; Neurotransmitter receptors and postsynaptic signal transmission; Highly calcium permeable postsynaptic nicotinic acetylcholine receptors; Presynaptic nicotinic acetylcholine receptors |
| fibromyalgia | mitochondrial_metabolic | 64 | 1 | METFORMIN: Mitochondrial complex I (NADH dehydrogenase) inhibitor (49 of 58 genes); GABAPENTIN: Voltage-gated calcium channel modulator (7 of 58 genes); PREGABALIN: Voltage-gated calcium channel modulator (7 of 58 genes) | Complex I biogenesis; Respiratory electron transport; Aerobic respiration and respiratory electron transport; Adrenaline,noradrenaline inhibits insulin secretion; Regulation of insulin secretion |
| ibs | neuronal | 49 | 1 | PREGABALIN: Voltage-gated calcium channel modulator (20 of 48 genes); LIDOCAINE: Sodium channel alpha subunit blocker (10 of 48 genes); MEXILETINE: Sodium channel alpha subunit blocker (10 of 48 genes) | Transmission across Chemical Synapses; Presynaptic depolarization and calcium channel opening; Neuronal System; Neurotransmitter receptors and postsynaptic signal transmission; Interaction between L1 and Ankyrins |
| fibromyalgia | neuronal | 48 | 1 | GABAPENTIN: Voltage-gated calcium channel modulator (20 of 48 genes); PREGABALIN: Voltage-gated calcium channel modulator (20 of 48 genes); ESLICARBAZEPINE ACETATE: Sodium channel alpha subunit blocker (10 of 48 genes) | Presynaptic depolarization and calcium channel opening; Transmission across Chemical Synapses; Neuronal System; Interaction between L1 and Ankyrins; Neurotransmitter receptors and postsynaptic signal transmission |
| pots | autonomic_cardiac | 17 | 1 | GABAPENTIN: Voltage-gated calcium channel modulator (8 of 8 genes) | Phase 2 - plateau phase; Phase 0 - rapid depolarisation; Cardiac conduction |
| endometriosis | mitochondrial_metabolic | 86 | 0.98 | METFORMIN HYDROCHLORIDE: Mitochondrial complex I (NADH dehydrogenase) inhibitor (48 of 49 genes) | Complex I biogenesis; Respiratory electron transport; Aerobic respiration and respiratory electron transport |
| endometriosis | autonomic_cardiac | 12 | 0.91 | LIDOCAINE: Sodium channel alpha subunit blocker (10 of 11 genes) | Phase 0 - rapid depolarisation; Cardiac conduction |
| migraine | vascular_endothelial | 39 | 0.55 | NITRIC OXIDE: Soluble guanylate cyclase activator (3 of 11 genes); NITROGLYCERIN: Soluble guanylate cyclase activator (3 of 11 genes); EPOPROSTENOL: Prostanoid IP receptor agonist (1 of 11 genes) | Nitric oxide stimulates guanylate cyclase; Platelet homeostasis |
| migraine | autonomic_cardiac | 37 | 0.42 | LACOSAMIDE: Sodium channel alpha subunit blocker (10 of 26 genes); LAMOTRIGINE: Sodium channel alpha subunit blocker (10 of 26 genes); LIDOCAINE HYDROCHLORIDE: Sodium channel alpha subunit blocker (10 of 26 genes) | Phase 0 - rapid depolarisation; Cardiac conduction; Phase 2 - plateau phase |
| eds_hsd | autonomic_cardiac | 28 | 0.36 | LIDOCAINE: Sodium channel alpha subunit blocker (10 of 22 genes) | Cardiac conduction; Phase 0 - rapid depolarisation; Ion homeostasis |
| migraine | neuronal | 121 | 0.32 | TOPIRAMATE: Carbonic anhydrase II inhibitor (31 of 113 genes); GABAPENTIN: Voltage-gated calcium channel modulator (20 of 113 genes); BUTALBITAL: GABA-A receptor (12 of 113 genes) | Neuronal System; Transmission across Chemical Synapses; Neurotransmitter receptors and postsynaptic signal transmission; Presynaptic depolarization and calcium channel opening; Interaction between L1 and Ankyrins |
| eds_hsd | neuronal | 62 | 0.22 | LIDOCAINE: Sodium channel alpha subunit blocker (10 of 36 genes) | NCAM signaling for neurite out-growth; NCAM1 interactions; Interaction between L1 and Ankyrins; L1CAM interactions; EPHA-mediated growth cone collapse |
| lyme_disease | connective_tissue | 26 | 0.12 | DOXYCYCLINE: Bacterial 70S ribosome inhibitor (4 of 26 genes) | Activation of Matrix Metalloproteinases; Syndecan interactions; Extracellular matrix organization |
| eds_hsd | vascular_endothelial | 49 | 0.07 | CELIPROLOL: Adrenergic receptor alpha-2 antagonist (3 of 27 genes) | Platelet degranulation ; Response to elevated platelet cytosolic Ca2+; Platelet activation, signaling and aggregation; GP1b-IX-V activation signalling; Platelet Aggregation (Plug Formation) |
| eds_hsd | other | 378 | 0.05 | LIDOCAINE: Sodium channel alpha subunit blocker (10 of 257 genes); CELIPROLOL: Adrenergic receptor alpha-2 antagonist (5 of 257 genes); MECASERMIN: Insulin-like growth factor I receptor agonist (1 of 257 genes) | Signaling by PDGF; Muscle contraction; Diseases of signal transduction by growth factor receptors and second messengers; Gastrulation; Diseases associated with glycosaminoglycan metabolism |
| lyme_disease | vascular_endothelial | 57 | 0.02 | DOXYCYCLINE: Bacterial 70S ribosome inhibitor (1 of 42 genes) | Cell surface interactions at the vascular wall; Kinesins; Platelet degranulation ; Response to elevated platelet cytosolic Ca2+; Coagulation pathway |
| lyme_disease | immune | 192 | 0.01 | DOXYCYCLINE: Bacterial 70S ribosome inhibitor (2 of 150 genes) | Neutrophil degranulation; Interleukin-10 signaling; Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell; Regulation of Complement cascade; Complement cascade |
| eds_hsd | connective_tissue | 76 | 0 | n/a | Extracellular matrix organization; ECM proteoglycans; Collagen formation; Assembly of collagen fibrils and other multimeric structures; Degradation of the extracellular matrix |

`top_drug_mechanisms_behind_overlap` (review addition) names the Open Targets drug records whose ChEMBL mechanism targets account for most of those genes. Drug mechanisms of the same class that supply the largest share of genes to two or more flagged condition-system links:

| top_drug_mechanism | n_links | condition_system_links |
|---|---|---|
| Voltage-gated calcium channel modulator | 6 | fibromyalgia:neuronal; ibs:neuronal; me_cfs:autonomic_cardiac; me_cfs:neuronal; pots:autonomic_cardiac; pots:neuronal |
| Sodium channel alpha subunit blocker | 3 | endometriosis:autonomic_cardiac; fibromyalgia:autonomic_cardiac; ibs:autonomic_cardiac |
| Mitochondrial complex I (NADH dehydrogenase) inhibitor | 3 | endometriosis:mitochondrial_metabolic; fibromyalgia:mitochondrial_metabolic; me_cfs:mitochondrial_metabolic |

One drug mechanism of this kind can therefore supply the same physiological system to several conditions; this is clinical precedence (what has been trialled), not disease biology.

Condition-system links by pre-specified support tier (rows) and the tier after removing flagged sources (columns; exploratory):

| support_status | literature_only | multi_source | no_molecular_data | not_supported | single_source | too_few_genes |
|---|---|---|---|---|---|---|
| literature_only | 23 | 0 | 0 | 0 | 0 | 0 |
| multi_source | 7 | 7 | 0 | 0 | 0 | 0 |
| no_molecular_data | 0 | 0 | 12 | 0 | 0 | 0 |
| not_supported | 0 | 0 | 0 | 18 | 0 | 0 |
| single_source | 0 | 0 | 0 | 9 | 2 | 0 |
| too_few_genes | 0 | 0 | 0 | 0 | 0 | 6 |

Post-hoc sensitivity S8 (ot_other without the drug-target datatype): supported condition-system pairs 38 (primary 48); gene/pathway permutation p-values are in the sensitivity table below.

## 6. Sensitivity analyses (plan section 9; S8 is post hoc)

Permutation p (same vs cross condition) per source pair and analysis:

| level | source_a | source_b | S1_ot_indirect | S2_gwas_genome_wide | S3_gwas_single_trait | S4_project_background | S5_without_pots | S6_mapmecfs_with_sex_strata | S7_top_level_anchors_only | S8_ot_other_without_drug_targets | primary |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gene | gwas_catalog | ot_literature | 0.0093 | 0.0030 | 0.0041 | 0.0033 | 0.0063 | 0.0063 | 0.0063 | 0.0063 | 0.0063 |
| gene | gwas_catalog | ot_other | 0.7925 | 0.5107 | 0.5698 | 0.6516 | 0.6272 | 0.6272 | 0.6272 | 0.2583 | 0.6272 |
| gene | ot_genetic | ot_literature | 1.0e-04 | 1.0e-04 | 1.0e-04 | 5.0e-05 | 5.0e-05 | 1.0e-04 | 1.0e-04 | 1.0e-04 | 1.0e-04 |
| gene | ot_genetic | ot_other | 4.0e-04 | 1.0e-04 | 1.0e-04 | 1.0e-04 | 1.0e-04 | 1.0e-04 | 1.0e-04 | 0.0026 | 1.0e-04 |
| gene | ot_genetic | gwas_catalog | 2.0e-04 | 0.0014 | 2.0e-04 | 2.0e-04 | 2.0e-04 | 2.0e-04 | 2.0e-04 | 2.0e-04 | 2.0e-04 |
| gene | ot_literature | ot_other | 5.0e-05 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 2.0e-04 | 5.0e-05 |
| pathway | gwas_catalog | ot_literature | 0.3764 | 0.2500 | 0.1833 | 0.3958 | 0.3958 | 0.3958 | 0.3958 | 0.3958 | 0.3958 |
| pathway | gwas_catalog | ot_other | 0.7403 | 1.0000 | 0.8875 | 0.8361 | 0.8361 | 0.8361 | 0.8361 | 1.0000 | 0.8361 |
| pathway | ot_genetic | ot_literature | 0.0142 | 4.7e-04 | 4.7e-04 | 4.7e-04 | 4.7e-04 | 4.7e-04 | 4.7e-04 | 4.7e-04 | 4.7e-04 |
| pathway | ot_genetic | ot_other | 5.0e-04 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 0.0083 | 5.0e-05 |
| pathway | ot_genetic | gwas_catalog | 0.0167 | 0.1667 | 0.0167 | 0.0167 | 0.0167 | 0.0167 | 0.0167 | 0.0167 | 0.0167 |
| pathway | ot_literature | ot_other | 2.0e-04 | 5.0e-05 | 5.0e-05 | 5.0e-05 | 1.0e-04 | 5.0e-05 | 5.0e-05 | 0.0014 | 5.0e-05 |

Supported condition-system pairs (excluding `other`) and flag concordance per analysis:

| analysis | n_conditions | n_supported_condition_system_pairs | flag_concordance_perm_p_any_source |
|---|---|---|---|
| primary | 11 | 48 | 0.339 |
| S1_ot_indirect | 12 | 54 | 0.587 |
| S2_gwas_genome_wide | 11 | 48 | 0.339 |
| S3_gwas_single_trait | 11 | 48 | 0.339 |
| S4_project_background | 11 | 48 | 0.339 |
| S5_without_pots | 10 | 44 | 0.390 |
| S6_mapmecfs_with_sex_strata | 11 | 48 | 0.339 |
| S7_top_level_anchors_only | 11 | 38 | 0.060 |
| S8_ot_other_without_drug_targets | 11 | 38 | 0.447 |

Flag-concordance permutation p for every definition in every analysis (review addition; the first two definitions are pre-specified, the rest exploratory):

| definition | S1_ot_indirect | S2_gwas_genome_wide | S3_gwas_single_trait | S4_project_background | S5_without_pots | S6_mapmecfs_with_sex_strata | S7_top_level_anchors_only | S8_ot_other_without_drug_targets | primary |
|---|---|---|---|---|---|---|---|---|---|
| any_source | 0.587 | 0.339 | 0.339 | 0.339 | 0.390 | 0.339 | 0.060 | 0.447 | 0.339 |
| condition_specific | 0.126 | 0.017 | 0.025 | 0.025 | 0.029 | 0.025 | 0.007 | 0.107 | 0.025 |
| non_literature_sources | 0.314 | 0.229 | 0.292 | 0.292 | 0.290 | 0.173 | 0.044 | 0.591 | 0.292 |
| condition_specific_non_literature | 0.495 | 0.113 | 0.159 | 0.159 | 0.179 | 0.159 | 0.044 | 0.782 | 0.159 |
| condition_specific_non_literature_unflagged | 0.411 | 0.450 | 0.450 | 0.450 | 0.374 | 0.450 | 0.211 | 0.712 | 0.450 |

## 7. Measurable biology (`measurable_biology`)

247 rows: 224 condition x system x measurement-class rows (condition -> system data-derived, system -> class curated) and 23 study-catalogue rows (omics classes already applied to the condition).

Support status per condition and system:

| condition_id | autonomic_cardiac | connective_tissue | immune | mitochondrial_metabolic | neuronal | vascular_endothelial |
|---|---|---|---|---|---|---|
| dysautonomia | multi_source | not_supported | literature_only | not_supported | literature_only | literature_only |
| eds_hsd | single_source | multi_source | literature_only | not_supported | single_source | multi_source |
| endometriosis | single_source | literature_only | literature_only | single_source | not_supported | literature_only |
| fibromyalgia | single_source | literature_only | literature_only | single_source | multi_source | literature_only |
| gastroparesis | not_supported | not_supported | multi_source | not_supported | single_source | not_supported |
| ibs | multi_source | literature_only | literature_only | not_supported | multi_source | literature_only |
| long_covid | not_supported | literature_only | literature_only | not_supported | not_supported | literature_only |
| lyme_disease | not_supported | multi_source | multi_source | not_supported | not_supported | multi_source |
| mcas | no_molecular_data | no_molecular_data | no_molecular_data | no_molecular_data | no_molecular_data | no_molecular_data |
| me_cfs | single_source | literature_only | literature_only | single_source | multi_source | literature_only |
| migraine | multi_source | literature_only | literature_only | not_supported | multi_source | multi_source |
| post_infectious_syndrome | no_molecular_data | no_molecular_data | no_molecular_data | no_molecular_data | no_molecular_data | no_molecular_data |
| pots | single_source | not_supported | literature_only | not_supported | single_source | literature_only |
| ptlds | too_few_genes | too_few_genes | too_few_genes | too_few_genes | too_few_genes | too_few_genes |

Demo cluster candidate measurement classes with molecular support:

| condition_id | physiological_system | measurement_class | support_status | sources_supporting | condition_specific | best_system_fdr |
|---|---|---|---|---|---|---|
| long_covid | immune | immune_assays | literature_only | ot_literature | True | 4.8e-38 |
| long_covid | immune | proteomics | literature_only | ot_literature | True | 4.8e-38 |
| long_covid | immune | blood_biomarkers | literature_only | ot_literature | True | 4.8e-38 |
| long_covid | vascular_endothelial | endothelial_function | literature_only | ot_literature | True | 1.0e-08 |
| long_covid | vascular_endothelial | microvascular_function | literature_only | ot_literature | True | 1.0e-08 |
| long_covid | vascular_endothelial | capillaroscopy | literature_only | ot_literature | True | 1.0e-08 |
| long_covid | vascular_endothelial | blood_biomarkers | literature_only | ot_literature | True | 1.0e-08 |
| long_covid | connective_tissue | vascular_imaging | literature_only | ot_literature | True | 0.016 |
| me_cfs | immune | immune_assays | literature_only | ot_literature | True | 6.2e-53 |
| me_cfs | immune | proteomics | literature_only | ot_literature | True | 6.2e-53 |
| me_cfs | immune | blood_biomarkers | literature_only | ot_literature | True | 6.2e-53 |
| me_cfs | vascular_endothelial | endothelial_function | literature_only | ot_literature | True | 1.6e-25 |
| me_cfs | vascular_endothelial | microvascular_function | literature_only | ot_literature | True | 1.6e-25 |
| me_cfs | vascular_endothelial | capillaroscopy | literature_only | ot_literature | True | 1.6e-25 |
| me_cfs | vascular_endothelial | blood_biomarkers | literature_only | ot_literature | True | 1.6e-25 |
| me_cfs | autonomic_cardiac | hrv | single_source | ot_other | True | 7.4e-14 |
| me_cfs | autonomic_cardiac | ecg_ambulatory | single_source | ot_other | True | 7.4e-14 |
| me_cfs | autonomic_cardiac | autonomic_testing | single_source | ot_other | True | 7.4e-14 |
| me_cfs | autonomic_cardiac | wearable_heart_rate | single_source | ot_other | True | 7.4e-14 |
| me_cfs | mitochondrial_metabolic | cpet | single_source | ot_other | True | 5.1e-16 |
| me_cfs | mitochondrial_metabolic | metabolomics | single_source | ot_other | True | 5.1e-16 |
| me_cfs | neuronal | small_fiber_testing | multi_source | ot_literature;ot_other | True | 1.1e-10 |
| me_cfs | neuronal | digital_cognitive | multi_source | ot_literature;ot_other | True | 1.1e-10 |
| me_cfs | connective_tissue | vascular_imaging | literature_only | ot_literature | True | 2.1e-13 |
| pots | immune | immune_assays | literature_only | ot_literature | False | 0.049 |
| pots | immune | proteomics | literature_only | ot_literature | False | 0.049 |
| pots | immune | blood_biomarkers | literature_only | ot_literature | False | 0.049 |
| pots | vascular_endothelial | endothelial_function | literature_only | ot_literature | True | 0.004 |
| pots | vascular_endothelial | microvascular_function | literature_only | ot_literature | True | 0.004 |
| pots | vascular_endothelial | capillaroscopy | literature_only | ot_literature | True | 0.004 |
| pots | vascular_endothelial | blood_biomarkers | literature_only | ot_literature | True | 0.004 |
| pots | autonomic_cardiac | hrv | single_source | ot_other | True | 2.7e-13 |
| pots | autonomic_cardiac | ecg_ambulatory | single_source | ot_other | True | 2.7e-13 |
| pots | autonomic_cardiac | autonomic_testing | single_source | ot_other | True | 2.7e-13 |
| pots | autonomic_cardiac | wearable_heart_rate | single_source | ot_other | True | 2.7e-13 |
| pots | neuronal | small_fiber_testing | single_source | ot_other | True | 4.0e-14 |
| pots | neuronal | digital_cognitive | single_source | ot_other | True | 4.0e-14 |

## 8. Figure 5 data

`results/tables/fig5_nodes.csv` (238 nodes: gene 168, pathway 37, measurement_class 15, condition 11, physiological_system 7) and `results/tables/fig5_edges.csv` (475 edges: condition_gene 193, gene_pathway 102, condition_pathway 79, condition_system 48, pathway_system 37, system_measurement 16). Draft figures: `results/figures/drafts/fig5_system_heatmap.png`, `results/figures/drafts/fig5_system_heatmap.svg`, `results/figures/drafts/fig5_graph_demo.png`, `results/figures/drafts/fig5_graph_demo.svg`.

## 9. Limitations

* Open Targets genetic evidence includes GWAS credible sets built from studies that are also in the GWAS Catalog, so their agreement is partly shared data, not independent replication.
* Open Targets literature evidence is Europe PMC co-mention; it tracks research attention and is biased towards well-studied (often immune) genes. It is kept as its own source and flagged `literature_only`.
* Agreement between the other source pairs is not independent replication either. Open Targets `genetic_association` (Orphanet, ClinVar/EVA, gene burden, GWAS credible sets) and `genetic_literature` (Gene2Phenotype, ClinGen, PanelApp, UniProt literature, counted in 'other') curate overlapping monogenic gene-disease knowledge. Europe PMC co-mention can be mined from the same GWAS and drug-trial papers that the genetic and drug evidence comes from; this was not measurable here, because only the top evidence items were retrieved (3 of 667 retrieved literature items cite a GWAS Catalog paper of the same condition).
* Open Targets 'other' evidence for long_covid and me_cfs is the clinical (drug-target) datatype only, and for pots almost only that. ChEMBL mechanisms list whole target families, so a single trialled drug (e.g. a gabapentinoid, metformin, a sodium-channel blocker) can make a physiological system look enriched, and the same drug record then does so for every condition it was trialled in.
* Confidence intervals are percentile bootstraps over 5-12 conditions and are approximate (likely too narrow); exact permutation p-values with 5-8 conditions have a floor of 1/n! (e.g. 0.0083 with 5).
* 'Condition-specific' in the system tables means only that the system's log2 fold is above the median of the other conditions for the same source; about half of all tests pass by construction.
* GWAS Catalog mapped genes are positional (nearest/overlapping genes) and multi-trait analyses attach the same locus to two conditions; sensitivity analyses S2/S3 restrict them.
* mapMECFS gene and protein lists are nominal p < 0.05 results from one small cohort with no FDR-significant row; they are included as published, with that label.
* POTS evidence comes from the NET-deficiency id only; MCAS has no molecular evidence in any source; post-infectious syndrome is a grouping node and PTLDS has too few genes.
* Physiological-system assignment uses Reactome's hierarchy plus curated sub-pathway anchors; Reactome has no autonomic top level, so autonomic/cardiac support depends entirely on the curated anchors (S7 shows the top-level-only result).
* System -> measurement-class links are curated assumptions with a stated rationale; they are not evidence that a measurement detects the condition.
* Enrichment is condition-level: it describes the biology public databases associate with a disease label, not the biology of any person, and it is not evidence of mechanism.

## Reproduce

`uv run python -m measure_it.omics.coherence && uv run python -m measure_it.omics.graph`
