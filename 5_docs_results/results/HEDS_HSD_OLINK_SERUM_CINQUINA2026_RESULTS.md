# Serum Olink proteome in hEDS/HSD: pre-specified case-control analysis

Dataset `heds_hsd_olink_serum_cinquina2026`: Cinquina et al., Clinical Proteomics 2026 (PMC13081554), Additional files 2 and 5. 460 Olink
Target 96 assays in serum from 88 hEDS, 88 HSD (2017 criteria) and 176 healthy controls. Plan: `docs/ANALYSIS_PLAN_HEDS_HSD_OLINK_SERUM_CINQUINA2026.md` (written
before any protein value was compared; amendment 1 added before computation). Analysis `heds_hsd_olink_serum_cinquina2026 2026-09-24.1`,
computed 2026-09-25T22:19:20+00:00.

## Answer first

**Primary (Italian, site-matched): 88 Italian hEDS/HSD vs 176 Italian controls,
458 assays, L2 logistic, 5-fold x 20 CV: AUROC 0.684 (0.624 to 0.747); label-permutation p <= 0.0010 (no shuffle reached the observed value)
(null mean 0.499, 95th percentile 0.571, 1000 shuffles).**
Decision by the locked rule: **the serum proteome discriminates Italian hEDS/HSD from Italian controls (pre-specified rule met)**.

* Sensitivity at 90% specificity: 0.30 (0.20-0.41).
* Per protein (S1, Italian sample): **24 of 458 assays have BH q < 0.05**.
* **Site negative control (S2):** American vs Italian *patients* (same diagnoses) are separated with AUROC 0.514 (0.441 to 0.586),
  permutation p = 0.4006. No site signal was detectable, so pooling sites is less of a concern than the design suggested.
* Whole cohort as in the paper (S3, site-confounded by design): AUROC 0.730 (0.684 to 0.776).
* hEDS vs HSD (S4, Italian): AUROC 0.480 (0.373 to 0.589), p = 0.6124 (paper: no DEPs between them).
* Probands only (S5, relatives excluded): AUROC 0.712 (0.636 to 0.786); plate indicators added (S7): AUROC 0.681 (0.621 to 0.744).
* 24 of the 24 proteins at q < 0.05 in the Italian site-matched sample are
  in the paper's list of 69 pooled DEPs. The strongest are among the paper's top-ranked proteins (MYOC and COMP lower;
  NUDT5 higher). The other 45 of the paper's DEPs do not reach
  q < 0.05 here. With half the patients this is partly a matter of power, so whether they depend on the American
  samples cannot be settled.
* **Size:** an AUROC of about 0.68, with 30% sensitivity at 90% specificity, is a modest group-level separation and
  not a usable test. The paper fitted its LASSO, random forest and classification tree on the full data and reports
  no held-out estimate. This one is modest.
* **Age is not controlled.** Controls are older on average (paper: 42 vs 36-38 y), and age is unavailable per control.
  Some of the separation could be age-related.

What it means: a serum protein profile distinguishes Italian clinic patients with hEDS/HSD from Italian healthy volunteers in held-out people. This is a group-level, case-control result: the controls are healthy volunteers (older on average, age/sex not released), not people with other causes of joint pain, fatigue or dysautonomia, so it is not evidence of a diagnostic test.

## Published claim vs our own computation

| | authors (published) | this analysis (from the released rows) |
|---|---|---|
| Sample | 88 hEDS + 88 HSD (44 Italian + 44 American each) vs 176 Italian controls, pooled | primary: Italian patients only (88) vs 176 controls; pooled cohort only as S3 |
| Method | Wilcoxon + BH per protein; LASSO -> RF / classification tree on the 69 DEPs, fitted on all data | pre-specified CV AUROC, person bootstrap, label-permutation null; per-protein Wilcoxon + BH |
| Per-protein result | 69 pooled DEPs (54 hEDS, 49 HSD, 0 hEDS vs HSD) | 24 Italian-only assays q < 0.05; overlap with the paper's pooled list: 24 of 69 |
| Site | "no DEPs" Italian vs American within hEDS or HSD | S2 CV AUROC 0.514 (0.441 to 0.586) for American vs Italian patients |

## Model results

| analysis | cases / controls | AUROC (95% CI) | permutation p | sensitivity at 90% specificity |
|---|---|---|---|---|
| primary: Italian hEDS+HSD vs controls | 88 / 176 | 0.684 (0.624 to 0.747) | <= 0.0010 (no shuffle reached the observed value) | 0.30 (0.20-0.41) |
| S2 site negative control: US vs Italian patients | 88 / 88 | 0.514 (0.441 to 0.586) | 0.4006 | 0.06 (0.03-0.12) |
| S3 whole cohort all patients vs controls (site-confounded by design) | 176 / 176 | 0.730 (0.684 to 0.776) | not run | 0.31 (0.20-0.42) |
| S4 hEDS vs HSD (Italian) | 44 / 44 | 0.480 (0.373 to 0.589) | 0.6124 | 0.12 (0.04-0.26) |
| S5 primary, probands only (relatives excluded) | 63 / 176 | 0.712 (0.636 to 0.786) | not run | 0.36 (0.23-0.48) |
| S7 primary + plate indicators | 88 / 176 | 0.681 (0.621 to 0.744) | not run | 0.29 (0.20-0.41) |

## Top 15 proteins (S1, Italian sample, sorted by p)

| protein | panel | n | Hedges g (95% CI) | P(case > control) | Mann-Whitney p | BH q |
|---|---|---|---|---|---|---|
| MYOC | DEV | 86 / 176 | -0.75 (-1.04 to -0.47) | 0.300 | 1.43e-07 | 0.0001 |
| CA14 | ODA | 86 / 171 | -0.61 (-0.90 to -0.36) | 0.320 | 2.61e-06 | 0.0006 |
| NUDT5 | DEV | 86 / 176 | 0.46 (0.19 to 0.75) | 0.655 | 4.65e-05 | 0.0057 |
| MAX | ODA | 86 / 171 | 0.42 (0.17 to 0.71) | 0.655 | 4.96e-05 | 0.0057 |
| FOXO1 | ODA | 86 / 171 | 0.49 (0.22 to 0.78) | 0.648 | 1.05e-04 | 0.0096 |
| RGMA | NEU | 87 / 172 | -0.50 (-0.75 to -0.26) | 0.358 | 1.82e-04 | 0.0139 |
| INPPL1 | ODA | 86 / 171 | 0.51 (0.23 to 0.82) | 0.641 | 2.34e-04 | 0.0153 |
| COMP | CAM | 88 / 176 | -0.47 (-0.75 to -0.21) | 0.367 | 4.13e-04 | 0.0217 |
| MIF | DEV | 86 / 176 | 0.46 (0.22 to 0.73) | 0.634 | 4.26e-04 | 0.0217 |
| BLVRB | DEV | 86 / 176 | 0.32 (0.08 to 0.63) | 0.629 | 6.75e-04 | 0.0238 |
| NUB1 | ODA | 86 / 171 | 0.39 (0.14 to 0.66) | 0.630 | 6.84e-04 | 0.0238 |
| MMP-10 | INF | 87 / 172 | -0.44 (-0.72 to -0.18) | 0.371 | 7.16e-04 | 0.0238 |
| CLEC10A | NEU | 87 / 172 | -0.46 (-0.73 to -0.20) | 0.372 | 7.58e-04 | 0.0238 |
| MAGED1 | ODA | 86 / 171 | 0.29 (0.05 to 0.56) | 0.629 | 7.74e-04 | 0.0238 |
| CA2 | DEV | 86 / 176 | 0.41 (0.12 to 0.70) | 0.628 | 7.78e-04 | 0.0238 |

Overlap with the paper's pooled DEP list (post hoc, as planned): ours only = ;
paper only (present in our assays) = 4E-BP1;AIFM1;AXIN1;BANK1;BMP-4;C1QTNF1;CA3;CAPG;CCL5;CD5;CD6;CD69;EN-RAGE;ERBB2IP;FCER2;FGR;GDF-8;HAGH;HAVCR2;IL7;LAT2;LXN;MAEA;MAP4K5;MCP-2;MESDC2;NTRK2;PARK7;PDGF-R-alpha;PEBP1;PLXDC1;PRKRA;PTPRJ;PXN;RASSF2;RRM2B;SIRT2;SMAD1;SOD1;SPOCK1;STAMBP;STXBP3;TIGAR;YES1;sFRP-3.

## Caveats carried with every number

Case-control; controls are healthy Italian volunteers (not people with other causes of joint pain), older on average than patients (paper: mean 42 vs 36-38 y) with no per-person age/sex released for controls; all controls are Italian while half the patients are American home draws shipped to Baltimore, so only the Italian-only contrast is site-matched; Olink NPX are relative log2 units. Bootstrap CIs hold the out-of-fold predictions fixed (person sampling only). Relatives of probands are
included in the primary (S5 excludes them).

## Outputs

`participants__heds_hsd_olink_serum_cinquina2026`, `participant_labs__heds_hsd_olink_serum_cinquina2026`, `phenotype_signatures__heds_hsd_olink_serum_cinquina2026` (459
rows); `results/tables/heds_hsd_olink_serum_cinquina2026_{model_performance, primary_decision, permutation_null, protein_effects,
paper_dep_overlap}.csv`, `heds_hsd_olink_serum_cinquina2026_run_metadata.json`; figure results/figures/drafts/heds_hsd_olink_serum_cinquina2026_auroc_and_null.png.

## Reproduce

```bash
uv run python -m measure_it.labs.heds_hsd_olink_serum_cinquina2026 --jobs 16      # deterministic, seed 20260923
uv run pytest tests/test_heds_hsd_olink_serum_cinquina2026.py
```
