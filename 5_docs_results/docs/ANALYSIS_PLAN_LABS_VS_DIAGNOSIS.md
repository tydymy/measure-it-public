# Analysis plan — do clinical labs / non-wearable measurements separate invisible-illness labels?

**Locked:** 2026-09-24, before any lab-by-label or metabolite-by-group statistic was computed. The file is written by
`measure_it.labs.plan.lock_plan()`; the analysis modules (`measure_it.labs.labs_vs_diagnosis`,
`measure_it.labs.appelman_metabolomics`) refuse to run if this file is missing or has been edited (sha256 check against
the plan text in the module; lock record `results/tables/labs_dx_plan_lock.json`). Anything added after the lock is
labelled **post hoc** in the results. Nulls are reported in the same tables, and with the same prominence, as positives.

**Question (from the project lead):** do clinical laboratory measurements, or other non-wearable measurements, separate
diseases under the invisible-illness umbrella from controls, and ideally from each other?

## What had been seen before this plan was written

* **NHANES 2011-2014.** The structure of the lab tables (`participant_labs__nhanes`, `participant_labs_extended__nhanes`:
  variables, layers, subsample weights, detection-limit flags, counts per layer) and label counts (number answering
  yes/no to each questionnaire item; rx reason-code group sizes) among adults aged >= 20 with the biochemistry profile
  (n = 10,106). The project's earlier wearable analysis used 32 labs inside a "clinical" feature block
  (results/WEARABLE_NHANES_RESULTS.md) but never reported a lab-by-label contrast; no lab-by-label statistic had been
  computed or printed in this work.
* **Appelman 2024 Source Data.** Sheet names, column names, `Group x Time` counts, that metabolite values are all
  positive with no missing cells at baseline, and the pathway-domain sheets. While printing sheet headers, the first
  1-3 data rows of `Supplementary Data Figure 6A/6B` (participants 1 and 2, both Long COVID; no control rows) were
  displayed; no group summary or comparison was seen. The paper's text (Results, Fig. 3 legend, Methods; Europe PMC
  PMC10766651) was read to record the published claims below **before** our computation.

## Part 1 — NHANES 2011-2014 labs vs diagnosis-like labels

### Population

NHANES 2011-2012 (G) and 2013-2014 (H) participants aged >= 20 at screening, not pregnant at the MEC exam, with the
standard biochemistry profile measured (`lab_albumin_g_dl` non-missing). Joined on SEQN only
(`participant_id = nhanes:<SEQN>`). The wearable-valid restriction used elsewhere is **not** applied (labs do not need
a wearable). Analyses are unweighted: internal comparisons, not US population estimates. Each label has its own
population (people for whom the label item was answered; rx labels: cycle H only, because reason-for-use codes exist
only in 2013-2014).

### Labels (fixed list; case = 1, control = 0)

| group | label_id | definition | controls |
|---|---|---|---|
| umbrella | `rx_fibromyalgia` | any prescription with reason-for-use ICD-10-CM M79.7 (H only) | other H participants in the population |
| umbrella | `rx_migraine` | G43 | same |
| umbrella | `rx_ibs` | K58 | same |
| umbrella | `rx_insomnia` | G47.0 / F51.0 | same |
| umbrella | `rx_myalgia` | M79.1 (statin-associated muscle pain is a likely contributor) | same |
| umbrella | `mecfs_like_proxy` | **PROXY, not ME/CFS**: DPQ040 >= 2 AND PFQ functional-limitation composite AND no affirmative answer to any exclusion diagnosis (`nhanes_cohort.PROXY_EXCLUSION_ITEMS`) | everyone else with DPQ040 and PFQ answered |
| umbrella | `fatigue_symptom` | DPQ040 >= 2 (symptom proxy) | DPQ040 answered |
| umbrella | `told_sleep_disorder` | SLQ060 = yes | SLQ060 = no |
| umbrella | `depression_phq9` | PHQ-9 >= 10 (screen) | PHQ-9 complete |
| comparator | `celiac_disease` | MCQ082 = yes | MCQ082 = no |
| comparator | `thyroid_problem` | MCQ160M = yes | = no |
| comparator | `rheumatoid_arthritis` | MCQ195 = rheumatoid | MCQ160A = no (no arthritis) |
| comparator | `osteoarthritis` | MCQ195 = osteoarthritis | MCQ160A = no |
| comparator | `psoriasis` | MCQ070 = yes | = no |
| comparator | `asthma` | MCQ010 = yes | = no |
| positive control | `diabetes_no_glycaemic_labs` | DIQ010 = yes; **glycaemic analytes removed** (HbA1c, serum glucose, fasting glucose, insulin, OGTT glucose) | DIQ010 = no (borderline excluded) |
| positive control (trivial) | `diabetes_with_glycaemic_labs` | same label, all analytes (sanity check; expected near-perfect) | same |
| positive control | `weak_failing_kidneys` | KIQ022 = yes (KIQ_U_G/H, downloaded for this analysis) | KIQ022 = no |
| positive control | `gout` | MCQ160N = yes | = no |
| positive control | `congestive_heart_failure` | MCQ160B = yes | = no |
| positive control | `anemia_treatment` | MCQ053 = yes (treated for anaemia in past 3 months) | = no |
| positive control | `liver_condition` | MCQ160L = yes | = no |

No analyte is removed for any other label (self-reported labels are not defined by a lab). The ME/CFS-like proxy's
exclusion list removes people reporting diabetes, thyroid problems, anaemia treatment etc., so labs linked to those
conditions are partly "defined out" of proxy cases; this is stated next to its results.

**Pairwise umbrella contrasts** (cases of one label vs cases of the other; people in both groups dropped):
`rx_fibromyalgia` vs `rx_migraine`; `rx_fibromyalgia` vs `rx_insomnia`; `rx_migraine` vs `rx_insomnia`;
`rx_fibromyalgia` vs `rheumatoid_arthritis` (H only); `mecfs_like_proxy` vs PHQ-9 >= 10 non-proxy.

### Analytes

* **Primary panel (tier `primary`)**: continuous analytes from full-sample layers measured in both cycles: CBC (20),
  biochemistry profile (`LBXSCH` dropped as a duplicate of total cholesterol), HbA1c, total and HDL cholesterol,
  25-OH vitamin D, vitamin B12 (G `LBXB12` and H Deming-adjusted `LBDB12` pooled), and the continuous analytes of the
  full-sample extended layers `urine_albumin_creatinine`, `nutrition_folate_b12_status`, `hormones_sex_steroids`,
  `tobacco_cotinine`. An analyte stays in the panel only if, in the population, it is >= 70% non-missing and < 50%
  below the detection limit. Urine flow-rate process variables are not analytes.
* **Extended tier (`extended`)**: every other continuous analyte in `participant_labs_extended__nhanes` (subsample
  layers analysed only among sampled people; special smoking-sample rows dropped; analytes >= 50% below LOD dropped)
  plus the subsample analytes of `participant_labs__nhanes` (fasting glucose, insulin, fasting triglycerides/LDL,
  2011-2012 thyroid panel). Binary results (serology, STI, HPV, celiac tTG/EMA, TB IGRA) are tested separately as
  `result_positive` (tier `extended_binary`). Extended-tier analytes enter the per-analyte screen only (no models):
  their subsample sizes differ per analyte.
* **Transforms (stateless, label-blind):** for an analyte with minimum >= 0 and skewness > 1 in the population,
  `log(x + c)` with `c` = half the smallest positive value (c = 0 when the minimum is > 0). Below-LOD values stay as
  released (NCHS fill LLOD/sqrt(2)).
* Urinary concentration analytes (component description contains "Urine"; not ratios, not creatinine itself) get log
  urinary creatinine (`URXUCR`, ALB_CR) as an additional covariate.

### Per-analyte contrasts

For each label x analyte with >= 10 measured cases and >= 10 measured controls: OLS of the z-scored (transformed)
analyte on case + (age - 50) + (age - 50)^2 + female, HC3 standard errors (the project's
`signatures.adjusted_std_difference`); effect = age/sex-adjusted standardised difference (SD units), 95% CI, p.
Sensitivity: the same + BMI (reported next to it, not used for the decision). Binary analytes: logistic regression on
the same covariates, effect = adjusted log odds ratio (>= 5 positives among cases and >= 5 among controls required).
**BH-FDR within label x tier.** Significance at q < 0.05. `underpowered` = fewer than 50 cases.

### Models (primary panel only; labels with >= 25 cases)

* Feature sets: `demographics` (age, female, 5 race/ethnicity indicators, income-poverty ratio, education ordinal;
  the project's `DEMOGRAPHIC_FEATURES`); `demographics_bmi` (+ BMI); `labs` (primary-panel analytes only);
  `demographics_labs`; `demographics_bmi_labs`.
* Model: L2 logistic regression, C chosen from {0.001, 0.01, 0.1, 1, 10} by inner 3-fold CV (log-loss) **inside each
  training fold**; median imputation with missingness indicators and standardisation fitted on the training fold only.
* CV: stratified 5-fold x 5 repeats (seed 20260923 + repeat), identical folds across feature sets.
* Metrics: AUROC (primary), AUPRC. Point estimate = mean over repeats; 95% CI = participant bootstrap (B = 1,000
  multinomial weights applied to the out-of-fold predictions of every repeat, averaged over repeats); deltas paired on
  the same resamples; two-sided bootstrap p for each delta.
* **Primary estimand:** Delta AUROC (`demographics_labs` - `demographics`). Also reported: labs-only AUROC and
  `demographics_bmi_labs` - `demographics_bmi`.
* **Label-permutation null** (labs-only model): 200 permutations, one CV repeat each, C fixed at the median C chosen in
  the observed labs-only run (fixed to keep the null affordable; a stated deviation from full re-tuning);
  p = (1 + #null >= observed) / 201.
* **Lab-block shuffle** (increment control): lab rows permuted across participants (demographics intact), 20 shuffles,
  repeat-1 folds; p = (1 + #shuffled Delta >= observed) / 21.
* **Decision rule:** "labs add to demographics" = Delta AUROC 95% CI excludes 0 **and** lab-block shuffle p <= 0.05.
  "labs alone discriminate" = labs-only AUROC 95% CI lower bound > 0.5 **and** permutation p <= 0.05 (labs-only models
  can learn age/sex from sex hormones and age-related analytes, so the increment is the primary reading). A Delta
  below +0.02 is described as negligible even when its CI excludes 0. BH across the umbrella labels' Delta bootstrap
  p-values is reported as a descriptive multiplicity check.
* **Positive-control gate:** the pipeline is judged able to detect lab signal if `diabetes_no_glycaemic_labs`,
  `weak_failing_kidneys` and `gout` all meet "labs add to demographics". If any fails, the nulls of this part are
  reported as uninterpretable for that reason.
* Pairwise contrasts: per-analyte contrasts as above; models (demographics vs demographics + labs) if both groups
  have >= 25 people.

### Interpretation fixed in advance

Rx reason codes identify **people taking a prescription for** the condition (treated cases), so lab differences can be
drug effects (statins and CK/lipids in `rx_myalgia`; antiepileptics, triptans, antidepressants, gabapentinoids in the
other rx groups). Self-reported "ever told" diagnoses have unknown recency and severity. Controls are the general
population, not healthy volunteers. NHANES 2011-2014 has no CRP/hs-CRP, ferritin or ESR. A per-analyte difference is
a group difference; only the cross-validated models speak to individual-level separation.

## Part 2 — MUSCLE-ME / Appelman 2024 Source Data (blood and muscle metabolomics)

### Data

`41467_2023_44432_MOESM4_ESM.xlsx` (Appelman et al., Nat Commun 2024;15:17, doi:10.1038/s41467-023-44432-3,
CC BY 4.0), registered as source `appelman_lc_pem_source`. Sheets `Supplementary Data Figure 6B` (venous blood /
plasma, 83 metabolites) and `6A` (vastus lateralis muscle, 116 metabolites); groups `Long COVID` and `Healthy`.
**No ME/CFS group exists in this file** (ME/CFS participants appear only in Charlton 2026, which released no
metabolomics), so an ME/CFS contrast is not possible. Sex and BMI from sheet `Table 1`; **age is not released** (the
sheet says so). Participant id `appelman_lc_pem_source:<ID>`. Same trial (NCT05225688) as Charlton 2026: not
independent evidence.

Timepoint for the primary analyses: **Baseline**. Blood: 25 Long COVID vs 21 healthy; muscle: 25 vs 19 (the authors
excluded 4 muscle samples at QC). Metabolite values are used as released (normalised by the authors to tissue weight /
plasma volume and internal standards), log2-transformed (all values > 0).

### Primary analyses (one per tissue)

* **Model:** in each training fold: standardise, keep the 10 metabolites with the largest |Welch t| **computed on the
  training fold only**, L2 logistic regression (C = 1). Stratified 5-fold CV x 50 repeats (seed 20260923 + repeat).
  AUROC = mean over repeats; 95% CI = stratified person bootstrap (2,000 resamples; weights applied to every repeat's
  out-of-fold predictions, averaged over repeats). **Permutation null:** 1,000 label permutations, each re-running the
  whole pipeline (selection included) on one CV repeat; p = (1 + #null >= observed) / 1,001.
* **Decision rule:** "metabolites separate LC from healthy" if the AUROC CI lower bound > 0.5 **and** permutation
  p <= 0.05; "strong" if, in addition, the lower bound >= 0.70; otherwise "null".

### Pre-specified secondary analyses (all reported)

* (a) Per metabolite: Hedges g (LC - healthy, log2) with a 2,000-bootstrap CI; Welch t p; Mann-Whitney p;
  sex-adjusted standardised difference (OLS log2 ~ LC + female, HC3), and sex + BMI adjusted as a sensitivity. BH-FDR
  across all metabolites of a tissue (primary family); BH within the authors' pathway domains also reported (their
  correction family; reproduction column).
* (b) Sex-only baseline and sex + metabolites (in-fold top-10 selection): Delta AUROC with the paired bootstrap.
* (c) All-metabolite L2 model (C from {0.01, 0.1, 1} by inner 3-fold CV) instead of top-10 selection.
* (d) Blood + muscle combined (people with both baseline samples), same pipeline as the primary.
* (e) VO2max (sheet `Figure1`) alone vs VO2max + blood metabolites (people with both), same CV; Delta AUROC.
* (f) PEM change (exploratory): per-metabolite Hedges g of the within-person change (Baseline -> 1 day after PEM),
  LC vs healthy, BH-FDR per tissue.
* (g) Direction vs the paper (table below): for each named metabolite, sign concordance of our Hedges g; for each named
  pathway, a per-person domain score (mean of z-scored log2 metabolites in the authors' domain) compared LC vs healthy
  (Hedges g, bootstrap CI), plus the share of domain metabolites whose g has the claimed sign.

### Published claims recorded before our computation (baseline, LC vs healthy; Appelman 2024 Results)

| tissue | claim (quoted from Results) | our check |
|---|---|---|
| muscle | "key metabolites of the tricarboxylic acid (TCA) cycle (including glutamate, FAD+, alpha-ketoglutarate and citric acid) were lower" | g < 0 for those 4; TCA domain score lower |
| muscle | "skeletal muscle glycolytic metabolites displayed few group differences" | Glycolysis domain score near 0 |
| muscle | "Skeletal muscle creatine concentrations were lower" | Creatine g < 0 |
| muscle | "lower S-adenosylmethionine (SAM)" | S-adenosyl methionine g < 0 |
| muscle | "lower hydroxyphenyl acetic acid" (the file has `Hydroxyphenyllactic acid`; mapped and flagged as a name mismatch) | g < 0 |
| muscle | "The ratio of citric acid to lactate in skeletal muscle was lower" | log2(citrate) - log2(lactate): g < 0 |
| muscle | "Many amino acids were not different between groups at rest" | Aminoacid domain score near 0 |
| blood | "Glycolytic metabolites in the venous blood were significantly higher" | Glycolysis domain score higher |
| blood | "pyruvate and other TCA cycle metabolites were lower" | Pyruvate g < 0; TCA domain score lower |
| blood | "Various blood metabolites within the purine pathway were lower" | Purines domain score lower |

The paper's statistics (generalised linear mixed model over three timepoints, Box-Cox, BH within pathway) are not
re-run exactly; the "reproduction" columns are our closest baseline-only analogue and are kept separate from our own
pre-specified analyses.

### Interpretation fixed in advance

46 people. With 25 vs 21, 80% power at two-sided alpha 0.05 needs |g| of about 0.85 for a single metabolite, larger
after FDR across 83-116 metabolites. Controls are healthy people who recovered from mild COVID-19, so separation is
against healthy people only (no disease, deconditioning or ME/CFS comparator). Patients able to undergo maximal CPET and
biopsies are the milder end. A cross-validated AUROC at n = 46 has a wide CI; the permutation null guards against
selection optimism.

## Outputs

`data/processed/phenotype_signatures__nhanes_labs`, `data/processed/phenotype_signatures__appelman_metabolomics`
(canonical phenotype_signatures columns), `results/tables/labs_dx_nhanes_*.csv`, `results/tables/labs_dx_appelman_*.csv`,
`results/LABS_VS_DIAGNOSIS_RESULTS.md`. Modules: `measure_it.labs.labs_vs_diagnosis` (Part 1),
`measure_it.labs.appelman_metabolomics` (Part 2), `measure_it.ingestion.appelman_lc_pem` (registration + audit),
`measure_it.labs.report` (results document).
