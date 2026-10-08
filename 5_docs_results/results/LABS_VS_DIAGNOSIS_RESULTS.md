# Do clinical labs and other non-wearable measurements separate invisible-illness labels?

**Plan (locked before any comparison):** `docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md` (sha256 `5572e783bb0ae949...`). **Modules:** `measure_it.labs.labs_vs_diagnosis` (Part 1), `measure_it.labs.appelman_metabolomics` (Part 2), `measure_it.ingestion.appelman_lc_pem`, `measure_it.labs.report`. **Tables:** `results/tables/labs_dx_nhanes_*.csv`, `results/tables/labs_dx_appelman_*.csv`. **Signature partitions:** `phenotype_signatures__nhanes_labs`, `phenotype_signatures__appelman_metabolomics`.

Reproduce: `uv run python -m measure_it.labs.plan && uv run python -m measure_it.labs.labs_vs_diagnosis && uv run python -m measure_it.ingestion.appelman_lc_pem && uv run python -m measure_it.labs.appelman_metabolomics && uv run python -m measure_it.labs.report`; tests: `uv run pytest tests/test_labs_vs_diagnosis.py`.

## Answer

**Part 1 (NHANES 2011-2014, 10,106 adults).** The positive-control gate passed: labs add clearly to demographics where they should diabetes_no_glycaemic_labs +0.147 (95% CI +0.136 to +0.159); weak_failing_kidneys +0.119 (95% CI +0.094 to +0.148); anemia_treatment +0.107 (95% CI +0.083 to +0.131); liver_condition +0.091 (95% CI +0.066 to +0.117); congestive_heart_failure +0.056 (95% CI +0.038 to +0.073); gout +0.048 (95% CI +0.031 to +0.065). So the pipeline can see lab signal, and the nulls below are informative within their sample sizes.
* **Labs-only models are not evidence of disease biology here.** They beat their permutation null for 24 of 25 modelled labels, but their AUROC is within 0.03 of (or below) the demographics-only AUROC for 16 of them: the lab panel re-learns age and sex (testosterone, creatinine, haemoglobin). The increment over demographics is the pre-specified reading.
* **Umbrella conditions identified by prescription reason codes: NULL.** Labs add nothing detectable to demographics for rx_fibromyalgia (n = 36) +0.022 (95% CI -0.027 to +0.072); rx_migraine (n = 81) -0.012 (95% CI -0.056 to +0.033); rx_insomnia (n = 186) +0.024 (95% CI -0.010 to +0.055); rx_myalgia (n = 104) -0.034 (95% CI -0.077 to +0.013). rx_ibs has 18 cases (no model). Primary-panel analytes with q < 0.05: rx_fibromyalgia 0, rx_migraine 2, rx_insomnia 12, rx_myalgia 3, rx_ibs 1.
* **Umbrella conditions against each other: NULL.** Rx fibromyalgia (case) vs Rx migraine (reference): 0 FDR analytes, Delta +0.057 (95% CI -0.050 to +0.156); Rx fibromyalgia (case) vs Rx insomnia (reference): 0 FDR analytes, no model; Rx migraine (case) vs Rx insomnia (reference): 0 FDR analytes, Delta +0.033 (95% CI -0.025 to +0.096); Rx fibromyalgia (case) vs rheumatoid arthritis (reference; 2013-2014 only): 2 FDR analytes, Delta +0.030 (95% CI -0.055 to +0.112). No routine lab tells fibromyalgia, migraine and insomnia (as treated conditions) apart.
* **ME/CFS-like proxy vs PHQ-9 >= 10 non-proxy: labs 'separate' them (+0.061 (95% CI +0.016 to +0.107)), by construction.** The separating analytes are lower glycaemia in the proxy (Glycohemoglobin -0.39; Glucose, serum -0.32; Mefox oxidation product -0.28): the proxy definition excludes diabetes and the other exclusion diagnoses, the depression comparison group does not. This is the case definition showing through, not a biological difference.
* **ME/CFS-like PROXY: a pre-specified 'labs add' that does not survive scrutiny.** Locked result +0.046 (95% CI +0.011 to +0.080) (n = 160); its strongest analytes are tobacco exposure and a LOWER HbA1c (Cotinine +0.47; Urinary Total NNAL +0.46; Platelet count +0.26; Mean cell volume +0.26; Glycohemoglobin -0.25; Blood urea nitrogen -0.24). Post hoc, with tobacco analytes removed and cotinine + BMI in the baseline the increment is +0.027 (-0.005 to +0.057); against controls without the proxy's exclusion diagnoses it is +0.005 (-0.020 to +0.032). The proxy's lab 'signal' is smoking plus the labs recognising the diseases the definition excludes, not an ME/CFS-like biology.
* **Symptom labels (fatigue, depressive symptoms, told sleep disorder): small increments, mostly BMI and smoking.** Locked: fatigue_symptom +0.027 (95% CI +0.017 to +0.037); depression_phq9 +0.028 (95% CI +0.016 to +0.039); told_sleep_disorder +0.048 (95% CI +0.032 to +0.063). Post hoc with cotinine + BMI in the baseline (tobacco analytes removed): fatigue_symptom +0.009 (+0.001 to +0.018); depression_phq9 +0.008 (-0.001 to +0.017); told_sleep_disorder -0.002 (-0.011 to +0.007).
* **Comparator diagnoses:** rheumatoid_arthritis +0.021 (95% CI +0.009 to +0.034); asthma +0.032 (95% CI +0.019 to +0.044); thyroid_problem +0.000 (95% CI -0.006 to +0.008); osteoarthritis +0.001 (95% CI -0.005 to +0.006); psoriasis +0.009 (95% CI -0.018 to +0.036); celiac_disease +0.020 (95% CI -0.047 to +0.087). Asthma (eosinophils) and rheumatoid arthritis carry small increments; thyroid problem, osteoarthritis, psoriasis and celiac disease (treated, mostly seronegative) do not.

**Part 2 (Appelman 2024 metabolomics, Long COVID vs healthy).** Plasma: CV AUROC 0.69 (0.56-0.82), permutation p 0.063 -> **null** under the locked rule; muscle: 0.67 (0.53-0.81), p 0.096 -> **null**. Secondary plasma + muscle: 0.76 (0.64-0.87), p 0.012 (meets the rule, but it is one of several secondaries at n = 25 vs 19). Per metabolite, q < 0.05 (BH per tissue): plasma 0/83, muscle 1/116 (Hydroxyphenyllactic acid). VO2max alone 0.76; adding plasma metabolites changes AUROC by -0.021 (-0.187 to +0.142). The paper's muscle TCA-lower direction is reproduced (individual CIs exclude 0, not FDR-significant at n = 44); its 'plasma glycolytic metabolites higher' is only partly reproduced (see the claims table). No ME/CFS group exists in this file.

**Bottom line.** Routine NHANES labs separate the conditions they are known to track (diabetes without glycaemic markers, kidney disease, anaemia, liver disease, heart failure, gout), and do not separate the invisible-illness labels available here from the general population or from each other beyond demographics, BMI and smoking. The one small open metabolomics set gives a borderline, underpowered Long COVID vs healthy signal that is not better than VO2max. These are nulls for these labels in these data (treated or self-reported cases, general-population or healthy controls), not evidence that no measurable biology exists.

## Part 1 — NHANES 2011-2014: laboratory analytes vs diagnosis-like labels

Population: 10,106 adults (>= 20, not pregnant, biochemistry profile measured; flow in `labs_dx_nhanes_population_flow.csv`). Analytes: 60 in the primary panel (models + screen), 380 extended-tier continuous analytes and 93 binary serology/STI/HPV results (screen only; subsample layers analysed among sampled people). Unweighted internal comparisons. Seed 20260923; 5-fold x 5 CV; 1000 participant bootstraps; 200 label permutations; 20 lab-block shuffles.

### Positive-control gate (can the pipeline see lab signal at all?)

Pre-specified: labs must add to demographics for `diabetes_no_glycaemic_labs`, `weak_failing_kidneys`, `gout`. Result: **PASSED** (diabetes_no_glycaemic_labs: True, weak_failing_kidneys: True, gout: True).

| label | cases / controls | AUROC demog. | AUROC labs only | AUROC demog.+labs | Delta (95% CI) | strongest FDR analytes (SD units, age/sex-adjusted) |
|---|---|---|---|---|---|---|
| diabetes_no_glycaemic_labs | 1278 / 8555 | 0.758 | 0.889 (0.880-0.899) | 0.904 (0.897-0.913) | +0.147 (+0.136 to +0.159) | Albumin creatinine ratio +0.72; Albumin, urine +0.63; Osmolality +0.54; Direct HDL-Cholesterol -0.49 |
| diabetes_with_glycaemic_labs | 1278 / 8555 | 0.758 | 0.940 (0.932-0.947) | 0.943 (0.936-0.949) | +0.185 (+0.172 to +0.197) | Glycohemoglobin +1.75; Glucose, serum +1.50; Albumin creatinine ratio +0.72; Albumin, urine +0.63 |
| weak_failing_kidneys | 340 / 9757 | 0.694 | 0.800 (0.773-0.827) | 0.813 (0.789-0.838) | +0.119 (+0.094 to +0.148) | Creatinine +1.45; Albumin creatinine ratio +1.06; Albumin, urine +0.92; Blood urea nitrogen +0.87 |
| gout | 418 / 9684 | 0.762 | 0.787 (0.766-0.810) | 0.810 (0.791-0.829) | +0.048 (+0.031 to +0.065) | Uric acid +0.58; Creatinine +0.42; Albumin creatinine ratio +0.39; Gamma glutamyl transferase +0.37 |
| congestive_heart_failure | 329 / 9759 | 0.810 | 0.847 (0.824-0.870) | 0.866 (0.846-0.884) | +0.056 (+0.038 to +0.073) | Creatinine +0.76; Albumin creatinine ratio +0.61; Uric acid +0.55; Albumin, urine +0.53 |
| anemia_treatment | 403 / 9694 | 0.699 | 0.804 (0.784-0.824) | 0.806 (0.786-0.825) | +0.107 (+0.083 to +0.131) | Red cell distribution width +0.92; Hemoglobin -0.72; Hematocrit -0.67; Mean cell hemoglobin -0.47 |
| liver_condition | 404 / 9689 | 0.656 | 0.725 (0.700-0.750) | 0.747 (0.722-0.772) | +0.091 (+0.066 to +0.117) | Aspartate aminotransferase AST +0.82; Gamma glutamyl transferase +0.71; Alanine aminotransferase ALT +0.68; Globulin +0.30 |

### Headline: every label, nulls included

Delta = AUROC(demographics + labs) - AUROC(demographics), paired participant bootstrap. 'labs only' p = label-permutation p. Reading applies the locked rule (CI > 0 and lab-block shuffle p <= 0.05; a Delta < 0.02 is called negligible). # FDR = primary-panel analytes with q < 0.05 / tested.

| label | group | cases / controls | AUROC demog. | AUROC labs only (perm p) | Delta demog.+labs (95% CI) | Delta over demog.+BMI | # FDR | reading |
|---|---|---|---|---|---|---|---|---|
| rx_fibromyalgia | umbrella | 36 / 5229 | 0.729 | 0.734 (0.005) | +0.022 (-0.027 to +0.072) | +0.003 (-0.040 to +0.048) | 0/60 | NULL: labs do not add to demographics |
| rx_migraine | umbrella | 81 / 5184 | 0.659 | 0.635 (0.005) | -0.012 (-0.056 to +0.033) | -0.012 (-0.053 to +0.030) | 2/60 | NULL: labs do not add to demographics |
| rx_ibs | umbrella | 18 / 5247 | - | - | - | - | 1/60 | no model (fewer than 25 cases or controls); per-analyte contrasts only |
| rx_insomnia | umbrella | 186 / 5079 | 0.677 | 0.688 (0.005) | +0.024 (-0.010 to +0.055) | +0.028 (-0.004 to +0.058) | 12/60 | NULL: labs do not add to demographics |
| rx_myalgia | umbrella | 104 / 5161 | 0.634 | 0.585 (0.030) | -0.034 (-0.077 to +0.013) | -0.046 (-0.089 to -0.002) | 3/60 | NULL: labs do not add to demographics |
| mecfs_like_proxy (PROXY) | umbrella | 160 / 8969 | 0.648 | 0.644 (0.005) | +0.046 (+0.011 to +0.080) | +0.050 (+0.014 to +0.084) | 9/60 | labs add to demographics |
| fatigue_symptom (PROXY) | umbrella | 1499 / 7630 | 0.651 | 0.658 (0.005) | +0.027 (+0.017 to +0.037) | +0.020 (+0.011 to +0.029) | 41/60 | labs add to demographics |
| told_sleep_disorder | umbrella | 951 / 9141 | 0.624 | 0.657 (0.005) | +0.048 (+0.032 to +0.063) | +0.001 (-0.008 to +0.011) | 37/60 | labs add to demographics |
| depression_phq9 (PROXY) | umbrella | 850 / 8246 | 0.712 | 0.697 (0.005) | +0.028 (+0.016 to +0.039) | +0.022 (+0.011 to +0.033) | 33/60 | labs add to demographics |
| celiac_disease | comparator | 42 / 10057 | 0.626 | 0.631 (0.010) | +0.020 (-0.047 to +0.087) | +0.016 (-0.045 to +0.079) | 4/60 | NULL: labs do not add to demographics |
| thyroid_problem | comparator | 999 / 9084 | 0.762 | 0.742 (0.005) | +0.000 (-0.006 to +0.008) | -0.001 (-0.008 to +0.006) | 24/60 | NULL: labs do not add to demographics |
| rheumatoid_arthritis | comparator | 460 / 7500 | 0.794 | 0.770 (0.005) | +0.021 (+0.009 to +0.034) | +0.017 (+0.006 to +0.028) | 23/60 | labs add to demographics |
| osteoarthritis | comparator | 1143 / 7500 | 0.827 | 0.764 (0.005) | +0.001 (-0.005 to +0.006) | -0.002 (-0.007 to +0.002) | 28/60 | NULL: labs do not add to demographics |
| psoriasis | comparator | 284 / 9811 | 0.580 | 0.581 (0.005) | +0.009 (-0.018 to +0.036) | -0.002 (-0.028 to +0.022) | 5/60 | NULL: labs do not add to demographics |
| asthma | comparator | 1519 / 8579 | 0.602 | 0.624 (0.005) | +0.032 (+0.019 to +0.044) | +0.025 (+0.013 to +0.038) | 24/60 | labs add to demographics |
| diabetes_no_glycaemic_labs | positive_control | 1278 / 8555 | 0.758 | 0.889 (0.005) | +0.147 (+0.136 to +0.159) | +0.112 (+0.101 to +0.122) | 44/58 | labs add to demographics |
| diabetes_with_glycaemic_labs | positive_control | 1278 / 8555 | 0.758 | 0.940 (0.005) | +0.185 (+0.172 to +0.197) | +0.148 (+0.136 to +0.159) | 46/60 | labs add to demographics |
| weak_failing_kidneys | positive_control | 340 / 9757 | 0.694 | 0.800 (0.005) | +0.119 (+0.094 to +0.148) | +0.115 (+0.089 to +0.141) | 37/60 | labs add to demographics |
| gout | positive_control | 418 / 9684 | 0.762 | 0.787 (0.005) | +0.048 (+0.031 to +0.065) | +0.031 (+0.016 to +0.046) | 26/60 | labs add to demographics |
| congestive_heart_failure | positive_control | 329 / 9759 | 0.810 | 0.847 (0.005) | +0.056 (+0.038 to +0.073) | +0.043 (+0.026 to +0.060) | 43/60 | labs add to demographics |
| anemia_treatment | positive_control | 403 / 9694 | 0.699 | 0.804 (0.005) | +0.107 (+0.083 to +0.131) | +0.108 (+0.083 to +0.132) | 27/60 | labs add to demographics |
| liver_condition | positive_control | 404 / 9689 | 0.656 | 0.725 (0.005) | +0.091 (+0.066 to +0.117) | +0.087 (+0.062 to +0.114) | 23/60 | labs add to demographics |
| pair_rx_fibromyalgia_vs_rx_migraine | pairwise | 32 / 77 | 0.487 | 0.548 (0.204) | +0.057 (-0.050 to +0.156) | - | 0/60 | NULL: labs do not add to demographics |
| pair_rx_fibromyalgia_vs_rx_insomnia | pairwise | 22 / 172 | - | - | - | - | 0/60 | no model (fewer than 25 cases or controls); per-analyte contrasts only |
| pair_rx_migraine_vs_rx_insomnia | pairwise | 69 / 174 | 0.610 | 0.638 (0.010) | +0.033 (-0.025 to +0.096) | - | 0/60 | NULL: labs do not add to demographics |
| pair_rx_fibromyalgia_vs_rheumatoid_arthritis | pairwise | 29 / 237 | 0.650 | 0.643 (0.005) | +0.030 (-0.055 to +0.112) | - | 2/60 | NULL: labs do not add to demographics |
| pair_mecfs_like_proxy_vs_depression_non_proxy (PROXY) | pairwise | 160 / 766 | 0.545 | 0.590 (0.005) | +0.061 (+0.016 to +0.107) | - | 15/60 | labs add to demographics |

### Umbrella labels: what the labs do and do not show

* **No increment over demographics (null or not modelled):** rx_fibromyalgia, rx_migraine, rx_insomnia, rx_myalgia.
* **Increment present but negligible (< 0.02 AUROC):** none.
* **Increment >= 0.02:** mecfs_like_proxy, fatigue_symptom, told_sleep_disorder, depression_phq9.
* BH across the umbrella labels' Delta bootstrap p-values (descriptive): rx_fibromyalgia 0.494, rx_migraine 0.632, rx_insomnia 0.230, rx_myalgia 0.232, mecfs_like_proxy 0.012, fatigue_symptom 0.003, told_sleep_disorder 0.003, depression_phq9 0.003.

### Which analytes differ, per label (primary panel, q < 0.05; age/sex-adjusted SD units, BMI-adjusted in brackets)

| label | cases | FDR primary / extended | strongest primary-panel differences |
|---|---|---|---|
| rx_fibromyalgia | 36 | 0 / 3 | none |
| rx_migraine | 81 | 2 / 2 | Glycohemoglobin -0.26 [-0.26]; Glucose, serum -0.24 [-0.26] |
| rx_ibs | 18 | 1 / 1 | 25OHD2+25OHD3 +0.79 [+0.75] |
| rx_insomnia | 186 | 12 / 8 | Urinary Total NNAL +0.37 [+0.37]; Folic acid +0.33 [+0.31]; Gamma glutamyl transferase +0.32 [+0.33]; Creatinine +0.30 [+0.31]; Cotinine +0.28 [+0.28]; Chloride -0.26 [-0.25] |
| rx_myalgia | 104 | 3 / 0 | Urinary Total NNAL +0.38 [+0.40]; Cotinine +0.34 [+0.36]; Triglycerides +0.30 [+0.24] |
| mecfs_like_proxy | 160 | 9 / 18 | Cotinine +0.47 [+0.47]; Urinary Total NNAL +0.46 [+0.47]; Platelet count +0.26 [+0.26]; Mean cell volume +0.26 [+0.27]; Glycohemoglobin -0.25 [-0.27]; Blood urea nitrogen -0.24 [-0.23] |
| fatigue_symptom | 1499 | 41 / 101 | Urinary Total NNAL +0.35 [+0.36]; Cotinine +0.34 [+0.36]; Gamma glutamyl transferase +0.23 [+0.17]; Albumin -0.23 [-0.14]; Direct HDL-Cholesterol -0.21 [-0.13]; Alkaline phosphotase +0.19 [+0.15] |
| told_sleep_disorder | 951 | 37 / 64 | Direct HDL-Cholesterol -0.29 [-0.09]; Albumin -0.23 [-0.03]; Total bilirubin -0.21 [-0.12]; White blood cell count +0.21 [+0.08]; Red cell distribution width +0.20 [+0.09]; Gamma glutamyl transferase +0.20 [+0.09] |
| depression_phq9 | 850 | 33 / 132 | Urinary Total NNAL +0.50 [+0.52]; Cotinine +0.47 [+0.49]; Triglycerides +0.27 [+0.20]; Alkaline phosphotase +0.27 [+0.22]; Gamma glutamyl transferase +0.27 [+0.21]; Direct HDL-Cholesterol -0.26 [-0.16] |
| celiac_disease | 42 | 4 / 30 | Total protein -0.66 [-0.67]; Globulin -0.56 [-0.53]; Mean cell volume +0.43 [+0.40]; Mean cell hemoglobin +0.41 [+0.39] |
| thyroid_problem | 999 | 24 / 41 | 25OHD2+25OHD3 +0.26 [+0.29]; RBC folate +0.20 [+0.20]; Alanine aminotransferase ALT +0.17 [+0.15]; Aspartate aminotransferase AST +0.15 [+0.15]; Monocyte number +0.13 [+0.10]; Red blood cell count +0.11 [+0.10] |
| rheumatoid_arthritis | 460 | 23 / 49 | Urinary Total NNAL +0.38 [+0.38]; Cotinine +0.35 [+0.36]; Gamma glutamyl transferase +0.29 [+0.21]; Alkaline phosphotase +0.26 [+0.18]; Folic acid +0.26 [+0.27]; Lactate dehydrogenase +0.26 [+0.21] |
| osteoarthritis | 1143 | 28 / 99 | Urinary Total NNAL +0.20 [+0.21]; RBC folate +0.20 [+0.16]; White blood cell count +0.20 [+0.09]; Total protein -0.19 [-0.17]; Segmented neutrophils num +0.19 [+0.10]; Cotinine +0.19 [+0.20] |
| psoriasis | 284 | 5 / 5 | Lymphocyte percent -0.33 [-0.33]; Segmented neutrophils percent +0.28 [+0.28]; Segmented neutrophils num +0.22 [+0.20]; Lymphocyte number -0.18 [-0.21]; Creatine Phosphokinase(CPK) -0.16 [-0.18] |
| asthma | 1519 | 24 / 90 | Eosinophils number +0.25 [+0.22]; Urinary Total NNAL +0.22 [+0.22]; Eosinophils percent +0.19 [+0.19]; Cotinine +0.18 [+0.18]; White blood cell count +0.16 [+0.10]; Segmented neutrophils num +0.15 [+0.10] |
| diabetes_no_glycaemic_labs | 1278 | 44 / 60 | Albumin creatinine ratio +0.72 [+0.71]; Albumin, urine +0.63 [+0.61]; Osmolality +0.54 [+0.51]; Direct HDL-Cholesterol -0.49 [-0.32]; Total Cholesterol -0.46 [-0.47]; Mefox oxidation product +0.45 [+0.40] |
| diabetes_with_glycaemic_labs | 1278 | 46 / 68 | Glycohemoglobin +1.75 [+1.69]; Glucose, serum +1.50 [+1.44]; Albumin creatinine ratio +0.72 [+0.71]; Albumin, urine +0.63 [+0.61]; Osmolality +0.54 [+0.51]; Direct HDL-Cholesterol -0.49 [-0.32] |
| weak_failing_kidneys | 340 | 37 / 31 | Creatinine +1.45 [+1.43]; Albumin creatinine ratio +1.06 [+1.05]; Albumin, urine +0.92 [+0.91]; Blood urea nitrogen +0.87 [+0.84]; Mefox oxidation product +0.69 [+0.66]; Methylmalonic Acid +0.64 [+0.65] |
| gout | 418 | 26 / 27 | Uric acid +0.58 [+0.48]; Creatinine +0.42 [+0.41]; Albumin creatinine ratio +0.39 [+0.37]; Gamma glutamyl transferase +0.37 [+0.31]; Albumin, urine +0.34 [+0.33]; Red cell distribution width +0.27 [+0.20] |
| congestive_heart_failure | 329 | 43 / 33 | Creatinine +0.76 [+0.71]; Albumin creatinine ratio +0.61 [+0.57]; Uric acid +0.55 [+0.44]; Albumin, urine +0.53 [+0.50]; Mefox oxidation product +0.50 [+0.44]; Total Cholesterol -0.49 [-0.50] |
| anemia_treatment | 403 | 27 / 39 | Red cell distribution width +0.92 [+0.88]; Hemoglobin -0.72 [-0.72]; Hematocrit -0.67 [-0.67]; Mean cell hemoglobin -0.47 [-0.43]; Creatinine +0.45 [+0.42]; MCHC -0.41 [-0.38] |
| liver_condition | 404 | 23 / 59 | Aspartate aminotransferase AST +0.82 [+0.81]; Gamma glutamyl transferase +0.71 [+0.69]; Alanine aminotransferase ALT +0.68 [+0.66]; Globulin +0.30 [+0.29]; Alkaline phosphotase +0.25 [+0.24]; Triglycerides +0.22 [+0.19] |
| pair_rx_fibromyalgia_vs_rx_migraine | 32 | 0 / 0 | none |
| pair_rx_fibromyalgia_vs_rx_insomnia | 22 | 0 / 0 | none |
| pair_rx_migraine_vs_rx_insomnia | 69 | 0 / 0 | none |
| pair_rx_fibromyalgia_vs_rheumatoid_arthritis | 29 | 2 / 0 | Albumin creatinine ratio -0.44 [-0.45]; Albumin, urine -0.39 [-0.39] |
| pair_mecfs_like_proxy_vs_depression_non_proxy | 160 | 15 / 4 | Glycohemoglobin -0.39 [-0.33]; Glucose, serum -0.32 [-0.27]; Mefox oxidation product -0.28 [-0.23]; Osmolality -0.27 [-0.24]; Mean platelet volume -0.24 [-0.22]; Mean cell volume +0.21 [+0.18] |

Full per-analyte table (all tiers, CIs, p, q, BMI-adjusted): `labs_dx_nhanes_effects.csv`; also `phenotype_signatures__nhanes_labs`.

### Umbrella conditions against each other (pairwise contrasts)

| contrast | n | # FDR primary | FDR analytes | AUROC demog. / demog.+labs | Delta (95% CI) |
|---|---|---|---|---|---|
| Rx fibromyalgia (case) vs Rx migraine (reference) | 32 / 77 | 0/60 | none | 0.487 / 0.544 | +0.057 (-0.050 to +0.156) |
| Rx fibromyalgia (case) vs Rx insomnia (reference) | 22 / 172 | 0/60 | none | - | no model (< 25 per group) |
| Rx migraine (case) vs Rx insomnia (reference) | 69 / 174 | 0/60 | none | 0.610 / 0.643 | +0.033 (-0.025 to +0.096) |
| Rx fibromyalgia (case) vs rheumatoid arthritis (reference; 2013-2014 only) | 29 / 237 | 2/60 | Albumin creatinine ratio -0.44; Albumin, urine -0.39 | 0.650 / 0.680 | +0.030 (-0.055 to +0.112) |
| ME/CFS-like PROXY (case) vs PHQ-9 >= 10 non-proxy (reference) | 160 / 766 | 15/60 | Glycohemoglobin -0.39; Glucose, serum -0.32; Mefox oxidation product -0.28; Osmolality -0.27 | 0.545 / 0.606 | +0.061 (+0.016 to +0.107) |

### Binary serology / STI / HPV (extended tier, q < 0.05)

Adjusted log odds ratios. Rows with fewer than 5 positive cases are fragile (a single person moves them).

| label_id | lab_name | layer | positive cases / cases | log OR (95% CI) | q |
|---|---|---|---|---|---|
| rx_fibromyalgia | HPV type 33 | infectious_hpv | 2 / 20 | +3.02 (+1.41 to +4.64) | 0.006 |
| rx_fibromyalgia | HPV type 35 | infectious_hpv | 2 / 20 | +2.45 (+0.91 to +3.99) | 0.018 |
| rx_fibromyalgia | HPV Type 51 | infectious_hpv | 1 / 29 | +3.77 (+1.36 to +6.17) | 0.018 |
| rx_migraine | Hepatitis A antibody | infectious_hepatitis | 22 / 81 | -0.84 (-1.34 to -0.35) | 0.035 |
| rx_ibs | HPV Type 66 | infectious_hpv | 1 / 12 | +3.03 (+0.86 to +5.20) | 0.037 |
| rx_insomnia | Hepatitis A antibody | infectious_hepatitis | 68 / 186 | -0.52 (-0.83 to -0.22) | 0.048 |
| rx_insomnia | HPV type 33 | infectious_hpv | 4 / 101 | +1.78 (+0.67 to +2.89) | 0.048 |
| fatigue_symptom | Hepatitis C antibody | infectious_hepatitis | 26 / 671 | +0.84 (+0.37 to +1.32) | 0.021 |
| fatigue_symptom | Herpes Simplex Virus II | infectious_herpes_cmv | 180 / 674 | +0.41 (+0.21 to +0.60) | 0.005 |
| told_sleep_disorder | Hepatitis A antibody | infectious_hepatitis | 382 / 951 | -0.42 (-0.56 to -0.28) | 2.8e-07 |
| told_sleep_disorder | Hepatitis B core antibody | infectious_hepatitis | 67 / 951 | -0.45 (-0.71 to -0.19) | 0.015 |
| told_sleep_disorder | Hepatitis C antibody | infectious_hepatitis | 21 / 411 | +0.94 (+0.43 to +1.44) | 0.008 |
| told_sleep_disorder | Herpes Simplex Virus II | infectious_herpes_cmv | 107 / 370 | +0.47 (+0.22 to +0.71) | 0.008 |
| told_sleep_disorder | TB coded result | infectious_tuberculosis_igra | 29 / 409 | -0.66 (-1.05 to -0.27) | 0.017 |
| depression_phq9 | Hepatitis C antibody | infectious_hepatitis | 17 / 387 | +0.88 (+0.32 to +1.43) | 0.031 |
| depression_phq9 | Roche HPV linear array summary result | infectious_hpv | 222 / 436 | +0.36 (+0.16 to +0.56) | 0.016 |
| depression_phq9 | Herpes Simplex Virus II | infectious_herpes_cmv | 117 / 380 | +0.58 (+0.34 to +0.82) | 2.4e-04 |
| depression_phq9 | Chlamydia Pgp3 MBA | infectious_sti_urine | 35 / 67 | +0.81 (+0.30 to +1.32) | 0.031 |
| depression_phq9 | Toxocara antibody result | infectious_parasite_serology_surplus | 72 / 785 | +0.42 (+0.16 to +0.68) | 0.031 |
| celiac_disease | HPV type 42 | infectious_hpv | 3 / 20 | +2.16 (+0.90 to +3.42) | 0.019 |
| celiac_disease | Tissue transglutaminase | immune_celiac_serology | 2 / 41 | +2.41 (+0.95 to +3.87) | 0.019 |
| rheumatoid_arthritis | Hepatitis C antibody | infectious_hepatitis | 11 / 215 | +1.14 (+0.45 to +1.83) | 0.029 |
| rheumatoid_arthritis | Hepatitis E IgG antibody | infectious_hepatitis | 23 / 459 | -0.74 (-1.18 to -0.31) | 0.029 |
| rheumatoid_arthritis | Toxocara antibody result | infectious_parasite_serology_surplus | 56 / 427 | +0.57 (+0.27 to +0.88) | 0.019 |
| osteoarthritis | Hepatitis A antibody | infectious_hepatitis | 495 / 1143 | -0.58 (-0.72 to -0.44) | 2.3e-14 |
| osteoarthritis | Hepatitis B core antibody | infectious_hepatitis | 91 / 1143 | -0.52 (-0.76 to -0.27) | 7.1e-04 |
| osteoarthritis | Hepatitis B Surface Antibody | infectious_hepatitis | 158 / 1142 | -0.46 (-0.65 to -0.27) | 6.1e-05 |
| osteoarthritis | Hepatitis C antibody | infectious_hepatitis | 18 / 508 | +0.93 (+0.34 to +1.52) | 0.035 |
| osteoarthritis | TB coded result | infectious_tuberculosis_igra | 45 / 504 | -0.52 (-0.87 to -0.18) | 0.046 |
| psoriasis | Hepatitis A antibody | infectious_hepatitis | 110 / 284 | -0.43 (-0.68 to -0.19) | 0.018 |
| psoriasis | Tissue transglutaminase | immune_celiac_serology | 6 / 283 | +1.58 (+0.72 to +2.44) | 0.018 |
| asthma | Hepatitis A antibody | infectious_hepatitis | 602 / 1519 | -0.35 (-0.46 to -0.24) | 1.2e-07 |
| anemia_treatment | HPV type 68 | infectious_hpv | 10 / 217 | +1.25 (+0.55 to +1.95) | 0.030 |
| anemia_treatment | Herpes Simplex Virus II | infectious_herpes_cmv | 67 / 188 | +0.54 (+0.22 to +0.87) | 0.035 |
| liver_condition | Hepatitis A antibody | infectious_hepatitis | 243 / 404 | +0.43 (+0.23 to +0.64) | 9.3e-04 |
| liver_condition | Hepatitis B core antibody | infectious_hepatitis | 96 / 404 | +1.03 (+0.79 to +1.28) | 4.2e-15 |
| liver_condition | Hepatitis B surface antigen | infectious_hepatitis | 10 / 404 | +1.29 (+0.61 to +1.97) | 0.003 |
| liver_condition | Hepatitis D antibody | infectious_hepatitis | 5 / 404 | +1.91 (+0.88 to +2.94) | 0.003 |
| liver_condition | Hepatitis B Surface Antibody | infectious_hepatitis | 111 / 404 | +0.44 (+0.21 to +0.68) | 0.003 |
| liver_condition | Hepatitis C antibody | infectious_hepatitis | 45 / 193 | +3.16 (+2.71 to +3.61) | 1.6e-40 |
| liver_condition | HPV Type 83 | infectious_hpv | 4 / 301 | +1.93 (+0.77 to +3.09) | 0.012 |

### POST HOC: is the symptom-label signal smoking? (not pre-specified)

Serum cotinine and urinary NNAL (tobacco exposure) are among the largest differences for several symptom labels, which was seen only after the locked analysis ran. Primary-panel contrasts re-estimated with log serum cotinine as an extra covariate (tobacco analytes omitted); number with q < 0.05:

| label | q < 0.05 age/sex (non-tobacco analytes) | q < 0.05 + cotinine | tested |
|---|---|---|---|
| rx_fibromyalgia | 0 | 0 | 58 |
| rx_migraine | 2 | 2 | 58 |
| rx_ibs | 1 | 1 | 58 |
| rx_insomnia | 10 | 5 | 58 |
| rx_myalgia | 1 | 0 | 58 |
| mecfs_like_proxy | 7 | 5 | 58 |
| fatigue_symptom | 39 | 34 | 58 |
| told_sleep_disorder | 35 | 33 | 58 |
| depression_phq9 | 31 | 30 | 58 |
| celiac_disease | 4 | 4 | 58 |
| thyroid_problem | 24 | 26 | 58 |
| rheumatoid_arthritis | 21 | 19 | 58 |
| osteoarthritis | 26 | 27 | 58 |
| psoriasis | 5 | 6 | 58 |
| asthma | 22 | 19 | 58 |
| diabetes_no_glycaemic_labs | 44 | 44 | 56 |
| diabetes_with_glycaemic_labs | 46 | 46 | 58 |
| weak_failing_kidneys | 37 | 38 | 58 |
| gout | 26 | 26 | 58 |
| congestive_heart_failure | 41 | 41 | 58 |
| anemia_treatment | 27 | 28 | 58 |
| liver_condition | 21 | 22 | 58 |
| pair_rx_fibromyalgia_vs_rx_migraine | 0 | 0 | 58 |
| pair_rx_fibromyalgia_vs_rx_insomnia | 0 | 0 | 58 |
| pair_rx_migraine_vs_rx_insomnia | 0 | 0 | 58 |
| pair_rx_fibromyalgia_vs_rheumatoid_arthritis | 2 | 2 | 58 |
| pair_mecfs_like_proxy_vs_depression_non_proxy | 15 | 11 | 58 |

### POST HOC: do the symptom-label increments survive smoking, BMI and the proxy's own exclusions?

Not pre-specified; run after the locked results were seen. Labs here exclude the tobacco analytes. (i) `smoking_adjusted`: baseline = demographics + log serum cotinine (+ BMI). (ii) `controls_without_exclusion_dx`: ME/CFS-like proxy cases vs controls with none of the proxy's exclusion diagnoses, so labs cannot score by recognising the excluded diseases (diabetes via HbA1c, anaemia, thyroid, liver).

| label_id | analysis | baseline | cases / n | AUROC base -> + labs | Delta (95% CI) |
|---|---|---|---|---|---|
| mecfs_like_proxy | smoking_adjusted | demographics + log cotinine | 160 / 9129 | 0.671 -> 0.694 | +0.022 (-0.008 to +0.053) |
| mecfs_like_proxy | smoking_adjusted | demographics + bmi + log cotinine | 160 / 9129 | 0.666 -> 0.693 | +0.027 (-0.005 to +0.057) |
| fatigue_symptom | smoking_adjusted | demographics + log cotinine | 1499 / 9129 | 0.660 -> 0.677 | +0.017 (+0.008 to +0.027) |
| fatigue_symptom | smoking_adjusted | demographics + bmi + log cotinine | 1499 / 9129 | 0.670 -> 0.680 | +0.009 (+0.001 to +0.018) |
| told_sleep_disorder | smoking_adjusted | demographics + log cotinine | 951 / 10092 | 0.625 -> 0.672 | +0.047 (+0.032 to +0.063) |
| told_sleep_disorder | smoking_adjusted | demographics + bmi + log cotinine | 951 / 10092 | 0.709 -> 0.706 | -0.002 (-0.011 to +0.007) |
| depression_phq9 | smoking_adjusted | demographics + log cotinine | 850 / 9096 | 0.724 -> 0.739 | +0.015 (+0.006 to +0.025) |
| depression_phq9 | smoking_adjusted | demographics + bmi + log cotinine | 850 / 9096 | 0.733 -> 0.741 | +0.008 (-0.001 to +0.017) |
| pair_mecfs_like_proxy_vs_depression_non_proxy | smoking_adjusted | demographics + log cotinine | 160 / 926 | 0.539 -> 0.608 | +0.069 (+0.025 to +0.115) |
| pair_mecfs_like_proxy_vs_depression_non_proxy | smoking_adjusted | demographics + bmi + log cotinine | 160 / 926 | 0.566 -> 0.609 | +0.043 (+0.004 to +0.083) |
| asthma | smoking_adjusted | demographics + log cotinine | 1519 / 10098 | 0.605 -> 0.631 | +0.026 (+0.014 to +0.038) |
| asthma | smoking_adjusted | demographics + bmi + log cotinine | 1519 / 10098 | 0.619 -> 0.638 | +0.019 (+0.008 to +0.029) |
| rheumatoid_arthritis | smoking_adjusted | demographics + log cotinine | 460 / 7960 | 0.797 -> 0.815 | +0.018 (+0.006 to +0.029) |
| rheumatoid_arthritis | smoking_adjusted | demographics + bmi + log cotinine | 460 / 7960 | 0.807 -> 0.819 | +0.012 (+0.002 to +0.022) |
| mecfs_like_proxy | controls_without_exclusion_dx | demographics | 160 / 5187 | 0.744 -> 0.749 | +0.005 (-0.020 to +0.032) |
| mecfs_like_proxy | controls_without_exclusion_dx | demographics + bmi | 160 / 5187 | 0.740 -> 0.750 | +0.009 (-0.018 to +0.037) |

### Deviations from the locked plan

* Primary-panel analytes failing the 70% coverage rule (estradiol, SHBG, serum hydroxycotinine) were screened in the extended tier rather than dropped.
* Same-analyte unit duplicates (urine albumin URXUMS = URXUMA; serum folate LBDFOT = LBDFOTSI) were dropped (same label, rank correlation > 0.999), like LBXSCH.
* NHANES reuses variable names across files (LBXHCT = hematocrit in CBC and hydroxycotinine in COT_H; LBXIN in GLU_G and INS_H): extended-tier collisions are namespaced `<var>__<layer>`; the INS_H insulin is also removed for the no-glycaemic diabetes label.
* Solver: newton-cholesky (same L2 objective as lbfgs, chosen for speed).
* The post-hoc cotinine-adjusted contrasts and post-hoc models above.
* The lab-block shuffle reached its floor (p = 1/21) for nearly every label, including labels whose Delta CI includes 0: shuffled lab blocks lower the penalised model's AUROC, so beating them is a weak test (as in results/WEARABLE_NHANES_RESULTS.md). The Delta CI is the operative part of the rule.

## Part 2 — Appelman 2024 (MUSCLE-ME trial) plasma and muscle metabolomics: Long COVID vs healthy

Source `appelman_lc_pem_source` (41467_2023_44432_MOESM4_ESM.xlsx sha256 726be4203131dad4848f7b7aef9bd3d038ad4c21d5ea55a022480a7de2d492db (Last-Modified Wed, 19 Feb 2025 20:38:37 GMT)). Baseline: plasma 25 Long COVID vs 21 healthy (the lead's '23 vs 21' is the VO2max-and-blood overlap), muscle 25 vs 19. **There is no ME/CFS group in this file**, so the ME/CFS contrast was not possible. Age is not released; sex (and BMI) only. 5-fold x 50 CV with in-fold top-10 selection; 2000 stratified bootstraps; 1000 permutations re-running the whole pipeline.

### Primary (locked): does each tissue's metabolome separate Long COVID from healthy?

| tissue | n LC / healthy | CV AUROC (95% bootstrap CI) | permutation null mean / 95th pct | perm p | locked decision |
|---|---|---|---|---|---|
| blood | 25 / 21 | 0.694 (0.555-0.819) | 0.496 / 0.674 | 0.063 | null |
| muscle | 25 / 19 | 0.674 (0.534-0.809) | 0.505 / 0.686 | 0.096 | null |

Reading: the bootstrap CI resamples people but not the model-fitting, so it can exclude 0.5 while the permutation test (which re-runs selection and fitting) does not reject; the locked rule needs both.

### Secondary analyses (pre-specified)

| model | tissue | n LC / healthy | CV AUROC (95% CI) | perm p | decision rule |
|---|---|---|---|---|---|
| blood_sex_only | blood | 25 / 21 | 0.335 (0.310-0.359) | - | - |
| blood_sex_plus_top10 | blood | 25 / 21 | 0.682 (0.542-0.810) | - | - |
| blood_all_l2cv | blood | 25 / 21 | 0.647 (0.500-0.784) | - | - |
| muscle_sex_only | muscle | 25 / 19 | 0.374 (0.287-0.464) | - | - |
| muscle_sex_plus_top10 | muscle | 25 / 19 | 0.662 (0.525-0.798) | - | - |
| muscle_all_l2cv | muscle | 25 / 19 | 0.696 (0.556-0.834) | - | - |
| blood_plus_muscle_top10 | blood+muscle | 25 / 19 | 0.762 (0.636-0.874) | 0.012 | metabolites separate LC from healthy |
| vo2max_only | blood (people with VO2max) | 23 / 21 | 0.762 (0.612-0.897) | - | - |
| vo2max_plus_blood_top10 | blood (people with VO2max) | 23 / 21 | 0.741 (0.610-0.854) | - | - |

| tissue | comparison | Delta AUROC (95% CI) | boot p |
|---|---|---|---|
| blood | blood_sex_plus_top10 minus blood_sex_only | +0.347 (+0.211 to +0.476) | 5.0e-04 |
| muscle | muscle_sex_plus_top10 minus muscle_sex_only | +0.288 (+0.135 to +0.451) | 5.0e-04 |
| blood (people with VO2max) | vo2max_plus_blood_top10 minus vo2max_only | -0.021 (-0.187 to +0.142) | 0.784 |

The sex-only models score below 0.5 (a known cross-validation artefact for a feature that carries no signal in a sex-matched sample: each training fold's sex imbalance is the mirror image of the test fold's), which inflates the sex + metabolites Delta; read the metabolite-only AUROCs instead. Adding plasma metabolites to VO2max does not improve on VO2max alone.

### (a) Per-metabolite differences

| tissue | analysis | tested | raw p < 0.05 | q < 0.05 (BH, tissue) | q < 0.05 sex-adjusted | q < 0.05 within domain (reproduction) | smallest p (Hedges g) |
|---|---|---|---|---|---|---|---|
| blood | baseline | 83 | 8 | 0 | 5 | 6 | Pyruvate -0.98; Dihydroxyacetone-P +0.95; Gluconate -0.91; Phosphoenolpyruvate +0.98; FAICAR -0.89; Glucose-6P +0.83 |
| muscle | baseline | 116 | 17 | 1 | 2 | 3 | Hydroxyphenyllactic acid -1.26; Citric acid -1.08; S-adenosyl methionine -0.96; Glycine -0.85; CoA-Glutathione -0.87; Alpha-Ketoglutarate -0.80 |
| blood | pem_change | 83 | 5 | 0 | 0 | 1 | Pyroglutamic acid +0.75; Asparagine -0.73; Allantoin +0.69; Gluconate +0.61; 2-PY -0.59; 2-Hydroxybutyric acid +0.55 |
| muscle | pem_change | 116 | 3 | 0 | 0 | 0 | Asparagine -0.90; Dihydroxyacetone-P -0.81; Methionine -0.63; Uric acid -0.58; Leucine -0.53; Phenylalanine -0.53 |

Expected false positives among raw p < 0.05 at the null: 4.2 of 83 (plasma), 5.8 of 116 (muscle).

### (g) Published claims vs our computation (kept separate)

Published = Appelman 2024 Results text (GLMM over three timepoints, Box-Cox, BH within pathway). Ours = baseline-only Hedges g on the released values with a stratified bootstrap CI; for pathway claims a per-person domain score (mean z-scored log2 metabolite in the authors' domain).

| tissue | published claim | item | our g (95% CI) | share of domain in claimed direction | concordance |
|---|---|---|---|---|---|
| blood | glycolytic metabolites in venous blood significantly higher | Glycolysis | +0.21 (-0.38 to +0.80) | 0.50 | same direction, CI includes 0 |
| blood | pyruvate lower | Pyruvate | -0.98 (-1.61 to -0.46) | - | same direction, CI excludes 0 |
| blood | other TCA-cycle metabolites lower | TCA | -0.32 (-0.93 to +0.24) | 0.71 | same direction, CI includes 0 |
| blood | various purine-pathway metabolites lower | Purines | -0.35 (-0.93 to +0.20) | 0.73 | same direction, CI includes 0 |
| muscle | TCA-cycle metabolites incl. glutamate lower | Glutamate | -0.57 (-1.05 to -0.07) | - | same direction, CI excludes 0 |
| muscle | TCA-cycle metabolites incl. FAD+ lower | FAD | -0.75 (-1.33 to -0.23) | - | same direction, CI excludes 0 |
| muscle | TCA-cycle metabolites incl. alpha-ketoglutarate lower | Alpha-Ketoglutarate | -0.80 (-1.43 to -0.25) | - | same direction, CI excludes 0 |
| muscle | TCA-cycle metabolites incl. citric acid lower | Citric acid | -1.08 (-1.81 to -0.56) | - | same direction, CI excludes 0 |
| muscle | key metabolites of the TCA cycle lower | TCA | -0.65 (-1.14 to -0.18) | 0.82 | same direction, CI excludes 0 |
| muscle | glycolytic metabolites displayed few group differences | Glycolysis | +0.25 (-0.34 to +0.83) | - | concordant (CI includes 0) |
| muscle | creatine concentrations lower | Creatine | -0.44 (-0.98 to +0.16) | - | same direction, CI includes 0 |
| muscle | lower S-adenosylmethionine (SAM) | S-adenosyl methionine | -0.96 (-1.84 to -0.33) | - | same direction, CI excludes 0 |
| muscle | lower hydroxyphenyl acetic acid (NAME MISMATCH: file has hydroxyphenyllactic acid) | Hydroxyphenyllactic acid | -1.26 (-1.93 to -0.74) | - | same direction, CI excludes 0 |
| muscle | ratio of citric acid to lactate lower | Citric acid/Lactate | -0.79 (-1.36 to -0.27) | - | same direction, CI excludes 0 |
| muscle | many amino acids not different at rest | Aminoacid | -0.29 (-0.89 to +0.25) | - | concordant (CI includes 0) |
