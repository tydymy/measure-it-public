# Harmonization audit (person_concepts)

Generated 2026-10-07T18:29:38+00:00 by `uv run python -m measure_it.harmonize.audit --with-synthetic --verify-codes`. Every number below is measured from the processed tables; contract: [HARMONIZATION_CONTRACT.md](HARMONIZATION_CONTRACT.md). Per-variable detail (aggregate, NHANES): `results/tables/harmonization_variables.csv`; per-dataset coverage: `results/tables/harmonization_coverage.csv`.

## 1. Rows and participants

| dataset | rows | participants | DEMOG | ICD10CM | LOINC | RXNORM | SURVEY | DEVICE | OMICS |
|---|---|---|---|---|---|---|---|---|---|
| nhanes | 1,640,811 | 19,931 | 39,862 (2 keys) | 294,884 (330 keys) | 1,154,313 (87 keys) | 28,994 (638 keys) | 122,758 (13 keys) | 0 (0 keys) | 0 (0 keys) |
| nhanes0306 | 956,417 | 20,470 | 40,940 (2 keys) | 187,238 (16 keys) | 604,235 (43 keys) | 25,967 (631 keys) | 98,037 (11 keys) | 0 (0 keys) | 0 (0 keys) |
| byod_audit_synth_ehr | 14,926 | 200 | 400 (2 keys) | 246 (8 keys) | 1,202 (5 keys) | 126 (8 keys) | 0 (0 keys) | 1,752 (9 keys) | 11,200 (80 keys) |
| byod_audit_synth_maestro | 26,971 | 280 | 560 (2 keys) | 613 (9 keys) | 2,530 (10 keys) | 0 (0 keys) | 3,157 (12 keys) | 4,960 (20 keys) | 15,151 (67 keys) |

Keys are feature keys at the level `feature_matrix` uses (ICD-10-CM 3-character categories).

## 2. Mapping coverage, NHANES

| measure | nhanes | nhanes0306 |
|---|---|---|
| Self-report condition items mapped to an ICD-10-CM category | 22/27 | 21/26 |
| 'Yes' answers on mapped items (share of all 'yes' answers) | 78.3% | 75.0% |
| Prescription reason-for-use ICD-10-CM rows (native) | 11,903 | 0 |
|   of which well-formed ICD-10-CM (kept) | 100.0% | n/a |
| Core lab variables mapped to LOINC (participant_labs) | 61/62 | 38/38 |
| Core lab result rows mapped | 98.3% | 100.0% |
| Core lab rows dropped for a unit outside the config | 0 | 0 |
| Extended lab variables mapped (participant_labs_extended) | 23/663 | n/a |
| Extended lab result rows mapped | 7.7% | n/a |
| Exam measures mapped (BP, pulse, anthropometry) | 7/7 | 7/7 |
| Distinct drug names (components) mapped to an RxNorm ingredient | 695/752 | 688/736 |
| Medication records (components) mapped | 99.2% | 98.1% |
| SURVEY items (PHQ-9, CDC Healthy Days, PFQ / DLQ / SLQ) | 13 | 11 |

Unmapped by design (counted, not guessed):

* nhanes condition: 5 variables / names, 45,485 rows; largest: `dx_cancer` (11,323: 'cancer or malignancy' spans C00-D49; the site item is not part of the); `dx_liver_condition` (11,312: 'any liver condition' spans K70-K77, B15-B19 and more); `dx_arthritis` (11,306: 'arthritis' spans M00-M25 and more; the type item (MCQ195) is mapped i); `dx_thyroid_problem` (11,302: 'thyroid problem' spans E00-E07 (hypo-, hyperthyroidism, goitre, nodul); `arthritis_type=other` (242: 'other arthritis' has no category)
* nhanes drug: 57 variables / names, 238 rows; largest: `INSULIN ISOPHANE` (39: ambiguous: 2 ingredients); `ANTI-INFECTIVES - UNSPECIFIED` (30: no exact or normalised RxNorm match); `INSULIN` (23: class name; the specific ingredient is not recorded); `PENICILLIN` (14: no exact or normalised RxNorm match); `ANTIHYPERTENSIVE AGENTS - UNSPECIFIED` (14: no exact or normalised RxNorm match); `ANTIHYPERLIPIDEMIC AGENTS - UNSPECIFIED` (10: no exact or normalised RxNorm match)
* nhanes measurement: 1043 variables / names, 3,807,877 rows; largest: `LBXHA` (15,868: categorical result); `LBXHBS` (15,853: categorical result); `LBDHEG` (14,722: categorical result); `LBXHBC` (14,722: categorical result); `LBDHEM` (14,722: categorical result); `LBDHBG` (14,717: categorical result)
* nhanes0306 condition: 5 variables / names, 40,322 rows; largest: `dx_cancer` (10,004: 'cancer or malignancy' spans C00-D49; the site item is not part of the); `dx_liver_condition` (9,999: 'any liver condition' spans K70-K77, B15-B19 and more); `dx_arthritis` (9,996: 'arthritis' spans M00-M25 and more; the type item (MCQ195) is mapped i); `dx_thyroid_problem` (9,995: 'thyroid problem' spans E00-E07 (hypo-, hyperthyroidism, goitre, nodul); `arthritis_type=other` (328: 'other arthritis' has no category)
* nhanes0306 drug: 48 variables / names, 521 rows; largest: `99999` (100: not a drug name (NHANES unknown / refused code)); `INSULIN ISOPHANE` (92: ambiguous: 2 ingredients); `PENICILLIN` (69: no exact or normalised RxNorm match); `ANTI-INFECTIVES - UNSPECIFIED` (48: no exact or normalised RxNorm match); `INSULIN` (35: class name; the specific ingredient is not recorded); `77777` (24: not a drug name (NHANES unknown / refused code))

## 3. User (BYOD) datasets

The two user datasets below are the SYNTHETIC examples (`examples/byod_synthetic/make_example.py`: `ehr` = ME/CFS with EHR labs, diagnoses and medications; `maestro` = MAESTRO-like layers without an EHR extract), ingested into a temporary folder for this audit and removed afterwards. They show the mechanics, not the coverage of a real study.

* `byod_audit_synth_ehr`: 14,926 rows, 200 participants; mapping methods {"native": 14726, "rule": 200}
* `byod_audit_synth_maestro`: 26,971 rows, 280 participants; mapping methods {"native": 23548, "curated_map": 2530, "self_report_map": 613, "rule": 280}

## 4. Overlap of feature keys (what a phenotype built in one dataset can be evaluated on)

| comparison | shared keys | by vocabulary | keys (user datasets) |
|---|---|---|---|
| nhanes AND nhanes0306 | 562 | {'RXNORM': 492, 'LOINC': 41, 'ICD10CM': 16, 'SURVEY': 11, 'DEMOG': 2} |  |
| byod_audit_synth_ehr AND nhanes | 19 of 23 non-device keys | {'RXNORM': 8, 'ICD10CM': 6, 'LOINC': 3, 'DEMOG': 2} | `DEMOG:age_mid`, `DEMOG:female`, `ICD10CM:E03`, `ICD10CM:F41`, `ICD10CM:I10`, `ICD10CM:J45`, `ICD10CM:K58`, `ICD10CM:R53`, `LOINC:2345-7`, `LOINC:3016-3`, `LOINC:718-7`, `RXNORM:10582`, `RXNORM:36437`, `RXNORM:4452`, `RXNORM:5640`, `RXNORM:6809`, `RXNORM:6918`, `RXNORM:6963`, `RXNORM:8787` |
| byod_audit_synth_ehr AND nhanes0306 | 16 of 23 non-device keys | {'RXNORM': 8, 'LOINC': 4, 'DEMOG': 2, 'ICD10CM': 2} | `DEMOG:age_mid`, `DEMOG:female`, `ICD10CM:I10`, `ICD10CM:J45`, `LOINC:1988-5`, `LOINC:2276-4`, `LOINC:2345-7`, `LOINC:718-7`, `RXNORM:10582`, `RXNORM:36437`, `RXNORM:4452`, `RXNORM:5640`, `RXNORM:6809`, `RXNORM:6918`, `RXNORM:6963`, `RXNORM:8787` |
| byod_audit_synth_ehr AND nhanes AND nhanes0306 | 14 of 23 non-device keys |  | `DEMOG:age_mid`, `DEMOG:female`, `ICD10CM:I10`, `ICD10CM:J45`, `LOINC:2345-7`, `LOINC:718-7`, `RXNORM:10582`, `RXNORM:36437`, `RXNORM:4452`, `RXNORM:5640`, `RXNORM:6809`, `RXNORM:6918`, `RXNORM:6963`, `RXNORM:8787` |
| byod_audit_synth_maestro AND nhanes | 16 of 33 non-device keys | {'LOINC': 7, 'ICD10CM': 7, 'DEMOG': 2} | `DEMOG:age_mid`, `DEMOG:female`, `ICD10CM:A69`, `ICD10CM:E03`, `ICD10CM:F41`, `ICD10CM:G43`, `ICD10CM:I73`, `ICD10CM:J45`, `ICD10CM:Z86`, `LOINC:1751-7`, `LOINC:1884-6`, `LOINC:2085-9`, `LOINC:2093-3`, `LOINC:2160-0`, `LOINC:2345-7`, `LOINC:2571-8` |
| byod_audit_synth_maestro AND nhanes0306 | 10 of 33 non-device keys | {'LOINC': 6, 'DEMOG': 2, 'ICD10CM': 2} | `DEMOG:age_mid`, `DEMOG:female`, `ICD10CM:J45`, `ICD10CM:Z86`, `LOINC:1751-7`, `LOINC:2085-9`, `LOINC:2093-3`, `LOINC:2160-0`, `LOINC:2345-7`, `LOINC:2571-8` |
| byod_audit_synth_maestro AND nhanes AND nhanes0306 | 10 of 33 non-device keys |  | `DEMOG:age_mid`, `DEMOG:female`, `ICD10CM:J45`, `ICD10CM:Z86`, `LOINC:1751-7`, `LOINC:2085-9`, `LOINC:2093-3`, `LOINC:2160-0`, `LOINC:2345-7`, `LOINC:2571-8` |

## 5. Code verification

LOINC (NLM Clinical Tables, exact LONG_COMMON_NAME match): 105/105 verified; 
ICD-10-CM: 64/64 curated codes found in the CMS FY2026 code list (data/processed/ontology_icd10cm_codes).

## 6. Known limits

* Self-reported conditions map to ICD-10-CM categories at most at confidence `medium`; several NHANES items (any thyroid problem, any liver condition, cancer, arthritis of unspecified type) have no single category and stay unmapped.
* Feature keys roll ICD-10-CM up to 3-character categories (contract), which can join unrelated codes: e.g. `ICD10CM:Z86` holds both Z86.16 (history of COVID-19, self-reported in a user dataset) and Z86.73 (history of stroke, NHANES self-report). Consumers that need precision use condition_level='code'.
* NHANES top-codes age (80 in 2011-2014, 85 in 2003-2006): their top band is `80+`, wider than the contract's 10-year bands.
* Prescription reason-for-use codes exist only in NHANES 2011-2014; categories that come only from them are missing (NaN), not absent, in 2003-2006.
* Environmental-chemical, infectious-serology (categorical) and microbiome laboratory layers are not mapped to LOINC: no clinical EHR counterpart is needed for the phenotype use case, and guessing codes for them is not allowed.
* Survey instruments of user datasets (COMPASS-31, DSQ-PEM, FUNCAP27, ...) have no NHANES counterpart except the weak crosswalks in configs/harmonize_survey_crosswalk.yaml; SURVEY overlap with NHANES is therefore small by construction.
* NMR-derived LOINC values (Nightingale) and routine clinical-chemistry values share codes after an explicit unit conversion but are not calibrated against each other.
