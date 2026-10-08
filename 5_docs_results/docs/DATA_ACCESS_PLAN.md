# Data access plan: people with wearables plus labs, EHR or omics, and invisible-illness labels

_Written 2026-09-24. Part A checks what the Snyder lab (Stanford) has released. Part B plans access to the controlled cohorts. Every count has a source. **verified** means we read or computed it on 2026-09-24 from the file or official page named. **unverified** means it comes from a secondary report, or we could not open the primary source. Any **estimate** is labelled as one and shows its arithmetic. Raw evidence is under `data/raw/snyder_candidates/` and `data/raw/controlled_cohorts/`._

## Summary

* **Only one open, no-registration release puts a wearable, clinical labs and multi-omics on the same named people: 22 people.** The 57 people in the Glucotypes study (Hall et al. 2018, PLoS Biol) wore a continuous glucose monitor (CGM), and the paper's supplement gives clinical labs for all of them. Of the 57, 22 have a subject ID (`1636-69-001` → `69-001`) that appears in the `SubjectID` column of the iPOP/HMP `clinical_tests.txt`. That column maps to the iPOP random ID (for example `ZOZOW1T`). All 22 also have clinical labs, plasma proteome, metabolome, cytokines, gut 16S and lipidome in the iPOP public bucket. We checked the crosswalk against the data: steady-state plasma glucose (SSPG, an insulin-resistance measure) is identical for 16 of the 17 people who have SSPG in both files (r = 0.97), and BMI agrees. The only wearable here is a CGM. There are no heart-rate or activity data for these people, and no invisible-illness labels.
* **The Basis smartwatch data (Li 2017; Dunn 2021) are openly downloadable.** The tar is 20.4 GB at `http://ipop-data.stanford.edu/wearable_data/Stanford_Wearables_data.tar`. It holds 43 people's second-level heart rate, accelerometry and skin temperature, plus one person (Subject 1, M.P. Snyder) with 73 dated clinical-lab visits. **The key that links the Basis files to iPOP participant IDs (and so to labs and omics) is not public.** Dunn 2021 joined the two inside a private `SECURE_data` directory. The only open wearable-plus-lab table for more than one person is `SSPG_model.xlsx`: 20 people, with daytime and night-time heart rate, steps, BMI and baseline SSPG.
* **dbGaP phs001719 (iPOP) holds no wearable data.** Its phenotype tables have 8 subject variables and 3 sample variables for 111 consented subjects.
* **The COVID wearable studies release heart rate, steps and labels only.** The phase-2 (Alavi 2022) zip has no surveys, even though 2,122 participants filled in at least one. Uwakwe 2025 (Stanford Digital Repository) releases self-reported symptoms and comorbidities for 126 people, of whom 31 have a Long COVID label. None of these releases has labs or omics.
* **Invisible-illness labels in the Snyder releases:** Long COVID in Uwakwe 2025 (31 of 126, self-reported). The same file has a self-reported "lyme" comorbidity flag, set for 1 of the 42 people who have a non-missing value. There are no ME/CFS, POTS, fibromyalgia, EDS or MCAS labels in any Snyder release. The iPOP releases carry insulin-resistance and diabetes status plus infection, immunisation and other event annotations only.
* **For controlled access, All of Us is the clear first choice.** Its public Data Browser (CDR "2025q4r5", 747,040 participants) reports 68,840 participants with Fitbit data. It also reports these EHR counts (rounded to 20): U09.9 3,200, G93.32 1,280, G93.3 (2022 parent code) 2,220, G90.A 1,300, M79.7 15,440, Q79.6 1,600, Q79.62 500, D89.4 660, A69.2 2,900, G43 51,700, K58 21,940, K31.84 5,800, N80 6,900. Counts of people who have both Fitbit data and one of these diagnoses are not public. They require the Registered Tier (Cohort Builder). The estimates below are labelled as such.

---

## Part A — Snyder lab / Stanford releases

### A.1 What is in the public iPOP Google Cloud bucket (verified)

We listed `gs://gbsc-gcp-project-ipop_public` anonymously through the JSON API: 14,291 objects, saved to `data/raw/snyder_candidates/ipop_public_bucket/gcs_bucket_full_listing.json`. No account or agreement was needed.

| Prefix | What it is (paper) | Per-person layers | Wearable? | People (computed from the files) | Verified |
|---|---|---|---|---|---|
| `HMP/` | iPOP / iHMP prediabetes cohort (Zhou et al. 2019 Nature) | `clinical_tests.txt` (51 clinical-lab analytes per visit), `proteome_abundance.txt`, `metabolome_abundance.txt`, `cytokine_abundance.txt` (62-plex), `gut_16s_abundance.txt`, `nares_16s_abundance.txt`, `RNAseq_abundance.txt` (70.8 MB, header not parsed), raw MS files (5.4 TB, not downloaded) | **no** | clinical 104 subjects / 969 visits; proteome 106 / 950; metabolome 102 / 981; cytokine 102 / 967; gut 16S 96 / 855; nares 16S 105 / 837; clinical ∩ proteome ∩ metabolome ∩ cytokine = **97**; union 106 | verified |
| `Lipidomics/` | iPOP lipidome (Hornburg et al. 2023) | demographics (SSPG, insulin-resistance status, BMI, sex, ethnicity, statin use), clinical upload, lipidome, cytokines | no | 112 subjects / 1,649 samples; 101 overlap HMP | verified |
| `Exercise/` | Stanford exercise (CPET) multi-omics (Contrepois et al. 2020 Cell) | proteome, metabolome, lipidome, transcriptome abundance | no | proteomics matrix 175 samples | verified (file header) |
| `PHD/PHD-paper-cohort-example-data.zip` | Personal Health Dashboard (Bahmani et al. 2021 Nat Commun) | 5-minute heart rate, steps, sleep for the first 90 days; sex; IR vs IS label (IR = insulin-resistant, IS = insulin-sensitive, set by SSPG) | **yes** | 33 participants (23 with an IR/IS label), 753,984 rows; the IDs are random 4-character codes with no link to iPOP | verified |
| `metabolic_subphenotype_db/…zip` | Metwally et al. 2024 Nat Biomed Eng | CGM (clinic and home) plus venous glucose; OGTT, HbA1c, fasting glucose; SSPG, disposition index, hepatic insulin resistance, GLP-1/GIP, HOMA, T2D polygenic score | **yes (CGM)** | 61 subjects; home CGM 29, clinic CGM 28; gold-standard metabolic phenotypes 32; the IDs `S01…` do not link to iPOP | verified |
| `NBME-22-1010C/` | earlier version of the same CGM dataset | same | yes (CGM) | 71 subjects in the demographics file | verified |
| `COVID-19/`, `COVID-19-Phase2/` | Mishra 2020, Alavi 2022 | heart rate, steps (phase-1 sleep) | yes | already ingested; see `data/raw/stanford_longcovid_wearables/DATA_AUDIT.md` | verified |

### A.2 Paper by paper

| Release | Data-availability statement (quoted or paraphrased from the paper) | What we found by following every link | Same-person layers and counts | Wearable in the open release? | Invisible-illness labels |
|---|---|---|---|---|---|
| **Li et al. 2017 PLoS Biol** "Digital health: tracking physiomes and activity…" | "These data are available at http://hmpdacc.org/data/wearable/stanford.tar" | That URL returns **404**. The Wayback CDX index shows it has been 404 since 2020-06-03, and the archived copy is truncated: the server replays a fixed 512-byte range. The **same file** (identical size of 20,377,763,840 bytes, Last-Modified 2016-11-30) is **live** at `ipop-data.stanford.edu/wearable_data/Stanford_Wearables_data.tar` (Dunn 2021 link). S1 Table (s009.xlsx) has only aggregates: 43 people with Basis data, 20 with SSPG. | Basis heart rate, skin temperature and steps for 43 people (S1 Table). Clinical labs only for Participant #1 (73 visits). SSPG for 20 people. | yes (in the tar) | Lyme disease detected in Participant #1 (case report, n = 1). Not a label column. |
| **Dunn et al. 2021 Nat Med** "Wearable sensors enable personalized predictions of clinical laboratory measurements" | "Intel Basis watch data are available on the Stanford iPOP site (…/Stanford_Wearables_data.tar) and in the Digital Health Data Repository (GitHub)" | We walked the tar headers with HTTP Range reads (158 tar entries) and saved the listing to `ipop_wearables_dunn2021/tar_listing_live_range.json`. It holds 48 `Basis_Watch_Data/Basis_NNN[b].csv` files (IDs 001–043, five of them with a second "b" file, 20.1 GB in total). The columns are `time, hr, accel_magnitude, skin_temp` at 1-second resolution; we checked the first 4 KB of Basis_001 and Basis_043. There is also `Subject1_rawdata.zip` (186 MB: Basis B1 and Peak, Masimo SpO₂, Scanadu, iHealth, Withings, MOVES, radiation, and `Subject1_Clinical_record.txt` = **73 visits × 38 lab analytes**, which we extracted and checked), per-figure source data, and a **template Data Access Agreement for "sensor data from wearable devices at higher resolution"** (we read it and did not sign it). The GitHub folder is now only a README that points to `Activity/Dataset_StanfordWearables`. | 54 iPOP participants wore Basis (paper). Openly: 43 Basis IDs. **Wearable-to-lab linkage is open for only 1 person** (Subject 1: 73 lab visits) plus the 20-person summary table `Fig6/SSPG_model.xlsx` (hrDayVal, hrNightVal, Steps, Delta HR, BMI, baseline SSPG). The paper's code (`jessilyn/wearables_vitalsigns/misc/scripts/load-data.R`) reads the joined wearable + lab file keyed by `iPOP_ID` from `../SECURE_data/`, so **that key is not public**. | yes | none |
| **Schüssler-Fiorenza Rose et al. 2019 Nat Med** (iPOP, n = 109) | "Raw omics data … are hosted on the NIH Human Microbiome 2 project site (portal.hmpdacc.org) under the T2D project along with clinical laboratory data through 2016. Data from participants who have not consented to make their data public are available on dbGAP (phs001719.v1.p1)." Wearables: "Basis watch … and a Fitbit Charge 2" (Methods). | The public omics and labs are the `HMP/` files above (verified). The HMP portal copy was not separately enumerated (unverified). **dbGaP phs001719.v1.p1**: 111 consented subjects. We downloaded the public data dictionaries. The phenotype table `pht008919` has 8 variables (random ID, insulin-resistance status, SSPG (64 non-missing), FPG, sex, race, age, BMI). The sample table `pht008920` has 3 variables. **No wearable variables.** Access is through a dbGaP Data Access Request. | Clinical labs + omics: 97–106 people (A.1). Wearable: none released as per-person files keyed to iPOP IDs. | no (Basis or Fitbit data not released with iPOP IDs) | none (reported health findings are case-level: atrial fibrillation, sleep apnoea, Lyme) |
| **Hall et al. 2018 PLoS Biol "Glucotypes"** | "Raw data of continuous glucose monitoring … are available as S1 Data" | S1 Data (s010.gz): 5-minute CGM for **57** people. S5 Data (s014.db, SQLite): per-person A1C, fasting glucose, OGTT 2 h, insulin, hsCRP, lipids, SSPG, height, weight, and diagnosis (38 non-diabetic, 14 prediabetic, 5 diabetic). IDs are `1636-69-NNN` / `1636-70-NNNN`. | **22 of 57 link to iPOP** random IDs (crosswalk saved to `glucotypes_hall2018/glucotypes_to_ipop_randomid_crosswalk.csv`). All 22 have iPOP clinical labs, proteome, metabolome, cytokines, gut 16S and lipidome. SSPG matches exactly for 16 of 17 (r = 0.97). | yes (CGM) | none |
| **Bahmani et al. 2021 Nat Commun (PHD)** | "The wearables data analyzed in this study are publicly available at …/PHD/PHD-paper-cohort-example-data.zip" | Verified contents (A.1) | 33 people: wearable + IR/IS label (from SSPG) for 23 | yes | none |
| **Metwally et al. 2024 Nat Biomed Eng** (metabolic subphenotypes) | "The de-identified glucose values from CGM and venous, along with other data types … can be downloaded from … metabolic_subphenotypes_db.zip" | Verified (A.1) | 61 people: CGM + OGTT + gold-standard insulin resistance / beta-cell function + PRS (32 fully phenotyped; 29 home CGM) | yes (CGM) | none |
| **Stanford COVID-19 wearables, phase 2** (Alavi 2022 Nat Med) | "De-identified raw heart rate and steps data … can be downloaded at …/COVID-19-Phase2-Wearables.zip. Source data are provided with this paper." | The paper reports that 2,122 participants filled in at least one survey. **No survey, symptom or lab file is released.** Labels exist only for the 84 positives (Source Data). See `DATA_AUDIT.md`. | wearable + COVID test date (84) | yes | acute COVID only |
| **Uwakwe et al. 2025 PLOS Digit Health** (Long COVID, SDR druid:cb174pb4851) | Stanford Digital Repository, ODC-By | Already ingested. `Demographic_Data.csv` has self-reported comorbidity flags. | 126 people: resting HR + demographics + symptoms; 31 Long COVID | yes | **Long COVID 31 of 126**; "lyme" = 1 of the 42 people with non-missing flags (both labels are self-reported) |
| Snyderome (Chen et al. 2012 Cell), `snyderome.stanford.edu` | Lab iPOP page, "DATA access" | One person's genome, CNVs, exome and RNA reads. The data page has no wearable links. | n = 1, omics only | no | none |

**Why we did not link the Basis files to iPOP IDs by matching.** SSPG and BMI in `SSPG_model.xlsx` do match single iPOP subjects for most of the 20 people. For example, paper ID 1 (SSPG 91.5, BMI 21.7) matches `ZOZOW1T`. But nothing public says that the paper IDs 1–20 are the Basis file numbers 001–043. Joining the files on quasi-identifiers would also be an attempt to re-link de-identified data, which the lab's template agreement forbids ("will not link clinical or phenotypic data to any Human Subject"). The project's guardrails treat such a join as fabricated. **We did not make that join. We record it as "linkage possible only with the lab's key".** The Glucotypes → iPOP crosswalk is different: it rests on the published subject-ID field (`69-001` in both releases). We confirmed it, not inferred it.

**What to ask the Snyder lab for** (email `mpsnyder@stanford.edu`; the iPOP portal says RNA-seq is "only available upon request"): (1) the `Basis_NNN` → iPOP random-ID key, which would put 43–54 people with minute-level wearable data together with about 100 people's worth of longitudinal labs and omics, (2) the Fitbit Charge 2 data from the later iPOP period, and (3) the higher-resolution sensor data under their template agreement. None of this gives invisible-illness labels. Its value is as a **multi-omics and wearable method-development cohort** of healthy and prediabetic people with infection and immunisation events.

---

## Part B — Access plan for controlled cohorts

### B.0 Ranked plan

The ranking is by the largest **verified** number of people who have a wearable, another personal data layer and an invisible-illness label, and by how quickly we can get access.

| Rank | Cohort | Wearable (verified n) | Other layers in the same people | Largest verified invisible-illness label counts | Access route / tier / cost / time | Why this rank |
|---|---|---|---|---|---|---|
| 1 | **All of Us** | Fitbit **68,840** (Data Browser; CDR v9 ">68,000") — minute-level HR, steps, sleep stages; no HRV/SpO₂ | EHR >481k, labs, surveys, WGS >535k, proteomics ~10k, RNA-seq ~9k | EHR: M79.7 15,440; U09.9 3,200; Q79.6 1,600; G90.A 1,300; G93.32 1,280; D89.4 660 (whole cohort; Fitbit overlap = **estimate** ≈9.2%) | institutional DURA plus training; Registered/Controlled Tier; free with $300 credits, then Google Cloud; about 2 weeks after the DURA | largest verified wearable n; every target ICD-10-CM code exists; open now |
| 2 | **RECOVER-Adult (BDC)** | raw Fitbit streams **6,569** (Terra/BDC v5.p4); sub-study 6,529 enrolled | study labs, symptom surveys, 935k-aliquot biospecimen inventory, EHR (recently added) | LCRI "likely Long COVID" 20% of 11,743 infected (**estimate** ≈2,349); study-defined ME/CFS 531 | dbGaP DAR via BDC (eRA Commons); free data, compute billed; 2–6 months | best Long COVID label quality; wearable + labs in the same people |
| 3 | **UK Biobank** | wrist accelerometry **103,568** (2013–2016, pre-COVID) | NMR 479,096; Olink 53,039; biochemistry, WGS, HES, GP | G93 first occurrence 7,716 (not specific); HES G93.3 1,462 (sub-agent) | AMS application, £3k–9k; **new applications closed; RAP re-opening "not before October" 2026** | largest omics × device overlap, but closed and no post-COVID device data |
| 4 | N3C / PCORnet | **none** | EHR, labs | N3C U09.9 20,532 (2022) | DUA plus DUR, free, about 10 business days after the DUA | EHR scale only; use it for label prevalence |
| 5 | Lifelines / HUNT4 / NAKO | thigh/hip accelerometry 15k–63k | questionnaires, biobank | self-reported CFS/FM/IBS; NAKO post-COVID questionnaire | European PI / ethics / fees | geography- or collaborator-gated |
| — | Snyder lab (open) | Basis 43; CGM 57 | labs + multi-omics for 22 (CGM) | none | open / email the lab for the Basis→iPOP key | method development only |

**Action order.**

1. Start the All of Us institutional DURA now; there is little else to wait for. In the first week of access, pull the Cohort Builder counts for Fitbit × each code to replace the estimates in this document.
2. Submit the RECOVER-Adult dbGaP DAR in parallel (phs003463; 2–6 months).
3. Prepare a UK Biobank application to submit when applications re-open in "late 2026".
4. Use N3C only for prevalence and comorbidity structure.

The pipeline code does not change. Each cohort becomes a `participant_wearable_features__<cohort>` partition computed inside its enclave, and only aggregates are exported.

### B.1 All of Us Research Program (NIH) — rank 1

**Access route.** Your institution first signs a Data Use and Registration Agreement (DURA) with the Vanderbilt data centre. Each researcher then sets up Google 2-step verification, verifies identity with Login.gov or ID.me, completes the Responsible Conduct of Research training and signs the Data User Code of Conduct. That gives Registered Tier access. Controlled Tier access (genomics, plus un-shifted dates and more detail) needs an extra Controlled Tier training. Eligibility: "researchers affiliated with academic institutions, health care institutions, not-for-profit organizations, or government entities" (researchallofus.org/register, verified 2026-09-24). **Cost:** new users get "$300 in computational credits"; after that, Google Cloud billing to your own account. **Time:** the institution is "contacted via email … within two weeks" of the request. The DURA itself takes as long as your signing office takes (not published; unverified). The individual steps are listed at about 5–10 + 10–20 + 30–90 minutes. Data never leave the Researcher Workbench. Aggregate results must follow the rule of no counts of 1–20.

**What exists (verified).**

| Quantity | Value | Source |
|---|---|---|
| Participants in the current CDR (Data Browser, `2025q4r5_combined`, used by `cdrversion-used`) | 747,040 | public.api.researchallofus.org `/v1/databrowser/cdrversion-used` (saved to `data/raw/controlled_cohorts/aou_databrowser/cdrversion-used.json`) |
| Participants with Fitbit data (Data Browser) | 68,840 | `/v1/databrowser/domain-totals` and `/fitbit-analysis-results` (every Fitbit data type returns the same 68,840) |
| CDR v9 (C2025Q4R6, data cutoff 2025-01-01): EHR, Fitbit, WGS, omics | EHR >481,000; Fitbit >68,000; WGS >535,000; RNA-seq ~9,000 samples; proteomics ~10,000 samples; clinical notes >99,000 | AoU support article 50653909888788 "Our Largest Genomic Dataset: CDR version 9" |
| CDR v8 (C2024Q3R4): Fitbit participants; share with EHR + physical measurements + genomics + surveys | 59,018; "Nearly half (46%)" | Patten et al. 2026 Nat Med (PMC13278962), abstract and Methods |
| Fitbit tables | `steps_intraday`, `heart_rate_minute_level`, `heart_rate_summary` (heart-rate zones), `activity_summary`, `sleep_level`, `sleep_daily_summary`, `device`; v9 adds `sleep_daily_summary_30dayavg / _counts / _ext`, `sleep_level_short` | AoU support "Data Types and Organization" (4619151535508) and the CDR v9 article |
| HRV / SpO₂ / skin temperature from Fitbit | **not in the 7 documented tables** (only heart rate, steps, activity, sleep and device). Raw accelerometry is not available ("raw accelerometry data are not available to researchers", Patten 2026) | same (verified absence in the table list; the column-level data dictionary was not opened) |
| Participants with any EHR condition / labs and measurements | 436,280 / 439,060 | Data Browser domain-totals |

**Invisible-illness EHR counts (Data Browser; all counts rounded to 20; "count" = participants with at least one record of the source code).** Saved to `data/raw/controlled_cohorts/aou_databrowser/searchConcepts_icd_codes.json` and `source_concepts.json`.

| Condition | Code | Participants | Note |
|---|---|---|---|
| Long COVID | ICD-10-CM U09.9 | 3,200 | standard concept "Post-acute COVID-19" (OMOP ext. 705076) = 4,080 including mapped codes; SNOMED 1119304009 = 620 |
| Long COVID (pre-Oct 2021 proxy) | B94.8 | 720 | non-specific |
| ME/CFS | G93.32 | 1,280 | specific code since Oct 2022 |
| Post-viral fatigue | G93.31 / G93.39 / G93.3 (parent) | 580 / 120 / 2,220 | |
| CFS (ICD-9) | 780.71 | 1,940 | |
| "Chronic fatigue, unspecified" | R53.82 | 14,820 | **Do not use as ME/CFS.** It maps to the SNOMED CFS concept 52702003 (16,980), which inflates that concept. |
| POTS | G90.A | 1,300 | SNOMED 371073003 = 1,980 |
| Historic POTS code | I49.8 | 8,500 | "other specified cardiac arrhythmias", so mostly not POTS |
| Dysautonomia, unspecified / other | G90.9 / G90.8 | 1,640 / 460 | |
| Fibromyalgia | M79.7 | 15,440 | PFHH survey: 19,300 answered the "still seeing a doctor for fibromyalgia" follow-up (self-report) |
| EDS (all) | Q79.6 (parent) | 1,600 | Q79.60 1,260; Q79.62 (hEDS) 500; ICD-9 756.83 1,060; SNOMED 398114001 2,260 |
| Hypermobility syndrome | M35.7 | 1,240 | |
| MCAS | D89.4 (parent) | 660 | D89.40 620; D89.42 60; D89.44 (HαT) 40 |
| Lyme | A69.2 | 2,900 | PTLDS has no specific code; A69.20 2,740 |
| Migraine | G43 | 51,700 | |
| IBS | K58 | 21,940 | PFHH survey follow-up 38,300 |
| Gastroparesis | K31.84 | 5,800 | |
| Endometriosis | N80 | 6,900 | PFHH survey follow-up 20,480 |
| Chronic fatigue (self-report, PFHH) | survey | 19,260 | follow-up respondents |

**Fitbit × diagnosis overlap.** The Data Browser does not publish intersections. The Cohort Builder in the Registered Tier returns them within minutes of access. Until then, here is a **labelled ESTIMATE** that assumes Fitbit sharing is independent of diagnosis:

`overlap ≈ EHR code count × (68,840 / 747,040) = code count × 0.0922`

This gives: U09.9 3,200 → **≈295**; G93.32 1,280 → **≈118**; G90.A 1,300 → **≈120**; M79.7 15,440 → **≈1,423**; Q79.6 1,600 → **≈147**; D89.4 660 → **≈61**; K31.84 5,800 → ≈534; N80 6,900 → ≈636.

The estimate could be off in either direction. Fitbit sharers are more often women and more often in the BYOD group, and these conditions are female-predominant, so the true overlap is probably higher. But Fitbit participants have EHR data less often than average (46% have every data type), which pushes the other way. The first action after access is to replace these estimates with real Cohort Builder counts.

**Exact tables and variables to request / use.** OMOP `condition_occurrence` (condition_source_value in U09.9, G93.32, G93.31, G93.39, G90.A, G90.9, M79.7, Q79.6x, M35.7, D89.4x, A69.2x, G43.x, K58.x, K31.84, N80.x; exclude R53.82 from ME/CFS). `measurement`: CBC, CMP, ferritin, TSH, CRP, tryptase, orthostatic vitals where coded (LOINC). `drug_exposure`: ivabradine, midodrine, fludrocortisone, beta-blockers, LDN, H1/H2 blockers, cromolyn. `observation`/`person`. PFHH survey items. Fitbit `heart_rate_minute_level`, `heart_rate_summary`, `steps_intraday`, `activity_summary`, `sleep_daily_summary`, `sleep_level`, `device`. Controlled Tier `srWGS` for the EDS/HSD genetics sub-analysis.

**Proposed first analysis.** This uses the same MeasurementAdapter pipeline as the existing `WearableHeartRateAdapter` and `WearableAdapter` (docs/ADAPTERS.md). Define the index date as the first qualifying code. For each case with at least 30 valid Fitbit days in the 12 months before and after the index date, and at least 600 minutes of heart rate per day, compute night-time resting HR (00:00–06:59 from `heart_rate_minute_level` with no steps in the previous 30 minutes), daily steps, sleep duration and efficiency, and the Fitbit resting-HR series. Controls: 4:1 matched on age, sex, BMI category, device model, enrolment year and pathway (BYOD vs WEAR), with no invisible-illness code. Primary contrasts:

* (a) POTS (G90.A) vs controls: night RHR and daytime HR increment (daytime − night RHR), a proxy for orthostatic tachycardia burden.
* (b) ME/CFS (G93.32) and Long COVID (U09.9) vs controls: steps and step-count variability (a crash/PEM proxy), and within-person post-index change against each person's pre-index baseline, reusing the Stanford placebo-onset negative control.
* (c) Negative-control outcomes (fracture, as in Patten 2026) and negative-control codes (R53.82).

Use participant-grouped CV and report sex-stratified results. Store outputs as `participant_wearable_features__aou` inside the Workbench only. Export only aggregates of 20 or more.

### Lab biomarkers in All of Us

_Added 2026-09-24. Companion to `docs/PUBLISHED_LAB_EVIDENCE.md` (what the literature claims for each lab test) and to B.1 above (EHR diagnosis and Fitbit counts). Every count here was read from the public Data Browser API on 2026-09-24 by `measure_it.labs.aou_lab_concepts` (through `measure_it.http`, disk-cached). Raw responses: `data/raw/controlled_cohorts/aou_databrowser/searchConcepts_labs.json`. The rows are appended to `results/tables/controlled_cohort_counts.csv` as **verified**._

**What was queried.** `POST https://public.api.researchallofus.org/v1/databrowser/searchConcepts` with `{"query": <term>, "domain": "MEASUREMENT"` (or `"PROCEDURE"`)`, "standardConceptFilter": "STANDARD_OR_CODE_ID_MATCH", "maxResults": 25, "minCount": 1}`, the endpoint the Data Browser page itself calls. CDR `2025q4r5_combined` (747,040 participants; 439,060 with any lab or measurement). We read every concept each search returned and kept a concept only when its code measures the target analyte. The match was by vocabulary and code, never by name similarity. The sex split for tryptase comes from `GET /v1/databrowser/concept-analysis-results?concept-ids=3019420&domain-id=Measurement`.

**How to read the counts.** "Participants" is the `countValue` field: people with at least one EHR record of the concept, including records mapped to it from other codes. "Source" is `sourceCountValue`: people whose record carried this exact code. Every value is a multiple of 20, and the smallest value the API returns is 20. We therefore report 20 as **≤20** (1–20 people), in line with the programme rule of no counts of 1–20 noted in B.1. These are people who **had the test**. They are not people with an abnormal result, and not people with a diagnosis. The Data Browser does not publish test × diagnosis intersections. Those need the Registered Tier.

| Target biomarker | Concept (vocabulary code) | Participants | Source code | Comment |
|---|---|---|---|---|
| **Serum tryptase** | LOINC 21582-2 Tryptase [Mass/volume] in Serum or Plasma | **4,420** | 3,340 | Sex at birth: female 3,460, male 960, other 40. Also SNOMED 121873004 "Tryptase measurement" 900 (standard-only) and LOINC 51834-0 [Moles/volume] 40 |
| Urinary leukotriene E4 | LOINC 33343-5 LTE4/creatinine, urine | 60 | ≤20 | LOINC 33344-3 (LTE4 mass/volume) ≤20; LOINC 101115-4 (24-h urine) ≤20 |
| Urinary 2,3-dinor-11β-PGF2α | LOINC 94381-1 /creatinine, 24-h urine | 60 | 40 | LOINC 97658-9 (random urine) 40; urinary PGD2 LOINC 14054-1 60 |
| Urinary N-methylhistamine | LOINC 44340-8 /creatinine, 24-h urine | 280 | 200 | LOINC 13781-0 (random urine) 80; 26053-9 and 12714-2 ≤20 each |
| Plasma histamine | LOINC 2416-6 | 260 | 160 | plasma histamine is not one of the consensus MCAS mediators; listed for completeness |
| KIT D816V | LOINC 88519-4 KIT c.2447A>T [Presence] | 100 | 100 | LOINC 55201-8 KIT targeted mutation analysis 100; CPT 81273 (KIT D816 analysis) 60. The EHR does not record the sample type (blood or marrow) or the assay sensitivity |
| TPSAB1 copy number (hereditary α-tryptasemia) | none | — | — | the search "TPSAB1" returns **no** concept |
| **Plasma norepinephrine** | LOINC 2666-6 | **1,080** | 840 | posture-specific codes: supine LOINC 1601-4 **≤20**, standing LOINC 17368-2 **≤20**; LOINC 14852-8 (moles) ≤20; standing epinephrine 95054-3 ≤20 |
| Plasma catecholamines (panel) | LOINC 2056-0 | 440 | 360 | 3-panel LOINC 34551-2 60; free fractionated LOINC 42493-7 40 |
| 24-h urine norepinephrine | LOINC 2668-2 | 1,740 | 1,260 | mostly ordered as a phaeochromocytoma work-up, not as an orthostatic test |
| **Cortisol** | LOINC 2143-6 serum/plasma | **20,160** | 16,700 | AM specimen LOINC 9813-7 6,680; 24-h urine free cortisol LOINC 2147-7 2,220; saliva LOINC 2142-8 920 |
| CGRP | none | — | — | the searches "calcitonin gene" and "CGRP" return **no** concept. "Calcitonin" (LOINC 1992-7, 800) is a different hormone, used in thyroid-cancer work-up |
| Adrenergic / muscarinic GPCR autoantibodies | none | — | — | "muscarinic" returns nothing. "adrenergic receptor" returns only an ADRA1B germline-variant concept and CPT 81401, and neither is an antibody assay |
| Ganglionic AChR antibody (autoimmune autonomic ganglionopathy) | LOINC 42233-7 | 920 | 640 | LOINC 94694-7 (immunoassay) 200. This is the only autonomic autoantibody with a concept, and it is not a POTS GPCR antibody |
| Blood / plasma volume | CPT 78110 / 78111 / 78122 (PROCEDURE) | ≤20 / 40 / 40 | same | radiolabelled volume-dilution procedures. CPT 78121 (red-cell volume) 60 |
| Orthostatic vitals and tests (context) | LOINC 69001-6 HR standing / 68999-2 HR supine; CPT 93660 tilt table; CPT 95924 autonomic testing | 500 / 460; 1,780; 820 | ≤20 / ≤20; 1,460; 740 | standing and supine HR are mostly mapped from other codes (source code ≤20) |

**What this means for each published biomarker.**

* **MCAS.** Tryptase (4,420 tested) is the only mast-cell test with a usable denominator. Only 660 people carry any D89.4 code (B.1), so the number of people who have both a D89.4 code and a tryptase result is at most 660, and probably far lower. The urinary mediators the consensus criteria rely on are almost absent: LTE4 60 people, 2,3-dinor-11β-PGF2α 60, N-methylhistamine 280 (24-h). Any subgroup of these (for example MCAS-coded × LTE4) will fall under the 20-person disclosure floor. They can be described, but they cannot be modelled. The **20% + 2 ng/mL** rule needs a baseline tryptase *and* an acute sample 1–4 h after a symptomatic episode. The EHR records the draw time but not the episode onset. Whether a pair meets the rule can therefore be approximated (two tryptase results close to an emergency or allergy encounter), but not verified. HαT (TPSAB1 copy number) is not an EHR measurement. The diagnosis code D89.44 exists, but only 40 people carry it (B.1). A genotype would have to come from the Controlled Tier srWGS. We have not checked whether short-read WGS can call TPSAB1 copy number reliably (unverified; published HαT studies use digital droplet PCR).
* **POTS.** Plasma norepinephrine exists for 1,080 people, but only ≤20 have a record coded as standing or supine. The hyperadrenergic criterion (upright NE ≥600 pg/mL) therefore **cannot be observed** in the EHR as coded. A position-unspecified NE cannot be read as an upright value. GPCR autoantibodies have no concept at all. Blood-volume measurement (CPT 78110/78111/78122) was done in ≤20–40 people per code. The POTS lab biomarkers in the literature are, in practice, not measurable in All of Us. The device route (Fitbit HR, B.1) and coded orthostatic vitals are more feasible.
* **Long COVID cortisol.** Serum cortisol is common (20,160; AM-specimen code 6,680). But cortisol follows a strong diurnal rhythm, is suppressed by glucocorticoid drugs, and is ordered mainly for suspected adrenal disease. The published low-cortisol finding needs time-of-draw control (use `measurement_datetime`, restricted to 06:00–10:00) and exclusion of steroid exposure (`drug_exposure`).
* **Migraine CGRP.** Not in the EHR (no concept). It would need new assays on biospecimens. We have not checked whether All of Us biospecimens can be requested for new assays (unverified).
* **ME/CFS, fibromyalgia, IBS, endometriosis.** The research biomarkers (NK cytotoxicity, metabolomic panels, anti-CdtB/anti-vinculin, salivary miRNA) are not routine clinical tests. We did not look them up in the Data Browser. CA-125 and faecal calprotectin are routine tests and could be added to `LAB_CONCEPTS` when needed. There is a precedent for routine labs: a UK Biobank study of ME/CFS blood traits replicated 9 of 14 traits in All of Us, with small effects (d 0.2–0.5; `docs/PUBLISHED_LAB_EVIDENCE.md` PLE-056). The All of Us proteomics release (~10,000 samples, B.1) is the only route to research-grade molecular markers on the same people. Which panel was used, and how many samples overlap each diagnosis code, must be checked in the Workbench.

**The ascertainment problem.** A mast-cell mediator test is ordered because a clinician already suspects mast-cell disease (or anaphylaxis, or mastocytosis). Testing is therefore not random. It depends on the suspected diagnosis, on symptoms, on access to an allergist, and on how much health care a person uses. Three consequences:

1. **Selection.** People who were tested are not representative of anyone with the condition. MCAS-coded people who were never tested cannot contribute a value. So "the tryptase distribution in MCAS" in the EHR really means "the distribution in MCAS patients whose clinician ordered tryptase".
2. **Incorporation.** Consensus MCAS criteria include a tryptase rise and urinary mediators. A D89.4x code may have been assigned *because* a test was positive, so comparing test values between coded and uncoded people partly measures the coding rule itself.
3. **Collider bias.** Being tested is caused both by the suspected condition and by clinical severity. Comparing only tested people can therefore create associations between look-alike diagnoses (POTS, EDS, IBS) and test values that do not exist in the population.

**Proposed design (Registered Tier, not yet run).**

* **Population.** All participants with ≥1 serum tryptase result (LOINC 21582-2 or its mapped concepts) = the *tested population*. The unit is the person. Baseline value = the lowest tryptase outside a ±24 h window of any emergency or anaphylaxis encounter (T78.2, T78.3, T88.6, T78.0x, T78.1).
* **Groups within the tested population:**
  * (A) MCAS-coded (D89.40–D89.43, D89.49; D89.44 HαT is analysed separately), split into tested *before* the first D89.4 code (primary; less incorporation) and tested only after;
  * (B) tested and never MCAS-coded;
  * (C) comparators among the tested: POTS (G90.A), EDS/HSD (Q79.6x, M35.7), IBS (K58.x), each without a D89.4 code;
  * (D) positive controls: systemic mastocytosis (D47.02, C96.2x) and anaphylaxis;
  * (E) negative-control codes: R53.82, fracture.
* **Outcomes.** Baseline tryptase (continuous, log scale). The share of people above 11.4 ng/mL, and above 8 ng/mL (the range that suggests HαT). Among people with a paired baseline and an acute result 1–4 h after an emergency or allergy encounter, the share whose pair meets **20% + 2**, labelled "approximate" because symptom onset is not recorded. Urinary mediators are reported descriptively, and only for cells of ≥20 people.
* **Handling ascertainment:**
  1. Report every contrast as *conditional on being tested*, never as population prevalence.
  2. Model testing itself. Fit the probability of ever being tested on age, sex, number of outpatient visits per year, allergy/immunology visits, emergency visits and the comparator codes. Use it (a) to compare **testing rates** across POTS, EDS, IBS and MCAS (clinician suspicion is itself an outcome), and (b) as inverse-probability-of-testing weights in a sensitivity analysis. The weights are only as good as the covariates, so both the weighted and the unweighted estimates are shown.
  3. Rank the comparisons by incorporation risk. "Before first code" is primary and "after first code" is secondary.
  4. Positive controls (mastocytosis must show high tryptase) and negative-control codes (R53.82 must not) check the pipeline.
  5. Report sex-stratified results. The tested population is 78% female: 3,460 of 4,420, rounded counts.
* **What it can and cannot show.** It can show whether coded MCAS patients who were tested differ in baseline tryptase from tested POTS, EDS and IBS patients, and how often clinicians test each group. It cannot estimate the sensitivity of any mediator for MCAS, because no reference standard independent of the tests exists in the EHR. It cannot validate the 20% + 2 rule, because symptom onset times are missing. Exports follow the ≥20 rule.

### B.2 NIH RECOVER-Adult on BioData Catalyst — rank 2 (largest *labelled* Long COVID cohort with wearables)

RECOVER-Adult is the only controlled resource we found that holds consumer-wearable streams, study-visit labs, a biospecimen inventory and a **research-grade Long COVID label** (the RECOVER Long COVID Research Index, LCRI, plus study-defined ME/CFS and POTS/autonomic criteria) for the same people. The label is defined at the participant level, not as an ICD code.

| Item | Value | Source | Status |
|---|---|---|---|
| Enrolled (NCT05172024) | 15,172 (dbGaP phs003463.v6.p5: 15,154 consented) | ClinicalTrials.gov API; https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003463.v6.p5 | verified |
| Digital-health sub-study | "A total of 6,529 participants enrolled … 4,122 (62%) connected a wearable device … Fitbits comprised the majority (83% of 4,122)". The study provided Fitbit Sense 2 or Charge 5 devices. The analysis includes HRV, resting HR, SpO₂, respiratory rate, steps and active minutes. | Vogel et al., JACC Adv 2026 (PMC13315523), L75–L78, L114 | verified |
| Participants with raw Fitbit Parquet streams on BDC (phs003463.v5.p4, released 2025-11-19) | 6,569 | Terra support "Resources for RECOVER" (article 16110566204443, updated 2026-03-10) | verified. This is higher than the 6,529 enrolled in the sub-study; the difference is undocumented. |
| LCRI 2024 update | 13,647 participants (11,743 with prior infection); "20%" of infected classified as likely Long COVID | "2024 Update of the RECOVER-Adult Long COVID Research Index", JAMA (PMC11862971) | verified (WebFetch of PMC page) |
| Likely Long COVID | **ESTIMATE** 0.20 × 11,743 ≈ 2,349 | arithmetic | estimate |
| Study-defined ME/CFS among infected | 531 of 11,785 (4.5%) vs 9 of 1,439 uninfected | PMC11968624 (Table: "Post-COVID-19 ME/CFS (n=531) … Infected (n=11,785)") | verified |
| Long COVID × Fitbit overlap | **ESTIMATE** if LCRI status is independent of sub-study participation: 6,569 × 0.20 ≈ 1,314. In the published Fitbit analysis, 498 of 1,475 (34%) had high LCRI. | PMC13315523 L62 | 498 of 1,475 **verified**; 1,314 is an estimate |
| Biospecimen inventory in the release | 935,596 aliquots (v3.p2) | BDC RECOVER release notes (bdcatalyst.gitbook.io) | unverified (sub-agent) |
| Tables | Answerdata (REDCap survey and visit forms, including symptoms, PASC index items, standard labs), Biospecimens, Concepts, Demographics, Visits, Fitbit (Parquet) | BDC release notes | unverified (sub-agent) |
| Omics or autoantibody results in dbGaP | not found. Pathobiology studies (20, launched March 2025) are receiving samples. | recovercovid.org news | unverified |

**Access.** Submit a dbGaP Data Access Request through BDC (phs003463 Adult, phs003461 Pediatrics, phs003768 Autopsy). You need an eRA Commons account, a Research Use Statement and a BDC Cloud Use Statement, and the PI must be at roughly tenure-track level (recovercovid.org/data). Aggregate counts are free in PIC-SURE without a login. **Time:** approval "often takes 2–6 months", then about a week for the data to appear on BDC. **Cost:** data and storage are free. Compute is billed through Terra or Seven Bridges; first-time non-commercial users can get $500 in pilot credits (not verified: search snippet only). Clinical-trial participant data (AUTONOMIC IVIG, ENERGIZE and SLEEP all record wearable outcomes) are listed as `ipdSharing: NO` on ClinicalTrials.gov and are not requestable yet.

**First analysis.** Parse the Fitbit Parquet files into the `participant_wearable_features__recover` partition (night RHR, HRV (Fitbit daily RMSSD), steps, SpO₂, sleep) using the same `WearableHeartRateAdapter` logic. Test LCRI-high vs LCRI-low, and infected vs uninfected, with participant-grouped CV. Then fit an **autonomic sub-phenotype** model: a study-defined POTS/orthostatic-intolerance symptom cluster against night-day HR increment and HRV, adjusted for age, sex, BMI and wear time. This replicates and extends Vogel 2026 as a check on the pipeline.

### B.3 N3C Enclave and PCORnet (RECOVER-EHR) — EHR scale, **no wearables**

* **N3C.** Access requires an institutional DUA, NIH IT-security and human-subjects training, and a Data Use Request per project, which the Data Access Committee decides "within 10 business days". Access is free, nothing can be downloaded, and non-US researchers get only de-identified or synthetic data. Size: "over 22.9 million patients from 84 contributing sites" (as of Aug 2024, JAMIA 2025). U09.9: 20,532 patients from 38 partners, among 5,434,528 with acute COVID-19 (as of 2022-08-10; Reese et al. eBioMedicine, PMC9769411 L79, **verified**). **No wearable data** (we checked the data overview, FAQ, Long COVID tenant and data-enhancement paper). Use it for code-level prevalence and co-occurrence of U09.9, G93.32, G90.A, M79.7, Q79.6x and D89.4x, and for labs such as tryptase and ferritin. It cannot be used for measurement validation. (Sub-agent report; the N3C and PMC pages were read by the sub-agent: **verified by sub-agent**, not re-read here.)
* **PCORnet Front Door.** Queries are distributed and results come back aggregated; simple queries are free, and turnaround is "weeks to months". No wearables.

### B.4 Other US cohorts with wearables and labels (smaller; from the sub-agent report, unverified here unless marked)

| Resource | Wearable n | Labels | Labs / omics in the same people | Access |
|---|---|---|---|---|
| Stanford Uwakwe 2025 (already ingested) | 126 | Long COVID 31 | none | open (SDR) — **verified** |
| NIH intramural PI-ME/CFS (Walitt 2024) | ActiGraph, 17 PI-ME/CFS + 21 healthy volunteers | adjudicated PI-ME/CFS | deep multi-omics, CSF, muscle | Pennsieve 356 (open, CC-BY; already in `device_candidates/pennsieve_356_mecfs_emg`) |
| OMF LIFT trial (NCT06366724, n = 160) | Garmin (steps, HRV) | ME/CFS with orthostatic intolerance; Long COVID | trial labs | participant data available under a data access agreement |
| Scripps DETECT | 37,146 enrolled; Long COVID analysis 279 vs 274 | self-report | none | no public route |
| TemPredict (Oura) | 63,153 | COVID self-report | antibody kits for 10,021 | application; Oura policy restricts third-party sharing |
| Yale LISTEN, CDC INSPIRE | no wearables | Long COVID (POTS 167 of 578 in LISTEN, self-report) | EHR via Hugo | not public / study team only |

### B.5 UK Biobank — rank 3 (largest multi-omics + accelerometry, but **closed to new applications**, and all accelerometry predates COVID)

**Status as of 2026-09-24 (verified, UKB Community Zendesk API):**

* FAQ "Can I apply to use UK Biobank?" (updated 2026-07-01): "Applications for new projects are currently not being accepted … We intend to start accepting new applications in late 2026."
* CEO post of 2026-09-18, "UKB-RAP: update on phased re-opening of access": "we do not now anticipate that re-opening access will start before October". Phased access will start first for institutions that have confirmed deletion of previously downloaded data.

Background, according to the sub-agent (UKB community posts plus press): UKB-RAP was closed on 22 April 2026 after "clear breaches of the UK Biobank Material Transfer Agreement". New data releases are paused until the "first half of 2027". Result downloads will initially need a manual output check.

**Access once it re-opens.** Register in the Access Management System, submit an application, sign the Material Transfer Agreement and pay. Analysis is on UKB-RAP only (DNAnexus). Eligibility: "all bona fide researchers … regardless of their location", commercial included. Fees before the pause (Wayback copy dated 2026-08-31, sub-agent): Tier 1 £3,000, Tier 2 £6,000, Tier 3 £9,000 for 3 years (ex VAT); students and lower-income countries £500. Compute on the RAP is billed per the 2026 rate card (for example £0.50/hr on-demand for a 16-core instance). Time: "an average of 15 weeks from application submission to data release" (pre-pause figure).

| Quantity | Value | Source | Status |
|---|---|---|---|
| Wrist accelerometry (Axivity AX3, 7 days), field 90001 | 103,568 participants / 112,330 items | https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=90001 | verified (re-read 2026-09-24) |
| Seasonal repeat accelerometry | 3,208 accepted; waves Nov 2017–Feb 2019; **no 2020s wave** | IJE 2026 (sub-agent) | unverified here |
| Olink proteomics (pharma proteomics pilot), field 30900 | 53,039 participants | showcase field 30900 | verified |
| NMR metabolomics (Nightingale), field 23474 | 479,096 participants | showcase field 23474 | verified |
| Biochemistry (e.g. cholesterol 30690) / haematology (WBC 30000) | 470,314 / 478,810 | showcase (sub-agent) | unverified here |
| WGS | 490,640 | PMC12443626 (sub-agent) | unverified here |
| Hospital-inpatient ICD-10 (field 41270): G93.3 (post-viral fatigue / ME) | 1,462 | showcase 41270 tree (sub-agent) | unverified here |
| … M79.7 fibromyalgia (sum over subcodes) | 5,072 (M79.7 alone 3,597) | same | unverified here |
| … G90 autonomic (sum) / Q79.6 EDS / A69.2 Lyme / U07.4 (post-COVID, UK usage) | 693 / 152 / 79 / 450 | same | unverified here |
| First-occurrence G93 (all sources, 3-character, not ME/CFS-specific), field 131114 | 7,716 | showcase 131114 | verified |
| Self-reported CFS / fibromyalgia (field 20002 tree) | 2,603 / 1,408 | showcase (sub-agent) | unverified here |
| U09.9, D89.4x, K31.84 | **not in WHO ICD-10**, so absent from HES coding | — | verified (coding system) |
| Health & Wellbeing (COVID) questionnaire completed | 201,684; Long COVID base cohort 8,668 with 2,751 cases | Nat Commun 2025, PMC12311030 (sub-agent) | unverified here |
| G93.3 inpatients who also have accelerometry | **ESTIMATE** 1,462 × 103,568 / ~502,000 ≈ 302 (upper bound: people who wore the accelerometer are healthier; 502k is the commonly cited cohort size) | arithmetic | estimate |
| Long COVID × accelerometry | **0 by design**: every accelerometer wear was before 2020 | — | inference from verified dates |

**First analysis (once access resumes).** Pre-diagnosis accelerometry as a *prodromal / risk* signal. Take participants with a later G93.3, M79.7 or G90.x code (hospital + primary care + self-report, requiring at least 2 sources for a case) and compare them with matched non-cases on wrist accelerometry features: daily acceleration, intensity distribution, sleep and the rest-activity rhythm. The features come from the existing `WearableAdapter` (the NHANES wrist-accelerometry pipeline, which uses the same device class and summary features). Then add NMR (n ≈ 479k) and Olink (≈53k) as a molecular-context layer. This is truly person-linked, so it may be called patient multi-omics.

Request these fields: 90001/90004/90012 plus the derived accelerometry fields (Category 1008/1009), 41270/41280, 42040 (GP), 20002, 131xxx first occurrence, 23400–23655 NMR, 30900 Olink, 30000–30890 biochemistry and haematology, and the COVID Health & Wellbeing questionnaire. This is Tier 2.

### B.6 Other cohorts (from the sub-agent report; verify before budgeting)

| Cohort | Device data | Relevant labels | Labs / omics | Access, cost and time | Rank note |
|---|---|---|---|---|---|
| **FinnGen** (Finland) | none cohort-wide | Risteys R13: fibromyalgia 4,291; G93.3 post-viral fatigue 334; G90 963; IBS 14,559; endometriosis 20,913; migraine 29,968. Finnish ICD-10 includes U09.9. | genotypes for all | individual data via FINBB Fingenious plus a Findata permit; summary statistics are public | labels and genomics only |
| **Lifelines** (NL) | ActiveLIFE thigh accelerometry in wave 4: about 15,382 of the 60,000 target collected so far (not yet in the catalogue) | Ballering 2022 COVID questionnaires, n = 76,422; self-reported CFS 1.3%, fibromyalgia 3.0%, IBS 9.7% of 94,516 (**ESTIMATE**: ≈1,229 / 2,835 / 9,168) | UGLI genotyping 38k, DAG3 microbiome/metabolomics 9.5k | application; about €3,635 for the academic Workspace (archived Apr 2026) | possible future rank 4 once ActiveLIFE is released |
| **HUNT4** (Norway) | Axivity on thigh and back: 32,644 valid | national registry linkage | biobank | Norwegian PI plus REK ethics approval | needs a Norwegian collaborator |
| **NAKO** (Germany) | hip ActiGraph: 63,236 valid | post-COVID questionnaire, 109,707 respondents | biobank | TransferHub; EU/EEA or GDPR-adequate countries only | |
| **Estonian Biobank** | none found | EHR (national health insurance) | NMR for everyone, WGS about 2.8k | ethics approval takes 3–6 months | no devices |
| **Evidation** | Homekit2020 (848 with dense Fitbit data), COVID Signals (target about 900; Garmin/Empatica) on Synapse | influenza-like illness / COVID only | none | Synapse "certified user" plus click-through Conditions for Use (**not accepted by us**) | no invisible-illness labels |
| Tohoku Medical Megabank | Fitbit Charge 5, "nearly 2,000", 1 year | records | multi-omics | not checked | watch |
| CLSA (Canada) | TicWatch plus ActiGraph | — | biobank | CAD $4,000 domestic / USD $10,000 international | |

### B.7 Guardrails for every cohort

* The only counts carried into this project are those read from a public page or computed from a downloaded file (marked verified). Every other figure is marked unverified or estimate in `results/tables/controlled_cohort_counts.csv`.
* No controlled cohort is joined to any public layer at the person level. The geographic and facility layers stay ecological.
* Where omics come from the same participant (All of Us proteomics, UK Biobank Olink/NMR, iPOP), results may be called patient multi-omics. Everywhere else they stay "condition-level molecular enrichment".
* R53.82 ("chronic fatigue, unspecified") and I49.8 are **not** used as ME/CFS or POTS labels. They serve as sensitivity or negative-control code sets.
* No accounts were created and no data-use agreements or click-through conditions were accepted for this plan. The Stanford template agreement found inside the Basis tar was read only.

### B.8 Evidence files

* `data/raw/snyder_candidates/`: iPOP bucket listing, HMP / Lipidomics / Exercise abundance and clinical files, PHD and CGM zips, Glucotypes S1/S5 data and the iPOP crosswalk, dbGaP phs001719 data dictionaries, the Basis tar listing (Range-walked), Subject 1 clinical record, SSPG tables and the template agreement (each with MANIFEST.json or PROVENANCE.json).
* `data/raw/controlled_cohorts/aou_databrowser/`: All of Us public API responses (`cdrversion-used`, `domain-totals`, `fitbit_analysis_results`, `searchConcepts_*`, `source_concepts`).
* The UK Biobank, RECOVER, N3C and European cohort figures marked "sub-agent" were read by a research sub-agent on 2026-09-24 from the URLs cited. We re-read the UK Biobank showcase fields 90001, 30900, 23474 and 131114, the UK Biobank status posts, and the RECOVER digital-health, LCRI, ME/CFS and N3C U09.9 sources ourselves.
