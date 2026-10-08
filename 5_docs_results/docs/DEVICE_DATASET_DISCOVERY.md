# Device-dataset discovery: can a public dataset show a device uncovering an invisible illness?

**Date:** 2026-09-24. **Module:** `src/measure_it/measurements/device_dataset_discovery.py`.
**Table:** `results/tables/device_dataset_candidates.csv` (90 rows, one per deduplicated candidate; every score keeps
its component columns). **Raw files:** `data/raw/device_candidates/<candidate_id>/` with `MANIFEST.json` (url, bytes,
sha256, retrieved_at); schema-level inspection output in `data/raw/device_candidates/inspection.json`.

**Reproduce:** `uv run python -m measure_it.measurements.device_dataset_discovery` (pipeline step
`measurement_device_dataset_discovery`; `xlrd` reads the two `.xls` files and `mat-io` decodes one MATLAB `table`
object; both are project dependencies in `uv.lock` since 2026-09-24).

## The short answer

* **No dataset processed here demonstrates that a device uncovers an invisible illness.** What exists is a small
  number of **open, clinically labelled case-control datasets with a device measurement** that are good enough to
  *test* the question with a pre-specified analysis. None has been tested yet.
* The best **open** candidate for the Long COVID / ME/CFS question is **MUSCLE-ME** (Charlton et al., Nat Commun
  2026). It has per-person **daily step counts** from an accelerometer and **maximal CPET**: 25 Long COVID patients
  who meet the Canadian Consensus Criteria, 26 ME/CFS patients and 30 matched healthy controls. Of these, 51 patients
  and 22 controls have step data.
* The next four are **fibromyalgia** clinic devices with clinical labels and 49-100 controls: infrared thermography
  (85 vs 92), force-plate posturography (80 vs 49) and a plantar-pressure platform (99 vs 100). A CFS bedside
  HR/BP provocation set (60 vs 20) completes the top five. The thermography authors' own published result is **null**.
* The best **wearable-HRV Long COVID** set (**COPERIA**, Polar H10 during rest, 6-minute walk and cold pressor;
  63 vs 66 post-COVID controls) needs one free Synapse login. The best-powered ME/CFS device sets (Cornell
  2-day CPET 84 vs 71 and wrist ActiGraph 58 vs 41) need one free mapMECFS login. Neither login was created.
* Every large wearable resource with a diagnosed or criteria-based label is **controlled access**. That includes
  RECOVER-Adult (Fitbit plus a 10-minute active stand in about 15k people), All of Us, UK Biobank (access
  currently suspended), Lifelines and HUNT4.
* **Omics (lead's question 1):** one open source links device and omics **within the same people**. The Appelman
  2024 Source Data has VO2max plus 83 blood and 116 muscle metabolites under shared study IDs: 23 Long COVID and
  21 healthy controls have both VO2max and baseline blood metabolomics (overlap counted, below). That is
  participant-linked omics under the project's guardrail. It is small.

Nulls come first, as the brief requires. Published nulls for devices in these conditions, found in or next to these
datasets:
* thermography in fibromyalgia (the authors' conclusion);
* resting HRV after COVID, where JMIR 2025 reuses Fit-COVID and finds PCC vs healthy not different;
* head-up tilt in LIINC Long COVID (negative in 23 of 26, no POTS on tilt);
* resting HR in the Stanford Long COVID set, this project's own null (AUROC 0.504).

## What was done, and what was not looked at

1. **Deduplicated** the 5 searchers' candidates into 90 entries. The same study reached through two routes is one
   row (for example LIINC tilt, Fit-COVID, PCH, Welltory and the Milan records). Where the same trial is published
   twice, the later record supersedes: Appelman 2024 is superseded by Charlton 2026 (both NCT05225688).
2. **Retrieved** listings, dictionaries and small per-person tables anonymously through `measure_it.download`,
   about 33 MB in total. Nothing was downloaded that was over 2 GB or was a bulk signal archive (COVISLEEP EDF,
   hEDS fMRI, mBrain21 and the ME/CFS imaging zips were checked by listing only). Access routes were probed through
   `measure_it.http`: the Synapse ACL and permissions endpoints, mapMECFS resource URLs, PhysioNet file URLs, the
   Zenodo record API, the Pennsieve download manifest and the Kaggle metadata.
3. **Inspected at schema level only.** The code counts rows, label values and people with non-missing device columns.
   It **never computed a device-by-group statistic** for any dataset, so the primary analyses below are specified
   before any outcome was seen. Things that were seen:
   * the label value counts;
   * the LIINC tilt interpretation categories, which are pooled and already published;
   * the thermography abstract's published null.

   A searcher reported that it printed PCH group mean steps. **The PCH steps analysis is therefore not blind** and
   is not recommended as a primary test.
4. **Did not** create any account, accept any data-use agreement, or request any restricted record.

Status counts (`access_status`): verified_open 51, controlled_access 23, verified_registration_only 6,
unreachable 6, not_as_described 4.

## Ranked verified-open datasets (case-control capable)

To be eligible for ranking, a dataset had to meet all of these:
* verified_open;
* a target-condition label whose code mapping is stated in the file;
* at least 20 cases and at least 20 controls with the device columns non-missing;
* not superseded by another row.

Score = label quality (0-3) + n (0-3, 1.5·log10 of the smaller group) + deployability (0-3) + completeness (0-1) +
controls (0-1) + device documented (0-1). The rubrics are in the module header.

| rank | dataset | condition | device | label | cases / controls with device | score |
|---|---|---|---|---|---|---|
| 1 | MUSCLE-ME, Charlton 2026 Source Data | Long COVID, ME/CFS | accelerometer daily steps; CPET | CCC (LC); diagnosed ME/CFS | 51 / 22 (steps); 49 / 30 (VO2) | 10.92 |
| 2 | Thermography in FM, PLOS ONE 0253281 | fibromyalgia | infrared camera, 6 regions | clinical (ACR per paper) | 85 / 92 | 10.89 |
| 3 | Force-plate posturography in FM, PLOS ONE 0196575 | fibromyalgia | force platform sway | clinical (ACR per searcher) | 80 / 49 | 10.54 |
| 4 | Foot baropodometry in FM, Mendeley pmvdx4dgxf | fibromyalgia | plantar pressure platform | clinical flag, criteria unstated | 99 / 100 | 9.99 |
| 5 | Neuromuscular strain in CFS, PLOS ONE 0159386 | ME/CFS | HR + BP during bedside provocation | Fukuda | 60 / 20 | 9.95 |
| 6 | PCH study, PLOS ONE 0352332 | Long COVID (self-report) | activPAL steps; 6MWT HR | self-report | 85 / 51 | 9.54 (not blind) |
| 7 | SLE/Sjogren Fitbit, Zenodo 8018238 | fibromyalgia inside SLE | Fitbit HR, sleep, steps | self-reported physician diagnosis | 27 / 63 (SLE controls) | 8.36 |
| 8 | EPOC, OpenNeuro ds007605 | Long COVID cognitive complaint | PVT reaction time; EEG | self-report | 60 / 37 | 8.34 |
| 9 | Fit-COVID, Mendeley t45r8yd6jd | prior COVID (not PCC) | short-term HRV | infection status | 20 / 20 | 7.95 |
| 10 | Shehadeh Long COVID cohort, Dataverse QYEBIN | Long COVID | HRV (metric unstated), VRI | clinic classification | 69 / 25 | 7.92 |
| 11-12 | migraine MRI ds005016; FM MRI ds001928 | migraine; FM | MRI | trial / clinical | 112 / 35; 20 / 20 | 7.32; 6.95 |

**Robustness of the ranking:** dropping any one score component keeps at least 4 of the 5 top datasets, and in five
of the six cases all 5. Charlton and thermography swap first place depending on weights. Removing the label-quality
component moves the self-report, non-blind PCH set to first, which is one reason label quality is weighted.

**Exact download URLs (top 5):**

1. Charlton 2026 Source Data:
   `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-75725-y/MediaObjects/41467_2026_75725_MOESM4_ESM.xlsx`
   (sheet `Main`, rows `Group == "MUSCLE-ME"`)
2. Thermography S1 Data: `https://ndownloader.figshare.com/files/28433660` (`pone.0253281.s001.xlsx`)
3. Force-plate S1 Dataset: `https://ndownloader.figshare.com/files/11317067` (S2, the strength reliability
   sheet, is at `https://ndownloader.figshare.com/files/11317070`)
4. Foot pressure `.sav`:
   `https://data.mendeley.com/public-files/datasets/pmvdx4dgxf/files/2d767323-8564-4b57-a2b4-126d91c4975f/file_downloaded`
   (needs a browser User-Agent)
5. Rowe CFS strain SPSS: `https://ndownloader.figshare.com/files/6101478` (`S1 File.SAV`)

## Pre-specified primary analyses for the top two (locked 2026-09-24, before any outcome was examined)

These are also stored verbatim in the `recommended_primary_analysis` column. No other feature, threshold or model
will be tried. Anything else will be labelled exploratory.

### 1. MUSCLE-ME (Long COVID + ME/CFS vs healthy): does a wearable step count separate patients from controls?

* **Data:** sheet `Main`, rows with `Group == "MUSCLE-ME"`. The 24 AGBRESA bed-rest volunteers are a separate
  pre-pandemic cohort and are excluded.
* **Primary contrast:** patients (Session `LC` or `ME`, n = 51 with `Steps`) vs healthy controls (Session `CON`,
  n = 22 with `Steps`).
* **Primary feature:** daily step count, one feature. The direction is fixed in advance (lower = patient) and nothing
  is fitted.
* **Statistic:** AUROC = P(control steps > patient steps), with a 95% CI from 2,000 stratified person bootstraps
  (seed 20260923).
* **Decision rule:**
  * "the wearable captures the illness-associated functional limitation" if the CI lower bound is at least 0.70;
  * "null" if the CI includes 0.50;
  * "weak" otherwise.
* **Missing data:** 8 of 30 controls have no steps. Report the complete-case result, plus worst-case and best-case
  bounds that assign those 8 to the lowest and the highest observed step value.
* **Pre-specified secondaries, all reported whatever they show:**
  * (a) LC vs CON and ME/CFS vs CON separately;
  * (b) `VO2_rel` (laboratory CPET) AUROC by the same method, as the lab comparator;
  * (c) LC vs ME/CFS by steps, a specificity check expected to be near 0.5.
* **Interpretation, fixed now:** activity limitation is part of the case definition (incorporation bias). A high
  AUROC would show that a consumer-grade sensor measures the defining limitation in *diagnosed* patients vs
  *healthy* people. It would not show that the sensor detects undiagnosed illness, or that it tells these illnesses
  apart from deconditioning or depression. The patients were also well enough to complete maximal CPET, so they sit
  at the milder end.

### 2. Thermography in fibromyalgia: does skin temperature add information beyond age and BMI?

* **Data:** S1 Data, 178 rows. Group 1 = FM (IDs `FM###`) and Group 2 = control (IDs `CG###`), as verified from the
  ID prefixes. Complete cases on the six regional **average** temperatures (`Neck_ave`, `Upper_back_ave`,
  `Lower_back_ave`, `Chest_ave`, `Knee_ave`, `Elbow_ave`): 85 FM vs 92 controls. The `*_min` and `*_max` columns
  are not used, which removes the option of choosing among 18 features.
* **Models:** (1) age + BMI; (2) the six `*_ave` temperatures; (3) both together. All are L2 logistic regressions
  with C = 1.0, standardised inside each fold.
* **Evaluation:** stratified 5-fold CV repeated 20 times (seed 20260923). The primary statistic is **Delta AUROC =
  model 3 minus model 1**, computed on pooled out-of-fold predictions, with a 95% CI from a paired person bootstrap
  (2,000 resamples). Device-only AUROC is tested against a 1,000-permutation label-shuffle null.
* **Decision rule:** "thermography adds information" only if the Delta AUROC CI lower bound is above 0 **and** the
  permutation p is below 0.05. Otherwise the result is reported as a replication of the authors' published null.
* **Caveat:** room temperature and season are not recorded, and that caveat goes with any result.

Suggested analyses for ranks 3-5 are in the CSV. They are not locked.

## What the verification changed (searcher claims that did not hold, or were sharpened)

| candidate | searcher claim | verified |
|---|---|---|
| al-Andalus FM ActiGraph (figshare 10013066) | 607 FM vs 39 non-FM controls | **not_as_described.** The 39 are recruits *not diagnosed by a rheumatologist*, and **all 39 have more than 10 tender points**. The authors' `.sps` excludes them as undiagnosed. There are no controls, so the set is case-only (536 FM with valid GT3X). |
| Robinson CFS HRV (PLOS ONE 0210394) | HRV for 51 CFS; controls unclear | **not_as_described.** The HRV sheet has no group or ID column, so rows cannot be assigned. |
| "Postural control ... Long COVID" (Zenodo 17563204) | likely simulated | **Confirmed.** FSS runs from -0.61 to 8.08 (scale 1-7), SF-36 PF reaches 105.41 and pain interference goes negative. Excluded. |
| FM proprioception (Zenodo 15188319 / 11085545) | likely simulated | **Confirmed and stronger.** The Symptom Severity score (an integer from 0 to 12) is non-integer and reaches 13.78 and 14.12; WPI is non-integer; the two records disagree. Excluded. |
| FM/CFS postural (PLOS ONE 0195111 `.mat`) | unreadable | **Decoded** with mat-io: 75 rows x 33 force-plate columns, 25 per `Code` 0/1/2. The code-to-group map is not in the file, so the set is not rankable. |
| PCH (PLOS ONE 0352332) | label column not found (one searcher) | `KlachtengerelateerdaanCOVID`: 87 PCC, 52 non-PCC, 1 missing. The set is **not blind** (see above). |
| Baroreflex LC (OSF hb6mn) | 16 folders (one searcher) / 52 (other) | The API lists **52** `SUB-xx` folders. No label file exists, so the authors' key is needed. |
| Shehadeh Long COVID (Dataverse QYEBIN) | HRV + VRI, 73 vs 41 | Counts confirmed (HRV for 69 cases and 25 controls). The README is an **unfilled template**, the HRV metric is undocumented and VRI's device is unnamed. The cohort is male-majority (82 M / 30 F). |
| Charlton 2026 | code only; data route unverified | **Source Data xlsx exists** (MOESM4), with per-person rows for LC, ME/CFS, controls and bed rest. |
| COPERIA (Synapse) | registration-only | Confirmed: AUTHENTICATED_USERS READ+DOWNLOAD, PUBLIC READ only, 0 access requirements, anonymous `canDownload=False`. |
| Homekit2020 / COVID Signals | qualified-researcher route | Confirmed by ACL: DOWNLOAD is limited to approved teams (3434178, 3464047). |
| Puebla EEG | 25-26 EDF | The record lists 24 EDF files plus Summary.xlsx. |
| Palombo ME/CFS uptime | DOI -3 (one searcher) | The figshare resource DOI is `10.1186/s12967-020-02583-7`. |
| mapMECFS packages | login-gated | Re-probed today. Every data file and data dictionary returned HTTP 403 (login page). |

## Question 1: omics with person-level signal

Condition-level enrichment is what the project has so far. Linked sources found in this sweep:

* **Open, participant-linked (counted):** the Appelman 2024 Source Data
  (`https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-44432-3/MediaObjects/41467_2023_44432_MOESM4_ESM.xlsx`).
  One file carries shared study IDs across sheets: `Figure1` (VO2max, peak power, GET), `Supplementary Data
  Figure2C` (accelerometer steps, 2 weeks), `Supplementary Data Figure 6A` (muscle metabolomics, 116 metabolites)
  and `6B` (blood metabolomics, 83 metabolites) at baseline, 1 day and 1 week after induced PEM. The ID overlaps are:

  | | VO2 | steps | baseline muscle metabolomics | baseline blood metabolomics | VO2 & muscle | VO2 & blood | steps & blood |
  |---|---|---|---|---|---|---|---|
  | Long COVID | 23 | 25 | 25 | 25 | 23 | 23 | 25 |
  | healthy | 21 | 13 | 19 | 21 | 19 | 21 | 13 |

  This can support a small person-level question. Two examples: do blood metabolites add to VO2max in separating
  Long COVID from controls (23 vs 21)? Does the post-PEM change in metabolites track a device measure? It is
  enough for a feasibility estimate, not a validation. The same trial feeds Charlton 2026, whose per-subject rows
  add muscle mitochondrial respiration and fibre typing (physiology, not omics) next to CPET and steps.
* **Registration-only:** Cornell ME/CFS CRC on mapMECFS. The `cor_id` key links wrist ActiGraph (58 vs 41) to
  exercise metabolomics, EV proteomics and PBMC scRNA-seq packages. The overlap cannot be counted without a login.
* **Controlled:** UK Biobank (accelerometry + NMR metabolomics + Olink, same people; access suspended in 2026),
  Lifelines (ECG + LifeLines-DEEP), and Hoel 2026 (ME/CFS proteome open on PRIDE, steps on request).

## Registration-only and controlled access (recorded, not requested)

* **One free account away:**
  * COPERIA (Synapse syn66697262; Polar H10 HRV, 63/66 valid);
  * mapMECFS Cornell Keller 2-day CPET (84/71) and Receno ActiGraph (58/41);
  * mapMECFS NIH PI-ME/CFS neurophysiology (HRV, tilt, Valsalva) and body-composition/exercise (accelerometry, CPET,
    grip);
  * mBrain21 on Kaggle (4 participants, 48.5 GB).
* **DUA / DAC / qualified-researcher review:**
  * Homekit2020, Evidation COVID Signals and Evidation ILI (Synapse teams + Conditions for Use);
  * CovIdentify (PhysioNet credentialed; files return 403);
  * JAX / Bateman Horne NASA lean test and MCAM (Data Use Certificate + NINDS DAC);
  * RECOVER-Adult and RECOVER-Pediatrics (dbGaP DAR + DUC);
  * All of Us (institutional DURA; Workbench only);
  * UK Biobank (application + fees; suspended after the April 2026 breach per press reports);
  * Lifelines (fees from EUR 4,500 per year);
  * HUNT4 (Norwegian PI + REK);
  * N3C.
* **Restricted on Zenodo (0 files visible, owner approval needed):** the Milan POTS records 4277100, 4419890,
  20327114, 14856129, 10559512, 14831596 and 18772019; BioICOPER 14282873; and 14865269, 8277005 and 5840869.
* **On request from authors:** Covid Collab, Scripps DETECT, Ruijgt / DeSportarts Firstbeat, CDC-FHCSD Fitbit, and
  Hoel steps.
* **No route found (unreachable):** RECOVER-AUTONOMIC (not yet released), STOP-PASC, Corona-Datenspende, the BHC
  ankle UpTime study, the Workwell Long COVID arm and Visible.

## Limitations

* The label criteria for ranks 2 and 3 (ACR) come from the papers as reported by the searchers. They were not
  re-read in full here. Rank 4's flag has no named criteria.
* Every open case-control set compares diagnosed patients with **healthy** volunteers. The clinical question is
  this illness versus other causes of the same symptoms, and case-control designs overstate accuracy for it. The
  only open set with symptomatic controls is Kedor 2022: post-COVID ME/CFS vs post-COVID non-ME/CFS on standing HR
  and grip, 17 vs 18, below the ranking floor.
* The largest open sets are fibromyalgia clinic devices (thermal camera, force plate, pressure platform), not
  consumer wearables. The only open, criteria-labelled consumer-sensor measure in Long COVID / ME/CFS with 20 or
  more controls is MUSCLE-ME steps.
* The ranking rubric is a transparent convenience, not an evidence grade. The components are in the CSV.

## Open items (for whoever owns shared files)

* Resolved 2026-09-24: `xlrd` and `mat-io` are project dependencies (`pyproject.toml`, `uv.lock`).
* Resolved 2026-09-24: the module runs as pipeline step `measurement_device_dataset_discovery`. It writes no processed
  parquet, so it does not need provenance columns under `write_table`.
* No `registry_entry.yaml` or `DATA_AUDIT.md` was written per candidate. These are discovery probes, not ingested
  sources. Whichever dataset is analysed first should get one.
