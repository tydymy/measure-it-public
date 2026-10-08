# Lab-dataset discovery: open per-person biomarker data with invisible-illness labels

**Date:** 2026-09-24. **Module:** `src/measure_it/labs/lab_dataset_discovery.py`.

**Table:** `results/tables/lab_dataset_candidates.csv`. It has 182 rows, one per searcher candidate. Duplicates are
flagged rather than dropped, and every score keeps its component columns. Rank robustness is in
`results/tables/lab_dataset_candidates_rank_robustness.csv`.

**Raw files:** `data/raw/lab_candidates/<candidate_id>/`, each with a `MANIFEST.json` (url, bytes, sha256,
retrieved_at). About 260 MB in total; no single file is over 70 MB. The schema-level inspection output is
`data/raw/lab_candidates/inspection.json`. The searchers' claims are kept verbatim as the discovery input in
`data/raw/lab_candidates/_searcher_inputs/`.

**Reproduce:** `uv run python -m measure_it.labs.lab_dataset_discovery [--no-download]` (pipeline step
`labs_dataset_discovery`). `xlrd` reads the legacy `.xls` supplements (a project dependency since 2026-09-24).

The method replicates `docs/DEVICE_DATASET_DISCOVERY.md`: parallel searchers, anonymous verification of every open
claim, schema-only inspection, a transparent ranking, and a locked plan written before any group comparison.

## The short answer

* **MCAS:** there is **no open per-person dataset with urinary LTE4, 2,3-dinor-11β-PGF2α, tetranor-PGDM or
  N-methylhistamine together with MCAS labels and controls.**
  * The only per-person MCAS lab file is Vysniauskaite 2015 (PLOS ONE S1): 238 MCAS and 19 systemic mastocytosis
    patients with plasma heparin and tryptase, and N-methylhistamine and chromogranin A as categories only.
  * It has **no controls**, and its tryptase is part of the SM criteria, so it cannot test a biomarker against an
    illness label.
  * Tryptase and KIT D816V appear only in case-only mastocytosis tables of 22-115 patients.
* **POTS:** there is **no open per-person dataset with supine/upright catecholamines, GPCR autoantibodies or blood
  volume together with controls.**
  * The largest autoantibody study (Hall 2022; 116 POTS vs 81 controls; published null) is available on request.
  * The only open catecholamine table is 30 paediatric orthostatic-intolerance patients in PDF, with no healthy
    controls.
* **The open, analysable lab sets with a criteria-based or clinic label and at least 20 people per group are:**
  * ME/CFS metabolomics: Metabolomics Workbench, journal supplements, PRIDE;
  * Long COVID immune profiling: Nature source-data and supplementary tables, Zenodo, Mendeley;
  * hEDS/HSD Olink;
  * endometriosis ELISAs;
  * Lyme metabolomics;
  * small fibromyalgia, IBS and migraine sets.

  31 datasets passed every eligibility rule. They are ranked below.
* **Three pre-specified analyses were run.** Rank 1 was analysed, and ranks 2 and 3 were tied within 0.01 of rank
  score, so both were analysed.

  | dataset | primary result (95% CI) | decision |
  |---|---|---|
  | **hEDS/HSD serum Olink** (Cinquina 2026), rank 1 | CV AUROC **0.684 (0.624-0.747)**, permutation p ≤ 0.001, 88 Italian patients vs 176 Italian controls | rule met |
  | **Endometriosis serum arginase-1** (RepOD PY1P9X), rank 2 | fit-free AUROC **0.846 (0.762-0.918)**, permutation p ≤ 0.0001, 120 vs 32 surgical controls | rule met |
  | **Long COVID cortisol** (MY-LC, Klein 2023), rank 3 | Delta AUROC over draw time + demographics **+0.253 (0.177-0.328)**; cortisol-only CV AUROC **0.963 (0.933-0.987)**, permutation p ≤ 0.001 | rule met |

  More on each:
  * **hEDS/HSD:** a modest group separation (30% sensitivity at 90% specificity). The site negative control was
    null.
  * **Endometriosis ARG1:** the result held within ELISA plates (0.854) and reproduces the paper's 0.848.
  * **Long COVID cortisol:** this reproduces the published 0.96 in the same released rows. The effect size (Hedges'
    g -2.8) is itself a caveat; see the results file.

  All three are **re-analyses of the authors' own released rows**, not independent replications, and all three use
  **healthy or surgical controls**, not the clinically relevant look-alike conditions.
* **One Long COVID run was invalid and was redone under a documented deviation.** The released MY-LC analytes are
  batch-integrated z-scores, so the plan's log10 transform silently dropped most Long COVID cases (7 vs 68 people).
  The invalid outputs are kept in `results/tables/_invalid_runs/`, and the deviation is in the Klein plan.

## What was done, and what was not

1. **Search.** Four parallel searchers (MCAS/POTS/EDS; Long COVID; ME/CFS; fibromyalgia/Lyme/migraine/IBS/
   endometriosis/gastroparesis) proposed 182 candidates. They used web search, the paperclip corpus, and repository
   APIs: Zenodo, figshare, Dryad, Mendeley, Dataverse, OSF, Metabolomics Workbench, MetaboLights, PRIDE, ImmPort,
   GEO, DataCite, Europe PMC and Pennsieve.

   mapMECFS was catalogued from its anonymous CKAN metadata API (119 packages). Every data resource answers HTTP 403
   with a login page; the clinical master labs file was re-probed here, 2026-09-24.
2. **Verification.** Every open claim was verified by retrieving a file or a listing anonymously through
   `measure_it.download` / `measure_it.http` into `data/raw/lab_candidates/<id>/`. There were no download errors.

   The downloads were small:
   * one per-person file, or the journal supplement bundle (Europe PMC `supplementaryFiles`);
   * MW REST datatables or factor tables;
   * MetaboLights / PRIDE sample sheets;
   * Zenodo, Dryad and PRIDE listings for large deposits.

   Gated routes were probed for their HTTP status only: mapMECFS 403, Dryad 401, RSU Dataverse.
3. **Inspection at schema level only.** The module counted:
   * rows;
   * label values;
   * people with the analyte present;
   * whether analytes share a person id.

   It never computed a biomarker-by-label contrast. The three analysis plans were written (21:13 UTC) before any
   group statistic was computed, and the Cinquina amendment at 21:14. The first outputs were written at 21:16
   (Cinquina ingest) and 21:21 (ARG1 results).

   What was seen before the plans:
   * label and covariate cross-tabs (for example, ELISA plate by group, and the Olink plate by group added as
     amendment 1);
   * published claims from the papers;
   * the names, but not the values, of the paper's top Olink proteins, visible when the supplement sheets were
     listed.
4. **Not done.** No account was created, no data-use agreement was accepted, and nothing behind a login was
   requested.

   Two records that appear to contain **direct personal identifiers** were **deliberately not downloaded**:
   * figshare 30127723 (Seibert 2026, PLOS ONE, post-COVID GPCR autoantibodies). The searcher saw column headers for
     surname, first name, date of birth, phone number and GP. It deleted its copies after seeing only the column
     headers and non-null counts.
   * Mendeley vzydrgtymt (Indonesian endometriosis vitamin D/zinc), which contains patient names and record numbers.

   Both are listed as `excluded_privacy`. Someone may want to notify the publishers.

**Status counts** (`access_status`): verified_open 109, unreachable 27, controlled 25, verified_registration_only 10,
duplicate 6, excluded_privacy 3, not_as_described 2.
* "unreachable" includes aggregate-only publications (recorded negatives) and Dryad files that refuse anonymous
  scripts.
* Of the 109 verified_open, 7 carry no lab analyte (`lab_values_in_file = no`), 4 have circular labels, and 9 have
  undocumented group codes.

## Ranked verified-open lab datasets (case-control capable)

**Eligibility.** All of the following:
* verified_open;
* label mapping stated in the file, or verified here;
* the label is not defined by the lab itself;
* the file carries lab values;
* at least 20 cases and at least 20 controls with the analyte, **counted here from the downloaded files**;
* not a duplicate.

**Score.** The rubric is in the module and was written before any analysis. The seven components are:
* label quality, 0-3;
* n, 0-3 (1.5·log10 of the smaller group);
* lab relevance, 0-3 (3 = a biomarker the brief names for that condition, such as cortisol for Long COVID);
* completeness, 0-1;
* controls, 0-1;
* multi-condition / disease comparators, 0-1;
* person linkage across analytes, 0-1.

| rank | dataset | condition | biomarkers | label | cases / controls with lab | score |
|---|---|---|---|---|---|---|
| 1 | Cinquina 2026, Europe PMC PMC13081554 | hEDS + HSD | 460 Olink serum proteins | 2017 hEDS criteria (per-patient checklists) | 176 / 176 | 12.00 |
| 2 | RepOD PY1P9X (Pliszkiewicz 2024) | endometriosis | serum ARG1, ARG2, arginase activity, glucose, creatinine, AST, ALT | laparoscopy + rASRM stage; surgical and healthy controls | 120 / 85 | 11.35 |
| 3 | Klein 2023 Nature, Suppl. Table 3 (MY-LC) | Long COVID | cortisol + 143 plasma mediators, flow, antibody reactivity (one row per person) | LC clinic cohort (symptoms > 6 wk after WHO-defined COVID); HC + CC | 98 / 79 | 11.34 |
| 4 | Serum MIF, BMC Womens Health 2020 (figshare) | endometriosis | serum MIF | surgical, rASRM | 150 / 150 | 11.00 |
| 5 | Baraniuk 2021 CSF metabolomics (figshare) | ME/CFS, Gulf War illness | ~187 targeted CSF metabolites | Fukuda + CCC | 45 / 21 (+20 GWI) | 10.98 |
| 6 | Che 2022, Metabolomics Workbench ST002000-3 | ME/CFS | 4 plasma platforms | CCC 2003 | 106 / 91 | 10.94 |
| 7 | Serum zinc/Fe/Mg/Cu, Zenodo 17773060 | migraine | Zn, Fe, Mg, Cu | ICHD-3, matched 1:1 | 53 / 53 | 10.59 |
| 8 | Nagy-Szakal 2018, MW ST000800-series | ME/CFS (± IBS) | targeted oxylipin panel | Fukuda/CCC | 50 / 50 | 10.55 |
| 9 | Naviaux 2016, MW ST000450 | ME/CFS | targeted plasma metabolomics | Fukuda / CCC / IOM | 45 / 39 | 10.39 |
| 10 | Hoel 2026, PRIDE PAD000026 | ME/CFS | SomaScan 7k | CCC | 50 / 29 | 10.19 |
| 11 | Mortensen 2017, ECM turnover, PLOS ONE 0185855 | IBS (+ Crohn's and UC comparators) | BGM, EL-NE, C5M, Pro-C5, CRP | clinical (Rome criteria not stated) | 22 / 24 | 10.01 |
| 12 | Talla 2023, Nat Commun | Long COVID | Olink Explore 1,463 | study-defined PASC (persistent symptoms; threshold not verified); recovered + uninfected | 55 / 46 | 9.99 |
| 13-31 | Caboni FM lysoPC; Lyme ST001223; BioMapAI JAX ME/CFS (clinical labs, transformed); endometriosis xenoestrogen; Germain 2022; MW ST003103 LC; Mexico LC metabolome; 2 endometriosis MW studies; ST000888 (FM arm); Cervia-Hasler SomaScan; Armstrong NMR; PTLDS ST001391; IBS MTBLS2774; 5 more LC sets | | | | see CSV | 7.95-9.98 |

**The rank-2/3 tie.** ARG1 (11.347) and MY-LC (11.340) differ by 0.007.
* The robustness table favours ARG1 for #2 when one of lab relevance, completeness, controls or linkage is dropped.
* It favours MY-LC when label quality is dropped.
* It favours the MIF or Baraniuk sets when the multi-condition or n component is dropped.

Rather than break the tie by a post hoc rule, **both were analysed**.

The brief prefers "LTE4/tryptase for MCAS or multiple umbrella conditions", and no MCAS dataset qualifies. Of the
top three, only Cinquina carries multiple umbrella labels: hEDS, HSD, and per-patient flags for dysautonomia/POTS,
chronic fatigue, MCAS/allergy and functional GI.

**Download URLs (top 3):**

1. Cinquina: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13081554/supplementaryFiles`, members
   `12014_2026_9588_MOESM4_ESM.xlsx` (NPX) and `..._MOESM6_ESM.xls` (clinical).
2. ARG1: `https://repod.icm.edu.pl/api/access/datafile/48824`.
3. MY-LC: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06651-y/MediaObjects/41586_2023_6651_MOESM4_ESM.xlsx`.
   It needs a browser User-Agent.

## Pre-specified analyses and results

The plans were locked before computation: `docs/ANALYSIS_PLAN_HEDS_HSD_OLINK_SERUM_CINQUINA2026.md`,
`docs/ANALYSIS_PLAN_ENDO_ARG1_REPOD.md` and `docs/ANALYSIS_PLAN_KLEIN2023_MYLC_ML_TABLE.md`. The results are in
`results/<ID>_RESULTS.md`. Published claims are kept in a separate column in every results file.

| dataset | primary (locked) | our result | published claim (separate) | key secondaries |
|---|---|---|---|---|
| hEDS/HSD Olink | Italian hEDS+HSD (88) vs Italian controls (176); 458 assays, L2 logistic, CV AUROC; rule: CI low > 0.5 and permutation p < 0.05 | **0.684 (0.624-0.747)**, perm. p ≤ 0.001; sensitivity at 90% specificity 0.30 | 69 DEPs pooled; LASSO/RF fitted on all data; no held-out estimate | site negative control, US vs Italian patients: 0.514 (0.441-0.586), p = 0.40, null; hEDS vs HSD 0.480 (null, as published); 24/458 proteins q < 0.05, all in the paper's list (MYOC, CA14, NUDT5, MAX, FOXO1, COMP ...); probands only 0.712 |
| Endometriosis ARG1 | endometriosis (120) vs surgical non-endometriosis controls (32), fit-free AUROC, higher = case | **0.846 (0.762-0.918)**, perm. p ≤ 0.0001; sensitivity at 90% specificity 0.66 | 0.848 (0.769-0.926), 105 vs 22 | plate-stratified 0.854 (0.775-0.928); vs healthy 0.918, batch unknown (no plate recorded for any healthy control); ARG2, activity and routine labs q > 0.2; no correlation with stage (rho ≈ 0) |
| MY-LC cortisol | Delta AUROC (draw time + age + sex + BMI + cortisol) minus covariates; rule: Delta CI > 0 and cortisol-only perm. p < 0.05 | **+0.253 (0.177-0.328)**; cortisol-only 0.963 (0.933-0.987), perm. p ≤ 0.001 | cortisol-only AUC 0.96 (0.93-0.99) in a matched subset | LC vs HC 0.973, vs CC 0.958; 144-analyte panel 0.947; 61/144 analytes q < 0.05; post hoc within draw-time tertiles 0.969; Hedges' g for cortisol -2.8 on z-scores; no batch column released |

## What the verification changed

| candidate | searcher claim | verified |
|---|---|---|
| MY-LC (Klein 2023) | cortisol etc. as concentrations (column names say ng/mL, pg/mL) | the values are **batch-integrated z-scores** (mean 0, SD 1; all 144 analytes have negatives). This caused Deviation 1; any concentration-scale re-use is wrong |
| MY-LC cohort codes | 3 = LC; 1/2 inferred | 1 = HC, 2 = CC, verified from the infection-test columns; `x0_Censor_Complete == 0` reproduces the paper's 40/39/99 |
| Cinquina P/H prefixes | P vs H unconfirmed | P = hEDS and H = HSD, confirmed against the MOESM6 sheets. **All 176 controls are Italian; half the patients are American** home draws (paper text). The primary was restricted to Italian patients; the site control came out null |
| Cinquina plates | not described | 4 Olink plates, every group on every plate (plan amendment 1) |
| RepOD ARG1 | 127 / 88 controls, types unclear | Ctrl = 25 benign-surgery controls, Ctrl2 = 9 myomas, Ctrl3 = 51 healthy + 3 other (`diagnosis` column); 120/85 with ARG1. **No ELISA plate for any healthy control**, so patient-vs-healthy cannot be checked for batch |
| Appelman 2024 | cortisol LC vs healthy 25/21 | baseline cortisol is present for 25/18. Below the ranking floor; now analysed for metabolomics by another module (`measure_it.labs.appelman_metabolomics`) |
| Kedor 2022 | labs for PCS vs PCS/CFS | the `Diagnose` label is written once per block (forward-filled here): 19 PCS/CFS vs 23 PCS, the only open symptomatic-control Long COVID / ME/CFS lab set. Below the floor |
| Lyme ST001223 | 329 vs 189 | those are samples. Subjects: 166 Lyme vs 95 healthy (duplicate runs per subject) |
| ST000888 fibromyalgia arm | FM comparator | 31 FM vs 89 healthy subjects, but the 95 features were **pre-selected as a Lyme signature**, so this is not an FM test |
| BioMapAI JAX ME/CFS | 48 clinical labs | values are MaAsLin-residualised / 0-1 scaled, **not clinical units**; the repo has no licence file. Ranked 15 |
| RepOD pediatric PTLDS MMP | `.tab` file | served as xlsx (fixed); control rows have a blank batch field (possible group-batch confounding) |
| MTBLS12737 (IBS urine), MTBLS13630 (migraine PFO) | metabolite data | **not_as_described**: MAFs have no abundance values |
| Yin 2024 LIINC Olink/CyTOF (Dryad, CC0) | open | **unreachable to scripts**: Dryad file API returns 401 anonymously (probed). A browser download may work |

## Registration-only, worth unlocking (one free account, no DUA reported)

* **mapMECFS** (RTI CKAN). Package metadata is anonymous; every file answers 403 with a login page.
  * The highest-value packages:
    * NIH PI-ME/CFS clinical master labs;
    * NIH blood and CSF flow cytometry;
    * NK chromium-release;
    * the Cornell CRC packages: plasma/EV cytokines, SomaLogic 7k, urine and plasma metabolomics, and anti-pathogen
      and autoantibody panels. `cor_id` links these to the wrist ActiGraph and 2-day CPET packages in
      `docs/DEVICE_DATASET_DISCOVERY.md`;
    * Columbia/Che 2025 multi-omics with exercise (biogenic amines, oxylipins, complex lipids, SomaLogic,
      cytokines).
  * **One account would unlock person-linked labs + device data in ME/CFS**, which no open source provides.
* **ImmPort** (free account + data-use terms): the IMPACC post-acute studies (Olink, metabolomics, cytokines) and
  LC studies with clinical labs; Lyme SDY1395 (mass spectrometry, CyTOF). SDY2270 (a mast-cell-disorder arm)
  reports `has_lab_test = N`, so it is not worth unlocking for labs.
* **Kaggle** Germain 2020 metabolomics mirror (the same data are open in the Europe PMC supplement, rank not
  computed).
* **Synapse**: the Mount Sinai replication cohort for Cervia-Hasler (SomaScan + clinical labs).

## Controlled / on request (recorded, not requested)

* **MCAS / mast cell:**
  * Vivli midostaurin trial: tryptase, KIT D816V burden.
  * All of Us EHR: LOINC tryptase, catecholamines and N-methylhistamine, if ordered.
  * PARIS tryptase CNV: on request.
* **POTS:**
  * Hall 2022 GPCR autoantibodies, 116 vs 81 (on request).
  * Kharraziha 2020 GPCR functional assay.
  * Malmö cohorts.
  * Tufvesson 2026 gut hormones.
  * Yao 2025, which has the ideal POTS + HSD + MCAS label structure (on request).
* **Long COVID:** RECOVER-Adult (25 routine labs, about 10k adults), PHOSP-COVID (Olink), Phetsouphanh 2022,
  Schultheiß 2022, Sotzny 2022 (GPCR autoantibodies), Greene 2024.
* **ME/CFS:** CDC MCAM (NK and clinical labs), UK ME/CFS Biobank (with an MS arm), the Chronic Fatigue Initiative
  (Hornig 51-plex cytokines), DecodeME (genotypes only), and the RSU Latvia viral/anti-cytokine set (restricted
  files).
* **Other:** NIDDK gastroparesis registry, BBMRI-NL migraine metabolome, restricted Zenodo migraine GFAP/NfL, and RSU
  fibromyalgia gut-barrier markers.

## Limitations

* **Recruitment design.** Every ranked dataset compares diagnosed people with healthy volunteers, or with surgical
  controls for ARG1. None tests a lab against the clinically relevant look-alike conditions. The only open lab sets
  with symptomatic comparators are all below the n floor or indirect:
  * Kedor: post-COVID with vs without ME/CFS, 19 vs 23;
  * the IBS ECM set: IBS vs Crohn's and UC;
  * Baraniuk: ME/CFS vs Gulf War illness;
  * ST000888: fibromyalgia inside a Lyme look-alike design.
* **The analyses re-use the authors' own released rows.** Agreement with the paper (MY-LC 0.96, ARG1 0.848) is
  **reproduction, not replication**.
* **Label-quality scores come from the papers' stated criteria as reported by the searchers.** They were re-read in
  full here only for the three analysed datasets.
* **The rubric is a convenience, not an evidence grade.** The rank-2/3 tie shows how sensitive the ordering is.
* **Many open ME/CFS and Long COVID omics sets are small** (17-60 per group) and have not been analysed. The ranking
  says where to look next, not what they show.

## Open items (for whoever owns shared files)

* Resolved 2026-09-24: `xlrd` is now a project dependency (`pyproject.toml`, `uv.lock`), so the Cinquina clinical
  sheet, the endometriosis MIF file and the `.xls` discovery inspections run under plain `uv run`.
* Resolved 2026-09-24: pipeline steps `labs_dataset_discovery` (offline: `download=False`, inspects the files on disk),
  `labs_heds_hsd_olink_cinquina2026`, `labs_klein2023_mylc` and `labs_endo_arg1_repod` run the modules; their
  `participants__<id>` / `participant_labs__<id>` partitions enter the `digital_person` unions (datasets kept
  separate, namespaced ids) and their `phenotype_signatures__<id>` partitions the canonical `phenotype_signatures`.
  The registry entries are merged into `SOURCE_REGISTRY.yaml` by `registry_merge`.
* The measurement classes used are `proteomics` (Cinquina) and `blood_biomarkers` (ARG1, MY-LC; MY-LC also lists
  `immune_assays`). Condition ids used: `eds_hsd`, `endometriosis`, `long_covid`.
