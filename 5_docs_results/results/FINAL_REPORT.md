# Measure It to Cure It: public-data alpha, final report

**Date:** 2026-09-24. **Build state:** last commit `209f666` (linked omics, GEO cohorts, published device evidence,
device dataset discovery and the two device analyses) plus an uncommitted working tree with the lab-dataset work
(labs vs diagnosis, lab dataset discovery and its three analyses, published lab evidence, NHANES lab layers, NHANES
2003-2006 replication, data access plan) and their pipeline integration (190 changed or untracked paths when this
header was written). Until that is committed, a clone does not match this report (§23).
**Brief:** `docs/SPEC.md`. **Conventions and guardrails:** `docs/CONVENTIONS.md`.

> **Update 2026-10-07.** Sections below describe the build of 2026-09-24 and are not
> rewritten. Since then: (1) a shared person-level vocabulary and device-subgroup phenotypes (`measure-it byod
> subgroup`), similar-people estimates and query packs (`measure-it byod similar`), BRFSS small-area Long COVID burden
> and measurement-specific clinician activity were added (docs/HARMONIZATION_CONTRACT.md, docs/MAESTRO_RUNBOOK.md;
> 40 sources); (2) the primary ranking now uses the BRFSS 2023 small-area estimate as Long COVID burden instead of the
> inherited Household Pulse state value, and Medicare billing activity as clinic capacity for measurements with a
> dedicated billing code (wearable bundle, autonomic testing, CPET; nailfold capillaroscopy has no code and keeps the
> density basis). Current numbers are in results/SCORING_RESULTS.md and results/STAGE_LAYERS_SCORING.md: for the
> primary query the old rules (S9) keep 2 of the new top-10; 4 of the top-10 are now in the top-10 in at least half
> of the Monte Carlo draws (0 before) and the median top-25 rank-interval width is 195 ranks (319 before); the
> clinic-location shuffle now moves clinic capacity (Spearman 0.487 vs 0.965). Rows 6 and 7 of "does not
> demonstrate" (§1) are therefore partly superseded; the within-state ordering of the Long COVID estimates is still
> modelled, not observed (results/BRFSS_SAE_RESULTS.md).

Numbers in this report are copied from files in this repository, and the file is named next to each (paths are
repository-relative; a bare `*.csv` name means `results/tables/<name>`). The 3,208-claim number audit
(`number_audit.csv`) covered the other documents, **not this report** (0 of its rows are for this file); the report's
numbers were checked afterwards by review, and those checks and corrections are folded in here. Where a component is
partial, not done, or rests on a proxy, the report says so in the same sentence.

---

## 1. Executive summary

**What was built.** A reproducible, public-data-only alpha of a modular *measurement-validation and deployment
engine* for invisible illnesses. It runs the SPEC chain (person-level phenotype -> objective wearable phenotype ->
ontology -> condition-level molecular context -> candidate measurement technology -> geographic burden and unmet need ->
clinic and research-site readiness -> candidate deployment opportunity) end to end for **Long COVID / ME/CFS x wearable
autonomic/activity monitoring**, over 37 registered public sources (`SOURCE_REGISTRY.yaml`, merged from
`data/raw/*/registry_entry.yaml`: 36 ingested, of which CDC PLACES at county and ZCTA level only; 1 partial) held
in four data layers that the code keeps apart. Beyond wearables, it re-analyses public person-level labs, omics and
other devices under plans locked before analysis, and curates published device and lab evidence (§14). It answers "what should we measure, where, and who could realistically
deploy it?" through 18 MCP tools, a FastAPI app, a five-page Streamlit dashboard and an LLM-free agent, and every
recommendation traces back to raw files with their sha256 (`docs/MCP.md`). No private CAPRIO/MAESTRO data are used;
capillaroscopy exists only as a stub adapter.

**What the alpha demonstrates**

| # | demonstrated | key number (file) |
|---|---|---|
| 1 | The whole chain rebuilds offline and deterministically | the 49-step DAG re-run offline on the final code (run `20260924T150209Z`): against the outputs it replaced, 279 tables identical, 4 changed only by this report's fixes, 3 in runtime fields only, and the canonical `phenotype_signatures` changed because the GEO-cohort and linked-omics signature partitions (§14) were added during the run. Earlier, two full runs of the then 47-step DAG (before the 2 figure steps were added) reproduced all 280 compared tables (`docs/REPRODUCIBILITY.md` §5.1, §5.7). Offline only: a cold build from an empty checkout has not been run. The §14 modules became 20 new steps of the now 75-step DAG afterwards; their offline run `20260924T215649Z` reproduced every output they had written as standalone modules (422 tables identical, 3 numerically equal, the 12 others the intended union, Data.gov and recommendation-text changes; §21) |
| 2 | Wrist accelerometry carries person-level information about function and symptoms | added to demographics it raises AUROC by +0.019 to +0.032 for 4 of 5 symptom/function targets (plus the mortality positive control), all 95% CIs > 0 (`nhanes_increments.csv`); the one robust result over clinical data, functional limitation, is partly definitional |
| 3 | The engine finds which objective measurements are registered for a condition, automatically, with an LLM-reviewer precision spot-check (not a domain-expert audit; recall unmeasured) | 10,646 text-mined mentions over 26 measurement classes; Long COVID: 85 trials use the wearable bundle objectively, 59 without its device-agnostic HRV class (`results/MEASUREMENT_DISCOVERY.md`) |
| 4 | Facility and research information changes where the engine would deploy (Test 5), partly by construction | adding desert, clinic capacity and research readiness to burden + vulnerability keeps 3 of the top-10 and 13 of the top-25 (Jaccard 0.35; contrast T5-2); the literal SPEC comparison, burden-only vs full, keeps 0 of 10 (tie-averaged 0.29: 7 counties tie at the burden-only top-10 boundary) and 4 of 25. Shuffling every facility location nationally still leaves composite Spearman 0.94 with the observed ranking (`test5_contrasts.csv`, `test6_scoring_clinic_shuffle.csv`) |
| 5 | Every recommendation is traceable to source rows and raw files | `trace_evidence` resolves all 104 object ids of a compact top-10 ranking (`docs/MCP.md` §5) |
| 6 | A new measurement modality runs through the same downstream code | the capillaroscopy stub bundle runs through the identical scoring path unchanged (S5; `test7_sensitivity_comparisons.csv`); what the swap changes is row 9 below |
| 7 | In public person-level data, a few objective measurements beyond wrist accelerometry separate **diagnosed** patients from **healthy or surgical** controls under plans locked before analysis (re-analyses of the authors' own released rows, not independent replication; spectrum and, for steps, incorporation bias; §14) | MUSCLE-ME hip-accelerometer daily steps, Long COVID + ME/CFS vs healthy: AUROC 0.831 (0.725-0.918), permutation p 0.0001, but 0.612 (0.459-0.763) if the 8 controls without steps had the lowest value, and no better than CPET VO2peak in the same people (difference -0.043, -0.192 to 0.092) (`results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md`); MY-LC serum cortisol, Long COVID vs controls: 0.963 (0.933-0.987), a finding two other published cohorts did not reproduce (PLE-042, PLE-043; `results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md`); endometriosis serum ARG1 vs surgical controls: 0.846 (0.762-0.918) (`results/ENDO_ARG1_REPOD_RESULTS.md`); hEDS/HSD serum Olink, Italian patients vs Italian controls: 0.684 (0.624-0.747), 30% sensitivity at 90% specificity (`results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md`) |
| 8 | The NHANES wearable conclusions largely carry over to another device and cohort | NHANES 2003-2006 hip ActiGraph counts (6,283 adults, never pooled with 2011-2014): 8 of 12 pre-specified Test 1 / Test 2 verdicts replicate (Test 1 4/6, Test 2 4/6); functional limitation keeps a small increment over clinical data (+0.0041, +0.0005 to +0.0075) that is fragile (CI > 0 in 2 of 3 model families; the 2003-04 -> 2005-06 temporal split CI includes 0) (`results/NHANES_2003_2006_RESULTS.md`) |
| 9 | The lab pipeline sees lab signal where it is known to exist, and published evidence is stored with verifiable quotes | NHANES positive-control gate passed (labs add +0.147, +0.136 to +0.159, for diabetes without glycaemic markers; +0.119 weak/failing kidneys; +0.048 gout; `results/LABS_VS_DIAGNOSIS_RESULTS.md`); 66 published device claims and 109 published lab claims, every quote re-fetched from Europe PMC and verified verbatim (`docs/PUBLISHED_DEVICE_EVIDENCE.md`, `docs/PUBLISHED_LAB_EVIDENCE.md`) |

**What the alpha does not demonstrate**

| # | not demonstrated | why (file) |
|---|---|---|
| 1 | An objective **autonomic** wearable phenotype for any invisible-illness label | NHANES 2011-2014 has no heart rate; the only Long COVID label with wearable heart rate among the public datasets processed here (Uwakwe 2025, 31 cases) gives a **null** wearable AUROC 0.504 (0.404-0.605) and increment -0.024 (-0.076 to 0.027); single-day acute RHR elevation is as common in placebo windows (62% vs 56%) (`results/WEARABLE_STANFORD_RESULTS.md`) |
| 2 | A validated ME/CFS phenotype | the NHANES group is a constructed **ME/CFS-like proxy** (132 cases, prevalence 0.017): 14 of 26 features differ after FDR, but the predictive increment is null (+0.0139, -0.0126 to +0.0387), wearable-only AUPRC is 0.024 (0.019-0.033), close to the base rate, and the difference is partly definitional (`results/WEARABLE_NHANES_RESULTS.md`, `nhanes_model_performance.csv`) |
| 3 | A useful wearable increment over clinical data | over a 50-variable clinical model the increment is robust only for functional limitation (+0.0069, +0.0040 to +0.0097) |
| 4 | Independent molecular support for autonomic measurement in the demo conditions | condition-level enrichment only; the ME/CFS and POTS autonomic/neuronal/mitochondrial system links are carried by drug-target gene families of trialled drugs (post-hoc flagged), and after removing them none of Long COVID, ME/CFS or POTS has a non-literature-supported system (`results/MOLECULAR_COHERENCE.md`). Until this report the `ask` agent printed these links as "supported physiological systems" without the flag, an overstatement in user-facing output; it now prints the support tiers (§16, §23 item 9) |
| 5 | Direct county-level burden for the demo conditions | Long COVID counties carry the **inherited state** HPS value; ME/CFS is a level-C claims proxy; POTS and dysautonomia are level D (8 of 14 conditions have no usable burden) (`results/GEOGRAPHY_RESULTS.md`) |
| 6 | A stable short list of regions | 0 of the top-10 counties stay in the top-10 in at least half of 1,000 Monte Carlo draws; median 5th-95th rank-interval width of the top-25 is 319 ranks (`results/SCORING_RESULTS.md`) |
| 7 | Measurement-specific clinic capability | clinic_capacity barely moves when facility locations are shuffled (Spearman 0.965): it measures facility and provider density. At composite level the shuffled ranking keeps Spearman 0.94 (0.936-0.945) with the observed one and 8.2 of the top-25, about the same retention as adding the facility components at all (T5-3: 8 of 25) (`test6_scoring_clinic_shuffle.csv`, `test5_contrasts.csv`) |
| 8 | Wearable x omics linkage within one person | mapMECFS per-participant physiology files are login-gated (29 of 29 anonymous probes HTTP 403); 0 wearable x omics pairs have comparable identifiers (`docs/MAP_MECFS_LINKAGE_AUDIT.md`) |
| 9 | Measurement-specific deployment targeting | swapping the wearable bundle for nailfold capillaroscopy (4 trials, 0 NIH core projects, no FDA product code) leaves the county ranking almost unchanged: composite Spearman 0.994, top-10 overlap 9 of 10, clinic_capacity Spearman 0.994, research_readiness 0.953. Burden, vulnerability and desert do not depend on the measurement, so under the equal weights the engine answers "given technology Y, where?" nearly the same way for any Y (`test7_sensitivity_comparisons.csv`). **Partly addressed 2026-09-25 (metric -> translation link, `docs/ANALYSIS_PLAN_METRIC_LINK.md`, `results/SCORING_RESULTS.md` §12):** each condition x measurement now carries a performance record from this project's own results (published claims flagged as tier 3); under the new `evidence_weighted` weights capillaroscopy and clinic autonomic testing have no performance record for Long COVID or ME/CFS and are not ranked (UNKNOWN, never scored 0 or 1), and the joint region x measurement top-25 becomes CPET 18 / wearable 7 (equal: 7/7/7/4). Within one measurement evidence still cannot reorder regions, and every record compares diagnosed patients with healthy controls |
| 10 | Person-level omics signal for ME/CFS or Long COVID that passes a pre-specified test or transfers between cohorts | mapMECFS person-linked omics: 0 of 5 layers pass (bar: one-sided permutation p < 0.01, Bonferroni over 5 layers; best CSF metabolomics AUROC 0.825, p 0.0130, which falls to 0.705 after in-fold age removal; `results/LINKED_OMICS_RESULTS.md`); GEO cohorts: 0 of 4 primary cross-cohort transfers replicate; the one near-perfect within-cohort classifier (GSE270045, 0.998) has cases and controls from different subject-id sources (post-hoc p 1.5e-08); the ME/CFS cfRNA classifier adds nothing to its covariates (-0.026, -0.077 to 0.095) (`geo_cohort_transfer.csv`, `geo_cohort_confounds.csv`, `geo_cohort_classifier.csv`) |
| 11 | Routine labs that separate the umbrella invisible-illness labels available in public data | NHANES 2011-2014: labs add nothing detectable to demographics for prescription-coded fibromyalgia (+0.022, -0.027 to +0.072), migraine, insomnia or myalgia, nor between them; the ME/CFS-like proxy's lab increment (+0.046) is smoking plus the definition's exclusion diagnoses (+0.005, -0.020 to +0.032, against controls without them); across 58 lab layers, 1 of 47 lab-block rows is supported (+0.0025) and 0 of 292 per-layer screens (`results/LABS_VS_DIAGNOSIS_RESULTS.md`, `results/NHANES_LAYERS_RESULTS.md`); Appelman 2024 metabolomics: plasma and muscle primary tests null (`results/LABS_VS_DIAGNOSIS_RESULTS.md` Part 2) |
| 12 | Any measurement tested against clinical look-alikes, and any open MCAS or POTS lab data | every dataset analysed in §14 uses healthy, surgical or general-population controls; no open person-level MCAS or POTS lab dataset with controls exists (`docs/LAB_DATASET_DISCOVERY.md`); fibromyalgia thermography fails its primary (+0.053, -0.028 to 0.129; `results/FM_THERMOGRAPHY_RESULTS.md`) |

**Bottom line.** The alpha shows that the *infrastructure* half of the thesis works on U.S. public data: fragmented
person-level, molecular, geographic and facility data can be normalised, kept in their own layers, combined into
transparent candidate deployment opportunities, stress-tested (Tests 1-7) and traced back to source. The *scientific*
half (that objective measurement reveals latent invisible-illness phenotypes) is **not** shown by the public data
processed here: the only Long COVID label with wearable heart rate is null, the NHANES ME/CFS evidence is a proxy, and
the MUSCLE-ME step counts (§14) do not survive the locked missing-control worst case. Nor does the
composite ranking yet depend much on the measurement or on where facilities really are (rows 7 and 9). Looking beyond
wearables (§14) adds candidates but not validation: hip-worn step counts, serum cortisol, serum ARG1 and a serum
proteome separate diagnosed patients from healthy or surgical controls in the authors' own released rows, while
person-linked and GEO omics, routine NHANES labs against umbrella labels, a small metabolomics set and thermography are
null or confounded, and no measurement has been tested against the clinical look-alikes it would have to tell apart. The engine reports this rather than
hiding it. Every ranked region is therefore a candidate deployment opportunity for **pilot evaluation**:
deployment there would generate the objective evidence that public data lack, not apply a validated diagnostic pathway.

---

## 2. Status at a glance

Status words: **done**, **partial** (the missing part is named), **not done (optional)**, **null / negative result**
(the component works; the evidence it produced is null or negative). Starting point: `docs/SPEC_COMPLIANCE.md`, a
**pre-report snapshot** of 2026-09-24 (it predates this report and the README, counts 1,177 tests and 38 scanned text
files); the rows below are updated for the final state and for the review of this report.

| component | status | note |
|---|---|---|
| 37 source ingestions with DATA_AUDIT.md and registry fragment (27 original + 10 added for §14) | done (36 ingested, **1 partial**) | mapMECFS partial: per-participant physiology gated behind a login; `SOURCE_REGISTRY.yaml` re-merged from the fragments: 37 sources |
| CDC PLACES | partial | county and ZCTA levels only; tract and place levels not ingested (nothing downstream uses them) |
| Optional public wearable benchmarks (PhysioNet, CovIdentify, Apple Watch sleep) | not done (optional) | the Stanford acute cohorts with placebo windows filled the signal-check role |
| NHANES raw 80-Hz mode; temporal/deep model | not done (optional / conditional) | SPEC allows both to be skipped |
| NHANES heart rate / HRV | not available in the data | never derived or imputed; the exam pulse is a clinical variable |
| Long COVID wearable label (Uwakwe 2025) | done; **null result** | AUROC 0.504 |
| Core data model (21 SPEC tables as DuckDB `spec.*` views; 11 required datasets) | done | `validate` (2026-09-24T22:27Z): PASS 6/6, 11/11 required datasets, 157/157 processed tables with provenance |
| Phases 1-6 | done | Phase 1 evidence is null (Long COVID) or proxy (ME/CFS) |
| MCP server (15 SPEC tools + 3), agent, API | done | minor signature supersets (§23 item 7); rankings are national only (no state filter) |
| Dashboard (5 SPEC pages) | done, **page 1 partial** | page 1 shows one area layer at a time (burden, diagnostic desert or opportunity), not independently toggled layers, and its facility overlay shows HRSA sites only (not NPPES or `clinic_registry`); the point overlays (HRSA, trial sites, NIH, specialist density) toggle independently |
| Adapter interface, WearableAdapter, CapillaroscopyAdapter stub | done | stub raises NotImplementedError by design |
| Validation Tests 1-7 | done | outcomes in §13 (several null or negative) |
| Figures 1-10, captions, dashboard screenshots | done | `results/figures/` |
| Primary demo story (SPEC, 9 steps) | **3 of 9 partial** | table below |
| `deployment_opportunities.parquet` member_components JSON | done | the 12,876 pre-fix `"nan"` strings for level-D set members were cleared by the final run `20260924T150209Z` (0 of 76,680 rows now hold one) |
| Person-linked omics, mapMECFS (5 layers, 21-42 people each) | done; **null result** | 0 of 5 layers pass the pre-specified test; no wearable is linkable (§14.4) |
| GEO case/control omics cohorts (7 subsets, 88 series screened) | done; **null / negative result**; write-up **partial** | 0 of 4 primary transfers replicate; GSE270045 source-confounded; cfRNA explained by covariates. No `results/GEO_COHORT_RESULTS.md` exists: numbers come from `results/tables/geo_cohort_*.csv` (§14.4) |
| Device dataset discovery (90 candidates) + 2 pre-specified device analyses | done | MUSCLE-ME steps: locked rule met on complete cases, **null** in the missing-control worst case; fibromyalgia thermography: **null** on its primary (§14.3-14.4) |
| Published device evidence (66 claims) and published lab evidence (109 claims) | done (targeted curation, not a systematic review) | every quote verified against Europe PMC; nothing reproduced by this project (§14.7) |
| Labs vs diagnosis (NHANES 2011-2014 + Appelman 2024 metabolomics) | done; **null result** for umbrella labels | positive-control gate passed (§14.4) |
| Lab dataset discovery (182 candidates) + 3 pre-specified lab analyses | done | 3 of 3 locked rules met (MY-LC cortisol, endometriosis ARG1, hEDS/HSD Olink); re-analyses of released rows against healthy or surgical controls (§14.3) |
| NHANES 2011-2014, 58 laboratory layers x accelerometry | done; **null result** for labs | 1 of 47 lab-block rows supported; accelerometry still adds for functional limitation over 6 of 7 lab bases (§14.4) |
| NHANES 2003-2006 hip-accelerometry replication | done | 8 of 12 pre-specified verdicts replicate (§7.5) |
| Data access plan (All of Us, RECOVER, UK Biobank, Snyder releases) | done (plan only) | counts of people with data or tests, read from public pages; nothing requested or analysed (§14.8) |
| Pipeline integration of the §14 modules | done | 20 steps added (74 in the DAG); offline re-run `20260924T215649Z` reproduced their outputs (§21). `measurement_phenotype_signal` / Phase 3 evidence still read only NHANES and Stanford (§23) |

**Primary demo story, step by step**

| step (SPEC) | status | note |
|---|---|---|
| 1 Public wearable data show an objectively measurable activity/autonomic phenotype | **partial** | activity differences for a constructed ME/CFS-like proxy (§7) and, in MUSCLE-ME, lower hip-accelerometer daily steps in diagnosed Long COVID / ME/CFS than in healthy controls, not robust to the 8 missing controls and part of the case definition (§14.3); no autonomic phenotype on any invisible-illness label |
| 2 Clinical and ontology data connect the phenotype to invisible illnesses | done | for the wearable datasets of §7 through proxy, self-report and prescription-code labels, not clinical diagnoses (§7, §19); MUSCLE-ME adds criteria-based Long COVID and diagnosed ME/CFS labels, normalised to the same ontology ids, against healthy controls only (§14.3) |
| 3 Molecular resources show pathways and candidate measurable biology | **partial** | condition-level enrichment exists, but no demo condition has non-literature, unflagged support for any system, so nothing supports autonomic measurement; the system -> measurement-class link is a curated assumption (§8) |
| 4 ClinicalTrials.gov, RePORTER and openFDA identify technologies under investigation or deployable | done | within the curated 26-class vocabulary; precision is an LLM-reviewer spot-check, recall unmeasured (§9) |
| 5 Federal population data locate burden and vulnerability | done | Long COVID county burden is an inherited state value, ME/CFS a claims proxy, POTS/dysautonomia level D (§10) |
| 6 NPPES, HRSA, ClinicalTrials.gov and RePORTER show where infrastructure is | done | locations are mostly ZCTA centroids; entity-resolution recall 0.528 (§11) |
| 7 The engine identifies measurement deserts and candidate deployment sites | **partial** | per-query absence notes and the diagnostic desert exist; the systematic nonmetro technology-experience measurement desert is **not supported** (§11) and the county desert partly measures small population (§10) |
| 8 An agent answers "what, where, who" | done | deterministic, LLM-free; national rankings only (§16) |
| 9 CAPRIO as a future device-specific implementation | done (interface level) | stub adapter and S5 swap; with today's data the ranking barely depends on the measurement (§1 row 9, §20) |

---

## 3. Exact datasets used

All 37 sources with a registry fragment (`data/raw/<source_id>/registry_entry.yaml`, from which
`SOURCE_REGISTRY.yaml` is merged: 37 sources): the 27 of the original build plus 10 added on 2026-09-24 for the
analyses of §14 (marked *new* below). Retrieval dates are the fragments' `retrieved_at` (UTC). Per-source audits:
`data/raw/<source_id>/DATA_AUDIT.md` (the 27 original audits contain all 17 SPEC fields, `docs/SPEC_COMPLIANCE.md` §7;
the 10 new audits were not part of that check). Discovery probes (`data/raw/device_candidates/`,
`data/raw/lab_candidates/`, `data/raw/snyder_candidates/`, `data/raw/controlled_cohorts/`) have no registry fragment:
they are catalogues of candidates, not ingested sources.

**Layer 1: person-level**

| source_id | source (publisher) | version / vintage | retrieved | unit of observation | status |
|---|---|---|---|---|---|
| `nhanes_2011_2014` | NHANES 2011-2012 (G) and 2013-2014 (H) (CDC NCHS) | public releases; PAM released Nov 2020; PAXMIN last-modified 2022-08-01; public-use Linked Mortality File 2019 | 2026-09-23 | participant (SEQN); participant-day and -minute for accelerometry | ingested |
| `stanford_longcovid_wearables` | Stanford Snyder Lab COVID-19 and Long COVID smartwatch datasets (Stanford University) | Mishra 2020 zip (updated 2020-10-19); Alavi 2022 zip (updated 2021-08-13); Uwakwe 2025 SDR druid:cb174pb4851 v3 (modified 2026-03-18) | 2026-09-23 | participant x minute / day / participant | ingested |
| `mapmecfs_nih_pi_mecfs` | NIH intramural post-infectious ME/CFS deep phenotyping, Walitt et al. 2024 (NINDS; mapMECFS/RTI; GEO/SRA; Pennsieve) | doi:10.1038/s41467-024-45107-3; GEO SuperSeries GSE251792 (last update 2024-03-06) | 2026-09-23 | participant (NIH study number) x omics modality x feature | **partial** |
| `nhanes_2003_2006` *new* | NHANES 2003-2004 (C) and 2005-2006 (D), hip accelerometry cohort (CDC NCHS) | public releases; PAXRAW_C/D.zip (2014-12-31); public-use Linked Mortality File 2019 | 2026-09-24 | participant (SEQN); participant-minute for hip accelerometry | ingested |
| `charlton_lc_mecfs_cpet_source` *new* | MUSCLE-ME Source Data, Charlton et al., Nat Commun 2026 (Springer Nature; Amsterdam UMC, Vrije Universiteit Amsterdam) | 41467_2026_75725_MOESM4_ESM.xlsx (Last-Modified 2026-07-28), rows `Group == 'MUSCLE-ME'` | 2026-09-24 | person (cross-sectional) | ingested |
| `appelman_lc_pem_source` *new* | Appelman et al. 2024 Long COVID PEM study Source Data, Nat Commun 15:17 (Springer Nature; Amsterdam UMC, Vrije Universiteit Amsterdam) | 41467_2023_44432_MOESM4_ESM.xlsx (Last-Modified 2025-02-19) | 2026-09-24 | person x timepoint (baseline, 1 day and 1 week after PEM induction) | ingested |
| `fm_thermography` *new* | Infrared thermography in fibromyalgia, PLOS ONE 2021 0253281 S1 Data (PLOS ONE / figshare; University of Valencia) | figshare file 28433660 (static supplement, 2021-06-24) | 2026-09-24 | participant (one resting thermography session) | ingested |
| `klein2023_mylc_ml_table` *new* | MY-LC per-participant table, Klein et al. Nature 2023, Suppl. Table 3 (Springer Nature; Yale, Icahn School of Medicine at Mount Sinai) | 41586_2023_6651_MOESM4_ESM.xlsx (static supplement, 2023) | 2026-09-24 | participant | ingested |
| `endo_arg1_repod` *new* | Serum arginase in endometriosis vs controls, RepOD PY1P9X (RepOD, University of Warsaw; Medical University of Warsaw) | endo_arg1.tab (static deposit); licence not stated in the file API | 2026-09-24 | participant | ingested |
| `heds_hsd_olink_serum_cinquina2026` *new* | Serum Olink proteome in hEDS/HSD, Cinquina et al., Clinical Proteomics 2026 (Europe PMC PMC13081554; Universita degli Studi di Brescia; The Ehlers-Danlos Society) | supplementary files (static, published 2026-03-07) | 2026-09-24 | participant (one serum sample) | ingested |
| `geo_cohorts` *new* | NCBI GEO public case/control omics cohorts, Long COVID and ME/CFS (NCBI / NLM; per-series submitters) | series as deposited on the retrieval date (sha256 per file in `MANIFEST.json`) | 2026-09-24 | person (one sample per subject) within one GEO cohort | ingested |

**Layer 2: condition-level molecular**

| source_id | source (publisher) | version / vintage | retrieved | unit of observation | status |
|---|---|---|---|---|---|
| `open_targets` | Open Targets Platform GraphQL API v4 | data release 26.06; API 26.6.3 | 2026-09-23 | target-disease association / evidence / drug / pathway | ingested |
| `gwas_catalog` | NHGRI-EBI GWAS Catalog REST API v2 | data release 2026-09-13; EFO v3.93.0; GRCh38.p14 | 2026-09-23 | curated SNP-trait association; study | ingested |
| `ncbi_geo_sra` | NCBI GEO / SRA E-utilities (metadata only) | gds last update 2026-09-22; sra 2026-09-23 | 2026-09-23 | GEO series / SRA study x condition | ingested |
| `reactome` | Reactome pathway knowledgebase | release 97 (versioned archive) | 2026-09-24 | pathway | ingested |
| `published_lab_evidence` *new* | Published lab and biomarker evidence for invisible illnesses, curated extraction (peer-reviewed journals via Europe PMC; curated by this project) | curated 2026-09-24 | 2026-09-24 | one published claim (publication x condition x biomarker x metric) | ingested |

**Layer 3: population / geographic**

| source_id | source (publisher) | version / vintage | retrieved | unit of observation | status |
|---|---|---|---|---|---|
| `cdc_long_covid` | Post-COVID Conditions, Household Pulse Survey (CDC NCHS / Census) | data.cdc.gov gsea-w83j, rows updated 2024-10-04; HPS phases 3.5-4.2 | 2026-09-23 | indicator x state/national x subgroup x period | ingested |
| `cdc_places` | CDC PLACES 2025 release (county + ZCTA) | BRFSS 2023 (35 measures), 2022 (5) | 2026-09-23 | county / ZCTA x measure | ingested |
| `cms_mmd` | CMS Mapping Medicare Disparities, Population View (CMS OMH) | claims years 2012-2022 final, 2023 preliminary | 2026-09-23 | geography x year x condition x measure (Medicare FFS) | ingested |
| `cdc_lyme` | Lyme disease cases by county of residence (CDC NCEZID) | cases 2001-2023 (file last-modified 2026-09-08) | 2026-09-23 | county x year | ingested |
| `cdc_svi` | CDC/ATSDR Social Vulnerability Index | SVI 2022 (ACS 2018-2022) | 2026-09-23 | county | ingested |
| `census_acs` | American Community Survey 5-year, Summary File (Census Bureau) | 2020-2024 | 2026-09-23 | state / county / ZCTA | ingested |
| `census_geography` | Cartographic boundaries, Gazetteer, ZCTA-county relationship (Census Bureau) | 2024 boundaries; 2020 relationship file | 2026-09-23 | state / county / ZCTA | ingested |
| `hrsa_hpsa` | HPSA and MUA/P designations (HRSA) | Data Warehouse extract 2026-09-23 | 2026-09-23 | shortage designation, summarised to county | ingested |
| `usda_rucc` | Rural-Urban Continuum Codes (USDA ERS) | RUCC 2023 | 2026-09-23 | county | ingested |

**Layer 4: facility / provider / research infrastructure**

| source_id | source (publisher) | version / vintage | retrieved | unit of observation | status |
|---|---|---|---|---|---|
| `nppes` | NPPES NPI downloadable file V.2 (CMS) | September 2026 V2 (npidata 20050523-20260913) | 2026-09-23 | NPI | ingested |
| `nucc_taxonomy` | Health Care Provider Taxonomy code set (NUCC) | 26.1 | 2026-09-23 | taxonomy code | ingested |
| `hrsa_health_centers` | Health Center Service Delivery and Look-Alike Sites (HRSA BPHC) | Data Warehouse extract 2026-09-23 | 2026-09-23 | health-center site | ingested |
| `clinicaltrials_gov` | ClinicalTrials.gov API v2 (NLM) | API 2.0.5; dataTimestamp 2026-09-23T09:00:05Z | 2026-09-23 | registered study; per-site rows | ingested |
| `nih_reporter` | NIH RePORTER Project API v2 (NIH OER) | latest date_added 2026-09-19 | 2026-09-23 | fiscal-year award record | ingested |
| `openfda_device` | openFDA device endpoints: classification, 510(k) incl. De Novo, PMA, registration and listing (FDA) | last_updated 2026-09-14; classification bulk 2026-09-23 | 2026-09-23 | product code; decision | ingested |

**Ontology, measurement and metadata**

| source_id | source (publisher) | version / vintage | retrieved | unit of observation | status |
|---|---|---|---|---|---|
| `published_device_evidence` *new* | Published device and objective-test evidence for invisible illnesses, curated extraction (peer-reviewed journals via Europe PMC; curated by this project; `data_layer = measurement`) | curated 2026-09-24 | 2026-09-24 | one published claim (publication x condition x device x metric) | ingested |
| `mondo_ontology` | Mondo Disease Ontology + Monarch KG (Monarch Initiative) | Mondo v2026-09-01; Monarch KG 2026-09-02 | 2026-09-23 | disease class; disease-phenotype association | ingested |
| `ebi_ols4` | EMBL-EBI Ontology Lookup Service | mondo 2026-09-01; efo 3.94.0; mesh 2025; hp 2026-09-01; doid 2026-08-31; ordo 4.9 | 2026-09-23 | ontology term | ingested |
| `cms_icd10cm` | ICD-10-CM code descriptions (CMS / NCHS) | FY2026 (April 1, 2026 update) | 2026-09-23 | code | ingested |
| `hgnc` | HGNC complete set | Last-Modified 2026-09-18 | 2026-09-23 | gene | ingested |
| `data_gov_catalog` | Data.gov catalog, metadata only (GSA) | Catalog API 0.1.0; DCAT-US | 2026-09-23 | catalog record per query | ingested |

---

## 4. Actual sample sizes

Row counts are from the processed tables (`docs/SPEC_COMPLIANCE.md` §2, §10; per-analysis flows from the result files named).

**Person-level (layer 1).** Namespaced ids (`nhanes:<SEQN>`, ...); no cross-dataset participant key exists by design.

| dataset | step | n | file |
|---|---|---|---|
| NHANES 2011-2014 | participants (DEMO) | 19,931 | `nhanes_cohort_flow.csv` |
| | with accelerometry data | 14,693 | |
| | passing the valid-wear rule (>= 4 valid days) | 12,955 | |
| | aged >= 18 | 8,954 | |
| | **analysis population** (excluding 88 pregnant) | **8,866** | |
| | targets (cases): functional limitation 8,423 (2,340); fatigue 8,230 (1,313); fair/poor health 8,295 (1,888); PHQ-9 >= 10 8,200 (737); ME/CFS-like proxy 7,796 (132); mortality 8,850 (770) | | `nhanes_targets.csv` |
| | prescription reason-code groups (2013-14 only): fibromyalgia 34, migraine 69, IBS 12, insomnia 146, myalgia 90 | | |
| Stanford Mishra 2020 | participants | 118 | `participants` |
| Stanford Alavi 2022 | participants | 2,123 | |
| Stanford Uwakwe 2025 | participants (Long COVID / not) | 126 (31 / 95) | `stanford_reproduction_benchmark.csv` |
| Stanford acute analyses | COVID-19-positive with a reference day / valid baseline / placebo-eligible | 116 / 107 / 84 (5,489 placebo onsets) | `results/WEARABLE_STANFORD_RESULTS.md` §B |
| mapMECFS | NIH study numbers with open omics | 47 (24 healthy volunteers, 23 PI-ME/CFS) | `docs/MAP_MECFS_LINKAGE_AUDIT.md` |
| | person-linked omics analysis, PI-ME/CFS / healthy volunteers per layer: CSF and serum SomaScan 21 / 21 (the same 42 people); PBMC RNA-seq 12 / 15; muscle RNA-seq 13 / 12; CSF metabolomics 10 / 11 | | `results/LINKED_OMICS_RESULTS.md` §2 |
| NHANES 2011-2014 labs vs diagnosis | adults >= 20, not pregnant, biochemistry measured | 10,106 | `results/LABS_VS_DIAGNOSIS_RESULTS.md` |
| | umbrella labels (cases): Rx fibromyalgia 36, Rx migraine 81, Rx IBS 18, Rx insomnia 186, Rx myalgia 104, ME/CFS-like proxy 160 | | |
| NHANES 2011-2014 lab layers | wearable-valid adults / with core labs and every earlier layer / with every full-sample lab layer and every earlier layer | 8,866 / 7,210 / 6,824 (58 lab layers) | `results/NHANES_LAYERS_RESULTS.md` §1 |
| NHANES 2003-2006 *new* | participants (DEMO) / with PAXRAW minutes / valid hip wear / aged >= 18 / **analysis population** (excluding 255 pregnant) | 20,470 / 14,631 / 9,691 / 6,538 / **6,283** | `results/NHANES_2003_2006_RESULTS.md` §2 |
| | targets (cases): functional limitation 5,825 (1,662); fatigue 2,885 (374); fair/poor health 5,916 (1,253); PHQ-9 >= 10 2,867 (152); ME/CFS-like proxy 2,673 (36); mortality 6,275 (1,448) | | |
| MUSCLE-ME (Charlton 2026) *new* | ingested: Long COVID / ME/CFS / healthy controls | 81: 25 / 26 / 30 | `data/raw/charlton_lc_mecfs_cpet_source/registry_entry.yaml` |
| | with daily steps (patients / controls); with VO2_rel | 73 (51 / 22); 79 | `results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md` |
| Appelman 2024 *new* | baseline plasma metabolomics Long COVID / healthy; baseline muscle; VO2max and plasma | 25 / 21; 25 / 19; 23 / 21 | `results/LABS_VS_DIAGNOSIS_RESULTS.md` Part 2 |
| Fibromyalgia thermography *new* | fibromyalgia / controls (all women); with the six regional averages | 86 / 92; 85 / 92 | `results/FM_THERMOGRAPHY_RESULTS.md` |
| MY-LC (Klein 2023) *new* | rows; analysis set Long COVID / healthy / convalescent controls; primary sample (complete cortisol and covariates) | 185; 99 / 40 / 39; 98 vs 77 | `results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md` |
| Endometriosis ARG1 (RepOD) *new* | rows; with pre-operative ARG1: endometriosis / surgical controls / healthy controls | 215; 120 / 32 / 53 | `results/ENDO_ARG1_REPOD_RESULTS.md`, registry fragment |
| hEDS/HSD Olink (Cinquina 2026) *new* | hEDS / HSD / healthy controls (controls all Italian; half the patients American); primary Italian patients vs controls | 88 / 88 / 176; 88 vs 176 | `results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md` |
| GEO cohorts *new* | series screened / passing rows / analysed subsets; cases + controls over the 7 subsets | 88 / 12 / 7; 252 + 234 | `data/raw/geo_cohorts/DATA_AUDIT.md` |
| | per subset, cases / controls: GSE226260_wb 17 / 20; GSE270045 19 / 17; GSE275334_lc 15 / 18; GSE293840 (cfRNA) 93 / 75; GSE156792 (methylation) 61 / 48; GSE16059 (twins) 32 / 32; GSE227375 (replication only) 15 / 24 | | `geo_cohort_classifier.csv` |
| All person-level tables | `participants` / `participant_wearable_features` / `participant_labs` | 22,345 / 15,315 / 735,587 in the build of `20260924T150209Z`; after the union with the §14 and NHANES 2003-2006 partitions: 43,826 / 30,019 / 1,414,960 (11 datasets in `participants`, never pooled) | |

**Condition-level molecular (layer 2).** `condition_molecular_evidence`: 70,295 rows = Open Targets 52,164 + GWAS
Catalog 2,772 + GEO 537 + SRA 257 + mapMECFS published analytes 14,565. `condition_registry`: 14 conditions; 298
ontology mappings.

**Geographic (layer 3).** 3,144 counties + 51 states/DC in the ranking universe (2024 vintage; Connecticut = 9 planning
regions). `geo_condition_burden` 558,246 rows; `geo_condition_features` 46,074 rows (3,291 geographies x 14
conditions); `geo_vulnerability` 3,144 counties.

**Facility (layer 4).** `providers` 1,456,346 NPIs; `facilities` 329,019; `clinic_registry` 321,536;
`research_site_registry` 8,722; HRSA sites 19,266; `clinical_trials` 7,910 (`trial_sites` 47,482, 27,828 U.S.);
`nih_projects` 9,342; FDA: 90 classification codes, 8,101 510(k), 132 PMA records.

**Derived.** `deployment_opportunities` 76,680 rows (48 condition-or-set x measurement x level combinations);
2,395 counties ranked for the primary query; `deployment_candidates` 10.

**Required final datasets** (`validate`: 11/11 present and non-empty):

| dataset | rows x columns |
|---|---|
| participant_wearable_features.parquet | 30,019 x 135 |
| participant_clinical_features.parquet | 40,401 x 174 |
| condition_registry.parquet | 14 x 50 |
| condition_molecular_evidence.parquet | 70,295 x 187 |
| condition_measurement_evidence.parquet | 420 x 98 |
| measurement_registry.parquet | 30 x 127 |
| geo_condition_features.parquet | 46,074 x 112 |
| clinic_registry.parquet | 321,536 x 57 |
| research_site_registry.parquet | 8,722 x 66 |
| deployment_opportunities.parquet | 76,680 x 140 |
| measure_it_public.duckdb | 1,410 MiB (1,478,504,448 bytes after run `20260924T222644Z`); 21 SPEC tables as `spec.*` views with the 8 provenance columns |

The two person-level rows include every partition of the digital-person unions (NHANES 2011-2014 and 2003-2006, the
Stanford releases, MUSCLE-ME for wearable features; NHANES 2011-2014 and 2003-2006 for clinical features); before the
§14 and NHANES 2003-2006 partitions they were 15,315 x 105 and 19,931 x 172 (run `20260924T150209Z`). The other rows are
unchanged. The datasets stay separate partitions with namespaced ids; none is pooled with another.

---

## 5. MapME/CFS linkage audit (partial)

Full audit: `docs/MAP_MECFS_LINKAGE_AUDIT.md` (generated by `measure_it.ingestion.mapmecfs`). Study: Walitt et al.,
Nat Commun 15, 907 (2024).

| question (SPEC) | answer |
|---|---|
| Downloaded datasets | GEO GSE251872 (PBMC RNA-seq), GSE245661 (muscle RNA-seq), GSE251790 (CSF SomaLogic), GSE254030 (serum SomaLogic); SRA PRJNA954397 metadata (stool metagenomics); Nature Source Data and Supplementary Data; the authors' GitHub repository; Pennsieve dataset 356 manifest |
| Not downloadable | every mapMECFS data file: 29 of 29 anonymous probes returned HTTP 403 "Login Required"; no account was created, so 85 mapMECFS files are catalogued as not automatable |
| Participant identifiers | NIH study number (1xx healthy volunteer, 3xx PI-ME/CFS), SomaLogic `Fatigue-NNN` (crosswalk published in GSE251790), microbiome `SID_NNN`, mapMECFS `map000NNN`, figure-local letters, group labels only |
| Subjects per modality (with ids) | PBMC RNA-seq 27; muscle RNA-seq 25; CSF SomaLogic 42; serum SomaLogic 42; CSF metabolomics 21; stool metagenomics 32 (own namespace) |
| Exact overlaps (study numbers) | PBMC x muscle 17; PBMC x CSF SomaLogic 27; muscle x CSF SomaLogic 21; CSF x serum SomaLogic 42; CSF SomaLogic x CSF metabolomics 19; 9 study numbers in all 5 modalities; union 47 |
| Wearable / actigraphy, HRV, tilt, CPET | in open files these carry **group labels only** (e.g. 24-h HRV: 33 rows, HV 19 / PI-ME/CFS 14, no id); the per-participant files exist only on mapMECFS (gated) |
| Can wearable and omics be linked? | **No, not with open data: 0 wearable x omics pairs have comparable identifiers.** `participant_wearable_features__mapmecfs` is not produced |
| What was built | a true omics-to-omics linked subset (`participant_omics_linked__mapmecfs`, 902,186 rows, 47 participants) and `condition_molecular_evidence__mapmecfs` (14,565 rows, labelled condition-level) |
| Pipeline without mapMECFS access | works: step `ingest_mapmecfs` runs offline with `fetch=False` |

Consistency checks recorded in the audit include: study-number prefix agrees with the deposited group in 27/27, 25/25
and 42/42 samples; SomaLogic crosswalk sex, age and group agree in 42/42; sex is concordant across deposits for 46 of
47 participants (MECFS_311 flagged). The paper's analytic cohort is 21 healthy volunteers + 17 PI-ME/CFS; the deposits
contain participants outside it and no open file flags which (this statement could not be checked against data:
`number_audit.csv` marks it unverifiable).

---

## 6. Wearable methods

Plans written before the analyses: `docs/ANALYSIS_PLAN_NHANES.md`, `docs/ANALYSIS_PLAN_STANFORD.md`; deviations are
listed at the end of each results file. Seed `20260923` throughout.

### 6.1 NHANES 2011-2014 wrist accelerometry (no heart rate exists)

| item | method |
|---|---|
| Input | released minute (PAXMIN), hour and day (PAXDAY) tables, MIMS units; joins only on SEQN (`participant_id = nhanes:<SEQN>`) |
| Valid wear | >= 4 valid days, a valid day being a full 1,440-minute recorded day with >= 600 valid wake-wear minutes |
| Features used in models (26) | activity volume (mean daily MIMS, active and sedentary minutes, sedentary fraction), day-to-day variability (SD, CV), weekend-weekday difference, wake/sleep-wear minutes, rest-activity rhythm (interdaily stability, intradaily variability, M10, L5, onsets, relative amplitude), fragmentation (active-to-sedentary and sedentary-to-active transition probabilities, bout lengths), active-period duration, algorithm sleep-proxy timing, regularity and duration (`nhanes_wearable_feature_dictionary.csv`) |
| Not derived | heart rate, HRV or any autonomic signal (none in NHANES 2011-2014); the 60-second exam pulse is a clinical variable, never a wearable one; raw 80-Hz mode not implemented (optional) |
| Feature sets | demographics (9), wearable (26), demographics + wearable (35), clinical (49-50: demographics, BMI, waist, BP, exam pulse, 32 labs, comorbidity and medication counts, smoking), clinical + wearable (75-76); target-defining items never used as predictors (`nhanes_predictor_exclusions.csv`) |
| Models | unpenalised logistic (`lr`); penalised logistic tuned inside training folds (`pen_lr`, **primary**); histogram gradient boosting (`hgb`) |
| Validation | participant-level repeated stratified 5-fold CV x 5 repeats; temporal split (train 2011-12, test 2013-14); 1,000 participant-bootstrap resamples, paired for increments; 200 label permutations; 20 wearable-row shuffles; AUROC, AUPRC, Brier, calibration intercept/slope, sensitivity/specificity at the training-fold prevalence |
| Unsupervised | PCA; UMAP for visualisation only; GMM by BIC; HDBSCAN (plus 4 sensitivity settings) |
| Signatures | age/age^2/sex-adjusted standardised differences per feature, HC3 CIs, BH-FDR over 26 features per phenotype |

### 6.2 Stanford consumer-wearable cohorts

| item | method |
|---|---|
| Releases | Mishra 2020 and Alavi 2022 (acute COVID-19 cohorts: minute heart rate and steps; **not** Long COVID labels); Uwakwe 2025 (the only Long COVID label among the Stanford releases processed: self-reported symptoms >= 12 weeks after confirmed COVID-19) |
| Data dictionary first | the DATA_AUDIT "Discovery" table records which release carries which label (`data/raw/stanford_longcovid_wearables/DATA_AUDIT.md`) |
| Uwakwe features (7 physiology + 3 wear-pattern) | whole-record resting-HR summaries of the authors' step-masked HR: mean, night mean, daily SD, IQR, day-to-day RMSSD of daily RHR (not beat-to-beat HRV), 30-day trend, CuSum alarms per 30 days; no infection date and no steps are released, so nothing is aligned to infection. Feature counts elsewhere: the models use the 7 physiology features (plus the 3 wear-pattern features as an artefact check: 10); the signature table tests each of the 7 unadjusted and age/sex-adjusted (14 tests: the "14 of 14" of the example JSON); the `WearableHeartRateAdapter` embeds a different set of 15 longitudinal features (Figure 10) |
| Uwakwe models | logistic, L2-logistic (primary), random forest; 10x5 repeated stratified CV; 2,000 participant-bootstrap resamples; 1,000 label permutations (full CV re-run each) and wearable-block permutations |
| Acute cohorts | personal baseline (>= 7 RHR days in days -42..-15), daily RHR and step z-scores; each participant's own **placebo onsets** before the infection as the negative control; participant-grouped CV for onset-vs-placebo windows; Kaplan-Meier recovery time |
| Leakage checks | 8 temporal-leakage and consistency checks, 0 violations (`stanford_temporal_leakage_checks.csv`) |

**Reproduction benchmark** (`stanford_reproduction_benchmark.csv`): the Mishra 2020 CuSum detector was ported from the
authors' R code and reproduced 220 of the 223 published alarms at the same timestamps; the pre-symptomatic detection
count is 14 of 24 vs the paper's 15 of 24 (the paper's alarm-to-episode rule is not stated; plausible rules give
14-18). Uwakwe 2025: cohort size and Table 1 agree exactly; the symptoms-only random forest gives ROC-AUC 0.862 vs the
paper's 0.907; the paper's heart-rate model (0.748) is **not reproducible** from the release because the
post-diagnosis 4-week window needs a diagnosis date that is not released. The Alavi 2022 NightSignal alerting was not
attempted.

### 6.3 Digital person and adapters

One Digital Phenotype Vector per dataset, never pooled (`docs/ADAPTERS.md` §8; `digital_person_summary.csv`): NHANES
adults 11,977 participants in 7 blocks; Stanford Alavi 2,123, Mishra 118, Uwakwe 126; mapMECFS 47 (omics-only). For
NHANES and Stanford the molecular route is labelled *condition-level molecular enrichment via the ontology*; mapMECFS is
labelled "participant-linked omics ... but NO wearable/physiology data is linkable in open data", so no
"Multimodal Person Representation" (wearable + linked omics) exists in this build. `WearableAdapter.embed` rebuilds all
38 stored accelerometry features from raw minutes for 523 participants (maximum absolute difference 2.9e-11;
`wearable_adapter_validation_summary.json`).

---

## 7. Phenotype modelling (Tests 1, 2 and 6 on person-level data)

### 7.1 NHANES: Test 1 (wearable vs demographic baseline) and Test 2 (increment over clinical)

`pen_lr`, 5x5 CV, AUROC with participant-bootstrap 95% CIs (`nhanes_model_performance.csv`, `nhanes_increments.csv`,
`nhanes_permutation_null.csv`; write-up `results/WEARABLE_NHANES_RESULTS.md` §1).

| target | n (cases) | wearable only | wearable - demographics | **Test 1**: demog.+wear - demog. | **Test 2**: clin.+wear - clin. | reading |
|---|---|---|---|---|---|---|
| Functional limitation | 8,423 (2,340) | 0.731 (0.718-0.743) | -0.0236 (-0.0364 to -0.0116) | +0.0320 (+0.0254 to +0.0383) | +0.0069 (+0.0040 to +0.0097) | T1 supported; T2 supported and robust (3/3 model families, temporal split, Holm p <= 0.012) |
| Fatigue (PHQ-9 item 4 >= 2) | 8,230 (1,313) | 0.614 (0.599-0.631) | -0.0252 (-0.0461 to -0.0037) | +0.0313 (+0.0219 to +0.0411) | +0.0046 (+0.0000 to +0.0093) | T1 supported; T2 **not robust** (Holm p 0.20; 1/3 families; temporal CI includes 0) |
| Fair/poor self-rated health | 8,295 (1,888) | 0.649 (0.636-0.664) | -0.0745 (-0.0908 to -0.0574) | +0.0192 (+0.0133 to +0.0251) | -0.0003 (-0.0017 to +0.0011) | T1 supported; T2 **null** |
| Depressive symptoms (PHQ-9 >= 10) | 8,200 (737) | 0.644 (0.623-0.663) | -0.0625 (-0.0890 to -0.0362) | +0.0279 (+0.0177 to +0.0377) | +0.0057 (+0.0005 to +0.0108) | T1 supported; T2 **not robust** (Holm p 0.17) |
| ME/CFS-like PROXY (not ME/CFS) | 7,796 (132) | 0.582 (0.536-0.627) | -0.0724 (-0.1352 to -0.0119) | +0.0139 (-0.0126 to +0.0387) | +0.0172 (-0.0037 to +0.0376) | T1 **null**; T2 **null** |
| All-cause mortality (positive control) | 8,850 (770) | 0.845 (0.830-0.860) | -0.0108 (-0.0243 to +0.0019) | +0.0193 (+0.0131 to +0.0254) | +0.0016 (-0.0012 to +0.0046) | control passed (the pipeline sees the known accelerometry-mortality signal over demographics); T2 AUROC null |

The other SPEC Phase 1 metrics for the same `pen_lr` 5x5 CV models (`nhanes_model_performance.csv`; 95% participant-bootstrap
CIs; AUPRC is compared with its no-skill baseline, the prevalence; sensitivity and specificity at the training-fold
prevalence as threshold; calibration intercept and slope of the cross-validated predictions):

| target | prevalence (AUPRC baseline) | AUPRC demographics | AUPRC wearable only | AUPRC demog.+wear | AUPRC clin.+wear | wearable-only sensitivity / specificity | clin.+wear calibration intercept / slope |
|---|---|---|---|---|---|---|---|
| Functional limitation | 0.278 | 0.512 (0.490-0.533) | 0.550 (0.529-0.571) | 0.604 (0.583-0.625) | 0.699 (0.680-0.718) | 0.645 (0.626-0.664) / 0.691 (0.679-0.703) | 0.002 / 0.990 (0.944-1.039) |
| Fatigue | 0.160 | 0.247 (0.230-0.267) | 0.232 (0.215-0.252) | 0.285 (0.264-0.309) | 0.356 (0.332-0.383) | 0.537 (0.512-0.562) / 0.622 (0.612-0.633) | -0.001 / 0.953 (0.888-1.023) |
| Fair/poor self-rated health | 0.228 | 0.435 (0.411-0.459) | 0.351 (0.332-0.373) | 0.455 (0.432-0.479) | 0.541 (0.519-0.567) | 0.593 (0.571-0.614) / 0.618 (0.607-0.629) | 0.002 / 0.975 (0.921-1.032) |
| Depressive symptoms | 0.090 | 0.183 (0.166-0.204) | 0.141 (0.127-0.157) | 0.209 (0.188-0.233) | 0.284 (0.254-0.318) | 0.590 (0.557-0.625) / 0.622 (0.612-0.632) | -0.002 / 0.979 (0.905-1.058) |
| **ME/CFS-like PROXY** | **0.017** | 0.029 (0.022-0.040) | **0.024 (0.019-0.033)** | 0.029 (0.023-0.039) | 0.035 (0.027-0.047) | 0.517 (0.440-0.593) / 0.593 (0.583-0.603) | -0.011 / 0.864 (0.683-1.040) |
| All-cause mortality (positive control) | 0.087 | 0.395 (0.360-0.431) | 0.424 (0.389-0.460) | 0.471 (0.434-0.510) | 0.522 (0.486-0.560) | 0.782 (0.753-0.809) / 0.776 (0.767-0.785) | 0.001 / 0.986 (0.931-1.044) |

For the one invisible-illness-adjacent label the precision is barely above chance: at a prevalence of 1.7%, the
wearable-only model's AUPRC is 0.024 and the best model's 0.035, and at the prevalence threshold it finds about half the
proxy cases (sensitivity 0.52) while flagging 41% of non-cases (specificity 0.59).

Readings:
* **Positive:** added to demographics, wearable features raise AUROC for 4 of 5 symptom/function targets and for the
  mortality positive control (+0.019 to +0.032).
* **Negative:** wearable features *alone* rank worse than a 9-variable demographic model for all 6 targets (point
  estimate); the 95% CI is below 0 for 5 of 6.
* **Small or absent over clinical data:** only functional limitation passes every robustness check, and that label is
  partly definitional (self-reported activity limitation, which device-measured activity is expected to track).
* Calibration of the clinical + wearable models is close to ideal in CV for the five common targets (slopes
  0.953-0.990); for the ME/CFS-like proxy the slope is 0.864 (0.683-1.040) and the wearable-only slope 0.811
  (0.294-1.288), i.e. poorly determined (`nhanes_model_performance.csv`; `results/WEARABLE_NHANES_RESULTS.md` §7).

**ME/CFS-like proxy signature** (132 cases vs 7,664; `nhanes_phenotype_signatures.csv`): 14 of 26 features pass FDR,
e.g. night-to-night SD of sleep midpoint +0.33 SD (+0.15 to +0.51), interdaily stability -0.32 (-0.50 to -0.14),
active-period duration -0.32 (-0.49 to -0.14). The sensitivity variants excluding depression (71 cases) or arthritis
(81 cases) have similar effect sizes but only 1 and 9 features pass FDR: they neither confirm nor rule out that the
signature is carried by co-occurring depressive symptoms or arthritis. 55.0% of proxy cases (of those with a PHQ-9)
have PHQ-9 >= 10 (`nhanes_proxy_characterisation.csv`).

**Unsupervised structure.** HDBSCAN finds **no cluster** (every participant labelled noise, also under 4 other
parameter settings); GMM BIC keeps improving to the k = 8 boundary with only moderately reproducible partitions
(bootstrap ARI mean 0.41). Activity patterns form a continuum, not discrete subtypes (`nhanes_hdbscan_sensitivity.csv`,
`nhanes_cluster_stability_summary.csv`).

**Prescription reason-code labels** (2013-14 only): fibromyalgia 13 of 26 FDR features (34 cases; drug-defined, may
reflect the medication); migraine 0 of 26 (69 cases) and IBS 0 of 26 (12 cases): **null**.

### 7.2 Stanford: the only Long COVID label with wearable heart rate among the public datasets processed (Uwakwe 2025)

`stanford_uwakwe_model_performance.csv`, `stanford_uwakwe_model_increments.csv`, `stanford_uwakwe_permutation_nulls.csv`.

| feature set (L2-logistic) | AUROC (95% CI) | label-permutation p | AUPRC (no-skill 0.246) | sensitivity / specificity | calibration intercept / slope |
|---|---|---|---|---|---|
| demographics / clinical baseline (age bin, sex, BMI, race, vaccination) | 0.702 (0.598-0.796) | < 0.001 | 0.449 (0.307-0.612) | 0.687 (0.536-0.823) / 0.648 (0.562-0.735) | -0.003 / 0.936 |
| wearable resting HR only | **0.504 (0.404-0.605)** | 0.300 | 0.288 (0.206-0.414) | 0.439 (0.300-0.583) / 0.635 (0.560-0.711) | 0.023 / -0.169 |
| demographics + wearable | 0.677 (0.573-0.769) | 0.008 | 0.420 (0.287-0.571) | 0.574 (0.432-0.711) / 0.672 (0.589-0.751) | 0.019 / 0.737 |
| increment (demographics + wearable - demographics) | **-0.024 (-0.076 to 0.027)** | wearable-block permutation p 0.469 | | | |

The wearable-only model's calibration slope is negative (-0.169): its predictions carry no usable ordering of risk.

**Test 1 and Test 2 are null on the invisible-illness label.** No model family rescues the wearable-only set
(logistic 0.523, random forest 0.557). 0 of 7 physiology features survive BH correction, and 0 of 7 reach p < 0.05
under the studentized permutation test. For scale, the 7 symptom items that *define* the label reach 0.857 (a
leakage-labelled reference, not a comparator for a measurement). With 31 cases the expected CI half-width was about
0.1 AUROC. The adapter's within-person RHR-elevation score against the label points the **opposite** way to the
hypothesis: AUROC 0.387, bootstrap CI 0.284-0.494 excludes 0.5, permutation p 0.069 (the file labels this
"inconclusive"); Long COVID participants also contribute fewer RHR days (median 355 vs 491), which the record-relative
score does not correct for (`wearable_hr_adapter_label_check.csv`).

### 7.3 Stanford acute COVID-19 cohorts (post-infection trajectory; NOT Long COVID)

`stanford_acute_endpoints.csv`, `stanford_acute_detection.csv`, `stanford_acute_recovery_km.csv`.

| endpoint (true onset vs the same people's placebo onsets) | estimate (95% CI) | reading |
|---|---|---|
| mean RHR z on day 0 | 1.03 (0.63 to 1.44) vs placebo band 0.01-0.13 | brief peak (descriptive) |
| day with RHR z > 2 in days -7..+14 | 62% true vs 56% placebo; difference 6.3 points (-5.3 to 19.7), p 0.332 | **null**: not specific to infection |
| mean RHR z, days 0..14 | difference +0.114 SD (-0.203 to 0.451), p 0.531 | **null** (pre-specified endpoint) |
| mean step z, days 0..14 | difference -1.103 SD (-1.334 to -0.880), p < 0.001 | supported: reduced walking during acute illness (expected illness behaviour) |
| onset-vs-placebo window A (days 0..7), model | AUROC 0.861 (0.792-0.917), p < 0.001 | carried by the step drop (mean RHR z alone 0.531, p 0.267) |
| pre-symptomatic window B (days -7..-1), model | AUROC 0.562 (0.473-0.646), p 0.115 | **null** |
| recovery to baseline band (from acute peak) | median 5 days (3-8) true vs 6 days (5-7) placebo "elevations" | not infection-specific |

### 7.4 Phase 3: what could make the phenotype visible

`results/PHASE3_MEASUREMENT_EVIDENCE.md` keeps five dimensions separate (phenotype_signal_strength,
measurement_evidence_strength, technology_maturity, regulatory_visibility, deployment_complexity) and computes no
composite. **18 of 25** phenotype x measurement-class cells have no public person-level signal result; at condition
level, 9 of 364 condition x class cells carry a verdict. None of the wearable datasets of §6 has a clinical diagnosis
label for any target condition, and no public dataset processed here has posture-linked heart rate, beat-to-beat HRV,
heart rate aligned to exertion, SpO2 or respiratory rate. This Phase 3 table predates §14 and does not include its
datasets, which do carry clinical labels (criteria-based Long COVID and diagnosed ME/CFS in MUSCLE-ME, ACR 2010
fibromyalgia with a clinic skin-temperature camera, among others).

| demo phenotype | public signal verdicts |
|---|---|
| Orthostatic / autonomic dysfunction | wearable_heart_rate **null** (Uwakwe Long COVID); 8 of 9 classes not observable |
| Activity intolerance / post-exertional pattern | accelerometry **mixed** (ME/CFS-like proxy: FDR features, null model increment; partly definitional); 4 of 5 classes not observable |
| Sleep / circadian disruption | accelerometry and sleep_objective **supported** on rx_insomnia (drug-defined label, signature only). The sleep_objective verdict rests on the NHANES wrist-accelerometry sleep **proxy** features: no processed dataset has PSG or wearable sleep staging, which is how the class is defined; wearable_heart_rate null |
| Post-infection physiological recovery | accelerometry **supported** (acute step drop); wearable_heart_rate **null** (acute RHR endpoints vs placebo) |

No wearable heart-rate result is `supported` for any condition or phenotype.

### 7.5 NHANES 2003-2006: replication with a different device and cohort

Plan `docs/ANALYSIS_PLAN_NHANES_2003_2006.md` (written before any 2003-2006 model was fitted); write-up
`results/NHANES_2003_2006_RESULTS.md`; tables `nhanes0306_replication.csv`, `nhanes0306_increments.csv`. The 2003-2006
cycles used a **hip-worn uniaxial ActiGraph 7164** (counts per minute, waking hours only), not the 2011-2014 wrist
MIMS; no feature is on a common scale and the cohorts are never pooled. What is compared is whether the same
pre-specified verdict holds. Analysis population 6,283 adults (§4); same models, CV, bootstrap and controls as §7.1.

| target | Test 1 (demog.+wear - demog.) 2011-2014 -> 2003-2006 | Test 2 (clin.+wear - clin.) 2011-2014 -> 2003-2006 | replicated (T1 / T2) |
|---|---|---|---|
| Functional limitation | +0.0320 -> +0.0193 (+0.0125 to +0.0266) | +0.0069 -> +0.0041 (+0.0005 to +0.0075) | yes / yes |
| Fatigue | +0.0313 -> +0.0286 (+0.0039 to +0.0509) | +0.0046 -> +0.0019 (-0.0102 to +0.0130) | yes / **no** |
| Fair/poor self-rated health | +0.0192 -> +0.0151 (+0.0086 to +0.0213) | -0.0003 -> +0.0012 (-0.0014 to +0.0035) | yes / yes (null both) |
| Depressive symptoms (PHQ-9 >= 10) | +0.0279 -> +0.0177 (-0.0144 to +0.0493) | +0.0057 -> +0.0035 (-0.0081 to +0.0147) | **no** / **no** |
| ME/CFS-like PROXY (36 cases in 2003-2006) | +0.0139 -> +0.0666 (+0.0220 to +0.1148) | +0.0172 -> +0.0560 (+0.0164 to +0.0987) | **no** / yes (null both: shuffle control fails) |
| All-cause mortality (positive control) | +0.0193 -> +0.0076 (+0.0047 to +0.0106) | +0.0016 -> +0.0009 (-0.0006 to +0.0024) | yes / yes (null both) |

**8 of 12 pre-specified verdicts replicate** (Test 1 4/6, Test 2 4/6). Test 1 holds for the well-powered targets.
The functional-limitation increment over clinical data reproduces in direction and size but is fragile: CI > 0 with
`pen_lr` and `hgb` but not `lr`, and the 2003-04 -> 2005-06 temporal Delta includes 0. Fatigue and depression Test 2
do not replicate (they were already not robust in 2011-2014). The ME/CFS-like proxy turns positive for Test 1 with 36
cases and a CI width of about 0.09; the better-powered 2011-2014 null (132 cases) stands, and the 2003-2006 signal is a
hypothesis, not a replication success or failure. Wearable-only models are again below or equal to demographics for
every target except the proxy; CRP adds nothing over the clinical set (all Delta AUROC within +/-0.002). Mortality
follow-up is about 6 years longer than for 2011-2014, so its AUROCs are not directly comparable.

---

## 8. Molecular enrichment (Test 3): condition-level, not patient multi-omics

Everything in this section is **condition-level molecular enrichment**: gene, variant, pathway and study lists that
public databases attach to a disease concept. Nothing is measured on, or linked to, any participant; no combined omics
score is computed. Write-up: `results/MOLECULAR_COHERENCE.md`; plan: `docs/ANALYSIS_PLAN_MOLECULAR.md`. Person-level
omics measured on the people of one cohort each (mapMECFS linked layers, GEO case/control cohorts, the MY-LC and hEDS/HSD
tables) are reported separately in §14; none is linked to a wearable participant.

**Sources** (`condition_molecular_evidence`, 70,295 rows): Open Targets 26.06 (target-disease 29,494; pathway 20,477;
known drug 1,082; text-mined 667; genetic association 444), GWAS Catalog 2026-09-13, curated GEO/SRA queries
(metadata only; `configs/geo_queries.yaml`), mapMECFS published differential analytes; Reactome 97 and HGNC for
over-representation analysis.

**Source counts for the demo conditions** (HGNC protein-coding genes; `test3_source_counts.csv`):

| condition | OT genetic | GWAS Catalog | OT literature (Europe PMC co-mention) | OT other (mostly drug targets) | mapMECFS | GEO series |
|---|---|---|---|---|---|---|
| long_covid | 2 | 0 | 300 | 46 | 0 | 52 |
| me_cfs | 70 | 9 | 1,587 | 109 | 1,387 | 28 |
| pots | 1 | 0 | 135 | 101 | 0 | 0 |
| dysautonomia | 28 | 1 | 313 | 42 | 0 | 49 |

MCAS and post-infectious syndrome have no protein-coding gene in any source (reported as UNKNOWN downstream, never imputed).

**Agreement across databases** (same-condition vs cross-condition overlap, condition-label permutation;
`test3_gene_agreement_control.csv`, `test3_pathway_agreement_control.csv`):

| source pair | gene level: difference (95% CI), perm. p | pathway level (Reactome): perm. p |
|---|---|---|
| OT genetic x GWAS Catalog (7 / 5 conditions) | 2.604 (1.498 to 3.560), 2.0e-04 | 0.017: above baseline |
| OT genetic x OT literature | 1.039 (0.428 to 1.615), 1.0e-04 | 4.7e-04: above |
| OT genetic x OT other | 1.228 (0.436 to 2.072), 1.0e-04 | 5.0e-05: above |
| GWAS Catalog x OT literature | 0.700 (**-0.072** to 1.327), 0.0063 | 0.396: **not** above |
| GWAS Catalog x OT other | -0.101 (-0.675 to 0.379), 0.627: **not** above | 0.836: **not** above |
| OT literature x OT other | 1.008 (0.571 to 1.473), 5.0e-05 | 5.0e-05: above |

Gene level: 5 of 6 source pairs agree more for the same condition than across conditions (for one of the five the
bootstrap CI includes 0). Pathway level: 4 of 6. These sources share curated knowledge and literature, so agreement is
not independent replication.

**Negative and cautionary results**
* The mapMECFS ME/CFS gene/protein list is **not** the best-matching condition in any of its 4 source comparisons
  (ranks 9/10, 7/8, 3/12, 4/12; `test3_mapmecfs_rank.csv`).
* Curated clinical flags vs data-supported systems: the pre-specified "any source" definition is not above chance
  (20 of 28 expected pairs, p 0.339); the "condition-specific" definition is nominal (18 of 28, Holm p 0.0495), but of
  its 18 concordant pairs 13 include Open Targets 'other' support (10 of them post-hoc flagged as mostly drug targets
  of trialled drugs; the 3 unflagged are eds_hsd:vascular_endothelial, lyme_disease:immune and migraine:neuronal), 5
  rest on literature co-mention alone and 0 involve genetic evidence (`test3_flag_concordance.csv`).
  Restricted to non-literature, unflagged support: those 3 of 28 (p 0.450). **Public molecular resources do not
  independently recover the curated clinical profile.**
* **Autonomic links are drug-target-driven.** For ME/CFS and POTS, the autonomic/cardiac, neuronal and
  mitochondrial/metabolic system links come from Open Targets 'other' evidence in which the genes are drug targets of
  trialled drugs (e.g. gabapentinoids: voltage-gated calcium channel family; metformin: complex I). One drug mechanism
  can supply the same system to several conditions (calcium-channel modulators: 6 condition-system links). After
  removing post-hoc-flagged sources, **none of Long COVID, ME/CFS or POTS has a system with non-literature support**; Long
  COVID's systems (immune, vascular/endothelial, connective tissue) are literature co-mention only. The `ask` agent
  and `get_molecular_context` expose this per system: the agent now prints each condition's systems in three tiers
  (non-literature unflagged / carried only by post-hoc-flagged sources / literature co-mention only; §16). Before
  this report it printed every enriched system as "supported", an overstatement in user-facing output.
* Open Targets literature evidence supports the immune system for 11 of 11 conditions tested (research-attention bias).

The `measurable_biology` graph (247 rows; `fig5_nodes.csv` 238 nodes, `fig5_edges.csv` 475 edges) links condition ->
system (data-derived) -> measurement class (curated assumption). A curated system -> class link is not evidence that
the measurement detects the condition.

---

## 9. Technology discovery

Write-up: `results/MEASUREMENT_DISCOVERY.md`; plan: `docs/ANALYSIS_PLAN_MEASUREMENT.md`. A trial using a measurement is
research activity, not evidence that it works; an FDA record is a deployment-readiness signal, not evidence that a
device diagnoses the target illness.

**Method.** Pattern text mining of ClinicalTrials.gov (7,910 trials: interventions, outcomes, eligibility, sites) and
NIH RePORTER (9,342 award records: titles and abstracts) for 26 curated measurement classes and 4 bundles
(`configs/measurements.yaml`, `configs/relevance.yaml`), with roles (outcome, intervention, eligibility), an objective
flag, a questionnaire filter and a non-human-context filter: 6,379 trial and 4,267 grant mentions after revision
(10,646; 16,355 with the patterns as first shipped). openFDA classification, 510(k)/De Novo, PMA and registration data
were mapped to classes with a mapping confidence (`measurement_regulatory_status`, 102 rows). The discovery runs
unchanged for any registry condition, but it can only find technology types inside the curated class vocabulary.

**Demo cluster** (objective-use trials / as outcome measure / NIH core projects; literal condition match):

| class | Long COVID | ME/CFS | POTS | dysautonomia | text-mining precision (record-clustered 95% CI; status) | tech maturity | deploy cx (1-5, curated) |
|---|---|---|---|---|---|---|---|
| blood_biomarkers | 110 / 91 / 38 | 32 / 27 / 26 | 12 / 9 / 0 | 38 / 34 / 22 | 0.83 (0.66-0.93; provisional) | established | 1 |
| autonomic_testing | 25 / 23 / 3 | 13 / 11 / 2 | 50 / 36 / 7 | 88 / 58 / 21 | 0.97 (0.83-0.99; meets) | established | 4 |
| hrv | 48 / 42 / 1 | 25 / 22 / 2 | 8 / 7 / 1 | 88 / 75 / 31 | 0.93 (0.78-0.98; provisional) | emerging | 2 |
| cpet | 53 / 47 / 2 | 11 / 10 / 5 | 7 / 5 / 0 | 11 / 6 / 3 | 1.00 (0.79-1.00; provisional) | established | 4 |
| accelerometry | 39 / 36 / 2 | 17 / 14 / 3 | 6 / 4 / 0 | 9 / 7 / 3 | 0.97 (0.78-1.00; provisional) | established | 1 |
| wearable_heart_rate | 14 / 9 / 0 | 5 / 5 / 0 | 1 / 1 / 0 | 9 / 7 / 2 | 0.85 (0.62-0.95; provisional) | established | 1 |
| ecg_ambulatory | 7 / 4 / 0 | 2 / 2 / 0 | 1 / 0 / 0 | 21 / 15 / 3 | 0.93 (0.67-0.99; provisional) | established | 2 |
| ppg | 5 / 5 / 0 | 1 / 1 / 0 | 1 / 1 / 0 | 3 / 3 / 1 | 1.00 (0.46-1.00; insufficient sample) | established | 1 |
| posture_detection | 3 / 2 / 0 | 0 / 0 / 0 | 1 / 1 / 0 | 3 / 3 / 0 | 0.00 (0.00-0.42; **below 0.8**) | no FDA category | 2 |
| capillaroscopy | 3 / 3 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 1 / 1 / 0 | 1.00 (0.43-1.00; insufficient sample) | no FDA category | 2 |

Precision status: "meets" only when the record-clustered lower bound is at least 0.8 (1 of the 10 classes shown:
autonomic_testing); "provisional" when the point estimate is at least 0.8 but the lower bound is not
(`measurement_precision_summary.csv`). Across all 26 classes, 16 have an adequate sample and a point estimate >= 0.8,
and only autonomic_testing clears 0.8 at the clustered lower bound.

Wearable autonomic/activity bundle: Long COVID 85 trials (13.0% of its 652 literal-match trials; 59, 9.0%, without the
device-agnostic HRV class), ME/CFS 36 (17.0%),
POTS 17 (12.2%), dysautonomia 103 (34.6%). Capillaroscopy (the future CAPRIO adapter) appears in 4 trials across all
14 conditions, 0 NIH core projects and has no FDA product code: a research-evidence gap, not an established deployment.

**Precision spot-check** (`measurement_precision_summary.csv`): a seeded sample of up to 30 mentions per class, from
demo-cluster records only, judged by an **AI reviewer (a language model, not a domain expert)** from the snippet alone.
These judgements are quality-assurance flags about the text mining, never data: no count or score uses them. Round 1:
15 of 26 classes below 0.8 strict precision. After one pattern revision and a fresh sample: 3 classes still below 0.8 (`immune_assays` 0.67,
`posture_detection` 0.00, `transcriptomics` 0.57) and 6 have too small a sample to call. The judgements are QA flags;
no score uses them. **Other cautions:** 144 of 298 dysautonomia trials link only through the broad term "autonomic
dysfunction" (through specific terms the bundle falls from 103 to 35 trials); the bundle's HRV member is device-agnostic
(laboratory, tilt table, Holter or wearable), so without HRV the bundle is 195 of 374 trials; the questionnaire filter
over-excludes (false-exclusion rate 0.77); recall is not measured, so a zero is not evidence of non-use.

**FDA context for the wearable bundle** (from the example recommendation's `technology.regulatory_context`):
accelerometry 3 mapped product codes (24 510(k) decisions), wearable_heart_rate 3 (343), hrv 2 (5; low mapping
confidence), ecg_ambulatory 10 (1,649), ppg 4 (59), posture_detection none found.

---

## 10. Geographic methods (Test 4)

Write-up: `results/GEOGRAPHY_RESULTS.md`; definitions pre-specified in `docs/BURDEN_DEFINITIONS.md`. All geographic
joins are ecological; nothing describes an individual.

**Burden evidence levels** (A direct measure; B closely matching coded condition; C symptom/comorbidity proxy; D no
usable burden estimate). Every burden row carries its level, its source resolution and an `inherited` flag.

| condition | level | primary measure | resolution |
|---|---|---|---|
| long_covid | A (state) | HPS "currently experiencing long COVID", % of adults, latest non-suppressed period | state; **county rows carry the state value, `inherited = true`** ("STATE value inherited as county context; not county prevalence") |
| lyme_disease | A | reported cases per 100k, 2023 | county (state = aggregated up) |
| fibromyalgia | B | CMS MMD 51 "Fibromyalgia, Chronic Pain and Fatigue", 2022 | county / state |
| migraine | B | CMS MMD 59 "Migraine and Other Chronic Headache", 2022 | county / state |
| me_cfs | **C** | CMS MMD 51 (contains R53.82 chronic fatigue, not G93.3x): a symptom proxy | county / state |
| ptlds | C | Lyme incidence (antecedent-exposure proxy) | county |
| pots, dysautonomia, eds_hsd, mcas, ibs, gastroparesis, endometriosis, post_infectious_syndrome | **D** | none; burden value null with a reason | - |

**Inherited state values.** No county long-COVID measure exists (HPS has none; PLACES 2025 has no COVID item). The
state value attached to 3,144 counties keeps `source_geographic_resolution = state`; no state estimate is downscaled.
It carries no within-state information (its Moran's I is unchanged by every within-state shuffle), and it is the burden
term of the long-COVID county diagnostic desert, flagged `desert_burden_resolution = 'state (inherited ...)'`. In
scoring its uncertainty is doubled and a per-county deviation is added (§12).

**Context and desert.** SVI 2022 (vulnerability, never burden); ACS 2020-2024; PLACES 2025; RUCC 2023; HPSA/MUA.
`diagnostic_desert = mean(pr(burden), pr(SVI), 1 - pr(relevant providers per 100k), 1 - pr(relevant trials within
50 km))`, percentile ranks within the 50 states + DC; an access-only variant is used for level-D conditions.

**Test 4: burden or population size?** (Spearman with population; `test4_population_correlations.csv`,
`test4_verdicts.csv`). County series show the i.i.d. 95% CI and, in brackets, the state-cluster bootstrap CI, which is
2.1-4.1 times wider; state series (n = 51) have i.i.d. CIs only.

| series | prevalence vs population | count vs population | reading |
|---|---|---|---|
| Long COVID, state (A) | -0.26 (-0.53 to 0.07; CI includes 0) | 0.95 (0.89 to 0.97) | prevalence not shown to depend on population (n = 51, CI includes 0); counts population-driven |
| Lyme, county (A) | 0.40 (0.37 to 0.43) [0.29 to 0.49] | 0.49 (0.47 to 0.52) [0.41 to 0.58] | mixed |
| CMS 51, county (B fibromyalgia / C ME/CFS) | 0.18 (0.14 to 0.21) [0.08 to 0.27] | no count (not defensible) | burden-driven |
| PLACES PHLTH county proxy (C) | -0.34 (-0.37 to -0.31) [-0.48 to -0.20] | 0.99 (0.99 to 1.00) [0.99 to 1.00] | counts population-driven; the prevalence top list is concentrated in small counties |
| diagnostic desert, 12 condition x level series | -0.36 to 0.31 | - | pre-specified verdicts: 10 burden-driven, 2 mixed (a rule that tests only the large-population tail; see below) |

**Reviewer correction (negative).** The pre-specified rule tests only the large-population tail. The county desert
top-25 lists contain 9-19 bottom-population-quartile counties (expected 6.3), significant after Holm for 3 of 6 series
(long COVID 19/25): **the county desert partly measures small population and rurality** through its access terms
(`test4_small_population_check.csv`). The desert's "burden-driven" verdicts therefore say only that large counties do
not dominate.

---

## 11. Facility matching

Write-up and plan: `docs/FACILITY_MATCHING.md`, `docs/ANALYSIS_PLAN_FACILITIES.md`.

**Registries.** `facilities` (329,019) resolves 369,339 NPPES organisation NPIs, 19,266 HRSA sites, 9,761 U.S.
ClinicalTrials.gov facility keys and 541 NIH organisations. `clinic_registry` (321,536; implementation-oriented) and
`research_site_registry` (8,722; research-oriented; 7,244 with a literal-match trial or title/abstract NIH match) keep
research readiness separate from burden. Coordinates: ZCTA centroid for 308,258 facilities, source point for 19,226,
city centroid for 308 (`docs/FACILITY_MATCHING.md`, geocode precision).

**Entity resolution quality.** V1 against HRSA's own site NPIs: precision **0.972**, recall **0.528** (precise, misses
about half of true links, so large academic centres can appear as several facilities). V2 location-shuffle false-link
share **0.0046** (a lower bound: it cannot see look-alike names in the same city).

**`find_candidate_clinics(geography, condition, measurement, radius_km=50)`.** Pool = geocoded facilities within the
radius of the county's Census internal point plus those inside the county. Six characteristics are **ranked
separately and never combined**: clinical_specialty_match, relevant_trial_history, NIH_research_activity,
technology_experience, community_access, distance_to_target_population. Candidates are drawn round-robin across the six
rankings, each with its raw values, ranks, plain-language reasons and the framing *"This facility has characteristics
suggesting it may be a viable implementation or study partner."* When no facility in the pool has trial history, NIH
activity or technology experience, the result carries an absence note ("a possible measurement/research desert ...
not proof that no such activity exists").

| example query (`facility_candidate_examples.csv`) | pool | eligible | with condition trials | with NIH (condition) | with technology experience |
|---|---|---|---|---|---|
| San Diego County x Long COVID x wearable monitoring | 2,615 | 1,664 | 6 | 5 | 3 |
| Cook County IL x ME/CFS x wearable monitoring | 8,135 | 5,159 | 2 | 3 | 0 |
| Harlan County KY x Long COVID x wearable monitoring | 216 | 184 | 0 | 0 | 0 |

**Matcher-level location shuffle** (201 counties, 200 national + 100 within-state draws; `test6_clinic_shuffle.csv`):
metro-nonmetro differences in Long COVID trials per candidate (0.145 vs null 0.093) and NIH share (0.031 vs 0.020)
survive BH-FDR and county resampling, but are small (about 0.05 trials per candidate). The earlier nonmetro
**technology-experience measurement desert is not supported** with the audited patterns (0.003 vs 0.005, p 0.26).
Primary care is a curated core specialty for every condition, so specialty match is a weak signal.

---

## 12. Scoring formula, weights and uncertainty

Method: `docs/SCORING.md`; plan and dated deviations: `docs/ANALYSIS_PLAN_SCORING.md`; results:
`results/SCORING_RESULTS.md`.

```
Opportunity_ws(g,c,m) = sum_k w_k * pr(component_k) / sum_k w_k      (over the available components)
components k: burden(g,c), vulnerability(g), diagnostic_desert(g,c), clinic_capacity(g,m), research_readiness(g,c,m)
pr() = percentile rank across all ranked geographies of the 50 states + DC (not within a state), per condition x
       measurement x level
```

| component | raw input |
|---|---|
| burden | primary burden measure (§10); a condition set = percentile of the mean of member percentiles, only where every member with a defined measure has a value |
| vulnerability | SVI 2022 overall |
| diagnostic_desert | the geography index (access-only for level D) |
| clinic_capacity | implementer-group NPIs and facilities per 100k, in the county and within 50 km (`configs/relevance.yaml`) |
| research_readiness | distinct condition trials, NIH core projects and technology-experience trials at facilities in the county or within 50 km |
| technology_saturation | reported only; used only in the `saturation_adjusted` weight set |

**Weight sets** (`configs/scoring.yaml`; default equal):

| weight set | burden | vulnerability | desert | clinic capacity | research readiness | 1 - saturation |
|---|---|---|---|---|---|---|
| equal (default) | 0.20 | 0.20 | 0.20 | 0.20 | 0.20 | 0 |
| burden_led | 0.40 | 0.15 | 0.15 | 0.15 | 0.15 | 0 |
| equity_led | 0.15 | 0.30 | 0.30 | 0.125 | 0.125 | 0 |
| capacity_led | 0.15 | 0.15 | 0.10 | 0.30 | 0.30 | 0 |
| study_partner | 0.25 | 0.10 | 0.10 | 0.20 | 0.35 | 0 |
| access_gap | 0.25 | 0.25 | 0.40 | 0.10 | 0 | 0 |
| burden_only | 1 | 0 | 0 | 0 | 0 | 0 |
| saturation_adjusted | 1/6 | 1/6 | 1/6 | 1/6 | 1/6 | 1/6 |

Since 2026-09-25 a ninth set, **`evidence_weighted`** (the pre-declared primary alternative; `equal` stays the default),
adds two absolute-scale components: `measurement_evidence` = tier factor x clip((AUROC lower 95% bound - 0.5) / 0.5, 0,
1) of the condition x measurement's performance record (`measurement_performance`), and `expected_yield` = burden
percentile x sensitivity at 0.90 specificity x implementation reach (share of the population within 50 km of an
implementer-group facility). Weights: burden, vulnerability, desert 0.15 each; clinic capacity, research readiness
0.10 each; measurement_evidence 0.15; expected_yield 0.20. A measurement whose performance is UNKNOWN or partial is not
ranked under it. Expected detectable cases and false positives are reported per region (a count only where a defensible
count exists: long COVID at state level; otherwise an index and an upper bound on false positives); they are planning
estimates for a candidate pilot, not predictions of diagnoses (`docs/SCORING.md` §9, `results/SCORING_RESULTS.md` §12).

**Rules.** Every component is kept raw and normalised in `deployment_opportunities` (204 columns since the 2026-09-25 metric link; 140 before). Level-D burden is
excluded and the other weights renormalised (the uncertainty says "the ranking carries no burden information").
Regions where a member with a defined burden has no value are **incomplete**: not ranked, listed separately (24
counties for ME/CFS-containing sets, e.g. the 9 Connecticut planning regions). Counties below 10,000 residents are
scored but not ranked.

**Uncertainty (Monte Carlo, 1,000 seeded draws).** Burden error from published CIs (HPS) or a binomial approximation
(CMS), times an evidence-level multiplier (A 1.0, B 1.25, C 2.0); an inherited state value's error x2 plus a per-county
deviation tau = 1.242 percentage points (the method-of-moments between-state SD of the HPS estimates, net of their
sampling error; the raw SD of the 51 state values is 1.78; `docs/SCORING.md`); weights ~ Dirichlet(20 x default).
Output per region: 5th-95th percentile rank interval and P(top 10). Sensitivity analyses S1 (tau = 0), S1b (inherited
multiplier 1 vs 2), S2 (PLACES proxy for long-COVID county burden), S4 (clinic capacity without primary care), S5
(measurement swap), S6 (specialist-only provider density).

---

## 13. Validation: SPEC Tests 1-7 (+ Test 8, temporal holdout)

| test | question (SPEC) | outcome | files |
|---|---|---|---|
| 1 Wearable signal | Does a wearable phenotype distinguish or stratify a relevant condition, compared with a demographic baseline? | **Mixed.** NHANES: added to demographics, +0.019 to +0.032 AUROC for 4 of 5 symptom/function targets and the mortality control (CIs > 0); the ME/CFS-like proxy is null (+0.0139, -0.0126 to +0.0387; wearable-only AUPRC 0.024 at prevalence 0.017); wearable-only AUROC below demographics for 6/6 (CI below 0 for 5/6). Stanford Long COVID: **null** (wearable AUROC 0.504 vs demographics 0.702). **Added (§7.5, §14):** NHANES 2003-2006 hip counts replicate Test 1 for 4 of 6 targets; MUSCLE-ME hip-accelerometer steps separate diagnosed Long COVID + ME/CFS from healthy controls (AUROC 0.831, 0.725-0.918) but not under the worst case for the 8 missing controls (0.612, 0.459-0.763), and the sex-only baseline is near chance (CV AUROC 0.509), so the increment over it (+0.496 pre-specified, +0.330 post-hoc per-repeat) is almost all steps; fibromyalgia thermography, a clinic camera rather than a wearable: device-only AUROC 0.630 (0.553-0.708) vs age + BMI 0.582 (0.497-0.663) | `nhanes_increments.csv`, `nhanes_model_performance.csv`, `stanford_uwakwe_model_performance.csv`, `nhanes0306_replication.csv`, `charlton_lc_mecfs_cpet_source_primary.csv`, `charlton_lc_mecfs_cpet_source_cv_delta_auroc.csv`, `fm_thermography_model_performance.csv`; Fig. 4 |
| 2 Multimodal increment | Does wearable physiology add over clinical variables? | **Small or null.** NHANES: robust only for functional limitation (+0.0069, +0.0040 to +0.0097); fatigue and depression meet the pre-specified rule but are not robust; fair/poor health, ME/CFS-like proxy and mortality null. Stanford: null (-0.024, -0.076 to 0.027). **Added (§7.5, §14):** over 58 NHANES lab layers, accelerometry still adds for functional limitation over 6 of 7 lab bases (e.g. +0.0054, +0.0024 to +0.0082, over all full-sample layers), while lab blocks add over clinical data in 1 of 47 rows; NHANES 2003-2006 replicates Test 2 for 4 of 6 targets (functional limitation +0.0041, fragile); thermography adds nothing demonstrable to age + BMI (+0.053, -0.028 to 0.129) | same; `nhanes_layers_headline.csv`, `nhanes0306_replication.csv`, `fm_thermography_delta_auroc.csv`; Fig. 4 |
| 3 Molecular coherence | Do disease-linked molecular resources return plausible systems, with source counts and agreement/disagreement? | **Partly.** Gene-level agreement above the cross-condition baseline for 5 of 6 source pairs, pathway level 4 of 6; curated flags not independently recovered (3 of 28 with non-literature, unflagged support); mapMECFS list not best-matching in 4 of 4 comparisons; demo-cluster autonomic links drug-target-driven. **Added (§14.4), person-level omics:** mapMECFS person-linked layers 0 of 5 pass; GEO cohorts 0 of 4 primary transfers replicate, GSE270045 source-confounded, ME/CFS cfRNA no better than its covariates. Test 3 at person level is therefore null in the public data processed | `test3_*.csv`, `linked_omics_model_performance.csv`, `geo_cohort_transfer.csv`, `geo_cohort_confounds.csv`; Fig. 5 |
| 4 Geography | Are hotspots driven by burden rather than population size? (prevalence and counts) | **Partly.** Prevalence series are not large-population-driven (Long COVID state rho -0.26, CI -0.53 to 0.07, which includes 0) while counts are (0.95); but the county diagnostic desert used in ranking is small-population-driven: its top-25 lists hold 9-19 bottom-quartile counties vs 6.3 expected (Holm-significant for 3 of 6 series; long COVID 19/25). The "10 burden-driven, 2 mixed" desert verdicts come from a rule that tests only the large-population tail | `test4_*.csv` |
| 5 Deployment | Does adding facility/research information change where the engine would deploy? | **Yes, partly by construction.** Cleanest contrast, T5-2 (burden + vulnerability, no facility or research information, vs full): top-10 overlap 3/10, top-25 13/25 (Jaccard 0.35), Kendall tau-b 0.736. Literal SPEC contrast, T5-1 (burden only vs full): top-10 0/10 (tie-averaged 0.29; 7 counties tie at the burden-only top-10 boundary), top-25 4/25 (Jaccard 0.087), tau-b 0.467 over 2,395 counties. Only the research-readiness part tracks real facility locations (Test 6) | `test5_*.csv`; Fig. 8 |
| 6 Negative controls | Shuffle phenotype labels, geography labels, clinic locations: does structure disappear? | **Phenotype labels:** yes (NHANES null means 0.496-0.500, p = 0.005; Uwakwe wearable-only already inside its null). **Geography labels:** the burden and SVI shuffles only confirm that the composite uses those components (the composite-burden Spearman falls from 0.657 to 0.098 because the composite is a weighted sum that contains the shuffled component: near-tautological); the substantive geography control is negative: no cross-source relationship survives spatial (CRH) correction, only the same-source positive control. **Clinic locations: largely no loss at composite level.** With every facility and provider location shuffled nationally the deployment ranking keeps Spearman 0.94 (0.936-0.945) with the observed ranking and 8.2 of the top-25 (1.7 of the top-10); clinic_capacity keeps 0.965 with its observed value (it measures facility density), research readiness falls to 0.465. **Added (§14), label permutations of the new analyses:** null means sit at chance throughout (MUSCLE-ME steps 0.501, p 0.0001; linked omics 0.47-0.51, none of 5 layers p < 0.01; thermography device-only 0.496, p 0.0050; hEDS/HSD Olink 0.499, p <= 0.001; NHANES 2003-2006 wearable-only 0.488-0.501); the hEDS/HSD site negative control (American vs Italian patients) is null (0.514, p 0.40); the NHANES lab positive-control gate passes | `nhanes_permutation_null.csv`, `stanford_uwakwe_permutation_nulls.csv`, `test6_geography_*.csv`, `test6_scoring_*.csv`, `test6_clinic_shuffle.csv`, `charlton_lc_mecfs_cpet_source_primary.csv`, `linked_omics_permutation_null.csv`, `fm_thermography_permutation_null.csv`, `heds_hsd_olink_serum_cinquina2026_permutation_null.csv`, `nhanes0306_permutation_null.csv` |
| 7 Ranking robustness | Re-rank under several weight sets; report unstable regions | **Order concordant, short list unstable.** Kendall's W 0.806 over 8 named sets (0.715 over 1,000 random weight vectors); 22 of the default top-25 meet the instability rule; median top-25 Monte Carlo interval width 319 ranks; 0 of the top-10 in the top-10 in >= 50% of draws | `test7_*.csv`, `test7_unstable_regions.csv`; Fig. 8 |
| 8 Temporal holdout (added 2026-09-25; not a SPEC test) | Would rankings built only from data dated <= 2021 have pointed to where Long COVID / ME/CFS research activity appeared in 2022-2026? (proxy outcome: research readiness, not clinical value) | **Null against population; readiness side partly supported; gap persistence holds.** Equal composite AUROC 0.563 (0.535-0.591) for any new relevant trial within 50 km, but population within 50 km alone 0.881 (0.863-0.897): Delta -0.318, and population-adjusted OR per SD 1.11 (0.96-1.27). evidence_weighted 0.493 (chance). Research readiness 0.819 (0.796-0.842), population-adjusted OR 2.30 (2.01-2.66), beats prior activity alone (+0.031, 0.015-0.047). Desert top decile: later-activity rate 0.062 vs 0.193 (lift 0.32, 0.18-0.51) | `TEMPORAL_HOLDOUT_RESULTS.md`, `temporal_holdout_*.csv`, `docs/ANALYSIS_PLAN_TEMPORAL_HOLDOUT.md` |

**Test 5 detail** (`test5_contrasts.csv`; `deployment_opportunities.parquet` `rank_burden_only`,
`composite_burden_only`): the burden-only top-10 is ten Mississippi counties with no condition trial or NIH project in
reach (research-readiness percentile at the 0.40 floor). Burden-only ranks 9-15 are a 7-way tie (composite 0.99551:
Clarke, Hancock, Lamar, Lauderdale, Newton, Perry and Walthall counties, MS; integer CMS percentages plus one inherited
state value per state), so which counties make the burden-only top-10 depends on tie-breaking: one of the tied
counties, Lauderdale County MS (burden-only rank 12), is rank 10 of the full ranking, and the tie-averaged top-10
overlap is 0.29 rather than 0. Burden + vulnerability alone (T5-2) keeps 3 of the full top-10 and 13 of the top-25;
burden + vulnerability + desert (T5-3) keeps 0 and 8 (Jaccard 0.190). Over the 16 county combinations that have a
burden-only ranking, the top-25 Jaccard with burden-only is 0.042-0.087 (the 8 level-D combinations, POTS and
dysautonomia x 4 bundles, have none; their 0 in `SCORING_RESULTS.md` is an empty list). Adding components that carry 40-80%
of the weight moves a short list by construction; the clinic-location shuffle (Test 6) shows that only the
research-readiness part tracks where real facilities are, and that random facility placement moves the top-25 about
as much as T5-3 (8.2 vs 8 of 25 kept).

**Test 6 detail.** Geography-label shuffles of burden and SVI (200 draws each, national and within-state) remove the
shuffled component's structure in every case but one (level-D conditions have no burden to shuffle). This part of the
control is close to tautological: it measures Spearman(composite, shuffled component), and the composite is a weighted
sum that still contains the shuffled component, so a fall is expected by arithmetic; it confirms the composite uses
the component, not that the structure is real. The informative results are two: for Long COVID the within-state
burden shuffle changes nothing, which shows that the inherited county value carries no within-state information
(`test6_scoring_label_shuffle.csv`); and across independent sources (claims vs BRFSS model vs HPS), no county or state
relationship survives the Clifford-Richardson-Hemon effective-sample-size correction plus Holm, while in simulation the
naive exchangeable shuffle rejects 48-72% of independent smooth-field pairs at nominal 5%
(`test6_geography_crh_calibration.csv`): the substantive geography negative control is a null. The clinic-location
shuffle at composite level (`test6_scoring_clinic_shuffle.csv`, 200 national draws): Spearman with the observed ranking
0.94 (0.936-0.945), top-25 overlap 8.2 (5.0-12.0), top-10 overlap 1.7; within-state draws give 0.94 and 8.8. 10 of 28
restated phenotype-level controls are null, including the Uwakwe wearable-only models
(`test6_scoring_phenotype_label_nulls.csv`).

**Test 7 detail** (`test7_summary.csv`, `test7_monte_carlo_decomposition.csv`). The instability comes from both
burden uncertainty and weights: for the primary query the median top-25 interval width is 319 with both, 233 with
burden noise only and 121 with weight noise only. Pairwise Kendall tau-b between named sets ranges 0.315
(capacity_led vs burden_only) to 0.904 (equal vs saturation_adjusted). S6 specialist-only (primary care removed from
provider density) keeps 5 of the top-10 and 17 of the top-25 (primary care is 86.1-91.9% of relevant NPIs;
`results/SCORING_RESULTS.md` §8, from `test7_sensitivity_specialist_only_qa.csv`). S2 (PLACES
proxy for long-COVID county burden) keeps 3 of the top-10.

**Test 8 — temporal holdout** (added 2026-09-25; plan `docs/ANALYSIS_PLAN_TEMPORAL_HOLDOUT.md`, written before any
holdout outcome was computed; results `results/TEMPORAL_HOLDOUT_RESULTS.md`, `results/tables/temporal_holdout_*.csv`,
table `temporal_holdout_rankings`). The unchanged scoring engine was run on inputs dated on or before 2021-12-31
(trials by first-posted date, NIH by fiscal year, NPIs by enumeration date, facilities by earliest dated source, CMS
MMD condition 51 of 2021; Long COVID has no burden measure before 2022, so the clean variant ranks on the ME/CFS
proxy); with no cutoff the same code reproduces the stored composites exactly (max abs diff 0). Outcome (primary): any
Long COVID or ME/CFS trial first posted after 2021 with a U.S. site in the county or within 50 km; 462 of 2,395 eligible
counties (0.193). This is a **proxy outcome** for research readiness; it says nothing about clinical value, yield or
patient benefit, which no public outcome measures.

| score (T = 2021, county, O1) | AUROC (95% CI) | Delta vs population within 50 km | population-adjusted OR per SD |
|---|---|---|---|
| equal composite (default) | 0.563 (0.535-0.591) | -0.318 (-0.347 to -0.288) | 1.11 (0.96-1.27) |
| evidence_weighted composite | 0.493 (0.467-0.523) | -0.388 (-0.417 to -0.355) | 0.94 (0.82-1.07) |
| study_partner composite | 0.744 (0.718-0.770) | -0.137 (-0.161 to -0.113) | 1.67 (1.45-1.93) |
| research_readiness | 0.819 (0.796-0.842) | -0.062 (-0.083 to -0.039) | 2.30 (2.01-2.66) |
| clinic_capacity | 0.677 (0.651-0.702) | -0.204 (-0.232 to -0.176) | 1.58 (1.37-1.85) |
| diagnostic_desert (expected negative) | 0.301 (0.273-0.326) | -0.580 | 0.60 (0.51-0.70) |
| burden (ME/CFS CMS proxy) | 0.484 (0.455-0.510) | -0.397 | 0.89 (0.78-1.02) |
| prior relevant trials within 50 km (baseline) | 0.787 (0.765-0.811) | -0.093 | |
| population within 50 km (baseline) | 0.881 (0.863-0.897) | | |

Reading (pre-specified rules): **no engine score beats population on AUROC**, so every composite is "explained by
population", and the evidence_weighted composite does not rank later activity above chance (1,000 random rankings:
0.471-0.530). The readiness side carries information beyond population and beyond persistence (research readiness
beats the prior-activity baseline by +0.031, 0.015-0.047; its OR stays 1.45, 1.18-1.77, with prior activity in the
model; study_partner is borderline there, 1.17, 1.00-1.36), but the equal composite dilutes it with components that
point the other way. The desert side shows **gap persistence**: the desert's top decile had a later-activity rate of
0.062 vs 0.193 (lift 0.32, 0.18-0.51). Robustness: T = 2019 / 2020 give the same pattern (equal 0.532 / 0.533;
research readiness 0.733 / 0.748); the new-facility outcome (sites without a prior relevant trial) gives equal 0.549
and research readiness 0.819; not date-filtering providers and facilities moves every AUROC by <= 0.012 (inside the
CI half-width), so that leakage path is small (SVI / ACS vintages cannot be bounded). The not-clean variant with the
post-period Long COVID burden gives equal 0.528. Leave-one-state-out normalisation: median within-state Spearman 0.946
(T = 2021) and 0.928 (current ranking; Delaware, Maine, North Dakota < 0.8); rolling T = 2021 vs current ranking:
Spearman 0.890, top-25 overlap 12 of 25.

---

## 14. Beyond wearables: labs, omics and other devices

The wearable results of §7 leave the scientific half of the thesis unshown. This section asks the same question of
every other open, person-level measurement found: do labs, omics or non-wrist devices separate an invisible-illness
label from controls in public data, under a test fixed in advance? It also curates what the literature claims
(§14.7) and where the controlled cohorts that could answer it better are (§14.8). No result here is linked to a
wearable participant of §7, to a geography or to a facility.

### 14.1 Why and how

* **Pre-specified plans, locked before any group comparison:** `docs/ANALYSIS_PLAN_CHARLTON_LC_MECFS_CPET_SOURCE.md`,
  `docs/ANALYSIS_PLAN_FM_THERMOGRAPHY.md`, `docs/ANALYSIS_PLAN_LINKED_OMICS.md`, `docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md`,
  `docs/ANALYSIS_PLAN_KLEIN2023_MYLC_ML_TABLE.md`, `docs/ANALYSIS_PLAN_ENDO_ARG1_REPOD.md`,
  `docs/ANALYSIS_PLAN_HEDS_HSD_OLINK_SERUM_CINQUINA2026.md`, `docs/ANALYSIS_PLAN_NHANES_LAYERS.md`,
  `docs/ANALYSIS_PLAN_NHANES_2003_2006.md`. The GEO-cohort plan is the `analysis_plan` block of `configs/geo_cohorts.yaml`
  (written before any case/control statistic; no separate plan document exists; module docstring of
  `measure_it.omics.geo_cohorts`). The two device primaries were first locked in `docs/DEVICE_DATASET_DISCOVERY.md`, and the three lab
  plans were written before any group statistic was computed (`docs/LAB_DATASET_DISCOVERY.md`); the discovery code
  inspected schemas and label counts only. The ordering is verified from file timestamps; what a person saw
  interactively cannot be verified from files. Every deviation is listed in its results file.
* **Datasets were found, not chosen for their results.** 90 device candidates and 182 lab candidates were verified
  anonymously and ranked by a transparent rubric (label quality, n, relevance, completeness, controls); the top-ranked
  open, criteria-labelled sets were analysed (both lab sets tied at rank 2/3 were analysed rather than break the tie
  post hoc). GEO series were screened by rule (88 series, 12 passing rows, 7 subsets analysed; `geo_cohort_screen.csv`).
* **Person-level re-analyses of public releases.** Each dataset is its own partition with namespaced ids
  (`participants__<source_id>`, `phenotype_signatures__<source_id>`); nothing is pooled across datasets, and omics are
  person-linked only where one file carries a shared id (mapMECFS study numbers, the MY-LC row per person).
* **Published claims are kept apart from our computation** in every results file, and nulls are reported first.

### 14.2 Summary

| dataset (source_id) | label | cases vs controls | measurement | primary pre-specified result (95% CI; p) | verdict | file |
|---|---|---|---|---|---|---|
| MUSCLE-ME, Charlton 2026 (`charlton_lc_mecfs_cpet_source`) | Long COVID (Canadian Consensus Criteria) + diagnosed ME/CFS vs healthy controls | 51 vs 22 (8 of 30 controls lack steps) | hip ActiGraph daily steps, fit-free | AUROC 0.831 (0.725-0.918); permutation p 0.0001; worst case for the missing controls 0.612 (0.459-0.763) | rule met on complete cases; **null** in the locked worst case | `results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md` |
| Klein 2023 MY-LC (`klein2023_mylc_ml_table`) | Long COVID clinic cohort vs healthy + convalescent controls | 98 vs 77 | serum cortisol (released batch-integrated z-scores) | Delta AUROC over draw time + age + sex + BMI +0.253 (0.177-0.328); cortisol alone 0.963 (0.933-0.987), p <= 0.001 | rule met | `results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md` |
| RepOD PY1P9X (`endo_arg1_repod`) | laparoscopy-confirmed endometriosis vs surgical non-endometriosis controls | 120 vs 32 | serum ARG1 ELISA, fit-free | AUROC 0.846 (0.762-0.918); p <= 0.0001 | rule met | `results/ENDO_ARG1_REPOD_RESULTS.md` |
| Cinquina 2026 (`heds_hsd_olink_serum_cinquina2026`) | Italian hEDS + HSD (2017 criteria) vs Italian healthy controls | 88 vs 176 | 458 serum Olink assays, L2 logistic | CV AUROC 0.684 (0.624-0.747); p <= 0.001 | rule met (modest) | `results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md` |
| mapMECFS person-linked omics (`mapmecfs_nih_pi_mecfs`) | PI-ME/CFS vs healthy volunteers | 10-21 vs 11-21 per layer | CSF and serum SomaScan, PBMC and muscle RNA-seq, CSF metabolomics | best: CSF metabolomics 0.825 (0.646-0.962), p 0.0130; bar p < 0.01, Bonferroni over 5 | **null**: 0 of 5 pass | `results/LINKED_OMICS_RESULTS.md` |
| GEO cohorts (`geo_cohorts`) | Long COVID or ME/CFS per study vs its controls | 7 subsets, 252 vs 234 | blood transcriptome, plasma cfRNA, T-cell methylation | cross-cohort transfer AUROC 0.356-0.561 over the 4 primary transfers, none with p < 0.05 | **null**: 0 of 4 replicate; within-cohort results confounded (§14.4) | `geo_cohort_transfer.csv` (no write-up file) |
| NHANES 2011-2014 labs (`nhanes_2011_2014`) | prescription-coded fibromyalgia, migraine, insomnia, myalgia vs general population | e.g. 36 vs 5,229 (fibromyalgia) | 60-analyte primary panel added to demographics | fibromyalgia +0.022 (-0.027 to +0.072); migraine -0.012 (-0.056 to +0.033); insomnia +0.024 (-0.010 to +0.055); myalgia -0.034 (-0.077 to +0.013) | **null**; positive controls pass | `results/LABS_VS_DIAGNOSIS_RESULTS.md` |
| NHANES 2011-2014, 58 lab layers | 6 symptom, function, proxy and mortality targets | e.g. 7,530 (2,007) | lab blocks over clinical; accelerometry over lab blocks | labs: 1 of 47 rows supported (+0.0025, +0.0002 to +0.0046); accelerometry over all full-sample layers, functional limitation +0.0054 (+0.0024 to +0.0082) | **null** for labs; accelerometry still adds for functional limitation | `results/NHANES_LAYERS_RESULTS.md` |
| Appelman 2024 (`appelman_lc_pem_source`) | Long COVID vs healthy | plasma 25 vs 21; muscle 25 vs 19 | 83 plasma / 116 muscle metabolites | plasma 0.694 (0.555-0.819), p 0.063; muscle 0.674 (0.534-0.809), p 0.096 | **null** (both tissues) | `results/LABS_VS_DIAGNOSIS_RESULTS.md` Part 2 |
| Sempere-Rubio 2021 (`fm_thermography`) | ACR 2010 fibromyalgia vs women without it | 85 vs 92 | infrared camera, six regional average skin temperatures | Delta AUROC over age + BMI +0.053 (-0.028 to 0.129); device-only 0.630 (0.553-0.708), p 0.0050 | **primary not met** (reproduces the authors' published null) | `results/FM_THERMOGRAPHY_RESULTS.md` |
| NHANES 2003-2006 (`nhanes_2003_2006`) | the §7 targets | 6,283 adults | hip ActiGraph counts | 8 of 12 Test 1 / Test 2 verdicts replicate | partly replicates §7.1 | §7.5 |

### 14.3 What separates cases from controls, and how far that goes

* **Hip-accelerometer steps and CPET in MUSCLE-ME (Long COVID, ME/CFS).** Daily steps separate 51 patients from 22
  healthy controls (0.831, 0.725-0.918). The group medians match the paper's Table 1 to within one step (7,153 / 4,718 /
  3,704 steps per day for controls / Long COVID / ME/CFS), so the rows are the paper's; the AUROCs are this project's
  own. The verdict does **not** survive the missing controls: with the 8 controls without steps set to the lowest
  observed value, AUROC falls to 0.612 (0.459-0.763), and those 8 are not missing at random (7 of 8 are women vs 8 of
  22 with steps; Fisher p 0.035, post hoc). In the same people, steps (0.824) and relative VO2peak from CPET (0.866)
  do not differ (-0.043, -0.192 to 0.092): a step counter is cheaper and needs no maximal exercise bout, but is not
  more informative. Steps do not tell Long COVID from ME/CFS (0.577, 0.414-0.726); Long COVID alone vs controls is
  "weak" under the locked rule (0.813, 0.675-0.918), ME/CFS alone meets it (0.848, 0.726-0.942). At a Youden threshold
  chosen inside each fold, the steps-only model finds 0.55 of patients at specificity 0.91, against healthy people
  only. Activity limitation is part of both case definitions, so this measures the defining limitation in diagnosed
  people, not undiagnosed illness (`results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md`).
* **Serum cortisol in MY-LC (Long COVID).** The released rows reproduce the published result: fit-free AUROC 0.966
  (0.937-0.988) against the paper's 0.96 (0.93-0.99); cortisol alone 0.963 (0.933-0.987); within draw-time tertiles
  (post hoc) 0.969. The size is itself a caveat: Hedges g -2.78 for one draw of a diurnal hormone, on batch-integrated
  z-scores with no batch or plate column, in a clinic-vs-advertisement recruitment design, so batch or recruitment
  differences can be neither checked nor excluded. **Two other published cohorts did not reproduce the finding:**
  cortisol did not differ in 91 PASC vs 87 controls (PLE-042), and a morning-cortisol clinic series found low values in
  1.2% of 86 (PLE-043); no meta-analysis was found (`docs/PUBLISHED_LAB_EVIDENCE.md`). A first run was invalid (the
  plan's log10 transform of z-scores dropped most cases, 7 vs 68 people) and was redone under a documented deviation;
  the invalid outputs are kept in `results/tables/_invalid_runs/` (`results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md`).
* **Serum ARG1 in endometriosis.** 0.846 (0.762-0.918) against women operated on for other benign conditions, the
  clinically closer comparator, but only 32 of them; within ELISA plates holding both groups 0.854 (0.775-0.928), so
  plate does not explain it; 0.66 sensitivity at 90% specificity. The paper reports 0.848 (0.769-0.926) on 105 vs 22.
  Against healthy controls 0.918 (0.863-0.964), but no ELISA plate is recorded for any healthy control, so that
  contrast cannot be checked for batch. ARG1 does not track rASRM stage (rho -0.00, p 0.991). Women referred with
  pelvic pain, the intended-use population, were not studied (`results/ENDO_ARG1_REPOD_RESULTS.md`).
* **Serum Olink proteome in hEDS/HSD.** 0.684 (0.624-0.747) for Italian patients vs Italian controls, 0.30 (0.20-0.41)
  sensitivity at 90% specificity: a modest group-level separation, not a usable test. The site negative control is
  null (American vs Italian patients 0.514, 0.441-0.586, p 0.40); probands only 0.712 (0.636-0.786); with plate
  indicators 0.681 (0.621-0.744); hEDS vs HSD 0.480 (0.373-0.589), null as published. 24 of 458 assays have q < 0.05,
  all in the paper's list of 69 pooled proteins. Controls are older on average and their age is not released, so part
  of the separation could be age (`results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md`).

### 14.4 What does not separate cases from controls

* **mapMECFS person-linked omics (PI-ME/CFS, 5 layers): 0 of 5 pass** (one-sided permutation p < 0.01, Bonferroni over
  5). CSF SomaScan 0.559 (0.427-0.686), p 0.2837; serum SomaScan 0.697 (0.564-0.815), p 0.0250; PBMC RNA-seq 0.703
  (0.546-0.851), p 0.0490; muscle RNA-seq 0.655 (0.454-0.840), p 0.1099; CSF metabolomics 0.825 (0.646-0.962), p 0.0130.
  Three layers are nominal (p < 0.05), a count that was not pre-specified in layers sharing people. CSF metabolomics is
  within Monte Carlo error of the pass line, but its participants are about 11 years younger and after in-fold age
  removal it falls to 0.705 (p 0.0619). Fusing CSF and serum SomaScan scores below serum alone (-0.049, -0.144 to
  0.029), and the two fluids do not agree on which proteins differ (rho of per-aptamer t 0.036 vs a joint-permutation
  null mean of 0.119). Agreement with the paper's directions is high but is a reproduction check on largely the same
  people. No wearable, HRV, tilt or CPET measurement is linkable in open data, so none of this is device evidence
  (`results/LINKED_OMICS_RESULTS.md`).
* **GEO case/control cohorts (Long COVID, ME/CFS): 0 of 4 primary transfers replicate.** No
  `results/GEO_COHORT_RESULTS.md` write-up exists; the numbers are from the tables. Within cohorts (`geo_cohort_classifier.csv`,
  1,000 permutations): GSE226260_wb 0.656 (0.433-0.804), p 0.088; GSE270045 0.998 (1.000-1.000), p 0.001 (Holm 0.007);
  GSE275334_lc (NanoString) 0.686 (0.511-0.886), p 0.039 (Holm 0.234); GSE293840 (ME/CFS cfRNA) 0.579 (0.503-0.671),
  p 0.107; GSE156792 (T-cell methylation) 0.574 (0.496-0.712), p 0.141; GSE16059 (discordant twins) 0.523
  (0.414-0.675), p 0.263; GSE227375 (replication only) 0.489 (0.315-0.688), p 0.514. Transfers between compatible
  primary cohorts (`geo_cohort_transfer.csv`): GSE226260_wb -> GSE270045 0.356 (0.180-0.559), GSE270045 -> GSE226260_wb
  0.409 (0.227-0.604), GSE16059 -> GSE227375 0.561 (0.367-0.743), GSE227375 -> GSE16059 0.467 (0.346-0.588); none
  replicated. Of 4 secondary transfers on NanoString panel genes, the one marked replicated (GSE275334_lc -> GSE270045,
  0.786, 0.612-0.921, p 0.004) targets the confounded cohort; the post-hoc panel-gene transfers between the two
  whole-blood cohorts do not replicate (0.663, p 0.048; 0.606, p 0.143; `geo_cohort_transfer_posthoc.csv`).
  * **GSE270045 has a source confound:** all 19 cases carry subject-id prefix `S`, the 17 controls `CCI` (13) and `HP`
    (4) (post-hoc chi-square p 1.5e-08; `geo_cohort_confounds.csv`), so collection source or batch cannot be separated
    from the label; 5,724 of 16,282 genes have q < 0.05 with genomic inflation 6.47 (`geo_cohort_de_summary.csv`).
    Its 50 signature rows in `phenotype_signatures__geo_cohorts` (the only cohort with rows) all carry
    `confound_flagged` and are labelled dataset-specific. Its "Long COVID" cases meet ME/CFS criteria
    (`data/raw/geo_cohorts/DATA_AUDIT.md`).
  * **The ME/CFS cfRNA signal is explained by covariates:** sex, batch and test site alone give 0.646 (0.544-0.711),
    and adding the cfRNA gives 0.619 (0.553-0.720), a change of -0.026 (-0.077 to 0.095); test site is associated with
    the label (p 0.0036), and no gene reaches q < 0.05 (genomic inflation 5.88).
  * GSE275334_lc's increment over its covariates (+0.307) is inflated by a covariate-only model below chance (0.367).
* **Routine NHANES labs vs umbrella labels: positive controls pass, umbrella labels null.** Labs add to demographics
  where they should (diabetes without glycaemic markers +0.147, +0.136 to +0.159; weak/failing kidneys +0.119; anaemia
  +0.107; liver +0.091; heart failure +0.056; gout +0.048). They add nothing detectable for the prescription-coded
  labels (§14.2), and the pairwise contrasts between them are null (e.g. fibromyalgia vs migraine +0.057, -0.050 to
  +0.156). The ME/CFS-like proxy's locked "labs add" (+0.046, +0.011 to +0.080) is carried by tobacco exposure and a
  lower HbA1c; post hoc it shrinks to +0.027 (-0.005 to +0.057) with cotinine and BMI in the baseline and +0.005
  (-0.020 to +0.032) against controls without the definition's exclusion diagnoses. Symptom-label increments are mostly
  BMI and smoking (post hoc: fatigue +0.009, depression +0.008, told sleep disorder -0.002). Labs-only models beat their
  permutation null for 24 of 25 labels, but for 16 they are within 0.03 of demographics: the panel re-learns age and sex
  (`results/LABS_VS_DIAGNOSIS_RESULTS.md` Part 1).
* **58 NHANES lab layers do not add robustly; accelerometry still adds for functional limitation.** Lab blocks over the
  clinical model: 1 of 47 estimable rows supported (functional limitation +0.0025, +0.0002 to +0.0046, Holm p 0.080),
  11 worse, 35 null; per-layer screen: 0 of 292 supported (251 null, 33 worse, 8 not robust). Accelerometry over the lab
  blocks: 9 of 47 rows supported, 6 of them functional limitation (e.g. +0.0054, +0.0024 to +0.0082, over all
  full-sample layers; +0.0069 over clinical data alone in §7.1). Every estimate is conditional on having the layer, and
  the subsample populations differ (`results/NHANES_LAYERS_RESULTS.md`).
* **Appelman 2024 metabolomics (Long COVID vs healthy): null.** Plasma and muscle primaries do not pass (above); the
  pre-specified secondary plasma + muscle model reaches 0.762 (0.636-0.874), p 0.012, one of several secondaries at
  25 vs 19. Per metabolite, q < 0.05 for 0 of 83 plasma and 1 of 116 muscle metabolites. VO2max alone gives 0.762
  (0.612-0.897) and adding plasma metabolites changes it by -0.021 (-0.187 to +0.142). No ME/CFS group exists in the file.
* **Fibromyalgia thermography: fails its primary.** Six regional skin temperatures add +0.053 (-0.028 to 0.129) to age
  + BMI; the device alone beats label shuffles (0.630, p 0.0050) but finds 25% (16-36%) of cases at 90% specificity.
  The difference runs opposite to the authors' autonomic hypothesis (warmer knee and lower back in fibromyalgia; knee
  mean +0.59 degC, BH q 0.010), as the paper reports. The one pre-specified sensitivity analysis with all 18
  temperatures clears 0 (+0.100, 0.011-0.185), weak evidence where the primary does not
  (`results/FM_THERMOGRAPHY_RESULTS.md`).
* **NHANES 2003-2006:** §7.5 (8 of 12 verdicts replicate; the ME/CFS-like proxy result with 36 cases is a hypothesis).

### 14.5 Caveats that apply to every result above

* **Re-analysis is not replication.** Every positive result re-analyses the authors' own released rows: MY-LC 0.966 vs
  the published 0.96, ARG1 0.846 vs 0.848, thermography 18 of 18 per-region significance calls, MUSCLE-ME Table 1
  medians, mapMECFS sign agreement. Appelman 2024 and MUSCLE-ME are the same trial (NCT05225688) and probably the
  same people. No result here was tested in an independent cohort.
* **Spectrum bias.** Controls are healthy or convalescent volunteers (MUSCLE-ME, MY-LC, Cinquina, thermography,
  mapMECFS), healthy, recovered or unaffected co-twin controls (GEO), surgical patients (ARG1) or the general population
  (NHANES), never people with the clinical look-alikes
  (deconditioning, depression, other chronic pain or fatigue). AUROCs against such controls overstate diagnostic use:
  in the one lab literature example that measured it, Lyme serology specificity falls from about 95% against healthy
  controls to about 80% in cross-sectional studies (PLE-091).
* **Incorporation bias.** Where the measurement is part of the label, accuracy is agreement with the definition: activity
  limitation in the Long COVID and ME/CFS case definitions (MUSCLE-ME steps); the tryptase 20% + 2 ng/mL rule and urinary
  mediators such as LTE4 that define MCAS (PLE-001, PLE-002); the orthostatic heart-rate rise that defines POTS (PDE-001);
  upright norepinephrine that defines hyperadrenergic POTS (PLE-025); gastric emptying that defines gastroparesis.
* **Small n, one site, recruitment channels.** 21-42 people per mapMECFS layer, 32 surgical controls for ARG1, 22
  controls with steps in MUSCLE-ME; MY-LC clinic patients vs advertised volunteers; thermography association members vs
  advertised controls; Cinquina's American patients were home draws while every control is Italian.

### 14.6 No open MCAS or POTS lab dataset

* **MCAS:** no open per-person dataset has urinary LTE4, 2,3-dinor-11β-PGF2α, tetranor-PGDM or N-methylhistamine with
  MCAS labels and controls. The only per-person MCAS lab file (Vysniauskaite 2015: 238 MCAS, 19 systemic mastocytosis)
  has no controls, and tryptase and KIT D816V appear only in case-only mastocytosis tables of 22-115 patients
  (`docs/LAB_DATASET_DISCOVERY.md`).
* **POTS:** no open per-person dataset has supine/upright catecholamines, GPCR autoantibodies or blood volume with
  controls; the largest autoantibody study (Hall 2022, 116 vs 81, published null) is available on request, and the only
  open catecholamine table is 30 paediatric orthostatic-intolerance patients in a PDF, without healthy controls. No
  public dataset processed has posture-linked heart rate (§7.4).
* **In All of Us these tests are rare** (people who had the test, not diagnoses; `docs/DATA_ACCESS_PLAN.md`,
  `controlled_cohort_counts.csv`): serum tryptase 4,420; urinary LTE4/creatinine 60; 2,3-dinor-11β-PGF2α (24-h) 60;
  N-methylhistamine (24-h) 280; plasma norepinephrine 1,080, of which <= 20 coded as standing or supine; no concept for
  GPCR autoantibodies. Any MCAS-coded (D89.4: 660 people) x LTE4 cell would fall under the 20-person disclosure floor.

### 14.7 Published evidence tables (claims, not results)

Two curated tables hold what the literature reports, one row per published claim, each with the verbatim quote it came
from. **Nothing in them was reproduced by this project**; they are targeted curations, not systematic reviews, so a
missing study is not evidence that none exists.

* **Devices** (`docs/PUBLISHED_DEVICE_EVIDENCE.md`; `published_device_evidence`, object ids
  `published_evidence:PDE-###`): 66 claims from 61 publications. The build re-fetches each abstract or open full text
  from Europe PMC: 66 of 66 quotes verified verbatim, 66 of 66 rows with their AUROC, sensitivity, specificity or stated
  counts found inside the quote, 66 of 66 DOIs matching. None of the 17 diagnostic-accuracy rows reports external
  validation (13 primary studies; 4 pooled reviews, where it does not apply). The best-supported device use is orthostatic
  heart rate with posture for POTS, including POTS inside Long COVID (objective orthostatic tachycardia in 31% of 467,
  PDE-045, and 20% of 221, PDE-046, with symptoms not picking out who has it), which is principled and partly circular
  because POTS is defined by that measurement. Endometriosis ultrasound has accuracy against surgery (deep disease 0.79 /
  0.94, PDE-059). ME/CFS 2-day CPET is contested (PDE-013 vs the null replication PDE-015); the pooled Long COVID HRV
  literature is null (PDE-038); RECOVER Fitbit shows group differences only (PDE-039); the one wearable Long COVID
  classifier (PDE-040) could not be reproduced from its release (§7.2). No MCAS device study was found.
* **Labs** (`docs/PUBLISHED_LAB_EVIDENCE.md`; `published_lab_evidence`, object ids `published_lab_evidence:PLE-###`;
  condition-level layer): 109 claims from 102 publications; all 109 quotes verified, and every AUROC, sensitivity,
  specificity and "as stated" count found inside its quote, on the last build. 46 rows use healthy comparators only, 19
  symptomatic or disease comparators and 23 both; 18 are flagged for incorporation bias. Apart from Lyme serology for
  the infection, no lab test for the target conditions has externally replicated accuracy against symptomatic
  look-alikes; positive discovery results that failed independent replication include POTS GPCR autoantibodies
  (PLE-031/032 vs PLE-033), Long COVID cortisol (PLE-040 vs PLE-042/043), migraine CGRP, IBS anti-CdtB/anti-vinculin and
  endometriosis salivary and serum miRNA; RECOVER found no clinically useful routine lab among 25 (1,880 vs 3,351,
  PLE-047).

### 14.8 Access routes that could answer the question better (counts of people with data or tests; nothing analysed)

Read from public pages on 2026-09-24 (`docs/DATA_ACCESS_PLAN.md`, `controlled_cohort_counts.csv`, discovery documents).
These count people who **have** a data type, a code or a test, not people analysed here, and not people with an abnormal
result. No account was created and no data-use agreement accepted.

| route | what it holds (verified counts) | access |
|---|---|---|
| All of Us (rank 1) | 747,040 participants (CDR 2025q4r5); Fitbit 68,840 (heart rate, steps, sleep; no HRV); EHR codes U09.9 3,200, G93.32 1,280, G90.A 1,300, M79.7 15,440, Q79.6 1,600, D89.4 660; serum cortisol 20,160, tryptase 4,420. Fitbit x diagnosis intersections are not public: the overlaps in the plan (e.g. U09.9 about 295) are labelled estimates | institutional DURA plus training; Registered Tier about two weeks after the DURA; aggregates only, no cells of 1-20 |
| RECOVER-Adult on BioData Catalyst (rank 2) | raw Fitbit streams for 6,569; digital-health sub-study 6,529 enrolled, 4,122 connected a device; study-defined post-COVID ME/CFS 531 of 11,785 infected; study labs and a biospecimen inventory in the same people | dbGaP Data Access Request; 2-6 months |
| mapMECFS | NIH PI-ME/CFS per-person HRV, tilt, actigraphy and CPET; Cornell 2-day CPET 84 vs 71 and wrist ActiGraph 58 vs 41, linked by `cor_id` to metabolomics and proteomics packages; every file returned HTTP 403 anonymously | one free account |
| COPERIA (Synapse) | Polar H10 HRV at rest, 6-minute walk and cold pressor, 63 Long COVID vs 66 post-COVID controls with valid data | one free Synapse login |
| ImmPort | IMPACC post-acute studies (Olink, metabolomics, cytokines), Long COVID studies with clinical labs, Lyme SDY1395; no counts recorded | free account plus data-use terms |

UK Biobank (wrist accelerometry 103,568, all before 2020, so 0 Long COVID by design) is closed to new applications
(`docs/DATA_ACCESS_PLAN.md` B.5). The only open release with a wearable, clinical labs and multi-omics on the same named
people is 22 Glucotypes participants with a CGM and no invisible-illness label (`docs/DATA_ACCESS_PLAN.md` Part A).

### 14.9 What this adds

Public person-level data now show four candidate measurable phenotypes beyond wrist accelerometry (hip steps, serum
cortisol, serum ARG1, a serum proteome), each evidence-supported only as a group difference between diagnosed people and
healthy or surgical controls in the authors' own released rows, and one of them (cortisol) not reproduced elsewhere.
Person-level omics, routine labs against umbrella labels, a small metabolomics set and thermography are null or
confounded. None of it involves autonomic measurement, none tests a measurement against clinical look-alikes, and none is
linked to the wearable participants or to the deployment ranking of §12, which is unchanged by this section.

---

## 15. Example deployment recommendation

Query (the SPEC example, dashboard page 5): *"Find 10 U.S. regions where wearable monitoring of autonomic/activity
abnormalities in Long COVID or ME/CFS would be useful to evaluate, and identify clinics/research sites that could
plausibly participate."* Call: `rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring",
"county", top_n=10, weight_set="equal")`. Full output: `results/example_deployment_recommendation.json` (SPEC schema per
region); also `deployment_candidates`, `demo_top10_regions.csv`, `demo_top10_candidate_sites.csv`,
`demo_top10_research_evidence.csv`.

| rank | county | composite | 5-95% rank interval | P(top 10) | burden level | burden / SVI / desert / capacity / readiness pct | condition trials, NIH core projects in reach |
|---|---|---|---|---|---|---|---|
| 1 | Oklahoma County, OK | 0.861 | 1-306 | 0.49 | C | 0.97 / 0.84 / 0.57 / 0.97 / 0.96 | 4, 1 |
| 2 | Wayne County, MI | 0.838 | 1-408 | 0.37 | C | 0.95 / 0.86 / 0.42 / 0.99 / 0.97 | 9, 4 |
| 3 | Franklin County, AL | 0.834 | 4-247 | 0.23 | C | 0.96 / 0.76 / 0.88 / 0.76 / 0.81 | 1, 0 |
| 4 | George County, MS | 0.830 | 2-288 | 0.30 | C | 0.99 / 0.84 / 0.97 / 0.50 / 0.85 | 2, 0 |
| 5 | Crittenden County, AR | 0.827 | 1-581 | 0.30 | C | 0.62 / 0.92 / 0.83 / 0.87 / 0.90 | 0, 1 |
| 6 | Bossier Parish, LA | 0.826 | 1-202 | 0.42 | C | 0.76 / 0.78 / 0.86 / 0.87 / 0.86 | 0, 2 |
| 7 | Webster Parish, LA | 0.825 | 1-302 | 0.45 | C | 0.69 / 0.88 / 0.76 / 0.93 / 0.86 | 0, 2 |
| 8 | Scotland County, NC | 0.823 | 1-289 | 0.26 | C | 0.92 / 1.00 / 0.85 / 0.95 / 0.40 | 0, 0 |
| 9 | Forrest County, MS | 0.818 | 4-323 | 0.21 | C | 0.99 / 0.95 / 0.86 / 0.89 / 0.40 | 0, 0 |
| 10 | Lauderdale County, MS | 0.815 | 5-315 | 0.18 | C | 1.00 / 0.96 / 0.84 / 0.89 / 0.40 | 0, 0 |

Rank 11, Muskogee County OK (composite 0.81491), trails rank 10 by 0.0001, so the two are tied at the three decimals
shown here and on the dashboard. Every burden level is C (the ME/CFS member is a claims proxy) and every long-COVID
member value is the state estimate inherited by the county. 6 of the 10 regions have no Long COVID or ME/CFS trial at a facility within reach, and their
candidate lists carry the matcher's absence note. 11 population-eligible regions are **not ranked** because their
burden is incomplete (9 Connecticut planning regions, North Slope Borough AK, Manassas Park VA).

Rank 1, abbreviated (values copied from the JSON; `...` marks cuts):

```jsonc
{"condition": {"condition_id": "long_covid_or_me_cfs", "kind": "condition_set", "members": ["long_covid", "me_cfs"]},
 "phenotype": [{"condition_id": "long_covid", "datasets": [
     {"dataset_id": "stanford_longcovid_uwakwe2025", "n_cases": 31, "n_controls": 95, "n_features_fdr_significant": 0,
      "null_results": "14 of 14 tested features did not pass BH-FDR q < 0.05"},
     {"dataset_id": "stanford_covid_mishra2020_alavi2022", "phenotype_label": "post-infection recovery trajectory (acute COVID-19 cohort; not a Long COVID label)", ...}]},
   {"condition_id": "me_cfs", "datasets": [{"dataset_id": "nhanes", "phenotype_id": "mecfs_like_proxy", "is_proxy": true,
      "n_cases": 132, "n_controls": 7664, "n_features_fdr_significant": 14, ...}]}],
 "measurement": {"measurement_id": "wearable_autonomic_activity_monitoring",
   "members": ["accelerometry", "wearable_heart_rate", "hrv", "ecg_ambulatory", "ppg", "posture_detection"],
   "adapter": "wearable", "adapter_status": "implemented", ...},
 "technology": {"regulatory_context": [{"measurement_id": "ecg_ambulatory", "n_mapped_product_codes": 10, "n_510k_total": 1649, ...}, ...],
   "guardrail": "Trials/grants show what researchers deploy, not that a measurement works; FDA records are deployment-readiness signals, not evidence that a device diagnoses the target illness."},
 "geography": {"name": "Oklahoma County, Oklahoma", "fips": "40109", "population": 806199,
   "burden": {"percentile": 0.9713, "members": {
       "long_covid": {"burden_value": 7.5, "burden_evidence_level": "A", "burden_inherited": true, "burden_source_resolution": "state", "burden_ci": [5.3, 10.2]},
       "me_cfs": {"burden_value": 28.0, "burden_evidence_level": "C", "burden_source_resolution": "county"}}},
   "burden_evidence_level": "C", "vulnerability": {"svi_overall": 0.8431, "percentile": 0.8432},
   "diagnostic_desert": {"index": 0.5243, "percentile": 0.5687}, "clinic_capacity": {"percentile": 0.9653},
   "research_readiness": {"percentile": 0.9585, "condition_trials_in_pool": 4, "condition_nih_core_projects_in_pool": 1, "technology_experience_trials_in_pool": 3},
   "composite": 0.8614, "rank": 1, "rank_interval_5_95": [1, 306], "p_top10_monte_carlo": 0.492,
   "ranks_under_weight_sets": {"equal": 1, "burden_led": 1, "equity_led": 56, "capacity_led": 1, "study_partner": 1,
                               "access_gap": 261, "burden_only": 72, "saturation_adjusted": 67}},
 "candidate_sites": [   // 20 facility records (about 15 distinct organisations), round-robin over six characteristics; 3 shown
   {"object_id": "facility:npi-1538213764", "facility_name": "MCBRIDE CLINIC ORTHOPEDIC HOSPITAL, LLC", "distance_km": 10.13,
    "selected_via": {"long_covid": "clinical_specialty_match", "me_cfs": "clinical_specialty_match"},
    "framing": "This facility has characteristics suggesting it may be a viable implementation or study partner."},
   {"object_id": "facility:ctgov-74f3a639213f", "facility_name": "University of Oklahoma Health Science Center - Oklahoma Clinical and Translational Science Institute - Appendix A & B",
    "selected_via": {"long_covid": "relevant_trial_history"}, ...},
   {"object_id": "facility:hrsa-BPS-H80-032732", "facility_name": "CITY RESCUE MISSION", "selected_via": {"long_covid": "community_access", "me_cfs": "community_access"}, ...}, ...],
 "research_evidence": {"condition_trials": {"n": 4, "object_ids": ["trial:NCT05172024", "trial:NCT05524532", "trial:NCT05595369", "trial:NCT06305780"]},
                       "nih_projects": {"n_core_projects": 1, "object_ids": ["nih:11451891"]}, ...},
 "molecular_context": [{"condition_id": "long_covid", "label": "condition-level molecular enrichment (public databases; not participant-linked; not patient multi-omics; not mechanism)", ...}, ...],
 "uncertainties": [
   "long_covid: burden is the state HPS estimate for OK inherited by every county of the state (source resolution: state). It carries no within-state variation and is not county prevalence; the Monte Carlo inflates the state estimate's uncertainty x2 (inherited-burden multiplier) and adds a per-county deviation (tau = 1.24 percentage points, assumed equal to the between-state SD).",
   "long_covid: the diagnostic desert's burden term is the inherited state value.",
   "me_cfs: burden is a level-C proxy (CMS MMD Fibromyalgia, Chronic Pain and Fatigue prevalence, unsmoothed actual, all ages); not the condition's prevalence; uncertainty band inflated x2.0.",
   "Rank interval (5th-95th percentile over 1000 draws of burden and weights): 1-306.",
   "me_cfs: No facility in the pool has relevant_trial_history, NIH_research_activity for this query: a possible measurement/research desert around this geography (absence in these public registries, not proof that no such activity exists)."],
 "provenance": {"object_ids": ["opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|40109", "geo:40109", ...]},
 "recommended_next_step": "Candidate deployment opportunity (rank 1 of 2395 under the 'equal' weights; 5th-95th percentile rank interval 1-306 over Monte Carlo draws of burden and weights). A pilot evaluation of wearable autonomic / activity monitoring for Long COVID or ME/CFS in Oklahoma County, Oklahoma could test feasibility and whether the measurement captures a measurable phenotype locally, working with candidate partners such as ... which have characteristics suggesting they may be viable implementation or study partners. This is a candidate deployment opportunity, not a validated diagnostic pathway; check the listed uncertainties first."}
```

Reading: Oklahoma County is rank 1 under four of the eight weight sets but 56 under `equity_led`, 72 under
`burden_only` and 261 under `access_gap`; its burden percentile rests on an inherited state value and a claims proxy.
It is a candidate for a pilot that would *measure* whether the phenotype is visible locally.

The 20 candidate sites are 20 facility records but about 15 distinct organisations: SAINTS MEDICAL GROUP, LLC appears
as three organisation NPIs (two in Oklahoma City, one in Shawnee), and the University of Oklahoma Health Sciences
Center and its Clinical and Translational Science Institute as four records (two NPPES NPIs, two ClinicalTrials.gov
facility keys). This follows from the entity resolution's recall of 0.528 (§11): many records of one organisation are
not merged. The picks are what each separate characteristic selects, not a judgement of fit: an orthopedic hospital is
chosen by clinical_specialty_match because its self-reported NPPES taxonomy falls in a curated specialty group, and
City Rescue Mission (an HRSA service-delivery site of the FQHC Community Health Centers of Oklahoma) by
community_access. Whether any
of them could run a wearable study is not assessed (`results/example_deployment_recommendation.json`,
`demo_top10_candidate_sites.csv`).

---

## 16. Example MCP queries (real outputs, abbreviated)

Captured 2026-09-24 through the in-memory FastMCP client (`docs/MCP.md` §6-7, which has the longer versions). Every tool
returns the same envelope (`tool`, `status`, `query`, `data`, `reason`, `caveats`, `provenance.object_ids`,
`truncation`, `meta`); a tool with no data returns `UNKNOWN / NOT AVAILABLE` with a reason.

```jsonc
// normalize_condition({"condition_or_code": "G90.A"})
{"status": "ok", "data": {"match_status": "matched", "input_type": "icd10cm",
  "matches": [{"canonical_condition_id": "pots", "primary_mondo_id": "MONDO:0001315", "object_id": "condition:pots", ...}],
  "cms_fy2026": {"code": "G90.A", "description": "Postural orthostatic tachycardia syndrome [POTS]", "billable": true}}}

// get_condition_burden({"condition": "POTS", "geography_level": "county"})          level D -> UNKNOWN, never a guess
{"status": "UNKNOWN / NOT AVAILABLE",
 "reason": "No public population measure: no CCW algorithm contains G90.A; PLACES has no orthostatic item.",
 "data": {"condition_id": "pots", "burden_evidence_level": "D", "value": "UNKNOWN / NOT AVAILABLE", "rows": []}}

// get_condition_burden({"condition": "Long COVID", "geography_level": "county", "geo_id": "40109"})
{"status": "ok", "data": {"burden_evidence_level": "A", "primary_source_resolution": "state",
  "rows": [{"geo_name": "Oklahoma County, Oklahoma", "value": 7.5, "ci_low": 5.3, "ci_high": 10.2, "inherited": true,
            "measure_label": "Currently experiencing long COVID, as a percentage of all adults [STATE value inherited as county context; not county prevalence]"}]}}

// trace_evidence({"object_id": "opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|40109"})
{"status": "ok", "data": {"components": {"burden": {"percentile": 0.9713, "evidence_level": "C", "inherited": true}, ...},
  "lineage": [..., {"relation": "burden row behind the burden component (geo_condition_burden)", "source_ids": ["cms_mmd"],
                    "row": {"source_record_id": "mmd:f:2022:v:51:county:40109:unsmoothed_actual:age=all", "value": 28.0, ...}}, ...],
  "sources": [{"source_id": "cms_mmd", "data_audit": [{"path": "data/raw/cms_mmd/DATA_AUDIT.md", "exists": true}],
               "raw": {"record_files": [{"path": "data/raw/cms_mmd/api/f_2022_v__prev_final_long_fltr12_racecat_all_sexcat_all_22_f.json",
                                          "bytes": 6952186, "sha256": "e215a50e7b7292e0322a2cc8f8ea29ba24eb3356a43c747a3d30883b537031f8"}]}}, ...]}}

// trace_evidence({"object_id": "participant:nhanes:62161"})                          person-level ids are refused
{"status": "UNKNOWN / NOT AVAILABLE", "data": {"refused": true,
  "reason": "refused: person-level rows (participants of NHANES, the Stanford wearable studies or MapMECFS) are summarised only. ..."}}

// search_us_open_data({"query": "long COVID"})                                        answered from the HTTP cache
{"status": "ok", "data": {"n_results": 20, "results": [{"title": "Long-term Care and COVID-19",
  "agency": "Centers for Disease Control and Prevention", "retrieve_from": "https://www.cdc.gov/nchs/covid19/npals.htm", ...}]},
 "caveats": ["Data.gov is a metadata catalog: it does not host the records; retrieve_from names the agency URL that does. ...",
             "Results are in Data.gov's own relevance order; this engine does not evaluate catalog relevance ..."]}
```

The top hit shown, "Long-term Care and COVID-19", is not about long COVID: the order is Data.gov's own relevance
ranking, which this engine does not evaluate. In the reproducible offline build the tool can answer only queries whose
responses are already in the HTTP cache (the 54 catalog queries run at ingestion, matched ignoring case and spacing
since this report); any other query returns UNKNOWN / NOT AVAILABLE with the list of cached queries.

**Deterministic agent** (`uv run measure-it ask "..."`; no language model; every sentence cites object ids). SPEC
demo question, abbreviated (`docs/MCP.md` §7). The question names four conditions, so the agent ranks the demo cluster
(Long COVID, ME/CFS, dysautonomia, POTS), not the two-condition set of §15; its top county therefore differs:

```markdown
## 1. What to measure
- Measurable phenotype in public data for Myalgic encephalomyelitis/chronic fatigue syndrome - nhanes: phenotype 'Unexplained fatigue with functional limitation (ME/CFS-LIKE PROXY; not ME/CFS)' (a PROXY label, not the condition itself), 132 cases vs 7,664 non-cases; 14 of 26 wearable features passed BH-FDR; ... [signature:nhanes|mecfs_like_proxy|interdaily_stability]
- Public person-level wearable signature for Postural orthostatic tachycardia syndrome: UNKNOWN / NOT AVAILABLE ... [condition:pots]
- Wearable autonomic / activity monitoring for Long COVID (post-COVID-19 condition): 85 registered trial(s) describe objective use (73 as an outcome measure) and 3 NIH core project(s) mention it; evidence tier high (research activity, not proof that it works). [measurement_evidence:long_covid|wearable_autonomic_activity_monitoring; ...]
- Condition-level molecular enrichment for Myalgic encephalomyelitis/chronic fatigue syndrome (public databases; not patient multi-omics, not mechanism): expression_study_metadata 30 records, genetic_association 38 records, known_drug 43 records, metadata_catalog 8 records; systems with non-literature support after removing post-hoc-flagged sources: none; carried only by post-hoc-flagged sources (mostly drug targets of trialled drugs, or one GWAS locus; not disease biology): autonomic_cardiac, mitochondrial_metabolic, neuronal; literature co-mention only (study attention): connective_tissue, immune, vascular_endothelial. [condition:me_cfs]
## 2. Where to measure it
- Burden of Postural orthostatic tachycardia syndrome at county level: UNKNOWN / NOT AVAILABLE (evidence level D) ... [condition:pots]
  1. Okmulgee County, Oklahoma: composite 0.888, rank interval 1-64 (P(top 10) 0.75); burden percentile 0.97 (evidence level C, includes a state value inherited by the county), ... [opportunity:autonomic_activity_invisible_illness|wearable_autonomic_activity_monitoring|40111; geo:40111]
## 3. Who could realistically deploy it
- For Long COVID (post-COVID-19 condition) x Wearable autonomic / activity monitoring [condition:long_covid], no facility in the pool has relevant_trial_history, NIH_research_activity, technology_experience: a possible measurement/research desert around this geography (absence in these public registries, not proof that no such activity exists). [geo:40111]
```

The molecular sentence is the one §8 describes: after removing post-hoc-flagged drug-target sources no system has
non-literature support. Before the review of this report the agent printed the same six systems as "supported
physiological systems", an overstatement that contradicted §8 (§23 item 9). For an unknown condition ("Zorblaxian
drift syndrome") the agent stops before ranking and reports the condition as UNKNOWN / NOT AVAILABLE. All 271 citations
(each answer's distinct object ids, summed; 164 distinct ids overall) of the six documented agent answers resolve with
`trace_evidence` (re-run after the review fixes; `docs/MCP.md` §7).

---

## 17. Figures and dashboard

Full captions and the code/data each figure is drawn from: [visuals/figures/CAPTIONS.md](../visuals/figures/CAPTIONS.md). PNG (300 dpi) and SVG in
`results/figures/`; interactive versions of Figures 6-9 in `results/maps/`. Every number on Figures 3-9 is read from
the processed or result tables when the figure is built.

**Figure 1. System architecture.** Person -> phenotype -> measurable biology -> technology -> geography -> clinic ->
deployment, grouped by join type: person-level joins on participant ids within one dataset (solid), concept-level
joins on ontology ids (dashed), ecological joins on place (dotted). The MeasurementAdapter is the only
modality-specific step; the MCP/agent layer reads every layer and returns UNKNOWN / NOT AVAILABLE when a tool has no data.

![Figure 1. System architecture](../visuals/figures/fig01_system_architecture.png)

**Figure 2. Public-data source architecture.** The sources of `SOURCE_REGISTRY.yaml` by data layer, with publisher,
unit, geographic levels, access gating and arrows to the canonical tables they write. The figure is drawn from the
registry when it is built: the file in `results/figures/` shows the registry of its last build
(37 sources: 36 ingested, 1 partial, 0 blocked; §3). To stay within the 247 mm page height with 37 sources, bands whose source rows would be taller than their table boxes show one line per source (name and publisher). mapMECFS is partial
(29 of 29 anonymous file probes returned HTTP 403); its open omics enter the molecular layer as condition-level evidence
and, for the omics-to-omics linked subset, as the person-level analysis of §14.4.

![Figure 2. Public-data source architecture](../visuals/figures/fig02_public_data_sources.png)

**Figure 3. Wearable phenotype differences.** (a) NHANES accelerometry, age/sex-adjusted differences (SD units, 95%
CI) for the ME/CFS-like PROXY (132 cases; 14 of 26 features pass FDR), Rx-defined fibromyalgia (34 cases), functional
limitation and the mortality positive control. (b) Stanford acute COVID-19 cohorts (NOT Long COVID): resting HR and
step z-scores around true vs placebo onsets; resting HR peaks at 1.03 SD on day 0 but the pre-specified days 0-14
contrast is null (+0.11 SD, -0.20 to +0.45); steps fall (-1.10 SD). (c) Uwakwe 2025 Long COVID: 0 of 7 adjusted
resting-HR features have p < 0.05; wearable-only AUROC 0.50: a null result.

![Figure 3. Wearable phenotype differences](../visuals/figures/fig03_wearable_phenotype_differences.png)

**Figure 4. Clinical-only versus clinical + wearable model performance.** (a) Cross-validated AUROC with 95% CIs for
five feature sets per NHANES target and three sets in the Uwakwe release. (b) Test 1 and Test 2 increments across model
families and the temporal split: Test 1 supported for 5 of 6 targets (4 of 5 symptom/function targets plus the
mortality control); Test 2 supported and robust only for functional limitation; Uwakwe increment -0.024 (null). AUPRC,
sensitivity, specificity and calibration are in §7.1-7.2, not on the figure. (c) Label-permutation controls.

![Figure 4. Clinical-only vs clinical + wearable](../visuals/figures/fig04_clinical_vs_wearable_models.png)

**Figure 5. Condition -> molecular evidence -> physiological system -> measurement class.** Condition-level molecular
enrichment from public databases, not participant-linked and not patient multi-omics. 9 of 18 drawn condition ->
pathway edges are post-hoc flagged (drug-target-only gene families); curated system -> measurement-class links are
dashed assumptions. (b) Support tier of every condition -> system link, with the tier after removing flagged sources.

![Figure 5. Molecular evidence to measurable physiology](../visuals/figures/fig05_molecular_measurable_physiology.png)

**Figure 6. Long COVID burden and unmet need.** (a, b) HPS "currently experiencing long COVID" by state (level A at
STATE resolution; 2.6% IL to 10.8% SD; most CIs overlap). (c) County diagnostic desert: its burden term is the STATE
estimate inherited by every county, not county prevalence; within a state counties differ only through SVI, provider
density and trial access.

![Figure 6. Long COVID burden and unmet need](../visuals/figures/fig06_long_covid_burden_unmet_need.png)

**Figure 7. Provider and research infrastructure.** (a) Relevant-provider density (long-COVID specialty groups, incl.
primary care); (b) specialists only (13.7% of relevant NPIs; 1,354 of 3,144 counties have none); (c) 17,818 active HRSA
service-delivery sites, 263 Long COVID or ME/CFS trials with a U.S. site (drawn at city centroids), 146 NIH-funded
organisations.

![Figure 7. Provider and research infrastructure overlay](../visuals/figures/fig07_infrastructure_overlay.png)

**Figure 8. Deployment opportunity map** (primary query, county, equal weights). (a) Composite for 2,395 ranked
counties; 725 small counties scored but not ranked; 24 counties with incomplete set burden hatched; top 10 numbered.
(b) Top-10 Monte Carlo rank intervals (widths 201-580; 0 of 10 in the top 10 in at least half the draws). (c) Ranks
under the 8 named weight sets (Kendall's W 0.806).

![Figure 8. Deployment opportunity map](../visuals/figures/fig08_deployment_opportunities.png)

**Figure 9. Geography drill-down: San Diego County, California.** Candidate pool (50 km), HRSA sites, trial sites and
NIH organisations; 18 candidate implementation or study partners with the characteristic that selected each; burden by
condition with evidence level (Long COVID = inherited California value; 8 conditions level D); the county's
components for the primary query (composite 0.54, rank 1,121 of 2,395; Monte Carlo 226-1,589).

![Figure 9. San Diego County drill-down](../visuals/figures/fig09_san_diego_drilldown.png)

**Figure 10. Modularity: wearables today, MeasurementAdapter, CAPRIO capillaroscopy tomorrow.** The two implemented
adapters (their own `describe_measurement()` output and reproduction checks: 38 features for 523 NHANES participants,
max difference 2.9e-11; 15 features for 126 Uwakwe participants, max difference 0), the four-method interface, and the
`CapillaroscopyAdapter` stub whose data methods raise NotImplementedError naming what CAPRIO must supply. The panel's
closing note reads the S5 swap result from `test7_sensitivity_comparisons.csv`: the code path is unchanged, and with
today's public data the candidate deployment opportunities differ only marginally (composite Spearman 0.994, top-10
overlap 9 of 10).

![Figure 10. Adapter modularity](../visuals/figures/fig10_adapter_modularity.png)

**Dashboard** (`uv run measure-it dashboard`; pages described in `docs/DASHBOARD.md` §3; screenshots at default
inputs, captured with headless Chromium):

| page | what it shows at the default input |
|---|---|
| Home | thesis, engine chain, the source registry by layer, the data-layer join rules |
| 1 National Opportunity Map (**partial** against the SPEC) | one choropleth at a time (opportunity, diagnostic desert or burden with evidence level: these three layers are mutually exclusive, not independently toggled as the SPEC asks) plus independently toggled point overlays for HRSA sites, trial sites, NIH organisations and specialist density; the facility overlay shows HRSA sites only, not NPPES or `clinic_registry` facilities; the burden-evidence box is always visible |
| 2 Condition Explorer ("Long COVID") | phenotype features, wearable abnormalities (the Long COVID label shown as a null result), condition-level molecular enrichment, candidate measurements with five separate dimensions, FDA technologies, burden sources |
| 3 Measurement Explorer ("wearable autonomic monitoring") | conditions where it is studied, signals and public datasets, trials, FDA device classes, regions with unmet need |
| 4 Geography Explorer ("San Diego County, California") | burden by condition (inherited value labelled), SVI, provider density, trials, NIH research, candidate technologies and facilities, nearby health centers |
| 5 Deployment Recommendation (SPEC query) | ranked table with component percentiles, rank interval and P(top 10), contribution chart, map, per-region evidence and uncertainties, unranked incomplete-burden regions, JSON download |

![Dashboard: home](../visuals/dashboard/home.png)

![Dashboard: National Opportunity Map](../visuals/dashboard/national_opportunity_map.png)

![Dashboard: Condition Explorer](../visuals/dashboard/condition_explorer.png)

![Dashboard: Measurement Explorer](../visuals/dashboard/measurement_explorer.png)

![Dashboard: Geography Explorer](../visuals/dashboard/geography_explorer.png)

![Dashboard: Deployment Recommendation](../visuals/dashboard/deployment_recommendation.png)

---

## 18. Negative and null results (listed together)

| area | result | file |
|---|---|---|
| Long COVID wearable (Uwakwe 2025) | wearable-only AUROC 0.504 (0.404-0.605), p 0.300; increment -0.024 (-0.076 to 0.027); 0 of 7 features after BH | `stanford_uwakwe_model_performance.csv` |
| Acute RHR | elevation flag 62% true vs 56% placebo (p 0.332); days 0-14 RHR z vs placebo p 0.531; pre-symptomatic window AUROC 0.562 (p 0.115) | `stanford_acute_endpoints.csv`, `stanford_acute_detection.csv` |
| NHANES wearable-only | below demographics for 6/6 targets (CI below 0 for 5/6) | `nhanes_increments.csv` |
| NHANES Test 2 | null for fair/poor health, ME/CFS-like proxy, mortality (AUROC); not robust for fatigue, depression | same |
| ME/CFS-like proxy | Test 1 and Test 2 null; wearable-only AUPRC 0.024 (0.019-0.033) against a base rate of 0.017 (best model 0.035); signature partly definitional; sensitivity variants underpowered | `results/WEARABLE_NHANES_RESULTS.md` §11, `nhanes_model_performance.csv` |
| Rx-defined migraine, IBS | 0 of 26 features pass FDR | `nhanes_phenotype_signatures.csv` |
| Unsupervised | HDBSCAN finds no cluster under 5 settings; GMM partitions moderately reproducible (ARI 0.41) | `nhanes_hdbscan_sensitivity.csv`, `nhanes_cluster_stability_summary.csv` |
| Temporal holdout (Test 8) | rankings built from data dated <= 2021 do not beat population within 50 km at predicting where Long COVID / ME/CFS trials appeared in 2022-2026 (equal composite AUROC 0.563 vs 0.881; population-adjusted OR 1.11, 0.96-1.27); evidence_weighted at chance (0.493); ME/CFS burden proxy unrelated (0.484) | `TEMPORAL_HOLDOUT_RESULTS.md`, `temporal_holdout_metrics.csv` |
| Phase 3 | 18 of 25 phenotype x class cells not observable in public person-level data; no wearable HR verdict is supported; the one sleep_objective verdict rests on an accelerometry sleep proxy, not PSG or sleep staging | `phase3_signal_results.csv` |
| Molecular | GWAS Catalog x OT other agreement not above baseline at gene or pathway level; mapMECFS list not best-matching in 4 of 4; curated flags not independently recovered; demo autonomic links drug-target-driven | `test3_*.csv` |
| Text mining | 3 classes below 0.8 precision after revision, 6 undetermined; only 1 of 26 classes (autonomic_testing) clears 0.8 at the record-clustered lower bound; precision judged by a language model; questionnaire filter false-exclusion 0.77 | `measurement_precision_summary.csv`, `measurement_pro_filter_summary.csv` |
| Geography | county desert top lists concentrated in small counties; no cross-source relationship survives spatial correction; the burden/SVI label shuffles are near-tautological (they confirm the composite uses those components) | `test4_small_population_check.csv`, `test6_geography_relationships.csv`, `test6_scoring_label_shuffle.csv` |
| Facility matcher | technology-experience measurement desert not supported (0.003 vs 0.005, p 0.26); entity-resolution recall 0.528 | `test6_clinic_shuffle.csv`, `facility_match_validation.csv` |
| Scoring | clinic_capacity not measurement-specific (shuffle Spearman 0.965; 0.665-0.994 across bundles); with all facility locations shuffled the composite ranking keeps Spearman 0.94 and 8.2 of the top-25 (as much change as T5-3's 8 of 25); the capillaroscopy swap leaves the ranking nearly unchanged (composite Spearman 0.994, top-10 9 of 10); top-10 unstable (0 of 10 with P(top 10) >= 0.5); 1,825 of 2,395 counties have no Long COVID or ME/CFS trial and no NIH core project in reach | `test6_scoring_clinic_shuffle.csv`, `test7_sensitivity_comparisons.csv`, `test7_summary.csv`, `test5_readiness_zero_counts.csv` |
| Reproduction | Uwakwe HR-model AUROC not reproducible from the release; symptoms-only model 0.862 vs 0.907 published | `stanford_reproduction_benchmark.csv` |
| Adapter | heart-rate adapter vs Long COVID label: point estimate in the opposite direction to the hypothesis (AUROC 0.387, CI 0.284-0.494 excludes 0.5; permutation p 0.069) | `wearable_hr_adapter_label_check.csv` |
| Deployment example | the 20 candidate sites of the rank-1 county are about 15 distinct organisations (entity-resolution recall 0.528) | `results/example_deployment_recommendation.json` |
| Person-linked omics (mapMECFS) | 0 of 5 layers pass (p < 0.01, Bonferroni); CSF metabolomics 0.825 (p 0.0130) falls to 0.705 (p 0.0619) after in-fold age removal; CSF + serum fusion below serum alone (-0.049, -0.144 to 0.029); CSF and serum disagree on which proteins differ (rho 0.036) | `results/LINKED_OMICS_RESULTS.md` |
| GEO cohorts | 0 of 4 primary cross-cohort transfers replicate (0.356-0.561); GSE270045's 0.998 confounded with subject-id source (p 1.5e-08); ME/CFS cfRNA below its covariates (-0.026, -0.077 to 0.095); 5 of 7 within-cohort classifiers not permutation-significant. No write-up file exists | `geo_cohort_transfer.csv`, `geo_cohort_confounds.csv`, `geo_cohort_classifier.csv` |
| MUSCLE-ME steps | worst case for the 8 missing controls 0.612 (0.459-0.763), null; steps no better than CPET (-0.043, -0.192 to 0.092); Long COVID vs ME/CFS 0.577 (0.414-0.726); Long COVID vs controls "weak" | `results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md` |
| Fibromyalgia thermography | primary Delta over age + BMI +0.053 (-0.028 to 0.129), rule not met; 25% sensitivity at 90% specificity | `results/FM_THERMOGRAPHY_RESULTS.md` |
| NHANES labs vs umbrella labels | Rx fibromyalgia, migraine, insomnia, myalgia: no increment over demographics; pairwise contrasts between them null; ME/CFS-like proxy increment explained by smoking and exclusion diagnoses (+0.005, -0.020 to +0.032); symptom-label increments mostly BMI and smoking | `results/LABS_VS_DIAGNOSIS_RESULTS.md` |
| NHANES 58 lab layers | lab blocks over clinical: 1 of 47 rows supported, 11 worse; per-layer screen 0 of 292 supported | `results/NHANES_LAYERS_RESULTS.md` |
| Appelman 2024 metabolomics | plasma 0.694 (p 0.063) and muscle 0.674 (p 0.096) null; 0 of 83 plasma metabolites at q < 0.05; adding plasma to VO2max -0.021 (-0.187 to +0.142) | `results/LABS_VS_DIAGNOSIS_RESULTS.md` Part 2 |
| NHANES 2003-2006 replication | 4 of 12 verdicts do not replicate: fatigue Test 2, depression Tests 1 and 2, ME/CFS-like proxy Test 1 (36 cases) | `results/NHANES_2003_2006_RESULTS.md` |
| MY-LC cortisol | reproduces in the released rows but not in two other published cohorts (PLE-042, PLE-043); batch not checkable | `results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md`, `docs/PUBLISHED_LAB_EVIDENCE.md` |
| Published evidence | 0 of 17 published device diagnostic-accuracy rows externally validated; RECOVER's 25 routine labs null (PLE-047); pooled Long COVID HRV null (PDE-038) | `docs/PUBLISHED_DEVICE_EVIDENCE.md`, `docs/PUBLISHED_LAB_EVIDENCE.md` |

---

## 19. Limitations

* **Labels.** No wearable dataset of §6-7 carries a clinical diagnosis of any target condition. Their strongest
  condition label is Uwakwe's self-reported Long COVID (126 participants, 31 cases); ME/CFS is a constructed NHANES
  proxy without post-exertional malaise, duration or clinical assessment; fibromyalgia, migraine and IBS are
  prescription reason codes (treated people, 2013-14 only). The datasets of §14 carry criteria-based or clinic labels,
  but only against healthy, surgical or general-population controls.
* **Re-analysis is not replication.** Every positive result of §14 re-analyses the authors' own released rows; agreement
  with the paper (MY-LC 0.96, ARG1 0.848, thermography 18 of 18 calls) is reproduction. Appelman 2024 and MUSCLE-ME
  are the same trial. No §14 result was tested in an independent cohort, and MY-LC cortisol was not reproduced by two
  other published cohorts (PLE-042, PLE-043).
* **Controls are not clinical look-alikes (spectrum bias).** No dataset processed compares a measurement between an
  invisible illness and the conditions it must be told apart from (deconditioning, depression, other chronic pain or
  fatigue); AUROCs against healthy, surgical or general-population controls overstate diagnostic use.
* **Incorporation bias.** Where the measurement is part of the label, "accuracy" is agreement with the definition:
  activity limitation in the Long COVID and ME/CFS case definitions (MUSCLE-ME steps), tryptase and urinary mediators
  such as LTE4 for MCAS, the orthostatic heart-rate rise for POTS (§14.5).
* **No open MCAS or POTS lab data.** No open per-person dataset has MCAS mediators or POTS catecholamines,
  autoantibodies or blood volume with controls, and in All of Us the posture-coded tests are held by <= 20 people
  (§14.6).
* **Small n.** 21-42 people per mapMECFS layer, 15-93 cases per GEO subset, 32 surgical controls for ARG1, 22 controls
  with steps in MUSCLE-ME, 36 proxy cases in NHANES 2003-2006; most §14 datasets come from one site.
* **Ethics: two public records with apparent patient identifiers were not used.** During lab discovery, searchers'
  schema probes found that figshare 30127723 (Seibert 2026, PLOS ONE, post-COVID GPCR autoantibodies, CellTrend ELISA,
  vs recovered controls) has populated columns for surname, first name, date of birth, telephone number and GP, and
  that Mendeley vzydrgtymt (Anwar 2023; serum and follicular vitamin D and zinc in infertile women with and without
  endometriosis, Indonesia) contains patient names and medical record numbers. For the figshare record the probe
  deleted its copies after seeing only column headers and non-null counts. Neither record was downloaded into
  `data/raw/lab_candidates/` or ingested, and both are marked `excluded_privacy` in `lab_dataset_candidates.csv`
  (`docs/LAB_DATASET_DISCOVERY.md` "What was done", which suggests notifying the publishers; no notification is
  recorded).
* **No autonomic signal in the population dataset.** NHANES 2011-2014 is wrist accelerometry only; nothing in it
  evaluates autonomic monitoring. Uwakwe heart rate is step-masked, whole-record and not aligned to infection.
* **Cross-sectional, unweighted models.** NHANES prediction metrics describe the unweighted wearable-valid sample, not
  the U.S. population; associations have no temporal direction (mortality excepted).
* **mapMECFS partial.** Per-participant actigraphy, HRV, tilt and CPET require an account; wearable x omics linkage is
  untestable with open data.
* **Molecular evidence is condition-level** and not independent across sources (shared curation, literature co-mention,
  drug-target families); it is not mechanism and not patient multi-omics.
* **Measurement discovery** is bounded by a curated 26-class vocabulary; precision was judged by a language model, not a
  domain expert, from snippets alone, on at most 30 mentions per class and only on demo-cluster records; the "meets 0.8"
  calls are provisional for every class but one because their lower confidence bounds fall below 0.8; recall is
  unmeasured.
* **Burden.** Long COVID county values are inherited state values (HPS is an experimental, low-response survey that
  ended 2024-09-16); CMS MMD covers Medicare fee-for-service only and publishes integer percentages; 8 of 14 conditions
  are level D, including POTS and dysautonomia, whose rankings contain no burden information.
* **Locations are centroids.** Facilities at ZIP/ZCTA or city centroids; the 50 km radius is measured from the county's
  Census internal point, a weak proxy for where people live in large counties.
* **Curated assumptions.** Condition specialty groups and measurement implementer groups (`configs/relevance.yaml`); an
  NPPES taxonomy does not mean a clinician evaluates or treats the condition. Primary care is 86-92% of the relevant
  NPIs, so provider density mostly measures primary-care supply.
* **Composite design.** The diagnostic desert contains burden and vulnerability, and its access terms oppose
  clinic_capacity and research_readiness, so the equal-weight composite partly trades them off.
* **Rankings are unstable** under burden uncertainty and weights (§13); the short list is a set of candidates, not an
  ordering to act on without the intervals.
* **Provenance of report tables.** Every processed table carries the provenance columns (`validate`: 157/157), but
  most CSV extracts in `results/tables/` do not: 131 of the 159 CSVs of the original build had no provenance column
  (133 lacked the full set of 9; `number_audit.csv` included; `docs/GUARDRAIL_AUDIT.md` §4); with the CSVs of §14:
  264 CSVs, 236 of them without a provenance column (238 without the full set of 9; recount of the CSV headers against
  `provenance.PROVENANCE_COLUMNS` on 2026-09-24). They are generated from processed tables that do.
* **Number audit coverage.** The 3,208-claim number audit (`number_audit.csv`) predates §14. The §14 result files
  record their own independent or adversarial reviews for linked omics, thermography and MUSCLE-ME only; §14 of this
  report is copied from the result files and tables and has not been audited separately.
* **Rankings are national.** `rank_deployment_opportunities` has no state filter; a question about one state gets the
  national ranking (the agent lists that state's counties among the national top 50 and reports a within-state
  ranking as UNKNOWN).
* **Measurement-agnostic rankings.** Burden, vulnerability and the diagnostic desert do not depend on the measurement,
  and clinic capacity is mostly facility density, so swapping the measurement changes the short list only marginally
  (§1 row 9).
* **Language check is lexical.** It catches the listed product-language terms, not every paraphrase.

---

## 20. How CAPRIO plugs in later

Interface (`src/measure_it/measurements/adapter.py`, unchanged since the scaffold):

```python
class MeasurementAdapter(ABC):
    def preprocess(self, raw_data) -> pd.DataFrame        # raw person-level data -> cleaned frame keyed by participant_id
    def embed(self, raw_data) -> pd.DataFrame             # -> one row per participant, numeric embedding columns
    def phenotype_score(self, embedding) -> pd.DataFrame  # -> one column per phenotype axis
    def describe_measurement(self) -> MeasurementDescription
```

| adapter | status | measurement ids | phenotype scores |
|---|---|---|---|
| `WearableAdapter` | implemented (NHANES reference) | accelerometry | low_activity, circadian_disruption, activity_fragmentation, sleep_proxy_irregularity |
| `WearableHeartRateAdapter` | implemented (personal baseline) | wearable_heart_rate, accelerometry | rhr_elevation_vs_personal_baseline, activity_reduction_vs_personal_baseline |
| `CapillaroscopyAdapter` | **stub**, `requires_private_data = True`; data methods raise NotImplementedError | capillaroscopy, microvascular_function | planned: capillary_rarefaction, abnormal_capillary_morphology, microhemorrhage_burden, reduced_capillary_flow |

**Path for CAPRIO:** nailfold image -> `preprocess` (QC'd capillary fields) -> `embed` (microvascular embedding) ->
`phenotype_score` (capillary phenotype scores against a healthy reference) -> `scores_to_long` -> the same ontology,
molecular, geography, clinic and scoring layers. `describe_measurement()` already works for the stub, so implementer
groups, curated deployment complexity, FDA context (capillaroscopy: no product code; microvascular_function: 3 codes)
and linked conditions resolve before any data exist. The `nailfold_capillaroscopy` bundle is already scored through the
identical code path (S5), and the contract test (`tests/test_adapters.py`) runs a synthetic adapter through every
downstream function and asserts that none names a modality: that is the modularity result. The same S5 run is also a
negative result for the product question "given technology Y, where should we deploy Y?": against the wearable bundle
the capillaroscopy ranking has composite Spearman 0.994 and top-10 overlap 9 of 10 (clinic_capacity 0.994,
research_readiness 0.953; `test7_sensitivity_comparisons.csv`), although capillaroscopy has 4 trials, 0 NIH core
projects and no FDA product code. Until a measurement's own person-level evidence, implementer network and research
activity differ enough to move the components, the deployment list barely depends on which measurement is deployed. Since
2026-09-25 the `evidence_weighted` ranking makes this explicit instead of hiding it: a measurement enters with its own
performance record, and the capillaroscopy bundle, which has none, is flagged UNKNOWN and not ranked rather than ranked
as if it were equivalent to the wearable bundle (`results/SCORING_RESULTS.md` §12). A CAPRIO validation study that
reports AUROC and sensitivity at a fixed specificity would enter `measurement_performance` as a tier-1 record.

**CAPRIO must supply** (`docs/ADAPTERS.md` §5): images and acquisition metadata with namespaced participant ids; the
versioned embedding model; phenotype definitions with a healthy reference cohort; a non-stub
`describe_measurement()`; ontology hooks (a `phenotype_axes` row with an HPO term, or a declared link to an existing
axis); `register_adapter("capillaroscopy", ...)`; its own reproduction validation. **What would not change:** the
measurement registry and discovery, ontology, molecular layer, geography, clinic and research-site registries and
matching, scoring and its sensitivity analyses, MCP tools, dashboard and the `participant_adapter_scores` schema.
Nothing about capillaroscopy performance is known or implied by the stub; in registered trials and grants it is
close to absent (4 trials, 0 NIH core projects).

---

## 21. Reproducibility

Full account: `docs/REPRODUCIBILITY.md`.

```bash
uv sync --frozen                                    # environment pinned in uv.lock (Python 3.12)
uv run measure-it pipeline --workers 16             # fresh clone: every source from the network (not yet run end to end)
uv run measure-it pipeline --offline --workers 16   # 75 steps from data/raw + data/_http_cache + data/interim only
uv run measure-it validate                          # provenance, required datasets, person-layer geography, ids, language, audits
uv run pytest -m "not network"                      # "not data" on a clone without data/processed
```

* **Offline guarantee.** `--offline` answers every HTTP call from the disk cache, refuses downloads not in a source's
  `MANIFEST.json` (url, bytes, sha256, retrieved_at) and installs a socket guard in every process. It reuses the derived
  caches in `data/interim` (rebuilt from `data/raw` when absent; `docs/REPRODUCIBILITY.md` §4).
* **Determinism.** The final code was re-run end to end offline (run `20260924T150209Z`, all 49 steps, 16 workers) and compared with the outputs it replaced (`measure_it.reproduce compare`; `results/pipeline_runs/20260924T150209Z_vs_pre_final.csv`): 279 processed and result tables identical (run timestamps normalised); 4 changed by this report's fixes only (`deployment_opportunities.member_components`: the `nan` strings cleared; `deployment_candidates.recommendation_json`, `demo_top10_candidate_sites.csv`, `facility_candidate_examples.csv`: the reworded absence note and technology-experience reason, and the phenotype-signal verdict); 3 in runtime fields only; and the canonical `phenotype_signatures` (441 -> 28,516 rows), because two contributors working in the same tree added `phenotype_signatures__geo_cohorts` and `__mapmecfs_omics` during the run (their 29 new tables were not pipeline outputs then; they are since the integration below). Every results/*.md write-up regenerated identical apart from generation times. 48 of the 49 steps ended ok; the final `validate` step failed on a registry entry of that work (fixed since; validate passes, see below). Before this, two full runs of the then 47-step DAG (`20260924T045326Z`, `20260924T062847Z`; the two figure steps were added later, so figure outputs were not among the 280 compared tables) matched in all 280 compared tables (231 identical, 46 identical after run-timestamp normalisation, 3 differing only in runtime fields), including Monte Carlo, permutation, bootstrap and shuffle outputs (`docs/REPRODUCIBILITY.md` §5)
* **Integration of §14 (2026-09-24, after the determinism run above).** The §14 modules are now 20 pipeline steps
  (74 in the DAG, all before `canonical_unions`; offline switches read saved raw files and texts). Run
  `20260924T215649Z` (the 20 steps plus `digital_person`, offline, 24.2 min) then `20260924T222227Z` and
  `20260924T223309Z` (Data.gov rules, `scoring_report`, `figures_data`, unions, registry merge, Figure 2, DuckDB,
  validation), compared with a copy taken before (`results/pipeline_runs/20260924_integration_vs_pre_integration.csv`):
  422 tables identical and 3 numerically equal, so every §14 output reproduced its standalone run; the 12 changed
  tables are the intended ones (the canonical person-level and `phenotype_signatures` unions now hold the §14 and NHANES
  2003-2006 partitions; the Data.gov log for the 10 new sources; the example recommendation's phenotype block; a caveat
  text; runtime fields). No deployment rank, score or site changed. No full 75-step run has been made since
  (`docs/REPRODUCIBILITY.md` §5.8).
* **Timing.** run `20260924T150209Z`: 104.4 min wall time with 16 workers on a host shared with other jobs; the reproduction run `20260924T062847Z` took 91.4 min. The NHANES and Stanford model steps are the long poles (84 and 91 min in the final run)
* **Inputs.** `data/raw` (40 GB, 37 sources), `data/_http_cache` (550 MB) and `data/interim` are not in git, and no data
  bundle is published. On a fresh clone the same pipeline without `--offline` fetches them from the agencies; this cold
  path has **not** been run end to end from an empty checkout, and sources that change (ClinicalTrials.gov, RePORTER,
  openFDA, Open Targets, ...) will return newer data (`docs/REPRODUCIBILITY.md` §8; README "Reproducing from scratch").
* **Checks run for this report (2026-09-24, after the §14 integration):** `uv run measure-it validate` (run
  `20260924T223309Z`, 22:34 UTC): **PASS**, 6/6 checks (provenance 157/157 processed tables; required datasets 11/11;
  person-layer tables without a geographic column 65/65; object ids 14/14; forbidden product language 0 unquoted uses in
  104 text files, this report included, and 66 MCP/API texts; DATA_AUDIT.md 37/37); `build_db` exit 0 with 21 of 21
  SPEC views; `uv run pytest -m "not network"`: **1,383 passed, 0 failed** in 233 s
  (`results/pipeline_runs/pytest_20260924_integration.txt`). The 4 failures of the earlier check (Data.gov rules,
  Figure 2 height, GEO-cohort trace source) are fixed (§23 item 12). `pytest -m "not network and not data"` selects 789
  tests, all passing in this tree; the earlier check that the no-data tests pass on a copy without any data (693 then)
  was not repeated.
* **Number audit.** 3,208 numeric claims in the other docs and results files, 2 dashboard claims and 3 claims of one of
  the 27 DATA_AUDIT files were re-derived from the tables: 3,170 matched, 14 were stale and fixed, 24 could not be
  checked (mostly latency timings on a shared host) (`results/tables/number_audit.csv`). **This report was not part of
  that audit** (0 rows); its numbers were checked by the 2026-09-24 review, whose corrections are included here. The
  §14 numbers were copied from the §14 result files and tables and were not part of either check.

---

## 22. Success criteria (SPEC 1-10)

Status from `docs/SPEC_COMPLIANCE.md` §13 (a pre-report snapshot), updated for the final state, with the qualifier
that the evidence requires.

| # | criterion | status | evidence and qualifier |
|---|---|---|---|
| 1 | At least one real public patient-level wearable dataset processed | **met** | NHANES 2011-2014 (12,955 wear-valid participants); three Stanford releases (2,360 participants with wearable features); added in §7.5 and §14: NHANES 2003-2006 hip accelerometry (6,283 adults in the analysis population) and MUSCLE-ME hip-accelerometer daily steps (73 people) |
| 2 | Wearable features linked to real clinical or disease information | **met** | NHANES questionnaires, labs, prescriptions and linked mortality on SEQN; Uwakwe Long COVID label. The disease labels of these wearable datasets are self-report, proxy or prescription-based, not clinical diagnoses. MUSCLE-ME links steps to criteria-based Long COVID and diagnosed ME/CFS labels in the same people, but only against 22 healthy controls with steps, and activity limitation is part of those labels (§14.3) |
| 3 | At least one invisible-illness example demonstrated | **met, with null / proxy evidence** | Long COVID / ME/CFS x wearable monitoring runs end to end; its wrist-wearable evidence is a null (Long COVID) and a proxy (ME/CFS). MUSCLE-ME hip steps separate diagnosed patients from healthy controls on complete cases (0.831) but not in the locked missing-control worst case (0.612); cortisol, ARG1 and a serum proteome are candidates from re-analyses against healthy or surgical controls; person-linked omics and routine labs are null (§14) |
| 4 | Ontology connects clinical concepts to molecular evidence | **met** | 14 conditions with MONDO ids (14/14), ICD-10-CM 13/14, EFO 13/14, MeSH 13/14; HPO terms only for 6 of 14; molecular rows carry the ontology id; the §14 datasets and published-evidence rows are normalised to the same condition ids |
| 5 | Candidate measurement technologies discovered automatically | **met within a curated vocabulary** | discovery runs for any registry condition; it cannot find a technology type outside the 26 classes; precision is an LLM-reviewer spot-check, recall unmeasured |
| 6 | U.S. public geographic data queried and normalised | **met** | 9 geographic sources on 2024 county/state FIPS with a ZCTA crosswalk; CDC PLACES at county and ZCTA only |
| 7 | Providers, facilities and trial sites mapped | **met** | Figure 7, dashboard overlays; coordinates are centroids for most facilities; the dashboard facility overlay shows HRSA sites only |
| 8 | Reproducible geographic deployment opportunities | **met (offline)** | offline 75-step pipeline; the final run on the final code reproduced every table except those changed by the review fixes, runtime fields and the §14 partitions added during the run (§21); a cold network build has not been run end to end. The rankings are reproducible but not stable under weights and burden uncertainty, and barely depend on the measurement |
| 9 | Every recommendation traceable to source data | **met** | `trace_evidence` over 20 namespaces; opportunity -> burden rows -> raw files with sha256; person-level ids refused |
| 10 | CAPRIO insertable as a new MeasurementAdapter without restructuring | **met (interface level)** | stub adapter, contract test and S5 swap; no CAPRIO data exist in this build, so nothing about capillaroscopy performance is shown |

---

## 23. Open items and state of the build

1. **Working tree not committed.** The last commit is `209f666`; the lab-dataset work, the §14 pipeline integration,
   the dependency additions (`xlrd`, `mat-io` in `pyproject.toml` / `uv.lock`) and this report are uncommitted working-tree
   changes on top of it. Commit the whole tree before anyone clones, then repeat the quickstart from the clone.
2. **Fresh-clone path.** A clone contains no data and no data bundle is published; the cold network build has not been
   run end to end from an empty checkout, and its memory peak and duration are unknown (README "Reproducing from
   scratch"). Publishing `data/processed` plus the DuckDB (about 2 GB) would let a new user run `ask`, the API, the
   dashboard, the MCP server and the tests without the 40 GB raw build.
3. **Scientific gap (data, not code):** no public invisible-illness label shows an autonomic wearable phenotype; demo
   story step 1 is partial. The next evidence has to come from a deployment that measures it (or from mapMECFS
   per-participant files with an account).
4. **mapMECFS per-participant physiology** needs a registered account; linkage within mapMECFS is untested.
5. **Optional or unprocessed inputs:** NHANES raw 80-Hz mode; PhysioNet/CovIdentify/Apple Watch benchmarks; a
   temporal/deep model; CDC PLACES tract and place levels; the Stanford phase-1 sleep files, which are downloaded but not
   processed (`docs/SPEC_COMPLIANCE.md` A2).
6. **Ontology coverage:** HPO terms empty for 8 of 14 conditions (no direct Monarch annotation); `ptlds` lacks an
   ICD-10-CM code, `mcas` an EFO id, `post_infectious_syndrome` a MeSH id.
7. **Interface deviations.** Compatible supersets: `get_condition_burden` defaults `geography_level` to "state";
   `find_candidate_clinics` makes `measurement` optional. Restrictions: dashboard page 1 shows one area layer at a time
   (burden, desert or opportunity are mutually exclusive, not independently toggled) and its facility overlay shows
   HRSA sites only (partial against the SPEC); `rank_deployment_opportunities` ranks nationally and has no state
   filter; a default top-10 ranking is about 90 kB of JSON, above Claude Code's default MCP output limit unless
   `MAX_MCP_OUTPUT_TOKENS` is raised (README).
8. **Report tables without provenance columns:** 236 of the 264 CSV extracts in `results/tables/` have no provenance
   column (§19); they are generated from processed tables that carry them.
9. **Agent molecular sentence (fixed 2026-09-24).** The `ask` agent listed every enriched system as "supported
   physiological systems" (for ME/CFS: autonomic_cardiac, connective_tissue, immune, mitochondrial_metabolic, neuronal,
   vascular_endothelial), without the post-hoc drug-target flag that `get_molecular_context` returns, which contradicted
   §8 in user-facing output. It now prints three tiers per condition (§16). The documented agent answers in
   `docs/MCP.md` §7 were re-run.
10. **`technology_evidence[*].phenotype_signal_strength` (fixed 2026-09-24).** It was null in all 140 rows of
    `results/example_deployment_recommendation.json` because the code read a key the dimension does not have; it now
    carries the verdict (after the final run: long_covid x wearable_heart_rate and x the bundle `null_result`, me_cfs x accelerometry and x the bundle `mixed`, the other members UNKNOWN / NOT AVAILABLE; test added).
11. **Stored `nan` strings (fixed by the final run).** The 12,876 `deployment_opportunities.member_components` rows that held the pre-fix `"nan"` string for level-D set members were rebuilt by run `20260924T150209Z`; 0 of 76,680 rows hold one now.
12. **§14 integration (done 2026-09-24).** The GEO-cohort, linked-omics, published-evidence, device and lab modules
    are pipeline steps (§21); Data.gov discovery rules exist for the 10 new sources (each query belongs to one
    source; topical hits on a journal deposit count as related, never as the dataset being indexed); Figure 2 fits the
    247 mm height with 37 sources (compact one-line rows in bands whose rows set the height); `trace_evidence` maps
    every signature dataset to its source (`nhanes0306` to NHANES 2003-2006, `geo_GSE*` to `geo_cohorts`, the lab and
    device releases to their own source ids) and refuses their participant ids; the digital-person unions include every
    person-level partition, datasets kept separate. **Still open:** `measurement_phenotype_signal` and the Phase 3
    table (`phase3_measurement_evidence`, `results/PHASE3_MEASUREMENT_EVIDENCE.md`) read only the NHANES 2011-2014 and
    Stanford results under their pre-specified plan, so the MUSCLE-ME accelerometry, thermography and lab results
    reach users through `phenotype_signatures` / `get_patient_phenotype_signature` and §14, not through
    `phenotype_signal_strength`; extending that plan would be a new, dated deviation. No `results/GEO_COHORT_RESULTS.md`
    write-up exists (§14.4 cites the tables). `docs/MCP.md` captured tool outputs predate the §14 signature partitions.
13. **Fixes to documented commands (2026-09-24, from a fresh-clone check).** `ask --json` and
    `python -m measure_it.tools` printed a `[scoring]` progress line before the JSON (now on stderr; test added);
    `measure-it dashboard` exited silently without a TTY on Streamlit's first-run prompt (now headless when stdin is
    not a terminal, bound to 127.0.0.1, and replaced by the Streamlit process so `kill <pid>` stops it); the `datagov`
    step crashed when every Data.gov query failed (now records the source `blocked`/`partial`; test added);
    `build-db` exited 0 with missing SPEC views and `validate` passed 0/0 checks on an empty checkout (both now fail);
    `ask` on an empty checkout printed a traceback (now one line); a failed pipeline step now prints its last error
    line and log path; run logs no longer contain absolute home paths; `SCORING_RESULTS.md` no longer assumes the
    checkout directory name; two tests that need processed tables are now marked `data`; offline
    `search_us_open_data` answers a cached query regardless of case and spacing; the matcher's absence note names its
    condition.

---

## Appendix: where everything is

| topic | file |
|---|---|
| Brief and conventions | `docs/SPEC.md`, `docs/CONVENTIONS.md` |
| Compliance and audits | `docs/SPEC_COMPLIANCE.md`, `docs/GUARDRAIL_AUDIT.md`, `results/tables/number_audit.csv` |
| Sources | `SOURCE_REGISTRY.yaml`, `data/raw/<source_id>/DATA_AUDIT.md` |
| Linkage audit | `docs/MAP_MECFS_LINKAGE_AUDIT.md` |
| Wearable results | `results/WEARABLE_NHANES_RESULTS.md`, `results/WEARABLE_STANFORD_RESULTS.md` |
| Measurement discovery, Phase 3 | `results/MEASUREMENT_DISCOVERY.md`, `results/PHASE3_MEASUREMENT_EVIDENCE.md` |
| Molecular coherence | `results/MOLECULAR_COHERENCE.md` |
| Geography, burden definitions | `results/GEOGRAPHY_RESULTS.md`, `docs/BURDEN_DEFINITIONS.md` |
| Facilities | `docs/FACILITY_MATCHING.md` |
| Scoring | `docs/SCORING.md`, `results/SCORING_RESULTS.md`, `results/example_deployment_recommendation.json` |
| MCP, API, dashboard | `docs/MCP.md`, `docs/DASHBOARD.md` |
| Adapters and CAPRIO | `docs/ADAPTERS.md` |
| Reproducibility | `docs/REPRODUCIBILITY.md`, `results/pipeline_runs/` |
| Analysis plans (written before each analysis) | `docs/ANALYSIS_PLAN_*.md`; GEO cohorts: `configs/geo_cohorts.yaml` (`analysis_plan` block) |
| NHANES 2003-2006 replication (§7.5) | `results/NHANES_2003_2006_RESULTS.md`, `nhanes0306_*.csv`, `data/raw/nhanes_2003_2006/DATA_AUDIT.md` |
| Person-linked omics, mapMECFS (§14.4) | `results/LINKED_OMICS_RESULTS.md`, `linked_omics_*.csv` |
| GEO case/control cohorts (§14.4) | `results/tables/geo_cohort_*.csv` (no `results/GEO_COHORT_RESULTS.md` exists), `data/raw/geo_cohorts/DATA_AUDIT.md`, module docstring of `src/measure_it/omics/geo_cohorts.py` |
| Device dataset discovery and device analyses (§14.3-14.4) | `docs/DEVICE_DATASET_DISCOVERY.md`, `device_dataset_candidates.csv`, `results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md`, `results/FM_THERMOGRAPHY_RESULTS.md` |
| Lab dataset discovery and lab analyses (§14.3) | `docs/LAB_DATASET_DISCOVERY.md`, `lab_dataset_candidates.csv`, `results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md`, `results/ENDO_ARG1_REPOD_RESULTS.md`, `results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md` |
| Labs vs diagnosis, NHANES lab layers (§14.4) | `results/LABS_VS_DIAGNOSIS_RESULTS.md`, `labs_dx_*.csv`, `results/NHANES_LAYERS_RESULTS.md`, `nhanes_layers_*.csv` |
| Published evidence (§14.7) | `docs/PUBLISHED_DEVICE_EVIDENCE.md`, `docs/PUBLISHED_LAB_EVIDENCE.md`, `data/raw/published_device_evidence/DATA_AUDIT.md`, `data/raw/published_lab_evidence/DATA_AUDIT.md` |
| Data access plan (§14.8) | `docs/DATA_ACCESS_PLAN.md`, `controlled_cohort_counts.csv`, `data/raw/controlled_cohorts/`, `data/raw/snyder_candidates/` |
| Figures | `visuals/figures/` (published copies, with `CAPTIONS.md`); the pipeline writes them to `results/figures/` and `results/maps/` |
