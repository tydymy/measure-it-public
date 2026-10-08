# SPEC compliance checklist (deliverables audit, 2026-09-24)

> **Snapshot note (added after the final report, 2026-09-24).** This checklist was written before
> `results/FINAL_REPORT.md` and `README.md` existed and before the final audit fixes. Its section 0 numbers are those of
> the snapshot (1,177 tests; the language check then scanned 38 text files); the current values are in
> `results/FINAL_REPORT.md` §20 and §22. Updated here after the report's review: FINAL_REPORT.md and README.md now
> exist; demo-story steps 3 and 7 and dashboard page 1 are **partial** (the report's §2 gives the reasons; rows
> updated below are marked). For the current status of every component read `results/FINAL_REPORT.md` §2 and §21; this file is the
> requirement-by-requirement evidence as audited at the time.
>
> **Update 2026-09-24 (integration of the labs, omics and device analyses).** Since this snapshot: 37 registered sources
> (36 ingested, 1 partial; 10 added for `results/FINAL_REPORT.md` §14, each with DATA_AUDIT.md and registry entry);
> the pipeline has 75 steps (5 NHANES layer / 2003-2006 steps and 20 steps for the §14 modules; `docs/REPRODUCIBILITY.md`
> §1, §5.8); `validate` PASS 6/6 (provenance 157/157 processed tables, required datasets 11/11, person-layer tables
> without geography 65/65, object ids 14/14, product language 0 unquoted uses in 104 text files, DATA_AUDIT.md 37/37);
> `pytest -m "not network"` 1,383 passed, 0 failed (789 of them need no data); `pyproject.toml` / `uv.lock` gained `xlrd`
> and `mat-io`. Rows A1-A4, B-C and the success criteria below are unchanged in status; the §14 person-level re-analyses
> (MUSCLE-ME hip steps + CPET, MY-LC cortisol, endometriosis ARG1, hEDS/HSD Olink, mapMECFS person-linked omics, GEO
> cohorts, NHANES labs vs labels, Appelman metabolomics, fibromyalgia thermography, NHANES 2003-2006 replication) add
> evidence to SPEC section A ("optional public wearable benchmarks" is still not done as specified) and are summarised
> in the report, not re-audited row by row here. Statements below that no processed dataset carries a clinical
> diagnosis label were true of the wearable datasets audited then; several §14 datasets carry clinical labels.

This file checks every concrete requirement in `docs/SPEC.md` against the repository as it stands on 2026-09-24.
Each row records a status, the evidence that was checked (file, table, test or command output) and the gaps. A row
is marked **done** only when the evidence was opened or run during this audit. Row counts are from the processed tables
on disk, not from other documents.

Status legend: **done**, **partial** (part of the requirement is met; the rest is named in "gaps"), **missing**,
**deliberately not done** (a documented decision), **not done (optional)** (the SPEC makes the item optional).

## 0. What was run for this audit

| check | command | result |
|---|---|---|
| repository validation | `uv run measure-it validate` | **PASS**, 6/6 checks: provenance 122/122 processed tables; required datasets 11/11; person-layer tables with no geographic column 32/32; object ids unique and namespaced in 14/14 row-identity tables; forbidden product language 0 unquoted uses in 38 text files; DATA_AUDIT.md for 27/27 registry sources |
| test suite | `uv run pytest -m "not network"` | 1,177 passed, 0 failed, 0 skipped (network tests deselected). The run printed progress dots only, because `-q` is also in `addopts`; the count is the number of dots. |
| MCP server, in-memory | FastMCP `Client(get_server())`, `list_tools()` | 18 tools: the 15 SPEC tools + 3 extras. Leading parameters match the SPEC signatures (section 4). |
| MCP server, stdio | `uv run --frozen measure-it serve-mcp` with `MEASURE_IT_OFFLINE=1`, via the FastMCP stdio client | 18 tools; `normalize_condition("G93.32")` returns `matched` -> `me_cfs` |
| deterministic agent | `MEASURE_IT_OFFLINE=1 uv run measure-it ask "What should we measure, where should we measure it, and who could realistically deploy it for Long COVID or ME/CFS with wearable autonomic monitoring?"` | 14.6 s. Returns sections for what to measure, where and who could deploy it, plus uncertainties and a provenance chain to raw files. Each sentence cites object ids. |
| generality | `tools.rank_deployment_opportunities("fibromyalgia", "accelerometry", "county")` and `("migraine", "retinal imaging", "state")` | both `ok`, "computed on the fly with measure_it.scoring.opportunity (same engine)". 2,395 counties are ranked. These combinations are not in the stored grid. |
| DuckDB | `information_schema` of `data/processed/measure_it_public.duckdb` | 21/21 SPEC tables present as `spec.*` views, each with the 8 SPEC provenance columns |

**Fixed during this audit:** `trace_evidence("opportunity:...")` listed 7 top-level sources. It omitted Census ACS,
whose B01003 total population sets rank eligibility (the population floor) and the per-100k clinic-capacity
denominators. `src/measure_it/mcp/trace.py` now adds `census_acs` (hint `acsdt5y2024-b01003.dat`). After the change,
`tests/test_trace.py`, `tests/test_mcp_server.py` and `tests/test_tools.py` all pass (106 passed).

## Summary

| SPEC section | done | partial | missing | not done (optional) / deliberately not done |
|---|---|---|---|---|
| Public data sources (22 rows) | 20 | 1 (CDC PLACES geographic levels) | 0 | 1 (A4 optional wearable benchmarks) |
| Core data model (21 tables + provenance) | 22 | 0 | 0 | 0 |
| Phases 1-6 | 6 | 0 | 0 | 0 (within Phase 1, the optional temporal/deep model was not attempted) |
| MCP tools (15) + UNKNOWN rule + output schema | 17 | 0 | 0 | 0 |
| Dashboard pages (5) | 4 | 1 (page 1: one area layer at a time; HRSA-only facility overlay; updated after the report review) | 0 | 0 |
| Adapter requirements (6) | 6 | 0 | 0 | 0 |
| DATA_AUDIT per source (27) + SOURCE_REGISTRY use | 28 | 0 | 0 | 0 |
| Reproducibility stack (13) | 13 | 0 | 0 | 0 |
| Validation tests 1-7 | 7 | 0 | 0 | 0 |
| Required final datasets (11) | 11 | 0 | 0 | 0 |
| Required figures (10) | 10 | 0 | 0 | 0 |
| Primary demo story (9) | 6 | 3 (steps 1, 3, 7; steps 3 and 7 updated after the report review) | 0 | 0 |
| Success criteria (10) | 10 | 0 | 0 | 0 |
| FINAL REPORT + README | 2 (written after this audit) | 0 | 0 | 0 |

At the time of this audit the largest open items were `results/FINAL_REPORT.md`, which did not yet exist, and
`README.md`, which was still the scaffold line; both were written afterwards. Section 15 lists every gap as audited.

## 1. Public data sources

| item | SPEC requirement | status | evidence | gaps |
|---|---|---|---|---|
| A1 NHANES 2011-2014 | demographics, medical/cardiovascular/diabetes/sleep/physical-functioning questionnaires, prescriptions, BP, BMI, labs, 7-day wrist accelerometry (minute/hour/day tables); joins only on official ids; one row per participant with derived wearable features; no invented HR/HRV | **done** | `data/raw/nhanes_2011_2014/`: DEMO, BPX, BMX, MCQ, BPQ, CDQ, DIQ, SLQ, PFQ, DLQ, HSQ, HUQ, DPQ, PAQ, SMQ, HIQ, INQ, RXQ_RX, CBC, BIOPRO, GHB, GLU, TCHOL, HDL, TRIGLY, VID, VITB12, THYROD, PAXHD/PAXDAY/PAXHR/PAXMIN (both cycles), linked mortality. `participant_wearable_features__nhanes`: 12,955 wear-valid participants x 54 columns. Features: mean/total daily MIMS, SD/CV of daily MIMS, sedentary_fraction, ASTP/SATP fragmentation, IS, IV, RA, M10/L5, active_period_duration_h, sedentary and active min/day, sleep_proxy_* window, weekday/weekend difference, valid wear time and rule. `participant_clinical_features__nhanes`: 19,931 x 170. `participant_id = nhanes:<SEQN>`. The DATA_AUDIT says "No heart-rate or HRV features exist or are derived"; the BPXPLS exam pulse is labelled as "NOT wearable heart rate". | The optional raw 80-Hz mode was not implemented (SPEC: "can be supported as an optional advanced mode"). |
| A2 Stanford Snyder Lab wearables | read the data dictionary first; HR baselines, RHR, deviation from baseline, steps, recovery trajectory, Long COVID discrimination; published analysis as a benchmark | **done** | Three releases are ingested: Mishra 2020 (118 participants), Alavi 2022 (2,123) and Uwakwe 2025 (126, of whom 31 have Long COVID). The DATA_AUDIT "Discovery" table records which release carries a genuine Long COVID label (Uwakwe only). `participant_wearable_daily__stanford_covid`. `results/WEARABLE_STANFORD_RESULTS.md`. `results/tables/stanford_reproduction_benchmark.csv`, `stanford_acute_recovery_km*.csv`, `stanford_temporal_leakage_checks.csv` (0 of 8 checks with violations). | On the only genuine Long COVID label, the wearable features give a **null** result (AUROC 0.504). The phase-1 sleep files are present but were not processed. |
| A3 NIH PI-ME/CFS / mapMECFS | attempt; `docs/MAP_MECFS_LINKAGE_AUDIT.md` with datasets, ids, subjects per modality, overlap counts and linkability; a linked subset only if true; pipeline works without it | **done** (registry status `partial`) | `docs/MAP_MECFS_LINKAGE_AUDIT.md`: 29 of 29 anonymous mapMECFS file downloads returned HTTP 403 "Login Required". Five open omics modalities link through the NIH study number (47 participants). **0 wearable x omics pairs** have comparable identifiers. Tables: `modality_inventory`, `modality_overlap`, `participant_omics_linked` (902,186 rows, omics-to-omics only) and `condition_molecular_evidence__mapmecfs` (14,565 rows, labelled condition-level). Pipeline step `ingest_mapmecfs` runs offline with `fetch=False`. | The per-participant actigraphy, HRV, tilt and CPET files need a mapMECFS account, so they are not used. |
| A4 Optional wearable benchmarks | PhysioNet, CovIdentify, Apple Watch sleep, "only if useful" | **not done (optional)** | No module, source or document uses them (a grep finds them only in SPEC.md). The Stanford acute-COVID releases (Mishra 2020, Alavi 2022) fill the signal-processing-check role through the true-onset vs placebo-window control (`stanford_acute_detection*.csv`). | Until this audit, no document recorded the decision. |
| B Disease/phenotype ontology | `condition_registry.parquet` with id, name, aliases, ICD-10-CM, MONDO, EFO, MeSH, HPO, parent, 8 flags, ontology provenance; the 12 listed conditions; public ontologies | **done** | `condition_registry` has 14 rows x 50 columns: long_covid, me_cfs, pots, dysautonomia, fibromyalgia, eds_hsd, mcas, lyme_disease, ptlds, post_infectious_syndrome, migraine, ibs, gastroparesis, endometriosis. Coverage: `icd10cm_codes` 13/14, `mondo_ids` 14/14, `efo_ids` 13/14, `mesh_ids` 13/14. It also has `hpo_terms`, `parent_condition`, the 8 flags (`infectious_or_post_infectious`, `autoimmune`, `autonomic`, `vascular`, `pain`, `fatigue_pem`, `gi`, `neurologic`) and `ontology_provenance`. Sources: `mondo_ontology` (Mondo v2026-09-01 + Monarch KG), `ebi_ols4`, `cms_icd10cm`. `condition_ontology_mappings` has 298 rows. | `hpo_terms` is filled for 6 of 14 conditions. The other 8 carry `hpo_annotation_status` = "no direct Monarch KG 2026-09-02 disease-phenotype annotation". `ptlds` has no ICD-10-CM code, `mcas` no EFO id, `post_infectious_syndrome` no MeSH id. The `icd10cm_blocks` column is empty for all 14 rows (not a SPEC field). |
| C Molecular evidence | `condition_molecular_evidence.parquet` from Open Targets, GWAS Catalog, GEO, SRA, MapME/CFS; per-edge fields; no arbitrary omics score | **done** | 70,295 rows. Open Targets 52,164 (target-disease 29,494; pathway 20,477; known drug 1,082; text-mined 667; genetic association 444). GWAS Catalog 2,772, GEO 537, SRA 257, MapMECFS 14,565. Columns: condition_id, ontology_id, evidence_type, entity_type/entity_id/entity_label, source_database, source_accession, evidence_direction, study_population, sample_size, date_retrieved, discovery_ancestry. `configs/geo_queries.yaml` holds the curated GEO/SRA queries. Graph: `measurable_biology` + `results/tables/fig5_nodes.csv` / `fig5_edges.csv`. MOLECULAR_COHERENCE.md: "No combined omics score is computed." | Open Targets is reached through GraphQL, not the Open Targets MCP (the SPEC allows either). |
| D Measurement registry | `measurement_registry.parquet`; the listed measurement classes; do not assume wearables only | **done** | 30 rows: 26 classes + 4 bundles. All 22 SPEC classes are present: wearable_heart_rate, ecg_ambulatory, ppg, hrv, accelerometry, posture_detection, continuous_spo2, respiratory_rate, sleep_objective, continuous_temperature, cpet, autonomic_testing, endothelial_function, vascular_imaging, capillaroscopy, retinal_imaging, digital_gait, blood_biomarkers, metabolomics, proteomics, immune_assays, microvascular_function. Extra classes: small_fiber_testing, digital_cognitive, transcriptomics, remote_patient_monitoring. | The class vocabulary is curated (`configs/measurements.yaml`); see success criterion 5. |
| D ClinicalTrials.gov | interventions, devices, outcome measures, endpoints, eligibility, phase/type, status, sponsors, facilities, city/state, technologies; not proof of effectiveness | **done** | `clinical_trials`: 7,910 (intervention_names, phases, study_type, overall_status, lead_sponsor, eligibility_criteria, ...). `trial_interventions` 13,021; `trial_outcomes` 54,022; `trial_sites` 47,482 (facility, city, state, zip, geopoint); `trial_conditions` 12,557. `measurement_mentions` has 10,646 rows (role, objective flag, snippet). The "not efficacy" caveat is in `results/MEASUREMENT_DISCOVERY.md` and in every trial tool. | - |
| D NIH RePORTER | grants for target conditions x wearable / biomarker / diagnostic / sensor / imaging / autonomic / microvascular / capillary / digital health / point of care; title, abstract, PI, institution, city/state/ZIP, funding, dates | **done** | `nih_projects`: 9,342 (project_title, abstract_text, pi_names, org_name/org_city/org_state/org_zipcode, award_amount, project_start_date/project_end_date). `nih_project_conditions`: 15,358. The augmentation terms are in `ingestion/nih_reporter.py` (lines 12-13, 93-99). | - |
| D openFDA device APIs | classification, 510(k), PMA, registration/listing; name, product code, class, specialty, clearance id, applicant, decision date; readiness signal only | **done** | `fda_device_classification` 90, `fda_510k` 8,101, `fda_pma` 132, `fda_registration_listing_counts` 95, `measurement_regulatory_status` 102 (with mapping_confidence and rationale). `get_regulatory_context` returns the caveat "deployment-readiness signal, not evidence of diagnosis". | - |
| E `geo_condition_features` | canonical geography; no finer than source resolution | **done** | 46,074 rows = 3,291 geographies x 14 conditions. Every row carries `burden_evidence_level`, `burden_source_resolution` and `burden_inherited`. | - |
| E CDC PLACES | county, place, tract, ZCTA measures: self-rated health, disability, depression, inactivity, chronic disease, obesity, sleep, access | **partial** | Uses the 2025 release at two levels: county (3,144 counties, 40 measures) and ZCTA (32,520 ZCTAs, 31 measures). `geo_context` holds `places_{GHLTH,DISABILITY,DEPRESSION,LPA,OBESITY,SLEEP,ACCESS2,...}_{crude,ageadj}` and `geo_context_zcta` holds the ZCTA values. | Census-tract and place levels were not ingested. Nothing downstream needs them, because the engine works at county/state level (ZCTA for context). Until this audit, no document recorded this choice. |
| E CDC Long COVID (Household Pulse) | respect source resolution; no downscaling | **done** | `geo_condition_burden__cdc_long_covid` holds only state- and national-resolution rows (15,912 rows: 10,404 state, 5,508 national). The 3,144 county rows of `geo_condition_features` carry `burden_inherited=True`, `burden_source_resolution="state"`. The tools and dashboard state "not county prevalence" (tested in `tests/test_dashboard.py`). | Long COVID has no within-state variation (this is a data limit, not a code gap). |
| E CMS Mapping Medicare Disparities | prevalence, utilization, hospitalization; `burden_evidence_level` A-D per condition | **done** | `geo_condition_burden__cms_mmd`: 351,936 rows (fibromyalgia B, migraine B, me_cfs C = the CCW "Fibromyalgia, Chronic Pain and Fatigue" proxy). `geo_context__cms_mmd`: 1,249,885 rows (prevalence, ED visit rate, hospitalization rate, costs). `geo_burden_definitions` has 42 rows; see `docs/BURDEN_DEFINITIONS.md`. Levels: A long_covid, lyme_disease; B fibromyalgia, migraine; C me_cfs, ptlds; D dysautonomia, pots, eds_hsd, mcas, gastroparesis, ibs, endometriosis, post_infectious_syndrome. | 8 of 14 conditions have no usable burden estimate (level D), and the tools return UNKNOWN for them. |
| E CDC disease-specific surveillance | Lyme/tickborne where applicable | **done** | `geo_condition_burden__cdc_lyme`: county case counts 2001-2023 (751,335 cases). The 8 legacy CT counties and 2 unmatched FIPS are reported in the DATA_AUDIT rather than force-mapped. | - |
| E CDC/ATSDR SVI | overall + socioeconomic, household, minority, housing/transport themes; at the right FIPS; not prevalence | **done** | `geo_vulnerability`: 3,144 counties, SVI 2022 (`rpl_themes`, `rpl_theme1`-`rpl_theme4`). `geo_vulnerability_puerto_rico` holds Puerto Rico. In scoring, SVI is its own component (`vulnerability_pct`) and is never used as burden. | - |
| E Census ACS | population, age, disability, insurance, poverty, income, transportation, broadband, rural/urban | **done** | `geo_context` (ACS 2020-2024 5-year): total_population, median_age, pct_age_18_64 / 65_plus, pct_with_disability, pct_uninsured, pct_below_poverty, median_household_income, pct_households_no_vehicle, pct_households_broadband (with MOEs). Rural/urban comes from USDA RUCC 2023 (`rucc_2023`, `is_metro`). | - |
| F NPPES | V2 file; NPI, type, name, taxonomy, specialty, address, city, state, ZIP, lat/lon; 13 specialty groups; NPI specialty is not treatment of the condition | **done** | `providers`: 1,456,346 rows (npi, entity_type, name, primary taxonomy, specialty_groups, address, city, state, zip5, lat/lon, county_fips, geocode_method), from NPPES September 2026 V.2. `provider_specialty_groups` has 30 groups, including all 13 SPEC groups (cardiology, electrophysiology, neurology, rheumatology, immunology_allergy, infectious_disease, pulmonology, gastroenterology, pain_medicine, pmr, primary_care, gynecology, vascular). Candidate reasons that rest on NPPES taxonomies carry the caveat "A taxonomy code does not mean the facility evaluates or treats ...". | Coordinates are ZCTA centroids, not street geocodes. The documents say so. |
| F HRSA | service-delivery sites: facility, location, type, metadata, geographic ids | **done** | `facilities__hrsa`: 19,266 sites (facility_name, health_center_type, site_type, location_setting, source lat/lon, county_fips, HPSA links). `hpsa_designations` and `mua_designations` come from HRSA data.hrsa.gov. | - |
| F `clinic_registry.parquet` | clinic registry | **done** | 321,536 facilities merged from NPPES organisations, HRSA sites, ClinicalTrials.gov sites and NIH organisations. Each has an `object_id = facility:<id>`, `match_methods` and `match_confidence_min`. See `docs/FACILITY_MATCHING.md`. | - |
| G `research_site_registry.parquet` | per facility: relevant trials, recruiting, conditions, technologies, recency, device/digital experience; NIH grants, active projects, funding, investigators, topics; separate from burden | **done** | 8,722 rows. Trial columns: n_trials, n_trials_site_recruiting, condition_ids_literal, measurement_classes_studied, latest_start_year, n_device_or_diagnostic_trials. NIH columns: n_nih_projects, n_nih_active_core_projects, nih_award_obligations_usd, n_nih_pis, nih_measurement_terms. `geo_condition_features` has per-geography `trials_in_geo_*` and `nih_*` counts. In scoring, research readiness is its own component (`research_readiness`). | NIH obligations are not expenditure; the tool caveat says so. |
| H Data.gov discovery | `search_us_open_data(query)` returning title, agency, description, access level, resource URLs, format, modified; route to the agency | **done** | `agents/datagov.py`. The offline call `search_us_open_data("long COVID")` returned 20 results with title, agency, description, access_level, resource_urls, resource_formats, modified, retrieve_from and the note "Data.gov does not host the records". `data_gov_search_results` and `data_gov_discovery_log` hold the discovery results. | In offline mode, a query that is not in the HTTP cache returns UNKNOWN. This is by design and documented in `docs/MCP.md`. |

## 2. Core data model (DuckDB + Parquet) and provenance fields

All 21 recommended tables exist in `data/processed/measure_it_public.duckdb` as `spec.<table>` views. Each view was
queried, and each has all 8 SPEC provenance columns: `source_name`, `source_record_id`, `source_version`,
`retrieved_at`, `source_geographic_resolution`, `evidence_type`, `evidence_level`, `provenance_notes`. Each also has
`data_layer`. `validate` confirms the same columns on 122/122 processed tables (the two boundary geoparquet files are
exempt).

| table | status | rows | note |
|---|---|---|---|
| participants | done | 22,345 | NHANES 19,931; Stanford Alavi 2,123, Uwakwe 126, Mishra 118; MapMECFS 47 (namespaced ids, no cross-dataset key) |
| participant_conditions | done | 46,421 | |
| participant_labs | done | 735,587 | NHANES |
| participant_wearable_features | done | 15,315 | NHANES 12,955 + Stanford 2,360 |
| participant_phenotype_embeddings | done | 4,152,767 | long format (participant x dimension) |
| conditions | done | 14 | view over `condition_registry` |
| condition_ontology_mappings | done | 298 | |
| condition_molecular_evidence | done | 70,295 | |
| measurements | done | 30 | view over `measurement_registry` |
| condition_measurement_evidence | done | 420 | 14 conditions x 30 measurements |
| measurement_regulatory_status | done | 102 | |
| geographies | done | 3,299 | 2024 counties + states (+ legacy CT flagged) |
| geo_condition_burden | done | 558,246 | |
| geo_context | done | 37,082 | |
| geo_vulnerability | done | 3,144 | |
| providers | done | 1,456,346 | |
| facilities | done | 329,019 | |
| clinical_trials | done | 7,910 | |
| trial_sites | done | 47,482 | |
| nih_projects | done | 9,342 | |
| deployment_candidates | done | 10 | the SPEC demo query's top 10, with `recommendation_json` and `recommended_next_step` |
| provenance fields on every table | done | - | `validate` provenance 122/122; `information_schema` check of the 21 views |

## 3. Phases 1-6

| phase | status | evidence | gaps |
|---|---|---|---|
| 1 Wearable phenotype | **done** | **NHANES**: models `lr`, `pen_lr` and `hgb` (`results/tables/nhanes_model_performance.csv`); 5-fold CV x 5 repeats at participant level; a 2011-12 -> 2013-14 temporal split; AUROC, AUPRC, sensitivity, specificity, calibration (`nhanes_calibration_bins.csv`) and participant-bootstrap CIs; 200 label permutations (`nhanes_permutation_null.csv`). Unsupervised: PCA, UMAP for visualisation only (`results/figures/drafts/nhanes_umap_clusters.png`), GMM by BIC and HDBSCAN. HDBSCAN labels everyone noise, and the report says "No discrete wearable subgroups". **Stanford**: logistic, L2-logistic and random forest; participant-grouped CV; 1,000 label permutations. Reports: `results/WEARABLE_NHANES_RESULTS.md`, `results/WEARABLE_STANFORD_RESULTS.md`. | The temporal/deep model was not attempted. The SPEC makes it conditional ("only add ... if it improves held-out performance"). No attempt or decision was recorded. |
| 2 Digital person | **done** | `digital_phenotype_datasets` (5 rows). Four rows are `digital_phenotype_vector` with `molecular_context_type = condition_level_molecular_enrichment_via_ontology`. MapMECFS is `omics_only_participant_linked` (47 participants; "NO wearable/physiology data is linkable in open data"). `digital_phenotype_blocks` has 391 feature rows; also `participant_phenotype_embeddings` and `docs/ANALYSIS_PLAN_DIGITAL_PERSON.md`. | No public dataset supports a wearable + linked-omics "Multimodal Person Representation". This is a data limit and is stated in the table. |
| 3 What could make the phenotype visible | **done** | `phase3_measurement_evidence` has 32 phenotype x measurement rows with observable_signals. It and `condition_measurement_evidence` (420 rows) both carry `phenotype_signal_strength`, `measurement_evidence_strength`, `technology_maturity`, `regulatory_visibility` and `deployment_complexity` as separate columns with `__` component columns. Neither has a composite column. `results/PHASE3_MEASUREMENT_EVIDENCE.md`. | - |
| 4 Geographic opportunity features | **done** | `deployment_opportunities`: 76,680 rows, 48 condition x measurement x level combinations. Component columns: raw `burden_*`, `vulnerability_pct`, `diagnostic_desert` (its raw desert components `desert_pct_*` are in `geo_condition_features`), `clinic_capacity`, `research_readiness`, `technology_saturation` (reported; not in the default composite), each with a percentile. 8 weight sets in `configs/scoring.yaml`; the default is equal (0.2 x 5). The Monte Carlo draws burden and Dirichlet weights 1,000 times (`rank_mc_p05` / `rank_mc_p95`). `burden_uncertainty_multiplier` widens proxy and inherited burden. See `docs/SCORING.md`. | - |
| 5 Clinic matching | **done** | `find_candidate_clinics("06073", "Long COVID", "wearable autonomic monitoring", radius_km=50)` returns six characteristics, each with a value and a separate rank: clinical_specialty_match, relevant_trial_history, NIH_research_activity, community_access, distance_to_target_population, technology_experience. Candidates are selected round-robin, with no combined score. The framing is verbatim: "This facility has characteristics suggesting it may be a viable implementation or study partner." Each ranked region in `rank_deployment_opportunities` carries `candidate_sites`; see `results/tables/demo_top10_candidate_sites.csv`. | For wearables, clinic capacity mostly measures primary-care density (Test 6b, section 9). |
| 6 Technology-to-clinic matching | **done** | One `rank_deployment_opportunities` record carries the whole chain: condition -> phenotype -> measurement -> technology (FDA context) -> geography -> candidate_sites. `recommended_next_step` opens with "Candidate deployment opportunity ..." (`results/example_deployment_recommendation.json`). | - |

## 4. Agentic MCP layer and output schema

The server is `src/measure_it/mcp/server.py` (FastMCP). It registers each tool from the facade
`src/measure_it/tools.py`; `.mcp.json` declares `uv run --frozen measure-it serve-mcp`. The signatures below are the
registered input schemas, read with `list_tools()`.

| SPEC tool | status | registered parameters (required in bold) | note |
|---|---|---|---|
| `search_condition(condition)` | done | **condition**, limit | |
| `normalize_condition(condition_or_code)` | done | **condition_or_code** | ICD-10-CM, MONDO/EFO/MeSH/HPO ids and aliases |
| `get_patient_phenotype_signature(condition)` | done | **condition**, top_n | |
| `get_molecular_context(condition)` | done | **condition**, top_n, include_broad | labelled condition-level molecular enrichment |
| `discover_candidate_measurements(condition, phenotype=None)` | done | **condition**, phenotype, top_n | |
| `get_measurement_evidence(measurement, condition)` | done | **measurement**, **condition**, max_object_ids | |
| `get_regulatory_context(technology)` | done | **technology**, max_records | |
| `get_condition_burden(condition, geography_level)` | done | **condition**, geography_level, geo_id, include_alternates, max_rows, order | `geography_level` has a default (`"state"`); the SPEC shows it without one. Callers that pass both arguments are unaffected. |
| `get_geographic_context(geography_id)` | done | **geography_id** | FIPS, `geo:`, `zcta:` or a name |
| `find_candidate_clinics(geography_id, condition, measurement)` | done | **geography_id**, **condition**, measurement, radius_km, max_sites | `measurement` is optional here; the SPEC shows it without a default. This is compatible. |
| `find_relevant_trials(condition, measurement=None)` | done | **condition**, measurement, literal_only, statuses, us_sites_only, max_trials | |
| `find_relevant_research_centers(condition)` | done | **condition**, measurement, top_n, include_expanded | |
| `rank_deployment_opportunities(condition, measurement, geography_level="county")` | done | **condition**, **measurement**, geography_level (default `county`), top_n, weight_set, include_small_population, max_sites_per_region, detail | |
| `search_us_open_data(query)` | done | **query**, rows | |
| `trace_evidence(object_id)` | done | **object_id**, max_rows | Person-level ids are refused. `trace_evidence("nhanes:62161")` returned UNKNOWN / refused. |
| Extras | done | `get_phenotype_measurement_evidence`, `get_measurable_biology`, `list_sources` | also 2 prompts and 4 resources (guardrails, data layers, object ids, tools) |
| No generated facts; UNKNOWN when no data | done | `tests/test_mcp_server.py::test_every_tool_ok_on_demo_and_unknown_on_nonsense`. The agent stops at UNKNOWN for "Zorblaxian drift syndrome" (docs/MCP.md section 7). The agent is template-only and uses no language model. | |

**Output schema.** Checked in `results/example_deployment_recommendation.json` (recommendations[0]); the result is **done**.
All SPEC keys are present: `condition`, `phenotype`, `measurement`, `technology` {`name`, `regulatory_context`,
`technology_evidence`}, `geography` {`name`, `fips`, `burden`, `burden_evidence_level`, `vulnerability`,
`diagnostic_desert`}, `candidate_sites` (20), `research_evidence` (5), `molecular_context`, `uncertainties`,
`provenance` and `recommended_next_step`. The record adds `candidate_site_query_status`, `technology.guardrail`,
and in `geography` the remaining components, rank, the Monte Carlo rank interval, P(top 10) and the ranks under each
weight set.

## 5. Dashboard (Streamlit, `dashboard/`; FastAPI in `src/measure_it/api/app.py`)

`tests/test_dashboard.py` runs every page in-process with Streamlit AppTest, using default and alternative inputs.
`tests/test_api.py` covers the API. Both passed in the full run. Screenshots are in `visuals/dashboard/`.
Page contents are described in `docs/DASHBOARD.md` section 3.

| page | SPEC contents | status | evidence | gaps |
|---|---|---|---|---|
| 1 National Opportunity Map | burden; diagnostic deserts; deployment opportunities; facilities; trial sites; layers toggled independently | **done** | `pages/1_National_Opportunity_Map.py`. In the screenshot (viewed): county choropleth, trial-site overlay, red outline on the top 10, always-visible burden-evidence box, ranked table. The choropleth radio picks Deployment opportunity / Diagnostic desert / Burden (with evidence level). Independent checkboxes: HRSA health-center sites, trial sites (+ active-only), NIH-funded organisations, specialist density. | The three area layers share one choropleth, so only one shows at a time; the point overlays toggle independently. The "facilities" overlay shows HRSA sites; `clinic_registry` facilities appear on pages 4-5, not on this map. |
| 2 Condition Explorer | input "Long COVID"; phenotype features, wearable abnormalities, molecular evidence, candidate measurements, technologies, burden sources | **done** | `pages/2_Condition_Explorer.py` has six sections, one per SPEC item. The test asserts the five dimensions as columns, "not patient multi-omics", and the Long COVID null result shown as a null. Screenshot `condition_explorer.png`. | - |
| 3 Measurement Explorer | input "wearable autonomic monitoring"; conditions studied, signals, public datasets, trials, FDA device classes, regions with unmet need | **done** | `pages/3_Measurement_Explorer.py`. The test asserts 10 unmet-need regions and the column "rank interval 5-95 (equal-weight MC)". Screenshot `measurement_explorer.png`. | - |
| 4 Geography Explorer | input "San Diego County, California"; burden, vulnerability, provider density, specialists, trials, NIH research, candidate technologies, nearby health centers | **done** | `pages/4_Geography_Explorer.py` (six sections). The test asserts "San Diego County, California" and "inherited by the county". Screenshot `geography_explorer.png`. | - |
| 5 Deployment Recommendation | the SPEC example query; ranked table, map, component scores, evidence, limitations | **done** | `pages/5_Deployment_Recommendation.py` shows the component percentiles, rank interval, P(top 10), a contribution chart, a map with candidate sites, per-region evidence, the unranked-burden list and a JSON download. Screenshot `deployment_recommendation.png`. | - |

## 6. CAPRIO / capillaroscopy adapter

| requirement | status | evidence |
|---|---|---|
| `MeasurementAdapter` with `preprocess`, `embed`, `phenotype_score`, `describe_measurement` | **done** | `src/measure_it/measurements/adapter.py`: an ABC with the four abstract methods and `MeasurementDescription` |
| `WearableAdapter` for the public alpha | **done** | `src/measure_it/wearables/adapter.py`: `WearableAdapter` (wrist accelerometry) + `WearableHeartRateAdapter` (Stanford resting HR). Outputs `participant_adapter_scores__nhanes` / `__stanford_covid` and `results/tables/wearable_adapter_validation.csv`. |
| stub `CapillaroscopyAdapter`, no private data | **done** | `src/measure_it/measurements/capillaroscopy_adapter.py`: `preprocess`, `embed` and `phenotype_score` raise `NotImplementedError`, naming what CAPRIO must supply. `describe_measurement()` returns status `stub` with `requires_private_data = True`. |
| documented image -> embedding -> capillary phenotype -> same pipeline | **done** | `docs/ADAPTERS.md` sections 3-6 ("What would NOT change") and Figure 10 (`fig10_adapter_modularity`) |
| downstream is modality-agnostic | **done** | `tests/test_adapters.py` runs a `SyntheticFakeAdapter` through the downstream functions unchanged and asserts that their source names no modality |
| capillaroscopy already flows through geography/clinic/scoring | **done** | `deployment_opportunities` has `nailfold_capillaroscopy` rows for all 6 condition(s)/sets (county + state). The agent answers "... nailfold capillaroscopy, where should we deploy it?" and reports the phenotype link as UNKNOWN (docs/MCP.md section 7). |

## 7. DATA_AUDIT per source and SOURCE_REGISTRY use

Each `data/raw/<source>/DATA_AUDIT.md` was checked for the 17 SPEC fields. The header table covers source, publishing
organisation, retrieval date, licence/access, update date, unit of observation, sample size, geography, person-level?,
geographic?, omics?, wearable? and true participant linkage?. The sections cover key variables, missingness, linkage
strategy and limitations. **27/27 audits contain all 17 fields**, and no header value is empty. The flags below are
the `SOURCE_REGISTRY.yaml` values.

| source | layer | status | person | geographic | omics | wearable | true linkage | fields | audit lines |
|---|---|---|---|---|---|---|---|---|---|
| `gwas_catalog` | condition_molecular | ingested | no | no | yes | no | no | 17/17 | 210 |
| `ncbi_geo_sra` | condition_molecular | ingested | no | no | yes | no | no | 17/17 | 209 |
| `open_targets` | condition_molecular | ingested | no | no | yes | no | no | 17/17 | 401 |
| `reactome` | condition_molecular | ingested | no | no | yes | no | no | 17/17 | 54 |
| `clinicaltrials_gov` | facility | ingested | no | yes | no | no | no | 17/17 | 406 |
| `hrsa_health_centers` | facility | ingested | no | yes | no | no | no | 17/17 | 108 |
| `nih_reporter` | facility | ingested | no | yes | no | no | no | 17/17 | 384 |
| `nppes` | facility | ingested | no | yes | no | no | no | 17/17 | 264 |
| `nucc_taxonomy` | facility | ingested | no | no | no | no | no | 17/17 | 59 |
| `openfda_device` | facility | ingested | no | no | no | no | no | 17/17 | 129 |
| `cdc_long_covid` | geographic | ingested | no | yes | no | no | no | 17/17 | 136 |
| `cdc_lyme` | geographic | ingested | no | yes | no | no | no | 17/17 | 238 |
| `cdc_places` | geographic | ingested | no | yes | no | no | no | 17/17 | 174 |
| `cdc_svi` | geographic | ingested | no | yes | no | no | no | 17/17 | 167 |
| `census_acs` | geographic | ingested | no | yes | no | no | no | 17/17 | 253 |
| `census_geography` | geographic | ingested | no | yes | no | no | no | 17/17 | 105 |
| `cms_mmd` | geographic | ingested | no | yes | no | no | no | 17/17 | 294 |
| `hrsa_hpsa` | geographic | ingested | no | yes | no | no | no | 17/17 | 156 |
| `usda_rucc` | geographic | ingested | no | yes | no | no | no | 17/17 | 112 |
| `data_gov_catalog` | metadata | ingested | no | no | no | no | no | 17/17 | 195 |
| `cms_icd10cm` | ontology | ingested | no | no | no | no | no | 17/17 | 85 |
| `ebi_ols4` | ontology | ingested | no | no | no | no | no | 17/17 | 134 |
| `hgnc` | ontology | ingested | no | no | no | no | no | 17/17 | 49 |
| `mondo_ontology` | ontology | ingested | no | no | no | no | no | 17/17 | 164 |
| `mapmecfs_nih_pi_mecfs` | person | partial | yes | no | yes | no | yes (omics-to-omics only; 0 wearable x omics) | 17/17 | 147 |
| `nhanes_2011_2014` | person | ingested | yes | no | no | yes | yes (SEQN across questionnaires, labs, accelerometry) | 17/17 | 535 |
| `stanford_longcovid_wearables` | person | ingested | yes | no | no | yes | no (wearable and labels share ids within each sub-dataset; nothing links across the three) | 17/17 | 252 |

**SOURCE_REGISTRY.yaml is read by the application (done).** It holds 27 sources, merged from the per-source
`registry_entry.yaml` files by pipeline step `registry_merge`. The readers found in the code, several also checked at
runtime:

| reader | what it uses the registry for |
|---|---|
| `mcp/trace.py` (`_registry`, `sources_for_table`) | maps each processed table to its sources; every `trace_evidence` result lists registry entries with version, retrieval date and the DATA_AUDIT path. Checked at runtime on the top opportunity id. |
| `tools.py` (`_provenance`, `_source_brief`, `list_sources`) | the `provenance.sources` block of every tool envelope; the `list_sources` tool |
| `api/app.py` | `GET /sources`, `GET /sources/{id}`, `/health` |
| `dashboard/app.py` | the home page's source table (via `list_sources`); the test asserts at least 20 sources |
| `agents/datagov.py` (`load_source_registry`) | searches Data.gov for every registry source |
| `figures/schematics.py` | Figure 2 is generated from the registry |
| `database.py` | the `_sources` table in DuckDB (27 rows) |
| `validate.py` | the data_audits check (27/27) |

## 8. Reproducibility stack

| item | status | evidence |
|---|---|---|
| Python | done | `requires-python >= 3.12` |
| uv | done | `uv.lock`, `uv sync --frozen`, `uv run measure-it ...` |
| DuckDB | done | `data/processed/measure_it_public.duckdb` (1,066 MB): 130 tables and views in `main` (including `_sources`), plus the `spec.*` views |
| Parquet | done | `data/processed/*.parquet` via `store.write_table` (+ `.meta.json`) |
| Polars or pandas | done | pandas (72 modules import it). Polars is a declared dependency, but nothing in `src/` imports it. |
| GeoPandas | done | 10 modules, among them `geography/crosswalk.py`, `geography/features.py`, the ClinicalTrials/NIH/HRSA/ACS ingestion, `scoring/recommend.py`, `api/geo.py` and `figures/data_maps.py`; 2024 county/state boundary geoparquet |
| scikit-learn | done | `wearables/models.py`, `wearables/stanford_models.py`, `wearables/unsupervised.py` |
| FastAPI | done | `src/measure_it/api/app.py`, `tests/test_api.py` |
| FastMCP | done | `src/measure_it/mcp/server.py`; stdio checked (section 0) |
| Plotly / Leaflet / MapLibre | done | Plotly choropleths in the dashboard; `results/maps/fig06-09_*.html` |
| pytest | done | 32 test files; 1,177 passed |
| cache API requests; no repeated federal hits | done | `data/_http_cache` (503 MB), `MANIFEST.json` per raw source, `--offline` mode with a socket guard (`docs/REPRODUCIBILITY.md` section 4) |
| timestamps and source versions | done | `retrieved_at` / `source_version` on every row; the registry; the MANIFEST sha256; determinism tables in `docs/REPRODUCIBILITY.md` section 5 |

## 9. Validation and anti-sycophancy tests

Every test has a result file and a stated outcome, and negative results are reported alongside positive ones.

| test | status | result files | stated outcome |
|---|---|---|---|
| 1 Wearable signal vs demographic baseline | **done** | `results/WEARABLE_NHANES_RESULTS.md` section 1, `results/tables/nhanes_increments.csv`; `results/WEARABLE_STANFORD_RESULTS.md` section A, `stanford_uwakwe_model_increments.csv`; Figure 4 | NHANES: adding wearable features to demographics raises AUROC for 5 of 6 targets (+0.019 to +0.032, 95% CIs > 0). The exception is the ME/CFS-like proxy, which is null. Wearable-only models score below demographics for 6/6 targets (point estimate; 95% CI below 0 for 5/6). Stanford Long COVID (Uwakwe, n = 126, 31 cases): **null**, wearable AUROC 0.504 vs demographics 0.702. |
| 2 Multimodal increment over clinical | **done** | same files | NHANES: the increment is small and robust only for functional limitation (+0.0069; Holm p <= 0.012). Fatigue and depression are unconfirmed; self-rated health, the ME/CFS-like proxy and mortality are null. Stanford: **null** (-0.024, 95% CI -0.076 to 0.027). |
| 3 Molecular coherence | **done** | `results/MOLECULAR_COHERENCE.md`, `results/tables/test3_*.csv` (source counts, gene/pathway agreement, controls, disagreements) | Gene-level agreement is above the cross-condition baseline for 5 of 6 source pairs by permutation; for one of the five the bootstrap CI includes 0. Pathway level: 4 of 6 pairs. The MapMECFS ME/CFS gene list is **not** the best-matching condition in any of its 4 source comparisons. |
| 4 Geography: burden vs population size | **done** | `results/GEOGRAPHY_RESULTS.md` section 4, `results/tables/test4_*.csv` | Prevalence and counts are both reported (e.g. Long COVID state prevalence rho with population -0.26; count 0.95). Pre-specified verdicts for the diagnostic desert: 10 of 12 series burden-driven, 2 mixed. The reviewer-added check found the desert top-25 over-represents bottom-population-quartile counties (9-19 vs 6.3 expected), so the desert partly measures small population and rurality. |
| 5 Deployment: does facility/research information change where we deploy? | **done** | `results/SCORING_RESULTS.md` section 5, `results/tables/test5_*.csv` | Yes, partly by construction. Burden-only vs full ranking: top-10 overlap 0/10, top-25 overlap 4/25 (Jaccard 0.087), Kendall tau-b 0.467. |
| 6 Negative controls (phenotype labels, geography labels, clinic locations) | **done** | phenotype labels: `nhanes_permutation_null.csv`, `stanford_uwakwe_permutation_nulls.csv`, `test6_scoring_phenotype_label_nulls.csv`. Geography labels: `test6_geography_*.csv`, `test6_scoring_label_shuffle*.csv`. Clinic locations: `test6_clinic_shuffle*.csv`, `test6_scoring_clinic_shuffle*.csv` | Burden-label shuffle: the Spearman of the composite with real burden drops from 0.657 to 0.098. Across-source geographic relationships: only the same-source positive control survives spatial-autocorrelation correction. **Clinic-location shuffle (negative result):** clinic_capacity barely moves (Spearman 0.965), because it measures facility density rather than measurement-specific capability. Research readiness does move (0.465). |
| 7 Ranking robustness | **done** | `results/SCORING_RESULTS.md` section 7, `results/tables/test7_*.csv`, `test7_unstable_regions.csv` | Kendall's W is 0.806 across the 8 named weight sets and 0.715 across 1,000 random weight vectors. 22 of the default top-25 meet the pre-specified instability rule. The median rank interval of the top-25 spans 319 ranks. None of the top-10 stays in the top-10 in at least half of the draws. |

## 10. Required final datasets

`validate` reports 11/11 present and non-empty. Row counts were read from the files:

| dataset | status | rows x columns |
|---|---|---|
| participant_wearable_features.parquet | done | 15,315 x 105 |
| participant_clinical_features.parquet | done | 19,931 x 172 (NHANES; Stanford demographics and comorbidities live in `participants` / `participant_conditions`) |
| condition_registry.parquet | done | 14 x 50 |
| condition_molecular_evidence.parquet | done | 70,295 x 187 |
| condition_measurement_evidence.parquet | done | 420 x 98 |
| measurement_registry.parquet | done | 30 x 127 |
| geo_condition_features.parquet | done | 46,074 x 112 |
| clinic_registry.parquet | done | 321,536 x 57 |
| research_site_registry.parquet | done | 8,722 x 66 |
| deployment_opportunities.parquet | done | 76,680 x 140 |
| measure_it_public.duckdb | done | `data/processed/measure_it_public.duckdb`, 1,066 MB |

## 11. Required figures

Each figure is in `results/figures/` as PNG (2161 px wide) and SVG, with a caption in `results/figures/CAPTIONS.md`.
Figures 1, 2 and 10 are generated by `figures.schematics` (Figure 2 from `SOURCE_REGISTRY.yaml`); Figures 3-9 by
`figures.data_figures`, from the processed tables. Figure 4 was viewed during this audit. The others were checked
for file, dimensions and caption, and by `tests/test_figures_data.py` / `tests/test_figures_schematics.py`.

| # | SPEC figure | file | status |
|---|---|---|---|
| 1 | System architecture | `fig01_system_architecture` | done |
| 2 | Public-data source architecture | `fig02_public_data_sources` | done |
| 3 | Wearable phenotype differences | `fig03_wearable_phenotype_differences` | done |
| 4 | Clinical-only vs clinical + wearable performance | `fig04_clinical_vs_wearable_models` (includes the Stanford null and the label-permutation nulls) | done |
| 5 | Condition -> molecular evidence -> measurable physiology | `fig05_molecular_measurable_physiology` | done |
| 6 | U.S. burden / unmet need, one condition | `fig06_long_covid_burden_unmet_need` (+ `results/maps/*.html`) | done |
| 7 | Provider / research infrastructure overlay | `fig07_infrastructure_overlay` | done |
| 8 | Deployment opportunity map | `fig08_deployment_opportunities` | done |
| 9 | Geography drill-down | `fig09_san_diego_drilldown` | done |
| 10 | Modularity: wearables <-> adapter <-> CAPRIO | `fig10_adapter_modularity` | done |

## 12. Primary demo story

| step | status | evidence | gaps |
|---|---|---|---|
| 1 Public person-level wearable data show an objectively measurable activity/autonomic phenotype | **partial** | **Activity:** NHANES shows measurable group differences for an ME/CFS-like proxy (14 of 26 wearable features pass FDR; e.g. interdaily stability -0.32 SD), and wearable data add over demographics for functional limitation, fatigue and fair/poor health. **Autonomic:** in acute COVID-19 (Stanford Mishra/Alavi), resting HR peaks at onset (mean z 1.03) and steps fall for longer. | No invisible-illness label shows an **autonomic** wearable phenotype. NHANES has no heart rate. On the only genuine Long COVID label (Uwakwe), whole-record resting-HR features are null. The acute single-day RHR elevation flag is not specific against placebo windows (62% vs 56%). The activity result is for a **proxy**, not diagnosed ME/CFS. The engine reports all of this, so the gap is in the evidence, not the code. |
| 2 Clinical and ontology data connect the phenotype to invisible illnesses | done | `phenotype_axes`, `condition_phenotypes` (Monarch/HPO), `measurement_phenotype_axis_map`, `phase3_measurement_evidence`, `normalize_condition` phenotype axes | |
| 3 Public molecular resources show pathways and candidate measurable biology | **partial** (updated after the report review) | `condition_molecular_evidence`, Test 3, `measurable_biology`, Figure 5, `get_molecular_context` (agent: systems per condition by support tier) | condition-level molecular enrichment only; after removing post-hoc-flagged drug-target sources no demo condition has a system with non-literature support, so nothing supports autonomic measurement; the system -> measurement-class link is a curated assumption |
| 4 ClinicalTrials.gov, RePORTER and openFDA identify technologies | done | `measurement_mentions`, `condition_measurement_evidence`, `measurement_regulatory_status`, `results/MEASUREMENT_DISCOVERY.md` | |
| 5 Federal population data show where burden and vulnerability concentrate | done | `geo_condition_features`, `get_condition_burden`, Figure 6 | Long COVID county values are inherited state values; ME/CFS burden is a level-C proxy; POTS and dysautonomia are level D. All are labelled. |
| 6 NPPES, HRSA, ClinicalTrials.gov and RePORTER show where infrastructure does or does not exist | done | `clinic_registry`, `research_site_registry`, Figure 7, dashboard page 1 overlays | |
| 7 Measurement deserts and candidate deployment sites | **partial** (updated after the report review) | `diagnostic_desert` (+ access-only variant), `find_candidate_clinics` with `absence_note` ("possible measurement/research desert"), `deployment_candidates` | the systematic nonmetro technology-experience measurement desert is not supported (`test6_clinic_shuffle.csv`); the county desert partly measures small population |
| 8 An agent answers "What should we measure, where, and who could deploy it?" | done | `uv run measure-it ask` (run in section 0); the MCP server for LLM clients, with the `deployment_question` prompt | |
| 9 CAPRIO presented as a future device-specific implementation | done | `docs/ADAPTERS.md`, Figure 10, the `CapillaroscopyAdapter` stub, `nailfold_capillaroscopy` rows in `deployment_opportunities`, `results/FINAL_REPORT.md` §19 | |

## 13. Success criteria

| # | criterion | status | evidence |
|---|---|---|---|
| 1 | At least one real public patient-level wearable dataset processed | **done** | NHANES 2011-2014 (12,955 wear-valid participants from PAXMIN/PAXDAY); three Stanford releases (2,360 participants with wearable features) |
| 2 | Wearable features linked to real clinical or disease information | **done** | NHANES: the same participant's questionnaires, labs, prescriptions and NCHS linked mortality, joined on SEQN. Uwakwe: the Long COVID label. |
| 3 | At least one invisible-illness example demonstrated | **done** | Long COVID / ME/CFS x wearable autonomic/activity monitoring runs end to end: `deployment_candidates`, `results/example_deployment_recommendation.json`, the agent answer. Its wearable evidence is a null (Long COVID) and a proxy (ME/CFS), stated as such. |
| 4 | Ontology connects clinical concepts to molecular evidence | **done** | `condition_registry` (MONDO/EFO ids) -> `condition_molecular_evidence.ontology_id` / `ontology_link_basis`; `get_molecular_context` |
| 5 | Candidate measurement technologies discovered automatically | **done** (within a curated vocabulary) | `measurements.discovery` mines ClinicalTrials.gov and RePORTER text for 26 curated measurement classes and counts objective use per condition. It runs unchanged for any registry condition, and an audit measured its precision per class. **Limit:** it cannot discover a technology type outside `configs/measurements.yaml`. The precision judgements came from a language-model reviewer and are disclosed as such (`results/MEASUREMENT_DISCOVERY.md`); they are QA flags, not score inputs. |
| 6 | U.S. public geographic data queried and normalised | **done** | 9 geographic sources on 2024 county/state FIPS, with a ZCTA crosswalk; `get_geographic_context("San Diego County, California")` |
| 7 | Providers, facilities and trial sites mapped | **done** | Figure 7, dashboard page 1 overlays, page 4 map; `providers`, `facilities__hrsa`, `trial_sites` with coordinates and geocode method |
| 8 | Reproducible geographic deployment opportunities | **done** | `uv run measure-it pipeline --offline` (49 steps at this audit; 74 since the 2026-09-24 integration); run-to-run determinism tables in `docs/REPRODUCIBILITY.md` section 5 |
| 9 | Every recommendation traceable to source data | **done** | `trace_evidence` over 20 namespaces (`docs/CONVENTIONS.md` section 7). `tests/test_trace.py::test_opportunity_lineage_reaches_features_burden_rows_and_raw_files` and `test_mcp_server.py::test_trace_round_trips_every_object_id_of_a_ranking` pass. This audit added Census ACS (the population denominator) to the opportunity trace. |
| 10 | CAPRIO insertable as a new MeasurementAdapter without restructuring | **done** | section 6; `tests/test_adapters.py` contract test |

## 14. Guardrails and final report

| requirement | status | evidence |
|---|---|---|
| Layer-1 participants never implied to be layer-3/4 people; no person-level geography | done | `validate` person_layer_no_geography 32/32; `trace_evidence` refuses participant ids; the agent states "those participants are not the people of any region, facility or trial" |
| Omics not linked to wearables are called condition-level molecular enrichment | done | `digital_phenotype_datasets.molecular_context_type`; the `get_molecular_context` label; MapMECFS linkage audit |
| Provenance on every derived table; every score exposes its components | done | sections 2 and 3 (Phases 3-4) |
| No LLM-generated numbers as data | done | the agent is deterministic and has no language model; tools return UNKNOWN for missing data. The one place a language model was used is the precision audit of text-mining snippets, which is disclosed and does not feed any score (success criterion 5). |
| Proxies labelled; no silent substitution | done | `burden_evidence_level` on every burden row; level C labelled "symptom/comorbidity proxy"; NHANES "ME/CFS-LIKE PROXY; not ME/CFS" |
| Product language | done | `validate` forbidden_language: 0 unquoted uses in 38 text files; `tests/test_dashboard.py::test_dashboard_source_text_uses_product_language` |
| `results/FINAL_REPORT.md` (executive summary, datasets, sample sizes, linkage audit, methods, scoring formula, validation, negative results, limitations, figures, MCP examples, example recommendation, CAPRIO) | **done after this audit** (was missing at the time) | Written after this audit. At the time the file did not exist. The material it needs does exist: `results/*_RESULTS.md`, `results/MOLECULAR_COHERENCE.md`, `results/MEASUREMENT_DISCOVERY.md`, `results/PHASE3_MEASUREMENT_EVIDENCE.md`, `docs/MAP_MECFS_LINKAGE_AUDIT.md`, `docs/SCORING.md`, `docs/MCP.md` sections 6-7, `docs/ADAPTERS.md`, `results/figures/CAPTIONS.md`, the dashboard screenshots and this checklist. |
| `README.md` | **done after this audit** (was missing at the time) | Written after this audit; at the time the file held one line: "measure-it-public (scaffold; README written at end of build)". |

## 15. Gaps (open issues), largest first

1. **`results/FINAL_REPORT.md` did not exist** at the time of this audit (written afterwards). It is a SPEC deliverable (section "FINAL REPORT").
2. **`README.md` was still the scaffold line** at the time of this audit (written afterwards).
3. **Demo step 1 is only partly supported by the data.** No invisible-illness label shows an autonomic wearable phenotype: NHANES has no heart rate, the Uwakwe Long COVID resting-HR result is null, and the acute-COVID RHR flag is not specific against placebo windows. The activity-pattern differences are for an ME/CFS-like proxy. The engine reports all of this correctly; the limit is in the public data.
4. **CDC PLACES census-tract and place levels were not ingested** (county and ZCTA only). Downstream needs neither.
5. **Optional items not done:** the NHANES raw 80-Hz mode, the PhysioNet/CovIdentify/Apple Watch benchmarks, and a temporal/deep model. The SPEC makes all three optional or conditional; no decision was recorded before this audit.
6. **MapMECFS per-participant physiology** (actigraphy, HRV, tilt, CPET) is behind a login, so wearable x omics linkage cannot be tested with open data.
7. **Ontology coverage:** HPO terms are empty for 8/14 conditions (no direct Monarch annotation, recorded per row). `ptlds` lacks an ICD-10-CM code, `mcas` an EFO id and `post_infectious_syndrome` a MeSH id. `icd10cm_blocks` is an empty column.
8. **Burden is weak for the demo conditions:** Long COVID counties inherit the state value, ME/CFS is a level-C proxy, and 8 of 14 conditions (including POTS and dysautonomia) are level D. All of this is labelled and propagated.
9. **Clinic capacity for wearables measures facility density.** The clinic-location shuffle barely changes it (Test 6b), and primary care is 86-92% of the relevant NPIs.
10. **Measurement discovery is bounded** by the 26-class curated vocabulary.
11. **Minor interface deviations:** in `get_condition_burden`, `geography_level` defaults to `"state"`; in `find_candidate_clinics`, `measurement` is optional. Both are compatible supersets of the SPEC signatures. On dashboard page 1, burden, desert and opportunity share one choropleth (the point overlays toggle independently), and the facilities overlay shows HRSA sites.
