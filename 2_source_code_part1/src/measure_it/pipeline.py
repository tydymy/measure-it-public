"""The reproducible pipeline: an explicit, ordered DAG of every step that produces the processed tables and results.

    uv run measure-it pipeline --list                    # the ordered steps, their entry points and dependencies
    uv run measure-it pipeline --offline --workers 16    # everything, from data/raw + data/_http_cache only
    uv run measure-it pipeline --only scoring_report     # one step (its inputs must exist)
    uv run measure-it pipeline --from measurement_discovery

Each step is one module entry point (the same function its `python -m` CLI calls). Steps run in fresh
subprocesses (`python -m measure_it.pipeline _step <name>`), so no module-level cache leaks between steps and
`MEASURE_IT_OFFLINE=1` reaches every process. With --workers N the scheduler runs independent steps concurrently
while the sum of their declared worker counts stays <= N; a step's own worker count is min(its module default, N)
(Step.default_workers; the module defaults except ingest_stanford_wearables, 8 instead of 12, so that the two
person-level ingestions fit a 16-worker budget together; results do not depend on worker counts because every
parallel task carries its own seed or is deterministic per participant). Long-pole steps (priority > 0: the NHANES
and Stanford ingestion and model steps) start first when ready, and while one of them waits for capacity no
lower-priority step is started. Otherwise STEPS order breaks ties.

Every step's start/end/duration/status is appended to results/pipeline_runs/<timestamp>.jsonl, and its stdout/stderr
to results/pipeline_runs/<timestamp>/<step>.log (git-ignored). A failed step stops the scheduling of its dependants;
the run exits non-zero.
"""
from __future__ import annotations

import argparse
import importlib
import json
import os
import re
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path

from .config import PROCESSED, PROJECT_ROOT, RESULTS, utc_now_iso

RUNS_DIR = RESULTS / "pipeline_runs"
OFFLINE_ENV = "MEASURE_IT_OFFLINE"
WORKERS_ENV = "MEASURE_IT_WORKERS"

GEO_INGEST = ["ingest_svi", "ingest_census_acs", "ingest_rucc", "ingest_hrsa", "ingest_cdc_places",
              "ingest_cdc_long_covid", "ingest_cdc_lyme", "ingest_cms_mmd"]


@dataclass(frozen=True)
class Step:
    name: str
    target: str                       # "module:function"
    deps: tuple[str, ...] = ()
    kwargs: dict = field(default_factory=dict)
    offline_kwargs: dict = field(default_factory=dict)   # merged over kwargs when MEASURE_IT_OFFLINE=1
    workers_kwarg: str | None = None  # the entry point's worker-count parameter
    default_workers: int = 1          # the module default (the value behind the committed results)
    description: str = ""
    priority: int = 0                 # long-pole steps start first when several are ready (scheduling only)

    def workers(self, budget: int) -> int:
        return max(1, min(self.default_workers, budget))


def _s(name, target, deps=(), description="", **kw) -> Step:
    return Step(name=name, target=target, deps=tuple(deps), description=description, **kw)


# ------------------------------------------------------------------------------------------------------------------
# The DAG, in topological order (the order used for --list, --from and ties)
# ------------------------------------------------------------------------------------------------------------------
STEPS: list[Step] = [
    # geography base + ontology
    _s("geography_crosswalk", "measure_it.geography.crosswalk:build", (),
       "2024 county/state/ZCTA geographies, ZCTA centroids, ZCTA-county crosswalk, boundary geoparquets"),
    _s("ontology", "measure_it.ontology.build_registry:run", (),
       "condition_registry, condition_ontology_mappings, condition_phenotypes, phenotype_axes, ICD-10-CM codes"),
    # geographic ingestion
    _s("ingest_svi", "measure_it.ingestion.svi:run", ("geography_crosswalk",), "CDC/ATSDR SVI -> geo_vulnerability"),
    _s("ingest_census_acs", "measure_it.ingestion.census_acs:run", ("geography_crosswalk", "ingest_svi"),
       "ACS 5-year summary file -> geo_context__acs, geo_context_zcta__acs"),
    _s("ingest_rucc", "measure_it.ingestion.rucc:run", ("geography_crosswalk",), "USDA RUCC -> geo_context__rucc"),
    _s("ingest_hrsa", "measure_it.ingestion.hrsa:run", ("geography_crosswalk",),
       "HRSA health-center sites, HPSA/MUA designations -> facilities__hrsa, geo_context__hpsa"),
    _s("ingest_cdc_places", "measure_it.ingestion.cdc_places:run", ("geography_crosswalk",),
       "CDC PLACES county + ZCTA -> geo_context__cdc_places, geo_context_zcta__cdc_places"),
    _s("ingest_cdc_long_covid", "measure_it.ingestion.cdc_long_covid:run", ("geography_crosswalk",),
       "HPS long-COVID state estimates -> geo_condition_burden__cdc_long_covid"),
    _s("ingest_cdc_lyme", "measure_it.ingestion.cdc_lyme:run", ("geography_crosswalk",),
       "CDC Lyme county case counts -> geo_condition_burden__cdc_lyme"),
    _s("ingest_cms_mmd", "measure_it.ingestion.cms_mmd:run", ("geography_crosswalk", "ontology"),
       "CMS Mapping Medicare Disparities -> geo_condition_burden__cms_mmd, geo_context__cms_mmd"),
    # facility / research / regulatory ingestion
    _s("ingest_cms_physician_service", "measure_it.ingestion.cms_physician_service:run", (),
       "Medicare Physician & Other Practitioners by Provider and Service (configured HCPCS) -> provider_measurement_activity"),
    _s("ingest_cms_partd_prescriber", "measure_it.ingestion.cms_partd_prescriber:run", (),
       "Medicare Part D Prescribers by Provider and Drug (configured drugs) -> provider_condition_rx_signals"),
    _s("ingest_nppes", "measure_it.ingestion.nppes:run", ("geography_crosswalk",),
       "NPPES V.2 + NUCC taxonomy -> providers, provider_specialty_groups, provider_density_county",
       default_workers=4),
    _s("ingest_clinicaltrials", "measure_it.ingestion.clinicaltrials:run", ("geography_crosswalk",),
       "ClinicalTrials.gov API v2 pages -> clinical_trials, trial_*, ctgov_facility_summary"),
    _s("ingest_nih_reporter", "measure_it.ingestion.nih_reporter:run", ("geography_crosswalk",),
       "NIH RePORTER -> nih_projects, nih_project_conditions, nih_project_query_hits",
       offline_kwargs={"offline": True}),
    _s("ingest_openfda", "measure_it.ingestion.openfda:run", (),
       "openFDA device APIs -> fda_*, measurement_regulatory_status"),
    # person-level ingestion
    _s("ingest_nhanes", "measure_it.ingestion.nhanes:run", (),
       "NHANES 2011-2014 clinical tables + wrist-accelerometry features (measure_it.wearables.nhanes_features)",
       default_workers=8, priority=10),
    _s("ingest_nhanes_labs_extended", "measure_it.ingestion.nhanes_labs_extended:run", ("ingest_nhanes",),
       "NHANES 2011-2014 full laboratory catalogue -> participant_labs_extended__nhanes, "
       "participant_layer_membership__nhanes, layer overlap map (results/tables/nhanes_lab_*, nhanes_layer_overlap)"),
    _s("ingest_nhanes_2003_2006", "measure_it.ingestion.nhanes_2003_2006:run", (),
       "NHANES 2003-2006 demographics/questionnaire/exam/labs (incl. CRP)/meds/linked mortality -> *__nhanes0306"),
    _s("wearable_nhanes0306_features", "measure_it.wearables.nhanes0306_features:run", ("ingest_nhanes_2003_2006",),
       "NHANES 2003-2006 hip ActiGraph PAXRAW -> participant_wearable_features/daily__nhanes0306 (Troiano rules)",
       workers_kwarg="workers", default_workers=4),
    _s("ingest_stanford_wearables", "measure_it.ingestion.stanford_wearables:run", (),
       "Stanford Snyder Lab wearable releases -> participants/conditions/features/daily __stanford_covid",
       workers_kwarg="workers", default_workers=8, priority=10),   # module default 12; per-participant tasks
    _s("ingest_mapmecfs", "measure_it.ingestion.mapmecfs:run", ("ontology",),
       "NIH PI-ME/CFS public deposits -> participants/omics __mapmecfs, modality inventory/overlap, molecular evidence",
       offline_kwargs={"fetch": False}),
    # person-level device / lab releases (journal supplements; one dataset each, namespaced ids, never pooled)
    _s("ingest_muscle_me_charlton", "measure_it.ingestion.muscle_me_charlton:run", (),
       "MUSCLE-ME Source Data (Charlton 2026: Long COVID, ME/CFS, healthy controls) -> participants / "
       "participant_wearable_features (hip steps) / participant_exam_features (CPET, biopsy) __charlton_lc_mecfs_cpet_source",
       offline_kwargs={"skip_download": True}),
    _s("ingest_fm_thermography", "measure_it.ingestion.fm_thermography:run", (),
       "Fibromyalgia infrared thermography S1 Data (PLOS ONE 2021) -> participants / participant_device_features "
       "__fm_thermography", offline_kwargs={"skip_download": True}),
    _s("ingest_appelman_lc_pem", "measure_it.ingestion.appelman_lc_pem:run", (),
       "Appelman 2024 Long COVID PEM Source Data (metabolomics, VO2max): download, registry entry, audit"),
    _s("datagov", "measure_it.agents.datagov:run", (),
       "Data.gov catalog discovery for every SOURCE_REGISTRY source + demo queries (metadata layer); reads the "
       "SOURCE_REGISTRY.yaml of the previous registry_merge, so a new source is searched from the next run on"),
    # condition-level molecular evidence
    _s("omics_open_targets", "measure_it.omics.open_targets:run", ("ontology",),
       "Open Targets GraphQL -> condition_molecular_evidence__open_targets (+ coverage)"),
    _s("omics_gwas_catalog", "measure_it.omics.gwas_catalog:run", ("ontology",),
       "GWAS Catalog REST v2 -> condition_molecular_evidence__gwas_catalog (+ coverage)"),
    _s("omics_geo_sra", "measure_it.omics.geo_sra:run", ("ontology",),
       "NCBI GEO/SRA metadata -> geo_study_catalog, sra_study_catalog, condition_molecular_evidence__geo"),
    _s("omics_geo_cohorts", "measure_it.omics.geo_cohorts:run", ("omics_geo_sra", "ontology"),
       "GEO case/control cohorts (Long COVID, ME/CFS): screen, within-cohort classifiers, cross-cohort transfer -> "
       "phenotype_signatures__geo_cohorts, results/tables/geo_cohort_*", workers_kwarg="n_jobs", default_workers=16,
       priority=5),
    _s("omics_person_linked", "measure_it.omics.person_linked:run", ("ingest_mapmecfs",),
       "mapMECFS person-linked omics (5 layers, PI-ME/CFS vs HV) -> phenotype_signatures__mapmecfs_omics, "
       "results/LINKED_OMICS_RESULTS.md, results/tables/linked_omics_*", workers_kwarg="workers", default_workers=16,
       priority=5),
    _s("omics_union", "measure_it.omics.query:build_union",
       ("omics_open_targets", "omics_gwas_catalog", "omics_geo_sra", "ingest_mapmecfs"),
       "canonical condition_molecular_evidence / condition_molecular_coverage (union of partitions)"),
    _s("omics_reference_data", "measure_it.omics.reference_data:run", (),
       "Reactome release 97 + HGNC reference files, registry entries and audits"),
    _s("omics_coherence", "measure_it.omics.coherence:run", ("omics_union", "omics_reference_data", "ontology"),
       "Test 3 molecular coherence -> results/tables/test3_*"),
    _s("omics_graph", "measure_it.omics.graph:run", ("omics_coherence", "omics_geo_sra"),
       "measurable_biology + Figure 5 graph data + results/MOLECULAR_COHERENCE.md"),
    # wearable models
    _s("wearable_models", "measure_it.wearables.models:run", ("ingest_nhanes",),
       "NHANES Tests 1/2 models, permutation and shuffle controls -> results/tables/nhanes_*",
       workers_kwarg="n_jobs", default_workers=6, priority=10),
    _s("wearable_signatures", "measure_it.wearables.signatures:build", ("ingest_nhanes", "wearable_models"),
       "phenotype_signatures__nhanes"),
    _s("wearable_unsupervised", "measure_it.wearables.unsupervised:run", ("ingest_nhanes",),
       "PCA/GMM/HDBSCAN structure -> participant_cluster_assignments__nhanes"),
    _s("wearable_nhanes_report", "measure_it.wearables.nhanes_report:write",
       ("wearable_models", "wearable_signatures", "wearable_unsupervised"), "results/WEARABLE_NHANES_RESULTS.md"),
    _s("wearable_nhanes_layers", "measure_it.wearables.nhanes_layers:run", ("ingest_nhanes_labs_extended", "wearable_models"),
       "NHANES lab-layer increments + accelerometry over all lab layers -> results/NHANES_LAYERS_RESULTS.md, "
       "results/tables/nhanes_layers_*", workers_kwarg="n_jobs", default_workers=16, priority=10),
    _s("wearable_nhanes0306_models", "measure_it.wearables.nhanes0306_models:run",
       ("wearable_nhanes0306_features", "wearable_models"),
       "NHANES 2003-2006 replication of Tests 1/2/6 -> results/NHANES_2003_2006_RESULTS.md, results/tables/nhanes0306_*",
       workers_kwarg="n_jobs", default_workers=6, priority=10),
    _s("stanford_models", "measure_it.wearables.stanford_models:run", ("ingest_stanford_wearables",),
       "Stanford Uwakwe/Mishra/Alavi analyses -> phenotype_signatures__stanford, results/WEARABLE_STANFORD_RESULTS.md",
       workers_kwarg="jobs", default_workers=6, priority=10),
    # device analyses (pre-specified plans docs/ANALYSIS_PLAN_CHARLTON_LC_MECFS_CPET_SOURCE.md, _FM_THERMOGRAPHY.md)
    _s("wearable_muscle_me_steps", "measure_it.wearables.muscle_me_steps_analysis:run", ("ingest_muscle_me_charlton",),
       "MUSCLE-ME hip-accelerometer steps and CPET vs healthy controls -> phenotype_signatures__charlton_lc_mecfs_cpet_source, "
       "results/CHARLTON_LC_MECFS_CPET_SOURCE_RESULTS.md"),
    _s("wearable_fm_thermography", "measure_it.wearables.fm_thermography_analysis:run", ("ingest_fm_thermography",),
       "FM thermography primary/secondary models -> phenotype_signatures__fm_thermography, "
       "results/FM_THERMOGRAPHY_RESULTS.md", workers_kwarg="jobs", default_workers=8),
    # labs and other non-wearable measurements (plan docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md, locked before analysis)
    _s("labs_plan", "measure_it.labs.plan:lock_plan", (),
       "lock (or verify) docs/ANALYSIS_PLAN_LABS_VS_DIAGNOSIS.md against the plan text in measure_it.labs.plan"),
    _s("labs_vs_diagnosis", "measure_it.labs.labs_vs_diagnosis:run",
       ("labs_plan", "ingest_nhanes", "ingest_nhanes_labs_extended"),
       "NHANES 2011-2014 labs vs invisible-illness labels (positive controls, umbrella labels) -> "
       "phenotype_signatures__nhanes_labs, results/tables/labs_dx_nhanes_*", workers_kwarg="n_jobs",
       default_workers=16),
    _s("labs_appelman_metabolomics", "measure_it.labs.appelman_metabolomics:run", ("labs_plan", "ingest_appelman_lc_pem"),
       "Appelman 2024 plasma/muscle metabolomics, Long COVID vs healthy -> phenotype_signatures__appelman_metabolomics, "
       "results/tables/labs_dx_appelman_*", workers_kwarg="n_jobs", default_workers=16),
    _s("labs_report", "measure_it.labs.report:write", ("labs_vs_diagnosis", "labs_appelman_metabolomics"),
       "results/LABS_VS_DIAGNOSIS_RESULTS.md"),
    _s("labs_klein2023_mylc", "measure_it.labs.klein2023_mylc_ml_table:run", (),
       "MY-LC per-participant table (Klein 2023): ingest + cortisol / panel analyses -> participants / participant_labs / "
       "phenotype_signatures __klein2023_mylc_ml_table, results/KLEIN2023_MYLC_ML_TABLE_RESULTS.md",
       workers_kwarg="jobs", default_workers=8),
    _s("labs_endo_arg1_repod", "measure_it.labs.endo_arg1_repod:run", (),
       "Endometriosis serum arginase (RepOD PY1P9X): ingest + fit-free AUROC -> participants / participant_labs / "
       "phenotype_signatures __endo_arg1_repod, results/ENDO_ARG1_REPOD_RESULTS.md"),
    _s("labs_heds_hsd_olink_cinquina2026", "measure_it.labs.heds_hsd_olink_serum_cinquina2026:run", (),
       "hEDS/HSD serum Olink (Cinquina 2026): ingest + CV AUROC -> participants / participant_labs / phenotype_signatures "
       "__heds_hsd_olink_serum_cinquina2026, results/HEDS_HSD_OLINK_SERUM_CINQUINA2026_RESULTS.md",
       workers_kwarg="jobs", default_workers=8),
    _s("labs_published_evidence", "measure_it.labs.published_lab_evidence:run", (),
       "published lab/biomarker claims with verbatim-quote verification -> published_lab_evidence",
       offline_kwargs={"fetch": False}),
    _s("labs_dataset_discovery", "measure_it.labs.lab_dataset_discovery:run", (),
       "public lab-dataset candidates: download, inspect, rank -> results/tables/lab_dataset_candidates*.csv",
       offline_kwargs={"download": False}),
    _s("labs_aou_lab_concepts", "measure_it.labs.aou_lab_concepts:run", (),
       "All of Us Data Browser lab-concept participant counts -> results/tables/controlled_cohort_counts.csv (appended, "
       "deduplicated); offline it reads the saved raw JSON only", offline_kwargs={"fetch": False}),
    # digital person + adapters
    _s("digital_person", "measure_it.wearables.digital_person:run",
       ("ingest_nhanes", "ingest_stanford_wearables", "ingest_mapmecfs", "ontology", "ingest_nhanes_2003_2006",
        "wearable_nhanes0306_features", "ingest_muscle_me_charlton", "ingest_fm_thermography", "labs_klein2023_mylc",
        "labs_endo_arg1_repod", "labs_heds_hsd_olink_cinquina2026"),
       "canonical participant unions + Digital Phenotype Vector tables"),
    _s("wearable_adapter", "measure_it.wearables.adapter:run", ("ingest_nhanes", "ingest_stanford_wearables"),
       "WearableAdapter reference, validation and participant_adapter_scores", workers_kwarg="workers",
       default_workers=4),
    # measurement discovery and evidence
    _s("measurement_discovery_v1", "measure_it.measurements.discovery:run",
       ("ingest_clinicaltrials", "ingest_nih_reporter"), "round-1 (as-shipped) pattern mentions -> data/interim",
       kwargs={"version": "v1"}),
    _s("measurement_discovery", "measure_it.measurements.discovery:run",
       ("ingest_clinicaltrials", "ingest_nih_reporter"), "measurement_mentions (audited patterns)"),
    _s("measurement_precision_audit", "measure_it.pipeline:precision_audit_summarize",
       ("measurement_discovery_v1", "measurement_discovery"),
       "precision + PRO-filter audit summaries from the judgement files"),
    _s("measurement_phenotype_evidence", "measure_it.measurements.phenotype_evidence:run",
       ("measurement_discovery", "measurement_precision_audit", "ingest_openfda", "ontology", "ingest_nhanes",
        "ingest_stanford_wearables", "ingest_mapmecfs", "omics_geo_sra", "omics_graph", "wearable_models",
        "wearable_signatures", "stanford_models"),
       "measurement_registry, condition_measurement_evidence, Phase 3 evidence, results/PHASE3_MEASUREMENT_EVIDENCE.md"),
    _s("measurement_report", "measure_it.measurements.report:run", ("measurement_phenotype_evidence",),
       "results/MEASUREMENT_DISCOVERY.md"),
    _s("measurement_published_evidence", "measure_it.measurements.published_evidence:run", (),
       "published device / objective-test claims with verbatim-quote verification -> published_device_evidence",
       offline_kwargs={"fetch": False}),
    _s("measurement_device_dataset_discovery", "measure_it.measurements.device_dataset_discovery:run", (),
       "public device-dataset candidates: download, inspect, rank -> results/tables/device_dataset_candidates.csv"),
    _s("measurement_adapters", "measure_it.measurements.adapters:run",
       ("measurement_phenotype_evidence", "ingest_openfda", "ontology"), "measurement_adapter_registry"),
    _s("ingest_brfss", "measure_it.ingestion.brfss:run", ("geography_crosswalk", "ingest_census_acs"),
       "BRFSS 2022/2023 Long COVID items -> brfss_long_covid_direct_estimates; MMWR 2023 comparison"),
    _s("brfss_sae", "measure_it.geography.sae:run", ("ingest_brfss", "ingest_svi", "ingest_rucc", "ingest_cdc_places"),
       "MRP county Long COVID estimates -> geo_condition_burden__brfss_long_covid_sae, results/BRFSS_SAE_RESULTS.md"),
    # geography features, facilities, scoring
    _s("geography_features", "measure_it.geography.features:run",
       tuple(GEO_INGEST) + ("geography_crosswalk", "ingest_clinicaltrials", "ingest_nih_reporter", "ingest_nppes",
                            "brfss_sae"),
       "geo_condition_burden, geo_context, geo_condition_features, Tests 4/6 -> results/GEOGRAPHY_RESULTS.md"),
    _s("facilities_registry", "measure_it.facilities.registry:run",
       ("ingest_nppes", "ingest_hrsa", "ingest_clinicaltrials", "ingest_nih_reporter", "ingest_rucc",
        "geography_crosswalk"),
       "facilities, facility_trials, facility_nih_projects, clinic_registry, research_site_registry"),
    _s("facilities_shuffle_test", "measure_it.facilities.shuffle_test:run", ("facilities_registry", "ingest_rucc"),
       "Test 6 clinic-location shuffle of the matcher -> results/tables/test6_clinic_shuffle*",
       workers_kwarg="workers", default_workers=6),
    _s("facilities_examples", "measure_it.facilities.query:write_examples", ("facilities_registry", "ontology"),
       "facility example outputs + draft overlay map"),
    # MAESTRO stages (docs/HARMONIZATION_CONTRACT.md): shared person vocabulary (B+C), BRFSS small-area Long COVID
    # burden and device-subgroup burden (E)
    _s("harmonize_nhanes", "measure_it.harmonize.run:run",
       ("ingest_nhanes", "ingest_nhanes_labs_extended", "ingest_nhanes_2003_2006"),
       "NHANES 2011-2014 / 2003-2006 -> person_concepts__nhanes, person_concepts__nhanes0306 (ICD-10-CM, LOINC, RxNorm)"),
    _s("geo_subgroup_burden", "measure_it.geography.subgroup:run",
       ("brfss_sae", "ingest_cdc_lyme", "geography_features"),
       "device-defined subgroup burden from subgroup_strata_fractions (byod subgroup) -> geo_subgroup_burden"),
    _s("similar_people", "measure_it.similar.cli:run_all", ("harmonize_nhanes", "ingest_nhanes", "ingest_nhanes_2003_2006"),
       "phenotypes from byod subgroup -> similar_reference_prevalence/_profile/_distance (NHANES, aggregate); "
       "OMOP + clinic query packs under results/byod/<id>/similar/ (no-op without a phenotype)"),
    _s("facilities_activity", "measure_it.facilities.activity:build",
       ("ingest_cms_physician_service", "ingest_cms_partd_prescriber", "ingest_nppes", "geography_features",
        "facilities_registry", "measurement_discovery"),
       "measurement-specific clinician activity -> geo_measurement_capacity, provider_measurement_summary, "
       "measurement_experience_sites"),
    _s("scoring_metric_link", "measure_it.scoring.metric_link:run",
       ("facilities_registry", "geography_features", "wearable_muscle_me_steps", "stanford_models", "wearable_models",
        "labs_appelman_metabolomics", "labs_klein2023_mylc", "measurement_published_evidence"),
       "measurement_performance (+ _records): per condition x measurement performance record (own computations "
       "first, published claims flagged), sensitivity at 0.90 specificity, implementation reach -> "
       "results/tables/measurement_performance*.csv, measurement_reach.csv (docs/ANALYSIS_PLAN_METRIC_LINK.md)"),
    _s("scoring_report", "measure_it.scoring.report:run",
       ("scoring_metric_link", "geography_features", "facilities_registry", "facilities_activity", "ingest_nppes",
        "measurement_phenotype_evidence",
        "measurement_adapters", "omics_union", "omics_graph", "wearable_models", "wearable_signatures",
        "stanford_models", "ingest_openfda", "wearable_adapter", "ontology"),
       "deployment_opportunities, Tests 5-7, sensitivity, demo recommendation -> results/SCORING_RESULTS.md",
       workers_kwarg="workers", default_workers=6),
    _s("facilities_activity_report", "measure_it.facilities.activity:report", ("facilities_activity", "scoring_report"),
       "new vs old clinic capacity (Spearman, location shuffle) -> results/tables/stage_f_capacity_*.csv",
       ),
    _s("scoring_stage_layers", "measure_it.scoring.stage_layers:run", ("scoring_report", "facilities_activity"),
       "S7-S9: BRFSS SAE burden / activity capacity (primary since 2026-10-07) vs the legacy rules -> "
       "results/STAGE_LAYERS_SCORING.md, results/tables/stage_layers_*.csv"),
    # publication data figures (after the scoring steps; read processed + result tables, write results/figures,
    # results/maps and the Figure 3-9 block of results/figures/CAPTIONS.md)
    _s("figures_data", "measure_it.figures.data_figures:run",
       ("scoring_report", "omics_graph", "wearable_models", "wearable_signatures", "stanford_models",
        "geography_features", "facilities_registry", "ingest_clinicaltrials", "ingest_hrsa", "ingest_nppes",
        "ontology"),
       "SPEC Figures 3-9 -> results/figures/fig03-09_* (PNG + SVG), results/maps/*.html, CAPTIONS.md block"),
    # Test 8: temporal holdout of the deployment hotspots (docs/ANALYSIS_PLAN_TEMPORAL_HOLDOUT.md)
    _s("scoring_temporal_holdout", "measure_it.scoring.temporal_holdout:run",
       ("scoring_report", "scoring_metric_link", "geography_features", "facilities_registry", "ingest_clinicaltrials",
        "ingest_nih_reporter", "ingest_nppes", "ingest_hrsa", "ingest_cms_mmd"),
       "temporal holdout (Test 8): rankings rebuilt from pre-cutoff inputs by the unchanged engine vs post-cutoff "
       "research activity -> temporal_holdout_rankings, results/tables/temporal_holdout_*.csv, "
       "results/TEMPORAL_HOLDOUT_RESULTS.md"),
    # assembly
    _s("canonical_unions", "measure_it.pipeline:canonical_unions",
       ("scoring_report", "scoring_temporal_holdout", "digital_person", "wearable_adapter", "omics_union", "facilities_registry",
        "geography_features", "wearable_unsupervised", "wearable_signatures", "stanford_models", "ingest_mapmecfs",
        "ingest_cms_mmd", "ingest_census_acs", "ingest_cdc_places", "ingest_nhanes", "ingest_nhanes_labs_extended",
        "wearable_nhanes0306_features", "omics_geo_cohorts", "omics_person_linked", "wearable_muscle_me_steps",
        "wearable_fm_thermography", "labs_vs_diagnosis", "labs_appelman_metabolomics", "labs_klein2023_mylc",
        "labs_endo_arg1_repod", "labs_heds_hsd_olink_cinquina2026", "labs_published_evidence",
        "measurement_published_evidence", "harmonize_nhanes", "brfss_sae", "geo_subgroup_burden", "similar_people", "facilities_activity"),
       "canonical <table> for every <table>__<source> family not built by its own module"),
    _s("registry_merge", "measure_it.registry:merge_registry_fragments",
       ("canonical_unions", "omics_reference_data", "datagov", "ingest_appelman_lc_pem"),
       "data/raw/*/registry_entry.yaml -> SOURCE_REGISTRY.yaml"),
    # schematics (Figures 1, 2, 10): Figure 2 is drawn from the SOURCE_REGISTRY.yaml that registry_merge writes
    _s("figures_schematics", "measure_it.figures.schematics:run",
       ("registry_merge", "wearable_adapter", "measurement_adapters"),
       "SPEC Figures 1, 2, 10 (schematics from SOURCE_REGISTRY.yaml and the adapter code) -> results/figures"),
    _s("build_db", "measure_it.database:build_db", ("canonical_unions", "registry_merge"),
       "data/processed/measure_it_public.duckdb (tables + SPEC views + _sources)"),
    _s("validate", "measure_it.validate:run_checks", ("build_db",),
       "provenance, required datasets, person-layer geography, object ids, product language, audits"),
]
STEP_BY_NAME = {s.name: s for s in STEPS}


def _check_dag() -> None:
    seen: set[str] = set()
    for s in STEPS:
        for d in s.deps:
            if d not in STEP_BY_NAME:
                raise ValueError(f"step {s.name}: unknown dependency {d}")
            if d not in seen:
                raise ValueError(f"step {s.name}: dependency {d} is listed after it (STEPS must be topological)")
        seen.add(s.name)


_check_dag()


# ------------------------------------------------------------------------------------------------------------------
# step helpers that are not a single module entry point
# ------------------------------------------------------------------------------------------------------------------
def precision_audit_summarize() -> None:
    """`python -m measure_it.measurements.precision_audit summarize` (the judgement files are inputs)."""
    from .config import TABLES
    from .measurements import precision_audit as PA
    TABLES.mkdir(parents=True, exist_ok=True)
    PA.summarize()
    PA.summarize_pro()


# Canonical tables that a module builds itself (not a plain union of its partitions); the union step skips them.
OWNED_CANONICALS = {
    "research_site_registry": "facilities_registry (per resolved facility; never union the ClinicalTrials.gov history)",
    "facilities": "facilities_registry (entity-resolved NPPES/HRSA/ClinicalTrials.gov/NIH facilities)",
    "geo_condition_burden": "geography_features (harmonised long burden table incl. PLACES proxies)",
    "geo_context": "geography_features (one wide row per state/county/ZCTA)",
    "condition_molecular_evidence": "omics_union (store.union_partitions)",
    "condition_molecular_coverage": "omics_union (store.union_partitions)",
    "participant_adapter_scores": "wearable_adapter (store.union_partitions)",
    **{t: "digital_person (union with dataset_id, source_partition and harmonised keys)" for t in (
        "participants", "participant_conditions", "participant_labs", "participant_medications",
        "participant_wearable_features", "participant_wearable_daily", "participant_clinical_features",
        "participant_omics_linked")},
}


def partition_families() -> dict[str, list[str]]:
    fam: dict[str, list[str]] = {}
    for p in sorted(PROCESSED.glob("*__*.parquet")):
        table, src = p.stem.split("__", 1)
        fam.setdefault(table, []).append(src)
    return fam


def _common_dtype(dtypes: list) -> str:
    """Target dtype for one column whose dtype differs across partitions: nullable Int64 when every partition holds
    integers (e.g. int64 vs Int64), Float64 when every partition is numeric, bool -> boolean, else string."""
    import pandas.api.types as pt

    if all(pt.is_integer_dtype(d) for d in dtypes):
        return "Int64"
    if all(pt.is_numeric_dtype(d) and not pt.is_bool_dtype(d) for d in dtypes):
        return "Float64"
    if all(pt.is_bool_dtype(d) for d in dtypes):
        return "boolean"
    return "string"


def canonical_unions() -> list[dict]:
    """Union every <table>__<source> family into <table> unless a module owns the canonical table.

    Columns whose dtype differs across partitions are cast to a common dtype (integer kinds -> Int64, numeric ->
    Float64, bool kinds -> boolean, anything else -> string; recorded); the result is written with
    store.write_table (provenance enforced)."""
    import pandas as pd

    from .store import partitions, write_table

    out = []
    for table, srcs in partition_families().items():
        if table in OWNED_CANONICALS:
            out.append({"table": table, "partitions": srcs, "action": "skipped", "reason": OWNED_CANONICALS[table]})
            continue
        frames = [pd.read_parquet(p) for p in partitions(table)]
        casts = []
        cols = sorted({c for f in frames for c in f.columns})
        for c in cols:
            dtypes = [f[c].dtype for f in frames if c in f.columns]
            if len({str(d) for d in dtypes}) > 1:
                target = _common_dtype(dtypes)
                casts.append(f"{c}->{target}")
                for f in frames:
                    if c in f.columns:
                        f[c] = f[c].astype(target)
        df = pd.concat(frames, ignore_index=True, sort=False)
        write_table(df, table, producer="measure_it.pipeline.canonical_unions",
                    description=f"union of {[f'{table}__{s}' for s in srcs]}; columns cast to a common dtype "
                                f"for dtype conflicts: {casts or 'none'}")
        out.append({"table": table, "partitions": srcs, "action": "unioned", "rows": int(len(df)),
                    "dtype_casts": casts})
        print(f"[unions] {table}: {len(df):,} rows from {srcs}" + (f" (cast {casts})" if casts else ""), flush=True)
    for r in out:
        if r["action"] == "skipped":
            print(f"[unions] {r['table']}: skipped ({r['reason']})", flush=True)
    return out


# ------------------------------------------------------------------------------------------------------------------
# selection and scheduling
# ------------------------------------------------------------------------------------------------------------------
def select_steps(only: list[str] | None = None, from_step: str | None = None) -> list[Step]:
    names = [s.name for s in STEPS]
    for n in (only or []) + ([from_step] if from_step else []):
        if n not in STEP_BY_NAME:
            raise SystemExit(f"unknown step {n!r}; run `measure-it pipeline --list`")
    if only:
        return [s for s in STEPS if s.name in set(only)]
    if from_step:
        return STEPS[names.index(from_step):]
    return list(STEPS)


def describe(steps: list[Step] | None = None) -> str:
    steps = steps or STEPS
    lines = [f"{'#':>3}  {'step':32s} {'entry point':58s} deps"]
    for i, s in enumerate(steps, 1):
        extra = []
        if s.kwargs:
            extra.append(f"kwargs={s.kwargs}")
        if s.offline_kwargs:
            extra.append(f"offline={s.offline_kwargs}")
        if s.workers_kwarg or s.default_workers > 1:
            extra.append(f"workers<={s.default_workers}")
        lines.append(f"{i:>3}  {s.name:32s} {s.target:58s} {', '.join(s.deps) or '-'}"
                     + (f"  [{'; '.join(extra)}]" if extra else ""))
        if s.description:
            lines.append(f"{'':5s}{s.description}")
    return "\n".join(lines)


def _log_event(path: Path, rec: dict) -> None:
    with open(path, "a") as fh:
        fh.write(json.dumps(rec, default=str) + "\n")


def _rel(path) -> str:
    """A path relative to the checkout, so run logs carry no home directory (and no user name)."""
    try:
        return Path(path).resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except (ValueError, OSError):
        return str(path)


def _argv_rel(argv: list[str]) -> list[str]:
    root = str(PROJECT_ROOT.resolve())
    return [_rel(a) if str(a).startswith(root) else a for a in argv]


ERROR_LINE_RE = re.compile(r"^(?:[A-Za-z_][\w.]*(?:Error|Exception|Miss|Exit|Interrupt)\b|[\w.]*Error:|\[step\] .*FAILED)")


def last_error_line(log_file: Path, max_lines: int = 400) -> str:
    """The most informative line of a failed step's log: its last exception line (e.g. 'OfflineCacheMiss: raw file
    ... needs the network'), else its last non-empty line."""
    try:
        lines = [ln.rstrip() for ln in Path(log_file).read_text(errors="replace").splitlines()[-max_lines:]]
    except OSError:
        return ""
    lines = [ln for ln in lines if ln.strip()]
    exc = [ln for ln in lines if ERROR_LINE_RE.match(ln.strip()) and not ln.strip().startswith("[step]")]
    line = (exc or [ln for ln in lines if not ln.strip().startswith("[step]")] or lines or [""])[-1].strip()
    return line[:300]


def failure_hint(err: str, offline: bool) -> str:
    if "OfflineCacheMiss" in err or (offline and "needs the network" in err):
        return ("offline mode answers only from data/raw and data/_http_cache; on a fresh clone run the pipeline "
                "without --offline first (README 'Reproducing from scratch')")
    if "not found; run the pipeline step" in err:
        return "an input table is missing: run the step that produces it (or the whole pipeline) first"
    return ""


def _reap(proc: subprocess.Popen) -> tuple[int | None, float | None]:
    """(return code, peak RSS in MB) of a finished step, or (None, None) while it runs. The peak RSS is that of the
    largest single process of the step (the step or one of its workers; os.wait4 rusage), not the sum over workers."""
    if not hasattr(os, "wait4"):
        return proc.poll(), None
    try:
        pid, wstatus, ru = os.wait4(proc.pid, os.WNOHANG)
    except ChildProcessError:
        return proc.poll(), None
    if pid == 0:
        return None, None
    rc = os.waitstatus_to_exitcode(wstatus)
    proc.returncode = rc
    kb = ru.ru_maxrss / (1024.0 if sys.platform == "darwin" else 1.0)   # bytes on macOS, kilobytes on Linux
    return rc, round(kb / 1024.0, 1)


def run(only: list[str] | None = None, from_step: str | None = None, offline: bool = False, workers: int = 1,
        run_id: str | None = None, echo=print) -> int:
    steps = select_steps(only, from_step)
    selected = {s.name for s in steps}
    run_id = run_id or time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = RUNS_DIR / f"{run_id}.jsonl"
    step_log_dir = RUNS_DIR / run_id
    step_log_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env[WORKERS_ENV] = str(workers)
    if offline:
        env[OFFLINE_ENV] = "1"
    else:
        env.pop(OFFLINE_ENV, None)
    env.setdefault("PYTHONUNBUFFERED", "1")
    t_run = time.time()
    _log_event(log_path, {"event": "run_start", "run_id": run_id, "at": utc_now_iso(), "offline": offline,
                          "workers": workers, "steps": [s.name for s in steps], "python": sys.version.split()[0],
                          "argv": _argv_rel(sys.argv), "cwd": "."})   # paths relative to the checkout
    echo(f"[pipeline] run {run_id}: {len(steps)} step(s), workers={workers}, offline={offline}; log {_rel(log_path)}")
    pending = list(steps)
    running: dict[str, tuple[subprocess.Popen, float, int, str, object]] = {}
    done: dict[str, str] = {}   # name -> ok | failed | skipped
    budget = max(1, workers)

    def deps_ok(s: Step) -> str | None:
        for d in s.deps:
            if d in selected:
                st = done.get(d)
                if st is None:
                    return None          # still waiting
                if st != "ok":
                    return "blocked"
        return "ready"

    while pending or running:
        # start whatever is ready and fits the worker budget (STEPS order breaks ties)
        used = sum(r[2] for r in running.values())
        order = {st.name: i for i, st in enumerate(STEPS)}
        reserved = False   # a ready long-pole step is waiting for capacity: start nothing of lower priority
        for s in sorted(pending, key=lambda x: (-x.priority, order[x.name])):
            if reserved and s.priority == 0:
                break
            st = deps_ok(s)
            if st == "blocked":
                pending.remove(s)
                done[s.name] = "skipped"
                _log_event(log_path, {"event": "step_skipped", "step": s.name, "at": utc_now_iso(),
                                      "reason": "a dependency failed or was skipped"})
                echo(f"[pipeline] SKIP {s.name} (dependency failed)")
                continue
            if st != "ready":
                continue
            w = s.workers(budget)
            if running and used + w > budget:
                reserved = reserved or s.priority > 0
                continue
            logf = open(step_log_dir / f"{s.name}.log", "w")
            cmd = [sys.executable, "-m", "measure_it.pipeline", "_step", s.name, "--workers", str(budget)]
            proc = subprocess.Popen(cmd, cwd=str(PROJECT_ROOT), env=env, stdout=logf, stderr=subprocess.STDOUT)
            started = utc_now_iso()
            running[s.name] = (proc, time.time(), w, started, logf)
            pending.remove(s)
            used += w
            _log_event(log_path, {"event": "step_start", "step": s.name, "at": started, "pid": proc.pid,
                                  "target": s.target, "workers": w, "log": _rel(step_log_dir / f"{s.name}.log")})
            echo(f"[pipeline] START {s.name} (pid {proc.pid}, workers {w})")
        # reap finished steps
        finished = False
        for name, (proc, t0, w, started, logf) in list(running.items()):
            rc, max_rss_mb = _reap(proc)
            if rc is None:
                continue
            logf.close()
            dur = time.time() - t0
            status = "ok" if rc == 0 else "failed"
            done[name] = status
            del running[name]
            finished = True
            rec = {"event": "step_end", "step": name, "start": started, "end": utc_now_iso(),
                   "duration_s": round(dur, 1), "status": status, "returncode": rc, "workers": w,
                   "max_rss_mb": max_rss_mb}
            if rc == 0:
                _log_event(log_path, rec)
                echo(f"[pipeline] DONE {name} in {dur:,.0f}s (rc {rc})")
            else:
                step_log = step_log_dir / f"{name}.log"
                err = last_error_line(step_log)
                hint = failure_hint(err, offline)
                _log_event(log_path, {**rec, "error": err, "log": _rel(step_log)})
                echo(f"[pipeline] FAIL {name} in {dur:,.0f}s (rc {rc}): {err or 'no error line in the log'}")
                echo(f"[pipeline]      log: {_rel(step_log)}" + (f"; hint: {hint}" if hint else ""))
        if not finished:
            time.sleep(1.0)
    failed = [n for n, st in done.items() if st != "ok"]
    _log_event(log_path, {"event": "run_end", "run_id": run_id, "at": utc_now_iso(),
                          "duration_s": round(time.time() - t_run, 1), "status": "ok" if not failed else "failed",
                          "failed_or_skipped": failed})
    echo(f"[pipeline] run {run_id} {'OK' if not failed else 'FAILED: ' + ', '.join(failed)} "
         f"in {time.time() - t_run:,.0f}s")
    return 0 if not failed else 1


# ------------------------------------------------------------------------------------------------------------------
# single-step execution (inside the subprocess)
# ------------------------------------------------------------------------------------------------------------------
def _exit(code: int) -> None:
    """Flush and leave without interpreter teardown (see omics.query.exit_cleanly: native-library teardown on this
    host sometimes aborts with exit 134 after all work is written)."""
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(code)


def run_step(name: str, budget: int = 1) -> None:
    s = STEP_BY_NAME[name]
    from .http import offline
    kwargs = dict(s.kwargs)
    if offline():
        kwargs.update(s.offline_kwargs)
    if s.workers_kwarg:
        kwargs[s.workers_kwarg] = s.workers(budget)
    mod, fn = s.target.split(":")
    func = getattr(importlib.import_module(mod), fn)
    print(f"[step] {name}: {s.target}({', '.join(f'{k}={v!r}' for k, v in kwargs.items())}) "
          f"offline={offline()} at {utc_now_iso()}", flush=True)
    t0 = time.time()
    try:
        out = func(**kwargs)
        if name == "validate" and isinstance(out, dict) and not out.get("ok", False):
            print(f"[step] {name}: validation FAILED", flush=True)
            _exit(3)
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
        print(f"[step] {name}: SystemExit({e.code}) after {time.time() - t0:,.1f}s", flush=True)
        _exit(code)
    except BaseException:  # noqa: BLE001 - report every failure with its traceback, then exit non-zero
        traceback.print_exc()
        print(f"[step] {name}: FAILED after {time.time() - t0:,.1f}s", flush=True)
        _exit(1)
    print(f"[step] {name}: ok in {time.time() - t0:,.1f}s", flush=True)
    _exit(0)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m measure_it.pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)
    st = sub.add_parser("_step")
    st.add_argument("name")
    st.add_argument("--workers", type=int, default=int(os.environ.get(WORKERS_ENV, "1")))
    a = ap.parse_args(argv)
    if a.cmd == "_step":
        run_step(a.name, a.workers)


if __name__ == "__main__":
    main()
