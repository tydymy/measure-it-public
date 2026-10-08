"""Shared person-level vocabulary (docs/HARMONIZATION_CONTRACT.md): every dataset -> `person_concepts__<dataset_id>`.

    uv run python -m measure_it.harmonize.run            # NHANES 2011-2014 + 2003-2006 -> person_concepts__nhanes*
    uv run python -m measure_it.harmonize.audit          # docs/HARMONIZATION_AUDIT.md (measured coverage)
    measure-it byod ingest <dir>                         # a user dataset -> person_concepts__byod_<id>

Modules
    schema       the person_concepts columns, contract checks, writer, age bands, EHR-evaluability per vocabulary
    nhanes       NHANES processed tables -> ICD10CM / LOINC / RXNORM / SURVEY / DEMOG rows (curated configs, RxNav)
    byod         a user dataset's partitions -> person_concepts (+ Nightingale NMR -> LOINC, self-report -> ICD-10-CM)
    rxnorm       drug names -> RxNorm ingredients (RxNav REST, cached)
    self_report  self-reported history -> ICD-10-CM; PRO instrument codes; survey crosswalk to reference items
    features     feature_matrix (wide participant x feature), icd10_category, feature_catalog
    phenotype    score_phenotype / phenotype_coverage / validate_phenotype for computable phenotypes
    audit        measured coverage and code verification -> docs/HARMONIZATION_AUDIT.md
Configs (curated): configs/harmonize_nhanes_conditions.yaml, harmonize_nhanes_loinc.yaml,
harmonize_self_report_icd10.yaml, harmonize_nightingale_loinc.yaml, harmonize_survey_crosswalk.yaml.
"""
