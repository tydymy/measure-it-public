# Runbook: a private device or omics dataset (first: MAESTRO capillaroscopy) through the engine

Everything below runs on the data owner's machine; person-level data, phenotypes, query packs and outreach lists stay
in local-only folders (`data/raw/byod_*`, `results/byod/`, `outreach/`). Interfaces: docs/HARMONIZATION_CONTRACT.md.
Data contract and privacy checks: docs/BRING_YOUR_OWN_DATA.md.

## 0. Prepare the folder

MAESTRO has no EHR extract, so the bridge to public data is demographics (`age_band`, `sex`), self-reported history
(`self_report_history.csv`, mapped to ICD-10-CM), questionnaire scores (`survey_<instrument>.csv`: COMPASS-31,
DSQ-PEM, FSS-9, PROMIS ...), the Nightingale panel (`omics_nightingale.csv`; clinical-chemistry measures are mapped to
LOINC by `configs/harmonize_nightingale_loinc.yaml`), ring summaries (`device_<name>.csv`) and the capillaroscopy
features (`device_capillaroscopy.csv`). Labels: one column per cohort (`long_covid`, `ptlds`, `lyme_disease`) against
the healthy volunteers. Start from `templates/byod/` and declare a `subgroup_analysis` (pre-specified: which label,
what counts as device-positive, which domains may explain it).

## 1. Ingest and evaluate the measurement

```bash
uv run measure-it byod validate  <folder>   # refuses identifiers, geography, dates; lists every reason
uv run measure-it byod ingest    <folder>   # person-layer partitions + person_concepts__byod_<id>
uv run measure-it byod evaluate  <folder>   # locked plan: does the device separate cases from controls?
```

## 2. Stage B+C: who are the device-positive patients, in shared terms?

```bash
uv run measure-it byod subgroup <folder>    # locked plan: device-positive vs -negative inside the label
```

Writes `results/byod/<id>/computable_phenotype.json` (features with codes, coefficients, `ehr_evaluable`),
`SUBGROUP.md`, and the aggregate tables `computable_phenotypes` and `subgroup_strata_fractions`.

## 3. Stage D: which people resemble them?

```bash
uv run measure-it byod similar <id>
```

NHANES survey-weighted prevalence of the phenotype pattern (a pre-pandemic resemblance profile, never Long COVID
prevalence; suppressed when feature coverage < 0.25), plus query packs under `results/byod/<id>/similar/`: All of Us
(BigQuery) and N3C (Spark) SQL you run inside those enclaves (All of Us forbids publishing cells of 1-20), and a clinic
export (code lists + README) a partner health system runs in its own EHR. The SQL is parse-checked, not executed; the
All of Us state concept and Fitbit column names must be checked against the CDR data dictionary first.

## 4. Stage E: where are they?

```bash
uv run measure-it pipeline --only geo_subgroup_burden
```

County subgroup counts = base burden (BRFSS 2023 small-area model for Long COVID; CDC county Lyme cases; post-Lyme as
13.7% of Lyme cases, Aucott 2022) x the stratum fractions from step 2. The Long COVID county values are modelled:
within-state differences come from county composition and covariates (results/BRFSS_SAE_RESULTS.md).

## 5. Stage F: who could use the device there?

```bash
uv run python -m measure_it.outreach --condition long_covid --measurement nailfold_capillaroscopy \
    --top 10 --per-county 20 --prefer-subgroup --phenotype-id byod_<id>:<name>
uv run python -m measure_it.outreach --condition lyme_disease_or_ptlds --measurement nailfold_capillaroscopy \
    --top 10 --per-county 20 --specialty-before-rx
```

Clinicians are ranked by measurement-specific Medicare billing, Part D signals, experience-site affiliation, then the
measurement's implementer specialties. Nailfold capillaroscopy has no CPT/HCPCS code, so its billing evidence is
related or proxy codes (configs/measurement_hcpcs.yaml); `--specialty-before-rx` puts rheumatology first.

## 6. Let the performance record move the map

```bash
uv run measure-it byod deploy <folder>      # evidence_weighted ranking with the device's own record
uv run python -m measure_it.scoring.stage_layers   # S7-S9: SAE burden / activity capacity vs the primary ranking
```

`byod remove <id>` deletes every person-level and derived row of the dataset and restores the rankings.
