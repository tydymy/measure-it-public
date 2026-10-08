# Launch agent contract: any dataset -> harmonized layers -> launch sites

Status: draft 2026-10-07. Binding for the three build parts below; change it only by
editing this file and saying so in the part's report.

Goal: one stand-alone command takes a dataset in (almost) any format, public or private (a device study such as
MAESTRO), harmonizes it to this engine's schema (`measure-it byod` folder contract + `person_concepts`,
docs/HARMONIZATION_CONTRACT.md) AND to the shape of the All of Us Curated Data Repository (OMOP CDM v5.4 plus the
All of Us survey and Fitbit tables), connects it to the phenotype, geography and clinician layers, and returns a
ranked, traceable list of candidate launch sites for the measurement device.

```
measure-it launch run <input path> --condition long_covid --measurement nailfold_capillaroscopy \
    [--out launch_runs/<id>] [--llm none|claude] [--approve] [--public] [--keep]
measure-it launch profile|map <input path> ...;  measure-it launch apply|report <run folder>
```

| part | package | owns |
|---|---|---|
| A | `measure_it.agent.read`, `measure_it.agent.profile`, `measure_it.agent.mapper`, `measure_it.agent.llm` | readers, profile, mapping proposals (rules, fuzzy, optional Claude) |
| B | `measure_it.agent.apply`, `measure_it.agent.omop`, `measure_it.agent.omop_run` | BYOD folder writer, de-identification transforms, OMOP / All of Us export, running the stage-D query packs locally against that export |
| C | `measure_it.agent.orchestrate`, `measure_it.agent.launch`, `measure_it.agent.report`, `measure_it.agent.cli` | the step plan, running the engine's layers, the launch ranking, the report, the CLI, end-to-end tests |

## 0. Rules (all parts)

1. **No language-model numbers.** A model may only PROPOSE how a column maps to a role and a code. Every proposal is
   validated deterministically (code format, code present in a local vocabulary or config, unit plausible) before it
   can be applied; a failed proposal becomes `unmapped` with the reason. All values, scores and rankings come from
   the data and the engine.
2. **What a model sees.** Only the profile's metadata: table and column names, inferred types, REDCap/data-dictionary
   labels and choices, units, numeric min/median/max, the top categorical levels with counts **only where the count is
   >= 11** (smaller levels are collapsed to `<11 other`), share missing. Never row-level records, participant ids,
   free text values, dates or anything matched by `measure_it.byod.checks` identifier or geography patterns. A run
   without `--llm claude` (the default) never contacts any model.
3. **Approval gate.** `map` writes `mapping.yaml` (human-editable). `apply` refuses an unapproved mapping unless
   `--approve` is given; the mapping records who approved it and when (`approved_by`, `approved_at`).
4. **Privacy.** Identifier and geography columns are dropped (never written), exact dates become per-participant
   `day_index`, exact ages become 10-year `age_band` (90+ top-coded); the BYOD validator runs on the written folder and
   must pass. Private outputs stay under the run folder, which carries a `.gitignore` (`*`).
5. Reuse, do not re-implement: `measure_it.byod` (validate/ingest/evaluate/subgroup/deploy), `measure_it.harmonize`,
   `measure_it.similar`, `measure_it.geography.subgroup`, `measure_it.outreach`, `measure_it.scoring`.

## 1. Profile (`profile.json`, part A -> parts A/C)

```json
{"input": "<path>", "format": "csv|tsv|excel|parquet|json|jsonl|xpt|sas7bdat|sav|dta|redcap|omop|fhir|zip-of-these",
 "tables": [{"name": "records", "source_file": "...", "sheet": null, "n_rows": 312, "shape": "wide|long|per_day",
   "id_candidates": ["record_id"], "columns": [
     {"name": "compass31_vaso", "label": "COMPASS-31 vasomotor score", "dtype": "float|int|bool|category|text|date|datetime",
      "n_missing": 4, "n_unique": 23, "numeric": {"min": 0, "median": 3.2, "max": 8.3}, "levels": [["0", 140], ["1", 98]],
      "unit": null, "choices": null, "code_pattern": "icd10|loinc|rxnorm|snomed|none", "flags": ["identifier?"]}]}],
 "dictionary": {"source": "redcap_data_dictionary.csv", "n_fields": 180} }
```

## 2. Mapping (`mapping.yaml`, part A -> part B)

One entry per source column (or per long-table column group):

```yaml
dataset_id: maestro_demo
condition_hint: long_covid
measurement_hint: nailfold_capillaroscopy
approved_by: null          # set by --approve or by hand
approved_at: null
tables:
  records:
    participant_id: record_id
    day_index_from: visit_date        # optional: a date column turned into day_index per participant
    columns:
      record_id:        {role: participant_id}
      cohort:           {role: label, condition: long_covid, value_map: {LC: 1, HC: 0}, label_basis: clinical_case_definition}
      age:              {role: age, transform: age_band}
      sex:              {role: sex, value_map: {F: female, M: male}}
      raynaud_hx:       {role: self_report_history, icd10cm: I73.0, value_map: {Yes: 1, No: 0}}
      glyca:            {role: lab, loinc: 82730-3, unit: mmol/L, platform: nmr_nightingale}
      compass31_vaso:   {role: survey, instrument: compass31, item: vasomotor}
      cap_density:      {role: device, block: capillaroscopy, measurement_class: capillaroscopy}
      olink_il6:        {role: omics, layer: proteomics, platform: olink}
      dx_code:          {role: diagnosis_code}            # a column of ICD-10-CM codes (long table)
      rx_code:          {role: medication_code}           # RxNorm codes (long table)
      zip:              {role: drop, reason: geography}
      name:             {role: drop, reason: identifier}
      notes:            {role: drop, reason: free_text}
      mystery_17:       {role: unmapped, reason: "no rule, fuzzy score 0.41 < 0.80"}
    proposed_by:        # per column: rule | dictionary | fuzzy | llm, with confidence 0-1 and a one-line rationale
      glyca: {method: fuzzy, confidence: 0.91, rationale: "matched configs/harmonize_nightingale_loinc.yaml 'GlycA'"}
subgroup_analysis: null    # optional, passed through to the BYOD manifest (device_positive definition)
```

Roles: `participant_id, label, age, age_band, sex, demographic, self_report_history, diagnosis_code, lab,
medication_code, medication_flag, survey, device, omics, date, drop, unmapped`.

Part A additions (2026-10-07, part A report; implemented in `measure_it.agent.read`, `.profile`, `.mapper`, `.llm`):

* Entry points: `read.read_input(path) -> Dataset` (`.format`, `.container`, `.tables: [Table(name, df, source_file,
  sheet, kind, labels, choices, field_types, forms, ...)]`, `.dictionary`, `.warnings`); `profile.profile_dataset(ds)
  -> dict`, `profile.write_profile(profile, out)`; `mapper.propose_mapping(profile, dataset=None, *, dataset_id,
  condition_hint, measurement_hint, threshold=0.80, llm=None, out=None) -> dict`, `mapper.write_mapping(mapping, out)`;
  `mapper.profile_cmd(path, out)`, `mapper.map_cmd(path, out, llm=None, condition_hint=None, measurement_hint=None)`.
  Both commands write `<out>/.gitignore` (`*`) if absent. `read` raises `ReadError` (unreadable / unsupported) and
  `MissingDependency` (names the missing reader package).
* Profile: `format` is the detected layout (`redcap`, `omop`, `fhir`) or the file format (`csv`, `tsv`, `excel`,
  `parquet`, `json`, `jsonl`, `xpt`, `sas7bdat`, `sav`, `dta`), or `mixed`; a ZIP or a folder is not a format but
  `container: zip|folder` (the contract's `zip-of-these`). Extra keys: top-level `container`, `warnings`,
  `small_cell_threshold` (11); per table `kind` (`redcap_records`, `redcap_repeating:<form>`, `omop:<table>`,
  `fhir:<resource>`), `n_columns`, `labels_from`; per column `share_missing`, `field_type` and `form` (REDCap). Flags:
  `identifier`, `geography`, `date`, `exact_age`, `free_text` (all from `byod.checks` name / value patterns, the REDCap
  `Identifier?` column, extra OMOP names such as `person_source_value`, `location_id`), weak `identifier?` /
  `geography?` / `date?` (matched in the label only), `structural` (REDCap event / repeat / `<form>_complete`),
  `id_values_look_identifying` (the first id candidate fails `byod.checks.id_problems`). A flagged column has
  `levels: null`, `numeric: null`, `code_pattern: none`. In `levels`, the collapsed small levels are
  `["<11 other", null]` unless their total is >= 11; levels beyond the top 20 are `["<other levels>", n]`.
  `code_pattern` uses check digits: LOINC mod-10, SNOMED CT Verhoeff (+ concept partition); RxNorm / SNOMED need a
  name or label hint (bare integers are never called codes).
* Mapping, top level extras: `condition_id` (the hint resolved), `input`, `format`, `generated_at`, `generated_by`,
  `fuzzy_threshold`, `llm` (null, or provider / model / counts sent, accepted, rejected / errors / log file),
  `summary` (`by_method`, `by_role`, `n_columns`). Table extras: `source_file`, `sheet`, `kind`, `shape`,
  `event_column`, `repeat_instance_column`, `long_lab: {code_column, value_column, unit_column, code_system}`.
  `proposed_by.method` is `rule | dictionary | fuzzy | llm | none` (`none` = unmapped).
* Column-spec extras: `value_map` keys are the source values as strings (integral numbers without `.0`); a `label`
  carries `level_conditions` (level -> condition id, `control` or null; other conditions' levels map to null) and
  `label_basis: null` (the data owner sets it; part B then defaults to `proxy` with a warning); `sex` may carry
  `unmapped_values` (codes without a value label: never guessed, `value_map` null); `participant_id` may carry
  `hash_ids: true`; `self_report_history` carries `label` (the config's condition label) or `from: text` (a long
  table of reported condition names, mapped per value with `configs/harmonize_self_report_icd10.yaml` on apply);
  `lab` carries `loinc_name`, and for Nightingale `platform: nmr_nightingale`, `unit` (source unit), `loinc_unit`,
  `factor` (copied from the config row, never computed), or `unit: null` + `unit_note` when the source gives no
  unit (units are never assumed); long-lab columns are `{role: lab, part: code|value|unit}` and a text code column
  carries `code_map: {local name: LOINC}` + `unmapped_levels`; `drop.reason` may also be `structural`, `date`,
  `llm_proposed`; `unmapped` rows decided by a rule (OMOP concept ids, SNOMED CT columns) are not sent to a model.
* LLM step (`llm="claude"` only): default model `claude-fable-5-1` (override `model=` / `MEASURE_IT_LLM_MODEL`),
  structured output (`output_config.format` json_schema), server-side refusal fallback (`fallbacks: "default"`).
  Roles `participant_id, label, age, age_band, sex, date` are never accepted from a model. Every request payload is
  appended to `<out>/llm_requests.jsonl` before it is sent.

## 3. Outputs of `apply` (part B)

`<out>/byod/` a folder that passes `measure-it byod validate` (manifest.yaml with data_use, labels, comparator,
measurement, primary_analysis, subgroup_analysis when given; participants.csv; device_*.csv; ehr_labs.csv;
ehr_diagnoses.csv; ehr_medications.csv; survey_*.csv; self_report_history.csv; omics_*.csv).
`<out>/omop/` All of Us-shaped OMOP CDM v5.4 tables as parquet: `person, observation_period, condition_occurrence,
measurement, drug_exposure, observation` (survey answers and self-reported history as All of Us PPI-style
observations), Fitbit-shaped `activity_summary, heart_rate_summary, sleep_daily_summary` when wearable daily data
exist, plus a local `concept` table for every code used (`vocabulary_id`, `concept_code`, local `concept_id` >=
2,000,000,000 as OMOP reserves for local concepts; standard concept ids are NOT invented). Dates in OMOP tables are
synthetic offsets from a fixed anchor (2000-01-01 + day_index), stated in `omop/README.md`.
`<out>/omop/query_pack_check.json` the result of running the stage-D query pack (transpiled to DuckDB with sqlglot)
against this export: row counts per step, any SQL error. This is the first execution of those packs.

Part B additions (2026-10-07, part B report; implemented in `measure_it.agent.apply`, `.omop`, `.omop_run`):

* Entry points: `apply_mapping(tables, mapping, out, *, approve, public=False, omop=True) -> dict` (writes
  `<out>/.gitignore`, `byod/`, `omop/`, `apply_report.json`, `mapping_applied.yaml`; the report carries the BYOD
  validator's full `validation` dict, `ok`, `data_use_confirmation_required`, `errors_other_than_data_use`, dropped /
  unmapped / privacy-guarded columns and dropped-value counts); `write_omop(harmonized, out_dir)`,
  `export_omop_from_byod(byod_dir, out_dir)`; `check_query_pack(omop_dir, phenotype_json_or_dict) -> dict`.
* `apply_mapping` RAISES `measure_it.agent.apply.ApplyRefused` (a `PermissionError`, `.reasons`) and writes nothing
  when the mapping is unapproved without `approve=True`, has no table with a participant_id, or the participant ids
  look identifying (patterns of `byod.checks.id_problems`, an identifying id-column name, or `id_hashing:
  salted_sha256` in the mapping) and `MEASURE_IT_BYOD_SALT` is unset. With the salt, ids become `h<16 hex>` (the
  `byod ingest` rule) and the manifest says `id_hashing: none` (already hashed).
* Optional mapping keys read by part B: `data_use` ({deidentified, use_permitted, statement, irb_or_dua_reference};
  the user's confirmation for private data — without it placeholders are written, validation fails on `data_use`,
  and `data_use_confirmation_required: true`); `source` (named in the public `statement`); `title`, `owner`,
  `description`, `licence`; `comparator` (`{type, description}` or a type string; default `healthy` + a caveat);
  `synthetic` (sets `demo` too); `primary_analysis` (overrides the default); `id_hashing`. Column specs may add
  `label_basis` and `definition` (labels; default basis `proxy` with a warning), `labels: [...]` (several labels
  from one column), `label`/`condition` (self-report name), `rxnorm` (medication_flag), `feature` (device/omics column
  name), `description` (device).
* Long lab tables: table-level `long_labs: {code: <col>, value: <col>, unit: <col>, name: <col>}` or column specs
  `{role: lab, part: code|value|unit|name}` (optional `code_map: {local: LOINC}` on the code column). Rows without a
  valid LOINC code or a numeric value are dropped and counted.
* `role: date` (or `day_index_from`): calendar dates -> day_index from the participant's first date over all
  tables; a numeric column is taken as an already-relative day number. Repeated rows + a day index make device /
  survey files per-day; repeated rows without one are averaged per participant (warning).
* A `lab` with `platform` nmr/nightingale whose LOINC and unit match a row of `configs/harmonize_nightingale_loinc.yaml`
  in Nightingale units is written to `omics_nightingale.csv` (manifest `omics.nightingale: {platform:
  nmr_nightingale, units: nightingale_standard}`) so the BYOD harmonizer keeps the platform flag; other platform labs
  go to `ehr_labs.csv` (platform kept only in OMOP `measurement_source_value`, warned).
* OMOP specifics: study-label cases are recorded as `condition_occurrence` with the condition's most specific
  ICD-10-CM base-population code and a local `study_label` type concept (so a query pack's base population is the
  study's case group); device features go to `measurement` under the local vocabulary `MEASURE_IT_DEVICE`;
  `person.year_of_birth` is null and the band is in the extension column `person.age_band_source_value`;
  `sex_at_birth_*` columns are on `person` (All of Us Workbench shape); `heart_rate_summary` adds
  `mean_heart_rate_device`. Every choice is listed in `omop/README.md`.
* `query_pack_check.json` runs two variants: `as_exported` (null year_of_birth: cohort 0 by construction) and
  `year_of_birth_from_age_band` (check-only year of birth from the band midpoint).

## 4. Orchestrated run (part C)

`plan.json` lists the steps actually run and why skipped ones were skipped (e.g. no label -> no evaluate; no device
positive definition -> no subgroup). Steps: profile -> map -> (approve) -> apply -> omop_export -> byod validate ->
data_use_gate -> byod ingest -> byod evaluate -> byod subgroup -> byod similar -> query_pack_check -> byod deploy ->
geo_subgroup_burden -> launch ranking -> outreach -> teardown -> report. `<out>/LAUNCH_REPORT.md` + `launch.json`: top
candidate launch counties for the condition x measurement with
every component (burden or device-subgroup burden, vulnerability, desert, measurement-specific clinic capacity,
research readiness, the device's own evidence where it exists), the Monte Carlo interval, the clinicians and sites to
contact (counts; names only in the outreach workbook), the similar-people estimates and query packs, and every caveat
the layers produced. "Greatest impact" is defined transparently: where a device subgroup exists, the launch ranking is
the engine's composite (same combination, same weight set) with the subgroup condition's burden member replaced by the
estimated device-subgroup rate per 100k adults per county (geo_subgroup_burden, percentile-ranked like the burden it
replaces; the estimated number of adults in the subgroup is reported next to it); else the engine's evidence-weighted
(if the device has a known performance record in the rankings) or equal-weight ranking; outside the scored grid the
burden-only fallback of `measure_it.outreach`; the report names which.

Amendments by part C (2026-10-07): (1) the subgroup re-score uses the subgroup RATE, not the count, as the burden
member, so the composite stays comparable with the engine's (burden is a prevalence percentile there); the count is
reported. (2) Extra recorded steps: `omop_export`, `data_use_gate` (private data whose data_use is not confirmed stop
here, before any ingest; confirmation lives in mapping.yaml's `data_use` block, which part B copies into the manifest,
or comes from `--confirm-deidentified --confirm-use-permitted --data-use-statement`), `query_pack_check` (after
`byod similar`, which writes the phenotype the pack is generated from) and `teardown` (`byod remove` unless `--keep`;
a run never modifies or removes a dataset id that was already ingested). (3) The device-subgroup burden is computed in
memory for the run's phenotype and written only under the run folder; a launch run never writes the shared
`geo_subgroup_burden` table. (4) Demo (synthetic) records are never deployed unless `--demo-rankings` (this run only,
incompatible with `--keep`; the teardown restores the rankings).
