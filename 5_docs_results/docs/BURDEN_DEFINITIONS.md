# Burden definitions and proxy map (geo_condition_burden / geo_condition_features)

Pre-specified 2026-09-23 before any result was computed (see `docs/ANALYSIS_PLAN_GEOGRAPHY.md`).
The same definitions live in code as `measure_it.geography.features.BURDEN_DEFINITIONS` and are written to the
table `geo_burden_definitions`; `tests/test_geo_features.py` checks that every target condition has a definition
at every level and that the table and this document agree on the primary measure ids and levels.

**Burden is not vulnerability, not access, not research activity.** Every burden number carries its
`burden_evidence_level`, its `source_geographic_resolution` and an `inherited` flag. A proxy is never substituted
silently: a condition without a usable measure gets level **D** and a null value with a reason.

## Evidence levels (configs/scoring.yaml)

| level | meaning | uncertainty multiplier |
|---|---|---|
| A | direct measure of the condition (survey self-report of the condition; surveillance case counts) | 1.0 |
| B | closely matching coded condition (claims algorithm that contains the condition's ICD-10-CM code) | 1.25 |
| C | symptom / comorbidity / antecedent proxy | 2.0 |
| D | no usable burden estimate (value null) | null |

`A` for an inherited value means "direct measure of the **state**", not of the county (see below).

## Long schema of `geo_condition_burden`

One row per published or derived estimate. Columns: `burden_row_id`, `geo_id`, `geo_level` (national/state/county),
`geo_name`, `state_fips`, `county_fips`, `geo_vintage`, `in_canonical_geographies`, `condition_id`, `measure_id`,
`measure_label`, `metric_type` (prevalence_pct / modeled_prevalence_pct / incidence_per_100k / case_count),
`value`, `value_unit`, `ci_low`, `ci_high`, `ci_level`, `numerator`, `denominator`, `denominator_description`,
`period`, `period_start`, `period_end`, `year`, `stratum`, `adjustment` (crude / age_standardized / age_adjusted),
`burden_evidence_level`, `measure_role` (primary / county_proxy / alternate / other), `is_primary_measure`,
`inherited`, `derivation` (as_published / computed_rate / inherited_from_state / aggregated_from_county /
copied_as_proxy), `proxy_for_condition`, `suppressed`, `estimate_kind` (survey / claims / surveillance / modeled),
`geo_match_status`, `source_table`, `source_measure_id`, plus the provenance columns
(`data_layer`, `source_name`, `source_record_id`, `source_version`, `retrieved_at`,
`source_geographic_resolution`, `evidence_type`, `evidence_level`, `provenance_notes`).

All rows of the three burden partitions are carried (HPS national subgroups, all Lyme years, all CMS years,
adjustments and age strata), then flagged; primary rows are a small subset.

## Primary measure per condition x geography level

| condition | level | primary measure_id | ev. level | resolution / inherited | rationale |
|---|---|---|---|---|---|
| long_covid | national | `lc_current_pct_all_adults`, latest period | A | national | HPS "currently experiencing long COVID", all adults; the only direct public measure |
| long_covid | state | `lc_current_pct_all_adults`, latest non-suppressed period per state | A | state | same; the latest period (Aug 20 - Sep 16, 2024) suppresses 6 states, which fall back to their most recent published period (`period` column records which) |
| long_covid | county | state primary value, `inherited=True` | A (of the state) | **state**, inherited | no county long-COVID measure exists (HPS has none; PLACES 2025 has no COVID item). The state value is context, never county prevalence |
| me_cfs | national/state/county | `cms_mmd_51_prev_actual_all` 2022 | C | as published | CCW "Fibromyalgia, Chronic Pain and Fatigue" contains R53.82 (chronic fatigue) but not G93.3x: a symptom proxy. 2022 = latest final claims year |
| fibromyalgia | national/state/county | `cms_mmd_51_prev_actual_all` 2022 | B | as published | the same CCW algorithm contains M79.7 (fibromyalgia) but is a broad pain/fatigue composite (24 % of FFS beneficiaries in 2023) |
| migraine | national/state/county | `cms_mmd_59_prev_actual_all` 2022 | B | as published | CCW "Migraine and Other Chronic Headache" (G43 + all G44 headache syndromes) |
| lyme_disease | county | `lyme_incidence_per_100k` 2023 | A | county | reported confirmed + probable cases, county of residence; 2023 is the latest year and the only year on the 2024 CT planning regions |
| lyme_disease | state / national | `lyme_incidence_per_100k_agg` 2023 | A | aggregated from county | sum of county cases / sum of county populations (aggregation up, not downscaling); excludes cases with unknown county (`state_cases_missing_from_county_file_pct`) |
| ptlds | county / state / national | Lyme 2023 incidence, copied with `proxy_for_condition` | C | as Lyme | antecedent-exposure proxy: PTLDS follows treated Lyme disease; no PTLDS surveillance or survey exists (Lyme audit). Lyme incidence is not PTLDS prevalence |
| pots | all | none | D | — | no public population measure; no CCW algorithm contains G90.A; PLACES has no orthostatic item |
| dysautonomia | all | none | D | — | no public population measure (G90.x absent from MMD-exposed CCW algorithms) |
| eds_hsd | all | none | D | — | no public population measure (Q79.6 / M35.7 absent) |
| mcas | all | none | D | — | no public population measure (D89.4 absent) |
| post_infectious_syndrome | all | none | D | — | grouping node; its members' measures cannot be unioned without overlap information |
| ibs | all | none | D | — | K58 absent from MMD-exposed CCW algorithms; no PLACES item |
| gastroparesis | all | none | D | — | K31.84 absent; no PLACES item |
| endometriosis | all | none | D | — | N80 absent; no PLACES item |

`measure_id` values: HPS ids are the ingestion's ids; CMS ids are `cms_mmd_<code>_prev_<actual|agestd>_<all|lt65>`;
Lyme ids are the ingestion's ids plus `_agg` for state/national aggregates; PLACES ids are
`places_<MEASURE>_<crude|ageadj>`.

## County C-proxy column (PLACES 2025, BRFSS 2023) — the proxy map

These are **level C** rows (`measure_role = county_proxy`), shown in a separate feature column
(`burden_proxy_*`), never as the primary burden:

| condition | PLACES measure | why this measure | why not others |
|---|---|---|---|
| long_covid | PHLTH (frequent physical distress, crude % of adults) | long COVID is defined by persisting symptoms and activity limitation; PHLTH (>= 14 poor-physical-health days of 30) is the closest symptom/function measure published at county level | COGNITION and DISABILITY are kept as alternates; GHLTH is broader self-rated health |
| me_cfs | PHLTH (crude) | fatigue/PEM-driven poor physical health days | as above |
| fibromyalgia | ARTHRITIS (crude) | the BRFSS item asks about "some form of arthritis, rheumatoid arthritis, gout, lupus, or fibromyalgia", so it literally contains the condition, though dominated by osteoarthritis | — |
| all others | none | no PLACES measure is specific enough | — |

PLACES values are **modeled small-area estimates** (multilevel models borrowing strength from state BRFSS and county
covariates); CDC cautions against using them for program evaluation. Kentucky and Pennsylvania have no 2023-BRFSS
measures in the 2025 release, so their counties have no PHLTH/ARTHRITIS value (not imputed).

State-level PLACES proxy rows (`measure_role = alternate`, `derivation = aggregated_from_county`) are
population(18+)-weighted means of county crude estimates; used only for the state-level proxy check (Test 6 R5).

## Pre-specified alternates (sensitivity only)

| condition | level | alternate measure_ids |
|---|---|---|
| long_covid | state | HPS mean of the 2024 non-suppressed periods (`lc_current_pct_all_adults_2024_mean`); strict latest period (suppressed -> null) |
| long_covid | county | PLACES PHLTH (proxy column), COGNITION, DISABILITY (crude and age-adjusted) |
| me_cfs, fibromyalgia, migraine | state/county | CMS 2022 age-standardized; CMS 2023 preliminary actual; CMS 2022 actual age < 65 |
| me_cfs | county | PLACES PHLTH |
| fibromyalgia | county | PLACES ARTHRITIS |
| lyme_disease, ptlds | county/state | pooled 2022-2023 incidence (`lyme_incidence_per_100k_2022_2023`); 2019 incidence (pre-2022 case definition) |

## Inheritance rule (CONVENTIONS §3.5)

A state value attached to a county keeps `source_geographic_resolution = state`, `inherited = True`,
`derivation = inherited_from_state`. It is written only for long_covid (the one condition whose only direct
measure is state-level). No other condition inherits: a county with a missing CMS value (e.g. Connecticut planning
regions, which CMS publishes only on legacy counties) stays null.

**Where the inherited value is used (review note, 2026-09-23).** As pre-specified in the analysis plan, the
long-COVID county `diagnostic_desert` uses the inherited state value as its burden component. Every county in a state
therefore gets the same burden percentile. That index ranks counties partly by their state's HPS estimate and is not
evidence of county long-COVID burden. Such rows carry `desert_burden_resolution = 'state (inherited; no within-state
variation)'` and an explicit `desert_definition`. Downstream scores should use `burden_inherited` /
`desert_burden_resolution` to widen uncertainty or exclude these rows from county-level burden rankings.

## Counts (burden_count) — only at the same resolution

| condition | level | count | kind |
|---|---|---|---|
| long_covid | state | HPS % x ACS 2020-2024 adults 18+ | `survey_prevalence_x_acs_adults` (period mismatch noted) |
| long_covid | county | none (an inherited rate times a county population would be a downscaled state estimate) | — |
| lyme_disease | county / state | reported cases 2023 | `surveillance_cases` |
| PLACES proxies | county | crude % x PLACES `total_pop_18plus` | `modeled` |
| CMS measures | any | none: prevalence applies to Medicare FFS beneficiaries whose counts are published only as size classes | — |
| ptlds, level D | any | none | — |

## Known traps carried into the table

* CMS MMD county values for Connecticut are on legacy counties (09001-09015): kept with
  `geo_match_status = matched_ct_legacy`; the 9 planning regions have no CMS value.
* 1:1 FIPS renames (bridged, 2026-09-24). CMS MMD still publishes two counties under codes retired in 2015:
  Oglala Lakota County SD as 46113 (Shannon County; now 46102) and Kusilvak Census Area AK as 02270 (Wade Hampton;
  now 02158). Both are renames with an unchanged boundary (Census "Substantial Changes to Counties and County
  Equivalent Entities"), so the CMS ingestion carries their rows to the 2024 code (layer choice: ingestion, where
  the source's county codes are matched to `geographies`; the harmonised burden table then needs no special case,
  and the CDC Lyme and HRSA ingestions apply the same two renames at their own ingestion). The code CMS published stays in
  `geo_condition_burden__cms_mmd.fips_as_published` and in `source_record_id`, `fips_bridge` names the rename, and
  the bridge sentence opens `provenance_notes` (carried into `geo_condition_burden`). Boundary changes are not
  bridged: 51515 (Bedford city VA, merged into Bedford County 2013) stays `unmatched`. Before the bridge both counties
  had no CMS value (me_cfs / fibromyalgia / migraine county burden null) although CMS publishes one. Other sources,
  checked 2026-09-24 over every processed table (string columns named like fips / geo_id / county): SVI 2022, ACS,
  PLACES, RUCC 2023 and the ZIP/ZCTA-geocoded NPPES, ClinicalTrials.gov and NIH tables are on the current codes
  (`46113` in the ZCTA tables is the Indiana ZCTA 46113); CDC Lyme carries 2001-2009 counts to the new codes
  (`cdc_lyme.FIPS_RENAMES_2000S`). HRSA published one in-effect MUA/P on a renamed code (MUA/P 00104, Wade Hampton
  Census Area, whole county, on 02270; no HPSA component uses 46113 or 02270). Until 2026-09-24 the HRSA ingestion
  treated it as a legacy-coded designation that only "may" cover 02158, so Kusilvak's MUA/P coverage read
  `undetermined_legacy_geography`; the HRSA ingestion now places components on the two 1:1 renames on the 2024 code
  (`hrsa.FIPS_RENAMES_1TO1`; `mua_designations.county_fips_bridge`, HPSA `county_fips_method = fips_rename_1to1`),
  and only boundary changes (02261 Valdez-Cordova split, 51515 Bedford city merger) and legacy Connecticut codes
  remain unplaced candidates. MUA/P coverage is context only: it enters no burden, desert or opportunity score.
* A defined primary measure without a value gets a measured, source-specific `burden_unknown_reason`
  (`features.not_published_reasons`): for CMS measures, a county absent from every MMD file ingested (2012-2023,
  prevalence and context measures; CMS omits cells with fewer than 11 beneficiaries and documents no other omission),
  a county without a 2022 cell, or a single omitted cell; for Lyme, a county code absent from the CDC county file.
  After the bridge the county gaps of the CMS measures are 13 Alaska county-equivalents (CMS publishes 17 Alaska
  codes, the 2024 vintage has 30), Kalawao County HI (15005) and Manassas Park city VA (51685), all absent from every
  MMD file, plus the Connecticut planning regions and the island areas (their own reasons).
* CMS publishes integer percentages; ties are common and small county differences are below resolution.
  2012-2015 vs 2016+ are different coding regimes (ICD-9 -> ICD-10); only 2022 is primary.
* me_cfs and fibromyalgia share the identical CMS 51 values; any correlation between the two is by construction.
* Lyme 2022+ uses a laboratory-only case definition in high-incidence states; not comparable with <= 2021.
  County of residence is not county of exposure; reported cases undercount infections.
* HPS is an experimental, low-response online survey; state CIs are wide (typically +/- 2 points around ~5 %).
