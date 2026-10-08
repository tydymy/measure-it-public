# Analysis plan — facility, clinic and research-site registries; clinic matching; Test 6 (clinic locations)

Written 2026-09-23, **before** any facility table, match or shuffle result was computed. Only schemas,
row counts and a handful of example rows of the input tables had been inspected. Deviations found
later are listed at the end of this file (section 9); nothing above that section is edited after results exist.

Scope: SPEC §F (provider and facility layer), §G (clinical research readiness), Phase 5 (clinic matching)
and the clinic-location part of validation Test 6. Modules: `measure_it.facilities.registry`,
`measure_it.facilities.matching`, `measure_it.facilities.query`, `measure_it.facilities.shuffle_test`.

## 1. Inputs (all already processed; DATA_AUDIT caveats that bind this work)

| table | unit | caveat carried into this work |
|---|---|---|
| `providers`, `provider_specialty_groups` (NPPES) | NPI | taxonomy is self-reported, not proof of treating any condition; lat/lon = ZIP5-as-ZCTA internal point; organisation NPIs under physician taxonomies are **group practices**, not facilities; subparts are separate NPIs |
| `facilities__hrsa` | HRSA site | point geocode from HRSA; 1,271 administrative-only sites do not see patients; no service-line fields |
| `trial_sites`, `clinical_trials`, `trial_conditions`, `trial_interventions`, `trial_outcomes` | NCT × site | geoPoint is a **city centroid** (labelled `point` upstream until 2026-09-23, now `city_centroid`); 22% of site names are sponsor placeholders; condition sets are MeSH-expanded (dysautonomia → MSA/Parkinson; mcas → mastocytosis/leukaemia), so `condition_literal_match` is the precision filter |
| `nih_projects`, `nih_project_conditions` | appl_id | one row per fiscal-year award; `award_amount` must be summed with `obligation_total` (once per appl_id, parent not subprojects); precision view = `match_tier == title_abstract` and not `likely_false_positive` (PASC false positives); VA awards have no amount |
| `geo_context__hpsa`, `geo_context__rucc`, `geographies` | 2024 county | HPSA is an administrative designation, not burden; CT = planning regions, legacy CT counties are not remapped |

## 2. Unified facility table and entity resolution

Source records:
* NPPES organisation NPIs (entity type 2) in any specialty/facility/research group (`providers`). Kind:
  `nppes_facility` (a facility-type taxonomy), `nppes_group_practice` (specialty taxonomies only), `nppes_research`.
* HRSA health-center sites (all; administrative-only sites flagged and not offered as clinic candidates).
* U.S. ClinicalTrials.gov facilities: `trial_sites` rows with `country_group == US` and `facility_identifiable`,
  additionally excluding non-facility locations (participants' homes, remote/virtual/online); keyed by the
  upstream `site_key` (normalised facility|city|state|country).
* NIH RePORTER awardee organisations in the U.S. (key: `org_ipf_code`, else organisation name).

Geocode precision (column `geocode_precision`, also `source_geographic_resolution`): `point` (HRSA address
geocode; RePORTER org point used when the ZIP has no ZCTA), `zcta_centroid` (NPPES; RePORTER ZIP; ClinicalTrials.gov
site ZIP when it resolves to a ZCTA in the registered state), `city_centroid` (ClinicalTrials.gov geoPoint when no
usable ZIP), `none`. A resolved facility takes the most precise member location (point > zcta_centroid > city_centroid).

Name normalisation: lower-case ASCII, `&`→and, punctuation removed, common abbreviations expanded (univ, ctr, hosp,
med, hlth, inst, natl, st→saint only as a leading token), legal suffixes removed (inc, llc, pc, pa, pllc, ltd, corp,
the …). *Informative tokens* = tokens that are not generic facility words (clinic, center, medical, health, hospital,
university, research, associates, group, …) or stop words.

Candidate links (never across states; never on name alone):

| method | rule | confidence |
|---|---|---|
| `shared_npi` | HRSA `site_npi` = NPPES NPI | 0.99 |
| `nppes_same_name_same_zip` | two NPPES organisations, identical normalised name, same ZIP5 | 0.95 |
| `exact_name_same_zip` | identical normalised name, same ZIP5, ≥1 informative token | 0.95 |
| `hrsa_org_name_same_zip` | HRSA health-center (organisation) name = other record's name, same ZIP5 | 0.85 |
| `fuzzy_name_same_zip` | token-sort similarity ≥ 92, shared informative token, same ZIP5 | 0.85 |
| `contained_name_same_zip` | informative tokens of one name ⊆ the other's, extra tokens only place names, same ZIP5 | 0.80 |
| `exact_name_same_city` | identical name, same city+state, ≤ 25 km apart, ≥ 2 informative tokens | 0.80 |
| `fuzzy_name_same_city` | similarity ≥ 95, same city+state, ≤ 25 km, ≥ 2 informative tokens | 0.75 |
| `contained_name_same_city` | containment as above with ≥ 2 informative tokens, same city+state, ≤ 25 km | 0.70 |

Clustering: edges processed by descending confidence (ties: distance, ids) with union-find under constraints —
at most one HRSA site and one NIH organisation per facility; no facility spans more than 25 km; for city-level
edges each record keeps only its best link per other source (prevents a chain such as one "Kaiser Permanente"
trial site absorbing every Kaiser practice in a city). NPPES-to-NPPES links are same-ZIP only.
`facility_id` = the anchor record's id (HRSA site > NPPES lowest NPI > NIH org > ClinicalTrials.gov key),
object id `facility:<facility_id>`. Every member keeps `match_method`, `match_confidence`, and the record it was
linked to; the facility keeps its weakest link (`match_confidence_min`).

### 2a. Validation of entity resolution (pre-specified)
* **V1 (ground truth).** HRSA sites with a valid `site_npi` present in `providers` give true HRSA↔NPPES links.
  Re-run the name/location rules with `shared_npi` disabled and report precision (strict: linked to the true NPI's
  NPPES record/group; lenient: linked to a record with the same normalised name as the true NPI) and recall.
* **V2 (location-shuffle false-link rate).** Permute the location fields (ZIP, city, state, coordinates) of the
  ClinicalTrials.gov and NIH records across records (seeded, 3 draws) and re-run cross-source link generation.
  Links that survive are chance matches; report them per method as an estimated false-link rate relative to the
  observed number of links.
* A seeded random sample of 25 links per method is written to `results/tables/facility_match_sample.csv` for
  human review; no numbers are claimed from it.

## 3. clinic_registry (implementation-oriented) and research_site_registry (research-oriented)

`clinic_registry`: resolved facilities with ≥1 NPPES or HRSA member. Specialty/facility groups (union of member
NPIs' self-reported taxonomies; HRSA membership adds `fqhc`), FQHC/look-alike flag and its source, HRSA site type,
county context (primary-care HPSA area/population designation, RUCC nonmetro), co-located individual NPIs by group in
the same ZIP (co-location, not affiliation), source object ids.

`research_site_registry`: resolved facilities with ≥1 U.S. target-condition trial site or NIH project. Trial
history by condition (literal-match view and MeSH-expanded view kept separately), recruiting (site status and
overall status), recency, interventional/observational, device/diagnostic-test trials, trials whose registered
text matches `configs/measurements.yaml` patterns (same fields as `trials_for_condition`: titles, summary,
keywords, interventions, outcomes; not eligibility), NIH projects (precision view and all tiers), active core
projects, award obligations via `obligation_total` (once per appl_id, no parent/subproject double count), distinct PI
profile ids, activity codes and condition topics. No burden variable enters this table; no composite score.

Long link tables: `facility_trials` (facility × NCT) and `facility_nih_projects` (facility × appl_id) so queries
filter by condition/measurement without re-deriving.

## 4. Clinic matching (Phase 5) — `find_candidate_clinics`

Inputs: geography (county FIPS, `geo:` id, state FIPS, or "County, State" name), condition (any alias; resolved
with `ontology.normalize.normalize_condition`, UNKNOWN unless status is `matched`), measurement (class id, bundle
id, label or alias from `configs/measurements.yaml` / `configs/relevance.yaml`; UNKNOWN if unresolved),
`radius_km` (default `scoring.yaml clinic_matching.default_radius_km` = 50), `max_sites` (default 15).

Candidate pool: geocoded facilities within `radius_km` of the geography's Census internal point, plus facilities
inside a county geography; HRSA administrative-only sites excluded. State geographies: facilities inside the state,
distance not ranked. UNKNOWN when the pool is empty.

Six characteristics, each ranked **separately** (never combined):

| characteristic | raw value(s) | rank key |
|---|---|---|
| clinical_specialty_match | condition groups (relevance.yaml `condition_specialties`) and measurement implementer groups (`measurement_implementers`, bundle members unioned) present in the facility's taxonomy groups | matches both > one; then number of distinct groups |
| relevant_trial_history | U.S. trials of the condition at the facility (literal view); recruiting/active; latest start year | count, then active, then recency |
| NIH_research_activity | NIH core projects matching the condition (precision view); active; obligations | core projects, active, obligations |
| community_access | flags: FQHC/look-alike (HRSA or NPPES FQHC taxonomy); community/rural-health/critical-access/public-health taxonomy; county primary-care HPSA area/population designation; county nonmetro (RUCC 4-9) | number of flags |
| distance_to_target_population | km to the internal point (+ inside flag) | ascending |
| technology_experience | trials at the facility (any target condition) whose registered text matches the measurement's classes | count, then count for this condition |

A facility is *eligible* only with a facility-intrinsic relevance signal (specialty match, trial history, NIH activity,
technology experience, or community-health facility type); county-level flags and distance alone never qualify it.
Selection of the returned list: round-robin over the six rankings (1st of each, then 2nd of each, …), recording the
characteristic(s) that selected each facility; ties broken by distance then facility_id. Each candidate carries raw
values, per-characteristic rank and number ranked, plain-language reasons, and the framing sentence "This facility has
characteristics suggesting it may be a viable implementation or study partner." Never "best".

## 5. Test 6 — clinic-location shuffle (negative control)

Question: do the matched sets depend on where facilities actually are, or would randomly placed facilities give
the same characteristics?

* Target geographies (fixed before running): 2024 counties in the 50 states + DC with 2020 population ≥ 10,000
  (`scoring.yaml ranking.min_population`); seeded (`config.SEED`) random sample of 100 metro (RUCC 1-3) and 100
  nonmetro (RUCC 4-9) counties, plus San Diego County (06073, the SPEC geography-explorer example).
* Pairs: primary (long_covid, wearable_autonomic_activity_monitoring); secondary (me_cfs, same bundle) and
  (pots, autonomic_function_testing).
* Null: permute the location tuple (lat, lon, county, state, geocode precision, and the county-level HPSA/RUCC
  flags that belong to the location) across all geocoded candidate facilities **nationally**, 200 seeded draws.
  Facility-intrinsic attributes (taxonomy groups, FQHC status, trials, NIH projects) stay with the facility.
  Secondary: 100 draws permuting within state.
* Invariance check (must hold exactly): the number of facilities within the radius of each geography is unchanged
  by the permutation (the set of occupied coordinates is identical).
* Metrics per geography × pair (on the returned matched set, `radius_km=50`, `max_sites=15`): mean condition trials
  per candidate; share with ≥1 condition trial; share with a condition-specialty group; share with ≥1 NIH condition
  project; share with technology experience; share FQHC/look-alike; share in a county with a primary-care HPSA area
  designation; mean distance; number of eligible facilities in the pool; plus km to the nearest facility with ≥1
  condition trial (unbounded) and Jaccard overlap of the observed vs shuffled matched set.
* Summaries, per stratum (all / metro / nonmetro): observed mean over geographies vs the null distribution of that
  mean (mean, 2.5–97.5 percentiles, two-sided empirical p = (1 + #|null − null mean| ≥ |obs − null mean|)/(1 + draws)),
  and the metro-minus-nonmetro contrast observed vs null.
* Pre-specified reading: the clinic layer carries location-specific structure if (a) observed matched sets are
  mostly replaced under shuffle (low Jaccard), and (b) research-history metrics show a metro–nonmetro contrast larger
  than the null (research infrastructure concentrated beyond facility density), and (c) community-access metrics of
  matched sets differ from the null (FQHCs located where shortage designations are). A metric whose observed value sits
  inside the null interval is reported as a **null result** for that characteristic.

No parameter (radius, thresholds, weights, sample) is tuned after seeing Test 6 output.

## 6. Research centres query — `find_relevant_research_centers(condition)`

Aggregates `facility_trials` (literal view by default) and `facility_nih_projects` (precision view by default) per
resolved facility for the condition: trial counts, recruiting/active, recency, device/diagnostic and measurement
trials, NIH projects/core/active/obligations/PIs, location, and object ids (`facility:`, `trial:`, `nih:`). UNKNOWN
when the condition does not resolve or has no trial or NIH record.

## 7. Outputs

Processed: `facilities`, `facility_source_links`, `facility_trials`, `facility_nih_projects`, `clinic_registry`,
`research_site_registry`. Results: `results/tables/facility_*.csv`, `results/tables/test6_clinic_shuffle*.csv`,
draft figure `results/figures/drafts/test6_clinic_shuffle.png`. Write-up: `docs/FACILITY_MATCHING.md`.

## 8. Language and guardrails

Candidate / may be a viable implementation or study partner / research readiness / measurement desert. A taxonomy
code does not mean a facility treats the condition; a trial registration is not evidence a measurement works; NIH
funding is a research-capability signal, not burden; ClinicalTrials.gov coordinates are city- or ZIP-level. No
patient location is inferred.

## 9. Deviations from this plan (added after the analysis)

1. **Entity-resolution rules were changed after the first build's human-review sample, and before any matching or
   Test 6 analysis.** The sample (kept as `results/tables/facility_match_sample_run1_before_rule_fixes.csv`) showed
   several defects. Numbered sites of one network linked to each other ("ClinEdge Site 146" / "Site 147"). Trial
   "facilities" that are only a city name ("Virginia Beach", "Akron") linked to local practices. Order-blind
   containment linked "Dept. of Gastroenterology, University of Utah Hospital" to "Utah Gastroenterology LLC".
   Added rules:
   * fuzzy and containment links require identical numeric tokens;
   * containment requires the same token order and at least one informative token that is not the city name;
   * trial-site names that are only a place, "Site Reference ID…" names and "contact sponsor for locations" names are
     excluded as placeholders;
   * the never-across-states rule is enforced in clustering for every link type, including same-ZIP NPPES grouping and
     shared NPI (32 edges rejected).

   V1 and V2 were re-run. First-run values (`results/tables/facility_build_stats_run1_before_rule_fixes.json`):
   328,496 facilities; V1 precision 0.970 strict / 0.988 lenient, recall 0.528; V2 false-link share 0.0053 of
   3,419 links. Final values: V1 0.972 / 0.990, recall 0.528; V2 0.0046 of 2,959 links.
2. Measurement resolution precedence (exact id > bundle label/alias > class name) was added because "CPET" matched
   both the `cpet` class id and a bundle alias.
3. Test 6 sample: San Diego County (the added demo county) is metro, so the metro stratum has 101 counties (201 in
   total).
4. The Test 6 summary adds an `outside_null_95` (percentile interval) flag next to the pre-specified empirical p. The
   two can disagree for metrics bounded at 0.
5. Descriptive additions that do not change any selection: `characteristics_absent_in_pool` / `absence_note` in
   matching results, a mobile-van note in reasons, and example outputs (`facility_candidate_examples.csv`,
   `facility_research_centers_*.csv`, draft overlay map).
6. `research_site_registry` is written as the canonical table next to the upstream partition
   `research_site_registry__clinicaltrials_gov`. `store.union_partitions("research_site_registry")` must not be run.
   (Integration fix 2026-09-23: the upstream table is now named `ctgov_facility_summary`, so the collision is gone.)
7. **Added at independent review (2026-09-23), after the Test 6 results existed.** No selection rule, radius, sample,
   draw count or seed was changed. Observed values, null intervals and empirical p are identical to the first run.
   * Multiplicity: BH-FDR q-values within each (pair, null) family of 40 tests (`q_bh_fdr_family`). The plan had no
     multiplicity rule, and the summary reports 240 p-values. The smallest attainable p with 200 draws is 1/201.
   * County-sampling uncertainty: the permutation p is conditional on the 201 sampled counties. A county bootstrap
     (4,000 draws, resampling within metro/nonmetro, seed `SEED+7`) of each county's observed-minus-null-mean
     difference now gives `geo_boot_diff` and `geo_boot_ci_lo/hi`. `robust_q05_and_geo_ci` marks differences with
     q ≤ 0.05 whose county interval excludes 0. The write-up treats only these as supported. Differences that pass
     the permutation test alone are reported as "sample-specific".
   * Criterion (a) (low Jaccard) holds by construction: a national permutation moves almost every facility away from
     any fixed point. It is reported as a sanity check, not as evidence.
   * `research_site_registry.has_precision_view_record` was added (literal-match trial or title/abstract NIH match).
     It changes no count.
8. **Integration (2026-09-23/24), after the Test 6 results existed; not pre-specified.** Trial tagging
   (`registry.trial_measurement_classes`, hence `facility_trials.measurement_classes`, technology experience and the
   research registry's measurement counts) now reads the audited class patterns that were moved into
   `configs/measurements.yaml` (integration item A1); before, it used the as-shipped patterns while the measurement
   text miner used the audited ones. Technology experience for the wearable bundle fell from 8.9% to 2.2% of matched
   Long COVID sites, so eligible pools, selected sites and every Test 6 metric were regenerated with the unchanged
   design, draws and seeds. The nonmetro technology-experience deficit that was supported before (0.024 vs 0.051) is
   not supported now (0.003 vs 0.005, p 0.26); `docs/FACILITY_MATCHING.md` §5 has the regenerated table and readings.
   The extended upstream placeholder flag (`trial_sites.facility_generic`, item A4) moves 532 excluded U.S. site rows
   between the exclusion categories; the kept 17,685 rows and every facility are identical.
