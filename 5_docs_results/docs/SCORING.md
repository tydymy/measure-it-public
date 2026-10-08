# Deployment-opportunity scoring: method reference

Modules: `measure_it.scoring.metric_link` (measurement performance, implementation reach, expected yield; table
`measurement_performance`; pipeline step `scoring_metric_link`, section 9),
`measure_it.scoring.opportunity` (components, composite, Monte Carlo, table `deployment_opportunities`),
`measure_it.scoring.controls` (Tests 5 and 6), `measure_it.scoring.sensitivity` (Test 7, S1-S6, burden-rule before/after),
`measure_it.scoring.recommend` (`rank_deployment_opportunities`, table `deployment_candidates`, demo outputs),
`measure_it.scoring.report` (runs everything, writes `results/SCORING_RESULTS.md`). Pipeline step: `scoring_report`
(`uv run measure-it pipeline --offline --only scoring_report`).
Pre-specified plan and dated deviations: `docs/ANALYSIS_PLAN_SCORING.md`. Results (all numbers):
`results/SCORING_RESULTS.md`. This file describes the method only.

```
uv run python -m measure_it.scoring.report            # build table + Tests 5-7 + S1-S6 + burden-rule before/after + demo + results markdown (~2 min with 6 workers)
uv run python -m measure_it.scoring.opportunity       # table only
uv run python -m measure_it.scoring.recommend         # demo query only (needs the table)
uv run pytest tests/test_scoring.py
```

## 1. What a row is

`deployment_opportunities` has one row per condition (or condition set) x measurement x geography, with
`object_id = opportunity:<condition>|<measurement>|<geo_id>` (CONVENTIONS section 7). Rows: conditions `long_covid`,
`me_cfs`, `pots`, `dysautonomia`; sets `autonomic_activity_invisible_illness` (the demo cluster of
`configs/conditions.yaml`) and `long_covid_or_me_cfs` (the SPEC demo query); measurement bundles
`wearable_autonomic_activity_monitoring`, `nailfold_capillaroscopy`, `autonomic_function_testing`,
`exercise_capacity_testing`; levels county (3,144 counties of the 50 states + DC, 2024 vintage) and state (51).
The row is **ecological**: it describes a place, the facilities around it and public registries, never a person.

`data_layer = derived`, `evidence_type = derived_score`, `evidence_level` = the burden evidence level.
`source_geographic_resolution` is `state` for state rows and for county rows whose burden is a state value
inherited by the county (CONVENTIONS section 3.5), otherwise `county`.

## 2. Components

Every component is kept raw and normalised. Normalisation: percentile rank `pr(x) = rank(x, average ties) / n` over
the non-null rows of the 50 states + DC, per condition x measurement x level (the same function the geography
module uses). A component that is an average of sub-percentiles is re-normalised with `pr()`.

| component | raw columns | normalised | source |
|---|---|---|---|
| burden(g,c) | `burden_value`, `burden_evidence_level`, `burden_inherited`, `burden_source_resolution`, CI, `burden_mc_sd` | `burden_pct` | primary measure of `geo_condition_features` (docs/BURDEN_DEFINITIONS.md) |
| vulnerability(g) | `svi_overall` | `vulnerability_pct` | CDC/ATSDR SVI 2022 (state: population-weighted county mean, derived upstream) |
| diagnostic_desert(g,c) | `diagnostic_desert`, `desert_pct_*`, `diagnostic_desert_variant`, `desert_burden_resolution` | `diagnostic_desert_pct` | the geography module's index; access-only variant for level D |
| clinic_capacity(g,m) | `cc_impl_providers_in_geo_*`, `cc_impl_providers_within_50km_*`, `cc_impl_facilities_in_geo_*`, `cc_impl_facilities_within_50km_*` (counts and per-100k rates, each with `_pct`) | `clinic_capacity_pct` | NPPES individual NPIs and `facilities` with >= 1 of the measurement's implementer groups (configs/relevance.yaml) |
| clinic_capacity, activity basis (since 2026-10-07) | `clinic_capacity_basis` = `measurement_activity`; density version kept in `clinic_capacity_density(_pct)` | `clinic_capacity_pct` | for measurements with a dedicated billing code (configs/measurement_hcpcs.yaml; wearable bundle, autonomic testing, CPET): mean percentile of active clinicians per 100k adults in the county and within 50 km, active = billed the code to >= 11 Medicare FFS beneficiaries in 2024 (`geo_measurement_capacity`, facilities.activity). Nailfold capillaroscopy has no code and keeps the density basis. `clinic_capacity_basis: implementer_density` restores the old rule (sensitivity S8) |
| research_readiness(g,c,m) | `rr_condition_trials_pool_n`, `rr_condition_nih_core_pool_n`, `rr_tech_trials_pool_n` (+ `_pct`), `rr_condition_trial_sites_pool_n`, `rr_tech_trials_condition_pool_n`, `pool_facilities_n` | `research_readiness_pct` | `facility_trials` (literal condition match), `facility_nih_projects` (precision view), registered measurement classes |
| technology_saturation(g,m) | `ts_tech_facilities_pool_n`, `technology_saturation` (per 100k within 50 km) | `technology_saturation_pct` | facilities with >= 1 trial mentioning a member class; reported only |

* **Pools.** "In the county" uses the facility's / provider's county FIPS; "within 50 km" is the great-circle distance
  from the county's Census internal point (`scoring.yaml clinic_matching.default_radius_km`). The research pool is
  "in the county OR within 50 km", which is exactly the `find_candidate_clinics` pool (tested). States use in-state
  counts only.
* **Distinct counts.** Research sub-measures count distinct NCT ids / NIH core projects over the pool, not
  facility-trial pairs, so one multi-site trial counts once per region.
* **Per-capita rates** use ACS population (county/state) and `population_within_50km` (ZCTA-based, geography
  module); a zero or missing denominator gives NaN, never infinity, and `clinic_capacity` then averages the available
  sub-percentiles.
* **Level D burden** (POTS, dysautonomia: no burden measure is defined at any level): `burden_value` null,
  `burden_excluded_level_D = True`, burden weight removed and the other weights renormalised
  (`burden_weight_renormalised`), `rank_burden_only` null, desert = access-only (labelled). Never ranked on burden.
* **Incomplete burden** (since 2026-09-24, plan deviation 7): a condition whose burden measure IS defined at the level
  but has no value in a region (me_cfs: the 24 counties without a CMS value, i.e. the 9 Connecticut planning regions,
  which CMS publishes only on the legacy counties, Manassas Park VA, 13 Alaska boroughs/census areas and Kalawao HI;
  26 before the 2026-09-24 FIPS bridge below) is not ranked there: `burden_incomplete = True`, `burden_incomplete_reason`,
  `burden_missing_source_reason`, `burden_members_missing`; burden, every composite and every rank null (also
  `rank_equal_incl_small`), no Monte Carlo; `burden_evidence_level = D`, `burden_excluded_level_D = False` (the burden
  weight is not renormalised away: the row is not ranked at all). Non-burden components are kept. Listed in
  `results/tables/scoring_incomplete_burden.csv` and in the `excluded_incomplete_burden` block of the query.
  Before 2026-09-24 me_cfs ranked these counties on the other four components with the burden weight renormalised.
  Two former gaps were crosswalk gaps, not missing publications: Oglala Lakota County SD (46102) and Kusilvak Census
  Area AK (02158) are the 2015 renamings of Shannon County (46113) and Wade Hampton Census Area (02270), which CMS
  still publishes. Since 2026-09-24 the CMS ingestion bridges these 1:1 renames (`cms_mmd.FIPS_RENAMES`;
  docs/BURDEN_DEFINITIONS.md), so both counties have a ME/CFS value and are complete (Oglala Lakota: equal-weight
  rank 443 of 2,395, burden-only rank 307 for the primary query; Kusilvak is below the population minimum);
  `report.bridged_fips_renames` reports the bridge in SCORING_RESULTS section 5. Every remaining gap is a county absent
  from every MMD file ingested, or a Connecticut planning region.
* **Long COVID burden since 2026-10-07** (`configs/scoring.yaml long_covid_burden: brfss_sae`): the BRFSS 2023
  multilevel-regression / post-stratification estimate (`geography.sae`, CDC MMWR definition,
  `geo_condition_burden__brfss_long_covid_sae`), at county, state and national level, `estimate_kind =
  modeled_small_area`, its own 95% interval, not inherited. Within-state differences come only from county
  composition and covariates (results/BRFSS_SAE_RESULTS.md: state agreement with CDC 52/52, leave-one-state-out MAE
  0.75 pp vs 0.92 for the national mean, metro-area check r = 0.16). Kentucky and Pennsylvania are absent from the
  2023 public file; their values are model predictions. The HPS rows stay in `geo_condition_burden` as alternates.
* **Inherited burden** (`long_covid_burden: hps_inherited`, the rule before 2026-10-07, sensitivity S7): the state HPS
  value, `burden_inherited = True`, source resolution `state`; the row's `uncertainties` says so, and the Monte Carlo
  adds a per-county deviation (section 4).

### Condition sets
* A member has a **defined** burden at a level when its burden measure exists there (`burden_defined_level` != D in
  `geo_condition_features`: long_covid and me_cfs at county and state; pots and dysautonomia never). Members without
  one are listed in `burden_members_excluded_level_D` and never enter the burden.
* **Completeness rule (plan deviation 7, 2026-09-24).** The raw set burden `burden_set_member_mean_pct` = mean of the
  burden percentiles of the members with a defined burden, computed only where EVERY such member has a value;
  `burden_pct = pr(burden_set_member_mean_pct)` over those complete regions, re-normalised like every other averaged
  component (reviewer correction 2026-09-23: the first build used the un-normalised mean, whose spread is about 25%
  smaller, so the set burden carried less than its nominal weight). A region where a defined member has no value is
  **incomplete** (above): no set burden, composite or rank, listed separately. Until 2026-09-23 the mean ran over the
  members with a value in the region, so the Connecticut planning regions, Oglala Lakota SD (then without a CMS value),
  North Slope AK and Manassas Park VA got a set burden from the inherited long-COVID state value alone, a different
  quantity ranked as if it were the same (Greater Bridgeport CT was rank 9 of the primary query and rank 10 of the
  demo cluster; Oglala Lakota rank 1 of the burden-only list). The old rule is kept only as `make_combo(..., burden_rule="available")` for
  the before/after tables (`results/tables/scoring_set_burden_rule_before_after.csv`, `scoring_set_burden_rule_top10.csv`).
* A single condition is a set of one: the same rule makes me_cfs's counties without a CMS value incomplete.
* `burden_value` is null for a set and the members' raw values, levels, CIs, deserts, `burden_defined_at_level` and
  `contributes_to_set_burden` are in the `member_components` JSON column; set evidence level = least direct
  contributing level (D in an incomplete region); `burden_inherited` = any member's value in the region inherited;
  `burden_members_contributing` lists the members whose percentile enters the set burden (empty when incomplete).
* `diagnostic_desert` = mean of the members' burden-informed deserts (access-only mean only if no member has one).
* `research_readiness` counts distinct trials / NIH core projects over the union of member conditions.
* `clinic_capacity` and `technology_saturation` do not depend on the condition.
* An ad hoc set (`"a+b"`) is scored on the fly by the same code (`rank_deployment_opportunities("A or B", ...)`).

### The adapter swap
A measurement enters only through `measurement_spec(m)` -> member classes + implementer groups (+ the registered
adapter's id/status from `measurements.adapters`). The capillaroscopy stub bundle runs through the identical
`make_combo -> facility_measures -> evaluate` path as the wearable bundle; nothing in the geography, clinic or
recommendation code knows which modality it is.

## 3. Composite and ranks

`composite_<ws> = sum_k w_k * normalised_k / sum_k w_k` over the AVAILABLE normalised components (burden_pct,
vulnerability_pct, diagnostic_desert_pct, clinic_capacity_pct, research_readiness_pct; for `evidence_weighted` also
the absolute-scale `measurement_evidence` and `expected_yield`, section 9), for every weight set of
`configs/scoring.yaml`, including `saturation_adjusted` (1/6 each on the five components and on
`1 - technology_saturation_pct`; moved from `opportunity.EXTRA_WEIGHT_SETS` into the config on 2026-09-23, values
unchanged).
`composite` / `rank` = the `equal` default. `rank_<ws>`: 1 = highest composite among ranking-eligible rows
(`rank_eligible`: counties with population >= `ranking.min_population`, all states, AND a complete burden); ties broken
by geo_id. Small counties keep their scores, are flagged (`small_population_flag`) and are ranked only in
`rank_equal_incl_small`. Incomplete-burden regions have no composite and are ranked nowhere.

## 4. Uncertainty: Monte Carlo rank intervals

1,000 seeded draws (`mc_n_draws`, `mc_seed`) per condition x measurement x level:
* burden per member: `value + e`, truncated at 0, `e ~ Normal(0, sd)` with
  `sd = CI half-width / 1.96 x evidence multiplier` (HPS) or `sd = sqrt(p(1-p)/n x 1e4 + 1/12) x multiplier`
  (CMS: n = geometric midpoint of the published Medicare FFS size class; 1/12 = integer rounding). Multipliers:
  `scoring.yaml evidence_level_uncertainty_multiplier` (A 1.0, B 1.25, C 2.0).
* inherited county burden: the state's error is drawn once per state per draw and shared by its counties, plus an
  independent per-county deviation `tau` = method-of-moments between-state SD of the HPS state estimates
  (`mc_tau_inherited`), an explicit assumption that counties vary around their state as much as states vary
  around the nation. Since 2026-09-23 the state's error is also multiplied by
  `scoring.yaml inherited_burden_uncertainty_multiplier` (2.0, the level-C value; column
  `burden_inherited_uncertainty_multiplier`, `mc_inherited_multiplier`): a state value carried to a county has no
  within-state information, so as a county burden it is treated like a proxy. Sensitivity S1b compares 1.0 (the
  earlier behaviour) with 2.0 for the primary query.
* desert recomputed with the drawn burden percentile; all percentiles recomputed per draw.
* weights ~ Dirichlet(20 x default weights); level-D rows drop the burden weight.
Columns: `rank_mc_p05`, `rank_mc_p50`, `rank_mc_p95`, `p_top10_mc`, `p_top25_mc`; decomposition with only burden
noise (`rank_mc_burden_p05/p95`) or only weight noise (`rank_mc_weights_p05/p95`).

## 5. Tests

* **Test 5** (`controls.test5_combo`): burden_only, burden+vulnerability, and burden+vulnerability+desert rankings vs
  the full ranking: top-10/25 overlap, Jaccard, tie-averaged overlap, expected overlap N^2/M, Kendall tau-b,
  Spearman; entering/leaving regions with component contributions `w_k (x_k - mean_k)` as the reason.
* **Test 6** (`controls.label_shuffle`, `controls.clinic_shuffle`, `controls.phenotype_label_nulls`): geography-label
  permutations of burden (among complete-burden regions only, `controls.permutation_among`, so the ranked set is the
  same in every draw; since 2026-09-24) and of SVI (national, within state); facility + individual-provider location permutations
  (national, within state) with clinic_capacity / research_readiness / technology_saturation recomputed for every
  county (the desert's provider/trial sub-percentiles come from the geography module and stay at their observed
  values; the permutation keeps the number of facilities per location slot, so it cannot move a pure density
  measure; since 2026-10-07, for an activity-basis clinic_capacity the activity clinicians' location slots are
  permuted the same way, `facilities.activity.capacity_pct_for`); phenotype-label permutation results restated from
  the wearable modules' tables. Test 8 (temporal holdout) keeps the density basis in every dated view, because the
  2024 billing data post-date every cutoff.
* **Test 7** (`sensitivity.test7_combo`): all named weight sets + 1,000 flat Dirichlet weight vectors; pairwise
  Kendall tau-b, Kendall's W, per-region rank range/IQR and top-N frequencies; unstable regions (top-25 under >= 1
  named set AND top-25 in < 50% of named sets or of random draws) in `results/tables/test7_unstable_regions.csv`.
* **Sensitivity** S1 tau = 0; S1b inherited-burden multiplier 1.0 vs 2.0 (added 2026-09-23, not pre-specified); S2
  PLACES PHLTH proxy for long-COVID county burden (counties without a PLACES value are incomplete under S2 and drop
  out of the comparison); S4 clinic capacity without primary care; S5 measurement swap; every comparison is over the
  regions ranked in both rankings (`n_compared`, Kendall tau-b and Spearman reported).
* **S6 `specialist_only`** (`sensitivity.s6_specialist_only`; declared 2026-09-24 before it was computed, plan
  deviation 8): primary care removed from relevant-provider density everywhere the composite uses it. (i) The
  diagnostic desert's provider term becomes `1 - pr(relevant specialists per 100k)`, the distinct individual NPIs with
  >= 1 taxonomy in the condition's core groups (configs/relevance.yaml `condition_specialties`) minus `primary_care`,
  per member condition, counted here from `providers` (`opportunity.relevant_provider_counts`: county = primary
  practice location in the county, as clinic_capacity; state = the NPI's state FIPS, which also counts the ~2% of NPIs
  without a county geocode, as the geography module does). The counts equal the geography module's
  `relevant_specialists_n` and, per group, `provider_density_county` (`test7_sensitivity_specialist_only_qa.csv`), and
  the resulting single-condition desert equals its `diagnostic_desert_sens_specialists` (tested). (ii) clinic_capacity
  uses the measurement's implementer groups minus `primary_care` (= S4; facility types such as `fqhc` and
  `community_health_clinic` stay). Reported against the primary equal-weight ranking: top-10/25 overlap, top-25
  Jaccard, Kendall tau-b, Spearman, for the primary query at county and state level and every condition x measurement
  at county level, with the two single swaps (desert provider term only; clinic capacity only) as decomposition
  (`test7_sensitivity_specialist_only.csv`, `_top25.csv`).
* **Burden-rule before/after** (`sensitivity.set_burden_rule_before_after`): every combination under the old
  (`available`) and the current (`complete`) rule: regions ranked, top-10/25 overlap, Kendall tau-b over the regions
  ranked under both, incomplete regions that were in the top-10/25 or the burden-only top-25, and the top-10 lists.

## 6. `rank_deployment_opportunities(condition, measurement, geography_level='county', top_n=10, weight_set='equal')`

Resolves aliases ('Long COVID or ME/CFS' -> `long_covid_or_me_cfs`; any "A or B" of resolvable conditions -> an ad hoc
set; 'wearable autonomic monitoring' -> the wearable bundle), reads `deployment_opportunities` (or scores the
combination on the fly with the same engine), and returns per region the SPEC output schema:

| key | content |
|---|---|
| condition | id, label, kind, members, `condition:` ids |
| phenotype | `get_patient_phenotype_signature` per member (public person-level cohorts; group differences; never these regions' people) |
| measurement | bundle, member classes, implementer groups, deployment complexity, adapter and status |
| technology | name; `regulatory_context` = `get_regulatory_context` per member class; `technology_evidence` = `get_measurement_evidence` per member condition; `measurement_performance` = the scored `measurement_performance` record of the queried condition (or set) x measurement (tier, AUROC and CI, sensitivity at 0.90 specificity, comparator, label basis, healthy-control caveat; UNKNOWN stays UNKNOWN); `member_condition_performance` for a set; `expected_yield_here` = the region's expected yield block |
| geography | name, fips, burden (+ members for a set), burden_evidence_level, vulnerability, diagnostic_desert, clinic_capacity, research_readiness, technology_saturation, composite, rank, `rank_basis` (weight set, rank column, universe), rank interval and P(top 10) with `monte_carlo_basis` (always the equal-weight Monte Carlo among population-eligible regions; under another weight set or with small counties included it says it is not an interval for the displayed rank), ranks under every weight set |
| candidate_sites | `find_candidate_clinics` per member condition, merged by facility, with lat / lon and `geocode_precision` of the facility, reasons, characteristics and the framing sentence |
| research_evidence | `trial:` and `nih:` ids at facilities in the region's pool (condition trials, technology-experience trials, NIH projects) |
| molecular_context | `get_molecular_context` summary per member, labelled condition-level molecular enrichment |
| uncertainties | row-specific (inherited, proxy, level D, rank interval, small population, matcher absence notes) and the measurement's performance tier, comparator and label basis, the healthy-control caveat, the yield basis and any evidence conflict |
| provenance | object ids, sources, versions, tools; `opportunity_row` = the region's `opportunity:` id with `stored` true/false. A precomputed combination cites its stored row first. A combination computed on the fly (an ad hoc condition set, a measurement class that is not one of the four scored bundles) has no stored row: its id is not cited in `object_ids` (it is marked "computed on the fly; not a stored row"), which list the traceable component ids instead (geo, member conditions, measurement, and each member condition's stored opportunity row where one exists) |
| recommended_next_step | a candidate deployment opportunity for pilot evaluation, never a validated diagnostic pathway; names the weight set of the rank and, when that is not the equal-weight rank among population-eligible regions, says the Monte Carlo interval belongs to the equal-weight rank (and gives it) |

`excluded_incomplete_burden` lists, with `excluded_incomplete_burden_note`, the regions of the combination that are
not ranked because their burden is incomplete (population-eligible ones unless `include_small_population=True`): name,
fips, opportunity (None for an on-the-fly combination) and geo object ids, the members without a value, the reason, the kept non-burden components, and
`UNKNOWN / NOT AVAILABLE` for burden, composite and rank.

Unresolvable inputs, an unknown weight set, a `top_n` that is not a positive integer, or a ranking that does not
exist (e.g. `burden_only` for a level-D condition) return `status = UNKNOWN / NOT AVAILABLE` with a reason and no
recommendations. `include_small_population=True` (county) ranks every scored county on the same composite for any
weight set (`rank_<ws>_incl_small`, ties by geo_id); small counties have no Monte Carlo interval. Measurements are
resolved by the one shared resolver `measure_it.measurements.resolve.resolve_measurement` (class ids / names / aliases
of configs/measurements.yaml, bundle ids / labels / aliases of configs/relevance.yaml, then the class patterns; the
same resolver backs every query tool and the agent), so the SPEC demo phrase 'wearable monitoring of
autonomic/activity abnormalities', 'wearable autonomic monitoring' and 'wearable' all give the wearable bundle. In the
geography block, `burden_evidence_level_note` states when
the level describes an inherited state measure rather than county prevalence; `deployment_candidates` rows with
inherited burden carry `source_geographic_resolution = state`.

## 7. Guardrails built into the code

* No county prevalence is created from a state value; inherited rows are flagged in the row, in `uncertainties`,
  in `source_geographic_resolution` and in the recommendation JSON.
* Proxy (C) burden gets a x2.0 uncertainty band and a row-level note; level D is excluded from burden-weighted
  ranking and flagged.
* A burden is never assembled from part of the members that have one: a region missing a defined member's value is
  incomplete, unranked and listed with the reason (no imputation, no mapping of legacy CT county values to planning
  regions).
* Every component is kept; composites are reported under every weight set with rank intervals.
* Person-level results (phenotype) and condition-level molecular evidence are attached as context only; they are
  not used to score places and are never joined to the geographic rows.
* Curated relevance maps (implementer groups, condition specialties) are assumptions (evidence_type
  `curated_config` upstream); the recommendation says so.

## 9. The metric -> translation link (added 2026-09-25)

Pre-specified in `docs/ANALYSIS_PLAN_METRIC_LINK.md` (before any number was computed); results in
`results/SCORING_RESULTS.md` section 12. Module `measure_it.scoring.metric_link`.

**Performance records.** One discrimination result per record: target condition (a condition set only when the
analysis pooled exactly its members), measurement class of the feature/device, comparator, label basis, n, AUROC with
95% CI, and the operating point at the pre-specified 0.90 specificity. Sources, closed list: this project's result
tables (MUSCLE-ME steps and VO2_rel, Uwakwe resting HR, NHANES ME/CFS-like proxy, Appelman VO2max and metabolomics,
MY-LC cortisol) and, flagged as tier 3, curated `published_device_evidence` diagnostic-accuracy rows. Quality tier:
1 own computation on a clinical/study case definition > 2 own computation on a proxy or self-reported label > 3
published claim (not reproduced) > 4 none. Sensitivity at 0.90 specificity: computed here from stored person-level
values where they exist (threshold = the ceil(0.9 n0)-th smallest control score; stratified bootstrap CI with the
threshold re-chosen); published only when the paper reports specificity >= 0.90; otherwise UNKNOWN (never
extrapolated). Records: table `measurement_performance_records`, `results/tables/measurement_performance_records.csv`.

**The scored record** of (condition, measurement class or bundle) is chosen by best tier, then status (known = AUROC
+ operating point > partial = AUROC only > UNKNOWN), then the most conservative evidence, then larger n. Rows without
any record are UNKNOWN; nothing is carried from a neighbouring condition, set, class or dataset. A lower-tier
own-computed record whose AUROC CI includes 0.5 next to a scored record whose CI excludes it sets `evidence_conflict`.
Table `measurement_performance` (`object_id = measurement_performance:<condition>|<measurement>`, trace namespace
`measurement_performance`): the grid of the scored conditions and sets x the four bundles, their member classes and
the extra classes with own records.

**Components (absolute 0-1 scales, not percentile ranks).**
* `measurement_evidence = tier_factor x clip((AUROC_lower95 - 0.5) / 0.5, 0, 1)`, tier factors 1.0 / 0.75 / 0.5.
  Constant within a condition x measurement: it cannot reorder regions of one measurement; it matters across
  measurements (joint ranking) and for whether a measurement is ranked at all.
* `expected_yield = pr(B) x op_sensitivity x reach`, B = burden percentile (prevalence-percentile path) or the
  percentile of the burden count (count path). Not re-normalised, so a less sensitive measurement contributes less.
* `reach(g, m)` = share of the region's ZCTA population (ACS) whose ZCTA internal point is within the default radius
  of a geocoded clinic-candidate facility with >= 1 implementer group of m. ZCTA -> county by the 2024 county of the
  internal point; regions without a populated ZCTA have reach UNKNOWN (independent cities).

**Expected yield and false positives** (columns of `deployment_opportunities`, every factor kept): count path only
where a defensible count exists at that resolution (a single condition, level A/B, not inherited, % of all adults:
long COVID at state level): `burden_count = value / 100 x adults 18+`, `expected_detectable_cases = count x
sensitivity x reach`, `expected_false_positives = (1 - specificity) x (adults - count) x reach`, `expected_ppv`,
`false_positives_per_detected_case`. Everywhere else (inherited county values, level-C proxies, condition sets): the
prevalence-percentile path, `expected_yield_index = burden_pct x sensitivity x reach` (an index, not a count) and
`expected_false_positives_upper = (1 - specificity) x reached adults`. Language fixed in every output: planning
estimates for a candidate pilot, not predictions of diagnoses; performance against healthy controls overstates
real-world performance.

**Weight set `evidence_weighted`** (configs/scoring.yaml; the pre-declared primary alternative, `equal` stays the
default): burden 0.15, vulnerability 0.15, diagnostic_desert 0.15, clinic_capacity 0.10, research_readiness 0.10,
measurement_evidence 0.15, expected_yield 0.20. **UNKNOWN rule:** performance UNKNOWN or partial -> no
evidence_weighted composite or rank anywhere in the combination (`evidence_weighted_status` says why; never scored 0
or 1; the equal-weight rank is kept and says it carries no performance information); performance known but burden
level D -> the yield weight is removed and the others renormalised (as for burden); known performance but reach
UNKNOWN in a region -> no composite there.

**Uncertainty.** `rank_mc_ew_*`, `p_top10_mc_ew`, `expected_yield_mc_p05/p50/p95`, `expected_false_positives_mc_*`:
1,000 draws of burden (as section 4), operating-point sensitivity and specificity (logit-normal from their 95% CIs)
and weights ~ Dirichlet(20 x evidence_weighted). measurement_evidence (a lower bound) and reach are not drawn. Under
`weight_set='evidence_weighted'` the recommendation's rank interval is this one and says so.

**Joint ranking.** `rank_equal_joint`, `rank_evidence_weighted_joint`: per condition x level, all ranking-eligible
region x measurement pairs of the four bundles ranked by the composite; bundles without known performance have no
evidence-weighted pairs. **T5-4** (`controls.test5_evidence_combo`): evidence_weighted vs equal per combination.
**S5 with evidence** (`sensitivity.s5_evidence`, `sensitivity.joint_composition`): within-bundle agreement under both
weight sets (or "not ranked" with the performance status) and the bundle composition of the joint top-10/25.
