# Analysis plan: deployment-opportunity scoring, recommendations and robustness (SPEC Phases 4-6, Tests 5-7)

Status: **pre-specified 2026-09-23, before any score, rank or test statistic was computed.** Only the schemas,
row counts, missingness and a few example rows of the input tables, the upstream docs
(`docs/BURDEN_DEFINITIONS.md`, `docs/FACILITY_MATCHING.md`, `results/GEOGRAPHY_RESULTS.md`) and the DATA_AUDIT
files of the consumed sources (cdc_long_covid, cms_mmd, cdc_svi, census_acs, cdc_places, nppes,
hrsa_health_centers, clinicaltrials_gov, nih_reporter) had been read. Deviations made later are listed at the
end, dated, with the reason; nothing above that section is edited after results exist.

Owner modules: `measure_it.scoring.opportunity` (components, composite, Monte Carlo, table),
`measure_it.scoring.controls` (Tests 5 and 6), `measure_it.scoring.sensitivity` (Test 7 and pre-specified
sensitivity analyses), `measure_it.scoring.recommend` (`rank_deployment_opportunities`, SPEC output schema,
`deployment_candidates`, demo query). Method reference: `docs/SCORING.md`. Results: `results/SCORING_RESULTS.md`.

## 1. Questions

* **Q1 (build).** For geography g, condition c and measurement m, compute the SPEC Phase 4 components
  burden(g,c), vulnerability(g), diagnostic_desert(g,c), clinic_capacity(g,m), research_readiness(g,c,m) and
  technology_saturation(g,m); keep every raw and normalised column; compose them under every weight set; and
  propagate proxy/inherited burden and weight uncertainty into a rank interval per row.
* **Q2 (Test 5).** Does adding facility and research information change where the engine would deploy a
  measurement, compared with ranking by burden alone? Reported plainly either way.
* **Q3 (Test 6, scoring).** Does the ranking lose structure when (a) the geography labels of burden or of
  vulnerability are shuffled, (b) clinic locations are shuffled? (c) Cite the phenotype-label permutation
  results of the wearable modules.
* **Q4 (Test 7).** How stable are the rankings across the configured weight sets and random weights? Which
  regions are unstable?
* **Q5 (recommendation).** For the SPEC demo query, return the SPEC output schema per region with every
  number traceable to an object id.

## 2. Units, universe, rows

* Geography levels: `county` (2024 vintage; Connecticut = 9 planning regions) and `state`. Universe = 50 states +
  DC (`geo_condition_features.in_analysis_universe`): 3,144 counties, 51 states. Every percentile and rank is
  computed within this universe, per condition x measurement x level.
* **Ranking-eligible** counties: `population_total >= 10,000` (`scoring.yaml ranking.min_population`); smaller
  counties are scored, reported and flagged `small_population_flag`, their `rank_*` is null, and a separate
  `rank_equal_incl_small` ranks all counties. All 51 states are eligible.
* Conditions (rows): `long_covid`, `me_cfs`, `pots`, `dysautonomia`; the demo cluster
  `autonomic_activity_invisible_illness` (members long_covid, me_cfs, pots, dysautonomia; `configs/conditions.yaml`);
  and the condition set named by the SPEC demo query, `long_covid_or_me_cfs` (members long_covid, me_cfs). Both
  sets use the same aggregation (section 3.7).
* Measurements (rows): bundles `wearable_autonomic_activity_monitoring` (primary), `nailfold_capillaroscopy`
  (the future CAPRIO adapter, run through the identical pipeline to demonstrate the swap),
  `autonomic_function_testing`, `exercise_capacity_testing`. Bundles resolve to member classes and implementer
  groups through `measure_it.facilities.matching.resolve_measurement` (configs/relevance.yaml, curated; since
  2026-09-24 a thin wrapper of the shared `measure_it.measurements.resolve.resolve_measurement`).
* Table `deployment_opportunities`: one row per (condition or set) x measurement x geography;
  `object_id = opportunity:<condition>|<measurement>|<geo_id>`.

**Primary query** for Tests 5-7: the SPEC demo query, `long_covid_or_me_cfs` x
`wearable_autonomic_activity_monitoring` x county, weight set `equal`. Every other combination is reported in
secondary tables.

## 3. Components (raw -> normalised)

Normalisation (configs/scoring.yaml): `pr(x) = rank(x, average ties) / n` over non-null universe rows, the same
function the geography module uses for the desert. A component that is itself an average of percentiles is
re-normalised with `pr()` before weighting. Every raw input, every sub-percentile, every component and every
normalised component is a column.

### 3.1 burden(g,c)
The primary measure of `geo_condition_features` (docs/BURDEN_DEFINITIONS.md), never substituted:
`burden_value`, `burden_evidence_level`, `burden_inherited`, `burden_source_resolution`, CI, multiplier.
* long_covid county: the state HPS value inherited (level A of the **state**, `inherited=True`). Flagged in
  `uncertainties` on every such row ("burden is the state estimate inherited; no within-state variation").
* me_cfs: CMS 51 2022 (level C symptom proxy); 24 counties have no value (level D for that row; 26 before the
  2026-09-24 FIPS-rename bridge of the CMS ingestion, deviation 7 update).
* pots, dysautonomia: level D everywhere.
* **Level D** rows: burden is null, `burden_excluded_level_D = True`, the burden weight is removed and the
  remaining weights renormalised to sum to 1 (`burden_weight_renormalised`), and the `burden_only` rank is null.
  They are never ranked on burden.
* normalised: `burden_pct = pr(burden_value)`.

### 3.2 vulnerability(g)
`svi_overall` (CDC/ATSDR SVI 2022; state value = population-weighted county mean, labelled derived upstream).
Normalised `vulnerability_pct = pr(svi_overall)`. Not burden.

### 3.3 diagnostic_desert(g,c)
The geography module's pre-specified index, taken as is: `mean(pr(burden), pr(SVI), 1-pr(relevant providers per
100k), 1-pr(relevant trials in the county or within 50 km))` with its four sub-percentiles kept. Level D:
`diagnostic_desert_access_only` is used and labelled `diagnostic_desert_variant = access_only (level D; not
burden-informed)`. For the long-COVID county index the burden component is the inherited state value
(`desert_burden_resolution` carried). Normalised `diagnostic_desert_pct = pr(diagnostic_desert)`.
Known overlap (reported, not corrected, because SPEC defines it this way): the desert re-uses burden and SVI, and
its access terms enter with the opposite sign to clinic_capacity and research_readiness.

### 3.4 clinic_capacity(g,m): implementation capacity for the measurement
Four raw sub-measures, each with a count and a per-capita rate:
* `cc_impl_providers_in_geo`: distinct individual NPIs (NPPES entity type 1) with >= 1 taxonomy in one of m's
  implementer specialty groups, primary practice in the county, per 100k ACS population;
* `cc_impl_providers_within_50km`: the same NPIs within 50 km (`clinic_matching.default_radius_km`) of the county
  internal point, per 100k population within 50 km (`population_within_50km`, geography module);
* `cc_impl_facilities_in_geo`: geocoded clinic-candidate facilities (`facilities`, `is_clinic_candidate`, the
  same base as `find_candidate_clinics`) with >= 1 implementer group (specialty or facility type; HRSA FQHC /
  look-alike sites count as `fqhc`), in the county, per 100k;
* `cc_impl_facilities_within_50km`: the same within 50 km, per 100k population within 50 km.
`clinic_capacity = mean(pr(4 rates))`; `clinic_capacity_pct = pr(clinic_capacity)`. States: the two in-state
rates only. Implementer groups are curated assumptions (configs/relevance.yaml) and say nothing about whether a
facility runs the measurement.

### 3.5 research_readiness(g,c,m): study capacity
Pool = geocoded facilities in the county or within 50 km of its internal point (the `find_candidate_clinics` pool);
state = facilities in the state. Raw sub-measures (distinct counts over the pool):
* `rr_condition_trials_pool_n`: distinct NCT ids registered at pool facilities whose record names the condition
  literally (`facility_trials.condition_ids_literal`; MeSH-expanded view not used);
* `rr_condition_nih_core_pool_n`: distinct NIH core projects of pool facilities, precision view
  (`facility_nih_projects.condition_ids_precision`);
* `rr_tech_trials_pool_n`: distinct NCT ids at pool facilities (any target condition) whose registered text
  mentions a member class of m (technology experience, as in `find_candidate_clinics`).
Reported, not scored: `rr_condition_trial_sites_pool_n` (facilities with >= 1 condition trial),
`rr_tech_trials_condition_pool_n`. `research_readiness = mean(pr(3 counts))`; normalised with `pr()`.
Counts are absolute (a study site is a site). Research readiness is not burden.

### 3.6 technology_saturation(g,m): reported, not in the default formula
Facilities in the pool with technology experience of m (>= 1 registered trial mentioning a member class), per
100k population within 50 km (state: per 100k state population). It is a research-deployment proxy for existing
use of the measurement, not market penetration (no public source measures device deployment by place).
`technology_saturation_pct = pr()`. Used only in the sensitivity weight set `saturation_adjusted` (section 6).

### 3.7 Condition sets (demo cluster; long_covid_or_me_cfs): documented aggregation
* burden: `burden_pct` = mean of the member `burden_pct` over members with a non-D value; `burden_value` is null
  for a set (members' raw values are kept in `member_components` JSON); set evidence level = the least direct
  contributing level (A < B < C); `burden_inherited` = any contributing member inherited; D members listed in
  `burden_members_excluded_level_D`.
* diagnostic_desert: mean of the member deserts over members with a burden-informed desert; access-only mean
  only if no member has one (labelled).
* clinic_capacity, technology_saturation: condition-independent (same as the single conditions).
* research_readiness: distinct trials / NIH core projects over the union of member conditions.
* Every member's own row stays in the table (the members' components are kept).

## 4. Composite and ranks

`composite_<ws> = sum_k w_k * normalised_k` over burden_pct, vulnerability_pct, diagnostic_desert_pct,
clinic_capacity_pct, research_readiness_pct, with the weight sets of configs/scoring.yaml (`equal` = default),
renormalised where burden is excluded. `rank_<ws>`: 1 = highest composite among eligible rows; ties broken by
geo_id (deterministic; tie counts at the top-N boundary reported where they matter). `composite`/`rank` = equal.

## 5. Uncertainty propagation (Monte Carlo)

1,000 seeded draws (`config.SEED`), per condition x measurement x level. In each draw:
1. **Burden** per member condition, `burden_draw = value + e` truncated at 0, `e ~ Normal(0, sd)`:
   * published CI (HPS, level A): `sd = (ci_high - ci_low) / (2 * 1.96) * multiplier(level)`;
   * no published CI (CMS 51, level C): `sd = sqrt(p(1-p)/n * 1e4 + 1/12) * multiplier` in percentage points,
     where n = the geometric midpoint of the published Medicare FFS size class (11-499 -> 74; 500-999 -> 707;
     1,000-4,999 -> 2,236; 5,000-9,999 -> 7,071; 10,000+ -> 10,000) and 1/12 is the variance of integer
     rounding; `multiplier` from `scoring.yaml evidence_level_uncertainty_multiplier` (A 1.0, B 1.25, C 2.0);
   * **inherited county values**: the state error is drawn once per state and shared by its counties, plus an
     independent per-county deviation `delta ~ Normal(0, tau)`, where `tau` = the method-of-moments between-state
     SD of the HPS state estimates (`sqrt(max(0, var(state values) - mean(se^2)))`). Assumption, stated in
     `uncertainties`: counties vary around their state value as much as states vary around the nation.
2. Recompute `burden_pct`, the desert (its burden percentile replaced by the drawn one; the other three
   sub-percentiles fixed) and `diagnostic_desert_pct`.
3. **Weights** `w ~ Dirichlet(20 * w_default)` (mean = default; sd about 0.087 per weight at the equal default);
   level D rows drop the burden weight and renormalise.
4. Composite and rank among eligible rows.
Per row: `rank_mc_p05`, `rank_mc_p50`, `rank_mc_p95` (the 5th-95th percentile rank interval), `p_top10_mc`,
`p_top25_mc`; plus two decompositions with one source switched off: burden noise only (default weights) and
weights only (observed burden): `rank_mc_burden_p05/p95`, `rank_mc_weights_p05/p95`.

## 6. Tests

### Test 5: does facility/research information change where we would deploy?
Primary query; secondary: every condition/set x measurement x level. Contrasts, each vs `equal` (full):
* **T5-1 (pre-specified primary, as in SPEC):** `burden_only` vs full;
* T5-2: need only (burden + vulnerability, 0.5/0.5) vs full: no facility or research information at all;
* T5-3: burden + vulnerability + desert (1/3 each) vs full: isolates adding clinic_capacity + research_readiness.
Metrics (eligible regions): top-10 and top-25 overlap, Jaccard, expected overlap under independence N^2/M,
Kendall tau-b and Spearman over all eligible regions; tie-averaged overlap (1,000 random tie-breaks) where the
contrast ranking has ties at the boundary; regions entering and leaving the top-N with their component
contributions `w_k * (x_k - mean_k)` as the reason. Reading (pre-specified, descriptive): facility/research
information **materially changes** the short list if the top-25 Jaccard of T5-1 (and of T5-3) is < 0.5. A change
is not evidence of a better list; Test 6 asks whether it is structure or noise.

### Test 6: negative controls for the scoring
(a) **Geography-label shuffles.** Permute the burden values across geographies (all members of a set with the
same permutation), nationally (200 draws) and within state (200 draws); recompute burden_pct, desert, composite,
rank. Separately permute SVI (vulnerability and the desert's SVI term) the same way. Statistics per draw (eligible
rows): Spearman of the shuffled composite with the observed composite; top-25 overlap with the observed list;
Spearman of the shuffled composite with the **real** (unshuffled) component (S). Reading: the ranking loses burden
(vulnerability) structure if the observed S exceeds the 97.5th percentile of the shuffled S. This is close to true
by construction when the component has weight; its size (how much of the list changes) is the informative part.
The within-state shuffle of an inherited long-COVID value is expected to change nothing (no within-state
information) and is reported as such.
(b) **Clinic-location shuffle.** The facilities contributor's shuffle (`measure_it.facilities.shuffle_test`)
permutes location tuples of geocoded candidate facilities and re-runs the matcher on 201 sampled counties; it does
not recompute scores. Here the same permutation scheme is applied to all geocoded clinic-candidate facilities,
and (for the provider sub-measures) to the locations of individual NPIs, nationally (200 draws) and within state
(100 draws); clinic_capacity, research_readiness and technology_saturation are recomputed for every county, then
composite and rank. Statistics: Spearman of shuffled with observed composite; top-25 overlap; Spearman of shuffled
with observed clinic_capacity and research_readiness (how much of each component's geography survives a shuffle
that keeps the set of occupied locations); S as in (a) for each component. Reading as in (a). Because the
permutation keeps location slots, per-capita facility density is preserved by design: if a component does not
change under the shuffle, it measures facility density rather than measurement-specific capability, and that is
reported as a finding.
(c) **Phenotype labels.** Cite, from the files that hold them, the label-permutation nulls of
`results/tables/nhanes_permutation_null.csv`, `nhanes_shuffle_control.csv`,
`stanford_uwakwe_permutation_nulls.csv` and `stanford_acute_detection.csv`.

### Test 7: ranking robustness
Per condition/set x measurement x level: ranks under the 7 configured weight sets plus `saturation_adjusted`
(equal weight 1/6 on the five components and on `1 - technology_saturation_pct`; defined in code because
`configs/scoring.yaml` is not ours to edit), and under 1,000 flat Dirichlet(1,1,1,1,1) random weight vectors.
Metrics: pairwise Kendall tau-b between named sets; Kendall's W over named sets (with and without `burden_only`)
and over random draws; per region min/max/range/IQR of rank across named sets and 5th-95th percentile across random
draws; fraction of named sets and of random draws in which the region is top-10 / top-25.
**Unstable region (pre-specified):** in the top-25 under at least one named set AND (top-25 in < 50% of named sets
OR top-25 in < 50% of random draws). Written to `results/tables/test7_unstable_regions.csv` with the reason.

## 7. Pre-specified sensitivity analyses (primary query unless stated)
* S1: Monte Carlo with tau = 0 (pure inheritance; no county deviation).
* S2: long_covid county burden replaced by the PLACES PHLTH C-proxy (`burden_proxy_value`, level C, multiplier
  2.0), for long_covid and the sets that contain it: Spearman and top-25 overlap with the primary ranking.
* S3: `saturation_adjusted` weight set (inside Test 7).
* S4: clinic_capacity without primary care (implementer groups minus `primary_care`), because primary care
  dominates NPPES counts.
* S5: the adapter swap: rank correlation and top-25 overlap between measurement bundles for the primary condition
  set (same pipeline, different measurement).

## 8. Recommendation output
`rank_deployment_opportunities(condition, measurement, geography_level='county', top_n=10, weight_set='equal')`
reads `deployment_opportunities` (computes on the fly with the same engine for combinations not in the table) and
returns per region the SPEC schema: condition; phenotype (`get_patient_phenotype_signature` per member);
measurement; technology {name, regulatory_context (`get_regulatory_context` per member class),
technology_evidence (`get_measurement_evidence`)}; geography {name, fips, burden, burden_evidence_level,
vulnerability, diagnostic_desert, + all components, rank interval}; candidate_sites (`find_candidate_clinics`
per member condition, with reasons and the framing sentence); research_evidence (trial:/nih: ids in the pool);
molecular_context (`get_molecular_context` summary, labelled condition-level molecular enrichment); uncertainties;
provenance (object ids + sources + versions); recommended_next_step (a candidate deployment opportunity / pilot
evaluation, never a validated diagnostic pathway). Aliases: 'Long COVID or ME/CFS' -> `long_covid_or_me_cfs`;
'wearable autonomic monitoring' -> the wearable bundle. Unresolvable inputs return UNKNOWN with a reason.

## 9. What is not done
No fitted model and no tuning on outcomes: all weights, radii, thresholds, draws and seeds are fixed above. No
county prevalence is created from a state value (the inherited value is carried with its flag; the Monte Carlo
deviation widens uncertainty, it is not an estimate). Nothing is joined to person-level data; phenotype and
molecular context are attached at the condition level only. No claim that a ranked region or site will produce
evidence; the output is a candidate deployment opportunity.

## Deviations (added after the analysis ran)

1. **2026-09-23, per-capita denominators (build QA, before any test statistic was read).** A few counties have
   `population_within_50km` of 0 or missing, which produced infinite per-100k rates. Rates with a non-positive or
   missing denominator are now NaN, and `clinic_capacity` is the mean of the AVAILABLE sub-percentiles (plan: mean of
   all four). Affects only counties without a within-50-km population value.
2. **2026-09-23, Test 6 scope widened (secondary, added).** The label shuffles (6a) were run for all six conditions/sets
   at county level and for the primary query at state level (plan: primary query). The clinic-location shuffle (6b)
   was run for the primary condition set with all four measurement bundles (plan: primary query). The primary rows
   are unchanged by this; the extra rows are secondary.
3. **2026-09-23, condition-set label.** The SPEC demo query set is labelled "Long COVID or ME/CFS" (the working label
   "Long COVID or ME/CFS (SPEC demo query)" read badly inside recommendation sentences). Presentation only.
4. **2026-09-23, descriptive tables added after the first results were read (not pre-specified; no test depends on
   them).** `test5_component_correlations.csv` (Spearman between the normalised components, to quantify the
   desert/capacity/readiness overlap stated in section 3.3) and `test5_readiness_zero_counts.csv` (how many regions
   have no condition trial and no NIH project in the pool, and the research_readiness percentile those regions
   still receive through technology-experience trials). Both were added because the first read of the demo output
   showed top-ranked counties with zero condition-specific research activity but a high research_readiness percentile.
5. **2026-09-23, Test 5 reason attribution for level D.** For level-D rows the burden contribution is undefined; it
   is now skipped when choosing the "largest contributions" in the reason text (the first run's text could list a
   NaN contribution). Statistics unchanged.
6. **2026-09-23, reviewer corrections (independent review, after all results existed).**
   * *Set burden re-normalised.* Section 3 requires every component that is an average of percentiles to be
     re-normalised with `pr()` before weighting (and `configs/scoring.yaml` rescales every component by percentile
     rank). The first build did this for desert, clinic_capacity and research_readiness but used the raw mean of the
     member burden percentiles as a set's `burden_pct`. Its spread among eligible counties was SD 0.213 against about
     0.28 for every other component, so under "equal" weights burden carried about 25% less than its nominal share.
     Now `burden_set_member_mean_pct` is the raw set burden (section 3.7) and `burden_pct = pr()` of it, in the
     observed ranking, in every Monte Carlo draw and in every shuffle. Single conditions are bit-identical
     (pr(pr(x)) = pr(x)). Effect on the equal-weight rankings (old vs new, eligible rows): primary query (county)
     top-10 membership unchanged (order changes: Hoke 8 -> 7, Wayne 7 -> 8), top-25 overlap 22/25, composite
     Spearman 0.994; demo cluster (county) top-10 overlap 8/10, top-25 22/25; state-level top-25 overlap 24-25/25.
     Every set number in `results/SCORING_RESULTS.md` was regenerated.
   * *Query-function edge cases.* `top_n` that is not a positive integer now returns UNKNOWN (a negative value used
     to return all but the last |n| regions); `include_small_population=True` now applies to every weight set (it was
     silently ignored except for `equal`); an all-level-D set no longer says "the set's burden uses the other
     members".
   * *Provenance.* `deployment_candidates` rows whose burden is the inherited state value carry
     `source_geographic_resolution = state` (CONVENTIONS 3.5; it was `county`); recommendation provenance lists the
     technology-experience trial ids as well as the condition-trial and NIH ids; the geography block gets
     `burden_evidence_level_note` when the evidence level describes an inherited state measure.
   * *Reporting.* Test 5's "leaving" reasons now say "smallest (most negative) contributions" (the second one could be
     positive) and mark rows whose set burden rests on one member; the results state that the 6(b) shuffle keeps the
     number of facilities per location (so it cannot move a density measure) and holds the desert's provider/trial
     terms at their observed values, so Test 5's change follows research readiness, not a measurement-specific
     clinic signal; the 12 eligible counties whose set burden rests on long COVID alone (no CMS value) are listed.
   * *Integration fixes (2026-09-23, after the review; not pre-specified).* (a) The `saturation_adjusted` weight set
     moved from `opportunity.EXTRA_WEIGHT_SETS` into `configs/scoring.yaml weight_sets` (same values). (b) An explicit
     `inherited_burden_uncertainty_multiplier` (2.0 = the level-C multiplier) now inflates the state-level error of an
     inherited county burden in the Monte Carlo, on top of the evidence-level multiplier and the per-county deviation
     tau: a state value carried to a county has no within-state information and is treated like a proxy for the
     county. New sensitivity S1b (multiplier 1.0 vs 2.0, primary query,
     `results/tables/test7_sensitivity_s1b_inherited_multiplier.csv`). Every Monte Carlo number of a combination with
     long COVID as a member was regenerated.
7. **2026-09-24, incomplete burden: the set burden is computed only where it is complete (deviation from section
   3.7; changed after results existed).** Section 3.7 averaged the member burden percentiles over the members with a
   non-D value *in each region*. Where a member that has a burden measure at the level had no value in a region
   (me_cfs at county level: no CMS value for the 9 Connecticut planning regions, which CMS publishes only on the
   legacy counties, Oglala Lakota SD, Manassas Park VA, 14 Alaska boroughs/census areas and Kalawao HI: 26 counties, 12 of them
   ranking-eligible), the set burden of that row rested on the other member alone: a different quantity from its
   neighbours', ranked as if it were the same. For `long_covid_or_me_cfs` it is the inherited long-COVID state value,
   which puts every county of a high-HPS state at the same high percentile; that made Greater Bridgeport Planning
   Region CT a top-10 county of the primary query (rank 9) and of the demo cluster (rank 10), and Oglala Lakota SD rank
   1 of the burden-only list. The positions were artefacts of the aggregation rule, which the first results already
   said. **Rule now:** a member has a *defined* burden at a level if its burden measure exists there
   (`burden_defined_level` != D in `geo_condition_features`; long_covid and me_cfs at county and state; pots and
   dysautonomia never). A condition's or set's burden is computed only in regions where EVERY member with a defined
   burden has a value (a single condition is a set of one); it is then the mean of those members' burden percentiles,
   re-normalised with `pr()` over the complete regions. Members with no defined burden (level D everywhere) are left
   out as before; if no member has a defined burden the burden weight is removed and the other weights renormalised,
   as before. Regions where a defined member has no value are **incomplete**: `burden_incomplete = True` with
   `burden_incomplete_reason` (why the region is not ranked) and `burden_missing_source_reason` (the source's reason for
   the missing value), burden, composite and every rank null, Monte Carlo not run; they
   are excluded from the condition's ranking and listed separately (`results/tables/scoring_incomplete_burden.csv`,
   the `excluded_incomplete_burden` block of `rank_deployment_opportunities`). Their non-burden components are kept.
   The geography-label burden shuffle (Test 6a) permutes burden values among complete regions only, so the ranked set
   is the same in every draw. Applying the rule to the single condition me_cfs (whose 26 counties without a CMS
   value were ranked on the other four components with the burden weight renormalised) is part of the same change:
   the same mixed-basis problem, the same rule. Why not impute or map: a CMS value cannot be carried from the legacy
   CT counties to the planning regions without a crosswalk that the geography module deliberately does not apply
   (CONVENTIONS section 4), and filling one member from another would be the same artefact. Not every gap is
   of that kind (reviewer note, 2026-09-24): Oglala Lakota County SD (46102) is the 2015 renaming of Shannon County
   (46113), and CMS MMD still publishes its value under 46113 (`geo_condition_burden`, `geo_vintage =
   cms_fips_not_in_2024`, 2012-2023). The geography module does not remap 1:1 renames, so its source reason ("no
   published value") is inaccurate and the county is incomplete because of a crosswalk gap; a 1:1 rename could be
   remapped without apportioning. That is a change for the geography module (reported to it); this ranking applies
   no mapping, and SCORING_RESULTS section 5 flags the case. The rule changes which
   regions are ranked, not the formula; before/after rankings are regenerated in every run from the old rule kept as
   an option (`make_combo(..., burden_rule="available")`) in `results/tables/scoring_set_burden_rule_before_after.csv`
   and `scoring_set_burden_rule_top10.csv`.
   **Effect (run 2026-09-24, all numbers regenerated in `results/SCORING_RESULTS.md`).** County rows ranked for me_cfs
   and the two sets: 2,394 (was 2,406); no state row is incomplete, and every combination without an incomplete region
   (long_covid, pots, dysautonomia; every state ranking) is identical under both rules. Primary query (equal weights):
   top-10 overlap 9 of 10 (Greater Bridgeport Planning Region CT, rank 9, is excluded; Lauderdale County MS moves from
   11 to 10), top-25 overlap 24 of 25, Kendall tau-b 0.998 over the 2,394 counties ranked under both rules. Demo
   cluster x wearable: 9 of 10 (Greater Bridgeport, rank 10, excluded; Forrest County MS 11 -> 10). Oglala Lakota County
   SD was rank 1 of the burden-only list of both sets. me_cfs alone: top-25 unchanged in all four measurements
   (tau 1.000; its incomplete counties ranked 238 and below x wearable, 189 and below over the four measurements).
   Across the 12 affected county combinations top-10 overlap is 9-10 and tau 0.998-1.000; in one (long_covid_or_me_cfs x nailfold_capillaroscopy) the re-normalisation of the set burden over the
   complete counties alone swaps the counties at ranks 10 and 11. Consequences elsewhere: Test 6a's burden shuffle
   permutes among complete counties only (long_covid alone reproduces its earlier draws exactly; me_cfs's and the sets' null means
   move by less than 0.004); S2 compares over the counties ranked under both rankings, and the 188 counties without a
   PLACES PHLTH value (KY 120, PA 67, TX 1) are incomplete under S2 instead of ranked with the burden weight
   renormalised.
   **Update (2026-09-24, later the same day; a data fix, not a method change).** The reviewer's rename note was acted on
   at the source layer: the CMS ingestion now bridges the two 1:1 FIPS renames CMS still publishes (46113 -> 46102
   Oglala Lakota SD, 02270 -> 02158 Kusilvak AK; `cms_mmd.FIPS_RENAMES`, docs/BURDEN_DEFINITIONS.md), and
   `geo_condition_features.burden_unknown_reason` now states the measured source reason (for the remaining CMS gaps:
   absent from every MMD file ingested) instead of "no published value ... (suppressed, cell omitted by the source,
   or entity absent from the source)". Both counties now have a ME/CFS value. Incomplete counties: 24 (was 26), 11
   ranking-eligible (was 12); county rows ranked for me_cfs and the two sets 2,395 (was 2,394). Oglala Lakota SD
   (population 13,491; CMS 51 prevalence 23 %, the county median) is complete and ranked: equal-weight rank 443,
   burden-only rank 307 of the primary query; Kusilvak AK is complete but below the population minimum. Primary
   query before/after the completeness rule: top-10 overlap 9 of 10, top-25 24 of 25, Kendall tau-b 0.999 over 2,395
   counties (0.998 over 2,394 before the bridge). Effect of the bridge on the stored rankings (equal weights, all 48
   combinations): top-10 membership unchanged everywhere (two adjacent swaps inside a top-10: demo cluster x autonomic
   testing ranks 6/7, me_cfs x wearable ranks 8/9); top-25 membership changes in 2 combinations (demo cluster x
   capillaroscopy: Wayne County TN out, Pickens County AL in; demo cluster x wearable: Palm Beach County FL out,
   Tippah County MS in); the primary query's top-25 is unchanged. The effect figures in the paragraph above are those
   of the first run of the rule.
8. **2026-09-24, sensitivity S6 `specialist_only` (declared before any of its statistics was computed; not part of
   the 2026-09-23 plan).** Primary care is 86-99% of the "relevant" individual NPIs of every condition (every core
   specialty set includes `primary_care`; results/GEOGRAPHY_RESULTS.md section 6), so relevant-provider density
   mostly measures primary-care supply. S6 removes `primary_care` from relevant-provider density everywhere the
   composite uses it: (i) the diagnostic desert's provider term becomes `1 - pr(relevant specialists per 100k)`, the
   distinct individual NPIs (NPPES entity type 1, primary practice location in the geography) with >= 1 taxonomy in
   the condition's core groups minus `primary_care` (per member condition for a set; the access-only variant for
   level D likewise), counted in this module from `providers` and checked against the geography module's
   `relevant_specialists_n` and, for the primary-care group, against `provider_density_county`; (ii)
   clinic_capacity uses the measurement's implementer groups minus `primary_care` (as S4; facility-type groups such as
   `fqhc` and `community_health_clinic` stay). Everything else is unchanged. Reported for the primary query (county;
   also state) and, secondarily, every condition x measurement at county level: top-10 and top-25 overlap, top-25
   Jaccard, Kendall tau-b and Spearman of the composite over eligible regions vs the primary (equal-weight) ranking,
   plus the two single swaps (desert only; clinic_capacity only = S4) to show which part moves the list. Reading
   (declared): descriptive; no threshold decides anything.
   **S6 result (deviation 8; SCORING_RESULTS section 8).** Primary query, county: top-10 overlap 5 of 10, top-25 17
   of 25 (Jaccard 0.515), Kendall tau-b 0.881; desert provider term alone 7 / 18 / 0.901; clinic capacity alone (= S4)
   5 / 17 / 0.854. State: 10 / 25 / 0.879. All 24 county combinations: top-10 overlap 4-8, tau 0.810-0.906. The
   specialist counts computed from `providers` equal the geography module's in every condition x geography row and
   `provider_density_county` in every group x county row checked.
