# Analysis plan: geographic burden, context and opportunity features (SPEC §E, Phase 4, Test 4, Test 6 geography labels)

Status: **pre-specified 2026-09-23, before any test statistic was computed.** Only schemas, row counts and
the DATA_AUDIT files of the inputs had been read when this was written. Deviations made later are listed
at the end, dated, with the reason; nothing above that section is edited after the analysis runs.

Owner modules: `measure_it.geography.features` (build + tests 4/6), `measure_it.geography.query` (MCP-facing queries).
Burden measure definitions (the proxy map): `docs/BURDEN_DEFINITIONS.md`, mirrored in code as
`features.BURDEN_DEFINITIONS` and written to the table `geo_burden_definitions`.

## 1. Questions

* **Q1 (build).** Put every public burden estimate for the 14 target conditions into one long table with an
  explicit evidence level (A/B/C/D), source resolution and inherited flag, and pick ONE primary measure per
  condition x geography level before looking at any result.
* **Q2 (Test 4).** Are geographic "hotspots" driven by burden or by population size? Compare rankings by
  prevalence/incidence, by absolute counts (only where a count is defensible at the same resolution) and by
  population alone; do the same for the diagnostic-desert index.
* **Q3 (Test 6, geography labels).** If burden values are shuffled across geographies (nationally, and within
  state), do (a) relationships with independent variables and (b) spatial autocorrelation collapse? Which
  relationships survive a within-state shuffle (i.e. are carried by between-state differences)?

## 2. Units and universes

* Geography levels in the feature table: `state` and `county` (2024 vintage; Connecticut = 9 planning regions).
  Rows are written for all 56 state-equivalents and 3,235 canonical counties; legacy CT counties are not
  feature rows (their CMS values stay in `geo_condition_burden`, flagged `matched_ct_legacy`).
* **Analysis universe** for every percentile, ranking and test: the 50 states + DC (51 states; 3,144 counties).
  Puerto Rico and the island areas are excluded from the universe because SVI ranks Puerto Rico separately
  (not poolable, per the SVI audit), PLACES and HPS do not cover them, and ACS/SVI do not cover the island areas.
  Their feature rows carry `diagnostic_desert = null` with a reason.
* Population: ACS 2020-2024 5-year `total_population` (per-100k denominators, Test 4 population);
  adult population = `n_age_18_64 + n_age_65_plus` (long-COVID counts).

## 3. Burden (primary measures, pre-specified; rationale in BURDEN_DEFINITIONS.md)

| condition | state primary (level) | county primary (level) | county C-proxy column |
|---|---|---|---|
| long_covid | HPS `lc_current_pct_all_adults`, latest non-suppressed period per state (A) | state value **inherited** (A, `inherited=True`, source resolution `state`) | PLACES PHLTH crude (C) |
| me_cfs | CMS MMD 51, 2022 final, unsmoothed actual, all ages (C) | same, county (C) | PLACES PHLTH crude (C) |
| fibromyalgia | CMS MMD 51 as above (B) | same (B) | PLACES ARTHRITIS crude (C) |
| migraine | CMS MMD 59 as above (B) | same (B) | none |
| lyme_disease | 2023 reported cases / population, aggregated from counties (A) | 2023 incidence per 100k (A) | none |
| ptlds | Lyme 2023 incidence as antecedent-exposure proxy (C) | same, county (C) | none |
| pots, dysautonomia, eds_hsd, mcas, post_infectious_syndrome, ibs, gastroparesis, endometriosis | none (D, null) | none (D, null) | none |

Pre-specified alternates (sensitivity only; never substituted silently): HPS pooled mean of 2024 periods; HPS strict
latest period; CMS 2022 age-standardized; CMS 2023 preliminary; Lyme pooled 2022-2023; Lyme 2019 (pre-2022 case
definition); PLACES age-adjusted versions; PLACES COGNITION and DISABILITY for long_covid.

Burden counts (only same-resolution, defensible): long_covid state = HPS % x ACS adults 18+; lyme county/state =
reported cases; PLACES proxy county counts = crude % x PLACES `total_pop_18plus` (labelled `modeled`). No counts for
CMS measures (FFS denominators are published only as size classes) and none for inherited county values.

Uncertainty multiplier from `configs/scoring.yaml` (A 1.0, B 1.25, C 2.0, D null) is attached to every burden value.

## 4. Features (one row per geo_id x condition; `object_id = geo:<fips>`)

* Vulnerability: SVI 2022 overall + 4 theme percentiles (county). State: population-weighted mean of county SVI
  percentiles, labelled as derived (no official state SVI exists).
* Relevant provider density: **primary variant = distinct individual NPIs (entity type 1) whose ANY taxonomy
  code falls in one of the condition's `condition_specialties.core` groups (configs/relevance.yaml), placed at the
  primary practice ZIP's county, per 100,000 ACS residents.** Variants kept as columns: specialists only (core groups
  minus primary_care); providers within 50 km of the county internal point per 100k population within 50 km (ZCTA
  internal points, ACS ZCTA population).
* HRSA: active service-delivery sites (incl. admin/service) per 100k; HPSA/MUA flags from `geo_context__hpsa`.
* Trials (literal-precision view: `trial_conditions.condition_literal_match = True`): distinct relevant trials with
  >= 1 U.S. site in the geography; recruiting (`overall_status = RECRUITING`); active (RECRUITING, NOT_YET_RECRUITING,
  ENROLLING_BY_INVITATION, ACTIVE_NOT_RECRUITING); completed; and for counties distinct relevant trials / recruiting
  trials with a U.S. site within 50 km and 100 km (haversine) of the county internal point. Site coordinates are
  ClinicalTrials.gov city centroids (audit), so radii are approximate.
* NIH (precision view: `match_tier = title_abstract` and not `likely_false_positive`): distinct appl_ids, core projects,
  active core projects and organizations geocoded to the geography.

### Diagnostic desert (pre-specified)
Within the analysis universe, per condition and level, with percentile rank `pr()` (average ties, in [0,1]):

    diagnostic_desert = mean( pr(burden), pr(SVI overall), 1 - pr(relevant providers per 100k),
                              1 - pr(relevant trials within 50 km) )            # county
    (state: trials with a U.S. site in the state instead of the 50 km count)

All four component percentiles and the raw components are kept as columns. Null (with `desert_unknown_reason`)
when any component is missing (level-D burden, suppressed/missing burden, outside the universe).
`diagnostic_desert_access_only` = mean of the three non-burden components (defined for level-D conditions too;
labelled as not burden-informed). Pre-specified sensitivity variants: specialists-only density; within-50-km
provider density; trials within 100 km; recruiting trials within 100 km; long_covid county burden = PLACES PHLTH
proxy instead of the inherited state value. Each variant is compared to the primary by Spearman rho and top-25 overlap.

## 5. Test 4 — hotspots vs population size

Series (all pre-specified): lyme_disease county (incidence, cases); lyme_disease state; long_covid state
(prevalence, adult count); long_covid county C-proxy PHLTH (crude prevalence, modeled count); CMS 51 county
(fibromyalgia B / me_cfs C, one series because the values are identical); CMS 59 migraine county; and the primary
diagnostic desert for every condition/level where it is defined.

Metrics per series:
1. Spearman rho of burden (prevalence/incidence) with total population; of count with population; of desert with
   population. 95% CI by percentile bootstrap (2,000 resamples, seed `config.SEED`).
2. Top-N overlap between rankings by prevalence, by count and by population (and desert vs population):
   N = 25 counties (`scoring.yaml ranking.top_n`), N = 10 states; overlap count, Jaccard, expected overlap under
   independence N^2/M, hypergeometric P(overlap >= observed). County rankings exclude counties below 10,000
   residents (`ranking.min_population`); sensitivity N = 100 and no population floor.
3. Partial Spearman (rank-based partial correlation) controlling for log population: count vs prevalence;
   desert vs burden; burden vs SVI.
4. OLS R^2 of log(count) on log(population) and on log(prevalence) (counts > 0 only).

Pre-specified reading: a ranking is called **population-driven** if rho with population > 0.7 and its top-N overlap
with the population ranking exceeds chance at P < 0.05; **burden-driven** if |rho| < 0.3 and the top-N overlap
is not above chance. Anything else is reported as mixed.

## 6. Test 6 — geography-label negative controls

Permutations: B = 1,000 (correlations) / 999 (Moran's I), seed `config.SEED`, two schemes: (i) national shuffle
of the burden variable across geographies; (ii) within-state shuffle (counties only). Two-sided permutation
P = (1 + #{|null| >= |observed|}) / (1 + B). "Collapse" = the permuted statistics centre on the no-structure value
(rho = 0; Moran's E[I] = -1/(n-1)) and the observed value lies outside the central 95% of the null.
Retained fraction under the within-state shuffle = mean(null) / observed.

Cross-source relationships (Spearman, counties/states in the universe with both values):
* R1 CMS 51 (claims) vs PLACES ARTHRITIS crude (BRFSS model; the BRFSS arthritis item names fibromyalgia).
* R2 CMS 51 vs PLACES PHLTH crude.
* R3 CMS 59 migraine vs PLACES PHLTH crude.
* R4 positive control (same source, different years): Lyme incidence 2023 vs 2019.
* R5 state: HPS long COVID (primary) vs PLACES PHLTH aggregated to states (pop 18+-weighted). This is also the
  only available validation of the long_covid county C-proxy.
* R6 state: HPS long COVID vs CMS 51 state prevalence.
(State relationships use the national scheme only.)

Spatial autocorrelation (counties): Moran's I with (a) contiguity weights from the 2024 1:5M county polygons
(polygons touching after a 100 m buffer in EPSG:5070), (b) a 100 km distance band between internal points; both
row-standardised; counties without neighbours dropped. Variables: log1p Lyme incidence 2023; CMS 51; CMS 59;
PLACES PHLTH; PLACES ARTHRITIS; SVI overall; log1p relevant-provider density (long_covid core groups); primary
diagnostic desert (me_cfs); and the inherited long_covid county value (expected to be unchanged by the within-state
shuffle, demonstrating that inherited values carry no within-state information).

## 7. Controls, leakage and what is not done

* No tuning on outcomes: every threshold, radius, year and weight above is fixed here; there is no fitted model.
* Ecological only: nothing is joined to person-level data; no county prevalence is created from a state estimate.
* Proxies are never promoted to primary after seeing results.
* Not done: small-area estimation/downscaling of HPS; smoothing of Lyme rates; any causal claim.

## 8. Outputs

`geo_condition_burden`, `geo_context`, `geo_condition_features`, `geo_burden_definitions` (data/processed);
`results/tables/test4_*.csv`, `test6_geography_*.csv`, `geo_*.csv`; draft figures in `results/figures/drafts/`;
`results/GEOGRAPHY_RESULTS.md` (every number read from those files).

## Deviations (added after the analysis ran)

1. **2026-09-23, desert trial component (build QA, before any test statistic was read).** The pre-specified county
   component "relevant trials within 50 km of the county internal point" undercounts trials sited inside large
   counties whose internal point is far from the population centre (San Diego County: 8 relevant long-COVID trials
   sited in the county, 3 within 50 km of its internal point, which lies ~48 km from the city). The primary component
   is now **distinct relevant trials with a U.S. site in the county OR within 50 km of its internal point**
   (`trials_in_geo_or_within_50km_n`); the 100 km and recruiting variants use the same union. The pre-specified
   internal-point-only version is kept as the sensitivity variant `diagnostic_desert_sens_trials_50km_internal_point_only`
   and compared to the primary in `results/tables/geo_desert_sensitivity.csv`; the number of counties affected per
   condition is in `results/tables/geo_trial_radius_qa.csv`.
2. **2026-09-23, operational definition of "collapse" (not pre-specified numerically).** A statistic is labelled
   collapsed when the observed value lies outside the central 95% of the permutation null AND the null mean is within
   0.05 (Spearman) / 0.02 (Moran's I) of the no-structure value. Retained fraction is computed as
   (null mean - no-structure value) / (observed - no-structure value), which equals mean(null)/observed for Spearman.
3. **2026-09-23, interpretation labels for the within-state shuffle (presentation only; the statistics are unchanged
   and were recomputed with the same seed).** The first run's single label set called R1-R3 "partially survives" even
   where the within-state-shuffle null mean exceeded the observed value. Labels are now scheme-specific: for the
   within-state shuffle, the null mean is reported as the between-state component and observed - null mean as the
   within-state component (`within_state_component`).
4. **Verdict wording for desert series.** The pre-specified label "burden-driven" means only "not population-driven"
   when applied to the diagnostic desert, whose components include provider and trial access.
5. **2026-09-23, de-duplication (presentation only).** The base Test 4 series (lyme_county, long_covid_state, ...)
   no longer repeat the diagnostic-desert rows, which are reported once per condition x level in the `desert_*` series.
6. **2026-09-23, bug fix in the primary-vs-alternate comparison table.** The CMS county comparisons first included
   legacy-CT and CMS `xx990` rows (n = 3,141 instead of the 3,118 canonical counties; 3,120 since the 2026-09-24 CMS FIPS
   rename bridge); the table is now restricted to
   canonical 2024 geographies. No test statistic depended on this table.

## Deviations added by the independent review (2026-09-23; statistics above are unchanged unless stated)

7. **Within-state permutation P.** The pre-specified P = (1 + #{|null| >= |observed|}) / (1 + B) is two-sided around
   the no-structure value. That is correct for the national shuffle, but not for the within-state shuffle, whose null
   is centred on the between-state component (R1 got P = 0.992, 0.994 in the current build, although its observed
   value lies below the within-state band; the null moved within Monte Carlo error when the 2026-09-24 CMS FIPS
   bridge added two county values to Test 6's shared random stream). Added `p_perm_vs_null_centre` = two-sided P around the null mean. Within-state rows of Test 6
   (relationships and Moran's I) are now read with it. `p_perm` is kept.
8. **Spatially adjusted inference for cross-source correlations (Test 6 R1-R6).** A national shuffle treats
   geographies as exchangeable, so its null is too narrow when both layers are spatially autocorrelated (all of them
   are; see Moran's I). Added the Clifford-Richardson-Hemon (1989) effective sample size on ranks (100 km distance
   classes up to 2,000 km for counties; 250 km up to 3,000 km for states), a t-test P at n_eff, Holm across R1-R6,
   and the minimum detectable |rho| at 80% power. Calibration: `crh_calibration` simulates independent smooth fields
   on real county points (`test6_geography_crh_calibration.csv`). Consequence: R3 (CMS 59 vs PLACES PHLTH) is
   significant only under the naive shuffle.
9. **Spatial dependence in Test 4 CIs.** County correlations also get a state-cluster bootstrap CI (resampling
   states with all their counties; 2,000 resamples, seed `config.SEED`). The i.i.d. CI is kept.
10. **Multiplicity.** Holm-adjusted P added for the Test 4 verdict family (`p_hypergeom_holm`), the small-population
    family, each Test 6 relationship scheme and each Moran's I scheme (`p_holm`).
11. **Small-population tail of Test 4.** The pre-specified top-N test compares only with the LARGEST-population
    ranking, so it cannot detect hotspots concentrated in small places. Added `test4_small_population_check.csv`:
    overlap of the top-N with the N smallest geographies, and the number of top-N geographies in the bottom
    population quartile (hypergeometric P, Holm). Also added `reviewer_reading` in `test4_verdicts.csv`. The
    pre-specified `verdict` column is unchanged.
12. **Ties in top-N lists.** CMS publishes integer percentages, so many counties tie at the top-N boundary (173 for
    CMS 59 at N = 25). The deterministic FIPS tie-break then decides membership. Added `overlap_tie_averaged` and
    `ties_at_boundary_a/_b` to `test4_topn_overlap.csv`, and tie-averaged bottom-quartile counts.
13. **Inherited burden inside the long-COVID county desert.** As pre-specified, the county long-COVID desert uses the
    inherited state value as its burden component. The first write-up nevertheless said this value is "never ranked
    as county burden", which was false. Added `desert_burden_resolution` ('state (inherited; no within-state
    variation)') and an explicit `desert_definition` for these rows. The results text is corrected.
    `query.get_geographic_context` attaches a `burden_value_note`.
14. **Query layer.** Condition resolution now also accepts a unique whole-word prefix ('Lyme'). Place names are
    accent-folded. 'US' resolves to the national level. A geography id at the wrong level, and a missing value
    for a specific geography, return the specific reason from `geo_condition_features.burden_unknown_reason`.
