# Analysis plan: temporal-holdout validation of the deployment hotspots (Test 8)

Status: **pre-specified 2026-09-25, before any holdout outcome was computed.** At the time of writing the following had
been read or counted, all of them inputs or data-availability facts, none of them a region-level outcome: the scoring
code (`measure_it.scoring.opportunity`, `metric_link`), the geography module's provider / trial / NIH counters
(`measure_it.geography.features.provider_counts`, `trial_counts`, `nih_counts`, `_desert_components`), the schemas of
`clinical_trials`, `trial_sites`, `trial_conditions`, `facility_trials`, `nih_projects`, `nih_project_conditions`,
`facility_nih_projects`, `providers`, `facilities`, `facilities__hrsa`, `geo_condition_burden(__cms_mmd)`, and three
national histograms (trial `start_year` over all 7,910 target-condition trials; NIH `fiscal_year` of the 9,342
RePORTER rows, FY2015-2026; NPPES `enumeration_date` year). No county or state was scored with past data and no
post-period count per region existed. Deviations made later are listed at the end, dated, with the reason; nothing
above that section is edited after results exist.

Owner module: `measure_it.scoring.temporal_holdout` (pipeline step `scoring_temporal_holdout`, after `scoring_report`,
before `canonical_unions`). Results: `results/TEMPORAL_HOLDOUT_RESULTS.md`, `results/tables/temporal_holdout_*.csv`,
processed table `temporal_holdout_rankings`, draft figure `results/figures/drafts/temporal_holdout_*.png`.

## 0. The question, and what it can and cannot answer

> Do regions the engine would have flagged using only data available at a past cutoff T turn out to be where
> invisible-illness research / clinical activity actually appeared after T?

There is no ground truth for "where deployment of a measurement is most valuable". What can be checked against later
public data is narrower:

1. **Research readiness (checkable).** Relevant trials, trial sites and NIH awards that appear after T are public,
   dated and geolocated. If the readiness side of the engine (research_readiness, clinic_capacity, and composites that
   weight them) points to places where study activity then appears, the engine's notion of "positioned to act" has
   external, temporal support. This is a **proxy outcome**: later research activity is evidence of research
   readiness, not of clinical value, diagnostic yield or patient benefit. Research also goes where money, prestige and
   large populations are; a hit does not show that deployment there would help anyone.
2. **Access gap (checkable only as gap persistence).** The desert side (diagnostic_desert, vulnerability) is built to
   point at places with LITTLE provider and trial activity. It should NOT predict new research activity; if it did,
   the desert would be mislabelled. The checkable statement is that flagged deserts remained under-served after T
   (their post-T activity rate is below the base rate). That says the gap did not close by itself; it says nothing
   about whether a deployment there would be valuable.
3. **Burden (checkable only as persistence of the input).** Whether the pre-T CMS proxy ranks counties the way the
   later CMS proxy does. This checks the stability of an input, not burden itself (the CMS measure is a level-C
   claims proxy for ME/CFS).
4. **Not checkable here:** whether a deployment would detect more cases, whether the measurement works in these
   places, patient benefit, cost, or whether the "right" places were chosen. No outcome in public data measures these.

Each component is therefore reported **separately with its pre-specified expected direction** (section 5.3), and the
composite results are read in the light of which components drive them.

## 1. Cutoff T

* **Primary T = 2021** (inputs dated on or before 2021-12-31; NIH fiscal year <= 2021). Reasons: (a) CMS MMD has a
  final 2021 claims year for the ME/CFS proxy (condition 51, "Fibromyalgia, Chronic Pain and Fatigue"), the same
  measure as the current primary burden (2022); (b) NIH RePORTER data here start at FY2015, so T = 2021 leaves 7
  pre-period fiscal years; (c) Long COVID research existed but was young (first trials registered 2020-2021; ICD-10
  U09.9 from October 2021), so the pre-period contains some Long COVID signal while most Long COVID activity is in
  the post-period; (d) the post-period (2022-01-01 to the retrieval date, 2026-09) is about 4.7 years, long enough for
  new trials and awards to appear.
* **Sensitivity T = 2019 and T = 2020** (same rules). T = 2019 is pre-COVID: there is no Long COVID activity before
  T, and the post-period contains the whole Long COVID research wave (RECOVER and others), which is a harder test.
* **Rolling check:** the T = 2019 and T = 2021 rankings are compared with each other and with the current ranking
  (section 7).

## 2. Pre-period inputs (only data dated <= T), and how each is filtered

The engine is **not forked**. `measure_it.scoring.opportunity.make_combo` / `combo_frame` / `evaluate` and
`metric_link.reach_for` run unchanged on a *dated input view*: a copy of the spatial base and of the
`geo_condition_features` rows in which every input dated after T is removed. The view replaces, for the duration of
the run, `opportunity.spatial()` (facility / provider base), `opportunity._features()` (burden and desert inputs) and
clears the engine's input caches (`_item_incidence`, `metric_link._zcta_covered`). A `T = None` view (no filtering)
must reproduce the stored `deployment_opportunities` composites of the primary query exactly (tested); that is the
guarantee that the holdout scores come from the same formula.

| input | pre-period rule | what stays current (vintage leakage) |
|---|---|---|
| Trials (research_readiness condition / technology trials; desert trial term) | `study_first_post_date <= T-12-31` (the registration was public at T). Condition matching = the literal view already in `trial_conditions` / `facility_trials`. | Site lists are today's records: a site added to a pre-T trial after T is counted as pre-T (not recoverable without registry history). |
| NIH (research_readiness NIH term) | `fiscal_year <= T` (rows of `facility_nih_projects`), precision view | Data start FY2015 (no earlier awards). Organisation geocodes are current. |
| CMS MMD burden (ME/CFS, level C) | CMS condition 51, unsmoothed actual, all ages, **claims year T** (matched to 2024 counties by the same ingestion) | none |
| Long COVID burden (HPS) | **Does not exist before mid-2022** (T = 2021 and earlier). Handled by two variants (section 3). | variant B only |
| Providers (clinic_capacity; desert provider term; reach) | individual NPIs with `enumeration_date <= T-12-31` | Current NPPES snapshot: taxonomy and practice address are today's; NPIs deactivated before the snapshot are absent (survivorship). |
| Facilities (clinic_capacity; reach) | a clinic-candidate facility exists at T if ANY of: an NPI in `npi_object_ids` enumerated <= T; an HRSA site with `site_added_to_scope_date <= T`; a trial at the facility first posted <= T; an NIH project at the facility with fiscal year <= T. Facilities with none of these dated sources are dropped (counted). | HRSA roster is active sites only (closed sites absent); facility attributes (specialty groups, FQHC flag) are current. |
| SVI (vulnerability, desert vulnerability term) | static: SVI 2022 | whole input (slow-moving; reported as leakage) |
| ACS population, adults, population within 50 km | static | whole input (denominators) |
| Measurement performance (evidence_weighted only) | static: the current `measurement_performance` record | the record comes from post-T data; it is constant within a condition x measurement, so it only rescales the expected-yield component, never reorders regions by itself |
| Curated configs (relevance.yaml, scoring.yaml) | static | written in 2026 |

**Leakage bound.** (i) Counts of what the date filters remove (share of current relevant NPIs enumerated after T,
share of facilities without a pre-T date, trials and NIH rows removed). (ii) A **static-context sensitivity**: the
holdout re-run with providers and facilities NOT date-filtered (today's snapshot). If the headline metrics move by
less than their bootstrap CI half-width, leakage through the provider / facility path is small relative to sampling
uncertainty; if they move more, the pre-period result is reported as sensitive to it. SVI / ACS leakage is not
bounded by a re-run (no earlier vintage is ingested); it is stated.

## 3. Rankings rebuilt with pre-period inputs

Scored with the unchanged engine, county and state level, universe = 50 states + DC (2024 counties), ranking-eligible
= population >= 10,000 (ACS) and complete pre-period burden (engine rule).

* **Q1 (primary): Long COVID or ME/CFS x wearable autonomic/activity monitoring, variant A ("clean").** Long COVID has
  no burden measure before T, so in the view its burden is undefined (level D) at every level: the set burden is the
  ME/CFS CMS proxy of year T alone (the engine's rule for a member without a defined burden); research readiness
  counts Long COVID and ME/CFS trials / NIH projects dated <= T; desert = the ME/CFS member's burden-informed desert.
  This is a clean holdout apart from the static context of section 2.
* **Q1-B (NOT a clean holdout):** as Q1 but Long COVID burden held at its current (post-period) HPS values
  (state estimate inherited by counties). Reported only to show how much the missing Long COVID burden matters,
  labelled "not a clean holdout" wherever it appears.
* **Q2: ME/CFS x wearable** (clean; CMS proxy of year T), with the ME/CFS-only outcome as its secondary outcome.
* Weight sets: **equal** (primary), **evidence_weighted** (pre-declared primary alternative; defined because the pooled
  and ME/CFS performance records are known), and every other named set of `configs/scoring.yaml` (reported,
  secondary). Single components: burden_pct, vulnerability_pct, diagnostic_desert_pct, clinic_capacity_pct,
  research_readiness_pct.
* **Baselines** (same eligible counties):
  * burden-only = the engine's `burden_only` rank (pre-period burden);
  * population-only = `population_within_50km` (primary: the outcome is defined within 50 km) and
    `population_total` (secondary);
  * random = uniform random ranking (AUROC 0.5 and precision = base rate in expectation; the empirical 95% range
    over 1,000 random rankings is reported);
  * prior activity (added here because persistence of research activity is the obvious explanation of any hit):
    distinct relevant trials first posted <= T with a U.S. site in the county or within 50 km.

## 4. Post-period outcomes (after T to the retrieval date)

Relevant = the literal-precision condition view (`trial_conditions.condition_literal_match`), conditions Long COVID
or ME/CFS for Q1 (ME/CFS only for the Q2 secondary outcome). NIH relevance = the precision view
(title/abstract match, likely false positives excluded) used by `geography.features.nih_counts`.

* **O1 (PRIMARY): any new relevant trial within reach** - a county is positive if >= 1 relevant trial first posted
  after T-12-31 has a U.S. site in the county or within 50 km of the county's Census internal point (the same
  "in county OR within 50 km" pool and the same counter, `geography.features.trial_counts`, as the desert's trial
  term). Binary.
* O1-count: number of such distinct trials (for the Spearman correlation).
* O1-in-county (sensitivity: stricter locality): a site inside the county.
* O1-start (sensitivity: date rule): trials partitioned by `start_date` instead of first-posted date.
* O2: ME/CFS-only version of O1 (secondary outcome for Q2 and for Q1).
* O3: new NIH award - a relevant NIH core project whose first fiscal year in the data is > T, at an organisation
  geocoded to the county (binary and count; `nih_counts` rule, in-county because the unit is the awardee's address).
* O4: new facility - a facility (from `facility_trials`) with >= 1 relevant post-T trial and no relevant trial first
  posted <= T, located in the county or within 50 km (binary and count). This is expansion to sites without prior
  relevant experience, i.e. the part of the activity that persistence alone cannot explain.
* O5 (burden input persistence): CMS condition 51 prevalence in 2022 (latest final year) vs year T; and the change
  2022 minus T.
* State level: O1-count, O3-count, O4-count per state and per 100k (Spearman); binary O1 reported but flagged as
  ceiling-limited if >= 80% of states are positive.

## 5. Metrics, inference and pre-specified readings

### 5.1 Metrics (county level; evaluation set = ranking-eligible counties of the pre-period ranking being evaluated)

For every score (composite under each weight set, each component, each baseline) against O1 (and, secondary, O2, O3,
O4, O1 sensitivities):
* AUROC of the score for the binary outcome (higher score = predicted positive; ties averaged);
* top-decile hit rate (share positive among the top 10% of eligible counties by rank) and lift = hit rate / base rate;
* precision@k for k = 10, 25, 100 (engine ranks; baselines ranked by score, ties by geo_id);
* Spearman correlation with the post-period count.

### 5.2 Inference
* **County bootstrap** (1,000 resamples of eligible counties with replacement, seed = config.SEED; rankings are fixed,
  top-k and top-decile re-selected within each resample) for a 95% percentile CI of every metric; **paired** bootstrap
  for the difference of AUROC (score minus each baseline).
* **State-cluster bootstrap** (resample states) for AUROC and Delta AUROC of the primary score: sensitivity for spatial
  dependence (overlapping 50 km pools make counties non-independent).
* **Permutation null** (1,000 permutations of the outcome across eligible counties; and within state): p-value for
  AUROC and precision@25 = (1 + #null >= observed) / 1,001.
* **Population control:** (a) population-only baseline and paired Delta AUROC; (b) population-stratified AUROC:
  AUROC within deciles of `population_within_50km`, pooled with weights n_pos x n_neg of the decile; (c) logistic
  regression O1 ~ log(1 + population_within_50km) + log(1 + population_total) + z(score): odds ratio per SD of the score
  with a county-bootstrap 95% CI and the likelihood-ratio test; (d) the same model plus log(1 + prior relevant trials
  within 50 km) (incremental over prior activity); (e) partial Spearman of the score with O1-count given
  log(1 + population_within_50km).

### 5.3 Expected directions (pre-specified) and readings

| score | expected relation to O1 / O3 / O4 | reading |
|---|---|---|
| research_readiness_pct | positive (persistence of research activity) | readiness supported if AUROC 95% CI > 0.5 AND it beats population-only (Delta AUROC CI > 0) |
| clinic_capacity_pct | positive (provider / facility density) | as above; expected to behave like population |
| burden_pct (ME/CFS CMS proxy) | no directional prediction | reported |
| vulnerability_pct | negative (research goes to less vulnerable places) | reported; negative is the expected, not a failure |
| diagnostic_desert_pct | **negative** (deserts have little activity by construction) | **gap persistence** if the desert top-decile hit rate is below the base rate (lift 95% CI upper bound < 1). An AUROC > 0.5 would mean the desert label points at places that then attracted research, i.e. the desert did not describe a lasting gap. |
| composite equal | mixes positive (readiness, capacity) and negative (desert, vulnerability) parts: expected weakly positive | **research-readiness validation of the composite passes** if all three hold: (a) AUROC 95% CI lower bound > 0.5; (b) permutation p < 0.05; (c) Delta AUROC vs population-only with 95% CI > 0. Passing (a)-(b) but not (c) = "explained by population". |
| composite evidence_weighted | as equal | same rule |
| study_partner / capacity_led | most positive of the named sets (readiness-weighted) | reported |
| access_gap | around or below 0.5 (no readiness weight; desert-led) | reported as a desert ranking (gap persistence rule) |

Adds over prior activity: the composite's odds ratio in model 5.2(d) with 95% CI excluding 1, and Delta AUROC vs the
prior-activity baseline. A composite that does not beat prior activity is reported as "no better than where activity
already was".

Every reading is reported whether positive or null; null and negative results sit in the same table as positive ones.

## 6. State level

Q1 at state level (51 states): Spearman of each score with O1-count, O3-count and O4-count, raw and per 100k
population, with state-bootstrap CIs, and the population-only baseline. Binary O1 AUROC only if fewer than 80% of
states are positive.

## 7. Stability

* **Leave-one-state-out (LOSO) normalisation.** For each state s (51), the county ranking of the primary query is
  recomputed with s **dropped from the normalisation reference**: every percentile rank the engine and the desert
  compute is taken against the counties of the other 50 states only (a county of s gets the percentile it would have
  if it were added to that reference: average rank among reference plus itself, divided by reference size plus one;
  reference counties keep the engine's definition within the reference). Reported per state: Spearman of s's
  counties' composites with and without s in the reference, and for the counties of the other states the Spearman of
  their composites and the top-25 overlap. Summary: median and minimum over states. Reading: median within-state
  Spearman >= 0.9 = stable; states with Spearman < 0.8 listed. Run for the current ranking (stored inputs) and for the
  T = 2021 ranking.
* **Rolling check.** Spearman / Kendall tau-b and top-10 / top-25 overlap between the T = 2019, T = 2020, T = 2021
  and current equal-weight county rankings (over counties eligible in both).

## 8. Guardrails

Ecological only: places, facilities and registries, never people. A county "hit" means a registered trial site or
NIH awardee was within reach; it does not describe any patient. Trial registration is a research-activity signal, not
evidence a measurement works. Every table carries provenance; `temporal_holdout_rankings` keeps every component of
every pre-period ranking. Language: "research-readiness validation", "proxy outcome", "candidate deployment
opportunity"; never a claim of clinical value.

## 9. What is not done

* No claim about clinical value, diagnostic yield, patient benefit or the "right" deployment places.
* No reconstruction of historical NPPES / HRSA / SVI vintages (not ingested); leakage is bounded, not removed.
* No Monte Carlo rank intervals for the pre-period rankings (the holdout evaluates the point ranking).
* No new data retrieval: everything comes from the processed tables already built (the step runs offline).

## Deviations (added after the analysis ran)

All dated 2026-09-25, after the first full run; none changes an outcome definition, a cutoff, a score or a reading
rule.

1. **Population-stratified AUROC CI:** 200 county-bootstrap resamples instead of 1,000 (it re-bins population deciles
   in every resample; computational). Every other CI uses 1,000.
2. **Random baseline:** as pre-specified it is summarised by 1,000 random rankings
   (`temporal_holdout_random_baseline.csv`). The per-score table also carries one random ranking scored like every
   other score (`baseline_random`); that row is a single draw and its CI / permutation p describe that draw only
   (in the first run it happened to fall at the upper edge of the random range).
3. **Burden-only baseline = burden component:** in the clean variant the set burden has one defined member (ME/CFS),
   so the engine's burden-only composite equals `burden_pct`; the two rows are identical by construction.
4. **O1-start (clarification):** only the outcome is partitioned by start date; the pre-period inputs stay
   partitioned by first-posted date.
5. **Trial dates:** `start_date` is month precision for most trials ("2003-12"); dates are parsed with mixed formats
   (first of the month). The first smoke run parsed only full dates; fixed before any reported number.
