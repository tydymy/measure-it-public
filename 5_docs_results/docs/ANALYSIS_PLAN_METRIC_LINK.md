# Analysis plan: the metric -> translation link (measurement performance, expected yield, evidence-weighted scoring)

Status: **pre-specified 2026-09-25, before any new number was computed.** At the time of writing, the inputs had been
read (the existing result tables and reports of the modules listed in section 2, `published_device_evidence`, the
schemas of `geo_condition_features`, `geo_context_zcta__acs`, `zcta_centroids` and `deployment_opportunities`), and
nothing below had been computed: no sensitivity at a fixed specificity, no reach, no expected yield, no
evidence-weighted composite or rank. Deviations made later are listed at the end, dated, with the reason; nothing above
that section is edited after results exist.

Owner module: `measure_it.scoring.metric_link` (performance table, reach, expected yield, joint ranking), wired into
`measure_it.scoring.opportunity` (components, composite, Monte Carlo), `measure_it.scoring.controls` (Test 5 contrast),
`measure_it.scoring.sensitivity` / `report` (Test 7, S5) and `measure_it.scoring.recommend` (output). Pipeline step
`scoring_metric_link` (before `scoring_report`). Method reference: `docs/SCORING.md` section 9. Results:
`results/SCORING_RESULTS.md` section 12.

## 0. The problem this answers

The deployment ranking did not depend on the measurement's own evidence. Phenotype and measurement evidence (e.g.
MUSCLE-ME hip steps AUROC 0.831, LC + ME/CFS vs healthy; Uwakwe whole-record resting HR for Long COVID AUROC 0.504)
was displayed next to recommendations but never entered the ranking: a measurement entered only through its implementer
groups (clinic capacity) and its trial mentions (research readiness), so swapping the wearable bundle for nailfold
capillaroscopy, for which no public performance evidence exists, gave almost the same top-10 (S5: Spearman 0.994,
9/10 overlap). This plan links the metric (how well the measurement separates cases from non-cases, in this project's
own computations first) to the translation decision (where, and whether, to deploy it).

Two facts fixed in advance, because they follow from arithmetic, not data:

1. **Within one condition x measurement, a per-measurement constant cannot reorder regions.** The measurement's
   evidence score is the same for every region, so it shifts every composite by the same amount. It can change which
   *measurement* is preferred (a ranking over region x measurement pairs), and whether a measurement is ranked at all
   (UNKNOWN rule), but not the order of regions for a fixed measurement.
2. **A percentile-normalised yield would lose the sensitivity.** pr(burden x sensitivity x reach) = pr(burden x reach)
   when sensitivity is constant within a combination. The yield component is therefore on an absolute scale
   (section 4.3), so that a less sensitive measurement contributes less, and region x measurement pairs are comparable.

## 1. Questions

* **Q1 (table).** For each (condition or condition set, measurement class or bundle) of the scored grid, what is the
  best available performance record, drawn only from this project's own computations and, separately flagged, from
  `published_device_evidence`? Where nothing exists the value is UNKNOWN.
* **Q2 (yield).** For each (geography, condition, measurement), what are the expected detectable cases (planning
  estimate) and the expected false positives at the pre-specified operating point, with every factor as a column and
  uncertainty from the performance CI and the burden Monte Carlo?
* **Q3 (scoring).** With `measurement_evidence` and `expected_yield` as components and a new pre-declared weight set
  `evidence_weighted`, how do the rankings change (T5-4: evidence_weighted vs equal), and does the adapter swap (S5)
  now differ because of evidence (joint ranking over region x measurement pairs), with UNKNOWN-performance
  measurements flagged instead of ranked as if equivalent?

## 2. Performance records (inputs; the list is closed)

A **record** is one discrimination result: a target condition (case group), a comparator, a measurement class (the
class of the feature or device), a metric, a CI, n, a label basis and source object ids. Records are never moved
between conditions or measurements: no value is imputed from a neighbouring condition, measurement or dataset.

### 2.1 Own-computed records (this project's result tables)

| id | target | class | analysis | source (read, not recomputed, unless marked NEW) |
|---|---|---|---|---|
| OWN-MM-STEPS-POOLED | long_covid_or_me_cfs | accelerometry | MUSCLE-ME daily steps, LC + ME/CFS vs healthy, complete case, fixed direction (lower = case) | `charlton_lc_mecfs_cpet_source_primary.csv` (analysis=primary); worst/best-case bounds from the same file as caveats |
| OWN-MM-STEPS-LC | long_covid | accelerometry | MUSCLE-ME steps, LC vs healthy | `..._secondary.csv` (a) |
| OWN-MM-STEPS-ME | me_cfs | accelerometry | MUSCLE-ME steps, ME/CFS vs healthy | `..._secondary.csv` (a) |
| OWN-MM-VO2-POOLED | long_covid_or_me_cfs | cpet | MUSCLE-ME relative VO2peak, LC + ME/CFS vs healthy | `..._secondary.csv` (b) |
| OWN-MM-VO2-LC | long_covid | cpet | MUSCLE-ME VO2_rel, LC vs healthy | **NEW**: same function (`muscle_me_steps_analysis.fixed_direction_auroc`, lower = case, 2,000 stratified bootstraps, seed 20260923) on the same person-level table |
| OWN-MM-VO2-ME | me_cfs | cpet | MUSCLE-ME VO2_rel, ME/CFS vs healthy | **NEW**, as above |
| OWN-UW-RHR | long_covid | wearable_heart_rate | Uwakwe 2025 whole-record resting-HR features, L2 logistic, 10x5 CV | `stanford_uwakwe_model_performance.csv` (wearable, l2_logistic) |
| OWN-NH-MECFS-PROXY | me_cfs | accelerometry | NHANES wearable-only pen_lr, ME/CFS-like PROXY vs rest of the wearable-valid sample | `nhanes_model_performance.csv` (mecfs_like_proxy, cv_5x5, wearable, pen_lr) |
| OWN-APP-VO2 | long_covid | cpet | Appelman 2024 VO2max-only CV model, LC vs healthy | `labs_dx_appelman_models.csv` (vo2max_only) |
| OWN-APP-METAB | long_covid | metabolomics | Appelman plasma metabolomics top-10 CV model (locked primary; null by its rule) | `labs_dx_appelman_models.csv` (blood_top10) |
| OWN-KL-CORT | long_covid | blood_biomarkers | MY-LC fit-free serum cortisol, LC vs HC + CC (primary sample) | `klein2023_mylc_ml_table_cortisol_fitfree.csv` (S6) |
| OWN-KL-CORT-CC | long_covid | blood_biomarkers | MY-LC fit-free cortisol, LC vs convalescent controls | same file (S2 LC vs CC) |

Not attached (listed so the omission is visible): Stanford acute-infection window detection (target is an acute
infection episode, not a condition in the grid); NHANES models of general-population labels (fatigue, functional
limitation, depression, fair/poor health, mortality: not grid conditions); fibromyalgia thermography, endometriosis
arginase and hEDS/HSD Olink (conditions outside the scored grid).

### 2.2 Published records (`published_device_evidence`; tier 3, never reproduced here)

Only `claim_type = diagnostic_accuracy` rows whose case group is a grid condition. The case group and the primary
device class are read from `condition_definition_in_study` / `device_or_test` (curated here, listed in full):

| row | target (case group) | primary class | metric available | attaches to bundle |
|---|---|---|---|---|
| PDE-006 | pots | ecg_ambulatory (chest ECG + accelerometer) | AUROC 0.900; se 0.95 / sp 0.85 | wearable |
| PDE-007 | pots (POTS within Long COVID) | ecg_ambulatory (24-h ECG rule) | se 0.9111 / sp 0.9636; sp is 53/55 vs Long COVID without POTS (the look-alike arm), so its CI uses n = 55 | wearable |
| PDE-008 | pots | ecg_ambulatory (Holter HRV) | AUROC 0.874; se 0.778 / sp 0.886 | wearable |
| PDE-002 | pots | autonomic_testing (10-min active stand) | se 0.87 / sp 0.67 | autonomic_function_testing |
| PDE-003 | pots | autonomic_testing (10-min tilt) | se 0.93 / sp 0.40 | autonomic_function_testing |
| PDE-022 | me_cfs | autonomic_testing (tilt + Doppler CBF) | sensitivity 0.90 only; specificity not reported | autonomic_function_testing |
| PDE-014 | me_cfs | cpet (2-day CPET classification) | 95.1% in-sample accuracy (not AUROC or sensitivity) | exercise_capacity_testing |
| PDE-040 | long_covid | wearable_heart_rate | AUROC 0.951 **with symptom features**: not a measurement-only metric; excluded from scoring (and not reproducible from the public release) | wearable (listed, excluded) |
| PDE-043 | long_covid | hrv | F1 56% (null); no AUROC or sensitivity | wearable and autonomic_function_testing (hrv is a member of both) |

PDE-049 (pupillometry) has no measurement class and attaches to nothing. Rows of conditions outside the grid are out
of scope.

### 2.3 Record fields and quality tier

`auroc`, `auroc_ci_low`, `auroc_ci_high`, `auroc_basis`, `op_sensitivity` (+ CI), `op_specificity` (+ CI),
`op_basis`, `n_cases`, `n_controls`, `label_basis`, `comparator`, `comparator_kind` (healthy / recovered after
infection / general population / look-alike or mixed), `quality_tier`, `source_object_ids`, `caveats`.

**Quality tier** (fixed): 1 = own-computed, clinical or study case definition applied by the investigators (MUSCLE-ME,
MY-LC, Appelman); 2 = own-computed, proxy or self-reported label (Uwakwe self-reported Long COVID; NHANES ME/CFS-like
proxy); 3 = published claim (not reproduced); 4 = none (UNKNOWN).

**AUROC lower bound.** Own records: the producer's 95% CI. Published: the published CI if one exists; otherwise the
Hanley-McNeil 95% CI computed here from the published AUROC and counts; a row with only a sensitivity/specificity pair
gets the single-operating-point bound AUROC >= (se + sp) / 2 (true for any concave ROC through that point), with lower
bound (se_lo + sp_lo) / 2 from Wilson 95% intervals of the published counts; `auroc_basis` says which. A sensitivity
without a specificity gives no AUROC (UNKNOWN).

### 2.4 Operating point (pre-specified specificity 0.90)

* **Own fixed-direction scores with person-level values** (MUSCLE-ME steps and VO2_rel; Uwakwe stored repeat-averaged
  out-of-fold probabilities): threshold c = the ceil(0.90 x n_controls)-th smallest control score (score oriented so
  that higher = case); a person is flagged if score > c. Sensitivity = share of cases flagged; achieved specificity =
  share of controls not flagged (>= 0.90 by construction; reported). 95% CI: 2,000 stratified bootstrap resamples
  (seed 20260923), c re-chosen in every resample. The threshold is chosen on the same controls (in-sample), which
  flatters it; the bootstrap re-choosing c only partly corrects this.
* **Own records without stored person-level scores** (NHANES proxy pen_lr, Appelman CV models): operating point
  UNKNOWN (not refitted here).
* **MY-LC cortisol**: the producer's stored sensitivity at 90% specificity; Wilson 95% CI from that proportion and
  n_cases (threshold-selection uncertainty not included; stated).
* **Published rows**: usable only when the reported specificity is >= 0.90: then (reported se, reported sp) is the
  operating point (sensitivity at exactly 0.90 specificity can only be equal or higher, so this is conservative), with
  Wilson CIs from the published counts. Reported specificity < 0.90 -> operating point UNKNOWN (never extrapolated).

### 2.5 Attachment and selection

* A record attaches to (condition c, class k) when its target is c (a condition set: only a record whose case group is
  exactly the set's members pooled) and its primary class is k; to (c, bundle b) when its primary class is a member of
  b. No record for a set is built from member records.
* `performance_status`: `known` = AUROC lower bound AND operating point; `partial` = AUROC lower bound only (no
  operating point at >= 0.90 specificity); `UNKNOWN` = no record, or no usable metric.
* **Scored record per (c, m)** (deterministic): among attached records not excluded as non-measurement-only, (1) best
  quality tier; (2) status known > partial > UNKNOWN; (3) the most conservative (lowest `measurement_evidence`);
  (4) larger n_cases; (5) record id. All attached records are kept in the row (`records` JSON) and in
  `results/tables/measurement_performance_records.csv`.
* `evidence_conflict = True` when another own-computed record of the same (c, m) has an AUROC CI that includes 0.5
  while the scored record's excludes it (e.g. a null resting-HR result next to a positive steps result); stated in the
  row and in the recommendation's uncertainties.

### 2.6 measurement_evidence (the evidence component)

`measurement_evidence = tier_factor x clip((AUROC_lower - 0.5) / 0.5, 0, 1)`: discrimination above chance at the lower
95% bound, on an absolute 0-1 scale (0 = no demonstrated discrimination; 1 = perfect), not percentile-normalised
(it is constant across regions). `tier_factor`: tier 1 = 1.0, tier 2 = 0.75, tier 3 = 0.5 (a published claim that
this project has not reproduced counts half). The undiscounted value is kept (`measurement_evidence_undiscounted`).
UNKNOWN when no AUROC lower bound exists.

Table `measurement_performance` (processed, provenance columns, `object_id =
measurement_performance:<condition>|<measurement>`): the grid of conditions long_covid, me_cfs, pots, dysautonomia,
long_covid_or_me_cfs, autonomic_activity_invisible_illness x measurements = the four scored bundles, their member
classes, and the two extra classes with own records (blood_biomarkers, metabolomics). One row per cell; UNKNOWN cells
are rows too.

## 3. Implementation reach (per measurement, per region)

`reach(g, m)` = share of the population of region g living in ZCTAs whose internal point is within the default radius
(50 km, `scoring.yaml clinic_matching.default_radius_km`) of at least one geocoded clinic-candidate facility (the
`find_candidate_clinics` / clinic_capacity base) with >= 1 of m's implementer groups (configs/relevance.yaml, curated).
ZCTA population: ACS (`geo_context_zcta__acs.total_population`); ZCTA location: `zcta_centroids`; ZCTA -> county:
`county_fips_2024_of_internal_point` (state: `state_fips`). Total population is used as the stand-in for the adult
share. A region with no populated ZCTA has reach UNKNOWN. Reach is condition-independent and is an access assumption
(an implementer-group facility exists nearby), not evidence that any facility runs the measurement.

## 4. Expected yield and false positives

### 4.1 Burden basis (per condition x level)

* **Count path**: only where a defensible count exists at that resolution: a single condition whose burden measure is
  level A or B, not inherited, in percent of all adults. burden_count = value / 100 x ACS adults 18+. In this build that
  is long_covid at state level (HPS, currently experiencing long COVID, % of adults).
* **Prevalence-percentile path**: everywhere else (long_covid county = inherited state value, not a county count;
  me_cfs = level-C claims proxy among Medicare FFS, not an ME/CFS count; condition sets). The burden percentile
  (`burden_pct`) replaces the count and no case count is produced; every row says so (`expected_yield_basis`).
* **Level D** (pots, dysautonomia; no burden measure): no yield (see the UNKNOWN rule, section 5).

### 4.2 Columns (every factor kept)

`reach`, `reach_population`, `reached_adults = adults18 x reach`, `op_sensitivity`, `op_specificity`,
`measurement_performance_status`, `measurement_performance_tier`, `expected_yield_basis`, and:

* count path: `expected_detectable_cases = burden_count x op_sensitivity x reach`;
  `expected_false_positives = (1 - op_specificity) x (adults18 - burden_count) x reach`;
  `expected_ppv = TP / (TP + FP)`; `false_positives_per_detected_case = FP / TP`;
* percentile path: `expected_yield_index = burden_pct x op_sensitivity x reach` (dimensionless, **not a count**);
  `expected_false_positives_upper = (1 - op_specificity) x reached_adults` (upper bound: every reached adult treated as
  a non-case; no prevalence is assumed). No PPV.

Language fixed in every output: expected yields are **planning estimates for a candidate pilot, not predictions of
diagnoses**; performance measured against healthy controls **overstates** real-world performance (the clinical
question is the illness versus other causes of the same symptoms); where prevalence is low most flags are false
positives.

### 4.3 The expected_yield component (normalised)

`expected_yield = pr(B) x op_sensitivity x reach`, where pr(B) = `burden_pct` (percentile path) or the percentile
rank of burden_count over the level's universe (count path). It is in [0, 1], varies with burden and reach across
regions, and scales with the measurement's sensitivity, so it is comparable across measurements for the same
condition and level. Not re-normalised with pr() (that would remove the sensitivity; section 0).

### 4.4 Uncertainty

Monte Carlo, 1,000 draws (seed 20260923 + 11) per combination: burden per member as in the existing Monte Carlo
(`opportunity.draw_burden`: CI x evidence multiplier, inherited-state error x the inherited multiplier + tau);
op_sensitivity and op_specificity each drawn on the logit scale, normal with SD = (logit(hi) - logit(lo)) / 3.92 from
their 95% CIs (a bound at 0 or 1 is clipped to [0.001, 0.999] first); weights ~ Dirichlet(20 x evidence_weighted
weights). Per draw: burden_pct, desert, yield, composite, rank. Outputs: 5th/50th/95th percentiles of
expected_detectable_cases (count path) or expected_yield_index, of expected false positives, and of the
evidence-weighted rank (`rank_mc_ew_p05/p50/p95`, `p_top10_mc_ew`, `p_top25_mc_ew`). measurement_evidence is not drawn
(it is already a lower 95% bound). Reach is not drawn (no sampling error; an access assumption).

## 5. Scoring

* New components `measurement_evidence` (section 2.6) and `expected_yield` (section 4.3), both absolute [0, 1], added
  to the component dictionary; the five original components, their normalisation and every existing weight set and
  composite are unchanged (continuity: `equal` stays the default and primary for Tests 5-7 as before).
* New weight set, pre-declared as the **primary alternative**: `evidence_weighted` = burden 0.15, vulnerability 0.15,
  diagnostic_desert 0.15, clinic_capacity 0.10, research_readiness 0.10, measurement_evidence 0.15, expected_yield 0.20
  (the measurement's own evidence carries 0.35; the rest keep their equal-weight proportions roughly). Added to
  `configs/scoring.yaml weight_sets` so every named-set analysis includes it.
* **UNKNOWN rule (fixed).**
  1. `performance_status` UNKNOWN or partial -> the combination gets **no** evidence_weighted composite or rank in any
     region: `evidence_weighted_status = "not ranked: measurement performance UNKNOWN"` (or `"... partial: no
     operating point at >= 90% specificity, so expected yield is UNKNOWN"`). It is never scored 0 (which would read as
     "shown not to work") nor ranked on the remaining components (which would read as "equivalent"). Its equal-weight
     rank is kept for continuity, and every output that shows it carries `measurement_performance_status` and the
     statement that the equal-weight rank contains no measurement-performance information.
  2. Performance known but burden level D (no burden measure): expected_yield is UNKNOWN because of the burden, not the
     measurement; the yield weight is removed and the other weights renormalised, as for level-D burden, and flagged
     (`expected_yield_excluded_level_D`).
  3. A region with known performance and complete burden but reach UNKNOWN: no evidence_weighted composite in that
     region (listed).
  4. Incomplete-burden regions stay unranked (deviation 7 of the scoring plan).
  5. Every recommendation's uncertainties name the performance tier (1-3, or UNKNOWN), the comparator, the label basis
     and the healthy-control caveat.
* **Joint ranking** (`rank_evidence_weighted_joint`, `rank_equal_joint`): per condition x level, all region x
  measurement pairs of the four scored bundles with a composite (ranking-eligible regions), ranked by the composite; for
  `evidence_weighted` the UNKNOWN-performance bundles have no pairs.

## 6. Analyses and pre-specified readings

* **T5-4** (Test 5 table): `evidence_weighted` vs `equal` for every combination where evidence_weighted is defined;
  top-10/25 overlap, Jaccard, Kendall tau-b, Spearman, entering/leaving regions with component contributions. Reading
  (descriptive): within one combination the change comes from expected_yield and the weight change only (section 0,
  fact 1).
* **Test 7**: `evidence_weighted` is one of the named weight sets (Kendall's W, pairwise tau, instability rule).
* **S5 with evidence** (primary condition set, county): (a) the within-bundle comparisons (wearable vs each other bundle)
  under `equal` (as before) and under `evidence_weighted` where both are defined, and "not ranked: performance UNKNOWN"
  where not; (b) the joint ranking: bundle composition of the joint top-10 and top-25 under `equal` and under
  `evidence_weighted`. Reading, fixed: the ranking "differs because of evidence" if the joint top-25 bundle composition
  differs between equal and evidence_weighted, or if a bundle ranked under equal is unranked under evidence_weighted
  for lack of performance evidence. Both are reported with the numbers, whichever way they fall.
* The new top-10 of the primary query under evidence_weighted and under equal, side by side, with the performance
  record used.

## 7. What is not done

No model is refitted; no performance value is imputed or carried between conditions, sets, classes or datasets; no
county prevalence is created from a state value (the county count path is not used for inherited burden); no expected
yield is presented as a forecast of diagnoses. Weights, tier factors, the 0.90 specificity, the radius and the seeds are
fixed above.

## Deviations (added after the analysis ran)

1. **2026-09-25, labelling fix after the first run (no score changed).** The first build derived the published rows'
   `comparator_kind` with a string rule, which labelled PDE-040 (infected participants without Long COVID) as
   "look-alike". It now uses the curated labels of the section 2.2 table (PDE-007: look-alike (Long COVID without POTS)
   + healthy; PDE-040: recovered after infection; the others healthy). `comparator_kind` enters no score.
2. **Reported, not a deviation:** because `evidence_weighted` is a named weight set (section 6), Test 7's
   Kendall's W across the named sets and Figure 8 panel (c) now include it; every equal-weight composite, rank and
   Monte Carlo interval is bit-identical to the build before this change (checked row by row over all 76,680 rows).
3. **2026-09-28, extension for user-supplied data (no number of this plan changed).** `measure_it.byod` adds a record
   source: user-supplied datasets evaluated by this engine under a plan locked before any outcome
   (docs/BRING_YOUR_OWN_DATA.md). Their records get the tier `user_supplied_own_computation` (numeric 2.5, evidence
   factor 0.6), ranked after tiers 1-2 and before tier 3 in the section 2.5 selection, and count as own-computed for
   `evidence_conflict`. A (condition, class) of a user record outside the grid of section 2.6 gets its own
   `measurement_performance` row (it is not scored). Records of `demo` datasets enter only with
   `MEASURE_IT_BYOD_DEMO=1`. With no user dataset evaluated, the records, the performance table, every composite and
   rank are identical to the build before this change (checked table by table after `measure-it byod remove`).
