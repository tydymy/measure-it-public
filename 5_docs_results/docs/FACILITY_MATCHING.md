# Facility, clinic and research-site registries; clinic matching; Test 6 (clinic locations)

Status 2026-09-23. Modules: `measure_it.facilities.registry` (build), `.matching` (`find_candidate_clinics`),
`.query` (`find_relevant_research_centers`), `.shuffle_test` (Test 6). Analysis plan written before the analysis:
`docs/ANALYSIS_PLAN_FACILITIES.md` (its §9 lists the deviations). Every number here is read from a file written by
this code; the file is named in brackets.

**Independent review (2026-09-23).** The pipeline was re-run end to end. Every table and results file reproduced
byte-for-byte apart from runtime stamps. The review changed the reading of Test 6, not its data:
* Multiplicity (BH-FDR) and county-sampling uncertainty (county bootstrap) were added to the Test 6 summary.
* Three contrasts that the first write-up reported as positive are not supported: metro−nonmetro technology
  experience, share with a Long COVID trial, and km to the nearest trial site.
* Criterion (a) (low Jaccard) is noted as true by construction.
* It is now stated that 1,478 of the 8,722 research sites have no precision-view record (new
  `has_precision_view_record` column), and that V2 is a lower bound on false links.

Reproduce, in order:

```
uv run python -m measure_it.facilities.registry        # ~4 min: facilities, registries, V1/V2 validation
uv run python -m measure_it.facilities.shuffle_test    # ~5 min, 6 workers: Test 6 (200 national + 100 within-state draws)
uv run python -m measure_it.facilities.query           # example outputs + draft overlay map
uv run pytest tests/test_facilities.py
```

## 1. What was built

| table (data/processed) | rows | unit |
|---|---|---|
| `facilities` | 329,019 | resolved facility (object id `facility:<facility_id>`) |
| `facility_source_links` | 398,907 | source record → facility, with `match_method`, `match_confidence`, linked record |
| `clinic_registry` | 321,536 | facilities with an NPPES organisation NPI or an HRSA site (implementation-oriented) |
| `research_site_registry` | 8,722 | facilities with ≥1 ingested U.S. trial site or NIH project (research-oriented); 7,244 have precision-view evidence (`has_precision_view_record`) |
| `facility_trials` | 17,598 | facility × NCT (literal and MeSH-expanded condition ids, status, registered measurement classes) |
| `facility_nih_projects` | 9,297 | facility × appl_id (precision-view and all-tier condition ids) |

[facility_resolution_summary.csv; .meta.json sidecars]

Source records: 369,339 NPPES organisation NPIs, 19,266 HRSA sites, 9,761 U.S. ClinicalTrials.gov facility keys,
541 NIH RePORTER organisations. Of 27,828 U.S. trial-site rows, 8,129 flagged upstream as placeholders
(`trial_sites.facility_generic`: sponsor placeholders such as "Research Site", and since the 2026-09-23 integration
also non-facility locations such as virtual/remote sites, "Contact <sponsor> for exact location", "Site Reference
ID/Investigator#", and names that are only the site's city), 1,995 unnamed sites and 19 further names that are only a
place were excluded; 17,685 rows were kept. (Before the upstream flag was extended, the same 8,148 excluded named
rows were counted as 7,597 placeholders, 137 non-facility locations and 414 place names; the kept rows are identical.)
[facility_resolution_summary.csv]

Geocode precision of the 329,019 facilities: `zcta_centroid` 308,258; `point` 19,226 (HRSA address geocodes and
RePORTER points); `city_centroid` 308 (ClinicalTrials.gov geoPoints where no usable ZIP); `none` 1,227.
ClinicalTrials.gov sites are placed at their ZIP's ZCTA centroid when the ZIP resolves in the registered state, else at
the geoPoint (a city centroid). [facility_resolution_summary.csv]

Kinds: 200,299 facilities are NPPES **group practices** (organisation NPIs under physician taxonomies only), 101,494
NPPES facility-type organisations, 19,266 HRSA sites, 7,180 trial-only sites, 368 NIH-only organisations and 412
research-taxonomy organisations (primary kind). [facility_resolution_summary.csv]

## 2. Entity resolution and how well it works

Links need a normalised-name agreement **and** the same ZIP5, or the same city+state within 25 km with at least two
informative name tokens. Never name only, never across states. Rules and confidences are in the plan §2. Clustering
allows at most one HRSA site and one NIH organisation per facility, no facility wider than 25 km, and no mixed
states. Edges rejected by these constraints: two HRSA sites 7,878; span > 25 km 250; cross-state 32 (NPPES
address-state typos with a correct ZIP, and one New Hampshire HRSA site whose HRSA-reported NPI belongs to its
Vermont parent); two NIH organisations 13. [facility_resolution_summary.csv]

Records attached by each method: single record 297,931; NPPES same name + ZIP 82,015; shared NPI 8,048;
HRSA-organisation name + ZIP 4,954; exact name + ZIP 3,823; contained name + ZIP 1,123; fuzzy name + ZIP 666;
exact name + city 229; contained name + city 97; fuzzy name + city 21. Facilities by weakest link: single record
297,931, high 26,915, medium 4,105, low 68. [facility_resolution_summary.csv]

**V1, ground truth (HRSA `site_npi`).** 5,758 HRSA sites carry an NPI that is an NPPES organisation in scope. With the
shared-NPI links switched off, the name/location rules linked 3,131 of them to some NPPES record: 3,043 to the true NPI
group, which is **precision 0.972**. 3,099 were linked to the true NPI's name (lenient precision 0.990). **Recall was
0.528**. The rules are precise and miss about half of true links, mostly because HRSA site names ("Clover Fork
Evarts") differ from the organisation's legal name. [facility_match_validation.csv]

**V2, location-shuffle false links.** Permuting the locations of ClinicalTrials.gov and NIH records (3 draws) left
13.7 links on average, against 2,959 observed cross-source links. That puts the estimated false-link share at
**0.0046** overall. By method it is 0.053 for contained-name-same-city (3.3 of 63), 0.018 for exact-name-same-city,
0.019 for fuzzy-name-same-city, and 0.001–0.004 for ZIP-level rules. [facility_match_shuffle_control.csv]

**V2 is a lower bound (review note).** A location shuffle can only catch a name that coincides with a record in some
other place. It cannot catch two different entities with similar names in the same city or ZIP. The human-review
sample has such links at the 0.70 tier: "North Houston Sleep Center" with "Houston Sleep Center", and "Puerto Rico
Health LLC" with "Community Health Foundation of Puerto Rico". In the second, the only informative tokens are the
state name. The city-level estimates also rest on 3 draws and 18–203 observed links, so they are imprecise. V1
checks HRSA↔NPPES links only. No ground truth exists for the ClinicalTrials.gov↔NPPES and NIH↔NPPES links that
feed the research registry.

**Negative for recall.** Only 1,239 of the 8,722 research sites were linked to an NPPES or HRSA record.
[facility_resolution_summary.csv] Large academic centres can therefore appear as several facilities: in the San
Diego example, "University of California San Diego" and "REGENTS OF THE UNIVERSITY OF CALIFORNIA" are separate.
A research site's missing clinical-taxonomy information often reflects an unlinked record, not the absence of
clinicians. Per-facility counts are lower bounds.

A seeded sample of 25 links per method is in `facility_match_sample.csv` for human review; no numbers are claimed
from it. The first run's sample (`facility_match_sample_run1_before_rule_fixes.csv`) exposed rule defects, and those
were fixed before any matching or Test 6 analysis (plan §9).

## 3. Registries

`clinic_registry` holds specialty/facility groups (the union of member NPIs' self-reported taxonomies; HRSA sites are
flagged FQHC/look-alike from the HRSA roster), the FQHC evidence source, HRSA site type/setting/location type, county
context (primary-care HPSA coverage, MUA/P coverage, RUCC 2023, nonmetro flag), and individual NPIs per SPEC group
in the same ZIP (**co-location, not affiliation**; not used in ranking). It also holds `npi:` and `hrsa_site:`
object ids.

`research_site_registry` holds, per facility: trials by condition in the literal view and in the MeSH-expanded view
(kept separate: the expanded dysautonomia set includes multiple system atrophy and Parkinson trials); site-recruiting,
overall-recruiting and open-status counts; first and latest start year; interventional/observational counts;
device/diagnostic-test trials; trials whose registered text matches `configs/measurements.yaml` classes and which
classes; NIH projects (precision view and all tiers); core and active core projects; award obligations counted once
per appl_id with `nih_reporter.obligation_total` (test-checked); distinct PI profile ids; activity codes; measurement
terms; `trial:` and `nih:` object ids. It contains no burden variable and no composite score.

**Entry is broad (review correction).** A facility enters `research_site_registry` with any ingested trial or NIH
record. That includes trials that match a target condition only through MeSH or synonym expansion, and NIH links that
are terms-field-only or likely false positives. 1,478 of the 8,722 rows (17%) have only such records. The remaining
7,244 have at least one literal-match trial or one title/abstract NIH match, flagged `has_precision_view_record`.
The aggregate recruiting, recency, interventional, device and measurement counts span all of a facility's ingested
trials, not only the literal-match ones. Per-condition literal counts are in `trials_by_condition_literal`.
`find_relevant_research_centers` filters to the precision view by default. [facility_resolution_summary.csv]

Research-centre counts per condition (`find_relevant_research_centers`) [facility_research_centers_summary.csv]:

| condition | view | facilities | with trials | with NIH | with both | states |
|---|---|---|---|---|---|---|
| long_covid | precision | 458 | 381 | 118 | 41 | 48 |
| me_cfs | precision | 126 | 73 | 63 | 10 | 31 |
| pots | precision | 141 | 134 | 13 | 6 | 32 |
| dysautonomia | precision | 274 | 165 | 132 | 23 | 46 |
| dysautonomia | expanded | 1,191 | 1,071 | 188 | 68 | 49 |

The dysautonomia expanded view is four times the precision view. That is the MeSH-expansion trap documented in the
ClinicalTrials.gov audit, and the literal view is the default.

## 4. Clinic matching — `find_candidate_clinics(geography_id, condition, measurement, radius_km=50, max_sites=15)`

**Inputs.** Geography: county FIPS, `geo:` id, state FIPS, abbreviation or name, or "County, State". Condition: any
alias, via `ontology.normalize.normalize_condition`; only status `matched` is accepted, otherwise UNKNOWN. Measurement:
class or bundle id, label or alias. Precedence is id > bundle label/alias > class name, so "CPET" resolves to the
class and "wearables" to the wearable bundle. Defaults come from `scoring.yaml`.

**Pool.** Geocoded facilities within the radius of the Census internal point, plus facilities inside the county.
HRSA administrative-only sites are excluded. A county with no geocoded facility in its pool returns UNKNOWN.
The pool is drawn from `facilities`, not only `clinic_registry`. It therefore includes research-only records: trial
sites with no NPPES/HRSA link and NIH awardee organisations, some of them companies (e.g. Advanced Brain Monitoring, a
for-profit R43/SBIR awardee, in the San Diego example). These are candidate *study* partners, not clinics. Their `primary_kind` says which kind they are.

**Six characteristics, ranked separately:** clinical_specialty_match, relevant_trial_history, NIH_research_activity,
technology_experience, community_access, distance_to_target_population. Definitions are in plan §4 and in
`CHARACTERISTIC_DEFINITIONS`. A facility is eligible only with a facility-level signal; county flags and distance
never qualify it.

**Selection.** The returned list is a round-robin over the six rankings, and each candidate records `selected_via`.
Each candidate also carries every raw value, its rank and the number ranked per characteristic, plain-language
reasons (with the guardrail caveats inline), `facility:`/`npi:`/`hrsa_site:`/`trial:`/`nih:` ids, and the framing
"This facility has characteristics suggesting it may be a viable implementation or study partner." No combined score
exists. When a pool has no facility with trial history, NIH activity or technology experience, the result lists
`characteristics_absent_in_pool` and an absence note (possible research/measurement desert, not proof of absence).

Examples [facility_candidate_examples.csv]:

| query | pool | eligible | with condition trials | with NIH (condition) | with technology experience |
|---|---|---|---|---|---|
| San Diego County × Long COVID × wearable autonomic monitoring | 2,615 | 1,664 | 6 | 5 | 3 |
| Cook County IL × ME/CFS × wearable autonomic monitoring | 8,135 | 5,159 | 2 | 3 | 0 |
| Harlan County KY × Long COVID × wearable autonomic monitoring | 216 | 184 | 0 | 0 | 0 |
| San Diego County × POTS × autonomic testing | 2,615 | 1,427 | 8 | 0 | 3 |

In Harlan County the 15 candidates were chosen by specialty match, community access and distance only (5 each). The
query returns the absence note. Technology experience for the wearable bundle is rare since the 2026-09-23 integration,
when trial tagging switched to the audited measurement patterns of `configs/measurements.yaml` (brand/device names now
need heart-rate context, bare "ECG" no longer counts, and so on): 3 San Diego facilities (41 with the as-shipped
patterns) and none in the Cook County pool (20 before), whose result now carries the absence note for
technology_experience. San Diego County's internal point lies in the east of the county: its matched
academic sites (UC San Diego, Scripps Research, VA San Diego) are 46–48 km away. **Distance to the internal point is a weak proxy for distance to the population in
large counties.**

## 5. Test 6 — clinic-location shuffle

Design (plan §5): 201 counties with population ≥ 10,000 in the 50 states + DC: 100 metro, 100 nonmetro, plus San
Diego, which is metro, so the metro stratum has 101. Radius 50 km, 15 sites. The location tuple of the 326,685 geocoded
candidate facilities was permuted nationally (200 draws) and within state (100 draws). Checks: pool sizes were
unchanged in every draw (0 violations), and the fast precomputed-slot path reproduced `matching.select` exactly (0
mismatches on 25 geographies × 3 pairs). [test6_clinic_shuffle_checks.csv]

**Statistics (review correction).** The summary holds 240 p-values: 3 pairs × 2 nulls × 10 metrics × 4 strata.
The smallest attainable p with 200 draws is 0.005. The permutation p is also conditional on the particular 201
counties sampled. Two columns were therefore added at review: BH-FDR q-values within each (pair, null) family of 40
tests, and a county bootstrap (4,000 draws, resampling counties within metro/nonmetro) of each county's
observed-minus-null difference. A difference is called **supported** only if q ≤ 0.05 AND the county 95% interval
excludes 0 (`robust_q05_and_geo_ci`). A difference that passes only the unadjusted permutation test is called
**sample-specific**. The review left the observed values, null intervals and p unchanged. The 2026-09-23 integration
did change them (below).

**Integration update (2026-09-23/24).** Trial tagging now uses the audited measurement patterns of
`configs/measurements.yaml` (they were previously applied only by the measurement text miner). Technology experience
for the wearable bundle therefore fell from 8.9% to 2.2% of matched Long COVID sites. Because technology experience is
one of the six characteristics that make a facility eligible and one of the six rankings the round-robin draws from,
the eligible pools and the selected sites changed too, and with them every metric. The design, draws and seeds are
unchanged. The table and readings below are from the regenerated `test6_clinic_shuffle.csv`. The main change in
reading: **the nonmetro technology-experience deficit (0.024 vs 0.051, supported) and the overall deficit (0.089 vs
0.106, supported) do not survive**. With the audited patterns, nonmetro technology experience is 0.003 vs 0.005
(p 0.26) and overall 0.022 vs 0.018 (p 0.13). FQHC/look-alike share and HPSA share also changed status (see (c)).

Primary pair, Long COVID × wearable autonomic/activity monitoring, national null. Each cell gives the observed value,
the null mean [2.5–97.5%], then p / q, then the county-bootstrap 95% CI of observed − null. ✔ = supported.
[test6_clinic_shuffle.csv]

| metric | all 201 | metro (101) | nonmetro (100) | metro − nonmetro |
|---|---|---|---|---|
| share with ≥1 Long COVID trial | 0.038 vs 0.038 [0.031, 0.045]; p 0.75 / q 0.79; CI [−0.006, +0.009] | 0.071 vs 0.063 [0.051, 0.077]; p 0.20 / q 0.24; CI [−0.006, +0.022] | 0.005 vs 0.012 [0.007, 0.017]; p 0.04 / q 0.061; CI [−0.010, −0.002] | 0.066 vs 0.052 [0.037, 0.066]; p 0.06 / q 0.085; CI [−0.000, +0.029] |
| mean Long COVID trials per candidate | 0.078 vs 0.065 [0.052, 0.081]; p 0.07 / q 0.093; CI [−0.004, +0.032] | 0.150 vs 0.111 [0.084, 0.139]; p 0.015 / q 0.033; CI [+0.006, +0.077] ✔ | 0.005 vs 0.018 [0.009, 0.029]; p 0.035 / q 0.058; CI [−0.017, −0.008] | 0.145 vs 0.093 [0.066, 0.121]; p 0.005 / q 0.013; CI [+0.018, +0.089] ✔ |
| share with NIH Long COVID project | 0.018 vs 0.013 [0.009, 0.017]; p 0.03 / q 0.052; CI [+0.000, +0.010] | 0.033 vs 0.023 [0.016, 0.030]; p 0.015 / q 0.033; CI [+0.002, +0.020] ✔ | 0.002 vs 0.003 [0.001, 0.006]; p 0.52 / q 0.56; CI [−0.003, +0.001] | 0.031 vs 0.020 [0.012, 0.027]; p 0.01 / q 0.025; CI [+0.003, +0.021] ✔ |
| share with technology experience | 0.022 vs 0.018 [0.013, 0.023]; p 0.13 / q 0.16; CI [−0.002, +0.010] | 0.040 vs 0.031 [0.022, 0.041]; p 0.065 / q 0.089; CI [−0.003, +0.022] | 0.003 vs 0.005 [0.002, 0.009]; p 0.26 / q 0.30; CI [−0.004, +0.001] | 0.038 vs 0.026 [0.018, 0.037]; p 0.02 / q 0.038; CI [−0.001, +0.024] |
| share FQHC/look-alike | 0.225 vs 0.247 [0.234, 0.261]; p 0.005 / q 0.013; CI [−0.043, −0.000] ✔ | 0.229 vs 0.286 [0.269, 0.304]; p 0.005 / q 0.013; CI [−0.086, −0.027] ✔ | 0.220 vs 0.207 [0.189, 0.226]; p 0.12 / q 0.16; CI [−0.019, +0.046] | 0.009 vs 0.079 [0.054, 0.100]; p 0.005 / q 0.013; CI [−0.116, −0.026] ✔ |
| share in primary-care HPSA county | 0.822 vs 0.812 [0.803, 0.820]; p 0.02 / q 0.038; CI [−0.001, +0.021] | 0.769 vs 0.750 [0.736, 0.766]; p 0.02 / q 0.038; CI [+0.002, +0.037] ✔ | 0.876 vs 0.875 [0.866, 0.884]; p 0.92 / q 0.92; CI [−0.012, +0.012] | −0.107 vs −0.126 [−0.142, −0.105]; p 0.04 / q 0.061; CI [−0.002, +0.041] |
| share with condition specialty group | 0.601 vs 0.643 [0.629, 0.655]; p 0.005 / q 0.013; CI [−0.067, −0.019] ✔ | 0.591 vs 0.594 [0.575, 0.613]; p 0.78 / q 0.80; CI [−0.030, +0.026] | 0.610 vs 0.692 [0.672, 0.711]; p 0.005 / q 0.013; CI [−0.121, −0.046] ✔ | −0.018 vs −0.098 [−0.123, −0.067]; p 0.005 / q 0.013; CI [+0.034, +0.128] ✔ |
| km to nearest Long COVID trial site | 108 vs 91 [79, 105]; p 0.025 / q 0.045; CI [+8, +26] ✔ | 81 vs 71 [59, 84]; p 0.10 / q 0.13; CI [+1, +21] | 136 vs 112 [99, 130]; p 0.005 / q 0.013; CI [+10, +39] ✔ | −55 vs −41 [−54, −28]; p 0.045 / q 0.066; CI [−32, +4] |

Jaccard overlap of the shuffled matched set with the observed one: 0.000 (national) and 0.022 (within state).

Reading against the pre-specified criteria:

* **(a) Location drives the output: yes, but trivially.** Shuffled matched sets share almost no facilities with the
  observed ones. A national permutation moves nearly every facility away from any fixed point, so this criterion holds
  by construction. It is a sanity check, not evidence.
* **(b) Research infrastructure is concentrated beyond facility density: supported for trial counts and NIH activity;
  not supported for technology experience.** The metro–nonmetro contrast in Long COVID trials per candidate (0.145 vs
  0.093) and in NIH share (0.031 vs 0.020) survives both the multiplicity adjustment and county resampling. Both
  effects are small in absolute terms: about 0.05 trials, or 1 NIH-active facility in 90, per candidate. The null
  contrast is itself positive because facility density follows population. Nonmetro targets are farther from the
  nearest Long COVID trial site than random placement would put them (136 km vs 112 km, supported). The
  technology-experience **measurement desert** reported before the integration is not supported with the audited
  patterns: nonmetro 0.003 vs 0.005 (p 0.26). Technology experience is now so rare (2.2% of matched sites) that the
  shuffle cannot separate it from random placement. Contrasts that are not supported:
  * technology experience, metro − nonmetro: 0.038 vs 0.026, p 0.02 / q 0.038, but the county CI includes 0
    (sample-specific; it is supported under the within-state null, below)
  * share with a Long COVID trial, metro − nonmetro: p 0.06 / q 0.085
  * km to the nearest trial site, metro − nonmetro: q 0.066, and the county CI includes 0
* **Null results.** The national average share of matched facilities with Long COVID trial history is unchanged by
  the shuffle (0.038 vs 0.038, p 0.75). The shuffle redistributes research presence between places; it does not
  remove it from the average. Also null at the metro level: share with a trial (p 0.20) and technology experience
  (p 0.065 / q 0.089). The nonmetro deficits in trial counts (p 0.035–0.04) do not survive the multiplicity
  adjustment (q 0.058–0.061), although their county intervals exclude 0.
* **(c) Community access: partly.** Observed matched sets hold fewer FQHC sites than random placement would give,
  overall (0.225 vs 0.247, supported, county CI upper bound −0.000) and in metro counties (0.229 vs 0.286, supported);
  the nonmetro FQHC share (p 0.12) and nonmetro HPSA share (p 0.92) are inside the null. The overall HPSA
  difference (0.822 vs 0.812, p 0.02 / q 0.038) is sample-specific: its county interval includes 0. (Before the
  integration the overall FQHC difference was sample-specific and the overall HPSA difference supported.)
* Observed matched sets carry condition-specialty groups less often than randomly placed facilities would, overall
  (0.601 vs 0.643) and in nonmetro counties (0.610 vs 0.692). Nonmetro matched sets are also closer to the internal
  point (21.7 vs 23.6 km). All are supported. Primary care is a curated "core" specialty for every condition
  (`configs/relevance.yaml`), so specialty match is a weak signal.
* **Secondary pairs** (national null; counts of supported metro−nonmetro contrasts) [test6_clinic_shuffle.csv]:
  * POTS × autonomic function testing: 10 of 10 supported (e.g. share with a POTS trial 0.048 vs 0.024, q 0.009).
  * ME/CFS × wearable: 8 of 10 (7 of 10 before the integration). Not supported: mean ME/CFS trials per candidate
    (p 0.07) and km to the nearest ME/CFS trial site (p 0.52). Technology experience (0.039 vs 0.027, q 0.042) and
    the share-with-trial contrast (0.024 vs 0.014, q 0.037) are supported only at their lower CI bound (+0.000).
    ME/CFS trial sites are rare: the overall share is 0.012 vs 0.009 (p 0.11, null).
  * Several of the 10 contrasts are not research metrics (pool eligibility, distance, FQHC, HPSA). "All 10" is not
    ten independent research findings.
* **Within-state null** (Long COVID): the same directions hold. Supported: metro−nonmetro trials per candidate
  (0.145 vs 0.097), the NIH contrast (0.031 vs 0.020), the share-with-trial contrast (0.066 vs 0.051), the
  technology-experience contrast (0.038 vs 0.025, q 0.019, CI [+0.002, +0.023]) and metro technology experience
  (0.040 vs 0.030), and nonmetro trials per candidate (0.005 vs 0.019, q 0.032). Not significant: the nonmetro
  technology-experience deficit (0.003 vs 0.005, p 0.29), which was supported before the integration.

Per-geography observed and null values: `test6_clinic_shuffle_geographies.csv`. Per-draw values:
`test6_clinic_shuffle_draws.csv`. Draft figure: `results/figures/drafts/test6_clinic_shuffle.png`. `outside_null_95`
(percentile interval) and the pre-specified p can disagree for metrics bounded at 0. Example: POTS nonmetro trials
per candidate is 0.000 against [0.001, 0.026], with p 0.06.

## 6. Limitations

* NPPES taxonomies are self-reported and many are stale. Organisation NPIs are billing entities; group practices are
  not facilities. A taxonomy does not mean a facility sees the condition.
* Locations are ZIP centroids (NPPES, RePORTER, most trial sites) or city centroids (308 trial facilities). Distances
  are approximate, and the county internal point is a poor stand-in for where people live in large counties.
* Entity resolution is precise but misses about half of true links (V1 recall 0.528). Research sites are mostly
  unlinked to clinical records, and spelling variants that could not be linked stay separate.
* Condition specialties and measurement implementer groups are curated assumptions (`configs/relevance.yaml`).
  Technology experience means a registered trial mentioned the measurement, not that the measurement works.
* Trial sets are high-precision samples (quoted-phrase search); NIH counts use the precision view. FY2026 RePORTER
  is incomplete.
* `research_site_registry` is written as the canonical table by `measure_it.facilities.registry`. The
  ClinicalTrials.gov per-site-key history it once sat next to (`research_site_registry__clinicaltrials_gov`) was
  renamed `ctgov_facility_summary` on 2026-09-23, so no `research_site_registry__*` partition exists and
  `store.union_partitions("research_site_registry")` can no longer overwrite the registry (the pipeline's union step
  also skips it).
* Entity-resolution precision (V1 0.972) is measured on HRSA↔NPPES links only. The ClinicalTrials.gov and NIH links
  that feed the research registry have only the V2 location-shuffle bound, which cannot see co-located look-alike
  names.
* The measurement resolver (since 2026-09-24 the shared `measure_it.measurements.resolve`, also used by the evidence
  tools and the agent) accepts a single word when it is the only label containing it at ≥ 50% of that label's
  words. For example, "monitoring" resolves to the wearable bundle. The match reason is returned with the result.
  Ambiguous words ("testing", "imaging", "heart rate") return UNKNOWN with the candidate ids.
* Test 6 conclusions are about the 201 sampled counties and one radius (50 km) and site count (15). The county
  bootstrap covers county sampling; radius and site count were not varied.
