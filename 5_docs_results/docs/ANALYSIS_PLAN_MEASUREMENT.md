# ANALYSIS PLAN — Measurement / technology discovery (SPEC §D, Phase 3)

Written 2026-09-23, **before** the text mining was run. Deviations made after this point are
listed at the end of the file, each with its date and reason.

Owner modules: `src/measure_it/measurements/{discovery,evidence,query,precision_audit}.py`.
Outputs: `measurement_mentions`, `measurement_registry`, `condition_measurement_evidence`,
`measurement_phenotype_axis_map` (data/processed), `results/tables/measurement_*.csv`,
`results/MEASUREMENT_DISCOVERY.md`, draft figures in `results/figures/drafts/`.

## 1. Questions

* **Q1 (discovery).** Which objective measurement classes (configs/measurements.yaml, 26 classes, 4
  bundles in configs/relevance.yaml) do registered trials and NIH-funded projects actually deploy for
  each target condition, and in which role (outcome measure, intervention, eligibility, description)?
  Primary view: the demo cluster (long_covid, me_cfs, pots, dysautonomia).
* **Q2 (validity of the miner).** How precise is the pattern-based miner, per class, when judged
  snippet by snippet? Where precision is below 0.8, does tightening the patterns fix it, judged on
  mentions that were not used to design the change?
* **Q3 (dimensions).** For each class and bundle: measurement_evidence_strength, technology_maturity,
  regulatory_visibility, deployment_complexity (curated), phenotype_signal_strength (left pending),
  each with its components, never combined into one number.

What these outputs are **not**: registration of a measurement in a trial, or a grant mentioning it,
is a research-activity / deployment-interest signal. It is not evidence that the measurement works,
detects the condition, or is valid for it. FDA records are deployment-readiness signals, not
evidence that a device diagnoses the target illness (CONVENTIONS §3.7).

## 2. Inputs and record universes

| Input | Use | Precision view (primary) | Sensitivity view |
|---|---|---|---|
| clinical_trials, trial_interventions, trial_outcomes | text to mine | — | — |
| trial_conditions | condition assignment | `condition_literal_match = True` (a searched term literally in conditions/keywords/titles/acronym) | all rows (adds MeSH/synonym expansion; for dysautonomia and mcas this is mostly other diseases, per the ctgov DATA_AUDIT) |
| nih_projects | text to mine | — | — |
| nih_project_conditions | condition assignment | `match_tier = title_abstract` and not `likely_false_positive` (the helper default; removes PASC = pancreatic stellate cells) | all tiers incl. terms-field-only |
| measurement_regulatory_status, fda_510k, fda_pma | FDA dimensions | curated product-code map, mapping_confidence kept | — |
| participant_* tables, modality_inventory__mapmecfs | which public person-level datasets observe each class | counts measured from the tables | — |
| geo_study_catalog, modality_inventory__mapmecfs | molecular assay link for omics classes | `study_type_in_scope = True` | — |
| phenotype_axes | phenotype filter in `discover_candidate_measurements` | — | — |

Grants: one RePORTER appl_id is one fiscal-year award record; multi-year projects repeat their
abstract. Counts are reported as distinct **core projects** (primary) and appl_ids (secondary).

## 3. Mining design

* **Fields and roles.** Trials: brief/official title, brief summary, keywords -> role `abstract`
  (the registry's study description); intervention name / description / other names -> `intervention`;
  outcome measure / description -> `outcome_measure` (outcome_type kept: primary/secondary/other);
  eligibility criteria -> `eligibility`. Grants: project title, abstract, public-health relevance ->
  `abstract`. RePORTER's thesaurus `terms` field is **not** mined (concept-level, low precision).
* **Mention.** One row per (record, field instance, measurement class) with >= 1 pattern hit; the
  snippet is +-150 characters around the first hit, with the hit in [[double brackets]]; number of
  hits and matched terms kept. Deterministic `mention_id` = `<source>|<record>|<field>|<sub_idx>|<measurement_id>`.
* **Objective vs questionnaire/PRO.** A mention is `objective_flag = False` when, within its local
  window (the whole field if <= 250 characters, else the sentence around the hit capped at +-200
  characters), (a) a class `exclude_pattern` matches, or (b) the PRO regex matches: questionnaire,
  survey, inventory, scale (not large-/small-/full-scale, scale-up, time scale), PROMIS, SF-36/12,
  RAND-36, EQ-5D, COMPASS-31, DePaul, Chalder, FSS, MFI, PHQ, GAD-7, HADS, BDI, PSQI, ISI, Epworth,
  self-report(ed), patient-reported, PRO/PROM, visual analog(ue), VAS, Likert, diary/diaries, rating,
  checklist, interview, and score-only phrases ("total/global/composite/symptom/severity/subscale/
  summary score", "score on/of the"). The flag and the term that fired are stored. This is
  deliberately conservative: a sentence that lists a device measure next to a questionnaire is marked
  non-objective. Its false-exclusion rate is audited (section 5).
* **Counting rules.** A trial "uses a class objectively" when it has >= 1 objective mention in role
  `outcome_measure`, `intervention` or `abstract`. Eligibility mentions are stored and counted
  separately (`n_trials_eligibility_only`) but do not enter headline counts: they are dominated by
  exclusion criteria ("abnormal ECG at screening"). "As outcome measure" = >= 1 objective mention in
  role `outcome_measure`. Recruiting = overall_status RECRUITING, NOT_YET_RECRUITING or
  ENROLLING_BY_INVITATION. Recent start = start_year >= 2021 and <= 2035 (placeholder years excluded).
  Grants: >= 1 objective mention in title/abstract/relevance.
* **Bundles** are resolved to member classes; bundle counts are the union of member records
  (never sums).

## 4. Dimensions (registry and condition x measurement table)

Every component is its own column. Each dimension also gets a tier with the rule written next to it;
tiers are within-dimension summaries and are never combined across dimensions.

* **measurement_evidence_strength** — components: n trials objective, n as outcome measure,
  n recruiting, n recent start (>= 2021), n distinct lead sponsors, n NIH core projects (and appl_ids).
  It measures research deployment activity, not efficacy or validity. Tier (registry, 14 conditions,
  literal view): high >= 50 objective trials, moderate 10-49, low 1-9, none 0. Tier (one condition):
  high >= 10, moderate 3-9, low 1-2, none 0.
* **technology_maturity** — components: n mapped product codes (by mapping confidence), n 510(k),
  n De Novo, n PMA (from measurement_regulatory_status, exact API counts), highest device class,
  first and most recent decision date, years since first clearance, n consumer-device clearances
  (OTC/wellness product codes, or consumer-electronics applicants: a curated list stated in code).
  Tier: established = a high/medium-confidence code with >= 10 decisions and a first decision >= 10
  years before 2026-09-23; emerging = decisions on record but not established; category_only =
  product code(s) with no decisions; no_fda_category = no product code (LDTs, research assays).
* **regulatory_visibility** — components: best regulatory_visibility value, best mapping confidence,
  n codes per confidence level, max US registered establishments on one code. Tier: visible = a
  high-confidence code with decisions; partial = decisions only under medium/low-confidence or
  filtered generic codes, or categories without decisions; not_visible = no product code found.
* **deployment_complexity** — curated ordinal 1-5 with reason and implementer groups from
  configs/relevance.yaml (evidence_type curated_config; labelled "curated assumption").
* **phenotype_signal_strength** — null, status `pending_phenotype_results` (filled later from the
  wearable results by another step).

## 5. Precision audit (Q2)

* **Population.** Objective mentions in counted roles (trial abstract/intervention/outcome; grant
  abstract) whose record is linked to a demo-cluster condition in the precision view. Grants are
  de-duplicated to the latest fiscal-year appl_id per core project before sampling.
* **Sample.** Per class, a seeded random sample (numpy default_rng, seed = config.SEED + stable class
  offset) of 30 mentions; classes with fewer than 30 such mentions are audited in full (census) and
  labelled as such; classes with 0 are reported as "no mentions".
* **Judging rubric** (applied to the snippet only, by an AI reviewer, not a domain expert):
  `yes` = the snippet describes objective use (performed, planned, or used as a test/criterion) of a
  measurement within the class **as scoped by its name in measurements.yaml** (e.g. a resting 12-lead
  ECG is outside "ECG (ambulatory/wearable/patch)"; a clinic spot-check vital sign is outside
  "continuous/overnight pulse oximetry" and "wearable heart rate"; a stopwatch sit-to-stand is outside
  "digital gait"); `no` = homonym/other meaning, questionnaire, drug/intervention name, background or
  negated statement, or a measurement outside the class scope; `unclear` = cannot tell from the snippet.
* **Metrics.** strict precision = yes / n (unclear counts against); lenient precision = yes / (yes+no);
  Wilson 95% CI for both. **Trigger:** strict precision < 0.8 -> tighten that class's patterns.
* **Tightening without tuning on the test set.** Changes are designed from the round-1 snippets, the
  whole corpus is re-mined, and the class is re-audited on a **fresh** seeded sample of 30 mentions
  from the tightened set that excludes every mention judged in round 1. The round-2 figure is the one
  reported. At most two rounds; a class still below 0.8 after round 2 is reported as such and flagged
  in the registry (`text_mining_precision_*` columns) rather than tuned further.
* **Tightened patterns** cannot be written into configs/measurements.yaml (shared file). They live in
  `discovery.PATTERN_REVISIONS` with a per-change rationale, are logged to
  `results/tables/measurement_pattern_changelog.csv`, and are reported back as a requested change to
  the shared config.
* **PRO-filter audit.** A seeded sample of 30 demo-cluster mentions flagged non-objective, judged for
  whether the PRO/exclude flag was wrong (false-exclusion rate).

## 6. Controls and sensitivity analyses

1. **Condition-assignment sensitivity:** every condition-level count is also computed with the
   MeSH/synonym-expanded trial set and with all RePORTER match tiers; the difference is reported,
   prominently for dysautonomia and mcas.
2. **Comparator conditions:** for each class, the share of demo-cluster literal trials that use it vs
   the share of other-condition literal trials (ratio, with Wilson intervals for each share;
   descriptive, not a hypothesis test; 26 classes, so multiplicity applies).
3. **Pattern-version sensitivity:** counts under the original (v1) and tightened patterns side by side.

## 7. Query layer

* `discover_candidate_measurements(condition, phenotype=None)`: condition resolved with the ontology
  module's `normalize_condition`; phenotype resolved against phenotype_axes (axis id, axis label, HPO
  id/label); unresolved input -> UNKNOWN with a reason, never a silent default. Candidates = classes
  with >= 1 objective trial or grant for the condition. **Ordering rule (lexicographic, transparent):**
  (1) n trials with objective use (literal view), (2) n trials using it as an outcome measure,
  (3) n NIH core projects, (4) technology_maturity tier, (5) lower deployment_complexity, (6) id.
  phenotype_signal_strength is not used while pending. With a phenotype, candidates are restricted to
  classes that observe that axis in `measurement_phenotype_axis_map`, whose basis is a text match
  between the axis's curated wearable_signals and the class patterns/signals/name, the axis's
  requires_accelerometry / requires_heart_rate flags, or a curated clinic-reference link (labelled).
* `get_measurement_evidence(measurement, condition)`: measurement resolved as class id/name, bundle
  id/label/alias, then the regulatory module's pattern/fuzzy matcher; returns the evidence row(s) with
  object_ids (trial:, nih:, fda:, geo_series:) and example snippets; UNKNOWN when nothing matches.

## 8. Figures (drafts)

1. condition x measurement heatmap of objective trials (literal view), demo cluster highlighted;
2. precision audit per class, round 1 and round 2, with Wilson CIs and the 0.8 line;
3. role mix (outcome / intervention / description) per class in the demo cluster.

## 9. Known limitations fixed in advance

Pattern mining misses measurements described in words not in the patterns (recall is not measured
except through the PRO-filter audit); registry text describes planned work; grants mention
measurements in background sentences; the AI reviewer is not a domain expert and judges snippets
without the full record.

## Deviations from this plan

1. **2026-09-23, after round 1 — non-human context flag (all classes).** The round-1 'no' judgements included
   animal and cell work (mostly in grant abstracts) across several classes. A mention whose local window names an
   animal or cell model (discovery.NONHUMAN_RE) is kept in measurement_mentions but does not count as evidence
   (`evidence_flag = objective_flag AND no non-human context`). Covered by the round-2 audit.
2. **2026-09-23, before round-2 sampling — QST rejected.** 'quantitative sensory testing' was drafted as a
   small_fiber_testing pattern and removed before the round-2 sample was drawn: its mentions were pain-threshold
   testing in fibromyalgia/migraine trials, not small-fiber nerve testing.
3. **2026-09-23, while judging round 2 — sentence-window bug fixed.** RePORTER abstracts are hard-wrapped (90% contain
   single line breaks), so a single newline was being treated as a sentence boundary and the PRO / non-human window for
   grants was one wrapped line. The plan specifies "the sentence around the hit"; a single newline is now not a
   boundary (a paragraph break or a newline followed by a list marker still is). The round-2 sample was not re-drawn;
   the summary also reports precision re-scored on the judged mentions still in the population
   (`precision_strict_current_population`).
4. **continuous_temperature** had no fresh post-revision mention to audit in round 2 (round 1 was a census). Its
   reported precision is the round-1 judgements of mentions surviving the revision, labelled as design data
   (`text_mining_precision_fresh_sample = False`).
5. **PRO filter not retuned.** Its false-exclusion rate was measured (results/tables/measurement_pro_filter_summary.csv)
   but the filter was not changed, because retuning it on the audited sample would have tuned on test data. A
   sensitivity count including PRO-flagged mentions is reported (`n_trials_any_counted_mention_sensitivity`).
6. **Additional outputs not in the plan:** `results/tables/measurement_manual_functional_tests.csv` (6MWT,
   sit-to-stand, TUG and gait-speed counts, which fall outside every class) and a `report.py` module that writes
   results/MEASUREMENT_DISCOVERY.md from the tables.
7. **2026-09-23, independent review (after all judging was complete; no pattern, sample or judgement changed).**
   * *Precision status and record-clustered intervals.* Judged mentions cluster within records (e.g. 8 of 30 round-2
     vascular_imaging mentions are one trial; 30 round-1 retinal_imaging mentions come from 9 records), so the
     mention-level Wilson interval is too narrow. `precision_audit.clustered_wilson` adds a conservative interval
     (worst-case Kish effective sample size (sum m)^2 / sum m^2 over records). A class now "meets 0.8" only with
     >= 10 judged mentions from >= 5 distinct records (`text_mining_precision_status`); respiratory_rate (1 fresh
     mention), small_fiber_testing (2), digital_gait (6), remote_patient_monitoring (8), ppg (10 from 4 trials) and
     capillaroscopy (10 from 4 trials) are `insufficient_sample`, not passes. The point-estimate flag
     `text_mining_precision_below_0_8` is unchanged.
   * *Dysautonomia broad-term sensitivity.* 144 of 298 literal-match dysautonomia trials link only through
     'autonomic dysfunction' (secondary autonomic dysfunction in diabetes, CKD, heart failure, sepsis, anaesthesia,
     Parkinson disease), the phrase HRV studies use for their outcome. `*_specific_terms` columns
     (evidence.BROAD_LITERAL_TERMS) report counts through the condition's other literal terms; the primary view is
     unchanged.
   * *HRV in the wearable bundle.* HRV is device-agnostic; bundles now also report `n_trials_objective_excluding_hrv`.
   * *Comparator.* Katz log intervals for the demo/other ratio, and the ratio with dysautonomia's broad-term-only
     trials dropped.
   * *Fixes.* Samsung Medison (hospital ultrasound) no longer counts as a consumer-device applicant; filtered subsets
     of a product code that a bundle also maps whole are no longer double counted in FDA decision totals; example
     object ids are written deterministically; `resolve_phenotype` returns UNKNOWN when several axes tie
     ('orthostatic', 'intolerance', 'sleep'); `resolve_measurement` gives an exact class id precedence over a bundle
     alias ('capillaroscopy' is the class) and returns UNKNOWN when only loose word matches point to several classes.
   * *Integration (2026-09-23, after the audit).* The tightened patterns were moved from `discovery.PATTERN_REVISIONS`
     into `configs/measurements.yaml` unchanged (per-class `audit_revision` blocks keep the as-shipped patterns as
     `patterns_v1` and the rationale per change), so the config is the single source of truth: the ClinicalTrials.gov
     trial tagging, the facility registry and the regulatory matcher (`get_regulatory_context`, which now also
     refuses questionnaire / patient-reported queries) use the audited patterns too. `discovery.PATTERN_REVISIONS`
     is derived from the config with the same content, so every mention's `pattern_version` (sha1 12e9f79f) and
     `measurement_mentions` are unchanged. The salivary-cortisol matches of `blood_biomarkers` (4 of 5 round-2 'no'
     judgements) are recorded as a `review_note`; the pattern was not changed after round 2.
