# Guardrail audit (2026-09-24)

An adversarial audit of the whole repository against the SPEC sections **NON-NEGOTIABLE SCIENTIFIC GUARDRAILS** and
**KEY PRODUCT LANGUAGE** (docs/SPEC.md) and docs/CONVENTIONS.md section 3. It covered code, processed tables, tool
outputs, MCP and API descriptions, agent answers, dashboard pages, figure text and docs. Every violation found was
fixed where it originated. This file lists each guardrail, where the code enforces it, how it was tested, the result,
and what remains a limitation.

**Summary.** No person-level row is joined to a geographic or facility row. No person-level table carries a location.
No omics are attributed to wearable participants. No level-D condition has a burden value. No state estimate is
presented as county prevalence. No forbidden phrase appears in any user-facing output tested. The audit found and
fixed eight problems:

| # | problem | where it came from | fix |
|---|---|---|---|
| 1 | `get_patient_phenotype_signature` guessed from partial or ambiguous condition text: `chronic` returned the migraine signatures, `headache` migraine, `syndrome` fibromyalgia + IBS, `chronic fatigue` the ME/CFS-like proxy | `wearables/signatures.py:_resolve_condition` used every match of a `partial`/`ambiguous` normaliser result; phenotype-text matching accepted generic words (`syndrome` -> "irritable bowel syndrome") | only an unambiguous `matched` result selects conditions, and the candidates are reported but not used (`candidates_not_used`). A query made only of generic label words (`GENERIC_LABEL_TOKENS`) does not text-match. Exact inputs (ME/CFS, PASC, depression, fatigue, ...) answer as before |
| 2 | `normalize_condition` returned envelope status `ok` for `partial`/`ambiguous` matches (e.g. `chronic` -> migraine), which a calling model could take as the answer. docs/MCP.md already said these "are not answers" | `tools.py:normalize_condition` | status `UNKNOWN / NOT AVAILABLE` with a reason. The candidates stay in `data.matches` / `data.other_candidates`, and none is chosen |
| 3 | The Condition Explorer showed a full page for a guessed first match (typing `chronic` showed migraine, `Lyme` showed Lyme disease). docs/DASHBOARD.md already said the page stops | `dashboard/pages/2_Condition_Explorer.py` | a partial or ambiguous match stops the page and lists the candidates. Nothing below the input is shown |
| 4 | `get_regulatory_context("wearable autonomic monitoring")` matched the bundle name to `autonomic_testing` (tilt table) by fuzzy word overlap (score 0.4). The bundle id `wearable_autonomic_activity_monitoring` returned UNKNOWN | `measurements/regulatory.py:match_measurement_classes` does not know bundles; the facade passed the text through | `tools.py:_regulatory_for_bundle`: when there is no exact or pattern class match, a bundle name resolves to its member classes (exact matches). When every match is fuzzy, a caveat says the classes are candidates, not an identification |
| 5 | For a level-D member of a condition set, `burden_source_resolution` and `burden_measure_id` read `"nan"` (e.g. POTS in the demo cluster). That is a string where the value is unknown | `scoring/opportunity.py:combo_frame` wrote `str(NaN)` | `UNKNOWN` / `None` at the source. `scoring/recommend.py:recommendation` normalises rows stored before the fix |
| 6 | Agent answers accepted out-of-scope premises without comment. For "Does the patient have POTS, and which clinic is best ...", the answer had no forbidden phrase but did not say that the engine does not assess individuals or rank clinics as best. For "multi-omics profile of NHANES participants living in San Diego County ... which gene causes it", it stopped without saying why the premise is impossible | `agents/deployment_agent.py:answer_question` | `scope_notes()` adds a fixed **Scope (what this engine does not do)** section when the question asks about an individual or diagnosis, which clinic is best, participant omics, causes or treatment, or where participants live. "Measurement: none named" became "none recognised" when a named device (e.g. "smartwatch") is not in the vocabulary |
| 7 | The geography Test 4 table listed `long_covid_county_proxy` and `cms51_county` with variable `prevalence` but without their measure or burden evidence level. A reader could take PLACES PHLTH (a level-C proxy) for long-COVID prevalence | `geography/features.py` report writer | `test4_series_legend()` prints each series' measure and level under the table. results/GEOGRAPHY_RESULTS.md now carries the same text the generator writes |
| 8 | Hand-typed numbers in dashboard text ("1,000 Monte Carlo draws", "x2 for an inherited value", "eight named weight sets") could go stale if the config changed | `dashboard/pages/5_Deployment_Recommendation.py` | read from the ranking output (`over N draws`) and `configs/scoring.yaml`. The rendered values are unchanged: 1,000 draws, x2, 8 weight sets |

Also changed: the `find_candidate_clinics` docstring. It said `never a "best clinic"`, which is excused in Markdown but
reaches the OpenAPI and MCP JSON with escaped quotes. It now says "never a quality ranking of facilities". The
`measure-it validate` language check now covers user-facing text it did not scan before (see guardrail 11).

---

## How the audit tested

| probe | what | result |
|---|---|---|
| Tool sweep | 708 calls to all 18 functions in `measure_it.tools`. Condition inputs were nonsense (`asdkjhqwe`, empty, SQL-like, 300 characters, `None`), level-D (POTS, endometriosis, IBS, MCAS, gastroparesis, Ehlers-Danlos, dysautonomia, post-infectious syndrome), ambiguous (fatigue, syndrome, chronic, pain, covid, Lyme, headache, tick) and aliases or codes (PASC, CFS, SEID, U09.9, G93.32, M79.7, MONDO/HPO ids, "long haul covid"). There were 12 measurement inputs, 17 geographies (Springfield, Washington County, Atlantis, 99999, zcta:92101, ...) and 19 object ids (participant ids, unknown NCT/NPI/K numbers, level-D opportunities). Each output was scanned with `validate.forbidden_hits` | before the fixes: 289 `ok`, 419 UNKNOWN, 0 errors, 0 forbidden terms. 12 of the `ok` answers were guesses: 8 `normalize_condition` partial/ambiguous matches (syndrome, chronic, covid, Lyme, post-infectious, headache, chronic fatigue, HP:0012378) and 4 phenotype signatures chosen through them (syndrome, chronic, headache, chronic fatigue). `wearable autonomic monitoring` got a fuzzy regulatory match, and the bundle id got UNKNOWN (problem 4). After: 278 `ok`, 430 UNKNOWN, 0 errors, 0 forbidden terms. All 12 guesses return UNKNOWN, and bundle names resolve to their member classes. The remaining `ok` answers to ambiguous text are `search_condition` candidate lists (6, by design), phenotype-text matches flagged `matched_by: "phenotype text"` (`fatigue`, `covid`), `get_phenotype_measurement_evidence("fatigue")` (fatigue is an axis of exactly one demo phenotype), and empty optional filters (5) |
| Agent | `uv run measure-it ask` on 5 questions: an individual-diagnosis / best-clinic / smartwatch question; capillaroscopy for two level-D conditions; the SPEC deployment question; a participant multi-omics + location + causation question; nonsense | no forbidden term. Level-D burden stays UNKNOWN (evidence level D). Inherited state values are labelled. Nothing is ranked without a resolved condition. The Scope section was added (problem 6) |
| API | 17 FastAPI endpoints via `TestClient`: `/`, `/health`, `/tools`, `/sources`, `/openapi.json`, one source, 11 tool calls | only hit: the quoted docstring phrase in `/openapi.json` and `/tools` (fixed). "Proof", "validated", "multi-omics" and "prevalence" appear only in negated caveats |
| Dashboard | 22 `AppTest` runs across the home page and 5 pages: defaults, the long-COVID county burden layer, POTS, and adversarial inputs (chronic, syndrome, Lyme, POTS, nonsense, smartwatch, qwerty, Springfield, Texas, Washington County). All text, dataframe cells and plotly specs were scanned | no exception, no forbidden term. The Condition Explorer guessed on partial input (problem 3, fixed) |
| Code | searched every module that reads a person-level table for merges or joins with geographic or facility tables. Scanned every person-layer parquet for location-like columns (`validate` regex plus tz, site, city, region, address, clinic, facility, npi, census, rural, ...) | person tables are joined only within a dataset on its own participant id (NHANES SEQN; Stanford native id within a sub-dataset; MapMECFS NIH study number, omics x omics). `measurements/evidence.py` reads person tables only for counts. No location columns |
| Tables | `geo_condition_burden` (558,246 rows), `geo_condition_features`, `deployment_opportunities` (76,680 rows), `deployment_candidates`, `geo_burden_definitions`, `digital_phenotype_datasets`, `SOURCE_REGISTRY.yaml` | every burden row has a level. 0 level-D rows carry a value. All 3,144 county rows at state resolution are `inherited = true` and say "not county prevalence" in `measure_label`. No level-D condition has a burden, burden count, proxy value or burden percentile in features or opportunities |
| Docs and figures | searched README, docs, results, CAPTIONS, data audits, configs and HTML maps for the forbidden terms and wider variants (proven, cure, patients have, definitive, diagnose, optimal, "the best", FDA-cleared for, validated for, diagnostic accuracy); read every "multi-omics" and "prevalence" mention in context; spot-checked hand-written numbers in docs/MCP.md, docs/DASHBOARD.md and docs/SCORING.md against the tables | all mentions are negations or quoted rules. The numbers checked match the current tables (Oklahoma County composite 0.861, interval 1-306; San Diego long-COVID 3.8 %; Oglala Lakota equal-weight rank 443, burden-only 307) |

Regression tests for every fix are in `tests/test_guardrail_audit.py` (13 tests). The full suite (`uv run pytest -m "not network"`, 1,190 tests) ran after the fixes with 1 failure. The language check caught a forbidden phrase in this document's first draft; it was reworded, and the affected test files (guardrail audit, tools, pipeline, MCP server: 90 tests) were re-run and pass. They cover the pure scope-note,
Test 4 legend and language-check logic, and data tests for the signature, normaliser, regulatory, level-D member,
agent-scope and Condition Explorer fixes.

---

## Guardrails, enforcement and results

### 1. Four data layers; a participant is never the same individual as anyone in the geographic or facility layers

* **Enforced in:** `provenance.py:add_provenance` / `check_provenance` (every processed row carries `data_layer`).
  `validate.py:check_person_geography`: a person-layer table may not carry a FIPS, ZIP, ZCTA, geo_id or lat/lon
  column. `mcp/trace.py:is_person_level_id` and `trace_evidence` refuse participant ids (`PERSON_REFUSAL`).
  `ingestion/nhanes.py` joins only on SEQN. `database.py:build_db` views add no cross-layer joins. Person-level results
  enter recommendations only as group statistics, labelled by `tools.py:get_patient_phenotype_signature` (caveat),
  `agents/deployment_agent.py:PERSON_NOTE` and the `shared_context.note` of `rank_deployment_opportunities`.
* **Tested:** code search for merges (above). `measure-it validate`: 32/32 person-layer tables have no geographic
  column. `trace_evidence` on `nhanes:62161`, `participant:nhanes:62161` and `stanford_covid:A0NVTRV` was refused. The
  trace of a `signature:` id returns one group-level row.
* **Result:** pass. No violation found.

### 2. Geographic joins are ecological; no individual location is inferred

* **Enforced in:** the geographic tables hold places only (`geography/features.py`, `geography/query.py` caveat
  "Ecological, place-level data: nothing here describes an individual or where any participant lives").
  `facilities/matching.py` locates facilities, never patients. `agents/deployment_agent.py:PERSON_LOCATION_NOTE`
  (new) answers "where do participants live" premises.
* **Tested:** location-like column scan of the 32 person-layer tables. The NHANES, Stanford and MapMECFS DATA_AUDITs
  record no geography below the nation. Agent question 4 asked for participants "living in San Diego County".
* **Result:** pass. The agent now states the premise is impossible (problem 6).

### 3. Omics belong to a wearable participant only with shared identifiers; otherwise "condition-level molecular enrichment"

* **Enforced in:** `ingestion/mapmecfs.py` (linkage audit -> docs/MAP_MECFS_LINKAGE_AUDIT.md,
  `modality_overlap__mapmecfs`): 0 wearable x omics pairs have comparable ids, so MapMECFS omics are linked omics x
  omics only. `wearables/digital_person.py` writes `molecular_context_statement` ("NOT individual multi-omics").
  `tools.py:MOLECULAR_LABEL` and `get_molecular_context`, `scoring/recommend.py` and `mcp/server.py:GUARDRAILS_MD`
  rule 5 label all molecular output.
* **Tested:** every "multi-omics" string in code, docs, results, the API and dashboard output is a negation.
  `digital_phenotype_datasets` gives molecular context type `condition_level_molecular_enrichment_via_ontology` for
  NHANES and Stanford, and `omics_only_participant_linked` (no wearable) for MapMECFS. `SOURCE_REGISTRY.yaml` gives
  MapMECFS `wearable: false` and linkage "omics x omics only".
* **Result:** pass.

### 4. Every derived table carries provenance fields

* **Enforced in:** `store.py:write_table` refuses a table without `PROVENANCE_COLUMNS`.
  `validate.py:check_provenance_all` checks them.
* **Tested:** `measure-it validate`: 122/122 processed tables pass. The 2 boundary geoparquets are exempt as geometry
  files.
* **Result:** pass for every table in `data/processed` and the DuckDB.
* **Limitation:** 130 of the 158 CSV extracts in `results/tables/` at the time of this audit (131 of 159 with
  `number_audit.csv`, which this audit added; report tables such as `test7_summary.csv` and
  `demo_top10_regions.csv`) have no provenance columns. They are generated by the named modules from processed tables
  that do. Rows that recommendations cite carry an `object_id` that `trace_evidence` resolves. They are not canonical
  data tables, but read literally, the SPEC's "every derived table" includes them. Adding the columns means changing
  about 15 report writers, which this audit did not do.

### 5. Every score exposes its component variables

* **Enforced in:** `scoring/opportunity.py:combo_frame` (raw components, `*_pct`, one composite and rank per weight
  set, `member_components` for sets, Monte Carlo interval columns). `geography/features.py` (`diagnostic_desert` with
  `desert_pct_*`). `measurements/*` (`measurement_evidence_strength__*`, `technology_maturity__*`,
  `regulatory_visibility`, the five Phase 3 dimensions kept separate). `facilities/matching.py:find_candidate_clinics`
  (six characteristics ranked separately, never one score). `participant_adapter_scores.components_json`.
* **Tested:** schema scan of every processed table for score, composite, index, tier, strength, rank, pct,
  readiness, capacity, saturation and desert columns; each has its component columns.
* **Result:** pass.

### 6. No LLM-generated number is data; tools return `UNKNOWN / NOT AVAILABLE` instead of guessing

* **Enforced in:** `tools.py:_run` / `_env` / `_status` (a tool never raises and never fills in), and
  `resolve_condition_arg` / `_cond_or_unknown` (only a `matched` condition is used). `normalize_condition` (fixed,
  problem 2). `wearables/signatures.py:_resolve_condition` (fixed, problem 1). `tools.py:get_regulatory_context` /
  `_regulatory_for_bundle` (fixed, problem 4). `geography/query.py:resolve_geography` returns UNKNOWN for ambiguous
  names ("Washington County": 30 matches). The agent is deterministic with no language model and uses only `matched`
  results. `mcp/server.py:GUARDRAILS_MD` / `INSTRUCTIONS` tell calling models to add nothing. The only numbers typed
  into code are curated configuration (`configs/*.yaml`, evidence type `curated_config`) and published paper values
  kept as labelled reproduction targets (`ingestion/stanford_wearables.py`).
* **Tested:** the 708-call tool sweep, the agent questions and the dashboard sweep. A search of `configs/` for data-like
  numbers found none.
* **Result:** fixed (problems 1-5 and 8). After the fixes, every nonsense input returns UNKNOWN. The `ok` answers to
  ambiguous text come from `search_condition` (a ranked candidate list by design), from text matches that the output
  flags with `matched_by: "phenotype text"`, and from empty optional filters (`measurement: ""` means no filter).
* **Limitations:** `get_geographic_context` gives an exact state name precedence over same-named counties ("New York" ->
  the state, "Washington" -> the state), and a unique county base name resolves without "County" ("Los Angeles" -> Los
  Angeles County). The resolved name is always returned. `get_regulatory_context("blood test")` still lists fuzzy
  candidates (autonomic_testing, blood_biomarkers, microvascular_function, all score 0.8), now with the fuzzy-only
  caveat.

### 7. No silent proxy substitution; `burden_evidence_level` is always exposed; level D has no burden

* **Enforced in:** `geography/features.py:_def` (pre-specified primary measure per condition x level,
  docs/BURDEN_DEFINITIONS.md) and `build_burden` (every row has `burden_evidence_level`, `measure_role`,
  `proxy_for_condition`). `geography/query.py:get_condition_burden` returns UNKNOWN plus a reason for level D.
  `tools.py:get_condition_burden` adds the level legend to every answer. `scoring/opportunity.py` excludes level D,
  renormalises the burden weight and says so in `uncertainties`. A condition set takes the level of its least direct
  member, and proxies widen the Monte Carlo band (C x2). `geography/features.py:test4_series_legend` is new
  (problem 7).
* **Tested:** all 558,246 burden rows have a level, and 0 level-D rows carry a value. Across 26,328 level-D feature
  rows and 25,560 level-D opportunity rows, no burden value, count, proxy, desert burden term or burden percentile is
  filled. The tool sweep ran POTS, endometriosis and six other level-D conditions through burden and ranking: burden
  UNKNOWN, `burden_excluded_level_D = true`, and the uncertainty says "the ranking carries no burden information".
* **Result:** pass after problems 5 and 7.
* **Limitation (resolved after this audit):** `deployment_opportunities.parquet` was not rebuilt in this audit, so its
  stored `member_components` JSON still held the pre-fix `"nan"` strings for level-D set members (12,876 rows). The full
  offline run `20260924T150209Z` on the final code rebuilt it: 0 of 76,680 rows hold one now.

### 8. State estimates are never downscaled or presented as county prevalence

* **Enforced in:** `geography/features.py:_inherited_rows` (`derivation = inherited_from_state`,
  `source_geographic_resolution = state`, `inherited = true`, and the `measure_label` suffix "[STATE value inherited
  as county context; not county prevalence]"). `geography/query.py:CAVEAT_INHERITED`. `tools.py:get_condition_burden`
  and `rank_deployment_opportunities` caveats. `scoring/recommend.py:recommendation` (`burden_evidence_level_note`).
  The Monte Carlo doubles an inherited value's uncertainty. Dashboard map warnings and colour-bar titles, and the
  Figure 6/8 captions, label it too.
* **Tested:** all 3,144 county rows with state resolution are inherited and labelled. The tool, agent, API and
  dashboard outputs were read in context: every "county prevalence" mention is a negation.
* **Result:** pass.
* **Note:** a long-COVID county row keeps `burden_evidence_level = A`, because the level describes the state HPS
  measure. Every output pairs it with `inherited = true` and resolution `state`.

### 9. A regulatory record is a deployment-readiness signal; trial use is not evidence of effectiveness or diagnostic validity

* **Enforced in:** `measurements/regulatory.py:REGULATORY_NOTE` and `ingestion/openfda.py:REGULATORY_NOTE` (a caveat
  on every regulatory answer). `ingestion/clinicaltrials.py:GUARDRAIL_NOTE` (on every trial answer).
  `agents/deployment_agent.py` sentence templates: "research activity, not proof that it works"; "a deployment-readiness
  signal, not evidence that it detects the condition". Candidate measurements are ordered "by registered research
  activity, not by validity".
* **Tested:** a search of docs, results, code and outputs for "FDA-cleared for", "approved for", "validated for",
  "diagnostic accuracy", "detects <condition>" and similar found only negations.
* **Result:** pass.

### 10. Research readiness stays separate from patient burden

* **Enforced in:** separate `research_readiness` and `burden` components in `scoring/opportunity.py`. The
  `find_relevant_research_centers` caveat says the order is "a display order ..., not a score or a quality ranking".
* **Result:** pass.

### 11. Product language

* **Enforced in:** `validate.py:check_language` / `forbidden_hits` (8 forbidden terms; quoted or rule-list mentions
  are excused). It now also scans README.md, `results/figures/*.md`, `results/*.json` string values,
  `data/raw/*/DATA_AUDIT.md`, and the MCP and API text returned by `service_texts()`: server instructions, guardrail
  and data-layer resources, every tool's title and description, and parameter docs. Also
  `figures/schematics.py:banned_phrases_in`, the product-language tests in `test_tools.py`, `test_dashboard.py`,
  `test_scoring.py` and `test_facilities.py`, and the framing in `scoring/recommend.py:next_step_text` ("candidate
  deployment opportunity, not a validated diagnostic pathway") and `facilities/matching.py:FRAMING` ("may be a viable
  implementation or study partner").
* **Tested:** repo-wide search, 708 tool outputs, 5 agent answers, 17 API responses, 22 dashboard runs (text,
  tables, plotly specs). `measure-it validate` now scans 70 files (this one included) and 66 MCP/API texts: 0 unquoted uses.
* **Result:** pass. "Measure It to Cure It" is the product name; it matches no forbidden term ("cures").
* **Limitation:** the check is lexical. It catches the listed phrases, not every paraphrase. Wider variants were
  searched by hand once, in this audit.

---

## Re-running this audit

```bash
uv run measure-it validate                                   # provenance, person-layer geography, language (incl. MCP/API text)
uv run pytest -m "not network" tests/test_guardrail_audit.py # regression tests for the fixes above
uv run measure-it ask "Does the patient have POTS, and which clinic is best in Texas?"   # Scope section expected
```
