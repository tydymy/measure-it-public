# Agentic MCP layer

The engine's query functions are exposed three ways, all through one facade:

| layer | where | what it is |
|---|---|---|
| tool facade | `src/measure_it/tools.py` | one plain-Python function per SPEC tool (+3 extras); each wraps an existing query function and returns the same JSON envelope |
| MCP server | `src/measure_it/mcp/server.py` | FastMCP 4 server (stdio) exposing every facade function as a typed, read-only MCP tool, plus 2 prompts and 4 resources |
| trace | `src/measure_it/mcp/trace.py` | `trace_evidence(object_id)`: any cited object id -> rows, provenance, SOURCE_REGISTRY entry, DATA_AUDIT.md, raw files, one level of lineage |
| deterministic agent | `src/measure_it/agents/deployment_agent.py` | `answer_question(question)` / `uv run measure-it ask "..."`: LLM-free, templated, every sentence cites object ids |

Nothing in this layer computes a new result: every number comes from the wrapped query function or from a stored table.
The layer only resolves inputs with the project normalisers, caps long lists (and says so), and adds caveats and
provenance.

## 1. Connecting

### Claude Code

`.mcp.json` at the repository root declares the server (project scope):

```json
{
  "mcpServers": {
    "measure-it-public": {
      "type": "stdio",
      "command": "uv",
      "args": ["run", "--frozen", "measure-it", "serve-mcp"],
      "env": {}
    }
  }
}
```

Start Claude Code in the repository root, approve the project server when asked (or run `claude mcp list` to check it),
and the 18 tools appear as `mcp__measure-it-public__<tool>`. Two settings are worth knowing:

* `MEASURE_IT_OFFLINE=1` in `env` runs the server with the project's offline guard: every tool answers from
  `data/processed`, `data/raw` and `data/_http_cache` only, and `search_us_open_data` returns `UNKNOWN / NOT AVAILABLE`
  for a query that is not in the HTTP cache instead of calling catalog.data.gov (the cached queries are the 54 run at
  ingestion, listed in `data_gov_discovery_log` and in the UNKNOWN reason; since 2026-09-24 a query that differs from a
  cached one only in case or spacing, e.g. "Long COVID" for "long COVID", is answered from that cached response and
  says so in a caveat). `search_us_open_data` is the only tool that can reach the network.
* Claude Code truncates MCP tool output above `MAX_MCP_OUTPUT_TOKENS` (default 25,000 tokens). The MCP text
  content is compact JSON. The largest default calls are close to that limit (measured 2026-09-24; tokens are the cl100k
  count, a proxy - Claude's tokenizer usually gives more on JSON):
  `rank_deployment_opportunities` top 10 compact ~90 kB / ~24k tokens (~85 kB / ~22k before the 2026-09-24 fixups
  added candidate-site coordinates and `provenance.opportunity_row`), `get_phenotype_measurement_evidence` ~83 kB /
  ~22k, `get_measurement_evidence` (wearable bundle) ~77 kB / ~20k, `discover_candidate_measurements` ~56 kB / ~14k
  (re-measured after the fixups through the in-memory MCP client, text content);
  `get_condition_burden` (county, 60 rows) is ~65 kB. Every other tool stays below ~50 kB. For these calls, set
  `MAX_MCP_OUTPUT_TOKENS=50000` in the environment Claude Code runs in (not in the server's `env`), or ask for fewer
  items (`top_n=5`, `max_object_ids`, `max_rows`).

### Claude Desktop

Claude Desktop does not start servers in a project directory, so point `uv` at the checkout
(`claude_desktop_config.json`, macOS `~/Library/Application Support/Claude/`, Windows `%APPDATA%\Claude\`):

```json
{
  "mcpServers": {
    "measure-it-public": {
      "command": "uv",
      "args": ["run", "--frozen", "--directory", "/absolute/path/to/measure-it-public", "measure-it", "serve-mcp"],
      "env": {"MEASURE_IT_OFFLINE": "1"}
    }
  }
}
```

### Anything else

```bash
uv run measure-it serve-mcp                 # stdio server (same as python -m measure_it.mcp.server)
uv run measure-it serve-mcp --transport http --port 8765   # streamable HTTP at /mcp, plus GET /health (docs/DEPLOYMENT.md)
uv run python -m measure_it.tools --benchmark                  # facade latency table (section 5)
uv run python -m measure_it.tools trace_evidence '{"object_id": "geo:06073"}'   # call one tool from the shell
```

In Python, without a subprocess (this is how the tests talk to the server):

```python
import asyncio
from fastmcp import Client
from measure_it.mcp.server import build_server

async def main():
    async with Client(build_server()) as c:
        r = await c.call_tool("normalize_condition", {"condition_or_code": "ME/CFS"})
        print(r.structured_content["data"]["matches"][0]["object_id"])     # condition:me_cfs
asyncio.run(main())
```

The server is safe on stdio: the MCP SDK points file descriptor 1 at stderr while it serves, so progress lines that
wrapped modules print (e.g. `[scoring] spatial base: ...`) go to stderr and never reach the JSON-RPC stream (checked
with a real `uv run measure-it serve-mcp` subprocess and the FastMCP stdio client). The scoring module now writes that
line to stderr itself, and `python -m measure_it.tools <tool> '<json>'` and `measure-it ask --json` redirect progress
output to stderr, so their stdout is JSON only.

## 2. The envelope every tool returns

```jsonc
{
  "tool": "get_condition_burden",
  "status": "ok",                         // "ok" | "UNKNOWN / NOT AVAILABLE" | "error"
  "query": {...},                         // the arguments
  "data": {...},                          // the wrapped function's result (NaN / inf -> null)
  "reason": "...",                        // only when status is not "ok": why there is no answer
  "caveats": ["...", "Tool output only: do not add facts, numbers or places that are not in it; report UNKNOWN / NOT AVAILABLE as unknown."],
  "provenance": {"producer": "measure_it.geography.query.get_condition_burden",
                 "tables": ["geo_condition_burden", "..."],
                 "sources": [{"source_id", "name", "source_version", "retrieved_at", "data_audit"}],
                 "object_ids": ["geo:40109", "..."], "n_object_ids": 1,
                 "trace_with": "trace_evidence(object_id) resolves any of these ids to rows, sources and raw files"},
  "truncation": [{"path": "data.rows", "returned": 60, "total": 3144}],   // every list that was capped
  "meta": {"elapsed_s": 0.0002, "cache_hit": true, "data_stamp": 1790243272.87}
}
```

Rules the facade enforces:

* **UNKNOWN stays UNKNOWN.** A wrapped function that finds nothing (unknown condition, level-D burden, no stored row,
  no trial) gives `status: "UNKNOWN / NOT AVAILABLE"` with its own `reason`; values it could not compute stay the
  sentinel string. An exception becomes `status: "error"` with the exception text, never a guess.
* **One condition resolver.** Condition arguments go through `normalize_condition` (the ontology normaliser); only an
  unambiguous `matched` result is used. `"fatigue"` returns UNKNOWN with the candidate conditions instead of silently
  becoming ME/CFS. `rank_deployment_opportunities` additionally accepts condition sets ("Long COVID or ME/CFS",
  "demo cluster", any "A or B").
* **Truncation is always reported.** Tool-specific caps (`max_rows`, `top_n`, `max_trials`, compact ranking) and a size
  guard (lists capped at 50, 25 and then 10 items while the JSON is above 60 kB; primary result lists such as
  `recommendations`, `rows`, `candidates` are never cut by the guard) each add a `truncation` entry.
* **Caching.** Results are memoised per (tool, arguments, data stamp); the stamp is the newest mtime of
  `data/processed`, `configs`, `results/tables` and `SOURCE_REGISTRY.yaml` (re-read at most every 2 s), so rebuilding a
  table invalidates every cached answer. `search_us_open_data` is not memoised (the HTTP layer caches catalog responses
  on disk, and an "unreachable" answer must not outlive an outage).

## 3. Tools

| tool | wraps | inputs | output (`data`) | caveats carried |
|---|---|---|---|---|
| `search_condition` | `ontology.normalize.search_condition` | `condition`, `limit` | ranked `candidates` with MONDO/EFO/MeSH/ICD-10-CM ids, `score`, `match_reason`, `object_id` | acronyms match exactly only |
| `normalize_condition` | `ontology.normalize.normalize_condition` | `condition_or_code` (text, alias, ICD-10-CM, MONDO/EFO/MeSH/HPO id) | `match_status` (matched / ambiguous / partial), `matches`, `other_candidates`, `cms_fy2026` for ICD codes, `phenotype_axis` when the text is also a phenotype | ambiguous / partial are not answers (envelope status UNKNOWN, candidates kept in `data`); HPO ids give phenotype associations, not equivalence |
| `get_patient_phenotype_signature` | `wearables.signatures.get_patient_phenotype_signature` | `condition` (condition or phenotype text), `top_n` | per dataset x phenotype: definition, label basis, `is_proxy`, n cases / controls, FDR-significant `top_features` (effect, CI, q, `signature:` id), null results | group statistics from public person-level data; proxies are not the condition; never linked to places |
| `get_molecular_context` | `omics.query.get_molecular_context` + `omics.coherence.get_molecular_coherence` | `condition`, `top_n`, `include_broad` | `label` (condition-level molecular enrichment), evidence by type with counts by source and top entities (each Ensembl gene, rsID, GWAS Catalog study and GEO series cited as `gene:` / `variant:` / `gwas_study:` / `geo_series:`), coverage, `coherence` (source counts, gene agreement, supported systems with `post_hoc_flagged`) | not participant-linked, not patient multi-omics, not mechanism; no combined omics score |
| `discover_candidate_measurements` | `measurements.query.discover_candidate_measurements` | `condition`, `phenotype` (optional), `top_n` | `candidates` in the documented research-activity order, each with the five separate dimensions, text-mining precision, example trial / NIH / FDA ids | the order is research activity, not validity; trials and grants are not efficacy |
| `get_measurement_evidence` | `measurements.query.get_measurement_evidence` | `measurement` (class / bundle / alias), `condition`, `max_object_ids` | per class/bundle: counts (objective use, outcome measure, recruiting, NIH core projects), trial / grant / FDA ids, snippets, phenotype-signal results | absence of a text-mined mention is not evidence of non-use |
| `get_phenotype_measurement_evidence` (extra) | `measurements.query.get_phenotype_measurement_evidence` | `phenotype` (demo phenotype or one of its axes) | Phase 3 table: observable signals, technologies, evidence per source, five dimensions, gap note | no composite score |
| `get_regulatory_context` | `measurements.regulatory.get_regulatory_context` | `technology`, `max_records` | matched classes (method, `measurement:` id), per-class summary (product codes, 510(k) / De Novo / PMA counts), `regulatory_status` rows with `fda:code:` ids, matching FDA records with `fda:` ids | a regulatory record is a deployment-readiness signal, not evidence of diagnosis; questionnaires refused |
| `get_condition_burden` | `geography.query.get_condition_burden` | `condition`, `geography_level` (national / state / county), `geo_id`, `include_alternates`, `max_rows`, `order` | `burden_evidence_level`, primary measure, `rows` (value, CI, `inherited`, source resolution, source), `summary` over all rows | level D -> UNKNOWN with the reason; inherited state values are never county prevalence; level C is a proxy |
| `get_geographic_context` | `geography.query.get_geographic_context` | `geography_id` (FIPS, `geo:`, `zcta:`, name) | population, burden by condition (with evidence levels), SVI, ACS, PLACES, rurality, HPSA, providers by specialty group, HRSA sites, trials, NIH projects | ecological; provider taxonomies are self-reported |
| `find_candidate_clinics` | `facilities.matching.find_candidate_clinics` | `geography_id`, `condition`, `measurement`, `radius_km`, `max_sites` | pool sizes, `n_with_characteristic`, `candidates` (round-robin over six characteristics, each with value and separate rank, reasons, `facility:` / `npi:` / `trial:` / `nih:` ids), `absence_note` | never one combined score, never a quality ranking; locations are centroids |
| `find_relevant_trials` | `ingestion.clinicaltrials.trials_for_condition` | `condition` (or condition set), `measurement` (class / bundle / alias, optional), `literal_only` (default true), `statuses`, `us_sites_only`, `max_trials` | `summary` (n, by status, recruiting), `trials` with `trial:` ids and the measurement match fields | registration is a research-activity signal |
| `find_relevant_research_centers` | `facilities.query.find_relevant_research_centers` | `condition`, `measurement`, `top_n`, `include_expanded` | `centers` with trial and NIH counts, obligations, object ids | display order, not a score; obligations are not expenditure |
| `rank_deployment_opportunities` | `scoring.recommend.rank_deployment_opportunities` | `condition` (condition or set), `measurement`, `geography_level` (county / state), `top_n`, `weight_set`, `include_small_population`, `max_sites_per_region`, `detail` (compact / full) | SPEC deployment-opportunity schema per region; compact mode shows the condition-level context once in `shared_context` and 3 candidate sites per region | candidate deployment opportunity, not a validated pathway; burden evidence level per region; rank interval |
| `search_us_open_data` | `agents.datagov.search_us_open_data` | `query`, `rows` | title, agency, description, access level, resource URLs / formats, modified, `retrieve_from` | Data.gov hosts metadata only |
| `trace_evidence` | `mcp.trace.trace_evidence` | `object_id`, `max_rows` | rows, provenance columns, sources (registry entry, DATA_AUDIT.md, MANIFEST, raw record files), lineage | person-level ids refused |
| `get_measurable_biology` (extra) | `omics.graph.get_measurable_biology` | `condition`, `measurement`, `include_unsupported` | measurement classes with molecular support (`molbio:` ids), system-level status | condition -> system data-derived, system -> class curated |
| `list_sources` (extra) | `SOURCE_REGISTRY.yaml` | `data_layer` (optional) | every source: publisher, version, retrieval date, layer, status, licence, outputs, audit path | layers are never joined as the same people |

Every tool description in the MCP listing starts from the function's docstring and adds: *use only facts in the
returned data; never add facts beyond tool output; if status is 'UNKNOWN / NOT AVAILABLE', report it as unknown with
the reason; cite provenance.object_ids*. Every parameter has a description; enumerations (`geography_level`, `detail`,
`order`) are JSON-schema enums, so an invalid value is rejected before the tool runs.

Prompts and resources:

| kind | name | content |
|---|---|---|
| prompt | `measure_it_guardrails` | the rules below + the data-layer table |
| prompt | `deployment_question(question)` | the rules + the documented tool order for "what / where / who" questions |
| resource | `measure-it://guardrails` | rules: no facts beyond tool output; UNKNOWN stays UNKNOWN; cite object ids; four layers never joined as the same people; condition-level molecular enrichment; burden evidence levels; research and regulatory signals are not efficacy; product language |
| resource | `measure-it://data-layers` | person / condition_molecular / geographic / facility (+ ontology, measurement, derived) and how each may be joined |
| resource | `measure-it://object-ids` | the namespaces `trace_evidence` resolves |
| resource | `measure-it://tools` | tool -> wrapped function |

Product language in every generated text: candidate, deployment opportunity, measurement desert, research readiness,
evidence-supported, measurable phenotype, implementation candidate. The guardrails tell the calling model never to
write "proves", "diagnosed by AI", "patient has", "definitive biomarker", "best clinic", "optimal treatment",
"confirms mechanism" or "cures"; the tests scan every tool output with the `validate` language check.

**User-supplied datasets (bring your own data).** The server also registers two read-only tools tagged `byod`,
outside the 18 above (so `measure_it.tools.TOOLS`, the HTTP API tool list and `/health` `n_tools` are unchanged):
`list_user_datasets()` (user datasets on this machine: title, owner, status, comparator, condition, class, group sizes,
plan hash) and `get_user_dataset_metric(dataset_id)` (locked plan, evaluation models, the tier-2.5 performance record
with its caveats, the last deploy report). Both return aggregate results only, never a person row or a participant id;
with nothing ingested they return UNKNOWN / NOT AVAILABLE. See docs/BRING_YOUR_OWN_DATA.md.

## 4. trace_evidence

| namespace | primary table | native id | one level of lineage |
|---|---|---|---|
| `condition` | condition_registry | canonical id, condition-set id (`long_covid_or_me_cfs`) or `a+b` | ontology mappings by ontology, burden definition per level (A-D), Monarch phenotype annotations; sets -> member conditions |
| `measurement` | measurement_registry | class or bundle id | member classes / containing bundles, adapter rows, example FDA records, configs |
| `trial` | clinical_trials | NCT id | matched conditions (literal flag), registered sites (U.S. examples), resolved facilities; raw: the `records/batch_*.json` page that holds the record |
| `nih` | nih_projects | appl_id | conditions (match tier), awardee facility; raw: `records.jsonl.gz` |
| `npi` | providers | NPI | resolved facility; raw: the NPPES V.2 ZIP named in `source_version` |
| `hrsa_site` | facilities__hrsa | BPHC site number | resolved facility; raw: the HRSA sites CSV |
| `facility` | facilities (+ clinic_registry / research_site_registry presence) | facility_id | source records resolved into it (`npi:`, `hrsa_site:`, ClinicalTrials.gov site keys), its trials and NIH projects |
| `fda` | fda_510k / fda_pma / fda_named_device_findings; `code:` -> fda_device_classification + measurement_regulatory_status | K / DEN / P number, `code:<product code>` | product code, mapped measurement classes; raw: `api/510k/<code>*.json` |
| `gene`, `variant`, `gwas_study` | condition_molecular_evidence | Ensembl id, rsID, GCST | conditions; summary by source / evidence type; raw: files named by the accession |
| `geo_series` | geo_study_catalog | GSE | conditions whose curated query found it |
| `geo` | geographies + geo_context (`US` -> national burden rows; `zcta:` -> geo_context_zcta) | FIPS | one geo_condition_features row per condition (burden value, evidence level, resolution, inherited, burden_row_id, desert) |
| `opportunity` | deployment_opportunities | `<condition>\|<measurement>\|<geo_id>` | component values (raw + percentile, burden evidence level) -> geo_condition_features rows -> geo_condition_burden rows and their sources and raw files; geo / condition / measurement ids |
| `signature` | phenotype_signatures | `<dataset>\|<phenotype>\|<feature>` | dataset (digital_phenotype_datasets) and label basis (definition, proxy flag, population, n); raw: the XPT / ZIP files named in its version |
| `measurement_evidence` | condition_measurement_evidence | `<condition>\|<measurement>` | condition, measurement, example trial / NIH / FDA / molecular ids, signal result |
| `measurement_signal` | measurement_phenotype_signal | `<dataset>\|<label>\|<class>\|<scope>` | measurement class, per-feature `signature:` rows |
| `phenotype_evidence` | phase3_measurement_evidence | `<phenotype>\|<measurement>` | measurement, signal result, example ids |
| `molbio` | measurable_biology | `<condition>\|<system>\|<class>[\|<catalogue>]` | condition, measurement class, Reactome pathways |
| `reactome` | results/tables/test3_reactome_ora.csv | R-HSA id | conditions whose gene lists overlap the pathway |
| `participant` (and `person`, `patient`, `seqn`, `subject`, any case) and bare participant ids (`nhanes:62161`, `stanford_...:<id>`, `mapmecfs:<id>`) | - | - | **refused**: person-level rows are summarised only |

Every source in a trace carries its SOURCE_REGISTRY fields, `data_audit` (path + exists), the MANIFEST path, and
`record_files` (path, URL, bytes, sha256, retrieval time) when the record's raw file can be singled out: a file whose path
contains the record's accession or code, a file named in its `source_version`, the ClinicalTrials.gov batch that holds
the NCT id (indexed once per process), or every file of a source with at most three. Otherwise it says so and points to
the MANIFEST. Lookups use `data/processed/measure_it_public.duckdb` when it was built after the table's Parquet file
(point lookups ~10 ms) and the Parquet file otherwise; both give the same rows (`lookup_backend` says which).

## 5. Latency

`uv run python -m measure_it.tools --benchmark` (MEASURE_IT_OFFLINE=1; aarch64, 20 cores; 2026-09-24). *Cold* = first
call of that tool in a fresh process with the tool cache empty (shared structures such as the ontology index may already
be loaded by the tools listed above it); *warm* = the same call again (memoised). kB = JSON size of the whole envelope
(`json.dumps` default separators; re-measured after the 2026-09-24 fixups, which added candidate-site coordinates and
`provenance.opportunity_row`; the timings are from the earlier benchmark run).

| tool | demo input | cold (s) | warm (s) | kB |
|---|---|---|---|---|
| search_condition | "long covid" | 0.29 | <0.001 | 2 |
| normalize_condition | "G93.32" | <0.001 | <0.001 | 2 |
| get_patient_phenotype_signature | "ME/CFS" | 0.04 | <0.001 | 33 |
| get_molecular_context | "Long COVID" | 0.62 | <0.001 | 33 |
| get_measurable_biology | "ME/CFS" | 0.01 | <0.001 | 30 |
| discover_candidate_measurements | "Long COVID", phenotype "orthostatic intolerance" | 0.06 | <0.001 | 58 |
| get_measurement_evidence | "wearable autonomic monitoring" x "Long COVID" | 0.10 | <0.001 | 80 |
| get_phenotype_measurement_evidence | "orthostatic/autonomic dysfunction" | 0.01 | <0.001 | 87 |
| get_regulatory_context | "wearable ECG patch" | 0.13 | <0.001 | 34 |
| get_condition_burden | "Long COVID", county | 0.25 | <0.001 | 67 |
| get_geographic_context | "San Diego County, California" | 0.09 | <0.001 | 48 |
| find_candidate_clinics | 06073, Long COVID, wearable autonomic monitoring | 0.47 | 0.001 | 47 |
| find_relevant_trials | Long COVID, wearable autonomic monitoring | 1.49 | <0.001 | 30 |
| find_relevant_research_centers | "ME/CFS" | 0.26 | <0.001 | 25 |
| rank_deployment_opportunities | Long COVID or ME/CFS x wearable autonomic monitoring, county, top 10 | 4.12 | 0.001 | 94 |
| search_us_open_data | "long COVID" (HTTP cache) | 0.003 | 0.002 | 28 |
| trace_evidence | the top opportunity id | 0.12 | <0.001 | 36 |
| list_sources | - | 0.001 | <0.001 | 50 |

The ranking's cold time is the spatial base of the scoring engine (1.5 s: 326,685 facilities, 1.07 M individual NPIs)
plus candidate sites for ten regions; over stdio the first call took 5.4 s and the repeat 0.007 s. Caching added in
this layer: memoised tool results (above), the DuckDB point-lookup connection, per-table schema / MANIFEST / raw-file
listings and the ClinicalTrials.gov batch index (lru caches keyed by file mtimes); the wrapped modules keep their own
table caches. A `trace_evidence` of every one of the 104 object ids in a compact top-10 ranking takes ~1.8 s in total (all resolve).

## 6. Real example calls (captured 2026-09-24 through the in-memory MCP client, abbreviated)

The outputs below are real; `...` marks cuts. The scoring tables were being rebuilt by another contributor on the
same day, so ranks can differ from `results/SCORING_RESULTS.md` at other times; the calls regenerate them.

### normalize_condition

```jsonc
// normalize_condition({"condition_or_code": "G90.A"})
{"tool": "normalize_condition", "status": "ok",
 "data": {"match_status": "matched", "input_type": "icd10cm",
          "matches": [{"canonical_condition_id": "pots", "preferred_name": "Postural orthostatic tachycardia syndrome",
                       "primary_mondo_id": "MONDO:0001315", "matched_id": "ICD10CM:G90.A", "predicate": "exact", "confidence": 1.0,
                       "match_reason": "exact ICD-10-CM code G90.A (CMS ICD-10-CM FY2026 description match (curated); predicate exact)",
                       "object_id": "condition:pots"}],
          "cms_fy2026": {"code": "G90.A", "description": "Postural orthostatic tachycardia syndrome [POTS]", "billable": true}, ...},
 "provenance": {"object_ids": ["condition:pots"], ...}, ...}
```

### get_condition_burden: a level-D condition and an inherited state value

```jsonc
// get_condition_burden({"condition": "POTS", "geography_level": "county"})
{"status": "UNKNOWN / NOT AVAILABLE",
 "reason": "No public population measure: no CCW algorithm contains G90.A; PLACES has no orthostatic item.",
 "data": {"condition_id": "pots", "burden_evidence_level": "D", "value": "UNKNOWN / NOT AVAILABLE", "rows": [], ...},
 "caveats": ["Ecological, place-level data: nothing here describes an individual or where any participant lives.",
             "burden_evidence_level D: A = direct condition measure, B = closely matching coded condition, C = symptom/comorbidity proxy, D = no usable burden estimate. Always state it.", ...]}

// get_condition_burden({"condition": "Long COVID", "geography_level": "county", "geo_id": "40109"})
{"status": "ok",
 "data": {"burden_evidence_level": "A", "primary_measure_id": "lc_current_pct_all_adults", "primary_source_resolution": "state",
          "rows": [{"geo_id": "40109", "geo_name": "Oklahoma County, Oklahoma", "value": 7.5, "value_unit": "percent",
                    "ci_low": 5.3, "ci_high": 10.2, "period": "Aug 20 - Sep 16, 2024", "burden_evidence_level": "A",
                    "source_geographic_resolution": "state", "inherited": true,
                    "measure_label": "Currently experiencing long COVID, as a percentage of all adults [STATE value inherited as county context; not county prevalence]",
                    "burden_row_id": "long_covid|inherited_from_state|county:40109|gsea-w83j:lc_current_pct_all_adults:Oklahoma:By State:Oklahoma:72",
                    "object_id": "geo:40109", ...}],
          "summary": {"n_rows": 1, "n_inherited_state_values": 1, "evidence_levels": ["A"], "source_resolutions": ["state"], ...}},
 "caveats": ["...", "Rows with inherited = true are STATE estimates attached to counties (source resolution 'state'); never describe them as county prevalence.", ...]}
```

### SPEC demo query: rank_deployment_opportunities

The Deployment Recommendation page's query ("Find 10 U.S. regions where wearable monitoring of autonomic/activity
abnormalities in Long COVID or ME/CFS would be useful to evaluate, and identify clinics/research sites that could
plausibly participate"):

```jsonc
// rank_deployment_opportunities({"condition": "Long COVID or ME/CFS", "measurement": "wearable autonomic monitoring",
//                                "geography_level": "county", "top_n": 10})           5.2 s cold, 90 kB
{"status": "ok",
 "data": {"condition_resolution": {"condition_id": "long_covid_or_me_cfs", "match_reason": "condition-set alias"},
          "weights": {"burden": 0.2, "vulnerability": 0.2, "diagnostic_desert": 0.2, "clinic_capacity": 0.2, "research_readiness": 0.2},
          "n_regions_ranked": 2395, "ranking_basis": "deployment_opportunities (precomputed)", "opportunity_rows_stored": true,
          "ranking_universe": "counties of the 50 states + DC with population >= 10,000",
          "shared_context": {"phenotype": [...], "molecular_context": [...], "technology": {"regulatory_context": [...], "technology_evidence": [...]}, ...},
          "recommendations": [
            {"condition": {"see": "data.shared_context.condition"}, "phenotype": {"see": "data.shared_context.phenotype"}, ...,
             "geography": {"name": "Oklahoma County, Oklahoma", "fips": "40109", "object_id": "geo:40109", "population": 806199.0,
                           "burden": {"percentile": 0.9713, "inherited": true,
                                      "members": {"long_covid": {"burden_value": 7.5, "burden_evidence_level": "A", "burden_inherited": true, "burden_source_resolution": "state"},
                                                  "me_cfs": {"burden_value": 28.0, "burden_evidence_level": "C", "burden_inherited": false, "burden_source_resolution": "county"}}},
                           "burden_evidence_level": "C",
                           "burden_evidence_level_note": "the set's evidence level is its least direct contributing member; the long-COVID member's burden is the state estimate inherited by the county (source resolution: state; ...)",
                           "vulnerability": {"svi_overall": 0.8431, "percentile": 0.8432}, "diagnostic_desert": {"index": 0.5243, "percentile": 0.5687},
                           "clinic_capacity": {"index": 0.8484, "percentile": 0.9653},
                           "research_readiness": {"index": 0.9376, "percentile": 0.9585, "condition_trials_in_pool": 4.0, "condition_nih_core_projects_in_pool": 1.0, "technology_experience_trials_in_pool": 3.0},
                           "composite": 0.8614, "rank": 1.0, "rank_interval_5_95": [1.0, 306.0], "p_top10_monte_carlo": 0.492, ...},
             "candidate_sites": [{"object_id": "facility:npi-1538213764", "facility_name": "MCBRIDE CLINIC ORTHOPEDIC HOSPITAL, LLC", "lat": 35.578131, "lon": -97.517635,
                                  "geocode_precision": "zcta_centroid", "distance_km": 10.13,
                                  "selected_via": {"long_covid": "clinical_specialty_match", "me_cfs": "clinical_specialty_match"}},
                                 {"object_id": "facility:ctgov-74f3a639213f", "facility_name": "University of Oklahoma Health Science Center - Oklahoma Clinical and Translational Science Institute - Appendix A & B",
                                  "selected_via": {"long_covid": "relevant_trial_history"}, ...},
                                 {"object_id": "facility:npi-1508144411", "facility_name": "University of Oklahoma Health Sciences Center",
                                  "selected_via": {"long_covid": "NIH_research_activity", "me_cfs": "technology_experience"}, ...}],
             "research_evidence": {"condition_trials": {"n": 4, "object_ids": ["trial:NCT05172024", "trial:NCT05524532"], "object_ids_note": "2 of 4 shown"},
                                   "technology_experience_trials": {"n": 3, ...}, "nih_projects": {"n_core_projects": 1, "object_ids": ["nih:11451891"], ...}},
             "uncertainties": ["long_covid: burden is the state HPS estimate for OK inherited by every county of the state (source resolution: state). It carries no within-state variation and is not county prevalence; ...",
                               "me_cfs: burden is a level-C proxy (CMS MMD Fibromyalgia, Chronic Pain and Fatigue prevalence, unsmoothed actual, all ages); not the condition's prevalence; uncertainty band inflated x2.0.",
                               "Rank interval (5th-95th percentile over 1000 draws of burden and weights): 1-306.", ...],
             "provenance": {"object_ids": ["opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|40109", "geo:40109", ...], "n_object_ids": 44,
                            "opportunity_row": {"object_id": "opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|40109", "stored": true, ...}},
             "recommended_next_step": "Candidate deployment opportunity (rank 1 of 2395 under the 'equal' weights; 5th-95th percentile rank interval 1-306 over Monte Carlo draws of burden and weights). A pilot evaluation of wearable autonomic / activity monitoring for Long COVID or ME/CFS in Oklahoma County, Oklahoma could test feasibility ... This is a candidate deployment opportunity, not a validated diagnostic pathway; check the listed uncertainties first."},
            ...9 more],
          "excluded_incomplete_burden": [...5 of 11 shown]},
 "truncation": [{"path": "data.recommendations[0].candidate_sites", "returned": 3, "total": 20}, ...],
 "provenance": {"n_object_ids": 104, "sources": [{"source_id": "cdc_long_covid", ...}, {"source_id": "cms_mmd", ...}, ...]}}
```

The ten regions of that call:

| rank | region | composite | rank interval 5-95 | P(top 10) | burden level | burden pct | SVI pct | desert pct | capacity pct | readiness pct |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Oklahoma County, OK | 0.861 | 1-306 | 0.49 | C | 0.97 | 0.84 | 0.57 | 0.97 | 0.96 |
| 2 | Wayne County, MI | 0.838 | 1-408 | 0.37 | C | 0.95 | 0.86 | 0.42 | 0.99 | 0.97 |
| 3 | Franklin County, AL | 0.834 | 4-247 | 0.23 | C | 0.96 | 0.76 | 0.88 | 0.76 | 0.81 |
| 4 | George County, MS | 0.830 | 2-288 | 0.30 | C | 0.99 | 0.84 | 0.97 | 0.50 | 0.85 |
| 5 | Crittenden County, AR | 0.827 | 1-581 | 0.30 | C | 0.62 | 0.92 | 0.83 | 0.87 | 0.90 |
| 6 | Bossier Parish, LA | 0.826 | 1-202 | 0.42 | C | 0.76 | 0.78 | 0.86 | 0.87 | 0.86 |
| 7 | Webster Parish, LA | 0.825 | 1-302 | 0.45 | C | 0.69 | 0.88 | 0.76 | 0.93 | 0.86 |
| 8 | Scotland County, NC | 0.823 | 1-289 | 0.26 | C | 0.92 | 1.00 | 0.85 | 0.95 | 0.40 |
| 9 | Forrest County, MS | 0.818 | 4-323 | 0.21 | C | 0.99 | 0.95 | 0.86 | 0.89 | 0.40 |
| 10 | Lauderdale County, MS | 0.815 | 5-315 | 0.18 | C | 1.00 | 0.96 | 0.84 | 0.89 | 0.40 |

Every region's burden evidence level is C (the least direct member: the ME/CFS burden is a CMS symptom/comorbidity
proxy) and every long-COVID burden is the state value inherited by the county; the wide rank intervals are the Monte
Carlo over burden uncertainty and weights.

### trace_evidence: from the recommendation back to raw source files

**Step 1 - the opportunity.** `trace_evidence({"object_id": "opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|40109"})`

```jsonc
{"status": "ok",
 "data": {"namespace": "opportunity", "tables": ["deployment_opportunities"], "rows": [{...the full 135-column row...}],
          "components": {"burden": {"raw": null, "percentile": 0.9713, "evidence_level": "C", "source_resolution": "county|state", "inherited": true,
                                    "members_contributing": "long_covid|me_cfs", "incomplete": false, "members_missing": [],
                                    "raw_note": "a condition set has no single raw burden value: its burden is the percentile of the mean of member burden percentiles; the member values are the geo_condition_features / geo_condition_burden rows in lineage"},
                         "vulnerability": {"raw": 0.8431, "percentile": 0.8432}, "diagnostic_desert": {"raw": 0.5243, "percentile": 0.5687},
                         "clinic_capacity": {"raw": 0.8484, "percentile": 0.9653}, "research_readiness": {"raw": 0.9376, "percentile": 0.9585},
                         "technology_saturation": {"raw": 0.1421, "percentile": 0.9556, "note": "reported only; not in the default formula"}},
          "composite_equal": 0.8614, "rank_equal": 1.0, "burden_incomplete": false, "burden_incomplete_reason": null, "burden_members_missing": [],
          "lineage": [{"object_id": "geo:40109", "relation": "geography"},
                      {"object_id": "measurement:wearable_autonomic_activity_monitoring", "relation": "measurement"},
                      {"object_id": "condition:long_covid", "relation": "member condition"}, {"object_id": "condition:me_cfs", "relation": "member condition"},
                      {"relation": "burden / vulnerability / desert inputs (geo_condition_features)", "table": "geo_condition_features",
                       "row": {"feature_id": "long_covid|40109", "burden_value": 7.5, "burden_evidence_level": "A", "burden_source_resolution": "state", "burden_inherited": true,
                               "burden_row_id": "long_covid|inherited_from_state|county:40109|gsea-w83j:lc_current_pct_all_adults:Oklahoma:By State:Oklahoma:72",
                               "svi_overall": 0.8431, "diagnostic_desert": 0.4608, ...}},
                      {"relation": "burden / vulnerability / desert inputs (geo_condition_features)", "table": "geo_condition_features",
                       "row": {"feature_id": "me_cfs|40109", "burden_value": 28.0, "burden_evidence_level": "C", "burden_source_resolution": "county",
                               "burden_row_id": "me_cfs|as_published|county:40109|mmd:f:2022:v:51:county:40109:unsmoothed_actual:age=all", "diagnostic_desert": 0.5877, ...}},
                      {"relation": "burden row behind the burden component (geo_condition_burden)", "table": "geo_condition_burden", "source_ids": ["cms_mmd"],
                       "row": {"measure_label": "CMS MMD Fibromyalgia, Chronic Pain and Fatigue prevalence, unsmoothed actual, all ages", "value": 28.0,
                               "value_unit": "percent of Medicare FFS beneficiaries (integer as published)", "burden_evidence_level": "C",
                               "source_table": "geo_condition_burden__cms_mmd", "source_record_id": "mmd:f:2022:v:51:county:40109:unsmoothed_actual:age=all", ...}},
                      {"relation": "burden row behind the burden component (geo_condition_burden)", "table": "geo_condition_burden", "source_ids": ["cdc_long_covid"],
                       "row": {"measure_label": "Currently experiencing long COVID, as a percentage of all adults [STATE value inherited as county context; not county prevalence]",
                               "value": 7.5, "inherited": true, "source_geographic_resolution": "state", "source_table": "geo_condition_burden__cdc_long_covid", ...}},
                      {"relation": "clinic_capacity / research_readiness inputs", "tables": ["providers (nppes)", "facilities", "facility_trials (clinicaltrials_gov)", "facility_nih_projects (nih_reporter)"], ...}],
          "sources": [
            {"source_id": "cms_mmd", "name": "CMS Mapping Medicare Disparities (MMD) Tool - Population View", "retrieved_at": "2026-09-23T20:39:38+00:00",
             "data_audit": [{"path": "data/raw/cms_mmd/DATA_AUDIT.md", "exists": true}],
             "raw": {"manifest": "data/raw/cms_mmd/MANIFEST.json", "n_manifest_files": 92,
                     "record_files": [{"path": "data/raw/cms_mmd/api/f_2022_v__prev_final_long_fltr12_racecat_all_sexcat_all_22_f.json",
                                       "url": "https://data.cms.gov/data-api/v1/mmd-tool/?_source=prev_final_long_fltr12_racecat_all_sexcat_all_22_f&year=22&... ...[277 characters]",
                                       "bytes": 6952186, "sha256": "e215a50e7b7292e0322a2cc8f8ea29ba24eb3356a43c747a3d30883b537031f8",
                                       "retrieved_at": "2026-09-23T20:45:04+00:00", "in_manifest": true}]}},
            {"source_id": "cdc_long_covid", "data_audit": [{"path": "data/raw/cdc_long_covid/DATA_AUDIT.md", "exists": true}],
             "raw": {"record_files": [{"path": "data/raw/cdc_long_covid/post_covid_conditions_gsea-w83j.csv",
                                       "url": "https://data.cdc.gov/api/views/gsea-w83j/rows.csv?accessType=DOWNLOAD", "bytes": 3596485,
                                       "sha256": "3356fa06bda30ef6bd01bb41f4ae349653ec6186e178b782fc745ef54811b923", ...}, ...]}},
            {"source_id": "cdc_svi", "raw": {"record_files": [{"path": "data/raw/cdc_svi/SVI_2022_US_COUNTY.csv", "bytes": 2301480, "sha256": "bc47d244153e...", ...}]}, ...},
            {"source_id": "nppes", "raw": {"record_files": [{"path": "data/raw/nppes/NPPES_Data_Dissemination_September_2026_V2.zip", "bytes": 1159503721, ...}, ...]}, ...},
            {"source_id": "hrsa_health_centers", ...}, {"source_id": "clinicaltrials_gov", ...}, {"source_id": "nih_reporter", ...}]}}
```

For a region whose burden is incomplete (a member with a defined burden measure has no value there, e.g. ME/CFS in
the Connecticut planning regions, `opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|09120`) the
trace carries `burden_incomplete: true`, `burden_incomplete_reason`, `burden_members_missing: ["me_cfs"]`, the
source's reason in `components.burden.missing_source_reason`, null `composite_equal` / `rank_equal`, and a
`composite_rank_note` (also in `notes`) that says why they are null; a complete county below the population minimum
gets a note that it is scored but not ranked, with its `rank_equal_incl_small`.

The burden component therefore rests on two raw files, each with its URL, byte count and sha256 from the source's
MANIFEST.json: the CMS MMD prevalence slice (level C proxy, county) and the CDC HPS long-COVID CSV (level A, but a
**state** value inherited by the county).

**Step 2 - a candidate site.** `trace_evidence({"object_id": "facility:npi-1538213764"})` gives the resolved facility
row and its source records (`npi:1538213764`, `npi:1659413649`, `npi:1932145505`, matched `nppes_same_name_same_zip`,
confidence 0.95). `trace_evidence({"object_id": "npi:1538213764"})` then gives the NPPES row:

```jsonc
{"status": "ok",
 "data": {"rows": [{"npi": "1538213764", "entity_type_label": "organization", "name": "MCBRIDE CLINIC ORTHOPEDIC HOSPITAL, LLC",
                    "primary_taxonomy_display_name": "Rehabilitation Hospital Unit", "city": "OKLAHOMA CITY", "state": "OK", "zip5": "73114",
                    "county_fips": "40109", "geocode_method": "zip5_as_zcta_internal_point",
                    "source_version": "NPPES_Data_Dissemination_September_2026_V2.zip (npidata_pfile_20050523-20260913.csv); NUCC 26.1", ...}],
          "lineage": [{"object_id": "facility:npi-1538213764", "relation": "resolved facility this NPI record was assigned to", "match_method": "nppes_same_name_same_zip", "match_confidence": 0.95}],
          "sources": [{"source_id": "nppes", "raw": {"record_files": [{"path": "data/raw/nppes/NPPES_Data_Dissemination_September_2026_V2.zip",
                                                                     "url": "https://download.cms.gov/nppes/NPPES_Data_Dissemination_September_2026_V2.zip",
                                                                     "bytes": 1159503721, "sha256": "6573e0d148c18a9a9d974248d77447744bd9666aa288c38d13ebf3d766e3cec4", ...}]}},
                      {"source_id": "nucc_taxonomy", "raw": {"record_files": [{"path": "data/raw/nucc_taxonomy/nucc_taxonomy_261.csv", ...}]}}]},
 "caveats": ["An NPPES taxonomy is self-reported; it does not mean the provider evaluates or treats any invisible illness.", ...]}
```

**Step 3 - research evidence.** `trace_evidence({"object_id": "trial:NCT05172024"})` gives the registry row ("NIH
RECOVER: A Multi-site Observational Study of Post-Acute Sequelae of SARS-CoV-2 Infection in Adults", COMPLETED), its
literal long-COVID match (terms "long COVID", "long haul COVID", "post-acute sequelae of SARS-CoV-2"), 86 sites (85 U.S.),
the resolved facilities, and the raw page that holds the record:
`data/raw/clinicaltrials_gov/records/batch_8e1ddc41df070933.json` (2,697,497 bytes, sha256 `b077034c0f16...`).

**Step 4 - the person-level evidence behind the phenotype context.**
`trace_evidence({"object_id": "signature:nhanes|mecfs_like_proxy|interdaily_stability"})` gives the aggregate row
(effect -0.32 SD, 95% CI -0.50 to -0.14), the label basis (`constructed_proxy`, `is_proxy: true`, 132 cases / 7,664
controls, the definition), the dataset row and the NHANES files named in its version (`DEMO_G.xpt`, `DEMO_H.xpt`,
`PAXDAY_G.xpt`, `PAXDAY_H.xpt`, `PAXHD_G.xpt`, ...). An individual participant cannot be traced:

```jsonc
// trace_evidence({"object_id": "participant:nhanes:62161"})
{"status": "UNKNOWN / NOT AVAILABLE",
 "data": {"refused": true, "namespace": "participant",
          "reason": "refused: person-level rows (participants of NHANES, the Stanford wearable studies or MapMECFS) are summarised only. This public tool never returns an individual's record; trace the group-level statistics instead (signature:, measurement_signal:, phenotype_evidence: object ids)."},
 "caveats": ["Person-level rows are summarised only; they are never joined to places or facilities and never exposed individually.", ...]}
```

### search_us_open_data

```jsonc
// search_us_open_data({"query": "long COVID"})       (answered from the HTTP cache)
{"status": "ok",
 "data": {"n_results": 20,
          "results": [{"rank": 1, "title": "Long-term Care and COVID-19", "agency": "Centers for Disease Control and Prevention",
                       "parent_organization": "U.S. Department of Health & Human Services", "access_level": "public",
                       "resource_urls": ["https://data.cdc.gov/api/v3/views/3j26-kg6d/export.csv?accessType=DOWNLOAD", ...],
                       "resource_formats": ["application/json", "application/xml", "text/csv"], "modified": "2026-09-10",
                       "retrieve_from": "https://www.cdc.gov/nchs/covid19/npals.htm", ...}, ...]},
 "caveats": ["Data.gov is a metadata catalog: it does not host the records; retrieve_from names the agency URL that does. A catalog hit is not evidence about a condition.",
             "Results are in Data.gov's own relevance order; this engine does not evaluate catalog relevance (the top hit can be unrelated to the query, e.g. 'Long-term Care and COVID-19' for 'long COVID').", ...]}
```

The top hit here is not about long COVID: catalog relevance is Data.gov's own and is not evaluated by this engine.

## 7. The deterministic agent: `uv run measure-it ask "..."`

`answer_question(question)` uses no language model. It

1. parses the question: condition spans against the ontology lexicon (preferred names, aliases, search terms, Mondo
   labels; acronyms such as POTS or ME/CFS case-sensitively; longest non-overlapping spans), each confirmed by
   `normalize_condition` (only a `matched` result is used); phenotype spans against the phenotype axes, their curated
   aliases (configs/relevance.yaml `phenotype_axis_aliases`, e.g. "PEM" -> post_exertional_malaise) and the Phase 3
   demo phenotypes; measurement spans against the vocabulary of the shared measurement resolver
   (`measure_it.measurements.resolve`: class ids / names / aliases, bundle ids / labels / aliases; a phrase that is
   exactly a bundle alias, e.g. "autonomic testing", gives the bundle); a county / state / FIPS geography; the level (`state` when the
   question asks for states) and the number of regions ("10 regions", default 5). All the parsed conditions together
   are one condition set (a named set when the members match, e.g. the demo cluster; otherwise "A or B");
2. calls the tools in this order: normalize_condition (search_condition when nothing matched) -> 
   get_phenotype_measurement_evidence (named phenotype) -> get_patient_phenotype_signature -> 
   discover_candidate_measurements -> get_measurement_evidence + get_regulatory_context -> get_molecular_context -> 
   get_condition_burden -> rank_deployment_opportunities -> (named geography: get_geographic_context + trace of its
   opportunity id) find_candidate_clinics for the named or top region -> find_relevant_research_centers + 
   find_relevant_trials -> trace_evidence of the top opportunity. The ranking is always national: for a named state
   at county level the agent asks for the national top 50, labels the list NATIONAL and lists that state's counties
   among the 50 separately (or says there are none; a within-state ranking is UNKNOWN because the tool has no state
   filter); a named state's stored row is reported as a rank among states, and its facility pool is the whole state;
3. fills fixed sentence templates from the tool output. Every factual sentence ends with the object ids it rests on
   (`[condition:pots]`, `[opportunity:...; geo:...]`, `[facility:...]`); UNKNOWN parts stay UNKNOWN and are listed; if no
   condition resolves it stops before ranking; opportunity ids are cited only for stored rankings (an on-the-fly ranking
   has no stored row and is cited by its geography, condition and measurement ids). `--json` prints the parse, every
   tool call (arguments, status, seconds) and the cited ids, and nothing else on stdout.

Template changes of the 2026-09-24 report review: the molecular sentence lists a condition's enriched systems in three
tiers (non-literature support after removing post-hoc-flagged sources / carried only by post-hoc-flagged sources,
mostly drug targets of trialled drugs / literature co-mention only) instead of "supported physiological systems";
person-level signatures from molecular datasets are called "molecular (not wearable) features", wearable datasets are
listed first, and a signature whose analysis flags a confound says so; the matcher's absence note names its condition
and measurement; the literal-match trial count explains how it differs from the "describe objective use" count.

Runs on 2026-09-24, re-run after the agent fixes of the report review (MEASURE_IT_OFFLINE=1; all 271 citations -
the sum over the six answers of each answer's distinct object ids; 164 distinct ids overall - resolve with
`trace_evidence`):

| question | parsed as | result |
|---|---|---|
| SPEC demo: "Where in the United States would deploying wearable-based autonomic or activity monitoring be most valuable for detecting objective abnormalities associated with Long COVID, ME/CFS, dysautonomia/POTS-like phenotypes, or related invisible illnesses, and which clinics or research sites are best positioned to implement the technology?" | 4 conditions = the demo cluster; wearable_autonomic_activity_monitoring; U.S., county | 42 tool calls, 88 cited ids; dysautonomia / POTS signatures and burden UNKNOWN (level D) |
| "Given phenotype orthostatic intolerance and candidate measurement technology nailfold capillaroscopy, where should we deploy it?" | phenotype axis orthostatic_intolerance; normalize_condition maps the text to dysautonomia (alias; other candidate POTS); nailfold_capillaroscopy bundle | 14 calls, 37 cited ids; capillaroscopy is not linked to the phenotype (UNKNOWN), burden level D |
| "Where should we deploy wearable monitoring for Zorblaxian drift syndrome?" (unknown condition) | no condition | stops after search_condition: condition UNKNOWN, nothing ranked |
| "Where should we deploy autonomic testing for POTS?" (level-D burden) | pots; autonomic_function_testing bundle | burden UNKNOWN (level D) in every region; the ranking carries no burden information; 14 calls, 36 cited ids |
| "Could we deploy wearable monitoring for Long COVID in San Diego County, California?" (named county) | long_covid; wearable bundle; geo:06073 | the county's stored row: rank 1,253 among counties in the county-level equal-weight ranking, burden level A but a state value inherited by the county; candidate facilities around it; 17 calls, 44 cited ids |
| "Where should we deploy wearable autonomic monitoring for Long COVID or ME/CFS in Texas?" (named state) | long_covid_or_me_cfs; wearable bundle; geo:48 | national county list labelled NATIONAL; Texas counties among the national top 50: Dallas County (national rank 37); Texas's stored state-level row: rank 4 among states; facility pool = the whole state; 26 calls, 65 cited ids |

SPEC demo question (abbreviated):

```markdown
## How the question was read
- Conditions: 'Long COVID' -> Long COVID (post-COVID-19 condition) (exact match to preferred_name 'Long COVID'); 'ME/CFS' -> Myalgic encephalomyelitis/chronic fatigue syndrome (exact match to alias 'ME/CFS'); 'dysautonomia' -> Dysautonomia (autonomic nervous system disorder) (...); 'POTS' -> Postural orthostatic tachycardia syndrome (exact match to alias 'POTS') [condition:long_covid; condition:me_cfs; condition:dysautonomia; condition:pots]
- Measurement: 'wearable-based autonomic or activity monitoring' -> wearable_autonomic_activity_monitoring (exact bundle id / label / alias 'wearable autonomic or activity monitoring' (bundles are preferred for deployment)) [measurement:wearable_autonomic_activity_monitoring]

## 1. What to measure
- Measurable phenotype in public data for Myalgic encephalomyelitis/chronic fatigue syndrome - nhanes: phenotype 'Unexplained fatigue with functional limitation (ME/CFS-LIKE PROXY; not ME/CFS)' (a PROXY label, not the condition itself), 132 cases vs 7,664 non-cases; 14 of 26 wearable features passed BH-FDR; e.g. sleep_proxy_midpoint_sd_h 0.33 (age_sex_adjusted_standardized_difference, 95% CI 0.15 to 0.51); interdaily_stability -0.32 (age_sex_adjusted_standardized_difference, 95% CI -0.50 to -0.14). [signature:nhanes|mecfs_like_proxy|sleep_proxy_midpoint_sd_h; signature:nhanes|mecfs_like_proxy|interdaily_stability]
- Public person-level wearable signature for Postural orthostatic tachycardia syndrome: UNKNOWN / NOT AVAILABLE - no wearable phenotype signature ... is labelled for ['pots']; ... [condition:pots]
- Person-level results are group statistics from public cohorts; those participants are not the people of any region, facility or trial named below.
- Wearable autonomic / activity monitoring for Long COVID (post-COVID-19 condition): 85 registered trial(s) describe objective use (73 as an outcome measure) and 3 NIH core project(s) mention it; evidence tier high (research activity, not proof that it works). [measurement_evidence:long_covid|wearable_autonomic_activity_monitoring; trial:NCT04768257; trial:NCT04806620]
- FDA context for ecg_ambulatory: decisions_on_record; 10 mapped product code(s), 1,649 510(k), 1 De Novo, 0 PMA decisions (a deployment-readiness signal, not evidence that it detects the condition). [measurement:ecg_ambulatory; fda:code:MLO; fda:code:MWJ]
- Condition-level molecular enrichment for Myalgic encephalomyelitis/chronic fatigue syndrome (public databases; not patient multi-omics, not mechanism): expression_study_metadata 30 records, genetic_association 38 records, known_drug 43 records, metadata_catalog 8 records; systems with non-literature support after removing post-hoc-flagged sources: none; carried only by post-hoc-flagged sources (mostly drug targets of trialled drugs, or one GWAS locus; not disease biology): autonomic_cardiac, mitochondrial_metabolic, neuronal; literature co-mention only (study attention): connective_tissue, immune, vascular_endothelial. [condition:me_cfs]
  ...

## 2. Where to measure it
- Burden of Long COVID (post-COVID-19 condition) at county level: evidence level A (direct condition measure), measure lc_current_pct_all_adults at state resolution; 3,144 of 3,144 geographies have a value (range 2.60-10.80 percent). These are STATE estimates inherited by counties (source resolution state), not county prevalence. [condition:long_covid]
- Burden of Postural orthostatic tachycardia syndrome at county level: UNKNOWN / NOT AVAILABLE (evidence level D) - No public population measure: no CCW algorithm contains G90.A; PLACES has no orthostatic item. [condition:pots]
- Candidate deployment opportunities for Long COVID / ME/CFS / dysautonomia-POTS-like phenotype x Wearable autonomic / activity monitoring (county level, 'equal' weights; 2,395 regions ranked; counties of the 50 states + DC with population >= 10,000): [...]
  1. Okmulgee County, Oklahoma: composite 0.888, rank interval 1-64 (P(top 10) 0.75); burden percentile 0.97 (evidence level C, includes a state value inherited by the county), vulnerability pct 0.92, diagnostic desert pct 0.99, clinic capacity pct 0.77, research readiness pct 0.79. [opportunity:autonomic_activity_invisible_illness|wearable_autonomic_activity_monitoring|40111; geo:40111]
  2. Oklahoma County, Oklahoma: composite 0.861, rank interval 1-324 (P(top 10) 0.47); ... [opportunity:autonomic_activity_invisible_illness|wearable_autonomic_activity_monitoring|40109; geo:40109]
  ...

## 3. Who could realistically deploy it
- Around Okmulgee County, Oklahoma for Long COVID (post-COVID-19 condition): 600 geocoded facilities within 50 km, 391 with a relevance signal. Facilities with characteristics suggesting they may be viable implementation or study partners (one per characteristic, round-robin; not a quality ranking): [geo:40111; condition:long_covid]
  - FAMILY CLINIC OF WELEETKA (WELEETKA, OK; 32.49 km; selected via community_access): Community access: community / rural-health / critical-access / public-health clinic taxonomy; its county has a primary-care HPSA designation (whole_county_population); its county is nonmetro (RUCC 8). [facility:npi-1699953687]
  ...
- For Long COVID (post-COVID-19 condition) x Wearable autonomic / activity monitoring [condition:long_covid], no facility in the pool has relevant_trial_history, NIH_research_activity, technology_experience: a possible measurement/research desert around this geography (absence in these public registries, not proof that no such activity exists). [geo:40111]
- Research readiness for Postural orthostatic tachycardia syndrome nationally (141 facilities with a trial site or NIH award; display order, not a ranking): University of Oklahoma Health Sciences Center (OKLAHOMA CITY, OK: 3 trials, 2 NIH core projects, 2 trials mentioning the measurement); ... [facility:npi-1508144411; ...]

## Provenance (top opportunity traced to sources and raw files)
- opportunity:autonomic_activity_invisible_illness|wearable_autonomic_activity_monitoring|40111 <- burden row long_covid|inherited_from_state|county:40111|gsea-w83j:lc_current_pct_all_adults:Oklahoma:By State:Oklahoma:72 (Currently experiencing long COVID, as a percentage of all adults [STATE value inherited as county context; not county prevalence], value 7.50 percent, evidence level A, source resolution state) from CDC NCHS Post-COVID Conditions (Household Pulse Survey). [...]
- Source cms_mmd: CMS Mapping Medicare Disparities (MMD) Tool - Population View (...); audit data/raw/cms_mmd/DATA_AUDIT.md; raw data/raw/cms_mmd/api/f_2022_v__prev_final_long_fltr12_racecat_all_sexcat_all_22_f.json. [...]
```

Capillaroscopy question (abbreviated):

```markdown
- Conditions: 'orthostatic intolerance' -> Dysautonomia (autonomic nervous system disorder) (exact match to alias 'orthostatic intolerance'; other candidates ['pots']) [condition:dysautonomia]
- Phenotype: 'orthostatic intolerance' (axis orthostatic_intolerance).
- nailfold_capillaroscopy is NOT linked to this phenotype in the phenotype-axis map: the Phase 3 table has no row for it, so its phenotype-level evidence is UNKNOWN / NOT AVAILABLE. [phenotype_evidence:orthostatic_autonomic|accelerometry; ...]
- Nailfold capillaroscopy (future CAPRIO adapter) for Dysautonomia (autonomic nervous system disorder): 10 registered trial(s) describe objective use (8 as an outcome measure) and 4 NIH core project(s) mention it; evidence tier high (research activity, not proof that it works). [measurement_evidence:dysautonomia|nailfold_capillaroscopy; trial:NCT01153581; trial:NCT03343574]
- FDA context for capillaroscopy: no_product_code_found; 0 mapped product code(s), UNKNOWN / NOT AVAILABLE 510(k), ... [measurement:capillaroscopy]
- Burden of Dysautonomia (autonomic nervous system disorder) at county level: UNKNOWN / NOT AVAILABLE (evidence level D) - No public population measure: G90.x absent from MMD-exposed CCW algorithms; no PLACES item. [condition:dysautonomia]
  1. Okmulgee County, Oklahoma: composite 0.861, rank interval 1-20 (P(top 10) 0.90); burden percentile UNKNOWN / NOT AVAILABLE (evidence level D), ... [opportunity:dysautonomia|nailfold_capillaroscopy|40111; geo:40111]
```

Unknown condition:

```markdown
**Stopped: no condition could be resolved, so no measurement, burden or deployment result can be attributed to it.** Nothing is ranked or attributed without a resolved input.

## UNKNOWN / NOT AVAILABLE
- Condition: UNKNOWN / NOT AVAILABLE - no registry condition matches 'Zorblaxian drift syndrome' (normalize_condition / search_condition).
```

## 8. Tests

| file | what |
|---|---|
| `tests/test_tools.py` | every SPEC tool exposed; each tool `ok` on a demo input and `UNKNOWN` on nonsense; envelope keys; JSON without NaN; no forbidden product language (the `validate` check); level-D burden UNKNOWN; inherited state values labelled; molecular output only ever called condition-level enrichment; ambiguous conditions not guessed; compact ranking keeps the SPEC schema; trials accept class / bundle / alias; memoisation; the agent on the capillaroscopy, unknown-condition, POTS, named-county and named-state questions (every sentence cites ids, every cited id traces; geography parsing; a named state is never presented as a within-state ranking); molecular entities and regulatory classes carry traceable object ids |
| `tests/test_mcp_server.py` | FastMCP in-memory client: tool list contains every SPEC tool with guardrail descriptions and described parameters; prompts / resources; every tool `ok` / `UNKNOWN` through MCP; `trace_evidence` round-trips every object id of a top-10 ranking; person-level refusal; schema rejects an invalid enum |
| `tests/test_trace.py` | every CONVENTIONS namespace resolves a real id with sources and existing DATA_AUDIT.md; unknown ids and namespaces UNKNOWN; person-level ids refused in every spelling (participant:, bare `nhanes:` / `stanford_...:` / `mapmecfs:` ids, any case); opportunity lineage reaches features, burden rows and raw files; signature lineage names dataset and label basis; trial trace finds its raw batch; facility trace lists source records |

`uv run pytest tests/test_tools.py tests/test_mcp_server.py tests/test_trace.py` (data tests need `data/processed`).

## 9. Limits

* The agent's parser is lexical: a condition, phenotype or measurement is found only when the question contains one of
  its registry terms (typos are not corrected); geography is found only as "<Name> County, <State>", "in <State>" or a
  FIPS / `geo:` id.
* "orthostatic intolerance" resolves to dysautonomia because the ontology lists it as a dysautonomia alias (POTS is the
  other candidate, reported); the Phase 3 phenotype is reported separately.
* Only precomputed condition x measurement x level combinations have stored `opportunity:` rows; other combinations are
  ranked on the fly by the same engine and cannot be traced by opportunity id.
* `rank_deployment_opportunities` ranks nationally; it cannot rank the counties of one state. The agent reports a named
  state's counties among the national top 50 and leaves a within-state ranking UNKNOWN.
* Trial, grant and Data.gov titles are quoted verbatim from their registries; Data.gov relevance ranking is not
  evaluated.
* Offline, `search_us_open_data` answers only the 54 catalog queries cached at ingestion (ignoring case and spacing).
* `get_measurement_evidence`, `get_phenotype_measurement_evidence` and a 10-region ranking are 80-90 kB of JSON (close
  to Claude Code's default MCP output limit, section 1).
