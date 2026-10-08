# API and dashboard

Two human- and machine-facing interfaces sit on the same tool facade as the MCP server (`docs/MCP.md`):

| interface | where | how it reaches the data |
|---|---|---|
| FastAPI app | `src/measure_it/api/app.py` (`measure_it.api.app:app`) | one GET and one POST route per tool, each calling the function of the same name in `measure_it.tools` and returning its JSON envelope unchanged |
| Streamlit dashboard | `dashboard/app.py` (home) + `dashboard/pages/1_...5_*.py`, helpers in `dashboard/components/` | imports `measure_it.tools` in-process (no HTTP); dense map layers read the same processed tables read-only (section 3) |
| display boundaries | `src/measure_it/api/geo.py` | simplified Census 2024 county / state polygons, shared by `/geojson/{level}` and the dashboard maps |

Neither interface computes a result. Every number on a page or in a response comes from a stored table through a
tool; when a tool has no data the answer is `UNKNOWN / NOT AVAILABLE` with the tool's reason, and nothing is filled in.

## 1. Running them

Both need the processed tables (`uv run measure-it pipeline --offline`, see `docs/REPRODUCIBILITY.md`).

```bash
uv run measure-it serve-api                    # http://127.0.0.1:8000, OpenAPI docs at /docs
uv run measure-it serve-api --offline --port 8080 --host 0.0.0.0
uv run uvicorn measure_it.api.app:app          # same app, plain uvicorn

uv run measure-it dashboard                    # http://localhost:8501 (opens a browser; headless without a TTY)
uv run measure-it dashboard --offline --headless --port 8502 --host 127.0.0.1
uv run streamlit run dashboard/app.py          # same dashboard, plain streamlit
```

`--offline` sets `MEASURE_IT_OFFLINE=1`: every tool answers from `data/processed`, `data/raw` and `data/_http_cache`
only, and `search_us_open_data` (the one tool that can reach the network) returns `UNKNOWN / NOT AVAILABLE` for a query
that is not in the HTTP cache (offline, a cached query is matched ignoring case and spacing). Both servers bind to
127.0.0.1 unless `--host` says otherwise. `measure-it dashboard` runs headless when stdin is not a terminal (Streamlit's
first-run e-mail prompt would otherwise stop it) and replaces itself with the Streamlit process, so Ctrl-C or
`kill <pid>` (the pid it prints) stops the server itself; `pkill -f` can match the calling shell.

Latency (aarch64, 20 cores, offline). The first call of the ranking in a process loads the scoring engine's spatial
base (326,685 facilities, 1.07 M individual NPIs; about 1.5 s) and takes 4-6 s; repeats are memoised by the tool facade
(< 0.01 s). A cold dashboard page renders in 2-8 s (Deployment Recommendation and Measurement Explorer are the slowest
because they rank), a warm one in about 1 s plus the browser's map drawing.

## 2. The API

| route | what it returns |
|---|---|
| `GET /` | index of the routes |
| `GET /health` | `status` (`ok`, or `degraded` when data/processed or SOURCE_REGISTRY.yaml is missing), tool count, processed-table count, data stamp, offline flag, the UNKNOWN sentinel |
| `GET /tools` | every tool: name, title, `spec` / `extra`, docstring, wrapped function, parameters (type, required, default, enum, description), routes |
| `GET /tools/{tool}` | the tool, arguments as query parameters (`?condition=Long%20COVID&geography_level=county`; lists repeat the key: `?statuses=RECRUITING&statuses=NOT_YET_RECRUITING`) |
| `POST /tools/{tool}` | the tool, arguments as one JSON object (`{"condition": "Long COVID or ME/CFS", "measurement": "wearable autonomic monitoring", "top_n": 10}`); unknown keys are rejected |
| `GET /sources` | `list_sources` envelope (SOURCE_REGISTRY.yaml summary); `?data_layer=person`; `?full=true` returns whole registry entries |
| `GET /sources/{source_id}` | one full registry entry and its DATA_AUDIT.md path; 404 with `status: UNKNOWN / NOT AVAILABLE` for an unknown id |
| `GET /geojson/{level}` | `level` = `county` (3,235 features, 1.8 MB) or `state` (56, 0.17 MB); `?state=06` for one state's counties; `?tolerance=` simplification in degrees (0-0.2; default 0.01 county, 0.02 state; 0 = none). `application/geo+json`, `ETag` + `Cache-Control: max-age=86400`, 304 on `If-None-Match` |

The 18 tools (15 SPEC tools + `get_phenotype_measurement_evidence`, `get_measurable_biology`, `list_sources`) and their
arguments are documented in `docs/MCP.md` section 3 and at `/docs`. Enumerations match the MCP server
(`geography_level`, `detail`, `order`); argument descriptions are the MCP server's.

**Response.** The body is the tool envelope `{tool, status, query, data, caveats, provenance, truncation, reason?, meta}`
(docs/MCP.md section 2), byte-for-byte the JSON of `measure_it.tools.<tool>(...)` apart from `meta` (timings).
HTTP status: 200 for `status: "ok"` and for `status: "UNKNOWN / NOT AVAILABLE"` (a valid answer: the data are missing
and `reason` says why); 500 for `status: "error"` (same body); 422 when an argument fails validation (missing, wrong
type, not in the enum, unknown JSON key) before the tool runs. The `X-Tool-Status` header repeats the envelope status.
CORS allows GET/POST from any origin (read-only public data).

```bash
$ curl -s "http://127.0.0.1:8000/tools/get_condition_burden?condition=POTS&geography_level=county"   # HTTP 200
{"tool": "get_condition_burden", "status": "UNKNOWN / NOT AVAILABLE",
 "reason": "No public population measure: no CCW algorithm contains G90.A; PLACES has no orthostatic item.", ...}

$ curl -s "http://127.0.0.1:8000/tools/get_condition_burden?condition=Long%20COVID&geography_level=county&geo_id=06073"
... "rows": [{"geo_id": "06073", "geo_name": "San Diego County, California", "value": 3.8, "value_unit": "percent",
              "burden_evidence_level": "A", "source_geographic_resolution": "state", "inherited": true,
              "object_id": "geo:06073", ...}] ...

$ curl -s -X POST http://127.0.0.1:8000/tools/rank_deployment_opportunities -H 'content-type: application/json' \
       -d '{"condition": "Long COVID or ME/CFS", "measurement": "wearable autonomic monitoring", "top_n": 10}'
```

The last value above is the California state Household Pulse Survey estimate carried by the county: `inherited: true`,
source resolution `state`. It is context, not San Diego County prevalence.

**Boundaries** (`measure_it.api.geo`): Census 2024 cartographic boundaries (1:5M; SOURCE_REGISTRY `census_geography`),
simplified per feature with shapely (`preserve_topology=True` within a feature, not across neighbours) and snapped to a
0.001-degree grid. Feature `id` = FIPS; properties `geo_id`, `name`, `state_abbr`, `object_id` (`geo:<fips>`) and
`in_ranking_universe` (50 states + DC). Neighbouring polygons can leave hairline gaps: the geometry is for drawing
only, never for spatial joins (the pipeline uses the unsimplified files). Built once per process and argument set
(about 0.5 s for counties) and served from memory.

## 3. The dashboard

Every page starts with its inputs and ends with a **Limitations and provenance** panel: page-specific limitations, the
four-data-layer rule, every UNKNOWN answer on the page with its reason, the caveats the tools returned, the tool calls
(tool, status, query, producer, tables, number of object ids, capped lists, time), the SOURCE_REGISTRY sources behind
them (version, retrieval date, DATA_AUDIT.md path), and the processed tables the page read directly. Pages that cite
object ids have a **trace evidence** expander: pick any id and `trace_evidence` resolves it to its row, sources,
DATA_AUDIT.md and one level of lineage.

Maps are Plotly geo maps (Albers USA for the nation, Mercator zoomed to a geography; no tile server). Choropleths use
one blue ramp (darker = higher); grey means no value, and hovering it gives the reason (level-D burden, a county below
10,000 residents that is not ranked, an incomplete member burden). Point layers use distinct colours and symbols
(orange circles = HRSA sites, aqua diamonds = trial sites, violet squares = NIH organisations or candidate
facilities, grey open circles = specialist density; the overlay colours pass the colour-vision-deficiency check of the
project's palette validator). Boundaries are the simplified geometry of section 2 (0.02 degrees for the national county
map, 1.4 MB) and are cached per process (`st.cache_resource`); table reads are cached with `st.cache_data` and keyed by
the tool facade's data stamp, so a rebuilt table invalidates them.

### Home (`dashboard/app.py`)

The thesis, the engine chain (person-level data -> objective phenotype -> ontology -> condition-level molecular
enrichment -> candidate measurement -> geographic burden -> clinic / research readiness -> candidate deployment
opportunity), links to the five pages, the source registry (`list_sources`: 27 sources by data layer) and the
data-layer table of the MCP server (which layers may be joined, and how).

### 1. National Opportunity Map (`pages/1_National_Opportunity_Map.py`)

* **Choropleth selector (one area layer at a time):** *Deployment opportunity* (the `composite_<weight set>` of `deployment_opportunities`, only
  for ranked regions), *Diagnostic desert* (the geography module's index; access-only variant for a level-D condition)
  or *Burden (with evidence level)* (the primary measure; a condition set shows the percentile of its members' mean
  burden percentile). County or state level.
* **Selectors:** condition or condition set (the six scored ones first; other registry conditions for burden and
  desert only), measurement bundle (the four scored bundles, with the adapter status; nailfold capillaroscopy is the
  stub adapter), weight set (the eight sets of `configs/scoring.yaml`).
* **Independent overlay toggles:** HRSA health-center service delivery sites (source coordinates); U.S. sites of
  trials registered for the condition (literal condition match; city-level geopoints aggregated per location; optional
  active/recruiting filter); NIH RePORTER awardee organisations with title/abstract matches (likely false positives
  excluded); specialist density (distinct individual NPIs in the condition's core specialty groups other than primary
  care per 100k, via `scoring.opportunity.provider_counts_for_groups`, the geography module's `relevant_specialists`
  definition).
* **Against the SPEC (partial):** the SPEC asks for independently toggled burden, diagnostic-desert, opportunity,
  facility and trial-site layers. Burden, desert and opportunity share the one choropleth above (mutually exclusive),
  and the facility overlay shows HRSA sites only, not NPPES or `clinic_registry` facilities.
* **Always visible:** the burden-evidence box (`get_condition_burden` per member): evidence level, measure and source,
  source resolution, the inherited-state-value caveat at county level, the proxy caveat for level C, the level-D
  consequence. The top 10 under the chosen weights (`rank_deployment_opportunities`) are outlined in red.
* **Table + CSV** of the chosen layer (ranked table with every component percentile right after the composite, then
  the equal-weight Monte Carlo rank interval and P(top 10), labelled as such for every weight set).
* **Inherited state values:** on the burden layer at county level, a condition whose counties carry their state
  estimate (Long COVID) is labelled on the colour bar itself ("inherited STATE value on counties") and a warning gives
  the share of counties affected; a condition set with such a member says "incl. inherited state value".
* Tools: `get_condition_burden`, `rank_deployment_opportunities`. Tables read directly: `deployment_opportunities`,
  `geo_condition_features`, `facilities__hrsa`, `trial_conditions` + `clinical_trials` + `trial_sites`,
  `nih_project_conditions` + `nih_projects`, `providers` (through the scoring engine).

### 2. Condition Explorer (`pages/2_Condition_Explorer.py`)

Input: condition text (default "Long COVID"; aliases, acronyms, ICD-10-CM codes and MONDO ids work), an optional
phenotype axis for the candidate measurements. `normalize_condition` resolves it (ambiguous or unknown text gives
UNKNOWN with the `search_condition` candidates and the page stops). Six sections:

1. **Phenotype features** (`get_patient_phenotype_signature`): per public dataset and label: definition, label basis,
   proxy flag, cases / controls, features tested, FDR-significant count, null results; per dataset the significant
   features (effect, CI, q, `signature:` id), the null features and model-level results.
2. **Wearable abnormalities**: forest plot of the FDR-significant wearable features (effect and 95% CI, one colour per
   dataset and label, proxy labels marked) and the null features per dataset. For Long COVID the only significant
   feature comes from an acute-infection cohort that is labelled "not a Long COVID label"; the Long COVID label itself
   (Uwakwe 2025, 31 cases) is a null result, shown as such.
3. **Molecular enrichment** (`get_molecular_context`), labelled *condition-level molecular enrichment, not patient
   multi-omics*: evidence types with record counts by source, top entities, Test 3 coherence (supported physiological
   systems with FDR and the post-hoc flag, gene-list agreement, cross-source genetic agreement).
4. **Candidate measurements** (`discover_candidate_measurements`): the five dimensions in separate columns next to
   the measurement (phenotype signal strength, measurement evidence strength, technology maturity, regulatory
   visibility, deployment complexity as a curated ordinal "n of 5"), then the counts behind them (objective-use
   trials, NIH core projects, FDA decisions), text-mining precision, the ordering note (research activity, not
   validity) and the reasons. A `null result` phenotype signal is shown as such.
5. **Technologies (FDA)** (`get_regulatory_context` per candidate class): one row per class (number of mapped product
   codes, device classes, best mapping confidence, 510(k) / De Novo / PMA totals, example codes); every product code
   (device, class, regulation, counts, mapping confidence) in an expander.
6. **Burden sources** (`get_condition_burden` at national, state and county level): evidence level, primary measure,
   source resolution, inherited flag, summary, rationale; the highest state values with their 95% CIs.

### 3. Measurement Explorer (`pages/3_Measurement_Explorer.py`)

Input: measurement class or bundle (default "wearable autonomic monitoring"), the condition used for trials and unmet
need, the unmet-need level and weight set (default `access_gap`: 0.4 on the diagnostic desert, 0 on research
readiness). Sections:

1. **Conditions where it is studied**: `get_measurement_evidence(measurement, condition)` for all 14 registry
   conditions: registered trials describing objective use (literal condition match), as outcome measure, recruiting,
   NIH core projects, evidence tier, phenotype-signal status; bar chart. A condition the tool answers with UNKNOWN /
   NOT AVAILABLE (the post-infectious grouping) is listed under the table with the tool's reason, not as a row of
   blanks.
2. **Signals & public datasets**: the registry's physiological signals and public person-level datasets
   (`trace_evidence("measurement:<id>")`) and the phenotype-signal results (with their role: an acute COVID-19 cohort
   is adjacent evidence, not the condition).
3. **Relevant trials** (`find_relevant_trials`): counts by status and type, the trials with `trial:` ids.
4. **FDA device classes** (`get_regulatory_context` per member class).
5. **Regions with unmet need** (`rank_deployment_opportunities`): top 10 with their components and a map of every
   ranked region's composite (precomputed combinations) or of the returned regions (combinations scored on the fly).
   The Monte Carlo rank interval is drawn around the equal weights only, so under `access_gap` (the default here) the
   column is "rank interval 5-95 (equal-weight MC)" and a caption says a rank can fall outside it.

### 4. Geography Explorer (`pages/4_Geography_Explorer.py`)

Input: county or state (default "San Diego County, California"; FIPS, `geo:` ids and state names work), the condition
and measurement used for facilities, and the radius (default 50 km). `get_geographic_context` answers most of it:

1. **Burden by condition** with value, unit (from `get_condition_burden` for the geography), evidence level, inherited
   flag, source resolution, CI, measure, diagnostic desert (and access-only variant), and the reason when unknown;
   usable levels first. An inherited value names itself in the value cell ("3.8 (CA STATE value, inherited)") and a
   warning counts the conditions whose county value is an inherited state estimate.
2. **Vulnerability & context**: SVI overall and the four themes (percentiles), ACS and PLACES context.
3. **Providers & specialists**: individual NPPES providers per 100k by specialty group; per condition the relevant
   specialists excluding primary care (count and per 100k) first, then the counts including primary care and the
   specialty groups.
4. **Trials & NIH research**: relevant trials with a site here (active, recruiting, completed, by condition) and
   NIH-funded research (award records, core projects, organisations).
5. **Candidate technologies & facilities**: `discover_candidate_measurements` for the condition (condition-level, the
   same list for every geography) and `find_candidate_clinics` (six characteristics ranked separately, reasons, absence
   note when a characteristic is missing in the pool).
6. **Nearby health centers (map)**: the geography's polygon, HRSA sites inside the county or within the radius of its
   internal point (inside the state for a state), and the candidate facilities; the HRSA table with distances. The
   view is the bounding box of the polygon and the points; neighbouring counties (or states) are drawn in grey from the
   same simplified Census boundaries, so no basemap coastline at another resolution is drawn over them (Alaska, which
   crosses the antimeridian, falls back to fitting its polygon).

### 5. Deployment Recommendation (`pages/5_Deployment_Recommendation.py`)

Default query: the SPEC example, "Find 10 U.S. regions where wearable monitoring of autonomic/activity abnormalities in
Long COVID or ME/CFS would be useful to evaluate, and identify clinics/research sites that could plausibly
participate." (`rank_deployment_opportunities("Long COVID or ME/CFS", "wearable autonomic monitoring", "county",
top_n=10, weight_set="equal", detail="full")`). Parameters: condition, measurement, level, weight set, number of
regions, candidate sites per region on the map, small counties.

* **Ranked table** with every component: rank, region, composite, then the five component percentiles (burden,
  vulnerability, desert, capacity, readiness), the rank interval (5th-95th percentile over 1,000 Monte Carlo draws of
  burden and weights), P(top 10), burden evidence level (and whether a member is an inherited state value), member
  burden values, state / FIPS / population, SVI, the desert / capacity / readiness indices, condition trials and NIH
  core projects in reach, technology saturation (reported only), the best and worst rank across the eight weight
  sets, and the `opportunity:` id (for a stored ranking; an on-the-fly ranking says "computed on the fly; not a stored
  row"). The Monte Carlo is drawn around the equal weights only: under another weight set the
  interval and P(top 10) columns are labelled "(equal-weight MC)" and the captions say the rank can fall outside the
  interval.
* **Charts:** each component's contribution (weight x percentile / sum of the weights of the available components;
  the bars add up to the composite) and the rank with its interval.
* **Map:** every ranked region's composite, the returned regions outlined, candidate sites (the lat / lon and geocode
  precision that `rank_deployment_opportunities` returns with each candidate site; before 2026-09-24 the page looked
  them up in `facilities` by facility id).
* **Evidence for one region:** the recommended next step (framed as a candidate deployment opportunity), candidate
  sites with their per-condition characteristic values and ranks and reasons, research evidence (trial / NIH ids),
  uncertainties, provenance with the trace_evidence expander.
* **Condition- and measurement-level context:** phenotype signatures (cases, controls, features tested,
  FDR-significant count and the null-result summary first, e.g. Long COVID 0 of 14), condition-level molecular
  enrichment, technology (regulatory context, technology evidence).
* **Regions not ranked** because their burden is incomplete, with the reason.
* **Recommendation JSON** download: the compact envelope (SPEC schema per region, shared context once) and the
  full-detail envelope (all sites and context), identical to the MCP tool and `/tools/rank_deployment_opportunities`.

### Page 6: Your data

`dashboard/pages/6_Your_Data.py` lists the user-supplied datasets added with `measure-it byod` (docs/BRING_YOUR_OWN_DATA.md)
and, for the selected one: SYNTHETIC / demo / comparator badges (a healthy-control comparator is shown as a warning), the
locked plan (sha256, time, plan JSON), the evaluation models, the tier-2.5 performance record and its caveats, the
ranking change of the last `byod deploy`, and the heatmap (condition x measurement evidence, user cells outlined) and
embedding (the dataset's own Digital Phenotype Vector PCA) images from `results/byod/<id>/`. It reads
`measure_it.byod.tools` (aggregate only); with no dataset ingested it says UNKNOWN / NOT AVAILABLE and shows the
commands of docs/WALKTHROUGH.md.

## 4. Caveats that apply to every page

* **Four data layers.** Person-level public cohorts (NHANES, Stanford wearable studies, MapMECFS) appear only as group
  statistics; they are never located and never joined to a place, facility or trial. Molecular evidence is
  condition-level molecular enrichment (public databases), never patient multi-omics. Geographic joins are ecological.
* **Burden evidence levels.** A direct, B coded condition, C symptom/comorbidity proxy, D none. Long COVID has no county
  measure: counties carry the state Household Pulse Survey estimate (inherited context, not county prevalence).
  ME/CFS burden is a level-C claims proxy (CMS "Fibromyalgia, Chronic Pain and Fatigue"); POTS and dysautonomia are
  level D, so their rankings contain no burden information.
* **Rankings are weight-sensitive and wide.** The rank intervals of the top 10 span hundreds of places
  (results/SCORING_RESULTS.md, Tests 5-7); the short list is a set of candidate deployment opportunities for pilot
  evaluation, not a validated diagnostic pathway.
* **Signals, not efficacy.** Registered trials and NIH grants are research-activity signals; FDA records are
  deployment-readiness signals; neither shows that a measurement works or detects a condition. Counts come from text
  mining with an audited precision per class.
* **Curated assumptions.** Condition specialty groups and measurement implementer groups (`configs/relevance.yaml`);
  an NPPES taxonomy does not mean a clinician evaluates or treats the condition.
* **Locations are centroids.** Facilities at ZIP/ZCTA or city centroids (HRSA sites at their source coordinates); the
  50 km radius is measured from the county's Census internal point.
* **Direct table reads are display-only.** The dense layers (choropleths, national point layers) read
  `data/processed` directly because the tools cap lists at 100 items; they select and aggregate rows
  and never compute a score. Albers USA does not draw Puerto Rico or the territories, which are outside the ranking
  universe (50 states + DC) anyway.

## 5. Tests and screenshots

```bash
uv run pytest tests/test_api.py tests/test_dashboard.py      # 48 + 11 tests; data tests need data/processed
```

`tests/test_api.py` calls every route with the FastAPI TestClient: GET and POST per tool with demo arguments (the body
must equal the facade's envelope apart from `meta`) and with inputs that must return UNKNOWN / NOT AVAILABLE with a
reason; 422 validation cases; the 500 mapping of an `error` envelope; `/health`, `/tools`, `/sources`,
`/sources/{id}`, `/geojson/{level}` (feature ids, grid precision, state filter, ETag 304). `tests/test_dashboard.py`
runs the home page and each page with `streamlit.testing.v1.AppTest` using the default and alternative inputs (other
levels, layers, weight sets, overlays, conditions, an ICD-10-CM code, a state, and nonsense input that must render
UNKNOWN) and checks that there is no exception, that the key elements and the Limitations and provenance panel
render, and that the rendered text passes the product-language check.

Screenshots of every page (default inputs, light theme, 1500 px wide; `dashboard/screenshots.py` writes full pages to `results/figures/dashboard/`, and the top of each page is published in `visuals/dashboard/`):
`home.png`, `national_opportunity_map.png`, `condition_explorer.png`, `measurement_explorer.png`,
`geography_explorer.png`, `deployment_recommendation.png`. They were captured on this aarch64 host (Ubuntu 24.04)
with Playwright's headless Chromium; to regenerate them:

```bash
uv run --with playwright playwright install chromium        # once: downloads Chromium headless shell (~115 MB)
uv run --with playwright python dashboard/screenshots.py    # starts streamlit offline on a free port, captures, stops
uv run --with playwright python dashboard/screenshots.py --only deployment_recommendation
uv run --with playwright python dashboard/screenshots.py --url http://127.0.0.1:8501 --out /tmp/shots   # a running server, e.g. the container
```

`dashboard/screenshots.py` waits until Streamlit's running indicator is gone, every expected Plotly chart has drawn
and the page height is stable, then grows the viewport to the content height (Streamlit scrolls inside its own
container) before the capture. If a headless browser cannot be installed, the tests above still exercise every page
headlessly; the screenshots can be taken from any machine with a browser that can reach `uv run measure-it dashboard`.
