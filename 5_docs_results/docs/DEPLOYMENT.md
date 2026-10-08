# Deployment

How to run the API, the dashboard and the MCP server locally, in containers from a prebuilt data bundle, and as a
static read-only snapshot for reviewers. Everything here is read-only: no service writes to the data, and with
`MEASURE_IT_OFFLINE=1` (always set in the container) no tool reaches the network.

| service | command | port | health |
|---|---|---|---|
| API (FastAPI) | `measure-it serve-api` | 8000 | `GET /health` (`status: ok`, or `degraded` without data) |
| dashboard (Streamlit) | `measure-it dashboard` | 8501 | `GET /_stcore/health` (`ok`) |
| MCP server (FastMCP) | `measure-it serve-mcp` (stdio) or `--transport http` | 8765 (HTTP only) | `GET /health` (HTTP transport only) |

The three share one tool facade (`measure_it.tools`), described in `docs/DASHBOARD.md` and `docs/MCP.md`.

## 1. Local run (a checkout with built data)

```bash
uv sync --frozen
uv run measure-it serve-api --offline                 # http://127.0.0.1:8000, OpenAPI at /docs
uv run measure-it dashboard --offline --headless      # http://127.0.0.1:8501
uv run measure-it serve-mcp                           # stdio (what .mcp.json starts)
uv run measure-it serve-mcp --transport http --offline --port 8765   # streamable HTTP at http://127.0.0.1:8765/mcp
```

All three bind to 127.0.0.1 unless `--host` says otherwise. The data must exist first (`uv run measure-it pipeline
--offline --workers 16`, README "Quickstart"), or be unpacked from a data bundle (section 3) at the project root.

## 2. Container run

One image (`Dockerfile`, multi-stage: `python:3.12-slim-bookworm` + uv, `uv sync --frozen --no-dev`) runs every
service; `docker-compose.yml` starts three containers from it. The image holds code, configs and docs, never data.

```bash
uv run measure-it bundle --unpack dist/bundle     # build the data bundle and unpack it (section 3)
docker compose up -d --build                      # API :8000, dashboard :8501, MCP (HTTP) :8765 on 127.0.0.1
docker compose ps                                 # all three should say (healthy) after about 10-60 s
uv run python deploy/smoke_test.py                # /health x3, every tool route (demo + UNKNOWN input), MCP over HTTP
uv run --with playwright python dashboard/screenshots.py --url http://127.0.0.1:8501 --out /tmp/shots   # headless pages
docker compose down
```

* **Data.** Compose mounts an unpacked bundle read-only: `$MEASURE_IT_BUNDLE/data`, `results`, `configs` and
  `SOURCE_REGISTRY.yaml` onto `/app/...` (`MEASURE_IT_BUNDLE` defaults to `./dist/bundle`). The configs and registry come
  from the bundle so that they always match the tables they were built with.
* **Hardening in compose.** Read-only root filesystem (`read_only: true`; only `/tmp` and the home directory are
  tmpfs), non-root user (uid 10001), `cap_drop: ALL`, `no-new-privileges`, ports published on 127.0.0.1 only,
  `restart: unless-stopped`, 6 GB memory limit per service.
* **Ports.** `MEASURE_IT_API_PORT`, `MEASURE_IT_DASHBOARD_PORT` and `MEASURE_IT_MCP_PORT` change the host ports, e.g.
  `MEASURE_IT_API_PORT=18000 docker compose up -d`.
* **MCP over stdio.** An MCP client that starts its server as a subprocess (Claude Desktop, Claude Code) can run the
  container with `-i` and no network at all:

  ```bash
  B=$PWD/dist/bundle
  docker run -i --rm --read-only --tmpfs /tmp --network none \
    -v "$B/data:/app/data:ro" -v "$B/results:/app/results:ro" -v "$B/configs:/app/configs:ro" \
    -v "$B/SOURCE_REGISTRY.yaml:/app/SOURCE_REGISTRY.yaml:ro" measure-it:local measure-it serve-mcp
  ```

  Put that command and its arguments in the client's `mcpServers` entry (`"command": "docker", "args": ["run", "-i", ...]`).
  The HTTP service of compose (`http://127.0.0.1:8765/mcp`, streamable HTTP) serves clients that connect by URL.
* **Baked image (optional).** To ship one self-contained image instead of image + bundle:
  `docker build --target baked --build-context bundle=dist/bundle -t measure-it:baked .` copies the unpacked bundle into
  `/app` (4.21 GB; built and checked offline with `--network none` on 2026-09-25). Run it without the volume mounts.
* **Check it was verified.** On 2026-09-25 (aarch64 host, Docker 29.2.1, Compose v5.0.2) the image was built, the stack
  came up with all three health checks `healthy`, `deploy/smoke_test.py` passed 36 of 36 tool calls (18 tools x demo
  input and UNKNOWN input) plus MCP over HTTP (18 tools listed, one call answered), the six dashboard pages rendered
  headlessly with no exception, and the stdio MCP server answered from a container with `--network none`. The image
  has not been built on x86_64.

## 3. The data bundle

`uv run measure-it bundle` (`src/measure_it/export/bundle.py`) writes
`dist/measure-it-data-<UTC date>.tar.zst` (`.tar.gz` when `zstd` is not installed), `...manifest.json` and `...sha256`.

| part | contents |
|---|---|
| `data/processed` | every `*.parquet`, `*.geoparquet`, `*.meta.json` and `measure_it_public.duckdb` (leave the DuckDB out with `--no-duckdb`: -1.4 GB; point lookups then read the Parquet tables) |
| `data/raw/<source>/` | only the provenance sidecars `trace_evidence` and the source summaries read (`MANIFEST.json`, `DATA_AUDIT.md`, `registry_entry.yaml`, `QUERY_LOG.json`, `queries.json`), plus `data/raw/data_gov_catalog` (stored Data.gov responses) and `data/raw/clinicaltrials_gov/records` (the registry batches a `trial:` id traces to). No raw data file |
| `data/_http_cache` | only the cached Data.gov catalog responses, so `search_us_open_data` answers its cached queries offline |
| `results` | `*.md` reports, `*.json` examples, `results/tables/**` (no figures, maps or pipeline run logs) |
| `configs`, `SOURCE_REGISTRY.yaml` | the configuration and registry the tables were built with |
| `BUNDLE_MANIFEST.json` | last archive member (also written beside the archive): per file size, SHA-256 of the archived bytes and row count (Parquet metadata, CSV records, rows per DuckDB table); git commit and dirty state; every source's version and retrieval date; what was included and excluded; files that changed while archiving |

The bundle of 2026-09-25 (commit `1c72464`, working tree dirty): 1,045 files, 2.48 GB unpacked, **1.00 GB** as
`.tar.zst` (0.997 GB), built in a few seconds on the 20-core host. Directory entries carry their source mtimes, so the
services' data stamp (`/health`) is the build's, not the unpack time.

```bash
uv run measure-it bundle                           # build (refuses when data/processed holds no table)
uv run measure-it bundle --unpack dist/bundle      # build and unpack
uv run measure-it bundle --verify dist/bundle      # re-hash every file against BUNDLE_MANIFEST.json (exit 1 on mismatch)
tar --zstd -xf dist/measure-it-data-<date>.tar.zst -C dist/bundle   # unpack by hand
```

**Refreshing the data.** Two routes:

* *Full offline rebuild* (a machine with `data/raw`, `data/_http_cache`, about 45 GB): `uv run measure-it pipeline
  --offline --workers 16` (91-104 min for the 49-step DAG measured; the 75-step DAG not timed), `uv run measure-it
  validate`, then `uv run measure-it bundle`. A cold build from the network fetches newer source versions and gives
  different numbers (docs/REPRODUCIBILITY.md sections 6-8).
* *Bundle only* (a deployment host): copy the new archive, check it against its `.sha256`, unpack into a fresh
  directory, `measure-it bundle --verify` it, point `MEASURE_IT_BUNDLE` at it and `docker compose up -d` (the services
  pick the new tables up on restart; the tool facade also invalidates its caches when the data stamp changes). Keep the
  previous directory until the new one is healthy, to roll back.

Build the bundle only when no pipeline is running: a file rewritten while it is archived is listed under
`files_changed_during_build` and the builder warns.

## 4. Static reviewer snapshot

`uv run measure-it export-static` (`src/measure_it/export/static_site.py`) writes `dist/static_snapshot/index.html`, one
self-contained page (about 1.8 MB: inline CSS, JS, data and simplified county boundaries) with the overview, the
primary deployment query (ranked table with components and rank intervals, evidence_weighted vs equal, ranks under
every weight set), a county choropleth with the top regions and candidate sites, evidence cards for the top 10, the
measurement performance table (UNKNOWN shown as UNKNOWN), results/FINAL_REPORT.md §14, the validation tests (§13) and
the limitations (§19). Every number is read when the page is built; nothing on it is live. Its only external request is
the pinned d3 script from cdn.jsdelivr.net for the map; without it the tables still render and the map says so. It
follows the page's light / dark colour scheme and has no horizontal page scroll at phone width.

```bash
uv run measure-it export-static
uv run --with playwright python scripts/screenshot_static_snapshot.py    # results/figures/static_snapshot_{desktop,mobile}.png
```

## 5. Resources

| item | measured |
|---|---|
| image | 1.73 GB (aarch64; the Python environment dominates: scipy, scikit-learn, numba, geopandas, streamlit) |
| bundle | 1.00 GB compressed, 2.48 GB unpacked (DuckDB 1.4 GB of it) |
| memory | after the smoke test and all dashboard pages: API 3.6 GB, dashboard 3.6 GB, MCP 0.5 GB (the scoring engine loads 326,685 facilities and 1.07 M individual NPIs on the first ranking call, about 1.5 s); compose caps each at 6 GB |
| CPU | idle when not answering; the first ranking call in a process takes 4-7 s, repeats are memoised |
| start-up | healthy within about 10 s; the first request to a tool loads its tables |

## 6. Security and privacy notes

* **No authentication, no authorisation, no TLS.** The services are meant for a laptop or a trusted network. Compose
  publishes them on 127.0.0.1 only. To expose them, put a reverse proxy in front (nginx, Caddy, Traefik, a cloud load
  balancer) that terminates TLS and requires authentication (SSO / OIDC or at least basic auth), and rate-limit the
  API: a ranking call costs several seconds of CPU. The API's CORS policy allows any origin (read-only public data);
  narrow it behind the proxy if the API is not meant to be called from other sites. The Streamlit dashboard and the
  MCP HTTP endpoint have no access control of their own either.
* **Read-only public data.** Every source is public (SOURCE_REGISTRY.yaml, `access_conditions`); person-level rows come
  from public, de-identified research releases and are summarised as group statistics; nothing in the bundle is
  protected health information, and no route accepts or stores user data. The two public lab records that appeared to
  contain patient identifiers (`PRIVACY_EXCLUDED` in `measure_it.labs.lab_dataset_discovery`) were never downloaded;
  the bundle builder refuses any path that names them, and the manifest lists them under `privacy_excluded_records`.
* **Offline.** `MEASURE_IT_OFFLINE=1` makes every tool answer from the bundle; `search_us_open_data` returns
  `UNKNOWN / NOT AVAILABLE` for a query that is not in the cached Data.gov responses. Streamlit's usage statistics and
  external-address lookup are switched off in the image.
* **Containers.** Non-root user, read-only root filesystem, no capabilities, data volumes mounted read-only.

## 7. What a production deployment would still need

* Authentication and TLS at a reverse proxy, rate limiting and request-size limits; an audit log of queries.
* A published, signed data bundle (object storage + checksum/signature) and a release process that ties bundle,
  image tag and git commit together; image scanning and an SBOM; builds for x86_64 as well as aarch64 (only aarch64 was
  built).
* A scheduled refresh: the offline pipeline cannot pick up new source versions, and the cold path has never been run
  end to end (README "Reproducing from scratch"). Time-varying sources (ClinicalTrials.gov, RePORTER, openFDA, Open
  Targets) need a refresh cadence and a changelog of how rankings moved.
* Monitoring beyond container health checks (latency, error rate, memory; the ranking engine's 3-4 GB per process
  limits how many workers fit on a host), log shipping, and horizontal scaling with the bundle on shared read-only
  storage.
* A tested backup / rollback of bundles, and a documented support and disclosure contact.
* A decision about the dashboard: Streamlit keeps one Python session per browser tab; for many concurrent users a
  cached static page (section 4) or a front end over the API scales better.
* Review of the scientific limits before any operational use: rankings are candidate deployment opportunities for pilot
  evaluation, weight-sensitive and wide (results/FINAL_REPORT.md §13, §19), not a validated diagnostic pathway.
