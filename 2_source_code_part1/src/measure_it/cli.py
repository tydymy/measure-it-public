"""Command-line entry point: `uv run measure-it <command>`.

    measure-it pipeline [--list] [--only STEP ...] [--from STEP] [--offline] [--workers N]
    measure-it build-db
    measure-it validate
    measure-it merge-registry
    measure-it serve-mcp [--transport stdio|http|sse] [--host 127.0.0.1] [--port 8765] [--offline]
    measure-it serve-api [--host 127.0.0.1] [--port 8000] [--offline] [--reload]
    measure-it dashboard [--host 127.0.0.1] [--port 8501] [--offline] [--headless]
    measure-it ask "QUESTION" [--json] [--top-n N] [--level county|state]
    measure-it bundle [--out dist] [--no-duckdb] [--unpack DIR] [--verify DIR]
    measure-it export-static [--out dist/static_snapshot/index.html]
    measure-it byod validate|ingest|evaluate|subgroup|deploy|remove|list|show ...   (bring your own data; docs/BRING_YOUR_OWN_DATA.md)
"""
from __future__ import annotations

import os
from typing import List, Optional

import typer

app = typer.Typer(no_args_is_help=True, add_completion=False)


def _register_byod() -> None:
    from .byod.cli import app as byod_app
    app.add_typer(byod_app, name="byod")


_register_byod()
app.add_typer(__import__("measure_it.agent.cli", fromlist=["app"]).app, name="launch")  # docs/LAUNCH_AGENT.md


@app.command("pipeline")
def pipeline(
    list_steps: bool = typer.Option(False, "--list", help="Print the ordered steps (entry point, dependencies) and exit."),
    only: Optional[List[str]] = typer.Option(None, "--only", help="Run only these steps (repeatable or comma-separated); "
                                                                  "their inputs must already exist."),
    from_step: Optional[str] = typer.Option(None, "--from", help="Run this step and every step after it."),
    offline: bool = typer.Option(False, "--offline", help="MEASURE_IT_OFFLINE=1: use data/raw and data/_http_cache "
                                                          "only; any cache miss or missing raw file fails the step."),
    workers: int = typer.Option(1, "--workers", min=1, help="Worker budget shared by concurrently running steps; a "
                                                            "step's own worker count is min(its default, N)."),
) -> None:
    """Run the explicit, ordered DAG of every step that produces data/processed and results/."""
    from . import pipeline as P
    only = [x.strip() for v in (only or []) for x in v.split(",") if x.strip()] or None   # --only a,b or repeated
    if list_steps:
        typer.echo(P.describe(P.select_steps(only, from_step) if (only or from_step) else None))
        raise typer.Exit(0)
    if offline:
        os.environ[P.OFFLINE_ENV] = "1"
    raise typer.Exit(P.run(only=only, from_step=from_step, offline=offline, workers=workers, echo=typer.echo))


@app.command("merge-registry")
def merge_registry() -> None:
    """Merge data/raw/*/registry_entry.yaml into SOURCE_REGISTRY.yaml."""
    from .registry import merge_registry_fragments
    reg = merge_registry_fragments()
    typer.echo(f"{len(reg['sources'])} sources in SOURCE_REGISTRY.yaml")


@app.command("build-db")
def build_db() -> None:
    """(Re)build data/processed/measure_it_public.duckdb: every processed table, SPEC core-data-model views
    (main.conditions, main.measurements and schema spec.*), _sources, _catalog, _spec_views."""
    from .database import build_db as _build
    typer.echo(str(_build()))   # exits non-zero (SystemExit) when a SPEC view has no underlying table


@app.command("validate")
def validate() -> None:
    """Provenance, required datasets, person-layer geography, object ids, product language, data audits.
    Exits non-zero when any check fails."""
    from .validate import run_checks
    report = run_checks()
    raise typer.Exit(0 if report["ok"] else 1)


@app.command("serve-mcp")
def serve_mcp(
    transport: str = typer.Option("stdio", "--transport", help="stdio (default; what .mcp.json starts), http "
                                                                "(streamable HTTP at /mcp, plus GET /health) or sse."),
    host: str = typer.Option("127.0.0.1", "--host", help="Interface for http/sse (0.0.0.0 exposes it)."),
    port: int = typer.Option(8765, "--port", help="Port for http/sse."),
    offline: bool = typer.Option(False, "--offline", help="MEASURE_IT_OFFLINE=1 for every tool call."),
) -> None:
    """Start the MCP server (FastMCP) exposing every agentic-layer tool; see docs/MCP.md, .mcp.json and
    docs/DEPLOYMENT.md (HTTP transport in the container)."""
    if transport not in ("stdio", "http", "streamable-http", "sse"):
        raise typer.BadParameter("transport must be stdio, http, streamable-http or sse")
    if offline:
        os.environ["MEASURE_IT_OFFLINE"] = "1"
    from .mcp.server import serve
    serve(transport=transport, host=host, port=port)


@app.command("serve-api")
def serve_api(
    host: str = typer.Option("127.0.0.1", "--host", help="Interface to bind (0.0.0.0 exposes it on the network)."),
    port: int = typer.Option(8000, "--port", help="Port."),
    offline: bool = typer.Option(False, "--offline", help="MEASURE_IT_OFFLINE=1: answer from data/processed and the "
                                                          "HTTP cache only (search_us_open_data never calls out)."),
    reload: bool = typer.Option(False, "--reload", help="Auto-reload on code changes (development)."),
) -> None:
    """Start the FastAPI app (measure_it.api.app:app): /health, /tools, GET+POST /tools/{tool}, /sources,
    /geojson/{level}; OpenAPI docs at /docs. See docs/DASHBOARD.md."""
    import uvicorn
    if offline:
        os.environ["MEASURE_IT_OFFLINE"] = "1"
    uvicorn.run("measure_it.api.app:app", host=host, port=port, reload=reload)


@app.command("dashboard")
def dashboard(
    port: int = typer.Option(8501, "--port", help="Port."),
    host: str = typer.Option("127.0.0.1", "--host", help="Interface to bind (0.0.0.0 exposes it on the network)."),
    offline: bool = typer.Option(False, "--offline", help="MEASURE_IT_OFFLINE=1 for every tool call the pages make."),
    headless: bool = typer.Option(False, "--headless", help="Do not open a browser (automatic when stdin is not a "
                                                            "terminal: nohup, CI, ssh without a TTY, a service)."),
) -> None:
    """Start the Streamlit dashboard (dashboard/app.py: home + 5 pages calling measure_it.tools). See docs/DASHBOARD.md.

    The process is replaced by Streamlit (exec), so Ctrl-C or `kill <pid>` stops the server itself."""
    import sys as _sys
    from .config import PROJECT_ROOT
    env = {**os.environ, **({"MEASURE_IT_OFFLINE": "1"} if offline else {})}
    # Streamlit's first-run e-mail prompt blocks, and exits without a message when stdin is not a TTY: run headless
    headless = headless or not _sys.stdin.isatty()
    cmd = [_sys.executable, "-m", "streamlit", "run", str(PROJECT_ROOT / "dashboard" / "app.py"),
           "--server.address", host, "--server.port", str(port), "--server.headless", "true" if headless else "false",
           "--browser.gatherUsageStats", "false", "--client.toolbarMode", "viewer"]
    typer.echo(f"measure-it dashboard: http://{'localhost' if host in ('127.0.0.1', 'localhost') else host}:{port} "
               f"(pid {os.getpid()}; stop with Ctrl-C or kill {os.getpid()})", err=True)
    _sys.stdout.flush()
    _sys.stderr.flush()
    os.chdir(PROJECT_ROOT)
    os.execve(_sys.executable, cmd, env)   # no orphaned child: this process becomes the Streamlit server


@app.command("ask")
def ask(
    question: str = typer.Argument(..., help="A deployment question, e.g. 'Given phenotype X and candidate "
                                             "measurement technology Y, where should we deploy Y?'"),
    as_json: bool = typer.Option(False, "--json", help="Print the full JSON result (parse, tool calls, sections)."),
    top_n: Optional[int] = typer.Option(None, "--top-n", min=1, max=25, help="Regions to list (default: from the "
                                                                              "question, else 5)."),
    level: Optional[str] = typer.Option(None, "--level", help="county or state (default: from the question)."),
) -> None:
    """Deterministic, LLM-free answer to 'what should we measure, where, and who could deploy it?' with cited object
    ids (measure_it.agents.deployment_agent)."""
    import contextlib
    import json as _json
    import sys as _sys
    from .agents.deployment_agent import answer_question
    try:
        # progress lines of the layers go to stderr so that stdout is only the answer (parseable with --json)
        with contextlib.redirect_stdout(_sys.stderr):
            res = answer_question(question, top_n=top_n, geography_level=level)
    except FileNotFoundError as e:
        typer.echo(f"measure-it ask: a processed table is missing ({str(e).splitlines()[0][:200]}). Build the data "
                   "first: `uv run measure-it pipeline --workers 16` on a fresh clone (network), or "
                   "`uv run measure-it pipeline --offline --workers 16` when data/raw and data/_http_cache exist "
                   "(README, 'Reproducing from scratch').", err=True)
        _sys.stdout.flush()
        _sys.stderr.flush()
        os._exit(1)
    typer.echo(_json.dumps(res, indent=1, default=str) if as_json else res["answer_markdown"])
    _sys.stdout.flush()
    _sys.stderr.flush()
    os._exit(0)   # skip native-library teardown (measure_it.omics.query.exit_cleanly explains the spurious abort)


@app.command("bundle")
def bundle(
    out: str = typer.Option("dist", "--out", help="Directory for the archive, its manifest and .sha256."),
    no_duckdb: bool = typer.Option(False, "--no-duckdb", help="Leave out data/processed/measure_it_public.duckdb "
                                                              "(the tools fall back to the Parquet tables)."),
    compression: Optional[str] = typer.Option(None, "--compression", help="zst (default when zstd is installed) "
                                                                          "or gz."),
    unpack: Optional[str] = typer.Option(None, "--unpack", help="After building, unpack the archive into this "
                                                                "directory (e.g. dist/bundle, what docker-compose "
                                                                "mounts)."),
    verify: Optional[str] = typer.Option(None, "--verify", help="Do not build: verify an unpacked bundle directory "
                                                                "against its BUNDLE_MANIFEST.json and exit."),
) -> None:
    """Assemble the read-only data bundle (processed tables, DuckDB, provenance sidecars, results tables, configs,
    SOURCE_REGISTRY.yaml) into dist/measure-it-data-<date>.tar.zst with BUNDLE_MANIFEST.json (sha256 and rows per
    file, git commit, source versions). See docs/DEPLOYMENT.md."""
    import json as _json
    from pathlib import Path
    from .export import bundle as B
    if verify:
        rep = B.verify_bundle(Path(verify))
        typer.echo(_json.dumps({k: (v[:20] if isinstance(v, list) else v) for k, v in rep.items()}, indent=1))
        raise typer.Exit(0 if rep["ok"] else 1)
    if compression not in (None, "zst", "gz"):
        raise typer.BadParameter("compression must be zst or gz")
    res = B.build_bundle(Path(out), include_duckdb=not no_duckdb, compression=compression, echo=typer.echo)
    if unpack:
        B.unpack_bundle(Path(res["archive"]), Path(unpack), echo=typer.echo)
    typer.echo(_json.dumps({k: v for k, v in res.items()}, indent=1))


@app.command("export-static")
def export_static(
    out: str = typer.Option("dist/static_snapshot/index.html", "--out", help="Output HTML file."),
) -> None:
    """Write the single-file, read-only reviewer snapshot (overview, primary deployment query, county map, top-10
    evidence cards, measurement performance, what works / what does not, validation tests, limitations). Every number
    is read from data/processed and results/ at build time."""
    from pathlib import Path
    from .export.static_site import build_static_site
    res = build_static_site(Path(out))
    typer.echo(f"{res['path']} ({res['bytes'] / 1e6:.2f} MB)")


if __name__ == "__main__":
    app()
