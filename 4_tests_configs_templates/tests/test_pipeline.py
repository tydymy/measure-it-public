"""Pipeline DAG, offline mode, build-db and validate (unit tests + data tests on the built outputs)."""
from __future__ import annotations

import importlib
import json
import os
import subprocess
import sys

import pytest

from measure_it import http, pipeline as P, validate as V
from measure_it.config import DUCKDB_PATH, PROCESSED, PROJECT_ROOT


# ------------------------------------------------------------------ DAG
def test_steps_are_topological_unique_and_importable():
    names = [s.name for s in P.STEPS]
    assert len(names) == len(set(names))
    seen = set()
    for s in P.STEPS:
        assert set(s.deps) <= seen, f"{s.name} depends on a later step"
        seen.add(s.name)
        mod, fn = s.target.split(":")
        assert callable(getattr(importlib.import_module(mod), fn)), s.target


def test_pipeline_covers_every_stage_named_in_the_brief():
    names = {s.name for s in P.STEPS}
    must = {"geography_crosswalk", "ontology", "omics_coherence", "omics_graph", "wearable_models",
            "wearable_signatures", "wearable_unsupervised", "wearable_nhanes_report", "stanford_models",
            "digital_person", "wearable_adapter", "measurement_adapters", "measurement_discovery",
            "measurement_report", "measurement_phenotype_evidence", "geography_features", "facilities_registry",
            "scoring_report", "canonical_unions", "registry_merge", "build_db", "validate"}
    assert must <= names
    ingest = {s.name for s in P.STEPS if s.name.startswith("ingest_")}
    mods = {p.stem for p in (PROJECT_ROOT / "src/measure_it/ingestion").glob("*.py") if p.stem != "__init__"}
    assert {f"ingest_{m}" for m in mods} <= ingest            # every ingestion module is a step
    assert P.STEPS[-1].name == "validate"


def test_select_steps_only_and_from():
    assert [s.name for s in P.select_steps(only=["validate", "ontology"])] == ["ontology", "validate"]
    names = [s.name for s in P.STEPS]
    tail = [s.name for s in P.select_steps(from_step="canonical_unions")]
    assert tail == names[names.index("canonical_unions"):]           # --from = that step and every later one
    # the real order after the unions: Figure 2 (figures_schematics) reads the SOURCE_REGISTRY.yaml that
    # registry_merge writes; figures_data runs before the unions (right after scoring_report)
    assert tail == ["canonical_unions", "registry_merge", "figures_schematics", "build_db", "validate"]
    assert names.index("scoring_report") < names.index("figures_data") < names.index("canonical_unions")
    with pytest.raises(SystemExit):
        P.select_steps(only=["no_such_step"])


def test_step_workers_capped_by_budget():
    s = P.STEP_BY_NAME["stanford_models"]
    assert s.workers(16) == s.default_workers and s.workers(2) == 2 and s.workers(0) == 1


def test_union_step_never_unions_module_owned_canonicals():
    for t in ("research_site_registry", "facilities", "geo_condition_burden", "geo_context", "participants"):
        assert t in P.OWNED_CANONICALS


# ------------------------------------------------------------------ offline mode
def test_offline_http_raises_on_cache_miss_and_serves_hits(tmp_path, monkeypatch):
    monkeypatch.setattr(http, "HTTP_CACHE", tmp_path)
    monkeypatch.setenv(http.OFFLINE_ENV, "1")
    url = "https://example.invalid/api"
    with pytest.raises(http.OfflineCacheMiss):
        http.get(url, params={"q": 1})
    key = http._cache_key("GET", url, {"q": 1}, None, None)
    meta, body = http._cache_paths(key)
    meta.parent.mkdir(parents=True)
    body.write_bytes(b'{"ok": true}')
    meta.write_text(json.dumps({"url": url, "status": 200, "headers": {}, "fetched_at": "t"}))
    r = http.get(url, params={"q": 1}, refresh=True)     # a refresh cannot happen offline: served from cache
    assert r.from_cache and r.json() == {"ok": True}


def test_offline_download_requires_manifested_raw_file(tmp_path, monkeypatch):
    from measure_it import download as D
    monkeypatch.setattr(D, "RAW", tmp_path)
    monkeypatch.setattr(D, "raw_dir", lambda sid: (tmp_path / sid).mkdir(parents=True, exist_ok=True) or tmp_path / sid)
    monkeypatch.setenv(http.OFFLINE_ENV, "1")
    with pytest.raises(http.OfflineCacheMiss):
        D.download_file("https://example.invalid/f.csv", "src_x", "f.csv")


def test_socket_guard_blocks_network_in_a_subprocess():
    code = ("import os, socket; os.environ['MEASURE_IT_OFFLINE']='1'\n"
            "import measure_it\n"
            "from measure_it.http import OfflineNetworkBlocked\n"
            "try:\n    socket.create_connection(('example.com', 80), timeout=2)\n    print('CONNECTED')\n"
            "except OfflineNetworkBlocked:\n    print('BLOCKED')\n")
    env = dict(os.environ, MEASURE_IT_OFFLINE="1")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, timeout=60)
    assert out.stdout.strip() == "BLOCKED", out.stdout + out.stderr


# ------------------------------------------------------------------ validate helpers
def test_forbidden_language_allows_quoted_rule_text():
    assert V.forbidden_hits('Avoid "proves", "cures" and "best clinic".') == []
    assert V.forbidden_hits("Never say `definitive biomarker`.") == []
    assert V.forbidden_hits("The terms `patient has` and `optimal treatment` are banned.") == []
    hits = V.forbidden_hits("This analysis proves the wearable works.")
    assert [h[1] for h in hits] == ["proves"]
    assert V.forbidden_hits("Exercise improves fitness; the county secures funds.") == []
    # a rule list excuses its bare items, as in SPEC.md ("Avoid (unless directly supported): diagnosed by AI; ...")
    assert V.forbidden_hits("Avoid (unless directly supported): diagnosed by AI; proves; cures; best clinic.") == []
    assert V.forbidden_hits("Never use: proves, cures and best clinic.") == []


def test_forbidden_language_rule_word_elsewhere_does_not_excuse_a_use():
    assert [h[1] for h in V.forbidden_hits("Wearables never prove anything, but this proves it.")] == ["proves"]
    assert [h[1] for h in V.forbidden_hits("We found the optimal treatment; we avoid overclaiming.")] == \
        ["optimal treatment"]
    assert [h[1] for h in V.forbidden_hits("Avoid claims: this proves nothing.")] == ["proves"]


def test_union_dtype_harmonisation_keeps_numbers_numeric():
    import numpy as np
    import pandas as pd
    assert P._common_dtype([np.dtype("int64"), pd.Int64Dtype()]) == "Int64"
    assert P._common_dtype([np.dtype("int64"), np.dtype("float64")]) == "Float64"
    assert P._common_dtype([np.dtype("bool"), pd.BooleanDtype()]) == "boolean"
    assert P._common_dtype([np.dtype("int64"), pd.StringDtype()]) == "string"


def test_geo_column_regex():
    for c in ("county_fips", "state_fips", "zip5", "zcta", "lat", "lon", "latitude", "geo_id", "county_fips_zip"):
        assert V.GEO_COLUMN_RE.search(c), c
    for c in ("participant_id", "sleep_proxy_onset_h", "population", "relative_amplitude", "translation"):
        assert not V.GEO_COLUMN_RE.search(c), c


# ------------------------------------------------------------------ built outputs
def test_failed_step_reports_its_error_line_log_path_and_hint(tmp_path):
    log = tmp_path / "geography_crosswalk.log"
    log.write_text("[step] geography_crosswalk: measure_it.geography.crosswalk:run() offline=True at x\n"
                   "Traceback (most recent call last):\n  File \"a.py\", line 3, in <module>\n"
                   "measure_it.http.OfflineCacheMiss: raw file data/raw/x.zip needs the network\n"
                   "[step] geography_crosswalk: FAILED after 1.0s\n")
    err = P.last_error_line(log)
    assert err.startswith("measure_it.http.OfflineCacheMiss") and "needs the network" in err
    assert "without --offline" in P.failure_hint(err, offline=True)
    assert P.last_error_line(tmp_path / "missing.log") == ""


def test_run_log_paths_are_relative_to_the_checkout():
    assert P._rel(PROJECT_ROOT / "results" / "pipeline_runs" / "x.log") == "results/pipeline_runs/x.log"
    assert P._argv_rel([str(PROJECT_ROOT / ".venv" / "bin" / "measure-it"), "pipeline"]) == [
        ".venv/bin/measure-it", "pipeline"]


def test_reap_reports_return_code_and_peak_rss():
    import time
    proc = subprocess.Popen([sys.executable, "-c", "b = bytearray(200 * 1024 * 1024); raise SystemExit(3)"])
    while True:
        rc, rss = P._reap(proc)
        if rc is not None:
            break
        time.sleep(0.05)
    assert rc == 3
    if hasattr(os, "wait4"):
        assert rss is not None and rss >= 150


def test_checks_with_nothing_to_check_fail_not_pass(tmp_path, monkeypatch):
    """An empty checkout must not PASS the provenance, person-layer or object-id checks on 0/0 tables."""
    monkeypatch.setattr(V, "PROCESSED", tmp_path)
    for check in (V.check_provenance_all, V.check_person_geography, V.check_object_ids):
        r = check()
        assert r["ok"] is False and r["details"], check.__name__


def test_build_db_exits_non_zero_when_spec_views_are_missing(tmp_path, monkeypatch):
    from measure_it import database as D
    monkeypatch.setattr(D, "PROCESSED", tmp_path)
    with pytest.raises(SystemExit) as e:
        D.build_db(path=tmp_path / "empty.duckdb")
    assert "SPEC views have no underlying table" in str(e.value.code)


@pytest.mark.data
def test_duckdb_has_tables_spec_views_and_sources():
    if not DUCKDB_PATH.exists():
        pytest.skip("run `measure-it build-db` first")
    import duckdb
    from measure_it.database import SPEC_VIEWS
    con = duckdb.connect(str(DUCKDB_PATH), read_only=True)
    views = {r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_schema = 'spec' AND table_type = 'VIEW'").fetchall()}
    assert set(SPEC_VIEWS) <= views
    assert con.execute("SELECT count(*) FROM main.conditions").fetchone()[0] == \
        con.execute("SELECT count(*) FROM condition_registry").fetchone()[0]
    assert con.execute("SELECT count(*) FROM main.measurements").fetchone()[0] > 0
    n_src = con.execute("SELECT count(*) FROM _sources").fetchone()[0]
    assert n_src >= 25
    tables = {r[0] for r in con.execute("SELECT name FROM _catalog WHERE kind = 'table'").fetchall()}
    assert {p.stem for p in PROCESSED.glob("*.parquet")} <= tables
    con.close()


@pytest.mark.data
def test_validate_checks_pass_on_the_build():
    for check in (V.check_provenance_all, V.check_person_geography, V.check_object_ids, V.check_language,
                  V.check_audits):
        r = check()
        assert r["ok"], (r["check"], r["summary"], r["details"][:5])


@pytest.mark.data
def test_no_research_site_registry_partition_remains():
    assert not list(PROCESSED.glob("research_site_registry__*.parquet"))
    assert (PROCESSED / "ctgov_facility_summary.parquet").exists()


@pytest.mark.data
def test_ask_json_stdout_is_parseable_json():
    """Progress lines of the scoring layer go to stderr, so `measure-it ask --json` prints JSON only."""
    env = dict(os.environ, MEASURE_IT_OFFLINE="1")
    out = subprocess.run([sys.executable, "-m", "measure_it.cli", "ask", "--json",
                          "Where should we deploy autonomic testing for POTS?"],
                         cwd=PROJECT_ROOT, env=env, capture_output=True, text=True, timeout=900)
    assert out.returncode == 0, out.stderr[-2000:]
    res = json.loads(out.stdout)
    assert res["answer_markdown"] and res["cited_object_ids"]
