"""The analytical DuckDB: every processed table, SPEC core-data-model views, and the source registry.

    uv run measure-it build-db        # -> data/processed/measure_it_public.duckdb

Contents
* one TABLE per data/processed/*.parquet (same name), plus the two boundary geoparquets
  (county_boundaries_2024, state_boundaries_2024; geometry as WKB BLOB);
* VIEWS for the SPEC "CORE DATA MODEL" names. Names that differ from a table name get a view in the main schema
  (``conditions`` -> condition_registry, ``measurements`` -> measurement_registry); every SPEC name also has a view in
  schema ``spec`` (``spec.participants``, ``spec.conditions``, ...), so the SPEC model can be queried uniformly;
* ``_sources``: one row per SOURCE_REGISTRY.yaml source (lists/dicts as JSON text);
* ``_catalog``: table/view name, kind, rows, producer, created_at (from the .meta.json sidecars);
* ``_spec_views``: SPEC name -> underlying table, with the table's row count.

The DB is rebuilt atomically (written to a temporary file, then renamed). Everything here stays ecological or
person-level exactly as in the parquet files; views add no joins across data layers.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import duckdb
import yaml

from .config import DUCKDB_PATH, PROCESSED, SOURCE_REGISTRY_PATH, utc_now_iso

# SPEC "CORE DATA MODEL" name -> processed table
SPEC_VIEWS: dict[str, str] = {
    "participants": "participants",
    "participant_conditions": "participant_conditions",
    "participant_labs": "participant_labs",
    "participant_wearable_features": "participant_wearable_features",
    "participant_phenotype_embeddings": "participant_phenotype_embeddings",
    "conditions": "condition_registry",
    "condition_ontology_mappings": "condition_ontology_mappings",
    "condition_molecular_evidence": "condition_molecular_evidence",
    "measurements": "measurement_registry",
    "condition_measurement_evidence": "condition_measurement_evidence",
    "measurement_regulatory_status": "measurement_regulatory_status",
    "geographies": "geographies",
    "geo_condition_burden": "geo_condition_burden",
    "geo_context": "geo_context",
    "geo_vulnerability": "geo_vulnerability",
    "providers": "providers",
    "facilities": "facilities",
    "clinical_trials": "clinical_trials",
    "trial_sites": "trial_sites",
    "nih_projects": "nih_projects",
    "deployment_candidates": "deployment_candidates",
}
GEO_FILES = ("county_boundaries_2024.geoparquet", "state_boundaries_2024.geoparquet")


def _q(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _meta(stem: str) -> dict:
    p = PROCESSED / f"{stem}.meta.json"
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return {}


def _sources_rows(path: Path = SOURCE_REGISTRY_PATH) -> tuple[list[str], list[list]]:
    reg = yaml.safe_load(path.read_text()) if path.exists() else {"sources": []}
    sources = reg.get("sources", []) or []
    cols: list[str] = []
    for s in sources:
        for k in s:
            if k not in cols:
                cols.append(k)
    rows = []
    for s in sources:
        row = []
        for c in cols:
            v = s.get(c)
            if isinstance(v, (list, dict)):
                v = json.dumps(v, default=str)
            elif v is not None and not isinstance(v, str):
                v = str(v)
            row.append(v)
        rows.append(row)
    return cols, rows


def build_db(path: Path = DUCKDB_PATH, threads: int | None = None) -> Path:
    """(Re)build the DuckDB catalogue; returns its path."""
    tmp = path.with_suffix(".building.duckdb")
    tmp.unlink(missing_ok=True)
    Path(str(tmp) + ".wal").unlink(missing_ok=True)
    con = duckdb.connect(str(tmp))
    if threads:
        con.execute(f"SET threads={int(threads)}")
    catalog = []
    for pq in sorted(PROCESSED.glob("*.parquet")):
        name = pq.stem
        con.execute(f"CREATE TABLE {_q(name)} AS SELECT * FROM read_parquet(?)", [str(pq)])
        n = con.execute(f"SELECT count(*) FROM {_q(name)}").fetchone()[0]
        m = _meta(name)
        catalog.append((name, "table", int(n), m.get("producer"), m.get("created_at"), str(pq.relative_to(PROCESSED.parent.parent))))
    for fn in GEO_FILES:
        p = PROCESSED / fn
        if p.exists():
            name = fn.split(".")[0]
            con.execute(f"CREATE TABLE {_q(name)} AS SELECT * FROM read_parquet(?)", [str(p)])
            n = con.execute(f"SELECT count(*) FROM {_q(name)}").fetchone()[0]
            catalog.append((name, "table (geoparquet; geometry as WKB)", int(n), "measure_it.geography.crosswalk", None,
                            str(p.relative_to(PROCESSED.parent.parent))))
    tables = {c[0] for c in catalog}
    # SPEC views
    con.execute("CREATE SCHEMA spec")
    spec_rows = []
    for view, table in SPEC_VIEWS.items():
        if table not in tables:
            spec_rows.append((view, table, None, "missing: underlying table not built"))
            continue
        con.execute(f"CREATE VIEW spec.{_q(view)} AS SELECT * FROM main.{_q(table)}")
        if view != table:
            con.execute(f"CREATE VIEW main.{_q(view)} AS SELECT * FROM main.{_q(table)}")
            catalog.append((view, f"view -> {table}", next(c[2] for c in catalog if c[0] == table), None, None, None))
        n = next(c[2] for c in catalog if c[0] == table)
        spec_rows.append((view, table, n, "ok"))
    con.execute("CREATE TABLE _spec_views (spec_name VARCHAR, table_name VARCHAR, n_rows BIGINT, status VARCHAR)")
    con.executemany("INSERT INTO _spec_views VALUES (?, ?, ?, ?)", spec_rows)
    # sources
    cols, rows = _sources_rows()
    if cols:
        con.execute("CREATE TABLE _sources (" + ", ".join(f"{_q(c)} VARCHAR" for c in cols) + ")")
        con.executemany("INSERT INTO _sources VALUES (" + ", ".join("?" for _ in cols) + ")", rows)
    # catalogue
    con.execute("CREATE TABLE _catalog (name VARCHAR, kind VARCHAR, n_rows BIGINT, producer VARCHAR, "
                "created_at VARCHAR, path VARCHAR)")
    if catalog:   # an empty checkout has no table (build-db then exits non-zero below)
        con.executemany("INSERT INTO _catalog VALUES (?, ?, ?, ?, ?, ?)", catalog)
    con.execute("CREATE TABLE _build (built_at VARCHAR, n_tables BIGINT, n_sources BIGINT)")
    con.execute("INSERT INTO _build VALUES (?, ?, ?)", [utc_now_iso(), len(tables), len(rows)])
    con.close()
    n_ok = sum(1 for r in spec_rows if r[3] == "ok")
    if n_ok < len(SPEC_VIEWS):
        # an incomplete catalogue is not installed (validate would count it as present); any earlier DuckDB is kept
        tmp.unlink(missing_ok=True)
        Path(str(tmp) + ".wal").unlink(missing_ok=True)
        missing = [f"{r[0]} -> {r[1]}" for r in spec_rows if r[3] != "ok"]
        raise SystemExit(f"[build-db] FAILED: {len(tables)} tables; {len(missing)} of {len(SPEC_VIEWS)} SPEC views have "
                         f"no underlying table ({', '.join(missing[:6])}{' ...' if len(missing) > 6 else ''}); "
                         f"{path.name} not written. Run the pipeline first")
    os.replace(tmp, path)
    Path(str(tmp) + ".wal").unlink(missing_ok=True)
    print(f"[build-db] {path}: {len(tables)} tables, {len(SPEC_VIEWS)} SPEC views ({n_ok} ok), {len(rows)} sources",
          flush=True)
    return path
