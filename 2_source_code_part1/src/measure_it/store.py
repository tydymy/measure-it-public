"""Processed-table storage: Parquet files under data/processed + a DuckDB catalog.

Naming convention
-----------------
Canonical tables are data/processed/<table>.parquet.
Source-specific partitions are data/processed/<table>__<source>.parquet
(e.g. participant_wearable_features__nhanes). `union_partitions(table)` stacks
the partitions into the canonical table. Each file gets a sidecar
<name>.meta.json with row count, columns, producer and creation time.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import duckdb
import pandas as pd

from .config import DUCKDB_PATH, PROCESSED, utc_now_iso
from .provenance import check_provenance


def processed_path(name: str) -> Path:
    return PROCESSED / f"{name}.parquet"


def write_table(df: pd.DataFrame, name: str, *, producer: str, require_provenance: bool = True,
                description: str = "") -> Path:
    if require_provenance:
        check_provenance(df, name)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    path = processed_path(name)
    fd, tmp = tempfile.mkstemp(dir=PROCESSED, prefix=f".{name}-", suffix=".parquet")
    os.close(fd)
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)
    os.chmod(path, 0o644)  # mkstemp creates 0600; processed tables must be readable by the dashboard/API
    meta = {
        "table": name,
        "rows": int(len(df)),
        "columns": {c: str(t) for c, t in df.dtypes.items()},
        "producer": producer,
        "created_at": utc_now_iso(),
        "description": description,
    }
    path.with_suffix(".meta.json").write_text(json.dumps(meta, indent=2))
    return path


def read_table(name: str, columns: list[str] | None = None) -> pd.DataFrame:
    path = processed_path(name)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; run the pipeline step that produces {name!r}")
    return pd.read_parquet(path, columns=columns)


def table_exists(name: str) -> bool:
    return processed_path(name).exists()


def partitions(table: str) -> list[Path]:
    return sorted(PROCESSED.glob(f"{table}__*.parquet"))


def union_partitions(table: str, *, producer: str = "store.union_partitions") -> Path:
    parts = partitions(table)
    if not parts:
        raise FileNotFoundError(f"no partitions for {table}")
    frames = [pd.read_parquet(p) for p in parts]
    df = pd.concat(frames, ignore_index=True, sort=False)
    return write_table(df, table, producer=producer,
                       description=f"union of {[p.stem for p in parts]}")


def build_duckdb(path: Path = DUCKDB_PATH) -> Path:
    """(Re)create the DuckDB catalog: delegates to measure_it.database.build_db (tables + SPEC views + _sources)."""
    from .database import build_db
    return build_db(path)


def _build_duckdb_tables_only(path: Path = DUCKDB_PATH) -> Path:
    """The original builder: one table per processed parquet file plus _catalog (kept for reference)."""
    tmp = path.with_suffix(".building.duckdb")
    tmp.unlink(missing_ok=True)
    con = duckdb.connect(str(tmp))
    for pq in sorted(PROCESSED.glob("*.parquet")):
        name = pq.stem
        con.execute(f'CREATE TABLE "{name}" AS SELECT * FROM read_parquet(?)', [str(pq)])
    con.execute("CREATE TABLE _catalog AS SELECT * FROM (VALUES (NULL::VARCHAR, NULL::BIGINT)) t(table_name, n_rows) WHERE false")
    for pq in sorted(PROCESSED.glob("*.parquet")):
        n = con.execute(f'SELECT count(*) FROM "{pq.stem}"').fetchone()[0]
        con.execute("INSERT INTO _catalog VALUES (?, ?)", [pq.stem, n])
    con.close()
    os.replace(tmp, path)
    return path


def connect(read_only: bool = True) -> duckdb.DuckDBPyConnection:
    return duckdb.connect(str(DUCKDB_PATH), read_only=read_only)
