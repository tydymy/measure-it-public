"""Compare a pipeline re-run with a baseline copy of data/processed and results/tables (determinism check).

    uv run python -m measure_it.reproduce compare --baseline data/interim/repro_backup/pre_integration \\
        --since 2026-09-24T05:00:00+00:00 --out results/tables/reproducibility_comparison.csv

For every processed table (*.parquet) and result table (results/tables/*.csv, *.json) in either tree:
row counts, column sets, and an order-independent content hash of the rows (a sum of per-row hashes in DuckDB, i.e.
the multiset of rows, equivalent to hashing the rows sorted by all columns).

Columns that legitimately change per run are those that carry the run's own timestamp (build/validation time
written by utc_now_iso(): created_at, built_at, retrieved_at of a derived table, "audit built <time>" in a
source_version, ...). They are detected, not assumed: a column is `run_stamped` when a value in the NEW table holds
an ISO-8601 timestamp inside the run window [--since, --until] (data dates outside the window, such as the shifted
future dates of a study release, never qualify). Such columns are hashed after replacing every ISO
timestamp with a placeholder, so a difference that is only the timestamp does not count, and any other difference
in the same column still does. When the hashes differ, each column is compared as a multiset to name the columns
that changed, and numeric columns are re-compared with a relative tolerance of 1e-9 (float noise).
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from .config import PROCESSED, PROJECT_ROOT, TABLES, utc_now_iso

TS_SQL = r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(\.\d+)?(\+00:00|Z|\+0000)?"
TS_RE = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:\+00:00|Z|\+0000)?")


def _q(c: str) -> str:
    return '"' + c.replace('"', '""') + '"'


def _reader(path: Path) -> str:
    p = str(path).replace("'", "''")
    if path.suffix == ".csv":
        return f"read_csv('{p}', all_varchar=true, header=true, sample_size=-1)"
    return f"read_parquet('{p}')"


def _columns(con, path: Path) -> list[tuple[str, str]]:
    return [(r[0], r[1]) for r in con.execute(f"DESCRIBE SELECT * FROM {_reader(path)}").fetchall()]


def _run_stamped(con, path: Path, cols: list[tuple[str, str]], since: str, until: str | None = None) -> list[str]:
    """Columns holding a timestamp inside the run window [since, until] (data dates outside it do not count)."""
    until = (until or "9999-12-31T23:59:59")[:19]
    out = []
    for c, t in cols:
        if not any(k in t.upper() for k in ("VARCHAR", "TIMESTAMP", "DATE")):
            continue
        n = con.execute(
            f"SELECT count(*) FROM (SELECT regexp_extract_all(CAST({_q(c)} AS VARCHAR), '{TS_SQL}') AS ts "
            f"FROM {_reader(path)}) WHERE len(list_filter(ts, x -> replace(x, ' ', 'T') >= '{since[:19]}' "
            f"AND replace(x, ' ', 'T') <= '{until}')) > 0"
        ).fetchone()[0]
        if n:
            out.append(c)
    return out


def _expr(c: str, t: str, normalize: bool) -> str:
    if normalize:
        return f"regexp_replace(CAST({_q(c)} AS VARCHAR), '{TS_SQL}', '<TS>', 'g')"
    return _q(c)


def _table_hash(con, path: Path, cols: list[tuple[str, str]], norm: set[str]) -> tuple[int, int]:
    if not cols:
        return 0, 0
    parts = ", ".join(_expr(c, t, c in norm) for c, t in cols)
    n, h = con.execute(f"SELECT count(*), sum(hash({parts})::HUGEINT) FROM {_reader(path)}").fetchone()
    return int(n), int(h or 0)


def _column_hashes(con, path: Path, cols: list[tuple[str, str]], norm: set[str]) -> dict[str, int]:
    if not cols:
        return {}
    sel = ", ".join(f"sum(hash({_expr(c, t, c in norm)})::HUGEINT)" for c, t in cols)
    row = con.execute(f"SELECT {sel} FROM {_reader(path)}").fetchone()
    return {c: int(v or 0) for (c, _), v in zip(cols, row)}


def _numeric_close(con, old: Path, new: Path, col: str) -> bool | None:
    """True when the column's sorted numeric values agree within rtol 1e-9 (None if not numeric)."""
    try:
        a = con.execute(f"SELECT TRY_CAST({_q(col)} AS DOUBLE) AS v FROM {_reader(old)} ORDER BY v NULLS LAST").df()["v"]
        b = con.execute(f"SELECT TRY_CAST({_q(col)} AS DOUBLE) AS v FROM {_reader(new)} ORDER BY v NULLS LAST").df()["v"]
    except Exception:  # noqa: BLE001
        return None
    if len(a) != len(b) or a.notna().sum() == 0 or a.notna().sum() != b.notna().sum():
        return None
    return bool(np.allclose(a.to_numpy(float), b.to_numpy(float), rtol=1e-9, atol=1e-12, equal_nan=True))


def _keyed_diff(con, old: Path, new: Path, common, tco: dict, stamped: set[str]) -> dict:
    """Row-level diff on a unique key (source_record_id, else object_id): rows only in old / only in new / changed."""
    names = [c for c, _ in common]
    for key in ("source_record_id", "object_id"):
        if key not in names:
            continue
        ok = all(con.execute(f"SELECT count(*) = count(DISTINCT {_q(key)}) FROM {_reader(p)}").fetchone()[0]
                 for p in (old, new))
        if not ok:
            continue
        rest = [(c, t) for c, t in common if c != key]
        ho = ", ".join(_expr(c, tco[c], c in stamped) for c, _ in rest) or "1"
        hn = ", ".join(_expr(c, t, c in stamped) for c, t in rest) or "1"
        r = con.execute(
            f"WITH a AS (SELECT CAST({_q(key)} AS VARCHAR) k, hash({ho}) h FROM {_reader(old)}), "
            f"b AS (SELECT CAST({_q(key)} AS VARCHAR) k, hash({hn}) h FROM {_reader(new)}) "
            "SELECT count(*) FILTER (WHERE b.k IS NULL), count(*) FILTER (WHERE a.k IS NULL), "
            "count(*) FILTER (WHERE a.k IS NOT NULL AND b.k IS NOT NULL AND a.h <> b.h) "
            "FROM a FULL OUTER JOIN b ON a.k = b.k").fetchone()
        return {"diff_key": key, "rows_only_old": int(r[0]), "rows_only_new": int(r[1]), "rows_changed": int(r[2])}
    return {}


def compare_file(con, old: Path | None, new: Path | None, since: str, until: str | None = None) -> dict:
    rec: dict = {"table": (new or old).name, "kind": "processed" if (new or old).suffix == ".parquet" else "result"}
    if old is None or not old.exists():
        rec.update(status="new_in_run")
        return rec
    if new is None or not new.exists():
        rec.update(status="missing_in_run")
        return rec
    if new.suffix == ".json":
        a, b = json.loads(old.read_text()), json.loads(new.read_text())
        na, nb = TS_RE.sub("<TS>", json.dumps(a, sort_keys=True)), TS_RE.sub("<TS>", json.dumps(b, sort_keys=True))
        if json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True):
            rec.update(status="identical")
        elif na == nb:
            rec.update(status="identical", run_stamped_columns="(timestamps only)")
        else:
            ka = set(a) if isinstance(a, dict) else set()
            diff = sorted(k for k in ka | (set(b) if isinstance(b, dict) else set())
                          if TS_RE.sub("<TS>", json.dumps(a.get(k), sort_keys=True, default=str)) !=
                          TS_RE.sub("<TS>", json.dumps(b.get(k), sort_keys=True, default=str))) if isinstance(a, dict) else []
            rec.update(status="changed", changed_columns=";".join(diff) or "(content)")
        return rec
    co, cn = _columns(con, old), _columns(con, new)
    so, sn = {c for c, _ in co}, {c for c, _ in cn}
    stamped = set(_run_stamped(con, new, cn, since, until))
    common = [(c, t) for c, t in cn if c in so]
    tco = dict(co)
    type_changed = sorted(c for c, t in common if tco[c] != t)
    no, ho = _table_hash(con, old, [(c, tco[c]) for c, _ in co], stamped)
    nn, hn = _table_hash(con, new, cn, stamped)
    rec.update(rows_old=no, rows_new=nn, n_cols_old=len(co), n_cols_new=len(cn),
               cols_added=";".join(sorted(sn - so)), cols_removed=";".join(sorted(so - sn)),
               dtype_changed=";".join(type_changed), run_stamped_columns=";".join(sorted(stamped)))
    if so == sn and not type_changed and no == nn and ho == hn:
        rec["status"] = "identical"
        return rec
    # column-level multiset comparison on the common columns
    ch_o = _column_hashes(con, old, [(c, tco[c]) for c, _ in common], stamped)
    ch_n = _column_hashes(con, new, common, stamped)
    changed = [c for c, _ in common if ch_o[c] != ch_n[c]]
    close = [c for c in changed if no == nn and _numeric_close(con, old, new, c)]
    rec["changed_columns"] = ";".join(changed)
    rec["float_noise_columns"] = ";".join(close)
    rec.update(_keyed_diff(con, old, new, common, tco, stamped))
    if so == sn and no == nn and changed and set(changed) == set(close):
        rec["status"] = "numerically_equal"
    else:
        rec["status"] = "changed"
    return rec


def compare(baseline: Path, since: str, until: str | None = None, processed: Path = PROCESSED,
            tables: Path = TABLES) -> pd.DataFrame:
    con = duckdb.connect()
    rows = []
    pairs = []
    for sub, cur, pat in (("processed", processed, "*.parquet"), ("tables", tables, "*.csv"),
                          ("tables", tables, "*.json")):
        names = {p.name for p in (baseline / sub).glob(pat)} | {p.name for p in cur.glob(pat)}
        for n in sorted(names):
            pairs.append((baseline / sub / n, cur / n))
    for old, new in pairs:
        try:
            rows.append(compare_file(con, old if old.exists() else None, new if new.exists() else None, since, until))
        except Exception as exc:  # noqa: BLE001
            rows.append({"table": new.name, "status": "error", "error": f"{type(exc).__name__}: {exc}"[:300]})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--baseline", type=Path, required=True, help="directory holding processed/ and tables/")
    c.add_argument("--since", required=True, help="run start (UTC ISO); timestamps in [since, until] mark run-stamped columns")
    c.add_argument("--until", default=None, help="run end (UTC ISO; default: now). Data dates outside the window "
                                                  "(e.g. shifted future study dates) never mark a column")
    c.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    base = a.baseline if a.baseline.is_absolute() else PROJECT_ROOT / a.baseline
    df = compare(base, a.since, a.until or utc_now_iso())
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(a.out, index=False)
    print(df["status"].value_counts().to_string())
    cols = [c for c in ("table", "status", "rows_old", "rows_new", "changed_columns", "cols_added", "cols_removed")
            if c in df]
    print(df[df["status"] != "identical"][cols].to_string(index=False, max_colwidth=80))


if __name__ == "__main__":
    main()
