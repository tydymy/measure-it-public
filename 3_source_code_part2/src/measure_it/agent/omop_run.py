"""Part B: execute the stage-D All of Us query pack (measure_it.similar.omop, BigQuery dialect) against the local
OMOP export, transpiled to DuckDB with sqlglot (docs/AGENT_CONTRACT.md §3).

    check_query_pack(omop_dir, phenotype_json) -> dict      also written to <omop_dir>/query_pack_check.json

sqlglot is NOT a project dependency: it is imported lazily (`uv run --with sqlglot ...`); without it the check
returns status `sqlglot_missing` and says how to run it.

What the check provides that a real All of Us CDR has and this export does not:
  * `concept_relationship`: one 'Maps to' row per local concept, mapping it to itself (the export's concept ids are
    local, so the "standard" target of a source concept is the concept itself);
  * `concept_ancestor`: one self-ancestry row per local concept (min/max levels 0), so an RxNorm ingredient's
    "descendants" are the ingredient concept itself (the export records ingredients);
  * the export's `year_of_birth` is null (only age bands exist), so every person fails the pack's adult filter; the
    check runs the pack twice: `as_exported`, and `year_of_birth_from_age_band` with a CHECK-ONLY year of birth
    = current year - floor(age-band midpoint) (40-49 -> age 44 -> age_mid 44.5 as the pack computes it; 90+ -> 92).
Per-step row counts: each CTE of the pack is counted on its own (WITH <ctes up to it> SELECT COUNT(*) FROM it); the
final aggregate rows are returned as the pack returns them (All of Us small-cell suppression applied by the SQL).
The counts are local, person-level-derived numbers about a private dataset: the JSON stays in the local-only run
folder.
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import pandas as pd

CDR_PLACEHOLDER = "`{CDR}`."
STUB_TABLES = ("concept_relationship", "concept_ancestor")


def _load_phenotype(phenotype_json) -> dict:
    from ..similar.phenotype_io import VOCAB_DOMAIN, load_phenotype
    if isinstance(phenotype_json, dict):
        tmp = json.loads(json.dumps(phenotype_json))
        for f in tmp.get("features", []):
            f.setdefault("vocabulary", f["feature_key"].split(":", 1)[0])
            f.setdefault("domain", VOCAB_DOMAIN.get(f["vocabulary"], "other"))
            f.setdefault("codes", [f["feature_key"].split(":", 1)[1]])
            f.setdefault("transform", "presence")
            f.setdefault("label", f["feature_key"])
        return tmp
    return load_phenotype(Path(phenotype_json))


def _write(omop_dir: Path, result: dict) -> dict:
    omop_dir.mkdir(parents=True, exist_ok=True)
    (omop_dir / "query_pack_check.json").write_text(json.dumps(result, indent=1, default=str))
    return result


def transpile_pack(sql_bq: str) -> tuple[str, list[str], list[tuple[str, str]]]:
    """BigQuery pack -> (DuckDB SQL, CDM tables referenced, [(cte name, DuckDB SQL of the CTE)])."""
    import sqlglot
    from sqlglot import exp
    ast = sqlglot.parse_one(sql_bq.replace(CDR_PLACEHOLDER, "`cdr`."), read="bigquery",
                            error_level=sqlglot.ErrorLevel.RAISE)
    cdm = set()
    for t in ast.find_all(exp.Table):
        if t.args.get("db") is not None and t.db == "cdr":          # map `{CDR}`.x -> the local table x
            cdm.add(t.name)
            t.set("db", None)
            t.set("catalog", None)
    duck = ast.sql(dialect="duckdb", pretty=True)
    with_ = ast.args.get("with") or ast.args.get("with_")
    ctes = [(c.alias, c.sql(dialect="duckdb")) for c in (with_.expressions if with_ else [])]
    return duck, sorted(cdm), ctes


def _connect(omop_dir: Path, *, yob_from_band: bool):
    import duckdb
    con = duckdb.connect(":memory:")
    present = []
    for p in sorted(omop_dir.glob("*.parquet")):
        name = p.stem
        if name == "person" and yob_from_band:
            from ..harmonize.schema import age_mid
            per = pd.read_parquet(p)
            mid = age_mid(per["age_band_source_value"].astype(object))
            per["year_of_birth"] = (dt.date.today().year - mid.apply(lambda v: int(v) if pd.notna(v) else None)
                                    ).astype("Int64")
            con.register("person_df", per)
            con.execute("CREATE TABLE person AS SELECT * FROM person_df")
            con.unregister("person_df")
        else:
            con.execute(f"CREATE TABLE {name} AS SELECT * FROM read_parquet(?)", [str(p)])
        present.append(name)
    if "concept" in present:
        con.execute("CREATE TABLE concept_relationship AS SELECT concept_id AS concept_id_1, concept_id AS "
                    "concept_id_2, 'Maps to' AS relationship_id, DATE '1970-01-01' AS valid_start_date, "
                    "DATE '2099-12-31' AS valid_end_date, CAST(NULL AS VARCHAR) AS invalid_reason FROM concept")
        con.execute("CREATE TABLE concept_ancestor AS SELECT concept_id AS ancestor_concept_id, concept_id AS "
                    "descendant_concept_id, 0 AS min_levels_of_separation, 0 AS max_levels_of_separation FROM concept")
    return con, present


def _run_variant(omop_dir: Path, duck_sql: str, ctes: list[tuple[str, str]], *, yob_from_band: bool) -> dict:
    con, present = _connect(omop_dir, yob_from_band=yob_from_band)
    steps, first_error = [], None
    try:
        for i, (name, _) in enumerate(ctes):
            q = "WITH " + ",\n".join(c for _, c in ctes[: i + 1]) + f"\nSELECT COUNT(*) FROM {name}"
            try:
                n = con.execute(q).fetchone()[0]
                steps.append({"step": name, "rows": int(n)})
            except Exception as exc:  # noqa: BLE001 - reported per step
                msg = f"{type(exc).__name__}: {str(exc).splitlines()[0][:400]}"
                steps.append({"step": name, "rows": None, "error": msg})
                first_error = first_error or {"step": name, "error": msg}
        try:
            res = con.execute(duck_sql).fetchdf()
            final = {"columns": list(res.columns), "rows": json.loads(res.to_json(orient="records")),
                     "error": None}
        except Exception as exc:  # noqa: BLE001
            final = {"columns": [], "rows": [], "error": f"{type(exc).__name__}: {str(exc).splitlines()[0][:400]}"}
    finally:
        con.close()
    return {"ok": first_error is None and final["error"] is None, "steps": steps, "first_error": first_error,
            "final": final, "tables_loaded": present + list(STUB_TABLES)}


def check_query_pack(omop_dir, phenotype_json) -> dict:
    """Run the stage-D All of Us BigQuery pack for `phenotype_json` (path or dict) against the export in `omop_dir`.
    Returns (and writes to <omop_dir>/query_pack_check.json) per-step row counts or the SQL error."""
    omop_dir = Path(omop_dir)
    result = {"checked_at": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
              "omop_dir": str(omop_dir), "source": "measure_it.similar.omop.build_sql(dialect='bigquery', "
                                                   "prefix='`{CDR}`.')",
              "dialect_in": "bigquery", "dialect_out": "duckdb"}
    try:
        import sqlglot
    except ImportError:
        result.update(status="sqlglot_missing", ok=False,
                      error="sqlglot is not installed (it is not a project dependency); run with "
                            "`uv run --with sqlglot ...` to execute the query pack")
        return _write(omop_dir, result)
    result["sqlglot_version"] = sqlglot.__version__
    try:
        ph = _load_phenotype(phenotype_json)
    except Exception as exc:  # noqa: BLE001
        result.update(status="phenotype_error", ok=False, error=f"{type(exc).__name__}: {exc}")
        return _write(omop_dir, result)
    result["phenotype_id"] = ph.get("phenotype_id")
    from ..similar.omop import build_sql
    try:
        sql_bq, info = build_sql(ph, dialect="bigquery", prefix=CDR_PLACEHOLDER)
    except Exception as exc:  # noqa: BLE001 - a stage-D generation failure is a finding, not a crash
        result.update(status="generation_error", ok=False, error=f"{type(exc).__name__}: {exc}")
        return _write(omop_dir, result)
    result["coverage"] = info["coverage"]
    result["evaluable"] = info["evaluable"]
    result["base_population"] = info["base_population"]
    result["threshold_lp"] = info["threshold_lp"]
    try:
        duck_sql, cdm_tables, ctes = transpile_pack(sql_bq)
    except Exception as exc:  # noqa: BLE001
        result.update(status="transpile_error", ok=False, error=f"{type(exc).__name__}: {str(exc)[:600]}")
        return _write(omop_dir, result)
    (omop_dir / "query_pack_check.bigquery.sql").write_text(sql_bq)
    (omop_dir / "query_pack_check.duckdb.sql").write_text(duck_sql + "\n")
    present = {p.stem for p in omop_dir.glob("*.parquet")} | set(STUB_TABLES)
    result["cdm_tables_referenced"] = cdm_tables
    result["cdm_tables_missing"] = sorted(set(cdm_tables) - present)
    result["stubs"] = {"concept_relationship": "'Maps to' self-map of every local concept",
                       "concept_ancestor": "self-ancestry of every local concept (levels 0)"}
    result["variants"] = {
        "as_exported": _run_variant(omop_dir, duck_sql, ctes, yob_from_band=False),
        "year_of_birth_from_age_band": _run_variant(omop_dir, duck_sql, ctes, yob_from_band=True),
    }
    ok = all(v["ok"] for v in result["variants"].values())
    result["ok"] = ok
    result["status"] = "executed" if ok else "sql_error"
    result["notes"] = [
        "First execution of a stage-D query pack: generated by measure_it.similar.omop, transpiled BigQuery -> "
        "DuckDB by sqlglot, `{CDR}`.<table> references mapped to the local export tables.",
        "as_exported: year_of_birth is null in the export, so the pack's adult filter keeps nobody (cohort = 0 is "
        "the expected result, not an error).",
        "year_of_birth_from_age_band: check-only year of birth from the age-band midpoint; ages 18-19 fall in the "
        "10-19 band (midpoint 14.5) and are excluded by the adult filter.",
        "State strata: the pack reads All of Us PPI state-of-residence answers (observation_source_concept_id "
        "1585249), which this export never holds (no geography), so every person is in stratum 'unknown'.",
        "Final rows apply All of Us suppression (counts 1-20 -> NULL), so a small synthetic run returns NULLs.",
    ]
    return _write(omop_dir, result)
