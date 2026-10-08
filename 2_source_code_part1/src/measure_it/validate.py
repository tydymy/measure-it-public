"""Build validation: `uv run measure-it validate` (exit code 0 = every check passed).

Checks
1. provenance    every data/processed/*.parquet passes provenance.check_provenance (the boundary geoparquets are
                 geometry support files, not data tables, and are listed as exempt);
2. required      every SPEC "REQUIRED FINAL DATASETS" file exists and is non-empty (the DuckDB must hold tables);
3. person_geo    no person-layer table (any row with data_layer = person) has a geographic identifier column
                 (county/state FIPS, ZIP, ZCTA, geo_id, latitude/longitude): person rows are never located;
4. object_ids    object_id is non-null/non-empty, namespaced ("<namespace>:<id>") and unique in the tables where it
                 identifies the row (OBJECT_ID_UNIQUE); link/long tables that repeat an object_id by design are listed
                 (an empty object_id there means "not addressable" and is counted, not failed);
5. language      the forbidden product terms do not appear in README.md, results/*.md, results/figures/*.md,
                 results/*.json (string values), docs/*.md, data/raw/*/DATA_AUDIT.md or dashboard/ text, nor in the
                 user-facing text of the MCP server and HTTP API (instructions, guardrails, tool descriptions, parameter
                 docs), except where quoted as a forbidden term (in quotes/backticks, or as a bare list item after a
                 rule word such as "Avoid:" on the same line);
6. audits        every SOURCE_REGISTRY.yaml source has data/raw/<source_id>/DATA_AUDIT.md, and every file its `audit`
                 field lists (';'-separated) exists.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pyarrow.parquet as pq

from .config import DOCS, DUCKDB_PATH, PROCESSED, PROJECT_ROOT, RAW, RESULTS, utc_now_iso
from .provenance import PROVENANCE_COLUMNS, check_provenance

REQUIRED_FINAL_DATASETS = [
    "participant_wearable_features", "participant_clinical_features", "condition_registry",
    "condition_molecular_evidence", "condition_measurement_evidence", "measurement_registry",
    "geo_condition_features", "clinic_registry", "research_site_registry", "deployment_opportunities",
]
REQUIRED_DB = DUCKDB_PATH

# Column names that locate a place. Person-layer tables must carry none of them.
GEO_COLUMN_RE = re.compile(
    r"^(?:geo_id|fips|county|county_fips|county_fips_\w+|state_fips|state_abbr|state_code|zip|zip5|zip_code|zipcode|"
    r"zcta|zcta5|postal_code|lat|lon|lng|latitude|longitude|geometry|census_tract|tract_fips|block_group)$"
    r"|(?:^|_)(?:fips|zcta|latitude|longitude)(?:_|$)|(?:^|_)(?:lat|lon)$",
    re.IGNORECASE)

# Tables whose object_id identifies the row (CONVENTIONS section 7). Other tables with an object_id column are link
# or long tables in which the id repeats by design (e.g. facility_trials: one row per facility x trial).
OBJECT_ID_UNIQUE = [
    "clinic_registry", "research_site_registry", "facilities", "condition_measurement_evidence", "measurement_registry",
    "deployment_opportunities", "deployment_candidates", "measurable_biology", "measurement_phenotype_signal",
    "phase3_measurement_evidence", "phenotype_signatures__nhanes", "phenotype_signatures__stanford",
    "phenotype_signatures", "geo_context", "measurement_performance",
]
OBJECT_ID_RE = re.compile(r"^[a-z_]+:\S")

FORBIDDEN_TERMS = ["proves", "diagnosed by AI", "patient has", "best clinic", "definitive biomarker",
                   "optimal treatment", "confirms mechanism", "cures"]
FORBIDDEN_RE = re.compile(r"\b(" + "|".join(re.escape(t) for t in FORBIDDEN_TERMS) + r")\b", re.IGNORECASE)
# a line that talks ABOUT the terms (a style rule), not one that uses them
MENTION_MARKER_RE = re.compile(r"\b(avoid|never|forbidden|do not (?:use|say|write)|don't (?:use|say)|banned|"
                               r"product language|not (?:say|write))\b", re.IGNORECASE)
QUOTES = "\"'`“”‘’"


def _check(name: str, ok: bool, details: list, summary: str) -> dict:
    return {"check": name, "ok": bool(ok), "summary": summary, "details": details}


def check_provenance_all() -> dict:
    bad, n = [], 0
    for p in sorted(PROCESSED.glob("*.parquet")):
        n += 1
        try:
            names = pq.read_schema(p).names
            cols = [c for c in PROVENANCE_COLUMNS if c in names]
            df = pq.read_table(p, columns=cols).to_pandas()
            check_provenance(df, p.stem)
        except Exception as exc:  # noqa: BLE001 - every failure is reported
            bad.append({"table": p.stem, "error": str(exc)[:300]})
    exempt = [p.name for p in sorted(PROCESSED.glob("*.geoparquet"))]
    if n == 0:   # nothing checked is not a pass
        bad.append({"problem": f"no processed tables under {PROCESSED.name}/ (run the pipeline first)"})
    return _check("provenance", not bad, bad,
                  f"{max(0, n - len(bad))}/{n} processed tables pass check_provenance; exempt geometry files: {exempt}")


def check_required() -> dict:
    missing = []
    for t in REQUIRED_FINAL_DATASETS:
        p = PROCESSED / f"{t}.parquet"
        if not p.exists():
            missing.append({"dataset": f"{t}.parquet", "problem": "missing"})
        elif pq.read_metadata(p).num_rows == 0:
            missing.append({"dataset": f"{t}.parquet", "problem": "empty"})
    if not REQUIRED_DB.exists():
        missing.append({"dataset": REQUIRED_DB.name, "problem": "missing"})
    else:
        import duckdb
        try:
            con = duckdb.connect(str(REQUIRED_DB), read_only=True)
            n = con.execute("SELECT count(*) FROM information_schema.tables WHERE table_type = 'BASE TABLE'").fetchone()[0]
            con.close()
            if n == 0:
                missing.append({"dataset": REQUIRED_DB.name, "problem": "no tables"})
        except Exception as exc:  # noqa: BLE001
            missing.append({"dataset": REQUIRED_DB.name, "problem": f"cannot open: {exc}"[:200]})
    return _check("required_datasets", not missing, missing,
                  f"{len(REQUIRED_FINAL_DATASETS) + 1 - len(missing)}/{len(REQUIRED_FINAL_DATASETS) + 1} SPEC required "
                  "final datasets present and non-empty")


def check_person_geography() -> dict:
    bad, person_tables = [], []
    for p in sorted(PROCESSED.glob("*.parquet")):
        names = pq.read_schema(p).names
        if "data_layer" not in names:
            continue
        layers = set(pq.read_table(p, columns=["data_layer"]).column("data_layer").unique().to_pylist())
        if "person" not in layers:
            continue
        person_tables.append(p.stem)
        geo_cols = [c for c in names if GEO_COLUMN_RE.search(c) and c not in PROVENANCE_COLUMNS]
        if geo_cols:
            bad.append({"table": p.stem, "geographic_columns": geo_cols})
    if not person_tables:   # nothing checked is not a pass
        bad.append({"problem": "no person-layer table found (run the pipeline first)"})
    return _check("person_layer_no_geography", not bad, bad,
                  f"{max(0, len(person_tables) - len(bad))}/{len(person_tables)} person-layer tables carry no geographic "
                  f"identifier column")


def check_object_ids() -> dict:
    import duckdb
    con = duckdb.connect()
    bad, checked, repeated = [], [], []
    for p in sorted(PROCESSED.glob("*.parquet")):
        if "object_id" not in pq.read_schema(p).names:
            continue
        # an empty object_id means "not addressable" (e.g. ClinicalTrials.gov / NIH rows of facility_source_links,
        # which have no CONVENTIONS namespace); it counts as missing, never as a malformed id
        n, nn, u, badfmt = con.execute(
            "SELECT count(*), count(NULLIF(object_id, '')), count(DISTINCT NULLIF(object_id, '')), "
            "count(*) FILTER (WHERE NULLIF(object_id, '') IS NOT NULL AND NOT regexp_matches(object_id, '^[a-z_]+:\\S')) "
            "FROM read_parquet(?)", [str(p)]).fetchone()
        rec = {"table": p.stem, "rows": n, "non_null": nn, "distinct": u, "bad_format": badfmt}
        if p.stem in OBJECT_ID_UNIQUE:
            checked.append(p.stem)
            if not (n == nn == u) or badfmt:
                bad.append(rec)
        else:
            repeated.append(f"{p.stem} ({u:,} distinct / {n:,} rows" + (f", {n - nn:,} not addressable)" if n > nn else ")"))
            if badfmt:
                bad.append(rec)
    missing = [t for t in OBJECT_ID_UNIQUE if not (PROCESSED / f"{t}.parquet").exists()]
    if not checked:   # nothing checked is not a pass
        bad.append({"problem": "no row-identity table with object_id found (run the pipeline first)"})
    return _check("object_ids", not bad, bad,
                  f"object_id unique, non-null and namespaced in {max(0, len(checked) - len(bad))}/{len(checked)} row-identity "
                  f"tables; repeated by design (link/long tables, format checked): {repeated}; "
                  f"declared but not built: {missing}")


def _text_files() -> list[Path]:
    files = sorted(RESULTS.glob("*.md")) + sorted(DOCS.glob("*.md"))
    files += [p for p in [PROJECT_ROOT / "README.md"] if p.exists()]
    files += sorted((RESULTS / "figures").glob("*.md")) + sorted(RESULTS.glob("*.json"))
    files += sorted(RAW.glob("*/DATA_AUDIT.md"))
    dash = PROJECT_ROOT / "dashboard"
    if dash.exists():
        files += sorted(p for p in dash.rglob("*") if p.is_file() and p.suffix.lower() in
                        {".md", ".py", ".html", ".txt", ".json", ".yaml", ".yml", ".js", ".ts", ".tsx", ".jsx"})
    return files


def _json_strings(obj) -> list[str]:
    """Every string value (and key) of a JSON document, so escaped quotes (\\") do not hide or fake a quoted term."""
    out, stack = [], [obj]
    while stack:
        v = stack.pop()
        if isinstance(v, dict):
            out += [k for k in v if isinstance(k, str)]
            stack += list(v.values())
        elif isinstance(v, list):
            stack += v
        elif isinstance(v, str):
            out.append(v)
    return out


def service_texts() -> dict[str, str]:
    """User-facing text of the MCP server and the HTTP API: instructions, guardrail prompt/resources, every tool's
    title, description and parameter docs (what a calling model or API client reads)."""
    from . import tools as T
    from .mcp import server as S
    out = {"mcp.INSTRUCTIONS": S.INSTRUCTIONS, "mcp.GUARDRAILS_MD": S.GUARDRAILS_MD,
           "mcp.DATA_LAYERS_MD": S.DATA_LAYERS_MD, "mcp.TOOL_PREAMBLE": S.TOOL_PREAMBLE}
    for name, fn in T.TOOLS.items():
        out[f"tool:{name}"] = "\n".join([str(S.TITLES.get(name, "")), S._description(name, fn)])
    for k, v in {**S.PARAM_DOCS, **{f"{a}.{b}": d for (a, b), d in S.PER_TOOL_DOCS.items()}}.items():
        out[f"param:{k}"] = str(v)
    try:
        from .api.app import DESCRIPTION
        out["api.DESCRIPTION"] = DESCRIPTION
    except Exception:  # noqa: BLE001 - the API extra may be absent; the MCP text is still checked
        pass
    out["mcp.object_ids"] = S._object_ids_md()
    out["mcp.tools"] = S._tools_md()
    return {k: v for k, v in out.items() if v}


_ITEM_OPEN = ":;,("
_ITEM_CLOSE = ":;,)."


def _listed_after_marker(line: str, start: int, end: int) -> bool:
    """True when the term at line[start:end] is a bare item of a list that a style rule introduces earlier on the
    line ("Avoid (unless directly supported): diagnosed by AI; proves; cures; ..."): a marker word precedes it and
    the text between the nearest list delimiters around it is the term alone ("a, b and c" lists included).
    A marker elsewhere on the line does not excuse a term used in a sentence ("this proves X, never Y")."""
    if not MENTION_MARKER_RE.search(line[:start]):
        return False
    lo = max(line.rfind(ch, 0, start) for ch in _ITEM_OPEN)
    his = [line.find(ch, end) for ch in _ITEM_CLOSE]
    hi = min([h for h in his if h >= 0], default=len(line))
    pieces = re.split(r"(?:^|\s+)(?:and|or)\s+", line[lo + 1:hi].strip(), flags=re.IGNORECASE)
    return line[start:end].lower() in {p.strip().lower() for p in pieces}


def forbidden_hits(text: str) -> list[tuple[int, str, str]]:
    """(line number, term, line) for every forbidden term not quoted as a forbidden term.

    A term is excused only when it is quoted ("proves", `cures`) or is a bare item of a list introduced by a rule
    word earlier on the same line (avoid / never / forbidden / banned / do not use ...)."""
    hits = []
    for i, line in enumerate(text.splitlines(), 1):
        for m in FORBIDDEN_RE.finditer(line):
            before = line[:m.start()].rstrip()[-1:] if line[:m.start()].rstrip() else ""
            after = line[m.end():].lstrip()[:1]
            quoted = before in QUOTES and after in QUOTES
            if quoted or _listed_after_marker(line, m.start(), m.end()):
                continue
            hits.append((i, m.group(0), line.strip()[:200]))
    return hits


def check_language() -> dict:
    bad, n = [], 0
    for f in _text_files():
        n += 1
        try:
            text = f.read_text(errors="replace")
            if f.suffix.lower() == ".json":
                text = "\n".join(_json_strings(json.loads(text)))
        except (OSError, ValueError):
            continue
        for ln, term, line in forbidden_hits(text):
            bad.append({"file": str(f.relative_to(PROJECT_ROOT)), "line": ln, "term": term, "text": line})
    try:
        svc = service_texts()
    except Exception as exc:  # noqa: BLE001 - reported, not raised
        svc = {}
        bad.append({"file": "MCP/API service text", "line": 0, "term": "", "text": f"could not load: {exc}"[:200]})
    for name, text in svc.items():
        for ln, term, line in forbidden_hits(text):
            bad.append({"file": name, "line": ln, "term": term, "text": line})
    return _check("forbidden_language", not bad, bad,
                  f"{n} text files and {len(svc)} MCP/API service texts scanned for {len(FORBIDDEN_TERMS)} forbidden "
                  f"terms; {len(bad)} unquoted use(s)")


def check_audits() -> dict:
    import yaml
    from .config import SOURCE_REGISTRY_PATH
    reg = yaml.safe_load(SOURCE_REGISTRY_PATH.read_text())
    bad, n = [], 0
    for s in reg.get("sources", []):
        n += 1
        sid = s["source_id"]
        # `audit` may list several files separated by ';' (e.g. MapMECFS: DATA_AUDIT.md; the linkage audit)
        listed = [PROJECT_ROOT / a.strip() for a in str(s.get("audit") or "").split(";") if a.strip()]
        own = RAW / sid / "DATA_AUDIT.md"
        if not (own.exists() and own.stat().st_size > 0):
            bad.append({"source_id": sid, "audit": s.get("audit"), "problem": f"{own.relative_to(PROJECT_ROOT)} missing"})
        missing = [str(a.relative_to(PROJECT_ROOT)) for a in listed if not a.exists()]
        if missing:
            bad.append({"source_id": sid, "audit": s.get("audit"), "problem": f"audit path(s) in registry do not exist: {missing}"})
    frag = sorted(p.parent.name for p in RAW.glob("*/registry_entry.yaml"))
    unmerged = sorted(set(frag) - {s["source_id"] for s in reg.get("sources", [])})
    for sid in unmerged:
        bad.append({"source_id": sid, "problem": "registry_entry.yaml not merged into SOURCE_REGISTRY.yaml"})
    return _check("data_audits", not bad, bad, f"{n - len(bad)}/{n} SOURCE_REGISTRY sources have a DATA_AUDIT.md"
                  + (f"; unmerged fragments: {unmerged}" if unmerged else ""))


CHECKS = [check_provenance_all, check_required, check_person_geography, check_object_ids, check_language, check_audits]


def run_checks(write: bool = True) -> dict:
    results = [c() for c in CHECKS]
    ok = all(r["ok"] for r in results)
    report = {"validated_at": utc_now_iso(), "ok": ok, "checks": results}
    lines = [f"measure-it validate ({report['validated_at']}): {'PASS' if ok else 'FAIL'}"]
    for r in results:
        lines.append(f"  [{'PASS' if r['ok'] else 'FAIL'}] {r['check']}: {r['summary']}")
        for d in r["details"][:25] if not r["ok"] else []:
            lines.append(f"         - {json.dumps(d, default=str)[:300]}")
        if not r["ok"] and len(r["details"]) > 25:
            lines.append(f"         ... {len(r['details']) - 25} more")
    text = "\n".join(lines)
    print(text, flush=True)
    if write:
        out = RESULTS / "pipeline_runs"
        out.mkdir(parents=True, exist_ok=True)
        (out / "validate_latest.json").write_text(json.dumps(report, indent=1, default=str))
    return report
