"""trace_evidence(object_id): resolve any cited object id back to its rows, provenance and raw source files.

Every row a recommendation cites is addressable as ``<namespace>:<native id>`` (docs/CONVENTIONS.md section 7).
``trace_evidence`` resolves every namespace of that table to

* ``rows``        the row(s) of the primary table (long text cut to ``max_text`` characters, rows capped),
* ``provenance``  the distinct provenance-column values of those rows (source_name, source_version, retrieved_at,
                  data_layer, evidence_type, evidence_level, source_geographic_resolution, provenance_notes),
* ``sources``     the SOURCE_REGISTRY.yaml entry of every source behind the rows (name, publisher, version,
                  retrieval time, licence, status), the path of its data/raw/<source>/DATA_AUDIT.md, its MANIFEST.json
                  and the raw files that hold (or were named for) this record when they can be identified,
* ``lineage``     one level of lineage: the object ids and rows this record was built from or points to
                  (opportunity -> component values -> geo_condition_features and burden rows and their sources;
                  signature -> dataset and label basis; facility -> the source records that were resolved into it; ...).

Unknown namespaces or ids return ``UNKNOWN / NOT AVAILABLE`` with the reason. Person-level ids (``participant:``) are
refused: person-level rows are summarised only (group statistics such as ``signature:`` rows) and are never exposed
through this public tool.

Lookups read the processed Parquet tables (the canonical store). When data/processed/measure_it_public.duckdb was
built after a table's Parquet file, the same table is read from the DuckDB (point lookups ~10 ms instead of a scan);
otherwise the Parquet file is scanned. Both give the same rows.
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import os
import re
import threading
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DOCS, DUCKDB_PATH, PROCESSED, PROJECT_ROOT, RAW, SOURCE_REGISTRY_PATH, TABLES, UNKNOWN
from ..provenance import PROVENANCE_COLUMNS

MAX_ROWS = 20
MAX_TEXT = 400
MAX_LINEAGE = 25

PERSON_LEVEL_NAMESPACES = {"participant", "person", "patient", "seqn", "subject"}
# participant_id is "<dataset_id>:<native id>" (CONVENTIONS section 2): a bare participant id (nhanes:62161,
# stanford_covid_mishra2020:A06L7KF, mapmecfs:MECFS_101) is refused too, whatever its case.
# The lab / device datasets added 2026-09-24 use their source_id as dataset_id (e.g. klein2023_mylc_ml_table:LC.MS.0001).
PERSON_DATASET_PREFIXES = ("nhanes", "stanford", "mapmecfs", "charlton_lc_mecfs_cpet_source", "fm_thermography",
                           "endo_arg1_repod", "heds_hsd_olink_serum_cinquina2026", "klein2023_mylc_ml_table",
                           "appelman_lc_pem_source")
_PREFIX_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_\-]*)\s*:")


def is_person_level_id(object_id: str) -> bool:
    m = _PREFIX_RE.match(str(object_id or ""))
    if not m:
        return False
    ns = m.group(1).lower()
    return ns in PERSON_LEVEL_NAMESPACES or ns.startswith(PERSON_DATASET_PREFIXES)


PERSON_REFUSAL =("refused: person-level rows (participants of NHANES, the Stanford wearable studies, MapMECFS or the "
                  "public lab / device datasets such as MUSCLE-ME, MY-LC, the fibromyalgia thermography, endometriosis "
                  "arginase and hEDS/HSD Olink releases) are "
                  "summarised only. This public tool never returns an individual's record; trace the group-level "
                  "statistics instead (signature:, measurement_signal:, phenotype_evidence: object ids).")

# CONVENTIONS section 7: namespace -> (primary table, native id description)
NAMESPACES: dict[str, tuple[str, str]] = {
    "condition": ("condition_registry", "canonical_condition_id (or a condition-set id)"),
    "measurement": ("measurement_registry", "measurement class id or bundle id"),
    "trial": ("clinical_trials", "NCT id"),
    "nih": ("nih_projects", "RePORTER appl_id"),
    "npi": ("providers", "NPI"),
    "hrsa_site": ("facilities__hrsa", "BPHC site number"),
    "facility": ("facilities", "facility_id"),
    "fda": ("fda_510k / fda_pma / fda_device_classification", "K/DEN/P number or code:<product code>"),
    "gene": ("condition_molecular_evidence", "Ensembl gene id"),
    "variant": ("condition_molecular_evidence", "rsID"),
    "gwas_study": ("condition_molecular_evidence", "GCST accession"),
    "geo_series": ("geo_study_catalog", "GSE accession"),
    "geo": ("geographies / geo_context", "FIPS (state 2, county 5), US, or zcta:<ZCTA>"),
    "opportunity": ("deployment_opportunities", "<condition>|<measurement>|<geo_id>"),
    "signature": ("phenotype_signatures", "<dataset>|<phenotype>|<feature>"),
    "measurement_evidence": ("condition_measurement_evidence", "<condition_id>|<measurement_id>"),
    "measurement_signal": ("measurement_phenotype_signal", "<dataset>|<label>|<class>|<scope>"),
    "phenotype_evidence": ("phase3_measurement_evidence", "<phenotype>|<measurement>"),
    "molbio": ("measurable_biology", "<condition>|<system>|<measurement_class>[|<catalogue>]"),
    "reactome": ("results/tables/test3_reactome_ora.csv", "Reactome stable id (R-HSA-...)"),
    "measurement_performance": ("measurement_performance", "<condition_id>|<measurement_id>"),
}

# Derived tables that no SOURCE_REGISTRY entry lists as a processed output: producer and the sources they draw on.
DERIVED_TABLES: dict[str, dict] = {
    "measurement_registry": {"producer": "measure_it.measurements.evidence",
                             "sources": ["clinicaltrials_gov", "nih_reporter", "openfda_device"],
                             "configs": ["configs/measurements.yaml", "configs/relevance.yaml",
                                         "configs/fda_product_code_map.yaml"]},
    "condition_measurement_evidence": {"producer": "measure_it.measurements.evidence",
                                       "sources": ["clinicaltrials_gov", "nih_reporter", "openfda_device"],
                                       "configs": ["configs/measurements.yaml", "configs/relevance.yaml"]},
    "measurement_phenotype_signal": {"producer": "measure_it.measurements.phenotype_evidence",
                                     "sources": ["nhanes_2011_2014", "stanford_longcovid_wearables"], "configs": []},
    "phase3_measurement_evidence": {"producer": "measure_it.measurements.phenotype_evidence",
                                    "sources": ["nhanes_2011_2014", "stanford_longcovid_wearables", "clinicaltrials_gov",
                                                "nih_reporter", "openfda_device"],
                                    "configs": ["configs/measurements.yaml"]},
    "measurable_biology": {"producer": "measure_it.omics.graph",
                           "sources": ["open_targets", "gwas_catalog", "mapmecfs_nih_pi_mecfs", "reactome", "hgnc"],
                           "configs": []},
    "facilities": {"producer": "measure_it.facilities.registry",
                   "sources": ["nppes", "nucc_taxonomy", "hrsa_health_centers", "clinicaltrials_gov", "nih_reporter",
                               "hrsa_hpsa", "usda_rucc"],
                   "configs": ["configs/specialty_groups.yaml"]},
    "facility_source_links": {"producer": "measure_it.facilities.registry",
                              "sources": ["nppes", "hrsa_health_centers", "clinicaltrials_gov", "nih_reporter"],
                              "configs": []},
    "facility_trials": {"producer": "measure_it.facilities.registry", "sources": ["clinicaltrials_gov"], "configs": []},
    "facility_nih_projects": {"producer": "measure_it.facilities.registry", "sources": ["nih_reporter"], "configs": []},
    "geo_context": {"producer": "measure_it.geography.features",
                    "sources": ["census_geography", "census_acs", "cdc_svi", "cdc_places", "usda_rucc", "hrsa_hpsa"],
                    "configs": []},
    "geo_context_zcta": {"producer": "measure_it.pipeline canonical_unions",
                         "sources": ["census_acs", "cdc_places", "census_geography"], "configs": []},
    "geo_condition_features": {"producer": "measure_it.geography.features",
                               "sources": ["cdc_long_covid", "cdc_brfss", "cms_mmd", "cdc_lyme", "cdc_places", "cdc_svi", "census_acs",
                                           "nppes", "clinicaltrials_gov", "nih_reporter", "hrsa_health_centers",
                                           "hrsa_hpsa", "usda_rucc"],
                               "configs": ["docs/BURDEN_DEFINITIONS.md", "configs/relevance.yaml"]},
    "geo_condition_burden": {"producer": "measure_it.pipeline canonical_unions (per-source partitions)",
                             "sources": ["cdc_long_covid", "cdc_brfss", "cdc_lyme", "cms_mmd", "cdc_places"], "configs": []},
    "deployment_opportunities": {"producer": "measure_it.scoring.opportunity",
                                 "sources": [], "configs": ["configs/scoring.yaml", "configs/relevance.yaml"]},
    "phenotype_signatures": {"producer": "measure_it.pipeline canonical_unions of phenotype_signatures__* "
                                         "(measure_it.wearables.signatures / stanford_models / muscle_me_steps_analysis "
                                         "/ fm_thermography_analysis, measure_it.omics.person_linked / geo_cohorts, "
                                         "measure_it.labs.*)",
                             "sources": ["nhanes_2011_2014", "stanford_longcovid_wearables", "mapmecfs_nih_pi_mecfs",
                                         "geo_cohorts", "charlton_lc_mecfs_cpet_source", "fm_thermography",
                                         "appelman_lc_pem_source", "endo_arg1_repod",
                                         "heds_hsd_olink_serum_cinquina2026", "klein2023_mylc_ml_table"],
                             "configs": []},
    "measurement_performance": {"producer": "measure_it.scoring.metric_link",
                                "sources": ["charlton_lc_mecfs_cpet_source", "stanford_longcovid_wearables",
                                            "nhanes_2011_2014", "appelman_lc_pem_source", "klein2023_mylc_ml_table",
                                            "published_device_evidence"],
                                "configs": ["docs/ANALYSIS_PLAN_METRIC_LINK.md", "configs/relevance.yaml"]},
    "measurement_performance_records": {"producer": "measure_it.scoring.metric_link",
                                        "sources": ["charlton_lc_mecfs_cpet_source", "stanford_longcovid_wearables",
                                                    "nhanes_2011_2014", "appelman_lc_pem_source",
                                                    "klein2023_mylc_ml_table", "published_device_evidence"],
                                        "configs": ["docs/ANALYSIS_PLAN_METRIC_LINK.md"]},
    "test3_reactome_ora": {"producer": "measure_it.omics.coherence",
                           "sources": ["reactome", "hgnc", "open_targets", "gwas_catalog", "mapmecfs_nih_pi_mecfs"],
                           "configs": []},
}

FACILITY_SOURCE_CODES = {"nppes": "nppes", "hrsa": "hrsa_health_centers", "ctgov": "clinicaltrials_gov",
                         "nih": "nih_reporter"}
MOLECULAR_SOURCE_DATABASES = {"GWAS Catalog": "gwas_catalog", "NCBI GEO": "ncbi_geo_sra", "NCBI SRA": "ncbi_geo_sra",
                              "Open Targets Platform": "open_targets",
                              "Open Targets Platform (Reactome annotation)": "open_targets"}


# ============================================================================================ JSON helpers

def to_jsonable(v, *, max_text: int | None = None):
    """Recursively convert to JSON-safe Python: numpy/pandas scalars -> python, NaN/inf/NA/NaT -> None,
    timestamps -> ISO strings, sets -> sorted lists, arrays/tuples -> lists. Strings longer than max_text are cut
    (with a marker) when max_text is given."""
    if v is None:
        return None
    if isinstance(v, dict):
        return {str(k): to_jsonable(x, max_text=max_text) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [to_jsonable(x, max_text=max_text) for x in v]
    if isinstance(v, (set, frozenset)):
        return [to_jsonable(x, max_text=max_text) for x in sorted(v, key=str)]
    if isinstance(v, np.ndarray):
        return [to_jsonable(x, max_text=max_text) for x in v.tolist()]
    if isinstance(v, pd.DataFrame):
        return [to_jsonable(r, max_text=max_text) for r in v.to_dict("records")]
    if isinstance(v, pd.Series):
        return to_jsonable(v.to_dict(), max_text=max_text)
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        f = float(v)
        return f if math.isfinite(f) else None
    if isinstance(v, str):
        if max_text and len(v) > max_text:
            return v[:max_text] + f" ...[cut at {max_text} of {len(v)} characters]"
        return v
    if v is pd.NA or v is pd.NaT:
        return None
    if isinstance(v, (pd.Timestamp, _dt.datetime, _dt.date)):
        return v.isoformat()
    if isinstance(v, Path):
        return str(v)
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(v, "item"):
        try:
            return to_jsonable(v.item(), max_text=max_text)
        except (TypeError, ValueError):
            pass
    return str(v)


def _maybe_json(v):
    """Decode a JSON list/object stored as text (e.g. uncertainties, example_*_object_ids)."""
    if isinstance(v, str) and v[:1] in "[{":
        try:
            return json.loads(v)
        except ValueError:
            return v
    return v


def _rel(p: Path) -> str:
    try:
        return str(Path(p).relative_to(PROJECT_ROOT))
    except ValueError:
        return str(p)


# ============================================================================================ lookups

def _q(name: str) -> str:
    return '"' + str(name).replace('"', '""') + '"'


class _Lookup:
    """Point lookups on processed tables: the DuckDB when it is newer than the table's Parquet, else the Parquet."""

    def __init__(self) -> None:
        self._lock = threading.RLock()
        self._file_con = None
        self._file_mtime: float | None = None
        self._file_tables: set[str] = set()
        self._mem = None

    def _file(self):
        p = DUCKDB_PATH
        if not p.exists():
            return None
        m = p.stat().st_mtime
        if self._file_con is None or m != self._file_mtime:
            import duckdb
            if self._file_con is not None:
                try:
                    self._file_con.close()
                except Exception:  # noqa: BLE001 - closing a stale handle
                    pass
                self._file_con = None
            try:
                con = duckdb.connect(str(p), read_only=True)
                self._file_tables = {r[0] for r in con.execute(
                    "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main' "
                    "AND table_type = 'BASE TABLE'").fetchall()}
            except Exception:  # noqa: BLE001 - locked or unreadable: fall back to Parquet
                self._file_con, self._file_mtime, self._file_tables = None, None, set()
                return None
            self._file_con, self._file_mtime = con, m
        return self._file_con

    def _memcon(self):
        if self._mem is None:
            import duckdb
            self._mem = duckdb.connect()
        return self._mem

    def backend(self, table: str) -> str:
        pq_path = PROCESSED / f"{table}.parquet"
        with self._lock:
            con = self._file()
            if con is not None and table in self._file_tables and self._file_mtime >= pq_path.stat().st_mtime:
                return "duckdb"
        return "parquet"

    def query(self, table: str, where: str, params: list | tuple = (), columns: list[str] | None = None,
              limit: int | None = None, order: str | None = None) -> pd.DataFrame:
        pq_path = PROCESSED / f"{table}.parquet"
        if not pq_path.exists():
            raise FileNotFoundError(f"table {table} is not built ({_rel(pq_path)} missing)")
        schema = _schema(table, pq_path.stat().st_mtime)
        cols = [c for c in (columns or schema) if c in schema]
        sel = ", ".join(_q(c) for c in cols) if columns else "*"
        tail = (f" ORDER BY {order}" if order else "") + (f" LIMIT {int(limit)}" if limit else "")
        with self._lock:
            con = self._file()
            if con is not None and table in self._file_tables and self._file_mtime >= pq_path.stat().st_mtime:
                cur = con.cursor()
                src = f"main.{_q(table)}"
            else:
                cur = self._memcon().cursor()
                src = "read_parquet(" + "'" + str(pq_path).replace("'", "''") + "')"
            try:
                return cur.execute(f"SELECT {sel} FROM {src} WHERE {where}{tail}", list(params)).fetchdf()
            finally:
                cur.close()

    def count(self, table: str, where: str, params: list | tuple = ()) -> int:
        pq_path = PROCESSED / f"{table}.parquet"
        if not pq_path.exists():
            return 0
        with self._lock:
            con = self._file()
            if con is not None and table in self._file_tables and self._file_mtime >= pq_path.stat().st_mtime:
                cur, src = con.cursor(), f"main.{_q(table)}"
            else:
                cur = self._memcon().cursor()
                src = "read_parquet(" + "'" + str(pq_path).replace("'", "''") + "')"
            try:
                return int(cur.execute(f"SELECT count(*) FROM {src} WHERE {where}", list(params)).fetchone()[0])
            finally:
                cur.close()


@lru_cache(maxsize=256)
def _schema(table: str, mtime: float) -> list[str]:
    import pyarrow.parquet as pq
    return list(pq.read_schema(PROCESSED / f"{table}.parquet").names)


LOOKUP = _Lookup()


def has_table(table: str) -> bool:
    return (PROCESSED / f"{table}.parquet").exists()


# ============================================================================================ source registry

def _mtime(p: Path) -> float:
    try:
        return p.stat().st_mtime
    except OSError:
        return 0.0


@lru_cache(maxsize=4)
def _registry(mtime: float) -> dict:
    import yaml
    reg = yaml.safe_load(SOURCE_REGISTRY_PATH.read_text()) if SOURCE_REGISTRY_PATH.exists() else {"sources": []}
    by_id = {s["source_id"]: s for s in reg.get("sources", []) or []}
    by_table: dict[str, list[str]] = {}
    for s in by_id.values():
        for t in s.get("processed_outputs") or []:
            by_table.setdefault(str(t).replace(".parquet", ""), []).append(s["source_id"])
    return {"by_id": by_id, "by_table": by_table, "generated_at": reg.get("generated_at")}


def registry() -> dict:
    return _registry(_mtime(SOURCE_REGISTRY_PATH))


def sources_for_table(table: str) -> list[str]:
    reg = registry()
    base = table.split("__")[0] if table not in reg["by_table"] else table
    if table in reg["by_table"]:
        return list(reg["by_table"][table])
    if table in DERIVED_TABLES:
        return list(DERIVED_TABLES[table]["sources"])
    if base in reg["by_table"]:
        return list(reg["by_table"][base])
    return list(DERIVED_TABLES.get(base, {}).get("sources", []))


@lru_cache(maxsize=64)
def _manifest(source_id: str, mtime: float) -> dict:
    p = RAW / source_id / "MANIFEST.json"
    try:
        return json.loads(p.read_text()).get("files", {})
    except (OSError, ValueError):
        return {}


def manifest(source_id: str) -> dict:
    return _manifest(source_id, _mtime(RAW / source_id / "MANIFEST.json"))


@lru_cache(maxsize=64)
def _raw_listing(source_id: str, mtime: float) -> tuple[str, ...]:
    base = RAW / source_id
    out = []
    for root, _dirs, files in os.walk(base):
        for f in files:
            if f in ("MANIFEST.json", "MANIFEST.lock", "registry_entry.yaml"):
                continue
            out.append(str(Path(root, f).relative_to(base)))
    return tuple(sorted(out))


def raw_listing(source_id: str) -> tuple[str, ...]:
    return _raw_listing(source_id, _mtime(RAW / source_id))


def raw_files(source_id: str, hints: list[str] | None = None, *, max_files: int = 5) -> dict:
    """Raw files of a source that hold (or were named for) a record: files whose path contains a hint (an accession,
    a product code, a file name quoted in source_version), else every file when the source has at most 3."""
    man = manifest(source_id)
    listing = raw_listing(source_id)
    hints = [h for h in (hints or []) if isinstance(h, str) and len(h) >= 3]
    hit: list[str] = []
    for h in hints:
        hl = h.lower()
        for f in listing:
            if hl in f.lower() and f not in hit:
                hit.append(f)
        for f in man:
            if Path(f).name.lower() in hl and f not in hit:   # file name quoted in a version string
                hit.append(f)
    basis = "file path contains the record's accession / code, or the file is named in source_version" if hit else ""
    if not hit and 0 < len(man) <= 3:
        hit, basis = list(man), "all files of this source (it has at most 3 raw files)"
    hit.sort(key=_raw_priority)
    files = []
    for f in hit[:max_files]:
        m = man.get(f, {})
        url = m.get("url")
        if isinstance(url, str) and len(url) > 200:
            url = url[:200] + f" ...[{len(url)} characters]"
        files.append({"path": _rel(RAW / source_id / f), "url": url, "bytes": m.get("bytes"),
                      "sha256": m.get("sha256"), "retrieved_at": m.get("retrieved_at"),
                      "in_manifest": f in man})
    return {"manifest": _rel(RAW / source_id / "MANIFEST.json") if (RAW / source_id / "MANIFEST.json").exists()
            else UNKNOWN, "n_manifest_files": len(man), "n_raw_files_on_disk": len(listing),
            "record_files": files, "n_record_files_matched": len(hit),
            "record_file_basis": basis or "the record's raw file cannot be singled out; see the MANIFEST (url, sha256, "
                                          "retrieved_at of every raw file)"}


DATA_EXT = (".csv", ".json", ".jsonl", ".gz", ".zip", ".xpt", ".parquet", ".tsv", ".txt", ".xlsx", ".yaml")


def _raw_priority(path: str) -> tuple:
    """Data files before documentation (pdf/html), then by name."""
    p = path.lower()
    ext = next((i for i, e in enumerate(DATA_EXT) if p.endswith(e)), len(DATA_EXT))
    return (ext >= len(DATA_EXT), "metadata" in p or "documentation" in p, p)


def source_summary(source_id: str, hints: list[str] | None = None) -> dict:
    s = registry()["by_id"].get(source_id)
    if s is None:
        return {"source_id": source_id, "status": UNKNOWN, "reason": "not in SOURCE_REGISTRY.yaml"}
    audits = [a.strip() for a in str(s.get("audit") or "").split(";") if a.strip()]
    own = RAW / source_id / "DATA_AUDIT.md"
    if _rel(own) not in audits:
        audits.insert(0, _rel(own))
    return {
        "source_id": source_id, "name": s.get("name"), "publisher": s.get("publisher"),
        "data_layer": s.get("data_layer"), "status": s.get("status"), "source_version": s.get("source_version"),
        "retrieved_at": s.get("retrieved_at"), "update_date": s.get("update_date"), "license": s.get("license"),
        "access_conditions": s.get("access_conditions"), "landing_url": s.get("landing_url"),
        "person_level": s.get("person_level"), "geographic": s.get("geographic"),
        "geographic_resolution": s.get("geographic_resolution"),
        "data_audit": [{"path": a, "exists": (PROJECT_ROOT / a).exists()} for a in audits],
        "ingestion_module": s.get("ingestion_module"),
        "raw": raw_files(source_id, hints),
    }


# ============================================================================================ result helpers

PROV_KEYS = [c for c in PROVENANCE_COLUMNS]


def _rows(df: pd.DataFrame, max_rows: int, max_text: int, drop: tuple = ()) -> tuple[list[dict], dict | None]:
    keep = [c for c in df.columns if c not in drop]
    recs = []
    for r in df[keep].head(max_rows).to_dict("records"):
        recs.append({k: to_jsonable(_maybe_json(v), max_text=max_text) for k, v in r.items()})
    trunc = ({"returned": min(len(df), max_rows), "total": int(len(df))} if len(df) > max_rows else None)
    return recs, trunc


def _prov(df: pd.DataFrame, table: str) -> list[dict]:
    cols = [c for c in PROV_KEYS if c in df.columns]
    if not cols or df.empty:
        return []
    d = df[cols].astype(str).drop_duplicates().head(10)
    return [{"table": table, **{k: (None if v in ("nan", "None", "<NA>") else v) for k, v in r.items()}}
            for r in d.to_dict("records")]


def _link(object_id: str, relation: str, **extra) -> dict:
    return {"object_id": object_id, "relation": relation, **{k: to_jsonable(v, max_text=MAX_TEXT) for k, v in extra.items()}}


class Trace:
    """Accumulates one trace result."""

    def __init__(self, object_id: str, namespace: str, native: str, max_rows: int, max_text: int):
        self.object_id, self.namespace, self.native = object_id, namespace, native
        self.max_rows, self.max_text = max_rows, max_text
        self.tables: list[str] = []
        self.rows: list[dict] = []
        self.n_rows = 0
        self.provenance: list[dict] = []
        self.source_ids: list[str] = []
        self.hints: dict[str, list[str]] = {}
        self.lineage: list[dict] = []
        self.lineage_truncated: dict | None = None
        self.notes: list[str] = []
        self.caveats: list[str] = []
        self.truncation: list[dict] = []
        self.extra: dict = {}

    def add_rows(self, table: str, df: pd.DataFrame, drop: tuple = ()) -> None:
        if df is None or df.empty:
            return
        recs, tr = _rows(df, self.max_rows, self.max_text, drop)
        self.tables.append(table)
        for r in recs:
            r["_table"] = table
        self.rows.extend(recs)
        self.n_rows += int(len(df))
        if tr:
            self.truncation.append({"path": f"data.rows[{table}]", **tr})
        self.provenance.extend(_prov(df, table))

    def add_sources(self, ids, hints: list[str] | None = None) -> None:
        for s in ids or []:
            if s and s not in self.source_ids:
                self.source_ids.append(s)
            if hints:
                self.hints.setdefault(s, [])
                self.hints[s] += [h for h in hints if h not in self.hints[s]]

    def add_link(self, item: dict) -> None:
        self.lineage.append(item)

    def result(self) -> dict:
        lin = self.lineage
        if len(lin) > MAX_LINEAGE:
            self.truncation.append({"path": "data.lineage", "returned": MAX_LINEAGE, "total": len(lin)})
            lin = lin[:MAX_LINEAGE]
        return {
            "status": "ok", "object_id": self.object_id, "namespace": self.namespace, "native_id": self.native,
            "tables": list(dict.fromkeys(self.tables)), "n_rows": self.n_rows, "rows": self.rows,
            "provenance": self.provenance,
            "sources": [source_summary(s, self.hints.get(s)) for s in self.source_ids],
            "lineage": lin, "notes": self.notes, "caveats": self.caveats, "truncation": self.truncation,
            "lookup_backend": {t: LOOKUP.backend(t) for t in dict.fromkeys(self.tables) if has_table(t)},
            **self.extra,
        }


def _unknown(object_id: str, reason: str, **extra) -> dict:
    return {"status": UNKNOWN, "object_id": object_id, "reason": reason, **extra}


def _eq(table: str, col: str, value, **kw) -> pd.DataFrame:
    return LOOKUP.query(table, f"{_q(col)} = ?", [value], **kw)


def _in(table: str, col: str, values: list, **kw) -> pd.DataFrame:
    values = [v for v in values if v is not None]
    if not values:
        return pd.DataFrame()
    ph = ", ".join("?" for _ in values)
    return LOOKUP.query(table, f"{_q(col)} IN ({ph})", values, **kw)


def _as_list(v) -> list[str]:
    """A JSON-list string, a '|'-joined string, a list/array or null -> list of strings."""
    v = _maybe_json(v)
    if isinstance(v, (list, tuple, np.ndarray)):
        return [str(x) for x in v if x is not None and str(x)]
    return _split(v)


def _split(s, sep: str = "|") -> list[str]:
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return []
    if isinstance(s, (list, tuple, np.ndarray)):
        return [str(x) for x in s if x is not None and str(x)]
    return [x for x in str(s).split(sep) if x]


# ============================================================================================ resolvers

def _t_condition(t: Trace) -> dict | None:
    cid = t.native
    from ..scoring.opportunity import condition_sets
    sets = condition_sets()
    if cid in sets or "+" in cid:
        members = sets[cid]["members"] if cid in sets else sorted(set(cid.split("+")))
        label = sets[cid]["label"] if cid in sets else " or ".join(members)
        reg = _in("condition_registry", "canonical_condition_id", members, columns=["canonical_condition_id"])
        found = set(reg["canonical_condition_id"]) if len(reg) else set()
        if found != set(members):
            return _unknown(t.object_id, f"condition set {cid!r}: not every member is a registry condition "
                                         f"({sorted(set(members) - found)} missing)")
        t.rows.append({"_table": "condition_sets (measure_it.scoring.opportunity.condition_sets; configs/conditions.yaml)",
                       "condition_set_id": cid, "label": label, "members": members,
                       "set_kind": sets[cid]["set_kind"] if cid in sets else "ad_hoc"})
        t.n_rows = 1
        t.tables.append("condition_registry")
        for m in members:
            t.add_link(_link(f"condition:{m}", "member condition"))
        t.notes.append("A condition set is scored as one burden aggregate of its members (docs/SCORING.md); it is a "
                       "grouping for deployment ranking, not a diagnosis.")
        t.add_sources(sources_for_table("condition_registry"))
        return None
    df = _eq("condition_registry", "canonical_condition_id", cid)
    if df.empty:
        return _unknown(t.object_id, f"{cid!r} is not a canonical_condition_id in condition_registry or a condition set")
    t.add_rows("condition_registry", df)
    maps = _eq("condition_ontology_mappings", "canonical_condition_id", cid,
               columns=["target_id", "target_ontology", "target_label", "predicate_condition", "mapping_source",
                        "in_registry_list"])
    if len(maps):
        by = maps["target_ontology"].value_counts().to_dict()
        t.add_link({"relation": "ontology mappings (condition_ontology_mappings)", "n_rows": int(len(maps)),
                    "by_ontology": to_jsonable(by),
                    "examples": to_jsonable(maps[maps["in_registry_list"].fillna(False).astype(bool)].head(8)
                                            if "in_registry_list" in maps else maps.head(8))})
    if has_table("geo_burden_definitions"):
        bd = _eq("geo_burden_definitions", "condition_id", cid,
                 columns=["geo_level", "burden_evidence_level", "primary_measure_id", "primary_source_resolution",
                          "inherited", "rationale"])
        if len(bd):
            t.add_link({"relation": "burden definition per geography level (geo_burden_definitions; "
                                    "docs/BURDEN_DEFINITIONS.md)", "rows": to_jsonable(bd, max_text=t.max_text)})
    n_ph = LOOKUP.count("condition_phenotypes", "canonical_condition_id = ?", [cid]) if has_table("condition_phenotypes") else 0
    t.add_link({"relation": "Monarch disease-phenotype annotations (condition_phenotypes)", "n_rows": n_ph})
    t.add_link({"relation": "curated names, aliases and flags", "config": "configs/conditions.yaml"})
    t.add_sources(sources_for_table("condition_registry") + ["cms_icd10cm"])
    return None


def _t_measurement(t: Trace) -> dict | None:
    df = _eq("measurement_registry", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, f"{t.native!r} is not a measurement class or bundle id in measurement_registry")
    t.add_rows("measurement_registry", df)
    r = df.iloc[0]
    for m in _as_list(r.get("member_ids")):
        if m != t.native:
            t.add_link(_link(f"measurement:{m}", "member measurement class"))
    for b in _as_list(r.get("bundle_ids")):
        t.add_link(_link(f"measurement:{b}", "bundle containing this class"))
    for f in (_maybe_json(r.get("fda_example_object_ids")) or [])[:5]:
        t.add_link(_link(f, "example FDA record of a mapped product code (deployment-readiness signal)"))
    if has_table("measurement_adapter_registry"):
        ad = _eq("measurement_adapter_registry", "object_id", t.object_id)
        if len(ad):
            keep = [c for c in ("adapter_id", "adapter_name", "status", "adapter_status", "modality", "bundles_naming_this_adapter")
                    if c in ad.columns]
            t.add_link({"relation": "measurement adapter(s) (measurement_adapter_registry)",
                        "rows": to_jsonable(ad[keep], max_text=t.max_text)})
    t.add_link({"relation": "curated measurement classes, patterns and bundles",
                "config": ["configs/measurements.yaml", "configs/relevance.yaml"]})
    t.add_sources(sources_for_table("measurement_registry"))
    t.caveats.append("Registered trials and funded grants show what researchers deploy, not that a measurement works; "
                     "FDA records are deployment-readiness signals.")
    return None


@lru_cache(maxsize=2)
def _ctgov_batch_index(mtime: float) -> dict[str, list[str]]:
    """NCT id -> raw record batch file(s) of data/raw/clinicaltrials_gov/records (built once per process)."""
    out: dict[str, list[str]] = {}
    base = RAW / "clinicaltrials_gov" / "records"
    if not base.exists():
        return out
    pat = re.compile(r'"nctId"\s*:\s*"(NCT\d{8})"')
    for f in sorted(base.glob("batch_*.json")):
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        for n in set(pat.findall(text)):
            out.setdefault(n, []).append(f"records/{f.name}")
    return out


def _t_trial(t: Trace) -> dict | None:
    nct = t.native.upper()
    df = _eq("clinical_trials", "nct_id", nct)
    if df.empty:
        return _unknown(t.object_id, f"{nct} is not in clinical_trials (the target-condition ClinicalTrials.gov pull)")
    t.add_rows("clinical_trials", df, drop=("eligibility_criteria", "arm_groups_json"))
    tc = _eq("trial_conditions", "nct_id", nct)
    for cid, g in tc.groupby("condition_id"):
        lit = bool(g["condition_literal_match"].astype(bool).any()) if "condition_literal_match" in g else None
        t.add_link(_link(f"condition:{cid}", "trial matched this condition (trial_conditions)",
                         literal_match=lit, matched_terms=sorted(set(g["matched_term"].astype(str)))[:5]
                         if "matched_term" in g else None))
    ts = _eq("trial_sites", "nct_id", nct, columns=["facility", "city", "state", "country_group", "county_fips",
                                                    "geocode_method", "source_geographic_resolution"])
    if len(ts):
        us = ts[ts["country_group"] == "US"]
        t.add_link({"relation": "registered sites (trial_sites)", "n_sites": int(len(ts)), "n_us_sites": int(len(us)),
                    "us_sites_example": to_jsonable(us.head(5))})
    if has_table("facility_trials"):
        ft = _eq("facility_trials", "nct_id", nct, columns=["facility_object_id", "measurement_classes"])
        for r in ft.head(10).itertuples():
            t.add_link(_link(r.facility_object_id, "resolved U.S. facility registered as a site of this trial",
                             measurement_classes=_split(r.measurement_classes)))
        if len(ft) > 10:
            t.notes.append(f"{len(ft)} resolved facilities list this trial; 10 shown in lineage.")
    idx = _ctgov_batch_index(_mtime(RAW / "clinicaltrials_gov" / "records"))
    t.add_sources(["clinicaltrials_gov"], idx.get(nct, []))
    t.caveats.append("A registered trial is a research-activity signal, not evidence that a measurement or "
                     "intervention works.")
    return None


def _t_nih(t: Trace) -> dict | None:
    try:
        appl = int(t.native)
    except ValueError:
        return _unknown(t.object_id, f"{t.native!r} is not a numeric RePORTER appl_id")
    df = _eq("nih_projects", "appl_id", appl)
    if df.empty:
        return _unknown(t.object_id, f"appl_id {appl} is not in nih_projects (the target-condition RePORTER pull)")
    t.add_rows("nih_projects", df, drop=("abstract_text", "public_health_relevance", "terms", "pref_terms",
                                         "spending_categories_desc", "ic_fundings", "principal_investigators"))
    pc = _eq("nih_project_conditions", "appl_id", appl)
    for r in pc.itertuples():
        t.add_link(_link(f"condition:{r.condition_id}", "project text mentions this condition (nih_project_conditions)",
                         match_tier=getattr(r, "match_tier", None),
                         likely_false_positive=getattr(r, "likely_false_positive", None)))
    if has_table("facility_nih_projects"):
        fn = _eq("facility_nih_projects", "appl_id", appl, columns=["facility_object_id", "org_name"])
        for r in fn.head(5).itertuples():
            t.add_link(_link(r.facility_object_id, "resolved awardee organisation (facility)", org_name=r.org_name))
    t.add_sources(["nih_reporter"], ["records.jsonl"])
    t.caveats.append("award_amount is a total-cost obligation for one fiscal year, not expenditure; a grant is a "
                     "research-capability signal, not evidence that a measurement works.")
    return None


def _t_npi(t: Trace) -> dict | None:
    df = _eq("providers", "npi", t.native)
    if df.empty:
        return _unknown(t.object_id, f"NPI {t.native} is not in providers (NPPES V.2, active U.S. records)")
    t.add_rows("providers", df, drop=("phone",))
    if has_table("facility_source_links"):
        fl = _eq("facility_source_links", "object_id", t.object_id,
                 columns=["facility_object_id", "match_method", "match_confidence"])
        for r in fl.head(5).itertuples():
            t.add_link(_link(r.facility_object_id, "resolved facility this NPI record was assigned to",
                             match_method=r.match_method, match_confidence=r.match_confidence))
    t.add_sources(["nppes", "nucc_taxonomy"], [str(df["source_version"].iloc[0])])
    t.caveats.append("An NPPES taxonomy is self-reported; it does not mean the provider evaluates or treats any "
                     "invisible illness.")
    return None


def _t_hrsa(t: Trace) -> dict | None:
    df = _eq("facilities__hrsa", "site_bphc_number", t.native)
    if df.empty:
        return _unknown(t.object_id, f"{t.native!r} is not a BPHC site number in facilities__hrsa")
    t.add_rows("facilities__hrsa", df, drop=("site_phone",))
    for fid in df["facility_id"].dropna().unique()[:3]:
        t.add_link(_link(f"facility:{fid}", "resolved facility of this HRSA site"))
    t.add_sources(["hrsa_health_centers"])
    t.caveats.append("HRSA lists health-center program sites; the source has no service-line or clinical-capability "
                     "detail.")
    return None


def _t_facility(t: Trace) -> dict | None:
    fid = t.native
    df = _eq("facilities", "facility_id", fid)
    if df.empty:
        return _unknown(t.object_id, f"{fid!r} is not a facility_id in facilities (clinic_registry / "
                                     "research_site_registry are built from it)")
    t.add_rows("facilities", df)
    reg_in = []
    for tab in ("clinic_registry", "research_site_registry"):
        if has_table(tab) and LOOKUP.count(tab, "facility_id = ?", [fid]):
            reg_in.append(tab)
    t.extra["registries_listing_this_facility"] = reg_in
    if has_table("facility_source_links"):
        fl = _eq("facility_source_links", "facility_id", fid,
                 columns=["record_id", "object_id", "source", "kind", "name", "match_method", "match_confidence"])
        t.add_link({"relation": "source records resolved into this facility (facility_source_links)",
                    "n_records": int(len(fl)), "by_source": to_jsonable(fl["source"].value_counts().to_dict())})
        for r in fl.head(15).itertuples():
            oid = r.object_id if isinstance(r.object_id, str) and r.object_id else None
            t.add_link({"object_id": oid, "relation": f"{r.source} source record ({r.kind})",
                        "record_id": r.record_id, "name": r.name, "match_method": r.match_method,
                        "match_confidence": to_jsonable(r.match_confidence),
                        "note": None if oid else "no CONVENTIONS namespace for this record kind (site key / org key)"})
    for tab, col in (("facility_trials", "nct_id"), ("facility_nih_projects", "appl_id")):
        if has_table(tab):
            d = _eq(tab, "facility_id", fid, columns=[col, "object_id"])
            if len(d):
                t.add_link({"relation": f"{tab}", "n": int(d[col].nunique()),
                            "object_ids": sorted(d["object_id"].dropna().unique().tolist())[:10]})
    srcs = [FACILITY_SOURCE_CODES.get(s, s) for s in _split(df["sources"].iloc[0])]
    t.add_sources(srcs + ["hrsa_hpsa", "usda_rucc"])
    t.caveats.append("A facility with these characteristics may be a viable implementation or study partner; this is "
                     "not a quality ranking. Locations are ZIP/city centroids or HRSA address points.")
    return None


def _t_fda(t: Trace) -> dict | None:
    nat = t.native
    if nat.lower().startswith("code:"):
        code = nat.split(":", 1)[1].upper()
        cls = _eq("fda_device_classification", "product_code", code)
        st = _eq("measurement_regulatory_status", "product_code", code)
        if cls.empty and st.empty:
            return _unknown(t.object_id, f"product code {code} is in neither fda_device_classification nor "
                                         "measurement_regulatory_status")
        t.add_rows("fda_device_classification", cls)
        t.add_rows("measurement_regulatory_status", st)
        for m in sorted(set(st["measurement_id"])) if len(st) else []:
            t.add_link(_link(f"measurement:{m}", "measurement class mapped to this product code "
                                                 "(configs/fda_product_code_map.yaml)"))
        n510 = LOOKUP.count("fda_510k", "product_code = ?", [code])
        t.add_link({"relation": "510(k)/De Novo decisions under this code (fda_510k)", "n": n510})
        t.add_sources(["openfda_device"], [code])
    else:
        key = nat.upper()
        if key.startswith("P"):
            df = _eq("fda_pma", "pma_number", key)
            table = "fda_pma"
        else:
            df = _eq("fda_510k", "clearance_id", key)
            table = "fda_510k"
            if df.empty and has_table("fda_named_device_findings"):
                df = _eq("fda_named_device_findings", "clearance_id", key)
                table = "fda_named_device_findings"
        if df.empty:
            return _unknown(t.object_id, f"{key} is not in fda_510k, fda_pma or fda_named_device_findings")
        t.add_rows(table, df)
        codes = sorted(set(df["product_code"].dropna().astype(str)))
        for c in codes:
            t.add_link(_link(f"fda:code:{c}", "product code of this decision"))
        for m in sorted({x for v in df["mapped_measurement_ids"] for x in (_maybe_json(v) or [])}):
            t.add_link(_link(f"measurement:{m}", "measurement class mapped to the product code"))
        t.add_sources(["openfda_device"], codes)
    t.caveats.append("A regulatory record is a deployment-readiness signal, not evidence that the technology diagnoses "
                     "the target illness.")
    return None


MOL_COLS = ["condition_id", "ontology_id", "evidence_type", "source_evidence_category", "entity_type", "entity_id",
            "entity_label", "source_database", "source_accession", "evidence_direction", "study_population",
            "sample_size", "source_score", "source_score_label", "p_value", "mapped_genes", "reported_trait",
            "data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
            "source_geographic_resolution", "evidence_level", "provenance_notes"]


def _molecular(t: Trace, where: str, params: list, what: str) -> dict | None:
    df = LOOKUP.query("condition_molecular_evidence", where, params, columns=MOL_COLS)
    if df.empty:
        return _unknown(t.object_id, f"no condition_molecular_evidence row for {what}")
    df = df.sort_values(["condition_id", "source_database", "evidence_type"]).reset_index(drop=True)
    t.add_rows("condition_molecular_evidence", df)
    t.extra["summary"] = {
        "n_rows": int(len(df)),
        "by_condition": to_jsonable(df["condition_id"].value_counts().to_dict()),
        "by_source_database": to_jsonable(df["source_database"].value_counts().to_dict()),
        "by_evidence_type": to_jsonable(df["evidence_type"].value_counts().to_dict()),
    }
    for c in sorted(df["condition_id"].dropna().unique()):
        t.add_link(_link(f"condition:{c}", "condition this molecular evidence was retrieved for"))
    hints = sorted(set(df["source_accession"].dropna().astype(str)) | set(df["ontology_id"].dropna().astype(str)
                                                                          .str.replace(":", "_")))[:20]
    for db in sorted(df["source_database"].dropna().unique()):
        t.add_sources([MOLECULAR_SOURCE_DATABASES.get(db, "mapmecfs_nih_pi_mecfs" if "mapMECFS" in db else db)], hints)
    t.caveats.append("Condition-level molecular enrichment from public databases: not measured on any participant, not "
                     "patient multi-omics and not mechanism.")
    return None


def _t_gene(t: Trace) -> dict | None:
    return _molecular(t, "entity_type = 'gene' AND entity_id = ?", [t.native], f"gene {t.native}")


def _t_variant(t: Trace) -> dict | None:
    return _molecular(t, "entity_type = 'variant' AND entity_id = ?", [t.native], f"variant {t.native}")


def _t_gwas_study(t: Trace) -> dict | None:
    r = _molecular(t, "source_database = 'GWAS Catalog' AND (entity_id = ? OR source_accession = ?)",
                   [t.native, t.native], f"GWAS Catalog study {t.native}")
    return r


def _t_geo_series(t: Trace) -> dict | None:
    df = _eq("geo_study_catalog", "gse", t.native.upper())
    if df.empty:
        return _unknown(t.object_id, f"{t.native} is not in geo_study_catalog (curated GEO text queries)")
    t.add_rows("geo_study_catalog", df)
    for c in sorted(df["condition_id"].dropna().unique()):
        t.add_link(_link(f"condition:{c}", "condition whose curated GEO query returned this series"))
    t.add_sources(["ncbi_geo_sra"], [t.native.upper()] + sorted(df["condition_id"].dropna().unique().tolist()))
    t.caveats.append("GEO study metadata only (no expression data analysed); condition-level, not participant-linked.")
    return None


FEATURE_COLS = ["feature_id", "object_id", "geo_id", "geo_level", "geo_name", "condition_id", "burden_value",
                "burden_measure_id", "burden_measure_label", "burden_value_unit", "burden_evidence_level",
                "burden_source_resolution", "burden_inherited", "burden_ci_low", "burden_ci_high", "burden_period",
                "burden_source_name", "burden_row_id", "burden_unknown_reason", "svi_overall", "diagnostic_desert",
                "diagnostic_desert_access_only", "desert_unknown_reason", "relevant_providers_per_100k",
                "trials_in_geo_or_within_50km_n", "nih_core_projects_n", "hrsa_sites_n", "source_name",
                "source_version", "retrieved_at", "source_geographic_resolution", "evidence_level"]
BURDEN_COLS = ["burden_row_id", "geo_id", "geo_level", "geo_name", "condition_id", "measure_id", "measure_label",
               "value", "value_unit", "ci_low", "ci_high", "period", "burden_evidence_level", "measure_role",
               "inherited", "derivation", "proxy_for_condition", "suppressed", "source_table", "source_measure_id",
               "data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
               "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes"]


def _burden_rows(t: Trace, row_ids: list[str], relation: str) -> None:
    ids = [r for r in row_ids if isinstance(r, str) and r]
    if not ids:
        return
    b = _in("geo_condition_burden", "burden_row_id", ids, columns=BURDEN_COLS)
    for r in b.to_dict("records"):
        src_ids = sources_for_table(str(r.get("source_table")))
        t.add_link({"relation": relation, "table": "geo_condition_burden", "row": to_jsonable(r, max_text=t.max_text),
                    "source_ids": src_ids})
        rid = str(r.get("source_record_id") or "")
        hints = [h for h in re.split(r"[:|]", rid) if len(h) >= 5 and re.search(r"\d", h) and "-" in h]
        mm = re.search(r"mmd:([a-z]):(\d{4}):([a-z]):", rid)
        if mm:   # CMS MMD API slice files are named f_<year>_<measure>__...
            hints.append(f"{mm.group(1)}_{mm.group(2)}_{mm.group(3)}__")
        if rid.startswith("brfss_sae:"):   # the small-area model is fitted on the BRFSS 2023 LLCP file
            hints += ["LLCP2023XPT.zip", "codebook23_llcp"]
        t.add_sources(src_ids, hints)


def _t_geo(t: Trace) -> dict | None:
    nat = t.native
    if nat.upper() in ("US", "USA"):
        b = LOOKUP.query("geo_condition_burden", "geo_level = 'national' AND measure_role = 'primary'", [],
                         columns=BURDEN_COLS)
        if b.empty:
            return _unknown(t.object_id, "no national primary burden rows in geo_condition_burden")
        t.add_rows("geo_condition_burden", b)
        for c in sorted(b["condition_id"].unique()):
            t.add_link(_link(f"condition:{c}", "national primary burden row"))
        for st in sorted(b["source_table"].dropna().unique()):
            t.add_sources(sources_for_table(st))
        t.notes.append("geo:US is the national level: the primary national burden rows of every condition are shown.")
        return None
    zm = re.fullmatch(r"zcta:(\d{5})", nat, flags=re.I)
    if zm:
        z = zm.group(1)
        df = LOOKUP.query("geo_context_zcta", "geo_id = ?", [z]) if has_table("geo_context_zcta") else pd.DataFrame()
        if df.empty:
            return _unknown(t.object_id, f"ZCTA {z} is not in geo_context_zcta")
        t.add_rows("geo_context_zcta", df)
        t.add_sources(sources_for_table("geo_context_zcta"))
        t.caveats.append("Ecological, place-level data; nothing describes an individual.")
        return None
    if not re.fullmatch(r"\d{2}|\d{5}", nat):
        return _unknown(t.object_id, f"{nat!r} is not a 2-digit state FIPS, 5-digit county FIPS, US or zcta:<ZCTA>")
    g = _eq("geographies", "geo_id", nat)
    if g.empty:
        return _unknown(t.object_id, f"no state/county with FIPS {nat} in geographies (2024 vintage)")
    t.add_rows("geographies", g)
    ctx = _eq("geo_context", "geo_id", nat) if has_table("geo_context") else pd.DataFrame()
    t.add_rows("geo_context", ctx)
    feats = _eq("geo_condition_features", "geo_id", nat, columns=FEATURE_COLS)
    for r in feats.sort_values("condition_id").to_dict("records"):
        t.add_link({"relation": "condition feature row (geo_condition_features)", "feature_id": r.get("feature_id"),
                    "condition_object_id": f"condition:{r.get('condition_id')}",
                    "burden_value": to_jsonable(r.get("burden_value")),
                    "burden_evidence_level": r.get("burden_evidence_level"),
                    "burden_source_resolution": r.get("burden_source_resolution"),
                    "burden_inherited": to_jsonable(r.get("burden_inherited")),
                    "burden_row_id": to_jsonable(r.get("burden_row_id")),
                    "burden_unknown_reason": to_jsonable(r.get("burden_unknown_reason"), max_text=t.max_text),
                    "diagnostic_desert": to_jsonable(r.get("diagnostic_desert"))})
    t.add_sources(["census_geography"] + (sources_for_table("geo_context") if len(ctx) else []))
    t.caveats.append("Ecological, place-level data: nothing here describes an individual or where any participant "
                     "lives. A state value attached to a county is inherited context, not county prevalence.")
    return None


def _isnull(v) -> bool:
    v = to_jsonable(v)
    return v is None or (isinstance(v, float) and math.isnan(v))


def null_score_note(r) -> str | None:
    """Why a stored opportunity row has a null composite and/or rank (None when both are set)."""
    comp, rank = r.get("composite_equal"), r.get("rank_equal")
    if not _isnull(comp) and not _isnull(rank):
        return None
    if bool(to_jsonable(r.get("burden_incomplete")) or False):
        members = ", ".join(_split(to_jsonable(r.get("burden_members_missing")))) or UNKNOWN
        return ("composite and rank are null because the burden is incomplete (member(s) with a defined burden "
                f"measure but no value in this region: {members}). {to_jsonable(r.get('burden_incomplete_reason'))} "
                "rank_deployment_opportunities lists the region under excluded_incomplete_burden; its vulnerability, "
                "desert, clinic_capacity and research_readiness components are kept in this row.")
    if _isnull(comp):
        return ("composite is null: no component of the weighted mean has a value for this region (see the component "
                "percentiles); it cannot be ranked.")
    if bool(to_jsonable(r.get("small_population_flag")) or False):
        pop = to_jsonable(r.get("population_total"))
        incl = to_jsonable(r.get("rank_equal_incl_small"))
        return (f"composite {float(comp):.3f} is computed but the county is not ranked: population "
                f"{'UNKNOWN' if _isnull(pop) else f'{int(pop):,}'} is below the ranking minimum (small denominators). "
                f"Rank among all counties including small ones (rank_equal_incl_small): "
                f"{'UNKNOWN' if _isnull(incl) else int(incl)}.")
    return (f"composite {float(comp):.3f} is computed but the region is not ranked (rank_eligible = "
            f"{to_jsonable(r.get('rank_eligible'))}).")


def _t_opportunity(t: Trace) -> dict | None:
    parts = t.native.split("|")
    if len(parts) != 3:
        return _unknown(t.object_id, "opportunity ids are opportunity:<condition>|<measurement>|<geo_id>")
    cid, mid, gid = parts
    df = _eq("deployment_opportunities", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, "not in deployment_opportunities (only the precomputed condition x measurement x "
                                     "level combinations are stored; others are computed on the fly by "
                                     "rank_deployment_opportunities and have no stored row)")
    t.add_rows("deployment_opportunities", df)
    r = df.iloc[0]
    from ..scoring.opportunity import COMPONENTS
    comps = {}
    for c in COMPONENTS:
        comps[c] = {"raw": to_jsonable(r.get(c if c not in ("burden", "vulnerability") else
                                             ("burden_value" if c == "burden" else "svi_overall"))),
                    "percentile": to_jsonable(r.get(f"{c}_pct"))}
    incomplete = bool(to_jsonable(r.get("burden_incomplete")) or False)
    missing = _split(to_jsonable(r.get("burden_members_missing")))
    comps["burden"].update({"evidence_level": r.get("burden_evidence_level"),
                            "source_resolution": to_jsonable(r.get("burden_source_resolution")),
                            "inherited": to_jsonable(r.get("burden_inherited")),
                            "modelled_small_area": ("small-area" in str(r.get("burden_source_name") or "")
                                                    or "small-area" in str(r.get("member_components") or "")),
                            "members_contributing": to_jsonable(r.get("burden_members_contributing")),
                            "members_excluded_level_D": to_jsonable(r.get("burden_members_excluded_level_D")),
                            "incomplete": incomplete,
                            "incomplete_reason": to_jsonable(r.get("burden_incomplete_reason")),
                            "members_missing": missing,
                            "missing_source_reason": to_jsonable(r.get("burden_missing_source_reason")),
                            "unknown_reason": to_jsonable(r.get("burden_unknown_reason"))})
    comps["technology_saturation"] = {"raw": to_jsonable(r.get("technology_saturation")),
                                      "percentile": to_jsonable(r.get("technology_saturation_pct")),
                                      "note": "reported only; not in the default formula"}
    if r.get("condition_kind") == "condition_set":
        comps["burden"]["raw_note"] = ("a condition set has no single raw burden value: its burden is the percentile of "
                                       "the mean of member burden percentiles; the member values are the "
                                       "geo_condition_features / geo_condition_burden rows in lineage")
    t.extra["components"] = comps
    t.extra["composite_equal"] = to_jsonable(r.get("composite_equal"))
    t.extra["rank_equal"] = to_jsonable(r.get("rank_equal"))
    # incomplete burden (docs/ANALYSIS_PLAN_SCORING.md deviation 7): why composite / rank are null, stated, not implied
    t.extra["burden_incomplete"] = incomplete
    t.extra["burden_incomplete_reason"] = to_jsonable(r.get("burden_incomplete_reason"))
    t.extra["burden_members_missing"] = missing
    note = null_score_note(r)
    if note:
        t.extra["composite_rank_note"] = note
        t.notes.append(note)
    t.add_link(_link(f"geo:{gid}", "geography"))
    t.add_link(_link(f"measurement:{mid}", "measurement"))
    members = _as_list(r.get("condition_members"))
    if not members:
        members = [cid]
    for m in members:
        t.add_link(_link(f"condition:{m}", "member condition"))
    level = r.get("geo_level")
    feats = LOOKUP.query("geo_condition_features", "geo_id = ? AND geo_level = ? AND condition_id IN (" +
                         ", ".join("?" for _ in members) + ")", [gid, level] + members, columns=FEATURE_COLS)
    burden_ids = []
    for f in feats.to_dict("records"):
        t.add_link({"relation": "burden / vulnerability / desert inputs (geo_condition_features)",
                    "table": "geo_condition_features", "row": to_jsonable(f, max_text=t.max_text),
                    "source_ids": sources_for_table("geo_condition_features")})
        burden_ids.append(f.get("burden_row_id"))
    _burden_rows(t, burden_ids, "burden row behind the burden component (geo_condition_burden)")
    t.add_link({"relation": "clinic_capacity / research_readiness inputs",
                "tables": ["providers (nppes)", "facilities", "facility_trials (clinicaltrials_gov)",
                           "facility_nih_projects (nih_reporter)"],
                "values": to_jsonable({k: r.get(k) for k in (
                    "cc_impl_providers_in_geo_per_100k", "cc_impl_providers_within_50km_per_100k",
                    "cc_impl_facilities_in_geo_per_100k", "cc_impl_facilities_within_50km_per_100k",
                    "rr_condition_trials_pool_n", "rr_condition_nih_core_pool_n", "rr_tech_trials_pool_n")}),
                "note": "object ids of the trials / grants / facilities in this pool: rank_deployment_opportunities "
                        "research_evidence and find_candidate_clinics"})
    t.add_sources(["cdc_svi"], ["PUERTORICO_COUNTY" if gid.startswith("72") else "US_COUNTY"])
    # ACS B01003 total population: rank eligibility (population floor) and the per-100k clinic-capacity denominators
    t.add_sources(["census_acs"], ["acsdt5y2024-b01003.dat"])
    t.add_sources(["nppes"], ["NPPES_Data_Dissemination"])
    t.add_sources(["hrsa_health_centers"], ["LookAlike_Sites.csv"])
    t.add_sources(["clinicaltrials_gov", "nih_reporter"])
    t.caveats.append("A ranked region is a candidate deployment opportunity for pilot evaluation, not a validated "
                     "diagnostic pathway; all joins are ecological.")
    return None


def _dataset_source(dataset_id: str) -> list[str]:
    """Source ids behind a person-level dataset id (phenotype_signatures.dataset_id / participant namespace)."""
    d = str(dataset_id)
    if d.startswith("nhanes0306"):
        return ["nhanes_2003_2006"]
    if d.startswith("nhanes"):                      # nhanes, nhanes_labs (NHANES 2011-2014 labs vs labels)
        return ["nhanes_2011_2014"]
    if d.startswith("stanford"):
        return ["stanford_longcovid_wearables"]
    if d.startswith("mapmecfs"):
        return ["mapmecfs_nih_pi_mecfs"]
    if re.match(r"^geo_(GSE\d+|cohort)", d):          # geo_GSE270045: one GEO case/control cohort
        return ["geo_cohorts"]
    if d in registry()["by_id"]:                    # datasets whose id is their source id (lab / device releases)
        return [d]
    return []


def _t_signature(t: Trace) -> dict | None:
    df = _eq("phenotype_signatures", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, "not in phenotype_signatures (object ids are signature:<dataset>|<phenotype>|"
                                     "<feature>)")
    t.add_rows("phenotype_signatures", df)
    r = df.iloc[0]
    ds = str(r.get("dataset_id"))
    t.extra["label_basis"] = {k: to_jsonable(r.get(k), max_text=t.max_text) for k in
                              ("phenotype_id", "phenotype_label", "phenotype_definition", "label_basis", "is_proxy",
                               "population", "n_cases", "n_controls", "evidence_level", "caveats")}
    if has_table("digital_phenotype_datasets"):
        dd = _eq("digital_phenotype_datasets", "dataset_id", ds, columns=["dataset_id", "representation_label",
                                                                          "n_participants_population",
                                                                          "population_definition"])
        if len(dd):
            t.add_link({"relation": "dataset (digital_phenotype_datasets)", "rows": to_jsonable(dd, max_text=t.max_text)})
    t.add_link({"relation": "dataset", "dataset_id": ds, "source_ids": _dataset_source(ds),
                "note": "aggregate statistics computed from person-level rows; no participant row is exposed"})
    if isinstance(r.get("condition_id"), str) and r.get("condition_id"):
        t.add_link(_link(f"condition:{r['condition_id']}", "condition the phenotype label refers to (proxy flag in "
                                                           "label_basis)"))
    # raw files named in the row's version / provenance text (e.g. DEMO_G.xpt, PAXDAY_H.xpt)
    t.add_sources(_dataset_source(ds), [str(r.get("source_version") or ""), str(r.get("provenance_notes") or "")])
    t.caveats.append("Group-level association in a public person-level dataset, never an individual finding or "
                     "diagnostic accuracy; proxy labels are not the named condition.")
    return None


def _t_measurement_evidence(t: Trace) -> dict | None:
    df = _eq("condition_measurement_evidence", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, "not in condition_measurement_evidence (<condition_id>|<measurement_id>)")
    t.add_rows("condition_measurement_evidence", df, drop=("example_snippets",))
    r = df.iloc[0]
    t.add_link(_link(r["condition_object_id"], "condition"))
    t.add_link(_link(r["measurement_object_id"], "measurement"))
    for k, rel in (("example_trial_object_ids", "example trial describing objective use"),
                   ("example_grant_object_ids", "example NIH project mentioning it"),
                   ("example_fda_object_ids", "example FDA record (deployment-readiness signal)"),
                   ("molecular_object_ids", "molecular context row")):
        for oid in (_maybe_json(r.get(k)) or [])[:5]:
            t.add_link(_link(oid, rel))
    sig = r.get("phenotype_signal_strength__signal_object_id")
    if isinstance(sig, str) and sig:
        t.add_link(_link(sig, "public person-level signal result behind phenotype_signal_strength"))
    t.add_sources(sources_for_table("condition_measurement_evidence"))
    t.caveats.append("Registered trials and funded grants show what researchers deploy, not that a measurement works "
                     "or detects the condition.")
    return None


def _t_measurement_signal(t: Trace) -> dict | None:
    df = _eq("measurement_phenotype_signal", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, "not in measurement_phenotype_signal")
    t.add_rows("measurement_phenotype_signal", df)
    r = df.iloc[0]
    t.add_link(_link(f"measurement:{r['measurement_id']}", "measurement class"))
    for oid in (_maybe_json(r.get("signature_object_ids")) or [])[:10]:
        t.add_link(_link(oid, "per-feature signature row"))
    t.add_sources(_dataset_source(r.get("dataset_id")))
    t.caveats.append("Group-level association in public data; label_class states how strong the label is.")
    return None


def _t_phenotype_evidence(t: Trace) -> dict | None:
    df = _eq("phase3_measurement_evidence", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, "not in phase3_measurement_evidence (<phenotype>|<measurement>)")
    t.add_rows("phase3_measurement_evidence", df)
    r = df.iloc[0]
    t.add_link(_link(r["measurement_object_id"], "measurement"))
    sig = r.get("phenotype_signal_strength__signal_object_id")
    if isinstance(sig, str) and sig:
        t.add_link(_link(sig, "public person-level signal result"))
    for k in ("ctgov__example_trial_object_ids", "nih__example_grant_object_ids", "fda__example_object_ids"):
        for oid in (_maybe_json(r.get(k)) or [])[:3]:
            t.add_link(_link(oid, k.split("__")[1].replace("_", " ")))
    t.add_sources(sources_for_table("phase3_measurement_evidence"))
    t.caveats.append("Five separate dimensions; no composite score. A phenotype-level signal on a reference label is "
                     "not a result for any specific condition unless stated.")
    return None


def _t_molbio(t: Trace) -> dict | None:
    df = _eq("measurable_biology", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, "not in measurable_biology (<condition>|<system>|<measurement_class>[|<catalogue>])")
    t.add_rows("measurable_biology", df)
    r = df.iloc[0]
    t.add_link(_link(f"condition:{r['condition_id']}", "condition (condition -> system link is data-derived)"))
    t.add_link(_link(f"measurement:{r['measurement_class']}", "measurement class (system -> class link is curated)"))
    for rid in re.findall(r"R-HSA-\d+", str(r.get("top_pathways") or ""))[:5]:
        t.add_link(_link(f"reactome:{rid}", "enriched Reactome pathway behind the system link"))
    t.add_sources(sources_for_table("measurable_biology"))
    t.caveats.append("Condition-level molecular enrichment; candidate measurable phenotype only, not evidence that the "
                     "measurement detects the condition.")
    return None


@lru_cache(maxsize=2)
def _reactome_ora(mtime: float) -> pd.DataFrame:
    return pd.read_csv(TABLES / "test3_reactome_ora.csv", dtype=str)


def _t_reactome(t: Trace) -> dict | None:
    p = TABLES / "test3_reactome_ora.csv"
    if not p.exists():
        return _unknown(t.object_id, "results/tables/test3_reactome_ora.csv not built (measure_it.omics.coherence)")
    ora = _reactome_ora(_mtime(p))
    df = ora[ora["reactome_id"] == t.native]
    if df.empty:
        return _unknown(t.object_id, f"{t.native} has no row (overlap > 0) in results/tables/test3_reactome_ora.csv")
    num = [c for c in ("pathway_size", "overlap", "p_value", "fdr", "log2_fold") if c in df.columns]
    df = df.copy()
    for c in num:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.sort_values(["fdr", "condition_id", "source"])
    t.add_rows("test3_reactome_ora", df)
    t.extra["summary"] = {"n_rows": int(len(df)), "n_enriched": int((df["enriched"].astype(str) == "True").sum()),
                          "conditions": sorted(df["condition_id"].unique().tolist())}
    for c in sorted(df["condition_id"].unique()):
        t.add_link(_link(f"condition:{c}", "condition whose gene list overlaps this pathway"))
    t.add_sources(DERIVED_TABLES["test3_reactome_ora"]["sources"][:2], [t.native])
    t.caveats.append("Over-representation analysis on condition-level gene lists; not mechanism and not person-linked.")
    return None


def _t_measurement_performance(t: Trace) -> dict | None:
    df = _eq("measurement_performance", "object_id", t.object_id)
    if df.empty:
        return _unknown(t.object_id, "not in measurement_performance (<condition_id>|<measurement_id>)")
    t.add_rows("measurement_performance", df)
    r = df.iloc[0]
    t.add_link(_link(f"condition:{r['condition_id']}", "condition or condition set"))
    t.add_link(_link(f"measurement:{r['measurement_id']}", "measurement class or bundle"))
    for oid in (_maybe_json(r.get("source_object_ids")) or [])[:10]:
        t.add_link(_link(oid, "result behind the scored performance record"
                         + (" (published claim, not reproduced)" if str(oid).startswith("published_evidence:") else "")))
    rec = _maybe_json(r.get("records")) or []
    ids = [x.get("record_id") for x in rec if isinstance(x, dict) and x.get("record_id")]
    if ids:
        t.add_rows("measurement_performance_records",
                   _in("measurement_performance_records", "record_id", ids))
    t.add_sources(sources_for_table("measurement_performance"))
    t.caveats.append("Performance against healthy or recovered controls overstates real-world performance; tier 3 rows "
                     "are published claims not reproduced by this project; UNKNOWN is never imputed from another "
                     "condition or measurement.")
    return None


RESOLVERS = {
    "condition": _t_condition, "measurement": _t_measurement, "trial": _t_trial, "nih": _t_nih, "npi": _t_npi,
    "hrsa_site": _t_hrsa, "facility": _t_facility, "fda": _t_fda, "gene": _t_gene, "variant": _t_variant,
    "gwas_study": _t_gwas_study, "geo_series": _t_geo_series, "geo": _t_geo, "opportunity": _t_opportunity,
    "signature": _t_signature, "measurement_evidence": _t_measurement_evidence,
    "measurement_signal": _t_measurement_signal, "phenotype_evidence": _t_phenotype_evidence, "molbio": _t_molbio,
    "reactome": _t_reactome, "measurement_performance": _t_measurement_performance,
}
assert set(RESOLVERS) == set(NAMESPACES)

OBJECT_ID_RE = re.compile(r"^([a-z_]+):(.+)$")


def parse_object_id(object_id: str) -> tuple[str | None, str | None]:
    m = OBJECT_ID_RE.match(str(object_id or "").strip())
    if not m:
        return None, None
    return m.group(1), m.group(2).strip()


def trace_evidence(object_id: str, *, max_rows: int = MAX_ROWS, max_text: int = MAX_TEXT) -> dict:
    """Resolve `<namespace>:<native id>` to rows, provenance, SOURCE_REGISTRY entries, DATA_AUDIT.md paths, raw files
    and one level of lineage. Returns status UNKNOWN (with a reason) for unknown namespaces or ids, and a refusal for
    person-level namespaces."""
    oid = str(object_id or "").strip()
    if is_person_level_id(oid):
        return _unknown(oid, PERSON_REFUSAL, refused=True, namespace=_PREFIX_RE.match(oid).group(1).lower())
    ns, nat = parse_object_id(oid)
    if ns is None:
        return _unknown(oid, "not an object id: expected '<namespace>:<native id>' (docs/CONVENTIONS.md section 7)",
                        known_namespaces=sorted(NAMESPACES))
    if ns in PERSON_LEVEL_NAMESPACES:
        return _unknown(oid, PERSON_REFUSAL, refused=True, namespace=ns)
    if ns not in RESOLVERS:
        return _unknown(oid, f"unknown namespace {ns!r}", known_namespaces=sorted(NAMESPACES), namespace=ns)
    t = Trace(oid, ns, nat, max(1, int(max_rows)), max(40, int(max_text)))
    try:
        out = RESOLVERS[ns](t)
    except FileNotFoundError as e:
        return _unknown(oid, f"table not built: {e}", namespace=ns)
    if out is not None:
        out.setdefault("namespace", ns)
        out.setdefault("native_id", nat)
        out.setdefault("primary_table", NAMESPACES[ns][0])
        return to_jsonable(out)
    res = t.result()
    res["primary_table"] = NAMESPACES[ns][0]
    return to_jsonable(res)


def namespaces_doc() -> list[dict]:
    return [{"namespace": k, "table": v[0], "native_id": v[1]} for k, v in NAMESPACES.items()] + [
        {"namespace": "participant", "table": "(person-level)", "native_id": "refused: person-level rows are summarised only"}]


def docs_path(name: str) -> str:
    return _rel(DOCS / name)
