"""Read-only data bundle: everything the API, dashboard and MCP server read at run time, in one archive.

    uv run measure-it bundle                          # dist/measure-it-data-<UTC date>.tar.zst (+ .manifest.json, .sha256)
    uv run measure-it bundle --no-duckdb              # without data/processed/measure_it_public.duckdb (-1.4 GB)
    uv run measure-it bundle --unpack dist/bundle     # build, then unpack into dist/bundle (what docker-compose mounts)
    uv run measure-it bundle --verify dist/bundle     # check an unpacked bundle against its BUNDLE_MANIFEST.json

What goes in (paths relative to the project root, unpacked in place at /app in the container):

* ``data/processed``: every ``*.parquet``, ``*.geoparquet`` and ``*.meta.json`` sidecar, and the DuckDB (optional);
* ``data/raw/<source>/``: only the provenance sidecars that ``trace_evidence`` and the source summaries read
  (MANIFEST.json, DATA_AUDIT.md, registry_entry.yaml, QUERY_LOG.json, queries.json; git keeps only their digests), plus
  ``data/raw/data_gov_catalog`` (the stored Data.gov responses) and ``data/raw/clinicaltrials_gov/records`` (the
  registry batches trace_evidence points a trial id to). No other raw file: raw data stay out of the bundle;
* ``data/_http_cache``: only the cached Data.gov catalog responses (hosts in HTTP_CACHE_HOSTS), so that
  ``search_us_open_data`` answers its cached queries offline;
* ``results``: the ``*.md`` reports, ``*.json`` examples and ``results/tables`` (figures, maps and pipeline run logs
  are left out);
* ``configs`` and ``SOURCE_REGISTRY.yaml``: the configuration the tables were built with.

``BUNDLE_MANIFEST.json`` (the last archive member, also written next to the archive) lists every file with its size,
SHA-256 (computed from the bytes written to the archive) and row count (Parquet metadata, CSV records, DuckDB tables),
the git commit and dirty state of the code, every source's version and retrieval date from SOURCE_REGISTRY.yaml, and
what was left out. The two lab records that appear to hold personal identifiers
(``measure_it.labs.lab_dataset_discovery.PRIVACY_EXCLUDED``) were never downloaded; the builder refuses to archive any
path that names them.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import platform
import shutil
import subprocess
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import yaml

from ..config import CONFIGS, DUCKDB_PATH, HTTP_CACHE, PROCESSED, PROJECT_ROOT, RAW, RESULTS, SOURCE_REGISTRY_PATH

BUNDLE_FORMAT = 1
MANIFEST_NAME = "BUNDLE_MANIFEST.json"
DIST = PROJECT_ROOT / "dist"
RAW_SIDECARS = ("MANIFEST.json", "DATA_AUDIT.md", "registry_entry.yaml", "QUERY_LOG.json", "queries.json")
RAW_WHOLE_DIRS = ("data_gov_catalog",)                       # search_us_open_data: stored responses + openapi
RAW_SUBDIRS = (("clinicaltrials_gov", "records"),)           # trace_evidence: NCT id -> raw registry batch file
HTTP_CACHE_HOSTS = ("catalog.data.gov", "api.gsa.gov")       # the Data.gov catalog API (both endpoints)
RESULTS_TOP_SUFFIXES = (".md", ".json")
EXCLUDED = {
    "data/raw": "raw source files (40 GB) except the provenance sidecars and the two directories listed in included",
    "data/interim": "derived caches rebuilt by the pipeline",
    "data/_http_cache": "every cached response except the Data.gov catalog hosts",
    "results/figures, results/maps, results/pipeline_runs": "figures, map files and pipeline run logs",
}


def _privacy_excluded() -> set[str]:
    try:
        from ..labs.lab_dataset_discovery import PRIVACY_EXCLUDED
        return set(PRIVACY_EXCLUDED)
    except Exception:  # noqa: BLE001 - keep the bundle usable if the module moves; the ids are also listed here
        return {"postcovid_gpcr_aab_seibert2026", "figshare_gpcr_autoab_postcovid_EXCLUDE", "endo_vitd_zinc_indonesia"}


# ------------------------------------------------------------------------------------------------ file selection

def _rel(p: Path) -> str:
    return p.relative_to(PROJECT_ROOT).as_posix()


def _walk(base: Path) -> Iterable[Path]:
    if not base.exists():
        return []
    return sorted(p for p in base.rglob("*") if p.is_file() and not p.name.startswith(".tmp-")
                  and p.name != "MANIFEST.lock")


def _http_cache_entries(hosts: tuple[str, ...] = HTTP_CACHE_HOSTS) -> list[Path]:
    from urllib.parse import urlparse
    out = []
    for meta in sorted(HTTP_CACHE.glob("*/*.meta.json")):
        try:
            host = urlparse(json.loads(meta.read_text()).get("url", "")).netloc
        except (OSError, ValueError):
            continue
        body = meta.with_name(meta.name[: -len(".meta.json")] + ".body")
        if host in hosts and body.exists():
            out += [meta, body]
    return out


def collect_files(include_duckdb: bool = True) -> list[Path]:
    """Absolute paths of every file that goes into the bundle (sorted, unique)."""
    files: list[Path] = []
    for p in sorted(PROCESSED.iterdir()) if PROCESSED.exists() else []:
        if p.is_file() and (p.suffix in (".parquet", ".geoparquet") or p.name.endswith(".meta.json")):
            files.append(p)
    if include_duckdb and DUCKDB_PATH.exists():
        files.append(DUCKDB_PATH)
    for src in sorted(RAW.iterdir()) if RAW.exists() else []:
        if not src.is_dir():
            continue
        if src.name in RAW_WHOLE_DIRS:
            files += list(_walk(src))
            continue
        files += [src / n for n in RAW_SIDECARS if (src / n).is_file()]
        for sid, sub in RAW_SUBDIRS:
            if src.name == sid:
                files += list(_walk(src / sub))
    files += _http_cache_entries()
    for p in sorted(RESULTS.iterdir()) if RESULTS.exists() else []:
        if p.is_file() and p.suffix in RESULTS_TOP_SUFFIXES:
            files.append(p)
    files += list(_walk(RESULTS / "tables"))
    files += list(_walk(CONFIGS))
    if SOURCE_REGISTRY_PATH.exists():
        files.append(SOURCE_REGISTRY_PATH)
    seen, out = set(), []
    for f in files:
        if f not in seen and "__pycache__" not in f.parts:
            seen.add(f)
            out.append(f)
    bad = [_rel(f) for f in out for x in _privacy_excluded() if x.lower() in _rel(f).lower()]
    if bad:
        raise RuntimeError(f"refusing to bundle files that name a privacy-excluded record: {bad[:5]}")
    return out


# ------------------------------------------------------------------------------------------------ row counts

def row_count(path: Path) -> int | dict | None:
    """Parquet / GeoParquet: rows from the file metadata; CSV: records (header excluded); DuckDB: rows per table."""
    try:
        if path.suffix in (".parquet", ".geoparquet"):
            import pyarrow.parquet as pq
            return int(pq.ParquetFile(path).metadata.num_rows)
        if path.suffix == ".csv":
            with open(path, newline="", encoding="utf-8", errors="replace") as fh:
                return max(sum(1 for _ in csv.reader(fh)) - 1, 0)
        if path.suffix == ".duckdb":
            import duckdb
            con = duckdb.connect(str(path), read_only=True)
            try:
                names = [r[0] for r in con.execute(
                    "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main' "
                    "AND table_type = 'BASE TABLE' ORDER BY 1").fetchall()]
                return {n: int(con.execute(f'SELECT count(*) FROM main."{n}"').fetchone()[0]) for n in names}
            finally:
                con.close()
    except Exception:  # noqa: BLE001 - an unreadable file gets no count, never a guessed one
        return None
    return None


# ------------------------------------------------------------------------------------------------ provenance

def git_state() -> dict:
    def run(*args: str) -> str | None:
        try:
            r = subprocess.run(["git", *args], cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=30)
            return r.stdout.strip() if r.returncode == 0 else None
        except (OSError, subprocess.SubprocessError):
            return None
    commit = run("rev-parse", "HEAD")
    status = run("status", "--porcelain", "--untracked-files=no")
    return {"commit": commit or "UNKNOWN / NOT AVAILABLE (not a git checkout)",
            "dirty_tracked_files": None if status is None else len([s for s in status.splitlines() if s.strip()]),
            "describe": run("describe", "--always", "--dirty")}


def source_versions() -> list[dict]:
    try:
        reg = yaml.safe_load(SOURCE_REGISTRY_PATH.read_text()) or {}
    except (OSError, yaml.YAMLError):
        return []
    keys = ("source_id", "name", "status", "source_version", "retrieved_at", "update_date", "data_layer")
    return [{k: (str(s.get(k)) if s.get(k) is not None else None) for k in keys} for s in reg.get("sources", [])]


# ------------------------------------------------------------------------------------------------ archive

class _HashingReader(io.RawIOBase):
    """File reader that hashes exactly the bytes tarfile copies into the archive."""

    def __init__(self, fh) -> None:
        self.fh, self.sha = fh, hashlib.sha256()

    def readable(self) -> bool:
        return True

    def read(self, n: int = -1) -> bytes:
        b = self.fh.read(n)
        self.sha.update(b)
        return b


def _open_sink(path: Path, compression: str):
    """(file object to write the tar stream to, closer) for zstd (external `zstd -T0`) or gzip (Python)."""
    if compression == "zst":
        exe = shutil.which("zstd")
        if exe is None:
            raise RuntimeError("zstd not found; use --compression gz")
        proc = subprocess.Popen([exe, "-T0", "-q", "-f", "-3", "-o", str(path)], stdin=subprocess.PIPE)

        def close() -> None:
            proc.stdin.close()
            if proc.wait() != 0:
                raise RuntimeError(f"zstd exited with {proc.returncode}")
        return proc.stdin, close
    import gzip
    gz = gzip.open(path, "wb", compresslevel=6)
    return gz, gz.close


def build_bundle(out_dir: Path | None = None, *, include_duckdb: bool = True, compression: str | None = None,
                 date: str | None = None, echo: Callable[[str], None] = print) -> dict:
    """Write dist/measure-it-data-<date>.tar.zst (or .tar.gz) with BUNDLE_MANIFEST.json; return a summary."""
    out_dir = Path(out_dir or DIST)
    out_dir.mkdir(parents=True, exist_ok=True)
    compression = compression or ("zst" if shutil.which("zstd") else "gz")
    date = date or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    archive = out_dir / f"measure-it-data-{date}.tar.{compression}"
    tmp = archive.with_name(archive.name + ".partial")
    t0 = time.time()
    files = collect_files(include_duckdb=include_duckdb)
    if not any(f.parent == PROCESSED and f.suffix == ".parquet" for f in files):
        raise RuntimeError("data/processed holds no Parquet table: run the pipeline first (README, Quickstart)")
    echo(f"bundle: {len(files)} files, {sum(f.stat().st_size for f in files) / 1e9:.2f} GB before compression")
    entries, changed = [], []
    sink, close = _open_sink(tmp, compression)
    try:
        with tarfile.open(fileobj=sink, mode="w|", format=tarfile.PAX_FORMAT) as tar:
            # directory entries first, with their source mtimes: the tools' data stamp (newest mtime in
            # data/processed, configs, results/tables) then matches the build, not the unpack time
            for d in sorted({p for f in files for p in f.relative_to(PROJECT_ROOT).parents if str(p) != "."}):
                info = tarfile.TarInfo(d.as_posix())
                info.type, info.mode = tarfile.DIRTYPE, 0o755
                info.mtime = int((PROJECT_ROOT / d).stat().st_mtime)
                info.uname = info.gname = "measure-it"
                tar.addfile(info)
            for i, f in enumerate(files, 1):
                st = f.stat()
                info = tarfile.TarInfo(_rel(f))
                info.size, info.mtime, info.mode = st.st_size, int(st.st_mtime), 0o644
                info.uname = info.gname = "measure-it"
                with open(f, "rb") as fh:
                    reader = _HashingReader(fh)
                    tar.addfile(info, reader)
                if f.stat().st_mtime != st.st_mtime:
                    changed.append(_rel(f))
                entries.append({"path": _rel(f), "bytes": st.st_size, "sha256": reader.sha.hexdigest(),
                                "mtime_utc": datetime.fromtimestamp(st.st_mtime, timezone.utc)
                                .replace(microsecond=0).isoformat(), "rows": row_count(f)})
                if i % 200 == 0:
                    echo(f"bundle: {i}/{len(files)} files")
            manifest = _manifest(entries, include_duckdb, compression, changed)
            raw = json.dumps(manifest, indent=1, default=str).encode()
            info = tarfile.TarInfo(MANIFEST_NAME)
            info.size, info.mtime, info.mode = len(raw), int(time.time()), 0o644
            tar.addfile(info, io.BytesIO(raw))
    finally:
        close()
    os.replace(tmp, archive)
    size = archive.stat().st_size
    sha = _sha256_file(archive)
    (out_dir / f"{archive.name}.sha256").write_text(f"{sha}  {archive.name}\n")
    man_path = out_dir / f"measure-it-data-{date}.manifest.json"
    man_path.write_bytes(raw)
    summary = {"archive": str(archive), "archive_bytes": size, "archive_sha256": sha, "manifest": str(man_path),
               "n_files": len(entries), "uncompressed_bytes": manifest["totals"]["bytes"],
               "files_changed_during_build": changed, "seconds": round(time.time() - t0, 1)}
    echo(f"bundle: {archive} {size / 1e9:.2f} GB ({manifest['totals']['bytes'] / 1e9:.2f} GB unpacked, "
         f"{len(entries)} files, {summary['seconds']} s)")
    if changed:
        echo(f"bundle: WARNING {len(changed)} file(s) changed while archiving (a pipeline was running?): "
             f"{changed[:5]}; rebuild the bundle when the pipeline has finished")
    return summary


def _manifest(entries: list[dict], include_duckdb: bool, compression: str, changed: list[str]) -> dict:
    groups: dict[str, dict] = {}
    for e in entries:
        parts = e["path"].split("/")
        g = "/".join(parts[:2]) if parts[0] in ("data", "results") and len(parts) > 2 else parts[0]
        d = groups.setdefault(g, {"files": 0, "bytes": 0})
        d["files"] += 1
        d["bytes"] += e["bytes"]
    tables = {Path(e["path"]).name.split(".")[0]: e["rows"] for e in entries
              if e["path"].startswith("data/processed/") and e["path"].endswith((".parquet", ".geoparquet"))}
    return {
        "bundle_format": BUNDLE_FORMAT,
        "name": "measure-it-public data bundle",
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "git": git_state(),
        "build_host": {"platform": platform.platform(), "machine": platform.machine(),
                       "python": platform.python_version()},
        "compression": compression,
        "includes_duckdb": include_duckdb,
        "offline_use": "unpack at the project root (or mount data/, results/, configs/, SOURCE_REGISTRY.yaml into "
                       "/app) and run with MEASURE_IT_OFFLINE=1; docs/DEPLOYMENT.md",
        "included": {"data/processed": "*.parquet, *.geoparquet, *.meta.json" + (", measure_it_public.duckdb"
                                                                                 if include_duckdb else ""),
                     "data/raw": "per-source " + ", ".join(RAW_SIDECARS) + "; whole: " + ", ".join(RAW_WHOLE_DIRS)
                                 + "; " + ", ".join("/".join(x) for x in RAW_SUBDIRS),
                     "data/_http_cache": "entries for hosts " + ", ".join(HTTP_CACHE_HOSTS),
                     "results": "*.md, *.json, tables/**", "configs": "all", "SOURCE_REGISTRY.yaml": "yes"},
        "excluded": EXCLUDED,
        "privacy_excluded_records": sorted(_privacy_excluded()),
        "privacy_note": "No personal health information: public, de-identified or aggregate sources only; the "
                        "privacy-excluded records were never downloaded and no bundled path names them.",
        "files_changed_during_build": changed,
        "totals": {"files": len(entries), "bytes": sum(e["bytes"] for e in entries), "groups": groups},
        "processed_table_rows": tables,
        "sources": source_versions(),
        "files": entries,
    }


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ unpack / verify

def unpack_bundle(archive: Path, dest: Path, echo: Callable[[str], None] = print) -> Path:
    """Unpack an archive into dest (created); returns dest. Uses the external zstd for .tar.zst."""
    archive, dest = Path(archive), Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    if archive.name.endswith(".zst"):
        proc = subprocess.Popen([shutil.which("zstd") or "zstd", "-dc", "-q", str(archive)], stdout=subprocess.PIPE)
        with tarfile.open(fileobj=proc.stdout, mode="r|") as tar:
            tar.extractall(dest, filter="data")
        if proc.wait() != 0:
            raise RuntimeError(f"zstd -d exited with {proc.returncode}")
    else:
        with tarfile.open(archive, mode="r:*") as tar:
            tar.extractall(dest, filter="data")
    echo(f"bundle: unpacked {archive.name} into {dest}")
    return dest


def verify_bundle(root: Path, *, check_hashes: bool = True) -> dict:
    """Compare an unpacked bundle with its BUNDLE_MANIFEST.json: missing files, size and SHA-256 mismatches."""
    root = Path(root)
    man = json.loads((root / MANIFEST_NAME).read_text())
    missing, size_bad, hash_bad = [], [], []
    for e in man["files"]:
        p = root / e["path"]
        if not p.is_file():
            missing.append(e["path"])
            continue
        if p.stat().st_size != e["bytes"]:
            size_bad.append(e["path"])
            continue
        if check_hashes and _sha256_file(p) != e["sha256"]:
            hash_bad.append(e["path"])
    ok = not (missing or size_bad or hash_bad)
    return {"ok": ok, "n_files": len(man["files"]), "missing": missing, "size_mismatch": size_bad,
            "sha256_mismatch": hash_bad, "git": man.get("git"), "created_at": man.get("created_at")}
