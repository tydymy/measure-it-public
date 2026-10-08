"""Bulk file downloads with a per-source MANIFEST.json (url, bytes, sha256, retrieved_at).

Files that already exist and are recorded in the manifest are not downloaded
again unless refresh=True. HTML error pages served as HTTP 200 are rejected.
Offline mode (MEASURE_IT_OFFLINE=1): a file must already exist and be recorded in
MANIFEST.json; anything else raises http.OfflineCacheMiss instead of downloading.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import time
from contextlib import contextmanager
from pathlib import Path

import requests

from .config import RAW, raw_dir, utc_now_iso
from .http import OFFLINE_ENV, HTMLInsteadOfData, OfflineCacheMiss, offline, user_agent_for


def manifest_path(source_id: str) -> Path:
    return raw_dir(source_id) / "MANIFEST.json"


@contextmanager
def _locked(path: Path):
    lock = path.with_suffix(".lock")
    with open(lock, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def load_manifest(source_id: str) -> dict:
    p = manifest_path(source_id)
    return json.loads(p.read_text()) if p.exists() else {"source_id": source_id, "files": {}}


def _record(source_id: str, rel: str, entry: dict) -> None:
    p = manifest_path(source_id)
    with _locked(p):
        m = load_manifest(source_id)
        m["files"][rel] = entry
        p.write_text(json.dumps(m, indent=2, sort_keys=True))


def record_file(source_id: str, rel_path: str, *, url: str, note: str = "", **extra) -> dict:
    """Record a file that was assembled locally (e.g. paged API responses) in MANIFEST.json."""
    path = raw_dir(source_id) / rel_path
    entry = {
        "url": url,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "retrieved_at": utc_now_iso(),
        "note": note,
        **extra,
    }
    _record(source_id, rel_path, entry)
    return entry


manifest_lock = _locked  # public name for modules that update MANIFEST.json themselves


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(chunk):
            h.update(b)
    return h.hexdigest()


def download_file(
    url: str,
    source_id: str,
    filename: str | None = None,
    *,
    refresh: bool = False,
    min_bytes: int = 1,
    allow_html: bool = False,
    timeout: float = 600,
    headers: dict | None = None,
    max_retries: int = 1,
) -> Path:
    """Download url into data/raw/<source_id>/<filename> and record it in the manifest.

    filename may contain subdirectories. max_retries > 1 retries 429/5xx and
    connection errors with exponential backoff.
    """
    dest_dir = raw_dir(source_id)
    filename = filename or url.rstrip("/").split("/")[-1].split("?")[0]
    dest = dest_dir / filename
    rel = str(dest.relative_to(RAW / source_id))
    manifest = load_manifest(source_id)
    if offline():
        # Offline reproduction: the raw file recorded in MANIFEST.json is the input; never re-download it.
        if dest.exists() and rel in manifest["files"]:
            return dest
        raise OfflineCacheMiss(f"{OFFLINE_ENV}=1: raw file data/raw/{source_id}/{rel} is missing or not in "
                               f"MANIFEST.json; downloading {url} needs the network")
    if dest.exists() and not refresh and rel in manifest["files"]:
        return dest

    hdrs = {"User-Agent": user_agent_for(url)}
    if headers:
        hdrs.update(headers)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(max(1, max_retries)):
        try:
            r = requests.get(url, stream=True, timeout=timeout, headers=hdrs)
        except requests.RequestException:
            if attempt + 1 >= max_retries:
                raise
            time.sleep(min(60, 2 ** (attempt + 1)))
            continue
        if r.status_code in (429, 500, 502, 503, 504) and attempt + 1 < max_retries:
            r.close()
            time.sleep(min(60, 2 ** (attempt + 1)))
            continue
        break
    with r:
        r.raise_for_status()
        first = True
        with open(tmp, "wb") as fh:
            for block in r.iter_content(chunk_size=1 << 20):
                if first:
                    head = block[:512].lstrip().lower()
                    if not allow_html and (head.startswith(b"<!doctype html") or head.startswith(b"<html")):
                        fh.close()
                        tmp.unlink(missing_ok=True)
                        raise HTMLInsteadOfData(f"{url} returned an HTML page, not a data file")
                    first = False
                fh.write(block)
        resp_headers = dict(r.headers)
        final_url = r.url.split("?", 1)[0]  # drop signed query strings (provider keys/tokens) before recording
    size = tmp.stat().st_size
    if size < min_bytes:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{url}: only {size} bytes (< {min_bytes})")
    os.replace(tmp, dest)
    _record(source_id, rel, {
        "url": url.split("?", 1)[0] if ("Signature=" in url or "X-Amz-" in url or "sig=" in url) else url,
        "final_url": final_url,
        "bytes": size,
        "sha256": sha256_file(dest),
        "retrieved_at": utc_now_iso(),
        "last_modified": resp_headers.get("Last-Modified"),
        "etag": resp_headers.get("ETag"),
        "content_type": resp_headers.get("Content-Type"),
    })
    return dest


def retrieved_at(source_id: str, filename: str) -> str | None:
    return load_manifest(source_id)["files"].get(filename, {}).get("retrieved_at")
