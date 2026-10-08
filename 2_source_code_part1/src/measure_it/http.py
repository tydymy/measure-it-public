"""Cached, polite HTTP for federal and biomedical APIs.

Every API response is cached on disk under data/_http_cache keyed by
method + URL + params + body, so reruns do not hit federal services again.
Pass refresh=True to bypass the cache for one call.

Guard: several federal hosts (CDC NHANES in particular) return an HTML error
page with HTTP 200 for moved files. `expect_json=True` and `reject_html=True`
turn that into an exception instead of silently caching garbage.

Offline mode: with the environment variable MEASURE_IT_OFFLINE=1 (set by
`measure-it pipeline --offline`), `request` answers only from the disk cache and
raises `OfflineCacheMiss` on a miss instead of touching the network; a
`refresh=True` call is served from the cache when an entry exists (it cannot be
refreshed offline) and raises otherwise. `download.download_file` and the few
modules that stream bulk files with `requests` directly call `ensure_online()`
first. `install_offline_guard()` additionally blocks every outbound TCP/UDP
connection at the socket level, so an offline run provably used cached/raw data
only (measure_it/__init__.py installs it whenever the variable is set).
"""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse

import requests

from .config import HTTP_CACHE, utc_now_iso

USER_AGENT = "measure-it-public/0.1 (public-data research alpha; contact: repository maintainer)"

# Minimum seconds between requests to the same host (per process).
HOST_MIN_INTERVAL = {
    "eutils.ncbi.nlm.nih.gov": 0.4,      # NCBI: 3 req/s without an API key
    "api.reporter.nih.gov": 1.0,         # RePORTER asks for <= 1 req/s
    "api.fda.gov": 0.3,                  # openFDA: 240 req/min without a key
    "clinicaltrials.gov": 0.2,
    "www.ebi.ac.uk": 0.2,
    "api.platform.opentargets.org": 0.2,
    "data.cdc.gov": 0.2,
    "api.census.gov": 0.2,
    "catalog.data.gov": 10.0,            # robots.txt Crawl-Delay: 10
}
DEFAULT_MIN_INTERVAL = 0.2

# www.cdc.gov (Akamai) answers HTTP 403 to non-default User-Agents; wwwn/ftp.cdc.gov do not.
HOST_USER_AGENT = {
    "www.cdc.gov": f"python-requests/{requests.__version__}",
}

# Request headers that change the response and therefore belong in the cache key.
# Only added to the key when present, so entries cached without them keep their keys.
VARY_HEADERS = {"range", "accept", "x-api-key", "authorization"}


def user_agent_for(url: str) -> str:
    return HOST_USER_AGENT.get(urlparse(url).netloc, USER_AGENT)

_last_call: dict[str, float] = {}
_lock = threading.Lock()


class HTMLInsteadOfData(RuntimeError):
    """A data endpoint answered with an HTML page (often an error page served as HTTP 200)."""


OFFLINE_ENV = "MEASURE_IT_OFFLINE"


class OfflineCacheMiss(RuntimeError):
    """Offline mode (MEASURE_IT_OFFLINE=1) needed the network: an HTTP cache miss or a missing raw file."""


class OfflineNetworkBlocked(OfflineCacheMiss):
    """A socket connection was attempted while the offline socket guard was installed."""


def offline() -> bool:
    """True when MEASURE_IT_OFFLINE is set to a true value (1/true/yes/on)."""
    return os.environ.get(OFFLINE_ENV, "").strip().lower() in {"1", "true", "yes", "on"}


def ensure_online(what: str) -> None:
    """Raise OfflineCacheMiss when offline mode forbids the network access described by `what`."""
    if offline():
        raise OfflineCacheMiss(f"{OFFLINE_ENV}=1: {what} would need the network (cache miss / missing raw file)")


_GUARD_INSTALLED = [False]
_LOOPBACK = {"localhost", "127.0.0.1", "::1", "0.0.0.0", ""}


def install_offline_guard() -> None:
    """Block outbound internet sockets in this process (loopback and AF_UNIX stay allowed).

    Idempotent. Every path that reaches the network (requests, urllib, raw sockets, DNS lookups) goes
    through socket.getaddrinfo / socket.connect, so a completed run with the guard installed is proof
    that no network I/O happened.
    """
    import socket

    if _GUARD_INSTALLED[0]:
        return
    real_connect, real_connect_ex = socket.socket.connect, socket.socket.connect_ex
    real_getaddrinfo = socket.getaddrinfo

    def _host(address) -> str | None:
        if isinstance(address, tuple) and address:
            return str(address[0])
        return None  # AF_UNIX path or other local address

    def _blocked(sock, address) -> bool:
        if getattr(sock, "family", None) not in (socket.AF_INET, socket.AF_INET6):
            return False
        h = _host(address)
        return h is not None and h not in _LOOPBACK and not h.startswith("127.")

    def connect(self, address):
        if _blocked(self, address):
            raise OfflineNetworkBlocked(f"{OFFLINE_ENV}=1: blocked socket connection to {address!r}")
        return real_connect(self, address)

    def connect_ex(self, address):
        if _blocked(self, address):
            raise OfflineNetworkBlocked(f"{OFFLINE_ENV}=1: blocked socket connection to {address!r}")
        return real_connect_ex(self, address)

    def getaddrinfo(host, *args, **kwargs):
        h = host.decode() if isinstance(host, bytes) else (str(host) if host is not None else "")
        if h not in _LOOPBACK and not h.startswith("127."):
            raise OfflineNetworkBlocked(f"{OFFLINE_ENV}=1: blocked DNS lookup of {h!r}")
        return real_getaddrinfo(host, *args, **kwargs)

    socket.socket.connect = connect
    socket.socket.connect_ex = connect_ex
    socket.getaddrinfo = getaddrinfo
    _GUARD_INSTALLED[0] = True


@dataclass
class CachedResponse:
    url: str
    status: int
    headers: dict
    content: bytes
    fetched_at: str
    from_cache: bool

    @property
    def text(self) -> str:
        return self.content.decode("utf-8", errors="replace")

    def json(self):
        return json.loads(self.content)

    def raise_for_status(self) -> None:
        if self.status >= 400:
            raise requests.HTTPError(f"HTTP {self.status} for {self.url}: {self.text[:300]}")


def _cache_key(method: str, url: str, params: dict | None, body, headers: dict | None = None) -> str:
    key = {"m": method.upper(), "u": url, "p": sorted((params or {}).items()), "b": body}
    vary = {k.lower(): v for k, v in (headers or {}).items() if k.lower() in VARY_HEADERS}
    if vary:
        key["h"] = sorted(vary.items())
    blob = json.dumps(key, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()


def _cache_paths(key: str) -> tuple[Path, Path]:
    d = HTTP_CACHE / key[:2]
    return d / f"{key}.meta.json", d / f"{key}.body"


def _looks_like_html(content: bytes, content_type: str) -> bool:
    head = content[:512].lstrip().lower()
    return "text/html" in content_type.lower() or head.startswith(b"<!doctype html") or head.startswith(b"<html")


def _throttle(url: str) -> None:
    host = urlparse(url).netloc
    interval = HOST_MIN_INTERVAL.get(host, DEFAULT_MIN_INTERVAL)
    with _lock:
        wait = _last_call.get(host, 0.0) + interval - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        _last_call[host] = time.monotonic()


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".tmp-")
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)
    os.replace(tmp, path)


def request(
    method: str,
    url: str,
    *,
    params: dict | None = None,
    json_body=None,
    data=None,
    headers: dict | None = None,
    refresh: bool = False,
    expect_json: bool = False,
    reject_html: bool = True,
    cache_errors: bool = False,
    timeout: float = 120,
    max_retries: int = 5,
) -> CachedResponse:
    body = json_body if json_body is not None else data
    key = _cache_key(method, url, params, body, headers)
    meta_path, body_path = _cache_paths(key)
    cached = meta_path.exists() and body_path.exists()
    if offline():
        if not cached:
            raise OfflineCacheMiss(f"{OFFLINE_ENV}=1: no HTTP cache entry for {method.upper()} {url} "
                                   f"params={params!r} (cache key {key[:12]})")
        refresh = False  # a refresh cannot happen offline; the cached response is the reproducible input
    if not refresh and cached:
        meta = json.loads(meta_path.read_text())
        return CachedResponse(
            url=meta["url"], status=meta["status"], headers=meta["headers"],
            content=body_path.read_bytes(), fetched_at=meta["fetched_at"], from_cache=True,
        )

    hdrs = {"User-Agent": user_agent_for(url)}
    if headers:
        hdrs.update(headers)
    last_exc: Exception | None = None
    for attempt in range(max_retries):
        _throttle(url)
        try:
            r = requests.request(method, url, params=params, json=json_body, data=data,
                                 headers=hdrs, timeout=timeout)
        except requests.RequestException as exc:
            last_exc = exc
            time.sleep(min(60, 2 ** attempt))
            continue
        if r.status_code in (429, 500, 502, 503, 504):
            retry_after = r.headers.get("Retry-After")
            time.sleep(float(retry_after) if retry_after and retry_after.isdigit() else min(60, 2 ** (attempt + 1)))
            last_exc = requests.HTTPError(f"HTTP {r.status_code} for {r.url}")
            continue
        break
    else:
        raise RuntimeError(f"giving up on {method} {url}: {last_exc}")

    ctype = r.headers.get("Content-Type", "")
    if (expect_json or reject_html) and r.status_code < 400 and _looks_like_html(r.content, ctype):
        if expect_json or reject_html:
            raise HTMLInsteadOfData(f"{url} returned HTML (HTTP {r.status_code}, {len(r.content)} bytes)")
    if expect_json and r.status_code < 400:
        json.loads(r.content)  # raises if not JSON

    resp = CachedResponse(url=r.url, status=r.status_code, headers=dict(r.headers),
                          content=r.content, fetched_at=utc_now_iso(), from_cache=False)
    if r.status_code < 400 or cache_errors:
        _atomic_write(body_path, r.content)
        _atomic_write(meta_path, json.dumps(
            {"url": resp.url, "status": resp.status, "headers": resp.headers,
             "fetched_at": resp.fetched_at, "method": method.upper(), "params": params, "body": body},
            default=str).encode())
    return resp


def get(url: str, **kw) -> CachedResponse:
    return request("GET", url, **kw)


def post(url: str, **kw) -> CachedResponse:
    return request("POST", url, **kw)


def get_json(url: str, **kw):
    r = request("GET", url, expect_json=True, **kw)
    r.raise_for_status()
    return r.json()


def post_json(url: str, **kw):
    r = request("POST", url, expect_json=True, **kw)
    r.raise_for_status()
    return r.json()
