"""Smoke test of a running deployment (docker compose or local servers): health, every tool route, MCP over HTTP.

    uv run python deploy/smoke_test.py                                   # defaults: API :8000, MCP :8765, dashboard :8501
    uv run python deploy/smoke_test.py --api http://127.0.0.1:18000 --mcp http://127.0.0.1:18765 \
        --dashboard http://127.0.0.1:18501

Every tool is called over GET with the demo arguments of tests/test_api.py (must answer `ok`) and with the input that
must come back `UNKNOWN / NOT AVAILABLE`; the MCP server must list 18 tools and answer one call. Exit 1 on any failure.
"""
from __future__ import annotations

import argparse
import asyncio
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from test_api import CASES  # noqa: E402  (demo and UNKNOWN arguments per tool)

UNKNOWN = "UNKNOWN / NOT AVAILABLE"


def get(url: str, timeout: float = 120) -> tuple[int, bytes, dict]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return r.status, r.read(), dict(r.headers)
    except urllib.error.HTTPError as e:
        return e.code, e.read(), dict(e.headers)


def query(args: dict) -> str:
    return urllib.parse.urlencode([(k, x) for k, v in args.items() if v is not None
                                   for x in (v if isinstance(v, list) else [v])])


async def mcp_check(url: str) -> dict:
    from fastmcp import Client
    async with Client(url) as c:
        tools = await c.list_tools()
        res = await c.call_tool("normalize_condition", {"condition_or_code": "G93.32"})
        env = json.loads(res.content[0].text)
        return {"n_tools": len(tools), "call_status": env.get("status"),
                "call_condition_id": (env.get("data") or {}).get("condition_id")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--api", default="http://127.0.0.1:8000")
    ap.add_argument("--mcp", default="http://127.0.0.1:8765")
    ap.add_argument("--dashboard", default="http://127.0.0.1:8501")
    a = ap.parse_args()
    fails, rows = [], []
    for name, url in (("api", f"{a.api}/health"), ("mcp", f"{a.mcp}/health"),
                      ("dashboard", f"{a.dashboard}/_stcore/health")):
        code, body, _ = get(url, 15)
        ok = code == 200 and (body == b"ok" or json.loads(body).get("status") == "ok")
        print(f"health {name:10s} {code} {body[:120].decode(errors='replace')}")
        if not ok:
            fails.append(f"health {name}")
    for tool, (demo, unknown) in CASES.items():
        for kind, args, want in (("demo", demo, "ok"), ("unknown", unknown, UNKNOWN)):
            t0 = time.time()
            code, body, hdr = get(f"{a.api}/tools/{tool}?{query(args)}")
            status = json.loads(body).get("status") if body[:1] == b"{" else None
            rows.append((tool, kind, code, status, round(time.time() - t0, 2)))
            if code != 200 or status != want:
                fails.append(f"{tool} {kind}: HTTP {code}, status {status!r} (want {want!r})")
    for r in rows:
        print(f"{r[0]:36s} {r[1]:8s} HTTP {r[2]}  {r[3]!s:26s} {r[4]:6.2f} s")
    m = asyncio.run(mcp_check(f"{a.mcp}/mcp"))
    print(f"mcp over HTTP: {m}")
    if m["n_tools"] != len(CASES) or m["call_status"] != "ok":
        fails.append(f"mcp: {m}")
    print(f"\n{len(rows)} tool calls, {len(fails)} failures")
    for f in fails:
        print("FAIL", f)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
