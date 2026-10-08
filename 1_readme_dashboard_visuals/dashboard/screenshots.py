"""Capture a full-page screenshot of the dashboard home page and each page to results/figures/dashboard/<page>.png.

    uv run --with playwright playwright install chromium          # once (downloads a headless Chromium)
    uv run --with playwright python dashboard/screenshots.py      # starts streamlit, captures, stops it

Options: --port 8599, --only national_opportunity_map (repeatable), --width 1500, --out results/figures/dashboard.
The server runs with MEASURE_IT_OFFLINE=1 (answers from data/processed and the HTTP cache only). A page counts as
rendered when Streamlit's "running" indicator has gone, every Plotly chart has drawn and the page height is stable.
"""
from __future__ import annotations

import argparse
import os
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "dashboard" / "app.py"
OUT = ROOT / "results" / "figures" / "dashboard"
PAGES = [  # (file name, url path, expected plotly charts at least)
    ("home", "", 0),
    ("national_opportunity_map", "National_Opportunity_Map", 1),
    ("condition_explorer", "Condition_Explorer", 1),
    ("measurement_explorer", "Measurement_Explorer", 1),
    ("geography_explorer", "Geography_Explorer", 1),
    ("deployment_recommendation", "Deployment_Recommendation", 2),
]


def _free_port(preferred: int) -> int:
    with socket.socket() as s:
        try:
            s.bind(("127.0.0.1", preferred))
            return preferred
        except OSError:
            s.bind(("127.0.0.1", 0))
            return s.getsockname()[1]


def _wait_health(port: int, timeout: float = 90) -> None:
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            with urllib.request.urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=2) as r:
                if r.status == 200:
                    return
        except OSError:
            time.sleep(0.5)
    raise RuntimeError("streamlit did not become healthy")


def _wait_rendered(page, n_plotly: int, timeout_s: float = 240) -> None:
    t0 = time.time()
    page.wait_for_selector('[data-testid="stAppViewContainer"]', timeout=60_000)
    last_h, stable = -1, 0
    while time.time() - t0 < timeout_s:
        running = page.locator('[data-testid="stStatusWidget"]').count() > 0
        skeleton = page.locator('[data-testid="stSkeleton"]').count() > 0
        plots = page.locator(".js-plotly-plot .main-svg").count()
        h = page.evaluate("document.body.scrollHeight")
        if not running and not skeleton and plots >= n_plotly:
            stable = stable + 1 if h == last_h else 0
            if stable >= 3:
                return
        last_h = h
        time.sleep(1.0)
    print(f"  warning: page not settled after {timeout_s}s (captured anyway)", flush=True)


def _capture(a, out: Path, base: str) -> int:
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        browser = p.chromium.launch()
        ctx = browser.new_context(viewport={"width": a.width, "height": 1000}, device_scale_factor=1,
                                  color_scheme="light")
        page = ctx.new_page()
        for name, path, n_plotly in PAGES:
            if a.only and name not in a.only:
                continue
            t0 = time.time()
            page.goto(f"{base}/{path}", wait_until="domcontentloaded")
            _wait_rendered(page, n_plotly)
            time.sleep(1.5)
            exc = page.locator('[data-testid="stException"]').count()
            # Streamlit scrolls inside its main container, so grow the viewport to the content height
            h = page.evaluate("""() => Math.max(...['[data-testid="stMain"]', '[data-testid="stAppViewContainer"]',
                'section.main'].map(s => document.querySelector(s)).filter(Boolean).map(e => e.scrollHeight),
                document.body.scrollHeight)""")
            page.set_viewport_size({"width": a.width, "height": int(min(h, 16000)) + 40})
            time.sleep(2.0)
            f = out / f"{name}.png"
            page.screenshot(path=str(f), full_page=True)
            page.set_viewport_size({"width": a.width, "height": 1000})
            shown = f.relative_to(ROOT) if f.resolve().is_relative_to(ROOT) else f
            print(f"{name}: {shown} ({time.time() - t0:.1f}s; exceptions on page: {exc})", flush=True)
        browser.close()
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8599)
    ap.add_argument("--only", action="append", default=None)
    ap.add_argument("--width", type=int, default=1500)
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--url", default=None, help="Capture a dashboard that is already running (e.g. the container, "
                                                 "http://127.0.0.1:8501) instead of starting one.")
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.url:
        return _capture(a, out, a.url.rstrip("/"))
    port = _free_port(a.port)
    env = {**os.environ, "MEASURE_IT_OFFLINE": "1"}
    cmd = [sys.executable, "-m", "streamlit", "run", str(APP), "--server.headless", "true", "--server.port",
           str(port), "--browser.gatherUsageStats", "false", "--theme.base", "light", "--server.runOnSave", "false",
           "--client.toolbarMode", "viewer"]
    log = tempfile.NamedTemporaryFile("w", prefix="streamlit_screenshot_", suffix=".log", delete=False)
    print(f"streamlit log: {log.name}", flush=True)
    proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    try:
        _wait_health(port)
        _capture(a, out, f"http://127.0.0.1:{port}")
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            proc.kill()
        log.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
