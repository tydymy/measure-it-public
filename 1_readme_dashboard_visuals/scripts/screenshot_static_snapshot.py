"""Open dist/static_snapshot/index.html headlessly (Playwright Chromium) at desktop and phone widths and capture
full-page screenshots to results/figures/static_snapshot_{desktop,mobile}.png.

    uv run --with playwright python scripts/screenshot_static_snapshot.py [--scheme light|dark]

Reports console errors, failed requests, requests to hosts other than cdn.jsdelivr.net, horizontal page overflow and
whether the county map drew."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "dist" / "static_snapshot" / "index.html"
OUT = ROOT / "results" / "figures"
VIEWS = {"desktop": {"width": 1440, "height": 1000}, "mobile": {"width": 390, "height": 844}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme", default="light", choices=["light", "dark"])
    ap.add_argument("--page", default=str(PAGE))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    from playwright.sync_api import sync_playwright
    bad = 0
    with sync_playwright() as p:
        b = p.chromium.launch()
        for name, vp in VIEWS.items():
            ctx = b.new_context(viewport=vp, device_scale_factor=1,
                                color_scheme=a.scheme, is_mobile=name == "mobile", has_touch=name == "mobile")
            pg = ctx.new_page()
            errors, hosts = [], set()
            pg.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
            pg.on("pageerror", lambda e: errors.append(str(e)))
            pg.on("request", lambda r: hosts.add(urlparse(r.url).netloc) if r.url.startswith("http") else None)
            pg.goto(Path(a.page).resolve().as_uri(), wait_until="networkidle")
            pg.wait_for_timeout(800)
            n_paths = pg.locator("#county-map svg path.county").count()
            overflow = pg.evaluate("document.documentElement.scrollWidth - document.documentElement.clientWidth")
            f = Path(a.out) / f"static_snapshot_{name}{a.suffix}.png"
            f.parent.mkdir(parents=True, exist_ok=True)
            pg.screenshot(path=str(f), full_page=True)
            other = sorted(h for h in hosts if h not in ("cdn.jsdelivr.net",))
            print(f"{name}: {f} | county paths {n_paths} | horizontal overflow {overflow}px | "
                  f"hosts {sorted(hosts)} | console errors {errors[:3]}")
            bad += bool(errors) + bool(other) + (overflow > 0) + (n_paths == 0)
            ctx.close()
        b.close()
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
