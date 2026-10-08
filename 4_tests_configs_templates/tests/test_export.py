"""Tests for the exports: the read-only data bundle (measure_it.export.bundle) and the static reviewer snapshot
(measure_it.export.static_site).

The bundle tests build a bundle from a small synthetic project tree (paths monkeypatched), so they need no data. The
snapshot tests are marked `data`: they build the page from data/processed and results/ into a temporary file and check
its size, product language, external hosts and that its numbers equal the processed tables.
"""
from __future__ import annotations

import hashlib
import json
import re
import tarfile
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import pytest

from measure_it.config import PROCESSED, RESULTS, UNKNOWN
from measure_it.export import bundle as B

HAS_DATA = (PROCESSED / "deployment_opportunities.parquet").exists() and (RESULTS / "FINAL_REPORT.md").exists()


# ------------------------------------------------------------------------------------------------ bundle

@pytest.fixture()
def fake_project(tmp_path, monkeypatch):
    root = tmp_path / "proj"
    processed, raw, cache = root / "data" / "processed", root / "data" / "raw", root / "data" / "_http_cache"
    results, configs = root / "results", root / "configs"
    for d in (processed, raw / "cdc_places", raw / "data_gov_catalog" / "search_responses",
              raw / "clinicaltrials_gov" / "records", raw / "lab_candidates", cache / "ab", results / "tables",
              results / "figures", results / "pipeline_runs", configs, root / "data" / "interim" / "x"):
        d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"a": range(7), "b": list("abcdefg")}).to_parquet(processed / "t1.parquet")
    (processed / "t1.meta.json").write_text("{}")
    (raw / "cdc_places" / "MANIFEST.json").write_text('{"files": {}}')
    (raw / "cdc_places" / "DATA_AUDIT.md").write_text("# audit\n")
    (raw / "cdc_places" / "big_raw_file.csv").write_text("x\n" * 100)                 # raw data: must stay out
    (raw / "data_gov_catalog" / "search_responses" / "q.json").write_text("{}")
    (raw / "clinicaltrials_gov" / "records" / "batch_0001.json").write_text("{}")
    (raw / "clinicaltrials_gov" / "big.zip").write_bytes(b"0" * 50)                    # raw data: must stay out
    for key, url in (("ab01", "https://catalog.data.gov/api/3/action/package_search"),
                     ("ab02", "https://api.fda.gov/device/510k.json")):
        (cache / "ab" / f"{key}.meta.json").write_text(json.dumps({"url": url}))
        (cache / "ab" / f"{key}.body").write_text("{}")
    (results / "FINAL_REPORT.md").write_text("# report\n")
    (results / "example.json").write_text("{}")
    (results / "tables" / "t.csv").write_text('a,b\n1,"multi\nline"\n2,x\n3,y\n')      # 3 records
    (results / "figures" / "f.png").write_bytes(b"png")
    (results / "pipeline_runs" / "run.jsonl").write_text("{}")
    (configs / "scoring.yaml").write_text("x: 1\n")
    (root / "SOURCE_REGISTRY.yaml").write_text(
        "sources:\n- source_id: cdc_places\n  source_version: '2024'\n  retrieved_at: '2026-09-01'\n")
    for name, val in (("PROJECT_ROOT", root), ("PROCESSED", processed), ("RAW", raw), ("HTTP_CACHE", cache),
                      ("RESULTS", results), ("CONFIGS", configs), ("SOURCE_REGISTRY_PATH", root / "SOURCE_REGISTRY.yaml"),
                      ("DUCKDB_PATH", processed / "measure_it_public.duckdb"), ("DIST", root / "dist")):
        monkeypatch.setattr(B, name, val)
    return root


def test_bundle_selects_only_runtime_files(fake_project):
    rel = {f.relative_to(fake_project).as_posix() for f in B.collect_files()}
    assert "data/processed/t1.parquet" in rel and "data/processed/t1.meta.json" in rel
    assert "data/raw/cdc_places/MANIFEST.json" in rel and "data/raw/cdc_places/DATA_AUDIT.md" in rel
    assert "data/raw/data_gov_catalog/search_responses/q.json" in rel
    assert "data/raw/clinicaltrials_gov/records/batch_0001.json" in rel
    assert "data/_http_cache/ab/ab01.body" in rel and "data/_http_cache/ab/ab02.body" not in rel   # data.gov only
    assert not any(p.endswith(("big_raw_file.csv", "big.zip")) for p in rel)                     # no raw data
    assert not any(p.startswith(("data/interim", "results/figures", "results/pipeline_runs")) for p in rel)
    assert {"results/FINAL_REPORT.md", "results/example.json", "results/tables/t.csv", "configs/scoring.yaml",
            "SOURCE_REGISTRY.yaml"} <= rel


def test_bundle_refuses_privacy_excluded_records(fake_project):
    bad = next(iter(B._privacy_excluded()))
    (fake_project / "data" / "raw" / "lab_candidates" / "MANIFEST.json").unlink(missing_ok=True)
    d = fake_project / "data" / "raw" / bad
    d.mkdir(parents=True)
    (d / "DATA_AUDIT.md").write_text("x")
    with pytest.raises(RuntimeError, match="privacy-excluded"):
        B.collect_files()


def test_bundle_manifest_is_correct_and_round_trips(fake_project, tmp_path):
    res = B.build_bundle(fake_project / "dist", compression="gz", date="2026-01-02", echo=lambda *_: None)
    archive = Path(res["archive"])
    assert archive.name == "measure-it-data-2026-01-02.tar.gz" and archive.stat().st_size == res["archive_bytes"]
    assert (archive.parent / f"{archive.name}.sha256").read_text().split()[0] == \
        hashlib.sha256(archive.read_bytes()).hexdigest()
    man = json.loads(Path(res["manifest"]).read_text())
    with tarfile.open(archive) as tar:
        names = tar.getnames()
        assert names[-1] == B.MANIFEST_NAME
        assert json.loads(tar.extractfile(B.MANIFEST_NAME).read()) == man
    files = {e["path"]: e for e in man["files"]}
    assert set(files) == {f.relative_to(fake_project).as_posix() for f in B.collect_files()}
    for p, e in files.items():   # sha256 and size of every entry equal the source bytes
        data = (fake_project / p).read_bytes()
        assert e["bytes"] == len(data) and e["sha256"] == hashlib.sha256(data).hexdigest(), p
    assert files["data/processed/t1.parquet"]["rows"] == 7
    assert files["results/tables/t.csv"]["rows"] == 3            # quoted newline is one record
    assert man["processed_table_rows"] == {"t1": 7}
    assert man["totals"]["files"] == len(files) and man["totals"]["bytes"] == sum(e["bytes"] for e in files.values())
    assert [s["source_id"] for s in man["sources"]] == ["cdc_places"] and man["sources"][0]["source_version"] == "2024"
    assert "commit" in man["git"] and man["bundle_format"] == B.BUNDLE_FORMAT
    assert set(man["privacy_excluded_records"]) == B._privacy_excluded()
    # unpack + verify, then tamper with one file
    dest = B.unpack_bundle(archive, tmp_path / "unpacked", echo=lambda *_: None)
    rep = B.verify_bundle(dest)
    assert rep["ok"] and rep["n_files"] == len(files)
    (dest / "results" / "tables" / "t.csv").write_text('a,b\n9,"multi\nline"\n2,x\n3,y\n')
    rep = B.verify_bundle(dest)
    assert not rep["ok"] and rep["sha256_mismatch"] == ["results/tables/t.csv"]


def test_bundle_without_processed_tables_fails_clearly(fake_project):
    (fake_project / "data" / "processed" / "t1.parquet").unlink()
    with pytest.raises(RuntimeError, match="run the pipeline first"):
        B.build_bundle(fake_project / "dist", compression="gz", echo=lambda *_: None)


@pytest.mark.data
@pytest.mark.skipif(not HAS_DATA, reason="needs data/processed")
def test_real_bundle_selection_has_runtime_tables_and_no_raw_data():
    rel = [f.relative_to(B.PROJECT_ROOT).as_posix() for f in B.collect_files()]
    assert "data/processed/deployment_opportunities.parquet" in rel
    assert "data/processed/county_boundaries_2024.geoparquet" in rel
    for p in rel:
        if p.startswith("data/raw/"):
            parts = p.split("/")
            assert parts[-1] in B.RAW_SIDECARS or parts[2] in B.RAW_WHOLE_DIRS or \
                "/".join(parts[2:4]) in {"/".join(x) for x in B.RAW_SUBDIRS}, p


# ------------------------------------------------------------------------------------------------ static snapshot

class _Page(HTMLParser):
    """Visible text, resource URLs (script src, link href, img/iframe/source src) and the <title>."""

    def __init__(self) -> None:
        super().__init__()
        self.text, self.resources, self.title, self._skip, self._in_title = [], [], "", 0, False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in ("script", "style"):
            self._skip += 1
        if tag == "title":
            self._in_title = True
        for key in ("src", "srcset", "data", "poster"):
            if a.get(key):
                self.resources.append(a[key])
        if tag == "link" and a.get("href"):
            self.resources.append(a["href"])

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self._skip -= 1
        if tag == "title":
            self._in_title = False

    def handle_data(self, data):
        if self._in_title:
            self.title += data
        elif not self._skip:
            self.text.append(data)


@pytest.fixture(scope="module")
def snapshot(tmp_path_factory):
    if not HAS_DATA:
        pytest.skip("needs data/processed and results/FINAL_REPORT.md")
    import os
    os.environ["MEASURE_IT_OFFLINE"] = "1"
    from measure_it.export import static_site as S
    out = tmp_path_factory.mktemp("static") / "index.html"
    res = S.build_static_site(out)
    html = out.read_text()
    page = _Page()
    page.feed(html)
    blob = re.search(r'<script type="application/json" id="snapshot-data">(.*?)</script>', html, re.S).group(1)
    return {"res": res, "html": html, "page": page, "data": json.loads(blob), "S": S}


@pytest.mark.data
def test_static_page_size_and_title(snapshot):
    assert snapshot["res"]["bytes"] == len(snapshot["html"].encode())
    assert snapshot["res"]["bytes"] <= 12 * 1024 * 1024
    assert 2 <= len(snapshot["page"].title.split()) <= 4


@pytest.mark.data
def test_static_page_requests_only_allowed_cdns(snapshot):
    html, S = snapshot["html"], snapshot["S"]
    for url in snapshot["page"].resources:
        host = urlparse(url).netloc
        assert host in S.ALLOWED_SCRIPT_HOSTS, url
    # no other way of loading anything: no CSS url()/@import, no fetch/XHR/websocket, no iframes or external anchors
    assert not re.search(r"url\(\s*['\"]?https?:", html) and "@import" not in html
    assert not re.search(r"\b(fetch|XMLHttpRequest|WebSocket|EventSource)\s*\(", html)
    assert "<iframe" not in html
    hosts = {urlparse(u).netloc for u in re.findall(r"https?://[^\s\"'<>)]+", html)}
    assert hosts <= set(S.ALLOWED_SCRIPT_HOSTS), hosts


@pytest.mark.data
def test_static_page_has_no_forbidden_language(snapshot):
    from measure_it.validate import forbidden_hits
    text = "\n".join(t.strip() for t in snapshot["page"].text if t.strip())
    assert forbidden_hits(text) == []


@pytest.mark.data
def test_static_page_has_every_view(snapshot):
    html, text = snapshot["html"], " ".join(snapshot["page"].text)
    for sid in ("overview", "query", "map", "cards", "performance", "works", "validation", "limits"):
        assert f'<section id="{sid}">' in html
    assert text.count("candidate sites") >= 10 and html.count('<details class="card"') == 10
    assert "four data layers" in text.lower() and "evidence_weighted" in text
    report = (RESULTS / "FINAL_REPORT.md").read_text()
    for h in (r"^## 13\. (.*)$", r"^### 14\.2 (.*)$", r"^### 14\.3 (.*)$", r"^### 14\.4 (.*)$", r"^## 19\. (.*)$"):
        assert re.search(h, report, re.M).group(1).strip() in text


@pytest.mark.data
def test_static_page_numbers_match_tables(snapshot):
    S, data, html = snapshot["S"], snapshot["data"], snapshot["html"]
    rows = S.opportunity_rows()
    top = rows[rows["rank_equal"].notna()].sort_values("rank_equal").head(10)
    assert [r["geo_id"] for r in data["top10_equal"]] == list(top["geo_id"])
    for r, t in zip(data["top10_equal"], top.to_dict("records")):
        for k in ("rank_equal", "composite_equal", "rank_mc_p05", "rank_mc_p95", "p_top10_mc", "burden_pct",
                  "vulnerability_pct", "diagnostic_desert_pct", "clinic_capacity_pct", "research_readiness_pct",
                  "rank_evidence_weighted"):
            assert r[k] == pytest.approx(float(t[k]), abs=1e-9), (t["geo_id"], k)
        row = re.search(rf'<tr data-geo="{t["geo_id"]}" data-rank="(\d+)" data-composite="([\d.]+)">(.*?)</tr>',
                        html, re.S)
        assert int(row.group(1)) == int(t["rank_equal"]) and row.group(2) == f"{t['composite_equal']:.3f}"
        assert f"{int(round(t['rank_mc_p05']))}-{int(round(t['rank_mc_p95']))}" in row.group(3)
    ew = rows[rows["rank_evidence_weighted"].notna()].sort_values("rank_evidence_weighted").head(10)
    assert [r["geo_id"] for r in data["top10_evidence_weighted"]] == list(ew["geo_id"])
    assert data["n_regions_ranked"] == int(rows["rank_equal"].notna().sum())


@pytest.mark.data
def test_static_page_performance_table_matches_and_keeps_unknown(snapshot):
    S, html = snapshot["S"], snapshot["html"]
    perf = S.performance_table()
    rows = re.findall(r'<tr data-cond="([^"]+)" data-meas="([^"]+)" class="[^"]*">(.*?)</tr>', html, re.S)
    assert len(rows) == len(perf)
    got = {(c, m): body for c, m, body in rows}
    for r in perf.to_dict("records"):
        body = got[(r["condition_id"], r["measurement_id"])]
        if r["performance_status"] == UNKNOWN:
            assert UNKNOWN in body and "UNKNOWN" in body.split("</td>")[4]    # AUROC cell says UNKNOWN
        elif pd.notna(r["auroc"]):
            assert f"{r['auroc']:.3f} ({r['auroc_ci_low']:.3f}-{r['auroc_ci_high']:.3f})" in body
