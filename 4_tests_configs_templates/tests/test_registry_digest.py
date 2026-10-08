"""registry.write_provenance_digest: per-download sidecars under data/raw -> two digest files."""
from __future__ import annotations

import json

from measure_it.registry import write_provenance_digest


def test_digest_gathers_manifests_and_audits_and_skips_uploads(tmp_path):
    raw = tmp_path / "data" / "raw"
    (raw / "src_a").mkdir(parents=True)
    (raw / "src_a" / "MANIFEST.json").write_text(json.dumps({"sha256": "abc", "url": "https://example.org/a.csv"}))
    (raw / "src_a" / "DATA_AUDIT.md").write_text("# DATA AUDIT — A\n\n| Field | Value |\n|---|---|\n")
    (raw / "src_b" / "nested").mkdir(parents=True)
    (raw / "src_b" / "nested" / "MANIFEST.json").write_text(json.dumps({"sha256": "def"}))
    (raw / "src_b" / "queries.json").write_text("not json")
    (raw / "byod_private").mkdir()
    (raw / "byod_private" / "MANIFEST.json").write_text(json.dumps({"sha256": "secret"}))
    (raw / "byod_private" / "DATA_AUDIT.md").write_text("# private\n")

    mpath, apath = write_provenance_digest(raw)

    files = json.loads(mpath.read_text())["files"]
    assert set(files) == {"data/raw/src_a/MANIFEST.json", "data/raw/src_b/nested/MANIFEST.json",
                          "data/raw/src_b/queries.json"}
    assert files["data/raw/src_a/MANIFEST.json"]["sha256"] == "abc"
    assert files["data/raw/src_b/queries.json"] == {"unreadable": True}
    md = apath.read_text()
    assert "## Source: src_a" in md and "## DATA AUDIT — A" in md   # audit heading demoted one level
    assert "byod_private" not in md and "secret" not in mpath.read_text()
