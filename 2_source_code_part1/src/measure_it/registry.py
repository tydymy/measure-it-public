"""SOURCE_REGISTRY.yaml: the application's index of every public data source.

Each ingestion module writes data/raw/<source_id>/registry_entry.yaml. This
module merges those fragments over the base entries in SOURCE_REGISTRY.yaml,
so parallel ingestion never edits one shared file. It also gathers the per-download
sidecars under data/raw (MANIFEST.json, QUERY_LOG.json, queries.json, DATA_AUDIT.md),
which git does not keep, into data/SOURCE_MANIFESTS.json and data/SOURCE_AUDITS.md.
"""
from __future__ import annotations

import json
from pathlib import Path

import yaml

from .config import RAW, SOURCE_REGISTRY_PATH, utc_now_iso

REQUIRED_FIELDS = [
    "source_id", "name", "publisher", "landing_url", "access_urls", "license", "access_conditions",
    "retrieved_at", "source_version", "update_date", "data_layer", "unit_of_observation",
    "sample_size", "geographic_resolution", "person_level", "geographic", "omics", "wearable",
    "participant_linkage", "true_participant_linkage_across_modalities", "status",
    "processed_outputs", "audit", "ingestion_module", "limitations",
]


def load_source_registry(path: Path = SOURCE_REGISTRY_PATH) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def source_entry(source_id: str) -> dict | None:
    for s in load_source_registry().get("sources", []):
        if s["source_id"] == source_id:
            return s
    return None


def write_registry_entry(entry: dict) -> Path:
    missing = [f for f in REQUIRED_FIELDS if f not in entry]
    if missing:
        raise ValueError(f"registry entry {entry.get('source_id')}: missing {missing}")
    p = RAW / entry["source_id"] / "registry_entry.yaml"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(yaml.safe_dump(entry, sort_keys=False, allow_unicode=True))
    return p


def merge_registry_fragments(path: Path = SOURCE_REGISTRY_PATH) -> dict:
    reg = load_source_registry(path)
    by_id = {s["source_id"]: s for s in reg.get("sources", [])}
    for frag in sorted(RAW.glob("*/registry_entry.yaml")):
        entry = yaml.safe_load(frag.read_text())
        base = by_id.get(entry["source_id"], {})
        base.update(entry)
        by_id[entry["source_id"]] = base
    reg["sources"] = sorted(by_id.values(), key=lambda s: (s.get("data_layer", ""), s["source_id"]))
    reg["generated_at"] = utc_now_iso()
    path.write_text(yaml.safe_dump(reg, sort_keys=False, allow_unicode=True, width=110))
    write_provenance_digest()
    return reg


MANIFEST_NAMES = ("MANIFEST.json", "QUERY_LOG.json", "queries.json")
MANIFESTS_DIGEST = RAW.parent / "SOURCE_MANIFESTS.json"
AUDITS_DIGEST = RAW.parent / "SOURCE_AUDITS.md"


def write_provenance_digest(raw: Path = RAW) -> tuple[Path, Path]:
    """One file per kind instead of one per download: every MANIFEST.json / QUERY_LOG.json / queries.json under
    data/raw keyed by its path, and every DATA_AUDIT.md in one document (one section per source). Uploaded
    datasets (data/raw/byod_*, local-only) are left out."""
    root = raw.parent.parent
    public = lambda q: not q.relative_to(raw).parts[0].startswith("byod_")  # noqa: E731
    manifests = {}
    for q in sorted(f for n in MANIFEST_NAMES for f in raw.rglob(n) if public(f)):
        try:
            manifests[q.relative_to(root).as_posix()] = json.loads(q.read_text())
        except json.JSONDecodeError:
            manifests[q.relative_to(root).as_posix()] = {"unreadable": True}
    mpath, apath = raw.parent / MANIFESTS_DIGEST.name, raw.parent / AUDITS_DIGEST.name
    mpath.write_text(json.dumps({"generated_at": utc_now_iso(), "files": dict(sorted(manifests.items()))},
                                indent=1, ensure_ascii=False) + "\n")
    audits = [a for a in sorted(raw.glob("*/DATA_AUDIT.md")) if public(a)]
    anchor = lambda a: "source-" + a.parent.name.replace("_", "-")  # noqa: E731
    parts = ["# Source data audits\n",
             f"_All {len(audits)} `data/raw/<source_id>/DATA_AUDIT.md` files in one document, written by "
             "`measure-it merge-registry` (pipeline step `registry_merge`). Download hashes and URLs are in "
             "`data/SOURCE_MANIFESTS.json`._\n",
             "\n".join(f"* [{a.parent.name}](#{anchor(a)})" for a in audits) + "\n"]
    for a in audits:
        # demote the audit's own headings one level under the per-source heading
        body = "\n".join("#" + ln if ln.startswith("#") else ln for ln in a.read_text().strip().splitlines())
        parts.append(f'<a id="{anchor(a)}"></a>\n\n## Source: {a.parent.name}\n\n{body}\n')
    apath.write_text("\n".join(parts))
    return mpath, apath
