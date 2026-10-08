"""Project paths and config loading. Every module resolves paths through here."""
from __future__ import annotations

from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIGS = PROJECT_ROOT / "configs"
DATA = PROJECT_ROOT / "data"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
PROCESSED = DATA / "processed"
HTTP_CACHE = DATA / "_http_cache"
RESULTS = PROJECT_ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
MAPS = RESULTS / "maps"
DOCS = PROJECT_ROOT / "docs"
SOURCE_REGISTRY_PATH = PROJECT_ROOT / "SOURCE_REGISTRY.yaml"
DUCKDB_PATH = PROCESSED / "measure_it_public.duckdb"

# Fixed seed for every stochastic step (splits, permutations, clustering).
SEED = 20260923

UNKNOWN = "UNKNOWN / NOT AVAILABLE"


def raw_dir(source_id: str) -> Path:
    p = RAW / source_id
    p.mkdir(parents=True, exist_ok=True)
    return p


def interim_dir(source_id: str) -> Path:
    p = INTERIM / source_id
    p.mkdir(parents=True, exist_ok=True)
    return p


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@lru_cache(maxsize=None)
def load_config(name: str) -> dict:
    """Load configs/<name>.yaml."""
    with open(CONFIGS / f"{name}.yaml") as fh:
        return yaml.safe_load(fh)
