"""BYOD shared pieces: identifiers, local-only paths, manifest loading and input-file discovery.

Identifiers
    manifest `dataset_id` (e.g. `synthetic_demo`)  ->  engine dataset id `byod_synthetic_demo`
    participant ids                                  ->  `byod_synthetic_demo:<native id or salted hash>`
    partitions                                       ->  `<table>__byod_synthetic_demo`

Local-only locations (each gets its own `.gitignore` containing `*`, so nothing in them can be committed even before
the repository's root .gitignore lists them):
    data/raw/byod_<id>/        registry entry, generated DATA_AUDIT.md, input manifest (file names, sizes, sha256)
    results/byod/<id>/         locked plan, evaluation tables, deploy report, heatmap / embedding data and figures
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from ..config import PROCESSED, RAW, RESULTS

PREFIX = "byod_"
DATASETS_TABLE = "byod_datasets"                  # one row per user dataset (metadata layer)
RECORDS_TABLE = "byod_performance_records"        # aggregate performance records read by the metric link
DEMO_ENV = "MEASURE_IT_BYOD_DEMO"                 # "1": demo-tagged records enter the metric link
SALT_ENV = "MEASURE_IT_BYOD_SALT"                 # salt for id_hashing: salted_sha256 (never written to disk)
PROVENANCE_NOTE = "user-supplied; not redistributed"
SOURCE_KIND = "user_supplied_own_computation"
RESULTS_ROOT = RESULTS / "byod"
LEDGER = RESULTS_ROOT / "LEDGER.jsonl"            # append-only record of every plan lock and evaluation
SCHEMA_VERSION = 1
MANIFEST_ID_RE = re.compile(r"^[a-z][a-z0-9_]{2,40}$")

# person-layer tables a user dataset may write (partition suffix __byod_<id>)
PERSON_TABLES = ["participants", "participant_conditions", "participant_labs", "participant_omics_linked",
                 "participant_wearable_features", "participant_wearable_daily", "participant_device_features",
                 "participant_adapter_scores", "phenotype_signatures", "participant_medications",
                 "participant_survey", "person_concepts"]
# aggregate (non-partition) tables written by `byod subgroup` (no person rows, no participant ids)
PHENOTYPES_TABLE = "computable_phenotypes"
STRATA_TABLE = "subgroup_strata_fractions"
COMPARATOR_TYPES = {"healthy": "healthy", "look_alike": "look-alike", "general_population": "general population"}
LABEL_BASES = {"clinical_case_definition", "clinical_diagnosis", "ehr_codes", "self_report", "proxy"}
BLOCK_KINDS = ("device", "ehr_labs", "ehr_diagnoses", "ehr_medications", "omics", "survey", "self_report_history")
BINARY_BLOCK_KINDS = ("ehr_diagnoses", "ehr_medications", "self_report_history")


def engine_id(manifest_id: str) -> str:
    """Manifest dataset_id -> engine dataset id (idempotent)."""
    return manifest_id if manifest_id.startswith(PREFIX) else PREFIX + manifest_id


def manifest_id_of(dataset_id: str) -> str:
    return dataset_id[len(PREFIX):] if dataset_id.startswith(PREFIX) else dataset_id


def local_only_dir(path: Path) -> Path:
    """mkdir + a `.gitignore` of `*`: git ignores everything below, whatever the root .gitignore says."""
    path.mkdir(parents=True, exist_ok=True)
    gi = path / ".gitignore"
    if not gi.exists():
        gi.write_text("# measure-it byod: local-only, user-supplied material; never commit\n*\n")
    return path


def raw_path(dataset_id: str) -> Path:
    return RAW / engine_id(dataset_id)


def results_path(dataset_id: str) -> Path:
    return RESULTS_ROOT / manifest_id_of(dataset_id)


def partition_name(table: str, dataset_id: str) -> str:
    return f"{table}__{engine_id(dataset_id)}"


def dataset_partitions(dataset_id: str) -> list[Path]:
    """Every processed file of this dataset (partitions and their sidecars)."""
    ds = engine_id(dataset_id)
    return sorted(p for p in PROCESSED.glob(f"*__{ds}.*") if p.name.split("__", 1)[1].split(".")[0] == ds)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_json(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


# ------------------------------------------------------------------------------------------------------------------
# input folder
# ------------------------------------------------------------------------------------------------------------------

@dataclass
class InputFile:
    path: Path
    kind: str            # participants | device | ehr_labs | ehr_diagnoses | ehr_medications | omics | survey |
                         # self_report_history
    block: str | None    # device_<name> | ehr_labs | ehr_diagnoses | ehr_medications | omics_<layer> |
                         # survey_<instrument> | self_report_history (None for participants)
    name: str            # device name / omics layer / survey instrument file name ('' otherwise)


@dataclass
class Inputs:
    root: Path
    manifest_path: Path
    files: list[InputFile] = field(default_factory=list)
    unknown_data_files: list[Path] = field(default_factory=list)
    ignored: list[Path] = field(default_factory=list)

    def by_kind(self, kind: str) -> list[InputFile]:
        return [f for f in self.files if f.kind == kind]

    def blocks(self) -> dict[str, InputFile]:
        return {f.block: f for f in self.files if f.block}

    @property
    def participants(self) -> InputFile | None:
        p = self.by_kind("participants")
        return p[0] if p else None


DATA_SUFFIXES = {".csv", ".parquet", ".tsv", ".xlsx", ".xls", ".json", ".txt", ".sav", ".dta", ".sas7bdat"}
NAME_RE = re.compile(r"^[a-z][a-z0-9_]{0,40}$")


def discover(root: Path) -> Inputs:
    """Classify the files of an input folder by the contract's file names."""
    root = Path(root)
    inp = Inputs(root=root, manifest_path=root / "manifest.yaml")
    for p in sorted(root.iterdir()):
        if p.is_dir() or p.name.startswith("."):
            continue
        stem, suf = p.stem.lower(), p.suffix.lower()
        if p.name in ("manifest.yaml", "manifest.yml"):
            inp.manifest_path = p
            continue
        if p.name in ("manifest.schema.json", "computable_phenotype.schema.json"):   # template docs, not data
            inp.ignored.append(p)
            continue
        if p.name == "participants.csv":
            inp.files.append(InputFile(p, "participants", None, ""))
        elif p.name == "ehr_labs.csv":
            inp.files.append(InputFile(p, "ehr_labs", "ehr_labs", ""))
        elif p.name == "ehr_diagnoses.csv":
            inp.files.append(InputFile(p, "ehr_diagnoses", "ehr_diagnoses", ""))
        elif p.name == "ehr_medications.csv":
            inp.files.append(InputFile(p, "ehr_medications", "ehr_medications", ""))
        elif p.name == "self_report_history.csv":
            inp.files.append(InputFile(p, "self_report_history", "self_report_history", ""))
        elif stem.startswith("survey_") and suf in (".csv", ".parquet"):
            name = stem[len("survey_"):]
            inp.files.append(InputFile(p, "survey", f"survey_{name}", name))
        elif stem.startswith("device_") and suf in (".csv", ".parquet"):
            name = stem[len("device_"):]
            inp.files.append(InputFile(p, "device", f"device_{name}", name))
        elif stem.startswith("omics_") and suf in (".csv", ".parquet"):
            name = stem[len("omics_"):]
            inp.files.append(InputFile(p, "omics", f"omics_{name}", name))
        elif suf in DATA_SUFFIXES:
            inp.unknown_data_files.append(p)
        else:
            inp.ignored.append(p)
    return inp


def load_manifest(path: Path) -> dict:
    with open(path) as fh:
        m = yaml.safe_load(fh)
    if not isinstance(m, dict):
        raise ValueError(f"{path.name}: not a YAML mapping")
    return m


def read_frame(path: Path):
    import pandas as pd
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    # everything as text first: the privacy checks see values exactly as supplied (leading zeros, dates)
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[""])


def analysis_defaults(pa: dict | None) -> dict:
    """The primary-analysis settings with engine defaults filled in (part of the locked plan)."""
    pa = dict(pa or {})
    pa.setdefault("combined_blocks", [])
    pa.setdefault("covariates", [])
    pa.setdefault("specificity", 0.90)
    pa.setdefault("n_splits", 5)
    pa.setdefault("n_repeats", 10)
    pa.setdefault("n_permutations", 200)
    pa.setdefault("n_bootstrap", 1000)
    pa.setdefault("C", 1.0)
    return pa


def subgroup_defaults(sa: dict | None) -> dict:
    """The subgroup-analysis settings with engine defaults filled in (part of the locked subgroup plan)."""
    sa = dict(sa or {})
    sa.setdefault("description", "")
    sa.setdefault("ehr_domains", ["condition", "measurement", "drug", "demographic", "survey"])
    sa.setdefault("max_features", 10)
    sa.setdefault("min_feature_participants", 5)
    sa.setdefault("exclude_feature_keys", [])
    sa.setdefault("n_splits", 5)
    sa.setdefault("n_repeats", 10)
    sa.setdefault("n_permutations", 200)
    sa.setdefault("n_bootstrap", 1000)
    sa.setdefault("C", 1.0)
    dp = dict(sa.get("device_positive") or {})
    if "block" in dp:
        dp.setdefault("direction", ">=")
    sa["device_positive"] = dp
    return sa


def demo_enabled() -> bool:
    import os
    return os.environ.get(DEMO_ENV, "") == "1"
