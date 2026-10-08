"""`measure-it byod remove <id>`: remove a user dataset and every record it produced, then restore the rankings.

Deletes the dataset's partitions (and sidecars), its rows in byod_datasets, byod_performance_records,
computable_phenotypes and subgroup_strata_fractions, its
data/raw/byod_<id>/ and results/byod/<id>/ folders; rebuilds the canonical unions and the Digital Phenotype Vectors;
rebuilds the metric link and re-scores every combination whose performance row changed (so the rankings return to
what they were without the dataset). results/byod/LEDGER.jsonl keeps a line saying the dataset was removed (the
evaluation history is not erased; `--purge-ledger` removes the dataset's lines too).
"""
from __future__ import annotations

import json
import shutil

import pandas as pd
import yaml

from ..config import SOURCE_REGISTRY_PATH
from ..store import read_table, table_exists
from . import common as K


def remove(dataset_id: str, *, rescore: bool = True, purge_ledger: bool = False, build_db: bool = False,
           echo=print) -> dict:
    ds = K.engine_id(dataset_id)
    from .person import dataset_row
    known = dataset_row(ds) is not None or bool(K.dataset_partitions(ds)) or K.raw_path(ds).exists()
    if not known:
        raise RuntimeError(f"no user dataset {ds} on this machine (see `measure-it byod list`)")
    removed = []
    for p in K.dataset_partitions(ds):
        p.unlink()
        removed.append(p.name)
    from .evaluate import write_records_table
    from .ingest import rebuild_person_layer, write_datasets_table
    if table_exists(K.RECORDS_TABLE):
        r = read_table(K.RECORDS_TABLE)
        write_records_table(r[r["dataset_id"] != ds].drop(columns=[c for c in r.columns if c in _PROV]))
    from .subgroup import drop_dataset
    drop_dataset(ds)                     # computable_phenotypes / subgroup_strata_fractions rows
    from ..similar.reference import drop_dataset as drop_similar
    removed += drop_similar(ds)          # similar_reference_* rows of its phenotypes (byod similar)
    if table_exists(K.DATASETS_TABLE):
        d = read_table(K.DATASETS_TABLE)
        write_datasets_table(d[d["dataset_id"] != ds].drop(columns=[c for c in d.columns if c in _PROV]))
    for folder in (K.raw_path(ds), K.results_path(ds)):
        if folder.exists():
            shutil.rmtree(folder)
            removed.append(str(folder.name) + "/")
    # never expected (the entry is not a SOURCE_REGISTRY fragment), but if someone merged it by hand, take it out
    reg = yaml.safe_load(SOURCE_REGISTRY_PATH.read_text()) if SOURCE_REGISTRY_PATH.exists() else None
    if reg and any(s.get("source_id") == ds for s in reg.get("sources", [])):
        reg["sources"] = [s for s in reg["sources"] if s.get("source_id") != ds]
        SOURCE_REGISTRY_PATH.write_text(yaml.safe_dump(reg, sort_keys=False, allow_unicode=True, width=110))
        removed.append("SOURCE_REGISTRY.yaml entry")
    if K.LEDGER.exists() and purge_ledger:
        lines = [ln for ln in K.LEDGER.read_text().splitlines() if json.loads(ln).get("dataset_id") != ds]
        K.LEDGER.write_text("\n".join(lines) + ("\n" if lines else ""))
    else:
        from .evaluate import _ledger
        _ledger({"event": "removed", "dataset_id": ds})
    echo(f"byod remove {ds}: deleted {len(removed)} item(s)")
    res = {"dataset_id": ds, "removed": removed}
    res["rebuild"] = rebuild_person_layer(echo=echo)
    _clean_canonicals(ds)
    if rescore:
        from .deploy import demo_state, refresh_rankings
        from .records import stored_records
        left = stored_records(include_demo=True)
        demo = demo_state() and bool(len(left)) and bool(left.get("demo", pd.Series(dtype=bool)).astype(bool).any())
        ref = refresh_rankings(demo, echo=echo)   # demo mode is kept only while a demo record remains
        res["rescored"] = ref["rescored"]
        res["changed_cells"] = ref["changed_cells"]
    if build_db:
        from ..database import build_db as _b
        _b()
    return res


_PROV = ("data_layer", "source_name", "source_record_id", "source_version", "retrieved_at",
         "source_geographic_resolution", "evidence_type", "evidence_level", "provenance_notes")


def _clean_canonicals(ds: str) -> None:
    """A canonical table whose family has no partition left keeps no row of the removed dataset."""
    from ..config import PROCESSED
    from ..store import partitions, write_table
    for table in K.PERSON_TABLES:
        p = PROCESSED / f"{table}.parquet"
        if partitions(table) or not p.exists():
            continue
        t = pd.read_parquet(p)
        if "dataset_id" in t.columns and (t["dataset_id"].astype(str) == ds).any():
            if (t["dataset_id"].astype(str) == ds).all():     # a family only user datasets create (participant_survey)
                for suf in (".parquet", ".meta.json"):
                    (PROCESSED / f"{table}{suf}").unlink(missing_ok=True)
                continue
            write_table(t[t["dataset_id"].astype(str) != ds], table, producer="measure_it.byod.remove",
                        description=f"{table} without the removed user dataset {ds}")
