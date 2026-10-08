"""Tool-facade functions for user datasets (same envelope as measure_it.tools; aggregate data only).

    list_user_datasets()                  -> the user datasets on this machine: manifest summary, status, plan hash
    get_user_dataset_metric(dataset_id)   -> the locked plan, the evaluation models and the performance record

Neither returns a person row or a participant id. When nothing is ingested the status is UNKNOWN / NOT AVAILABLE.
"""
from __future__ import annotations

import json

import pandas as pd

from .. import tools as T
from ..config import UNKNOWN
from ..store import read_table, table_exists
from . import common as K

CAVEAT = ("User-supplied, local-only data: not public and not re-analysable by others; the performance record is "
          "tier user_supplied_own_computation (2.5). Against healthy controls, performance is optimistic.")
KEEP = ("dataset_id", "title", "owner", "licence", "synthetic", "demo", "status", "comparator_type",
        "primary_label", "primary_condition_id", "measurement_class", "measurement_bundle", "n_participants",
        "n_cases_primary", "n_controls_primary", "plan_sha256", "ingested_at", "evaluated_at", "deployed_at",
        "deployed_with_demo")


def _datasets() -> pd.DataFrame:
    return read_table(K.DATASETS_TABLE) if table_exists(K.DATASETS_TABLE) else pd.DataFrame()


def list_user_datasets() -> dict:
    """User-supplied (bring-your-own-data) datasets on this machine: title, owner, status, comparator, condition,
    measurement class, group sizes and the locked plan's hash. UNKNOWN / NOT AVAILABLE when none is ingested."""
    def fn():
        d = _datasets()
        if d.empty:
            return T._env("list_user_datasets", UNKNOWN, {"datasets": []}, caveats=[CAVEAT],
                          reason="no user dataset has been ingested (measure-it byod ingest <dir>)")
        rows = [{k: r.get(k) for k in KEEP} for r in d.to_dict("records")]
        return T._env("list_user_datasets", T.OK, {"datasets": rows, "n": len(rows)}, caveats=[CAVEAT])
    return T._run("list_user_datasets", {}, fn)


def get_user_dataset_metric(dataset_id: str) -> dict:
    """One user dataset's locked analysis plan, evaluation models (AUROC / AUPRC with CIs, sensitivity at the
    pre-specified specificity, permutation p), its tier-2.5 performance record and the last deploy report. Aggregate
    results only: no person rows, no participant ids."""
    ds = K.engine_id(str(dataset_id))

    def fn():
        d = _datasets()
        if d.empty or not (d["dataset_id"] == ds).any():
            return T._env("get_user_dataset_metric", UNKNOWN, None, caveats=[CAVEAT],
                          reason=f"no user dataset {ds!r} (list_user_datasets shows what exists)")
        meta = {k: v for k, v in d[d["dataset_id"] == ds].iloc[0].to_dict().items() if k in KEEP}
        out = {"dataset": meta}
        lock = K.results_path(ds) / "plan_lock.json"
        out["locked_plan"] = json.loads(lock.read_text()) if lock.exists() else UNKNOWN
        mp = K.results_path(ds) / "evaluation_models.csv"
        out["models"] = pd.read_csv(mp).to_dict("records") if mp.exists() else UNKNOWN
        rec = read_table(K.RECORDS_TABLE) if table_exists(K.RECORDS_TABLE) else pd.DataFrame()
        rec = rec[rec["dataset_id"] == ds] if len(rec) else rec
        if len(rec):
            r = rec.iloc[0].to_dict()
            for k in ("caveats", "source_object_ids", "target_members", "source_tables"):
                if isinstance(r.get(k), str) and r[k].startswith("["):
                    r[k] = json.loads(r[k])
            out["performance_record"] = r
        else:
            out["performance_record"] = UNKNOWN
        rp = K.results_path(ds) / "deploy_report.json"
        out["deploy_report"] = json.loads(rp.read_text()) if rp.exists() else UNKNOWN
        status = T.OK if len(rec) else UNKNOWN
        return T._env("get_user_dataset_metric", status, out, caveats=[CAVEAT],
                      reason=None if len(rec) else "dataset ingested but not evaluated (measure-it byod evaluate)")
    return T._run("get_user_dataset_metric", {"dataset_id": ds}, fn)
