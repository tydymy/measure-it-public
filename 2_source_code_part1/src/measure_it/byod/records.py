"""User-supplied performance records for the metric link (measure_it.scoring.metric_link.build_records).

One aggregate record per evaluated user dataset (its locked primary analysis), stored in `byod_performance_records`.
Tier `user_supplied_own_computation` (numeric 2.5): computed by this engine under a plan locked before any outcome was
computed, with CV, a participant bootstrap and a label-permutation null, on data that are not public and cannot be
re-analysed by anyone else. It ranks below this project's own computations on public data (tiers 1-2) and above
published claims that were never re-analysed (tier 3); its evidence factor is 0.6 (tier 1: 1.0, 2: 0.75, 3: 0.5).

Records tagged `demo` (the SYNTHETIC tutorial dataset) are left out unless MEASURE_IT_BYOD_DEMO=1
(`measure-it byod deploy --demo`), so they never reach a default ranking.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from ..store import read_table, table_exists
from . import common as K

USER_TIER = 2.5
USER_TIER_FACTOR = 0.6
RECORD_FIELDS = ("record_id", "source_kind", "quality_tier", "target_condition", "target_members", "primary_class",
                 "dataset_id", "analysis", "label_basis", "comparator", "comparator_kind", "n_cases", "n_controls",
                 "auroc", "auroc_ci_low", "auroc_ci_high", "auroc_basis", "op_sensitivity", "op_sensitivity_ci_low",
                 "op_sensitivity_ci_high", "op_specificity", "op_specificity_ci_low", "op_specificity_ci_high",
                 "op_threshold", "op_basis", "measurement_only", "exclusion_reason", "source_object_ids",
                 "source_tables", "caveats", "computed_here")
EXTRA_FIELDS = ("user_supplied", "demo", "synthetic", "plan_sha256")
JSON_FIELDS = ("target_members", "source_object_ids", "source_tables", "caveats")


def stored_records(include_demo: bool | None = None) -> pd.DataFrame:
    """Rows of byod_performance_records (demo rows only when enabled)."""
    if not table_exists(K.RECORDS_TABLE):
        return pd.DataFrame()
    df = read_table(K.RECORDS_TABLE)
    include_demo = K.demo_enabled() if include_demo is None else include_demo
    if not include_demo and "demo" in df.columns:
        df = df[~df["demo"].astype(bool)]
    return df


def performance_records(include_demo: bool | None = None) -> list[dict]:
    """Keyword dicts for metric_link._rec (list-valued fields decoded from JSON)."""
    out = []
    for r in stored_records(include_demo).to_dict("records"):
        kw = {k: r.get(k) for k in RECORD_FIELDS if k in r}
        for k in JSON_FIELDS:
            v = kw.get(k)
            kw[k] = json.loads(v) if isinstance(v, str) and v.startswith("[") else (v if isinstance(v, list) else [])
        for k in ("n_cases", "n_controls", "auroc", "auroc_ci_low", "auroc_ci_high", "op_sensitivity",
                  "op_sensitivity_ci_low", "op_sensitivity_ci_high", "op_specificity", "op_specificity_ci_low",
                  "op_specificity_ci_high", "op_threshold"):
            v = kw.get(k)
            kw[k] = float(v) if v is not None and pd.notna(v) else np.nan
        kw["quality_tier"] = USER_TIER
        kw["measurement_only"] = bool(kw.get("measurement_only", True))
        if not kw.get("exclusion_reason"):
            kw["exclusion_reason"] = None
        for k in EXTRA_FIELDS:
            if k in r:
                kw[k] = bool(r[k]) if k != "plan_sha256" else str(r[k])
        out.append(kw)
    return out


def extra_grid_cells(existing: set[tuple[str, str]]) -> list[tuple[str, str]]:
    """(condition, class) cells of user records outside the scored grid, so measurement_performance shows them."""
    cells = []
    for r in performance_records():
        cell = (str(r["target_condition"]), str(r["primary_class"]))
        if cell not in existing and cell not in cells:
            cells.append(cell)
    return cells
