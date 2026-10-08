"""Tests for measure_it.labs.endo_arg1_repod (ingest + pre-specified analysis)."""
from __future__ import annotations

import pandas as pd
import pytest

from measure_it.config import PROCESSED, TABLES
from measure_it.labs import endo_arg1_repod as M
from measure_it.provenance import PROVENANCE_COLUMNS


def test_num_parses_na_and_commas():
    s = M.num(pd.Series(["1.5", "n/a", "N/A", "2,5", None]))
    assert s.iloc[0] == 1.5 and s.iloc[3] == 2.5 and s.iloc[[1, 2, 4]].isna().all()


def _t(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p)


@pytest.mark.data
def test_participants_and_labs():
    p = _t(f"participants__{M.SOURCE_ID}")
    assert p.participant_id.str.startswith(f"{M.SOURCE_ID}:").all() and p.participant_id.is_unique
    assert set(PROVENANCE_COLUMNS) <= set(p.columns)
    assert (p.data_layer == "person").all() and (p.source_geographic_resolution == "none").all()
    assert p.source_group.value_counts().to_dict() == {"Patient": 127, "Ctrl3": 54, "Ctrl": 25, "Ctrl2": 9}
    labs = _t(f"participant_labs__{M.SOURCE_ID}")
    assert set(labs.participant_id) <= set(p.participant_id)


@pytest.mark.data
def test_primary_decision_follows_rule():
    p = TABLES / f"{M.SOURCE_ID}_primary_decision.csv"
    if not p.exists():
        pytest.skip("analysis not run")
    d = pd.read_csv(p).iloc[0]
    assert (d.n_cases, d.n_controls) == (120, 32)
    met = d.ci_low > 0.5 and d.perm_p < 0.05
    assert d.decision.startswith("ARG1 separates") == met
    s = _t(f"phenotype_signatures__{M.SOURCE_ID}")
    assert s.object_id.str.startswith(f"signature:{M.SOURCE_ID}|endometriosis|").all() and s.object_id.is_unique
