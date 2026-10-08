"""Tests for measure_it.labs.klein2023_mylc_ml_table (ingest + pre-specified analysis)."""
from __future__ import annotations

import pandas as pd
import pytest

from measure_it.config import PROCESSED, TABLES
from measure_it.labs import klein2023_mylc_ml_table as M
from measure_it.provenance import PROVENANCE_COLUMNS


def test_analyte_name_and_unit_parsing():
    assert M.analyte_name(M.CORT) == "Cortisol" and M.analyte_unit(M.CORT) == "ng/mL"
    assert M.participant_id("LC.MS.0001") == "klein2023_mylc_ml_table:LC.MS.0001"


def test_floor_is_half_smallest_positive():
    assert M._floor([0, 2.0, 4.0, None]) == 1.0 and M._floor([0, 0]) == 1e-9


def _t(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p)


@pytest.mark.data
def test_analysis_set_matches_paper():
    p = _t(f"participants__{M.SOURCE_ID}")
    assert len(p) == 185 and p.participant_id.is_unique and set(PROVENANCE_COLUMNS) <= set(p.columns)
    a = p[p.in_analysis_set]
    assert a.group_label.value_counts().to_dict() == {"LC": 99, "HC": 40, "CC": 39}
    assert (p.source_geographic_resolution == "none").all()
    labs = _t(f"participant_labs__{M.SOURCE_ID}")
    assert labs.lab_variable.nunique() == 144


@pytest.mark.data
def test_primary_decision_follows_rule():
    p = TABLES / f"{M.SOURCE_ID}_primary_decision.csv"
    if not p.exists():
        pytest.skip("analysis not run")
    d = pd.read_csv(p).iloc[0]
    met = d.delta_ci_low > 0 and d.cortisol_only_perm_p < 0.05
    assert d.decision.startswith("cortisol adds") == met
    s = _t(f"phenotype_signatures__{M.SOURCE_ID}")
    assert s.object_id.is_unique and s.object_id.str.startswith(f"signature:{M.SOURCE_ID}|long_covid|").all()
