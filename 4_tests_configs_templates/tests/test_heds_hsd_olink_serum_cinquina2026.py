"""Tests for measure_it.labs.heds_hsd_olink_serum_cinquina2026 (ingest + pre-specified analysis)."""
from __future__ import annotations

import pandas as pd
import pytest

from measure_it.config import PROCESSED, TABLES
from measure_it.labs import heds_hsd_olink_serum_cinquina2026 as M
from measure_it.provenance import PROVENANCE_COLUMNS


def test_assay_filter_drops_duplicates_and_sparse():
    vals = pd.DataFrame({"A|1": [1.0, 2, 3, 4], "A|2": [1.0, 2, 3, 4], "B|3": [1.0, None, None, 4]})
    at = pd.DataFrame({"assay_key": ["A|1", "A|2", "B|3"], "assay": ["A", "A", "B"]})
    out = M.assay_filter(vals, at).set_index("assay_key")
    assert out.retained.to_dict() == {"A|1": True, "A|2": False, "B|3": False}


def _t(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p)


@pytest.mark.data
def test_participants_groups_and_sites():
    p = _t(f"participants__{M.SOURCE_ID}")
    assert len(p) == 352 and p.participant_id.is_unique and set(PROVENANCE_COLUMNS) <= set(p.columns)
    c = p.groupby(["group_label", "site"]).size().to_dict()
    assert c == {("HSD", "Italy"): 44, ("HSD", "USA"): 44, ("control", "Italy"): 176, ("hEDS", "Italy"): 44,
                 ("hEDS", "USA"): 44}
    assert p.loc[p.group_label == "control", "age_years"].isna().all()     # not released for controls
    assert (p.source_geographic_resolution == "none").all()


@pytest.mark.data
def test_primary_is_italian_only_and_rule_applied():
    p = TABLES / f"{M.SOURCE_ID}_primary_decision.csv"
    if not p.exists():
        pytest.skip("analysis not run")
    d = pd.read_csv(p).iloc[0]
    assert (d.n_cases, d.n_controls) == (88, 176)
    met = d.ci_low > 0.5 and d.perm_p < 0.05
    assert d.decision.startswith("the serum proteome discriminates") == met
    s = _t(f"phenotype_signatures__{M.SOURCE_ID}")
    assert s.object_id.is_unique and s.object_id.str.startswith(f"signature:{M.SOURCE_ID}|eds_hsd|").all()
