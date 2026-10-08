"""Tests for measure_it.labs.lab_dataset_discovery: candidate table integrity and ranking rules."""
from __future__ import annotations

import math

import pandas as pd
import pytest

from measure_it.config import TABLES
from measure_it.labs import lab_dataset_discovery as D

STATUSES = {"verified_open", "verified_registration_only", "controlled", "not_as_described", "unreachable", "duplicate",
            "excluded_privacy"}


def test_privacy_excluded_records_are_never_downloaded():
    for cid in D.PRIVACY_EXCLUDED:
        assert not D.FILES.get(cid), cid


def test_counts_helper_counts_people_not_values():
    lab = pd.Series(["a", "a", "b", "c", None])
    has = pd.Series([True, False, True, True, True])
    c = D._counts(lab, ["a"], ["b", "c"], has)
    assert (c["n_cases_file"], c["n_controls_file"], c["n_cases_with_lab"], c["n_controls_with_lab"]) == (2, 2, 1, 2)


@pytest.mark.data
def test_candidate_table_contract():
    p = TABLES / "lab_dataset_candidates.csv"
    if not p.exists():
        pytest.skip("run measure_it.labs.lab_dataset_discovery first")
    df = pd.read_csv(p)
    assert df.candidate_id.is_unique
    assert set(df.access_status) <= STATUSES
    el = df[df.rank_eligible == True]  # noqa: E712
    assert (el.access_status == "verified_open").all()
    assert (el.n_cases_with_lab >= D.MIN_PER_GROUP).all() and (el.n_controls_with_lab >= D.MIN_PER_GROUP).all()
    assert el.label_mapping_unverified.isna().all() and el.circular_label.isna().all()
    comp = el[D.SCORE_COLS].sum(axis=1)
    assert ((comp - el.rank_score).abs() < 1e-6).all()
    for r in el.itertuples():
        assert math.isclose(r.score_n, round(min(3.0, 1.5 * math.log10(min(r.n_cases_with_lab, r.n_controls_with_lab))), 3))
    excl = df[df.access_status == "excluded_privacy"]
    assert excl.raw_dir.isna().all() or (excl.raw_dir == "").all()
