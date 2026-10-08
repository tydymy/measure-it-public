"""Tests for measure_it.omics.geo_cohorts (GEO case/control omics cohorts for Long COVID / ME/CFS).

Pure-function tests run anywhere; @pytest.mark.data tests check the outputs of
    uv run python -m measure_it.omics.geo_cohorts
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from measure_it.config import RAW, TABLES, load_config
from measure_it.omics import geo_cohorts as G
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.store import read_table, table_exists

GSM_TEXT = """^SAMPLE = GSM1
!Sample_title = S01.100.PASC
!Sample_source_name_ch1 = Whole blood
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = status: PASC
!Sample_characteristics_ch1 = draw: 3_6M
!Sample_platform_id = GPL1
^SAMPLE = GSM2
!Sample_title = S01.200.PASC
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = status: PASC
!Sample_characteristics_ch1 = draw: 6_12M
!Sample_platform_id = GPL1
^SAMPLE = GSM3
!Sample_title = S02.50.Recovered
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = status: Recovered
!Sample_characteristics_ch1 = draw: Acute
!Sample_platform_id = GPL1
^SAMPLE = GSM4
!Sample_title = S03.90.Recovered
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = status: Recovered
!Sample_characteristics_ch1 = draw: 3_6M
!Sample_platform_id = GPL2
"""


def test_parse_gsm_text_characteristics_become_columns():
    df = G.parse_gsm_text(GSM_TEXT)
    assert df["gsm"].tolist() == ["GSM1", "GSM2", "GSM3", "GSM4"]
    assert df.loc[0, "c:status"] == "PASC" and df.loc[2, "c:draw"] == "Acute"
    assert df.loc[3, "platform_id"] == "GPL2"


def test_label_rule_filters_excludes_acute_and_keeps_earliest_per_subject():
    meta = G.parse_gsm_text(GSM_TEXT)
    rule = {"filter": {"platform_id": ["GPL1"]}, "label": {"field": "c:status", "case": ["PASC"], "control": ["Recovered"]},
            "exclude_timepoints": {"field": "c:draw", "values": ["Acute"]},
            "subject": {"from": "title", "regex": r"^([^.]+)\."},
            "order": {"from": "title", "regex": r"^[^.]+\.(\d+)\.", "numeric": True}}
    out = G.apply_label_rule(meta, rule)
    # GSM4 is on another platform, GSM3 is acute; S01 keeps its earliest draw (day 100)
    assert out["gsm"].tolist() == ["GSM1"] and out["label"].tolist() == [1]


def test_label_rule_twin_pairs_are_completed():
    meta = pd.DataFrame({"gsm": ["a", "b", "c", "d", "e"], "c:dx": ["CFS", "unaffected", "ICF", "unaffected", "CFS"],
                         "c:pair": ["1", "1", "2", "2", "3"]})
    rule = {"label": {"field": "c:dx", "case": ["CFS"], "control": ["unaffected"]}, "pair": {"field": "c:pair"}}
    out = G.apply_label_rule(meta, rule)
    # pair 2 has no CFS case, pair 3 has no control: only pair 1 survives
    assert out["gsm"].tolist() == ["a", "b"] and out["pair"].nunique() == 1


def test_permute_labels_within_pairs_keeps_one_case_per_pair():
    y = np.array([1, 0, 0, 1, 1, 0])
    g = np.array(["p1", "p1", "p2", "p2", "p3", "p3"])
    rng = np.random.default_rng(0)
    for _ in range(20):
        yp = G.permute_labels(y, g, rng)
        assert all(yp[g == k].sum() == 1 for k in np.unique(g))


def test_permute_labels_unpaired_is_a_permutation():
    y = np.array([1, 1, 0, 0, 0])
    yp = G.permute_labels(y, np.array(list("abcde")), np.random.default_rng(1))
    assert sorted(yp) == sorted(y)


def test_welch_t_sign_and_zero_variance():
    X = np.array([[2.0, 1.0], [3.0, 1.0], [0.0, 1.0], [1.0, 1.0]])
    t = G.welch_t(X, np.array([1, 1, 0, 0]))
    assert t[0] > 0 and t[1] == 0


def test_sex_chromosome_features_are_flagged():
    assert G._is_autosomal("7") and G._is_autosomal(None) and G._is_autosomal("Un")
    assert not G._is_autosomal("X") and not G._is_autosomal("Y") and not G._is_autosomal("X|Y") and not G._is_autosomal("MT")


def test_parse_rcc_code_summary():
    txt = "<Code_Summary>\nCodeClass,Name,Accession,Count\nEndogenous,CD3E,NM_1,120\nPositive,POS_A,ERCC,5000\n</Code_Summary>"
    df = G.parse_rcc(txt)
    assert df["Count"].tolist() == [120.0, 5000.0] and df["CodeClass"].tolist() == ["Endogenous", "Positive"]


def test_plan_is_fixed_and_single_specification():
    plan = load_config("geo_cohorts")["analysis_plan"]
    c = plan["classifier"]
    assert (c["k_features"], c["C"], c["variance_filter_top"], c["n_permutations"]) == (50, 0.1, 5000, 1000)
    tiers = {p["tier"] for p in plan["transfer"]["pairs"]}
    assert tiers == {"primary", "secondary"}


# --------------------------------------------------------------------------- outputs
@pytest.mark.data
def test_every_catalog_series_is_screened_with_a_reason():
    scr = pd.read_csv(TABLES / G.SCREEN_CSV)
    uni = G.catalog_universe()
    assert set(uni["gse"]) <= set(scr["gse"])
    assert "NOT_IN_CONFIG" not in set(scr["fail_reasons"].dropna())
    fails = scr[scr["decision"] == "fail"]
    assert fails["fail_reasons"].fillna("").str.len().gt(0).all()
    passes = scr[scr["decision"] == "pass"]
    assert (passes[["n_case_subjects", "n_control_subjects"]].min(axis=1) >= 15).all()


@pytest.mark.data
def test_selection_follows_the_rule():
    sel = pd.read_csv(TABLES / "geo_cohort_selection.csv")
    prim = sel[sel["role"] == "primary"]
    assert prim.groupby("condition_id").size().max() <= 3
    assert sel[sel["role"] == "replication_only"].groupby("condition_id").size().max() <= 1
    for cond, d in prim.groupby("condition_id"):
        assert d["n_analysable"].is_monotonic_decreasing


@pytest.mark.data
def test_classifier_results_are_complete_and_bounded():
    cls = pd.read_csv(TABLES / "geo_cohort_classifier.csv")
    sel = pd.read_csv(TABLES / "geo_cohort_selection.csv")
    assert set(cls["subset_id"]) == set(sel["subset_id"])
    assert cls["perm_p"].between(1 / 1001 - 1e-12, 1).all() and (cls["n_permutations"] == 1000).all()
    assert (cls["auroc_ci_low"] <= cls["auroc_ci_high"]).all()
    null = pd.read_csv(TABLES / "geo_cohort_permutation_null.csv")
    assert null.groupby("subset_id").size().eq(1000).all()


@pytest.mark.data
def test_transfer_replication_flag_matches_rule():
    tr = pd.read_csv(TABLES / "geo_cohort_transfer.csv")
    assert (tr["tier"] == "primary").sum() == 4
    ok = tr.dropna(subset=["transfer_auroc"])
    assert (ok["replicated"] == ((ok["ci_low"] > 0.5) & (ok["perm_p"] < 0.05))).all()


@pytest.mark.data
@pytest.mark.skipif(not table_exists(G.SIGNATURE_TABLE), reason="run the module first")
def test_signatures_only_for_permutation_passing_subsets():
    sig = read_table(G.SIGNATURE_TABLE)
    cls = pd.read_csv(TABLES / "geo_cohort_classifier.csv")
    passing = {f"geo_{s}" for s in cls.loc[cls["perm_p"] < 0.05, "subset_id"]}
    assert set(sig["dataset_id"]) <= passing
    assert all(c in sig.columns for c in PROVENANCE_COLUMNS)
    if len(sig):
        assert sig["object_id"].str.startswith("signature:").all()
        assert (sig["q_value_bh"] < 0.05).all()
        assert set(sig["replication_status"]) <= {"dataset_specific", "replicated_transfer"}
        assert (sig["data_layer"] == "person").all()


@pytest.mark.data
def test_registry_entry_and_audit_exist():
    assert (RAW / G.SOURCE_ID / "registry_entry.yaml").exists()
    audit = (RAW / G.SOURCE_ID / "DATA_AUDIT.md").read_text()
    assert "True participant linkage across modalities? | no" in audit


@pytest.mark.data
def test_ci_is_reported_with_the_estimator_it_brackets():
    """The pre-specified bootstrap CI is of the repeat-averaged-score AUROC, not of the mean per-repeat AUROC."""
    cls = pd.read_csv(TABLES / "geo_cohort_classifier.csv")
    assert "auroc_avg_score" in cls.columns
    assert (cls["auroc_avg_score"] >= cls["auroc_ci_low"] - 1e-9).all()
    assert (cls["auroc_avg_score"] <= cls["auroc_ci_high"] + 1e-9).all()
    cov = cls.dropna(subset=["cov_only_auroc"])
    # a below-chance covariates-only CV AUROC is flagged, so a delta against it is not read as an omics increment
    assert (cov["cov_only_below_chance"] == (cov["cov_only_auroc"] < 0.5)).all()


@pytest.mark.data
def test_confound_rows_name_their_test():
    conf = pd.read_csv(TABLES / "geo_cohort_confounds.csv")
    assert conf["test"].notna().all()
    # a 2x2 table is Fisher-exact; anything larger (e.g. the 3-prefix GSE270045 table) is chi-square, not Fisher
    n_rows = conf["table"].map(lambda t: len(json.loads(t)) - ("concordance_with_recorded_sex" in t))
    fisher = conf["test"] == "fisher_exact"
    assert (n_rows[fisher] == 2).all()
    assert conf.loc[conf["test"].str.startswith("chi_square"), "table"].map(lambda t: len(json.loads(t))).gt(2).all()
