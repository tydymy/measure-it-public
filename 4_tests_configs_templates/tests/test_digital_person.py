"""Tests for the digital person representation (canonical unions + per-dataset Digital Phenotype Vector).

Unit tests use small synthetic frames created inside this module only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from measure_it import config
from measure_it.config import UNKNOWN
from measure_it.store import partitions, read_table, table_exists
from measure_it.wearables import digital_person as DP

MAPMECFS_LABEL = ("participant-linked omics exist for 47 participants but NO wearable/physiology data is linkable in "
                  "open data, so this is not a wearable+omics multimodal person")


# ------------------------------------------------------------------------------------------ unit
def test_transform_feature_kinds():
    dims, t = DP.transform_feature(pd.Series(["yes", "no", None]), DP.Feat("b", "x", "binary", "t"))
    assert dims["b"].tolist()[:2] == [1.0, 0.0] and np.isnan(dims["b"].iloc[2])
    dims, t = DP.transform_feature(pd.Series(["a", "b", None, UNKNOWN]), DP.Feat("c", "x", "categorical", "t"))
    assert set(dims) == {"c=a", "c=b"} and dims["c=a"].isna().tolist() == [False, False, True, True]
    dims, t = DP.transform_feature(pd.Series([0.0, 6.0, 12.0]), DP.Feat("h", "x", "clock_hours", "t"))
    assert dims["h__sin"].round(6).tolist() == [0.0, 1.0, 0.0] and dims["h__cos"].round(6).tolist() == [1.0, 0.0, -1.0]
    skewed = pd.Series(np.r_[np.zeros(50), [1000.0]])
    dims, t = DP.transform_feature(skewed, DP.Feat("s", "x", "numeric", "t"))
    assert t.startswith("log1p") and dims["s"].max() == pytest.approx(np.log1p(1000.0))
    dims, t = DP.transform_feature(pd.Series([1.0, 2.0, 3.0]), DP.Feat("n", "x", "numeric", "t"))
    assert t == "identity"


def _synthetic_spec(n=60, seed=1):
    rng = np.random.default_rng(seed)
    idx = pd.Index([f"synthdp:{i}" for i in range(n)], name="participant_id")
    x1 = rng.normal(10, 2, n)
    x1[:6] = np.nan  # 10 % missing
    frame = pd.DataFrame({
        "x1": x1, "x2": rng.normal(size=n), "b1": rng.integers(0, 2, n).astype(float),
        "c1": rng.choice(["p", "q", "r"], n), "lab1": rng.integers(0, 2, n).astype(float),
        "mostly_missing": np.r_[[1.0, 2.0], np.full(n - 2, np.nan)], "constant": np.ones(n),
        "w1": rng.normal(size=n), "w2": rng.normal(size=n)}, index=idx)
    feats = [DP.Feat("x1", "labs", "numeric", "synthetic"), DP.Feat("x2", "labs", "numeric", "synthetic"),
             DP.Feat("b1", "symptoms_conditions", "binary", "synthetic"),
             DP.Feat("c1", "demographics", "categorical", "synthetic"),
             DP.Feat("lab1", "symptoms_conditions", "binary", "synthetic", label_defining=True),
             DP.Feat("mostly_missing", "labs", "numeric", "synthetic"),
             DP.Feat("constant", "labs", "numeric", "synthetic"),
             DP.Feat("w1", "wearable", "numeric", "synthetic"), DP.Feat("w2", "wearable", "numeric", "synthetic"),
             DP.Feat("not_there", "labs", "numeric", "synthetic")]
    return DP.DatasetSpec(dataset_id="synthdp", frame=frame, feats=feats, population="synthetic test population",
                          source_tables=["synthetic"], unavailable_blocks={"medications": "not in synthetic data"})


def test_build_dpv_synthetic():
    emb, blocks, pca = DP.build_dpv(_synthetic_spec())
    assert set(emb["dataset_id"]) == {"synthdp"}
    z = emb[emb["dimension"] == "z:x1"]["value"]
    assert len(z) == 54 and z.mean() == pytest.approx(0, abs=1e-9) and z.std(ddof=0) == pytest.approx(1)
    miss = emb[emb["dimension"] == "missing:x1"]
    assert len(miss) == 60 and miss["value"].sum() == 6
    assert not (emb["dimension"] == "missing:x2").any()  # no indicator without missingness
    b = blocks.set_index("feature")
    assert not b.loc["mostly_missing", "included"] and "missing" in b.loc["mostly_missing", "exclusion_reason"]
    assert not b.loc["constant", "included"] and "zero variance" in b.loc["constant", "exclusion_reason"]
    assert not b.loc["not_there", "included"]
    assert b.loc["(block not available)", "block"] == "medications"
    assert set(b.loc["c1", "dimensions"].split(";")) == {"c1=p", "c1=q", "c1=r"}
    # label-defining features are dimensions but excluded from the PCA fit
    sc = pca[(pca["block"] == "symptoms_conditions") & (pca["component"] == 1)].iloc[0]
    assert sc["n_features"] == 1 and sc["n_label_defining_excluded_from_fit"] == 1
    assert (emb["dimension"] == "z:lab1").any()
    assert {"pca:labs:pc1", "pca:wearable:pc1", "available:wearable"} <= set(emb["dimension"])
    assert set(emb["dimension_type"]) == {"zscore", "missing_indicator", "block_available", "pca"}
    # the embeddings themselves carry the label-defining flag (not only the blocks table)
    assert emb.loc[emb["dimension"] == "z:lab1", "label_defining"].all()
    assert not emb.loc[emb["dimension"] != "z:lab1", "label_defining"].any()


def test_build_dpv_non_spec_block_gets_availability_and_pca():
    spec = _synthetic_spec()
    rng = np.random.default_rng(3)
    for c in ("x3", "x4"):
        spec.frame[c] = rng.normal(size=len(spec.frame))
        spec.feats.append(DP.Feat(c, "omics_availability", "numeric", "synthetic"))
    emb, blocks, pca = DP.build_dpv(spec)
    assert (emb["dimension"] == "available:omics_availability").sum() == len(spec.frame)
    assert (pca["block"] == "omics_availability").any()


def test_harmonize_participant_keys():
    df = pd.DataFrame({"participant_id": ["nhanes:1", "stanford_longcovid_uwakwe2025:0", "mapmecfs:MECFS_101", "x:2"],
                       "dataset_id": ["nhanes", "stanford_longcovid_uwakwe2025", "mapmecfs", "x"],
                       "sex": ["male", "female", "Female", UNKNOWN], "age_years": [85.0, np.nan, np.nan, np.nan],
                       "age_range": [np.nan, "60-69", np.nan, UNKNOWN], "age_years_somalogic": [np.nan, np.nan, 44.0, np.nan]})
    h = DP.harmonize_participant_keys(df)
    assert h["sex_harmonized"].tolist() == ["male", "female", "female", UNKNOWN]
    assert h["age_band_harmonized"].tolist() == ["80+", "60-69", "40-49", UNKNOWN]
    assert np.isnan(h["age_years_numeric"].iloc[1]) and h["age_years_numeric"].iloc[2] == 44.0


# ------------------------------------------------------------------------------------------ data
UNIONS_BUILT = pytest.mark.skipif(not table_exists("participants"), reason="run digital_person first")


@pytest.mark.data
@UNIONS_BUILT
def test_unions_keep_dataset_and_namespaced_ids():
    for t in DP.UNION_TABLES:
        df = read_table(t, columns=["participant_id", "dataset_id", "source_partition"])
        expected = sum(len(pd.read_parquet(p, columns=["participant_id"])) for p in partitions(t))
        assert len(df) == expected, t
        assert (df["participant_id"].str.split(":", n=1).str[0] == df["dataset_id"]).all(), t
    p = read_table("participants", columns=["participant_id", "dataset_id", "sex_harmonized"])
    assert not p["participant_id"].duplicated().any()
    assert {"nhanes", "mapmecfs", "stanford_longcovid_uwakwe2025"} <= set(p["dataset_id"])
    chk = pd.read_csv(config.TABLES / "digital_person_union_checks.csv")
    assert (chk["n_prefix_mismatch"] == 0).all() and (chk["rows_expected"] == chk["rows_written"]).all()
    assert (chk["n_participants_not_in_participants"] == 0).all()


EMB_BUILT = pytest.mark.skipif(not table_exists("participant_phenotype_embeddings"), reason="run digital_person first")


@pytest.mark.data
@EMB_BUILT
def test_embeddings_are_per_dataset_and_never_pooled():
    e = read_table("participant_phenotype_embeddings", columns=["participant_id", "dataset_id", "dimension",
                                                                "dimension_type", "value"])
    assert set(e["dataset_id"].astype(str)) == {"nhanes", "stanford_covid_mishra2020", "stanford_covid_alavi2022",
                                                "stanford_longcovid_uwakwe2025", "mapmecfs"}
    per = e.groupby("participant_id", observed=True)["dataset_id"].nunique()
    assert (per == 1).all()
    assert (e["participant_id"].str.split(":", n=1).str[0] == e["dataset_id"].astype(str)).all()
    # within-dataset z-scores are centred per dimension
    z = e[e["dimension_type"].astype(str) == "zscore"]
    means = z.groupby([z["dataset_id"].astype(str), "dimension"], observed=True)["value"].mean()
    assert means.abs().max() < 1e-6
    # outcome labels are never dimensions
    dims = " ".join(e["dimension"].astype(str).unique()).lower()
    for bad in ("long_covid_label", "cohort_group", "group_label", "mort_", "seqn", "wtmec", "reference_date",
                "rhr_acute", "days_to_rhr_recovery"):
        assert bad not in dims, bad


@pytest.mark.data
@EMB_BUILT
def test_dataset_labels_and_molecular_context_statements():
    d = read_table("digital_phenotype_datasets").set_index("dataset_id")
    assert d.loc["mapmecfs", "representation_label"] == MAPMECFS_LABEL
    assert d.loc["mapmecfs", "representation_type"] == "omics_only_participant_linked"
    for ds in ("nhanes", "stanford_covid_mishra2020", "stanford_covid_alavi2022", "stanford_longcovid_uwakwe2025"):
        s = d.loc[ds, "molecular_context_statement"]
        assert "condition-level molecular enrichment" in s and "NOT individual multi-omics" in s
        assert d.loc[ds, "molecular_context_type"] == "condition_level_molecular_enrichment_via_ontology"
    assert not d["pooled_across_datasets"].astype(bool).any()
    assert "wearable" not in d.loc["mapmecfs", "blocks_available"]


@pytest.mark.data
@EMB_BUILT
def test_blocks_table_flags_label_defining_and_unavailable_blocks():
    b = read_table("digital_phenotype_blocks")
    uw = b[b["dataset_id"] == "stanford_longcovid_uwakwe2025"]
    chronic = uw[uw["feature"].str.startswith("symptom[chronic]")]
    assert len(chronic) > 0 and chronic["label_defining"].astype(bool).all()
    assert not uw[uw["feature"].str.startswith("symptom[acute]")]["label_defining"].astype(bool).any()
    mm = b[(b["dataset_id"] == "mapmecfs") & (b["block"] == "wearable")]
    assert len(mm) == 1 and not mm["included"].astype(bool).iloc[0]
    nh = b[(b["dataset_id"] == "nhanes") & b["included"].astype(bool)]
    assert set(DP.BLOCKS) <= set(nh["block"])


@pytest.mark.data
@EMB_BUILT
def test_query_functions_return_unknown_instead_of_guessing():
    assert DP.get_digital_phenotype("nhanes:0")["status"] == UNKNOWN
    assert DP.describe_dataset_representation("not_a_dataset")["status"] == UNKNOWN
    assert DP.molecular_context_route("not_a_dataset:1")["status"] == UNKNOWN
    rep = DP.describe_dataset_representation("mapmecfs")
    assert rep["status"] == "ok" and "wearable" in rep["blocks_unavailable"]
    mm = DP.molecular_context_route("mapmecfs:MECFS_101")
    assert mm["molecular_context_type"] == "participant_linked_omics_omics_only"
    assert mm["wearable_or_physiology_linkable"] is False
    adult = read_table("participant_clinical_features__nhanes", columns=["participant_id", "age_years"])
    pid = adult.loc[adult["age_years"] >= 18, "participant_id"].iloc[0]
    g = DP.get_digital_phenotype(pid)
    assert g["status"] == "ok" and g["dataset_id"] == "nhanes" and "demographics" in g["blocks"]
    child = adult.loc[adult["age_years"] < 18, "participant_id"].iloc[0]
    assert DP.get_digital_phenotype(child)["status"] == UNKNOWN
    # a block with z dimensions is reported available (omics_availability was missed before review)
    mm = DP.get_digital_phenotype("mapmecfs:MECFS_101")
    assert mm["blocks"]["omics_availability"]["available"] is True


@pytest.mark.data
@EMB_BUILT
def test_label_defining_flag_on_embeddings():
    e = read_table("participant_phenotype_embeddings", columns=["dataset_id", "dimension", "label_defining"])
    uw = e[e["dataset_id"].astype(str) == "stanford_longcovid_uwakwe2025"]
    chronic = uw[uw["dimension"].astype(str).str.startswith("z:symptom[chronic]")]
    assert len(chronic) and chronic["label_defining"].all()
    assert not e.loc[e["dimension"].astype(str).str.startswith("pca:"), "label_defining"].any()
    assert not e.loc[e["dataset_id"].astype(str) == "nhanes", "label_defining"].any()


@pytest.mark.data
def test_molecular_route_uses_exact_matches_only():
    """NHANES 'Headache' prescription reasons (ICD-10-CM R51) must not be routed to migraine."""
    c = read_table("participant_conditions__nhanes", columns=["participant_id", "icd10cm_code", "condition_label"])
    codes = c.groupby("participant_id")["icd10cm_code"].apply(lambda s: set(s.dropna().astype(str)))
    labels = c.groupby("participant_id")["condition_label"].apply(lambda s: set(s.astype(str)))
    headache_only = [p for p in codes.index if codes[p] & {"R51", "R51.P"} and not any(
        x.startswith("G43") for x in codes[p]) and not any("igraine" in x for x in labels[p])]
    assert headache_only
    r = DP.molecular_context_route(headache_only[0])
    routed = [x["canonical_condition_id"] for x in r["conditions_via_ontology"]] \
        if isinstance(r["conditions_via_ontology"], list) else []
    assert "migraine" not in routed
    mig = [p for p in codes.index if "G43" in codes[p]]
    rm = DP.molecular_context_route(mig[0])
    assert "migraine" in [x["canonical_condition_id"] for x in rm["conditions_via_ontology"]]
