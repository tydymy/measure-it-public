"""Tests for measure_it.mcp.trace (trace_evidence).

Pure tests run anywhere; tests marked `data` need data/processed (the pipeline outputs).
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from measure_it.config import PROCESSED, PROJECT_ROOT, TABLES, UNKNOWN
from measure_it.mcp import trace as TR


# ------------------------------------------------------------------------------------------------ pure

def test_parse_object_id():
    assert TR.parse_object_id("condition:long_covid") == ("condition", "long_covid")
    assert TR.parse_object_id("fda:code:DQA") == ("fda", "code:DQA")
    assert TR.parse_object_id("opportunity:a|b|06073") == ("opportunity", "a|b|06073")
    assert TR.parse_object_id("nonsense") == (None, None)
    assert TR.parse_object_id("") == (None, None)


def test_every_conventions_namespace_has_a_resolver():
    expected = {"condition", "measurement", "trial", "nih", "npi", "hrsa_site", "facility", "fda", "gene", "variant",
                "gwas_study", "geo_series", "geo", "opportunity", "signature", "measurement_evidence",
                "measurement_signal", "phenotype_evidence", "molbio", "reactome", "measurement_performance"}
    assert set(TR.RESOLVERS) == expected == set(TR.NAMESPACES)


@pytest.mark.parametrize("oid", ["participant:nhanes:62161", "participant:stanford_covid:ABC", "person:1",
                                 "nhanes:62161", "stanford_covid_mishra2020:A06L7KF", "mapmecfs:MECFS_101",
                                 "Participant:nhanes:62161", " SEQN:62161", "nhanes0306:21005",
                                 "klein2023_mylc_ml_table:LC.MS.0001", "charlton_lc_mecfs_cpet_source:P25",
                                 "fm_thermography:FM001", "endo_arg1_repod:Patient-1",
                                 "heds_hsd_olink_serum_cinquina2026:C54"])
def test_person_level_ids_are_refused(oid):
    r = TR.trace_evidence(oid)
    assert r["status"] == UNKNOWN and r["refused"] is True
    assert "summarised only" in r["reason"]
    assert "rows" not in r


@pytest.mark.parametrize("oid", ["nonsense", "foo:bar", "", ":x"])
def test_malformed_or_unknown_namespace_is_unknown(oid):
    r = TR.trace_evidence(oid)
    assert r["status"] == UNKNOWN and r["reason"]
    assert "known_namespaces" in r


def test_to_jsonable_nan_inf_numpy():
    import numpy as np
    x = {"a": float("nan"), "b": float("inf"), "c": np.int64(3), "d": np.array([1.0, np.nan]), "e": pd.NA,
         "f": {1, 2}, "g": pd.Timestamp("2026-09-24"), "h": "x" * 50}
    out = TR.to_jsonable(x, max_text=10)
    json.dumps(out, allow_nan=False)
    assert out["a"] is None and out["b"] is None and out["c"] == 3 and out["d"] == [1.0, None] and out["e"] is None
    assert out["f"] == [1, 2] and out["g"].startswith("2026-09-24") and out["h"].startswith("x" * 10 + " ...[cut")


# ------------------------------------------------------------------------------------------------ data

def _first(table: str, col: str, where=None):
    df = pd.read_parquet(PROCESSED / f"{table}.parquet", columns=[col] + ([where[0]] if where else []))
    if where:
        df = df[df[where[0]] == where[1]]
    return str(df[col].dropna().iloc[0])


def sample_ids() -> list[str]:
    ids = [
        f"condition:{_first('condition_registry', 'canonical_condition_id')}",
        "condition:long_covid_or_me_cfs",
        _first("measurement_registry", "object_id"),
        f"trial:{_first('clinical_trials', 'nct_id')}",
        f"nih:{_first('nih_projects', 'appl_id')}",
        f"npi:{_first('providers', 'npi')}",
        f"hrsa_site:{_first('facilities__hrsa', 'site_bphc_number')}",
        _first("facilities", "object_id"),
        f"fda:{_first('fda_510k', 'clearance_id')}",
        f"fda:code:{_first('fda_device_classification', 'product_code')}",
        f"fda:{_first('fda_pma', 'pma_number')}",
        f"gene:{_first('condition_molecular_evidence', 'entity_id', ('entity_type', 'gene'))}",
        f"variant:{_first('condition_molecular_evidence', 'entity_id', ('entity_type', 'variant'))}",
        f"gwas_study:{_first('condition_molecular_evidence', 'entity_id', ('source_database', 'GWAS Catalog'))}",
        f"geo_series:{_first('geo_study_catalog', 'gse')}",
        f"geo:{_first('geographies', 'geo_id', ('geo_level', 'county'))}",
        f"geo:{_first('geographies', 'geo_id', ('geo_level', 'state'))}",
        "geo:US",
        _first("deployment_opportunities", "object_id"),
        _first("phenotype_signatures", "object_id"),
        _first("condition_measurement_evidence", "object_id"),
        _first("measurement_phenotype_signal", "object_id"),
        _first("phase3_measurement_evidence", "object_id"),
        _first("measurable_biology", "object_id"),
        _first("measurement_performance", "object_id"),
        f"reactome:{pd.read_csv(TABLES / 'test3_reactome_ora.csv', usecols=['reactome_id'])['reactome_id'].iloc[0]}",
    ]
    return ids


@pytest.mark.data
def test_every_namespace_resolves_a_real_id():
    ids = sample_ids()
    assert {TR.parse_object_id(i)[0] for i in ids} == set(TR.NAMESPACES)
    for oid in ids:
        r = TR.trace_evidence(oid)
        assert r["status"] == "ok", (oid, r.get("reason"))
        assert r["object_id"] == oid
        assert r["n_rows"] >= 1 and r["rows"], oid
        json.dumps(r, allow_nan=False)
        assert r["sources"], f"{oid}: no source"
        for s in r["sources"]:
            if s.get("status") == UNKNOWN:
                continue
            assert s["data_audit"], oid
            assert any((PROJECT_ROOT / a["path"]).exists() for a in s["data_audit"]), (oid, s["source_id"])
            assert "raw" in s and "manifest" in s["raw"]


@pytest.mark.data
@pytest.mark.parametrize("oid", ["condition:zzz_not_a_condition", "measurement:no_such_class", "trial:NCT00000000",
                                 "nih:1", "npi:0000000000", "hrsa_site:NOPE", "facility:nope", "fda:K000000",
                                 "fda:code:ZZZ", "gene:ENSG00000000000", "variant:rs0", "gwas_study:GCST0",
                                 "geo_series:GSE0", "geo:99999", "geo:abc", "opportunity:a|b|c", "opportunity:bad",
                                 "signature:x|y|z", "measurement_evidence:x|y", "measurement_signal:x",
                                 "phenotype_evidence:x|y", "molbio:x|y|z", "reactome:R-HSA-0"])
def test_unknown_ids_return_unknown_with_reason(oid):
    r = TR.trace_evidence(oid)
    assert r["status"] == UNKNOWN, oid
    assert r["reason"]


@pytest.mark.data
def test_opportunity_lineage_reaches_features_burden_rows_and_raw_files():
    oid = _first("deployment_opportunities", "object_id", ("condition_id", "long_covid_or_me_cfs"))
    r = TR.trace_evidence(oid)
    assert r["status"] == "ok"
    comps = r["components"]
    assert {"burden", "vulnerability", "diagnostic_desert", "clinic_capacity", "research_readiness"} <= set(comps)
    assert comps["burden"]["evidence_level"] in ("A", "B", "C", "D")
    kinds = {lk.get("table") for lk in r["lineage"]}
    assert "geo_condition_features" in kinds and "geo_condition_burden" in kinds
    linked = {lk.get("object_id") for lk in r["lineage"] if lk.get("object_id")}
    geo_id = oid.split("|")[-1]
    assert f"geo:{geo_id}" in linked and "condition:long_covid" in linked and "condition:me_cfs" in linked
    burden = [lk for lk in r["lineage"] if lk.get("table") == "geo_condition_burden"]
    for lk in burden:
        assert lk["source_ids"], lk
        assert lk["row"]["burden_evidence_level"] in ("A", "B", "C")
        if lk["row"]["condition_id"] == "long_covid" and lk["row"]["geo_level"] == "county":
            # inherited HPS state value (hps_inherited) or the labelled BRFSS small-area model (brfss_sae, default)
            assert ((lk["row"]["inherited"] is True and lk["row"]["source_geographic_resolution"] == "state")
                    or lk["row"]["derivation"] == "mrp_small_area_estimate")
    src = {s["source_id"]: s for s in r["sources"]}
    lc_src = "cdc_brfss" if "cdc_brfss" in src else "cdc_long_covid"
    assert lc_src in src and "cms_mmd" in src and "cdc_svi" in src
    files = [f["path"] for f in src[lc_src]["raw"]["record_files"]]
    assert files and any(f.endswith((".csv", ".zip")) for f in files)
    assert all((PROJECT_ROOT / f).exists() for f in files)


@pytest.mark.data
def test_opportunity_trace_explains_incomplete_burden_and_null_rank():
    """An incomplete-burden row (Greater Bridgeport CT: no CMS ME/CFS value for the planning regions) surfaces
    burden_incomplete, its reason, the missing members and why composite/rank are null; a small complete county says
    it is scored but not ranked; a ranked row has no note."""
    oid = "opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|09120"
    r = TR.trace_evidence(oid)
    assert r["status"] == "ok"
    assert r["burden_incomplete"] is True and r["burden_members_missing"] == ["me_cfs"]
    assert "me_cfs" in r["burden_incomplete_reason"]
    b = r["components"]["burden"]
    assert b["incomplete"] is True and b["members_missing"] == ["me_cfs"] and b["incomplete_reason"]
    assert "Connecticut" in b["missing_source_reason"]
    assert r["composite_equal"] is None and r["rank_equal"] is None
    assert "incomplete" in r["composite_rank_note"] and r["composite_rank_note"] in r["notes"]
    # Kusilvak AK (02158): complete since the CMS FIPS bridge, population < 10,000 -> scored, not ranked
    s = TR.trace_evidence("opportunity:long_covid_or_me_cfs|wearable_autonomic_activity_monitoring|02158")
    assert s["burden_incomplete"] is False and s["burden_members_missing"] == []
    assert s["composite_equal"] is not None and s["rank_equal"] is None
    assert "not ranked" in s["composite_rank_note"] and "rank_equal_incl_small" in s["composite_rank_note"]
    d = pd.read_parquet(PROCESSED / "deployment_opportunities.parquet", columns=["object_id", "rank_equal"])
    ok = TR.trace_evidence(str(d.loc[d["rank_equal"] == 1, "object_id"].iloc[0]))
    assert ok["burden_incomplete"] is False and "composite_rank_note" not in ok


@pytest.mark.data
def test_signature_lineage_names_dataset_and_label_basis():
    oid = _first("phenotype_signatures", "object_id", ("phenotype_id", "mecfs_like_proxy"))
    r = TR.trace_evidence(oid)
    assert r["status"] == "ok"
    lb = r["label_basis"]
    assert lb["is_proxy"] is True and lb["label_basis"] and lb["phenotype_definition"]
    assert any(lk.get("relation") == "dataset" for lk in r["lineage"])
    assert [s["source_id"] for s in r["sources"]] == ["nhanes_2011_2014"]
    # aggregate row only: no participant identifier anywhere in the trace
    assert "participant_id" not in json.dumps(r)


@pytest.mark.data
def test_trial_trace_finds_raw_record_batch():
    oid = f"trial:{_first('clinical_trials', 'nct_id')}"
    r = TR.trace_evidence(oid)
    ct = next(s for s in r["sources"] if s["source_id"] == "clinicaltrials_gov")
    files = ct["raw"]["record_files"]
    assert files and files[0]["path"].startswith("data/raw/clinicaltrials_gov/records/batch_")
    assert oid.split(":")[1] in (PROJECT_ROOT / files[0]["path"]).read_text()


@pytest.mark.data
def test_facility_trace_lists_source_records():
    oid = _first("facilities", "object_id", ("primary_kind", "hrsa_health_center_site"))
    r = TR.trace_evidence(oid)
    assert r["status"] == "ok"
    ids = [lk.get("object_id") for lk in r["lineage"] if lk.get("object_id")]
    assert any(i.startswith("hrsa_site:") for i in ids)
    assert "hrsa_health_centers" in [s["source_id"] for s in r["sources"]]


@pytest.mark.data
def test_every_signature_dataset_traces_to_its_source():
    """One signature per dataset_id (NHANES 2003-2006, GEO cohorts, the lab and device releases, ...) resolves, and
    its dataset link names a registered source with a DATA_AUDIT."""
    df = pd.read_parquet(PROCESSED / "phenotype_signatures.parquet", columns=["dataset_id", "object_id"])
    firsts = df.dropna().groupby("dataset_id")["object_id"].first()
    assert len(firsts) >= 10
    for ds, oid in firsts.items():
        r = TR.trace_evidence(oid)
        assert r["status"] == "ok", (oid, r.get("reason"))
        link = next(lk for lk in r["lineage"] if lk.get("relation") == "dataset")
        assert link["source_ids"], (ds, oid)
        assert all((PROJECT_ROOT / "data/raw" / s / "DATA_AUDIT.md").exists() for s in link["source_ids"]), ds


def test_dataset_source_mapping():
    assert TR._dataset_source("nhanes0306") == ["nhanes_2003_2006"]
    assert TR._dataset_source("nhanes_labs") == ["nhanes_2011_2014"]
    assert TR._dataset_source("geo_GSE270045") == ["geo_cohorts"]
    assert TR._dataset_source("no_such_dataset") == []
