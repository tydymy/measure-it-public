"""Tests for stage F: configs/measurement_hcpcs.yaml, configs/condition_rx_signals.yaml, the CMS ingestion modules and
measure_it.facilities.activity."""
from __future__ import annotations

import re

import numpy as np
import pytest

from measure_it.config import load_config
from measure_it.facilities import activity as A
from measure_it.ingestion import cms_partd_prescriber as PD
from measure_it.ingestion import cms_physician_service as PS
from measure_it.provenance import check_provenance
from measure_it.store import read_table

CODE_RE = re.compile(r"^(?:\d{5}|\d{4}[FTU]|[A-V]\d{4})$")


# --------------------------------------------------------------------------- configs (pure)

def test_hcpcs_config_covers_every_measurement_id():
    cfg = load_config("measurement_hcpcs")["measurements"]
    classes = {c["id"] for c in load_config("measurements")["measurement_classes"]}
    assert classes <= set(cfg), f"classes without an entry: {classes - set(cfg)}"
    bundles = set(load_config("relevance")["measurement_bundles"])
    assert bundles <= set(load_config("measurement_hcpcs").get("bundles", {}))


def test_hcpcs_codes_well_formed_and_sourced():
    cfg = load_config("measurement_hcpcs")["measurements"]
    for mid, spec in cfg.items():
        for key in ("codes", "related_codes"):
            for c in spec.get(key) or []:
                assert CODE_RE.match(str(c["code"])), (mid, c["code"])
                assert c.get("rationale") and c.get("confidence") in {"high", "medium", "low"}, (mid, c["code"])
                assert c.get("source", "").startswith("http"), (mid, c["code"])
        assert spec.get("summary"), mid


def test_capillaroscopy_has_no_dedicated_code_only_labelled_proxies():
    spec = load_config("measurement_hcpcs")["measurements"]["capillaroscopy"]
    assert not spec.get("codes")
    assert spec.get("related_codes") and all(c.get("proxy") is True for c in spec["related_codes"])
    assert "no" in spec["summary"].lower() and "cpt" in spec["summary"].lower()


def test_code_map_roles():
    cm = PS.code_map({"measurements": {
        "a": {"codes": [{"code": "11111"}], "related_codes": [{"code": "22222"}, {"code": "33333", "proxy": True}]},
        "b": {"codes": [{"code": "22222"}], "related_codes": [{"code": "22222"}]}}})
    r = cm.set_index(["measurement_class", "hcpcs_cd"])["code_role"]
    assert r[("a", "11111")] == "dedicated" and r[("a", "22222")] == "related" and r[("a", "33333")] == "proxy"
    assert r[("b", "22222")] == "dedicated" and len(cm) == 4  # duplicate keeps the strongest role


def test_rx_config_has_no_lyme_or_long_covid_direct_signal():
    cfg = load_config("condition_rx_signals")
    rows = PD.drug_rows(cfg)
    rel = {}
    for s in cfg["signals"].values():
        for c, v in s["relation"].items():
            rel.setdefault(c, set()).add(v)
    assert "lyme_disease" not in rel and "ptlds" not in rel
    assert rel.get("long_covid", set()) <= {"related_phenotype"} and rel.get("me_cfs", set()) <= {"related_phenotype"}
    assert set(rows["specificity"]) <= {"low", "medium", "high"}


def test_capacity_basis_order():
    roles = {"capillaroscopy": {"proxy"}, "microvascular_function": set(), "autonomic_testing": {"dedicated"},
             "hrv": {"related"}}
    assert A.capacity_basis("capillaroscopy", roles) == "proxy"
    assert A.capacity_basis("autonomic_function_testing", roles) == "dedicated"
    assert A.capacity_basis("hrv", roles) == "related"
    assert A.capacity_basis("metabolomics", roles) == "none"


def test_measurement_members_are_generic():
    mm = A.measurement_members()
    assert mm["nailfold_capillaroscopy"]["members"] == ["capillaroscopy", "microvascular_function"]
    assert mm["cpet"] == {"kind": "class", "members": ["cpet"]}
    with pytest.raises(KeyError):
        A.members_of("not_a_measurement")


def test_percentile_rank_matches_scoring():
    from measure_it.scoring.opportunity import percentile_rank
    x = np.array([3.0, 1.0, np.nan, 1.0, 7.0])
    np.testing.assert_allclose(A.percentile_rank(x), percentile_rank(x), equal_nan=True)


# --------------------------------------------------------------------------- data tests

@pytest.mark.data
def test_activity_table_contract_and_suppression():
    t = read_table("provider_measurement_activity")
    check_provenance(t, "provider_measurement_activity")
    assert {"npi", "measurement_class", "hcpcs_cd", "tot_services", "data_year", "payer_scope", "code_role"} <= set(t)
    assert (t["payer_scope"] == "Medicare FFS Part B").all()
    assert t["tot_benes_max_pos"].min() >= 11  # CMS suppression rule holds in what we ingested
    classes = {c["id"] for c in load_config("measurements")["measurement_classes"]}
    assert set(t["measurement_class"]) <= classes
    assert (t.loc[t["measurement_class"] == "capillaroscopy", "code_role"] == "proxy").all()


@pytest.mark.data
def test_rx_table_contract():
    t = read_table("provider_condition_rx_signals")
    check_provenance(t, "provider_condition_rx_signals")
    assert t["tot_claims"].min() >= 11
    assert set(t["drug_generic"]) <= set(PD.drug_rows(PD.load_signals())["gnrc_name"])


@pytest.mark.data
def test_geo_measurement_capacity():
    cap = read_table("geo_measurement_capacity")
    check_provenance(cap, "geo_measurement_capacity")
    ids = set(A.measurement_members())
    assert set(cap["measurement_id"]) == ids
    c = cap[cap["geo_level"] == "county"]
    assert c.groupby("measurement_id").size().nunique() == 1  # every measurement on the same county universe
    assert (c["active_clinicians_within_50km_n"].fillna(0) >= c["active_clinicians_n"].fillna(0)).mean() > 0.99
    assert (cap.loc[cap["measurement_id"] == "capillaroscopy", "capacity_basis"] == "proxy").all()
    assert (cap.loc[cap["measurement_id"] == "nailfold_capillaroscopy", "capacity_basis"] != "dedicated").all()
    none = cap[cap["capacity_basis"] == "none"]
    assert none["measurement_capacity"].isna().all() and {"posture_detection", "metabolomics"} <= set(none["measurement_id"])
    # counts reconcile with the provider summary (county = NPPES / CMS-ZIP county of individual NPIs)
    s = read_table("provider_measurement_summary")
    s = s[(s["measurement_id"] == "autonomic_function_testing") & (s["entity_code"] == "I")
          & (s["services_dedicated"] > 0)]
    k = c[c["measurement_id"] == "autonomic_function_testing"].set_index("geo_id")["active_clinicians_n"]
    exp = s.groupby("county_fips")["npi"].nunique()
    common = exp.index.intersection(k.index)
    assert (k.loc[common].astype(int) == exp.loc[common]).all()


@pytest.mark.data
def test_experience_sites_us_and_generic():
    e = read_table("measurement_experience_sites")
    check_provenance(e, "measurement_experience_sites")
    assert e["state_fips"].notna().all()
    assert ((e["n_trials"] > 0) | (e["n_nih_projects"] > 0)).all()
    for r in e.sample(min(200, len(e)), random_state=0).itertuples():
        assert set(r.member_classes_matched.split("|")) <= set(A.members_of(r.measurement_id))
