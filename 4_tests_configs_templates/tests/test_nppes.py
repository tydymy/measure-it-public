"""Tests for NPPES + NUCC ingestion (measure_it.ingestion.nppes) and specialty helpers.

Unit tests need no data. @pytest.mark.data tests read data/processed outputs and the
raw NUCC CSV; run `uv run python -m measure_it.ingestion.nppes` first.
"""
from __future__ import annotations

import json
from datetime import date

import pandas as pd
import pytest

from measure_it.config import CONFIGS, PROCESSED, RAW
from measure_it.facilities import specialties as sp
from measure_it.ingestion import nppes
from measure_it.provenance import PROVENANCE_COLUMNS

# ------------------------------------------------------------------ unit: specialties helpers

TOY = {
    "groups": {
        **{g: {"group_type": "specialty", "spec_listed": True,
               "codes": [{"code": "207RC0000X", "display_name": "x", "role": "physician"}]}
           for g in sp.SPEC_SPECIALTY_GROUPS},
        "electrophysiology": {"group_type": "specialty", "spec_listed": True,
                              "codes": [{"code": "207RC0001X", "display_name": "EP", "role": "physician"}]},
        "cardiology": {"group_type": "specialty", "spec_listed": True,
                       "codes": [{"code": "207RC0000X", "display_name": "CV", "role": "physician"},
                                 {"code": "207RC0001X", "display_name": "EP", "role": "physician"}]},
        "fqhc": {"group_type": "facility", "spec_listed": False,
                 "codes": [{"code": "261QF0400X", "display_name": "FQHC", "role": "facility"}]},
    }
}


def test_normalize_taxonomy_code():
    assert sp.normalize_taxonomy_code(" 207rc0000x ") == "207RC0000X"
    assert sp.normalize_taxonomy_code("207RC0000") is None       # 9 chars
    assert sp.normalize_taxonomy_code("207RC0000Y") is None      # must end in X
    assert sp.normalize_taxonomy_code("") is None
    assert sp.normalize_taxonomy_code(None) is None
    assert sp.normalize_taxonomy_code(float("nan")) is None


def test_code_to_groups_overlap_and_union():
    m = sp.code_to_groups(TOY)
    assert m["207RC0001X"] == ("cardiology", "electrophysiology")
    assert "fqhc" in m["261QF0400X"]
    assert sp.groups_for_codes(["207RC0001X", "261QF0400X", None, "junk"], m) == \
        ["cardiology", "electrophysiology", "fqhc"]
    assert sp.groups_for_codes(["101Y00000X"], m) == []
    assert sp.join_groups(["b", "a", "b", ""]) == "a|b"


def test_group_code_frame_shape():
    f = sp.group_code_frame(TOY)
    assert set(f.columns) >= {"specialty_group", "group_type", "taxonomy_code", "role", "spec_listed"}
    assert len(f[f.taxonomy_code == "207RC0001X"]) == 2


def test_select_primary_taxonomy():
    assert sp.select_primary_taxonomy(["A00000000X", "B00000000X"], ["N", "Y"]) == ("B00000000X", "primary_switch_Y")
    assert sp.select_primary_taxonomy(["A00000000X"], [None]) == ("A00000000X", "single_code")
    assert sp.select_primary_taxonomy(["A00000000X", "A00000000X"], ["N", "N"]) == ("A00000000X", "single_code")
    assert sp.select_primary_taxonomy(["A00000000X", "B00000000X"], ["N", "N"]) == \
        ("A00000000X", "first_listed_no_primary_flag")
    assert sp.select_primary_taxonomy(["A00000000X", "B00000000X"], ["Y", "Y"]) == ("A00000000X", "first_Y_of_multiple")
    assert sp.select_primary_taxonomy([None, ""], ["Y", "N"]) == (None, "no_taxonomy")


def test_validate_group_config_catches_errors():
    assert sp.validate_group_config(TOY) == []
    bad = json.loads(json.dumps(TOY))
    del bad["groups"]["vascular"]
    bad["groups"]["fqhc"]["codes"].append({"code": "261QF0400X", "display_name": "dup", "role": "facility"})
    bad["groups"]["cardiology"]["codes"].append({"code": "BADCODE", "display_name": "x", "role": "physician"})
    bad["groups"]["rheumatology"]["codes"][0]["role"] = "wizard"
    bad["groups"]["pmr"]["excluded_considered"] = [{"code": "207RC0000X", "reason": "x"}]
    probs = " ".join(sp.validate_group_config(bad))
    for frag in ["'vascular' missing", "duplicate code 261QF0400X", "malformed code 'BADCODE'", "bad role",
                 "both included and excluded"]:
        assert frag in probs


def test_role_from_nucc():
    assert sp.role_from_nucc("Allopathic & Osteopathic Physicians", "Internal Medicine", "Individual") == "physician"
    assert sp.role_from_nucc("Physician Assistants & Advanced Practice Nursing Providers", "Nurse Practitioner",
                             "Individual") == "nurse_practitioner"
    assert sp.role_from_nucc("Hospitals", "General Acute Care Hospital", "Non-Individual") == "facility"
    assert sp.role_from_nucc("Other Service Providers", "Specialist", "Individual") == "other_individual"


def test_real_config_structure():
    cfg = sp.load_specialty_groups(CONFIGS / "specialty_groups.yaml")
    g = cfg["groups"]
    # exactly the 13 SPEC groups are flagged spec_listed
    assert {k for k, v in g.items() if v["spec_listed"]} == set(sp.SPEC_SPECIALTY_GROUPS)
    codes = {k: {c["code"] for c in v["codes"]} for k, v in g.items()}
    assert codes["electrophysiology"] == {"207RC0001X"}                  # cardiac EP only
    assert {"2084N0600X", "204R00000X"} <= codes["neurology"]             # neurophysiology -> neurology
    assert "207RC0001X" in codes["cardiology"]
    assert "363A00000X" not in codes["primary_care"]                      # PAs excluded
    assert "363LF0000X" in codes["primary_care"]                          # family NPs included
    assert "208M00000X" not in codes["primary_care"]                      # hospitalists excluded
    assert codes["fqhc"] == {"261QF0400X"}
    assert codes["rural_health_clinic"] == {"261QR1300X"}
    assert "282N00000X" in codes["general_acute_care_hospital"]
    for v in g.values():   # facility groups contain only non-individual codes except the research group
        if v["group_type"] == "facility":
            assert {c["role"] for c in v["codes"]} == {"facility"}


# ------------------------------------------------------------------ unit: nppes pure helpers

@pytest.mark.parametrize("raw,expected", [
    ("CA", "CA"), ("ca", "CA"), ("CALIFORNIA", "CA"), ("California", "CA"), ("CA - CALIFORNIA", "CA"),
    ("MD-MARYLAND", "MD"), ("PUERTO RICO (PR)", "PR"), ("P.R.", "PR"), ("FL.", "FL"), ("US VIRGIN ISLANDS", "VI"),
    ("AE", "AE"), ("FM", "FM"), ("CLARK COUNTY", None), ("N/A", None), ("10025", None), (None, None),
    ("225400000X", None), ("PUERTO RRICO", None),
])
def test_normalize_state(raw, expected):
    assert nppes.normalize_state(raw) == expected


@pytest.mark.parametrize("country,state,expected", [
    ("US", "NY", ("NY", "us_state_dc")),
    ("US", "DC", ("DC", "us_state_dc")),
    (None, "TX", ("TX", "us_state_dc")),
    ("US", "PR", ("PR", "us_territory")),
    ("US", "GU", ("GU", "us_territory")),
    ("US", "AE", ("AE", "dropped_military_apo_fpo")),
    ("US", "FM", ("FM", "dropped_freely_associated_state")),
    ("US", "N/A", (None, "dropped_unresolved_state")),
    ("CA", "ON", (None, "dropped_non_us_country")),
    ("MX", "CA", (None, "dropped_non_us_country")),
    ("UM", "PUERTO RICO", ("PR", "us_territory")),
    ("UM", "CALIFORNIA", ("CA", "us_state_dc")),
    ("UM", "ISRAEL", (None, "dropped_unresolved_state")),
])
def test_classify_us_location(country, state, expected):
    assert nppes.classify_us_location(country, state) == expected


def test_zip5_and_names_and_deactivation():
    assert nppes.zip5("021151234") == "02115"
    assert nppes.zip5("96350 1821") == "96350"
    assert nppes.zip5("AE") is None
    assert nppes.zip5(None) is None
    assert nppes.zip5("1234") is None
    assert nppes.display_name(1, "JANE", "DOE", None) == "JANE DOE"
    assert nppes.display_name(2, None, None, "ACME CLINIC") == "ACME CLINIC"
    assert nppes.display_name("2", "X", "Y", "") is None
    assert nppes.is_deactivated(date(2020, 1, 1), None) is True
    assert nppes.is_deactivated(date(2020, 1, 1), date(2020, 2, 1)) is False
    assert nppes.is_deactivated(date(2020, 3, 1), date(2020, 2, 1)) is True
    assert nppes.is_deactivated(None, date(2020, 2, 1)) is False


def test_membership_resolution():
    gm = pd.Series(["zip5_as_zcta_internal_point", "not_geocoded"], index=[5, 9])
    out = nppes.membership_resolution(gm)
    assert out.tolist() == ["zcta_centroid", "none"] and out.index.tolist() == [5, 9]


def test_top_ungeocoded_zips():
    ng = pd.DataFrame({"zip5": ["44195", "44195", "94143"], "city": ["Cleveland", "CLEVELAND", "SAN FRANCISCO"],
                       "state": ["OH", "OH", "CA"]})
    assert nppes._top_ungeocoded_zips(ng) == "44195: 2 (CLEVELAND OH); 94143: 1 (SAN FRANCISCO CA)"


def test_connecticut_summary():
    prov = pd.DataFrame({
        "state": ["CT", "CT", "CT", "CT", "NY"],
        "county_fips": ["09110", "09003", "36027", None, "09190"],
    })
    s = nppes.connecticut_summary(prov)
    assert s == {"ct_address_providers": 4, "ct_geocoded_to_2024_planning_region": 1,
                 "ct_geocoded_to_legacy_county": 1, "ct_geocoded_to_other_state_county": 1,
                 "ct_not_geocoded": 1, "any_provider_on_legacy_ct_county": True}


# ------------------------------------------------------------------ data tests

def _need(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{p} missing; run `uv run python -m measure_it.ingestion.nppes`")
    return pd.read_parquet(p)


@pytest.fixture(scope="module")
def providers():
    return _need("providers")


@pytest.fixture(scope="module")
def membership():
    return _need("provider_specialty_groups")


@pytest.fixture(scope="module")
def density():
    return _need("provider_density_county")


@pytest.fixture(scope="module")
def counts():
    return _need("provider_taxonomy_counts")


@pytest.fixture(scope="module")
def stats():
    p = RAW / "nppes" / "ingest_stats.json"
    if not p.exists():
        pytest.skip("ingest_stats.json missing")
    return json.loads(p.read_text())


@pytest.mark.data
def test_config_matches_nucc_file():
    files = sorted((RAW / "nucc_taxonomy").glob("nucc_taxonomy_*.csv"))
    if not files:
        pytest.skip("NUCC CSV not downloaded")
    nucc = pd.read_csv(files[-1], dtype=str).fillna("")
    cfg = sp.load_specialty_groups()
    assert sp.validate_against_nucc(cfg, nucc) == []
    # roles in the YAML agree with NUCC grouping
    idx = nucc.set_index("Code")
    for gid, g in cfg["groups"].items():
        for c in g["codes"]:
            r = idx.loc[c["code"]]
            assert sp.role_from_nucc(r["Grouping"], r["Classification"], r["Section"]) == c["role"], (gid, c)


@pytest.mark.data
def test_providers_basic_integrity(providers):
    assert providers["npi"].is_unique
    assert providers["npi"].str.fullmatch(r"\d{10}").all()
    assert set(providers["entity_type"].unique()) <= {1, 2}
    assert (providers["specialty_groups"].str.len() > 0).all()
    assert providers["primary_taxonomy_code"].notna().all()
    assert providers["us_location_category"].isin(["us_state_dc", "us_territory"]).all()
    assert providers["state"].isin(nppes.US_STATES_DC | nppes.US_TERRITORIES).all()
    z = providers["zip5"].dropna()
    assert z.str.fullmatch(r"\d{5}").all()
    # organisations never carry an individual credential
    assert providers.loc[providers["entity_type"] == 2, "credential"].isna().all()


@pytest.mark.data
def test_providers_no_deactivated(providers, stats):
    d = providers[providers["deactivation_date"].notna()]
    assert (d["reactivation_date"].notna() & (d["reactivation_date"] >= d["deactivation_date"])).all()
    assert stats["deactivated_excluded"] > 0
    x = stats["deactivation_crosscheck"]
    # every NPI in the CMS deactivation report is flagged deactivated in the full file, or is newer than it
    assert x["active_in_full_file"] + x["not_in_full_file"] <= 5


@pytest.mark.data
def test_primary_groups_subset_of_any(providers):
    for anyg, prim in zip(providers["specialty_groups"], providers["primary_specialty_groups"]):
        if prim:
            assert set(prim.split("|")) <= set(anyg.split("|"))


@pytest.mark.data
def test_group_membership_consistent_with_codes(providers):
    cfg = sp.load_specialty_groups()
    m = sp.code_to_groups(cfg)
    sample = providers.sample(n=min(5000, len(providers)), random_state=20260923)
    for codes, groups, prim, pgroups in zip(sample["all_taxonomy_codes"], sample["specialty_groups"],
                                            sample["primary_taxonomy_code"], sample["primary_specialty_groups"]):
        assert sp.groups_for_codes(codes.split("|"), m) == groups.split("|")
        assert sp.join_groups(m.get(prim, ())) == (pgroups or "")


@pytest.mark.data
def test_geocoding_fields(providers):
    g = providers[providers["geocode_method"] == "zip5_as_zcta_internal_point"]
    assert len(g) / len(providers) > 0.9
    assert g["lat"].between(-15, 72).all() and g["lon"].between(-180, 180).all()
    assert g["county_fips"].str.fullmatch(r"\d{5}").all()   # incl. relationship-file fallback for water points
    assert set(g["county_assignment_method"]) <= {"zcta_point_in_2024_county", "zcta_largest_land_share_rel2020"}
    # the fallback is rare and never assigns legacy Connecticut counties
    fb = g[g["county_assignment_method"] == "zcta_largest_land_share_rel2020"]
    assert len(fb) < 0.01 * len(g) and not fb["county_fips"].str.startswith("09").any()
    assert (g["geo_state_matches_address_state"] == True).mean() > 0.99  # noqa: E712
    assert (g["source_geographic_resolution"] == "zcta_centroid").all()
    ng = providers[providers["geocode_method"] == "not_geocoded"]
    assert ng["lat"].isna().all() and ng["county_fips"].isna().all()
    assert (ng["source_geographic_resolution"] == "state").all()


@pytest.mark.data
@pytest.mark.parametrize("name", ["providers", "provider_specialty_groups", "provider_density_county",
                                  "provider_taxonomy_counts", "nucc_taxonomy"])
def test_provenance_and_caveat(name):
    df = _need(name)
    assert all(c in df.columns for c in PROVENANCE_COLUMNS)
    assert (df["data_layer"] == "facility").all()
    if name != "nucc_taxonomy":
        assert df["provenance_notes"].str.contains("does NOT mean the provider evaluates or treats Long COVID").all()
        assert df["provenance_notes"].str.contains("billing").all()
        assert df["provenance_notes"].str.contains("ZIP centroids").all()


@pytest.mark.data
def test_membership_matches_providers(providers, membership):
    assert set(membership["npi"]) == set(providers["npi"])
    per = membership.groupby("npi")["specialty_group"].apply(lambda s: "|".join(sorted(s)))
    assert (per.reindex(providers["npi"]).values == providers["specialty_groups"].values).all()
    assert not membership.duplicated(["npi", "specialty_group"]).any()
    # county_fips is ZIP-centroid derived: provenance must say so row by row
    has_cty = membership["county_fips"].notna()
    assert (membership.loc[has_cty, "source_geographic_resolution"] == "zcta_centroid").all()
    assert (membership.loc[~has_cty, "source_geographic_resolution"] == "none").all()


@pytest.mark.data
def test_connecticut_vintage(providers, stats):
    legacy = {f"090{i:02d}" for i in range(1, 16, 2)}
    assert not providers["county_fips"].isin(legacy).any()
    ct = providers[(providers["state"] == "CT") & providers["county_fips"].str.startswith("09", na=False)]
    assert len(ct) > 0 and ct["county_fips"].str.fullmatch(r"091[1-9]0").all()
    assert stats["connecticut"]["ct_geocoded_to_legacy_county"] == 0
    assert "Connecticut vintage" in (RAW / "nppes" / "DATA_AUDIT.md").read_text()


@pytest.mark.data
def test_density_counts(density, membership):
    cols = ["n_individual_providers", "n_individual_providers_primary", "n_individual_physicians",
            "n_organizations", "n_organizations_primary"]
    assert (density[cols] >= 0).all().all()
    assert (density["n_individual_providers_primary"] <= density["n_individual_providers"]).all()
    assert (density["n_individual_physicians"] <= density["n_individual_providers"]).all()
    assert (density["n_organizations_primary"] <= density["n_organizations"]).all()
    assert not density.duplicated(["county_fips", "specialty_group"]).any()
    # zeros are kept: every county has a row for every group
    assert density.groupby("specialty_group")["county_fips"].nunique().nunique() == 1
    # totals equal the geocoded membership rows (counties outside the 2024 list would break this)
    geo = membership[membership["county_fips"].notna()]
    ind = geo[geo["entity_type"] == 1].groupby("specialty_group").size()
    tot = density.groupby("specialty_group")["n_individual_providers"].sum()
    assert (tot.reindex(ind.index) == ind).all()


@pytest.mark.data
def test_taxonomy_counts(counts, providers):
    grp = counts[counts["count_level"] == "specialty_group"].set_index("specialty_group")
    for g in sp.SPEC_SPECIALTY_GROUPS:
        assert grp.loc[g, "n_npis_any"] > 0, g
        n = providers["specialty_groups"].str.split("|").apply(lambda s: g in s).sum()
        assert grp.loc[g, "n_npis_any"] == n
    assert (grp["n_npis_primary"] <= grp["n_npis_any"]).all()
    assert (grp["n_individual_any"] + grp["n_organizations_any"] == grp["n_npis_any"]).all()
    codes = counts[counts["count_level"] == "taxonomy_code"]
    assert set(codes["config_status"]) == {"included", "excluded_considered"}
    assert codes["taxonomy_code"].is_unique and codes["taxonomy_code"].is_monotonic_increasing  # deterministic order
    # sanity: the big specialties are big, and cardiac EP is a subset of cardiology
    assert grp.loc["primary_care", "n_npis_any"] > grp.loc["cardiology", "n_npis_any"] > \
        grp.loc["electrophysiology", "n_npis_any"]


@pytest.mark.data
def test_nucc_table():
    t = _need("nucc_taxonomy")
    assert t["taxonomy_code"].is_unique
    assert t["taxonomy_code"].str.fullmatch(r"[0-9A-Z]{9}X").all()
    assert (t.loc[t["taxonomy_code"] == "207RC0001X", "specialty_groups"].iloc[0]
            == "cardiology|electrophysiology")


@pytest.mark.data
def test_audit_and_registry_written():
    for sid in ("nppes", "nucc_taxonomy"):
        assert (RAW / sid / "DATA_AUDIT.md").exists()
        assert (RAW / sid / "registry_entry.yaml").exists()
        assert (RAW / sid / "MANIFEST.json").exists()
    txt = (RAW / "nppes" / "DATA_AUDIT.md").read_text()
    assert "does NOT mean the provider evaluates or treats Long COVID" in txt
    assert "Version 2" in txt
