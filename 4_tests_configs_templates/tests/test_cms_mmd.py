"""Tests for measure_it.ingestion.cms_mmd (CMS Mapping Medicare Disparities).

Pure-function tests run anywhere; @pytest.mark.data tests read the processed outputs
(run `uv run python -m measure_it.ingestion.cms_mmd` first).
"""
from __future__ import annotations

import json
import math
import shutil
from urllib.parse import parse_qs, urlparse

import pandas as pd
import pytest
import yaml

from measure_it.config import RAW, load_config
from measure_it.ingestion import cms_mmd as m
from measure_it.provenance import PROVENANCE_COLUMNS, check_provenance
from measure_it.registry import REQUIRED_FIELDS
from measure_it.store import processed_path, read_table

BURDEN = "geo_condition_burden__cms_mmd"
CONTEXT = "geo_context__cms_mmd"
COVERAGE = "condition_burden_coverage__cms_mmd"


# ------------------------------------------------------------------ pure functions
def test_normalize_fips_pads_and_labels():
    assert m.normalize_fips("4", "s") == "04"
    assert m.normalize_fips("04", "s") == "04"
    assert m.normalize_fips("1001", "c") == "01001"
    assert m.normalize_fips("56045", "c") == "56045"
    assert m.normalize_fips("", "n") == "US"
    with pytest.raises(ValueError):
        m.normalize_fips("123456", "c")
    with pytest.raises(ValueError):
        m.normalize_fips("AB", "s")


def test_fips_rename_bridge_is_one_to_one_and_county_only():
    g = pd.Series(["46113", "02270", "51515", "46", "01001", "02990"])
    lvl = pd.Series(["county", "county", "county", "state", "county", "county"])
    new, published = m.bridge_fips_renames(g, lvl)
    assert new.tolist() == ["46102", "02158", "51515", "46", "01001", "02990"]   # 51515 (a merger) is not bridged
    assert published.tolist() == g.tolist()
    assert set(m.FIPS_RENAMES) == {"46113", "02270"}
    assert "46102" in m.fips_bridge_note("46113") and m.fips_bridge_note("51515") is None


def test_normalize_stratum_treats_blank_dot_null_as_all():
    for v in ("", ".", None, float("nan"), " . "):
        assert m.normalize_stratum(v) == "all"
    assert m.normalize_stratum("0") == "0"


def test_decode_dencat():
    assert m.decode_dencat("3") == ("1,000-4,999", 1000.0, 4999.0)
    label, lo, hi = m.decode_dencat("5")
    assert label == "10,000+" and lo == 10000.0 and math.isnan(hi)
    label, lo, hi = m.decode_dencat(".")
    assert label is None and math.isnan(lo)


def test_resolve_source_file_replicates_app_first_match():
    xw = [
        {"population": "f", "measure": "v", "year": "2", "elig": ".", "race_code": "", "sex_code": ".",
         "adjust": "", "dual": ".", "url": "prev_year_2"},
        {"population": "f", "measure": "v", "year": "22", "elig": "", "race_code": ".", "sex_code": ".",
         "adjust": "12", "dual": "", "url": "prev_fltr12_22_f"},
        {"population": "f", "measure": "v", "year": "22", "elig": "", "race_code": ".", "sex_code": ".",
         "adjust": "34", "dual": "", "url": "prev_fltr34_22_f"},
    ]
    opts = {"population": "f", "measure": "v", "year": "22", "elig": ".", "race_code": ".", "sex_code": ".",
            "adjust": "1", "dual": "."}
    src, crit = m.resolve_source_file(xw, opts)
    assert src == "prev_fltr12_22_f" and "adjust" in crit
    src, _ = m.resolve_source_file(xw, {**opts, "adjust": "3"})
    assert src == "prev_fltr34_22_f"
    src, _ = m.resolve_source_file(xw, {**opts, "measure": "h"})
    assert src is None


def test_plan_requests_resolves_each_adjustment_separately():
    # Regression: principal cost 2012-2018 lives in one file per adjustment. Resolving once with
    # adjust=1 and asking that file for fltr 1|2 silently dropped every age-standardized row.
    xw = [
        {"population": "f", "measure": "p", "year": "2", "elig": ".", "race_code": "", "sex_code": ".",
         "adjust": "1", "dual": ".", "url": "pc_fltr_1"},
        {"population": "f", "measure": "p", "year": "2", "elig": ".", "race_code": "", "sex_code": "",
         "adjust": "2", "dual": ".", "url": "pc_fltr_2"},
        {"population": "f", "measure": "t", "year": "2", "elig": "", "race_code": "", "sex_code": "",
         "adjust": "", "dual": ".", "url": "tc_all"},
    ]
    plan = [p for p in m.plan_requests(xw) if p["year_code"] == "2" and p["measure"] in ("p", "t")]
    by = {(p["measure"], p["source_file"]): p for p in plan}
    assert by[("p", "pc_fltr_1")]["fltrs"] == ["1"] and by[("p", "pc_fltr_2")]["fltrs"] == ["2"]
    assert parse_qs(urlparse(by[("p", "pc_fltr_2")]["url"]).query)["fltr"] == ["2"]
    assert by[("t", "tc_all")]["fltrs"] == ["1", "2"]
    assert len({p["raw_file"] for p in plan}) == len(plan)


def test_icd10_text_search_and_registry_merge():
    txt = "Anemia 2 years C94.6, D46.0, D50.9 (any DX) ... Kidney M35.04, M35.0A, N01.0"
    codes = m.icd10_codes_in_text(txt)
    assert codes == ["C94.6", "D46.0", "D50.9", "M35.04", "M35.0A", "N01.0"]
    assert m.codes_matching(["M35.7"], codes) == []          # EDS hypermobility is not Sjogren M35.0x
    assert m.codes_matching(["D46"], codes) == ["D46.0"]
    # sourced registry codes already covered by a prefix are not duplicated; new ones are added
    assert m.merge_search_codes(["G93.3"], ["G93.32", "R53.82"]) == ["G93.3", "R53.82"]
    assert m.merge_search_codes(["K58"], None) == ["K58"]


def test_build_query_filters_and_2019_ed_exception():
    url, params = m.build_query("f1", year_code="23", measure="v", conditions=["51", "59"],
                                agecats=[".", "0"], fltrs=["1", "2"])
    q = parse_qs(urlparse(url).query)
    assert q["condition"] == ["51|59"]
    assert q["agecat"] == [".|0|IS NULL"]
    assert q["sexcat"] == [".|IS NULL"] and q["racecat"] == [".|IS NULL"] and q["dual"] == [".|IS NULL"]
    assert "eligcat" in q and "elig" not in q
    assert q["geography"] == ["c|s|n"]
    url, _ = m.build_query("f2", year_code="9", measure="e", conditions=["10"], agecats=["."], fltrs=["1"])
    q = parse_qs(urlparse(url).query)
    assert "elig" in q and "eligcat" not in q


def test_parse_ccw_algorithms_and_code_search():
    page = """<table><tr><th>Algorithms</th><th>Reference Period</th><th>ICD9</th><th>ICD10</th><th>Claims</th></tr>
    <tr><td>Fibromyalgia, Chronic Pain and Fatigue</td><td>2 Years</td><td>DX 729.1</td>
    <td>DX G89.29, M79.7, R53.82 (any DX on the claim)</td><td>At least 1 inpatient claim</td></tr>
    <tr><td>Chronic Pain 4</td><td>1 Year</td><td>DX 338</td><td>DX G43.001, M79.7</td><td>2 claims</td></tr></table>"""
    df = m.parse_ccw_algorithms(page)
    assert list(df["ccw_algorithm_name"]) == ["Fibromyalgia, Chronic Pain and Fatigue", "Chronic Pain"]
    assert df.iloc[0]["icd10_codes"] == ["G89.29", "M79.7", "R53.82"]
    assert m.codes_matching(["G93.3"], df.iloc[0]["icd10_codes"]) == []
    assert m.codes_matching(["R53.82"], df.iloc[0]["icd10_codes"]) == ["R53.82"]
    assert m.codes_matching(["G43"], df.iloc[1]["icd10_codes"]) == ["G43.001"]


def test_evidence_level_mapping_is_conservative():
    assert m.evidence_level_for("fibromyalgia") == "B"
    assert m.evidence_level_for("migraine") == "B"
    assert m.evidence_level_for("me_cfs") == "C"
    for cid in ("pots", "long_covid", "eds_hsd", "ibs", "endometriosis", "lyme_disease", "dysautonomia",
                "mcas", "gastroparesis", "ptlds"):
        assert m.evidence_level_for(cid) == "D"
    assert all(level != "A" for _, _, level in m.BURDEN_MAPPINGS)
    # every mapping target is a real configured condition
    ids = {c["id"] for c in load_config("conditions")["conditions"]}
    assert {cid for cid, _, _ in m.BURDEN_MAPPINGS} <= ids
    assert set(m.TARGET_SEARCH_CODES) <= ids | {"post_infectious_syndrome"}


def test_limitation_text_is_the_mandated_wording():
    assert "Medicare FFS beneficiaries (mostly >= 65 plus disability-entitled)" in m.LIMITATION
    assert "not representative of working-age invisible-illness populations" in m.LIMITATION
    assert "coded-condition prevalence depends on diagnosis/coding practice" in m.LIMITATION


# ------------------------------------------------------------------ processed outputs
def _need(name):
    if not processed_path(name).exists():
        pytest.skip(f"{name} not built; run uv run python -m measure_it.ingestion.cms_mmd")
    return read_table(name)


@pytest.fixture(scope="module")
def burden():
    return _need(BURDEN)


@pytest.fixture(scope="module")
def context():
    return _need(CONTEXT)


@pytest.fixture(scope="module")
def coverage():
    return _need(COVERAGE)


@pytest.mark.data
def test_burden_provenance_and_layer(burden):
    check_provenance(burden, BURDEN)
    assert set(PROVENANCE_COLUMNS) <= set(burden.columns)
    assert set(burden["data_layer"]) == {"geographic"}
    assert set(burden["evidence_type"]) == {"administrative_claims"}
    assert (burden["evidence_level"] == burden["burden_evidence_level"]).all()
    assert burden["provenance_notes"].str.contains(m.LIMITATION, regex=False).all()


@pytest.mark.data
def test_burden_levels_and_condition_mapping(burden):
    assert set(burden["burden_evidence_level"]) <= {"B", "C"}
    pairs = set(zip(burden["condition_id"], burden["source_condition_code"], burden["burden_evidence_level"]))
    assert pairs == {("fibromyalgia", "51", "B"), ("me_cfs", "51", "C"), ("migraine", "59", "B")}
    assert (burden.loc[burden["source_condition_code"] == "51", "ccw_algorithm_name"]
            == "Fibromyalgia, Chronic Pain and Fatigue").all()
    # the C proxy and the B mapping share identical numbers (same source measure, two labels)
    f = burden[burden["condition_id"] == "fibromyalgia"].set_index("source_record_id")["prevalence_pct"]
    c = burden[burden["condition_id"] == "me_cfs"].set_index("source_record_id")["prevalence_pct"]
    assert f.sort_index().equals(c.sort_index())


@pytest.mark.data
def test_burden_values_and_keys(burden):
    v = burden["prevalence_pct"].dropna()
    assert len(v) > 0 and v.between(0, 100).all()
    key = ["geo_id", "year", "condition_id", "adjustment", "stratum_age"]
    assert not burden.duplicated(key).any()
    assert set(burden["year"]) == set(range(2012, 2024))
    assert set(burden["adjustment"]) == {"unsmoothed_actual", "unsmoothed_age_standardized"}
    assert set(burden["stratum_age"]) == {"all", "<65"}


@pytest.mark.data
def test_geography_resolution_is_never_downscaled(burden, context):
    for df in (burden, context):
        assert (df["source_geographic_resolution"] == df["geo_level"]).all()
        cty, st = df[df["geo_level"] == "county"], df[df["geo_level"] == "state"]
        assert cty["geo_id"].str.fullmatch(r"\d{5}").all()
        assert st["geo_id"].str.fullmatch(r"\d{2}").all()
        assert st["county_fips"].isna().all()
        assert (df.loc[df["geo_level"] == "national", "geo_id"] == "US").all()


@pytest.mark.data
def test_connecticut_legacy_not_forced_onto_planning_regions(burden):
    ct = burden[(burden["geo_level"] == "county") & (burden["state_fips"] == "09")]
    assert len(ct) > 0
    legacy = ct["geo_id"].str.match(r"^090(0[1-9]|1[0-5])$")
    assert (ct.loc[legacy, "geo_match_status"] == "matched_ct_legacy").all()
    regions = ct["geo_id"].str.match(r"^091[1-9]0$")
    # a planning-region row may only exist if CMS itself published that FIPS
    assert (ct.loc[regions, "geo_match_status"] == "matched_2024").all()
    assert set(burden["geo_match_status"]) <= {"matched_2024", "matched_ct_legacy", "unmatched", "national"}


@pytest.mark.data
def test_unmatched_rows_are_explained_not_forced(burden, context):
    geos = read_table("geographies")
    known = set(geos["geo_id"])
    for df in (burden, context):
        um = df[df["geo_match_status"] == "unmatched"]
        assert um["unmatched_reason"].notna().all()
        assert not um["geo_id"].isin(known).any()
        assert df.loc[df["geo_match_status"] != "unmatched", "unmatched_reason"].isna().all()
        c990 = um["geo_id"].str.fullmatch(r"\d{2}990")
        assert um.loc[c990, "unmatched_reason"].str.startswith("cms_code_ending_990").all()
        assert (um.loc[~c990, "unmatched_reason"] == "fips_not_in_2024_census_vintage").all()
        matched = df[df["geo_match_status"].isin(["matched_2024", "matched_ct_legacy"])]
        assert matched["geo_id"].isin(known).all()


@pytest.mark.data
def test_fips_renames_carried_to_2024_codes(burden, context):
    """CMS publishes Oglala Lakota SD / Kusilvak AK under 46113 / 02270 (1:1 renames of 2015); ingestion bridges them."""
    for df in (burden, context):
        cty = df[df["geo_level"] == "county"]
        assert not cty["geo_id"].isin(["46113", "02270"]).any()
        br = cty[cty["fips_bridge"].notna()]
        assert set(zip(br["fips_as_published"], br["geo_id"])) == {("46113", "46102"), ("02270", "02158")}
        assert (br["geo_match_status"] == "matched_2024").all() and br["unmatched_reason"].isna().all()
        assert br["provenance_notes"].str.startswith("FIPS bridge:").all()
        assert br.apply(lambda r: f":{r['fips_as_published']}:" in r["source_record_id"], axis=1).all()
        rest = cty[cty["fips_bridge"].isna()]
        assert (rest["fips_as_published"] == rest["geo_id"]).all()
    b = burden[(burden["geo_id"] == "46102") & burden["default_view"] & (burden["condition_id"] == "me_cfs")]
    assert set(b["year"]) == set(range(2012, 2024)) and b["prevalence_pct"].notna().all()
    assert "51515" in set(burden.loc[burden["geo_match_status"] == "unmatched", "geo_id"])   # merger: not bridged


@pytest.mark.data
def test_county_match_rate_is_high(burden):
    c = burden[(burden["geo_level"] == "county") & burden["default_view"] & (burden["condition_id"] == "migraine")]
    for year, g in c.groupby("year"):
        overall = (g["geo_match_status"] != "unmatched").mean()
        assert overall > 0.98, (year, overall)
        census_style = g[~g["geo_id"].str.fullmatch(r"\d{2}990")]
        rate = (census_style["geo_match_status"] != "unmatched").mean()
        assert rate > 0.995, (year, rate)


@pytest.mark.data
def test_zero_values_are_flagged(burden):
    z = burden["prevalence_pct"] == 0
    assert (burden.loc[z, "value_zero_possibly_suppressed"]).all()
    assert not burden.loc[~z, "value_zero_possibly_suppressed"].any()


@pytest.mark.data
def test_context_is_not_burden(context):
    check_provenance(context, CONTEXT)
    assert "burden_evidence_level" not in context.columns
    assert context["evidence_level"].str.startswith("context").all()
    assert set(context["context_category"]) <= {"utilization", "spending", "population_composition",
                                                "infection_exposure"}
    covid = context[context["source_condition_code"] == "134"]
    assert len(covid) > 0
    assert covid["provenance_notes"].str.contains("NOT Long COVID").all()
    assert covid["related_condition_ids"].isna().all()
    assert not context["measure"].eq("prevalence").where(context["source_condition_code"].isin(["51", "59"]),
                                                          False).any()


@pytest.mark.data
def test_every_measure_year_has_both_adjustments(context, burden):
    # every (measure, year, MMD condition) the tool serves at county level has actual AND age-standardized rows
    for df, code_col in ((context, "measure_code"), (burden, "measure")):
        c = df[df["geo_level"] == "county"]
        adj = c.groupby([code_col, "year", "source_condition_code"])["adjustment"].agg(frozenset)
        bad = adj[adj != frozenset({"unsmoothed_actual", "unsmoothed_age_standardized"})]
        assert bad.empty, bad.head()
    pc = context[(context["measure"] == "average_principal_cost") & (context["year"] <= 2018)]
    assert (pc["adjustment"] == "unsmoothed_age_standardized").sum() > 5000


@pytest.mark.data
def test_raw_slices_match_their_manifest_query(context, burden):
    # the fltr list in each slice's MANIFEST url is exactly the set of adjustments parsed from it
    manifest = json.loads((RAW / m.SOURCE_ID / "MANIFEST.json").read_text())["files"]
    inv = {v: k for k, v in m.FLTR_CODES.items()}
    both = pd.concat([context[["raw_file", "adjustment"]], burden[["raw_file", "adjustment"]]])
    got = both.groupby("raw_file")["adjustment"].agg(lambda s: {inv[a] for a in s})
    for raw_file, fl in got.items():
        q = parse_qs(urlparse(manifest[raw_file]["url"]).query)
        assert set(q["fltr"][0].split("|")) == fl, raw_file


@pytest.mark.data
def test_coverage_table_lists_every_configured_condition(coverage, burden):
    check_provenance(coverage, COVERAGE)
    ids = {c["id"] for c in load_config("conditions")["conditions"]}
    assert set(coverage["condition_id"]) == ids
    lv = coverage.set_index("condition_id")["burden_evidence_level"]
    assert lv["fibromyalgia"] == "B" and lv["migraine"] == "B" and lv["me_cfs"] == "C"
    for cid in ("pots", "long_covid", "eds_hsd", "ibs", "endometriosis", "lyme_disease"):
        assert lv[cid] == "D"
    d = coverage[coverage["burden_evidence_level"] == "D"]
    assert (d["n_burden_rows"] == 0).all()
    assert not burden["condition_id"].isin(d["condition_id"]).any()
    fib = coverage.set_index("condition_id").loc["fibromyalgia"]
    assert "M79.7" in fib["ccw_icd10_codes"] and "R53.82" in fib["ccw_icd10_codes"]
    me = coverage.set_index("condition_id").loc["me_cfs"]
    assert "G93.3" not in me["ccw_icd10_codes"]
    assert coverage["rationale"].str.len().gt(40).all()
    # the 30 CCW chronic conditions (MMD-exposed, PDF-only) were searched, or the gap is stated
    if shutil.which("pdftotext"):
        assert (coverage["search_hits_in_ccw_30_chronic_conditions"] != "not searched").all()
        assert coverage["ccw_algorithm_sets_searched"].str.contains("CCW 30 chronic-condition algorithms (the set", regex=False).all()
        # no target code occurs there, so no D condition needs review
        assert not coverage["rationale"].str.contains("REVIEW").any()
    else:
        d_rat = coverage.loc[coverage["burden_evidence_level"] == "D", "rationale"]
        assert d_rat.str.contains("not searched|grouping", regex=True).all()


@pytest.mark.data
def test_search_codes_cover_the_sourced_condition_registry(coverage):
    if not processed_path("condition_registry").exists():
        pytest.skip("condition_registry not built")
    reg = read_table("condition_registry", columns=["canonical_condition_id", "icd10cm_codes"])
    searched = coverage.set_index("condition_id")["icd10_codes_searched"].str.split(", ")
    for r in reg.itertuples():
        if r.canonical_condition_id not in searched.index:
            continue
        keys = searched[r.canonical_condition_id]
        for code in r.icd10cm_codes:
            assert any(code == k or code.startswith(k) for k in keys), (r.canonical_condition_id, code)


@pytest.mark.data
def test_registry_entry_and_audit_exist():
    d = RAW / m.SOURCE_ID
    entry = yaml.safe_load((d / "registry_entry.yaml").read_text())
    assert not [f for f in REQUIRED_FIELDS if f not in entry]
    assert entry["person_level"] is False and entry["geographic"] is True
    assert any("not representative of working-age" in x for x in entry["limitations"])
    audit = (d / "DATA_AUDIT.md").read_text()
    for needle in ("source_id", "Missingness", "Linkage strategy", "Limitations", "burden_evidence_level",
                   "Connecticut"):
        assert needle in audit
    manifest = json.loads((d / "MANIFEST.json").read_text())
    api = [k for k in manifest["files"] if k.startswith("api/")]
    assert len(api) >= 60
    assert all(len(manifest["files"][k]["sha256"]) == 64 for k in api)
