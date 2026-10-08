"""Tests for measure_it.harmonize: person_concepts contract, curated configs, feature matrices, computable-phenotype
scoring. Pure tests run on a fresh clone; tests marked `data` need the harmonized NHANES partitions
(`uv run python -m measure_it.harmonize.run`)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from measure_it.config import PROCESSED, load_config
from measure_it.harmonize import schema as S
from measure_it.harmonize.features import feature_catalog, feature_matrix, icd10_category
from measure_it.harmonize.phenotype import phenotype_coverage, score_phenotype, validate_phenotype
from measure_it.harmonize.self_report import crosswalk_for, instrument_code, map_history

HAS_NHANES_CONCEPTS = (PROCESSED / "person_concepts__nhanes.parquet").exists()


def loinc_check_digit_ok(code: str) -> bool:
    """LOINC Mod 10 check digit (the algorithm published by Regenstrief)."""
    body, check = code.split("-")
    digits = body.zfill(7)[::-1]
    odd = "".join(d for i, d in enumerate(digits) if i % 2 == 0)
    even = "".join(d for i, d in enumerate(digits) if i % 2 == 1)
    total = sum(int(c) for c in str(int(odd[::-1]) * 2)) + sum(int(c) for c in even)
    return (10 - total % 10) % 10 == int(check)


# ------------------------------------------------------------------------------------------------ pure helpers

@pytest.mark.parametrize("code,cat", [("I73.00", "I73"), ("G9332", "G93"), ("u09.9", "U09"), ("E11", "E11"),
                                      ("", None), (None, None)])
def test_icd10_category(code, cat):
    assert icd10_category(code) == cat


def test_age_bands():
    a = pd.Series([3, 44, 79.9, 80, 85, 95, np.nan])
    assert S.age_band_from_years(a).tolist() == ["0-9", "40-49", "70-79", "80-89", "80-89", "90+", None]
    assert S.age_band_from_years(a, topcode=80).tolist()[3:6] == ["80+", "80+", "80+"]
    assert S.age_band_from_supplied(pd.Series(["40-44", "45-54", "90+", "85+", None])).tolist() == [
        "40-49", "45-54", "90+", "80+", None]
    assert S.age_mid(pd.Series(["40-49", "90+"])).tolist() == [44.5, 92.5]


def _concepts(rows):
    df = pd.DataFrame(rows)
    df["dataset_id"] = df["participant_id"].str.split(":").str[0]
    for c, v in (("mapping_method", "native"), ("mapping_confidence", "high"), ("source_variable", "x")):
        if c not in df.columns:
            df[c] = v
    return S.coerce(df)


def toy_concepts() -> pd.DataFrame:
    rows = []
    for ds, ids in (("a", ["a:1", "a:2", "a:3"]), ("b", ["b:1", "b:2"])):
        for i, p in enumerate(ids):
            rows.append({"participant_id": p, "domain": "demographic", "vocabulary": "DEMOG", "concept_code": "sex",
                         "value_as_string": "female" if i % 2 == 0 else "male"})
            rows.append({"participant_id": p, "domain": "demographic", "vocabulary": "DEMOG",
                         "concept_code": "age_band", "value_as_string": "40-49"})
    # dataset a: EHR-like codes (presence only); dataset b: survey-like (explicit 0 answers)
    rows += [{"participant_id": "a:1", "domain": "condition", "vocabulary": "ICD10CM", "concept_code": "I73.00",
              "value_as_number": 1.0},
             {"participant_id": "a:1", "domain": "condition", "vocabulary": "ICD10CM", "concept_code": "I73.9",
              "value_as_number": 1.0},
             {"participant_id": "b:1", "domain": "condition", "vocabulary": "ICD10CM", "concept_code": "I10",
              "value_as_number": 1.0},
             {"participant_id": "b:2", "domain": "condition", "vocabulary": "ICD10CM", "concept_code": "I10",
              "value_as_number": 0.0},
             {"participant_id": "a:2", "domain": "condition", "vocabulary": "ICD10CM", "concept_code": "I10",
              "value_as_number": 1.0},
             {"participant_id": "a:1", "domain": "measurement", "vocabulary": "LOINC", "concept_code": "2093-3",
              "value_as_number": 180.0, "unit": "mg/dL"},
             {"participant_id": "a:1", "domain": "measurement", "vocabulary": "LOINC", "concept_code": "2093-3",
              "value_as_number": 200.0, "unit": "mg/dL"},
             {"participant_id": "a:2", "domain": "drug", "vocabulary": "RXNORM", "concept_code": "6809",
              "value_as_number": 1.0},
             {"participant_id": "b:1", "domain": "survey", "vocabulary": "SURVEY", "concept_code": "PHQ9:total",
              "value_as_number": 12.0}]
    return _concepts(rows)


def test_concepts_contract_checks():
    c = toy_concepts()
    assert S.check_concepts(c) == []
    bad = c.copy()
    bad.loc[0, "vocabulary"] = "LOINC"            # demographic domain with a LOINC vocabulary
    assert any("domain/vocabulary" in e for e in S.check_concepts(bad))
    bad = c.copy()
    bad.loc[bad["vocabulary"] == "ICD10CM", "mapping_method"] = "self_report_map"
    assert any("capped at medium" in e for e in S.check_concepts(bad))
    bad = c.copy()
    bad.loc[bad["concept_code"] == "2093-3", "unit"] = ["mg/dL", "mmol/L"]
    assert any("more than one unit" in e for e in S.check_concepts(bad))
    bad = c.assign(zip_code="02139")
    assert any("geographic" in e for e in S.check_concepts(bad))


def test_feature_matrix_semantics():
    fm = feature_matrix(None, concepts=toy_concepts())
    assert list(fm.index) == ["a:1", "a:2", "a:3", "b:1", "b:2"]
    assert fm.loc["a:1", "LOINC:2093-3"] == 190.0 and np.isnan(fm.loc["a:2", "LOINC:2093-3"])
    assert fm.loc["a:1", "ICD10CM:I73"] == 1.0                    # I73.00 and I73.9 -> one category
    # dataset a records I73 as presence only: everyone else in a is 0; dataset b never records I73: NaN
    assert fm.loc["a:3", "ICD10CM:I73"] == 0.0 and np.isnan(fm.loc["b:1", "ICD10CM:I73"])
    # I10 in dataset b has explicit 0 answers: a b participant without a row would be NaN; a's are 0-filled
    assert fm.loc["b:2", "ICD10CM:I10"] == 0.0 and fm.loc["a:3", "ICD10CM:I10"] == 0.0
    assert fm.loc["a:2", "RXNORM:6809"] == 1.0 and fm.loc["a:1", "RXNORM:6809"] == 0.0
    assert fm.loc["a:1", "DEMOG:female"] == 1.0 and fm.loc["a:2", "DEMOG:female"] == 0.0
    assert fm.loc["a:1", "DEMOG:age_mid"] == 44.5
    by_code = feature_matrix(None, concepts=toy_concepts(), condition_level="code")
    assert "ICD10CM:I73.00" in by_code.columns and "ICD10CM:I73" not in by_code.columns
    sel = feature_matrix(None, ["ICD10CM:I73", "SURVEY:PHQ9:total", "LOINC:9999-9"], concepts=toy_concepts())
    assert list(sel.columns) == ["ICD10CM:I73", "SURVEY:PHQ9:total", "LOINC:9999-9"]
    assert sel["LOINC:9999-9"].isna().all()
    nan = feature_matrix(None, concepts=toy_concepts(), absent="nan")
    assert np.isnan(nan.loc["a:3", "ICD10CM:I73"])
    dom = feature_matrix(None, concepts=toy_concepts(), domains=["survey"])
    assert list(dom.columns) == ["SURVEY:PHQ9:total"] and len(dom) == 5


def test_feature_catalog_ehr_evaluable():
    cat = feature_catalog(toy_concepts()).set_index("feature_key")
    assert cat.loc["ICD10CM:I73", "codes"] == ["I73.00", "I73.9"]
    assert cat.loc["ICD10CM:I73", "ehr_evaluable"] and cat.loc["LOINC:2093-3", "ehr_evaluable"]
    assert not cat.loc["SURVEY:PHQ9:total", "ehr_evaluable"]
    assert cat.loc["LOINC:2093-3", "unit"] == "mg/dL"


# ------------------------------------------------------------------------------------------------ curated configs

def test_icd10_configs_are_well_formed_and_capped():
    nh = load_config("harmonize_nhanes_conditions")
    sr = load_config("harmonize_self_report_icd10")
    for r in nh["mappings"] + sr["mappings"]:
        assert S.ICD10CM_RE.fullmatch(r["icd10cm"]), r
        assert r["confidence"] in ("medium", "low"), r           # self-report: never high
        assert r["rationale"]
    for u in nh["unmapped"]:
        assert u["reason"]


@pytest.mark.skipif(not (PROCESSED / "ontology_icd10cm_codes.parquet").exists(), reason="needs ontology tables")
@pytest.mark.data
def test_icd10_config_codes_exist_in_cms_list():
    from measure_it.store import read_table
    known = set(read_table("ontology_icd10cm_codes", columns=["code"])["code"])
    for name in ("harmonize_nhanes_conditions", "harmonize_self_report_icd10"):
        for r in load_config(name)["mappings"]:
            assert r["icd10cm"] in known, (name, r["icd10cm"])


def test_loinc_configs_check_digits_and_units():
    nh = load_config("harmonize_nhanes_loinc")
    ng = load_config("harmonize_nightingale_loinc")
    for r in nh["mappings"] + ng["mappings"]:
        assert loinc_check_digit_ok(r["loinc"]), r["loinc"]
        assert r["unit"] and r["confidence"] in ("high", "medium", "low") and r["loinc_name"]
    for r in ng["mappings"]:
        assert r["confidence"] in ("medium", "low") and r["factor"] > 0 and r["source_unit"]
    # one unit per LOINC code (contract): rows sharing a code agree on the unit
    units = {}
    for r in nh["mappings"]:
        assert units.setdefault(r["loinc"], r["unit"]) == r["unit"], r["loinc"]
    # a Nightingale row that meets an NHANES LOINC must convert into the NHANES unit
    for r in ng["mappings"]:
        if r["loinc"] in units:
            assert units[r["loinc"]] == r["unit"], r["measure"]


def test_self_report_mapping():
    r = map_history("Raynaud's phenomenon")
    assert r == {"status": "mapped", "icd10cm": "I73.0", "label": "Raynaud's syndrome", "confidence": "medium",
                 "how": "curated alias"}
    assert map_history("  POTS ")["icd10cm"] == "G90.A"
    assert map_history("chronic lyme")["confidence"] == "low"
    assert map_history("seasonal allergies")["status"] == "unmapped"          # never fuzzy-matched
    assert map_history("Raynauds-like fingers")["status"] == "unmapped"
    sup = map_history("whatever the participant wrote", "M797")
    assert sup["icd10cm"] == "M79.7" and sup["confidence"] == "medium"


def test_survey_instruments_and_crosswalk():
    assert instrument_code("compass31") == "COMPASS31" and instrument_code("fss") == "FSS9"
    assert instrument_code("my_new_scale") == "MY_NEW_SCALE"
    cw = crosswalk_for("FSS9:total")
    assert cw and cw[0]["reference"] == "PHQ9:item4_tired" and cw[0]["match_quality"] == "related_construct"
    assert crosswalk_for("COMPASS31:vasomotor")[0]["match_quality"] == "none"
    same = crosswalk_for("PHQ9:total")
    assert same[0]["reference"] == "PHQ9:total" and same[0]["match_quality"] == "same_item"
    assert crosswalk_for("UNKNOWN:item") == []


def test_rxnorm_components_and_mocked_lookup(monkeypatch):
    from measure_it.harmonize import rxnorm as R
    assert R.split_components("HYDROCHLOROTHIAZIDE; LISINOPRIL") == ["HYDROCHLOROTHIAZIDE", "LISINOPRIL"]
    calls = {"rxcui.json": {"idGroup": {"rxnormId": ["221124"]}},
             "rxcui/221124/properties.json": {"properties": {"tty": "PIN", "name": "metoprolol succinate"}},
             "rxcui/221124/related.json": {"relatedGroup": {"conceptGroup": [
                 {"tty": "IN", "conceptProperties": [{"rxcui": "6918", "name": "metoprolol", "tty": "IN"}]}]}}}
    monkeypatch.setattr(R, "_get", lambda path, params=None: calls[path])
    R.ingredient_for_name.cache_clear()
    r = R.ingredient_for_name("METOPROLOL SUCCINATE")
    assert r["status"] == "mapped" and r["rxcui"] == "6918" and r["tty"] == "IN" and "via" in r
    assert R.ingredient_for_name("99999")["status"] == "unmapped"
    monkeypatch.setattr(R, "_get", lambda path, params=None: {"idGroup": {}})
    R.ingredient_for_name.cache_clear()
    assert R.ingredient_for_name("NOT A DRUG")["status"] == "unmapped"
    R.ingredient_for_name.cache_clear()


# ------------------------------------------------------------------------------------------------ phenotype scoring

def toy_phenotype() -> dict:
    return {
        "phenotype_id": "byod_toy_set:sub", "version": 1, "condition_id": "long_covid",
        "base_population": {"description": "d", "icd10cm": ["U09.9"], "label_column": "lc"},
        "subgroup_definition": "device-positive = score >= 2", "model": "l1_logistic", "intercept": -1.0,
        "features": [
            {"feature_key": "ICD10CM:I73", "domain": "condition", "vocabulary": "ICD10CM", "codes": ["I73.0"],
             "transform": "presence", "center": 0.3, "scale": 0.45, "coefficient": 1.0,
             "coefficient_standardized": 0.45, "label": "Raynaud", "ehr_evaluable": True},
            {"feature_key": "SURVEY:COMPASS31:vasomotor", "domain": "survey", "vocabulary": "SURVEY", "codes": [],
             "transform": "standardized_value", "center": 2.0, "scale": 1.0, "coefficient": 0.5,
             "coefficient_standardized": 0.5, "label": "vasomotor", "ehr_evaluable": False}],
        "decision_threshold": 0.5,
        "performance": {"auroc": 0.7, "auroc_ci_low": 0.6, "auroc_ci_high": 0.8, "n_positive": 30, "n_negative": 40,
                        "comparator": "cases without the finding"},
        "strata_fractions": [{"age_band": "all", "sex": "all", "n_cases": 70, "n_positive": 30, "fraction": 0.43,
                              "ci_low": 0.32, "ci_high": 0.54, "suppressed": False}],
        "plan_sha256": "0" * 64, "caveats": []}


def test_score_phenotype_formula_and_coverage():
    ph = toy_phenotype()
    fm = pd.DataFrame({"ICD10CM:I73": [1.0, 0.0, np.nan], "SURVEY:COMPASS31:vasomotor": [3.0, 2.0, 4.0]},
                      index=["p1", "p2", "p3"])
    p = score_phenotype(ph, fm)
    lp = np.array([-1 + 1.0 + 0.5, -1 + 0.0 + 0.0, -1 + 0.3 + 1.0])     # p3: I73 missing -> center 0.3
    assert np.allclose(p.to_numpy(), 1 / (1 + np.exp(-lp)))
    assert p.attrs["coverage"]["ICD10CM:I73"] == pytest.approx(2 / 3)
    assert p.attrs["n_rows_with_any_missing"] == 1 and p.attrs["missing_features"] == []
    strict = score_phenotype(ph, fm, missing="nan")
    assert strict.isna().tolist() == [False, False, True]
    # a reference cohort that never asked the survey item: reported, never silently dropped
    ref = pd.DataFrame({"ICD10CM:I73": [1.0, 0.0]}, index=["r1", "r2"])
    q = score_phenotype(ph, ref)
    assert q.attrs["missing_features"] == ["SURVEY:COMPASS31:vasomotor"]
    assert q.attrs["coefficient_mass_available"] == pytest.approx(0.45 / 0.95)
    cov = phenotype_coverage(ph, ref)
    assert cov["n_features_observed"].tolist() == [1, 1]
    assert np.allclose(cov["coefficient_mass_observed"], 0.45 / 0.95)
    assert np.allclose(score_phenotype(ph, fm, output="logit").to_numpy(), lp)


def test_phenotype_schema():
    ph = toy_phenotype()
    assert validate_phenotype(ph) == []
    bad = json.loads(json.dumps(ph))
    bad["features"][0]["vocabulary"] = "DEVICE"           # device features are never part of a phenotype
    bad["features"][0]["feature_key"] = "DEVICE:capillaroscopy:score"
    assert validate_phenotype(bad)
    bad = json.loads(json.dumps(ph))
    del bad["features"][1]["center"]                       # standardized features need center and scale
    assert validate_phenotype(bad)
    bad = json.loads(json.dumps(ph))
    bad["strata_fractions"] = [{"age_band": "40-49", "sex": "female", "n_cases": 12, "n_positive": 5,
                                "fraction": 0.4, "ci_low": 0.2, "ci_high": 0.6, "suppressed": False}]
    assert validate_phenotype(bad)                         # the pooled 'all' row is required
    del ph["features"][0]["ehr_evaluable"]
    assert validate_phenotype(ph)


# ------------------------------------------------------------------------------------------------ NHANES (data)

@pytest.mark.data
@pytest.mark.skipif(not HAS_NHANES_CONCEPTS, reason="needs person_concepts__nhanes (measure_it.harmonize.run)")
@pytest.mark.parametrize("ds", ["nhanes", "nhanes0306"])
def test_nhanes_person_concepts_contract(ds):
    from measure_it.provenance import check_provenance
    from measure_it.store import read_table
    t = read_table(f"person_concepts__{ds}")
    check_provenance(t, f"person_concepts__{ds}")
    assert S.check_concepts(t) == []
    assert set(t["data_layer"]) == {"person"} and set(t["source_geographic_resolution"]) == {"none"}
    P = read_table(f"participants__{ds}", columns=["participant_id"])
    sex = t[(t["vocabulary"] == "DEMOG") & (t["concept_code"] == "sex")]
    assert sex["participant_id"].is_unique and set(sex["participant_id"]) == set(P["participant_id"])
    assert set(t["domain"]) >= {"demographic", "condition", "measurement", "drug", "survey"}
    assert t.loc[t["concept_code"] == "2093-3", "unit"].unique().tolist() == ["mg/dL"]
    assert t.loc[t["vocabulary"] == "DEMOG", "value_as_string"].isin(
        ["female", "male", "other", "unknown", "80+"] + [f"{a}-{a + 9}" for a in range(0, 80, 10)]).all()
    fm = feature_matrix([ds], ["ICD10CM:I10", "LOINC:2093-3", "DEMOG:age_mid", "RXNORM:6809"])
    assert fm["ICD10CM:I10"].dropna().isin([0.0, 1.0]).all() and fm["LOINC:2093-3"].median() > 150
    assert fm["RXNORM:6809"].mean() > 0.01                  # metformin users exist in both cycles
