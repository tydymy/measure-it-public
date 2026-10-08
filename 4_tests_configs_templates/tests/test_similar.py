"""Tests for stage D (measure_it.similar): the survey estimator against hand computations, coverage, phenotype
application on a synthetic reference with a known weighted prevalence, OMOP query generation (+ sqlglot parse when
installed: `uv run --with sqlglot pytest tests/test_similar.py`), the clinic export, and an end-to-end run against
a temporary processed/results tree. Tests marked `data` use the real NHANES partitions."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from measure_it.config import PROCESSED
from measure_it.similar import omop, phenotype_io, reference, survey


# --------------------------------------------------------------------------------------------- fixtures


def _phenotype(**over) -> dict:
    ph = {
        "phenotype_id": "byod_synth:capillary_positive", "version": 1, "condition_id": "long_covid",
        "base_population": {"description": "Long COVID (U09.9)", "icd10cm": ["U09.9"], "label_column": "lc"},
        "model": "l1_logistic", "intercept": -1.0,
        "features": [
            {"feature_key": "ICD10CM:I73", "domain": "condition", "vocabulary": "ICD10CM", "codes": ["I73"],
             "transform": "presence", "coefficient": 2.0, "label": "Raynaud / peripheral vascular"},
            {"feature_key": "LOINC:2093-3", "domain": "measurement", "vocabulary": "LOINC", "codes": ["2093-3"],
             "transform": "standardized_value", "center": 190.0, "scale": 40.0, "unit": "mg/dL",
             "coefficient": 0.5, "label": "Total cholesterol", "platform": "nmr_nightingale"},
            {"feature_key": "RXNORM:6809", "domain": "drug", "vocabulary": "RXNORM", "codes": ["6809"],
             "transform": "presence", "coefficient": -0.5, "label": "metformin"},
            {"feature_key": "DEMOG:female", "domain": "demographic", "vocabulary": "DEMOG", "codes": ["female"],
             "transform": "presence", "coefficient": 0.5},
            {"feature_key": "SURVEY:COMPASS31:total", "domain": "survey", "vocabulary": "SURVEY",
             "codes": ["COMPASS31:total"], "transform": "standardized_value", "center": 30, "scale": 15,
             "coefficient": 1.0},
            {"feature_key": "DEVICE:evie:resting_hr_mean", "domain": "device", "vocabulary": "DEVICE",
             "codes": ["evie:resting_hr_mean"], "transform": "standardized_value", "center": 65, "scale": 8,
             "coefficient": 0.5, "label": "ring resting heart rate"},
            {"feature_key": "DEVICE:evie:hrv_rmssd", "domain": "device", "vocabulary": "DEVICE",
             "codes": ["evie:hrv_rmssd"], "transform": "standardized_value", "center": 40, "scale": 15,
             "coefficient": -0.5},
        ],
        "decision_threshold": 0.5,
        "performance": {"auroc": 0.71, "n_positive": 40, "n_negative": 60},
        "strata_fractions": [{"age_band": "all", "sex": "all", "n_cases": 100, "n_positive": 40, "fraction": 0.4,
                              "suppressed": False}],
        "device_positive_centroid": {"ICD10CM:I73": 0.6, "LOINC:2093-3": 200.0, "DEMOG:female": 0.8},
        "caveats": ["synthetic test phenotype"],
    }
    ph.update(over)
    for f in ph["features"]:
        f.setdefault("label", f["feature_key"])
    return ph


def _local_scorer(ph, X):
    return reference.linear_predictor(ph, X) >= phenotype_io.threshold_lp(ph)


# --------------------------------------------------------------------------------------------- estimator


TOY = pd.DataFrame({"h": [1, 1, 1, 2, 2, 2], "j": [1, 1, 2, 1, 2, 2], "w": [1, 2, 1, 2, 2, 1],
                    "y": [1, 0, 1, 1, 0, 0]})


def test_taylor_proportion_matches_hand_computation():
    e = survey.domain_estimate(TOY.y, TOY.w, TOY.h, TOY.j)
    assert e.estimate == pytest.approx(4 / 9)
    # PSU totals of z = w(y-p)/W: S1 (-3, 5)/81, S2 (10, -12)/81 -> v = 2*(32 + 242)/6561 = 548/6561
    assert e.se == pytest.approx(math.sqrt(548) / 81)
    assert e.df == 2 and e.n_unweighted == 6 and e.sum_weights == 9
    from scipy import stats
    q = stats.t.ppf(0.975, 2)
    assert e.ci_low == pytest.approx(max(0, 4 / 9 - q * e.se)) and e.ci_high == pytest.approx(min(1, 4 / 9 + q * e.se))


def test_domain_estimate_keeps_full_design():
    dom = np.array([1, 1, 1, 1, 1, 0], dtype=bool)   # drop the last unit (S2, PSU 2)
    e = survey.domain_estimate(TOY.y, TOY.w, TOY.h, TOY.j, dom)
    assert e.estimate == pytest.approx(0.5)
    # z: S1 (-0.5, 0.5)/8, S2 (1, -1)/8 -> v = 2*(0.25 + 0.25)/64 + 2*(1 + 1)/64 = 5/64
    assert e.se == pytest.approx(math.sqrt(5) / 8)
    # subsetting the data first (wrong for variance) would give a different SE here only by chance; the domain
    # estimate must equal the estimate computed with y=0, w kept and z=0 outside the domain
    y2 = TOY.y.astype(float).where(dom)
    assert survey.domain_estimate(y2, TOY.w, TOY.h, TOY.j).se == pytest.approx(e.se)


def test_mean_estimate_and_lonely_psu():
    e = survey.domain_estimate([2.0, 4.0, 6.0], [1, 1, 2], [1, 1, 2], [1, 2, 1], proportion=False)
    assert e.estimate == pytest.approx(4.5)
    # stratum 2 has one PSU: removed (contributes 0) by default
    z = np.array([1 * (2 - 4.5), 1 * (4 - 4.5), 2 * (6 - 4.5)]) / 4
    assert e.se == pytest.approx(math.sqrt(2 * ((z[0] - z[:2].mean()) ** 2 + (z[1] - z[:2].mean()) ** 2)))
    e2 = survey.domain_estimate([2.0, 4.0, 6.0], [1, 1, 2], [1, 1, 2], [1, 2, 1], proportion=False,
                                lonely_psu="adjust")
    assert e2.se > e.se


def test_large_sample_recovers_known_weighted_prevalence():
    rng = np.random.default_rng(1)
    n = 20000
    strata = rng.integers(0, 15, n)
    psu = rng.integers(0, 2, n)
    grp = rng.random(n) < 0.3                          # 30% of the sample, weight 1; others weight 4
    w = np.where(grp, 1.0, 4.0)
    y = np.where(grp, rng.random(n) < 0.6, rng.random(n) < 0.1).astype(float)
    true_p = (0.3 * 1 * 0.6 + 0.7 * 4 * 0.1) / (0.3 * 1 + 0.7 * 4)
    e = survey.domain_estimate(y, w, strata, psu)
    assert abs(e.estimate - true_p) < 3 * e.se + 1e-3
    assert e.ci_low < true_p < e.ci_high
    assert e.nchs_reliability == "reliable"


def test_korn_graubard_and_nchs_rules():
    lo, hi, neff = survey.korn_graubard(0.1, 0.01, 2000, 30)
    assert lo < 0.1 < hi and neff == pytest.approx(900)
    assert survey.nchs_reliability(0.1, 0.0, 0.4, 500, 30).startswith("suppress: Korn-Graubard")
    assert survey.nchs_reliability(0.1, 0.05, 0.15, 20, 30).startswith("suppress: effective")
    assert survey.nchs_reliability(0.1, 0.08, 0.12, 500, 5).startswith("review")


def test_weighted_quantiles():
    q = survey.weighted_quantiles([1, 2, 3, 4], [1, 1, 1, 1], [0.5])
    assert q[0] == pytest.approx(2.5)
    q = survey.weighted_quantiles([1, 2], [1, 3], [0.5])
    assert 1.5 < q[0] <= 2


# --------------------------------------------------------------------------------------------- coverage


def test_coverage_overall_and_by_domain():
    ph = _phenotype()
    ev = {"ICD10CM:I73": True, "LOINC:2093-3": True, "DEMOG:female": True}
    cov = phenotype_io.coverage(ph, ev)
    total = 2 + 0.5 + 0.5 + 0.5 + 1 + 0.5 + 0.5
    assert cov["coverage"] == pytest.approx(3.0 / total)
    assert cov["by_domain"]["survey"]["evaluable_share_of_domain"] == 0
    assert cov["by_domain"]["condition"]["evaluable_share_of_domain"] == 1
    assert cov["n_evaluable"] == 3 and cov["n_features"] == 7


def test_threshold_and_base_population():
    ph = _phenotype(decision_threshold=0.25)
    assert phenotype_io.threshold_lp(ph) == pytest.approx(math.log(1 / 3))
    bp = phenotype_io.base_population(ph)
    assert bp["icd10cm"] == ["U09.9"] and not bp["proxy"]
    ptlds = phenotype_io.base_population({"condition_id": "ptlds", "base_population": {}})
    assert ptlds["proxy"] and ptlds["icd10cm"] == ["A69.2"] and "PROXY" in ptlds["label"]
    lyme = phenotype_io.base_population({"condition_id": "lyme_disease", "base_population": {}})
    assert lyme["icd10cm"] == ["A69.2"]


# --------------------------------------------------------------------------------------------- reference


def _synthetic_reference(seed=7, n_per_psu=60):
    """2 strata x 2 PSUs. Features chosen so the phenotype rule is positive exactly when I73 present
    (lp = -1 + 2*I73 + 0.5*z_chol + 0.5*female with |0.5 z| < 0.5 by construction)."""
    rng = np.random.default_rng(seed)
    rows = []
    sid = 0
    for h in (1, 2):
        for j in (1, 2):
            for _ in range(n_per_psu):
                sid += 1
                rows.append({"participant_id": f"ref:{sid}", "seqn": sid, "cycle_code": "G",
                             "age_years": float(rng.integers(18, 80)), "sex": "male",
                             "mec_examined": True, "sdmvstra": h, "sdmvpsu": j,
                             "wtmec4yr_pooled": float(rng.choice([1000.0, 3000.0]))})
    d = pd.DataFrame(rows)
    i73 = (rng.random(len(d)) < 0.25).astype(float)
    X = pd.DataFrame({"ICD10CM:I73": i73,
                      "LOINC:2093-3": 190 + 40 * rng.uniform(-0.9, 0.9, len(d)),
                      "DEMOG:female": 0.0}, index=d.participant_id)
    X.loc[X.index[::10], "LOINC:2093-3"] = np.nan          # item missingness -> centre
    return d, X


def test_reference_prevalence_equals_weighted_share_with_known_rule():
    d, X = _synthetic_reference()
    spec = reference.REFERENCES["nhanes"]
    res = reference.analyse_reference(_phenotype(), X, d, spec, scorer=_local_scorer)
    prev = res["prevalence"]
    nat = prev.query("variant == 'all_available_features' and age_band == 'all' and sex == 'all'").iloc[0]
    w = d.set_index("participant_id")["wtmec4yr_pooled"]
    expected = float((w * X["ICD10CM:I73"]).sum() / w.sum())
    assert nat["prevalence"] == pytest.approx(expected)
    e = survey.domain_estimate(X["ICD10CM:I73"].reindex(d.participant_id).values, w.values, d.sdmvstra, d.sdmvpsu)
    assert nat["se"] == pytest.approx(e.se)
    assert nat["n_unweighted"] == len(d) and nat["n_positive_unweighted"] == int(X["ICD10CM:I73"].sum())
    total = 5.5
    assert nat["feature_coverage"] == pytest.approx(3.0 / total)
    cbd = json.loads(nat["coverage_by_domain"])
    assert cbd["survey"] == 0 and cbd["device"] == 0 and cbd["measurement"] == 1
    assert nat["weighted_share_complete_features"] < 1
    # NMR-mapped cholesterol present -> a second variant without it
    assert set(prev["variant"]) == {"all_available_features", "excluding_nmr_mapped"}
    nmr = prev.query("variant == 'excluding_nmr_mapped' and age_band == 'all' and sex == 'all'").iloc[0]
    assert nmr["feature_coverage"] == pytest.approx(2.5 / total)
    # strata: female rows are empty in this all-male synthetic sample
    fem = prev.query("variant == 'all_available_features' and age_band == 'all' and sex == 'female'").iloc[0]
    assert fem["n_unweighted"] == 0 and np.isnan(fem["prevalence"])


def test_reference_profile_and_distance_are_aggregate():
    d, X = _synthetic_reference()
    res = reference.analyse_reference(_phenotype(), X, d, reference.REFERENCES["nhanes"], scorer=_local_scorer)
    prof = res["profile"].set_index("feature_key")
    assert prof.loc["ICD10CM:I73", "mean_phenotype_positive"] == pytest.approx(1.0)
    assert prof.loc["ICD10CM:I73", "mean_others"] == pytest.approx(0.0)
    assert not prof.loc["SURVEY:COMPASS31:total", "available_in_reference"]
    assert prof.loc["LOINC:2093-3", "nmr_platform_mismatch"]
    dist = res["distance"].set_index("group")
    assert dist.loc["all_adults", "n_features_shared"] == 3
    assert 0 <= dist.loc["all_adults", "q50"] <= 1
    # positives (I73 = 1, centroid 0.6) sit closer on that feature than others (I73 = 0)
    assert dist.loc["phenotype_positive", "q50"] < dist.loc["others", "q50"]
    for df in (res["prevalence"], res["profile"], res["distance"]):
        assert "participant_id" not in df.columns and "seqn" not in df.columns


def test_age_band_topcode():
    b = reference.age_band(pd.Series([18, 25, 49, 79, 80, 85]), 80)
    assert list(b) == ["18-19", "20-29", "40-49", "70-79", "80+", "80+"]
    b = reference.age_band(pd.Series([81, 85]), 85)
    assert list(b) == ["80+", "80+"]


def test_gower_binary_and_numeric():
    X = pd.DataFrame({"a": [0.0, 1.0, np.nan], "b": [0.0, 10.0, 5.0]})
    dist, used = reference.gower_to_centroid(X, {"a": 1.0, "b": 0.0, "zz": 3.0}, {"a"})
    assert used == ["a", "b"]
    assert dist.iloc[0] > dist.iloc[1]
    assert not np.isnan(dist.iloc[2])                   # missing 'a' skipped, 'b' used


def test_reduced_intercept_and_standardized_weights():
    ph = _phenotype()
    ph["features"][0]["center"] = 0.2                      # I73 presence, training prevalence 0.2
    ph["features"][0]["coefficient_standardized"] = 4.0     # coverage weight overrides |coefficient| = 2
    ev = {"LOINC:2093-3": True}
    assert phenotype_io.reduced_intercept(ph, ev) == pytest.approx(-1.0 + 2.0 * 0.2)
    cov = phenotype_io.coverage(ph, ev)
    assert cov["coverage"] == pytest.approx(0.5 / (4 + 0.5 + 0.5 + 0.5 + 1 + 0.5 + 0.5))
    # the local rule holds a missing presence feature at its centre, like harmonize.phenotype.score_phenotype
    X = pd.DataFrame({"LOINC:2093-3": [190.0]})
    assert reference.linear_predictor(ph, X).iloc[0] == pytest.approx(-1.0 + 2.0 * 0.2)


def test_same_item_survey_crosswalk_only():
    ph = _phenotype(features=[
        {"feature_key": "SURVEY:HRQOL:gh", "vocabulary": "SURVEY", "domain": "survey", "codes": ["HRQOL:gh"],
         "transform": "standardized_value", "center": 2, "scale": 1, "coefficient": 1.0,
         "reference_crosswalk": [{"reference": "CDC_HRQOL4:general_health", "match_quality": "same_item"}]},
        {"feature_key": "SURVEY:FSS9:total", "vocabulary": "SURVEY", "domain": "survey", "codes": ["FSS9:total"],
         "transform": "standardized_value", "center": 4, "scale": 1, "coefficient": 1.0,
         "reference_crosswalk": [{"reference": "PHQ9:item4_tired", "match_quality": "related_construct"}]}])
    assert reference.same_item_crosswalk(ph) == {"SURVEY:HRQOL:gh": "SURVEY:CDC_HRQOL4:general_health"}


# --------------------------------------------------------------------------------------------- OMOP


def test_icd_category_prefix_from_code_match():
    f = {"feature_key": "ICD10CM:G93", "codes": ["G93.32", "G93.31"], "code_match": "category_prefix"}
    assert omop.icd_codes(f) == ["G93"]
    assert omop.icd_codes(dict(f, code_match="exact")) == ["G93.32", "G93.31"]
    assert omop._code_predicate("c", ["G93"]) == "(c LIKE 'G93%')"


def test_omop_sql_generation_and_coverage():
    ph = _phenotype()
    sql, info = omop.build_sql(ph, dialect="bigquery", prefix="`{CDR}`.")
    assert "vocabulary_id = 'ICD10CM'" in sql and "'Maps to'" in sql and "concept_ancestor" in sql
    assert "LIKE 'I73%'" in sql and "'2093-3'" in sql and "concept_code = 'mg/dL'" in sql
    assert "BETWEEN 1 AND 20" in sql and "complementary" in sql
    assert "COMPASS31" not in sql.split("WITH", 1)[1]           # survey not evaluable in OMOP
    assert info["coverage"]["coverage"] == pytest.approx(3.5 / 5.5)
    assert not info["evaluable"]["SURVEY:COMPASS31:total"] and not info["evaluable"]["DEVICE:evie:hrv_rmssd"]
    spark, _ = omop.build_sql(ph, dialect="spark")
    assert "percentile_approx" in spark and "BETWEEN 1 AND 19" in spark and "location" in spark
    fb, finfo = omop.build_sql(ph, dialect="bigquery", prefix="`{CDR}`.", fitbit=True)
    assert [x["fitbit_proxy"] for x in finfo["fitbit"]] == ["resting_hr"]   # HRV has no Fitbit table
    assert finfo["coverage"]["coverage"] == pytest.approx(4.0 / 5.5)


def test_omop_ptlds_uses_labelled_proxy():
    ph = _phenotype(condition_id="ptlds", base_population={})
    sql, info = omop.build_sql(ph, dialect="spark")
    assert info["base_population"]["proxy"] and "LIKE 'A69.2%'" not in sql  # A69.2 is an exact code here
    assert "'A69.2'" in sql and "PROXY" in sql


def test_sql_parses_in_both_dialects(tmp_path):
    pytest.importorskip("sqlglot")
    ph = _phenotype()
    m = omop.write_query_pack(ph, tmp_path)
    for name in m["files"].values():
        if name.endswith(".sql"):
            omop.parse_check((tmp_path / name).read_text(), "spark" if name.startswith("n3c") else "bigquery")


def test_clinic_export(tmp_path):
    ph = _phenotype()
    out = omop.clinic_export(ph, tmp_path / "clinic")
    icd = pd.read_csv(tmp_path / "clinic" / "icd10cm.csv")
    loinc = pd.read_csv(tmp_path / "clinic" / "loinc.csv")
    ne = pd.read_csv(tmp_path / "clinic" / "not_evaluable_in_ehr.csv")
    assert icd["code"].tolist() == ["I73"] and loinc.loc[0, "unit"] == "mg/dL"
    assert "NMR" in loinc.loc[0, "platform_caveat"]
    assert set(ne["feature_key"]) == {"SURVEY:COMPASS31:total", "DEVICE:evie:resting_hr_mean", "DEVICE:evie:hrv_rmssd"}
    readme = (tmp_path / "clinic" / "README.md").read_text()
    assert "Scoring rule" in readme and f"{out['coverage']['coverage']:.0%}" in readme
    assert "diagnosed by AI" not in readme and "patient has" not in readme


# --------------------------------------------------------------------------------------------- end to end


def test_run_end_to_end_on_temporary_tree(tmp_path, monkeypatch):
    import measure_it.store as store
    from measure_it.similar import cli
    proc, res_root = tmp_path / "processed", tmp_path / "byod"
    proc.mkdir()
    d, X = _synthetic_reference()
    d.assign(data_layer="person").to_parquet(proc / "participants__nhanes.parquet")
    pd.DataFrame({"participant_id": ["x"], "vocabulary": ["LOINC"], "concept_code": ["2093-3"],
                  "source_variable": ["LBXTC"]}).to_parquet(proc / "person_concepts__nhanes.parquet")
    (res_root / "synth").mkdir(parents=True)
    (res_root / "synth" / "computable_phenotype.json").write_text(json.dumps(_phenotype()))
    for mod in (store, reference, cli):
        monkeypatch.setattr(mod, "PROCESSED", proc, raising=False)
    monkeypatch.setattr(phenotype_io, "RESULTS_ROOT", res_root)
    monkeypatch.setattr(reference, "_feature_matrix", lambda ds, keys, xw=None: X[[k for k in keys if k in X.columns]])
    monkeypatch.setattr(reference, "score", _local_scorer)
    s = cli.run("byod_synth", references=("nhanes", "nhanes0306"))
    ent = s["phenotypes"][0]
    assert ent["references"]["nhanes"]["status"] == "ok"
    assert ent["references"]["nhanes0306"]["status"] == "skipped"
    prev = pd.read_parquet(proc / "similar_reference_prevalence.parquet")
    assert set(prev["data_layer"]) == {"derived"} and set(prev["evidence_type"]) == {"survey_estimate"}
    assert prev["object_id"].is_unique and prev["reference_population"].str.contains("no Long COVID item").all()
    assert "participant_id" not in prev.columns
    out = res_root / "synth" / "similar"
    assert (out / ".gitignore").read_text().strip().endswith("*")
    assert any(p.name.startswith("aou_bigquery__") for p in out.iterdir())
    # re-run replaces this phenotype's rows instead of appending
    cli.run("byod_synth", references=("nhanes",))
    assert len(pd.read_parquet(proc / "similar_reference_prevalence.parquet")) == len(prev)


# --------------------------------------------------------------------------------------------- real data


@pytest.mark.data
@pytest.mark.skipif(not (PROCESSED / "participants__nhanes.parquet").exists(), reason="NHANES not ingested")
def test_nhanes_design_columns_and_estimator_on_real_data():
    spec = reference.REFERENCES["nhanes"]
    d = reference.load_design(spec)
    c = pd.read_parquet(PROCESSED / "participant_clinical_features__nhanes.parquet",
                        columns=["participant_id", "dx_hypertension"])
    d = d.merge(c, on="participant_id")
    dom = (d.age_years >= 20) & d.sex.eq("female") & d.mec_examined
    e = survey.domain_estimate(d.dx_hypertension, d[spec.mec_weight], d.sdmvstra, d.sdmvpsu, dom)
    # self-reported hypertension, women 20+, 2011-2014: about a third (NCHS Data Brief 220: 29-33% measured/aware)
    assert 0.25 < e.estimate < 0.40 and e.df >= 25 and e.nchs_reliability == "reliable"


@pytest.mark.data
@pytest.mark.skipif(not (PROCESSED / "person_concepts__nhanes.parquet").exists(), reason="harmonize not built")
def test_real_score_phenotype_agrees_with_local_rule():
    """The stage-BC scorer and this module's local re-implementation (used by the synthetic tests) agree on NHANES."""
    from measure_it.harmonize.features import feature_matrix
    keys = ["DEMOG:female", "DEMOG:age_mid", "SURVEY:PHQ9:total", "ICD10CM:I10"]
    X = feature_matrix(["nhanes"], keys)
    ph = _phenotype(intercept=-0.7, features=[
        {"feature_key": "DEMOG:female", "vocabulary": "DEMOG", "domain": "demographic", "codes": [],
         "transform": "presence", "center": 0.5, "coefficient": 0.6},
        {"feature_key": "DEMOG:age_mid", "vocabulary": "DEMOG", "domain": "demographic", "codes": [],
         "transform": "standardized_value", "center": 45.0, "scale": 15.0, "coefficient": 0.4},
        {"feature_key": "SURVEY:PHQ9:total", "vocabulary": "SURVEY", "domain": "survey", "codes": ["PHQ9:total"],
         "transform": "standardized_value", "center": 5.0, "scale": 5.0, "coefficient": 0.8},
        {"feature_key": "ICD10CM:I10", "vocabulary": "ICD10CM", "domain": "condition", "codes": ["I10"],
         "transform": "presence", "center": 0.3, "coefficient": 0.5},
        {"feature_key": "SURVEY:COMPASS31:total", "vocabulary": "SURVEY", "domain": "survey",
         "codes": ["COMPASS31:total"], "transform": "standardized_value", "center": 30, "scale": 15,
         "coefficient": 1.0}])
    real = reference.score(ph, X)
    local = _local_scorer(ph, X)
    assert len(real) == len(X) and (real.values == local.values).all()


def test_low_coverage_suppresses_prevalence_and_warns_in_sql():
    d, X = _synthetic_reference()
    ph = _phenotype()
    ph["features"][4]["coefficient"] = 20.0                 # survey item (absent from NHANES) dominates the mass
    res = reference.analyse_reference(ph, X, d, reference.REFERENCES["nhanes"], scorer=_local_scorer)
    nat = res["prevalence"].query("variant == 'all_available_features' and age_band == 'all' and sex == 'all'").iloc[0]
    assert nat["coverage_flag"].startswith("insufficient") and np.isnan(nat["prevalence"])
    assert nat["nchs_reliability"].startswith("suppress: feature coverage")
    assert phenotype_io.coverage_flag(0.4).startswith("low") and phenotype_io.coverage_flag(0.9) == "adequate"
    sql, info = omop.build_sql(ph, dialect="spark")
    assert sql.startswith("-- WARNING: feature coverage") and info["coverage"]["flag"].startswith("insufficient")


def test_suppression_covers_the_phenotype_negative_difference():
    """Regression (found by the first local execution, agent.omop_run): with n_base 30 and n_positive 28 the pack used
    to publish both, so 2 phenotype-negative people were recoverable by subtraction. Execute the pack's own
    suppression CTEs and final SELECT (DuckDB) on hand-made counts."""
    import duckdb
    ph = _phenotype()
    for dialect in ("spark", "bigquery"):
        sql, _ = omop.build_sql(ph, dialect=dialect)
        tail = sql[sql.index("flagged AS ("):]
        counts = ("by_level AS (SELECT * FROM (VALUES ('all','all',30,28), ('state','SC',25,25), "
                  "('state','GA',5,3), ('state','TX',400,150)) t(level, stratum, n_base, n_positive))")
        rows = duckdb.connect().sql("WITH " + counts + ",\n" + tail).fetchall()
        got = {(r[0], r[1]): (r[2], r[3]) for r in rows}
        assert got[("all", "all")][1] is None                 # 30 - 28 = 2 negatives: positive count suppressed
        assert got[("state", "GA")] == (None, None)            # base 5
        assert got[("state", "TX")] == (400, 150)              # large cells unchanged
