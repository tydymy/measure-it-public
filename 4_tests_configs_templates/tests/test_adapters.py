"""Tests for the measurement adapters, the adapter registry and the downstream adapter contract.

The contract test uses `SyntheticFakeAdapter`, defined here. Its data are synthetic and exist ONLY inside
this test module; nothing synthetic is written to data/processed.
"""
from __future__ import annotations

import inspect

import numpy as np
import pandas as pd
import pytest

from measure_it.config import UNKNOWN
from measure_it.measurements import adapters as A
from measure_it.measurements.adapter import MeasurementAdapter, MeasurementDescription
from measure_it.measurements.capillaroscopy_adapter import CapillaroscopyAdapter
from measure_it.store import read_table, table_exists
from measure_it.wearables import adapter as W
from measure_it.wearables import nhanes_features as NF


# ------------------------------------------------------------------------------------------ helpers
def test_weighted_quantiles_simple():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.ones(4)
    assert W.weighted_quantiles(x, w, [0.5])[0] == pytest.approx(2.5)
    # doubling a weight moves the median towards that value
    assert W.weighted_quantiles(x, np.array([1, 1, 1, 5.0]), [0.5])[0] > 3.0
    assert np.isnan(W.weighted_quantiles([np.nan], [1.0], [0.5])[0])


def test_classify_result_rules():
    assert W.classify_result({"auroc": 0.55, "ci_low": 0.49, "ci_high": 0.6, "perm_p_two_sided": 0.01}).startswith("null")
    assert W.classify_result({"auroc": 0.39, "ci_low": 0.28, "ci_high": 0.49, "perm_p_two_sided": 0.07}).startswith("inconclusive")
    assert W.classify_result({"auroc": 0.6, "ci_low": 0.58, "ci_high": 0.62, "perm_p_two_sided": 0.001}).startswith("association")


def test_age_band_edges():
    assert W.age_band(6) == "6-11" and W.age_band(11.5) == "6-11" and W.age_band(17.9) == "12-17"
    assert W.age_band(18) == "18-29" and W.age_band(79.9) == "70-79"
    assert W.age_band(80) == "80+" and W.age_band(5) is None and W.age_band(np.nan) is None


# ------------------------------------------------------------------------------------------ WearableAdapter (unit)
def _toy_minutes(n_days=5, pid="toyset:1"):
    """Tiny hand-checkable minute frame: full days, wake 06:00-22:00 with MIMS 20, sleep otherwise MIMS 1."""
    rows = []
    for d in range(n_days):
        for m in range(1440):
            wake = 360 <= m < 1320
            rows.append((pid, d * 1440 + m, 20.0 if wake else 1.0, 1 if wake else 2, 0, (d % 7) + 1))
    return pd.DataFrame(rows, columns=["participant_id", "abs_min", "mims", "pred", "qf", "dow"])


def test_preprocess_requires_columns_and_index():
    ad = W.WearableAdapter()
    with pytest.raises(ValueError):
        ad.preprocess(pd.DataFrame({"participant_id": ["a:1"], "mims": [1.0]}))
    with pytest.raises(ValueError):
        ad.preprocess(pd.DataFrame({"participant_id": ["a:1"], "mims": [1.0], "pred": [1]}))


def test_preprocess_timestamp_derives_abs_min_and_nhanes_dow():
    ts = pd.to_datetime(["2024-06-02 23:59", "2024-06-03 00:00"])  # a Sunday then a Monday
    df = pd.DataFrame({"participant_id": ["a:1", "a:1"], "timestamp": ts, "mims": [1.0, 2.0], "pred": [1, 1]})
    c = W.WearableAdapter().preprocess(df)
    assert c["abs_min"].tolist() == [1439, 1440]
    assert c["dow"].tolist() == [1, 2]  # NHANES coding: Sunday = 1, Monday = 2
    assert c["day_index"].tolist() == [1, 2]


def test_preprocess_restores_mims_precision_at_threshold():
    df = pd.DataFrame({"participant_id": ["a:1"], "abs_min": [0], "mims": [np.float32(10.558)], "pred": [1]})
    c = W.WearableAdapter().preprocess(df)
    assert bool(c["active"].iloc[0])  # float32(10.558) < 10.558 unless restored


def test_day_table_definitions():
    raw = _toy_minutes(n_days=2)
    raw.loc[5, "qf"] = 3            # a QC-flagged wake? no: minute 5 is sleep -> removed from sleep count
    raw.loc[400, "mims"] = -0.01    # MIMS not computable in a wake minute: counted as wake wear, not in MIMS sum
    d = W.WearableAdapter.day_table(W.WearableAdapter().preprocess(raw))
    d1 = d[d["day_index"] == 1].iloc[0]
    assert d1["minutes_with_data"] == 1440
    assert d1["wake_wear_min"] == 960
    assert d1["sleep_wear_min"] == 480 - 1
    assert d1["qc_flag_count"] == 3
    assert d1["total_mims"] == pytest.approx(959 * 20.0 + 479 * 1.0)
    assert bool(d1["full_day"]) and bool(d1["valid_day"])


def test_embed_toy_participant_passes_rule_and_has_all_features():
    emb = W.WearableAdapter().embed(_toy_minutes(n_days=5))
    assert list(emb.columns[:2]) == ["participant_id", "passes_valid_wear_rule"]
    assert set(W.FEATURE_COLUMNS) <= set(emb.columns)
    r = emb.iloc[0]
    assert bool(r["passes_valid_wear_rule"]) and r["n_valid_days"] == 5
    assert r["active_min_per_day"] == pytest.approx(960)  # MIMS 20 >= 10.558 for every wake minute
    assert r["sedentary_fraction"] == pytest.approx(0.0)
    assert r["interdaily_stability"] == pytest.approx(1.0)  # identical days


def test_embed_short_record_fails_rule():
    emb = W.WearableAdapter().embed(_toy_minutes(n_days=3))
    assert not bool(emb.iloc[0]["passes_valid_wear_rule"])
    assert np.isnan(emb.iloc[0]["interdaily_stability"])


def test_phenotype_score_against_synthetic_reference():
    """A participant at every stratum median scores 0 on every axis (reference built in-test)."""
    feats = sorted({f for comps in W.AXES.values() for f, _ in comps})
    rows = []
    for sex in ("female", "male", "all"):
        for _, _, band in W.AGE_BANDS + [(0, 0, "all_6plus")]:
            for f in feats:
                rows.append({"sex": sex, "age_band": band, "feature": f, "feature_kind": "feature",
                             "w_median": 1.0, "robust_scale": 0.5, "stratum_fallback": ""})
            for ax in W.AXES:
                rows.append({"sex": sex, "age_band": band, "feature": f"axis:{ax}", "feature_kind": "axis_composite",
                             "w_median": 0.0, "robust_scale": 1.0, "stratum_fallback": ""})
    ref = pd.DataFrame(rows)
    emb = pd.DataFrame({"participant_id": ["x:1", "x:2", "x:3", "x:4"], "passes_valid_wear_rule": [True, True, False, True],
                        "age_years": [45, np.nan, 45, 4], "sex": ["female", "male", "male", "female"],
                        **{f: [1.0, 1.5, 1.0, 1.0] for f in feats}})
    sc = W.WearableAdapter(reference=ref).phenotype_score(emb)
    assert sc.loc[0, list(W.AXES)].abs().max() == pytest.approx(0.0)
    assert sc.loc[1, "reference_stratum"] == "all|all_6plus"         # missing age -> pooled, flagged
    assert sc.loc[1, "activity_fragmentation"] == pytest.approx(1.0)  # astp (+): (1.5 - 1) / 0.5
    assert sc.loc[1, "low_activity"] == pytest.approx((1 - 1 - 1) / 3)  # mean of +1, -1, -1 component z
    assert sc.loc[2, "score_status"] == "insufficient_valid_wear" and sc.loc[2, list(W.AXES)].isna().all()
    assert sc.loc[3, "score_status"] == "outside_reference_age_range" and sc.loc[3, list(W.AXES)].isna().all()
    for ax, comps in W.AXES.items():
        for f, _ in comps:
            assert f"{ax}__z__{f}" in sc.columns  # every score keeps its components


def test_describe_measurement_ids_are_registry_ids():
    d = W.WearableAdapter().describe_measurement()
    assert d.status == "implemented" and not d.requires_private_data
    assert d.measurement_ids == ["accelerometry"]
    assert d.embedding_features == list(NF.FEATURE_DEFINITIONS)
    assert A.validate_description(d) == []
    assert A.validate_description(W.WearableHeartRateAdapter().describe_measurement()) == []


def test_hr_personal_baseline_features_synthetic():
    dates = pd.date_range("2030-01-01", periods=60, freq="D")
    rhr = np.r_[np.full(50, 60.0) + np.tile([-1.0, 1.0], 25), np.full(10, 70.0)]
    daily = pd.DataFrame({"participant_id": "s:1", "date": dates, "rhr_mean": rhr, "valid_wear_day": True,
                          "steps": np.nan, "rhr_night_mean": np.nan, "hr_minutes": 1000})
    f = W.personal_baseline_features(daily)
    assert f["personal_rhr_median"] == pytest.approx(60.0, abs=1.0)
    assert f["rhr_p95_7d_personal_z"] > 2.5  # the last 10 days are far above the personal baseline
    assert np.isnan(f["steps_p5_7d_personal_z"])
    emb = W.WearableHeartRateAdapter().embed(daily)
    sc = W.WearableHeartRateAdapter().phenotype_score(emb)
    assert sc.loc[0, "rhr_elevation_vs_personal_baseline"] > 2.5
    assert np.isnan(sc.loc[0, "activity_reduction_vs_personal_baseline"])


def test_pooled_fallback_reference_is_consistent_with_pooled_scoring():
    """A participant without age/sex is scored against the pooled reference; the pooled composite reference
    must be built the same way, so scoring a whole reference population pooled gives median ~0, scale ~1."""
    rng = np.random.default_rng(0)
    n = 6000
    age = rng.integers(6, 86, n).astype(float)
    sex = rng.choice(["female", "male"], n)
    feats = sorted({f for comps in W.AXES.values() for f, _ in comps})
    pop = pd.DataFrame({"participant_id": [f"synthref:{i}" for i in range(n)], "age_years": age, "sex": sex,
                        W.REFERENCE_WEIGHT: np.ones(n), "passes_valid_wear_rule": True})
    for f in feats + [f for f in W.REFERENCE_EXTRA_FEATURES if f not in feats]:
        pop[f] = 5.0 * age / 80.0 + rng.normal(size=n)  # strong age dependence: pooled z wider than stratum z
    pop["age_band"] = pop["age_years"].map(W.age_band)
    ref = W.build_reference(pop)
    anon = pop.assign(age_years=np.nan, sex=None)
    sc = W.WearableAdapter(reference=ref).phenotype_score(anon)
    assert (sc["reference_stratum"] == "all|all_6plus").all()
    for ax in W.AXES:
        v = sc[ax].dropna()
        q25, q50, q75 = np.percentile(v, [25, 50, 75])
        assert abs(q50) < 0.05, ax
        assert (q75 - q25) / W.IQR_TO_SD == pytest.approx(1.0, abs=0.05), ax


def test_weighted_auroc_and_design_bootstrap():
    rng = np.random.default_rng(1)
    y = rng.integers(0, 2, 400).astype(bool)
    s = rng.normal(size=400) + 0.8 * y
    assert W.weighted_auroc(y, s, np.ones(400)) == pytest.approx(A.auroc(y, s))
    w = np.where(y, 1.0, 1.0)
    w[:200] = 3.0  # weights only on part of the sample: still a valid AUROC in [0, 1]
    assert 0 <= W.weighted_auroc(y, s, w) <= 1
    stratum = np.repeat(np.arange(20), 20)
    psu = np.tile(np.repeat([1, 2], 10), 20)
    lo, hi = W.design_bootstrap_auroc_ci(y, s, stratum, psu, n_boot=200)
    assert lo < A.auroc(y, s) < hi


# ------------------------------------------------------------------------------------------ stub
def test_capillaroscopy_stub_raises_and_describes():
    cap = CapillaroscopyAdapter()
    for fn in (cap.preprocess, cap.embed, cap.phenotype_score):
        with pytest.raises(NotImplementedError, match="CAPRIO must supply"):
            fn(pd.DataFrame())
    d = cap.describe_measurement()
    assert d.status == "stub" and d.requires_private_data is True
    assert d.measurement_ids == ["capillaroscopy", "microvascular_function"]
    assert {"capillary density", "capillary morphology", "microhemorrhages", "capillary flow"} == set(d.signals)
    assert d.public_training_data == [] and d.embedding_features == []
    assert len(d.phenotypes) == 4
    assert A.validate_description(d) == []


def test_stub_has_no_participant_scores():
    d = CapillaroscopyAdapter().describe_measurement()
    with pytest.raises(ValueError, match="stub"):
        A.scores_to_long(d, pd.DataFrame({"participant_id": []}), dataset_id="caprio")


# ------------------------------------------------------------------------------------------ registry
def test_registry_lists_and_gets_adapters():
    names = {a["name"] for a in A.list_adapters()}
    assert {"wearable", "wearable_heart_rate", "capillaroscopy"} <= names
    assert isinstance(A.get_adapter("wearable"), W.WearableAdapter)
    assert isinstance(A.get_adapter("wearable_accelerometry_v1"), W.WearableAdapter)
    assert isinstance(A.get_adapter("CapillaroscopyAdapter"), CapillaroscopyAdapter)
    with pytest.raises(KeyError):
        A.get_adapter("no_such_adapter")


def test_bundle_resolution_from_relevance_config():
    assert isinstance(A.adapter_for_bundle("wearable_autonomic_activity_monitoring"), W.WearableAdapter)
    assert isinstance(A.adapter_for_bundle("wearables"), W.WearableAdapter)        # alias
    assert isinstance(A.adapter_for_bundle("CAPRIO"), CapillaroscopyAdapter)       # alias
    r = A.resolve_bundle("CPET")
    assert r["status"] == "resolved" and r["adapter_name"] == UNKNOWN and A.adapter_for_bundle("CPET") is None
    assert A.resolve_bundle("not a bundle")["status"] == UNKNOWN
    w = A.resolve_bundle("wearable monitoring")
    assert "wearable_heart_rate" in w["adapters_covering_members"] and w["object_id"].startswith("measurement:")
    # SPEC demo wording is an alias of the wearable bundle
    assert A.resolve_bundle("wearable monitoring of autonomic/activity abnormalities")["bundle_id"] == \
        "wearable_autonomic_activity_monitoring"


def test_bundle_may_declare_several_adapters():
    w = A.resolve_bundle("wearable_autonomic_activity_monitoring")
    assert w["adapter_names"] == ["wearable", "wearable_heart_rate"] and w["adapter_name"] == "wearable"
    assert w["adapters_registered"] == ["wearable", "wearable_heart_rate"]
    got = A.adapters_for_bundle("wearables")
    assert [type(a) for a in got] == [W.WearableAdapter, W.WearableHeartRateAdapter]
    assert "wearable_autonomic_activity_monitoring" in next(
        a for a in A.list_adapters() if a["name"] == "wearable_heart_rate")["bundles_naming_this_adapter"]
    # both config forms, and none
    assert A.bundle_adapter_names({"adapter": "capillaroscopy"}) == ["capillaroscopy"]
    assert A.bundle_adapter_names({"adapters": ["wearable", "wearable_heart_rate"]}) == ["wearable", "wearable_heart_rate"]
    assert A.bundle_adapter_names({"adapters": ["wearable"], "adapter": "capillaroscopy"}) == ["wearable", "capillaroscopy"]
    assert A.bundle_adapter_names({"adapter": None}) == [] and A.bundle_adapter_names({}) == []
    assert A.adapters_for_bundle("CPET") == [] and A.adapters_for_bundle("not a bundle") == []


def test_axis_list_parsing_is_exact():
    assert A._axis_list(np.array(["fatigue", "sleep_disturbance"])) == ["fatigue", "sleep_disturbance"]
    assert A._axis_list('["chronic_fatigue"]') == ["chronic_fatigue"]
    assert A._axis_list("chronic_fatigue; tachycardia") == ["chronic_fatigue", "tachycardia"]
    assert A._axis_list(None) == [] and A._axis_list(np.nan) == []
    cp = pd.DataFrame({"canonical_condition_id": ["a", "b"], "hpo_id": ["HP:1", "HP:2"],
                       "association_scope": ["direct", "direct"], "phenotype_axes": ["chronic_fatigue", "fatigue"]})
    assert [c["canonical_condition_id"] for c in A._axis_conditions(cp, "fatigue")] == ["b"]


@pytest.mark.data
def test_condition_measurement_evidence_route_is_live_and_skips_empty_cells():
    if not table_exists("condition_measurement_evidence"):
        pytest.skip("condition_measurement_evidence not built")
    cme = read_table("condition_measurement_evidence")
    c = A.conditions_for_description(W.WearableAdapter().describe_measurement())
    route = c[c["route"] == "condition_measurement_evidence"]
    assert len(route) and (route["canonical_condition_id"] != UNKNOWN).all()
    assert route["object_id"].str.startswith("measurement_evidence:").all()
    acc = cme[cme["measurement_id"] == "accelerometry"]
    linked = set(acc.loc[acc["measurement_evidence_strength"] != "none", "condition_id"])
    empty = set(acc.loc[acc["measurement_evidence_strength"] == "none", "condition_id"])
    assert set(route["canonical_condition_id"]) == linked and not (set(route["canonical_condition_id"]) & empty)


def test_validate_description_flags_problems():
    bad = MeasurementDescription(adapter_id="x", modality="m", measurement_ids=["not_a_class"], signals=[],
                                 input_unit_of_observation="", embedding_features=[], phenotypes=[],
                                 public_training_data=[], status="beta")
    probs = " | ".join(A.validate_description(bad))
    assert "not_a_class" in probs and "phenotypes is empty" in probs and "status" in probs


# ------------------------------------------------------------------------------------------ contract test
class SyntheticFakeAdapter(MeasurementAdapter):
    """A made-up modality whose data are SYNTHETIC and exist only inside this test.

    It deliberately uses a measurement id and phenotype the wearable adapter does not, to show that the
    downstream functions depend only on the MeasurementDescription and the participant scores.
    """

    DATASET = "synthetic_fake_dataset"

    def preprocess(self, raw_data):
        return raw_data.dropna()

    def embed(self, raw_data):
        d = self.preprocess(raw_data)
        return d.groupby("participant_id", as_index=False)["signal"].agg(["mean", "std"]).rename(
            columns={"mean": "emb_mean", "std": "emb_sd"})

    def phenotype_score(self, embedding):
        z = (embedding["emb_mean"] - embedding["emb_mean"].median()) / embedding["emb_mean"].std()
        return pd.DataFrame({"participant_id": embedding["participant_id"],
                             "abnormal_circadian_rhythm": z,
                             "abnormal_circadian_rhythm__z__emb_mean": z,
                             "synthetic_unmapped_phenotype": -z,
                             "score_status": "scored"})

    def describe_measurement(self):
        return MeasurementDescription(
            adapter_id="synthetic_fake_v0", modality="synthetic test modality", measurement_ids=["retinal_imaging"],
            signals=["synthetic signal"], input_unit_of_observation="participant x synthetic sample",
            embedding_features=["emb_mean", "emb_sd"],
            phenotypes=["abnormal_circadian_rhythm", "synthetic_unmapped_phenotype"],
            public_training_data=[], requires_private_data=True, status="implemented",
            limitations=["synthetic data inside tests only"])


@pytest.fixture
def fake_registered():
    A.register_adapter("synthetic_fake", SyntheticFakeAdapter, overwrite=True)
    yield
    A.unregister_adapter("synthetic_fake")


def _synthetic_raw(n=40, seed=0):
    rng = np.random.default_rng(seed)
    pid = np.repeat([f"{SyntheticFakeAdapter.DATASET}:{i}" for i in range(n)], 10)
    return pd.DataFrame({"participant_id": pid, "signal": rng.normal(size=n * 10) + np.repeat(np.arange(n) / n, 10)})


@pytest.mark.data   # measurement_context reads the processed ontology, measurement and FDA tables
def test_contract_fake_adapter_flows_through_unchanged_downstream(fake_registered):
    fake = A.get_adapter("synthetic_fake")
    assert "synthetic_fake" in {a["name"] for a in A.list_adapters()}
    desc = fake.describe_measurement()
    assert A.validate_description(desc) == []
    # measurement layer: ids resolve against configs/measurements.yaml + relevance.yaml (curated)
    mids = A.resolve_measurement_ids(desc).set_index("measurement_id")
    assert mids.loc["retinal_imaging", "status"] == "resolved"
    assert mids.loc["retinal_imaging", "deployment_complexity"] == 2
    assert set(mids.loc["retinal_imaging", "implementer_groups"]) == {"primary_care", "fqhc",
                                                                      "general_acute_care_hospital"}
    # ontology layer: an axis-id phenotype resolves, an unknown one returns UNKNOWN with a reason
    ph = A.resolve_phenotypes(desc).set_index("phenotype")
    assert ph.loc["abnormal_circadian_rhythm", "link_type"] == "exact_axis_id"
    assert ph.loc["abnormal_circadian_rhythm", "hpo_id"] == "HP:0006979"
    assert ph.loc["synthetic_unmapped_phenotype", "status"] == UNKNOWN
    assert ph.loc["synthetic_unmapped_phenotype", "reason"]
    # condition layer + full context (same call the scoring / clinic layers make)
    ctx = A.measurement_context(desc)
    assert ctx["contract_problems"] == [] and ctx["implementer_groups"]
    assert isinstance(ctx["conditions"], list)
    # participant scores -> standard long rows; summary with (synthetic) labels
    scores = fake.phenotype_score(fake.embed(_synthetic_raw()))
    long = A.scores_to_long(desc, scores, dataset_id=SyntheticFakeAdapter.DATASET)
    assert list(long.columns) == A.LONG_COLUMNS
    assert set(long["phenotype"]) == set(desc.phenotypes) and len(long) == 2 * len(scores)
    assert long["components_json"].str.contains("emb_mean").any()
    labels = pd.Series((np.arange(40) >= 20).astype(float), index=[f"{SyntheticFakeAdapter.DATASET}:{i}" for i in range(40)])
    summ = A.summarize_scores(long, labels=labels).set_index("phenotype")
    assert summ.loc["abnormal_circadian_rhythm", "auroc"] > 0.5 > summ.loc["synthetic_unmapped_phenotype", "auroc"]


def test_contract_same_code_path_for_every_adapter(fake_registered):
    """measurement_context returns the same structure for the fake, the wearable and the stub adapters."""
    keys = None
    for name in ("synthetic_fake", "wearable", "wearable_heart_rate", "capillaroscopy"):
        ctx = A.measurement_context(A.get_adapter(name).describe_measurement())
        keys = keys or set(ctx)
        assert set(ctx) == keys
    cap = A.measurement_context(CapillaroscopyAdapter().describe_measurement())
    assert any("stub" in c for c in cap["caveats"])
    assert all(p["status"] == UNKNOWN for p in cap["phenotypes"])  # capillary phenotypes not yet in the ontology
    assert {"rheumatology", "vascular"} <= set(cap["implementer_groups"])


def test_contract_downstream_functions_are_modality_agnostic():
    """The downstream functions never mention a modality: they only read the description and scores."""
    for fn in (A.validate_description, A.resolve_measurement_ids, A.resolve_phenotypes, A.conditions_for_description,
               A.measurement_context, A.scores_to_long, A.summarize_scores):
        src = inspect.getsource(fn).lower()
        for token in ("wearableadapter", "capillaroscopyadapter", "mims", "nhanes", "stanford", "capillar"):
            assert token not in src, f"{fn.__name__} references {token!r}"


def test_contract_rejects_unnamespaced_or_pooled_ids(fake_registered):
    fake = A.get_adapter("synthetic_fake")
    scores = fake.phenotype_score(fake.embed(_synthetic_raw()))
    with pytest.raises(ValueError, match="namespaced"):
        A.scores_to_long(fake.describe_measurement(), scores, dataset_id="some_other_dataset")
    with pytest.raises(ValueError, match="missing"):
        A.scores_to_long(fake.describe_measurement(), scores.drop(columns=["synthetic_unmapped_phenotype"]),
                         dataset_id=SyntheticFakeAdapter.DATASET)


# ------------------------------------------------------------------------------------------ real data
@pytest.mark.data
def test_wearable_adapter_reproduces_nhanes_features():
    cmp, summ = W.validate_nhanes_reproduction(n_pass_per_cycle=6, n_fail_per_cycle=2, workers=1, n_short_per_cycle=2)
    assert summ["all_pass"], cmp[~cmp["pass"]]
    assert summ["n_features"] == len(NF.FEATURE_DEFINITIONS)
    # the short-record edge case (< 9 recorded days: different minute-grid length) is exercised too
    assert summ["n_short_record_checked"] == 4 and summ["n_short_record_embedded_passing"] == 4
    assert set(cmp["sample"]) == {"random_rule_passing", "all_rule_passing_with_fewer_than_9_recorded_days"}


@pytest.mark.data
def test_hr_adapter_reproduces_stanford_longitudinal_features():
    # stored daily table (wrapper check) + the adapter's own minute path on released Uwakwe minute files
    cmp, summ = W.validate_stanford_reproduction(native_ids=["0", "1", "2", "3", "4"])
    assert summ["all_pass"], cmp[~cmp["pass"]]
    raw = cmp[cmp["check"] == "uwakwe_raw_minutes"]
    assert len(raw) == len(W.LONGITUDINAL_COLUMNS) and (raw["n_participants"] == 5).all()


@pytest.mark.data
def test_reference_table():
    if not table_exists(W.REFERENCE_TABLE):
        pytest.skip("run `uv run python -m measure_it.wearables.adapter` first")
    ref = read_table(W.REFERENCE_TABLE)
    strata = ref[ref["feature_kind"] == "feature"].drop_duplicates(["sex", "age_band"])
    assert len(strata) == 2 * len(W.AGE_BANDS) + 1
    comp = ref[ref["feature"].isin({f for c in W.AXES.values() for f, _ in c}) | ref["feature"].str.startswith("axis:")]
    assert (comp["robust_scale"] > 0).all() and (comp["n"] >= W.MIN_STRATUM_N).all()
    assert set(ref["data_layer"]) == {"derived"} and (ref["weight_variable"] == "wtmec4yr_pooled").all()
    # sanity: activity volume declines with age in adults (reference medians)
    m = ref[(ref["feature"] == "mean_daily_mims") & (ref["sex"] == "female")].set_index("age_band")["w_median"]
    assert m["18-29"] > m["70-79"]


@pytest.mark.data
def test_adapter_scores_and_registry_tables():
    if not table_exists("participant_adapter_scores"):
        pytest.skip("run `uv run python -m measure_it.wearables.adapter` first")
    s = read_table("participant_adapter_scores")
    prefix = s["participant_id"].str.split(":", n=1).str[0]
    assert (prefix == s["dataset_id"]).all()
    assert set(s["adapter_status"]) == {"implemented"}
    assert not s.duplicated(["adapter_id", "participant_id", "phenotype"]).any()
    n = s[s["dataset_id"] == "nhanes"]
    scored = n[n["score_status"] == "scored"]
    # reference participants: median axis score close to 0 by construction (robust z against own strata)
    assert scored.groupby("phenotype")["score"].median().abs().max() < 0.25
    if table_exists("measurement_adapter_registry"):
        r = read_table("measurement_adapter_registry")
        assert r["object_id"].str.startswith("measurement:").all()
        assert (r.loc[r["adapter_name"] == "capillaroscopy", "status"] == "stub").all()
