"""Tests for the NHANES 2011-2014 extended laboratory layers (measure_it.ingestion.nhanes_labs_extended) and the
layered increment analysis (measure_it.wearables.nhanes_layers).

Unit tests need no data; @pytest.mark.data tests check processed outputs
(run `uv run python -m measure_it.ingestion.nhanes_labs_extended` and `... -m measure_it.wearables.nhanes_layers`).
"""
import numpy as np
import pandas as pd
import pytest

from measure_it.config import PROCESSED, TABLES
from measure_it.ingestion import nhanes_labs_extended as E
from measure_it.provenance import PROVENANCE_COLUMNS
from measure_it.validate import GEO_COLUMN_RE

# ------------------------------------------------------------------ unit tests


def test_component_of_handles_underscored_names():
    assert E.component_of("HEPB_S_G") == ("HEPB_S", "G")
    assert E.component_of("UHMS_G") == ("UHMS", "G")
    assert E.component_of("SSHCV_E") == ("SSHCV", "E")


@pytest.mark.parametrize("lc,analytes,labels,expected", [
    ("LBDBPBLC", ["LBXBPB", "LBDBPBSI"], {}, "LBXBPB"),                                   # name rule
    ("URDUA3LC", ["URXUAS", "URXUAS3", "URXUAS5"],
     {"URDUA3LC": "Urinary Arsenous acid comment code", "URXUAS3": "Urinary Arsenous acid (ug/L)",
      "URXUAS5": "Urinary Arsenic acid (ug/L)", "URXUAS": "Urinary arsenic, total (ug/L)"}, "URXUAS3"),  # label rule
    ("URDMEALC", ["URXNDEA", "URXNMEA"],
     {"URDMEALC": "NMEA Comment Code", "URXNMEA": "N-Nitrosoethylmethylamine (NMEA) (ng/L)",
      "URXNDEA": "N-Nitrosodiethylamine (NDEA) (ng/L)"}, "URXNMEA"),                      # abbreviation
    ("LBDMMALC", ["LBXMMASI"], {}, "LBXMMASI"),                                          # SI-only analyte
    ("SSBFOAL", ["SSNFOA", "SSBFOA"], {}, "SSBFOA"),                                      # surplus naming
    ("LBDXYZLC", ["LBXABC"], {"LBDXYZLC": "Something comment code", "LBXABC": "Other (ug/L)"}, None),
])
def test_match_lc(lc, analytes, labels, expected):
    assert E.match_lc(lc, analytes, labels) == expected


def test_positive_code_map_requires_both_directions():
    cb = {"V": {"codes": {"1": "Positive", "2": "Negative", "3": "Indeterminate", ".": "Missing"}},
          "W": {"codes": {"1": "Genotype 1a", "2": "Genotype 1b"}}}
    assert E.positive_code_map("V", cb) == {1.0: 1.0, 2.0: 0.0}
    assert E.positive_code_map("W", cb) == {}


def test_decide_rules():
    assert E.decide("BFRPOL_G", "BFRs - Pooled Samples", "Nov 2020", "u")[0] == "excluded"
    assert E.decide("SSEVD_G", "EV-D68", "Withdrawn", "")[0] == "excluded"
    assert E.decide("CBC_H", "CBC", "2016", "u")[0] == "core"
    d = E.decide("UHMS_G", "Metals - Urine - Special Sample", "2013", "u")
    assert d[0] == "ingested" and d[1] == "env_metals_urine" and "special" in d[2]
    assert E.decide("OMP_all_years", "Oral Microbiome Project", "Oct 2022", "")[1] == "oral_microbiome_16s"
    with pytest.raises(KeyError):
        E.decide("NEWTHING_H", "A new panel", "2030", "u")


def test_layer_family():
    assert E.layer_family("env_pfas_serum") == "environmental chemicals"
    assert E.layer_family("infectious_hpv") == "infectious serology / molecular"
    assert E.layer_family("oral_microbiome_16s") == "oral microbiome (16S)"


def test_bh_and_holm():
    from measure_it.wearables import nhanes_layers as NL
    p = np.array([0.01, 0.04, 0.03, 0.5])
    q = NL.bh(p)
    assert np.all(q >= p) and q[0] == pytest.approx(0.04) and q[3] == pytest.approx(0.5)
    h = NL.holm(p)
    assert h[0] == pytest.approx(0.04) and np.all(np.diff(h[np.argsort(p)]) >= 0)


def test_transform_continuous_logs_only_skewed_nonnegative():
    from measure_it.wearables import nhanes_layers as NL
    rng = np.random.default_rng(0)
    skewed = pd.Series(rng.lognormal(0, 1.5, 500))
    sym = pd.Series(rng.normal(0, 1, 500))
    assert NL._transform_continuous(skewed)[1].startswith("log")
    assert NL._transform_continuous(sym)[1] == "none"


# ------------------------------------------------------------------ data tests

needs_data = pytest.mark.data


def _t(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built")
    return pd.read_parquet(p)


@needs_data
def test_extended_table_contract():
    ext = _t("participant_labs_extended__nhanes")
    assert set(PROVENANCE_COLUMNS) <= set(ext.columns)
    assert ext["participant_id"].str.match(r"^nhanes:\d+$").all()
    assert (ext["data_layer"] == "person").all()
    assert not [c for c in ext.columns if GEO_COLUMN_RE.search(c) and c not in PROVENANCE_COLUMNS]
    assert not ext.duplicated(["participant_id", "component_file", "lab_variable"]).any()
    assert ext["weight_variable"].notna().all() and (ext["weight_variable"] != "").all()
    # below_lod is exactly lc_code == 1 where a comment code exists
    has = ext["lc_code"].notna()
    assert (ext.loc[has, "below_lod"].astype(bool) == (ext.loc[has, "lc_code"] == 1)).all()
    assert ext.loc[~has, "below_lod"].isna().all()
    assert ext["lab_layer"].nunique() >= 40
    # participants all exist in the core participants table (join on SEQN only)
    parts = _t("participants__nhanes")
    assert ext["participant_id"].isin(parts["participant_id"]).all()


@needs_data
def test_membership_matches_layer_table():
    m = _t("participant_layer_membership__nhanes")
    lt = pd.read_csv(TABLES / "nhanes_lab_layers.csv")
    assert m["participant_id"].is_unique and len(m) == 19931
    assert int(m["analysis_population"].sum()) == 8866
    for r in lt.itertuples():
        col = r.lab_layer if r.lab_layer in E.NONLAB_LAYERS else f"lab__{r.lab_layer}"
        assert int(m.loc[m["analysis_population"], col].sum()) == r.n_analysis_population
    ext = _t("participant_labs_extended__nhanes")
    a = ext[ext["analyte_role"] == "analyte"]
    n_pfas = a.loc[a["lab_layer"] == "env_pfas_serum", "participant_id"].nunique()
    assert int(m["lab__env_pfas_serum"].sum()) == n_pfas


@needs_data
def test_inventory_covers_every_listed_file():
    inv = pd.read_csv(TABLES / "nhanes_lab_inventory.csv")
    assert set(inv["decision"]) <= {"ingested", "core", "excluded"}
    assert (inv.loc[inv["decision"] == "excluded", "reason"].fillna("") != "").all()
    assert {"G", "H"} == set(inv["cycle_code"])
    assert not inv.loc[inv["file"].str.contains("POL_"), "decision"].eq("ingested").any()   # pooled never ingested


@needs_data
def test_overlap_table_consistent():
    ov = pd.read_csv(TABLES / "nhanes_layer_overlap.csv")
    lay = ov[ov["row_type"] == "layer"].set_index("set")
    assert lay.loc["wearable_valid", "n_analysis_population"] == 8866
    cum = ov[ov["row_type"] == "cumulative__analysis_population"]["n_analysis_population"].to_numpy()
    assert (np.diff(cum) <= 0).all()          # adding layers can only shrink the intersection
    fam = ov[ov["row_type"] == "family_combination__analysis_population"]
    assert fam["n_analysis_population"].sum() == 8866   # exact intersections partition the population


@needs_data
def test_layers_results_if_built():
    p = TABLES / "nhanes_layers_headline.csv"
    if not p.exists():
        pytest.skip("layer analysis not run")
    h = pd.read_csv(p)
    est = h[h["estimable"]]
    assert (est["n_cases"] >= 30).all()
    assert est["verdict"].notna().all()
    assert {"labs_full", "wear_over_full"} <= set(est["comparison"])
    s = pd.read_csv(TABLES / "nhanes_layers_screen.csv")
    e = s[s["estimable"]]
    assert e["bh_q"].between(0, 1).all()
    # restricted populations never exceed the target population of the full cohort
    assert (e["n"] <= 8866).all()
