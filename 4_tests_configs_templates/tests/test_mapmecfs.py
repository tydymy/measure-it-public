"""Tests for measure_it.ingestion.mapmecfs (NIH PI-ME/CFS deposits + linkage audit).

Unit tests exercise pure functions on synthetic inputs. Tests marked `data`
check the processed outputs of `uv run python -m measure_it.ingestion.mapmecfs`.
"""
from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd
import pytest
import yaml

from measure_it.config import PROCESSED, RAW, UNKNOWN
from measure_it.ingestion import mapmecfs as mm
from measure_it.provenance import check_provenance
from measure_it.registry import REQUIRED_FIELDS

SOURCE = mm.SOURCE_ID


# --------------------------------------------------------------------------- unit tests
@pytest.mark.parametrize("label,expected", [
    ("Peripheral Blood Mononuclear Cells, HV, Baseline MECFS_103", "103"),
    ("Vastus Lateralis Muscle, PI-ME/CFS, Baseline [MECFS_325]", "325"),
    ("Patient #312", "312"),
    ("patient #322", "322"),
    ("Control #113", "113"),
    ("GSM7989349_Walitt101", "101"),
    ("GSM7846603_FibroWalitt110", "110"),
    ("NT_FIBRO_WALITT_101", "101"),
    ("S101", "101"),
    ("103", "103"),
    (103.0, "103"),
])
def test_study_number_recognises_nih_formats(label, expected):
    assert mm.study_number(label) == expected


@pytest.mark.parametrize("label", [
    "SID_340", "SID_101", "map000288-00-02", "Fatigue-126", "HV1", "MECFS1", "HV A", "PI-ME/CFS A",
    "NINR-00418", "HV", "PI-ME/CFS", "", None, float("nan"), "250", "Patient #512",
])
def test_study_number_rejects_other_namespaces(label):
    # other namespaces must never be coerced into NIH study numbers
    assert mm.study_number(label) is None


def test_participant_id_is_namespaced():
    assert mm.participant_id("101") == "mapmecfs:MECFS_101"
    assert mm.native_id("325") == "MECFS_325"


@pytest.mark.parametrize("text,expected", [
    ("healthy volunteer (HV)", "HV"), ("Control", "HV"), ("Peripheral Blood Mononuclear Cells, HV, Baseline", "HV"),
    ("post-infectious Myalgic encephalomyelitis/chronic fatigue syndrome (PI-ME/CFS) volunteer", "PI-ME/CFS"),
    ("Patient", "PI-ME/CFS"), ("Vastus Lateralis Muscle, PI-ME/CFS, Baseline", "PI-ME/CFS"), ("", ""), ("Serum", ""),
])
def test_normalise_group(text, expected):
    assert mm.normalise_group(text) == expected


def test_blocks_splits_on_empty_header_cells():
    header = ["a", "b", "", "c", "d", "e", "", "", "f"]
    assert mm._blocks(header) == [(0, 1), (3, 5), (8, 8)]


def test_parse_adat_synthetic(tmp_path):
    rf = ["PlateId", "SampleId", "SampleType", "SampleDescription"]
    lines = ["!Checksum\tx", "^HEADER", "!Title\tT", "^COL_DATA", "^ROW_DATA", "^TABLE_BEGIN",
             "\t" * len(rf) + "SeqId\t1-1_1\t2-2_2",
             "\t" * len(rf) + "EntrezGeneSymbol\tAAA\tBBB",
             "\t" * len(rf) + "Target\tA\tB",
             "\t".join(rf) + "\t\t\t",
             "P1\tFatigue-101\tSample\t\t\t1.5\t2.5",
             "P1\tQC1\tQC\t\t\t9\t9",
             "P1\tFatigue-102\tSample\tSerum and CSF from separate visits\t\t3\t"]
    p = tmp_path / "x.adat.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    meta, colmeta, vals = mm.parse_adat(p)
    assert list(colmeta["SeqId"]) == ["1-1_1", "2-2_2"]
    assert list(meta["SampleId"]) == ["Fatigue-101", "QC1", "Fatigue-102"]
    assert vals.shape == (3, 2)
    assert vals[0, 1] == 2.5 and np.isnan(vals[2, 1])


def _mod(mid, ns, ids=(), sns=(), obtained=True, wearable=False, omics=False):
    return mm.Modality(mid, mid, "d", "a", "u", obtained, obtained, "f", "k", wearable, omics, "", ns,
                       native_ids=list(ids), study_numbers=list(sns), n_subjects=len(ids) or None,
                       blocked_reason="" if obtained else "login required")


def test_overlap_table_rules():
    a = _mod("pbmc", mm.STUDY_NUMBER_NS, ["MECFS_101", "MECFS_302"], ["101", "302"], omics=True)
    b = _mod("serum", "somalogic_fatigue_code", ["Fatigue-1", "Fatigue-2"], ["101", "303"], omics=True)
    c = _mod("stool", "microbiome_subject_sid", ["SID_101"], omics=True)
    d = _mod("hrv", "group_label_only", wearable=True)
    e = _mod("accel_gated", "not_obtained", obtained=False, wearable=True)
    t = mm.overlap_table([a, b, c, d, e]).set_index(["modality_a", "modality_b"])
    assert t.loc[("pbmc", "serum"), "comparable"] and t.loc[("pbmc", "serum"), "n_overlap"] == 1
    assert "crosswalk" in t.loc[("pbmc", "serum"), "reason"]
    # SID_101 must not be matched to study number 101; the overlap is unknown (NA), not zero
    assert not t.loc[("pbmc", "stool"), "comparable"] and pd.isna(t.loc[("pbmc", "stool"), "n_overlap"])
    assert not t.loc[("pbmc", "hrv"), "comparable"] and pd.isna(t.loc[("pbmc", "hrv"), "n_overlap"])
    assert not t.loc[("pbmc", "accel_gated"), "comparable"] and pd.isna(t.loc[("pbmc", "accel_gated"), "n_overlap"])


def test_adj_p_consistency_detects_non_adjustments():
    p = [0.001, 0.01, 0.02, 0.04]
    assert mm.adj_p_consistent_with_p(p, [0.004, 0.02, 0.0267, 0.04]) is True          # BH-like
    assert mm.adj_p_consistent_with_p(p, [0.0005, 0.02, 0.03, 0.04]) is False          # adj < p
    assert mm.adj_p_consistent_with_p(p, [0.03, 0.02, 0.035, 0.04]) is False           # not monotone in p
    assert mm.adj_p_consistent_with_p(p, [np.nan] * 4) is None


def test_somalogic_value_scale_against_adat():
    rng = np.random.default_rng(0)
    ref = rng.uniform(50, 5000, size=(20, 3))                       # 20 samples x 3 SOMAmers (RFU)
    seq = ["a", "b", "c"]
    raw = pd.DataFrame(ref[:10].T)                                  # table rows = SOMAmers, columns = samples
    logged = np.log2(raw)
    assert mm.somalogic_value_scale(seq, raw, seq, ref)[0] == "RFU"
    assert mm.somalogic_value_scale(seq, logged, seq, ref)[0] == "log2_RFU"
    assert mm.somalogic_value_scale(seq, raw * 3.7, seq, ref)[0] == "unknown"
    assert mm.somalogic_value_scale(seq, logged, None, None)[0] == "log2_RFU"     # heuristic without ADAT


def test_birth_sex_value_alignment():
    # androgen high in Birth Sex 0 rows, low in Birth Sex 1 rows, in both blocks
    df = pd.DataFrame({"ID": ["HV", "101", "HV", "102", "PI-ME/CFS", "301", "PI-ME/CFS", "302"],
                       "Birth Sex": [0, 1, 0, 1, 0, 1, 0, 1],
                       "androsterone sulfate": [2.0, 0.3, 1.8, 0.2, 2.2, 0.25, 1.5, 0.3]})
    df["study_number"] = df["ID"].map(mm.study_number)
    al = mm.birth_sex_value_alignment(df)
    assert al["HV"]["auc"] == 1.0 and al["PI-ME/CFS"]["auc"] == 1.0
    shifted = df.assign(**{"androsterone sulfate": np.roll(df["androsterone sulfate"].to_numpy(), 1)})
    assert mm.birth_sex_value_alignment(shifted)["HV"]["auc"] < 0.5


def test_scan_github_r_objects(tmp_path, monkeypatch):
    monkeypatch.setattr(mm, "GITHUB_DIR", tmp_path)
    (tmp_path / "x.RData").write_bytes(gzip.compress(b"RDX3\nX\x00Fatigue-101\x00Fatigue-102\x00SID_7\x00HV"))
    (tmp_path / "y.rds").write_bytes(gzip.compress(b"X\n\x00\x00MECFS_301\x00"))
    rows = {r["file"]: r["identifier_types"] for r in mm.scan_github_r_objects()}
    assert rows["github: x.RData"] == "somalogic_fatigue_code (2), microbiome_subject_sid (1)"
    assert rows["github: y.rds"] == f"{mm.STUDY_NUMBER_NS} (1)"


def test_render_linkage_audit_works_with_nothing_downloaded():
    """The pipeline must produce an audit even when no deposit could be fetched."""
    ctx = {"mods": [], "ovl": pd.DataFrame(), "parts": pd.DataFrame(), "log": [
        mm.AccessRecord("mapMECFS", "https://www.mapmecfs.org/", "GET", 403, "login_required", "Login Required.")],
        "checks": [], "written": {}, "source_scan": []}
    text = mm.render_linkage_audit(ctx)
    assert "No participant-level identifiers could be obtained" in text
    assert "## 6. Can wearable/actigraphy and omics be linked?" in text
    assert "Login Required" in text


def test_run_succeeds_when_nothing_is_downloadable(tmp_path, monkeypatch):
    """run() must still write the audit docs and a registry entry when no deposit is available.

    Every write target is redirected into tmp_path so real outputs are untouched.
    """
    raw = tmp_path / "raw"
    raw.mkdir()
    for name, sub in (("RAW", ""), ("EXTRACT", "extracted"), ("SOURCE_DATA_DIR", "extracted/source_data/x"),
                      ("SUPP_DIR", "extracted/supp_data/x"), ("GITHUB_DIR", "extracted/github/x")):
        monkeypatch.setattr(mm, name, raw / sub if sub else raw)
    monkeypatch.setattr(mm, "LINKAGE_AUDIT_PATH", tmp_path / "MAP_MECFS_LINKAGE_AUDIT.md")
    monkeypatch.setattr(mm, "PROCESSED", tmp_path / "processed")
    monkeypatch.setattr(mm, "load_manifest", lambda sid: {"source_id": sid, "files": {}})
    tables, registry = {}, {}
    monkeypatch.setattr(mm, "write_table", lambda df, name, **kw: tables.setdefault(name, df))
    monkeypatch.setattr(mm, "write_registry_entry", lambda e: registry.update(e))

    ctx = mm.run(fetch=False)

    assert ctx["status"] == "blocked"
    assert "participants__mapmecfs" not in tables and "condition_molecular_evidence__mapmecfs" not in tables
    assert not [f for f in REQUIRED_FIELDS if f not in registry]
    assert registry["wearable_omics_linkage"] is False
    assert registry["true_participant_linkage_across_modalities"] is False
    text = (tmp_path / "MAP_MECFS_LINKAGE_AUDIT.md").read_text()
    assert "No participant-level identifiers could be obtained" in text
    assert (raw / "DATA_AUDIT.md").exists()


# --------------------------------------------------------------------------- data tests
def _read(name):
    p = PROCESSED / f"{name}.parquet"
    if not p.exists():
        pytest.skip(f"{name} not built; run `uv run python -m measure_it.ingestion.mapmecfs`")
    return pd.read_parquet(p)


def _registry():
    p = RAW / SOURCE / "registry_entry.yaml"
    if not p.exists():
        pytest.skip("registry entry not written")
    return yaml.safe_load(p.read_text())


@pytest.mark.data
def test_registry_entry_complete_and_honest():
    e = _registry()
    assert not [f for f in REQUIRED_FIELDS if f not in e]
    assert e["status"] in {"partial", "blocked"}
    assert e["wearable_omics_linkage"] is False
    assert e["wearable"] is False
    assert "login" in e["access_conditions"].lower()


@pytest.mark.data
def test_participants_table():
    p = _read("participants__mapmecfs")
    check_provenance(p, "participants__mapmecfs")
    assert p["participant_id"].is_unique
    assert p["participant_id"].str.fullmatch(r"mapmecfs:MECFS_[13]\d{2}").all()
    assert (p["data_layer"] == "person").all()
    assert not p["has_wearable_or_hrv_data_open"].any()
    # prefix rule 1xx HV / 3xx PI-ME/CFS agrees with every deposited group label
    known = p[p["group_concordant_across_sources"] & (p["group_label"] != UNKNOWN)]
    assert len(known) > 0
    assert (known["group_label"] == known["group_from_study_number_prefix"]).all()
    assert (p["n_open_omics_modalities"] >= 1).all()


@pytest.mark.data
def test_linked_omics_only_joins_on_published_ids():
    o = _read("participant_omics_linked__mapmecfs")
    p = _read("participants__mapmecfs")
    check_provenance(o.head(1000), "participant_omics_linked__mapmecfs")
    assert set(o["participant_id"]) == set(p["participant_id"])
    assert set(o["linkage_method"].astype(str)) <= {"direct_study_number", "crosswalk_fatigue_code_via_GSE251790"}
    assert set(o["modality"].astype(str)) <= {"pbmc_rnaseq", "muscle_rnaseq", "csf_somalogic", "serum_somalogic",
                                               "csf_metabolomics"}
    # one sample per participant per modality
    per = o.groupby(["participant_id", "modality"], observed=True)["sample_accession"].nunique()
    assert (per == 1).all()
    assert o["source_record_id"].is_unique
    # every participant flag matches presence in the omics table
    for m in ("pbmc_rnaseq", "muscle_rnaseq", "csf_somalogic", "serum_somalogic", "csf_metabolomics"):
        have = set(o.loc[o["modality"] == m, "participant_id"])
        assert have == set(p.loc[p[f"has_{m}"], "participant_id"])


@pytest.mark.data
def test_no_wearable_partition_when_not_linked():
    e = _registry()
    if not e["wearable_omics_linkage"]:
        assert not (PROCESSED / "participant_wearable_features__mapmecfs.parquet").exists()


@pytest.mark.data
def test_overlap_table_consistent_with_inventory():
    inv = _read("modality_inventory__mapmecfs")
    ov = _read("modality_overlap__mapmecfs")
    check_provenance(inv, "modality_inventory__mapmecfs")
    check_provenance(ov, "modality_overlap__mapmecfs")
    n = len(inv)
    assert len(ov) == n * (n - 1) // 2
    # no wearable x omics pair is comparable, and no count is reported for non-comparable pairs
    wo = ov[(ov["a_is_wearable"] & ov["b_is_omics"]) | (ov["b_is_wearable"] & ov["a_is_omics"])]
    assert len(wo) > 0 and not wo["comparable"].any()
    assert ov.loc[~ov["comparable"], "n_overlap"].isna().all()
    # comparable pairs never exceed either side
    c = ov[ov["comparable"]]
    assert (c["n_overlap"] <= c[["n_a", "n_b"]].min(axis=1)).all()
    # overlap ids listed == count
    assert (c["overlap_ids"].fillna("").map(lambda s: len([x for x in s.split(",") if x])) == c["n_overlap"]).all()
    # gated mapMECFS files are recorded as not obtained
    gated = inv[inv["id_namespace"] == "not_obtained"]
    assert (~gated["obtained"]).all() and (gated["blocked_reason"].str.contains("Login Required")).all()


@pytest.mark.data
def test_condition_evidence():
    ev = _read("condition_molecular_evidence__mapmecfs")
    check_provenance(ev, "condition_molecular_evidence__mapmecfs")
    assert (ev["data_layer"] == "condition_molecular").all()
    assert (ev["evidence_type"] == "published_biomarker").all()
    assert (ev["condition_id"] == "me_cfs").all()
    assert (ev["pmid"] == "38383456").all() and (ev["doi"] == "10.1038/s41467-024-45107-3").all()
    assert ev["p_value"].between(0, 1).all()
    assert set(ev["direction"]) <= {"up_in_PI-ME/CFS", "down_in_PI-ME/CFS", "no_change",
                                    "no_difference_in_medians", UNKNOWN}
    assert "participant_id" not in ev.columns
    # direction of RNA-seq rows must follow the sign of the published logFC
    rn = ev[ev["effect_measure"].str.startswith("limma")]
    assert len(rn) > 0
    assert ((rn["effect_value"] > 0) == (rn["direction"] == "up_in_PI-ME/CFS")).all()
    # sample sizes present for every row
    assert ev["n_cases"].notna().all() and ev["n_controls"].notna().all()
    # an adjusted-p column that is not an adjustment of p is never used as FDR evidence (SD14A)
    bad = ev[ev["adj_p_consistent_with_p"] == False]  # noqa: E712
    assert bad["fdr_significant"].isna().all()
    assert not (ev["evidence_level"] == "single_cohort_fdr_significant")[ev["fdr_significant"].isna()].any()
    assert (ev.loc[ev["table_id"] == "14A", "table_scope"] == "significant_results_only").all()
    # SD17 effects respect the published value scale (17B is log2 RFU: effect = difference of medians)
    for tid, scale in (("17A", "RFU"), ("17B", "log2_RFU")):
        t = ev[ev["table_id"] == tid]
        if t.empty:
            continue
        assert (t["value_scale"] == scale).all()
        want = (t["median_pi_mecfs"] - t["median_hv"]) if scale == "log2_RFU" else \
            np.log2(t["median_pi_mecfs"] / t["median_hv"])
        assert np.allclose(t["effect_value"], want, equal_nan=True)
    # ontology id comes from the condition registry, never typed in
    assert ev["condition_ontology_id"].nunique() == 1
    assert ev["condition_ontology_id"].iloc[0] == UNKNOWN or ev["condition_ontology_id"].iloc[0].startswith("MONDO:")


@pytest.mark.data
def test_linkage_audit_document():
    p = mm.LINKAGE_AUDIT_PATH
    if not p.exists():
        pytest.skip("linkage audit not written")
    text = p.read_text()
    for section in ("## Bottom line", "## 1. Public deposits", "## 2. Participant identifier schemes",
                    "## 3. Subjects per modality", "## 4. Exact identifier overlap", "## 6. Can wearable/actigraphy",
                    "## 8. What is blocked and why", "## 10. Other ME/CFS datasets"):
        assert section in text
    assert "cannot truthfully be linked" in text
    assert (RAW / SOURCE / "DATA_AUDIT.md").exists()
    log = json.loads((RAW / SOURCE / "access_log.json").read_text())
    assert not any("X-Amz-Security-Token" in (r.get("url") or "") for r in log)
