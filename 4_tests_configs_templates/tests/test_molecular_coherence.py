"""Tests for measure_it.omics.coherence (Test 3) and measure_it.omics.graph (measurable_biology, Fig 5).

Unit tests use synthetic inputs. Tests marked `data` check the outputs of
`uv run python -m measure_it.omics.coherence && uv run python -m measure_it.omics.graph`.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from measure_it.config import PROCESSED, TABLES, UNKNOWN, load_config
from measure_it.omics import coherence as C
from measure_it.omics import graph as G
from measure_it.provenance import EVIDENCE_TYPES, PROVENANCE_COLUMNS, check_provenance


# =========================================================================== unit: statistics
def test_bh_fdr_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = np.array([0.001, 0.02, 0.03, 0.04, 0.5, 0.9, np.nan])
    got = C.bh_fdr(p)
    want = multipletests(p[:-1], method="fdr_bh")[1]
    assert np.allclose(got[:-1], want) and np.isnan(got[-1])


def test_hypergeom_p_matches_scipy_and_edge_cases():
    assert C.hypergeom_p(0, 100, 10, 10) == 1.0
    assert math.isclose(C.hypergeom_p(3, 100, 10, 10), stats.hypergeom.sf(2, 100, 10, 10))
    # complete overlap of two 5-gene sets out of 20000 is extreme
    assert C.hypergeom_p(5, 20000, 5, 5) < 1e-15


def test_overlap_stats():
    a, b = {"A", "B", "C", "D"}, {"C", "D", "E"}
    s = C.overlap_stats(a, b, 100)
    assert s["overlap"] == 2 and s["jaccard"] == pytest.approx(2 / 5)
    assert s["expected_overlap"] == pytest.approx(4 * 3 / 100)
    assert s["log2_fold"] == pytest.approx(np.log2(2.5 / (0.12 + 0.5)))
    assert s["hypergeom_p"] == pytest.approx(stats.hypergeom.sf(1, 100, 4, 3))


def test_label_permutation_exact_and_monte_carlo():
    S = np.full((4, 4), 0.0)
    np.fill_diagonal(S, 1.0)                     # perfect same-condition agreement
    r = C.label_permutation_test(S)
    assert r["exact"] and r["n_perm"] == 24
    assert r["p"] == pytest.approx(1 / 24)       # only the identity reaches T_obs
    rng = np.random.default_rng(0)
    N = rng.normal(size=(9, 9))                  # 9! > 50,000 -> Monte Carlo
    r2 = C.label_permutation_test(N, n_perm=2000)
    assert not r2["exact"] and r2["n_perm"] == 2000 and 0 < r2["p"] <= 1
    assert np.isnan(C.label_permutation_test(np.ones((1, 1)))["p"])


def test_same_vs_cross_empirical_p():
    S = pd.DataFrame([[5.0, 1.0, 0.0], [1.0, 4.0, 2.0], [0.0, 2.0, 0.5]], index=list("abc"), columns=list("abc"))
    per = C.same_vs_cross(S).set_index("condition_id")
    assert per.loc["a", "n_cross"] == 4 and per.loc["a", "empirical_p"] == pytest.approx(1 / 5)
    assert per.loc["c", "same"] == 0.5 and per.loc["c", "empirical_p"] > 0.5


def test_bootstrap_ci_brackets_mean():
    lo, hi = C.bootstrap_mean_ci([1.0, 2.0, 3.0, 4.0, 5.0])
    assert lo < 3.0 < hi
    assert all(np.isnan(C.bootstrap_mean_ci([1.0])))


def test_ora_recovers_planted_pathway():
    bg = frozenset(f"G{i}" for i in range(1000))
    lib = {"P1": frozenset(f"G{i}" for i in range(20)), "P2": frozenset(f"G{i}" for i in range(500, 540)),
           "TOO_SMALL": frozenset({"G1", "G2"})}
    res = C.ora({f"G{i}" for i in range(10)} | {"G700", "NOT_IN_BG"}, lib, bg).set_index("reactome_id")
    assert "TOO_SMALL" not in res.index
    assert res.loc["P1", "overlap"] == 10 and res.loc["P1", "enriched"]
    assert res.loc["P2", "overlap"] == 0 and not res.loc["P2", "enriched"]
    assert res.loc["P1", "set_size_in_background"] == 11      # NOT_IN_BG dropped


# =========================================================================== unit: curated tables + hierarchy
def test_plan_tables_parse_and_use_known_vocabularies():
    anchors = C.system_anchor_table()
    assert set(anchors["level"]) <= {"top", "sub"}
    assert anchors["reactome_id"].str.match(r"^R-HSA-\d+$").all()
    assert set(anchors["physiological_system"]) <= set(C.SYSTEMS)
    cmap = C.system_measurement_map()
    known = {m["id"] for m in load_config("measurements")["measurement_classes"]}
    assert set(cmap["measurement_class"]) <= known
    assert (cmap["rationale"].str.len() > 20).all()          # every curated row carries a rationale
    assert "other" not in set(cmap["physiological_system"])  # 'other' maps to no measurement


def test_parse_markdown_table_errors():
    with pytest.raises(ValueError):
        C.parse_markdown_table("no table here", "X")
    txt = "<!-- X:BEGIN -->\n| a | b |\n|---|---|\n| 1 | 2 |\n<!-- X:END -->"
    assert C.parse_markdown_table(txt, "X").to_dict("records") == [{"a": "1", "b": "2"}]


def test_classify_pathways_anchor_rule():
    names = {"T1": "Immune System", "T2": "Signal Transduction", "T3": "Metabolism", "A": "child of immune",
             "B": "adrenoceptors", "B1": "child of adrenoceptors", "C": "two parents", "D": "catecholamine"}
    parents = {"A": ["T1"], "B": ["T2"], "B1": ["B"], "C": ["T1", "T3"], "D": ["T3"]}
    tops = {"T1", "T2", "T3"}
    anchors = pd.DataFrame([
        {"reactome_id": "T1", "level": "top", "physiological_system": "immune"},
        {"reactome_id": "T3", "level": "top", "physiological_system": "mitochondrial_metabolic"},
        {"reactome_id": "B", "level": "sub", "physiological_system": "autonomic_cardiac"},
        {"reactome_id": "D", "level": "sub", "physiological_system": "autonomic_cardiac"}])
    cl = C.classify_pathways(names, parents, tops, anchors).set_index("reactome_id")["physiological_systems"]
    assert cl["A"] == "immune"
    assert cl["B1"] == "autonomic_cardiac"          # sub-anchor inherited through the hierarchy
    assert cl["T2"] == "other"                      # top level without an anchor
    assert cl["C"] == "immune;mitochondrial_metabolic"   # several parents -> several systems
    assert cl["D"] == "autonomic_cardiac"           # sub-anchor overrides its Metabolism top level
    top_only = C.classify_pathways(names, parents, tops, anchors, use_sub_anchors=False).set_index("reactome_id")
    assert top_only.loc["D", "physiological_systems"] == "mitochondrial_metabolic"
    assert top_only.loc["B1", "physiological_systems"] == "other"


def _hgnc():
    return pd.DataFrame([
        {"hgnc_id": "HGNC:1", "symbol": "MARCHF3", "locus_group": "protein-coding gene", "status": "Approved",
         "prev_symbol": "MARCH3", "alias_symbol": np.nan, "ensembl_gene_id": "ENSG00000173926"},
        {"hgnc_id": "HGNC:2", "symbol": "C18orf21", "locus_group": "protein-coding gene", "status": "Approved",
         "prev_symbol": np.nan, "alias_symbol": "SHARED", "ensembl_gene_id": "ENSG2"},
        {"hgnc_id": "HGNC:3", "symbol": "GENE3", "locus_group": "protein-coding gene", "status": "Approved",
         "prev_symbol": np.nan, "alias_symbol": "SHARED", "ensembl_gene_id": "ENSG3"},
        {"hgnc_id": "HGNC:4", "symbol": "LINC1", "locus_group": "non-coding RNA", "status": "Approved",
         "prev_symbol": np.nan, "alias_symbol": np.nan, "ensembl_gene_id": "ENSG4"},
    ])


def test_symbol_mapper_rules():
    m = C.SymbolMapper(_hgnc())
    assert m.map_symbol("C18ORF21") == ("C18orf21", "approved")          # case-insensitive
    assert m.map_symbol("MARCH3") == ("MARCHF3", "previous_symbol")
    assert m.map_symbol("SHARED") == (None, "ambiguous_alias")            # never guessed
    assert m.map_symbol("2022-03-03 00:00:00") == (None, "date_mangled")  # spreadsheet date, not a gene
    assert m.map_symbol("NOPE") == (None, "unmapped")
    assert m.map_symbol(np.nan) == (None, "empty")
    assert m.map_ensembl("ENSG00000173926.5") == ("MARCHF3", "ensembl")   # version stripped
    assert m.map_ensembl("P12345") == (None, "no_ensembl")
    assert m.is_protein_coding("MARCHF3") and not m.is_protein_coding("LINC1")


def test_gene_sets_scope_rules():
    """Open Targets primary scope drops broad ids and indirect-only targets of ids with descendants."""
    m = C.SymbolMapper(_hgnc())
    base = {"source_database": C.OT_SOURCE_DB, "source_evidence_category": "overall_association",
            "ot_datatype_score__genetic_association": 0.5, "ot_datatype_score__literature": np.nan,
            "source_name": "ot", "ontology_id": "MONDO:1"}
    ev = pd.DataFrame([
        {**base, "condition_id": "x", "entity_id": "ENSG00000173926", "entity_label": "MARCHF3",
         "ontology_match": "exact", "n_descendants_in_source": 0, "has_direct_association": True},
        {**base, "condition_id": "x", "entity_id": "ENSG2", "entity_label": "C18orf21",
         "ontology_match": "exact", "n_descendants_in_source": 5, "has_direct_association": False},  # indirect only
        {**base, "condition_id": "y", "entity_id": "ENSG3", "entity_label": "GENE3",
         "ontology_match": "broad", "n_descendants_in_source": 0, "has_direct_association": True},   # broad id
        {**base, "condition_id": "x", "entity_id": "ENSG4", "entity_label": "LINC1",
         "ontology_match": "exact", "n_descendants_in_source": 0, "has_direct_association": True},   # non-coding
    ])
    for col in ("entity_type", "mapped_genes", "p_value", "reported_trait", "source_accession", "table_id",
                "analyte_class", "gene_symbol", "ensembl_gene_id"):
        ev[col] = np.nan
    gs, _ = C.build_gene_sets(ev, m, C.PRIMARY)
    assert set(gs.loc[gs["source"] == "ot_genetic", "gene_symbol"]) == {"MARCHF3"}
    gs2, _ = C.build_gene_sets(ev, m, C.Settings("s", ot_scope="indirect"))
    assert set(gs2.loc[gs2["source"] == "ot_genetic", "gene_symbol"]) == {"MARCHF3", "C18orf21", "GENE3"}


def test_expected_systems_from_flags():
    reg = pd.DataFrame([{"canonical_condition_id": "a", "infectious_or_post_infectious": True, "autoimmune": False,
                         "autonomic": True, "vascular": False, "neurologic": False, "fatigue_pem": True}])
    assert C.expected_systems(reg)["a"] == {"immune", "autonomic_cardiac", "mitochondrial_metabolic"}


def test_support_tier():
    assert G.support_tier([], 0, 0) == "no_molecular_data"
    assert G.support_tier([], 0, 3) == "too_few_genes"
    assert G.support_tier([], 2, 50) == "not_supported"
    assert G.support_tier(["ot_literature"], 3, 50) == "literature_only"
    assert G.support_tier(["ot_genetic", "gwas_catalog"], 3, 50) == "single_source"   # one 'genetic' family
    assert G.support_tier(["ot_genetic", "ot_literature"], 3, 50) == "multi_source"


# =========================================================================== data tests
def _csv(name):
    p = TABLES / f"{name}.csv"
    if not p.exists():
        pytest.skip(f"{name} not built; run `uv run python -m measure_it.omics.coherence`")
    return pd.read_csv(p)


def _mb():
    p = PROCESSED / "measurable_biology.parquet"
    if not p.exists():
        pytest.skip("measurable_biology not built; run `uv run python -m measure_it.omics.graph`")
    return pd.read_parquet(p)


@pytest.mark.data
@pytest.mark.parametrize("name", ["test3_gene_sets", "test3_source_counts", "test3_gene_agreement_pairs",
                                  "test3_gene_agreement_control", "test3_reactome_ora", "test3_system_enrichment",
                                  "test3_flag_concordance", "test3_disagreements", "test3_sensitivity",
                                  "fig5_nodes", "fig5_edges"])
def test_result_tables_carry_provenance(name):
    df = _csv(name)
    check_provenance(df, name)
    assert set(df["data_layer"]) == {"derived"}
    assert set(df["evidence_type"]) <= EVIDENCE_TYPES
    assert not [c for c in df.columns if "participant" in c.lower() or "omics_score" in c.lower()]


@pytest.mark.data
def test_gene_sets_are_hgnc_protein_coding_and_scoped():
    gs = _csv("test3_gene_sets")
    assert set(gs["source"]) <= set(C.SOURCES)
    assert not gs.duplicated(["condition_id", "source", "gene_symbol"]).any()
    assert set(gs.loc[gs["source"] == "mapmecfs", "condition_id"]) == {"me_cfs"}
    # the grouping node has no primary-scope evidence and MCAS has none at all
    assert "post_infectious_syndrome" not in set(gs["condition_id"]) and "mcas" not in set(gs["condition_id"])
    has_ens = gs["ensembl_gene_id"].fillna("") != ""
    assert (gs.loc[has_ens, "object_id"] == "gene:" + gs.loc[has_ens, "ensembl_gene_id"]).all()


@pytest.mark.data
def test_source_counts_cover_registry_and_report_known_gaps():
    counts = _csv("test3_source_counts").set_index("condition_id")
    from measure_it.store import read_table
    assert set(counts.index) == set(read_table("condition_registry")["canonical_condition_id"])
    assert counts.loc["mcas", "n_gene_sources_nonempty"] == 0
    # fibromyalgia: GWAS Catalog genes but no Open Targets genetic evidence (reported as a disagreement)
    assert counts.loc["fibromyalgia", "genes_gwas_catalog"] > 0 and counts.loc["fibromyalgia", "genes_ot_genetic"] == 0
    dis = _csv("test3_disagreements")
    fib = dis[(dis["condition_id"] == "fibromyalgia") & (dis["kind"] == "presence_absence")
              & (dis["source_a"] == "ot_genetic") & (dis["source_b"] == "gwas_catalog")]
    assert len(fib) == 1
    # no FDR-significant mapMECFS gene-level row exists; the table must not claim otherwise
    assert counts.loc["me_cfs", "mapmecfs_gene_level_rows_fdr_significant"] == 0


@pytest.mark.data
def test_agreement_tables_consistent():
    pairs = _csv("test3_gene_agreement_pairs")
    assert (pairs["overlap"] <= pairs[["n_a", "n_b"]].min(axis=1)).all()
    assert pairs["hypergeom_p"].between(0, 1).all()
    assert pairs["same_condition"].any() and (~pairs["same_condition"]).any()   # control rows present
    ctrl = _csv("test3_gene_agreement_control")
    ok = ctrl[ctrl["status"] == "ok"]
    assert len(ok) > 0 and ok["perm_p"].between(0, 1).all()
    assert (ok["n_conditions"] >= 2).all()
    # every pair lacking >= 2 shared conditions is UNKNOWN, never a p-value
    assert ctrl.loc[ctrl["status"] != "ok", "perm_p"].isna().all()
    assert ctrl.loc[ctrl["status"] != "ok", "status"].str.startswith(UNKNOWN).all()


@pytest.mark.data
def test_system_enrichment_rules():
    se = _csv("test3_system_enrichment")
    assert set(se["physiological_system"]) == set(C.SYSTEMS)
    sup = se[se["supported"]]
    assert (sup["fdr"] < C.FDR_ALPHA).all() and (sup["n_enriched_pathways"] >= 1).all()
    assert (se["set_size_in_background"] >= C.MIN_ORA_GENES).all()
    ora = _csv("test3_reactome_ora")
    enr = ora[ora["enriched"]]
    assert (enr["fdr"] < C.FDR_ALPHA).all() and (enr["overlap"] >= C.MIN_OVERLAP).all()
    assert ora["pathway_size"].between(C.PW_MIN, C.PW_MAX).all()


@pytest.mark.data
def test_measurable_biology_table():
    mb = _mb()
    check_provenance(mb, "measurable_biology")
    assert mb["object_id"].is_unique and mb["object_id"].str.startswith("molbio:").all()
    known = {m["id"] for m in load_config("measurements")["measurement_classes"]}
    assert set(mb["measurement_class"]) <= known
    assert set(mb["data_layer"]) == {"derived"}
    sysrows = mb[mb["physiological_system"] != G.ASSAY_PLATFORM]
    assert set(sysrows["physiological_system"]) <= set(C.SYSTEMS) - {"other"}
    assert set(mb["support_status"]) <= {"multi_source", "single_source", "literature_only", "not_supported",
                                         "no_molecular_data", "too_few_genes", "study_catalogue_evidence"}
    # every system row says which link is data-derived and which curated
    assert sysrows["system_measurement_link"].str.startswith("curated").all()
    has = ~sysrows["support_status"].isin(["no_molecular_data", "too_few_genes"])
    assert sysrows.loc[has, "condition_system_link"].str.startswith("data_derived").all()
    assert sysrows.loc[~has, "condition_system_link"].str.startswith(UNKNOWN).all()
    assert (sysrows.loc[~has, "evidence_type"] == "curated_config").all()
    # MCAS: every row explicit, none supported
    mc = sysrows[sysrows["condition_id"] == "mcas"]
    assert len(mc) > 0 and (mc["support_status"] == "no_molecular_data").all()
    # supported rows agree with the Test 3 system table
    se = _csv("test3_system_enrichment")
    s = sysrows[sysrows["n_sources_supporting"] > 0].drop_duplicates(["condition_id", "physiological_system"])
    for r in s.itertuples():
        t = se[(se["condition_id"] == r.condition_id) & (se["physiological_system"] == r.physiological_system)
               & se["supported"]]
        assert set(t["source"]) == set(r.sources_supporting.split(";"))


@pytest.mark.data
def test_fig5_graph_is_closed_and_typed():
    nodes, edges = _csv("fig5_nodes"), _csv("fig5_edges")
    ids = set(nodes["node_id"])
    assert nodes["node_id"].is_unique
    assert set(edges["source_node"]) <= ids and set(edges["target_node"]) <= ids
    assert set(nodes["node_type"]) <= {"condition", "gene", "pathway", "physiological_system", "measurement_class"}
    cur = edges[edges["edge_type"] == "system_measurement"]
    assert len(cur) > 0 and (cur["link_type"] == "curated").all()
    assert (edges.loc[edges["edge_type"].isin(["condition_gene", "condition_pathway", "condition_system"]),
                      "link_type"] == "data_derived").all()
    for ext in ("png", "svg"):
        assert (G.DRAFTS / f"fig5_system_heatmap.{ext}").exists()
        assert (G.DRAFTS / f"fig5_graph_demo.{ext}").exists()


@pytest.mark.data
def test_mapmecfs_partition_has_common_columns():
    p = PROCESSED / "condition_molecular_evidence__mapmecfs.parquet"
    if not p.exists():
        pytest.skip("mapmecfs partition not built")
    ev = pd.read_parquet(p)
    for c in ("entity_type", "entity_id", "entity_label", "ontology_id", "source_database", "source_accession"):
        assert c in ev.columns, c
    from measure_it.store import read_table
    reg = read_table("condition_registry").set_index("canonical_condition_id")
    assert (ev["ontology_id"] == reg.loc["me_cfs", "primary_mondo_id"]).all()
    assert set(ev["entity_type"]) <= {"gene", "analyte"}


@pytest.mark.data
def test_query_functions_return_unknown_not_guesses():
    if not (TABLES / "test3_source_counts.csv").exists() or not (PROCESSED / "measurable_biology.parquet").exists():
        pytest.skip("outputs not built")
    r = C.get_molecular_coherence("mast cell activation syndrome")
    assert r["status"] == UNKNOWN and r.get("condition_id") == "mcas"
    r = C.get_molecular_coherence("zzz not a condition qqq")
    assert r["status"] == UNKNOWN
    r = C.get_molecular_coherence("ME/CFS")
    assert r["status"] == "available" and r["condition_id"] == "me_cfs"
    assert "not measured on any person" in " ".join(r["caveats"])
    m = G.get_measurable_biology("mcas")
    assert m["status"] == UNKNOWN
    m = G.get_measurable_biology("long_covid", "wearable autonomic monitoring")
    assert m["resolved_classes"] and set(m["resolved_classes"]) >= {"hrv", "wearable_heart_rate"}
    assert m["status"] in ("available", "not_supported", UNKNOWN)
    json.dumps(m, default=str)   # JSON-serialisable for the MCP layer
    # a tested negative is not reported as missing data
    mb = _mb()
    lc = mb[(mb["condition_id"] == "long_covid") & (mb["physiological_system"] == "autonomic_cardiac")]
    if len(lc) and (lc["support_status"] == "not_supported").all():
        assert m["status"] == "not_supported"
    # a condition with only study-catalogue rows says its system layer is UNKNOWN
    p = G.get_measurable_biology("ptlds")
    assert p["system_level"]["status"] == UNKNOWN or p["status"] == UNKNOWN
    r = C.get_molecular_coherence("migraine")
    for a in r["gene_agreement"]:
        if a["empirical_p_vs_cross_condition"] != UNKNOWN:
            assert a["min_attainable_empirical_p"] <= a["empirical_p_vs_cross_condition"] + 1e-12


def test_independent_loci_merges_nearby_positions():
    # HLA-DQA1/DQB1 (6:32.66 Mb) and HLA-DPB1 (6:33.09 Mb) are one locus at a 1 Mb window
    assert C.independent_loci(["6:32658495", "6:33087546"]) == 1
    assert C.independent_loci(["6:32658495", "6:35087546", "10:52220867"]) == 3
    assert C.independent_loci(["bad", "1:x"]) == 0


def test_overlap_diagnostics():
    info = pd.DataFrame([
        {"condition_id": "c", "source": "ot_other", "gene_symbol": "A", "ot_datatypes": "clinical", "gwas_loci": ""},
        {"condition_id": "c", "source": "ot_other", "gene_symbol": "B", "ot_datatypes": "animal_model;clinical;literature",
         "gwas_loci": ""},
        {"condition_id": "c", "source": "gwas_catalog", "gene_symbol": "H1", "ot_datatypes": "", "gwas_loci": "6:32658495"},
        {"condition_id": "c", "source": "gwas_catalog", "gene_symbol": "H2", "ot_datatypes": "", "gwas_loci": "6:33087546"},
    ]).set_index(["condition_id", "source", "gene_symbol"])
    assert C.overlap_diagnostics(["A", "B"], "c", "ot_other", info)["share_overlap_drug_target_only"] == 0.5
    assert C.overlap_diagnostics(["H1", "H2"], "c", "gwas_catalog", info)["n_independent_gwas_loci"] == 1


# =========================================================================== review additions (2026-09-23)
def test_condition_specific_empirical_uses_attainable_threshold():
    # 16 cross pairs: smallest attainable p = 1/17 = 0.0588 > 0.05 -> above-all-cross must count as specific
    assert C.condition_specific_empirical(1 / 17, 16)
    assert not C.condition_specific_empirical(2 / 17, 16)
    # 22 cross pairs: 1/23 = 0.043 < 0.05 -> usual threshold
    assert C.condition_specific_empirical(1 / 23, 22)
    assert not C.condition_specific_empirical(2 / 23, 22)
    assert not C.condition_specific_empirical(np.nan, 5) and not C.condition_specific_empirical(0.01, 0)


def test_same_vs_cross_reports_attainable_p():
    S = pd.DataFrame([[5.0, 1.0, 0.0], [1.0, 4.0, 2.0], [0.0, 2.0, 0.5]], index=list("abc"), columns=list("abc"))
    per = C.same_vs_cross(S).set_index("condition_id")
    assert per.loc["a", "min_attainable_empirical_p"] == pytest.approx(1 / 5)
    assert per.loc["a", "same_above_all_cross"] and not per.loc["c", "same_above_all_cross"]


def test_holm():
    from statsmodels.stats.multitest import multipletests
    p = np.array([0.0247, 0.3391])
    assert np.allclose(C.holm(p), multipletests(p, method="holm")[1])


@pytest.mark.data
def test_disagreements_never_list_above_all_cross_pairs():
    dis = _csv("test3_disagreements")
    per = _csv("test3_gene_agreement_per_condition")
    per = per[per["statistic"] == "log2_fold"]
    d = dis[dis["kind"] == "overlap_not_condition_specific"].merge(
        per[["condition_id", "source_a", "source_b", "same_above_all_cross"]],
        on=["condition_id", "source_a", "source_b"], how="left")
    assert not d["same_above_all_cross"].fillna(False).astype(bool).any()
    assert {"n_cross_pairs", "min_attainable_empirical_p"} <= set(dis.columns)


@pytest.mark.data
def test_flag_concordance_definitions_labelled_and_adjusted():
    fc = _csv("test3_flag_concordance").set_index("observed_definition")
    assert set(C.FLAG_DEFINITIONS_PRESPECIFIED) <= set(fc.index)
    assert (fc.loc[list(C.FLAG_DEFINITIONS_PRESPECIFIED), "definition_status"] == "pre-specified").all()
    assert fc.loc[list(C.FLAG_DEFINITIONS_PRESPECIFIED), "perm_p_holm_prespecified"].notna().all()
    others = fc.drop(index=list(C.FLAG_DEFINITIONS_PRESPECIFIED))
    assert others["definition_status"].str.startswith("exploratory").all()
    assert others["perm_p_holm_prespecified"].isna().all()
    sens = _csv("test3_sensitivity")
    assert (sens["level"] == "system_flag_concordance").any()


@pytest.mark.data
def test_fig5_edges_carry_post_hoc_flags():
    edges = _csv("fig5_edges")
    assert "post_hoc_flag" in edges.columns
    se = _csv("test3_system_enrichment")
    assert "top_drug_mechanisms_behind_overlap" in se.columns
    flagged = se[(se["source"] == "ot_other") & se["supported"] & (se["share_overlap_drug_target_only"] >= 0.5)]
    # a drug record is not always retrieved for the id carrying the clinical evidence (e.g. dysautonomia
    # component ids), so require names for most, not all, flagged rows
    assert flagged["top_drug_mechanisms_behind_overlap"].fillna("").str.len().gt(0).mean() >= 0.5
