"""Tests for measure_it.omics (open_targets, gwas_catalog, geo_sra, query).

Unit tests are pure functions (no network). @pytest.mark.data tests check the outputs of
    uv run python -m measure_it.omics.open_targets
    uv run python -m measure_it.omics.gwas_catalog
    uv run python -m measure_it.omics.geo_sra
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import pytest
import yaml

from measure_it.config import RAW, UNKNOWN, load_config
from measure_it.omics import geo_sra as GS
from measure_it.omics import gwas_catalog as GW
from measure_it.omics import open_targets as OT
from measure_it.omics import query as Q
from measure_it.provenance import EVIDENCE_TYPES, PROVENANCE_COLUMNS
from measure_it.registry import REQUIRED_FIELDS
from measure_it.store import read_table, table_exists

CTX = {"condition_id": "me_cfs", "ontology_id": "MONDO:0005404", "ontology_id_role": "primary",
       "ontology_match": "exact", "source_disease_id": "MONDO_0005404",
       "source_disease_label": "myalgic encephalomeyelitis/chronic fatigue syndrome", "date_retrieved": "2026-09-23"}


# --------------------------------------------------------------------------- query helpers
def test_curie_short_roundtrip():
    assert Q.curie_to_short("MONDO:0005404") == "MONDO_0005404"
    assert Q.short_to_curie("MONDO_0005404") == "MONDO:0005404"
    assert Q.short_to_curie("MONDO:0005404") == "MONDO:0005404"
    assert Q.curie_to_short("EFO:0004540") == "EFO_0004540"


def test_finalize_columns_orders_schema_first_and_keeps_extras():
    df = pd.DataFrame([{"extra": 1, "condition_id": "x", "source_score": "0.5", "data_layer": "condition_molecular",
                        "evidence_type": "known_drug"}])
    out = Q.finalize_columns(df)
    assert list(out.columns[: len(Q.EVIDENCE_COLUMNS)]) == Q.EVIDENCE_COLUMNS
    assert "extra" in out.columns and out["source_score"].dtype == "float64"
    assert out.columns.tolist().count("evidence_type") == 1


def test_stringify_object_columns_makes_lists_json():
    df = pd.DataFrame({"a": [[1, 2], None], "b": [1, "x"]})
    out = Q.stringify_object_columns(df)
    assert out.loc[0, "a"] == "[1, 2]"
    assert out["b"].tolist() == ["1", "x"]


def test_layer_label_is_enrichment_not_patient_omics():
    assert "condition-level molecular enrichment" in Q.LAYER_LABEL
    assert "not patient multi-omics" in Q.LAYER_LABEL


# --------------------------------------------------------------------------- Open Targets (pure)
ASSOC = [{"score": 0.8, "target": {"id": "ENSG1", "approvedSymbol": "A", "approvedName": "a", "biotype": "protein_coding"},
          "datatypeScores": [{"id": "genetic_association", "score": 0.7}, {"id": "literature", "score": 0.2}],
          "datasourceScores": [{"id": "gwas_credible_sets", "score": 0.7}, {"id": "europepmc", "score": 0.2}]},
         {"score": 0.1, "target": {"id": "ENSG2", "approvedSymbol": "B", "approvedName": "b", "biotype": "protein_coding"},
          "datatypeScores": [{"id": "literature", "score": 0.1}], "datasourceScores": [{"id": "europepmc", "score": 0.1}]}]


def test_association_rows_keep_source_scores_exactly():
    rows = OT.association_rows(CTX, ASSOC, {"ENSG1": 0.75}, count=2)
    assert rows[0]["source_score"] == 0.8 and rows[0]["source_score_label"] == OT.SCORE_LABEL
    assert rows[0]["ot_datatype_score__genetic_association"] == 0.7
    assert "ot_datatype_score__genetic_association" not in rows[1]  # absent stays absent (null), not 0
    assert rows[0]["direct_association_score"] == 0.75 and rows[0]["has_direct_association"]
    assert np.isnan(rows[1]["direct_association_score"]) and not rows[1]["has_direct_association"]
    assert json.loads(rows[0]["ot_datasource_scores_json"]) == {"europepmc": 0.2, "gwas_credible_sets": 0.7}
    assert rows[0]["association_rank"] == 1 and rows[0]["evidence_type"] == "target_disease_association"


def _ev(**kw):
    base = {"id": "ev1", "datasourceId": "gwas_credible_sets", "datatypeId": "genetic_association", "score": 0.9,
            "resourceScore": None, "target": {"id": "ENSG1", "approvedSymbol": "A"},
            "disease": {"id": "MONDO_0005404", "name": "x"}, "diseaseFromSource": None, "studyId": None,
            "studySampleSize": None, "cohortShortName": None, "variantRsId": None, "variant": None,
            "credibleSet": {"studyLocusId": "sl", "pValueMantissa": 3.0, "pValueExponent": -9, "beta": 0.1,
                            "variant": {"id": "1_100_A_G", "rsIds": ["rs1"]},
                            "study": {"id": "GCST1", "traitFromSource": "CFS", "nSamples": 1000, "nCases": 10,
                                      "nControls": 990, "initialSampleSize": "10 cases, 990 controls",
                                      "pubmedId": "1", "publicationFirstAuthor": "X",
                                      "discoverySamples": [{"ancestry": "European", "sampleSize": 1000}]}},
            "literature": [], "publicationYear": None, "pValueMantissa": None, "pValueExponent": None,
            "oddsRatio": None, "beta": None, "confidence": None, "drug": None, "clinicalStage": None,
            "directionOnTrait": None, "directionOnTarget": None, "reactionId": None, "reactionName": None,
            "pathways": None, "contrast": None, "log2FoldChangeValue": None, "log2FoldChangePercentileRank": None,
            "biologicalModelId": None, "biologicalModelAllelicComposition": None,
            "diseaseModelAssociatedModelPhenotypes": None, "releaseVersion": "26.06"}
    base.update(kw)
    return base


def test_evidence_row_genetic_uses_credible_set_study():
    r = OT.evidence_row(CTX, _ev(), "genetic_association")
    assert r["evidence_type"] == "genetic_association" and r["source_accession"] == "GCST1"
    assert r["variant_rsids"] == "rs1" and r["sample_size"] == 1000
    assert r["p_value"] == pytest.approx(3e-9)
    assert r["study_population"].startswith("European (n=1000)")
    assert r["evidence_direction"] is None  # not defined by the source for this item
    assert r["evidence_from_descendant"] is False


def test_evidence_row_clinical_maps_to_known_drug_with_direction():
    e = _ev(datasourceId="clinical_precedence", datatypeId="clinical", credibleSet=None,
            drug={"id": "CHEMBL1", "name": "D"}, clinicalStage="PHASE_2", directionOnTarget="LoF",
            directionOnTrait="protect", disease={"id": "MONDO_9", "name": "child"})
    r = OT.evidence_row(CTX, e, "known_drug")
    assert r["evidence_type"] == "known_drug" and r["source_accession"] == "CHEMBL1"
    assert r["evidence_direction"] == "target:LoF, trait:protect"
    assert r["evidence_from_descendant"] is True
    assert r["evidence_type"] in EVIDENCE_TYPES


def test_every_datatype_mapping_is_controlled_vocabulary():
    assert set(OT.DATATYPE_TO_EVIDENCE_TYPE.values()) | {OT.DEFAULT_EVIDENCE_TYPE} <= EVIDENCE_TYPES


def test_select_top_evidence_caps_per_target():
    rows = [{"score": s, "target": {"id": t}} for t in ("a", "b") for s in (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3)]
    out = OT.select_top_evidence(rows, per_target=2)
    assert len(out) == 4 and all(r["score"] >= 0.8 for r in out)


def test_drug_rows_summarise_mechanism_and_trials():
    drugs = [{"id": "h", "maxClinicalStage": "PHASE_2",
              "drug": {"id": "CHEMBL19019", "name": "NALTREXONE", "drugType": "Small molecule",
                       "maximumClinicalStage": "APPROVAL",
                       "mechanismsOfAction": {"rows": [{"mechanismOfAction": "Opioid antagonist", "actionType": "ANTAGONIST",
                                                        "targetName": "x", "targets": [{"id": "ENSG3", "approvedSymbol": "OPRM1"}]}]}},
              "clinicalReports": [{"id": "nct1", "source": "AACT", "type": "CLINICAL_TRIAL", "trialOverallStatus": "RECRUITING",
                                   "diseases": [{"disease": {"id": "MONDO_0005404"}}]},
                                  {"id": "nct2", "source": "AACT", "type": "CLINICAL_TRIAL", "trialOverallStatus": None,
                                   "diseases": []}]}]
    r = OT.drug_rows(CTX, drugs)[0]
    assert r["entity_type"] == "drug" and r["entity_id"] == "CHEMBL19019"
    assert r["mechanism_target_symbols"] == "OPRM1" and r["action_types"] == "ANTAGONIST"
    assert json.loads(r["trial_status_counts"]) == {"RECRUITING": 1, "not reported": 1}
    assert r["reports_cite_queried_id"] is True and np.isnan(r["source_score"])


def test_pathway_rows_count_annotated_targets_without_scoring():
    tp = {"ENSG1": {"pathways": [{"pathwayId": "R-HSA-1", "pathway": "P1", "topLevelTerm": "T"}]},
          "ENSG2": {"pathways": [{"pathwayId": "R-HSA-1", "pathway": "P1", "topLevelTerm": "T"},
                                 {"pathwayId": "R-HSA-2", "pathway": "P2", "topLevelTerm": "T"}]}}
    rows = {r["entity_id"]: r for r in OT.pathway_rows(CTX, ASSOC, tp)}
    assert rows["R-HSA-1"]["n_associated_targets_annotated"] == 2
    assert rows["R-HSA-1"]["n_genetic_association_targets_annotated"] == 1
    assert rows["R-HSA-1"]["top_target_symbols"] == "A;B"
    assert np.isnan(rows["R-HSA-2"]["source_score"])


def test_md_table_handles_lists_and_pipes():
    t = OT._md_table(pd.DataFrame([{"a": ["x", "y"], "b": "p|q", "c": np.nan}]))
    assert "x; y" in t and "p/q" in t


# --------------------------------------------------------------------------- GWAS Catalog (pure)
def test_reported_genes_from_loci():
    payload = {"_embedded": {"loci": [{"author_reported_genes": [{"gene_name": "LRP1"}, "TRPM8", {"gene_name": "LRP1"}]},
                                      {"author_reported_genes": []}]}}
    assert GW.reported_genes_from_loci(payload) == ["LRP1", "TRPM8"]
    assert GW.reported_genes_from_loci(None) == []


def test_ancestry_sample_size_sums_initial_only():
    payload = {"_embedded": {"ancestries": [
        {"type": "initial", "number_of_individuals": 100, "ancestral_groups": [{"ancestral_group": "European"}]},
        {"type": "initial", "number_of_individuals": 50, "ancestral_groups": [{"ancestral_group": "East Asian"}]},
        {"type": "replication", "number_of_individuals": 999, "ancestral_groups": []}]}}
    n, pop = GW.ancestry_sample_size(payload)
    assert n == 150 and pop == "European (n=100); East Asian (n=50)"
    n2, pop2 = GW.ancestry_sample_size({"_embedded": {"ancestries": []}})
    assert np.isnan(n2) and pop2 == ""


def test_beta_direction_and_float_parsing():
    assert GW.parse_beta_direction({"beta": "1.826 unit decrease"}) == "decrease"
    assert GW.parse_beta_direction({"beta": "-"}) is None
    assert GW.parse_beta_direction({"beta_direction": "increase"}) == "increase"
    assert np.isnan(GW.to_float("-")) and np.isnan(GW.to_float(None)) and GW.to_float("0.964") == 0.964


def test_gwas_association_row():
    a = {"association_id": 7, "risk_frequency": "0.58", "pvalue_mantissa": 5, "pvalue_exponent": -8, "p_value": 5e-08,
         "or_value": "0.964", "beta": "-", "range": "[0.951-0.977]", "efo_traits": [{"efo_id": "MONDO_0005404", "efo_trait": "x"}],
         "reported_trait": ["CFS"], "accession_id": "GCST9", "mapped_genes": ["G1", "G2"],
         "snp_effect_allele": ["rs325506-G"], "snp_allele": [{"rs_id": "rs325506", "effect_allele": "G"}],
         "pubmed_id": "123", "locations": ["1:1"]}
    r = GW.association_row(CTX, a, ["G1"], {"initial_sample_size": "10 cases"}, 100.0, "European (n=100)", "efo_id")
    assert r["entity_type"] == "variant" and r["entity_id"] == "rs325506" and r["risk_allele"] == "G"
    assert r["or_value"] == 0.964 and np.isnan(r["beta_value"]) and r["beta_text"] is None
    assert r["beta_unit"] is None and r["beta_direction"] is None
    assert r["mapped_genes"] == "G1;G2" and r["reported_genes"] == "G1"
    assert r["association_trait_is_queried_id"] and np.isnan(r["source_score"])
    assert r["evidence_direction"] is None and r["sample_size"] == 100.0
    assert r["reported_genes_status"] == "reported"
    unknown = GW.association_row(CTX, a, None, None, np.nan, "", "efo_id")
    assert unknown["reported_genes"] is None and unknown["reported_genes_status"].startswith("loci request failed")
    assert GW.association_row(CTX, a, [], None, np.nan, "", "efo_id")["reported_genes_status"] == "none reported"


def test_parse_beta_text_all_v2_patterns():
    # the three textual forms the REST v2 API uses (no numeric beta field exists in v2)
    assert GW.parse_beta_text("0.0420454 unit decrease") == (0.0420454, "unit", "decrease")
    assert GW.parse_beta_text("0.07189597 unit increase") == (0.07189597, "unit", "increase")
    assert GW.parse_beta_text("0.030532 increase") == (0.030532, None, "increase")
    v, u, d = GW.parse_beta_text("-")
    assert np.isnan(v) and u is None and d is None
    assert np.isnan(GW.parse_beta_text(None)[0])


def test_gwas_association_row_keeps_beta_from_text():
    a = {"association_id": 8, "beta": "1.826 unit decrease", "or_value": None, "p_value": 2e-9,
         "efo_traits": [], "snp_allele": [{"rs_id": "rs1", "effect_allele": "T"}], "snp_effect_allele": ["rs1-T"],
         "accession_id": "GCST1"}
    r = GW.association_row(CTX, a, [], None, np.nan, "", "efo_id")
    assert r["beta_text"] == "1.826 unit decrease" and r["beta_value"] == 1.826
    assert r["beta_unit"] == "unit" and r["beta_direction"] == "decrease"
    assert r["evidence_direction"] == "beta decrease per effect allele (as reported)"
    assert np.isnan(r["or_value"]) and np.isnan(r["source_score"])


# --------------------------------------------------------------------------- query tool (pure)
def test_excluding_broad_ids_keeps_text_query_rows():
    d = pd.DataFrame({
        "ontology_match": ["broad", "broad", "narrow", None],
        "source_disease_id": ["MONDO_0021669", "", "MONDO_0011479", None],  # OT id row, GEO text row, OT, mapMECFS
    })
    keep = Q._keep_when_excluding_broad(d).tolist()
    assert keep == [False, True, True, True]


def test_entity_counts_ranks_gwas_variants_by_smallest_p():
    df = pd.DataFrame({
        "entity_id": ["rs1", "rs2", "rs2", "rs3"], "entity_type": "variant", "entity_label": ["a", "b", "b", "c"],
        "source_score": np.nan, "sample_size": [10.0, 5.0, 5.0, 1000.0], "source_database": "GWAS Catalog",
        "source_accession": ["G1", "G2", "G3", "G4"], "source_score_label": None,
        "p_value": [1e-8, 1e-20, 1e-6, 1e-7], "mentions_condition": np.nan,
    })
    out = Q._entity_counts(df, 3)
    assert [e["entity_id"] for e in out] == ["rs2", "rs1", "rs3"]  # 1e-20 < 1e-8 < 1e-7
    assert out[0]["min_reported_p"] == 1e-20 and "p-value" in out[0]["rank_basis"]


def test_descendant_evidence_drivers_counts_descendant_share():
    base = {"condition_id": "dysautonomia", "ontology_id": "MONDO:0001292", "source_disease_label": "ANS disorder",
            "n_descendants_in_source": 41, "entity_type": "gene"}
    rows = [
        {**base, "source_evidence_category": "genetic_association", "evidence_from_descendant": True,
         "evidence_disease_label": "paraganglioma"},
        {**base, "source_evidence_category": "genetic_association", "evidence_from_descendant": True,
         "evidence_disease_label": "paraganglioma"},
        {**base, "source_evidence_category": "literature", "evidence_from_descendant": False,
         "evidence_disease_label": "ANS disorder"},
        {**base, "source_evidence_category": OT.ASSOC_CATEGORY, "evidence_from_descendant": None,
         "evidence_disease_label": None},  # association rows are not evidence rows
        {**base, "ontology_id": "MONDO:0005404", "n_descendants_in_source": 0,  # no descendants -> not listed
         "source_evidence_category": "literature", "evidence_from_descendant": False, "evidence_disease_label": "x"},
    ]
    out = OT.descendant_evidence_drivers(pd.DataFrame(rows))
    assert len(out) == 1
    r = out.iloc[0]
    assert r["evidence_rows"] == 3 and r["evidence_rows_from_descendants"] == 2
    assert r["genetic_evidence_rows"] == 2 and r["genetic_rows_from_descendants"] == 2
    assert r["top_genetic_evidence_diseases"] == "paraganglioma (2)"


# --------------------------------------------------------------------------- GEO / SRA (pure)
EXPXML = ('<Summary><Title>t</Title><Platform instrument_model="NovaSeq">ILLUMINA</Platform></Summary>'
          '<Submitter acc="ERA1" center_name="QIB"/><Experiment acc="ERX1"/><Study acc="ERP146749" name="IgG-Seq ME"/>'
          '<Organism taxid="9606" ScientificName="Homo sapiens"/><Sample acc="ERS1"/>'
          '<Library_descriptor><LIBRARY_STRATEGY>WGA</LIBRARY_STRATEGY><LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>'
          '<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION></Library_descriptor><Bioproject>PRJEB61661</Bioproject>'
          '<Biosample>SAMEA1</Biosample>')


def test_parse_sra_expxml():
    p = GS.parse_sra_expxml(EXPXML, '<Run acc="ERR1"/><Run acc="ERR2"/>')
    assert p["study_acc"] == "ERP146749" and p["bioproject"] == "PRJEB61661"
    assert p["library_strategy"] == "WGA" and p["platform"] == "ILLUMINA" and p["run_accs"] == ["ERR1", "ERR2"]


def test_mentions_and_build_term():
    assert GS.mentions("Long COVID in PASC patients", [r"long[ -]?covid", r"\bPASC\b"]) == [r"long[ -]?covid", r"\bPASC\b"]
    assert GS.mentions("pancreatic PASCAL", [r"\bPASC\b"]) == []
    cfg = {"geo_filter": "F", "sra_filter": "S"}
    assert GS.build_term("x", "gds", cfg) == "(x) AND F" and GS.build_term("x", "sra", cfg) == "(x) AND S"


def test_geo_record():
    r = GS.geo_record("2001", {"accession": "GSE1", "gpl": "24676;21290", "pubmedids": ["38383456"], "pdat": "2023/12/22",
                               "gdstype": "Expression profiling by array", "entrytype": "GSE", "n_samples": 27})
    assert r["platform"] == "GPL24676;GPL21290" and r["release_date"] == "2023-12-22" and r["pubmed_ids"] == "38383456"


def test_precision_sample_is_seeded_and_frozen():
    geo = pd.DataFrame({"condition_id": ["long_covid"] * 15 + ["me_cfs"] * 3,
                        "gse": [f"GSE{i}" for i in range(18)], "study_type_in_scope": True})
    a = GS.precision_sample(geo)
    b = GS.precision_sample(geo)
    assert a["gse"].tolist() == b["gse"].tolist()
    assert (a["condition_id"] == "long_covid").sum() == 10 and (a["condition_id"] == "me_cfs").sum() == 3
    grown = pd.concat([geo, pd.DataFrame({"condition_id": ["long_covid"] * 5, "gse": [f"GSE9{i}" for i in range(5)],
                                          "study_type_in_scope": True})])
    c = GS.precision_sample(grown, previous=a)
    assert set(c.loc[c.condition_id == "long_covid", "gse"]) == set(a.loc[a.condition_id == "long_covid", "gse"])


def test_geo_query_config_rules():
    cfg = load_config("geo_queries")
    reg = set(read_table("condition_registry")["canonical_condition_id"]) if table_exists("condition_registry") else set(cfg["conditions"])
    assert set(cfg["conditions"]) == reg, "every registry condition needs at least one curated query"
    for cid, qs in cfg["conditions"].items():
        assert qs, cid
        for q in qs:
            t = q["term"]
            for acr in ("POTS[", "MCAS[", "IBS[", "EDS["):
                assert acr not in t, f"bare acronym {acr} in {q['id']}"
            if "PASC[" in t:
                assert "COVID" in t, "PASC must be paired with a COVID term"
            assert set(q["db"]) <= {"gds", "sra"} and q.get("rationale")
        for pat in cfg["qc_patterns"][cid]:
            re.compile(pat)
    assert "Homo sapiens" in cfg["geo_filter"] and "gse" in cfg["geo_filter"]


# --------------------------------------------------------------------------- data tests
PARTS = {
    "open_targets": "condition_molecular_evidence__open_targets",
    "gwas_catalog": "condition_molecular_evidence__gwas_catalog",
    "ncbi_geo_sra": "condition_molecular_evidence__geo",
}


def _need(name):
    if not table_exists(name):
        pytest.skip(f"{name} not built")
    return read_table(name)


@pytest.mark.data
@pytest.mark.parametrize("part", list(PARTS.values()))
def test_partition_schema_and_provenance(part):
    df = _need(part)
    assert len(df) > 0
    for c in Q.EVIDENCE_COLUMNS + PROVENANCE_COLUMNS:
        assert c in df.columns, c
    assert set(df["data_layer"]) == {"condition_molecular"}
    assert set(df["evidence_type"]) <= EVIDENCE_TYPES
    assert set(df["entity_type"]) <= {"gene", "variant", "pathway", "drug", "study", "analyte"}
    assert df["condition_id"].notna().all() and df["source_database"].notna().all()
    assert not [c for c in df.columns if "omics_score" in c.lower() or "combined_score" in c.lower()]
    reg = set(read_table("condition_registry")["canonical_condition_id"])
    assert set(df["condition_id"]) <= reg


@pytest.mark.data
def test_open_targets_outputs():
    df = _need(PARTS["open_targets"])
    a = df[df["source_evidence_category"] == OT.ASSOC_CATEGORY]
    assert (a["evidence_type"] == "target_disease_association").all()
    assert not a.duplicated(["condition_id", "ontology_id", "entity_id"]).any()  # one row per id x target
    assert (a["source_score_label"] == OT.SCORE_LABEL).all()
    assert a["source_score"].between(0, 1).all()
    dt_cols = [c for c in a.columns if c.startswith("ot_datatype_score__")]
    assert "ot_datatype_score__genetic_association" in dt_cols and "ot_datatype_score__literature" in dt_cols
    vals = a[dt_cols].to_numpy(dtype=float).ravel()
    vals = vals[~np.isnan(vals)]
    assert len(vals) and ((vals >= 0) & (vals <= 1)).all()
    # only source-provided scores: drug and pathway rows carry none
    assert df.loc[df["entity_type"].isin(["drug", "pathway"]), "source_score"].isna().all()
    assert {"gene", "drug", "pathway"} <= set(df["entity_type"])
    # descendants: the API default equals enableIndirect:true for every queried id
    log = json.loads((RAW / "open_targets" / "run_log.json").read_text())
    for sid, c in log["default_includes_descendants"].items():
        assert c["default"] == c["enableIndirect_true"] >= c["direct"], sid
    # every association call retrieved each target exactly once (pagination defect guarded)
    for key, pg in log["association_paging"].items():
        assert pg["count"] == pg["rows"] == pg["unique_targets"], key
    n_by_id = a.groupby(["condition_id", "ontology_id"])["entity_id"].nunique()
    api = {k.split("|")[0]: v["count"] for k, v in log["association_paging"].items() if k.endswith("|default")}
    for (cid, oid), n in n_by_id.items():
        assert n == api[Q.curie_to_short(oid)], (cid, oid)
    assert "Open Targets Platform data release" in df["source_version"].iloc[0]
    cov = _need("condition_molecular_coverage__open_targets")
    assert set(cov["condition_id"]) == set(read_table("condition_registry")["canonical_condition_id"])
    # a condition none of whose ids exist in Open Targets is recorded, not dropped
    none_ids = cov[cov["coverage_status"] == "no identifier"]["condition_id"].tolist()
    for cid in none_ids:
        assert cid not in set(a["condition_id"])


@pytest.mark.data
def test_gwas_outputs():
    df = _need(PARTS["gwas_catalog"])
    v = df[df["entity_type"] == "variant"]
    s = df[df["entity_type"] == "study"]
    assert len(v) > 0 and len(s) > 0
    assert v["entity_id"].str.contains(r"rs\d+").mean() > 0.9
    assert v["p_value"].dropna().between(0, 1).all()
    assert s["source_accession"].str.startswith("GCST").all()
    assert df["source_score"].isna().all()
    # every curated beta text is kept and split (REST v2 has no numeric beta field)
    has_beta = v["beta_text"].notna()
    assert has_beta.any() and v.loc[has_beta, "beta_value"].notna().all()
    assert set(v.loc[has_beta, "beta_direction"]) <= {"increase", "decrease"}
    assert v.loc[~has_beta, "beta_value"].isna().all()
    cov = _need("condition_molecular_coverage__gwas_catalog")
    assert set(cov["condition_id"]) == set(read_table("condition_registry")["canonical_condition_id"])
    assert {"zero hits", "identifier not present in source", "queried"} & set(cov["coverage_status"])
    # a 404 trait id never produces rows
    absent = cov[cov["id_status_in_source"].str.startswith("not a GWAS")]["ontology_id"]
    assert not set(absent) & set(df["ontology_id"])


@pytest.mark.data
def test_geo_outputs_and_precision():
    geo = _need("geo_study_catalog")
    sra = _need("sra_study_catalog")
    ev = _need(PARTS["ncbi_geo_sra"])
    assert geo["gse"].str.match(r"^GSE\d+$").all()
    assert (geo["taxon"].str.contains("Homo sapiens")).all()
    assert geo["n_samples"].ge(1).all()
    assert sra["study_acc"].str.match(r"^[SED]RP\d+$").all()
    assert (ev["entity_type"] == "study").all()
    assert set(ev["sample_size_basis"].dropna().str.contains("not participants")) == {True}
    assert len(ev) == len(geo) + len(sra)
    prec = json.loads((RAW / "ncbi_geo_sra" / "precision_check.json").read_text())
    for cid in Q.DEMO_CLUSTER:
        p = prec[cid]
        assert p["n_sampled"] == min(10, p["n_series_in_pool"]) or p["n_series_in_pool"] == 0 or cid == "pots"
        assert p["n_unreviewed"] == 0, f"{cid}: sampled GEO series not reviewed"
        assert p["precision"] is not None
    qlog = json.loads((RAW / "ncbi_geo_sra" / "query_log.json").read_text())
    assert {q["condition_id"] for q in qlog} == set(read_table("condition_registry")["canonical_condition_id"])


@pytest.mark.data
def test_union_contains_all_partitions():
    u = _need("condition_molecular_evidence")
    parts = sorted(p for p in (RAW.parent / "processed").glob("condition_molecular_evidence__*.parquet"))
    n = sum(len(pd.read_parquet(p, columns=["data_layer"])) for p in parts)
    assert len(u) == n
    for name in PARTS.values():
        if table_exists(name):
            src = read_table(name, columns=["source_name"])["source_name"].iloc[0]
            assert (u["source_name"] == src).any()


@pytest.mark.data
@pytest.mark.parametrize("source_id", list(PARTS))
def test_audit_and_registry_entry(source_id):
    d = RAW / source_id
    if not (d / "registry_entry.yaml").exists():
        pytest.skip(f"{source_id} not run")
    entry = yaml.safe_load((d / "registry_entry.yaml").read_text())
    assert not [f for f in REQUIRED_FIELDS if f not in entry]
    assert entry["data_layer"] == "condition_molecular" and entry["person_level"] is False
    assert entry["true_participant_linkage_across_modalities"] is False
    audit = (d / "DATA_AUDIT.md").read_text()
    for field in ("source_id", "Retrieval date", "Sample size", "Missingness", "Linkage strategy", "Limitations"):
        assert field in audit
    assert (d / "MANIFEST.json").exists()
    assert "patient multi-omics" not in audit.replace("never patient multi-omics", "").replace(
        "not patient multi-omics", "")


@pytest.mark.data
def test_get_molecular_context():
    if not table_exists("condition_molecular_evidence"):
        pytest.skip("union not built")
    r = Q.get_molecular_context("Long COVID", top_n=3)
    assert r["condition_id"] == "long_covid" and r["status"] == "available"
    assert "not patient multi-omics" in r["layer"]
    assert r["evidence"] and all("counts_by_source" in v for v in r["evidence"].values())
    assert r["provenance"]
    u = Q.get_molecular_context("definitely not a condition qqq")
    assert u["status"] == UNKNOWN and u["evidence"] == UNKNOWN
    by_id = Q.get_molecular_context("MONDO:0005404")
    assert by_id["condition_id"] == "me_cfs"
    # excluding broad ids must not drop POTS's text-query GEO/SRA rows (its primary id is a broad match)
    pots = Q.get_molecular_context("pots", include_broad=False)
    assert "expression_study_metadata" in pots["evidence"] or "metadata_catalog" in pots["evidence"]


@pytest.mark.data
def test_molecular_context_reads_tables_once_per_version(monkeypatch):
    """Table-level cache: the evidence / coverage tables and the all-condition genetic agreement are read and computed
    once per table version, not on every call; the result is identical to an uncached call."""
    if not table_exists("condition_molecular_evidence"):
        pytest.skip("union not built")
    Q._evidence_tables.cache_clear()
    first = Q.get_molecular_context("me_cfs", top_n=3)
    reads = []
    real = Q.read_table
    monkeypatch.setattr(Q, "read_table", lambda name, **kw: reads.append(name) or real(name, **kw))
    again = Q.get_molecular_context("me_cfs", top_n=3)
    Q.get_molecular_context("long_covid", top_n=3)
    assert reads == [] and json.dumps(first, default=str) == json.dumps(again, default=str)
    assert Q._evidence_tables.cache_info().hits >= 2
    Q._evidence_tables.cache_clear()
