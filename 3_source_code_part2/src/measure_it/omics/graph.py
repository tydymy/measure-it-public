"""Phenotype -> measurable-biology graph (SPEC section C / Phase 2 naming, Figure 5 data).

Builds, from the Test 3 outputs written by ``measure_it.omics.coherence``:

* processed table ``measurable_biology`` — condition x physiological system x candidate
  measurement class. The condition -> system link is **data-derived** (Reactome enrichment of
  public condition-level gene sets); the system -> measurement-class link is **curated** (table in
  docs/ANALYSIS_PLAN_MOLECULAR.md, with a rationale per row). A second row family records which
  omics measurement classes public study catalogues have already applied to a condition
  (research-readiness context, not molecular support).
* ``results/tables/fig5_nodes.csv`` / ``fig5_edges.csv`` and draft figures
  ``results/figures/drafts/fig5_*``.
* ``results/MOLECULAR_COHERENCE.md`` — the Test 3 write-up; every number is read from a file
  written by this package.

Everything is condition-level molecular enrichment; nothing describes a person.

Reproduce: ``uv run python -m measure_it.omics.coherence && uv run python -m measure_it.omics.graph``
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import FIGURES, RESULTS, TABLES, UNKNOWN, load_config, utc_now_iso
from ..provenance import PROVENANCE_COLUMNS, add_provenance
from ..store import read_table, table_exists, write_table
from . import coherence as C

PRODUCER = "measure_it.omics.graph"
TABLE = "measurable_biology"
DRAFTS = FIGURES / "drafts"
REPORT_PATH = RESULTS / "MOLECULAR_COHERENCE.md"
DEMO = ("long_covid", "me_cfs", "pots")
ASSAY_PLATFORM = "cross_system_assay_platform"
TOP_GENES, TOP_PATHWAYS = 5, 3

CONDITION_CAVEATS = {
    "pots": "Open Targets evidence comes only from MONDO:0011479, which Mondo labels 'POTS due to NET deficiency' "
            "(a monogenic form); Open Targets labels it plain POTS, so literature co-mentions may concern POTS in general.",
    "dysautonomia": "broad concept: Open Targets ids include autonomic nervous system disorder (direct targets only kept) "
                    "and GEO series include familial dysautonomia / MSA / Parkinson's studies.",
    "eds_hsd": "Ehlers-Danlos syndrome id has 42 descendant subtypes in Open Targets; only direct associations kept, "
               "which still include monogenic EDS genes.",
    "me_cfs": "mapMECFS gene lists are nominal p < 0.05 results from one small cohort with no FDR-significant rows.",
    "mcas": "no identifier usable by Open Targets or the GWAS Catalog and no GEO/SRA hits.",
    "post_infectious_syndrome": "grouping node with a broad id only; excluded from the primary analysis.",
    "ptlds": "only 4 Open Targets targets; too few genes for enrichment.",
}


# =========================================================================== inputs
def _csv(name: str) -> pd.DataFrame:
    p = TABLES / f"{name}.csv"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing; run `uv run python -m measure_it.omics.coherence` first")
    return pd.read_csv(p)


def _meta() -> dict:
    p = TABLES / "test3_run_metadata.json"
    return json.loads(p.read_text()) if p.exists() else {}


def measurement_names() -> dict[str, str]:
    return {m["id"]: m["name"] for m in load_config("measurements")["measurement_classes"]}


# =========================================================================== measurable_biology
def support_tier(sources_supporting: list[str], n_tested: int, n_genes_any: int) -> str:
    """Support tier of a condition -> system link (plan section 6)."""
    if n_tested == 0:
        return "too_few_genes" if n_genes_any > 0 else "no_molecular_data"
    fams = {C.SOURCE_FAMILY[s] for s in sources_supporting}
    if len(fams) >= 2:
        return "multi_source"
    if fams == {"literature"}:
        return "literature_only"
    if fams:
        return "single_source"
    return "not_supported"


FLAG_TEXT = {
    "gwas_single_locus": "post-hoc: GWAS Catalog support comes from genes of a single GWAS locus (positional mapping of "
                         "one region, e.g. the MHC), not from independent loci",
    "ot_other_mostly_drug_targets": "post-hoc: Open Targets 'other' support is mostly drug-target genes of trialled drugs "
                                    "(clinical precedence; one mechanism can list a whole gene family), not disease biology",
}


def post_hoc_flags(sup: pd.DataFrame) -> tuple[list[str], set[str], list[str]]:
    """Post-hoc diagnostic flags for the sources supporting one condition-system link (not pre-specified).

    Returns (flags, flagged sources, detail strings naming the drug records behind drug-target support).
    """
    flags, srcs, details = [], set(), []
    for x in sup.itertuples():
        if x.source == "gwas_catalog" and getattr(x, "n_independent_gwas_loci", np.nan) == 1:
            flags.append("gwas_single_locus")
            srcs.add(x.source)
        if x.source == "ot_other" and getattr(x, "share_overlap_drug_target_only", np.nan) >= 0.5:
            flags.append("ot_other_mostly_drug_targets")
            srcs.add(x.source)
            mech = getattr(x, "top_drug_mechanisms_behind_overlap", "")
            if isinstance(mech, str) and mech:
                details.append(f"drug records behind the Open Targets 'other' genes: {mech}")
    return sorted(set(flags)), srcs, details


def _system_rows(reg: pd.DataFrame, se: pd.DataFrame, counts: pd.DataFrame, cmap: pd.DataFrame,
                 names: dict) -> pd.DataFrame:
    rows = []
    for _, r in reg.iterrows():
        c = r["canonical_condition_id"]
        crow = counts[counts["condition_id"] == c]
        n_genes_any = int(crow[[f"genes_{s}" for s in C.SOURCES]].sum(axis=1).iloc[0]) if len(crow) else 0
        sc = se[se["condition_id"] == c]
        tested = sorted(sc["source"].unique())
        for sysid, g in cmap.groupby("physiological_system", sort=False):
            ss = sc[sc["physiological_system"] == sysid]
            sup = ss[ss["supported"] == True]  # noqa: E712
            srcs = sorted(sup["source"])
            tier = support_tier(srcs, len(tested), n_genes_any)
            spec = sup["condition_specific"]
            spec_srcs = sorted(sup.loc[spec == True, "source"])  # noqa: E712
            by_src = {x.source: {"n_enriched_pathways": int(x.n_enriched_pathways), "system_fdr": float(x.fdr),
                                 "system_log2_fold": round(float(x.log2_fold), 4),
                                 "condition_specific": (None if pd.isna(x.condition_specific) else bool(x.condition_specific))}
                      for x in ss.itertuples()}
            best = sup.sort_values("fdr").head(1)
            top_pw = "; ".join(dict.fromkeys(
                p for x in sup.sort_values("fdr")["top_enriched_pathways"].fillna("") for p in x.split("; ") if p))
            genes = sorted({gname for x in sup["genes_in_system"].fillna("") for gname in x.split(";") if gname})
            flags, flagged_srcs, flag_details = post_hoc_flags(sup)
            robust_tier = support_tier([x for x in srcs if x not in flagged_srcs], len(tested), n_genes_any)
            caveats = [CONDITION_CAVEATS[c]] if c in CONDITION_CAVEATS else []
            caveats += [FLAG_TEXT[f] for f in flags] + flag_details
            if tier == "literature_only":
                caveats.append("supported only by Europe PMC co-mention (study attention), not by genetic or experimental evidence")
            if tier in ("multi_source", "single_source", "literature_only") and not spec_srcs:
                caveats.append("no supporting source shows this system above the other conditions' median "
                               "(not condition-specific)")
            for m in g.itertuples():
                has_data = tier not in ("no_molecular_data", "too_few_genes")
                rows.append({
                    "object_id": f"molbio:{c}|{sysid}|{m.measurement_class}",
                    "condition_id": c, "condition_label": r["preferred_name"], "ontology_id": r["primary_mondo_id"],
                    "physiological_system": sysid, "physiological_system_label": C.SYSTEM_LABEL[sysid],
                    "measurement_class": m.measurement_class,
                    "measurement_name": names.get(m.measurement_class, UNKNOWN),
                    "support_status": tier,
                    "post_hoc_flags": ";".join(flags),
                    "support_status_without_flagged_sources": robust_tier,
                    "condition_system_link": ("data_derived: Reactome enrichment of public condition-level gene sets"
                                              if has_data else f"{UNKNOWN}: {tier.replace('_', ' ')}"),
                    "system_measurement_link": "curated (docs/ANALYSIS_PLAN_MOLECULAR.md section 7)",
                    "link_evidence": ("data_derived_condition_to_system + curated_system_to_measurement"
                                      if has_data else "curated_system_to_measurement_only (no molecular data)"),
                    "n_sources_tested": len(tested), "sources_tested": ";".join(tested),
                    "n_sources_supporting": len(srcs), "sources_supporting": ";".join(srcs),
                    "n_source_families_supporting": len({C.SOURCE_FAMILY[s] for s in srcs}),
                    "source_families_supporting": ";".join(sorted({C.SOURCE_FAMILY[s] for s in srcs})),
                    "condition_specific": (bool(spec_srcs) if len(srcs) else pd.NA),
                    "condition_specific_sources": ";".join(spec_srcs),
                    "n_enriched_pathways_supporting": int(sup["n_enriched_pathways"].sum()),
                    "evidence_by_source": json.dumps(by_src, sort_keys=True),
                    "best_system_fdr": float(best["fdr"].iloc[0]) if len(best) else np.nan,
                    "max_system_log2_fold": float(sup["log2_fold"].max()) if len(sup) else np.nan,
                    "top_pathways": top_pw[:600], "top_genes": ";".join(genes[:25]),
                    "linking_rule": ("system supported by a source when its gene set is enriched for the system's "
                                     "Reactome gene universe (BH FDR < 0.05) and >= 1 Reactome pathway of the system is "
                                     "enriched (FDR < 0.05, overlap >= 2); systems from the Reactome hierarchy anchor "
                                     "table; measurement class from the curated system map"),
                    "mapping_rationale": m.rationale,
                    "caveats": " | ".join(caveats),
                    "_etype": "pathway_annotation" if has_data else "curated_config",
                })
    return pd.DataFrame(rows)


def _platform_rows(reg: pd.DataFrame, names: dict) -> pd.DataFrame:
    """Omics measurement classes already applied to a condition in public study catalogues."""
    rows = []
    geo = read_table("geo_study_catalog") if table_exists("geo_study_catalog") else pd.DataFrame()
    sra = read_table("sra_study_catalog") if table_exists("sra_study_catalog") else pd.DataFrame()
    ev = C.load_evidence()
    mm = ev[ev["source_name"].astype(str).str.contains(C.MAPMECFS_SOURCE_TAG, regex=False)]
    geo_map = {"transcriptomics": r"Expression profiling|Non-coding RNA profiling", "proteomics": r"Protein profiling"}
    mm_map = {"transcriptomics": ["transcript", "transposable_element"], "proteomics": ["protein_aptamer"],
              "metabolomics": ["metabolite", "stool_metabolite", "lipid_class"], "immune_assays": ["immune_cell_population"]}
    for _, r in reg.iterrows():
        c = r["canonical_condition_id"]
        g = geo[(geo["condition_id"] == c) & geo["study_type_in_scope"].fillna(False).astype(bool)] if len(geo) else geo
        s = sra[sra["condition_id"] == c] if len(sra) else sra
        for mclass, pat in geo_map.items():
            gg = g[g["study_type"].astype(str).str.contains(pat)] if len(g) else g
            n_named = int(gg["mentions_condition"].fillna(False).astype(bool).sum()) if len(gg) else 0
            n_sra = 0
            if mclass == "transcriptomics" and len(s):
                n_sra = int((s["library_strategies"].astype(str).str.contains("RNA-Seq|miRNA-Seq|ncRNA-Seq")
                             & s["mentions_condition"].fillna(False).astype(bool)).sum())
            if len(gg) == 0 and n_sra == 0:
                continue
            rows.append({"condition_id": c, "measurement_class": mclass, "catalogue": "NCBI GEO / SRA",
                         "catalogue_key": "geo_sra",
                         "n_studies": int(len(gg)), "n_studies_naming_condition": n_named,
                         "n_sra_studies_naming_condition": n_sra,
                         "detail": f"{len(gg)} in-scope GEO series of this type ({n_named} name the condition in "
                                   f"title/summary); {n_sra} SRA RNA studies name it in the title",
                         "accessions": ";".join(gg.loc[gg["mentions_condition"].fillna(False).astype(bool), "gse"]
                                                .astype(str).head(10)) if len(gg) else "",
                         "_etype": "expression_study_metadata"})
        m = mm[mm["condition_id"] == c]
        for mclass, classes in mm_map.items():
            mt = m[m["analyte_class"].isin(classes)]
            if mt.empty:
                continue
            rows.append({"condition_id": c, "measurement_class": mclass,
                         "catalogue": "NIH PI-ME/CFS (Walitt 2024) published supplementary tables",
                         "catalogue_key": "mapmecfs",
                         "n_studies": 1, "n_studies_naming_condition": 1, "n_sra_studies_naming_condition": 0,
                         "detail": f"{len(mt)} published analyte rows ({int((pd.to_numeric(mt['p_value'], errors='coerce') < 0.05).sum())} "
                                   f"nominal p < 0.05; FDR-significant: {int((mt['fdr_significant'] == True).sum())}; "  # noqa: E712
                                   f"FDR status UNKNOWN: {int(mt['fdr_significant'].isna().sum())})",
                         "accessions": ";".join(sorted(set(mt["supplementary_table"].astype(str)))),
                         "_etype": "published_biomarker"})
    out = []
    for x in rows:
        c = x["condition_id"]
        rr = reg[reg["canonical_condition_id"] == c].iloc[0]
        out.append({
            "object_id": f"molbio:{c}|{ASSAY_PLATFORM}|{x['measurement_class']}|{x['catalogue_key']}", "condition_id": c,
            "condition_label": rr["preferred_name"], "ontology_id": rr["primary_mondo_id"],
            "physiological_system": ASSAY_PLATFORM,
            "physiological_system_label": "cross-system assay platform (study catalogue)",
            "measurement_class": x["measurement_class"], "measurement_name": names.get(x["measurement_class"], UNKNOWN),
            "support_status": "study_catalogue_evidence", "post_hoc_flags": "",
            "support_status_without_flagged_sources": "study_catalogue_evidence",
            "condition_system_link": "not applicable (assay platform, not a physiological system)",
            "system_measurement_link": "data_derived: the class has been applied to this condition in public studies",
            "link_evidence": "data_derived_study_catalogue",
            "n_sources_tested": np.nan, "sources_tested": x["catalogue"], "n_sources_supporting": np.nan,
            "sources_supporting": x["catalogue"], "n_source_families_supporting": np.nan,
            "source_families_supporting": "", "condition_specific": pd.NA, "condition_specific_sources": "",
            "n_enriched_pathways_supporting": np.nan, "evidence_by_source": json.dumps(
                {k: x[k] for k in ("n_studies", "n_studies_naming_condition", "n_sra_studies_naming_condition")}),
            "best_system_fdr": np.nan, "max_system_log2_fold": np.nan, "top_pathways": "",
            "top_genes": "", "linking_rule": "GEO in-scope series type / SRA library strategy / mapMECFS analyte class -> "
                                             "omics measurement class (study existence = research readiness, not effect)",
            "mapping_rationale": x["detail"], "caveats": ("study metadata from text queries; presence of a study is "
                                                          "not evidence of an effect") + (f" | accessions: {x['accessions']}"
                                                                                          if x["accessions"] else ""),
            "_etype": x["_etype"]})
    return pd.DataFrame(out)


def build_measurable_biology() -> pd.DataFrame:
    reg = read_table("condition_registry")
    se = _csv("test3_system_enrichment")
    counts = _csv("test3_source_counts")
    cmap = C.system_measurement_map()
    names = measurement_names()
    missing = set(cmap["measurement_class"]) - set(names)
    if missing:
        raise ValueError(f"curated map uses measurement classes absent from configs/measurements.yaml: {missing}")
    df = pd.concat([_system_rows(reg, se, counts, cmap, names), _platform_rows(reg, names)], ignore_index=True)
    meta = _meta()
    vers = C.version_string(meta.get("versions", {})) or UNKNOWN
    built = utc_now_iso()
    parts = []
    for et, g in df.groupby("_etype", sort=False):
        parts.append(add_provenance(
            g.drop(columns=["_etype"]), data_layer="derived",
            source_name=f"{PRODUCER} (measurable biology: {C.LAYER_LABEL} -> curated measurement classes)",
            source_version=vers, retrieved_at=built, evidence_type=et, source_record_id="object_id",
            source_geographic_resolution="none", evidence_level=g["support_status"],
            provenance_notes=("condition -> physiological system is data-derived from Test 3 (results/tables/"
                              "test3_system_enrichment.csv); system -> measurement class is a curated assumption "
                              "(docs/ANALYSIS_PLAN_MOLECULAR.md); study-catalogue rows are research-readiness context"
                              if et != "curated_config" else
                              "no molecular data for this condition: only the curated system -> measurement map applies")))
    out = pd.concat(parts, ignore_index=True)
    out["condition_specific"] = out["condition_specific"].astype("boolean")
    for c in ("n_sources_tested", "n_sources_supporting", "n_source_families_supporting",
              "n_enriched_pathways_supporting"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("float64")
    if out["object_id"].duplicated().any():
        raise ValueError("measurable_biology object_id not unique")
    return out


# =========================================================================== Figure 5 data
def build_fig5(mb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    gs = _csv("test3_gene_sets")
    ora = _csv("test3_reactome_ora")
    se = _csv("test3_system_enrichment")
    reg = read_table("condition_registry")
    cmap = C.system_measurement_map()
    names = measurement_names()
    nodes: dict[str, dict] = {}
    edges: list[dict] = []

    def node(nid, ntype, label, layer, **attr):
        nodes.setdefault(nid, {"node_id": nid, "node_type": ntype, "label": label, "layer": layer, **attr})

    info = gs.set_index(["condition_id", "source", "gene_symbol"])[["ot_datatypes", "gwas_loci"]]
    other_dt = set(C.OT_OTHER_DATATYPES)

    def edge_flag(c: str, s: str, genes: list[str]) -> str:
        """Post-hoc flag of the genes behind one condition->pathway edge (same rule as test3_system_enrichment)."""
        genes = [g for g in genes if g]
        if not genes:
            return ""
        d = C.overlap_diagnostics(genes, c, s, info)
        if s == "ot_other" and d["share_overlap_drug_target_only"] >= 0.5:
            return "ot_other_mostly_drug_targets"
        if s == "gwas_catalog" and d["n_independent_gwas_loci"] == 1:
            return "gwas_single_locus"
        return ""

    conds = sorted(set(se["condition_id"]))
    for c in conds:
        r = reg[reg["canonical_condition_id"] == c].iloc[0]
        node(f"condition:{c}", "condition", r["preferred_name"], 0, object_id=f"condition:{c}",
             demo_cluster=c in DEMO)
        # top genes per source
        for s in C.SOURCES:
            g = gs[(gs["condition_id"] == c) & (gs["source"] == s)]
            if g.empty:
                continue
            asc = s in ("gwas_catalog", "mapmecfs")
            top = g.sort_values(["rank_stat", "gene_symbol"], ascending=[asc, True]).head(TOP_GENES)
            for x in top.itertuples():
                nid = f"gene:{x.gene_symbol}"
                node(nid, "gene", x.gene_symbol, 1, object_id=(f"gene:{x.ensembl_gene_id}"
                                                                if isinstance(x.ensembl_gene_id, str) and x.ensembl_gene_id else ""))
                dts = set(str(x.ot_datatypes).split(";")) & other_dt if s == "ot_other" else set()
                edges.append({"source_node": f"condition:{c}", "target_node": nid, "edge_type": "condition_gene",
                              "evidence_source": s, "link_type": "data_derived", "weight": 1.0,
                              "statistic": float(x.rank_stat), "statistic_label": x.rank_stat_label,
                              "condition_id": c, "detail": f"top {TOP_GENES} {C.SOURCE_LABEL[s]} gene",
                              "post_hoc_flag": "drug_target_only" if dts == {"clinical"} else ""})
            # top enriched pathways per source
            o = ora[(ora["condition_id"] == c) & (ora["source"] == s) & (ora["enriched"] == True)]  # noqa: E712
            for x in o.sort_values(["fdr", "p_value"]).head(TOP_PATHWAYS).itertuples():
                pid = f"pathway:{x.reactome_id}"
                node(pid, "pathway", x.reactome_name, 2, object_id=f"reactome:{x.reactome_id}",
                     physiological_systems=x.physiological_systems)
                edges.append({"source_node": f"condition:{c}", "target_node": pid, "edge_type": "condition_pathway",
                              "evidence_source": s, "link_type": "data_derived",
                              "weight": float(-np.log10(max(x.fdr, 1e-300))), "statistic": float(x.fdr),
                              "statistic_label": "Reactome ORA BH FDR", "condition_id": c,
                              "detail": f"overlap {x.overlap} of {x.pathway_size}",
                              "post_hoc_flag": edge_flag(c, s, str(x.overlap_genes).split(";"))})
                top_genes_here = set(top["gene_symbol"])
                for gname in str(x.overlap_genes).split(";"):
                    if gname in top_genes_here:
                        edges.append({"source_node": f"gene:{gname}", "target_node": pid, "edge_type": "gene_pathway",
                                      "evidence_source": "reactome", "link_type": "data_derived", "weight": 1.0,
                                      "statistic": np.nan, "statistic_label": "Reactome membership",
                                      "condition_id": c, "detail": "gene annotated to pathway in Reactome"})
    # pathway -> system (Reactome hierarchy)
    for nid, n in list(nodes.items()):
        if n["node_type"] != "pathway":
            continue
        for sysid in str(n.get("physiological_systems", "other")).split(";"):
            node(f"system:{sysid}", "physiological_system", C.SYSTEM_LABEL.get(sysid, sysid), 3)
            edges.append({"source_node": nid, "target_node": f"system:{sysid}", "edge_type": "pathway_system",
                          "evidence_source": "reactome_hierarchy", "link_type": "data_derived (Reactome hierarchy + "
                                                                                  "curated anchor table)",
                          "weight": 1.0, "statistic": np.nan, "statistic_label": "", "condition_id": "",
                          "detail": "system from the anchor rule of plan 5.3"})
    # condition -> system (support summary) and system -> measurement (curated)
    sysrows = mb[mb["physiological_system"] != ASSAY_PLATFORM].drop_duplicates(["condition_id", "physiological_system"])
    for x in sysrows.itertuples():
        if x.condition_id not in conds or x.support_status in ("not_supported", "no_molecular_data", "too_few_genes"):
            continue
        node(f"system:{x.physiological_system}", "physiological_system", x.physiological_system_label, 3)
        edges.append({"source_node": f"condition:{x.condition_id}", "target_node": f"system:{x.physiological_system}",
                      "edge_type": "condition_system", "evidence_source": x.sources_supporting,
                      "link_type": "data_derived", "weight": float(x.n_source_families_supporting),
                      "statistic": float(x.best_system_fdr), "statistic_label": "best system-level BH FDR",
                      "condition_id": x.condition_id,
                      "detail": f"{x.support_status}; above other conditions' median={x.condition_specific}; "
                                f"tier without post-hoc-flagged sources={x.support_status_without_flagged_sources}",
                      "post_hoc_flag": x.post_hoc_flags if isinstance(x.post_hoc_flags, str) else ""})
    for m in cmap.itertuples():
        node(f"system:{m.physiological_system}", "physiological_system", C.SYSTEM_LABEL[m.physiological_system], 3)
        node(f"measurement:{m.measurement_class}", "measurement_class", names.get(m.measurement_class, m.measurement_class),
             4, object_id=f"measurement:{m.measurement_class}")
        edges.append({"source_node": f"system:{m.physiological_system}", "target_node": f"measurement:{m.measurement_class}",
                      "edge_type": "system_measurement", "evidence_source": "curated_config", "link_type": "curated",
                      "weight": 1.0, "statistic": np.nan, "statistic_label": "", "condition_id": "",
                      "detail": m.rationale})
    nd = pd.DataFrame(list(nodes.values()))
    ed = pd.DataFrame(edges).drop_duplicates(["source_node", "target_node", "edge_type", "evidence_source", "condition_id"])
    ed["post_hoc_flag"] = ed["post_hoc_flag"].fillna("")
    return nd, ed


# =========================================================================== draft figures
def draft_figures(mb: pd.DataFrame, nodes: pd.DataFrame, edges: pd.DataFrame) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    DRAFTS.mkdir(parents=True, exist_ok=True)
    out = []
    # ---- heatmap: conditions x systems, number of supporting source families
    sysrows = mb[mb["physiological_system"] != ASSAY_PLATFORM].drop_duplicates(["condition_id", "physiological_system"])
    systems = [s for s in C.SYSTEMS if s != "other"]
    conds = list(dict.fromkeys(read_table("condition_registry")["canonical_condition_id"]))
    mat = np.full((len(conds), len(systems)), np.nan)
    spec = np.zeros_like(mat, dtype=bool)
    lit = np.zeros_like(mat, dtype=bool)
    robust = np.full(mat.shape, -1)
    flag_src = {"ot_other_mostly_drug_targets": "ot_other", "gwas_single_locus": "gwas_catalog"}
    for i, c in enumerate(conds):
        for j, s in enumerate(systems):
            r = sysrows[(sysrows["condition_id"] == c) & (sysrows["physiological_system"] == s)]
            if r.empty or r["support_status"].iloc[0] in ("no_molecular_data", "too_few_genes"):
                continue
            mat[i, j] = float(r["n_source_families_supporting"].iloc[0])
            spec[i, j] = bool(r["condition_specific"].fillna(False).iloc[0])
            lit[i, j] = r["support_status"].iloc[0] == "literature_only"
            fl = [flag_src[f] for f in str(r["post_hoc_flags"].iloc[0] or "").split(";") if f in flag_src]
            if fl:
                srcs = [x for x in str(r["sources_supporting"].iloc[0]).split(";") if x and x not in fl]
                robust[i, j] = len({C.SOURCE_FAMILY[x] for x in srcs})
    fig, ax = plt.subplots(figsize=(9, 8))
    cmap = plt.get_cmap("Blues", 5)
    cmap.set_bad("#e6e6e6")
    im = ax.imshow(np.ma.masked_invalid(mat), cmap=cmap, vmin=-0.5, vmax=4.5, aspect="auto")
    for i in range(len(conds)):
        for j in range(len(systems)):
            if np.isnan(mat[i, j]):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="#777777")
                continue
            txt = f"{int(mat[i, j])}" + ("*" if spec[i, j] else "") + (" L" if lit[i, j] else "")
            if robust[i, j] >= 0:
                txt += f"\n(D->{int(robust[i, j])})"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.5, color="black" if mat[i, j] < 3 else "white")
    ax.set_xticks(range(len(systems)), [C.SYSTEM_LABEL[s] for s in systems], rotation=35, ha="right")
    ax.set_yticks(range(len(conds)), conds)
    cb = fig.colorbar(im, ax=ax, ticks=range(5))
    cb.set_label("source families supporting the system (shared data possible)")
    ax.set_title("Condition -> physiological system (Reactome enrichment of public gene sets)\n"
                 "* above the other conditions' median in >= 1 source;  L = literature co-mention only;  "
                 "n/a = no/too few genes\n(D->k) = post-hoc flag (drug-target gene family or single GWAS locus); "
                 "k families remain without it", fontsize=8.5)
    fig.tight_layout()
    for ext in ("png", "svg"):
        p = DRAFTS / f"fig5_system_heatmap.{ext}"
        fig.savefig(p, dpi=200)
        out.append(p)
    plt.close(fig)

    # ---- layered graph for the demo cluster: condition -> pathway -> system -> measurement class
    e = edges
    colors = {"ot_genetic": "#1b9e77", "gwas_catalog": "#d95f02", "ot_literature": "#7570b3",
              "ot_other": "#e7298a", "mapmecfs": "#66a61e"}
    cp = e[(e["edge_type"] == "condition_pathway") & e["condition_id"].isin(DEMO)]
    cs = e[(e["edge_type"] == "condition_system") & e["condition_id"].isin(DEMO)]
    pw_nodes = list(dict.fromkeys(cp["target_node"]))
    sys_nodes = [f"system:{s}" for s in C.SYSTEMS if f"system:{s}" in set(cs["target_node"])
                 | set(e[(e["edge_type"] == "pathway_system") & e["source_node"].isin(pw_nodes)]["target_node"])]
    sm = e[(e["edge_type"] == "system_measurement") & e["source_node"].isin(sys_nodes)]
    ms_nodes = list(dict.fromkeys(sm["target_node"]))
    cols = {0: [f"condition:{c}" for c in DEMO], 1: pw_nodes, 2: sys_nodes, 3: ms_nodes}
    pos = {}
    for x, ids in cols.items():
        for k, nid in enumerate(ids):
            pos[nid] = (x, 1 - (k + 0.5) / max(len(ids), 1))
    lab = dict(zip(nodes["node_id"], nodes["label"]))
    fig, ax = plt.subplots(figsize=(15, max(7, 0.28 * max(len(v) for v in cols.values()) + 2)))
    for r in cp.itertuples():
        (x0, y0), (x1, y1) = pos[r.source_node], pos[r.target_node]
        flagged = isinstance(r.post_hoc_flag, str) and bool(r.post_hoc_flag)
        ax.plot([x0, x1], [y0, y1], color=colors.get(r.evidence_source, "grey"), lw=0.9, alpha=0.7,
                ls=":" if flagged else "-")
    for r in e[(e["edge_type"] == "pathway_system") & e["source_node"].isin(pw_nodes)].itertuples():
        if r.target_node in pos:
            (x0, y0), (x1, y1) = pos[r.source_node], pos[r.target_node]
            ax.plot([x0, x1], [y0, y1], color="#999999", lw=0.6, alpha=0.6)
    for r in sm.itertuples():
        (x0, y0), (x1, y1) = pos[r.source_node], pos[r.target_node]
        ax.plot([x0, x1], [y0, y1], color="#444444", lw=0.8, ls="--", alpha=0.7)
    for nid, (x, y) in pos.items():
        t = lab.get(nid, nid)
        t = t if len(t) <= 55 else t[:52] + "..."
        ax.text(x, y, t, ha="center", va="center", fontsize=7 if x == 1 else 8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#bbbbbb", lw=0.5))
    for x, title in enumerate(["condition", "top enriched Reactome pathways\n(colour = evidence source)",
                               "physiological system\n(Reactome hierarchy)", "candidate measurement class\n(curated, dashed)"]):
        ax.text(x, 1.02, title, ha="center", va="bottom", fontsize=9, weight="bold")
    short = {"ot_genetic": "Open Targets genetic association", "gwas_catalog": "GWAS Catalog mapped genes",
             "ot_literature": "Open Targets literature co-mention (Europe PMC)",
             "ot_other": "Open Targets other datatypes (mostly drug targets here)",
             "mapmecfs": "mapMECFS published analytes (nominal p)"}
    for s, col in colors.items():
        ax.plot([], [], color=col, label=short[s])
    ax.plot([], [], color="black", ls=":", label="post-hoc flag: >= 50 % of the pathway genes are drug-target-only "
                                                 "(one drug mechanism can list a gene family)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=7, frameon=False)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(-0.02, 1.08)
    ax.axis("off")
    ax.set_title("Figure 5 (draft): condition -> molecular evidence -> physiological system -> measurement class "
                 "(demo cluster; condition-level molecular enrichment)", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "svg"):
        p = DRAFTS / f"fig5_graph_demo.{ext}"
        fig.savefig(p, dpi=200)
        out.append(p)
    plt.close(fig)
    return out


# =========================================================================== write-up
def _fmt(x, nd=3) -> str:
    if x is None or x is pd.NA or (isinstance(x, float) and np.isnan(x)):
        return "n/a"
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, float) and x.is_integer() and abs(x) < 1e7 and nd <= 3:
        return str(int(x))
    if isinstance(x, float) and abs(x) < 1e-3 and x != 0:
        return f"{x:.1e}"
    return f"{x:.{nd}f}" if isinstance(x, float) else str(x)


def _md(df: pd.DataFrame, cols: list[str] | None = None, nd: int = 3) -> str:
    d = df[cols] if cols else df
    if d.empty:
        return "_none_\n"
    head = "| " + " | ".join(d.columns) + " |\n|" + "---|" * len(d.columns) + "\n"
    body = "".join("| " + " | ".join(_fmt(v, nd) if not isinstance(v, str) else v.replace("|", "/") for v in row)
                   + " |\n" for row in d.itertuples(index=False))
    return head + body


def write_report(mb: pd.DataFrame, nodes: pd.DataFrame, edges: pd.DataFrame, figs: list[Path]) -> Path:
    meta = _meta()
    drop = PROVENANCE_COLUMNS
    counts = _csv("test3_source_counts").drop(columns=drop)
    qc = _csv("test3_gene_mapping_qc").drop(columns=drop)
    gc = _csv("test3_gene_agreement_control").drop(columns=drop)
    gper = _csv("test3_gene_agreement_per_condition").drop(columns=drop)
    mrank = _csv("test3_mapmecfs_rank").drop(columns=drop)
    pc = _csv("test3_pathway_agreement_control").drop(columns=drop)
    se = _csv("test3_system_enrichment").drop(columns=drop)
    fc = _csv("test3_flag_concordance").drop(columns=drop)
    dis = _csv("test3_disagreements").drop(columns=drop)
    sens = _csv("test3_sensitivity").drop(columns=drop)
    top = _csv("test3_toplevel_crosscheck").drop(columns=drop)
    svc_p = TABLES / "test3_reactome_service_crosscheck.csv"
    svc = pd.read_csv(svc_p).drop(columns=drop) if svc_p.exists() else pd.DataFrame()
    ora = _csv("test3_reactome_ora")
    gpairs = _csv("test3_gene_agreement_pairs")
    L: list[str] = []
    A = L.append
    v = meta.get("versions", {})
    A("# Molecular coherence (Test 3) and the phenotype -> measurable-biology graph\n\n")
    A(f"Generated {utc_now_iso()} by `{PRODUCER}.write_report` from files written by `measure_it.omics.coherence` "
      f"(built {meta.get('built_at', UNKNOWN)}) and `{PRODUCER}`. Plan (written before the analysis): "
      "`docs/ANALYSIS_PLAN_MOLECULAR.md`. Every number below is read from `results/tables/test3_*.csv`, "
      "`results/tables/fig5_*.csv` or `data/processed/measurable_biology.parquet`.\n\n")
    A("**Scope.** Condition-level molecular enrichment only: gene lists that public databases attach to a disease "
      "concept. Nothing here is measured on, or linked to, any participant, and no result describes a person. "
      "No combined omics score is computed.\n\n")
    A("Versions: " + "; ".join(f"{k}: {x}" for k, x in v.items()) + ".\n\n")
    A("**Independent review (2026-09-23).** The gene-level and pathway-level statistics were recomputed "
      "independently and match. Corrections made in review: (1) the per-condition disagreement list used a "
      "fixed empirical-p threshold of 0.05 that cannot be reached with 16-18 cross-condition pairs, so pairs whose "
      "same-condition overlap is above every cross pair were listed as disagreements; the criterion is now "
      "p <= max(0.05, 1/(1 + n_cross)). (2) The flag-concordance result is now reported with a Holm adjustment "
      "over the two pre-specified definitions, the sources behind each concordant pair, and exploratory definitions "
      "without literature and post-hoc-flagged evidence. (3) The drug records behind Open Targets 'other' support "
      "are named. (4) Post-hoc flags are drawn in the Figure 5 drafts. (5) The query tool separates a tested "
      "negative from missing data. Details are in the plan's 'Deviations' section.\n\n")

    # ---------------- headline
    g = gc[(gc["statistic"] == "log2_fold") & (gc["status"] == "ok")].copy()
    p = pc[(pc["statistic"] == "spearman_rho_neglog10p") & (pc["status"] == "ok")].copy()
    A("## Headline results\n\n")
    A("Positive and negative results are listed together. Permutation p-values test whether sources agree more "
      "for the **same** condition than for **different** conditions (condition-label shuffle).\n\n")
    for r in g.itertuples():
        verdict = ("above the cross-condition baseline" if r.perm_p < 0.05 and r.diff_ci95_low > 0 else
                   "above the cross-condition baseline by the permutation test, but the bootstrap CI of the "
                   "difference includes 0" if r.perm_p < 0.05 else "NOT above the cross-condition baseline")
        A(f"* Gene level, {r.source_a} x {r.source_b} ({r.n_conditions} conditions): same-condition mean log2 "
          f"fold {_fmt(r.same_mean)} vs cross-condition {_fmt(r.cross_mean_all_pairs)}; difference "
          f"{_fmt(r.mean_diff_same_minus_cross)} (95% CI {_fmt(r.diff_ci95_low)} to {_fmt(r.diff_ci95_high)}); "
          f"permutation p = {_fmt(r.perm_p, 4)} (BH across pairs {_fmt(r.perm_p_bh_across_pairs, 4)}) -> agreement "
          f"**{verdict}**.\n")
    for r in p.itertuples():
        verdict = ("above the cross-condition baseline" if r.perm_p < 0.05 and r.diff_ci95_low > 0 else
                   "above the cross-condition baseline by the permutation test, but the bootstrap CI of the "
                   "difference includes 0" if r.perm_p < 0.05 else "NOT above the cross-condition baseline")
        A(f"* Pathway level, {r.source_a} x {r.source_b} ({r.n_conditions} conditions): same-condition Spearman rho "
          f"{_fmt(r.same_mean)} vs {_fmt(r.cross_mean_all_pairs)} (difference 95% CI {_fmt(r.diff_ci95_low)} to "
          f"{_fmt(r.diff_ci95_high)}); permutation p = {_fmt(r.perm_p, 4)} (BH across pairs "
          f"{_fmt(r.perm_p_bh_across_pairs, 4)}) -> **{verdict}**.\n")
    for r in mrank.itertuples():
        if r.status == "ok":
            A(f"* mapMECFS (me_cfs) x {r.other_source}: overlap {int(r.me_cfs_overlap)} genes, log2 fold "
              f"{_fmt(r.me_cfs_log2_fold)}; me_cfs ranks {int(r.me_cfs_rank_of_n)} of {int(r.n_conditions_compared)} "
              f"conditions (share of conditions matching at least as well = {_fmt(r.rank_p)}; best match: "
              f"{r.best_condition}) -> me_cfs is "
              f"{'the best-matching condition' if int(r.me_cfs_rank_of_n) == 1 else '**not** the best-matching condition'}.\n")
    for r in fc.itertuples():
        holm_txt = (f"; Holm-adjusted over the two pre-specified definitions {_fmt(r.perm_p_holm_prespecified, 4)}"
                    if not pd.isna(r.perm_p_holm_prespecified) else "")
        A(f"* Systems vs curated clinical flags ({r.observed_definition.replace('_', ' ')}; {r.definition_status}): "
          f"{int(r.n_concordant)} of {int(r.n_expected_condition_system_pairs)} expected condition-system pairs observed "
          f"(null mean {_fmt(r.perm_null_mean_sum, 1)}); permutation p = {_fmt(r.perm_p, 4)}{holm_txt}.\n")
    cs_row = fc[fc["observed_definition"] == "condition_specific"]
    if len(cs_row):
        cp_ = str(cs_row["concordant_pairs"].iloc[0]).split("; ")
        n_lit_only = sum(1 for x in cp_ if x and "(ot_literature)" in x)
        n_flag = sum(1 for x in cp_ if "[post-hoc flag]" in x)
        n_gen = sum(1 for x in cp_ if "ot_genetic" in x or "gwas_catalog" in x)
        A(f"* **Reading of the flag-concordance result (review).** The only nominally significant pre-specified "
          f"definition is the 'condition-specific' one, and 'condition-specific' here means only that the system's "
          f"fold is above the median of the other conditions (about half of all tests pass by construction). Of its "
          f"{len([x for x in cp_ if x])} concordant pairs, {n_lit_only} rest on Europe PMC co-mention alone, "
          f"{n_flag} include Open Targets 'other' support that is mostly drug targets of trialled drugs (post-hoc "
          f"flag), and {n_gen} involve genetic evidence. Literature co-mention and trialled-drug targets both reflect "
          "what is already believed clinically about a condition, which is what the curated flags also encode, so "
          "this concordance is partly circular.")
        rb = fc[fc["observed_definition"] == "condition_specific_non_literature_unflagged"]
        if len(rb):
            rb = rb.iloc[0]
            ok = rb["perm_p"] < 0.05
            A(f" Restricted to condition-specific, non-literature, unflagged support: {int(rb['n_concordant'])} of "
              f"{int(rb['n_expected_condition_system_pairs'])} expected pairs observed (null mean "
              f"{_fmt(rb['perm_null_mean_sum'], 1)}; p = {_fmt(rb['perm_p'], 4)}; exploratory, added in review)."
              + ("" if ok else " The public molecular resources therefore do **not** independently recover the "
                               "curated clinical profile."))
        fcs = sens[(sens["level"] == "system_flag_concordance") & (sens["source_a"] == "condition_specific")]
        lost = fcs[(fcs["analysis"] != "primary") & (fcs["perm_p"] >= 0.05)]
        if len(fcs) > 1:
            A(" Across the sensitivity analyses the condition-specific concordance p ranges "
              f"{_fmt(fcs['perm_p'].min(), 3)}-{_fmt(fcs['perm_p'].max(), 3)}; it is >= 0.05 in "
              + (", ".join(f"{a} (p = {_fmt(p_, 3)})" for a, p_ in zip(lost["analysis"], lost["perm_p"])) or "none")
              + " (S8 is post hoc).")
        A("\n")
    n_enr = ora[ora["enriched"] == True].groupby("source")["reactome_id"].count()  # noqa: E712
    n_sets = ora.groupby("source")["condition_id"].nunique()
    A("* Enriched Reactome pathways (FDR < 0.05, overlap >= 2) summed over conditions, by source: " +
      "; ".join(f"{s}: {int(n_enr.get(s, 0))} across {int(n_sets.get(s, 0))} condition sets" for s in C.SOURCES) + ".\n")
    sysmb_h = mb[mb["physiological_system"] != ASSAY_PLATFORM].drop_duplicates(["condition_id", "physiological_system"])
    for c in DEMO:
        d = sysmb_h[sysmb_h["condition_id"] == c]
        parts = []
        for tier in ("multi_source", "single_source", "literature_only"):
            xs = d[d["support_status"] == tier]
            if len(xs):
                parts.append(f"{tier}: " + ", ".join(
                    f"{C.SYSTEM_LABEL[x.physiological_system]}"
                    + (f" [{x.post_hoc_flags}]" if x.post_hoc_flags else "") for x in xs.itertuples()))
        robust = d[d["support_status_without_flagged_sources"].isin(["multi_source", "single_source"])]
        A(f"* Demo cluster, {c}: " + ("; ".join(parts) if parts else "no supported system") +
          f". Systems with non-literature support after removing post-hoc-flagged evidence: "
          f"{', '.join(C.SYSTEM_LABEL[x] for x in robust['physiological_system']) or 'none'}.\n")
    A("\n")

    # ---------------- source counts
    A("## 1. Source counts per condition\n\n")
    A("Gene counts are HGNC protein-coding genes in the primary scope (plan section 3). GEO/SRA counts are study "
      "metadata from text queries and are never used as gene evidence.\n\n")
    A(_md(counts, ["condition_id", "ot_unique_targets_all_ids", "ot_unique_targets_direct", "gwas_associations",
                   "gwas_studies", "geo_series_in_scope", "sra_studies", "mapmecfs_rows", "genes_ot_genetic",
                   "genes_gwas_catalog", "genes_ot_literature", "genes_ot_other", "genes_mapmecfs",
                   "n_gene_sources_nonempty"]))
    A("\nIdentifier harmonisation (distinct raw identifiers per source and their outcome):\n\n")
    A(_md(qc))
    none = counts[counts["n_gene_sources_nonempty"] == 0]["condition_id"].tolist()
    A(f"\nConditions with no protein-coding gene from any source: {', '.join(none) or 'none'} "
      "(reported as UNKNOWN / NOT AVAILABLE downstream, never imputed).\n\n")

    # ---------------- gene level
    A("## 2. Gene-level agreement (Q1)\n\n")
    A(f"Background N = {meta.get('background_size_genes', UNKNOWN)} HGNC approved protein-coding genes. "
      "log2 fold = log2((k + 0.5)/(E + 0.5)).\n\n")
    A(_md(gc[gc["statistic"] == "log2_fold"],
          ["source_a", "source_b", "n_conditions", "same_mean", "cross_mean_all_pairs", "mean_diff_same_minus_cross",
           "diff_ci95_low", "diff_ci95_high", "perm_p", "n_perm", "perm_exact", "perm_p_bh_across_pairs", "status"]))
    A("\nJaccard (size-dependent; for reading only):\n\n")
    A(_md(gc[gc["statistic"] == "jaccard"], ["source_a", "source_b", "n_conditions", "same_mean",
                                             "cross_mean_all_pairs", "perm_p"], nd=4))
    same = gpairs[gpairs["same_condition"] == True]  # noqa: E712
    A("\nSame-condition overlaps (every condition with both sets non-empty):\n\n")
    A(_md(same.sort_values(["source_a", "source_b", "condition_a"]),
          ["source_a", "source_b", "condition_a", "n_a", "n_b", "overlap", "jaccard", "log2_fold", "hypergeom_p",
           "hypergeom_fdr_same_condition"], nd=3))
    A("\nmapMECFS (me_cfs only) against each condition's set of another source:\n\n")
    A(_md(mrank))
    A("\n")

    # ---------------- pathway level
    A("## 3. Pathway-level agreement (Q2, Reactome)\n\n")
    A(f"Local over-representation analysis over {meta.get('reactome_pathways_tested', UNKNOWN)} Reactome pathways "
      f"({C.PW_MIN}-{C.PW_MAX} genes), background {meta.get('reactome_background_size', UNKNOWN)} Reactome-annotated "
      "protein-coding genes, BH-FDR per gene set.\n\n")
    tab = ora.groupby(["condition_id", "source"]).agg(genes_in_background=("set_size_in_background", "first"),
                                                     enriched_pathways=("enriched", "sum")).reset_index()
    A(_md(tab))
    A("\n")
    A(_md(pc, ["source_a", "source_b", "statistic", "n_conditions", "same_mean", "cross_mean_all_pairs",
               "mean_diff_same_minus_cross", "diff_ci95_low", "diff_ci95_high", "perm_p", "n_perm", "status"]))
    A("\nCross-checks: ")
    if len(top):
        t = top.iloc[0]
        A(f"hierarchy-derived top level agrees with the Open Targets `reactome_top_level_term` for "
          f"{int(t.n_top_level_agrees)} of {int(t.n_found_in_reactome_hierarchy)} pathway ids "
          f"({_fmt(float(t.share_agrees_of_found))}). ")
    if len(svc):
        A("Reactome AnalysisService vs local ORA:\n\n")
        cols = [c for c in ["condition_id", "source", "status", "n_genes_submitted", "n_pathways_compared",
                            "spearman_rho_p_values", "n_enriched_local", "n_fdr05_service", "n_in_both", "jaccard"]
                if c in svc.columns]
        A(_md(svc, cols))
    A("\n")

    # ---------------- systems
    A("## 4. Physiological systems (Q3)\n\n")
    A("A system is supported by a source when its gene set is enriched for the system's Reactome gene universe "
      "(BH FDR < 0.05) and at least one pathway of the system is enriched. `*` = condition-specific (above the "
      "median of the other conditions for the same source).\n\n")
    tbl = se.assign(cell=pd.Series(np.where(se["supported"], "S", "."), index=se.index)
                    + pd.Series(np.where(se["condition_specific"] == True, "*", ""), index=se.index))  # noqa: E712
    piv = tbl.pivot_table(index=["condition_id", "source"], columns="physiological_system", values="cell",
                          aggfunc="first").reindex(columns=C.SYSTEMS).reset_index()
    A(_md(piv))
    A("\nFlag concordance (expected systems from curated `condition_registry` flags vs data-supported systems):\n\n")
    A(_md(fc, ["observed_definition", "definition_status", "n_conditions", "n_expected_condition_system_pairs",
               "n_observed_condition_system_pairs", "n_concordant", "perm_null_mean_sum", "perm_p",
               "perm_p_holm_prespecified", "per_system"]))
    A("\nConcordant pairs and the sources behind them (`[post-hoc flag]` = drug-target gene family or single GWAS "
      "locus):\n\n")
    A(_md(fc, ["observed_definition", "concordant_pairs"]))
    lit = se[se["supported"] & (se["source"] == "ot_literature")]
    n_lit_imm = int(((lit["physiological_system"] == "immune")).sum())
    n_lit_sets = int(se[se["source"] == "ot_literature"]["condition_id"].nunique())
    A(f"\nSource-level bias check: Open Targets literature sets support the immune system for {n_lit_imm} of "
      f"{n_lit_sets} conditions tested.\n\n")

    # ---------------- disagreements
    A("## 5. Disagreement between sources (explicit)\n\n")
    kinds = dis["kind"].value_counts()
    A("Counts by kind: " + "; ".join(f"{k}: {int(n)}" for k, n in kinds.items()) + ".\n\n")
    pa = dis[dis["kind"] == "presence_absence"]
    A("Presence/absence disagreements (one source names genes, the other none):\n\n")
    A(_md(pa, ["condition_id", "source_a", "source_b", "n_a", "n_b", "detail"]))
    fib = dis[(dis["condition_id"] == "fibromyalgia") & (dis["kind"] == "presence_absence")
              & dis["source_a"].isin(["ot_genetic"]) & dis["source_b"].isin(["gwas_catalog"])]
    if len(fib):
        A(f"\nFibromyalgia: {fib['detail'].iloc[0]}. The GWAS Catalog audit lists possible reasons (release timing, "
          "Open Targets credible-set/L2G assignment vs the Catalog's positional mapping, sub-threshold Catalog "
          "associations); which one applies was not verified here.\n")
    gl = gper[gper["statistic"] == "log2_fold"]
    ncr = gl["n_cross"][gl["n_cross"] > 1] if "n_cross" in gl else pd.Series(dtype=float)
    rng = (f"With {int(ncr.min())}-{int(ncr.max())} cross pairs per condition the smallest attainable empirical p is "
           f"{_fmt(1 / (1 + ncr.max()), 3)}-{_fmt(1 / (1 + ncr.min()), 3)}, so a" if len(ncr) else "A")
    A("\nPer-condition specificity: the same-condition fold is compared with the cross-condition pairs that share "
      f"one side with it. {rng} pair counts as condition-specific when its p is <= max(0.05, 1/(1 + n_cross)), i.e. "
      "when the same-condition fold is above every cross pair whenever 0.05 cannot be reached (review correction: "
      "the pre-review build used a fixed 0.05 threshold, which listed pairs such as migraine and IBS Open Targets "
      "genetic x GWAS Catalog as disagreements although their same-condition fold was the highest of all pairs).\n\n")
    A("Gene-level pairs with shared genes that are not condition-specific, or with no shared genes:\n\n")
    A(_md(dis[dis["kind"].isin(["overlap_not_condition_specific", "no_shared_genes"])],
          ["condition_id", "source_a", "source_b", "n_a", "n_b", "overlap", "empirical_p_vs_cross",
           "min_attainable_empirical_p", "kind"]))
    A("\nGene-level pairs whose same-condition overlap is above every cross-condition pair:\n\n")
    A(_md(dis[dis["kind"] == "agreement"], ["condition_id", "source_a", "source_b", "n_a", "n_b", "overlap",
                                            "empirical_p_vs_cross", "min_attainable_empirical_p"]))
    A("\n")

    # ---------------- post-hoc diagnostics
    A("## 5b. Post-hoc diagnostics (NOT pre-specified; added after inspecting the results)\n\n")
    A("Two features of the evidence were noticed after the primary analysis ran: (i) positional GWAS mapping can "
      "place several genes of one locus into the same pathway; (ii) the Open Targets clinical (drug-target) datatype "
      "lists every target of a trialled drug mechanism, so one mechanism can contribute a whole gene family. The "
      "pre-specified support rule is unchanged; these diagnostics are reported beside it.\n\n")
    sup = se[se["supported"] == True]  # noqa: E712
    gw = sup[sup["source"] == "gwas_catalog"]
    if len(gw):
        A("GWAS-supported systems and the number of independent GWAS loci (positions > 1 Mb apart) behind the "
          "genes of the system's enriched pathways:\n\n")
        A(_md(gw, ["condition_id", "physiological_system", "overlap", "n_independent_gwas_loci",
                   "genes_in_enriched_pathways", "top_enriched_pathways"]))
    oo = sup[sup["source"] == "ot_other"]
    if len(oo):
        n_dt = int((oo["share_overlap_drug_target_only"] >= 0.5).sum())
        A(f"\nOpen Targets 'other'-supported condition-system links: {len(oo)}; in {n_dt} of them at least half of "
          "the genes of the system's enriched pathways carry only drug-target (clinical) evidence among the "
          "non-literature, non-genetic datatypes:\n\n")
        A(_md(oo.sort_values("share_overlap_drug_target_only", ascending=False),
              ["condition_id", "physiological_system", "overlap", "share_overlap_drug_target_only",
               "top_drug_mechanisms_behind_overlap", "top_enriched_pathways"], nd=2))
        A("\n`top_drug_mechanisms_behind_overlap` (review addition) names the Open Targets drug records whose ChEMBL "
          "mechanism targets account for most of those genes. Drug mechanisms of the same class that supply the "
          "largest share of genes to two or more flagged condition-system links:\n\n")
        fl = oo[oo["share_overlap_drug_target_only"] >= 0.5]
        rec = {}
        for x in fl.itertuples():
            first = str(x.top_drug_mechanisms_behind_overlap or "").split("; ")[0]
            if not first or first == "nan":
                continue
            mech = first.split(": ", 1)[-1].rsplit(" (", 1)[0]
            rec.setdefault(mech, []).append(f"{x.condition_id}:{x.physiological_system}")
        rt = pd.DataFrame([{"top_drug_mechanism": k, "n_links": len(v), "condition_system_links": "; ".join(sorted(v))}
                           for k, v in rec.items() if len(v) >= 2]).sort_values("n_links", ascending=False) \
            if any(len(v) >= 2 for v in rec.values()) else pd.DataFrame()
        A(_md(rt) if len(rt) else "_none_\n")
        A("\nOne drug mechanism of this kind can therefore supply the same physiological system to several "
          "conditions; this is clinical precedence (what has been trialled), not disease biology.\n")
    sysmb0 = mb[mb["physiological_system"] != ASSAY_PLATFORM].drop_duplicates(["condition_id", "physiological_system"])
    cmp_ = pd.crosstab(sysmb0["support_status"], sysmb0["support_status_without_flagged_sources"]).reset_index()
    A("\nCondition-system links by pre-specified support tier (rows) and the tier after removing flagged sources "
      "(columns; exploratory):\n\n")
    A(_md(cmp_))
    s8 = sens[sens["analysis"] == "S8_ot_other_without_drug_targets"]
    if len(s8):
        s8s = s8[s8["level"] == "system"]
        A(f"\nPost-hoc sensitivity S8 (ot_other without the drug-target datatype): supported condition-system pairs "
          f"{_fmt(float(s8s['same_mean'].iloc[0]) if len(s8s) else np.nan)} (primary "
          f"{_fmt(float(sens[(sens['analysis'] == 'primary') & (sens['level'] == 'system')]['same_mean'].iloc[0]))}); "
          "gene/pathway permutation p-values are in the sensitivity table below.\n\n")

    # ---------------- sensitivity
    A("## 6. Sensitivity analyses (plan section 9; S8 is post hoc)\n\n")
    s_ok = sens[sens["level"].isin(["gene", "pathway"])]
    piv = s_ok.pivot_table(index=["level", "source_a", "source_b"], columns="analysis", values="perm_p",
                           aggfunc="first").reset_index()
    A("Permutation p (same vs cross condition) per source pair and analysis:\n\n")
    A(_md(piv, nd=4))
    ss = sens[sens["level"] == "system"][["analysis", "n_conditions", "same_mean", "perm_p"]].rename(
        columns={"same_mean": "n_supported_condition_system_pairs", "perm_p": "flag_concordance_perm_p_any_source"})
    A("\nSupported condition-system pairs (excluding `other`) and flag concordance per analysis:\n\n")
    A(_md(ss))
    fcs = sens[sens["level"] == "system_flag_concordance"]
    if len(fcs):
        A("\nFlag-concordance permutation p for every definition in every analysis (review addition; the first two "
          "definitions are pre-specified, the rest exploratory):\n\n")
        pv = fcs.pivot_table(index="source_a", columns="analysis", values="perm_p", aggfunc="first")
        order = [d for d in ["any_source", "condition_specific", "non_literature_sources",
                             "condition_specific_non_literature", "condition_specific_non_literature_unflagged"]
                 if d in pv.index]
        A(_md(pv.reindex(order).reset_index().rename(columns={"source_a": "definition"}), nd=3))
    A("\n")

    # ---------------- measurable biology
    A("## 7. Measurable biology (`measurable_biology`)\n\n")
    sysmb = mb[mb["physiological_system"] != ASSAY_PLATFORM]
    A(f"{len(mb)} rows: {len(sysmb)} condition x system x measurement-class rows (condition -> system data-derived, "
      f"system -> class curated) and {len(mb) - len(sysmb)} study-catalogue rows (omics classes already applied "
      "to the condition).\n\n")
    st = sysmb.drop_duplicates(["condition_id", "physiological_system"]).pivot_table(
        index="condition_id", columns="physiological_system", values="support_status", aggfunc="first").reset_index()
    A("Support status per condition and system:\n\n")
    A(_md(st))
    dm = sysmb[sysmb["condition_id"].isin(DEMO) & sysmb["support_status"].isin(
        ["multi_source", "single_source", "literature_only"])]
    A("\nDemo cluster candidate measurement classes with molecular support:\n\n")
    A(_md(dm, ["condition_id", "physiological_system", "measurement_class", "support_status", "sources_supporting",
               "condition_specific", "best_system_fdr"]))
    A("\n")

    # ---------------- figure
    A("## 8. Figure 5 data\n\n")
    A(f"`results/tables/fig5_nodes.csv` ({len(nodes)} nodes: " +
      ", ".join(f"{k} {int(n)}" for k, n in nodes["node_type"].value_counts().items()) +
      f") and `results/tables/fig5_edges.csv` ({len(edges)} edges: " +
      ", ".join(f"{k} {int(n)}" for k, n in edges["edge_type"].value_counts().items()) + "). Draft figures: " +
      ", ".join(f"`{p.relative_to(RESULTS.parent)}`" for p in figs) + ".\n\n")

    # ---------------- limitations
    A("## 9. Limitations\n\n")
    for x in [
        "Open Targets genetic evidence includes GWAS credible sets built from studies that are also in the GWAS "
        "Catalog, so their agreement is partly shared data, not independent replication.",
        "Open Targets literature evidence is Europe PMC co-mention; it tracks research attention and is biased "
        "towards well-studied (often immune) genes. It is kept as its own source and flagged `literature_only`.",
        "Agreement between the other source pairs is not independent replication either. Open Targets "
        "`genetic_association` (Orphanet, ClinVar/EVA, gene burden, GWAS credible sets) and `genetic_literature` "
        "(Gene2Phenotype, ClinGen, PanelApp, UniProt literature, counted in 'other') curate overlapping monogenic "
        "gene-disease knowledge. Europe PMC co-mention can be mined from the same GWAS and drug-trial papers that "
        "the genetic and drug evidence comes from; this was not measurable here, because only the top evidence items "
        "were retrieved (3 of 667 retrieved literature items cite a GWAS Catalog paper of the same condition).",
        "Open Targets 'other' evidence for long_covid and me_cfs is the clinical (drug-target) datatype only, and "
        "for pots almost only that. ChEMBL mechanisms list whole target families, so a single trialled drug (e.g. a "
        "gabapentinoid, metformin, a sodium-channel blocker) can make a physiological system look enriched, and the "
        "same drug record then does so for every condition it was trialled in.",
        "Confidence intervals are percentile bootstraps over 5-12 conditions and are approximate (likely too "
        "narrow); exact permutation p-values with 5-8 conditions have a floor of 1/n! (e.g. 0.0083 with 5).",
        "'Condition-specific' in the system tables means only that the system's log2 fold is above the median of "
        "the other conditions for the same source; about half of all tests pass by construction.",
        "GWAS Catalog mapped genes are positional (nearest/overlapping genes) and multi-trait analyses attach the same "
        "locus to two conditions; sensitivity analyses S2/S3 restrict them.",
        "mapMECFS gene and protein lists are nominal p < 0.05 results from one small cohort with no FDR-significant "
        "row; they are included as published, with that label.",
        "POTS evidence comes from the NET-deficiency id only; MCAS has no molecular evidence in any source; "
        "post-infectious syndrome is a grouping node and PTLDS has too few genes.",
        "Physiological-system assignment uses Reactome's hierarchy plus curated sub-pathway anchors; Reactome has no "
        "autonomic top level, so autonomic/cardiac support depends entirely on the curated anchors (S7 shows the "
        "top-level-only result).",
        "System -> measurement-class links are curated assumptions with a stated rationale; they are not evidence that "
        "a measurement detects the condition.",
        "Enrichment is condition-level: it describes the biology public databases associate with a disease label, "
        "not the biology of any person, and it is not evidence of mechanism.",
    ]:
        A(f"* {x}\n")
    A("\n## Reproduce\n\n`uv run python -m measure_it.omics.coherence && uv run python -m measure_it.omics.graph`\n")
    REPORT_PATH.write_text("".join(L))
    return REPORT_PATH


# =========================================================================== run + query
def _prov_table(df: pd.DataFrame, name: str, notes: str, record) -> pd.DataFrame:
    meta = _meta()
    return add_provenance(df.reset_index(drop=True), data_layer="derived",
                          source_name=f"{PRODUCER} (Figure 5 graph data; {C.LAYER_LABEL})",
                          source_version=C.version_string(meta.get("versions", {})) or UNKNOWN,
                          retrieved_at=utc_now_iso(), evidence_type="derived_score", source_record_id=record,
                          source_geographic_resolution="none", evidence_level="condition_level_molecular_enrichment",
                          provenance_notes=notes)


def run() -> dict:
    mb = build_measurable_biology()
    write_table(mb, TABLE, producer=PRODUCER,
                description="condition x physiological system x candidate measurement class; condition->system "
                            "data-derived (Reactome enrichment of public gene sets), system->class curated")
    nodes, edges = build_fig5(mb)
    nodes = _prov_table(nodes, "fig5_nodes", "Figure 5 nodes (conditions, top genes, top enriched pathways, systems, "
                                             "measurement classes)", "node_id")
    edges = _prov_table(edges, "fig5_edges", "Figure 5 edges; link_type says data-derived vs curated",
                        lambda d: d["source_node"] + "->" + d["target_node"] + "|" + d["edge_type"] + "|"
                        + d["evidence_source"].astype(str) + "|" + d["condition_id"].astype(str))
    TABLES.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(TABLES / "fig5_nodes.csv", index=False)
    edges.to_csv(TABLES / "fig5_edges.csv", index=False)
    figs = draft_figures(mb, nodes, edges)
    rep = write_report(mb, nodes, edges, figs)
    return {"measurable_biology": len(mb), "fig5_nodes": len(nodes), "fig5_edges": len(edges),
            "figures": [str(p) for p in figs], "report": str(rep)}


def get_measurable_biology(condition: str, measurement: str | None = None, *, include_unsupported: bool = False) -> dict:
    """Candidate measurement classes with condition-level molecular support for one condition.

    Returns rows of ``measurable_biology`` (condition -> system data-derived, system -> class curated).
    `measurement` filters to one class or a bundle from configs/relevance.yaml. Returns the UNKNOWN
    sentinel with a reason when the table is missing, the condition is unknown, or no molecular
    data exist for it. Never a statement about a person.
    """
    from .query import _resolve_condition
    cid, how = _resolve_condition(condition)
    if cid is None:
        return {"query": condition, "status": UNKNOWN, "reason": how.get("reason", "condition not recognised")}
    if not table_exists(TABLE):
        return {"query": condition, "condition_id": cid, "status": UNKNOWN,
                "reason": f"{TABLE} not built; run `uv run python -m measure_it.omics.graph`"}
    mb = read_table(TABLE)
    d = mb[mb["condition_id"] == cid]
    classes = None
    if measurement:
        from ..measurements.resolve import resolve_measurement
        rm = resolve_measurement(measurement)   # the shared class / bundle / alias resolver
        if rm["status"] != "matched":
            return {"query": condition, "condition_id": cid, "measurement": measurement, "status": UNKNOWN,
                    "reason": f"measurement: {rm.get('reason')}", "candidates": rm.get("candidates", [])}
        classes = list(rm["member_ids"])
        d = d[d["measurement_class"].isin(classes)]
    if d.empty:
        return {"query": condition, "condition_id": cid, "measurement": measurement, "status": UNKNOWN,
                "reason": "no measurable_biology row for this condition/measurement"}
    sysd = d[d["physiological_system"] != ASSAY_PLATFORM]
    if len(sysd) and sysd["support_status"].isin(["no_molecular_data", "too_few_genes"]).all() and \
            (d["physiological_system"] == ASSAY_PLATFORM).sum() == 0:
        return {"query": condition, "condition_id": cid, "measurement": measurement, "status": UNKNOWN,
                "reason": f"no condition-level molecular data ({sysd['support_status'].iloc[0]})",
                "caveats": sysd["caveats"].iloc[0]}
    keep = d if include_unsupported else d[~d["support_status"].isin(
        ["not_supported", "no_molecular_data", "too_few_genes"])]
    cols = ["object_id", "physiological_system", "measurement_class", "measurement_name", "support_status",
            "post_hoc_flags", "support_status_without_flagged_sources", "sources_supporting", "condition_specific",
            "best_system_fdr", "top_pathways", "link_evidence", "mapping_rationale", "caveats"]
    recs = [{k: (None if (isinstance(v, float) and np.isnan(v)) or v is pd.NA else v) for k, v in r.items()}
            for r in keep[cols].to_dict("records")]
    # physiological-system layer status, stated separately from the study-catalogue rows (review fix:
    # a condition with only catalogue rows must not read as having molecular support)
    sys_all = mb[(mb["condition_id"] == cid) & (mb["physiological_system"] != ASSAY_PLATFORM)]
    if len(sys_all) and sys_all["support_status"].isin(["no_molecular_data", "too_few_genes"]).all():
        system_level = {"status": UNKNOWN, "reason": f"no condition -> system link can be tested "
                                                     f"({sys_all['support_status'].iloc[0].replace('_', ' ')})"}
    else:
        system_level = {"status": "tested",
                        "supported_systems": sorted(set(sys_all.loc[sys_all["support_status"].isin(
                            ["multi_source", "single_source", "literature_only"]), "physiological_system"])),
                        "not_supported_systems": sorted(set(sys_all.loc[sys_all["support_status"] == "not_supported",
                                                                        "physiological_system"]))}
    tested_negative = (not recs and len(sysd) > 0 and (sysd["support_status"] == "not_supported").any())
    status = "available" if recs else ("not_supported" if tested_negative else UNKNOWN)
    reason = "" if recs else (
        "molecular data exist and were tested, but no source supports a physiological system mapped to this "
        "measurement (a negative result, not missing data)" if tested_negative
        else "no molecularly supported system maps to this measurement")
    return {"query": condition, "condition_id": cid, "measurement": measurement, "resolved_classes": classes,
            "layer": C.LAYER_LABEL, "status": status, "reason": reason, "system_level": system_level,
            "rows": recs or UNKNOWN,
            "not_supported": sorted(set(d.loc[d["support_status"] == "not_supported", "physiological_system"])),
            "caveats": ["condition -> system links are data-derived condition-level enrichment; system -> measurement "
                        "links are curated assumptions (docs/ANALYSIS_PLAN_MOLECULAR.md).",
                        "Candidate measurable phenotype only; not evidence that the measurement detects the condition.",
                        "Rows with post_hoc_flags rest partly on drug-target gene families of trialled drugs or on one "
                        "GWAS locus; read support_status_without_flagged_sources before citing them."]}


if __name__ == "__main__":  # pragma: no cover
    from .query import exit_cleanly
    print(json.dumps(run(), indent=1, default=str))
    exit_cleanly(0)
