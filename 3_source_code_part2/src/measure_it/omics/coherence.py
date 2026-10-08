"""Molecular coherence across public sources (SPEC validation test 3).

Question: do disease-linked public molecular resources return coherent, condition-specific
biology? For every registry condition, gene lists are built from four source families
(Open Targets genetic / literature / other datatypes, GWAS Catalog mapped genes, mapMECFS
published analytes), harmonised to HGNC protein-coding symbols, and compared

* gene level   — overlap, Jaccard, hypergeometric enrichment against a stated background;
* pathway level — Reactome over-representation (local, FDR) and agreement of the enriched
  pathways; the Reactome AnalysisService is used as a cross-check;
* system level  — enriched pathways are classified into physiological systems by the Reactome
  hierarchy (anchor table in docs/ANALYSIS_PLAN_MOLECULAR.md).

Every agreement statistic is set against a cross-condition baseline (the same two sources for
*different* conditions) and a condition-label permutation test. Everything is condition-level
molecular enrichment; nothing here is measured on, or linked to, a person.

Reproduce: ``uv run python -m measure_it.omics.coherence`` (plan: docs/ANALYSIS_PLAN_MOLECULAR.md).
"""
from __future__ import annotations

import io
import itertools
import json
import math
import re
import zipfile
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .. import http
from ..config import DOCS, SEED, TABLES, UNKNOWN, utc_now_iso
from ..provenance import add_provenance
from ..store import read_table, table_exists

PRODUCER = "measure_it.omics.coherence"
PLAN_PATH = DOCS / "ANALYSIS_PLAN_MOLECULAR.md"
LAYER_LABEL = "condition-level molecular enrichment (not participant-linked; not patient multi-omics)"

# --------------------------------------------------------------------------- reference data
# Reactome is pinned to a versioned release archive (not /download/current/) and HGNC is read from data/raw/hgnc;
# both are registered sources (measure_it.omics.reference_data: MANIFEST.json, DATA_AUDIT.md, registry entry).
from . import reference_data as REF  # noqa: E402

REACTOME_RELEASE = REF.REACTOME_RELEASE
REACTOME_DL = REF.REACTOME_DL
REACTOME_GMT = REF.reactome_url("gmt")
REACTOME_PATHWAYS = REF.reactome_url("pathways")
REACTOME_RELATION = REF.reactome_url("relation")
REACTOME_ANALYSIS = "https://reactome.org/AnalysisService/identifiers/projection"  # live cross-check only (cached)
HGNC_URL = REF.HGNC_URL

# --------------------------------------------------------------------------- analysis constants (plan sections 3-6)
SOURCES = ["ot_genetic", "gwas_catalog", "ot_literature", "ot_other", "mapmecfs"]
SOURCE_LABEL = {
    "ot_genetic": "Open Targets genetic_association targets",
    "gwas_catalog": "GWAS Catalog mapped genes",
    "ot_literature": "Open Targets literature (Europe PMC co-mention) targets",
    "ot_other": "Open Targets other datatypes (drug targets, animal models, expression, genetic literature, somatic, pathway)",
    "mapmecfs": "mapMECFS / NIH PI-ME/CFS published differential genes and proteins (nominal p < 0.05)",
}
# ot_genetic and gwas_catalog share GWAS data (OT gwas_credible_sets) -> one family for "multi-source"
SOURCE_FAMILY = {"ot_genetic": "genetic", "gwas_catalog": "genetic", "ot_literature": "literature",
                 "ot_other": "ot_other", "mapmecfs": "mapmecfs"}
OT_OTHER_DATATYPES = ["clinical", "animal_model", "rna_expression", "genetic_literature", "somatic_mutation",
                      "affected_pathway"]
SYSTEMS = ["immune", "neuronal", "vascular_endothelial", "autonomic_cardiac", "mitochondrial_metabolic",
           "connective_tissue", "other"]
SYSTEM_LABEL = {"immune": "immune", "neuronal": "neuronal", "vascular_endothelial": "vascular/endothelial",
                "autonomic_cardiac": "autonomic/cardiac", "mitochondrial_metabolic": "mitochondrial/metabolic",
                "connective_tissue": "connective tissue", "other": "other"}
FLAG_TO_SYSTEM = {  # plan section 6, fixed before results
    "immune": ["infectious_or_post_infectious", "autoimmune"],
    "autonomic_cardiac": ["autonomic"],
    "vascular_endothelial": ["vascular"],
    "neuronal": ["neurologic"],
    "mitochondrial_metabolic": ["fatigue_pem"],
}
MIN_ORA_GENES = 5
PW_MIN, PW_MAX = 10, 500
FDR_ALPHA = 0.05
MIN_OVERLAP = 2
N_PERM = 20_000
EXACT_PERM_LIMIT = 50_000
N_BOOT = 2_000
MAPMECFS_PRIMARY_TABLES = ("16A", "19A", "17A", "17B")
MAPMECFS_SEX_TABLES = ("16B", "16C", "19B", "19C")
MAPMECFS_GENE_CLASSES = ("transcript", "protein_aptamer")
CROSSCHECK_SETS = [("me_cfs", "ot_genetic"), ("migraine", "gwas_catalog"), ("long_covid", "ot_literature")]
MULTI_TRAIT_RE = re.compile(r"pleiotropy|\bMTAG\b|\bor\b|and/or", re.I)
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}|^\d{1,2}:\d{2}:\d{2}$")
OT_SOURCE_DB = "Open Targets Platform"
GWAS_SOURCE_DB = "GWAS Catalog"
MAPMECFS_SOURCE_TAG = "mapmecfs_nih_pi_mecfs"


@dataclass(frozen=True)
class Settings:
    """One analysis configuration: the primary analysis or a pre-specified sensitivity (plan 9)."""
    name: str = "primary"
    ot_scope: str = "primary"               # primary: exact/narrow ids, direct targets for ids with descendants; indirect: all
    gwas_filter: str = "all"                # all | genome_wide (p < 5e-8) | single_trait
    background: str = "hgnc_protein_coding"  # | project_universe
    exclude_conditions: tuple = ()
    mapmecfs_tables: tuple = MAPMECFS_PRIMARY_TABLES
    anchors: str = "full"                   # full | top_only
    ot_other_exclude_clinical: bool = False  # post-hoc S8: drop the clinical (drug-target) datatype from ot_other
    description: str = "primary analysis (plan sections 3-6)"


PRIMARY = Settings()
SENSITIVITY = [
    Settings("S1_ot_indirect", ot_scope="indirect",
             description="Open Targets: all ids incl. broad, descendants included (API default)"),
    Settings("S2_gwas_genome_wide", gwas_filter="genome_wide", description="GWAS Catalog associations p < 5e-8 only"),
    Settings("S3_gwas_single_trait", gwas_filter="single_trait",
             description="GWAS Catalog without pleiotropy/MTAG/'or' multi-trait reported traits"),
    Settings("S4_project_background", background="project_universe",
             description="gene background = protein-coding genes named by any source for any condition"),
    Settings("S5_without_pots", exclude_conditions=("pots",), description="pots (NET-deficiency id only) excluded"),
    Settings("S6_mapmecfs_with_sex_strata", mapmecfs_tables=MAPMECFS_PRIMARY_TABLES + MAPMECFS_SEX_TABLES,
             description="mapMECFS gene set adds the sex-stratified tables SD16B/C, SD19B/C"),
    Settings("S7_top_level_anchors_only", anchors="top_only",
             description="system classification from Reactome top-level anchors only (no curated sub-pathways)"),
    Settings("S8_ot_other_without_drug_targets", ot_other_exclude_clinical=True,
             description="POST HOC (not in the plan): ot_other without the clinical (drug-target) datatype, added "
                         "after drug-mechanism target families were seen to drive ot_other support"),
]
LOCUS_WINDOW_BP = 1_000_000   # post-hoc diagnostic: GWAS positions closer than this form one locus


# =========================================================================== pure statistics
def bh_fdr(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN kept as NaN)."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return out
    q = p[ok]
    n = len(q)
    order = np.argsort(q)
    ranked = q[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def hypergeom_p(k: int, M: int, K: int, n: int) -> float:
    """One-sided P(X >= k), X ~ Hypergeometric(population M, K successes, n draws)."""
    if k <= 0:
        return 1.0
    return float(stats.hypergeom.sf(k - 1, M, K, n))


def log2_fold(k: int, expected: float) -> float:
    """Continuity-corrected log2 enrichment, log2((k + 0.5) / (E + 0.5))."""
    return float(np.log2((k + 0.5) / (expected + 0.5)))


def overlap_stats(a: set, b: set, background_size: int) -> dict:
    """Overlap statistics of two gene sets drawn from a background of the stated size."""
    na, nb = len(a), len(b)
    k = len(a & b)
    union = len(a | b)
    expected = na * nb / background_size if background_size else np.nan
    return {"n_a": na, "n_b": nb, "overlap": k, "jaccard": (k / union) if union else np.nan,
            "expected_overlap": expected, "log2_fold": log2_fold(k, expected),
            "hypergeom_p": hypergeom_p(k, background_size, na, nb) if na and nb else np.nan}


def label_permutation_test(S: np.ndarray, *, n_perm: int = N_PERM, seed: int = SEED,
                           exact_limit: int = EXACT_PERM_LIMIT) -> dict:
    """Condition-label permutation test on a square agreement matrix.

    S[i, j] = agreement of source A for condition i with source B for condition j. Observed
    statistic = mean of the diagonal (same condition); null = mean of S[i, pi(i)] over
    permutations pi of B's labels. Exact enumeration when n! <= exact_limit (p = share of
    permutations with T >= T_obs, identity included); otherwise Monte Carlo with
    p = (1 + #{T >= T_obs}) / (1 + n_perm).
    """
    S = np.asarray(S, dtype=float)
    n = S.shape[0]
    if n < 2:
        return {"n_conditions": n, "t_obs": float(np.nanmean(np.diag(S))) if n else np.nan,
                "null_mean": np.nan, "null_sd": np.nan, "p": np.nan, "n_perm": 0, "exact": False}
    t_obs = float(np.nanmean(np.diag(S)))
    rows = np.arange(n)
    if math.factorial(n) <= exact_limit:
        perms = np.array(list(itertools.permutations(range(n))))
        t = np.nanmean(S[rows, perms], axis=1)
        p = float(np.mean(t >= t_obs - 1e-12))
        exact = True
    else:
        rng = np.random.default_rng(seed)
        perms = np.array([rng.permutation(n) for _ in range(n_perm)])
        t = np.nanmean(S[rows, perms], axis=1)
        p = float((1 + np.sum(t >= t_obs - 1e-12)) / (1 + n_perm))
        exact = False
    return {"n_conditions": n, "t_obs": t_obs, "null_mean": float(np.mean(t)), "null_sd": float(np.std(t)),
            "p": p, "n_perm": int(len(perms)), "exact": exact}


def bootstrap_mean_ci(d, *, n_boot: int = N_BOOT, seed: int = SEED) -> tuple[float, float]:
    """Percentile 95 % CI of the mean of d by resampling its elements (conditions)."""
    d = np.asarray([x for x in d if not np.isnan(x)], dtype=float)
    if len(d) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    means = rng.choice(d, size=(n_boot, len(d)), replace=True).mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def same_vs_cross(S: pd.DataFrame) -> pd.DataFrame:
    """Per condition: same-condition value vs the cross-condition values that share one side.

    S is indexed by condition of source A (rows) and condition of source B (columns). For a
    condition c present on both sides: cross values = S[c, c'] and S[c', c] for c' != c.
    Returns same, cross mean, difference and the empirical p = (1 + #cross >= same)/(1 + n_cross).
    """
    rows = []
    for c in [x for x in S.index if x in S.columns]:
        same = S.loc[c, c]
        cross = [S.loc[c, x] for x in S.columns if x != c] + [S.loc[x, c] for x in S.index if x != c]
        cross = np.array([v for v in cross if not pd.isna(v)], dtype=float)
        if pd.isna(same) or not len(cross):
            continue
        emp = float((1 + np.sum(cross >= same)) / (1 + len(cross)))
        rows.append({"condition_id": c, "same": float(same), "cross_mean": float(cross.mean()),
                     "cross_median": float(np.median(cross)), "n_cross": len(cross),
                     "diff_same_minus_cross": float(same - cross.mean()),
                     "empirical_p": emp, "min_attainable_empirical_p": 1.0 / (1 + len(cross)),
                     "same_above_all_cross": bool(np.all(cross < same))})
    return pd.DataFrame(rows)


def condition_specific_empirical(empirical_p: float, n_cross: int, alpha: float = 0.05) -> bool:
    """Per-condition specificity from the empirical p of `same_vs_cross`.

    With n_cross cross-condition pairs the smallest attainable empirical p is 1/(1 + n_cross); with
    16-18 cross pairs that is 0.053-0.059, so a fixed 0.05 threshold can never be met (review
    correction, 2026-09-23). A condition counts as specific when p <= max(alpha, 1/(1 + n_cross)),
    i.e. when the same-condition value is above every cross pair sharing a side whenever alpha is
    unattainable.
    """
    if empirical_p is None or np.isnan(empirical_p) or not n_cross:
        return False
    return bool(empirical_p <= max(alpha, 1.0 / (1 + n_cross)) + 1e-12)


def ora(genes: set, library: dict[str, frozenset], background: frozenset,
        min_size: int = PW_MIN, max_size: int = PW_MAX) -> pd.DataFrame:
    """Over-representation of `genes` in each library set sized [min_size, max_size] (hypergeometric, BH)."""
    g = set(genes) & background
    M, n = len(background), len(g)
    rows = []
    for pid, members in library.items():
        K = len(members)
        if K < min_size or K > max_size:
            continue
        hit = g & members
        k = len(hit)
        rows.append({"reactome_id": pid, "pathway_size": K, "overlap": k, "set_size_in_background": n,
                     "background_size": M, "expected_overlap": K * n / M if M else np.nan,
                     "p_value": hypergeom_p(k, M, K, n), "overlap_genes": ";".join(sorted(hit))})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["fdr"] = bh_fdr(df["p_value"].to_numpy())
    df["log2_fold"] = [log2_fold(k, e) for k, e in zip(df["overlap"], df["expected_overlap"])]
    df["enriched"] = (df["fdr"] < FDR_ALPHA) & (df["overlap"] >= MIN_OVERLAP)
    return df


# =========================================================================== curated tables in the plan
def parse_markdown_table(text: str, marker: str) -> pd.DataFrame:
    """Parse the markdown table between <!-- MARKER:BEGIN --> and <!-- MARKER:END -->."""
    m = re.search(rf"<!--\s*{marker}:BEGIN\s*-->(.*?)<!--\s*{marker}:END\s*-->", text, re.S)
    if not m:
        raise ValueError(f"table {marker} not found")
    lines = [ln.strip() for ln in m.group(1).strip().splitlines() if ln.strip().startswith("|")]
    header = [c.strip() for c in lines[0].strip("|").split("|")]
    body = [ln for ln in lines[1:] if not re.fullmatch(r"\|[\s\-:|]+\|", ln)]
    rows = [[c.strip() for c in ln.strip("|").split("|")] for ln in body]
    bad = [r for r in rows if len(r) != len(header)]
    if bad:
        raise ValueError(f"{marker}: malformed rows {bad[:2]}")
    return pd.DataFrame(rows, columns=header)


@lru_cache(maxsize=1)
def _plan_text() -> str:
    return PLAN_PATH.read_text()


def system_anchor_table() -> pd.DataFrame:
    """Curated Reactome anchor -> physiological system table (plan section 5.3)."""
    t = parse_markdown_table(_plan_text(), "REACTOME_SYSTEM_ANCHORS")
    unknown = set(t["physiological_system"]) - set(SYSTEMS)
    if unknown:
        raise ValueError(f"unknown systems in anchor table: {unknown}")
    return t


def system_measurement_map() -> pd.DataFrame:
    """Curated physiological system -> candidate measurement class table (plan section 7)."""
    t = parse_markdown_table(_plan_text(), "SYSTEM_MEASUREMENT_MAP")
    unknown = set(t["physiological_system"]) - set(SYSTEMS)
    if unknown:
        raise ValueError(f"unknown systems in measurement map: {unknown}")
    return t


# =========================================================================== reference data (HGNC, Reactome)
class SymbolMapper:
    """HGNC-based harmonisation of gene symbols and Ensembl ids to approved symbols."""

    def __init__(self, hgnc: pd.DataFrame):
        h = hgnc[hgnc["status"].fillna("Approved") == "Approved"].copy()
        self.approved = {s.upper(): s for s in h["symbol"]}
        self.locus_group = dict(zip(h["symbol"], h["locus_group"]))
        self.ensembl_of = dict(zip(h["symbol"], h["ensembl_gene_id"].fillna("")))
        self.hgnc_id_of = dict(zip(h["symbol"], h["hgnc_id"]))
        self.by_ensembl = {e: s for e, s in zip(h["ensembl_gene_id"], h["symbol"]) if isinstance(e, str) and e}
        self.prev: dict[str, set] = {}
        self.alias: dict[str, set] = {}
        for sym, prev, alias in zip(h["symbol"], h["prev_symbol"], h["alias_symbol"]):
            for col, target in ((prev, self.prev), (alias, self.alias)):
                if isinstance(col, str):
                    for s in col.split("|"):
                        if s.strip():
                            target.setdefault(s.strip().upper(), set()).add(sym)
        self.protein_coding = frozenset(s for s, g in self.locus_group.items() if g == "protein-coding gene")

    def map_symbol(self, raw) -> tuple[str | None, str]:
        s = "" if raw is None or (isinstance(raw, float) and np.isnan(raw)) else str(raw).strip()
        if not s:
            return None, "empty"
        if DATE_RE.match(s):
            return None, "date_mangled"
        u = s.upper()
        if u in self.approved:
            return self.approved[u], "approved"
        for table, how in ((self.prev, "previous_symbol"), (self.alias, "alias")):
            hits = table.get(u)
            if hits:
                if len(hits) == 1:
                    return next(iter(hits)), how
                return None, f"ambiguous_{how}"
        return None, "unmapped"

    def map_ensembl(self, raw) -> tuple[str | None, str]:
        s = "" if raw is None or (isinstance(raw, float) and np.isnan(raw)) else str(raw).strip()
        if not s.startswith("ENSG"):
            return None, "no_ensembl"
        sym = self.by_ensembl.get(s.split(".")[0])
        return (sym, "ensembl") if sym else (None, "ensembl_unmapped")

    def is_protein_coding(self, sym: str | None) -> bool:
        return sym is not None and sym in self.protein_coding


@dataclass
class Reference:
    mapper: SymbolMapper
    hgnc_meta: dict
    reactome_version: str
    reactome_meta: dict
    names: dict[str, str]
    parents: dict[str, list[str]]
    top_levels: set[str]
    gmt_raw: dict[str, list[str]]
    library: dict[str, frozenset] = field(default_factory=dict)
    background: frozenset = frozenset()


def _file_text(path: Path, source_id: str) -> tuple[str, dict]:
    m = REF.file_meta(source_id, path.name)
    return path.read_text(), {"url": m["url"], "fetched_at": m["fetched_at"], "last_modified": m["last_modified"]}


@lru_cache(maxsize=1)
def load_reference() -> Reference:
    """HGNC + Reactome (pinned release) reference data from data/raw/{hgnc,reactome} (reference_data.ensure_files)."""
    paths = REF.ensure_files()
    hgnc_txt, hgnc_meta = _file_text(paths["hgnc"], REF.HGNC_SOURCE)
    hgnc = pd.read_csv(io.StringIO(hgnc_txt), sep="\t", dtype=str,
                       usecols=["hgnc_id", "symbol", "locus_group", "status", "prev_symbol", "alias_symbol",
                                "ensembl_gene_id"])
    mapper = SymbolMapper(hgnc)
    ver = str(REACTOME_RELEASE)  # pinned by the archive URL
    pw_txt, pw_meta = _file_text(paths["pathways"], REF.REACTOME_SOURCE)
    rel_txt, rel_meta = _file_text(paths["relation"], REF.REACTOME_SOURCE)
    gmt_meta = REF.file_meta(REF.REACTOME_SOURCE, paths["gmt"].name)
    z = zipfile.ZipFile(paths["gmt"])
    gmt_txt = z.read(z.namelist()[0]).decode()
    pw = pd.read_csv(io.StringIO(pw_txt), sep="\t", header=None, names=["stid", "name", "species"], dtype=str)
    pw = pw[pw["species"] == "Homo sapiens"]
    names = dict(zip(pw["stid"], pw["name"]))
    rel = pd.read_csv(io.StringIO(rel_txt), sep="\t", header=None, names=["parent", "child"], dtype=str)
    rel = rel[rel["parent"].str.startswith("R-HSA") & rel["child"].str.startswith("R-HSA")]
    parents: dict[str, list[str]] = {}
    for p, c in zip(rel["parent"], rel["child"]):
        parents.setdefault(c, []).append(p)
    top_levels = set(names) - set(parents)
    gmt_raw = {}
    for line in gmt_txt.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[1].startswith("R-HSA"):
            gmt_raw[parts[1]] = parts[2:]
    ref = Reference(mapper=mapper, hgnc_meta=hgnc_meta, reactome_version=str(ver).strip(),
                    reactome_meta={"gmt": REACTOME_GMT, "fetched_at": gmt_meta["fetched_at"], "pathways": pw_meta,
                                   "relation": rel_meta},
                    names=names, parents=parents, top_levels=top_levels, gmt_raw=gmt_raw)
    lib = {}
    for pid, genes in gmt_raw.items():
        mapped = set()
        for g in genes:
            s, _ = mapper.map_symbol(g)
            if mapper.is_protein_coding(s):
                mapped.add(s)
        if mapped:
            lib[pid] = frozenset(mapped)
    ref.library = lib
    ref.background = frozenset(set().union(*lib.values()))
    return ref


def ancestors(pid: str, parents: dict[str, list[str]], _memo: dict | None = None) -> set[str]:
    """All ancestors of a pathway in the Reactome hierarchy (not including itself)."""
    memo = _memo if _memo is not None else {}
    if pid in memo:
        return memo[pid]
    out: set[str] = set()
    for p in parents.get(pid, []):
        out.add(p)
        out |= ancestors(p, parents, memo)
    memo[pid] = out
    return out


def classify_pathways(names: dict[str, str], parents: dict[str, list[str]], top_levels: set[str],
                      anchors: pd.DataFrame, *, use_sub_anchors: bool = True) -> pd.DataFrame:
    """Assign each pathway to physiological system(s) by walking the Reactome hierarchy.

    Rule (plan 5.3): sub-pathway anchors among ancestors-or-self win; otherwise the systems of the
    top-level anchors among its ancestors-or-self; otherwise 'other'.
    """
    top_anchor = {r.reactome_id: r.physiological_system for r in anchors.itertuples() if r.level == "top"}
    sub_anchor = {r.reactome_id: r.physiological_system for r in anchors.itertuples() if r.level == "sub"}
    memo: dict = {}
    rows = []
    for pid, name in names.items():
        anc = ancestors(pid, parents, memo) | {pid}
        tops = sorted(names.get(a, a) for a in anc if a in top_levels)
        sub_hits = sorted(a for a in anc if a in sub_anchor) if use_sub_anchors else []
        if sub_hits:
            systems = sorted({sub_anchor[a] for a in sub_hits})
            rule = "sub-pathway anchor: " + "; ".join(names.get(a, a) for a in sub_hits)
        else:
            top_hits = sorted(a for a in anc if a in top_anchor)
            if top_hits:
                systems = sorted({top_anchor[a] for a in top_hits})
                rule = "top-level anchor: " + "; ".join(names.get(a, a) for a in top_hits)
            else:
                systems = ["other"]
                rule = "no anchor among ancestors (top level: " + "; ".join(tops) + ")"
        rows.append({"reactome_id": pid, "reactome_name": name, "top_levels": ";".join(tops),
                     "physiological_systems": ";".join(systems), "classification_rule": rule})
    return pd.DataFrame(rows)


def system_gene_universe(library: dict[str, frozenset], classes: pd.DataFrame) -> dict[str, frozenset]:
    """Union of genes of all pathways assigned to each system."""
    sys_of = dict(zip(classes["reactome_id"], classes["physiological_systems"].str.split(";")))
    out: dict[str, set] = {s: set() for s in SYSTEMS}
    for pid, genes in library.items():
        for s in sys_of.get(pid, ["other"]):
            out[s] |= genes
    return {s: frozenset(g) for s, g in out.items()}


# =========================================================================== gene sets
def _num(s) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def load_evidence() -> pd.DataFrame:
    return read_table("condition_molecular_evidence")


def _ot_rows(ev: pd.DataFrame, settings: Settings) -> pd.DataFrame:
    a = ev[(ev["source_database"] == OT_SOURCE_DB) & (ev["source_evidence_category"] == "overall_association")]
    if settings.ot_scope == "primary":
        scope = a["ontology_match"].isin(["exact", "narrow"])
        desc = _num(a["n_descendants_in_source"]).fillna(0) > 0
        direct = a["has_direct_association"].map(lambda v: v is True or str(v).lower() == "true")
        a = a[scope & (~desc | direct)]
    return a


def _gwas_rows(ev: pd.DataFrame, settings: Settings) -> pd.DataFrame:
    v = ev[(ev["source_database"] == GWAS_SOURCE_DB) & (ev["entity_type"] == "variant")]
    if settings.gwas_filter == "genome_wide":
        v = v[_num(v["p_value"]) < 5e-8]
    elif settings.gwas_filter == "single_trait":
        v = v[~v["reported_trait"].fillna("").str.contains(MULTI_TRAIT_RE)]
    return v


def _mapmecfs_rows(ev: pd.DataFrame, settings: Settings) -> pd.DataFrame:
    m = ev[ev["source_name"].astype(str).str.contains(MAPMECFS_SOURCE_TAG, regex=False)]
    return m[m["table_id"].isin(settings.mapmecfs_tables) & m["analyte_class"].isin(MAPMECFS_GENE_CLASSES)
             & (_num(m["p_value"]) < 0.05)]


def build_gene_sets(ev: pd.DataFrame, mapper: SymbolMapper, settings: Settings = PRIMARY
                    ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Long table (condition_id, source, gene_symbol, rank_stat, ...) + identifier-mapping QC per source.

    rank_stat is the source's own statistic used only to order genes for display (Open Targets
    datatype score, higher = stronger; GWAS / mapMECFS minimum p, lower = stronger).
    """
    rows: list[dict] = []
    qc: dict[str, dict] = {s: {} for s in SOURCES}

    def note(source, key, how):
        qc[source].setdefault(key, how)

    # ---- Open Targets
    a = _ot_rows(ev, settings)
    other = [d for d in OT_OTHER_DATATYPES if not (settings.ot_other_exclude_clinical and d == "clinical")]
    dt = {"ot_genetic": ["genetic_association"], "ot_literature": ["literature"], "ot_other": other}
    all_dt = [c for c in a.columns if c.startswith("ot_datatype_score__")]
    present_mat = a[all_dt].apply(_num).fillna(0) > 0
    names_dt = np.array([c.replace("ot_datatype_score__", "") for c in all_dt])
    a = a.assign(_dts=[";".join(names_dt[row]) for row in present_mat.to_numpy()])
    for src, dts in dt.items():
        cols = [f"ot_datatype_score__{d}" for d in dts if f"ot_datatype_score__{d}" in a.columns]
        score = a[cols].apply(_num).max(axis=1)
        sub = a[score.fillna(0) > 0].assign(_score=score)
        for (cid, eid, lab), g in sub.groupby(["condition_id", "entity_id", "entity_label"], dropna=False):
            sym, how = mapper.map_ensembl(eid)
            if sym is None:
                sym, how = mapper.map_symbol(lab)
            note(src, f"{eid}|{lab}", how if mapper.is_protein_coding(sym) or sym is None else "non_protein_coding")
            if not mapper.is_protein_coding(sym):
                continue
            present = sorted({d for x in g["_dts"] for d in x.split(";") if d})
            rows.append({"condition_id": cid, "source": src, "gene_symbol": sym, "rank_stat": float(g["_score"].max()),
                         "rank_stat_label": f"max Open Targets {'/'.join(dts)} datatype score (higher = stronger)",
                         "n_records": int(len(g)), "source_ids": ";".join(sorted(set(g["ontology_id"].astype(str)))),
                         "ot_datatypes": ";".join(sorted(present)), "gwas_loci": ""})
    # ---- GWAS Catalog
    v = _gwas_rows(ev, settings)
    for cid, g in v.groupby("condition_id"):
        best: dict[str, list] = {}
        loci: dict[str, set] = {}
        for genes, p, acc, loc in zip(g["mapped_genes"], _num(g["p_value"]), g["source_accession"], g["locations"]):
            if not isinstance(genes, str):
                continue
            for raw in {x.strip() for x in re.split(r"[;,]", genes) if x.strip()}:
                sym, how = mapper.map_symbol(raw)
                note("gwas_catalog", raw, how if mapper.is_protein_coding(sym) or sym is None else "non_protein_coding")
                if not mapper.is_protein_coding(sym):
                    continue
                b = best.setdefault(sym, [np.inf, 0, set()])
                b[0] = min(b[0], p if not np.isnan(p) else np.inf)
                b[1] += 1
                b[2].add(str(acc))
                if isinstance(loc, str):
                    loci.setdefault(sym, set()).update(x.strip() for x in re.split(r"[;,]", loc) if ":" in x)
        for sym, (p, n, accs) in best.items():
            rows.append({"condition_id": cid, "source": "gwas_catalog", "gene_symbol": sym, "rank_stat": p,
                         "rank_stat_label": "minimum curated association p-value (lower = stronger)",
                         "n_records": n, "source_ids": ";".join(sorted(accs)[:10]), "ot_datatypes": "",
                         "gwas_loci": ";".join(sorted(loci.get(sym, set())))})
    # ---- mapMECFS
    m = _mapmecfs_rows(ev, settings)
    best = {}
    for ens, symtxt, p, tid in zip(m["ensembl_gene_id"], m["gene_symbol"], _num(m["p_value"]), m["table_id"]):
        sym, how = mapper.map_ensembl(ens)
        cands = []
        if sym is not None:
            cands = [(sym, how)]
            if isinstance(symtxt, str) and DATE_RE.match(symtxt.strip()):
                note("mapmecfs", f"{symtxt}|{ens}", "date_mangled_recovered_via_ensembl")
        else:
            parts = str(symtxt).split() if isinstance(symtxt, str) and not DATE_RE.match(symtxt.strip()) else [symtxt]
            cands = [mapper.map_symbol(x) for x in parts]
        for sym, how in cands:
            key = f"{symtxt}|{ens}|{sym}"
            note("mapmecfs", key, how if mapper.is_protein_coding(sym) or sym is None else "non_protein_coding")
            if not mapper.is_protein_coding(sym):
                continue
            b = best.setdefault(sym, [np.inf, 0, set()])
            b[0] = min(b[0], p)
            b[1] += 1
            b[2].add(str(tid))
    for sym, (p, n, tids) in best.items():
        rows.append({"condition_id": "me_cfs", "source": "mapmecfs", "gene_symbol": sym, "rank_stat": p,
                     "rank_stat_label": "minimum published nominal p-value (lower = stronger; no FDR-significant rows)",
                     "n_records": n, "source_ids": ";".join("SD" + t for t in sorted(tids)), "ot_datatypes": "",
                     "gwas_loci": ""})
    gs = pd.DataFrame(rows, columns=["condition_id", "source", "gene_symbol", "rank_stat", "rank_stat_label",
                                     "n_records", "source_ids", "ot_datatypes", "gwas_loci"])
    if settings.exclude_conditions:
        gs = gs[~gs["condition_id"].isin(settings.exclude_conditions)]
    gs = gs.drop_duplicates(["condition_id", "source", "gene_symbol"])
    # deterministic row order (set iteration order is randomised per process)
    gs = gs.assign(_o=gs["source"].map({s: i for i, s in enumerate(SOURCES)})).sort_values(
        ["condition_id", "_o", "gene_symbol"]).drop(columns="_o")
    gs["ensembl_gene_id"] = gs["gene_symbol"].map(mapper.ensembl_of).fillna("")
    gs["hgnc_id"] = gs["gene_symbol"].map(mapper.hgnc_id_of).fillna("")
    qrows = []
    for src, d in qc.items():
        counts = pd.Series(list(d.values()), dtype=object).value_counts().to_dict()
        qrows.append({"source": src, "n_distinct_raw_identifiers": len(d), **{f"n_{k}": int(v) for k, v in counts.items()}})
    return gs.reset_index(drop=True), pd.DataFrame(qrows).fillna(0)


def sets_by(gs: pd.DataFrame, source: str) -> dict[str, set]:
    d = gs[gs["source"] == source]
    return {c: set(g["gene_symbol"]) for c, g in d.groupby("condition_id")}


# =========================================================================== analysis
@dataclass
class Results:
    settings: Settings
    gene_sets: pd.DataFrame
    mapping_qc: pd.DataFrame
    background_size: int
    gene_pairs: pd.DataFrame          # all (c_a, c_b) pairs per source pair
    gene_control: pd.DataFrame        # per source pair: same vs cross + permutation
    gene_per_condition: pd.DataFrame  # per condition x source pair: same, cross, empirical p
    mapmecfs_rank: pd.DataFrame
    ora: pd.DataFrame
    pathway_pairs: pd.DataFrame
    pathway_control: pd.DataFrame
    pathway_per_condition: pd.DataFrame
    classes: pd.DataFrame
    system_enrichment: pd.DataFrame
    flag_concordance: pd.DataFrame
    ora_sets: list


def _pair_matrix(sets_a: dict, sets_b: dict, N: int, universe: frozenset | None) -> pd.DataFrame:
    rows = []
    for ca, a in sets_a.items():
        for cb, b in sets_b.items():
            aa, bb = (a & universe, b & universe) if universe is not None else (a, b)
            if not aa or not bb:
                continue
            rows.append({"condition_a": ca, "condition_b": cb, "same_condition": ca == cb, **overlap_stats(aa, bb, N)})
    return pd.DataFrame(rows)


def _control_rows(pairs: pd.DataFrame, stat: str, level: str, sa: str, sb: str
                  ) -> tuple[dict, pd.DataFrame]:
    """Same-vs-cross summary + permutation test for one source pair and one statistic."""
    if pairs.empty:
        return {"level": level, "source_a": sa, "source_b": sb, "statistic": stat, "n_conditions": 0,
                "status": UNKNOWN + ": no condition has both sets"}, pd.DataFrame()
    S = pairs.pivot_table(index="condition_a", columns="condition_b", values=stat, aggfunc="first")
    both = sorted(set(S.index) & set(S.columns))
    per = same_vs_cross(S)
    if per.empty or len(both) < 2:
        base = {"level": level, "source_a": sa, "source_b": sb, "statistic": stat, "n_conditions": len(both),
                "conditions": ";".join(both),
                "status": UNKNOWN + f": {len(both)} condition(s) with both sets; a cross-condition control needs >= 2"}
        if len(both) == 1:
            c = both[0]
            base["same_mean"] = float(S.loc[c, c])
        return base, per.assign(source_a=sa, source_b=sb, statistic=stat, level=level)
    sq = S.loc[both, both].to_numpy(dtype=float)
    perm = label_permutation_test(sq)
    lo, hi = bootstrap_mean_ci(per["diff_same_minus_cross"].to_numpy())
    off = pairs[~pairs["same_condition"]][stat].astype(float)
    row = {"level": level, "source_a": sa, "source_b": sb, "statistic": stat, "n_conditions": len(both),
           "conditions": ";".join(both), "same_mean": float(np.nanmean(np.diag(sq))),
           "cross_mean_all_pairs": float(off.mean()) if len(off) else np.nan,
           "mean_diff_same_minus_cross": float(per["diff_same_minus_cross"].mean()),
           "diff_ci95_low": lo, "diff_ci95_high": hi, "perm_null_mean": perm["null_mean"],
           "perm_null_sd": perm["null_sd"], "perm_p": perm["p"], "n_perm": perm["n_perm"],
           "perm_exact": perm["exact"], "n_conditions_same_above_all_cross":
               int((per["empirical_p"] <= 1 / (1 + per["n_cross"])).sum()), "status": "ok"}
    return row, per.assign(source_a=sa, source_b=sb, statistic=stat, level=level)


def independent_loci(positions: list[str], window: int = LOCUS_WINDOW_BP) -> int:
    """Number of loci among 'chr:pos' strings, merging positions closer than `window` on a chromosome."""
    by_chr: dict[str, list[int]] = {}
    for x in positions:
        try:
            ch, pos = x.split(":")[:2]
            bp = int(float(pos))
        except ValueError:
            continue
        by_chr.setdefault(ch, []).append(bp)
    n = 0
    for pos in by_chr.values():
        pos = sorted(pos)
        n += 1 + sum(1 for a, b in zip(pos, pos[1:]) if b - a > window)
    return n


def overlap_diagnostics(genes: list[str], condition: str, source: str, info: pd.DataFrame) -> dict:
    """Post-hoc diagnostics of the genes behind a system's enriched pathways (not part of the support rule).

    gwas_catalog: number of independent GWAS loci behind the overlapping genes (positional mapping can
    put several genes of one locus, e.g. the MHC, into one pathway). ot_other: share of overlapping genes
    whose only non-literature, non-genetic Open Targets evidence is the clinical (drug-target) datatype.
    """
    out = {"n_independent_gwas_loci": np.nan, "share_overlap_drug_target_only": np.nan}
    if not genes:
        return out
    if source == "gwas_catalog":
        pos = []
        for gname in genes:
            key = (condition, source, gname)
            if key in info.index:
                pos += [x for x in str(info.loc[key, "gwas_loci"]).split(";") if ":" in x]
        out["n_independent_gwas_loci"] = independent_loci(pos) if pos else np.nan
    elif source == "ot_other":
        other = set(OT_OTHER_DATATYPES)
        n = 0
        for gname in genes:
            key = (condition, source, gname)
            if key in info.index:
                dts = set(str(info.loc[key, "ot_datatypes"]).split(";")) & other
                n += dts == {"clinical"}
        out["share_overlap_drug_target_only"] = n / len(genes)
    return out


def drug_mechanism_index(ev: pd.DataFrame) -> dict[tuple[str, str], set[str]]:
    """(condition_id, target symbol) -> {'DRUG: mechanism'} from the Open Targets drug records.

    Used only to name the trialled-drug mechanisms behind Open Targets 'other' (clinical) support
    (review addition): ChEMBL mechanisms can list a whole protein family (e.g. every voltage-gated
    calcium-channel subunit for a gabapentinoid), so one drug record can supply a system's genes.
    """
    d = ev[(ev["source_database"] == OT_SOURCE_DB) & (ev["entity_type"] == "drug")]
    out: dict[tuple[str, str], set[str]] = {}
    if d.empty or "mechanism_target_symbols" not in d:
        return out
    for c, lab, eid, moa, tg in zip(d["condition_id"], d["entity_label"], d["entity_id"], d["mechanism_of_action"],
                                    d["mechanism_target_symbols"]):
        name = lab if isinstance(lab, str) and lab else str(eid)
        m = str(moa).split(";")[0].strip() if isinstance(moa, str) else UNKNOWN
        for t in str(tg).split(";") if isinstance(tg, str) else []:
            if t.strip():
                out.setdefault((c, t.strip()), set()).add(f"{name}: {m}")
    return out


def top_drug_mechanisms(genes: list[str], condition: str, index: dict, k: int = 3) -> str:
    """The k drug records that account for most of `genes` (count of genes each record targets)."""
    cnt: dict[str, int] = {}
    for g in genes:
        for m in index.get((condition, g), ()):
            cnt[m] = cnt.get(m, 0) + 1
    top = sorted(cnt.items(), key=lambda t: (-t[1], t[0]))[:k]
    return "; ".join(f"{m} ({n} of {len(genes)} genes)" for m, n in top)


def analyze(ev: pd.DataFrame, ref: Reference, settings: Settings = PRIMARY) -> Results:
    mapper = ref.mapper
    drug_index = drug_mechanism_index(ev)
    gs, qc = build_gene_sets(ev, mapper, settings)
    if settings.background == "project_universe":
        universe = frozenset(gs["gene_symbol"])
    else:
        universe = mapper.protein_coding
    N = len(universe)
    sets = {s: sets_by(gs, s) for s in SOURCES}

    # ---------------- Q1 gene level
    pair_frames, ctrl_rows, per_frames = [], [], []
    for sa, sb in itertools.combinations(SOURCES, 2):
        pm = _pair_matrix(sets[sa], sets[sb], N, universe)
        if not pm.empty:
            pm.insert(0, "source_a", sa)
            pm.insert(1, "source_b", sb)
            pair_frames.append(pm)
        for stat in ("log2_fold", "jaccard"):
            row, per = _control_rows(pm, stat, "gene", sa, sb)
            ctrl_rows.append(row)
            if not per.empty:
                per_frames.append(per)
    gene_pairs = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()
    gene_control = pd.DataFrame(ctrl_rows)
    ok = (gene_control["status"] == "ok") & (gene_control["statistic"] == "log2_fold")
    gene_control["perm_p_bh_across_pairs"] = np.nan
    gene_control.loc[ok, "perm_p_bh_across_pairs"] = bh_fdr(gene_control.loc[ok, "perm_p"].to_numpy())
    gene_per = pd.concat(per_frames, ignore_index=True) if per_frames else pd.DataFrame()
    if not gene_pairs.empty:
        same = gene_pairs["same_condition"]
        gene_pairs["hypergeom_fdr_same_condition"] = np.nan
        gene_pairs.loc[same, "hypergeom_fdr_same_condition"] = bh_fdr(gene_pairs.loc[same, "hypergeom_p"].to_numpy())

    # mapMECFS: one condition -> rank of me_cfs among conditions of the other source
    mrows = []
    mm = sets["mapmecfs"].get("me_cfs", set())
    for sb in SOURCES[:-1]:
        vals = {c: overlap_stats(mm & universe, b & universe, N) for c, b in sets[sb].items() if b & universe}
        if not mm or "me_cfs" not in vals:
            mrows.append({"other_source": sb, "status": UNKNOWN + ": no me_cfs set in one of the sources",
                          "n_conditions_compared": len(vals)})
            continue
        f = pd.Series({c: v["log2_fold"] for c, v in vals.items()})
        mrows.append({"other_source": sb, "status": "ok", "n_conditions_compared": len(f),
                      "me_cfs_overlap": vals["me_cfs"]["overlap"], "me_cfs_log2_fold": f["me_cfs"],
                      "me_cfs_hypergeom_p": vals["me_cfs"]["hypergeom_p"],
                      "me_cfs_rank_of_n": int((f > f["me_cfs"]).sum() + 1),
                      "rank_p": float((f >= f["me_cfs"]).sum() / len(f)),
                      "other_conditions_log2_fold_median": float(f.drop("me_cfs").median()) if len(f) > 1 else np.nan,
                      "best_condition": f.idxmax()})
    mapmecfs_rank = pd.DataFrame(mrows)

    # ---------------- Q2 pathway level (ORA)
    lib, bg = ref.library, ref.background
    anchors = system_anchor_table()
    classes = classify_pathways(ref.names, ref.parents, ref.top_levels, anchors,
                                use_sub_anchors=(settings.anchors == "full"))
    ora_frames, ora_sets = [], []
    for s in SOURCES:
        for c, g in sets[s].items():
            n_bg = len(g & bg)
            if n_bg < MIN_ORA_GENES:
                continue
            o = ora(g, lib, bg)
            o.insert(0, "condition_id", c)
            o.insert(1, "source", s)
            ora_frames.append(o)
            ora_sets.append((c, s))
    ora_df = pd.concat(ora_frames, ignore_index=True) if ora_frames else pd.DataFrame()
    sys_of = dict(zip(classes["reactome_id"], classes["physiological_systems"]))
    if not ora_df.empty:
        ora_df["reactome_name"] = ora_df["reactome_id"].map(ref.names)
        ora_df["physiological_systems"] = ora_df["reactome_id"].map(sys_of).fillna("other")
    # agreement matrices
    pw_pairs, pw_ctrl, pw_per = [], [], []
    if not ora_df.empty:
        mlogp = ora_df.pivot_table(index="reactome_id", columns=["condition_id", "source"], values="p_value",
                                   aggfunc="first").apply(lambda x: -np.log10(x))
        enr = ora_df[ora_df["enriched"]].groupby(["condition_id", "source"])["reactome_id"].apply(set).to_dict()
        tested_n = mlogp.shape[0]
        for sa, sb in itertools.combinations(SOURCES, 2):
            ca = [c for (c, s) in mlogp.columns if s == sa]
            cb = [c for (c, s) in mlogp.columns if s == sb]
            rows = []
            for x in ca:
                for y in cb:
                    va, vb = mlogp[(x, sa)], mlogp[(y, sb)]
                    rho = stats.spearmanr(va, vb).statistic if va.nunique() > 1 and vb.nunique() > 1 else np.nan
                    ea, eb = enr.get((x, sa), set()), enr.get((y, sb), set())
                    k, un = len(ea & eb), len(ea | eb)
                    rows.append({"source_a": sa, "source_b": sb, "condition_a": x, "condition_b": y,
                                 "same_condition": x == y, "spearman_rho_neglog10p": rho,
                                 "n_enriched_a": len(ea), "n_enriched_b": len(eb), "n_enriched_both": k,
                                 "jaccard_enriched": k / un if un else np.nan,
                                 "hypergeom_p_enriched_overlap_anticonservative":
                                     hypergeom_p(k, tested_n, len(ea), len(eb)) if ea and eb else np.nan})
            pp = pd.DataFrame(rows)
            if not pp.empty:
                pw_pairs.append(pp)
            for stat in ("spearman_rho_neglog10p", "jaccard_enriched"):
                row, per = _control_rows(pp, stat, "pathway", sa, sb)
                pw_ctrl.append(row)
                if not per.empty:
                    pw_per.append(per)
    pathway_pairs = pd.concat(pw_pairs, ignore_index=True) if pw_pairs else pd.DataFrame()
    pathway_control = pd.DataFrame(pw_ctrl)
    if not pathway_control.empty:
        okp = (pathway_control["status"] == "ok") & (pathway_control["statistic"] == "spearman_rho_neglog10p")
        pathway_control["perm_p_bh_across_pairs"] = np.nan
        pathway_control.loc[okp, "perm_p_bh_across_pairs"] = bh_fdr(pathway_control.loc[okp, "perm_p"].to_numpy())
    pathway_per = pd.concat(pw_per, ignore_index=True) if pw_per else pd.DataFrame()

    # ---------------- Q3 system level
    universes = system_gene_universe(lib, classes)
    gene_info = gs.set_index(["condition_id", "source", "gene_symbol"])[["ot_datatypes", "gwas_loci"]]
    srows = []
    M = len(bg)
    for (c, s) in ora_sets:
        g = sets[s][c] & bg
        o = ora_df[(ora_df["condition_id"] == c) & (ora_df["source"] == s)]
        for sysid in SYSTEMS:
            G = universes[sysid] & bg
            k = len(g & G)
            E = len(G) * len(g) / M
            in_sys = o[o["physiological_systems"].str.split(";").map(lambda xs: sysid in xs)]
            enr_sys = in_sys[in_sys["enriched"]].sort_values("fdr")
            srows.append({"condition_id": c, "source": s, "source_family": SOURCE_FAMILY[s],
                          "physiological_system": sysid, "set_size_in_background": len(g),
                          "system_universe_size": len(G), "background_size": M, "overlap": k,
                          "expected_overlap": E, "log2_fold": log2_fold(k, E), "p_value": hypergeom_p(k, M, len(G), len(g)),
                          "n_pathways_tested_in_system": int(len(in_sys)), "n_enriched_pathways": int(len(enr_sys)),
                          "top_enriched_pathways": "; ".join(enr_sys["reactome_name"].astype(str).head(5)),
                          "top_enriched_pathway_ids": ";".join(enr_sys["reactome_id"].head(5)),
                          "best_pathway_fdr": float(enr_sys["fdr"].min()) if len(enr_sys) else np.nan,
                          "genes_in_system": ";".join(sorted(g & G)[:25]),
                          "genes_in_enriched_pathways": ";".join(sorted(
                              {x for og in enr_sys["overlap_genes"] for x in str(og).split(";") if x})[:40]),
                          **overlap_diagnostics(sorted({x for og in enr_sys["overlap_genes"]
                                                        for x in str(og).split(";") if x}), c, s, gene_info),
                          "top_drug_mechanisms_behind_overlap": top_drug_mechanisms(sorted(
                              {x for og in enr_sys["overlap_genes"] for x in str(og).split(";") if x}), c, drug_index)
                          if s == "ot_other" else ""})
    system_enrichment = pd.DataFrame(srows)
    if not system_enrichment.empty:
        system_enrichment["fdr"] = bh_fdr(system_enrichment["p_value"].to_numpy())
        system_enrichment["supported"] = (system_enrichment["fdr"] < FDR_ALPHA) & \
            (system_enrichment["n_enriched_pathways"] >= 1)
        pct, spec, med = [], [], []
        for r in system_enrichment.itertuples():
            others = system_enrichment[(system_enrichment["source"] == r.source)
                                       & (system_enrichment["physiological_system"] == r.physiological_system)
                                       & (system_enrichment["condition_id"] != r.condition_id)]["log2_fold"]
            if len(others) == 0:
                pct.append(np.nan)
                spec.append(pd.NA)
                med.append(np.nan)
            else:
                pct.append(float((others < r.log2_fold).mean()))
                spec.append(bool(r.log2_fold > others.median()))
                med.append(float(others.median()))
        system_enrichment["specificity_percentile"] = pct
        system_enrichment["other_conditions_median_log2_fold"] = med
        system_enrichment["condition_specific"] = pd.array(spec, dtype="boolean")

    flag = flag_concordance(system_enrichment)
    return Results(settings=settings, gene_sets=gs, mapping_qc=qc, background_size=N, gene_pairs=gene_pairs,
                   gene_control=gene_control, gene_per_condition=gene_per, mapmecfs_rank=mapmecfs_rank,
                   ora=ora_df, pathway_pairs=pathway_pairs, pathway_control=pathway_control,
                   pathway_per_condition=pathway_per, classes=classes, system_enrichment=system_enrichment,
                   flag_concordance=flag, ora_sets=ora_sets)


def expected_systems(registry: pd.DataFrame) -> dict[str, set]:
    """Systems expected from the curated clinical flags (plan section 6)."""
    out = {}
    for _, r in registry.iterrows():
        out[r["canonical_condition_id"]] = {s for s, flags in FLAG_TO_SYSTEM.items()
                                            if any(bool(r.get(f, False)) for f in flags)}
    return out


def post_hoc_flagged(se: pd.DataFrame) -> pd.Series:
    """Rows whose support carries a post-hoc flag (single GWAS locus; >= 50 % drug-target-only genes)."""
    def col(name):
        return pd.to_numeric(se[name], errors="coerce") if name in se else pd.Series(np.nan, index=se.index)
    ot = (se["source"] == "ot_other") & (col("share_overlap_drug_target_only") >= 0.5)
    gw = (se["source"] == "gwas_catalog") & (col("n_independent_gwas_loci") == 1)
    return (ot | gw).fillna(False).astype(bool)


# Flag-concordance definitions: the first two are pre-specified (plan section 6); the rest are exploratory.
FLAG_DEFINITIONS_PRESPECIFIED = ("any_source", "condition_specific")


def flag_concordance(system_enrichment: pd.DataFrame, registry: pd.DataFrame | None = None) -> pd.DataFrame:
    """Concordance of data-supported systems with curated-flag systems vs a shuffled-profile null.

    `any_source` and `condition_specific` are pre-specified. `non_literature_sources` was added by the
    analyst after the results (exploratory). `condition_specific_non_literature` and
    `condition_specific_non_literature_unflagged` were added in review (exploratory): Europe PMC
    co-mention and the drug targets of trialled drugs both reflect what is already believed clinically
    about a condition, which is also what the curated flags encode, so concordance driven by them is
    partly circular. Holm adjustment is applied across the two pre-specified definitions.
    """
    if system_enrichment is None or system_enrichment.empty:
        return pd.DataFrame([{"observed_definition": "any", "status": UNKNOWN + ": no system enrichment"}])
    reg = registry if registry is not None else read_table("condition_registry")
    exp = expected_systems(reg)
    flagged = list(FLAG_TO_SYSTEM)
    conds = sorted(system_enrichment["condition_id"].unique())
    se = system_enrichment
    sup_ = se["supported"].astype(bool)
    spec_ = se["condition_specific"].fillna(False).astype(bool)
    nonlit = se["source_family"] != "literature"
    unflag = ~post_hoc_flagged(se)
    rows = []
    for label, mask in (("any_source", sup_),
                        ("condition_specific", sup_ & spec_),
                        ("non_literature_sources", sup_ & nonlit),
                        ("condition_specific_non_literature", sup_ & spec_ & nonlit),
                        ("condition_specific_non_literature_unflagged", sup_ & spec_ & nonlit & unflag)):
        sup = se[mask]
        obs = {c: set(sup.loc[sup["condition_id"] == c, "physiological_system"]) & set(flagged) for c in conds}
        E = np.array([[s in exp.get(c, set()) for s in flagged] for c in conds], dtype=int)
        O = np.array([[s in obs[c] for s in flagged] for c in conds], dtype=int)
        # statistic: sum_c |E_c ∩ O_c|; null: permute rows of O (condition labels)
        S = E @ O.T   # S[i, j] = |E_i ∩ O_j|
        perm = label_permutation_test(S.astype(float))
        tp = int((E & O).sum())
        conc = []
        for c in conds:
            for s in flagged:
                if s in exp.get(c, set()) and s in obs[c]:
                    x = sup[(sup["condition_id"] == c) & (sup["physiological_system"] == s)]
                    tags = [f"{r.source}{'[post-hoc flag]' if not unflag.loc[r.Index] else ''}" for r in x.itertuples()]
                    conc.append(f"{c}:{s}({'+'.join(sorted(tags))})")
        rows.append({"observed_definition": label,
                     "definition_status": "pre-specified" if label in FLAG_DEFINITIONS_PRESPECIFIED else
                     ("exploratory (analyst, post hoc)" if label == "non_literature_sources"
                      else "exploratory (added in review, post hoc)"),
                     "n_conditions": len(conds), "conditions": ";".join(conds),
                     "n_expected_condition_system_pairs": int(E.sum()),
                     "n_observed_condition_system_pairs": int(O.sum()), "n_concordant": tp,
                     "concordance_share_of_expected": tp / E.sum() if E.sum() else np.nan,
                     "concordance_share_of_observed": tp / O.sum() if O.sum() else np.nan,
                     "stat_sum_concordant": perm["t_obs"] * len(conds) if len(conds) else np.nan,
                     "perm_null_mean_sum": perm["null_mean"] * len(conds) if len(conds) else np.nan,
                     "perm_p": perm["p"], "n_perm": perm["n_perm"], "perm_exact": perm["exact"],
                     "per_system": json.dumps({s: {"expected": int(E[:, i].sum()), "observed": int(O[:, i].sum()),
                                                   "both": int((E[:, i] & O[:, i]).sum())}
                                               for i, s in enumerate(flagged)}),
                     "concordant_pairs": "; ".join(conc),
                     "status": "ok" if len(conds) >= 2 else UNKNOWN})
    out = pd.DataFrame(rows)
    pre = out["observed_definition"].isin(FLAG_DEFINITIONS_PRESPECIFIED) & out["perm_p"].notna()
    out["perm_p_holm_prespecified"] = np.nan
    if pre.any():
        out.loc[pre, "perm_p_holm_prespecified"] = holm(out.loc[pre, "perm_p"].to_numpy())
    return out


def holm(p) -> np.ndarray:
    """Holm step-down adjusted p-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (n - rank) * p[i]))
        adj[i] = running
    return adj


# =========================================================================== cross-checks
def reactome_service_crosscheck(res: Results, ref: Reference) -> pd.DataFrame:
    """Local ORA vs the Reactome AnalysisService for the pre-specified gene sets (plan 5.4)."""
    rows = []
    for c, s in CROSSCHECK_SETS:
        genes = sorted(res.gene_sets[(res.gene_sets["condition_id"] == c) & (res.gene_sets["source"] == s)]["gene_symbol"])
        base = {"condition_id": c, "source": s, "n_genes_submitted": len(genes)}
        local = res.ora[(res.ora["condition_id"] == c) & (res.ora["source"] == s)]
        if len(genes) < MIN_ORA_GENES or local.empty:
            rows.append({**base, "status": UNKNOWN + ": set too small for ORA"})
            continue
        try:
            r = http.post(REACTOME_ANALYSIS, params={"interactors": "false", "pageSize": "10000", "page": "1",
                                                     "sortBy": "ENTITIES_PVALUE", "order": "ASC",
                                                     "resource": "TOTAL", "pValue": "1", "includeDisease": "true"},
                          data="\n".join(genes).encode(), headers={"Content-Type": "text/plain"}, expect_json=True,
                          timeout=300)
            r.raise_for_status()
            js = r.json()
        except Exception as exc:  # service down -> recorded, never invented
            rows.append({**base, "status": f"{UNKNOWN}: AnalysisService error {type(exc).__name__}: {str(exc)[:150]}"})
            continue
        svc = pd.DataFrame([{"reactome_id": p["stId"], "svc_p": p["entities"]["pValue"], "svc_fdr": p["entities"]["fdr"],
                             "svc_found": p["entities"]["found"]} for p in js.get("pathways", [])])
        m = local.merge(svc, on="reactome_id", how="inner")
        rho = stats.spearmanr(m["p_value"], m["svc_p"]).statistic if len(m) > 2 else np.nan
        loc_sig = set(local.loc[local["enriched"], "reactome_id"]) & set(svc["reactome_id"])
        svc_sig = set(svc.loc[svc["svc_fdr"] < FDR_ALPHA, "reactome_id"]) & set(local["reactome_id"])
        un = loc_sig | svc_sig
        rows.append({**base, "status": "ok", "service_summary_token": js.get("summary", {}).get("token", ""),
                     "service_identifiers_not_found": js.get("identifiersNotFound", np.nan),
                     "n_pathways_local_tested": len(local), "n_pathways_service": len(svc),
                     "n_pathways_compared": len(m), "spearman_rho_p_values": rho,
                     "n_enriched_local": len(loc_sig), "n_fdr05_service": len(svc_sig),
                     "n_in_both": len(loc_sig & svc_sig), "jaccard": len(loc_sig & svc_sig) / len(un) if un else np.nan,
                     "reactome_version": ref.reactome_version})
    return pd.DataFrame(rows)


def ot_toplevel_crosscheck(ev: pd.DataFrame, classes: pd.DataFrame) -> pd.DataFrame:
    """Hierarchy-derived Reactome top level vs the top-level term Open Targets reports per pathway."""
    ot = ev[(ev["entity_type"] == "pathway") & ev["reactome_top_level_term"].notna()][
        ["entity_id", "reactome_top_level_term"]].drop_duplicates("entity_id")
    m = ot.merge(classes[["reactome_id", "top_levels"]], left_on="entity_id", right_on="reactome_id", how="left")
    m["found_in_hierarchy"] = m["top_levels"].notna()
    m["agrees"] = [bool(isinstance(t, str) and o in t.split(";")) for o, t in zip(m["reactome_top_level_term"],
                                                                                 m["top_levels"])]
    return pd.DataFrame([{"n_ot_pathway_ids": len(m), "n_found_in_reactome_hierarchy": int(m["found_in_hierarchy"].sum()),
                          "n_top_level_agrees": int(m["agrees"].sum()),
                          "share_agrees_of_found": float(m.loc[m["found_in_hierarchy"], "agrees"].mean())
                          if m["found_in_hierarchy"].any() else np.nan,
                          "examples_disagree": "; ".join(
                              f"{a}: OT '{b}' vs '{c}'" for a, b, c in
                              m.loc[m["found_in_hierarchy"] & ~m["agrees"],
                                    ["entity_id", "reactome_top_level_term", "top_levels"]].head(5).itertuples(index=False))}])


# =========================================================================== descriptive tables
def source_counts(ev: pd.DataFrame, gs: pd.DataFrame) -> pd.DataFrame:
    """Per condition: record counts of every source (raw + primary-scope gene counts)."""
    reg = read_table("condition_registry")
    geo = read_table("geo_study_catalog") if table_exists("geo_study_catalog") else pd.DataFrame()
    sra = read_table("sra_study_catalog") if table_exists("sra_study_catalog") else pd.DataFrame()
    cov = read_table("condition_molecular_coverage") if table_exists("condition_molecular_coverage") else pd.DataFrame()
    ot = ev[ev["source_database"] == OT_SOURCE_DB]
    ota = ot[ot["source_evidence_category"] == "overall_association"]
    gw = ev[ev["source_database"] == GWAS_SOURCE_DB]
    mm = ev[ev["source_name"].astype(str).str.contains(MAPMECFS_SOURCE_TAG, regex=False)]
    rows = []
    for _, r in reg.iterrows():
        c = r["canonical_condition_id"]
        a = ota[ota["condition_id"] == c]
        g = gw[gw["condition_id"] == c]
        gv = g[g["entity_type"] == "variant"]
        ge = geo[geo["condition_id"] == c] if len(geo) else geo
        sr = sra[sra["condition_id"] == c] if len(sra) else sra
        m = mm[mm["condition_id"] == c]
        cv = cov[cov["condition_id"] == c] if len(cov) else cov
        n_gene = {s: int(((gs["condition_id"] == c) & (gs["source"] == s)).sum()) for s in SOURCES}
        rows.append({
            "condition_id": c, "preferred_name": r["preferred_name"], "primary_mondo_id": r["primary_mondo_id"],
            "grouping_only": bool(r["grouping_only"]),
            "ot_association_rows_all_ids": int(len(a)), "ot_unique_targets_all_ids": int(a["entity_id"].nunique()),
            "ot_unique_targets_direct": int(a.loc[a["has_direct_association"].map(lambda v: v is True or str(v) == "True"),
                                                  "entity_id"].nunique()),
            "ot_known_drug_rows": int(((ot["condition_id"] == c) & (ot["entity_type"] == "drug")).sum()),
            "gwas_associations": int(len(gv)),
            "gwas_associations_p_lt_5e-8": int((_num(gv["p_value"]) < 5e-8).sum()),
            "gwas_studies": int((g["entity_type"] == "study").sum()),
            "geo_series_in_scope": int(ge["study_type_in_scope"].fillna(False).astype(bool).sum()) if len(ge) else 0,
            "geo_series_in_scope_naming_condition": int((ge["study_type_in_scope"].fillna(False).astype(bool)
                                                         & ge["mentions_condition"].fillna(False).astype(bool)).sum())
            if len(ge) else 0,
            "sra_studies": int(len(sr)),
            "sra_studies_naming_condition": int(sr["mentions_condition"].fillna(False).astype(bool).sum()) if len(sr) else 0,
            "mapmecfs_rows": int(len(m)),
            "mapmecfs_gene_level_rows_p_lt_05_unstratified_tables": int(
                (m["table_id"].isin(MAPMECFS_PRIMARY_TABLES) & m["analyte_class"].isin(MAPMECFS_GENE_CLASSES)
                 & (_num(m["p_value"]) < 0.05)).sum()) if len(m) else 0,
            "mapmecfs_gene_level_rows_fdr_significant": int(
                (m["analyte_class"].isin(MAPMECFS_GENE_CLASSES) & (m["fdr_significant"] == True)).sum())  # noqa: E712
            if len(m) else 0,
            **{f"genes_{s}": n_gene[s] for s in SOURCES},
            "n_gene_sources_nonempty": int(sum(v > 0 for v in n_gene.values())),
            "coverage_notes": "; ".join(sorted({f"{x.source_database}: {x.coverage_status}" for x in cv.itertuples()
                                                if x.coverage_status not in ("queried",)}))[:500] if len(cv) else "",
        })
    return pd.DataFrame(rows)


def disagreements(counts: pd.DataFrame, gene_pairs: pd.DataFrame, gene_per: pd.DataFrame,
                  system: pd.DataFrame) -> pd.DataFrame:
    """Explicit disagreement list: presence/absence, zero overlap, non-specific overlap, system discordance."""
    rows = []
    per = gene_per[(gene_per["statistic"] == "log2_fold")] if not gene_per.empty else gene_per
    for r in counts.itertuples():
        c = r.condition_id
        for sa, sb in itertools.combinations(SOURCES, 2):
            na, nb = getattr(r, f"genes_{sa}"), getattr(r, f"genes_{sb}")
            if sa == "mapmecfs" or sb == "mapmecfs":
                if c != "me_cfs":
                    continue
            if na == 0 and nb == 0:
                continue
            if (na == 0) != (nb == 0):
                have, lack = (sa, sb) if na else (sb, sa)
                rows.append({"condition_id": c, "level": "gene", "source_a": sa, "source_b": sb,
                             "kind": "presence_absence", "n_a": na, "n_b": nb,
                             "detail": f"{SOURCE_LABEL[have]} name {max(na, nb)} protein-coding genes; "
                                       f"{SOURCE_LABEL[lack]} name none"})
                continue
            pr = gene_pairs[(gene_pairs["source_a"] == sa) & (gene_pairs["source_b"] == sb)
                            & (gene_pairs["condition_a"] == c) & (gene_pairs["condition_b"] == c)]
            if pr.empty:
                continue
            k = int(pr["overlap"].iloc[0])
            pc = per[(per["source_a"] == sa) & (per["source_b"] == sb) & (per["condition_id"] == c)]
            emp = float(pc["empirical_p"].iloc[0]) if len(pc) else np.nan
            n_cross = int(pc["n_cross"].iloc[0]) if len(pc) else 0
            min_p = 1.0 / (1 + n_cross) if n_cross else np.nan
            specific = condition_specific_empirical(emp, n_cross)
            if k == 0:
                kind, detail = "no_shared_genes", f"{na} vs {nb} genes, none shared"
            elif not np.isnan(emp) and not specific:
                kind, detail = ("overlap_not_condition_specific",
                                f"{k} shared genes but the same-condition fold is not above every cross-condition "
                                f"pair sharing a side (empirical p = {emp:.3f}; smallest attainable "
                                f"{min_p:.3f} with {n_cross} cross pairs)")
            elif not np.isnan(emp):
                kind, detail = ("agreement", f"{k} shared genes; same-condition fold above all {n_cross} cross-condition "
                                             f"pairs sharing a side (empirical p = {emp:.3f}, the smallest attainable)")
            else:
                kind, detail = "agreement_uncontrolled", f"{k} shared genes; no cross-condition control (single condition)"
            rows.append({"condition_id": c, "level": "gene", "source_a": sa, "source_b": sb, "kind": kind,
                         "n_a": na, "n_b": nb, "overlap": k, "hypergeom_p": float(pr["hypergeom_p"].iloc[0]),
                         "empirical_p_vs_cross": emp, "n_cross_pairs": n_cross, "min_attainable_empirical_p": min_p,
                         "detail": detail})
    # system-level: supported by one family, tested but not supported by another
    if system is not None and not system.empty:
        for (c, sysid), g in system.groupby(["condition_id", "physiological_system"]):
            if sysid == "other":
                continue
            sup = sorted(set(g.loc[g["supported"], "source"]))
            nsup = sorted(set(g.loc[~g["supported"], "source"]))
            if sup and nsup:
                rows.append({"condition_id": c, "level": "system", "source_a": ";".join(sup), "source_b": ";".join(nsup),
                             "kind": "system_supported_by_some_sources_only",
                             "detail": f"{SYSTEM_LABEL[sysid]}: supported by {', '.join(sup)}; tested but not "
                                       f"supported by {', '.join(nsup)}"})
    return pd.DataFrame(rows)


# =========================================================================== provenance + writing
def versions(ev: pd.DataFrame, ref: Reference) -> dict:
    def first(db):
        v = ev.loc[ev["source_database"] == db, "source_version"].dropna()
        return str(v.iloc[0]) if len(v) else UNKNOWN
    mm = ev.loc[ev["source_name"].astype(str).str.contains(MAPMECFS_SOURCE_TAG, regex=False), "source_version"]
    return {"open_targets": first(OT_SOURCE_DB), "gwas_catalog": first(GWAS_SOURCE_DB),
            "mapmecfs": str(mm.iloc[0]) if len(mm) else UNKNOWN,
            "reactome": f"Reactome release {ref.reactome_version} ({REACTOME_GMT}, fetched {ref.reactome_meta['fetched_at']})",
            "hgnc": f"HGNC complete set (Last-Modified {ref.hgnc_meta['last_modified']}, fetched {ref.hgnc_meta['fetched_at']})"}


def version_string(v: dict) -> str:
    """Compact per-row source_version; the full strings live in results/tables/test3_run_metadata.json."""
    def short(k, x):
        x = str(x)
        if k == "open_targets":
            return x.split(";")[0].replace("Open Targets Platform data release", "Open Targets")
        if k == "gwas_catalog":
            return x.split(";")[0].replace("GWAS Catalog data release", "GWAS Catalog")
        if k == "mapmecfs":
            return "Walitt 2024 Nat Commun suppl. (PMID 38383456)" if x != UNKNOWN else x
        if k == "reactome":
            return x.split(" (")[0]
        if k == "hgnc":
            m = re.search(r"Last-Modified ([^,]+, [^)]*?)(?:, fetched|\))", x)
            return "HGNC " + (m.group(1) if m else x)
        return x
    return "; ".join(short(k, x) for k, x in v.items())


def with_prov(df: pd.DataFrame, *, evidence_type: str, record_id, notes: str, vers: dict,
              built_at: str) -> pd.DataFrame:
    out = add_provenance(df.reset_index(drop=True), data_layer="derived",
                         source_name=f"{PRODUCER} (Test 3; condition-level molecular enrichment)",
                         source_version=version_string(vers), retrieved_at=built_at, evidence_type=evidence_type,
                         source_record_id=record_id, source_geographic_resolution="none",
                         evidence_level="condition_level_molecular_enrichment",
                         provenance_notes=notes + " Plan: docs/ANALYSIS_PLAN_MOLECULAR.md; versions: "
                                                  "test3_run_metadata.json.")
    return out


def _write_csv(df: pd.DataFrame, name: str) -> Path:
    TABLES.mkdir(parents=True, exist_ok=True)
    p = TABLES / f"{name}.csv"
    df.to_csv(p, index=False)
    return p


def sensitivity_summary(res: Results) -> list[dict]:
    rows = []
    for level, ctrl, stat in (("gene", res.gene_control, "log2_fold"),
                              ("pathway", res.pathway_control, "spearman_rho_neglog10p")):
        if ctrl is None or ctrl.empty:
            continue
        for r in ctrl[ctrl["statistic"] == stat].itertuples():
            rows.append({"analysis": res.settings.name, "description": res.settings.description, "level": level,
                         "source_a": r.source_a, "source_b": r.source_b, "statistic": stat,
                         "n_conditions": r.n_conditions, "same_mean": getattr(r, "same_mean", np.nan),
                         "mean_diff_same_minus_cross": getattr(r, "mean_diff_same_minus_cross", np.nan),
                         "diff_ci95_low": getattr(r, "diff_ci95_low", np.nan),
                         "diff_ci95_high": getattr(r, "diff_ci95_high", np.nan),
                         "perm_p": getattr(r, "perm_p", np.nan), "status": r.status})
    se = res.system_enrichment
    fc = res.flag_concordance
    any_row = fc[fc["observed_definition"] == "any_source"]
    rows.append({"analysis": res.settings.name, "description": res.settings.description, "level": "system",
                 "statistic": "n_supported_condition_system_pairs (excluding 'other')",
                 "n_conditions": int(se["condition_id"].nunique()) if len(se) else 0,
                 "same_mean": int(se[se["supported"] & (se["physiological_system"] != "other")]
                                  .drop_duplicates(["condition_id", "physiological_system"]).shape[0]) if len(se) else 0,
                 "perm_p": float(any_row["perm_p"].iloc[0]) if len(any_row) and "perm_p" in any_row else np.nan,
                 "status": "flag-concordance permutation p in perm_p"})
    # every flag-concordance definition per analysis (review addition: robustness of the p = 0.025 claim)
    for r in fc.itertuples():
        if getattr(r, "status", "") != "ok":
            continue
        rows.append({"analysis": res.settings.name, "description": res.settings.description,
                     "level": "system_flag_concordance", "source_a": r.observed_definition, "source_b": "",
                     "statistic": "sum over conditions of |expected systems ∩ observed systems|",
                     "n_conditions": r.n_conditions, "same_mean": r.n_concordant,
                     "mean_diff_same_minus_cross": r.n_concordant - r.perm_null_mean_sum,
                     "perm_p": r.perm_p, "status": r.definition_status})
    return rows


def run(*, sensitivity: bool = True, crosscheck: bool = True) -> dict:
    built_at = utc_now_iso()
    ev = load_evidence()
    ref = load_reference()
    vers = versions(ev, ref)
    res = analyze(ev, ref, PRIMARY)
    written = {}

    def w(df, name, **kw):
        out = with_prov(df, vers=vers, built_at=built_at, **kw)
        _write_csv(out, name)
        written[name] = len(out)

    gs = res.gene_sets.copy()
    gs["object_id"] = np.where(gs["ensembl_gene_id"] != "", "gene:" + gs["ensembl_gene_id"], "")
    gs["condition_object_id"] = "condition:" + gs["condition_id"]
    w(gs, "test3_gene_sets", evidence_type="derived_score",
      record_id=lambda d: d["condition_id"] + "|" + d["source"] + "|" + d["gene_symbol"],
      notes="Test 3 gene set member (HGNC protein-coding); rank_stat orders genes for display only.")
    w(res.mapping_qc, "test3_gene_mapping_qc", evidence_type="derived_score", record_id="source",
      notes="Identifier harmonisation outcome per distinct raw identifier and source (HGNC).")
    counts = source_counts(ev, res.gene_sets)
    counts["object_id"] = "condition:" + counts["condition_id"]
    w(counts, "test3_source_counts", evidence_type="derived_score", record_id="condition_id",
      notes="Per-condition source counts; GEO/SRA counts are study metadata (never gene evidence).")
    gp = res.gene_pairs.copy()
    w(gp, "test3_gene_agreement_pairs", evidence_type="derived_score",
      record_id=lambda d: d["source_a"] + "|" + d["source_b"] + "|" + d["condition_a"] + "|" + d["condition_b"],
      notes=f"Gene-level overlap for every (source A condition, source B condition) pair; background N = "
            f"{res.background_size} HGNC approved protein-coding genes; same_condition=False rows are the "
            "cross-condition baseline.")
    w(res.gene_control, "test3_gene_agreement_control", evidence_type="derived_score",
      record_id=lambda d: d["level"] + "|" + d["source_a"] + "|" + d["source_b"] + "|" + d["statistic"],
      notes="Same-condition vs cross-condition gene agreement per source pair with condition-label permutation test "
            "and bootstrap CI over conditions.")
    w(res.gene_per_condition, "test3_gene_agreement_per_condition", evidence_type="derived_score",
      record_id=lambda d: d["source_a"] + "|" + d["source_b"] + "|" + d["condition_id"] + "|" + d["statistic"],
      notes="Per condition: same-condition agreement vs the cross-condition pairs sharing one side.")
    w(res.mapmecfs_rank, "test3_mapmecfs_rank", evidence_type="derived_score", record_id="other_source",
      notes="mapMECFS (me_cfs only) gene set compared with every condition's set of another source; rank of me_cfs.")
    ora_out = res.ora[res.ora["overlap"] > 0].copy()
    ora_out["object_id"] = "condition:" + ora_out["condition_id"]
    w(ora_out, "test3_reactome_ora", evidence_type="pathway_annotation",
      record_id=lambda d: d["condition_id"] + "|" + d["source"] + "|" + d["reactome_id"],
      notes=f"Reactome ORA, pathways {PW_MIN}-{PW_MAX} genes, background N = {len(ref.background)}; overlap > 0 rows "
            "only (overlap-0 pathways, p = 1, were in the FDR and Spearman calculations).")
    w(res.pathway_pairs, "test3_pathway_agreement_pairs", evidence_type="pathway_annotation",
      record_id=lambda d: d["source_a"] + "|" + d["source_b"] + "|" + d["condition_a"] + "|" + d["condition_b"],
      notes="Pathway-level agreement for every source-pair x condition-pair; the enriched-overlap hypergeometric p "
            "is anti-conservative (nested pathways).")
    w(res.pathway_control, "test3_pathway_agreement_control", evidence_type="pathway_annotation",
      record_id=lambda d: d["level"] + "|" + d["source_a"] + "|" + d["source_b"] + "|" + d["statistic"],
      notes="Same-condition vs cross-condition pathway agreement with condition-label permutation test.")
    w(res.pathway_per_condition, "test3_pathway_agreement_per_condition", evidence_type="pathway_annotation",
      record_id=lambda d: d["source_a"] + "|" + d["source_b"] + "|" + d["condition_id"] + "|" + d["statistic"],
      notes="Per condition pathway agreement vs cross-condition pairs.")
    w(res.classes, "test3_reactome_system_classification", evidence_type="pathway_annotation",
      record_id="reactome_id",
      notes="Every human Reactome pathway with its top-level ancestors and physiological system(s) under the anchor "
            "rule of plan 5.3 (curated anchor table).")
    se = res.system_enrichment.copy()
    se["object_id"] = "condition:" + se["condition_id"]
    w(se, "test3_system_enrichment", evidence_type="pathway_annotation",
      record_id=lambda d: d["condition_id"] + "|" + d["source"] + "|" + d["physiological_system"],
      notes="System-level hypergeometric enrichment (system gene universe within the Reactome background), BH over "
            "all tests; supported = FDR < 0.05 and >= 1 enriched pathway in the system.")
    w(res.flag_concordance, "test3_flag_concordance", evidence_type="derived_score", record_id="observed_definition",
      notes="Data-supported systems vs systems expected from curated condition_registry flags; condition-profile "
            "permutation null.")
    dis = disagreements(counts, res.gene_pairs, res.gene_per_condition, res.system_enrichment)
    w(dis, "test3_disagreements", evidence_type="derived_score",
      record_id=lambda d: d["condition_id"] + "|" + d["level"] + "|" + d["source_a"].astype(str) + "|" + d["source_b"].astype(str),
      notes="Explicit cross-source disagreements (presence/absence, zero overlap, non-specific overlap, system).")
    w(ot_toplevel_crosscheck(ev, res.classes), "test3_toplevel_crosscheck", evidence_type="pathway_annotation",
      record_id=lambda d: pd.Series(["ot_vs_reactome_hierarchy"] * len(d)),
      notes="Reactome top level from the hierarchy walk vs Open Targets reactome_top_level_term.")
    if crosscheck:
        w(reactome_service_crosscheck(res, ref), "test3_reactome_service_crosscheck", evidence_type="pathway_annotation",
          record_id=lambda d: d["condition_id"] + "|" + d["source"],
          notes="Local ORA vs Reactome AnalysisService (projection, no interactors) for the pre-specified sets.")
    sens_rows = sensitivity_summary(res)
    sens_results = {}
    if sensitivity:
        for st in SENSITIVITY:
            r2 = analyze(ev, ref, st)
            sens_results[st.name] = r2
            sens_rows += sensitivity_summary(r2)
    sens = pd.DataFrame(sens_rows)
    w(sens, "test3_sensitivity", evidence_type="derived_score",
      record_id=lambda d: d["analysis"] + "|" + d["level"] + "|" + d["source_a"].fillna("").astype(str) + "|"
      + d["source_b"].fillna("").astype(str),
      notes="Primary and pre-specified sensitivity analyses (plan section 9): permutation p per source pair and "
            "number of supported condition-system pairs.")
    meta = {"built_at": built_at, "versions": vers, "background_size_genes": res.background_size,
            "reactome_background_size": len(ref.background), "reactome_pathways_in_library": len(ref.library),
            "reactome_pathways_tested": int(res.ora["reactome_id"].nunique()) if len(res.ora) else 0,
            "ora_sets": [f"{c}|{s}" for c, s in res.ora_sets], "written": written}
    (TABLES / "test3_run_metadata.json").write_text(json.dumps(meta, indent=1, default=str))
    return {"results": res, "sensitivity": sens_results, "meta": meta}


# =========================================================================== query tool
def _read_csv(name: str) -> pd.DataFrame | None:
    p = TABLES / f"{name}.csv"
    return pd.read_csv(p) if p.exists() else None


def get_molecular_coherence(condition: str) -> dict:
    """Test 3 summary for one condition: source counts, gene agreement vs control, supported systems.

    Condition-level molecular enrichment only. Returns the UNKNOWN sentinel with a reason when the
    condition is not in the registry, the Test 3 outputs are not built, or no source has genes.
    """
    from .query import _resolve_condition
    cid, how = _resolve_condition(condition)
    if cid is None:
        return {"query": condition, "status": UNKNOWN, "reason": how.get("reason", "condition not recognised")}
    counts = _read_csv("test3_source_counts")
    if counts is None:
        return {"query": condition, "condition_id": cid, "status": UNKNOWN,
                "reason": "Test 3 outputs not built; run `uv run python -m measure_it.omics.coherence`"}
    row = counts[counts["condition_id"] == cid]
    if row.empty:
        return {"query": condition, "condition_id": cid, "status": UNKNOWN, "reason": "condition absent from Test 3"}
    r = row.iloc[0]
    out = {"query": condition, "condition_id": cid, "layer": LAYER_LABEL, "resolution": how,
           "source_counts": {k: (v.item() if hasattr(v, "item") else v) for k, v in r.items()
                             if k.startswith(("ot_", "gwas_", "geo_", "sra_", "mapmecfs_", "genes_", "n_gene"))}}
    if int(r["n_gene_sources_nonempty"]) == 0:
        out.update({"status": UNKNOWN, "reason": "no source names any protein-coding gene for this condition",
                    "coverage_notes": r.get("coverage_notes", "")})
        return out
    per = _read_csv("test3_gene_agreement_per_condition")
    pairs = _read_csv("test3_gene_agreement_pairs")
    agree = []
    if pairs is not None:
        same = pairs[(pairs["condition_a"] == cid) & (pairs["condition_b"] == cid)]
        for x in same.itertuples():
            e = per[(per["condition_id"] == cid) & (per["source_a"] == x.source_a) & (per["source_b"] == x.source_b)
                    & (per["statistic"] == "log2_fold")] if per is not None else pd.DataFrame()
            agree.append({"source_a": x.source_a, "source_b": x.source_b, "n_a": int(x.n_a), "n_b": int(x.n_b),
                          "overlap": int(x.overlap), "jaccard": float(x.jaccard), "log2_fold": float(x.log2_fold),
                          "hypergeom_p": float(x.hypergeom_p),
                          "empirical_p_vs_cross_condition": float(e["empirical_p"].iloc[0]) if len(e) else UNKNOWN,
                          "min_attainable_empirical_p": (float(1 / (1 + e["n_cross"].iloc[0])) if len(e) else UNKNOWN),
                          "same_condition_above_all_cross_pairs": (
                              condition_specific_empirical(float(e["empirical_p"].iloc[0]), int(e["n_cross"].iloc[0]))
                              if len(e) else UNKNOWN)})
    se = _read_csv("test3_system_enrichment")
    systems = UNKNOWN
    if se is not None:
        s = se[(se["condition_id"] == cid) & se["supported"]]
        flagged = post_hoc_flagged(s)
        systems = [{"physiological_system": x.physiological_system, "source": x.source, "fdr": float(x.fdr),
                    "n_enriched_pathways": int(x.n_enriched_pathways),
                    "condition_specific": (None if pd.isna(x.condition_specific) else bool(x.condition_specific)),
                    "condition_specific_definition": "system log2 fold above the median of the other conditions "
                                                     "for the same source (about half pass by construction)",
                    "post_hoc_flagged": bool(flagged.loc[x.Index]),
                    "top_drug_mechanisms_behind_overlap": (getattr(x, "top_drug_mechanisms_behind_overlap", "")
                                                           if isinstance(getattr(x, "top_drug_mechanisms_behind_overlap",
                                                                                 ""), str) else ""),
                    "top_enriched_pathways": x.top_enriched_pathways} for x in s.itertuples()] or UNKNOWN
    dis = _read_csv("test3_disagreements")
    out.update({"status": "available", "gene_agreement": agree or UNKNOWN, "supported_systems": systems,
                "disagreements": dis[(dis["condition_id"] == cid) & (dis["kind"] != "agreement")]["detail"].tolist()
                if dis is not None else UNKNOWN,
                "caveats": ["Condition-level molecular enrichment from public databases; not measured on any person.",
                            "Open Targets genetic evidence and GWAS Catalog share underlying GWAS data.",
                            "Open Targets literature evidence is text-mined co-mention (study attention).",
                            "Open Targets 'other' evidence for several conditions is mostly the targets of trialled "
                            "drugs; one drug mechanism can list a whole gene family (see post-hoc flags).",
                            "mapMECFS gene lists are nominal p < 0.05 results with no FDR-significant rows."]})
    return out


if __name__ == "__main__":  # pragma: no cover
    import sys
    from .query import exit_cleanly
    out = run(sensitivity="--no-sensitivity" not in sys.argv, crosscheck="--no-crosscheck" not in sys.argv)
    print(json.dumps(out["meta"], indent=1, default=str)[:4000])
    exit_cleanly(0)
