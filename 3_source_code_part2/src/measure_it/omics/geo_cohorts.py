"""Public GEO case/control omics cohorts for Long COVID and ME/CFS: screen, person-level classifiers, replication.

Question
--------
Do public case/control omics cohorts carry *replicable person-level* signal for Long COVID or ME/CFS?
This module answers it with a pre-specified protocol: the ``analysis_plan`` block of
configs/geo_cohorts.yaml (file last written 2026-09-24 09:33:10 -07:00, before any case/control statistic;
its ``fixed_at`` field reads 09:35+00:00, which is the local clock mislabelled as UTC). The separate
docs/ANALYSIS_PLAN_GEO_COHORTS.md was NOT written (the harness refused .md writes from the subagent).
Nulls are reported as prominently as positives.

Steps
-----
1. Screen (metadata only): every GEO series that ``geo_study_catalog`` lists for long_covid or
   me_cfs. Sample characteristics come from GEO ``acc.cgi`` (GSM text) and the file inventory from
   the GEO FTP directory listings. Eligibility rules and per-series label rules live in
   ``configs/geo_cohorts.yaml``; the code evaluates them (subject counts per group in the analysable
   subset, matrix availability) -> results/tables/geo_cohort_screen.csv.
2. Selection (pre-specified): per condition, the three largest passing blood-derived subsets
   (+ at most one replication-only subset when no two primary subsets are compatible).
3. For each selected subset: penalised logistic classifier with in-fold feature selection
   (repeated stratified, subject/pair-grouped CV), subject-level bootstrap CI, 1,000-permutation
   null, covariates-only baseline where covariates exist; per-feature differential analysis with
   BH-FDR; confound audit (label x batch/site/sex; expression-inferred sex).
4. Replication: cross-dataset transfer (train on one cohort, test on the other) for compatible
   pairs and top-100 signature overlap against a random-gene null.
5. Outputs: results/tables/geo_cohort_*.csv (results/GEO_COHORT_RESULTS.md is still to be written by hand
   from these tables), ``phenotype_signatures__geo_cohorts`` (only for cohorts whose classifier passes the
   permutation test; every row labelled dataset-specific), DATA_AUDIT.md and registry entry.

Guardrails
----------
* Person-level omics here are measured on the people of ONE GEO cohort each. They are not linked
  to any wearable participant or any other dataset: there is no cross-dataset participant key.
* A within-dataset AUROC is not a diagnostic claim. Case and control samples in public cohorts are
  often collected or processed in different batches or sites; the confound audit reports this.
* No variant search: one classifier specification, one DE model per cohort, fixed in the plan.

Reproduce: ``uv run python -m measure_it.omics.geo_cohorts`` (``--screen-only`` for step 1;
``--offline`` to use data/raw + the HTTP cache only).
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import tarfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import RAW, SEED, TABLES, UNKNOWN, load_config, raw_dir, utc_now_iso
from ..download import download_file, record_file
from ..http import get

SOURCE_ID = "geo_cohorts"
CONFIG_NAME = "geo_cohorts"
MODULE = "measure_it.omics.geo_cohorts"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_PLATFORM_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/platforms"
GEO_DOWNLOAD = "https://www.ncbi.nlm.nih.gov/geo/download/"
NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
CONDITIONS = ("long_covid", "me_cfs")
SIGNATURE_TABLE = "phenotype_signatures__geo_cohorts"
SCREEN_CSV = "geo_cohort_screen.csv"


# =========================================================================== metadata transport
def _stub(gse: str) -> str:
    return gse[:-3] + "nnn"


def gsm_text(gse: str) -> str:
    """GSM-level metadata of a series (GEO acc.cgi text, brief view); cached + archived under data/raw."""
    params = {"acc": gse, "targ": "gsm", "form": "text", "view": "brief"}
    r = get(GEO_ACC, params=params, reject_html=False)
    r.raise_for_status()
    _archive_text(f"metadata/{gse}_gsm_brief.txt", r.content, url=f"{GEO_ACC}?acc={gse}&targ=gsm&form=text&view=brief",
                  fetched_at=r.fetched_at)
    return r.text


def series_text(gse: str) -> str:
    params = {"acc": gse, "targ": "self", "form": "text", "view": "brief"}
    r = get(GEO_ACC, params=params, reject_html=False)
    r.raise_for_status()
    _archive_text(f"metadata/{gse}_series_brief.txt", r.content, url=f"{GEO_ACC}?acc={gse}&targ=self&form=text&view=brief",
                  fetched_at=r.fetched_at)
    return r.text


def _archive_text(rel: str, content: bytes, *, url: str, fetched_at: str) -> None:
    dest = raw_dir(SOURCE_ID) / rel
    if dest.exists() and dest.read_bytes() == content:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(content)
    record_file(SOURCE_ID, rel, url=url, note=f"GEO metadata text (HTTP cache fetched_at {fetched_at})")


def listing(gse: str, sub: str) -> list[dict]:
    """Files in a GEO FTP series sub-directory ('suppl' or 'matrix'): name + size string (HTML index)."""
    url = f"{GEO_FTP}/{_stub(gse)}/{gse}/{sub}/"
    r = get(url, reject_html=False, cache_errors=True)
    if r.status >= 400:
        return []
    rows = re.findall(r'<a href="([^"?/][^"]*)">[^<]*</a>\s+(\S+ \S+)\s+(\S+)', r.text)
    return [{"file": a, "modified": m, "size": s} for a, m, s in rows]


def ncbi_counts_files(gse: str) -> list[str]:
    """NCBI-generated RNA-seq count files offered for a series (empty when GEO has none)."""
    r = get(GEO_DOWNLOAD, params={"type": "rnaseq_counts", "acc": gse}, reject_html=False, cache_errors=True)
    if r.status >= 400:
        return []
    return sorted(set(re.findall(r'file=([^"&]+\.tsv\.gz)', r.text)))


# =========================================================================== metadata parsing
_KEEP = {"title", "source_name_ch1", "organism_ch1", "platform_id", "type", "library_strategy", "molecule_ch1"}


def parse_gsm_text(text: str) -> pd.DataFrame:
    """GEO GSM text -> one row per sample; characteristics become 'c:<key>' columns."""
    rows: list[dict] = []
    cur: dict | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=", 1)[1].strip()}
            rows.append(cur)
        elif cur is not None and line.startswith("!Sample_"):
            k, _, v = line[len("!Sample_"):].partition(" = ")
            if k == "characteristics_ch1":
                kk, sep, vv = v.partition(": ")
                key = f"c:{kk.strip()}" if sep else f"c:{v.strip()}"
                val = vv.strip() if sep else ""
                cur[key] = f"{cur[key]}; {val}" if key in cur else val
            elif k in _KEEP:
                cur[k] = f"{cur[k]}|{v}" if k in cur else v
    return pd.DataFrame(rows)


def _extract(df: pd.DataFrame, spec: dict) -> pd.Series:
    col = spec.get("from", "gsm")
    s = df[col].astype(str)
    if "regex" in spec:
        s = s.str.extract(spec["regex"], expand=False)
    return s


def apply_label_rule(meta: pd.DataFrame, rule: dict) -> pd.DataFrame:
    """Apply one subset rule -> analysable samples (gsm, subject, label 1=case 0=control, pair, covariates).

    Deterministic and metadata-only. Steps: subset filter -> timepoint exclusion -> label -> subject ->
    one sample per subject (earliest by `order`, else first GSM) -> twin-pair completion.
    """
    df = meta.copy()
    for col, vals in (rule.get("filter") or {}).items():
        if col not in df:
            return df.iloc[0:0].assign(subject="", label=0)
        df = df[df[col].astype(str).isin([str(v) for v in vals])]
    ex = rule.get("exclude_timepoints")
    if ex and ex["field"] in df:
        df = df[~df[ex["field"]].astype(str).isin(ex["values"])]
    lab = rule["label"]
    field = lab["field"]
    if field not in df:
        return df.iloc[0:0].assign(subject="", label=0)
    val = df[field].astype(str)
    if "case_regex" in lab:
        is_case = val.str.contains(lab["case_regex"], regex=True)
        is_ctrl = val.str.contains(lab["control_regex"], regex=True)
    else:
        is_case = val.isin(lab["case"])
        cfield = lab.get("control_field", field)
        is_ctrl = df[cfield].astype(str).isin(lab["control"]) if cfield in df else pd.Series(False, index=df.index)
    df = df[is_case ^ is_ctrl].copy()
    df["label"] = is_case.loc[df.index].astype(int)
    df["subject"] = _extract(df, rule.get("subject", {"from": "gsm"}))
    df = df[df["subject"].notna() & (df["subject"] != "nan")]
    if "order" in rule:
        o = _extract(df, rule["order"])
        df["_order"] = pd.to_numeric(o, errors="coerce") if rule["order"].get("numeric") else o
    else:
        df["_order"] = df["gsm"]
    # a subject with conflicting labels is dropped (never silently assigned)
    conflict = df.groupby("subject")["label"].nunique()
    df = df[~df["subject"].isin(conflict[conflict > 1].index)]
    df = df.sort_values(["subject", "_order", "gsm"]).groupby("subject", as_index=False).head(1)
    if "pair" in rule:
        df["pair"] = df[rule["pair"]["field"]].astype(str)
        ok = df.groupby("pair")["label"].agg(lambda x: sorted(x.tolist()) == [0, 1])
        df = df[df["pair"].isin(ok[ok].index)]
    else:
        df["pair"] = df["subject"]
    return df.drop(columns=["_order"]).sort_values("gsm").reset_index(drop=True)


# =========================================================================== screen
def catalog_universe() -> pd.DataFrame:
    from ..store import read_table
    cat = read_table("geo_study_catalog")
    cat = cat[cat["condition_id"].isin(CONDITIONS)].copy()
    return cat[["condition_id", "gse", "title", "n_samples", "study_type", "study_type_in_scope", "platform",
                "supplementary_file_types"]].sort_values(["condition_id", "gse"]).reset_index(drop=True)


def _matrix_check(gse: str, rule: dict, suppl: list[dict], mtx: list[dict], ncbi: list[str]) -> tuple[bool | None, str]:
    m = rule.get("matrix", {}) or {}
    kind = m.get("kind", "none_checked")
    names = {f["file"]: f["size"] for f in suppl + mtx}
    if kind == "none_checked":
        return None, "not checked (subset failed an earlier criterion)"
    if kind == "label_aware":
        return False, f"{m.get('file')} present but label-aware: {m.get('note', '')}"
    if kind == "suppl_listed":
        return None, m.get("note", "")
    if kind == "ncbi_counts":
        ok = m["file"] in ncbi
        return ok, f"NCBI-generated counts {'available' if ok else 'NOT available'}: {m['file']}"
    files = m.get("files") or [m.get("file")]
    missing = [f for f in files if f not in names]
    if missing:
        return False, f"missing on GEO FTP: {missing}"
    return True, "; ".join(f"{f} ({names[f]})" for f in files)


def screen(offline: bool = False) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Evaluate every catalog series against configs/geo_cohorts.yaml. Returns (screen table, analysable samples)."""
    cfg = load_config(CONFIG_NAME)
    series_cfg = cfg["series"]
    min_n = int(cfg["min_per_group"])
    blood = set(cfg["blood_compartments"])
    uni = catalog_universe()
    rows, samples = [], {}
    cat_conds = uni.groupby("gse")["condition_id"].agg(lambda x: sorted(set(x))).to_dict()
    for rec in uni.itertuples(index=False):
        gse, cond = rec.gse, rec.condition_id
        base = {"condition_id": cond, "gse": gse, "catalog_title": rec.title, "catalog_n_samples": int(rec.n_samples),
                "catalog_study_type": rec.study_type, "study_type_in_scope": bool(rec.study_type_in_scope)}
        entry = series_cfg.get(gse)
        if entry is None:
            rows.append({**base, "subset_id": gse, "decision": "fail", "fail_reasons": "NOT_IN_CONFIG",
                         "note": "series missing from configs/geo_cohorts.yaml"})
            continue
        meta = parse_gsm_text(gsm_text(gse))
        organisms = sorted(set(meta.get("organism_ch1", pd.Series(dtype=str)).dropna()))
        platforms = sorted(set(meta.get("platform_id", pd.Series(dtype=str)).dropna()))
        base.update({"geo_n_samples": len(meta), "organisms": "; ".join(organisms), "platforms": "; ".join(platforms)})
        if "exclude" in entry:
            rows.append({**base, "subset_id": gse, "decision": "fail", "fail_reasons": entry["exclude"],
                         "note": entry.get("note", "")})
            continue
        suppl, mtx = listing(gse, "suppl"), listing(gse, "matrix")
        ncbi = ncbi_counts_files(gse)
        for rule in entry["subsets"]:
            conds = {rule["condition"], rule.get("also_condition")}
            note = ""
            if cond not in conds:
                # a subset for a condition the catalog did not list this series under: report it once
                if rule["condition"] in cat_conds[gse] or cond != cat_conds[gse][0]:
                    continue
                note = f"series catalogued under {cat_conds[gse]}; subset label is {rule['condition']}"
            df = apply_label_rule(meta, rule)
            n_case = int((df["label"] == 1).sum()) if len(df) else 0
            n_ctrl = int((df["label"] == 0).sum()) if len(df) else 0
            human = organisms == ["Homo sapiens"]
            reasons = []
            if not human:
                reasons.append("NON_HUMAN_OR_MIXED_ORGANISM")
            if rule.get("design_exclusion"):
                reasons.append(rule["design_exclusion"])
            if n_case < min_n or n_ctrl < min_n:
                reasons.append("TOO_FEW")
            mat_ok, mat_note = _matrix_check(gse, rule, suppl, mtx, ncbi)
            if not reasons and mat_ok is not True:
                reasons.append("NO_USABLE_MATRIX" if mat_ok is False or mat_ok is None else "")
            is_blood = rule["compartment"] in blood
            analysed_under = rule["condition"]
            decision = "pass" if not reasons else "fail"
            row_cond = cond if cond in conds else analysed_under
            if decision == "pass" and row_cond != analysed_under:
                decision = "pass_analysed_under_other_condition"
            rows.append({**base, "condition_id": row_cond, "subset_id": rule["id"], "subset_description": rule.get("description", ""),
                         "label_rule": json.dumps(rule["label"]), "filter": json.dumps(rule.get("filter", {})),
                         "n_case_subjects": n_case, "n_control_subjects": n_ctrl, "n_analysable": n_case + n_ctrl,
                         "compartment": rule["compartment"], "blood_derived": is_blood, "assay": rule["assay"],
                         "compat_class": rule["compat_class"], "matrix_available": mat_ok, "matrix_note": mat_note,
                         "ncbi_generated_counts": "; ".join(ncbi), "covariates": "; ".join(rule.get("covariates", [])),
                         "timing": rule.get("timing", ""), "decision": decision,
                         "fail_reasons": "; ".join(r for r in reasons if r), "note": note})
            if decision == "pass":
                samples[rule["id"]] = df
    out = pd.DataFrame(rows)
    return out, samples


def select_subsets(scr: pd.DataFrame) -> pd.DataFrame:
    """Pre-specified selection: top-3 blood-derived passing subsets per condition by analysable n, + replication partner."""
    cfg = load_config(CONFIG_NAME)
    k = int(cfg["max_primary_per_condition"])
    out = []
    for cond in CONDITIONS:
        p = scr[(scr["condition_id"] == cond) & (scr["decision"] == "pass")].copy()
        p["min_group"] = p[["n_case_subjects", "n_control_subjects"]].min(axis=1)
        p["blood_derived"] = p["blood_derived"].astype(bool)
        blood = p[p["blood_derived"]].sort_values(["n_analysable", "min_group", "gse"], ascending=[False, False, True])
        nonblood = p[~p["blood_derived"]].sort_values(["n_analysable", "min_group", "gse"], ascending=[False, False, True])
        ranked = pd.concat([blood, nonblood]) if len(blood) < k else blood
        prim = ranked.head(k).assign(role="primary")
        rest = ranked.iloc[k:]
        classes = prim["compat_class"].value_counts()
        rep = []
        if not (classes > 1).any():
            for cls in prim["compat_class"]:
                cand = rest[rest["compat_class"] == cls]
                if len(cand):
                    rep.append(cand.head(1).assign(role="replication_only"))
                    break
        sel = pd.concat([prim, *rep])
        sel["selection_rank"] = range(1, len(sel) + 1)
        out.append(sel)
    return pd.concat(out, ignore_index=True)


# =========================================================================== downloads
def fetch_matrix_files(rule: dict) -> list[Path]:
    """Download the processed matrix file(s) a selected subset needs (data/raw/geo_cohorts/<GSE>/...)."""
    gse = rule["id"].split("_")[0]
    m = rule["matrix"]
    files = m.get("files") or [m.get("file")]
    paths = []
    for f in files:
        if m["kind"] == "ncbi_counts":
            url = f"{GEO_DOWNLOAD}?type=rnaseq_counts&acc={gse}&format=file&file={f}"
        elif m["kind"] == "series_matrix":
            url = f"{GEO_FTP}/{_stub(gse)}/{gse}/matrix/{f}"
        else:
            url = f"{GEO_FTP}/{_stub(gse)}/{gse}/suppl/{f}"
        paths.append(download_file(url, SOURCE_ID, f"{gse}/{f}", max_retries=4, timeout=1800))
    if m.get("platform_annot") == "GPL13534_manifest":
        paths.append(fetch_450k_manifest())
    elif m.get("platform_annot"):
        gpl = m["platform_annot"]
        url = f"{GEO_PLATFORM_FTP}/{gpl[:-3]}nnn/{gpl}/annot/{gpl}.annot.gz"
        paths.append(download_file(url, SOURCE_ID, f"platforms/{gpl}.annot.gz", max_retries=4))
    return paths


def fetch_gene_annotation() -> Path:
    """NCBI Gene human gene_info (GeneID, Symbol, chromosome, Ensembl xrefs).

    GEO's own Human.GRCh38.p13.annot.tsv.gz link (/geo/download/?format=file&type=rnaseq_counts&file=...) answered
    an automated client with a reCAPTCHA challenge page on 2026-09-24, so the NCBI Gene FTP file is used instead
    (same GeneID namespace as the NCBI-generated counts).
    """
    return download_file(NCBI_GENE_INFO_URL, SOURCE_ID, "annotation/Homo_sapiens.gene_info.gz", max_retries=4)


def fetch_450k_manifest() -> Path:
    f = "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
    url = f"{GEO_PLATFORM_FTP}/GPL13nnn/GPL13534/suppl/{f}"
    return download_file(url, SOURCE_ID, f"platforms/{f}", max_retries=4)


def rule_by_id(subset_id: str) -> dict:
    for gse, entry in load_config(CONFIG_NAME)["series"].items():
        for r in entry.get("subsets", []) or []:
            if r["id"] == subset_id:
                return r
    raise KeyError(subset_id)


def run_screen(offline: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    scr, samples = screen(offline)
    sel = select_subsets(scr)
    TABLES.mkdir(parents=True, exist_ok=True)
    scr.to_csv(TABLES / SCREEN_CSV, index=False)
    sel.to_csv(TABLES / "geo_cohort_selection.csv", index=False)
    for sid in sel["subset_id"]:
        samples[sid].to_csv(raw_dir(SOURCE_ID) / f"samples_{sid}.csv", index=False)
    return scr, sel, samples


# =========================================================================== annotation
SEX_Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y"]
AUTOSOMES = {str(i) for i in range(1, 23)}


@lru_cache(maxsize=1)
def gene_info() -> pd.DataFrame:
    p = RAW / SOURCE_ID / "annotation" / "Homo_sapiens.gene_info.gz"
    gi = pd.read_csv(p, sep="\t", usecols=["GeneID", "Symbol", "dbXrefs", "chromosome", "type_of_gene"], dtype=str)
    gi["ensembl"] = gi["dbXrefs"].str.extract(r"Ensembl:(ENSG\d+)", expand=False)
    return gi


def ensembl_map() -> pd.DataFrame:
    gi = gene_info()
    e = gi.dropna(subset=["ensembl"]).drop_duplicates("ensembl")
    return e.set_index("ensembl")[["Symbol", "chromosome"]]


def symbol_chrom() -> dict:
    gi = gene_info().drop_duplicates("Symbol")
    return dict(zip(gi["Symbol"], gi["chromosome"]))


def _is_autosomal(chrom: str | float | None) -> bool:
    """True unless the annotation places the feature on X, Y or MT (unknown chromosome is kept, per plan)."""
    if chrom is None or (isinstance(chrom, float) and np.isnan(chrom)) or str(chrom) in {"", "-", "nan", "Un"}:
        return True
    parts = set(re.split(r"[|;, ]+", str(chrom)))
    return not (parts & {"X", "Y", "MT", "chrX", "chrY", "chrM", "XY"})


# =========================================================================== dataset loaders
class Dataset:
    """One analysable subset: X (samples x features, log scale, sex-chromosome features removed), labels, groups."""

    def __init__(self, sid: str, X: pd.DataFrame, samples: pd.DataFrame, feat: pd.DataFrame, sexmat: pd.DataFrame | None,
                 kind: str):
        self.sid, self.X, self.samples, self.feat, self.sexmat, self.kind = sid, X, samples, feat, sexmat, kind
        assert list(X.index) == list(samples.index)

    @property
    def y(self) -> np.ndarray:
        return self.samples["label"].to_numpy(int)

    @property
    def groups(self) -> np.ndarray:
        return self.samples["pair"].astype(str).to_numpy()

    @property
    def paired(self) -> bool:
        return self.samples["pair"].nunique() < len(self.samples)


MATRIX_DROPS: dict[str, list[str]] = {}


def _restrict_to_matrix(sid: str, s: pd.DataFrame, present: set) -> pd.DataFrame:
    """Keep samples that the deposited matrix contains; record the ones it lacks (never imputed)."""
    missing = [g for g in s.index if g not in present]
    MATRIX_DROPS[sid] = missing
    return s.loc[[g for g in s.index if g in present]]


def _samples(sid: str) -> pd.DataFrame:
    s = pd.read_csv(raw_dir(SOURCE_ID) / f"samples_{sid}.csv", dtype=str)
    s["label"] = s["label"].astype(int)
    return s.set_index("gsm")


def _counts_to_logcpm(counts: pd.DataFrame, plan: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """counts: features x samples (non-negative). Returns (log2 CPM+1 of kept features, samples x features; CPM filter mask)."""
    counts = counts.astype(float).clip(lower=0)
    lib = counts.sum(axis=0)
    cpm = counts / lib * 1e6
    keep = (cpm >= plan["cpm_min"]).mean(axis=1) >= plan["cpm_min_fraction"]
    return np.log2(cpm.loc[keep] + 1).T, np.log2(cpm + 1).T


def _finish(sid, logx_all: pd.DataFrame, logx: pd.DataFrame, chrom: pd.Series, samples, kind) -> Dataset:
    sex_cols = [g for g in ["XIST", *SEX_Y_GENES] if g in logx_all.columns]
    sexmat = logx_all[sex_cols] if sex_cols else None
    auto = [c for c in logx.columns if _is_autosomal(chrom.get(c))]
    feat = pd.DataFrame({"feature": auto, "chromosome": [chrom.get(c) for c in auto]})
    return Dataset(sid, logx[auto], samples, feat, sexmat, kind)


def _rnaseq_by_symbol(counts: pd.DataFrame, id_kind: str) -> tuple[pd.DataFrame, pd.Series]:
    """Aggregate a count matrix (rows = ids) to gene symbols; returns (counts by symbol, chromosome per symbol)."""
    if id_kind == "ensembl":
        em = ensembl_map()
        ids = counts.index.str.replace(r"\.\d+$", "", regex=True)
        sym = pd.Series(ids, index=counts.index).map(em["Symbol"])
        sym = sym.fillna(pd.Series(ids, index=counts.index))
    elif id_kind == "geneid":
        gi = gene_info().set_index("GeneID")
        sym = pd.Series(counts.index.astype(str), index=counts.index).map(gi["Symbol"])
        sym = sym.fillna(pd.Series("GeneID:" + counts.index.astype(str), index=counts.index))
    else:
        sym = pd.Series(counts.index, index=counts.index)
    agg = counts.groupby(sym.values).sum()
    sc = symbol_chrom()
    em = ensembl_map()
    chrom = pd.Series({s: sc.get(s, em["chromosome"].get(s)) for s in agg.index})
    return agg, chrom


def load_dataset(sid: str) -> Dataset:
    plan = load_config(CONFIG_NAME)["analysis_plan"]
    s = _samples(sid)
    base = raw_dir(SOURCE_ID)
    if sid == "GSE226260_wb":
        c = pd.read_csv(base / "GSE226260/GSE226260_AdditionalSamples.rawCounts.csv.gz", index_col=0)
        c = c[s["title"].tolist()]
        c.columns = s.index
        agg, chrom = _rnaseq_by_symbol(c, "ensembl")
        logx, logall = _counts_to_logcpm(agg, plan)
        return _finish(sid, logall, logx, chrom, s, "rnaseq")
    if sid == "GSE270045":
        c = pd.read_csv(base / "GSE270045/GSE270045_LC_counts.tsv.gz", sep="\t", index_col=0)
        c = c[s["c:sample id"].tolist()]
        c.columns = s.index
        agg, chrom = _rnaseq_by_symbol(c, "symbol")
        logx, logall = _counts_to_logcpm(agg, plan)
        return _finish(sid, logall, logx, chrom, s, "rnaseq")
    if sid == "GSE293840":
        c = pd.read_csv(base / "GSE293840/GSE293840_raw_counts_all.csv.gz", index_col=0)
        num = s["subject"].str.extract(r"(\d+)$", expand=False).astype(int)
        c = c[[f"cfs_cfrna_{n}" for n in num]]
        c.columns = s.index
        agg, chrom = _rnaseq_by_symbol(c, "ensembl")
        logx, logall = _counts_to_logcpm(agg, plan)
        return _finish(sid, logall, logx, chrom, s, "rnaseq")
    if sid == "GSE227375":
        c = pd.read_csv(base / "GSE227375/GSE227375_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="\t", index_col=0)
        s = _restrict_to_matrix(sid, s, set(c.columns))
        c = c[s.index.tolist()]
        agg, chrom = _rnaseq_by_symbol(c, "geneid")
        logx, logall = _counts_to_logcpm(agg, plan)
        return _finish(sid, logall, logx, chrom, s, "rnaseq")
    if sid == "GSE275334_lc":
        x = nanostring_matrix(base / "GSE275334/GSE275334_RAW.tar", s.index.tolist())
        sc = symbol_chrom()
        chrom = pd.Series({g: sc.get(g) for g in x.columns})
        return _finish(sid, x, x, chrom, s, "nanostring")
    if sid == "GSE16059":
        x, chrom = affy_matrix(base / "GSE16059/GSE16059_series_matrix.txt.gz", base / "platforms/GPL570.annot.gz",
                               s.index.tolist())
        return _finish(sid, x, x, chrom, s, "array")
    if sid == "GSE156792":
        x, chrom = methylation_matrix(base / "GSE156792/GSE156792_series_matrix.txt.gz", s.index.tolist())
        return _finish(sid, x, x, chrom, s, "methylation")
    raise KeyError(sid)


def parse_rcc(text: str) -> pd.DataFrame:
    m = re.search(r"<Code_Summary>\s*\n(.*?)</Code_Summary>", text, re.S)
    rows = [l.split(",") for l in m.group(1).strip().splitlines()]
    df = pd.DataFrame(rows[1:], columns=rows[0])
    df["Count"] = df["Count"].astype(float)
    return df


def nanostring_matrix(tar_path: Path, gsms: list[str]) -> pd.DataFrame:
    """RCC -> positive-control then housekeeping geometric-mean scaling -> log2(count+1), endogenous genes."""
    per = {}
    with tarfile.open(tar_path) as tf:
        for mem in tf.getmembers():
            gsm = mem.name.split("_")[0]
            if gsm in gsms:
                per[gsm] = parse_rcc(gzip.decompress(tf.extractfile(mem).read()).decode("utf-8", "replace"))
    counts = pd.DataFrame({g: d.set_index("Name")["Count"] for g, d in per.items()})[gsms]
    cls = per[gsms[0]].set_index("Name")["CodeClass"]
    gm = lambda v: np.exp(np.log(np.clip(v, 1, None)).mean(axis=0))
    pos = counts.loc[cls[cls == "Positive"].index]
    f_pos = gm(pos).mean() / gm(pos)
    c1 = counts * f_pos
    hk = c1.loc[cls[cls == "Housekeeping"].index]
    f_hk = gm(hk).mean() / gm(hk)
    c2 = c1 * f_hk
    endo = cls[cls == "Endogenous"].index
    return np.log2(c2.loc[endo] + 1).T


def _read_series_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        lines = [l for l in fh if not l.startswith("!") and l.strip()]
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip('"')
    return df


def affy_matrix(sm_path: Path, annot_path: Path, gsms: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    v = _read_series_matrix(sm_path)[gsms].astype(float)
    if np.nanpercentile(v.values, 99) > 100:
        v = np.log2(v.clip(lower=1))
    with gzip.open(annot_path, "rt", errors="replace") as fh:
        lines = [l for l in fh if not l.startswith(("#", "!", "^"))]
    ann = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    ann = ann[["ID", "Gene symbol", "Chromosome annotation"]].dropna(subset=["Gene symbol"])
    ann = ann[~ann["Gene symbol"].str.contains("///")]
    ann["chrom"] = ann["Chromosome annotation"].str.extract(r"Chromosome ([0-9XYMT]+)", expand=False)
    ann = ann.set_index("ID")
    v = v.loc[v.index.intersection(ann.index)].dropna()
    v["_sym"] = ann.loc[v.index, "Gene symbol"]
    v["_mean"] = v[gsms].mean(axis=1)
    best = v.sort_values("_mean", ascending=False).drop_duplicates("_sym")
    x = best.set_index("_sym")[gsms].T
    chrom = ann.drop_duplicates("Gene symbol").set_index("Gene symbol")["chrom"]
    return x, chrom


def methylation_matrix(sm_path: Path, gsms: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    b = _read_series_matrix(sm_path)[gsms].astype(float)
    b = b[~b.index.str.startswith("rs")].dropna()
    b = b.clip(0.001, 0.999)
    m = np.log2(b / (1 - b))
    man = RAW / SOURCE_ID / "platforms" / "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
    mf = pd.read_csv(man, skiprows=7, usecols=["IlmnID", "CHR"], dtype=str, low_memory=False).dropna()
    chrom = mf.set_index("IlmnID")["CHR"]
    return m.T.astype(np.float32), chrom


# =========================================================================== statistics
def welch_t(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    a, b = X[y == 1], X[y == 0]
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    se = np.sqrt(va / len(a) + vb / len(b))
    se[se == 0] = np.inf
    return (a.mean(0) - b.mean(0)) / se


def fit_score(Xtr, ytr, Xte, k: int, C: float, Ctr=None, Cte=None):
    """In-fold z-score -> Welch |t| top-k -> L2 logistic. Optional covariate blocks are prepended unselected."""
    from sklearn.linear_model import LogisticRegression
    mu, sd = Xtr.mean(0), Xtr.std(0, ddof=1)
    sd[sd == 0] = 1.0
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    idx = np.argsort(-np.abs(welch_t(Ztr, ytr)), kind="stable")[:k]
    A, B = Ztr[:, idx], Zte[:, idx]
    if Ctr is not None:
        A, B = np.hstack([Ctr, A]), np.hstack([Cte, B])
    m = LogisticRegression(C=C, class_weight="balanced", max_iter=5000).fit(A, ytr)
    return m.decision_function(B), idx


def fit_cov(Ctr, ytr, Cte, C: float):
    from sklearn.linear_model import LogisticRegression
    m = LogisticRegression(C=C, class_weight="balanced", max_iter=5000).fit(Ctr, ytr)
    return m.decision_function(Cte)


def cv_splits(y, groups, n_folds: int, seed: int):
    from sklearn.model_selection import StratifiedGroupKFold
    return list(StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed).split(np.zeros(len(y)), y, groups))


def cv_scores(X, y, groups, repeats: list[int], p: dict, C_mat=None, mode: str = "omics") -> np.ndarray:
    """OOF decision scores, shape (len(repeats), n). mode: omics | cov | omics_cov."""
    out = np.full((len(repeats), len(y)), np.nan)
    for i, r in enumerate(repeats):
        for tr, te in cv_splits(y, groups, p["n_folds"], SEED + r):
            if mode == "cov":
                mu, sd = C_mat[tr].mean(0), C_mat[tr].std(0, ddof=1)
                sd[sd == 0] = 1
                out[i, te] = fit_cov((C_mat[tr] - mu) / sd, y[tr], (C_mat[te] - mu) / sd, p["covariate_C"])
            elif mode == "omics_cov":
                mu, sd = C_mat[tr].mean(0), C_mat[tr].std(0, ddof=1)
                sd[sd == 0] = 1
                out[i, te], _ = fit_score(X[tr], y[tr], X[te], p["k_features"], p["C"],
                                          (C_mat[tr] - mu) / sd, (C_mat[te] - mu) / sd)
            else:
                out[i, te], _ = fit_score(X[tr], y[tr], X[te], p["k_features"], p["C"])
    return out


def auroc(y, s) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s))


def auprc(y, s) -> float:
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(y, s))


def _rank_norm(S: np.ndarray) -> np.ndarray:
    """Per-repeat rank-normalised scores so that repeats are averaged on a common scale."""
    from scipy.stats import rankdata
    return np.vstack([rankdata(s) / len(s) for s in S])


def cluster_bootstrap_auc(y, s, groups, n: int, seed: int, s2=None) -> tuple[float, float, np.ndarray]:
    rng = np.random.default_rng(seed)
    ug = np.unique(groups)
    members = {g: np.where(groups == g)[0] for g in ug}
    vals, diffs = [], []
    for _ in range(n):
        pick = np.concatenate([members[g] for g in rng.choice(ug, len(ug), replace=True)])
        yy = y[pick]
        if yy.min() == yy.max():
            continue
        a = auroc(yy, s[pick])
        vals.append(a)
        if s2 is not None:
            diffs.append(a - auroc(yy, s2[pick]))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi), np.array(diffs)


def permute_labels(y, groups, rng) -> np.ndarray:
    """Label permutation; for twin pairs (one case + one control per group) swap within pairs with probability 0.5."""
    if len(np.unique(groups)) < len(y):
        yp = y.copy()
        for g in np.unique(groups):
            ix = np.where(groups == g)[0]
            if rng.random() < 0.5:
                yp[ix] = yp[ix][::-1]
        return yp
    return rng.permutation(y)


def _perm_stat(X, yp, groups, p) -> float:
    reps = list(range(p["permutation_repeats"]))
    S = cv_scores(X, yp, groups, reps, p)
    return float(np.mean([auroc(yp, s) for s in S]))


def variance_filter(X: pd.DataFrame, top: int) -> pd.DataFrame:
    v = X.var(axis=0)
    return X[v.sort_values(ascending=False, kind="stable").index[:top]]


def covariate_matrix(ds: Dataset, cols: list[str], numeric: list[str]) -> np.ndarray | None:
    if not cols:
        return None
    parts = []
    for c in cols:
        v = ds.samples[c]
        if c in numeric:
            x = pd.to_numeric(v, errors="coerce")
            parts.append(x.fillna(x.median()).to_frame(c))
        else:
            parts.append(pd.get_dummies(v.fillna("missing").astype(str), prefix=c, drop_first=True).astype(float))
    return pd.concat(parts, axis=1).to_numpy(float)


def classify(ds: Dataset, n_jobs: int = 16) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    from joblib import Parallel, delayed
    plan = load_config(CONFIG_NAME)["analysis_plan"]
    p = plan["classifier"]
    Xf = variance_filter(ds.X, p["variance_filter_top"])
    X, y, g = Xf.to_numpy(float), ds.y, ds.groups
    reps = list(range(p["n_repeats"]))
    S = cv_scores(X, y, g, reps, p)
    aucs = np.array([auroc(y, s) for s in S])
    aps = np.array([auprc(y, s) for s in S])
    savg = _rank_norm(S).mean(0)
    # NOTE: the pre-specified CI is a cluster bootstrap of the AUROC of the repeat-AVERAGED score, which is a
    # different (ensemble) estimator from auroc_mean (mean of per-repeat AUROCs). Both point values are reported
    # so that a CI can be read against the estimator it brackets (auroc_avg_score).
    lo, hi, _ = cluster_bootstrap_auc(y, savg, g, p["n_bootstrap"], SEED)
    obs_perm = float(aucs[: p["permutation_repeats"]].mean())
    rng = np.random.default_rng(SEED + 7)
    perms = [permute_labels(y, g, rng) for _ in range(p["n_permutations"])]
    null = np.array(Parallel(n_jobs=n_jobs)(delayed(_perm_stat)(X, yp, g, p) for yp in perms))
    pval = (1 + int((null >= obs_perm).sum())) / (1 + len(null))
    res = {"subset_id": ds.sid, "n_cases": int(y.sum()), "n_controls": int((1 - y).sum()), "n_features_input": ds.X.shape[1],
           "n_features_after_variance_filter": X.shape[1], "auroc_mean": float(aucs.mean()), "auroc_sd_repeats": float(aucs.std(ddof=1)),
           "auroc_avg_score": auroc(y, savg), "auroc_ci_low": lo, "auroc_ci_high": hi, "auprc_mean": float(aps.mean()), "prevalence": float(y.mean()),
           "perm_statistic_obs": obs_perm, "perm_null_mean": float(null.mean()), "perm_null_q95": float(np.quantile(null, 0.95)),
           "perm_p": pval, "n_permutations": len(null)}
    cov_cols = plan["covariates"].get(ds.sid, [])
    Cm = covariate_matrix(ds, cov_cols, plan["numeric_covariates"])
    if Cm is not None:
        Sc = cv_scores(X, y, g, reps, p, Cm, mode="cov")
        So = cv_scores(X, y, g, reps, p, Cm, mode="omics_cov")
        ca = np.array([auroc(y, s) for s in Sc])
        oa = np.array([auroc(y, s) for s in So])
        sc_avg, so_avg = _rank_norm(Sc).mean(0), _rank_norm(So).mean(0)
        clo, chi, _ = cluster_bootstrap_auc(y, sc_avg, g, p["n_bootstrap"], SEED)
        olo, ohi, diffs = cluster_bootstrap_auc(y, so_avg, g, p["n_bootstrap"], SEED, s2=sc_avg)
        res.update({"covariates": "; ".join(cov_cols), "cov_only_auroc": float(ca.mean()),
                    "cov_only_auroc_avg_score": auroc(y, sc_avg), "omics_plus_cov_auroc_avg_score": auroc(y, so_avg),
                    "cov_only_below_chance": bool(ca.mean() < 0.5), "cov_only_ci_low": clo,
                    "cov_only_ci_high": chi, "omics_plus_cov_auroc": float(oa.mean()), "omics_plus_cov_ci_low": olo,
                    "omics_plus_cov_ci_high": ohi, "delta_auroc_over_cov": float(oa.mean() - ca.mean()),
                    "delta_ci_low": float(np.percentile(diffs, 2.5)), "delta_ci_high": float(np.percentile(diffs, 97.5))})
    else:
        res.update({"covariates": UNKNOWN})
    rep = pd.DataFrame({"subset_id": ds.sid, "repeat": reps, "auroc": aucs, "auprc": aps})
    nul = pd.DataFrame({"subset_id": ds.sid, "permutation": range(len(null)), "null_statistic": null})
    return res, rep, nul


def bh(p: np.ndarray) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests
    return multipletests(p, method="fdr_bh")[1]


def differential(ds: Dataset) -> tuple[pd.DataFrame, dict]:
    from scipy import stats
    plan = load_config(CONFIG_NAME)["analysis_plan"]
    X = ds.X.to_numpy(float)
    y = ds.y
    if ds.paired:
        s = ds.samples.assign(_i=range(len(y)))
        case = s[s["label"] == 1].set_index("pair")["_i"]
        ctrl = s[s["label"] == 0].set_index("pair")["_i"].loc[case.index]
        D = X[case.to_numpy()] - X[ctrl.to_numpy()]
        n = D.shape[0]
        est = D.mean(0)
        se = D.std(0, ddof=1) / np.sqrt(n)
        dof = n - 1
        model = "paired t-test (case minus co-twin)"
    else:
        cov_cols = plan["covariates"].get(ds.sid, [])
        Cm = covariate_matrix(ds, cov_cols, plan["numeric_covariates"])
        Dm = np.column_stack([np.ones(len(y)), y] + ([Cm] if Cm is not None else []))
        beta, *_ = np.linalg.lstsq(Dm, X, rcond=None)
        resid = X - Dm @ beta
        dof = len(y) - np.linalg.matrix_rank(Dm)
        s2 = (resid ** 2).sum(0) / dof
        XtX_inv = np.linalg.pinv(Dm.T @ Dm)
        est = beta[1]
        se = np.sqrt(s2 * XtX_inv[1, 1])
        model = "OLS feature ~ label" + (" + " + " + ".join(cov_cols) if cov_cols else "")
    with np.errstate(divide="ignore", invalid="ignore"):
        t = est / se
    pv = 2 * stats.t.sf(np.abs(t), dof)
    ok = np.isfinite(pv)
    q = np.full(len(pv), np.nan)
    q[ok] = bh(pv[ok])
    tcrit = stats.t.ppf(0.975, dof)
    de = pd.DataFrame({"subset_id": ds.sid, "feature": ds.X.columns, "estimate": est, "se": se,
                       "ci_low": est - tcrit * se, "ci_high": est + tcrit * se, "t": t, "p_value": pv, "q_value_bh": q,
                       "mean_cases": X[y == 1].mean(0), "mean_controls": X[y == 0].mean(0)})
    chi = stats.chi2.isf(pv[ok], 1)
    summ = {"subset_id": ds.sid, "model": model, "n_features_tested": int(ok.sum()), "n_q_lt_0_05": int((q < 0.05).sum()),
            "n_q_lt_0_10": int((q < 0.10).sum()), "frac_p_lt_0_05": float((pv[ok] < 0.05).mean()),
            "genomic_inflation_lambda": float(np.median(chi) / stats.chi2.ppf(0.5, 1)), "dof": int(dof)}
    return de.sort_values("p_value").reset_index(drop=True), summ


def confounds(ds: Dataset) -> list[dict]:
    from scipy import stats
    rows = []
    s = ds.samples
    for col in [c for c in s.columns if c in ("c:batch", "c:test site", "c:Sex", "c:gender", "platform_id", "c:time of_blood_draw")]:
        tab = pd.crosstab(s[col].fillna("missing"), s["label"])
        if tab.shape[0] < 2:
            continue
        p = stats.fisher_exact(tab.values)[1] if tab.shape == (2, 2) else stats.chi2_contingency(tab.values)[1]
        rows.append({"subset_id": ds.sid, "check": f"label x {col}", "table": json.dumps({str(k): v for k, v in tab.to_dict("index").items()}),
                     "p_value": float(p), "test": _test_name(tab)})
    for col in ("c:age (years)", "c:age"):
        if col in s:
            a = pd.to_numeric(s[col], errors="coerce")
            p = stats.mannwhitneyu(a[s["label"] == 1].dropna(), a[s["label"] == 0].dropna()).pvalue
            rows.append({"subset_id": ds.sid, "check": f"label x {col} (Mann-Whitney)",
                         "table": json.dumps({"median_cases": float(a[s["label"] == 1].median()), "median_controls": float(a[s["label"] == 0].median())}),
                         "p_value": float(p), "test": "mann_whitney"})
    # POST HOC (added after the first results, not in the plan): does the alphabetic prefix of the sample /
    # subject identifier (a proxy for collection source) align with the label?
    pref = s["subject"].astype(str).str.extract(r"^([A-Za-z]+)", expand=False).fillna("")
    if pref.nunique() > 1:
        tab = pd.crosstab(pref, s["label"])
        p = stats.fisher_exact(tab.values)[1] if tab.shape == (2, 2) else stats.chi2_contingency(tab.values)[1]
        rows.append({"subset_id": ds.sid, "check": "POST HOC: label x subject-id alphabetic prefix",
                     "table": json.dumps({str(k): v for k, v in tab.to_dict("index").items()}), "p_value": float(p),
                     "test": _test_name(tab)})
    proxy = batch_proxy(ds)
    if proxy is not None and proxy.nunique() > 1:
        tab = pd.crosstab(proxy, s["label"])
        p = stats.chi2_contingency(tab.values)[1] if tab.shape != (2, 2) else stats.fisher_exact(tab.values)[1]
        rows.append({"subset_id": ds.sid, "check": f"POST HOC: label x {proxy.name}",
                     "table": json.dumps({str(k): v for k, v in tab.to_dict("index").items()}), "p_value": float(p),
                     "test": _test_name(tab)})
    sm = ds.sexmat
    if sm is not None and "XIST" in sm and sum(g in sm for g in SEX_Y_GENES) >= 2:
        z = (sm - sm.mean()) / sm.std(ddof=1).replace(0, 1)
        score = z[[g for g in SEX_Y_GENES if g in z]].mean(axis=1) - z["XIST"]
        male = (score > 0).map({True: "male", False: "female"})
        tab = pd.crosstab(male, s["label"])
        p = stats.fisher_exact(tab.values)[1] if tab.shape == (2, 2) else np.nan
        row = {"subset_id": ds.sid, "check": "label x expression-inferred sex", "table": json.dumps({str(k): v for k, v in tab.to_dict("index").items()}),
               "p_value": float(p), "test": "fisher_exact" if tab.shape == (2, 2) else "not_testable"}
        rec = next((c for c in ("c:Sex", "c:gender", "c:sex") if c in s and s[c].notna().any()), None)
        if rec:
            ok = s[rec].notna()
            conc = (s.loc[ok, rec].str.lower().str.startswith("m") == (male[ok] == "male")).mean()
            row["table"] = json.dumps({**json.loads(row["table"]), "concordance_with_recorded_sex": round(float(conc), 3)})
        rows.append(row)
    return rows


def _test_name(tab: pd.DataFrame) -> str:
    return "fisher_exact" if tab.shape == (2, 2) else f"chi_square_{tab.shape[0]}x{tab.shape[1]}"


def batch_proxy(ds: Dataset) -> pd.Series | None:
    """POST HOC processing-batch proxies read from deposited metadata: 450K chip (Sentrix) id from the GSM title,
    NanoString cartridge id from each RCC file."""
    s = ds.samples
    if ds.kind == "methylation":
        chip = s["title"].str.extract(r"(\d{10,12})_R\d\dC\d\d", expand=False)
        return chip.rename("450K chip (Sentrix) id")
    if ds.kind == "nanostring":
        tar_path = raw_dir(SOURCE_ID) / "GSE275334/GSE275334_RAW.tar"
        cart = {}
        with tarfile.open(tar_path) as tf:
            for mem in tf.getmembers():
                gsm = mem.name.split("_")[0]
                if gsm in s.index:
                    txt = gzip.decompress(tf.extractfile(mem).read()).decode("utf-8", "replace")
                    m = re.search(r"CartridgeID,(\S+)", txt)
                    cart[gsm] = m.group(1) if m else "unknown"
        return pd.Series(cart).reindex(s.index).rename("NanoString cartridge id")
    return None


# =========================================================================== replication
def transfer(src: Dataset, tgt: Dataset, de_src: pd.DataFrame, de_tgt: pd.DataFrame, tier: str, n_jobs: int = 16,
             restrict: set | None = None) -> tuple[dict, dict]:
    from joblib import Parallel, delayed
    from scipy import stats
    plan = load_config(CONFIG_NAME)["analysis_plan"]
    p, tp = plan["classifier"], plan["transfer"]
    uni = set(src.X.columns) & set(tgt.X.columns)
    if restrict is not None:
        uni &= set(restrict)
    uni = sorted(uni)
    base = {"source": src.sid, "target": tgt.sid, "tier": tier, "n_shared_genes": len(uni)}
    if len(uni) < p["k_features"]:
        return {**base, "note": "too few shared genes"}, {**base}
    Xs_raw = src.X[uni]
    keep = variance_filter(Xs_raw, p["variance_filter_top"]).columns
    zs = lambda X: ((X - X.mean()) / X.std(ddof=1).replace(0, 1)).fillna(0)
    Xs, Xt = zs(src.X[uni])[keep].to_numpy(float), zs(tgt.X[uni])[keep].to_numpy(float)
    ys, yt, gt = src.y, tgt.y, tgt.groups
    score, idx = fit_score(Xs, ys, Xt, p["k_features"], p["C"])
    a = auroc(yt, score)
    lo, hi, _ = cluster_bootstrap_auc(yt, score, gt, tp["n_bootstrap"], SEED)
    rng = np.random.default_rng(SEED + 11)
    null = np.array([auroc(permute_labels(yt, gt, rng), score) for _ in range(tp["n_permutations"])])
    pv = (1 + int((null >= a).sum())) / (1 + len(null))
    res = {**base, "n_features_used": int(len(keep)), "target_n_cases": int(yt.sum()), "target_n_controls": int((1 - yt).sum()),
           "transfer_auroc": a, "ci_low": lo, "ci_high": hi, "perm_p": pv, "perm_null_q95": float(np.quantile(null, 0.95)),
           "replicated": bool(lo > 0.5 and pv < 0.05), "selected_genes": ";".join(np.array(keep)[idx][:50])}
    # signature overlap
    o = plan["overlap"]
    ds_, dt_ = de_src.set_index("feature").loc[uni], de_tgt.set_index("feature").loc[uni]
    top_s = set(ds_["p_value"].nsmallest(o["top_n"]).index)
    top_t = set(dt_["p_value"].nsmallest(o["top_n"]).index)
    ov = top_s & top_t
    rng = np.random.default_rng(SEED + 13)
    uarr = np.array(uni)
    null_ov = np.array([len(set(rng.choice(uarr, o["top_n"], replace=False)) & top_t) for _ in range(o["n_null_draws"])])
    conc = int(sum(np.sign(ds_.loc[g, "estimate"]) == np.sign(dt_.loc[g, "estimate"]) for g in ov))
    rho = stats.spearmanr(ds_["t"], dt_["t"], nan_policy="omit").statistic
    ovr = {**base, "top_n": o["top_n"], "overlap": len(ov), "expected_overlap_random": float(null_ov.mean()),
           "overlap_perm_p": (1 + int((null_ov >= len(ov)).sum())) / (1 + len(null_ov)),
           "hypergeom_p": float(stats.hypergeom.sf(len(ov) - 1, len(uni), o["top_n"], o["top_n"])),
           "direction_concordant": conc, "direction_concordance_binom_p": float(stats.binomtest(conc, len(ov), 0.5, alternative="greater").pvalue) if ov else np.nan,
           "spearman_de_t_descriptive": float(rho), "overlap_genes": ";".join(sorted(ov)),
           "n_fdr05_source": int((ds_["q_value_bh"] < 0.05).sum()), "n_fdr05_target": int((dt_["q_value_bh"] < 0.05).sum())}
    return res, ovr


# =========================================================================== run
def holm(p: list[float]) -> list[float]:
    from statsmodels.stats.multitest import multipletests
    return list(multipletests(p, method="holm")[1])


def run(n_jobs: int = 16) -> dict:
    scr, sel, _ = run_screen()
    plan = load_config(CONFIG_NAME)["analysis_plan"]
    sids = sel["subset_id"].tolist()
    data, res, reps, nulls, desum, detop, conf, de_all = {}, [], [], [], [], [], [], {}
    for sid in sids:
        ds = load_dataset(sid)
        data[sid] = ds
        r, rp, nl = classify(ds, n_jobs)
        row = sel.set_index("subset_id").loc[sid]
        r.update({"condition_id": row["condition_id"], "role": row["role"], "compartment": row["compartment"], "assay": row["assay"]})
        res.append(r), reps.append(rp), nulls.append(nl)
        de, summ = differential(ds)
        de_all[sid] = de
        desum.append(summ)
        detop.append(de.head(25))
        conf.extend(confounds(ds))
        print(sid, {k: r[k] for k in ("auroc_mean", "auroc_ci_low", "auroc_ci_high", "perm_p")}, summ["n_q_lt_0_05"], flush=True)
    cls = pd.DataFrame(res)
    cls["perm_p_holm"] = holm(cls["perm_p"].tolist())
    cls["within_dataset_signal"] = cls["perm_p"] < plan["classifier"]["alpha"]
    tr, ov = [], []
    for pr in plan["transfer"]["pairs"]:
        a, b = transfer(data[pr["source"]], data[pr["target"]], de_all[pr["source"]], de_all[pr["target"]], pr["tier"], n_jobs)
        tr.append(a), ov.append(b)
        print("transfer", pr, {k: a.get(k) for k in ("transfer_auroc", "ci_low", "ci_high", "perm_p")}, flush=True)
    # POST HOC (added after the pre-specified results were seen; cannot change any call): the two whole-blood
    # Long COVID cohorts transferred on the genes of the NanoString Immune Exhaustion panel only, to ask whether the
    # panel-restricted agreement seen in the secondary transfers also holds between the whole-blood cohorts.
    panel = set(data["GSE275334_lc"].X.columns) if "GSE275334_lc" in data else set()
    ph = []
    for a_, b_ in (("GSE226260_wb", "GSE270045"), ("GSE270045", "GSE226260_wb")):
        if panel and a_ in data and b_ in data:
            x, _ = transfer(data[a_], data[b_], de_all[a_], de_all[b_], "post_hoc_panel_genes", n_jobs, restrict=panel)
            ph.append(x)
    pd.DataFrame(ph).to_csv(TABLES / "geo_cohort_transfer_posthoc.csv", index=False)
    tr = pd.DataFrame(tr)
    prim = tr["tier"] == "primary"
    tr.loc[prim, "perm_p_holm_primary"] = holm(tr.loc[prim, "perm_p"].tolist())
    out = {"classifier": cls, "classifier_repeats": pd.concat(reps), "permutation_null": pd.concat(nulls),
           "de_summary": pd.DataFrame(desum), "de_top": pd.concat(detop), "confounds": pd.DataFrame(conf),
           "transfer": tr, "overlap": pd.DataFrame(ov)}
    for k, v in out.items():
        v.to_csv(TABLES / f"geo_cohort_{k}.csv", index=False)
    write_signatures(cls, tr, de_all, data, sel, out["confounds"])
    write_registry_and_audit(scr, sel, data)
    return out


def write_signatures(cls: pd.DataFrame, tr: pd.DataFrame, de_all: dict, data: dict, sel: pd.DataFrame,
                     conf: pd.DataFrame | None = None) -> None:
    from ..provenance import add_provenance
    from ..store import write_table
    plan = load_config(CONFIG_NAME)["analysis_plan"]
    rows = []
    for r in cls.itertuples():
        if not r.within_dataset_signal:
            continue
        rep = tr[(tr["tier"] == "primary") & ((tr["source"] == r.subset_id) | (tr["target"] == r.subset_id))]
        level = "replicated_transfer" if (len(rep) and rep["replicated"].any()) else "dataset_specific"
        de = de_all[r.subset_id]
        de = de[de["q_value_bh"] < 0.05].head(50)
        flags = []
        if conf is not None and len(conf):
            c = conf[(conf["subset_id"] == r.subset_id) & (conf["p_value"] < 0.05)]
            flags = c["check"].tolist()
        caveat = ("dataset-specific; public cohort, batch/site structure largely unrecorded; not participant-linked to any "
                  "other dataset" + (f"; CONFOUND FLAGGED: {'; '.join(flags)}" if flags else ""))
        gse = r.subset_id.split("_")[0]
        for d in de.itertuples():
            rows.append({"object_id": f"signature:geo_{r.subset_id}|{r.condition_id}_case_control|{d.feature}",
                         "dataset_id": f"geo_{r.subset_id}", "phenotype_id": f"{r.condition_id}_case_control",
                         "phenotype_label": f"{r.condition_id} case vs control ({gse})", "condition_id": r.condition_id,
                         "is_proxy": False, "label_basis": "study case definition recorded in GEO sample metadata",
                         "feature": d.feature, "feature_label": d.feature, "effect_measure": "difference_in_log2_units" if data[r.subset_id].kind != "methylation" else "difference_in_M_value",
                         "effect_size": float(d.estimate), "ci_low": float(d.ci_low), "ci_high": float(d.ci_high),
                         "null_value": 0.0, "se": float(d.se), "p_value": float(d.p_value), "q_value_bh": float(d.q_value_bh),
                         "fdr_significant": True,
                         "direction": "higher_in_cases" if d.estimate > 0 else "lower_in_cases", "n_cases": int(r.n_cases),
                         "n_controls": int(r.n_controls), "classifier_auroc": float(r.auroc_mean), "classifier_perm_p": float(r.perm_p),
                         "replication_status": level, "compartment": r.compartment, "assay": r.assay,
                         "confound_flagged": bool(flags), "caveats": caveat})
    df = pd.DataFrame(rows, columns=["object_id", "dataset_id", "phenotype_id", "phenotype_label", "condition_id", "is_proxy",
                                     "label_basis", "feature", "feature_label", "effect_measure", "effect_size", "ci_low",
                                     "ci_high", "null_value", "se", "p_value", "q_value_bh", "fdr_significant", "direction", "n_cases", "n_controls", "classifier_auroc", "classifier_perm_p",
                                     "replication_status", "compartment", "assay", "confound_flagged", "caveats"])
    df = add_provenance(df, data_layer="person", source_name="NCBI GEO public case/control series (processed matrices)",
                        source_version=f"GEO FTP/acc.cgi retrieved {utc_now_iso()[:10]}", retrieved_at=utc_now_iso(),
                        evidence_type="differential_expression", source_record_id="object_id",
                        evidence_level=(df["replication_status"] + np.where(df["confound_flagged"], "_confound_flagged", "")
                                        if len(df) else ""),
                        provenance_notes="Cohort-level contrast computed by measure_it.omics.geo_cohorts from person-level GEO profiles; "
                                         "written only for subsets whose classifier passed the pre-specified permutation test")
    write_table(df, SIGNATURE_TABLE, producer=MODULE,
                description="GEO cohort DE features for subsets with a permutation-significant classifier (dataset-specific unless replicated)")


def _md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        out.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(out)


def _retrieval_window() -> str:
    p = RAW / SOURCE_ID / "MANIFEST.json"
    times = sorted(v.get("retrieved_at") for v in json.loads(p.read_text())["files"].values() if v.get("retrieved_at"))
    return f"{times[0]} to {times[-1]}" if times else UNKNOWN


def write_registry_and_audit(scr: pd.DataFrame, sel: pd.DataFrame, data: dict) -> None:
    """DATA_AUDIT.md (generated from the measured tables) + registry_entry.yaml via write_registry_entry."""
    from ..registry import write_registry_entry
    window = _retrieval_window()
    size = {sid: {"cases": int(ds.y.sum()), "controls": int((1 - ds.y).sum()), "features_analysed": int(ds.X.shape[1]),
                  "missing_values": int(np.isnan(ds.X.to_numpy(float)).sum())} for sid, ds in data.items()}
    reasons = scr["fail_reasons"].fillna("").str.split("; ").explode()
    reasons = reasons[reasons != ""].value_counts().rename_axis("fail_reason").reset_index(name="rows")
    selt = sel[["condition_id", "subset_id", "role", "n_case_subjects", "n_control_subjects", "compartment", "assay",
                "compat_class", "matrix_note"]].copy()
    miss = pd.DataFrame([{"subset_id": k, **v} for k, v in size.items()])
    L = ["# DATA AUDIT — NCBI GEO public case/control omics cohorts (Long COVID, ME/CFS)\n",
         "Generated by `measure_it.omics.geo_cohorts` from the files it downloaded; numbers are measured.\n",
         "| Field | Value |", "|---|---|",
         f"| source_id | {SOURCE_ID} |",
         f"| Source (dataset/API name, exact files/endpoints) | GEO sample metadata {GEO_ACC}?acc=<GSE>&targ=gsm&form=text&view=brief; "
         f"FTP series directories {GEO_FTP}/<GSEnnn>/<GSE>/(suppl,matrix)/; NCBI-generated counts {GEO_DOWNLOAD}?type=rnaseq_counts; "
         f"NCBI Gene {NCBI_GENE_INFO_URL}; GPL570 annot; Illumina 450K v1.1 manifest (GPL13534 suppl) |",
         "| Publishing organization | NCBI / National Library of Medicine (GEO); data submitted by the study authors of each series |",
         f"| Retrieval date (UTC) | {window} |",
         "| Source version / release | GEO series as deposited on the retrieval date (see MANIFEST.json sha256 per file) |",
         "| Source update date / cadence | per series (submitter updates) |",
         "| License / access conditions | open GEO records, no login; submitter terms apply; no controlled-access data touched |",
         "| Unit of observation | one sample = one person of one GEO cohort (one sample per subject kept) |",
         f"| Sample size (actual, as ingested) | {sum(v['cases'] for v in size.values())} cases + {sum(v['controls'] for v in size.values())} controls over {len(size)} analysed subsets (table below) |",
         "| Geography (resolution, vintage) | none |",
         "| Person-level? | yes: person-level omics profiles within each cohort |", "| Geographic? | no |",
         "| Omics? | yes (whole-blood / PBMC / leukocyte transcriptome, plasma cfRNA, T-cell DNA methylation) |",
         "| Wearable? | no |",
         "| True participant linkage across modalities? | no: each cohort is separate; no identifier links any GEO sample to any wearable or other dataset |\n",
         "## Screen\n",
         f"{len(scr)} screened rows ({scr['gse'].nunique()} series) from `geo_study_catalog` (long_covid, me_cfs); "
         f"{int((scr['decision'] == 'pass').sum())} passing rows. Every row with its reason: `results/tables/{SCREEN_CSV}`. Failure reasons:\n",
         _md_table(reasons),
         "\n## Selected subsets (pre-specified rule; `results/tables/geo_cohort_selection.csv`)\n", _md_table(selt),
         "\n## Files / endpoints retrieved\n",
         f"All files and GEO metadata texts are under `data/raw/{SOURCE_ID}/` with URL, bytes, sha256 and retrieval time in "
         "MANIFEST.json. Analysable sample lists (label, subject, pair, covariates): `samples_<subset>.csv`.\n",
         "## Key variables\n",
         "`label` (1 = case per the study's own definition in GEO metadata, 0 = control), `subject`, `pair` (twin pair for "
         "GSE16059), covariates as recorded (sex, age, batch, test site where present). Features: log2 CPM (RNA-seq), "
         "log2 normalised counts (NanoString), log2 intensity (Affymetrix), M-values (450K).\n",
         "## Missingness (measured, analysed matrices)\n", _md_table(miss),
         "\nCovariates recorded per subset: GSE293840 sex/batch/test site/age bin (age bin mixes bins 1-4 and years across "
         "batches; not used), GSE275334 sex/age, GSE227375 sex, GSE16059 sex (matched in pairs); none for GSE226260_wb, GSE270045, GSE156792.\n",
         "## Access problems recorded\n",
         "* GEO's `Human.GRCh38.p13.annot.tsv.gz` link (`/geo/download/?format=file&type=rnaseq_counts&file=...`) served a "
         "reCAPTCHA challenge page (HTTP 200, text/html) to an automated client on 2026-09-24; NCBI Gene `gene_info` used instead. "
         "The GSE227375 NCBI-generated count file downloaded normally.\n"
         "* GSE156792 M/U intensity files contain no probe identifiers (485,512 unlabelled rows); the series-matrix beta values were used.\n"
         "* GSE339049 (monocyte ATAC) raw files withheld by the submitters; its processed peak matrix was built with the labels (not used).\n",
         "## Linkage strategy\n",
         "Within a cohort, matrix columns are matched to GSM samples by title / sample id / subject number (see loader). "
         "Nothing is joined across cohorts at the person level; cross-cohort analysis is model transfer on shared gene symbols.\n",
         "## Limitations and caveats\n",
         "Small cohorts; case/control processing batches mostly unrecorded; test site partly confounded with label in "
         "GSE293840; GSE270045 'Long COVID' cases meet ME/CFS criteria; timing post-infection varies (GSE226260_wb 3 to >12 "
         "months); medication unrecorded; GEO-deposited processed matrices differ in pipeline across cohorts.\n",
         "## Processed outputs\n",
         f"`{SIGNATURE_TABLE}` (rows only for subsets with a permutation-significant classifier) and "
         "`results/tables/geo_cohort_*.csv`.\n", "## Reproduce\n", f"`uv run python -m {MODULE}`\n"]
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text("\n".join(L))
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "NCBI GEO public case/control omics cohorts (Long COVID, ME/CFS)",
        "publisher": "NCBI / National Library of Medicine (GEO); per-series submitters",
        "landing_url": "https://www.ncbi.nlm.nih.gov/geo/",
        "access_urls": [GEO_ACC, GEO_FTP, GEO_DOWNLOAD, NCBI_GENE_INFO_URL],
        "license": "open GEO records; submitter terms apply to individual series",
        "access_conditions": "open, no login; GEO's /geo/download annotation link served a reCAPTCHA page to an automated client (worked around with NCBI Gene FTP)",
        "retrieved_at": window.split(" to ")[-1],
        "source_version": "GEO series as deposited on the retrieval date; sha256 per file in MANIFEST.json",
        "update_date": "per series",
        "data_layer": "person",
        "unit_of_observation": "person (one sample per subject) within one GEO cohort",
        "sample_size": {"screened_rows": int(len(scr)), "screened_series": int(scr["gse"].nunique()),
                        "passing_rows": int((scr["decision"] == "pass").sum()), "analysed_subsets": size},
        "geographic_resolution": "none",
        "person_level": True, "geographic": False, "omics": True, "wearable": False,
        "participant_linkage": "within-cohort only; no cross-dataset participant key",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": [SIGNATURE_TABLE, "results/tables/geo_cohort_*.csv"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": MODULE,
        "limitations": ["small cohorts (33-168 subjects)", "case/control batch structure mostly unrecorded",
                        "replication testable only for two blood-transcriptome pairs",
                        "person-level within each cohort; never linked to wearable participants"],
    })


def main() -> None:  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--screen-only", action="store_true")
    ap.add_argument("--download", action="store_true", help="download matrices of the selected subsets and stop")
    ap.add_argument("--jobs", type=int, default=16)
    args = ap.parse_args()
    if args.screen_only or args.download:
        scr, sel, _ = run_screen()
        print(scr["decision"].value_counts().to_string())
        print(sel[["condition_id", "subset_id", "role", "n_case_subjects", "n_control_subjects", "compat_class"]].to_string())
        if args.download:
            fetch_gene_annotation()
            for sid in sel["subset_id"]:
                print(sid, [str(p) for p in fetch_matrix_files(rule_by_id(sid))])
        return
    run(args.jobs)


if __name__ == "__main__":  # pragma: no cover
    main()
