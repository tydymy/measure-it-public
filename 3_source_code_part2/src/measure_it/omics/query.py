"""Condition-level molecular enrichment/evidence: shared schema, helpers and the query tool.

This module is the single interface to the ``condition_molecular`` layer.

* Shared by the three ingestion modules (``open_targets``, ``gwas_catalog``,
  ``geo_sra``): the common row schema (``EVIDENCE_COLUMNS``), the condition ->
  ontology-identifier resolution from the already-built ontology tables
  (``condition_query_ids``), raw-response archiving with MANIFEST entries
  (``save_raw_json``) and the partition union (``build_union``).
* Query tool: ``get_molecular_context(condition)`` -> dict grouped by
  evidence_type with counts per source and top entities, plus provenance, or
  ``UNKNOWN / NOT AVAILABLE`` when nothing is recorded.

Guardrails (SPEC Phase 2, CONVENTIONS section 3)
------------------------------------------------
Everything here is **condition-level molecular enrichment/evidence** gathered
from public databases about a disease concept. None of it is measured on, or
linked to, any participant in the person-level layer, and it must never be
described as "patient multi-omics". No combined "omics score" is computed:
each row keeps only the score its source publishes (``source_score`` +
``source_score_label``) or none at all.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import RAW, UNKNOWN, raw_dir, utc_now_iso
from ..download import _record as _manifest_record  # shared, locked MANIFEST writer
from ..provenance import PROVENANCE_COLUMNS
from ..store import read_table, table_exists, union_partitions

EVIDENCE_TABLE = "condition_molecular_evidence"
COVERAGE_TABLE = "condition_molecular_coverage"
LAYER_LABEL = "condition-level molecular enrichment/evidence (not participant-linked; not patient multi-omics)"

# Common columns every partition written by the omics modules carries (plus provenance).
EVIDENCE_COLUMNS = [
    "condition_id",            # registry canonical_condition_id
    "ontology_id",             # CURIE of the registry identifier that was queried (MONDO:0005404)
    "ontology_id_role",        # role of that id for the condition in condition_ontology_mappings
    "ontology_match",          # exact / narrow / broad / related (registry predicate for the id)
    "source_disease_id",       # the id as the source spells it (MONDO_0005404), or '' for text queries
    "source_disease_label",    # the source's own label for that id
    "evidence_type",           # provenance EVIDENCE_TYPES vocabulary
    "source_evidence_category",  # the source's own category (OT datatype id, GEO series type, ...)
    "entity_type",             # gene / variant / pathway / drug / study / analyte
    "entity_id",
    "entity_label",
    "source_database",
    "source_accession",
    "evidence_direction",      # only when the source defines it
    "study_population",
    "sample_size",
    "sample_size_basis",
    "source_score",            # source-provided only; never computed here
    "source_score_label",
    "date_retrieved",
]

COVERAGE_COLUMNS = [
    "condition_id", "ontology_id", "ontology_id_role", "ontology_match", "source_database",
    "id_status_in_source", "source_disease_label", "coverage_status", "n_records", "detail",
]

DEMO_CLUSTER = ("long_covid", "me_cfs", "pots")

# Ordering of the source's own clinical-stage labels (used only to sort drugs for display).
CLINICAL_STAGE_ORDER = {"APPROVAL": 10, "PHASE_4": 9, "PHASE_3": 8, "PHASE_2_3": 7, "PHASE_2": 6, "PHASE_1_2": 5,
                        "PHASE_1": 4, "EARLY_PHASE_1": 3, "IND": 2, "PRECLINICAL": 1, "UNKNOWN": 0}


# --------------------------------------------------------------------------- condition ids
def curie_to_short(curie: str) -> str:
    """'MONDO:0005404' -> 'MONDO_0005404' (the form Open Targets and the GWAS Catalog use)."""
    return str(curie).replace(":", "_", 1)


def short_to_curie(short: str) -> str:
    """'MONDO_0005404' -> 'MONDO:0005404'."""
    s = str(short)
    return s.replace("_", ":", 1) if "_" in s and ":" not in s else s


def _as_list(v) -> list:
    if v is None:
        return []
    if isinstance(v, float) and np.isnan(v):
        return []
    if isinstance(v, (list, tuple, np.ndarray)):
        return [x for x in v if x is not None and str(x) != ""]
    return [v]


@lru_cache(maxsize=1)
def _registry() -> pd.DataFrame:
    return read_table("condition_registry")


@lru_cache(maxsize=1)
def _mappings() -> pd.DataFrame:
    return read_table("condition_ontology_mappings")


def condition_query_ids(include_obsolete_efo: bool = False) -> pd.DataFrame:
    """Every registry MONDO/EFO identifier per condition, with its role and match predicate.

    Identifiers come only from ``condition_registry.mondo_ids`` and ``efo_ids`` (built and
    verified by the ontology module); ``related_mondo_ids`` (context such as a broader parent)
    are not the condition and are not queried. With ``include_obsolete_efo`` the EFO ids the
    ontology module recorded as obsolete (replaced by a Mondo class) are added with
    ``id_status='obsolete_in_OLS4'`` so a source still using them is not silently missed.
    Conditions with no identifier get one row with ``ontology_id=''``.
    """
    reg = _registry()
    maps = _mappings()
    mm = maps[maps["target_ontology"].isin(["MONDO", "EFO"])]
    rows = []
    for _, r in reg.iterrows():
        cid = r["canonical_condition_id"]
        ids = []
        for i in _as_list(r["mondo_ids"]) + _as_list(r["efo_ids"]):
            if i not in ids:
                ids.append(i)
        sub = mm[mm["canonical_condition_id"] == cid]
        for i in ids:
            m = sub[(sub["target_id"] == i) & (sub["target_ontology"] == "MONDO")]
            if m.empty:
                m = sub[sub["target_id"] == i]
            m0 = m.iloc[0] if len(m) else None
            rows.append({
                "condition_id": cid,
                "condition_label": r["preferred_name"],
                "ontology_id": i,
                "source_query_id": curie_to_short(i),
                "ontology_id_role": (m0["mondo_role"] if m0 is not None else "registry_list"),
                "ontology_match": (m0["predicate_condition"] if m0 is not None else ""),
                "registry_label": (m0["target_label"] if m0 is not None else ""),
                "id_status": "current",
                "is_primary": i == r.get("primary_mondo_id"),
                "grouping_only": bool(r.get("grouping_only", False)),
            })
        if include_obsolete_efo:
            obs = sub[(sub["target_ontology"] == "EFO") & sub["target_status"].astype(str).str.startswith("obsolete")]
            for _, o in obs.iterrows():
                if o["target_id"] in ids:
                    continue
                ids.append(o["target_id"])
                rows.append({
                    "condition_id": cid, "condition_label": r["preferred_name"],
                    "ontology_id": o["target_id"], "source_query_id": curie_to_short(o["target_id"]),
                    "ontology_id_role": f"obsolete EFO id ({o['mondo_role']} class {o['mondo_id']})",
                    "ontology_match": o["predicate_condition"], "registry_label": o["target_label"],
                    "id_status": "obsolete_in_OLS4", "is_primary": False,
                    "grouping_only": bool(r.get("grouping_only", False)),
                })
        if not ids:
            rows.append({"condition_id": cid, "condition_label": r["preferred_name"], "ontology_id": "",
                         "source_query_id": "", "ontology_id_role": "", "ontology_match": "",
                         "registry_label": "", "id_status": "no identifier", "is_primary": False,
                         "grouping_only": bool(r.get("grouping_only", False))})
    return pd.DataFrame(rows)


def condition_search_terms() -> dict[str, list[str]]:
    reg = _registry()
    return {r["canonical_condition_id"]: _as_list(r["search_terms"]) for _, r in reg.iterrows()}


# --------------------------------------------------------------------------- raw archiving
def save_raw_json(source_id: str, relpath: str, payload, *, url: str, request: dict | None = None,
                  retrieved_at: str | None = None) -> Path:
    """Write an API response under data/raw/<source_id>/<relpath> and record it in MANIFEST.json."""
    base = raw_dir(source_id)
    dest = base / relpath
    dest.parent.mkdir(parents=True, exist_ok=True)
    blob = json.dumps(payload, indent=1, sort_keys=True, default=str).encode()
    dest.write_bytes(blob)
    _manifest_record(source_id, str(dest.relative_to(RAW / source_id)), {
        "url": url,
        "request": request or {},
        "bytes": len(blob),
        "sha256": hashlib.sha256(blob).hexdigest(),
        "retrieved_at": retrieved_at or utc_now_iso(),
        "content_type": "application/json",
    })
    return dest


def retrieval_window(source_id: str) -> str:
    """First and last response time recorded in data/raw/<source_id>/MANIFEST.json."""
    p = RAW / source_id / "MANIFEST.json"
    if not p.exists():
        return UNKNOWN
    times = sorted(v.get("retrieved_at") for v in json.loads(p.read_text()).get("files", {}).values()
                   if v.get("retrieved_at"))
    return f"{times[0]} to {times[-1]} (API responses, from MANIFEST.json)" if times else UNKNOWN


def safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(s))[:120]


# --------------------------------------------------------------------------- frame helpers
def finalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Put EVIDENCE_COLUMNS first (creating missing ones as null) and provenance last."""
    out = df.copy()
    for c in EVIDENCE_COLUMNS:
        if c not in out.columns:
            out[c] = None
    out["sample_size"] = pd.to_numeric(out["sample_size"], errors="coerce").astype("float64")
    out["source_score"] = pd.to_numeric(out["source_score"], errors="coerce").astype("float64")
    extra = [c for c in out.columns if c not in EVIDENCE_COLUMNS and c not in PROVENANCE_COLUMNS]
    prov = [c for c in PROVENANCE_COLUMNS if c in out.columns and c != "evidence_type"]
    return out[EVIDENCE_COLUMNS + extra + prov]


def stringify_object_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Parquet-safe: lists/dicts in object columns become JSON strings; mixed scalars become str."""
    out = df.copy()
    for c in out.columns:
        if out[c].dtype == object:
            vals = out[c]
            if vals.map(lambda v: isinstance(v, (list, dict, tuple, np.ndarray))).any():
                out[c] = vals.map(lambda v: json.dumps(list(v) if isinstance(v, (tuple, np.ndarray)) else v,
                                                       default=str)
                                  if isinstance(v, (list, dict, tuple, np.ndarray)) else v)
            kinds = {type(v) for v in out[c].dropna()}
            if len(kinds) > 1:
                out[c] = out[c].map(lambda v: None if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
    return out


def build_union() -> dict:
    """Union every condition_molecular_evidence__* (and coverage) partition with store.union_partitions."""
    out = {}
    for table in (EVIDENCE_TABLE, COVERAGE_TABLE):
        try:
            p = union_partitions(table, producer="measure_it.omics.query.build_union")
            out[table] = len(pd.read_parquet(p, columns=["data_layer"]))
        except FileNotFoundError:
            out[table] = 0
    return out


def exit_cleanly(code: int = 0) -> None:
    """Flush and leave without interpreter teardown.

    Measured on this aarch64 host: roughly 1 in 6 short processes that import the store stack
    abort *after* all work is done with "terminate called without an active exception"
    (exit 134) during native-library teardown. Every table is already written atomically
    (os.replace) by then, so skipping teardown only removes a spurious failure code.
    """
    import os
    import sys
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(code)


# --------------------------------------------------------------------------- cross-source agreement
def _split_symbols(series: pd.Series) -> set[str]:
    out: set[str] = set()
    for v in series.dropna():
        out.update(x.strip().upper() for x in str(v).split(";") if x.strip())
    return out


def genetic_gene_agreement(ev: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per condition: genes with Open Targets genetic-association evidence vs GWAS Catalog mapped genes.

    Counts only (SPEC validation 3: report agreement/disagreement across databases); nothing is
    combined into a score. Open Targets genes = targets whose `genetic_association` datatype score
    is > 0 on any registry id of the condition; GWAS genes = the Catalog's mapped genes of every
    curated association of any registry id. Compared by upper-case gene symbol.
    """
    if ev is None:
        ev = read_table(EVIDENCE_TABLE) if table_exists(EVIDENCE_TABLE) else pd.DataFrame()
    rows = []
    for cid in _registry()["canonical_condition_id"]:
        d = ev[ev["condition_id"] == cid] if len(ev) else ev
        ot = d[(d.get("source_evidence_category") == "overall_association")] if len(d) else d
        col = "ot_datatype_score__genetic_association"
        ot_genes = _split_symbols(ot.loc[ot[col].fillna(0) > 0, "entity_label"]) if len(ot) and col in ot else set()
        gw = d[(d.get("source_database") == "GWAS Catalog") & (d.get("entity_type") == "variant")] if len(d) else d
        gw_genes = _split_symbols(gw["mapped_genes"]) if len(gw) and "mapped_genes" in gw else set()
        both = ot_genes & gw_genes
        rows.append({"condition_id": cid, "open_targets_genetic_genes": len(ot_genes),
                     "gwas_catalog_mapped_genes": len(gw_genes), "in_both": len(both),
                     "open_targets_only": len(ot_genes - gw_genes), "gwas_catalog_only": len(gw_genes - ot_genes),
                     "shared_examples": ";".join(sorted(both)[:15])})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- query tool
def _resolve_condition(condition: str) -> tuple[str | None, dict]:
    reg = _registry()
    c = str(condition or "").strip()
    if c in set(reg["canonical_condition_id"]):
        return c, {"status": "matched", "input_type": "canonical_condition_id"}
    try:
        from ..ontology.normalize import normalize_condition
        res = normalize_condition(c)
    except Exception as exc:  # ontology tables missing
        return None, {"status": UNKNOWN, "reason": f"normalization unavailable: {exc}"}
    if res.get("status") in ("matched", "partial", "ambiguous") and res.get("matches"):
        return res["matches"][0]["canonical_condition_id"], {
            "status": res["status"], "input_type": res.get("input_type"),
            "match_reason": res["matches"][0].get("match_reason", ""),
            "other_matches": [m["canonical_condition_id"] for m in res["matches"][1:]]}
    return None, {"status": UNKNOWN, "reason": res.get("reason", "no registry condition matches")}


def _entity_counts(df: pd.DataFrame, top_n: int) -> list[dict]:
    """Top entities of ONE entity type, ranked by the source-provided score when the source gives one;
    otherwise (study metadata) by the condition-mention QC flag, then record count, then sample size.
    The ranking basis is reported with every entity."""
    if df.empty:
        return []
    d = df.copy()
    d["_mention"] = d["mentions_condition"].fillna(False).astype(bool) if "mentions_condition" in d else False
    g = d.groupby("entity_id", dropna=False)
    agg = g.agg(entity_type=("entity_type", "first"), entity_label=("entity_label", "first"),
                n_records=("entity_id", "size"), max_source_score=("source_score", "max"),
                max_sample_size=("sample_size", "max"), mentions_condition=("_mention", "max"),
                sources=("source_database", lambda s: sorted(set(map(str, s)))),
                accessions=("source_accession", lambda s: sorted(set(map(str, s.dropna())))[:5])).reset_index()
    if agg["max_source_score"].notna().any():
        agg = agg.sort_values(["max_source_score", "n_records"], ascending=[False, False])
        basis = "max source-provided score (" + "; ".join(sorted(set(df["source_score_label"].dropna().astype(str)))) + ")"
    elif "max_clinical_stage_for_indication" in d and d["max_clinical_stage_for_indication"].notna().any():
        stage = d.groupby("entity_id")["max_clinical_stage_for_indication"].agg(
            lambda s: max((CLINICAL_STAGE_ORDER.get(str(x), 0) for x in s), default=0))
        agg["_stage"] = agg["entity_id"].map(stage).fillna(0)
        agg["max_clinical_stage"] = agg["entity_id"].map(
            d.groupby("entity_id")["max_clinical_stage_for_indication"].first())
        agg = agg.sort_values(["_stage", "n_records", "entity_label"], ascending=[False, False, True])
        basis = "source-reported maximum clinical stage for the indication (Open Targets), then records"
    elif "adj_p_value" in d and d["adj_p_value"].notna().any():
        agg["_p"] = agg["entity_id"].map(d.groupby("entity_id")["adj_p_value"].min())
        agg["min_published_adj_p"] = agg["_p"]
        agg = agg.sort_values(["_p", "n_records"], ascending=[True, False], na_position="last")
        basis = "smallest published adjusted p-value (source-reported statistic, not a score)"
    elif "p_value" in d and d["p_value"].notna().any():
        # e.g. GWAS Catalog variants: rank by the smallest curated p-value (a statistic, not a score)
        agg["_p"] = agg["entity_id"].map(d.groupby("entity_id")["p_value"].min())
        agg["min_reported_p"] = agg["_p"]
        agg = agg.sort_values(["_p", "n_records", "entity_id"], ascending=[True, False, True], na_position="last")
        basis = "smallest reported p-value (source-reported statistic, not a score)"
    else:
        has_mention = "mentions_condition" in df and df["mentions_condition"].notna().any()
        agg = agg.sort_values(["mentions_condition", "n_records", "max_sample_size", "entity_id"],
                              ascending=[False, False, False, True])
        basis = ("no source score: condition named in title/summary first, then records, then sample size"
                 if has_mention else "no source score: number of records, then sample size")
    out = []
    for _, r in agg.head(top_n).iterrows():
        out.append({"entity_type": r["entity_type"], "entity_id": r["entity_id"],
                    "entity_label": r["entity_label"], "n_records": int(r["n_records"]),
                    "max_source_score": (None if pd.isna(r["max_source_score"]) else float(r["max_source_score"])),
                    "max_sample_size": (None if pd.isna(r["max_sample_size"]) else float(r["max_sample_size"])),
                    **({"max_clinical_stage": r["max_clinical_stage"]} if "max_clinical_stage" in agg else {}),
                    **({"min_published_adj_p": float(r["min_published_adj_p"])}
                       if "min_published_adj_p" in agg and pd.notna(r["min_published_adj_p"]) else {}),
                    **({"min_reported_p": float(r["min_reported_p"])}
                       if "min_reported_p" in agg and pd.notna(r["min_reported_p"]) else {}),
                    "sources": r["sources"], "example_accessions": r["accessions"], "rank_basis": basis})
    return out


def _rank_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Open Targets: rank genes by their association rows (overall association score), not by evidence items."""
    if "source_evidence_category" in df:
        assoc = df[df["source_evidence_category"] == "overall_association"]
        if len(assoc):
            return assoc
    return df


def _keep_when_excluding_broad(d: pd.DataFrame) -> pd.Series:
    """Rows to keep when broad identifiers are excluded.

    Only rows retrieved *by* a broad identifier (Open Targets / GWAS Catalog queries on an id whose
    registry predicate is 'broad') are dropped. GEO/SRA rows come from curated text queries
    (source_disease_id == ''); their ontology_id is only the condition's primary id as a label, so the
    predicate of that id says nothing about the match and they are kept.
    """
    broad = d["ontology_match"].fillna("").astype(str) == "broad"
    if "source_disease_id" in d:
        id_based = d["source_disease_id"].fillna("").astype(str) != ""
        return ~(broad & id_based)
    return ~broad


def _evidence_stamp() -> tuple:
    from ..store import processed_path
    return tuple(processed_path(n).stat().st_mtime_ns if processed_path(n).exists() else 0
                 for n in (EVIDENCE_TABLE, COVERAGE_TABLE))


@lru_cache(maxsize=2)
def _evidence_tables(stamp: tuple) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """(condition_molecular_evidence, condition_molecular_coverage, genetic_gene_agreement over every condition), read
    and computed once per table version; get_molecular_context only filters them. The agreement is None when it cannot
    be computed (reported as UNKNOWN, never invented)."""
    ev = read_table(EVIDENCE_TABLE)
    cov = read_table(COVERAGE_TABLE) if table_exists(COVERAGE_TABLE) else pd.DataFrame(columns=COVERAGE_COLUMNS)
    try:
        agr = genetic_gene_agreement(ev)
    except Exception:  # a partition missing -> agreement unknown, not invented
        agr = None
    return ev, cov, agr


def get_molecular_context(condition: str, *, top_n: int = 10, include_broad: bool = True) -> dict:
    """Condition-level molecular enrichment/evidence for one condition.

    Returns {"condition": ..., "status": ..., "evidence": {evidence_type: {"n_records", "counts_by_source",
    "counts_by_entity_type", "top_entities"}}, "coverage": [...], "provenance": [...], "caveats": [...]}.
    Missing data -> the UNKNOWN sentinel, never an invented value.
    """
    cid, how = _resolve_condition(condition)
    if cid is None:
        return {"query": condition, "status": UNKNOWN, "reason": how.get("reason", ""), "evidence": UNKNOWN}
    if not table_exists(EVIDENCE_TABLE):
        return {"query": condition, "condition_id": cid, "status": UNKNOWN,
                "reason": f"{EVIDENCE_TABLE} has not been built", "evidence": UNKNOWN}
    ev, cov, agr = _evidence_tables(_evidence_stamp())
    cond_col = "condition_id"
    # every partition, mapMECFS included, carries the common entity/source columns (no per-source view needed)
    d = ev[ev[cond_col] == cid].copy()
    if not include_broad and "ontology_match" in d:
        d = d[_keep_when_excluding_broad(d)]
    cov = cov[cov["condition_id"] == cid]
    reg = _registry().set_index("canonical_condition_id")
    result = {
        "query": condition,
        "condition_id": cid,
        "condition_label": reg.loc[cid, "preferred_name"],
        "resolution": how,
        "layer": LAYER_LABEL,
        "status": "available" if len(d) else UNKNOWN,
        "evidence": {},
        "coverage": [{k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in r.items()}
                     for r in cov[COVERAGE_COLUMNS].to_dict("records")] if len(cov) else UNKNOWN,
        "provenance": [],
        "caveats": [
            "Condition-level molecular enrichment/evidence from public databases; not measured on any "
            "participant and not linkable to person-level wearable data.",
            "No combined omics score is computed; source_score is whatever the source publishes "
            "(Open Targets association/evidence scores) and is not comparable across sources.",
            "Open Targets associations are indirect by default (descendant diseases included).",
            "GEO/SRA rows are study metadata from curated text queries, not analysed data; sample counts "
            "are GEO samples or SRA runs, not participants.",
        ],
    }
    if d.empty:
        result["evidence"] = UNKNOWN
        result["reason"] = "no molecular evidence rows for this condition in any source"
        return result
    for et, grp in d.groupby("evidence_type"):
        result["evidence"][et] = {
            "n_records": int(len(grp)),
            "counts_by_source": {str(k): int(v) for k, v in grp["source_database"].value_counts().items()},
            "counts_by_entity_type": {str(k): int(v) for k, v in grp["entity_type"].value_counts().items()},
            "counts_by_source_category": ({str(k): int(v) for k, v in grp["source_evidence_category"].value_counts().items()}
                                          if "source_evidence_category" in grp else {}),
            "n_distinct_entities": int(grp["entity_id"].nunique()),
            "top_entities": {str(t): _entity_counts(_rank_rows(sub), top_n) for t, sub in grp.groupby("entity_type")},
        }
        if "mentions_condition" in grp and grp["mentions_condition"].notna().any():
            result["evidence"][et]["n_records_naming_condition_in_title_or_summary"] = int(
                grp["mentions_condition"].fillna(False).astype(bool).sum())
    try:
        result["cross_source_genetic_agreement"] = agr[agr["condition_id"] == cid].iloc[0].to_dict()
    except Exception:  # agreement not computable (None) or no row -> unknown, not invented
        result["cross_source_genetic_agreement"] = UNKNOWN
    prov_cols = ["source_name", "source_version", "retrieved_at", "data_layer"]
    prov = d[prov_cols].drop_duplicates().to_dict("records")
    result["provenance"] = prov
    return result


if __name__ == "__main__":  # pragma: no cover
    import sys
    print(json.dumps(get_molecular_context(sys.argv[1] if len(sys.argv) > 1 else "long_covid"), indent=1, default=str)[:20000])
