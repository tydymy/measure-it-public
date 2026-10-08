"""NHGRI-EBI GWAS Catalog REST API v2 -> condition-level genetic-association evidence.

What this module does
---------------------
For every MONDO/EFO identifier the ontology layer recorded for a registry condition
(``omics.query.condition_query_ids``, including EFO ids the ontology module flagged as
obsolete so a Catalog still using them is not silently missed):

1. ``GET /v2/efo-traits/{id}`` — is the id a GWAS Catalog trait, and its label (404 = not a trait).
2. ``GET /v2/associations?efo_id=`` (default) and again with ``show_child_trait=true``;
   paginated. rsID, effect/risk allele, p-value, OR, beta (+ unit/direction), CI,
   risk-allele frequency, mapped genes, reported trait, study accession, PMID.
3. ``GET /v2/associations/{id}/loci`` for every association — author-reported genes.
4. ``GET /v2/studies?efo_id=`` (default + ``show_child_trait=true``) and
   ``GET /v2/studies/{accession}`` for any study referenced by an association but not
   listed; ``GET /v2/studies/{accession}/ancestries`` for structured sample sizes.
5. ``GET /v2/metadata`` — data release date, EFO version, API version.

Rows: one per (condition, queried id, association) with ``entity_type='variant'`` and
one per (condition, queried id, study) with ``entity_type='study'`` — studies with no
curated association are kept, so zero-association results stay visible. Zero-hit
conditions are recorded explicitly in ``condition_molecular_coverage__gwas_catalog``.

Guardrails: condition-level molecular evidence, not participant data; no score is
computed (``source_score`` is null — the Catalog publishes p-values and effect sizes,
kept in their own columns, not a score).

Reproduce: ``uv run python -m measure_it.omics.gwas_catalog``
"""
from __future__ import annotations

import argparse
import json
import re
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode

import numpy as np
import pandas as pd

from ..config import raw_dir
from ..http import get
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from . import query as Q

SOURCE_ID = "gwas_catalog"
SOURCE_DB = "GWAS Catalog"
BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"
LANDING_URL = "https://www.ebi.ac.uk/gwas/"
DOCS_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/docs"
SPEC_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/rest-api-doc.yaml"
MODULE = "measure_it.omics.gwas_catalog"
PARTITION = "condition_molecular_evidence__gwas_catalog"
COVERAGE = "condition_molecular_coverage__gwas_catalog"
PAGE_SIZE = 200
WORKERS = 4


# --------------------------------------------------------------------------- transport
def _get(path: str, params: dict | None = None, *, raw_rel: str | None = None, allow_404: bool = False,
         max_retries: int = 5):
    url = f"{BASE}/{path.lstrip('/')}"
    r = get(url, params=params, expect_json=True, cache_errors=allow_404, max_retries=max_retries)
    if allow_404 and r.status == 404:
        return None, r
    r.raise_for_status()
    payload = r.json()
    if raw_rel:
        Q.save_raw_json(SOURCE_ID, raw_rel, payload,
                        url=url + (("?" + urlencode(params)) if params else ""),
                        request={"method": "GET", "params": params or {}}, retrieved_at=r.fetched_at)
    return payload, r


def fetch_metadata() -> tuple[dict, str]:
    payload, r = _get("metadata", raw_rel="metadata.json")
    return payload, r.fetched_at


def fetch_trait(short_id: str) -> tuple[dict | None, int]:
    payload, r = _get(f"efo-traits/{short_id}", allow_404=True, raw_rel=None)
    Q.save_raw_json(SOURCE_ID, f"traits/{Q.safe_name(short_id)}.json",
                    {"http_status": r.status, "body": payload if payload is not None else r.text[:2000]},
                    url=f"{BASE}/efo-traits/{short_id}", request={"method": "GET"}, retrieved_at=r.fetched_at)
    return payload, r.status


def _paged(path: str, params: dict, embedded_key: str, raw_prefix: str) -> tuple[int, list[dict]]:
    out, page, total = [], 0, 0
    while True:
        p = {**params, "size": PAGE_SIZE, "page": page}
        payload, _ = _get(path, p, raw_rel=f"{raw_prefix}__page{page}.json")
        items = (payload.get("_embedded") or {}).get(embedded_key, [])
        out.extend(items)
        pg = payload.get("page") or {}
        total = pg.get("totalElements", len(out))
        page += 1
        if page >= pg.get("totalPages", 0) or not items:
            break
    return total, out


def fetch_associations(short_id: str, child: bool) -> tuple[int, list[dict]]:
    params = {"efo_id": short_id}
    if child:
        params["show_child_trait"] = "true"
    tag = "child" if child else "default"
    return _paged("associations", params, "associations", f"associations/{Q.safe_name(short_id)}__{tag}")


def fetch_studies(short_id: str, child: bool) -> tuple[int, list[dict]]:
    params = {"efo_id": short_id}
    if child:
        params["show_child_trait"] = "true"
    tag = "child" if child else "default"
    return _paged("studies", params, "studies", f"studies/{Q.safe_name(short_id)}__{tag}")


FETCH_ERRORS: dict[str, str] = {}


def _tolerant(kind: str, key: str, path: str, raw_rel: str) -> dict | None:
    """Per-record detail calls (loci, ancestries): a failure is recorded, never filled in."""
    try:
        payload, _ = _get(path, raw_rel=raw_rel, max_retries=3)
        return payload
    except Exception as exc:  # HTTP 5xx after retries, timeouts
        FETCH_ERRORS[f"{kind}:{key}"] = str(exc)[:200]
        return None


def fetch_loci(association_id) -> dict | None:
    return _tolerant("loci", str(association_id), f"associations/{association_id}/loci", f"loci/{association_id}.json")


def fetch_study(acc: str) -> dict:
    payload, _ = _get(f"studies/{acc}", raw_rel=f"study_records/{Q.safe_name(acc)}.json")
    return payload


def fetch_ancestries(acc: str) -> dict | None:
    return _tolerant("ancestries", acc, f"studies/{acc}/ancestries", f"ancestries/{Q.safe_name(acc)}.json")


# --------------------------------------------------------------------------- parsing (pure)
def reported_genes_from_loci(loci_payload: dict | None) -> list[str]:
    genes = []
    for loc in ((loci_payload or {}).get("_embedded") or {}).get("loci", []):
        for g in loc.get("author_reported_genes") or []:
            name = g.get("gene_name") if isinstance(g, dict) else str(g)
            if name and name not in genes:
                genes.append(name)
    return genes


def ancestry_sample_size(anc_payload: dict | None) -> tuple[float, str]:
    """Sum of number_of_individuals over the study's 'initial' (discovery) ancestry entries."""
    items = ((anc_payload or {}).get("_embedded") or {}).get("ancestries", [])
    init = [a for a in items if (a.get("type") or "").lower() == "initial" and a.get("number_of_individuals") is not None]
    if not init:
        return np.nan, ""
    n = float(sum(a["number_of_individuals"] for a in init))
    groups = []
    for a in init:
        g = ", ".join(x.get("ancestral_group", "") for x in a.get("ancestral_groups") or [])
        groups.append(f"{g or 'NR'} (n={a['number_of_individuals']})")
    return n, "; ".join(groups)


def parse_beta_direction(assoc: dict) -> str | None:
    """Direction of a beta effect as the Catalog reports it ('increase'/'decrease'), else None."""
    d = assoc.get("beta_direction")
    if d:
        return str(d)
    b = assoc.get("beta")
    if isinstance(b, str):
        m = re.search(r"\b(increase|decrease)\b", b)
        if m:
            return m.group(1)
    return None


# REST API v2 association records carry the beta only as text: "<magnitude> [<unit>] increase|decrease"
# (measured on every association retrieved 2026-09-23: 'N unit decrease', 'N unit increase', 'N increase').
# There is no numeric beta / beta_unit / beta_direction field in v2.
_BETA_TEXT = re.compile(r"^\s*(?P<value>[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(?P<unit>.*?)\s*"
                        r"\b(?P<direction>increase|decrease)\b\s*$", re.IGNORECASE)


def parse_beta_text(beta) -> tuple[float, str | None, str | None]:
    """Split the Catalog's beta text into (magnitude as written, unit, direction); NaN/None when absent.

    The magnitude is the number exactly as the Catalog writes it (unsigned; the sign is carried by
    the direction word). Nothing is recomputed or harmonised.
    """
    if not isinstance(beta, str):
        return (to_float(beta), None, None)
    m = _BETA_TEXT.match(beta)
    if not m:
        return (np.nan, None, None)
    return (float(m.group("value")), (m.group("unit") or None), m.group("direction").lower())


def to_float(v) -> float:
    try:
        if v is None or (isinstance(v, str) and v.strip() in ("", "-", "NR", "NA")):
            return np.nan
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def association_row(ctx: dict, a: dict, reported_genes: list[str] | None, study: dict | None,
                    study_n: float, study_pop: str, found_via: str) -> dict:
    """reported_genes=None means the /loci call failed (unknown), [] means the Catalog lists none."""
    snps = a.get("snp_allele") or []
    rsids = [s.get("rs_id") for s in snps if s.get("rs_id")]
    alleles = [s.get("effect_allele") for s in snps if s.get("effect_allele")]
    traits = a.get("efo_traits") or []
    trait_ids = [t.get("efo_id") for t in traits]
    direction = parse_beta_direction(a)
    beta_value, beta_unit, _ = parse_beta_text(a.get("beta"))
    return {
        **ctx,
        "evidence_type": "genetic_association",
        "source_evidence_category": "GWAS association",
        "entity_type": "variant",
        "entity_id": ";".join(rsids) if rsids else ";".join(a.get("snp_effect_allele") or []) or str(a["association_id"]),
        "entity_label": ";".join(a.get("snp_effect_allele") or []) or None,
        "source_database": SOURCE_DB,
        "source_accession": a.get("accession_id"),
        "gwas_association_id": str(a["association_id"]),
        "evidence_direction": (f"beta {direction} per effect allele (as reported)" if direction else None),
        "study_population": study_pop or ("; ".join((study or {}).get("discovery_ancestry") or []) or None),
        "sample_size": study_n,
        "sample_size_basis": ("sum of number_of_individuals over the study's initial ancestry entries (GWAS Catalog)"
                              if not np.isnan(study_n) else None),
        "initial_sample_size_text": (study or {}).get("initial_sample_size"),
        "source_score": np.nan, "source_score_label": None,
        "rsids": ";".join(rsids) or None,
        "risk_allele": ";".join(alleles) or None,
        "snp_effect_allele": ";".join(a.get("snp_effect_allele") or []) or None,
        "risk_allele_frequency": to_float(a.get("risk_frequency")),
        "p_value": to_float(a.get("p_value")),
        "pvalue_mantissa": to_float(a.get("pvalue_mantissa")),
        "pvalue_exponent": to_float(a.get("pvalue_exponent")),
        "pvalue_description": a.get("pvalue_description") or None,
        "or_value": to_float(a.get("or_value")),
        "beta_text": a.get("beta") if a.get("beta") not in (None, "-") else None,
        "beta_value": beta_value,            # magnitude as written in beta_text (unsigned)
        "beta_unit": beta_unit,              # unit word from beta_text (e.g. 'unit'), None if not stated
        "beta_direction": direction,         # 'increase' / 'decrease' per effect allele, from beta_text
        "ci_text": a.get("range"),
        "ci_lower": to_float(a.get("ci_lower")),
        "ci_upper": to_float(a.get("ci_upper")),
        "mapped_genes": ";".join(a.get("mapped_genes") or []) or None,
        "reported_genes": ";".join(reported_genes or []) or None,
        "reported_genes_status": ("loci request failed (unknown)" if reported_genes is None
                                  else ("none reported" if not reported_genes else "reported")),
        "reported_trait": ";".join(a.get("reported_trait") or []) or None,
        "gwas_efo_trait_ids": ";".join(trait_ids) or None,
        "gwas_efo_trait_labels": ";".join(t.get("efo_trait", "") for t in traits) or None,
        "association_trait_is_queried_id": ctx["source_disease_id"] in trait_ids,
        "locations": ";".join(a.get("locations") or []) or None,
        "multi_snp_haplotype": a.get("multi_snp_haplotype"),
        "snp_interaction": a.get("snp_interaction"),
        "pubmed_id": str(a.get("pubmed_id")) if a.get("pubmed_id") else None,
        "first_author": a.get("first_author"),
        "found_via": found_via,
    }


def study_row(ctx: dict, s: dict, study_n: float, study_pop: str, n_assoc: int, found_via: str) -> dict:
    traits = s.get("efo_traits") or []
    return {
        **ctx,
        "evidence_type": "genetic_association",
        "source_evidence_category": "GWAS study",
        "entity_type": "study",
        "entity_id": s.get("accession_id"),
        "entity_label": s.get("disease_trait"),
        "source_database": SOURCE_DB,
        "source_accession": s.get("accession_id"),
        "evidence_direction": None,
        "study_population": study_pop or ("; ".join(s.get("discovery_ancestry") or []) or None),
        "sample_size": study_n,
        "sample_size_basis": ("sum of number_of_individuals over the study's initial ancestry entries (GWAS Catalog)"
                              if not np.isnan(study_n) else None),
        "initial_sample_size_text": s.get("initial_sample_size"),
        "replication_sample_size_text": s.get("replication_sample_size"),
        "discovery_ancestry": "; ".join(s.get("discovery_ancestry") or []) or None,
        "replication_ancestry": "; ".join(s.get("replication_ancestry") or []) or None,
        "cohort": ";".join(s.get("cohort") or []) or None,
        "pubmed_id": str(s.get("pubmed_id")) if s.get("pubmed_id") else None,
        "reported_trait": s.get("disease_trait"),
        "gwas_efo_trait_ids": ";".join(t.get("efo_id", "") for t in traits) or None,
        "gwas_efo_trait_labels": ";".join(t.get("efo_trait", "") for t in traits) or None,
        "association_trait_is_queried_id": ctx["source_disease_id"] in [t.get("efo_id") for t in traits],
        "full_summary_stats_available": s.get("full_summary_stats_available"),
        "n_associations_for_query": n_assoc,
        "genotyping_technologies": ";".join(s.get("genotyping_technologies") or []) or None,
        "source_score": np.nan, "source_score_label": None,
        "found_via": found_via,
    }


# --------------------------------------------------------------------------- pipeline
def run(*, limit_ids: int | None = None) -> dict:
    raw_dir(SOURCE_ID)
    meta, retrieved = fetch_metadata()
    version = (f"GWAS Catalog data release {meta.get('data_release_date')}; REST API {meta.get('version')} "
               f"(api release {meta.get('api_release_date')}); EFO {meta.get('efo_version')}; "
               f"gene build {meta.get('gene_build')}")
    ids = Q.condition_query_ids(include_obsolete_efo=True)
    if limit_ids:
        ids = ids.head(limit_ids)

    per_id: dict[str, dict] = {}
    for sid in sorted(set(ids["source_query_id"]) - {""}):
        trait, status = fetch_trait(sid)
        rec = {"trait_status": status, "trait_label": (trait or {}).get("efo_trait"),
               "assoc": {}, "studies": {}}
        for child in (False, True):
            tag = "child" if child else "default"
            rec["assoc"][tag] = fetch_associations(sid, child)
            rec["studies"][tag] = fetch_studies(sid, child)
        per_id[sid] = rec

    # union associations/studies per id (default first), remember how each was found
    assoc_by_id: dict[str, dict[str, tuple[dict, str]]] = {}
    studies_by_id: dict[str, dict[str, tuple[dict, str]]] = {}
    for sid, rec in per_id.items():
        a_map, s_map = {}, {}
        for tag in ("default", "child"):
            for a in rec["assoc"][tag][1]:
                key = str(a["association_id"])
                if key not in a_map:
                    a_map[key] = (a, "efo_id" if tag == "default" else "efo_id+show_child_trait")
            for s in rec["studies"][tag][1]:
                key = s["accession_id"]
                if key not in s_map:
                    s_map[key] = (s, "efo_id" if tag == "default" else "efo_id+show_child_trait")
        assoc_by_id[sid], studies_by_id[sid] = a_map, s_map

    all_assoc_ids = sorted({k for m in assoc_by_id.values() for k in m})
    study_records: dict[str, dict] = {}
    for m in studies_by_id.values():
        for acc, (s, _) in m.items():
            study_records.setdefault(acc, s)
    referenced = {a["accession_id"] for m in assoc_by_id.values() for a, _ in m.values() if a.get("accession_id")}
    missing_studies = sorted(referenced - set(study_records))

    with ThreadPoolExecutor(WORKERS) as ex:
        loci = dict(zip(all_assoc_ids, ex.map(fetch_loci, all_assoc_ids)))
        for acc, rec in zip(missing_studies, ex.map(fetch_study, missing_studies)):
            study_records[acc] = rec
        accs = sorted(study_records)
        anc = dict(zip(accs, ex.map(fetch_ancestries, accs)))
    study_n = {acc: ancestry_sample_size(anc.get(acc)) for acc in study_records}

    rows, coverage = [], []
    for _, idr in ids.iterrows():
        sid = idr["source_query_id"]
        base = {"condition_id": idr["condition_id"], "ontology_id": idr["ontology_id"],
                "ontology_id_role": idr["ontology_id_role"], "ontology_match": idr["ontology_match"],
                "source_database": SOURCE_DB}
        if not sid:
            coverage.append({**base, "id_status_in_source": "no identifier", "source_disease_label": None,
                             "coverage_status": "no identifier", "n_records": 0, "detail": "no MONDO/EFO id"})
            continue
        rec = per_id[sid]
        ctx = {"condition_id": idr["condition_id"], "ontology_id": idr["ontology_id"],
               "ontology_id_role": idr["ontology_id_role"], "ontology_match": idr["ontology_match"],
               "ontology_id_status": idr["id_status"],
               "source_disease_id": sid, "source_disease_label": rec["trait_label"],
               "registry_label": idr["registry_label"], "date_retrieved": retrieved[:10]}
        a_map, s_map = assoc_by_id[sid], studies_by_id[sid]
        n_assoc_by_study: dict[str, int] = {}
        for key, (a, via) in a_map.items():
            acc = a.get("accession_id")
            n_assoc_by_study[acc] = n_assoc_by_study.get(acc, 0) + 1
            n, pop = study_n.get(acc, (np.nan, ""))
            lp = loci.get(key)
            rows.append(association_row(ctx, a, reported_genes_from_loci(lp) if lp is not None else None,
                                        study_records.get(acc), n, pop, via))
        study_accs = dict(s_map)
        for acc in n_assoc_by_study:  # studies reached only through an association
            if acc and acc not in study_accs and acc in study_records:
                study_accs[acc] = (study_records[acc], "referenced by association")
        for acc, (s, via) in study_accs.items():
            n, pop = study_n.get(acc, (np.nan, ""))
            rows.append(study_row(ctx, s, n, pop, n_assoc_by_study.get(acc, 0), via))
        d_cnt = rec["assoc"]["default"][0]
        c_cnt = rec["assoc"]["child"][0]
        s_d, s_c = rec["studies"]["default"][0], rec["studies"]["child"][0]
        if rec["trait_status"] == 404:
            status, id_status = "identifier not present in source", "not a GWAS Catalog trait (HTTP 404)"
        elif len(a_map) == 0 and len(s_map) == 0:
            status, id_status = "zero hits", "present"
        elif len(a_map) == 0:
            status, id_status = "studies only, zero associations", "present"
        else:
            status, id_status = "queried", "present"
        coverage.append({**base, "id_status_in_source": id_status, "source_disease_label": rec["trait_label"],
                         "coverage_status": status, "n_records": len(a_map),
                         "detail": json.dumps({"trait_http_status": rec["trait_status"],
                                               "associations_default": d_cnt, "associations_show_child_trait": c_cnt,
                                               "associations_union": len(a_map),
                                               "studies_default": s_d, "studies_show_child_trait": s_c,
                                               "studies_union": len(s_map)})})

    df = pd.DataFrame(rows)
    notes = ("Condition-level genetic-association evidence from the GWAS Catalog; not participant-linked. "
             "p-values/OR/beta are as curated by the Catalog; no score is computed. sample_size = sum of the "
             "study's initial-ancestry number_of_individuals.")
    if len(df):
        df = add_provenance(df.drop(columns=["evidence_type"]), data_layer="condition_molecular",
                            source_name=f"{SOURCE_DB} REST API v2", source_version=version, retrieved_at=retrieved,
                            evidence_type="genetic_association",
                            source_record_id=lambda x: x["source_disease_id"] + "|" + x["entity_type"] + "|"
                            + x.get("gwas_association_id", pd.Series("", index=x.index)).fillna(x["source_accession"]).astype(str),
                            evidence_level="curated_gwas_catalog_record", provenance_notes=notes)
        out = Q.stringify_object_columns(Q.finalize_columns(df))
    else:
        out = Q.finalize_columns(pd.DataFrame(columns=Q.EVIDENCE_COLUMNS))
    write_table(out, PARTITION, producer=MODULE,
                description="GWAS Catalog associations and studies per registry condition id")
    cov = pd.DataFrame(coverage)
    cov = _condition_level(cov)
    cov = add_provenance(cov, data_layer="condition_molecular", source_name=f"{SOURCE_DB} REST API v2",
                         source_version=version, retrieved_at=retrieved, evidence_type="metadata_catalog",
                         source_record_id=lambda x: x["condition_id"] + "|" + x["ontology_id"].astype(str),
                         provenance_notes="coverage of registry identifiers in the GWAS Catalog")
    write_table(Q.stringify_object_columns(cov), COVERAGE, producer=MODULE,
                description="GWAS Catalog coverage per condition x registry identifier (zero hits explicit)")
    run_log = {"metadata": meta, "version": version, "retrieved_at": retrieved,
               "n_ids_queried": len(per_id), "n_associations_unique": len(all_assoc_ids),
               "n_loci_calls": len(loci), "n_studies": len(study_records),
               "n_studies_fetched_individually": len(missing_studies),
               "detail_call_failures": dict(FETCH_ERRORS),
               "per_id": {sid: {"trait_status": r["trait_status"], "trait_label": r["trait_label"],
                                "associations_default": r["assoc"]["default"][0],
                                "associations_child": r["assoc"]["child"][0],
                                "studies_default": r["studies"]["default"][0],
                                "studies_child": r["studies"]["child"][0]} for sid, r in per_id.items()}}
    (raw_dir(SOURCE_ID) / "run_log.json").write_text(json.dumps(run_log, indent=1, default=str))
    summary = summarize(out, cov)
    write_audit(out, cov, run_log, summary)
    write_registry(out, cov, run_log, summary)
    union = Q.build_union()
    return {"rows": len(out), "version": version, "union": union, "zero_association_conditions": summary["zero_assoc"],
            "zero_evidence_conditions": summary["zero_any"]}


def _condition_level(cov: pd.DataFrame) -> pd.DataFrame:
    extra = []
    for cid, g in cov.groupby("condition_id"):
        if (g["id_status_in_source"] == "present").sum() == 0:
            extra.append({"condition_id": cid, "ontology_id": "(all registry ids)", "ontology_id_role": "condition",
                          "ontology_match": "", "source_database": SOURCE_DB,
                          "id_status_in_source": "no identifier usable by source", "source_disease_label": None,
                          "coverage_status": "no identifier", "n_records": 0,
                          "detail": "no registry id is a GWAS Catalog trait: " + ";".join(g["ontology_id"].astype(str))})
    return pd.concat([cov, pd.DataFrame(extra)], ignore_index=True) if extra else cov


def summarize(df: pd.DataFrame, cov: pd.DataFrame) -> dict:
    per = []
    for cid in Q._registry()["canonical_condition_id"]:
        d = df[df["condition_id"] == cid] if len(df) else df
        a = d[d["entity_type"] == "variant"] if len(d) else d
        s = d[d["entity_type"] == "study"] if len(d) else d
        genes = set()
        for g in a.get("mapped_genes", pd.Series(dtype=str)).dropna():
            genes.update(x for x in str(g).split(";") if x)
        rgenes = set()
        for g in a.get("reported_genes", pd.Series(dtype=str)).dropna():
            rgenes.update(x for x in str(g).split(";") if x)
        c = cov[cov["condition_id"] == cid]
        per.append({
            "condition_id": cid,
            "gwas_trait_ids": ";".join(sorted(c[c["id_status_in_source"] == "present"]["ontology_id"])) or "none",
            "associations": int(a["gwas_association_id"].nunique()) if len(a) else 0,
            "associations_p_lt_5e-8": int(a[a["p_value"] < 5e-8]["gwas_association_id"].nunique()) if len(a) else 0,
            "unique_rsids": int(a["entity_id"].nunique()) if len(a) else 0,
            "mapped_genes": len(genes), "reported_genes": len(rgenes),
            "studies": int(s["entity_id"].nunique()) if len(s) else 0,
            "studies_with_associations": int(s[s["n_associations_for_query"] > 0]["entity_id"].nunique()) if len(s) else 0,
            "pubmed_ids": int(d["pubmed_id"].dropna().nunique()) if len(d) else 0,
        })
    per = pd.DataFrame(per)
    return {"per_condition": per,
            "zero_assoc": per[per["associations"] == 0]["condition_id"].tolist(),
            "zero_any": per[(per["associations"] == 0) & (per["studies"] == 0)]["condition_id"].tolist()}


def write_audit(df: pd.DataFrame, cov: pd.DataFrame, run_log: dict, summary: dict) -> None:
    from .open_targets import _md_table, _missingness
    per = summary["per_condition"]
    a = df[df["entity_type"] == "variant"]
    s = df[df["entity_type"] == "study"]
    meta = run_log["metadata"]
    L = ["# DATA AUDIT — NHGRI-EBI GWAS Catalog (REST API v2)\n", "| Field | Value |", "|---|---|",
         f"| source_id | {SOURCE_ID} |",
         f"| Source (dataset/API name, exact files/endpoints) | GWAS Catalog REST API v2: GET {BASE}/metadata, "
         "/efo-traits/{id}, /associations?efo_id=[&show_child_trait=true], /associations/{id}/loci, "
         "/studies?efo_id=[&show_child_trait=true], /studies/{acc}, /studies/{acc}/ancestries (OpenAPI: "
         f"{SPEC_URL}) |",
         "| Publishing organization | NHGRI-EBI GWAS Catalog (EMBL-EBI / NHGRI) |",
         f"| Retrieval date (UTC) | {Q.retrieval_window(SOURCE_ID)} |",
         f"| Source version / release | {run_log['version']} |",
         f"| Source update date / cadence | data release {meta.get('data_release_date')} (weekly Catalog releases) |",
         "| License / access conditions | open API, no key; EMBL-EBI terms of use; API code Apache-2.0 |",
         "| Unit of observation | curated SNP-trait association (variant row) and GWAS study (study row) per registry id |",
         f"| Sample size (actual, as ingested) | {len(df)} rows: {len(a)} association rows "
         f"({a['gwas_association_id'].nunique() if len(a) else 0} unique associations), {len(s)} study rows "
         f"({s['entity_id'].nunique() if len(s) else 0} unique GCST accessions) |",
         "| Geography (resolution, vintage) | none (concept-level); study ancestry/country metadata only |",
         "| Person-level? | no |", "| Geographic? | no |",
         "| Omics? | yes — condition-level genetic association evidence only |", "| Wearable? | no |",
         "| True participant linkage across modalities? | no — published summary associations; no participants |\n",
         "## Files / endpoints retrieved\n",
         f"Every response is saved under `data/raw/{SOURCE_ID}/` (metadata.json, traits/, associations/, studies/, "
         "loci/, study_records/, ancestries/) with URL, parameters and sha256 in MANIFEST.json; `run_log.json` has "
         f"per-id counts. Page size {PAGE_SIZE}. {run_log['n_loci_calls']} loci calls, {run_log['n_studies']} study "
         f"ancestry calls ({run_log['n_studies_fetched_individually']} studies fetched individually because only an "
         f"association referenced them). Detail calls that failed after retries (recorded, not filled in): "
         f"{len(run_log['detail_call_failures'])} — "
         + (", ".join(sorted(run_log["detail_call_failures"])[:20]) or "none")
         + ". Rows whose /loci call failed carry `reported_genes_status = 'loci request failed (unknown)'`.\n",
         "## Per-id results (which registry ids are GWAS Catalog traits)\n",
         _md_table(pd.DataFrame([{"source_id": k, **v} for k, v in run_log["per_id"].items()])),
         "\n`show_child_trait=true` was run alongside the default query; the union is kept and `found_via` records "
         "which query returned each row. Measured: for MONDO_0001292 the child-trait query returned fewer "
         "associations than the default, so neither query alone is complete.\n",
         "## Per-condition counts\n",
         _md_table(per),
         f"\nConditions with zero GWAS associations: {', '.join(summary['zero_assoc']) or 'none'}.",
         f"\nConditions with zero GWAS associations AND zero studies: {', '.join(summary['zero_any']) or 'none'}.\n",
         "## Agreement with Open Targets genetic-association targets (counts only)\n",
         "Genes with a non-zero Open Targets `genetic_association` score vs GWAS Catalog mapped genes, per condition, "
         "compared by symbol (`omics.query.genetic_gene_agreement`). Differences reflect release timing ("
         f"{_ot_version()} vs this Catalog release), Open Targets' credible-set/L2G gene assignment vs the Catalog's "
         "positional mapping, and the Catalog's inclusion of sub-threshold associations.\n",
         _md_table(_agreement()),
         "\n## Coverage table (every registry id, including obsolete EFO ids)\n",
         _md_table(cov[["condition_id", "ontology_id", "ontology_id_role", "ontology_match", "id_status_in_source",
                        "source_disease_label", "coverage_status", "n_records"]]),
         "\n## Key variables\n",
         "* variant rows: `entity_id` = rsID(s); `risk_allele` = effect allele as curated; `p_value` (+ mantissa/"
         "exponent, `pvalue_description`); `or_value`; `beta_text` (the Catalog's text, e.g. '0.042 unit decrease'; "
         "REST v2 has no numeric beta field) split into `beta_value` (magnitude as written, unsigned), `beta_unit` "
         "and `beta_direction`; `ci_text`/`ci_lower`/`ci_upper`; `risk_allele_frequency`; `mapped_genes` (Ensembl-mapped by the Catalog); `reported_genes` "
         "(author-reported, from /loci); `reported_trait`; `source_accession` = GCST study; `pubmed_id`; "
         "`evidence_direction` only when the Catalog reports a beta direction (increase/decrease per effect allele).",
         "* study rows: `entity_id` = GCST accession; `initial_sample_size_text`; `discovery_ancestry`; `cohort`; "
         "`n_associations_for_query` (0 = the study is indexed to the trait but has no curated association); "
         "`full_summary_stats_available`.",
         "* `sample_size` = sum of `number_of_individuals` over the study's initial-ancestry entries (structured); "
         "`initial_sample_size_text` keeps the curated free text (cases/controls).",
         "* `source_score` is null: the Catalog publishes statistics, not a score.\n",
         "## Missingness\n", "Association rows:\n",
         _md_table(_missingness(a, ["rsids", "risk_allele", "p_value", "or_value", "beta_value", "beta_text",
                                    "ci_text", "risk_allele_frequency", "mapped_genes", "reported_genes",
                                    "sample_size", "study_population", "evidence_direction", "pubmed_id"])) if len(a) else "(none)",
         "\nStudy rows:\n",
         _md_table(_missingness(s, ["sample_size", "initial_sample_size_text", "discovery_ancestry", "cohort",
                                    "pubmed_id"])) if len(s) else "(none)",
         "\n## Linkage strategy\n",
         "Registry condition -> MONDO/EFO ids (ontology module) -> GWAS Catalog `efo_id` (same CURIE with '_'). "
         "Obsolete EFO ids recorded by the ontology module were also tried; all are reported in the coverage table. "
         "No free-text trait search is used, so 'zero hits' means no curated association/study is mapped to the id. "
         "Variants join to Open Targets evidence by rsID/GCST accession; nothing joins to participants.\n",
         "## Limitations and caveats\n",
         "* Curated top associations only (Catalog inclusion criteria, typically p < 1e-5 as reported); absence "
         "of associations is not evidence of no genetic effect — several studies are indexed with summary "
         "statistics but no curated association (e.g. Long COVID).",
         "* Large PheWAS-style studies (e.g. Million Veteran Program PheCodes) contribute many associations with "
         "no author-reported genes.",
         "* ORs and betas are per the curated effect allele; strand/allele harmonisation is not performed.",
         "* Ancestry composition is predominantly European in most studies (see `study_population`).",
         "* Condition-level molecular evidence only — never patient multi-omics.\n",
         "## Processed outputs\n",
         f"* `{PARTITION}`: {len(df)} rows ({len(a)} variant, {len(s)} study)",
         f"* `{COVERAGE}`: {len(cov)} rows",
         "* `condition_molecular_evidence` / `condition_molecular_coverage`: union of all partitions\n",
         "## Reproduce\n", f"`uv run python -m {MODULE}`\n"]
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text("\n".join(L))


def _agreement() -> pd.DataFrame:
    """Cross-source gene agreement computed from the current partitions (OT + this GWAS partition)."""
    from ..store import read_table, table_exists
    parts = [read_table(n) for n in ("condition_molecular_evidence__open_targets", PARTITION) if table_exists(n)]
    if len(parts) < 2:
        return pd.DataFrame([{"note": "Open Targets partition not built; agreement not computed"}])
    return Q.genetic_gene_agreement(pd.concat(parts, ignore_index=True, sort=False))


def _ot_version() -> str:
    from ..store import read_table, table_exists
    n = "condition_molecular_evidence__open_targets"
    return str(read_table(n, columns=["source_version"])["source_version"].iloc[0]) if table_exists(n) else "(not built)"


def write_registry(df: pd.DataFrame, cov: pd.DataFrame, run_log: dict, summary: dict) -> None:
    per = summary["per_condition"].set_index("condition_id")
    meta = run_log["metadata"]
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "NHGRI-EBI GWAS Catalog REST API v2",
        "publisher": "NHGRI-EBI GWAS Catalog (EMBL-EBI / NHGRI)",
        "landing_url": LANDING_URL,
        "access_urls": [BASE, DOCS_URL, SPEC_URL],
        "license": "EMBL-EBI terms of use (data); Apache-2.0 (API code)",
        "access_conditions": "open API, no key",
        "retrieved_at": run_log["retrieved_at"],
        "source_version": run_log["version"],
        "update_date": str(meta.get("data_release_date")),
        "data_layer": "condition_molecular",
        "unit_of_observation": "curated SNP-trait association; GWAS study",
        "sample_size": {
            "rows": int(len(df)),
            "unique_associations": int(df.loc[df["entity_type"] == "variant", "gwas_association_id"].nunique()) if len(df) else 0,
            "unique_studies": int(df.loc[df["entity_type"] == "study", "entity_id"].nunique()) if len(df) else 0,
            "per_condition_associations": {k: int(v) for k, v in per["associations"].items()},
            "per_condition_studies": {k: int(v) for k, v in per["studies"].items()},
        },
        "geographic_resolution": "none",
        "person_level": False, "geographic": False, "omics": True, "wearable": False,
        "participant_linkage": "not applicable (published summary associations)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": [PARTITION, COVERAGE, "condition_molecular_evidence", "condition_molecular_coverage"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": MODULE,
        "limitations": [
            "Curated top associations only; zero associations is not evidence of no genetic effect",
            "Zero-association conditions: " + (", ".join(summary["zero_assoc"]) or "none"),
            "Zero association and zero study conditions: " + (", ".join(summary["zero_any"]) or "none"),
            "Effect sizes per curated effect allele; no harmonisation",
            "Condition-level molecular evidence only; not participant-linked",
        ],
    })


def main() -> None:  # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--limit-ids", type=int, default=None)
    args = ap.parse_args()
    print(json.dumps(run(limit_ids=args.limit_ids), indent=1, default=str))
    Q.exit_cleanly()


if __name__ == "__main__":  # pragma: no cover
    main()
