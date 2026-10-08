"""NCBI GEO / SRA study catalogue for the registry conditions (E-utilities; metadata only).

What this module does
---------------------
1. Reads the curated, reviewed queries in ``configs/geo_queries.yaml`` (one or more per
   condition) and runs each through Entrez ``esearch`` on ``db=gds`` (Homo sapiens GEO
   series) and ``db=sra`` (Homo sapiens SRA experiment packages). The exact term sent, the
   hit count, Entrez's query translation and any quoted phrase Entrez reports as not found
   are logged (``data/raw/ncbi_geo_sra/query_log.json``).
2. ``esummary`` for every hit. No study data are downloaded.
3. ``geo_study_catalog``: one row per condition x GSE (title, summary, platform(s), number of
   samples, series type, PubMed ids, release date, BioProject, matched condition, the
   queries that found it, whether the series type is in scope, and an automatic flag
   whether title/summary mention the condition).
   ``sra_study_catalog``: SRA experiments aggregated to condition x study (SRP/ERP/DRP):
   BioProject, study title, library strategies/sources, platforms, number of matched
   experiments and runs, and whether the BioProject is also a GEO series' BioProject.
4. ``condition_molecular_evidence__geo``: study-level metadata rows (entity_type 'study')
   from both catalogues in the common condition_molecular schema.
5. Precision check for the demo cluster (long_covid, me_cfs, pots): a seeded sample of up
   to 10 in-scope GEO series per condition is written to ``precision_sample.csv``; the
   manual review (``precision_review.yaml``, one verdict + reason per GSE) is read back and
   precision is reported. Series not yet reviewed are counted as unreviewed, never guessed.

Guardrails: study metadata are condition-level molecular *evidence availability*, not
analysed data and not participant-linked. ``sample_size`` is the GEO sample (GSM) count
or the SRA run count, not a count of participants.

Reproduce: ``uv run python -m measure_it.omics.geo_sra``
"""
from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np
import pandas as pd
import yaml

from ..config import SEED, load_config, raw_dir
from ..http import get
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from . import query as Q

SOURCE_ID = "ncbi_geo_sra"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
LANDING_URL = "https://www.ncbi.nlm.nih.gov/geo/"
MODULE = "measure_it.omics.geo_sra"
PARTITION = "condition_molecular_evidence__geo"
COVERAGE = "condition_molecular_coverage__geo"
GEO_TABLE = "geo_study_catalog"
SRA_TABLE = "sra_study_catalog"
CONFIG_NAME = "geo_queries"
RETMAX = 10000
SUMMARY_BATCH = 200
PRECISION_N = 10


# --------------------------------------------------------------------------- transport
def eutil(endpoint: str, params: dict, raw_rel: str | None = None) -> tuple[dict, str]:
    url = f"{EUTILS}/{endpoint}.fcgi"
    r = get(url, params={**params, "retmode": "json"}, expect_json=True)
    r.raise_for_status()
    payload = r.json()
    if raw_rel:
        Q.save_raw_json(SOURCE_ID, raw_rel, payload, url=url, request={"method": "GET", "params": params},
                        retrieved_at=r.fetched_at)
    return payload, r.fetched_at


def einfo(db: str) -> dict:
    payload, _ = eutil("einfo", {"db": db}, raw_rel=f"einfo_{db}.json")
    info = payload["einforesult"]["dbinfo"][0]
    return {"db": db, "count": info.get("count"), "lastupdate": info.get("lastupdate"),
            "dbbuild": info.get("dbbuild")}


def esearch(db: str, term: str, tag: str) -> dict:
    payload, ts = eutil("esearch", {"db": db, "term": term, "retmax": RETMAX}, raw_rel=f"esearch/{db}__{tag}.json")
    res = payload["esearchresult"]
    warn = res.get("warninglist") or {}
    return {"db": db, "term": term, "count": int(res.get("count", 0)), "ids": res.get("idlist", []),
            "querytranslation": res.get("querytranslation"),
            "phrases_not_found": warn.get("quotedphrasesnotfound") or [],
            "phrases_ignored": warn.get("phrasesignored") or [],
            "retrieved_at": ts}


def esummary(db: str, ids: list[str]) -> dict[str, dict]:
    out = {}
    ids = sorted(set(ids), key=int)
    for i in range(0, len(ids), SUMMARY_BATCH):
        chunk = ids[i:i + SUMMARY_BATCH]
        payload, _ = eutil("esummary", {"db": db, "id": ",".join(chunk)},
                           raw_rel=f"esummary/{db}__batch{i // SUMMARY_BATCH:04d}.json")
        res = payload.get("result", {})
        for uid in res.get("uids", []):
            out[uid] = res[uid]
    return out


# --------------------------------------------------------------------------- pure helpers
def build_term(term: str, db: str, cfg: dict) -> str:
    filt = cfg["geo_filter"] if db == "gds" else cfg["sra_filter"]
    return f"({term}) AND {filt}"


def mentions(text: str, patterns: list[str]) -> list[str]:
    """Which QC patterns (case-insensitive regex) occur in text."""
    t = str(text or "")
    return [p for p in patterns if re.search(p, t, flags=re.IGNORECASE)]


def parse_sra_expxml(expxml: str, runs_xml: str) -> dict:
    """Parse the XML fragments of one SRA esummary record."""
    root = ET.fromstring(f"<root>{expxml}</root>")
    def attr(path, name):
        el = root.find(path)
        return el.get(name) if el is not None else None
    def text(path):
        el = root.find(path)
        return el.text.strip() if el is not None and el.text else None
    runs_root = ET.fromstring(f"<root>{runs_xml or ''}</root>")
    run_accs = [r.get("acc") for r in runs_root.findall("Run") if r.get("acc")]
    return {
        "experiment_acc": attr("Experiment", "acc"),
        "experiment_title": text("Summary/Title"),
        "study_acc": attr("Study", "acc"),
        "study_title": attr("Study", "name"),
        "bioproject": text("Bioproject"),
        "biosample": text("Biosample"),
        "organism": attr("Organism", "ScientificName"),
        "platform": text("Summary/Platform"),
        "instrument": attr("Summary/Platform", "instrument_model"),
        "library_strategy": text("Library_descriptor/LIBRARY_STRATEGY"),
        "library_source": text("Library_descriptor/LIBRARY_SOURCE"),
        "library_selection": text("Library_descriptor/LIBRARY_SELECTION"),
        "center_name": attr("Submitter", "center_name"),
        "run_accs": run_accs,
    }


def geo_record(uid: str, s: dict) -> dict:
    gpl = [f"GPL{x}" for x in str(s.get("gpl") or "").split(";") if x]
    pmids = [str(x) for x in (s.get("pubmedids") or [])]
    return {
        "geo_uid": uid,
        "gse": s.get("accession"),
        "title": s.get("title"),
        "summary": s.get("summary"),
        "platform": ";".join(gpl) or None,
        "n_samples": s.get("n_samples"),
        "study_type": s.get("gdstype"),
        "pubmed_ids": ";".join(pmids) or None,
        "release_date": (str(s.get("pdat") or "").replace("/", "-") or None),
        "taxon": s.get("taxon"),
        "entry_type": s.get("entrytype"),
        "bioproject": s.get("bioproject") or None,
        "supplementary_file_types": s.get("suppfile") or None,
        "gds_accession": s.get("gds") or None,
    }


# --------------------------------------------------------------------------- pipeline
def run() -> dict:
    raw_dir(SOURCE_ID)
    cfg = load_config(CONFIG_NAME)
    study_types = set(cfg["geo_study_types"])
    info = {db: einfo(db) for db in ("gds", "sra")}
    version = (f"NCBI Entrez db=gds last update {info['gds']['lastupdate']} ({info['gds']['count']} records); "
               f"db=sra last update {info['sra']['lastupdate']} ({info['sra']['count']} records)")
    reg_ids = set(Q._registry()["canonical_condition_id"])

    qlog, hits = [], defaultdict(lambda: defaultdict(list))   # hits[db][(cid, uid)] -> [query ids]
    retrieved = None
    for cid, queries in cfg["conditions"].items():
        if cid not in reg_ids:
            raise ValueError(f"{CONFIG_NAME}.yaml condition {cid!r} is not in condition_registry")
        for q in queries:
            for db in q["db"]:
                term = build_term(q["term"], db, cfg)
                res = esearch(db, term, f"{cid}__{q['id']}")
                retrieved = retrieved or res["retrieved_at"]
                if res["count"] > RETMAX:
                    raise RuntimeError(f"{q['id']} {db}: {res['count']} hits exceed RETMAX {RETMAX}")
                qlog.append({"condition_id": cid, "query_id": q["id"], "db": db, "term": term,
                             "count": res["count"], "querytranslation": res["querytranslation"],
                             "phrases_not_found": res["phrases_not_found"], "rationale": q.get("rationale", "")})
                for uid in res["ids"]:
                    hits[db][(cid, uid)].append(q["id"])
    missing_conditions = sorted(reg_ids - set(cfg["conditions"]))

    geo_sum = esummary("gds", [u for (_, u) in hits["gds"]])
    sra_sum = esummary("sra", [u for (_, u) in hits["sra"]])
    qterms = {(x["query_id"], x["db"]): x["term"] for x in qlog}

    # ---- GEO catalogue
    geo_rows = []
    for (cid, uid), qids in hits["gds"].items():
        s = geo_sum.get(uid)
        if not s:
            continue
        rec = geo_record(uid, s)
        if rec["entry_type"] != "GSE":
            continue
        found = mentions(f"{rec['title']} {rec['summary']}", cfg["qc_patterns"].get(cid, []))
        geo_rows.append({"condition_id": cid, **rec,
                         "study_type_in_scope": any(t.strip() in study_types for t in str(rec["study_type"]).split(";")),
                         "study_type_all_in_scope": all(t.strip() in study_types for t in str(rec["study_type"]).split(";")),
                         "mentions_condition": bool(found), "qc_patterns_matched": ";".join(found) or None,
                         "query_ids": ";".join(sorted(set(qids))),
                         "query_terms": " || ".join(qterms[(q, "gds")] for q in sorted(set(qids)))})
    geo = pd.DataFrame(geo_rows)

    # ---- SRA catalogue (experiments -> study)
    exp_rows = []
    for (cid, uid), qids in hits["sra"].items():
        s = sra_sum.get(uid)
        if not s:
            continue
        p = parse_sra_expxml(s.get("expxml", ""), s.get("runs", ""))
        exp_rows.append({"condition_id": cid, "sra_uid": uid, **p, "createdate": s.get("createdate"),
                         "query_ids": sorted(set(qids))})
    exps = pd.DataFrame(exp_rows)
    geo_bioprojects = defaultdict(set)
    for _, r in geo.iterrows():
        if r["bioproject"]:
            geo_bioprojects[r["bioproject"]].add(r["gse"])
    sra_rows = []
    if len(exps):
        for (cid, study), g in exps.groupby(["condition_id", "study_acc"], dropna=False):
            runs = sorted({a for lst in g["run_accs"] for a in lst})
            bps = sorted(set(g["bioproject"].dropna()))
            qids = sorted({q for lst in g["query_ids"] for q in lst})
            title = g["study_title"].dropna().iloc[0] if g["study_title"].notna().any() else None
            found = mentions(title or "", cfg["qc_patterns"].get(cid, []))
            sra_rows.append({
                "condition_id": cid, "study_acc": study, "bioproject": ";".join(bps) or None,
                "study_title": title,
                "library_strategies": ";".join(sorted(set(g["library_strategy"].dropna()))) or None,
                "library_sources": ";".join(sorted(set(g["library_source"].dropna()))) or None,
                "platforms": ";".join(sorted(set(g["platform"].dropna()))) or None,
                "center_names": ";".join(sorted(set(g["center_name"].dropna())))[:500] or None,
                "n_experiments_matched": int(len(g)),
                "n_runs_matched": len(runs),
                "n_biosamples_matched": int(g["biosample"].nunique()),
                "first_createdate": min(g["createdate"].dropna()) if g["createdate"].notna().any() else None,
                "last_createdate": max(g["createdate"].dropna()) if g["createdate"].notna().any() else None,
                "linked_geo_series": ";".join(sorted({x for b in bps for x in geo_bioprojects.get(b, ())})) or None,
                "mentions_condition": bool(found), "qc_patterns_matched": ";".join(found) or None,
                "query_ids": ";".join(qids),
                "query_terms": " || ".join(qterms[(q, "sra")] for q in qids),
            })
    sra = pd.DataFrame(sra_rows)

    # ---- precision check (demo cluster)
    precision = precision_check(geo)

    # ---- write catalogues
    retrieved = retrieved or pd.Timestamp.utcnow().isoformat()
    notes = ("GEO/SRA study metadata found by curated Entrez text queries (configs/geo_queries.yaml); metadata only, "
             "no data downloaded; condition-level evidence availability, not participant-linked.")
    geo_out = add_provenance(geo, data_layer="condition_molecular", source_name="NCBI GEO (Entrez E-utilities db=gds)",
                             source_version=version, retrieved_at=retrieved, evidence_type="expression_study_metadata",
                             source_record_id=lambda x: x["condition_id"] + "|" + x["gse"],
                             evidence_level="text_query_match", provenance_notes=notes)
    write_table(Q.stringify_object_columns(geo_out), GEO_TABLE, producer=MODULE,
                description="GEO series per condition from curated Entrez queries (metadata only)")
    sra_out = add_provenance(sra if len(sra) else pd.DataFrame(columns=["condition_id", "study_acc"]),
                             data_layer="condition_molecular", source_name="NCBI SRA (Entrez E-utilities db=sra)",
                             source_version=version, retrieved_at=retrieved, evidence_type="metadata_catalog",
                             source_record_id=lambda x: x["condition_id"] + "|" + x["study_acc"].astype(str),
                             evidence_level="text_query_match", provenance_notes=notes)
    write_table(Q.stringify_object_columns(sra_out), SRA_TABLE, producer=MODULE,
                description="SRA studies per condition (experiments aggregated) from curated Entrez queries")

    ev = evidence_rows(geo, sra, retrieved)
    ev_frames = []
    for (et, src), g in ev.groupby(["evidence_type", "source_name_"]):
        ev_frames.append(add_provenance(
            g.drop(columns=["evidence_type", "source_name_"]), data_layer="condition_molecular", source_name=src,
            source_version=version, retrieved_at=retrieved, evidence_type=et,
            source_record_id=lambda x: x["condition_id"] + "|" + x["source_accession"].astype(str),
            evidence_level="text_query_match", provenance_notes=notes))
    ev_out = Q.stringify_object_columns(Q.finalize_columns(pd.concat(ev_frames, ignore_index=True)))
    write_table(ev_out, PARTITION, producer=MODULE,
                description="GEO series and SRA studies as study-level condition_molecular evidence rows")

    cov = coverage_rows(geo, sra, qlog, missing_conditions)
    cov = add_provenance(cov, data_layer="condition_molecular", source_name="NCBI GEO/SRA (Entrez E-utilities)",
                         source_version=version, retrieved_at=retrieved, evidence_type="metadata_catalog",
                         source_record_id=lambda x: x["condition_id"] + "|" + x["source_database"],
                         provenance_notes="per-condition GEO/SRA study counts from curated queries")
    write_table(Q.stringify_object_columns(cov), COVERAGE, producer=MODULE,
                description="GEO/SRA coverage per condition (zero hits explicit)")

    (raw_dir(SOURCE_ID) / "query_log.json").write_text(json.dumps(qlog, indent=1))
    run_log = {"einfo": info, "version": version, "retrieved_at": retrieved, "precision": precision,
               "n_geo_uids": len(geo_sum), "n_sra_uids": len(sra_sum), "missing_conditions": missing_conditions}
    (raw_dir(SOURCE_ID) / "run_log.json").write_text(json.dumps(run_log, indent=1, default=str))
    summary = summarize(geo, sra)
    write_audit(geo, sra, ev_out, cov, qlog, run_log, summary)
    write_registry(geo, sra, ev_out, run_log, summary)
    union = Q.build_union()
    return {"geo_rows": len(geo), "sra_rows": len(sra), "evidence_rows": len(ev_out), "union": union,
            "precision": {k: v.get("precision") for k, v in precision.items()}}


def evidence_rows(geo: pd.DataFrame, sra: pd.DataFrame, retrieved: str) -> pd.DataFrame:
    ids = Q.condition_query_ids()
    prim = ids[ids["is_primary"]].set_index("condition_id")
    def onto(cid):
        if cid in prim.index:
            r = prim.loc[cid]
            return r["ontology_id"], r["ontology_id_role"], r["ontology_match"]
        return None, None, None
    rows = []
    for _, r in geo.iterrows():
        oid, role, match = onto(r["condition_id"])
        rows.append({
            "condition_id": r["condition_id"], "ontology_id": oid, "ontology_id_role": role, "ontology_match": match,
            "ontology_link_basis": "registry primary id of the condition; the GEO match itself is a text query",
            "source_disease_id": "", "source_disease_label": None,
            "evidence_type": "expression_study_metadata", "source_name_": "NCBI GEO (Entrez E-utilities db=gds)",
            "source_evidence_category": r["study_type"],
            "entity_type": "study", "entity_id": r["gse"], "entity_label": r["title"],
            "source_database": "NCBI GEO", "source_accession": r["gse"],
            "evidence_direction": None, "study_population": None,
            "sample_size": r["n_samples"], "sample_size_basis": "GEO samples (GSM) in the series; not participants",
            "source_score": np.nan, "source_score_label": None,
            "date_retrieved": retrieved[:10],
            "study_type_in_scope": r["study_type_in_scope"], "mentions_condition": r["mentions_condition"],
            "pubmed_ids": r["pubmed_ids"], "release_date": r["release_date"], "platform": r["platform"],
            "bioproject": r["bioproject"], "query_ids": r["query_ids"],
        })
    for _, r in sra.iterrows():
        oid, role, match = onto(r["condition_id"])
        rows.append({
            "condition_id": r["condition_id"], "ontology_id": oid, "ontology_id_role": role, "ontology_match": match,
            "ontology_link_basis": "registry primary id of the condition; the SRA match itself is a text query",
            "source_disease_id": "", "source_disease_label": None,
            "evidence_type": "metadata_catalog", "source_name_": "NCBI SRA (Entrez E-utilities db=sra)",
            "source_evidence_category": r["library_strategies"],
            "entity_type": "study", "entity_id": r["study_acc"], "entity_label": r["study_title"],
            "source_database": "NCBI SRA", "source_accession": r["study_acc"],
            "evidence_direction": None, "study_population": None,
            "sample_size": r["n_runs_matched"], "sample_size_basis": "SRA runs matched by the query; not participants",
            "source_score": np.nan, "source_score_label": None,
            "date_retrieved": retrieved[:10],
            "mentions_condition": r["mentions_condition"], "bioproject": r["bioproject"],
            "linked_geo_series": r["linked_geo_series"], "query_ids": r["query_ids"],
        })
    return pd.DataFrame(rows)


def coverage_rows(geo: pd.DataFrame, sra: pd.DataFrame, qlog: list[dict], missing: list[str]) -> pd.DataFrame:
    rows = []
    ids = Q.condition_query_ids()
    prim = ids[ids["is_primary"]].set_index("condition_id")
    for cid in Q._registry()["canonical_condition_id"]:
        oid = prim.loc[cid, "ontology_id"] if cid in prim.index else ""
        g = geo[geo["condition_id"] == cid] if len(geo) else geo
        s = sra[sra["condition_id"] == cid] if len(sra) else sra
        ql = [q for q in qlog if q["condition_id"] == cid]
        for db, df, n in (("NCBI GEO", g, int(g["study_type_in_scope"].sum()) if len(g) else 0),
                          ("NCBI SRA", s, len(s))):
            if cid in missing:
                status = "no curated query"
            elif n > 0:
                status = "queried"
            elif db == "NCBI GEO" and len(df):
                status = f"zero in-scope series ({len(df)} series of other types)"
            else:
                status = "zero hits"
            rows.append({
                "condition_id": cid, "ontology_id": oid, "ontology_id_role": "primary (text query, not id-based)",
                "ontology_match": "", "source_database": db, "id_status_in_source": "text query",
                "source_disease_label": None, "coverage_status": status, "n_records": n,
                "detail": json.dumps({
                    "queries": [q["query_id"] for q in ql if (q["db"] == "gds") == (db == "NCBI GEO")],
                    "esearch_counts": {q["query_id"]: q["count"] for q in ql if (q["db"] == "gds") == (db == "NCBI GEO")},
                    **({"series_all_types": len(g), "series_in_scope_types": n,
                        "series_mentioning_condition": int(g["mentions_condition"].sum()) if len(g) else 0}
                       if db == "NCBI GEO" else
                       {"studies": len(s), "bioprojects": int(s["bioproject"].nunique()) if len(s) else 0,
                        "experiments": int(s["n_experiments_matched"].sum()) if len(s) else 0,
                        "runs": int(s["n_runs_matched"].sum()) if len(s) else 0})})})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- precision check
def precision_sample(geo: pd.DataFrame, n: int = PRECISION_N, previous: pd.DataFrame | None = None) -> pd.DataFrame:
    """Seeded sample of up to n in-scope GEO series per demo-cluster condition.

    If a previous (already reviewed) sample exists and all its series are still in the pool, it is
    kept, so a later re-run against a grown GEO does not silently replace reviewed series.
    """
    out = []
    for cid in Q.DEMO_CLUSTER:
        g = geo[(geo["condition_id"] == cid) & (geo["study_type_in_scope"])] if len(geo) else geo
        if g.empty:
            g = geo[geo["condition_id"] == cid] if len(geo) else geo  # fall back to any series type
        g = g.sort_values("gse")
        prev = set(previous.loc[previous["condition_id"] == cid, "gse"]) if previous is not None and len(previous) else set()
        if prev and prev <= set(g["gse"]) and len(prev) == min(n, len(g)):
            take = g[g["gse"].isin(prev)]
        else:
            take = g if len(g) <= n else g.sample(n=n, random_state=SEED)
        out.append(take.assign(sampled_from=len(g)))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def precision_check(geo: pd.DataFrame) -> dict:
    d = raw_dir(SOURCE_ID)
    prev = pd.read_csv(d / "precision_sample.csv") if (d / "precision_sample.csv").exists() else None
    samp = precision_sample(geo, previous=prev)
    if len(samp):
        samp[["condition_id", "gse", "study_type", "n_samples", "title", "summary", "query_ids", "mentions_condition",
              "sampled_from"]].to_csv(d / "precision_sample.csv", index=False)
    review_path = d / "precision_review.yaml"
    review = yaml.safe_load(review_path.read_text()) if review_path.exists() else {}
    verdicts = {(r["condition_id"], r["gse"]): r for r in (review or {}).get("reviews", [])}
    out = {}
    for cid in Q.DEMO_CLUSTER:
        s = samp[samp["condition_id"] == cid] if len(samp) else samp
        rev = [verdicts.get((cid, g)) for g in s["gse"]] if len(s) else []
        done = [v for v in rev if v is not None]
        rel = sum(1 for v in done if v.get("relevant") is True)
        flag = dict(zip(s["gse"], s["mentions_condition"])) if len(s) else {}
        done_m = [v for v in done if flag.get(v["gse"])]
        rel_m = sum(1 for v in done_m if v.get("relevant") is True)
        out[cid] = {"n_series_in_pool": int(s["sampled_from"].iloc[0]) if len(s) else 0,
                    "n_sampled": int(len(s)), "n_reviewed": len(done), "n_relevant": rel,
                    "n_unreviewed": int(len(s)) - len(done),
                    "precision": (round(rel / len(done), 3) if done else None),
                    "auto_mentions_condition": int(s["mentions_condition"].sum()) if len(s) else 0,
                    "n_reviewed_with_auto_flag": len(done_m), "n_relevant_with_auto_flag": rel_m,
                    "precision_with_auto_flag": (round(rel_m / len(done_m), 3) if done_m else None),
                    "false_positives": [v["gse"] for v in done if v.get("relevant") is not True],
                    "reviewer": (review or {}).get("reviewer")}
    (d / "precision_check.json").write_text(json.dumps(out, indent=1))
    return out


# --------------------------------------------------------------------------- audit + registry
def summarize(geo: pd.DataFrame, sra: pd.DataFrame) -> dict:
    per = []
    for cid in Q._registry()["canonical_condition_id"]:
        g = geo[geo["condition_id"] == cid] if len(geo) else geo
        s = sra[sra["condition_id"] == cid] if len(sra) else sra
        per.append({
            "condition_id": cid,
            "geo_series_all_types": len(g),
            "geo_series_in_scope": int(g["study_type_in_scope"].sum()) if len(g) else 0,
            "geo_in_scope_mentioning_condition": int((g["study_type_in_scope"] & g["mentions_condition"]).sum()) if len(g) else 0,
            "geo_samples_in_scope": int(g.loc[g["study_type_in_scope"], "n_samples"].sum()) if len(g) else 0,
            "sra_studies": len(s),
            "sra_bioprojects": int(s["bioproject"].nunique()) if len(s) else 0,
            "sra_studies_mentioning_condition": int(s["mentions_condition"].sum()) if len(s) else 0,
            "sra_runs_matched": int(s["n_runs_matched"].sum()) if len(s) else 0,
            "sra_studies_linked_to_geo": int(s["linked_geo_series"].notna().sum()) if len(s) else 0,
        })
    per = pd.DataFrame(per)
    zero = per[(per["geo_series_in_scope"] == 0) & (per["sra_studies"] == 0)]["condition_id"].tolist()
    return {"per_condition": per, "zero_any": zero,
            "zero_geo": per[per["geo_series_in_scope"] == 0]["condition_id"].tolist()}


def write_audit(geo, sra, ev, cov, qlog, run_log, summary) -> None:
    from .open_targets import _md_table, _missingness
    per = summary["per_condition"]
    prec = pd.DataFrame([{"condition_id": k, **v} for k, v in run_log["precision"].items()])
    ql = pd.DataFrame(qlog)[["condition_id", "query_id", "db", "count", "phrases_not_found", "term"]]
    ql["phrases_not_found"] = ql["phrases_not_found"].map(lambda v: "; ".join(v) if v else "")
    types = geo["study_type"].value_counts().reset_index() if len(geo) else pd.DataFrame()
    dys = geo[(geo["condition_id"] == "dysautonomia") & geo["study_type_in_scope"]] if len(geo) else geo
    dys_n = len(dys)
    dys_fd = int(dys["title"].str.contains("familial dysautonomia|multiple system atrophy|parkinson", case=False,
                                            regex=True).sum()) if len(dys) else 0
    strat = (sra["library_strategies"].str.split(";").explode().value_counts().reset_index()
             if len(sra) else pd.DataFrame())
    L = ["# DATA AUDIT — NCBI GEO / SRA (Entrez E-utilities, metadata only)\n", "| Field | Value |", "|---|---|",
         f"| source_id | {SOURCE_ID} |",
         f"| Source (dataset/API name, exact files/endpoints) | NCBI Entrez E-utilities: {EUTILS}/einfo.fcgi, "
         "esearch.fcgi and esummary.fcgi on db=gds (GEO) and db=sra (SRA); curated terms in "
         "configs/geo_queries.yaml |",
         "| Publishing organization | NCBI / National Library of Medicine |",
         f"| Retrieval date (UTC) | {Q.retrieval_window(SOURCE_ID)} |",
         f"| Source version / release | {run_log['version']} |",
         "| Source update date / cadence | continuous (daily database updates) |",
         "| License / access conditions | open; NCBI E-utilities usage policy (<= 3 requests/s without API key); "
         "individual GEO/SRA submissions carry their submitters' terms; controlled-access (dbGaP) data not touched |",
         "| Unit of observation | GEO series (GSE) x condition; SRA study (SRP/ERP/DRP, experiments aggregated) x condition |",
         f"| Sample size (actual, as ingested) | {len(geo)} condition x GEO-series rows "
         f"({geo['gse'].nunique() if len(geo) else 0} unique GSE; {int(geo['study_type_in_scope'].sum()) if len(geo) else 0} "
         f"in-scope types); {len(sra)} condition x SRA-study rows "
         f"({sra['study_acc'].nunique() if len(sra) else 0} unique studies, "
         f"{int(sra['n_experiments_matched'].sum()) if len(sra) else 0} experiments) |",
         "| Geography (resolution, vintage) | none |",
         "| Person-level? | no (study metadata only; nothing downloaded) |", "| Geographic? | no |",
         "| Omics? | yes — catalogue of public omics studies (condition-level evidence availability) |",
         "| Wearable? | no |",
         "| True participant linkage across modalities? | no — study-level metadata; sample counts are GSM/runs |\n",
         "## Files / endpoints retrieved\n",
         f"Every esearch/esummary/einfo response is saved under `data/raw/{SOURCE_ID}/` (esearch/, esummary/) with "
         "URL, parameters and sha256 in MANIFEST.json. `query_log.json` records every term exactly as sent, its count, "
         "Entrez's query translation and any quoted phrase Entrez reported as not found. `precision_sample.csv`, "
         "`precision_review.yaml` and `precision_check.json` hold the manual precision check.\n",
         "## Queries (term as sent; phrases Entrez did not find are silently dropped by Entrez)\n",
         _md_table(ql),
         "\n## Per-condition counts\n",
         "`geo_series_in_scope` = Homo sapiens GSE whose series type is expression / methylation / non-coding RNA / "
         "protein profiling; `*_mentioning_condition` = the automatic regex QC (title+summary for GEO, study title "
         "for SRA) found a condition term. SRA runs/experiments are those matched by the query, not whole projects.\n",
         _md_table(per),
         f"\nConditions with zero in-scope GEO series AND zero SRA studies: {', '.join(summary['zero_any']) or 'none'}.",
         f"\nConditions with zero in-scope GEO series: {', '.join(summary['zero_geo']) or 'none'}.\n",
         "## Precision check (demo cluster: long_covid, me_cfs, pots)\n",
         "A seeded random sample (seed 20260923) of up to 10 in-scope GEO series per condition "
         "(`precision_sample.csv`) was reviewed by reading each series title and summary; verdicts and reasons are in "
         "`precision_review.yaml`. Criterion: the series studies people with the condition (or their samples) "
         "as a study group, not merely mentions it. The reviewer field records who made the judgement. "
         "`precision_with_auto_flag` restricts to sampled series whose title/summary match the condition QC regex "
         "(the tighter view a downstream user can apply with `mentions_condition`).\n",
         _md_table(prec),
         "\n## GEO series types retrieved\n", _md_table(types) if len(types) else "(none)",
         "\n## SRA library strategies (study level)\n", _md_table(strat) if len(strat) else "(none)",
         "\n## Key variables\n",
         "* geo_study_catalog: `gse`, `title`, `summary`, `platform` (GPL ids), `n_samples` (GSM count), `study_type` "
         "(GEO series type), `pubmed_ids`, `release_date`, `bioproject`, `condition_id`, `query_ids`/`query_terms`, "
         "`study_type_in_scope`, `mentions_condition`.",
         "* sra_study_catalog: `study_acc`, `bioproject`, `study_title`, `library_strategies`, `library_sources`, "
         "`platforms`, `n_experiments_matched`, `n_runs_matched`, `n_biosamples_matched`, `linked_geo_series`, "
         "`mentions_condition`.",
         "* condition_molecular_evidence__geo: one study row per catalogue row (entity_type 'study'); `sample_size` = "
         "GSM count (GEO) or matched runs (SRA), never participants; `ontology_id` = the condition's registry primary "
         "id (the match itself is a text query).\n",
         "## Missingness\n", "GEO catalogue:\n",
         _md_table(_missingness(geo, ["title", "summary", "platform", "n_samples", "study_type", "pubmed_ids",
                                      "release_date", "bioproject"])) if len(geo) else "(none)",
         "\nSRA catalogue:\n",
         _md_table(_missingness(sra, ["study_title", "bioproject", "library_strategies", "platforms",
                                      "first_createdate"])) if len(sra) else "(none)",
         "\n## Linkage strategy\n",
         "Condition -> curated Entrez text query (not an ontology id; GEO/SRA are not consistently ontology-indexed). "
         "GEO series link to SRA through BioProject (`linked_geo_series`). No sample- or participant-level record is "
         "retrieved; nothing here joins to the person layer. The MapME/CFS partition covers GSE251872 etc. separately.\n",
         "## Limitations and caveats\n",
         "* Text queries over [All Fields] also match sample characteristics and related-work mentions; precision "
         "is measured on the demo cluster only and the automatic `mentions_condition` flag is a coarse QC.",
         "* Entrez drops quoted phrases absent from its phrase index without failing; the config avoids them and "
         "the query log lists any that remain.",
         "* GEO sample counts include controls and technical replicates; SRA counts are runs, not people.",
         "* SRA human-organism filter excludes metagenome-labelled microbiome submissions (e.g. 'human gut "
         "metagenome'), which matters for IBS/ME/CFS microbiome work.",
         "* Metadata catalogue only: presence of a study is research-readiness context, not evidence of an effect.",
         f"* dysautonomia is a broad concept: of its {dys_n} in-scope GEO series, {dys_fd} have 'familial dysautonomia', "
         f"'multiple system atrophy' or 'Parkinson' in the title (neurodegenerative/monogenic autonomic failure, not "
         "POTS-like dysautonomia).\n",
         "## Processed outputs\n",
         f"* `{GEO_TABLE}`: {len(geo)} rows", f"* `{SRA_TABLE}`: {len(sra)} rows",
         f"* `{PARTITION}`: {len(ev)} rows", f"* `{COVERAGE}`: {len(cov)} rows",
         "* `condition_molecular_evidence` / `condition_molecular_coverage`: union of all partitions\n",
         "## Reproduce\n", f"`uv run python -m {MODULE}`\n"]
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text("\n".join(L))


def write_registry(geo, sra, ev, run_log, summary) -> None:
    per = summary["per_condition"].set_index("condition_id")
    prec = run_log["precision"]
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "NCBI GEO / SRA (Entrez E-utilities, metadata only)",
        "publisher": "NCBI / National Library of Medicine",
        "landing_url": LANDING_URL,
        "access_urls": [f"{EUTILS}/esearch.fcgi", f"{EUTILS}/esummary.fcgi", f"{EUTILS}/einfo.fcgi"],
        "license": "open NCBI records; submitter terms apply to individual datasets",
        "access_conditions": "open API, no key (<= 3 requests/s)",
        "retrieved_at": run_log["retrieved_at"],
        "source_version": run_log["version"],
        "update_date": f"gds {run_log['einfo']['gds']['lastupdate']}; sra {run_log['einfo']['sra']['lastupdate']}",
        "data_layer": "condition_molecular",
        "unit_of_observation": "GEO series x condition; SRA study x condition (metadata)",
        "sample_size": {
            "geo_condition_series_rows": int(len(geo)),
            "geo_unique_series": int(geo["gse"].nunique()) if len(geo) else 0,
            "sra_condition_study_rows": int(len(sra)),
            "sra_unique_studies": int(sra["study_acc"].nunique()) if len(sra) else 0,
            "per_condition_geo_series_in_scope": {k: int(v) for k, v in per["geo_series_in_scope"].items()},
            "per_condition_sra_studies": {k: int(v) for k, v in per["sra_studies"].items()},
            "precision_demo_cluster": {k: v.get("precision") for k, v in prec.items()},
        },
        "geographic_resolution": "none",
        "person_level": False, "geographic": False, "omics": True, "wearable": False,
        "participant_linkage": "not applicable (study metadata)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",
        "processed_outputs": [GEO_TABLE, SRA_TABLE, PARTITION, COVERAGE, "condition_molecular_evidence",
                              "condition_molecular_coverage"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": MODULE,
        "limitations": [
            "Text-query matches; precision measured only for the demo cluster (see DATA_AUDIT.md)",
            "Metadata only; sample counts are GSM/runs, not participants",
            "SRA human-organism filter excludes metagenome-labelled submissions",
            "Zero in-scope GEO series: " + (", ".join(summary["zero_geo"]) or "none"),
        ],
    })


def main() -> None:  # pragma: no cover
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args()
    print(json.dumps(run(), indent=1, default=str))
    Q.exit_cleanly()


if __name__ == "__main__":  # pragma: no cover
    main()
