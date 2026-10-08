"""Open Targets Platform (GraphQL API v4) -> condition-level molecular evidence.

What this module does
---------------------
For every MONDO/EFO identifier the ontology layer recorded for a registry condition
(``omics.query.condition_query_ids``):

1. ``disease(efoId)`` metadata: Open Targets label, descendants, obsolete terms. An id
   the Platform does not know is recorded as "identifier not present in source".
2. ``associatedTargets`` with the API default (which, measured at run time, is the
   *indirect* association set: descendant diseases are included) and again with
   ``enableIndirect: false`` (direct only). Paginated at the API maximum page size.
   The overall association score and every per-datatype / per-datasource score are
   kept exactly as the API reports them and labelled
   "Open Targets association score (source-provided)".
3. Top evidence by datatype: for each requested datatype (genetic_association,
   known_drug [the 26.x API calls this datatype ``clinical``], literature,
   rna_expression, animal_model, affected_pathway, plus genetic_literature), the
   10 targets with the highest datatype score and their evidence items
   (``disease.evidences``, indirect, filtered to that datatype's datasources).
4. Known drugs / clinical candidates: ``disease.drugAndClinicalCandidates`` (the 26.x
   replacement of ``knownDrugs``) with drug, mechanism of action, maximum clinical
   stage for the indication and the trial statuses of the linked clinical reports.
5. Reactome pathways of the associated targets as the API returns them
   (``Target.pathways``), aggregated to one row per condition id x pathway with
   the number of associated targets annotated to it. No enrichment statistic is
   computed; these are annotation counts.

Outputs
-------
* ``condition_molecular_evidence__open_targets`` (data_layer condition_molecular)
* ``condition_molecular_coverage__open_targets`` (per condition x id: present / absent / counts)
* raw GraphQL responses under ``data/raw/open_targets/`` + MANIFEST.json, DATA_AUDIT.md,
  registry_entry.yaml

Guardrail: condition-level molecular enrichment/evidence only. No combined omics score.

Reproduce: ``uv run python -m measure_it.omics.open_targets``
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from ..config import raw_dir
from ..http import post
from ..provenance import add_provenance
from ..registry import write_registry_entry
from ..store import write_table
from . import query as Q

SOURCE_ID = "open_targets"
SOURCE_DB = "Open Targets Platform"
API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
LANDING_URL = "https://platform.opentargets.org/"
DOCS_URL = "https://platform-docs.opentargets.org/data-access/graphql-api"
MODULE = "measure_it.omics.open_targets"
PARTITION = "condition_molecular_evidence__open_targets"
COVERAGE = "condition_molecular_coverage__open_targets"
SCORE_LABEL = "Open Targets association score (source-provided)"
EVIDENCE_SCORE_LABEL = "Open Targets evidence score (source-provided)"

PAGE_SIZE = 3000            # API maximum (a size of 5000 is rejected: "size must be between 0 and 3000")
TOP_TARGETS_PER_DATATYPE = 10
EVIDENCE_PAGE_SIZE = 200
MAX_EVIDENCE_PAGES = 3
TOP_EVIDENCE_PER_TARGET = 5
PATHWAY_BATCH = 500
TOP_SYMBOLS_PER_PATHWAY = 25

# requested label -> datatype ids the API may use for it (first one observed wins)
REQUESTED_DATATYPES = {
    "genetic_association": ["genetic_association"],
    "genetic_literature": ["genetic_literature"],
    "known_drug": ["known_drug", "clinical"],
    "literature": ["literature"],
    "rna_expression": ["rna_expression"],
    "animal_model": ["animal_model"],
    "affected_pathway": ["affected_pathway"],
}

# OT datatype -> provenance evidence_type vocabulary
DATATYPE_TO_EVIDENCE_TYPE = {
    "genetic_association": "genetic_association",
    "genetic_literature": "genetic_association",
    "known_drug": "known_drug",
    "clinical": "known_drug",
    "literature": "text_mined",
    "affected_pathway": "pathway_annotation",
}
DEFAULT_EVIDENCE_TYPE = "target_disease_association"
# Association rows (one per target) are identified by this category; evidence items of datatypes without
# a more specific evidence_type (animal_model, rna_expression, somatic_mutation) share evidence_type
# 'target_disease_association' but carry their Open Targets datatype as source_evidence_category.
ASSOC_CATEGORY = "overall_association"

# --------------------------------------------------------------------------- GraphQL documents
META_Q = """{ meta { name product apiVersion { x y z suffix } dataVersion { year month iteration } dataPrefix } }"""

DISEASE_Q = """query disease($id: String!) {
  disease(efoId: $id) {
    id name description dbXRefs obsoleteTerms descendants
    therapeuticAreas { id name }
    parents { id name }
  }
}"""

_ASSOC_BODY = """count
      rows {
        score
        target { id approvedSymbol approvedName biotype }
        datatypeScores { id score }
        datasourceScores { id score }
      }"""
ASSOC_DEFAULT_Q = ("query assoc($id: String!, $index: Int!, $size: Int!) { disease(efoId: $id) { id "
                   "associatedTargets(page: {index: $index, size: $size}) { " + _ASSOC_BODY + " } } }")
ASSOC_DIRECT_Q = ("query assocDirect($id: String!, $index: Int!, $size: Int!) { disease(efoId: $id) { id "
                  "associatedTargets(page: {index: $index, size: $size}, enableIndirect: false) { "
                  + _ASSOC_BODY + " } } }")
# Multi-page variants. Measured 2026-09-23 (data release 26.06): paging the default (score-descending)
# order across the 3000-row page boundary repeats some targets and skips others, reproducibly
# (endometriosis: 3561 rows, 3521 unique; autonomic nervous system disorder: 3246 rows, 3225 unique;
# post-infectious disorder: 8240 rows, 8239 unique), while orderByScore "score asc" returned every
# target exactly once for all of them. Ids with more than PAGE_SIZE associations are therefore paged
# in ascending order and re-sorted by score client-side; completeness is asserted.
ASSOC_DEFAULT_ASC_Q = ("query assocAsc($id: String!, $index: Int!, $size: Int!) { disease(efoId: $id) { id "
                       "associatedTargets(page: {index: $index, size: $size}, orderByScore: \"score asc\") { "
                       + _ASSOC_BODY + " } } }")
ASSOC_DIRECT_ASC_Q = ("query assocDirectAsc($id: String!, $index: Int!, $size: Int!) { disease(efoId: $id) { id "
                      "associatedTargets(page: {index: $index, size: $size}, enableIndirect: false, "
                      "orderByScore: \"score asc\") { " + _ASSOC_BODY + " } } }")
PAGING_LOG: dict[str, dict] = {}


class IncompleteAssociations(RuntimeError):
    pass
ASSOC_INDIRECT_COUNT_Q = ("query assocInd($id: String!) { disease(efoId: $id) { id "
                          "associatedTargets(page: {index: 0, size: 0}, enableIndirect: true) { count } } }")

EVIDENCE_Q = """query ev($id: String!, $ens: [String!]!, $ds: [String!], $size: Int, $cursor: String) {
  disease(efoId: $id) {
    evidences(ensemblIds: $ens, enableIndirect: true, datasourceIds: $ds, size: $size, cursor: $cursor) {
      count cursor
      rows {
        id datasourceId datatypeId score resourceScore
        target { id approvedSymbol }
        disease { id name }
        diseaseFromSource diseaseFromSourceMappedId
        studyId studySampleSize studyCases cohortShortName
        variantRsId variant { id rsIds }
        credibleSet {
          studyLocusId pValueMantissa pValueExponent beta
          variant { id rsIds }
          study { id traitFromSource nSamples nCases nControls initialSampleSize pubmedId
                  publicationFirstAuthor discoverySamples { ancestry sampleSize } }
        }
        literature publicationYear
        pValueMantissa pValueExponent oddsRatio beta confidence
        drug { id name } clinicalStage
        directionOnTrait directionOnTarget
        reactionId reactionName pathways { id name }
        contrast log2FoldChangeValue log2FoldChangePercentileRank
        biologicalModelId biologicalModelAllelicComposition
        diseaseModelAssociatedModelPhenotypes { id label }
        releaseVersion
      }
    }
  }
}"""

DATASOURCE_PROBE_Q = """query probe($id: String!, $ens: [String!]!, $ds: [String!]) {
  disease(efoId: $id) { evidences(ensemblIds: $ens, enableIndirect: true, datasourceIds: $ds, size: 1) {
    count rows { datasourceId datatypeId } } }
}"""

DRUGS_Q = """query drugs($id: String!) {
  disease(efoId: $id) {
    id
    drugAndClinicalCandidates {
      count
      rows {
        id maxClinicalStage
        drug {
          id name drugType maximumClinicalStage
          mechanismsOfAction { rows { mechanismOfAction actionType targetName targets { id approvedSymbol } } }
        }
        clinicalReports {
          id source type clinicalStage trialPhase trialOverallStatus year url
          diseases { diseaseFromSource disease { id } }
        }
      }
    }
  }
}"""

PATHWAYS_Q = """query tp($ids: [String!]!) {
  targets(ensemblIds: $ids) { id approvedSymbol pathways { pathwayId pathway topLevelTerm } }
}"""


# --------------------------------------------------------------------------- transport
class GraphQLError(RuntimeError):
    pass


def gql(query: str, variables: dict | None = None, *, raw_rel: str | None = None) -> tuple[dict, str]:
    """POST a GraphQL document through the shared cached HTTP layer. Returns (data, fetched_at)."""
    body = {"query": query, "variables": variables or {}}
    r = post(API_URL, json_body=body, expect_json=True, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    payload = r.json()
    if payload.get("errors"):
        raise GraphQLError(f"{variables}: {payload['errors'][0].get('message', '')[:300]}")
    if raw_rel:
        Q.save_raw_json(SOURCE_ID, raw_rel, payload, url=API_URL,
                        request={"method": "POST", "graphql_variables": variables,
                                 "graphql_operation": query.split("{", 1)[0].strip()[:80]},
                        retrieved_at=r.fetched_at)
    return payload["data"], r.fetched_at


# --------------------------------------------------------------------------- fetchers
def fetch_meta() -> tuple[dict, str]:
    data, ts = gql(META_Q, raw_rel="meta.json")
    return data["meta"], ts


def release_string(meta: dict) -> str:
    dv, av = meta["dataVersion"], meta["apiVersion"]
    it = f".{dv['iteration']}" if dv.get("iteration") else ""
    return (f"Open Targets Platform data release {dv['year']}.{dv['month']}{it}; "
            f"API {av['x']}.{av['y']}.{av['z']}{('-' + av['suffix']) if av.get('suffix') else ''}")


def fetch_disease(short_id: str) -> dict | None:
    data, _ = gql(DISEASE_Q, {"id": short_id}, raw_rel=f"diseases/{Q.safe_name(short_id)}.json")
    return data["disease"]


def _page_all(q: str, short_id: str, tag: str) -> tuple[int, list[dict]]:
    rows, index = [], 0
    while True:
        data, _ = gql(q, {"id": short_id, "index": index, "size": PAGE_SIZE},
                      raw_rel=f"associations/{Q.safe_name(short_id)}__{tag}__page{index}.json")
        at = data["disease"]["associatedTargets"]
        count = at["count"]
        rows.extend(at["rows"])
        index += 1
        if len(rows) >= count or not at["rows"]:
            return count, rows


def fetch_associations(short_id: str, *, direct: bool) -> tuple[int, list[dict]]:
    """All associated targets of one disease id, sorted by overall score (descending)."""
    tag = "direct" if direct else "default"
    count, rows = _page_all(ASSOC_DIRECT_Q if direct else ASSOC_DEFAULT_Q, short_id, tag)
    strategy = "single page (API default order)"
    if count > PAGE_SIZE:
        count, rows = _page_all(ASSOC_DIRECT_ASC_Q if direct else ASSOC_DEFAULT_ASC_Q, short_id, f"{tag}_asc")
        strategy = "paged with orderByScore 'score asc', re-sorted by score client-side"
    uniq = {r["target"]["id"] for r in rows}
    PAGING_LOG[f"{short_id}|{tag}"] = {"count": count, "rows": len(rows), "unique_targets": len(uniq),
                                       "strategy": strategy}
    if len(uniq) != count or len(rows) != count:
        raise IncompleteAssociations(f"{short_id} {tag}: API count {count}, rows {len(rows)}, unique {len(uniq)}")
    rows = sorted(rows, key=lambda r: -r["score"])
    return count, rows


def fetch_indirect_true_count(short_id: str) -> int:
    data, _ = gql(ASSOC_INDIRECT_COUNT_Q, {"id": short_id},
                  raw_rel=f"associations/{Q.safe_name(short_id)}__indirect_true_count.json")
    return data["disease"]["associatedTargets"]["count"]


def fetch_evidence(short_id: str, target_ids: list[str], datasource_ids: list[str], tag: str) -> tuple[int, list[dict]]:
    rows, cursor, total = [], None, 0
    for page in range(MAX_EVIDENCE_PAGES):
        data, _ = gql(EVIDENCE_Q, {"id": short_id, "ens": target_ids, "ds": datasource_ids,
                                   "size": EVIDENCE_PAGE_SIZE, "cursor": cursor},
                      raw_rel=f"evidence/{Q.safe_name(short_id)}__{tag}__page{page}.json")
        ev = data["disease"]["evidences"]
        total = ev["count"]
        rows.extend(ev["rows"])
        cursor = ev.get("cursor")
        if not cursor or len(rows) >= total or not ev["rows"]:
            break
    return total, rows


def fetch_drugs(short_id: str) -> tuple[int, list[dict]]:
    data, _ = gql(DRUGS_Q, {"id": short_id}, raw_rel=f"drugs/{Q.safe_name(short_id)}.json")
    d = data["disease"]["drugAndClinicalCandidates"]
    return d["count"], d["rows"]


def fetch_target_pathways(target_ids: list[str]) -> dict[str, dict]:
    out = {}
    ids = sorted(set(target_ids))
    for i in range(0, len(ids), PATHWAY_BATCH):
        chunk = ids[i:i + PATHWAY_BATCH]
        data, _ = gql(PATHWAYS_Q, {"ids": chunk}, raw_rel=f"target_pathways/batch_{i // PATHWAY_BATCH:04d}.json")
        for t in data["targets"]:
            out[t["id"]] = t
    return out


def datasource_datatype_map(examples: dict[str, tuple[str, str]]) -> dict[str, str]:
    """Measure which datatype each datasource belongs to by fetching one evidence item for it.

    examples: datasource_id -> (disease short id, target id) where that datasource scored.
    """
    out = {}
    for ds, (sid, tid) in sorted(examples.items()):
        data, _ = gql(DATASOURCE_PROBE_Q, {"id": sid, "ens": [tid], "ds": [ds]},
                      raw_rel=f"datasource_probe/{Q.safe_name(ds)}.json")
        rows = data["disease"]["evidences"]["rows"]
        if rows:
            out[ds] = rows[0]["datatypeId"]
    return out


# --------------------------------------------------------------------------- row builders (pure)
def association_rows(ctx: dict, assoc_rows: list[dict], direct_scores: dict[str, float], count: int) -> list[dict]:
    """One row per associated target with overall + per-datatype scores exactly as reported."""
    out = []
    for rank, r in enumerate(assoc_rows, start=1):
        t = r["target"]
        row = {**ctx,
               "evidence_type": "target_disease_association",
               "source_evidence_category": ASSOC_CATEGORY,
               "entity_type": "gene", "entity_id": t["id"], "entity_label": t["approvedSymbol"],
               "target_name": t.get("approvedName"), "target_biotype": t.get("biotype"),
               "source_database": SOURCE_DB,
               "source_accession": f"{ctx['source_disease_id']}--{t['id']}",
               "source_score": r["score"], "source_score_label": SCORE_LABEL,
               "association_mode": "API default (indirect: includes descendant diseases)",
               "association_rank": rank,
               "n_associated_targets_for_id": count,
               "direct_association_score": direct_scores.get(t["id"], np.nan),
               "has_direct_association": t["id"] in direct_scores,
               "ot_datasource_scores_json": json.dumps({x["id"]: x["score"] for x in r["datasourceScores"]},
                                                       sort_keys=True)}
        for x in r["datatypeScores"]:
            row[f"ot_datatype_score__{x['id']}"] = x["score"]
        out.append(row)
    return out


def _pvalue(m, e):
    if m is None or e is None:
        return np.nan
    try:
        return float(m) * 10.0 ** int(e)
    except (TypeError, ValueError):
        return np.nan


def evidence_row(ctx: dict, e: dict, requested: str) -> dict:
    """Flatten one Open Targets evidence item (source values only)."""
    dt = e["datatypeId"]
    cs = e.get("credibleSet") or {}
    st = cs.get("study") or {}
    variant = e.get("variant") or cs.get("variant") or {}
    rsids = variant.get("rsIds") or ([e["variantRsId"]] if e.get("variantRsId") else [])
    study_id = e.get("studyId") or st.get("id")
    lit = e.get("literature") or []
    drug = e.get("drug") or {}
    pws = e.get("pathways") or []
    sample_size = st.get("nSamples") if st.get("nSamples") is not None else e.get("studySampleSize")
    if st.get("discoverySamples"):
        pop = "; ".join(f"{s['ancestry']} (n={s['sampleSize']})" for s in st["discoverySamples"])
        if st.get("initialSampleSize"):
            pop = f"{pop} | initial sample: {st['initialSampleSize']}"
    elif e.get("cohortShortName"):
        pop = e["cohortShortName"]
    else:
        pop = None
    if dt in ("clinical", "known_drug"):
        acc = drug.get("id")
    elif dt == "affected_pathway":
        acc = e.get("reactionId") or (pws[0]["id"] if pws else None)
    elif dt == "literature":
        acc = ";".join(lit[:10]) if lit else None
    elif dt == "animal_model":
        acc = e.get("biologicalModelId")
    else:
        acc = study_id or (";".join(lit[:10]) if lit else None)
    dirs = []
    if e.get("directionOnTarget"):
        dirs.append(f"target:{e['directionOnTarget']}")
    if e.get("directionOnTrait"):
        dirs.append(f"trait:{e['directionOnTrait']}")
    return {
        **ctx,
        "evidence_type": DATATYPE_TO_EVIDENCE_TYPE.get(dt, DEFAULT_EVIDENCE_TYPE),
        "source_evidence_category": dt,
        "requested_datatype": requested,
        "ot_datasource": e["datasourceId"],
        "entity_type": "gene", "entity_id": e["target"]["id"], "entity_label": e["target"]["approvedSymbol"],
        "source_database": SOURCE_DB,
        "source_accession": acc,
        "ot_evidence_id": e["id"],
        "evidence_disease_id": e["disease"]["id"], "evidence_disease_label": e["disease"]["name"],
        "evidence_from_descendant": e["disease"]["id"] != ctx["source_disease_id"],
        "disease_from_source": e.get("diseaseFromSource"),
        "evidence_direction": ", ".join(dirs) if dirs else None,
        "study_population": pop,
        "sample_size": sample_size,
        "sample_size_basis": ("GWAS study nSamples (Open Targets study index)" if st.get("nSamples") is not None
                              else ("studySampleSize as reported by the datasource" if e.get("studySampleSize") else None)),
        "source_score": e["score"], "source_score_label": EVIDENCE_SCORE_LABEL,
        "resource_score": e.get("resourceScore"),
        "study_id": study_id,
        "variant_id": variant.get("id"),
        "variant_rsids": ";".join(rsids) if rsids else None,
        "p_value": (_pvalue(e.get("pValueMantissa"), e.get("pValueExponent"))
                    if e.get("pValueMantissa") is not None else _pvalue(cs.get("pValueMantissa"), cs.get("pValueExponent"))),
        "odds_ratio": e.get("oddsRatio"),
        "beta": e.get("beta") if e.get("beta") is not None else cs.get("beta"),
        "gwas_trait_from_source": st.get("traitFromSource"),
        "literature_pmids": ";".join(lit) if lit else None,
        "publication_year": e.get("publicationYear"),
        "drug_id": drug.get("id"), "drug_name": drug.get("name"),
        "clinical_stage": e.get("clinicalStage"),
        "reactome_id": e.get("reactionId"), "reactome_name": e.get("reactionName"),
        "expression_contrast": e.get("contrast"), "log2_fold_change": e.get("log2FoldChangeValue"),
        "animal_model_id": e.get("biologicalModelId"),
        "animal_model_phenotypes": ";".join(p["label"] for p in (e.get("diseaseModelAssociatedModelPhenotypes") or [])[:10]) or None,
        "ot_release_version": e.get("releaseVersion"),
    }


def select_top_evidence(rows: list[dict], per_target: int = TOP_EVIDENCE_PER_TARGET) -> list[dict]:
    """Keep the highest-scoring evidence items per target (rows arrive sorted by score, descending)."""
    seen: Counter = Counter()
    out = []
    for r in sorted(rows, key=lambda x: -(x.get("score") or 0)):
        tid = r["target"]["id"]
        if seen[tid] < per_target:
            out.append(r)
            seen[tid] += 1
    return out


def drug_rows(ctx: dict, drugs: list[dict]) -> list[dict]:
    out = []
    for d in drugs:
        drug = d.get("drug") or {}
        moa = ((drug.get("mechanismsOfAction") or {}).get("rows")) or []
        mech = sorted({m["mechanismOfAction"] for m in moa})
        acts = sorted({m["actionType"] for m in moa if m.get("actionType")})
        tsym = sorted({t["approvedSymbol"] for m in moa for t in (m.get("targets") or [])})
        tids = sorted({t["id"] for m in moa for t in (m.get("targets") or [])})
        reps = d.get("clinicalReports") or []
        rep_dis = {x["disease"]["id"] for r in reps for x in (r.get("diseases") or []) if x.get("disease")}
        statuses = Counter(r.get("trialOverallStatus") or "not reported" for r in reps)
        out.append({
            **ctx,
            "evidence_type": "known_drug",
            "source_evidence_category": "drugAndClinicalCandidates",
            "entity_type": "drug", "entity_id": drug.get("id"), "entity_label": drug.get("name"),
            "source_database": SOURCE_DB, "source_accession": drug.get("id"),
            "ot_indication_id": d["id"],
            "drug_type": drug.get("drugType"),
            "max_clinical_stage_for_indication": d.get("maxClinicalStage"),
            "drug_maximum_clinical_stage_any_indication": drug.get("maximumClinicalStage"),
            "mechanism_of_action": "; ".join(mech) or None,
            "action_types": "; ".join(acts) or None,
            "mechanism_target_symbols": ";".join(tsym) or None,
            "mechanism_target_ids": ";".join(tids) or None,
            "n_clinical_reports": len(reps),
            "clinical_report_ids": ";".join(sorted(r["id"] for r in reps)[:50]) or None,
            "clinical_report_sources": ";".join(sorted({r.get("source") or "" for r in reps})) or None,
            "clinical_report_types": ";".join(sorted({r.get("type") or "" for r in reps})) or None,
            "trial_status_counts": json.dumps(dict(statuses), sort_keys=True),
            "report_disease_ids": ";".join(sorted(rep_dis)) or None,
            "reports_cite_queried_id": ctx["source_disease_id"] in rep_dis,
            "source_score": np.nan, "source_score_label": None,
        })
    return out


def pathway_rows(ctx: dict, assoc_rows: list[dict], target_pathways: dict[str, dict]) -> list[dict]:
    """Aggregate Reactome annotations of the associated targets to condition-id x pathway counts."""
    by_pw: dict[str, dict] = {}
    members: dict[str, list[tuple[float, str]]] = defaultdict(list)
    genetic: Counter = Counter()
    for r in assoc_rows:
        tid = r["target"]["id"]
        has_gen = any(x["id"] == "genetic_association" and x["score"] > 0 for x in r["datatypeScores"])
        for pw in (target_pathways.get(tid) or {}).get("pathways") or []:
            pid = pw["pathwayId"]
            by_pw.setdefault(pid, pw)
            members[pid].append((r["score"], r["target"]["approvedSymbol"]))
            if has_gen:
                genetic[pid] += 1
    out = []
    for pid, pw in by_pw.items():
        mem = sorted(members[pid], key=lambda x: -x[0])
        out.append({
            **ctx,
            "evidence_type": "pathway_annotation",
            "source_evidence_category": "Target.pathways (Reactome) of associated targets",
            "entity_type": "pathway", "entity_id": pid, "entity_label": pw["pathway"],
            "reactome_top_level_term": pw.get("topLevelTerm"),
            "source_database": f"{SOURCE_DB} (Reactome annotation)", "source_accession": pid,
            "n_associated_targets_annotated": len(mem),
            "n_genetic_association_targets_annotated": genetic[pid],
            "top_target_symbols": ";".join(s for _, s in mem[:TOP_SYMBOLS_PER_PATHWAY]),
            "source_score": np.nan, "source_score_label": None,
        })
    return out


# --------------------------------------------------------------------------- pipeline
def run(*, limit_ids: int | None = None) -> dict:
    raw_dir(SOURCE_ID)
    meta, meta_ts = fetch_meta()
    version = release_string(meta)
    retrieved = meta_ts
    ids = Q.condition_query_ids(include_obsolete_efo=False)
    if limit_ids:
        ids = ids.head(limit_ids)

    # 1. resolve every id in the Platform
    diseases: dict[str, dict | None] = {}
    for sid in sorted(set(ids["source_query_id"]) - {""}):
        diseases[sid] = fetch_disease(sid)

    # 2. associations (default = indirect; direct separately) for present ids
    assoc: dict[str, tuple[int, list]] = {}
    direct: dict[str, tuple[int, list]] = {}
    indirect_true_count: dict[str, int] = {}
    for sid, d in diseases.items():
        if d is None:
            continue
        assoc[sid] = fetch_associations(sid, direct=False)
        direct[sid] = fetch_associations(sid, direct=True)
        indirect_true_count[sid] = fetch_indirect_true_count(sid)

    # datasource -> datatype, measured from one evidence item per datasource
    examples = {}
    for sid, (_, rows) in assoc.items():
        for r in rows:
            for x in r["datasourceScores"]:
                examples.setdefault(x["id"], (sid, r["target"]["id"]))
    ds_map = datasource_datatype_map(examples)
    datatype_sources: dict[str, list[str]] = defaultdict(list)
    for ds, dt in ds_map.items():
        datatype_sources[dt].append(ds)
    observed_dt = set(datatype_sources)
    requested_to_api = {}
    for req, cands in REQUESTED_DATATYPES.items():
        hit = next((c for c in cands if c in observed_dt), None)
        requested_to_api[req] = hit

    all_rows: list[dict] = []
    coverage: list[dict] = []
    evidence_log: list[dict] = []
    drug_log: dict[str, dict] = {}
    all_targets: set[str] = set()
    for sid, (cnt, rows) in assoc.items():
        all_targets.update(r["target"]["id"] for r in rows)
    target_pw = fetch_target_pathways(sorted(all_targets)) if all_targets else {}

    for _, idr in ids.iterrows():
        sid = idr["source_query_id"]
        base_cov = {"condition_id": idr["condition_id"], "ontology_id": idr["ontology_id"],
                    "ontology_id_role": idr["ontology_id_role"], "ontology_match": idr["ontology_match"],
                    "source_database": SOURCE_DB}
        if not sid:
            coverage.append({**base_cov, "id_status_in_source": "no identifier", "source_disease_label": None,
                             "coverage_status": "no identifier", "n_records": 0,
                             "detail": "condition has no MONDO/EFO id in the registry"})
            continue
        d = diseases.get(sid)
        if d is None:
            coverage.append({**base_cov, "id_status_in_source": "not present in source", "source_disease_label": None,
                             "coverage_status": "identifier not present in source", "n_records": 0,
                             "detail": f"disease(efoId: {sid}) returned null in {version}"})
            continue
        ctx = {"condition_id": idr["condition_id"], "ontology_id": idr["ontology_id"],
               "ontology_id_role": idr["ontology_id_role"], "ontology_match": idr["ontology_match"],
               "source_disease_id": d["id"], "source_disease_label": d["name"],
               "registry_label": idr["registry_label"],
               "label_differs_from_registry": str(d["name"]).strip().lower() != str(idr["registry_label"]).strip().lower(),
               "n_descendants_in_source": len(d["descendants"]),
               "date_retrieved": retrieved[:10]}
        cnt, rows = assoc[sid]
        dcnt, drows = direct[sid]
        dscores = {r["target"]["id"]: r["score"] for r in drows}
        a_rows = association_rows(ctx, rows, dscores, cnt)
        all_rows.extend(a_rows)

        # top evidence by datatype
        n_ev = 0
        for req, api_dt in requested_to_api.items():
            if api_dt is None:
                evidence_log.append({"source_disease_id": sid, "requested": req, "api_datatype": None,
                                     "status": "datatype not observed in the associations of any queried id", "n_targets": 0,
                                     "evidence_total": 0, "evidence_retrieved": 0, "evidence_kept": 0})
                continue
            scored = [(r[f"ot_datatype_score__{api_dt}"], r["entity_id"]) for r in a_rows
                      if r.get(f"ot_datatype_score__{api_dt}", 0) and r.get(f"ot_datatype_score__{api_dt}") > 0]
            scored.sort(key=lambda x: -x[0])
            top = [t for _, t in scored[:TOP_TARGETS_PER_DATATYPE]]
            if not top:
                evidence_log.append({"source_disease_id": sid, "requested": req, "api_datatype": api_dt,
                                     "status": "no target with this datatype", "n_targets": 0,
                                     "evidence_total": 0, "evidence_retrieved": 0, "evidence_kept": 0})
                continue
            total, ev = fetch_evidence(sid, top, sorted(datatype_sources[api_dt]), tag=req)
            kept = select_top_evidence(ev)
            all_rows.extend(evidence_row(ctx, e, req) for e in kept)
            n_ev += len(kept)
            evidence_log.append({"source_disease_id": sid, "requested": req, "api_datatype": api_dt,
                                 "status": "ok", "n_targets": len(top), "n_targets_with_datatype": len(scored),
                                 "evidence_total": total, "evidence_retrieved": len(ev), "evidence_kept": len(kept)})

        # known drugs / clinical candidates
        dc, drugs = fetch_drugs(sid)
        d_rows = drug_rows(ctx, drugs)
        all_rows.extend(d_rows)
        drug_log[sid] = {"count": dc, "rows": len(drugs),
                         "rows_citing_queried_id": sum(r["reports_cite_queried_id"] for r in d_rows)}

        # pathways of associated targets
        p_rows = pathway_rows(ctx, rows, target_pw)
        all_rows.extend(p_rows)

        n_gen = sum(1 for r in a_rows if (r.get("ot_datatype_score__genetic_association") or 0) > 0)
        coverage.append({**base_cov, "id_status_in_source": "present", "source_disease_label": d["name"],
                         "coverage_status": "queried" if cnt else "zero hits", "n_records": cnt,
                         "detail": json.dumps({
                             "n_associated_targets_default": cnt, "n_associated_targets_direct": dcnt,
                             "n_associated_targets_enableIndirect_true": indirect_true_count[sid],
                             "n_genetic_association_targets": n_gen,
                             "n_known_drug_rows": len(d_rows), "n_pathways": len(p_rows),
                             "n_top_evidence_rows": n_ev, "n_descendants": len(d["descendants"])})})

    # ---- frames + provenance
    df = pd.DataFrame(all_rows)
    notes = ("Condition-level molecular enrichment/evidence from Open Targets; not participant-linked. "
             "source_score is the Open Targets score as published (associations: overall association score; "
             "evidence rows: evidence score); no combined omics score is computed. Associations use the API "
             "default, which includes descendant diseases (indirect).")
    frames = []
    for et, g in df.groupby("evidence_type"):
        frames.append(add_provenance(
            g.drop(columns=["evidence_type"]), data_layer="condition_molecular", source_name=f"{SOURCE_DB} GraphQL API v4",
            source_version=version, retrieved_at=retrieved, evidence_type=et,
            source_record_id=lambda x: x["source_disease_id"].astype(str) + "|" + x["entity_type"].astype(str) + "|"
            + x["entity_id"].astype(str) + "|" + x["source_evidence_category"].astype(str) + "|"
            + x.get("ot_evidence_id", pd.Series("", index=x.index)).fillna("").astype(str),
            evidence_level=g["source_evidence_category"].map(
                lambda c: "source_association_score" if c == ASSOC_CATEGORY else "source_record"),
            provenance_notes=notes))
    out = Q.stringify_object_columns(Q.finalize_columns(pd.concat(frames, ignore_index=True)))
    write_table(out, PARTITION, producer=MODULE,
                description="Open Targets condition-level molecular evidence (associations, top evidence, known drugs, "
                            "Reactome pathways of associated targets)")

    cov = pd.DataFrame(coverage)
    cov = _condition_level_coverage(cov)
    cov = add_provenance(cov, data_layer="condition_molecular", source_name=f"{SOURCE_DB} GraphQL API v4",
                         source_version=version, retrieved_at=retrieved, evidence_type="metadata_catalog",
                         source_record_id=lambda x: x["condition_id"] + "|" + x["ontology_id"].astype(str),
                         provenance_notes="coverage of registry identifiers in Open Targets")
    write_table(Q.stringify_object_columns(cov), COVERAGE, producer=MODULE,
                description="Open Targets coverage per condition x registry identifier")

    run_log = {"meta": meta, "version": version, "retrieved_at": retrieved,
               "default_includes_descendants": {sid: {"default": assoc[sid][0], "enableIndirect_true": indirect_true_count[sid],
                                                       "direct": direct[sid][0]} for sid in assoc},
               "association_paging": PAGING_LOG,
               "datasource_to_datatype": ds_map, "requested_to_api_datatype": requested_to_api,
               "evidence_log": evidence_log, "drug_log": drug_log,
               "n_unique_targets_pathway_lookup": len(all_targets),
               "n_targets_with_pathways": sum(1 for t in target_pw.values() if t.get("pathways"))}
    (raw_dir(SOURCE_ID) / "run_log.json").write_text(json.dumps(run_log, indent=1, default=str))
    summary = summarize(out, cov, run_log)
    write_audit(out, cov, run_log, summary)
    write_registry(out, cov, run_log, summary)
    union = Q.build_union()
    return {"rows": len(out), "coverage_rows": len(cov), "version": version, "union": union, **summary["totals"]}


def _condition_level_coverage(cov: pd.DataFrame) -> pd.DataFrame:
    """Add one condition-level row per condition: 'no identifier usable by source' when no id resolved."""
    extra = []
    for cid, g in cov.groupby("condition_id"):
        present = g[g["id_status_in_source"] == "present"]
        if present.empty:
            extra.append({"condition_id": cid, "ontology_id": "(all registry ids)", "ontology_id_role": "condition",
                          "ontology_match": "", "source_database": SOURCE_DB,
                          "id_status_in_source": "no identifier usable by source", "source_disease_label": None,
                          "coverage_status": "no identifier", "n_records": 0,
                          "detail": "none of the registry MONDO/EFO ids exists in the Open Targets disease index: "
                                    + ";".join(g["ontology_id"].astype(str))})
    return pd.concat([cov, pd.DataFrame(extra)], ignore_index=True) if extra else cov


def summarize(df: pd.DataFrame, cov: pd.DataFrame, run_log: dict) -> dict:
    per = []
    reg = Q._registry()
    for cid in reg["canonical_condition_id"]:
        d = df[df["condition_id"] == cid]
        a = d[d["source_evidence_category"] == ASSOC_CATEGORY]
        gen = a[a.get("ot_datatype_score__genetic_association", pd.Series(dtype=float)).fillna(0) > 0] if len(a) else a
        a_ex = a[a["ontology_match"].isin(["exact", "narrow"])]
        a_pr = a[a["ontology_id_role"] == "primary"]
        gen_pr = gen[gen["ontology_id_role"] == "primary"] if len(gen) else gen
        per.append({
            "condition_id": cid,
            "ids_present": ";".join(sorted(set(d["ontology_id"]))) or "",
            "unique_targets": int(a["entity_id"].nunique()),
            "unique_targets_exact_or_narrow_ids": int(a_ex["entity_id"].nunique()),
            "unique_targets_primary_id": int(a_pr["entity_id"].nunique()),
            "unique_genetic_association_targets": int(gen["entity_id"].nunique()),
            "genetic_association_targets_primary_id": int(gen_pr["entity_id"].nunique()),
            "direct_targets": int(a[a["has_direct_association"].fillna(False).astype(bool)]["entity_id"].nunique()) if len(a) else 0,
            "known_drugs": int(d[d["entity_type"] == "drug"]["entity_id"].nunique()),
            "pathways": int(d[d["entity_type"] == "pathway"]["entity_id"].nunique()),
            "top_evidence_rows": int(len(d[(d["entity_type"] == "gene") & (d["source_evidence_category"] != ASSOC_CATEGORY)])),
        })
    per = pd.DataFrame(per)
    zero = per[per["unique_targets"] == 0]["condition_id"].tolist()
    return {"per_condition": per, "zero_evidence": zero,
            "totals": {"unique_targets": int(df[df["source_evidence_category"] == ASSOC_CATEGORY]["entity_id"].nunique()),
                       "association_rows": int((df["source_evidence_category"] == ASSOC_CATEGORY).sum()),
                       "conditions_zero_targets": zero}}


def descendant_evidence_drivers(df: pd.DataFrame, top_k: int = 5) -> pd.DataFrame:
    """For ids with descendant diseases: which diseases the retrieved top evidence is annotated to.

    Counts evidence rows only (as retrieved: top targets per datatype); it shows whether a
    condition-level target list is driven by the queried disease itself or by descendants
    (e.g. autonomic nervous system disorder -> paraganglioma genetics). No score is computed.
    """
    cols = ["condition_id", "ontology_id", "source_disease_label", "n_descendants_in_source",
            "evidence_rows", "evidence_rows_from_descendants", "pct_from_descendants",
            "genetic_evidence_rows", "genetic_rows_from_descendants", "top_evidence_diseases",
            "top_genetic_evidence_diseases"]
    if df.empty or "evidence_disease_label" not in df:
        return pd.DataFrame(columns=cols)
    ev = df[(df["entity_type"] == "gene") & (df["source_evidence_category"] != ASSOC_CATEGORY)]
    ev = ev[pd.to_numeric(ev["n_descendants_in_source"], errors="coerce").fillna(0) > 0]
    out = []
    gen_cats = {"genetic_association", "genetic_literature"}

    def top(s: pd.Series) -> str:
        return "; ".join(f"{k} ({v})" for k, v in s.value_counts().head(top_k).items())

    for (cid, oid, lab), g in ev.groupby(["condition_id", "ontology_id", "source_disease_label"], sort=False):
        desc = g["evidence_from_descendant"].astype(str).str.lower() == "true"
        gg = g[g["source_evidence_category"].isin(gen_cats)]
        gdesc = gg["evidence_from_descendant"].astype(str).str.lower() == "true"
        out.append({"condition_id": cid, "ontology_id": oid, "source_disease_label": lab,
                    "n_descendants_in_source": int(g["n_descendants_in_source"].iloc[0]),
                    "evidence_rows": int(len(g)), "evidence_rows_from_descendants": int(desc.sum()),
                    "pct_from_descendants": round(100 * desc.mean(), 1) if len(g) else np.nan,
                    "genetic_evidence_rows": int(len(gg)), "genetic_rows_from_descendants": int(gdesc.sum()),
                    "top_evidence_diseases": top(g["evidence_disease_label"]),
                    "top_genetic_evidence_diseases": top(gg["evidence_disease_label"]) if len(gg) else ""})
    return pd.DataFrame(out, columns=cols)


# --------------------------------------------------------------------------- audit + registry
def _md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    def cell(v) -> str:
        if isinstance(v, (list, tuple, set, np.ndarray)):
            return "; ".join(map(str, v))
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return ""
        return str(v).replace("|", "/").replace("\n", " ")
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(cell(v) for v in r.values) + " |")
    return "\n".join(lines)


def _missingness(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    for c in cols:
        if c not in df:
            continue
        s = df[c]
        miss = s.isna() | (s.astype(str).isin(["", "None", "nan"]))
        rows.append({"column": c, "rows": len(s), "missing": int(miss.sum()), "missing_pct": round(100 * miss.mean(), 1)})
    return pd.DataFrame(rows)


def write_audit(df: pd.DataFrame, cov: pd.DataFrame, run_log: dict, summary: dict) -> None:
    per = summary["per_condition"]
    meta = run_log["meta"]
    counts = df.groupby(["evidence_type", "entity_type"]).size().reset_index(name="rows")
    a = df[df["source_evidence_category"] == ASSOC_CATEGORY]
    ev = df[(df["entity_type"] == "gene") & (df["source_evidence_category"] != ASSOC_CATEGORY)]
    drugs = df[df["entity_type"] == "drug"]
    miss_a = _missingness(a, ["source_score", "direct_association_score", "ot_datatype_score__genetic_association",
                              "ot_datatype_score__literature", "ot_datatype_score__animal_model",
                              "ot_datatype_score__clinical", "ot_datatype_score__rna_expression",
                              "ot_datatype_score__affected_pathway"])
    miss_e = _missingness(ev, ["study_population", "sample_size", "evidence_direction", "variant_rsids", "p_value",
                               "literature_pmids", "drug_id", "reactome_id"])
    miss_d = _missingness(drugs, ["mechanism_of_action", "action_types", "max_clinical_stage_for_indication",
                                  "clinical_report_ids"])
    desc = pd.DataFrame([{"source_disease_id": k, **v} for k, v in run_log["default_includes_descendants"].items()])
    ev_log = pd.DataFrame(run_log["evidence_log"])
    cov_view = cov[["condition_id", "ontology_id", "ontology_id_role", "ontology_match", "id_status_in_source",
                    "source_disease_label", "coverage_status", "n_records"]]
    lab = df.drop_duplicates(["ontology_id"])[["condition_id", "ontology_id", "registry_label", "source_disease_label",
                                                "label_differs_from_registry", "n_descendants_in_source"]]
    L = []
    L.append("# DATA AUDIT — Open Targets Platform (GraphQL API v4)\n")
    L.append("| Field | Value |\n|---|---|")
    L.append(f"| source_id | {SOURCE_ID} |")
    L.append(f"| Source (dataset/API name, exact files/endpoints) | Open Targets Platform GraphQL API, POST {API_URL} "
             "(queries: meta, disease, disease.associatedTargets [default and enableIndirect:false], disease.evidences, "
             "disease.drugAndClinicalCandidates, targets.pathways) |")
    L.append("| Publishing organization | Open Targets (EMBL-EBI, Wellcome Sanger Institute, partners) |")
    L.append(f"| Retrieval date (UTC) | {Q.retrieval_window(SOURCE_ID)} |")
    L.append(f"| Source version / release | {run_log['version']} (meta.name: {meta['name']}) |")
    L.append("| Source update date / cadence | quarterly Platform data releases (release above is the one served at retrieval) |")
    L.append("| License / access conditions | open API, no key; Open Targets Platform data are released under CC0 1.0 |")
    L.append("| Unit of observation | target (gene) x disease association; evidence item; drug x disease indication; "
             "Reactome pathway annotation of associated targets |")
    L.append(f"| Sample size (actual, as ingested) | {len(df)} rows: {int(len(a))} association rows "
             f"({summary['totals']['unique_targets']} unique targets), {len(ev)} top-evidence rows, "
             f"{len(drugs)} drug rows, {int((df['entity_type'] == 'pathway').sum())} pathway rows |")
    L.append("| Geography (resolution, vintage) | none (concept-level) |")
    L.append("| Person-level? | no |\n| Geographic? | no |\n| Omics? | yes — condition-level molecular evidence only |")
    L.append("| Wearable? | no |")
    L.append("| True participant linkage across modalities? | no — disease-concept level; no participants |\n")
    L.append("## Files / endpoints retrieved\n")
    L.append(f"Every GraphQL response is saved under `data/raw/{SOURCE_ID}/` (meta.json, diseases/, associations/, "
             "evidence/, drugs/, target_pathways/, datasource_probe/) and listed with sha256, bytes and the GraphQL "
             "variables in MANIFEST.json. `run_log.json` holds the measured datasource->datatype map, the evidence "
             "retrieval log and the descendant check. Responses are also cached by `measure_it.http`.\n")
    L.append("## Descendant diseases: does the association call include them by default?\n")
    L.append("Measured per id: the association count with the API default equals the count with "
             "`enableIndirect: true` and exceeds the direct-only count wherever the disease has descendants. "
             "The default call is therefore **indirect (descendants included)**; the partition keeps the default "
             "scores and adds `direct_association_score` / `has_direct_association` from the direct-only call.\n")
    L.append(_md_table(desc))
    L.append("\nPagination (measured): paging the API's default score-descending order across the 3000-row page "
             "boundary repeated some targets and skipped others, reproducibly (e.g. endometriosis 3561 rows but 3521 "
             "unique targets). Ids with more than 3000 associations are paged with `orderByScore: \"score asc\"`, which "
             "returned every target exactly once, and re-sorted by score; the module fails if the retrieved unique "
             "targets differ from the API count.\n")
    L.append(_md_table(pd.DataFrame([{"id_call": k, **v} for k, v in run_log.get("association_paging", {}).items()])))
    L.append("\n## Identifier coverage (registry ids -> Open Targets disease index)\n")
    L.append(_md_table(cov_view))
    L.append("\nLabel check (Open Targets label vs the registry/Mondo label of the same id):\n")
    L.append(_md_table(lab))
    L.append("\n## Per-condition counts\n")
    L.append("`unique_targets` counts every target associated with any registry id of the condition (including "
             "descendants and broad/component ids); `unique_targets_exact_or_narrow_ids` restricts to ids whose registry "
             "predicate is exact or narrow (a broad primary id such as post-infectious disorder is excluded). "
             "`unique_targets_primary_id` / `genetic_association_targets_primary_id` restrict to the condition's "
             "registry primary id. "
             "`direct_targets` = targets with a direct (non-propagated) association. Genetic-association targets = "
             "targets with a non-zero "
             "`genetic_association` datatype score.\n")
    L.append(_md_table(per))
    L.append(f"\nConditions with zero Open Targets associations: {', '.join(summary['zero_evidence']) or 'none'}.\n")
    drivers = descendant_evidence_drivers(df)
    L.append("## Which diseases drive the evidence of ids with descendants (measured)\n")
    L.append("Retrieved top-evidence rows of every id that has descendant diseases in Open Targets, by the disease "
             "each evidence item is annotated to. Where most rows come from descendants, the condition-level "
             "target counts above describe those descendants, not the condition itself; use "
             "`unique_targets_primary_id` or `has_direct_association` for the condition itself.\n")
    L.append(_md_table(drivers) if len(drivers) else "(no id with descendants has evidence rows)")
    widened = per[(per["unique_targets"] > 0) & (per["unique_targets_primary_id"] < 0.5 * per["unique_targets"])]
    if len(widened):
        L.append("\nConditions whose condition-level target count is mostly from non-primary ids (component/"
                 "narrower/broad ids and their descendants): "
                 + "; ".join(f"{r.condition_id} {r.unique_targets_primary_id} of {r.unique_targets} targets "
                             f"({r.genetic_association_targets_primary_id} of {r.unique_genetic_association_targets} "
                             "genetic) on the primary id" for r in widened.itertuples()) + ".")
    desc_driven = drivers[pd.to_numeric(drivers["pct_from_descendants"], errors="coerce") >= 80]
    if len(desc_driven):
        L.append("\nIds whose retrieved evidence is at least 80% annotated to descendant diseases (their target "
                 "lists describe the descendants): "
                 + "; ".join(f"{r.condition_id} {r.ontology_id} {r.source_disease_label} "
                             f"({r.pct_from_descendants}% of {r.evidence_rows} evidence rows; direct-only targets: "
                             f"{run_log['default_includes_descendants'].get(Q.curie_to_short(r.ontology_id), {}).get('direct', 'n/a')})"
                             for r in desc_driven.itertuples()) + ".")
    L.append("")
    L.append("## Datatypes and datasources (measured)\n")
    L.append("Datasource -> datatype map measured from one evidence item per datasource: "
             + ", ".join(f"`{k}`->`{v}`" for k, v in sorted(run_log["datasource_to_datatype"].items())) + ".\n")
    L.append("Requested datatype -> API datatype used: "
             + ", ".join(f"{k}->{v}" for k, v in run_log["requested_to_api_datatype"].items())
             + ". In this release the drug datatype is named `clinical` (datasource `clinical_precedence`); the "
             "`knownDrugs` field no longer exists and `drugAndClinicalCandidates` replaces it.\n")
    L.append("Top-evidence retrieval log (10 highest-scoring targets per datatype; up to "
             f"{EVIDENCE_PAGE_SIZE * MAX_EVIDENCE_PAGES} evidence items fetched; top {TOP_EVIDENCE_PER_TARGET} per "
             "target kept):\n")
    L.append(_md_table(ev_log) if len(ev_log) else "(none)")
    L.append("\n## Key variables\n")
    L.append("* `source_score` — association rows: Open Targets overall association score (0-1, source-provided, "
             "label 'Open Targets association score (source-provided)'); evidence rows: the evidence item score "
             "(source-provided). Drug and pathway rows have no source score (null).")
    L.append("* `ot_datatype_score__<datatype>` — per-datatype association scores exactly as returned "
             "(genetic_association, genetic_literature, somatic_mutation, clinical, affected_pathway, rna_expression, "
             "literature, animal_model); null = the datatype contributed no evidence.")
    L.append("* `ot_datasource_scores_json` — per-datasource scores as returned.")
    L.append("* `direct_association_score`, `has_direct_association` — from the `enableIndirect:false` call.")
    L.append("* evidence rows: `ot_datasource`, `evidence_disease_id` (the disease the evidence is annotated to; "
             "`evidence_from_descendant`), `study_id`, `variant_rsids`, `p_value`, `odds_ratio`, `beta`, "
             "`literature_pmids`, `drug_id`, `clinical_stage`, `reactome_id`, `evidence_direction` "
             "(`target:GoF/LoF`, `trait:risk/protect`, only where Open Targets defines it).")
    L.append("* drug rows: `max_clinical_stage_for_indication`, `drug_maximum_clinical_stage_any_indication`, "
             "`mechanism_of_action`, `action_types`, `mechanism_target_symbols`, `trial_status_counts`, "
             "`reports_cite_queried_id` (false = the indication reaches the queried id through a related disease in "
             "the clinical reports).")
    L.append("* pathway rows: `n_associated_targets_annotated`, `n_genetic_association_targets_annotated`, "
             "`top_target_symbols`, `reactome_top_level_term` (annotation counts, no enrichment test).\n")
    L.append("## Missingness\n")
    L.append("Association rows:\n")
    L.append(_md_table(miss_a))
    L.append("\nTop-evidence rows:\n")
    L.append(_md_table(miss_e))
    L.append("\nDrug rows:\n")
    L.append(_md_table(miss_d) if len(miss_d) else "(none)")
    L.append("\n## Linkage strategy\n")
    L.append("Registry condition -> MONDO/EFO ids from `condition_registry.mondo_ids`/`efo_ids` (ontology module) "
             "-> Open Targets disease id (same CURIE with '_' separator). No text matching and no guessed ids: an "
             "id absent from the Open Targets disease index is recorded as 'identifier not present in source'. Targets "
             "join to other tables by Ensembl gene id only. Nothing here is joinable to participants.\n")
    L.append("## Limitations and caveats\n")
    L.append("* Associations are indirect by default: a broad id (e.g. autonomic nervous system disorder, "
             "Ehlers-Danlos syndrome, post-infectious disorder) pulls in every descendant disease, including "
             "monogenic subtypes and neoplasms (see the measured table of evidence diseases above). "
             "`ontology_match`/`ontology_id_role`, `unique_targets_primary_id` and `has_direct_association` must be "
             "used before interpreting a condition's target list.")
    L.append("* Literature (Europe PMC text mining) and animal-model (IMPC) evidence dominate target counts for "
             "poorly characterised conditions; they are co-mention/phenotype-similarity signals, not mechanism.")
    L.append("* Open Targets scores are relative, source-specific heuristics; they are not probabilities and "
             "must not be combined with GWAS/GEO counts into a single score.")
    L.append("* For POTS the only Open Targets id is MONDO_0011479, which the registry/Mondo labels 'POTS due to "
             "NET deficiency' (a monogenic form, predicate narrow) while Open Targets labels it 'postural "
             "orthostatic tachycardia syndrome'.")
    L.append("* `drugAndClinicalCandidates` is not documented as ontology-propagated; `reports_cite_queried_id` "
             "records whether the linked clinical reports name the queried id itself.")
    L.append("* Condition-level molecular enrichment/evidence only — never patient multi-omics.\n")
    L.append("## Processed outputs\n")
    L.append(f"* `{PARTITION}`: {len(df)} rows\n" + _md_table(counts))
    L.append(f"* `{COVERAGE}`: {len(cov)} rows")
    L.append("* `condition_molecular_evidence` / `condition_molecular_coverage`: union of all partitions "
             "(store.union_partitions)\n")
    L.append("## Reproduce\n")
    L.append(f"`uv run python -m {MODULE}`\n")
    (raw_dir(SOURCE_ID) / "DATA_AUDIT.md").write_text("\n".join(L))


def write_registry(df: pd.DataFrame, cov: pd.DataFrame, run_log: dict, summary: dict) -> None:
    per = summary["per_condition"].set_index("condition_id")
    write_registry_entry({
        "source_id": SOURCE_ID,
        "name": "Open Targets Platform GraphQL API v4",
        "publisher": "Open Targets (EMBL-EBI, Wellcome Sanger Institute and partners)",
        "landing_url": LANDING_URL,
        "access_urls": [API_URL, DOCS_URL],
        "license": "CC0 1.0 (Open Targets Platform data)",
        "access_conditions": "open API, no key",
        "retrieved_at": run_log["retrieved_at"],
        "source_version": run_log["version"],
        "update_date": f"data release {run_log['meta']['dataVersion']['year']}.{run_log['meta']['dataVersion']['month']}",
        "data_layer": "condition_molecular",
        "unit_of_observation": "target-disease association / evidence item / drug indication / pathway annotation",
        "sample_size": {
            "rows": int(len(df)),
            "association_rows": int((df["source_evidence_category"] == ASSOC_CATEGORY).sum()),
            "unique_targets": summary["totals"]["unique_targets"],
            "per_condition_unique_targets": {k: int(v) for k, v in per["unique_targets"].items()},
            "per_condition_unique_targets_primary_id": {k: int(v) for k, v in per["unique_targets_primary_id"].items()},
            "per_condition_genetic_association_targets": {k: int(v) for k, v in per["unique_genetic_association_targets"].items()},
            "per_condition_known_drugs": {k: int(v) for k, v in per["known_drugs"].items()},
        },
        "geographic_resolution": "none",
        "person_level": False,
        "geographic": False,
        "omics": True,
        "wearable": False,
        "participant_linkage": "not applicable (condition-level evidence, no participants)",
        "true_participant_linkage_across_modalities": False,
        "status": "ingested",  # coverage gaps (e.g. no Open Targets id for a condition) are listed in limitations
        "processed_outputs": [PARTITION, COVERAGE, "condition_molecular_evidence", "condition_molecular_coverage"],
        "audit": f"data/raw/{SOURCE_ID}/DATA_AUDIT.md",
        "ingestion_module": MODULE,
        "limitations": [
            "Associations are indirect by default (descendant diseases included); broad/component ids pull in "
            "monogenic subtypes and neoplasms (per-condition counts are also given for the primary id only)",
            "Scores are Open Targets heuristics (source-provided), not probabilities; no combined omics score",
            "Conditions with no Open Targets identifier: " + (", ".join(summary["zero_evidence"]) or "none"),
            "Literature/animal-model evidence dominates target counts for poorly characterised conditions",
            "Condition-level molecular enrichment/evidence only; not participant-linked",
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
