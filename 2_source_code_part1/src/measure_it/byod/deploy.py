"""`measure-it byod deploy <dir>`: let a user dataset's performance record reach the metric link and the ranking.

    1. metric link rebuilt (measure_it.scoring.metric_link.build(write_csv=False): the committed results/tables CSVs
       are not touched) with the user records; a demo record only with --demo (MEASURE_IT_BYOD_DEMO=1);
    2. every scored (condition, measurement) combination whose performance row changed is re-scored at county and
       state level with the unchanged engine (opportunity.make_combo / combo_frame: same seeds, same Monte Carlo),
       spliced into deployment_opportunities, and the joint (region x measurement) ranks recomputed;
       deployment_candidates is refreshed when the primary demo combination changed;
    3. report: performance before / after, evidence_weighted top-10 before / after, equal-weight ranks unchanged,
       joint top-25 bundle composition before / after, the recommendation under evidence_weighted;
    4. heatmap data (condition x measurement evidence, before / after) and the dataset's own embedding (its Digital
       Phenotype Vector PCA coordinates, within-dataset only), CSV + PNG, under results/byod/<id>/;
    5. guards: validate.check_person_geography() must pass, and no derived table may contain one of the dataset's
       participant ids. Only the aggregate record crosses into the ranking; scoring never reads person tables.

Everything describes candidate deployment opportunities for pilot evaluation, not validated diagnostic pathways.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import UNKNOWN, utc_now_iso
from ..store import read_table, table_exists, write_table
from . import common as K
from . import person as PER

STATE = K.RESULTS_ROOT / "state.json"
OPP_KEYS = ["object_id", "condition_id", "measurement_id", "geo_level", "geo_id", "geo_name", "state_abbr",
            "rank_eligible", "rank_equal", "rank_evidence_weighted", "rank_equal_joint", "rank_evidence_weighted_joint",
            "composite_equal", "composite_evidence_weighted", "measurement_performance_status",
            "measurement_performance_record_id", "evidence_weighted_status"]
PERF_KEYS = ["performance_status", "quality_tier", "selected_record_id", "measurement_evidence", "op_sensitivity",
             "op_specificity", "auroc", "auroc_ci_low"]


def demo_state() -> bool:
    try:
        return bool(json.loads(STATE.read_text()).get("demo"))
    except (OSError, ValueError):
        return False


def _set_demo(demo: bool) -> None:
    os.environ[K.DEMO_ENV] = "1" if demo else "0"
    K.local_only_dir(K.RESULTS_ROOT)
    STATE.write_text(json.dumps({"demo": bool(demo), "updated_at": utc_now_iso()}))


def _same(a, b) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        try:
            fa, fb = float(a), float(b)
        except (TypeError, ValueError):
            return a == b
        return (np.isnan(fa) and np.isnan(fb)) or abs(fa - fb) < 1e-12
    return (a if a == a else None) == (b if b == b else None)


def changed_cells(before: pd.DataFrame, after: pd.DataFrame) -> list[tuple[str, str]]:
    b = before.set_index(["condition_id", "measurement_id"])
    out = []
    for r in after.to_dict("records"):
        key = (r["condition_id"], r["measurement_id"])
        if key not in b.index:
            out.append(key)
            continue
        old = b.loc[key]
        if any(not _same(old.get(k), r.get(k)) for k in PERF_KEYS if k in r):
            out.append(key)
    return out


def rescore(cells: list[tuple[str, str]], echo=print) -> list[tuple[str, str, str]]:
    """Re-score the given (condition, measurement) combinations at every level and splice them in."""
    from ..scoring import opportunity as O
    df = read_table(O.TABLE)
    S = O.spatial()
    done = []
    key = list(zip(df["condition_id"], df["measurement_id"], df["geo_level"]))
    pos = pd.Series(range(len(df)), index=pd.MultiIndex.from_tuples(key)).sort_index()
    new_parts = []
    for c, m in cells:
        for level in O.LEVELS:
            old_idx = np.sort(pos.loc[(c, m, level)].to_numpy()) if (c, m, level) in pos.index else np.array([], int)
            cb = O.make_combo(c, m, level)
            f = O._provenance(O.combo_frame(cb), S)
            if len(old_idx):
                if list(df.iloc[old_idx]["geo_id"]) != list(f["geo_id"]):
                    raise RuntimeError(f"re-scored {c} x {m} x {level}: region order differs from the stored rows")
                f.index = df.index[old_idx]
            else:
                f.index = pd.RangeIndex(len(df) + sum(len(p) for p in new_parts),
                                        len(df) + sum(len(p) for p in new_parts) + len(f))
            new_parts.append(f)
            done.append((c, m, level))
            echo(f"  re-scored {c} x {m} x {level}: {len(f)} rows")
    if not new_parts:
        return done
    new = pd.concat(new_parts)
    keep = df.drop(index=[i for i in new.index if i in df.index])
    out = pd.concat([keep, new.reindex(columns=df.columns)]).sort_index()
    t = pd.to_numeric(out["measurement_performance_tier"], errors="coerce")
    if t.notna().all() and (t % 1 == 0).all():   # integer tiers again once no 2.5 (user-supplied) record is scored
        out["measurement_performance_tier"] = t.astype("int64")
    O.add_joint_ranks(out)
    assert out["object_id"].is_unique
    write_table(out.reset_index(drop=True), O.TABLE, producer=f"{O.PRODUCER} (measure_it.byod re-score)",
                description="Deployment-opportunity components, composites under every weight set, and Monte Carlo "
                            f"rank intervals; combinations re-scored by measure_it.byod: {sorted(set(cells))}. "
                            "docs/SCORING.md")
    return done


def refresh_rankings(demo: bool, echo=print) -> dict:
    """Metric link + re-scoring of the combinations whose performance row changed. Returns before/after frames."""
    from ..scoring import metric_link as ML
    from ..scoring import opportunity as O
    _set_demo(demo)
    before_perf = read_table(ML.TABLE) if table_exists(ML.TABLE) else pd.DataFrame(columns=["condition_id",
                                                                                            "measurement_id"])
    before_opp = read_table(O.TABLE, columns=OPP_KEYS)
    res = ML.build(write=True, write_csv=False)
    ML.records.cache_clear()
    ML._performance_cached.cache_clear()
    after_perf = res["performance"]
    scored = {(c, m) for c in O.all_condition_ids() for m in O.MEASUREMENTS}
    cells = [c for c in changed_cells(before_perf, after_perf) if c in scored]
    echo(f"  metric link rebuilt ({len(after_perf)} rows); scored combinations whose performance changed: "
         f"{cells or 'none'}")
    done = rescore(cells, echo=echo)
    if (O.PRIMARY["condition"], O.PRIMARY["measurement"]) in cells:
        from ..scoring import recommend as REC
        r = REC.rank_deployment_opportunities(**REC.DEMO_ARGS)
        r["query_text"] = REC.DEMO_QUERY
        write_table(REC.deployment_candidates_table(r, REC.DEMO_QUERY), "deployment_candidates", producer=REC.PRODUCER,
                    description="One row per recommendation of the SPEC demo query, with the full SPEC-schema JSON.")
        echo("  deployment_candidates refreshed (the primary demo combination changed)")
    after_opp = read_table(O.TABLE, columns=OPP_KEYS)
    return {"before_perf": before_perf, "after_perf": after_perf, "before_opp": before_opp, "after_opp": after_opp,
            "changed_cells": cells, "rescored": done, "demo": demo}


# ------------------------------------------------------------------------------------------------------------------
# report
# ------------------------------------------------------------------------------------------------------------------

def _perf_brief(p: pd.DataFrame, c: str, m: str) -> dict:
    r = p[(p["condition_id"] == c) & (p["measurement_id"] == m)]
    if not len(r):
        return {"performance_status": UNKNOWN, "note": "no row"}
    r = r.iloc[0]
    out = {k: (None if isinstance(r.get(k), float) and not np.isfinite(r.get(k)) else r.get(k)) for k in PERF_KEYS}
    out["quality_tier_label"] = r.get("quality_tier_label")
    return out


def _top(g: pd.DataFrame, col: str, n: int = 10) -> list[dict]:
    t = g[g[col].notna()].sort_values(col).head(n)
    return [{"rank": int(r[col]), "geo_id": r["geo_id"], "name": r["geo_name"]} for r in t.to_dict("records")]


def combo_change(before: pd.DataFrame, after: pd.DataFrame, c: str, m: str, level: str) -> dict:
    sel = lambda d: d[(d["condition_id"] == c) & (d["measurement_id"] == m) & (d["geo_level"] == level)]  # noqa
    b, a = sel(before).set_index("geo_id"), sel(after).set_index("geo_id")
    ew = "rank_evidence_weighted"
    out = {"condition_id": c, "measurement_id": m, "geo_level": level,
           "n_ranked_evidence_weighted_before": int(b[ew].notna().sum()),
           "n_ranked_evidence_weighted_after": int(a[ew].notna().sum()),
           "evidence_weighted_status_after": str(a["evidence_weighted_status"].dropna().iloc[0])
           if a["evidence_weighted_status"].notna().any() else "",
           "equal_weight_ranks_unchanged": bool((b["rank_equal"].reindex(a.index).fillna(-1)
                                                == a["rank_equal"].fillna(-1)).all()),
           "top10_evidence_weighted_before": _top(b.reset_index(), ew),
           "top10_evidence_weighted_after": _top(a.reset_index(), ew)}
    tb = {x["geo_id"] for x in out["top10_evidence_weighted_before"]}
    ta = {x["geo_id"] for x in out["top10_evidence_weighted_after"]}
    out["top10_overlap"] = len(tb & ta) if tb and ta else None
    common = b[ew].notna() & a[ew].reindex(b.index).notna()
    if common.sum() > 2:
        from scipy.stats import spearmanr
        out["spearman_evidence_weighted_common"] = float(spearmanr(b.loc[common, ew],
                                                                   a[ew].reindex(b.index)[common]).statistic)
    # joint ranking: where does this measurement enter the region x measurement list for the condition?
    jb = before[(before["condition_id"] == c) & (before["geo_level"] == level)]
    ja = after[(after["condition_id"] == c) & (after["geo_level"] == level)]
    col = "rank_evidence_weighted_joint"
    for tag, j in (("before", jb), ("after", ja)):
        top = j[j[col] <= 25]
        out[f"joint_top25_composition_{tag}"] = top["measurement_id"].value_counts().to_dict()
        mm = j[(j["measurement_id"] == m) & j[col].notna()]
        out[f"best_joint_rank_{tag}"] = int(mm[col].min()) if len(mm) else None
        out[f"n_joint_pairs_{tag}"] = int(j[col].notna().sum())
    return out


def relevant_cells(rec: dict) -> list[tuple[str, str]]:
    from ..scoring import opportunity as O
    from ..scoring.metric_link import _measurement
    cls = rec["primary_class"]
    cells = [(rec["target_condition"], b) for b in O.MEASUREMENTS if cls in _measurement(b)["members"]]
    return cells + [(rec["target_condition"], cls)]


def deploy_report(ds: str, rec: dict | None, ref: dict) -> dict:
    rows = []
    cells = relevant_cells(rec) if rec else []
    for c, m in cells:
        b, a = _perf_brief(ref["before_perf"], c, m), _perf_brief(ref["after_perf"], c, m)
        uses = str(a.get("selected_record_id") or "").startswith("BYOD-")
        why = ("the user record is the selected record" if uses else
               f"the selected record stays {a.get('selected_record_id')} (tier order: a better tier is present)"
               if a.get("selected_record_id") else "no usable record")
        if rec and rec.get("demo") and not ref["demo"]:
            why = "demo record excluded (run with --demo to include it)"
        rows.append({"condition_id": c, "measurement_id": m, "before": b, "after": a, "uses_user_record": uses,
                     "why": why})
    changes = [combo_change(ref["before_opp"], ref["after_opp"], c, m, lv) for c, m, lv in ref["rescored"]]
    return {"dataset_id": ds, "at": utc_now_iso(), "demo_mode": ref["demo"],
            "record_id": rec.get("record_id") if rec else None, "relevant_cells": rows,
            "changed_scored_combinations": [list(x) for x in ref["changed_cells"]], "ranking_changes": changes,
            "framing": "Candidate deployment opportunities for pilot evaluation, not validated diagnostic pathways. "
                       "Only the aggregate performance record crossed into the ranking."}


def report_markdown(rep: dict, rec: dict | None) -> str:
    L = [f"# Deploy report: {rep['dataset_id']}", "", f"{rep['at']}; demo mode {rep['demo_mode']}. {rep['framing']}"
         + (f" (Nothing changed in this run; the comparison below is from the deploy of "
            f"{rep['ranking_changes_from_deploy_at']}.)" if rep.get("ranking_changes_from_deploy_at") else ""), ""]
    if rec:
        L.append(f"Record `{rec['record_id']}`: {rec['target_condition']} x {rec['primary_class']}, AUROC "
                 f"{rec['auroc']:.3f} ({rec['auroc_ci_low']:.3f}-{rec['auroc_ci_high']:.3f}), tier "
                 "user_supplied_own_computation" + (" (demo)" if rec.get("demo") else "") + ".")
    L += ["", "## Performance cells", "", "| condition | measurement | before | after | why |", "|---|---|---|---|---|"]
    for r in rep["relevant_cells"]:
        f = lambda p: (f"{p.get('performance_status')} / {p.get('selected_record_id') or '-'} / evidence "  # noqa
                       f"{p['measurement_evidence']:.3f}" if isinstance(p.get("measurement_evidence"), float)
                       else f"{p.get('performance_status')} / {p.get('selected_record_id') or '-'}")
        L.append(f"| {r['condition_id']} | {r['measurement_id']} | {f(r['before'])} | {f(r['after'])} | {r['why']} |")
    for key, r in (rep.get("recommendations_evidence_weighted") or {}).items():
        L += ["", f"## Recommendation under evidence_weighted: {key.replace('|', ' x ')} (county, top 3)", ""]
        for i, x in enumerate(r.get("top3") or [], 1):
            L.append(f"{i}. {x['name']}: candidate partners {', '.join(x['candidate_sites']) or UNKNOWN}")
        if r.get("top3"):
            L += ["", "Uncertainty carried by every recommendation: " + "; ".join(
                str(u) for u in r["top3"][0]["uncertainties"])[:900]]
    L += ["", "## Ranking changes (re-scored combinations)", ""]
    if not rep["ranking_changes"]:
        L.append("No scored combination changed.")
    for ch in rep["ranking_changes"]:
        L += [f"### {ch['condition_id']} x {ch['measurement_id']} ({ch['geo_level']})", "",
              f"* regions ranked under evidence_weighted: {ch['n_ranked_evidence_weighted_before']} -> "
              f"{ch['n_ranked_evidence_weighted_after']}; equal-weight ranks unchanged: "
              f"{ch['equal_weight_ranks_unchanged']}",
              f"* evidence_weighted top-10 before: {' | '.join(x['name'] for x in ch['top10_evidence_weighted_before']) or 'not ranked'}",
              f"* evidence_weighted top-10 after: {' | '.join(x['name'] for x in ch['top10_evidence_weighted_after']) or 'not ranked'}",
              f"* joint (region x measurement) top-25 under evidence_weighted, by measurement: before "
              f"{ch['joint_top25_composition_before']}, after {ch['joint_top25_composition_after']}; best joint rank "
              f"of this measurement {ch['best_joint_rank_before']} -> {ch['best_joint_rank_after']}", ""]
    return "\n".join(L) + "\n"


# ------------------------------------------------------------------------------------------------------------------
# heatmap / embedding data
# ------------------------------------------------------------------------------------------------------------------

def heatmap_data(before: pd.DataFrame, after: pd.DataFrame) -> pd.DataFrame:
    from ..scoring import opportunity as O
    keep = lambda d: d[d["condition_id"].isin(O.all_condition_ids())]  # noqa
    a = keep(after)[["condition_id", "measurement_id", "measurement_evidence", "performance_status",
                     "selected_record_id", "quality_tier"]]
    b = keep(before)[["condition_id", "measurement_id", "measurement_evidence", "selected_record_id"]].rename(
        columns={"measurement_evidence": "measurement_evidence_before", "selected_record_id": "selected_record_before"})
    h = a.merge(b, on=["condition_id", "measurement_id"], how="left")
    h["user_record"] = h["selected_record_id"].astype(str).str.startswith("BYOD-")
    h["changed"] = [not _same(x, y) for x, y in zip(h["measurement_evidence"], h["measurement_evidence_before"])]
    return h


def embedding_data(ds: str, label: str) -> pd.DataFrame:
    import pyarrow.dataset as pads

    from ..config import PROCESSED
    t = pads.dataset(PROCESSED / "participant_phenotype_embeddings.parquet").to_table(
        columns=["participant_id", "dataset_id", "block", "dimension", "dimension_type", "value"],
        filter=(pads.field("dataset_id") == ds) & (pads.field("dimension_type") == "pca")).to_pandas()
    if t.empty:
        return pd.DataFrame()
    t["block"] = t["block"].astype(str)
    w = t.pivot_table(index=["participant_id", "block"], columns="dimension", values="value").reset_index()
    rows = []
    for b, g in w.groupby("block"):
        pc1, pc2 = f"pca:{b}:pc1", f"pca:{b}:pc2"
        if pc1 in g and pc2 in g:
            rows.append(pd.DataFrame({"participant_id": g["participant_id"], "block": b, "pc1": g[pc1],
                                      "pc2": g[pc2]}))
    e = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(e):
        P = read_table(K.partition_name("participants", ds), columns=["participant_id", label])
        e = e.merge(P.rename(columns={label: "label"}), on="participant_id", how="left")
        e = e.dropna(subset=["pc1", "pc2"])
    return e


def figures(out_dir: Path, heat: pd.DataFrame, emb: pd.DataFrame, ds: str) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paths = []
    from ..scoring import opportunity as O
    if len(heat):
        cols = [m for m in O.MEASUREMENTS if m in set(heat["measurement_id"])] + sorted(
            set(heat["measurement_id"]) - set(O.MEASUREMENTS))
        rows = [c for c in O.all_condition_ids() if c in set(heat["condition_id"])]
        M = heat.pivot_table(index="condition_id", columns="measurement_id", values="measurement_evidence").reindex(
            index=rows, columns=cols)
        U = heat.pivot_table(index="condition_id", columns="measurement_id", values="user_record",
                             aggfunc="max").reindex(index=rows, columns=cols).fillna(False)
        fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(cols)), 0.5 * len(rows) + 2.2))
        cmap = matplotlib.colormaps["viridis"].with_extremes(bad="#d9d9d9")
        im = ax.imshow(np.ma.masked_invalid(M.to_numpy(float)), cmap=cmap, vmin=0, vmax=1, aspect="auto")
        for i in range(len(rows)):
            for j in range(len(cols)):
                v = M.iat[i, j]
                ax.text(j, i, "UNK" if not np.isfinite(v) else f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if (np.isfinite(v) and v < 0.5) else "black")
                if bool(U.iat[i, j]):
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="#e4572e", lw=2.2))
        ax.set_xticks(range(len(cols)), [c.replace("_", " ") for c in cols], rotation=60, ha="right", fontsize=7)
        ax.set_yticks(range(len(rows)), [r.replace("_", " ") for r in rows], fontsize=7)
        ax.set_title(f"measurement_evidence after adding {ds} (red outline: the user record is selected; UNK = "
                     "UNKNOWN, never 0)", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.02, label="measurement_evidence (tier-discounted)")
        fig.tight_layout()
        p = out_dir / "heatmap_measurement_evidence.png"
        fig.savefig(p, dpi=130)
        plt.close(fig)
        paths.append(str(p))
    if len(emb):
        present = list(dict.fromkeys(emb["block"]))
        pref = ["wearable", "clinical_exam"] + [b for b in present if b.startswith("omics_")] + ["labs",
                                                                                                "symptoms_conditions"]
        blocks = [b for b in pref if b in present][:4] or present[:4]
        fig, axes = plt.subplots(1, len(blocks), figsize=(3.6 * len(blocks), 3.4), squeeze=False)
        for ax, b in zip(axes[0], blocks):
            g = emb[emb["block"] == b]
            for lab, col, name in ((0, "#4c78a8", "control"), (1, "#e4572e", "case")):
                gg = g[g["label"] == lab]
                ax.scatter(gg["pc1"], gg["pc2"], s=9, alpha=0.7, c=col, label=name)
            ax.set_title(f"{b} block (within {ds})", fontsize=8)
            ax.set_xlabel("PC1", fontsize=7)
            ax.set_ylabel("PC2", fontsize=7)
        axes[0][0].legend(fontsize=7)
        fig.suptitle("Digital Phenotype Vector PCA (descriptive; labels colour points only, never fitted)", fontsize=8)
        fig.tight_layout()
        p = out_dir / "embedding_dpv_pca.png"
        fig.savefig(p, dpi=130)
        plt.close(fig)
        paths.append(str(p))
    return paths


# ------------------------------------------------------------------------------------------------------------------
# guards
# ------------------------------------------------------------------------------------------------------------------

def no_geography_guard(ds: str) -> dict:
    """person_layer_no_geography passes and no derived / geographic / facility table holds a participant id of ds."""
    from ..scoring import metric_link as ML
    from ..scoring import opportunity as O
    from ..validate import check_person_geography
    chk = check_person_geography()
    leaks = []
    needle = f"{ds}:"
    for name in (O.TABLE, ML.TABLE, ML.RECORDS_TABLE, "deployment_candidates", K.RECORDS_TABLE):
        if not table_exists(name):
            continue
        t = read_table(name)
        for c in t.columns:
            if t[c].dtype == object and t[c].astype(str).str.contains(needle, regex=False).any():
                leaks.append(f"{name}.{c}")
    ok = chk["ok"] and not leaks
    if not ok:
        raise RuntimeError(f"no-geography guard failed: {chk['details'][:3]} participant ids in {leaks}")
    return {"person_layer_no_geography": chk["summary"], "participant_ids_in_derived_tables": leaks}


def deploy(root_or_id: str | Path, *, demo: bool = False, export_static: bool = False, build_db: bool = False,
           echo=print) -> dict:
    p = Path(str(root_or_id))
    ds = K.engine_id(K.load_manifest(K.discover(p).manifest_path)["dataset_id"]) if p.is_dir() else \
        K.engine_id(str(root_or_id))
    meta = PER.dataset_row(ds)
    if meta is None:
        raise RuntimeError(f"{ds} is not ingested")
    from .records import stored_records
    recs = stored_records(include_demo=True)
    rec = recs[recs["dataset_id"] == ds].iloc[0].to_dict() if len(recs) and (recs["dataset_id"] == ds).any() else None
    if rec is None:
        raise RuntimeError(f"{ds} has no performance record: run `measure-it byod evaluate` first")
    if rec.get("demo") and not demo:
        echo(f"byod deploy {ds}: this is a demo record; it stays EXCLUDED from the rankings (pass --demo to include "
             "it). Rebuilding the metric link without it.")
    else:
        echo(f"byod deploy {ds}: record {rec['record_id']} enters the metric link" + (" (demo mode)" if demo else ""))
    ref = refresh_rankings(demo, echo=echo)
    rep = deploy_report(ds, rec, ref)
    out = K.local_only_dir(K.results_path(ds))
    old_path = out / "deploy_report.json"
    if not ref["changed_cells"] and old_path.exists():
        # a repeat deploy that changed nothing keeps the comparison of the deploy that did (same record, same mode)
        old = json.loads(old_path.read_text())
        if old.get("record_id") == rep["record_id"] and old.get("demo_mode") == rep["demo_mode"] \
                and old.get("ranking_changes"):
            rep["ranking_changes"] = old["ranking_changes"]
            rep["relevant_cells"] = old["relevant_cells"]
            rep["changed_scored_combinations"] = old["changed_scored_combinations"]
            rep["ranking_changes_from_deploy_at"] = old.get("ranking_changes_from_deploy_at", old.get("at"))
    heat = heatmap_data(ref["before_perf"], ref["after_perf"])
    heat.to_csv(out / "heatmap_measurement_evidence.csv", index=False)
    emb = embedding_data(ds, meta["primary_label"])
    if len(emb):
        emb.to_csv(out / "embedding_dpv_pca.csv", index=False)
    rep["figures"] = figures(out, heat, emb, ds)
    rep["guard"] = no_geography_guard(ds)
    # the recommendation for the user's cell(s) under evidence_weighted (where ranked)
    from ..scoring import recommend as REC
    recos = {}
    for ch in rep["ranking_changes"]:
        if ch["geo_level"] != "county":
            continue
        r = REC.rank_deployment_opportunities(ch["condition_id"], ch["measurement_id"], "county", top_n=3,
                                              weight_set="evidence_weighted", with_context=True)
        recos[f"{ch['condition_id']}|{ch['measurement_id']}"] = {
            "status": r.get("status"), "reason": r.get("reason"),
            "measurement_performance": {k: (r.get("measurement_performance") or {}).get(k) for k in (
                "object_id", "performance_status", "quality_tier", "selected_record_id", "tier_text")},
            "top3": [{"name": x["geography"]["name"], "fips": x["geography"]["fips"],
                      "candidate_sites": [c.get("facility_name") for c in x.get("candidate_sites", [])][:3],
                      "uncertainties": [u for u in x.get("uncertainties", []) if "Measurement performance" in str(u)],
                      "recommended_next_step": x.get("recommended_next_step")}
                     for x in r.get("recommendations", [])]}
    rep["recommendations_evidence_weighted"] = recos
    (out / "deploy_report.json").write_text(json.dumps(rep, indent=1, default=str))
    (out / "DEPLOY_REPORT.md").write_text(report_markdown(rep, rec))
    from .ingest import upsert_dataset_row
    m = {k: v for k, v in meta.items() if k not in ("data_layer", "source_name", "source_record_id", "source_version",
                                                    "retrieved_at", "source_geographic_resolution", "evidence_type",
                                                    "evidence_level", "provenance_notes")}
    included = not (rec.get("demo") and not demo)
    m.update({"status": "deployed" if included else "evaluated (demo record excluded)",
              "deployed_at": utc_now_iso(), "deployed_with_demo": bool(demo)})
    upsert_dataset_row(m)
    if export_static:
        from ..export.static_site import build_static_site
        rep["static_snapshot"] = build_static_site(out / "static_snapshot" / "index.html")["path"]
    if build_db:
        from ..database import build_db as _b
        _b()
    for ch in rep["ranking_changes"]:
        echo(f"  {ch['condition_id']} x {ch['measurement_id']} ({ch['geo_level']}): ranked under evidence_weighted "
             f"{ch['n_ranked_evidence_weighted_before']} -> {ch['n_ranked_evidence_weighted_after']} regions; "
             f"best joint rank {ch['best_joint_rank_before']} -> {ch['best_joint_rank_after']}; equal-weight ranks "
             f"unchanged {ch['equal_weight_ranks_unchanged']}")
    echo(f"  guard: {rep['guard']['person_layer_no_geography']}; participant ids in derived tables: none")
    echo(f"  report: {out / 'DEPLOY_REPORT.md'}")
    return rep
