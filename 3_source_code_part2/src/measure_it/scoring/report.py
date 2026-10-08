"""Run the scoring analyses (Tests 5-7, sensitivity S1-S6, burden-rule before/after, demo) -> results/SCORING_RESULTS.md.

    uv run python -m measure_it.scoring.report            # everything (table build + tests + demo + report)
    uv run python -m measure_it.scoring.report --no-build # reuse data/processed/deployment_opportunities

Every number in results/SCORING_RESULTS.md is read from a table this module (or measure_it.scoring.*) wrote in the
same run. Plan: docs/ANALYSIS_PLAN_SCORING.md.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from ..config import FIGURES, PROJECT_ROOT, RESULTS, TABLES, UNKNOWN, utc_now_iso
from ..store import read_table, table_exists
from . import controls as C
from . import metric_link as ML
from . import opportunity as O
from . import recommend as REC
from . import sensitivity as SENS

PKEY = (O.PRIMARY["condition"], O.PRIMARY["measurement"], O.PRIMARY["level"])
T = {}  # tables written in this run (name -> DataFrame)


def _w(df: pd.DataFrame, name: str) -> pd.DataFrame:
    TABLES.mkdir(parents=True, exist_ok=True)
    df.to_csv(TABLES / f"{name}.csv", index=False)
    T[name] = df
    return df


def _combos(levels=O.LEVELS):
    for level in levels:
        for mid in O.MEASUREMENTS:
            for cid in O.all_condition_ids():
                yield cid, mid, level


# ------------------------------------------------------------------------------------------------------------------
# analyses
# ------------------------------------------------------------------------------------------------------------------

def run_test5() -> None:
    summ, moves = [], []
    for cid, mid, level in _combos():
        cb = O.make_combo(cid, mid, level)
        ev = O.evaluate(cb)
        s, m = C.test5_combo(cb, ev)
        summ += s
        moves += m
        s, m = C.test5_evidence_combo(cb, ev)   # T5-4 (docs/ANALYSIS_PLAN_METRIC_LINK.md section 6)
        summ += s
        moves += m
    _w(pd.DataFrame(summ), "test5_contrasts")
    _w(pd.DataFrame(moves), "test5_movements")
    # descriptive (added after the plan): how the components co-vary (the desert re-uses burden and SVI, and its
    # access terms oppose capacity/readiness), and how many regions have zero condition trials / NIH projects
    df = read_table(O.TABLE)
    comps = [O.PCT[k] for k in O.COMPONENTS] + ["technology_saturation_pct"]
    rows = []
    for (cid, mid, lev), g in df[df["rank_eligible"]].groupby(["condition_id", "measurement_id", "geo_level"]):
        for i, a in enumerate(comps):
            for b in comps[i + 1:]:
                rows.append({"condition_id": cid, "measurement_id": mid, "geo_level": lev, "component_a": a,
                             "component_b": b, "spearman": C.spearman(g[a].to_numpy(float), g[b].to_numpy(float),
                                                                      np.ones(len(g), bool)),
                             "n": int(g[[a, b]].notna().all(axis=1).sum())})
    _w(pd.DataFrame(rows), "test5_component_correlations")
    z = []
    for (cid, mid, lev), g in df[df["rank_eligible"]].groupby(["condition_id", "measurement_id", "geo_level"]):
        zero = (g["rr_condition_trials_pool_n"] == 0) & (g["rr_condition_nih_core_pool_n"] == 0)
        for n in (10, 25, None):
            sub = g if n is None else g[g["rank_equal"] <= n]
            zs = zero.loc[sub.index]
            z.append({"condition_id": cid, "measurement_id": mid, "geo_level": lev,
                      "subset": "all eligible" if n is None else f"equal top-{n}", "n_regions": len(sub),
                      "n_zero_condition_trials_and_nih": int(zs.sum()),
                      "share_zero": float(zs.mean()) if len(sub) else np.nan,
                      "research_readiness_pct_min_among_zero": float(sub.loc[zs, "research_readiness_pct"].min())
                      if zs.any() else np.nan,
                      "research_readiness_pct_max_among_zero": float(sub.loc[zs, "research_readiness_pct"].max())
                      if zs.any() else np.nan})
    _w(pd.DataFrame(z), "test5_readiness_zero_counts")


def run_test6(workers: int = 6) -> None:
    summ, draws = [], []
    for cid in O.all_condition_ids():
        for level in O.LEVELS:
            if level == "state" and cid != O.PRIMARY["condition"]:
                continue
            cb = O.make_combo(cid, O.PRIMARY["measurement"], level)
            for what in ("burden", "vulnerability"):
                for scheme in ("national", "within_state"):
                    if level == "state" and scheme == "within_state":
                        continue
                    if what == "burden" and all(m.level_d.all() for m in cb.members):
                        summ.append({**C.combo_key(cb), "shuffled": what, "scheme": scheme, "n_draws": 0,
                                     "note": "level D: no burden to shuffle"})
                        continue
                    s, d = C.label_shuffle(cb, what, scheme)
                    summ.append(s)
                    draws += d
    _w(pd.DataFrame(summ), "test6_scoring_label_shuffle")
    _w(pd.DataFrame(draws), "test6_scoring_label_shuffle_draws")
    cs, cd = [], []
    for mid in O.MEASUREMENTS:
        s, d = C.clinic_shuffle(O.PRIMARY["condition"], mid, "county", workers=workers)
        cs += s
        cd += d
    _w(pd.DataFrame(cs), "test6_scoring_clinic_shuffle")
    _w(pd.DataFrame(cd), "test6_scoring_clinic_shuffle_draws")
    _w(C.phenotype_label_nulls(), "test6_scoring_phenotype_label_nulls")


def run_test7() -> None:
    summ, pairs, regs = [], [], []
    for cid, mid, level in _combos():
        cb = O.make_combo(cid, mid, level)
        s, p, r = SENS.test7_combo(cb)
        summ.append(s)
        pairs += p
        if (cid, mid, level) == PKEY:
            _w(r[r["eligible"]].sort_values("rank_equal"), "test7_region_stability_primary")
        regs.append(r[r["in_top25_any_named"]])
    _w(pd.DataFrame(summ), "test7_summary")
    _w(pd.DataFrame(pairs), "test7_pairwise_tau")
    reg = pd.concat(regs, ignore_index=True)
    _w(reg, "test7_region_stability_top25_any")
    _w(reg[reg["unstable"]].sort_values(["condition_id", "measurement_id", "geo_level", "rank_equal"]),
       "test7_unstable_regions")


def run_sensitivity() -> None:
    cb = O.make_combo(*PKEY)
    _w(SENS.s1_tau(cb), "test7_sensitivity_s1_tau")
    _w(SENS.s1b_inherited_multiplier(cb), "test7_sensitivity_s1b_inherited_multiplier")
    rows = []
    for cid in ("long_covid", O.QUERY_SET_ID, "autonomic_activity_invisible_illness"):
        a = O.make_combo(cid, O.PRIMARY["measurement"], "county")
        b = O.make_combo(cid, O.PRIMARY["measurement"], "county", burden_override="proxy")
        rows.append(SENS.compare(a, b, "S2 long-COVID county burden = PLACES PHLTH C-proxy (b) vs inherited state (a)"))
    for mid in O.MEASUREMENTS:
        a = O.make_combo(O.PRIMARY["condition"], mid, "county")
        b = O.make_combo(O.PRIMARY["condition"], mid, "county", exclude_groups=("primary_care",))
        rows.append(SENS.compare(a, b, "S4 clinic_capacity without primary care (b) vs with (a)"))
    base = O.make_combo(*PKEY)
    for mid in O.MEASUREMENTS[1:]:
        b = O.make_combo(O.PRIMARY["condition"], mid, "county")
        rows.append(SENS.compare(base, b, "S5 measurement swap: wearable bundle (a) vs other bundle (b)",
                                 {"adapter_b": b.meas.get("adapter_name"), "adapter_status_b": b.meas.get("adapter_status")}))
    _w(pd.DataFrame(rows), "test7_sensitivity_comparisons")
    # S5 with evidence (docs/ANALYSIS_PLAN_METRIC_LINK.md section 6): within-bundle agreement under equal and
    # evidence_weighted, and the joint region x measurement ranking
    _w(SENS.s5_evidence(O.PRIMARY["condition"], "county"), "test7_sensitivity_s5_evidence")
    # S6 specialist_only (declared 2026-09-24, plan deviation 8): primary query at county and state; secondary: every
    # condition x measurement at county level
    s6, reg = [], None
    for cid, mid, level in _combos():
        if level == "state" and (cid, mid) != PKEY[:2]:
            continue
        r, g = SENS.s6_specialist_only(cid, mid, level, regions=(cid, mid, level) == PKEY)
        s6 += r
        reg = g if g is not None else reg
    _w(pd.DataFrame(s6), "test7_sensitivity_specialist_only")
    _w(reg, "test7_sensitivity_specialist_only_top25")
    _w(SENS.s6_provider_qa(), "test7_sensitivity_specialist_only_qa")


def run_burden_rule() -> None:
    """Plan deviation 7: incomplete-burden regions (listed separately) and the before/after of the rule change."""
    summ, tops = [], []
    for cid, mid, level in _combos():
        s, t = SENS.set_burden_rule_before_after(cid, mid, level)
        summ.append(s)
        if s["n_incomplete"]:
            tops.append(t)
    _w(pd.DataFrame(summ), "scoring_set_burden_rule_before_after")
    _w(pd.concat(tops, ignore_index=True) if tops else pd.DataFrame(), "scoring_set_burden_rule_top10")
    df = read_table(O.TABLE)
    inc = df[df["burden_incomplete"].astype(bool)]
    cols = ["condition_id", "measurement_id", "geo_level", "geo_id", "geo_name", "state_abbr", "population_total",
            "small_population_flag", "burden_members_missing", "burden_missing_source_reason",
            "burden_incomplete_reason", "vulnerability_pct",
            "diagnostic_desert_pct", "clinic_capacity_pct", "research_readiness_pct", "object_id"]
    _w(inc[cols].sort_values(["condition_id", "measurement_id", "geo_level", "small_population_flag", "geo_id"]),
       "scoring_incomplete_burden")


def run_metric_link() -> None:
    """Metric -> translation link outputs read by section 12 (plan docs/ANALYSIS_PLAN_METRIC_LINK.md)."""
    df = read_table(O.TABLE)
    _w(SENS.joint_composition(df), "scoring_joint_ranking_composition")
    prim_all = df[(df["condition_id"] == PKEY[0]) & (df["geo_level"] == PKEY[2])]
    jt = prim_all[prim_all[f"rank_{O.EVIDENCE_WEIGHT_SET}_joint"] <= 25].sort_values(f"rank_{O.EVIDENCE_WEIGHT_SET}_joint")
    je = prim_all[prim_all["rank_equal_joint"] <= 25].sort_values("rank_equal_joint")
    cols = ["measurement_id", "geo_id", "geo_name", "state_abbr", "composite_equal", "composite_evidence_weighted",
            "rank_equal", "rank_evidence_weighted", "rank_equal_joint", "rank_evidence_weighted_joint",
            "measurement_performance_status", "measurement_performance_tier", "measurement_evidence",
            "expected_yield", "reach", "object_id"]
    _w(pd.concat([je[cols].assign(list="equal joint top-25"), jt[cols].assign(list="evidence_weighted joint top-25")]),
       "scoring_joint_ranking_top25_primary")
    prim = prim_all[prim_all["measurement_id"] == PKEY[1]]
    cols = ["geo_id", "geo_name", "state_abbr", "population_total", "rank_equal", "rank_evidence_weighted",
            "rank_mc_p05", "rank_mc_p95", "rank_mc_ew_p05", "rank_mc_ew_p95", "p_top10_mc", "p_top10_mc_ew",
            "composite_equal", "composite_evidence_weighted", "burden_pct", "vulnerability_pct",
            "diagnostic_desert_pct", "clinic_capacity_pct", "research_readiness_pct", "measurement_evidence",
            "expected_yield", "reach", "op_sensitivity", "op_specificity", "expected_yield_index",
            "expected_yield_mc_p05", "expected_yield_mc_p95", "reached_adults", "expected_false_positives_upper",
            "expected_false_positives_mc_p05", "expected_false_positives_mc_p95", "object_id"]
    top = prim[(prim["rank_equal"] <= 10) | (prim[f"rank_{O.EVIDENCE_WEIGHT_SET}"] <= 10)]
    _w(top[cols].sort_values(f"rank_{O.EVIDENCE_WEIGHT_SET}"), "metric_link_top10_primary")
    # count path (long COVID, state): expected detectable cases and false positives, every factor a column
    lc = df[(df["condition_id"] == "long_covid") & (df["geo_level"] == "state") &
            (df["measurement_performance_status"] == "known")]
    _w(lc[["measurement_id", "geo_id", "geo_name", "burden_value", "burden_ci_low", "burden_ci_high", "adults_18plus",
           "burden_count", "reach", "op_sensitivity", "op_specificity", "expected_detectable_cases",
           "expected_yield_mc_p05", "expected_yield_mc_p95", "expected_false_positives", "expected_false_positives_mc_p05",
           "expected_false_positives_mc_p95", "expected_ppv", "false_positives_per_detected_case",
           "rank_evidence_weighted", "rank_equal", "measurement_performance_record_id"]]
       .sort_values(["measurement_id", "rank_evidence_weighted"]), "metric_link_long_covid_state_yield")
    # which combinations are ranked under evidence_weighted, and why not
    st = df.groupby(["condition_id", "measurement_id", "geo_level"]).agg(
        performance_status=("measurement_performance_status", "first"),
        tier=("measurement_performance_tier", "first"), record=("measurement_performance_record_id", "first"),
        measurement_evidence=("measurement_evidence", "first"), op_sensitivity=("op_sensitivity", "first"),
        expected_yield_basis=("expected_yield_basis", "first"),
        n_ranked_equal=("rank_equal", lambda x: int(x.notna().sum())),
        n_ranked_evidence_weighted=(f"rank_{O.EVIDENCE_WEIGHT_SET}", lambda x: int(x.notna().sum()))).reset_index()
    _w(st, "metric_link_combination_status")
    # eligible regions dropped from the evidence-weighted ranking for lack of a reach value (plan rule 5.3)
    drop = df[df["rank_eligible"].astype(bool) & (df["measurement_performance_status"] == "known")
              & df[f"composite_{O.EVIDENCE_WEIGHT_SET}"].isna() & (df["expected_yield_basis"] != "level_D")]
    _w(drop[["condition_id", "measurement_id", "geo_level", "geo_id", "geo_name", "reach", "reach_zcta_n",
             "evidence_weighted_status"]], "metric_link_unranked_no_reach")
    _w(pd.read_csv(TABLES / "measurement_performance.csv"), "measurement_performance")
    _w(pd.read_csv(TABLES / "measurement_performance_records.csv"), "measurement_performance_records")


def mc_decomposition() -> None:
    df = read_table(O.TABLE, columns=["condition_id", "measurement_id", "geo_level", "rank_equal", "rank_mc_p05",
                                      "rank_mc_p95", "rank_mc_burden_p05", "rank_mc_burden_p95",
                                      "rank_mc_weights_p05", "rank_mc_weights_p95", "p_top10_mc", "p_top25_mc",
                                      "mc_n_draws", "mc_tau_inherited"])
    rows = []
    for (cid, mid, lev), g in df.groupby(["condition_id", "measurement_id", "geo_level"]):
        top = g[g["rank_equal"] <= 25]
        top10 = g[g["rank_equal"] <= 10]
        rows.append({"condition_id": cid, "measurement_id": mid, "geo_level": lev, "n_top25": len(top),
                     "median_width_full_top25": float((top["rank_mc_p95"] - top["rank_mc_p05"]).median()),
                     "median_width_burden_only_top25": float((top["rank_mc_burden_p95"] - top["rank_mc_burden_p05"]).median()),
                     "median_width_weights_only_top25": float((top["rank_mc_weights_p95"] - top["rank_mc_weights_p05"]).median()),
                     "median_p95_full_top10": float(top10["rank_mc_p95"].median()),
                     "n_top10_with_p_top10_ge_0_5": int((top10["p_top10_mc"] >= 0.5).sum()),
                     "mean_p_top10_of_top10": float(top10["p_top10_mc"].mean()),
                     "mc_n_draws": int(g["mc_n_draws"].iloc[0]), "tau_inherited_pp": float(g["mc_tau_inherited"].iloc[0])})
    _w(pd.DataFrame(rows), "test7_monte_carlo_decomposition")


# ------------------------------------------------------------------------------------------------------------------
# figures (drafts)
# ------------------------------------------------------------------------------------------------------------------

def draft_figures() -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    D = FIGURES / "drafts"
    D.mkdir(parents=True, exist_ok=True)
    out = []
    df = read_table(O.TABLE)
    p = df[(df["condition_id"] == PKEY[0]) & (df["measurement_id"] == PKEY[1]) & (df["geo_level"] == PKEY[2])]
    top = p[p["rank_equal"] <= 25].sort_values("rank_equal")
    fig, axes = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw={"width_ratios": [1.2, 1]})
    comps = ["burden_pct", "vulnerability_pct", "diagnostic_desert_pct", "clinic_capacity_pct", "research_readiness_pct"]
    im = axes[0].imshow(top[comps].to_numpy(float), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    axes[0].set_yticks(range(len(top)), [f"{int(r)}. {n} ({s})" for r, n, s in
                                         zip(top["rank_equal"], top["geo_name"].str.replace(" County", ""),
                                             top["state_abbr"])], fontsize=7)
    axes[0].set_xticks(range(len(comps)), [c.replace("_pct", "").replace("_", "\n") for c in comps], fontsize=8)
    fig.colorbar(im, ax=axes[0], shrink=0.6, label="percentile (normalised component)")
    y = np.arange(len(top))
    axes[1].errorbar(top["rank_mc_p50"], y, xerr=[top["rank_mc_p50"] - top["rank_mc_p05"],
                                                   top["rank_mc_p95"] - top["rank_mc_p50"]],
                     fmt="o", color="#1f5fa8", capsize=2, label="Monte Carlo median, 5th-95th pct")
    axes[1].plot(top["rank_equal"], y, "D", color="black", ms=4, label="observed rank (equal weights)")
    axes[1].plot(top["rank_burden_only"], y, "x", color="#e31a1c", label="burden-only rank")
    axes[1].set_xscale("log")
    axes[1].invert_yaxis()
    axes[1].set_yticks([])
    axes[1].set_xlabel("rank (log scale)")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)
    fig.suptitle("DRAFT - Top-25 counties: Long COVID or ME/CFS x wearable autonomic/activity monitoring. Components "
                 "(left) and rank uncertainty (right)", fontsize=10)
    fig.tight_layout()
    f = D / "scoring_primary_top25_components.png"
    fig.savefig(f, dpi=120)
    plt.close(fig)
    out.append(str(f))
    if "test7_region_stability_primary" in T:
        r = T["test7_region_stability_primary"]
        r = r[r["rank_equal"] <= 25]
        names = [c for c in r.columns if c.startswith("rank_") and c.replace("rank_", "") in O.weight_sets()]
        fig, ax = plt.subplots(figsize=(10, 8))
        yy = np.arange(len(r))
        for c in names:
            ax.plot(r[c], yy, "o", ms=4, label=c.replace("rank_", ""), alpha=0.8)
        ax.set_xscale("log")
        ax.set_yticks(yy, [f"{int(a)}. {b}" for a, b in zip(r["rank_equal"], r["geo_name"])], fontsize=7)
        ax.invert_yaxis()
        ax.axvline(25, color="grey", ls="--", lw=0.8)
        ax.legend(fontsize=7, loc="lower right")
        ax.set_xlabel("rank under each named weight set (log)")
        ax.set_title("DRAFT - Test 7: ranks of the equal-weight top-25 under every named weight set", fontsize=10)
        fig.tight_layout()
        f = D / "scoring_test7_rank_by_weight_set.png"
        fig.savefig(f, dpi=120)
        plt.close(fig)
        out.append(str(f))
    if "test6_scoring_label_shuffle_draws" in T:
        d = T["test6_scoring_label_shuffle_draws"]
        s = T["test6_scoring_label_shuffle"]
        d = d[(d["condition_id"] == PKEY[0]) & (d["geo_level"] == "county")]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        for ax, (what, col) in zip(axes[:2], (("burden", "S_rho_with_real_burden_pct"),
                                              ("vulnerability", "S_rho_with_real_vulnerability_pct"))):
            for scheme, color in (("national", "#9aa0a6"), ("within_state", "#fdbf6f")):
                x = d[(d["shuffled"] == what) & (d["scheme"] == scheme)][col]
                ax.hist(x, bins=30, color=color, alpha=0.8, label=f"{scheme} shuffle")
            so = s[(s["condition_id"] == PKEY[0]) & (s["geo_level"] == "county") & (s["shuffled"] == what)]
            ax.axvline(float(so["S_observed"].iloc[0]), color="black", lw=2, label="observed")
            ax.set_title(f"Spearman(composite, real {what} percentile)", fontsize=9)
            ax.legend(fontsize=7)
        if "test6_scoring_clinic_shuffle_draws" in T:
            cd = T["test6_scoring_clinic_shuffle_draws"]
            cd = cd[cd["measurement_id"] == PKEY[1]]
            for scheme, color in (("national", "#9aa0a6"), ("within_state", "#fdbf6f")):
                axes[2].hist(cd[cd["scheme"] == scheme]["rho_capacity_shuffled_vs_observed"], bins=30, color=color,
                             alpha=0.8, label=f"clinic_capacity, {scheme}")
                axes[2].hist(cd[cd["scheme"] == scheme]["rho_readiness_shuffled_vs_observed"], bins=30,
                             color="#1f78b4" if scheme == "national" else "#a6cee3", alpha=0.8,
                             label=f"research_readiness, {scheme}")
            axes[2].set_title("clinic-location shuffle: Spearman(shuffled, observed component)", fontsize=9)
            axes[2].legend(fontsize=7)
        fig.suptitle("DRAFT - Test 6 (scoring) negative controls, Long COVID or ME/CFS x wearable monitoring, county",
                     fontsize=10)
        fig.tight_layout()
        f = D / "scoring_test6_controls.png"
        fig.savefig(f, dpi=120)
        plt.close(fig)
        out.append(str(f))
    return out


# ------------------------------------------------------------------------------------------------------------------
# markdown
# ------------------------------------------------------------------------------------------------------------------

def _state_abbr() -> dict:
    from ..geography.crosswalk import STATE_ABBR_TO_FIPS
    return STATE_ABBR_TO_FIPS


# 1:1 county renames with a new FIPS since 2010 that CMS MMD still publishes under the legacy code; bridged at CMS
# ingestion (measure_it.ingestion.cms_mmd.FIPS_RENAMES, the single list). Mergers and splits are not 1:1 and not
# bridged.
def bridged_fips_renames(condition_id: str = "me_cfs") -> list[dict]:
    """The 1:1 FIPS renames bridged at CMS ingestion and the primary-measure rows each 2024 county now carries
    (fips_as_published = the legacy code). Descriptive: it reports the bridge, it does not apply one."""
    from ..ingestion.cms_mmd import FIPS_RENAMES
    cms = read_table("geo_condition_burden__cms_mmd", columns=["geo_id", "fips_as_published", "condition_id", "year",
                                                              "adjustment", "stratum_age"])
    b = read_table("geo_condition_burden", columns=["geo_id", "condition_id", "measure_id", "is_primary_measure"])
    prim = b[(b["condition_id"] == condition_id) & b["is_primary_measure"].fillna(False).astype(bool)]
    names = O.spatial().levels["county"].geo.set_index("geo_id")["geo_name"]
    out = []
    for legacy, (new, what) in FIPS_RENAMES.items():
        rows = cms[(cms["geo_id"] == new) & (cms["fips_as_published"] == legacy)
                   & (cms["condition_id"] == condition_id)]
        yrs = rows["year"].dropna().astype(int)
        out.append({"geo_id": new, "geo_name": names.get(new, new), "legacy_fips": legacy, "rename": what,
                    "n_cms_rows_bridged": int(len(rows)),
                    "years": f"{yrs.min()}-{yrs.max()}" if len(yrs) else UNKNOWN,
                    "has_primary_value": bool((prim["geo_id"] == new).any()),
                    "legacy_code_left_in_burden": bool((b["geo_id"] == legacy).any())})
    return out


def _fmt(v, nd=3):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}"
    if isinstance(v, (float, np.floating)):
        if float(v).is_integer() and abs(v) >= 1:
            return f"{int(v):,}"
        return f"{v:.{nd}f}"
    return str(v)


def md_table(df: pd.DataFrame, cols: list[str] | None = None, nd: int = 3) -> str:
    cols = cols or list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df[cols].itertuples(index=False):
        lines.append("| " + " | ".join(_fmt(v, nd) for v in r) + " |")
    return "\n".join(lines)


def _p(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).to_numpy()
    return df[m]


def write_markdown(demo: dict, figs: list[str], maps: list[str]) -> str:
    df = read_table(O.TABLE)
    meta = json.loads((O.PROCESSED / f"{O.TABLE}.meta.json").read_text())
    ws = O.weight_sets()
    prim = _p(df, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2])
    L = []
    A = L.append
    A("# Scoring results: deployment opportunities, recommendations, Tests 5-7")
    A("")
    A(f"_Generated by `measure_it.scoring.report` at {utc_now_iso()} from tables it wrote in the same run "
      f"(`results/tables/test5_*`, `test6_scoring_*`, `test7_*`, `scoring_*`, `demo_*`; "
      f"`data/processed/deployment_opportunities` "
      f"built {meta['created_at']}). Plan (pre-specified): `docs/ANALYSIS_PLAN_SCORING.md`; method: `docs/SCORING.md`._")
    A("")
    A("Every ranked region below is a **candidate deployment opportunity** for pilot evaluation. The ranking is "
      "ecological (places, facilities, registries), rests on curated assumptions and proxy or inherited burden, and "
      "is not a validated diagnostic pathway.")
    A("")
    # ---- headline
    t5 = T["test5_contrasts"]
    t5p = _p(t5, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2])
    bo25 = _p(t5p, contrast="burden_only", top_n=25).iloc[0]
    bo10 = _p(t5p, contrast="burden_only", top_n=10).iloc[0]
    nd25 = _p(t5p, contrast="need_plus_desert", top_n=25).iloc[0]
    t7 = T["test7_summary"]
    t7p = _p(t7, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2]).iloc[0]
    mcd = T["test7_monte_carlo_decomposition"]
    mcp = _p(mcd, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2]).iloc[0]
    ls = T["test6_scoring_label_shuffle"]
    cs = T["test6_scoring_clinic_shuffle"]
    csp = _p(cs, measurement_id=PKEY[1], scheme="national")
    cap = csp[csp["component"] == "clinic_capacity_pct"].iloc[0]
    rdy = csp[csp["component"] == "research_readiness_pct"].iloc[0]
    A("## 0. Headline (primary query: Long COVID or ME/CFS x wearable autonomic/activity monitoring x county)")
    A("")
    A(f"* **Test 5, positive (partly by construction): facility and research information changes where the engine "
      f"would deploy.** Burden-only "
      f"vs full ranking: top-10 overlap {_fmt(bo10['overlap'])} of 10, top-25 overlap {_fmt(bo25['overlap'])} of 25 "
      f"(Jaccard {_fmt(bo25['jaccard'])}; tie-averaged overlap {_fmt(bo25['overlap_tie_averaged'], 1)}), Kendall "
      f"tau-b {_fmt(bo25['kendall_tau_b'])} over {_fmt(bo25['n_rank_corr'])} eligible counties. Adding "
      f"clinic_capacity + research_readiness to burden + vulnerability + desert alone keeps "
      f"{_fmt(nd25['overlap'])} of the top-25 (Jaccard {_fmt(nd25['jaccard'])}). Pre-specified reading "
      f"(Jaccard < 0.5): material change = {_fmt(bool(bo25['material_change_jaccard_lt_0_5']))}. The part of the change "
      f"that follows real locations is research readiness (where trial- and NIH-active facilities are); clinic "
      f"capacity is facility density and does not respond to the location shuffle (Test 6b).")
    A(f"* **The short list is not stable.** Monte Carlo (burden CI x evidence multiplier, "
      f"x{_fmt(float(prim['mc_inherited_multiplier'].iloc[0]))} for a state value inherited by a county where one is "
      f"still used, inherited-burden deviation, "
      f"Dirichlet weights; {_fmt(mcp['mc_n_draws'])} draws): median 5th-95th percentile rank-interval width of the "
      f"top-25 = {_fmt(mcp['median_width_full_top25'], 0)} ranks; only {_fmt(mcp['n_top10_with_p_top10_ge_0_5'])} of "
      f"the top-10 are in the top-10 in at least half of the draws. Test 7: Kendall's W across the "
      f"{len(ws)} named weight sets = {_fmt(t7p['kendalls_w_named'])} "
      f"(without burden_only {_fmt(t7p['kendalls_w_named_excl_burden_only'])}), across 1,000 random weight vectors "
      f"= {_fmt(t7p['kendalls_w_random'])}; {_fmt(t7p['n_default_top25_unstable'])} of the default top-25 meet the "
      f"pre-specified instability rule. The overall ordering is concordant; which counties make the top 10-25 "
      f"depends on the weights.")
    rho_cc = float(cap["component_rho_shuffled_vs_observed_mean"])
    if "clinic_capacity_basis" in prim.columns and prim["clinic_capacity_basis"].iloc[0] == "measurement_activity":
        A(f"* **Test 6(b): clinic capacity now follows where the measurement is performed** (Spearman of shuffled "
          f"with observed clinic_capacity {_fmt(rho_cc)} on average, national shuffle; it was 0.94-0.97 for the "
          f"implementer-density capacity used before 2026-10-07). Since 2026-10-07 the primary query's clinic_capacity "
          f"counts clinicians who billed the measurement's dedicated codes to Medicare FFS beneficiaries in 2024 "
          f"(facilities.activity), and the shuffle permutes those clinicians' locations. Most counties have few or no "
          f"such clinicians, and sparsity alone also lets a shuffle move a measure. Research readiness also moves "
          f"(Spearman {_fmt(rdy['component_rho_shuffled_vs_observed_mean'])}).")
    else:
        A(f"* **Test 6(b), negative for clinic capacity: the clinic-location shuffle barely moves clinic_capacity** "
          f"(Spearman of shuffled with observed clinic_capacity {_fmt(rho_cc)} "
          f"on average, national shuffle), because the shuffle keeps the set of occupied locations and wearable "
          f"implementer groups include primary care, which is most facilities. Clinic capacity for this measurement "
          f"measures facility/provider density, not measurement-specific capability. Research readiness does move "
          f"(Spearman {_fmt(rdy['component_rho_shuffled_vs_observed_mean'])}).")
    lsb = _p(ls, condition_id=PKEY[0], geo_level="county", shuffled="burden", scheme="national").iloc[0]
    A(f"* **Test 6(a): shuffling burden labels removes the ranking's burden structure** (Spearman of the composite "
      f"with real burden {_fmt(lsb['S_observed'])} observed vs {_fmt(lsb['S_null_mean'])} "
      f"[{_fmt(lsb['S_null_p2_5'])}, {_fmt(lsb['S_null_p97_5'])}] shuffled; top-25 overlap with the observed list "
      f"{_fmt(lsb['top25_overlap_mean'], 1)} on average). Long-COVID county burden is the inherited state value: a "
      f"within-state shuffle of it changes nothing (section 5).")
    A("* **Burden is weak everywhere it matters.** Long-COVID county burden is the state HPS value inherited by every "
      "county (no within-state information); ME/CFS burden is a level-C claims proxy (CMS 'Fibromyalgia, Chronic Pain "
      "and Fatigue'); POTS and dysautonomia have no usable burden (level D), so their rankings contain no burden "
      "information at all.")
    bra = T["scoring_set_burden_rule_before_after"]
    brp = _p(bra, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2]).iloc[0]
    brt = _p(T["scoring_set_burden_rule_top10"], condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2])
    entered = "; ".join(f"{r.geo_name} ({r.rank_after:.0f})" for r in brt[brt["status"] == "entered"].itertuples())
    A(f"* **Incomplete burden is no longer ranked (rule changed 2026-09-24, plan deviation 7).** A set's burden is "
      f"computed only where every member with a defined burden measure has a value. The "
      f"{_fmt(brp['n_incomplete'])} counties where ME/CFS has no CMS value "
      f"({_fmt(brp['n_incomplete_population_eligible'])} with population >= {O.min_population():,}, among them the 9 "
      f"Connecticut planning regions) used to be ranked in the two sets on the inherited long-COVID state value alone "
      f"(and for ME/CFS alone on the other components without burden); they are now listed separately, not ranked "
      f"(section 5). Primary query before/after: top-10 overlap "
      f"{_fmt(brp['top10_overlap_before_after'])} of 10 (left: {brp['incomplete_in_top10_before'] or 'none'}; "
      f"entered: {entered or 'none'}), top-25 overlap {_fmt(brp['top25_overlap_before_after'])} of 25, Kendall tau-b "
      f"{_fmt(brp['kendall_tau_b_common'])} over the {_fmt(brp['n_common'])} counties ranked under both rules.")
    s6 = T["test7_sensitivity_specialist_only"]
    s6p = _p(s6, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2], variant="specialist_only").iloc[0]
    qa = T["test7_sensitivity_specialist_only_qa"]
    shares = qa[qa["condition_id"].isin(O.condition_spec(PKEY[0])["members"]) & (qa["geo_level"] == "county")
                & qa["share_primary_care"].notna()]
    A(f"* **S6 specialist_only (declared 2026-09-24 before it was computed): counting relevant providers without "
      f"primary care keeps {_fmt(s6p['top10_overlap'])} of the top-10 and {_fmt(s6p['top25_overlap'])} of the "
      f"top-25** (Kendall tau-b {_fmt(s6p['kendall_tau_b'])} over {_fmt(s6p['n_compared'])} counties). Primary care is "
      f"{100 * shares['share_primary_care'].min():.0f}-{100 * shares['share_primary_care'].max():.0f}% of the relevant "
      "NPIs of the query's conditions, so the desert's provider term and clinic capacity mostly measure primary-care "
      "supply; which counties make the top-10 depends on that choice (section 8).")
    A(metric_link_headline())
    A("")
    # ---- build
    A("## 1. What was built")
    A("")
    n_combo = df.groupby(["condition_id", "measurement_id", "geo_level"]).ngroups
    A(f"`data/processed/deployment_opportunities.parquet`: {_fmt(len(df))} rows, {df.shape[1]} columns, {n_combo} "
      f"condition (or set) x measurement x level combinations; `object_id = opportunity:<condition>|<measurement>|<geo_id>`. "
      f"Conditions: {', '.join(O.SINGLE_CONDITIONS)}; sets: {', '.join(O.condition_sets())}. Measurements: "
      f"{', '.join(O.MEASUREMENTS)}. Levels: county (3,144 counties of the 50 states + DC; "
      f"{_fmt(int((~prim['small_population_flag'].astype(bool)).sum()))} with population >= {O.min_population():,}, of "
      f"which {_fmt(int(prim['rank_eligible'].sum()))} are ranked for the primary query; the others have an incomplete "
      f"burden) and state (51).")
    A("")
    lev = df[df["geo_level"] == "county"].groupby(["condition_id"]).agg(
        rows=("geo_id", "size"), level_A=("burden_evidence_level", lambda s: int((s == "A").sum())),
        level_C=("burden_evidence_level", lambda s: int((s == "C").sum())),
        level_D=("burden_evidence_level", lambda s: int((s == "D").sum())),
        inherited=("burden_inherited", "sum"), burden_excluded_D=("burden_excluded_level_D", "sum"),
        incomplete=("burden_incomplete", "sum")).reset_index()
    lev = lev.assign(**{c: lev[c] // len(O.MEASUREMENTS) for c in ("rows", "level_A", "level_C", "level_D", "inherited",
                                                                   "burden_excluded_D", "incomplete")})
    A("County rows per condition by burden evidence level (per measurement; a set's level is its least direct "
      "contributing member; `burden_excluded_D` = no member has a defined burden, so the burden weight is removed and "
      "the other weights renormalised; `incomplete` = a member with a defined burden has no value in the county, so "
      "there is no burden, composite or rank, level D):")
    A("")
    A(md_table(lev))
    A("")
    pinc = prim[prim["burden_incomplete"].astype(bool)]
    A(f"Incomplete burden (plan deviation 7): {len(pinc)} counties have no CMS value for ME/CFS "
      f"({int((~pinc['small_population_flag']).sum())} with population >= {O.min_population():,}); for me_cfs, "
      "long_covid_or_me_cfs and the demo cluster they are scored on the other components but not ranked, and are "
      "listed in `results/tables/scoring_incomplete_burden.csv` with the reason (section 5). Incomplete state rows "
      f"(all combinations): {int(df.loc[df['geo_level'] == 'state', 'burden_incomplete'].astype(bool).sum())}.")
    A("")
    A(f"Inherited-burden deviation used in the Monte Carlo: tau = {_fmt(float(prim['mc_tau_inherited'].iloc[0]))} "
      "percentage points (method-of-moments between-state SD of the HPS state estimates). The state estimate's own "
      f"error is inflated x{_fmt(float(prim['mc_inherited_multiplier'].iloc[0]))} when it is carried to a county "
      "(`inherited_burden_uncertainty_multiplier` in configs/scoring.yaml: a state value has no within-state "
      "information, so as a county burden it is treated like a level-C proxy; added 2026-09-23, sensitivity S1b "
      "in section 8).")
    A("")
    # ---- formula
    A("## 2. Formula and weights")
    A("")
    A("`composite_ws(g,c,m) = sum_k w_k * pr(component_k)` over burden, vulnerability, diagnostic_desert, "
      "clinic_capacity, research_readiness (weighted mean over available components; level-D burden removed and the "
      "other weights renormalised). `pr()` = percentile rank within the 50 states + DC. Components: burden = primary "
      "measure of `geo_condition_features` (condition set: pr of the mean of the burden percentiles of the members "
      "with a defined burden measure, computed only where every such member has a value (plan deviation 7, "
      "2026-09-24) and re-normalised like every other averaged component (reviewer correction 2026-09-23)); "
      "vulnerability = SVI overall; diagnostic_desert = the geography module's "
      "index (access-only for level D); clinic_capacity = mean pr of implementer-group individual providers and "
      "implementer facilities per 100k, in the county and within 50 km; research_readiness = mean pr of distinct "
      "condition trials, NIH core projects and technology-experience trials at facilities in the county or within "
      "50 km. technology_saturation is reported and used only in `saturation_adjusted`.")
    A("")
    wt = pd.DataFrame([{"weight_set": k, **{c: v.get(c, 0.0) for c in O.COMPONENTS + ["low_saturation"] +
                                            O.EVIDENCE_COMPONENTS}} for k, v in ws.items()])
    A(md_table(wt, nd=3))
    A("")
    A("`evidence_weighted` (pre-declared 2026-09-25 as the primary alternative; `equal` stays the default for "
      "continuity) adds two components on an absolute 0-1 scale: `measurement_evidence` (tier factor x "
      "clip((AUROC lower 95% bound - 0.5) / 0.5, 0, 1) of the scored performance record) and `expected_yield` (burden "
      "percentile x sensitivity at 0.90 specificity x implementation reach). A combination whose measurement "
      "performance is UNKNOWN or partial has no evidence_weighted composite or rank (section 12; "
      "docs/ANALYSIS_PLAN_METRIC_LINK.md).")
    A("")
    # ---- top regions
    A("## 3. Top regions (primary query, equal weights)")
    A("")
    top = prim[prim["rank_equal"] <= 15].sort_values("rank_equal")
    cols = ["rank_equal", "geo_name", "state_abbr", "composite_equal", "burden_pct", "burden_evidence_level",
            "vulnerability_pct", "diagnostic_desert_pct", "clinic_capacity_pct", "research_readiness_pct",
            "rank_mc_p05", "rank_mc_p95", "p_top10_mc", "rank_burden_only"]
    A(md_table(top, cols, nd=2))
    A("")
    n_inh = int(top["burden_inherited"].sum())
    n_me = int(top["burden_members_contributing"].str.contains("me_cfs").sum())
    A(f"{n_inh} of these {len(top)} counties carry inherited long-COVID burden (the state HPS estimate; see "
      f"`uncertainties` in the table and the JSON); in {n_me} of them the ME/CFS level-C claims proxy also contributes"
      + (f" (the other {len(top) - n_me}: "
         + "; ".join(top.loc[~top["burden_members_contributing"].str.contains("me_cfs"), "geo_name"]) + ")"
         if n_me < len(top) else " (every ranked county has both members since the 2026-09-24 completeness rule)") +
      ". `rank_mc_p05/p95` is the 5th-95th percentile rank over the Monte Carlo draws.")
    A("")
    zc = T["test5_readiness_zero_counts"]
    zc10 = _p(zc, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2], subset="equal top-10").iloc[0]
    zca = _p(zc, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2], subset="all eligible").iloc[0]
    # the research-readiness floor: the tie of counties with no condition trial, NIH core project or
    # technology-experience trial in reach (reviewer correction 2026-09-24: the sentence used to say the zero-count
    # top-10 counties were 'carried by technology-experience trials' even when they sat at this floor)
    rr_floor = float(prim.loc[prim["rank_eligible"].astype(bool), "research_readiness_pct"].min())
    z_lo, z_hi = zc10["research_readiness_pct_min_among_zero"], zc10["research_readiness_pct_max_among_zero"]
    n_z10 = int(zc10["n_zero_condition_trials_and_nih"])
    if n_z10 == 0:
        z_txt = "In the top-10, every county has a Long COVID or ME/CFS trial or NIH core project in reach."
    elif abs(z_hi - rr_floor) < 1e-9:
        z_txt = (f"In the top-10, {n_z10} counties have neither and sit at the research-readiness floor "
                 f"({_fmt(rr_floor, 2)}, the tie of counties with no condition trial, no NIH core project and no "
                 "technology-experience trial in reach): they reach the top-10 on the other components.")
    else:
        z_txt = (f"In the top-10, {n_z10} counties have neither; their research_readiness percentile is "
                 f"{_fmt(z_lo, 2)}-{_fmt(z_hi, 2)} (floor {_fmt(rr_floor, 2)}): above the floor it is carried by "
                 "technology-experience trials of any target condition.")
    A(f"**Research readiness is compressed by zeros.** {_fmt(zca['n_zero_condition_trials_and_nih'])} of "
      f"{_fmt(zca['n_regions'])} eligible counties ({_fmt(zca['share_zero'])}) have no Long COVID or ME/CFS trial and "
      f"no NIH core project at any facility in the county or within 50 km, and all of them tie at the bottom of those "
      f"two sub-percentiles. {z_txt} A high research_readiness percentile does not by itself mean condition-specific "
      "study experience (technology-experience trials of any target condition count too); the recommendation JSON "
      "lists the trial and NIH ids so the reader can see which it is.")
    A("")
    ptop_state = _p(df, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level="state")
    A("State level (same query):")
    A("")
    A(md_table(ptop_state[ptop_state["rank_equal"] <= 10].sort_values("rank_equal"),
               ["rank_equal", "geo_name", "composite_equal", "burden_pct", "vulnerability_pct", "diagnostic_desert_pct",
                "clinic_capacity_pct", "research_readiness_pct", "rank_mc_p05", "rank_mc_p95", "rank_burden_only"], nd=2))
    A("")
    A("Monte Carlo decomposition (median rank-interval width of each combination's top-25; `burden_only` = burden "
      "noise with default weights, `weights_only` = Dirichlet weights with observed burden):")
    A("")
    A(md_table(mcd[mcd["geo_level"] == "county"].sort_values(["measurement_id", "condition_id"]),
               ["condition_id", "measurement_id", "median_width_full_top25", "median_width_burden_only_top25",
                "median_width_weights_only_top25", "n_top10_with_p_top10_ge_0_5", "mean_p_top10_of_top10"], nd=2))
    A("")
    # ---- demo
    A("## 4. Demo recommendation (SPEC demo query)")
    A("")
    A("The same query under the primary alternative `evidence_weighted` is in "
      "`results/example_deployment_recommendation_evidence_weighted.json` and "
      "`results/tables/demo_top10_regions_evidence_weighted.csv` (section 12). Every recommendation's technology block "
      "now carries the measurement's performance record and, per region, the expected yield and false positives.")
    A("")
    A(f"> {REC.DEMO_QUERY}")
    A("")
    A(f"`rank_deployment_opportunities{tuple(REC.DEMO_ARGS.values())}` resolved the condition to "
      f"`{demo['condition_resolution']['condition_id']}` ({demo['condition_resolution']['match_reason']}) and the "
      f"measurement to `{demo['measurement_resolution']['measurement_id']}`. Output: "
      "`results/example_deployment_recommendation.json` (SPEC schema per region), `data/processed/deployment_candidates`, "
      "`results/tables/demo_top10_regions.csv`, `demo_top10_candidate_sites.csv`, `demo_top10_research_evidence.csv`.")
    A("")
    dr = T.get("demo_top10_regions")
    if dr is not None and len(dr):
        A(md_table(dr, ["rank", "geo_name", "state", "rank_interval_p05", "rank_interval_p95", "p_top10_monte_carlo",
                        "burden_evidence_level", "condition_trials_in_pool", "condition_nih_core_projects_in_pool",
                        "technology_experience_trials_in_pool", "n_candidate_sites"], nd=2))
        A("")
        dsites = T["demo_top10_candidate_sites"]
        n_research = int(dr["condition_trials_in_pool"].fillna(0).gt(0).sum())
        exl = demo.get("excluded_incomplete_burden") or []
        if exl:
            A(f"Not ranked (incomplete burden; the JSON's `excluded_incomplete_burden` block, with the reason and the "
              f"kept components): {'; '.join(e['name'] for e in exl)}.")
            A("")
        A(f"Candidate sites returned: {_fmt(len(dsites))} facility rows over {len(dr)} regions (median "
          f"{_fmt(float(dr['n_candidate_sites'].median()), 0)} per region). {n_research} of the {len(dr)} regions have "
          "at least one Long COVID or ME/CFS trial registered at a facility in the county or within 50 km; in the "
          "others the candidate list rests on specialty match, community access and distance, and the JSON carries "
          "the matcher's absence note (a possible measurement/research desert).")
        A("")
    # ---- Test 5
    A("## 5. Test 5: does facility/research information change where we would deploy?")
    A("")
    A("Primary query (eligible counties):")
    A("")
    A(md_table(t5p, ["contrast", "top_n", "overlap", "jaccard", "overlap_tie_averaged", "n_tied_at_boundary_contrast",
                     "expected_overlap_independent", "kendall_tau_b", "spearman", "material_change_jaccard_lt_0_5"]))
    A("")
    cc = T["test5_component_correlations"]
    ccp = _p(cc, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2])
    A("Spearman correlations between the normalised components (primary query, eligible counties; "
      "`test5_component_correlations.csv` has every combination). The desert's access terms oppose clinic_capacity "
      "and research_readiness, so the full composite partly trades them off:")
    A("")
    A(md_table(ccp, ["component_a", "component_b", "spearman", "n"]))
    A("")
    rs = ccp[(ccp["component_a"] == "research_readiness_pct") & (ccp["component_b"] == "technology_saturation_pct")]
    if len(rs):
        A(f"research_readiness and technology_saturation correlate {_fmt(float(rs['spearman'].iloc[0]))} (both count "
          "technology-experience trials), so the `saturation_adjusted` weight set, which enters (1 - saturation), "
          "largely cancels research readiness rather than adding separate information.")
        A("")
    mv = T["test5_movements"]
    mvp = _p(mv, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2], contrast="burden_only", top_n=10)
    A("Regions entering / leaving the top-10 (burden-only vs full), with the component contributions "
      "`w_k (x_k - mean_k)` as the reason:")
    A("")
    A(md_table(mvp, ["movement", "geo_name", "state_abbr", "rank_full", "rank_contrast", "reason"], nd=3))
    A("")
    agg = t5[(t5["top_n"] == 25)].groupby(["contrast", "geo_level"]).agg(
        combos=("jaccard", "size"), median_jaccard=("jaccard", "median"), min_jaccard=("jaccard", "min"),
        max_jaccard=("jaccard", "max"), n_material_change=("material_change_jaccard_lt_0_5", "sum"),
        median_tau=("kendall_tau_b", "median")).reset_index()
    A("All combinations (top-25):")
    A("")
    A(md_table(agg))
    A("")
    csw = _p(T["test6_scoring_clinic_shuffle"], measurement_id=PKEY[1], scheme="national",
             component="clinic_capacity_pct").iloc[0]
    rsw = _p(T["test6_scoring_clinic_shuffle"], measurement_id=PKEY[1], scheme="national",
             component="research_readiness_pct").iloc[0]
    n_mat = int(_p(t5, geo_level="county", contrast="burden_only", top_n=25)["material_change_jaccard_lt_0_5"].sum())
    n_cc = len(_p(t5, geo_level="county", contrast="burden_only", top_n=25))
    A(f"Reading. The pre-specified criterion (top-25 Jaccard < 0.5) is met in {n_mat} of {n_cc} county "
      f"combinations: adding clinic_capacity and research_readiness (40% of the default weight) re-orders the short "
      f"list. That is partly by construction (any two added components with 40% weight would move a list), so the "
      f"useful questions are whether the change follows where facilities actually are, and whether the new list is "
      f"stable. On the first, the answer is narrower than 'facility geography': placing the same facilities and "
      f"providers at shuffled locations keeps {_fmt(csw['top25_overlap_mean'], 1)} of the observed top-25 on average "
      f"(Test 6b), but that shuffle leaves clinic_capacity almost unchanged (Spearman "
      f"{_fmt(csw['component_rho_shuffled_vs_observed_mean'])} with its observed value, because it keeps the number of "
      f"facilities per location) and holds the desert's provider/trial terms at their observed values. What it does "
      f"move is research_readiness (Spearman {_fmt(rsw['component_rho_shuffled_vs_observed_mean'])}). So the change "
      f"tracks where trial- and NIH-active facilities are; the test says nothing about whether clinic_capacity "
      f"locates measurement-capable clinics. On the second: the list is not stable (Test 7 and the Monte Carlo "
      f"intervals in section 3).")
    A("")
    # incomplete burden (plan deviation 7): regions excluded from the ranking, before/after of the rule change
    bra = T["scoring_set_burden_rule_before_after"]
    brt = T["scoring_set_burden_rule_top10"]
    aff = bra[bra["n_incomplete"] > 0]
    pinc = prim[prim["burden_incomplete"].astype(bool) & ~prim["small_population_flag"].astype(bool)]
    A("**Incomplete burden (rule changed 2026-09-24; plan deviation 7).** Until 2026-09-23 a set's burden was the mean "
      "of the member percentiles over the members with a value in the region, so where ME/CFS has no CMS value the "
      "set burden rested on the inherited long-COVID state value alone, a different quantity ranked as if it were the "
      "same (in the build of 2026-09-23 it put Greater Bridgeport Planning Region CT in the top-10 and Oglala Lakota "
      "County SD, then without a CMS value, at burden-only rank 1). Now the burden is computed only where every member with a defined burden measure has a value; the other "
      "regions keep their non-burden components, get no burden, composite or rank, and are listed separately. The "
      "same rule applies to a single condition (a set of one): me_cfs's counties without a CMS value were ranked on "
      "the other four components with the burden weight renormalised and are now incomplete too. "
      f"{len(pinc)} population-eligible counties are incomplete for the primary query:")
    A("")
    A(md_table(pinc.sort_values("geo_id"), ["geo_name", "state_abbr", "population_total", "burden_members_missing",
                                             "burden_missing_source_reason"]))
    A("")
    # 1:1 FIPS renames (reviewer note 2026-09-24): bridged at CMS ingestion since 2026-09-24, so not incomplete
    br = bridged_fips_renames()
    pr = prim.set_index("geo_id")
    parts = []
    for r in br:
        g = r["geo_id"]
        rk = pr["rank_equal"].get(g) if g in pr.index else None
        bo = pr["rank_burden_only"].get(g) if g in pr.index and "rank_burden_only" in pr else None
        state = ("still incomplete" if g in set(pinc["geo_id"]) else
                 f"now complete (equal-weight rank {_fmt(rk, 0)}, burden-only rank {_fmt(bo, 0)} for the primary "
                 "query)" if rk is not None and rk == rk else "complete but not ranked (small population)")
        parts.append(f"{r['geo_name']} ({g}): {r['rename']}; CMS MMD still publishes it under {r['legacy_fips']} "
                     f"({r['n_cms_rows_bridged']} ME/CFS-mapped rows, {r['years']}), bridged to {g} at ingestion "
                     f"(`fips_as_published = {r['legacy_fips']}`); {state}")
    A("**1:1 FIPS renames are bridged, not gaps.** " + "; ".join(parts) + ". Until 2026-09-24 the CMS ingestion left "
      "these rows on the legacy codes, so both counties were listed as incomplete with a 'no published value' reason; "
      "the bridge is a code identity (same boundary), unlike the Connecticut planning regions, whose values CMS "
      "publishes only on the legacy counties and which are not apportioned. The regions listed above are absent from "
      "every MMD file ingested (reason per region in `burden_missing_source_reason`).")
    A("")
    A("Before/after (old rule `available` vs current rule `complete`, equal weights, counties ranked under both rules "
      "for tau; `results/tables/scoring_set_burden_rule_before_after.csv` has every combination, 0 changes wherever "
      "no region is incomplete):")
    A("")
    A(md_table(aff[aff["geo_level"] == "county"].sort_values(["measurement_id", "condition_id"]),
               ["condition_id", "measurement_id", "n_incomplete_population_eligible", "top10_overlap_before_after",
                "top25_overlap_before_after", "kendall_tau_b_common", "incomplete_in_top10_before",
                "incomplete_in_burden_only_top25_before"]))
    A("")
    for cid in (PKEY[0], "autonomic_activity_invisible_illness"):
        t = _p(brt, condition_id=cid, measurement_id=PKEY[1], geo_level=PKEY[2])
        if len(t):
            A(f"Top-10 before and after, {O.condition_spec(cid)['label']} x {PKEY[1]} (county):")
            A("")
            A(md_table(t, ["geo_name", "state_abbr", "rank_before", "rank_after", "rank_burden_only_before",
                           "rank_burden_only_after", "status"]))
            A("")
    other = aff[(aff["top10_overlap_before_after"] < 10) & (aff["incomplete_in_top10_before"].fillna("") == "")]
    A("The composites of the complete counties move only through the re-normalisation of the set burden over the "
      "complete counties (Kendall tau-b above). Where an incomplete county was in the old top-10, the change is its "
      "removal and the next county moving up"
      + (f"; in {len(other)} combination(s) ({'; '.join(f'{r.condition_id} x {r.measurement_id}' for r in other.itertuples())}) "
         "the re-normalisation alone swaps a county at the top-10 boundary" if len(other) else "")
      + ". These incomplete counties may still be deployment opportunities; this ranking cannot place them, and "
      "`rank_deployment_opportunities` returns them in `excluded_incomplete_burden`.")
    A("")
    lvd = t5[(t5["top_n"] == 25) & (t5["contrast"] == "burden_only") & t5["condition_id"].isin(["pots", "dysautonomia"])]
    A(f"POTS and dysautonomia have no burden-only ranking (level D): `n_eligible_contrast` = "
      f"{_fmt(int(lvd['n_eligible_contrast'].max()) if len(lvd) else 0)}. Their full rankings rest on vulnerability, "
      "the access-only desert, clinic capacity and research readiness only.")
    A("")
    lc = _p(t5, condition_id="long_covid", measurement_id=PKEY[1], geo_level="county", contrast="burden_only", top_n=25)
    if len(lc):
        A(f"Long COVID alone at county level: the burden-only ranking is the inherited state value, so "
          f"{_fmt(lc['n_tied_at_boundary_contrast'].iloc[0])} counties are tied at the top-25 boundary; the FIPS "
          f"tie-break overlap is {_fmt(lc['overlap'].iloc[0])} and the tie-averaged overlap "
          f"{_fmt(lc['overlap_tie_averaged'].iloc[0], 1)}. A county burden-only list for long COVID is a list of "
          "one state's counties, not a burden ranking.")
        A("")
    # ---- Test 6
    A("## 6. Test 6: negative controls")
    A("")
    A("(a) Geography-label shuffles (200 permutations each; S = Spearman of the composite with the REAL, unshuffled "
      "component; the ranking loses that component's structure if the observed S exceeds the 97.5th percentile of "
      "the shuffled S):")
    A("")
    A(md_table(ls, ["condition_id", "geo_level", "shuffled", "scheme", "n_draws", "S_observed", "S_null_mean",
                    "S_null_p97_5", "p_one_sided", "loses_structure", "rho_with_observed_mean", "top25_overlap_mean",
                    "ranking_unchanged_all_draws"]))
    A("")
    A("The S criterion is close to true by construction whenever a component carries weight; the informative part is "
      "how far the list moves (`rho_with_observed_mean`, `top25_overlap_mean`). Rows with "
      "`ranking_unchanged_all_draws = yes` are shuffles that cannot change anything: for long COVID alone, the "
      "within-state shuffle permutes identical inherited state values among a state's counties, which demonstrates "
      "that the county burden carries no within-state information.")
    A("")
    A("(b) Clinic-location shuffle (facility and individual-provider location tuples permuted; 200 national, 100 "
      "within-state draws; clinic_capacity, research_readiness and technology_saturation recomputed for every county; "
      "the diagnostic desert's provider and trial sub-percentiles come from the geography module and are held at "
      "their observed values, so this control tests the two scoring components only):")
    A("")
    A(md_table(cs, ["measurement_id", "scheme", "component", "n_draws", "S_observed", "S_null_mean", "S_null_p97_5",
                    "loses_structure", "component_rho_shuffled_vs_observed_mean", "rho_with_observed_mean",
                    "top25_overlap_mean", "pool_size_invariance_violations"]))
    A("")
    ccn = cs[(cs["scheme"] == "national") & (cs["component"] == "clinic_capacity_pct")].set_index("measurement_id")
    rrn = cs[(cs["scheme"] == "national") & (cs["component"] == "research_readiness_pct")]
    ccr = ccn["component_rho_shuffled_vs_observed_mean"]
    A("Reading: the permutation keeps the number of facilities (and of providers) at every location slot and moves "
      "only which facility sits there, so it cannot change a pure density measure by design. A component that stays "
      "highly correlated with its observed value therefore measures facility/provider density rather than where "
      "measurement-specific capability sits. Since 2026-10-07 the bundles with a dedicated billing code (wearable, "
      "autonomic testing, CPET) use activity capacity, whose clinicians' locations are permuted instead, and nailfold "
      "capillaroscopy keeps the density basis (mean Spearman shuffled vs observed "
      + ", ".join(f"{k} {_fmt(v)}" for k, v in ccr.items() if k != "exercise_capacity_testing")
      + f"; exercise capacity testing {_fmt(ccr.get('exercise_capacity_testing', np.nan))}). research_readiness, which counts condition trials, NIH "
      f"projects and technology-experience trials, moves for every bundle "
      f"({_fmt(rrn['component_rho_shuffled_vs_observed_mean'].min())}-"
      f"{_fmt(rrn['component_rho_shuffled_vs_observed_mean'].max())}). Where clinic_capacity stays near 1 it measures "
      "density (the reason it was replaced where a billing code exists), not a failure of the rest of the pipeline. The facilities "
      "contributor's matcher-level shuffle (`results/tables/test6_clinic_shuffle.csv`, docs/FACILITY_MATCHING.md "
      "section 5) is a separate test and is not re-stated here.")
    A("")
    ph = T["test6_scoring_phenotype_label_nulls"]
    A("(c) Phenotype-label permutations, restated from the wearable modules' own tables:")
    A("")
    A(md_table(ph, ["source_file", "target", "model_or_score", "control", "observed", "null_mean", "null_upper", "p",
                    "reading"], nd=3))
    A("")
    n_not = int((~ph["observed_above_null_p_le_0_05"]).sum())
    A(f"{n_not} of {len(ph)} phenotype-level controls are NOT distinguishable from their null (p > 0.05); the "
      "Long COVID wearable-only models of the Stanford/Uwakwe cohort are among them. The deployment layer does not "
      "use these person-level results to score regions; they appear in each recommendation's `phenotype` block as "
      "condition-level context.")
    A("")
    # ---- Test 7
    A("## 7. Test 7: ranking robustness")
    A("")
    A(md_table(t7[t7["geo_level"] == "county"].sort_values(["measurement_id", "condition_id"]),
               ["condition_id", "measurement_id", "kendalls_w_named", "kendalls_w_named_excl_burden_only",
                "kendalls_w_random", "pairwise_tau_median", "pairwise_tau_min_excl_burden_only",
                "n_regions_top25_any_named", "n_unstable", "n_default_top25_unstable",
                "median_named_range_default_top25"]))
    A("")
    A("State level:")
    A("")
    A(md_table(t7[t7["geo_level"] == "state"].sort_values(["measurement_id", "condition_id"]),
               ["condition_id", "measurement_id", "kendalls_w_named", "kendalls_w_random", "pairwise_tau_median",
                "n_regions_top25_any_named", "n_unstable", "n_default_top25_unstable"]))
    A("")
    un = T["test7_unstable_regions"]
    unp = _p(un, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2]).sort_values("rank_equal").head(20)
    A(f"Unstable regions, primary query (first 20 of {_fmt(len(_p(un, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2])))} "
      "by equal-weight rank; all combinations in `results/tables/test7_unstable_regions.csv`):")
    A("")
    A(md_table(unp, ["geo_name", "state_abbr", "rank_equal", "rank_named_min", "rank_named_max", "frac_named_top25",
                     "frac_random_top25", "unstable_reason"], nd=2))
    A("")
    pt = T["test7_pairwise_tau"]
    ptp = _p(pt, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2])
    A("Pairwise Kendall tau-b between named weight sets (primary query):")
    A("")
    A(md_table(ptp, ["weight_set_a", "weight_set_b", "kendall_tau_b", "top25_overlap"]))
    A("")
    # ---- sensitivity
    A("## 8. Sensitivity analyses")
    A("")
    s1 = T["test7_sensitivity_s1_tau"]
    A(f"S1, tau = 0 (pure inheritance) vs tau = {_fmt(float(s1['tau_pp'].iloc[0]))}: median top-25 interval width "
      f"{_fmt(float((s1['rank_mc_p95_with_tau'] - s1['rank_mc_p05_with_tau']).median()), 0)} with tau vs "
      f"{_fmt(float((s1['rank_mc_p95_tau0'] - s1['rank_mc_p05_tau0']).median()), 0)} without; mean P(top-25) of the "
      f"top-25 {_fmt(float(s1['p_top25_with_tau'].mean()))} vs {_fmt(float(s1['p_top25_tau0'].mean()))} "
      "(`test7_sensitivity_s1_tau.csv`).")
    A("")
    s1b = T["test7_sensitivity_s1b_inherited_multiplier"]
    A(f"S1b, inherited-burden multiplier {_fmt(float(s1b['inherited_multiplier_configured'].iloc[0]))} (configured; "
      f"not pre-specified) vs {_fmt(float(s1b['inherited_multiplier_alt'].iloc[0]))} (the value before 2026-09-23): "
      f"median top-25 interval width "
      f"{_fmt(float((s1b['rank_mc_p95_configured'] - s1b['rank_mc_p05_configured']).median()), 0)} vs "
      f"{_fmt(float((s1b['rank_mc_p95_alt'] - s1b['rank_mc_p05_alt']).median()), 0)}; mean P(top-25) of the top-25 "
      f"{_fmt(float(s1b['p_top25_configured'].mean()))} vs {_fmt(float(s1b['p_top25_alt'].mean()))}; top-10 counties "
      f"with P(top-10) >= 0.5: {int((s1b.loc[s1b['rank_equal'] <= 10, 'p_top10_configured'] >= 0.5).sum())} vs "
      f"{int((s1b.loc[s1b['rank_equal'] <= 10, 'p_top10_alt'] >= 0.5).sum())} "
      "(`test7_sensitivity_s1b_inherited_multiplier.csv`).")
    A("")
    sc = T["test7_sensitivity_comparisons"]
    A(md_table(sc, ["analysis", "condition_id_a", "measurement_id_b", "n_compared", "spearman_composite",
                    "kendall_tau_b", "top10_overlap", "top25_overlap", "spearman_clinic_capacity",
                    "spearman_research_readiness"]))
    A("")
    fx = O._features()
    fx = fx[(fx["condition_id"] == "long_covid") & (fx["geo_level"] == "county") & fx["burden_proxy_value"].isna()]
    miss_st = fx["geo_id"].str[:2].map(lambda x: next((a for a, f in _state_abbr().items() if f == x), x))
    A("Each comparison is over the counties ranked in both (a) and (b). Under S2 the "
      f"{len(fx)} counties without a PLACES PHLTH value ("
      + ", ".join(f"{k} {v}" for k, v in miss_st.value_counts().items())
      + "; PLACES lacks the 2023 BRFSS measures for Kentucky and Pennsylvania) are incomplete and drop out, by the "
      "same completeness rule (before 2026-09-24 they were ranked with the burden weight renormalised).")
    A("")
    s6 = T["test7_sensitivity_specialist_only"]
    s6q = T["test7_sensitivity_specialist_only_qa"]
    s6p = _p(s6, condition_id=PKEY[0], measurement_id=PKEY[1])
    A("**S6 specialist_only** (declared 2026-09-24 before it was computed; plan deviation 8). Primary care is removed "
      "from relevant-provider density wherever the composite uses it: the diagnostic desert's provider term becomes "
      "1 - pr(relevant specialists per 100k) (the condition's core groups minus primary_care, per member condition) and "
      "clinic_capacity uses the implementer groups minus primary_care (facility types such as FQHCs stay). The two "
      "decomposition rows swap one term at a time (clinic capacity only = S4).")
    A("")
    A(md_table(s6p.sort_values(["geo_level", "declared"], ascending=[True, False]),
               ["geo_level", "variant", "n_compared", "top10_overlap", "top25_overlap", "top25_jaccard", "kendall_tau_b",
                "spearman_composite", "spearman_diagnostic_desert_pct", "spearman_clinic_capacity_pct",
                "entering_top10"]))
    A("")
    s6c = s6[(s6["geo_level"] == "county") & (s6["variant"] == "specialist_only")]
    A(f"All {len(s6c)} county combinations (specialist_only vs primary): top-10 overlap "
      f"{_fmt(int(s6c['top10_overlap'].min()))}-{_fmt(int(s6c['top10_overlap'].max()))} (median "
      f"{_fmt(float(s6c['top10_overlap'].median()), 1)}), top-25 overlap {_fmt(int(s6c['top25_overlap'].min()))}-"
      f"{_fmt(int(s6c['top25_overlap'].max()))}, Kendall tau-b {_fmt(s6c['kendall_tau_b'].min())}-"
      f"{_fmt(s6c['kendall_tau_b'].max())} (`test7_sensitivity_specialist_only.csv`; the primary query's regions in "
      "`test7_sensitivity_specialist_only_top25.csv`).")
    A("")
    cq = s6q[s6q["check"].str.contains("this module")]
    sh = s6q[s6q["share_primary_care"].notna()]
    A(f"Provider counts (`test7_sensitivity_specialist_only_qa.csv`): the relevant-NPI and specialist counts computed "
      f"here from `providers` equal the geography module's `relevant_providers_n` / `relevant_specialists_n` in "
      f"{_fmt(int(cq.loc[cq['check'].str.contains('geo_condition_features'), 'n_identical'].sum()))} of "
      f"{_fmt(int(cq.loc[cq['check'].str.contains('geo_condition_features'), 'n_geographies'].sum()))} condition x "
      f"geography rows, and the per-group county counts equal `provider_density_county` in "
      f"{_fmt(int(cq.loc[cq['check'].str.contains('provider_density_county'), 'n_identical'].sum()))} of "
      f"{_fmt(int(cq.loc[cq['check'].str.contains('provider_density_county'), 'n_geographies'].sum()))} group x county "
      f"rows. Primary care is {100 * sh['share_primary_care'].min():.1f}-{100 * sh['share_primary_care'].max():.1f}% of "
      "the relevant NPIs (long_covid, me_cfs, pots, dysautonomia; county and state sums).")
    A("")
    s5 = sc[sc["analysis"].str.startswith("S5")]
    s4 = sc[sc["analysis"].str.startswith("S4")]
    s5c = s5[s5["measurement_id_b"] == "nailfold_capillaroscopy"].iloc[0]
    A(f"S5 (the adapter swap) runs the nailfold-capillaroscopy stub bundle and the two clinic-test bundles through the "
      f"identical code path (`make_combo` -> `facility_measures` -> `evaluate`); only the member classes and curated "
      f"implementer groups change. Wearable vs capillaroscopy: composite Spearman {_fmt(s5c['spearman_composite'])}, "
      f"clinic_capacity Spearman {_fmt(s5c['spearman_clinic_capacity'])}, top-10 overlap "
      f"{_fmt(s5c['top10_overlap'])}, top-25 overlap {_fmt(s5c['top25_overlap'])}. Clinic capacity correlates "
      f"{_fmt(s5['spearman_clinic_capacity'].min())}-{_fmt(s5['spearman_clinic_capacity'].max())} across bundles; the "
      f"lowest is exercise capacity testing, the only bundle without primary care among its implementer groups. "
      f"Removing primary care (S4) changes clinic capacity (Spearman {_fmt(s4['spearman_clinic_capacity'].min())}-"
      f"{_fmt(s4['spearman_clinic_capacity'].max())} with the original) but keeps "
      f"{_fmt(int(s4['top25_overlap'].min()))}-{_fmt(int(s4['top25_overlap'].max()))} of the top-25. These are "
      "equal-weight rankings, which contain no measurement-performance information; S5 with the measurement's own "
      "evidence is in section 12.")
    A("")
    # ---- negative / null
    A("## 9. Negative and null results (as prominent as the positive ones)")
    A("")
    s5 = T["test7_sensitivity_comparisons"]
    s5 = s5[s5["analysis"].str.startswith("S5")]
    A(f"* Clinic capacity is not measurement-specific: the location shuffle leaves it almost unchanged "
      f"(mean Spearman {_fmt(cap['component_rho_shuffled_vs_observed_mean'])} for the wearable bundle), and it "
      f"correlates {_fmt(s5['spearman_clinic_capacity'].min())}-{_fmt(s5['spearman_clinic_capacity'].max())} "
      "across measurement bundles (S5).")
    zcc = T["test5_readiness_zero_counts"]
    zcc = zcc[(zcc["geo_level"] == "county") & (zcc["subset"] == "equal top-25")]
    A("* Research readiness is not condition-specific study experience: counties with no condition trial and no NIH "
      "project in reach tie at the bottom of those sub-percentiles, and technology-experience trials of any target "
      "condition can lift their research_readiness percentile (up to "
      f"{_fmt(float(zcc['research_readiness_pct_max_among_zero'].max()), 2)} among the county top-25 lists; section 3).")
    A(f"* The top-10 is fragile: {_fmt(mcp['n_top10_with_p_top10_ge_0_5'])} of 10 primary-query counties stay in the "
      "top-10 in half the Monte Carlo draws; see the per-region intervals in section 3.")
    A("* Long-COVID county burden has no within-state information (inherited); county rankings for long COVID "
      "differentiate counties within a state only through vulnerability, desert access terms, capacity and readiness.")
    A("* POTS and dysautonomia rankings contain no burden information (level D everywhere).")
    s6pp = _p(T["test7_sensitivity_specialist_only"], condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2],
              variant="specialist_only").iloc[0]
    A(f"* Relevant-provider density is mostly primary-care supply: counted without primary care (S6), "
      f"{_fmt(10 - int(s6pp['top10_overlap']))} of the primary query's top-10 counties change.")
    A(f"* {n_not} of {len(ph)} restated phenotype-level controls are null, including the Long COVID wearable-only models.")
    A("* The composite double-counts: the diagnostic desert contains burden and vulnerability, and its access terms "
      "oppose clinic_capacity and research_readiness (as SPEC defines them); the equal-weight composite therefore "
      "partly cancels access against capacity.")
    A("")
    # ---- limitations
    A("## 10. Limitations")
    A("")
    A("* Measurement performance (section 12) is measured against healthy or recovered controls, which overstates "
      "real-world performance; the only look-alike comparators in the table are MY-LC convalescent controls (cortisol, "
      "not a scored bundle) and PDE-007's Long COVID without POTS (a published claim). Expected yields are planning "
      "estimates for a candidate pilot, not predictions of diagnoses.")
    A("")
    for u in REC.GENERAL_UNCERTAINTIES:
        A(f"* {u}")
    A("* HPS is an experimental, low-response online survey; its state estimates have wide CIs and the series ended "
      "2024-09-16. CMS 51 covers Medicare FFS beneficiaries only, is published as integer percentages with size-class "
      "denominators, and is a symptom composite for ME/CFS (level C).")
    A("* The Monte Carlo sd for CMS values (binomial SE at the size-class midpoint + rounding) and the inherited-burden "
      "deviation tau are pre-specified modelling assumptions, not published uncertainty; the inherited-burden "
      "multiplier (x2, the level-C value) is a later modelling assumption (sensitivity S1b).")
    A("* technology_saturation is a research-deployment proxy (facilities with registered trials mentioning the "
      "measurement), not market penetration; no public source measures device deployment by place.")
    A("* Research readiness counts facilities linked by entity resolution (V1 recall 0.53 for HRSA-NPPES links); "
      "large academic centres may appear as several facilities, and ClinicalTrials.gov sites without a usable ZIP are "
      "city centroids.")
    A("* The ranking universe excludes counties below 10,000 residents (scored and flagged, not ranked), and, for "
      "me_cfs and the sets that contain it, the counties without a CMS value (incomplete burden: the Connecticut "
      "planning regions, which CMS publishes only on the legacy counties, and a few others; listed separately).")
    A("")
    L.extend(metric_link_section())
    A("## 13. Files")
    A("")
    for name in sorted(T):
        A(f"* `results/tables/{name}.csv`")
    for f in figs + maps:
        A(f"* `{_repo_relative(f)}`")
    A("* `results/example_deployment_recommendation.json`")
    A("* `data/processed/deployment_opportunities.parquet`, `data/processed/deployment_candidates.parquet`")
    A("")
    text = "\n".join(L) + "\n"
    (RESULTS / "SCORING_RESULTS.md").write_text(text)
    return text


def _mid_short(m: str) -> str:
    return {"wearable_autonomic_activity_monitoring": "wearable", "nailfold_capillaroscopy": "capillaroscopy",
            "autonomic_function_testing": "autonomic testing", "exercise_capacity_testing": "CPET"}.get(m, m)


def _comp_txt(row) -> str:
    return ", ".join(f"{_mid_short(m)} {int(row[f'n_{m}'])}" for m in O.MEASUREMENTS if row.get(f"ranked_{m}")) or "none"


def metric_link_headline() -> str:
    top = T["metric_link_top10_primary"]
    ew10 = set(top.loc[top["rank_evidence_weighted"] <= 10, "geo_id"])
    eq10 = set(top.loc[top["rank_equal"] <= 10, "geo_id"])
    t5 = T["test5_contrasts"]
    t5e = _p(t5, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2], contrast=O.EVIDENCE_WEIGHT_SET)
    r25 = t5e[t5e["top_n"] == 25].iloc[0]
    jc = _p(T["scoring_joint_ranking_composition"], condition_id=PKEY[0], geo_level=PKEY[2])
    j_eq = jc[(jc["weight_set"] == "equal") & (jc["top_n"] == 25)].iloc[0]
    j_ew = jc[(jc["weight_set"] == O.EVIDENCE_WEIGHT_SET) & (jc["top_n"] == 25)].iloc[0]
    s5 = T["test7_sensitivity_s5_evidence"]
    unr = [_mid_short(m) for m in O.MEASUREMENTS if not j_ew.get(f"ranked_{m}")]
    cap = s5[(s5["measurement_id_b"] == "nailfold_capillaroscopy") & (s5["weight_set"] == "equal")].iloc[0]
    return (f"* **Metric -> translation link (added 2026-09-25; docs/ANALYSIS_PLAN_METRIC_LINK.md, pre-specified): the "
            f"measurement's own evidence now enters the ranking where it exists, and only there.** Under the new "
            f"primary-alternative weight set `evidence_weighted`, the primary query keeps {len(ew10 & eq10)} of the "
            f"equal-weight top-10 (top-25 overlap {_fmt(r25['overlap'])} of 25, Kendall tau-b "
            f"{_fmt(r25['kendall_tau_b'])}); within one measurement the change comes only from expected_yield and the "
            "weights, because the evidence score is the same for every region. Across measurements it decides: "
            f"{', '.join(unr) or 'no bundle'} have no usable performance record for Long COVID or ME/CFS (UNKNOWN) and "
            f"are not ranked under evidence_weighted, where under equal weights capillaroscopy ranked almost like the "
            f"wearable bundle (Spearman {_fmt(cap['spearman_composite'])}, top-10 overlap {_fmt(cap['top10_overlap'])}). "
            f"Joint top-25 over region x measurement pairs: equal = {_comp_txt(j_eq)}; evidence_weighted = "
            f"{_comp_txt(j_ew)} (section 12). Performance is against healthy controls (overstated) and the expected "
            "yields are planning estimates for a candidate pilot, not predictions of diagnoses.")


def evidence_pair_text(jp: pd.DataFrame) -> str:
    """Data-driven comparison of the bundles that do have known performance for the primary condition set."""
    perf = T["measurement_performance"]
    pp = perf[(perf["condition_id"] == PKEY[0]) & perf["measurement_id"].isin(O.MEASUREMENTS)
              & (perf["performance_status"] == "known")].sort_values("measurement_evidence", ascending=False)
    if len(pp) < 2:
        return ""
    j = jp[(jp["weight_set"] == O.EVIDENCE_WEIGHT_SET) & (jp["top_n"] == 25)].iloc[0]
    parts = [f"{_mid_short(r.measurement_id)} ({r.selected_record_id}: AUROC lower bound {r.auroc_ci_low:.3f}, "
             f"measurement_evidence {r.measurement_evidence:.3f}, sensitivity {r.op_sensitivity:.2f}; "
             f"{int(j[f'n_{r.measurement_id}'])} of the joint evidence-weighted top-25)" for r in pp.itertuples()]
    txt = ("Among the bundles with known performance: " + "; ".join(parts) + ". ")
    if {"exercise_capacity_testing", "wearable_autonomic_activity_monitoring"} <= set(pp["measurement_id"]):
        txt += ("Both records come from the same MUSCLE-ME people (diagnosed patients vs healthy controls), so the "
                "difference between them is a statement about that case-control sample, not about deployability: "
                "CPET is a maximal exercise test (deployment complexity 4; post-exertional risk) with narrower "
                "implementation reach, and the composite does not weigh burden on the patient.")
    return txt


def metric_link_section() -> list[str]:
    L = []
    A = L.append
    perf = T["measurement_performance"]
    rec = T["measurement_performance_records"]
    st = T["metric_link_combination_status"]
    top = T["metric_link_top10_primary"]
    s5 = T["test7_sensitivity_s5_evidence"]
    jc = T["scoring_joint_ranking_composition"]
    jt = T["scoring_joint_ranking_top25_primary"]
    lcy = T["metric_link_long_covid_state_yield"]
    nr = T["metric_link_unranked_no_reach"]
    t5 = T["test5_contrasts"]
    A("## 12. The metric -> translation link: measurement performance, expected yield, evidence-weighted ranking")
    A("")
    A("Plan: `docs/ANALYSIS_PLAN_METRIC_LINK.md` (pre-specified 2026-09-25, before any number below was computed). "
      "Module `measure_it.scoring.metric_link` (pipeline step `scoring_metric_link`); tables "
      "`data/processed/measurement_performance` (`object_id = measurement_performance:<condition>|<measurement>`), "
      "`measurement_performance_records`, and the new columns of `deployment_opportunities`.")
    A("")
    A("**Why this was needed.** Until 2026-09-25 a measurement entered the ranking only through its curated implementer "
      "groups (clinic capacity) and its trial mentions (research readiness). The phenotype evidence was displayed next "
      "to recommendations but never scored, so a measurement with no public performance evidence (nailfold "
      "capillaroscopy) ranked almost like the wearable bundle. Two facts were fixed before computing: a per-measurement "
      "constant cannot reorder regions within one measurement, and a percentile-normalised yield would lose the "
      "sensitivity; so `measurement_evidence` and `expected_yield` are on absolute 0-1 scales.")
    A("")
    A("### 12.1 Performance records (tiers: 1 own computation on a case definition; 2 own computation on a proxy or "
      "self-reported label; 3 published claim, not reproduced)")
    A("")
    A(md_table(rec.sort_values(["quality_tier", "target_condition", "record_id"]),
               ["record_id", "target_condition", "primary_class", "quality_tier", "comparator_kind", "n_cases",
                "n_controls", "auroc", "auroc_ci_low", "auroc_ci_high", "op_sensitivity", "op_sensitivity_ci_low",
                "op_sensitivity_ci_high", "op_specificity", "performance_status", "measurement_evidence"], nd=3))
    A("")
    A("`op_sensitivity` = sensitivity at >= 0.90 specificity: for MUSCLE-ME and Uwakwe computed here from the stored "
      "person-level values (threshold on the controls, in-sample; stratified bootstrap CI with the threshold "
      "re-chosen); MY-LC is the producer's value; published rows only when the paper reports specificity >= 0.90 "
      "(PDE-007). OWN-MM-VO2-LC and OWN-MM-VO2-ME are the two new AUROCs (same function as the producer). Records "
      "without person-level predictions (NHANES proxy, Appelman CV models) have no operating point: partial.")
    A("")
    used = perf[perf["performance_status"] != UNKNOWN]
    A(f"### 12.2 The scored record per condition x measurement ({len(perf)} grid cells; "
      f"{int((perf['performance_status'] == 'known').sum())} known, {int((perf['performance_status'] == 'partial').sum())} "
      f"partial, {int((perf['performance_status'] == UNKNOWN).sum())} UNKNOWN)")
    A("")
    A(md_table(used.sort_values(["condition_id", "measurement_id"]),
               ["condition_id", "measurement_id", "performance_status", "quality_tier", "selected_record_id", "auroc",
                "auroc_ci_low", "op_sensitivity", "op_specificity", "measurement_evidence", "evidence_conflict"], nd=3))
    A("")
    sc4 = st[(st["geo_level"] == "county")].copy()
    A("Scored combinations (county): performance status and how many counties each ranking covers.")
    A("")
    A(md_table(sc4.sort_values(["condition_id", "measurement_id"]),
               ["condition_id", "measurement_id", "performance_status", "tier", "record", "measurement_evidence",
                "op_sensitivity", "expected_yield_basis", "n_ranked_equal", "n_ranked_evidence_weighted"], nd=3))
    A("")
    conf = perf[perf["evidence_conflict"].astype(bool)]
    for r in conf.itertuples():
        A(f"* Evidence conflict, {r.condition_id} x {r.measurement_id}: {r.evidence_conflict_note}. The scored record is "
          f"{r.selected_record_id} (higher tier); the conflict is stated in every recommendation's uncertainties.")
    A("* Every UNKNOWN cell is UNKNOWN: nothing is imputed from a neighbouring condition, set or measurement. The demo "
      "cluster has no record at all (no analysis pools its four members). ME/CFS x autonomic testing has one attached "
      "published row (PDE-022, sensitivity without specificity) and so no usable metric.")
    A("")
    # ---- top-10
    A("### 12.3 Primary query: top-10 under evidence_weighted and under equal")
    A("")
    ew = top.sort_values("rank_evidence_weighted")
    A(md_table(ew[ew["rank_evidence_weighted"] <= 10],
               ["rank_evidence_weighted", "rank_equal", "geo_name", "state_abbr", "composite_evidence_weighted",
                "rank_mc_ew_p05", "rank_mc_ew_p95", "p_top10_mc_ew", "burden_pct", "expected_yield", "reach",
                "expected_yield_index", "expected_false_positives_upper"], nd=2))
    A("")
    eq = top.sort_values("rank_equal")
    A("Equal weights (continuity):")
    A("")
    A(md_table(eq[eq["rank_equal"] <= 10], ["rank_equal", "rank_evidence_weighted", "geo_name", "state_abbr",
                                              "composite_equal", "rank_mc_p05", "rank_mc_p95", "p_top10_mc"], nd=2))
    A("")
    t5e = t5[(t5["contrast"] == O.EVIDENCE_WEIGHT_SET)]
    t5p = _p(t5e, condition_id=PKEY[0], measurement_id=PKEY[1], geo_level=PKEY[2])
    r10, r25 = t5p[t5p["top_n"] == 10].iloc[0], t5p[t5p["top_n"] == 25].iloc[0]
    p0 = ew.iloc[0]
    A(f"T5-4 (evidence_weighted vs equal, primary query): top-10 overlap {_fmt(r10['overlap'])} of 10, top-25 "
      f"{_fmt(r25['overlap'])} of 25 (Jaccard {_fmt(r25['jaccard'])}), Kendall tau-b {_fmt(r25['kendall_tau_b'])}, "
      f"Spearman {_fmt(r25['spearman'])}. The wearable record used is OWN-MM-STEPS-POOLED (MUSCLE-ME steps, LC + "
      f"ME/CFS vs healthy, tier 1): measurement_evidence {_fmt(float(p0['measurement_evidence']))} for every county, "
      f"sensitivity {_fmt(float(p0['op_sensitivity']))} at specificity {_fmt(float(p0['op_specificity']))}; reach is "
      f"1.0 in {int((top['reach'] == 1).sum())} of these {len(top)} counties, so expected_yield is essentially the "
      "burden percentile x sensitivity and the list moves towards higher-burden counties. The evidence-weighted "
      "Monte Carlo intervals (burden, sensitivity/specificity CI, Dirichlet weights) are wider than the equal-weight "
      "ones: the short list is no more stable than before.")
    A("")
    t5c = t5e[(t5e["top_n"] == 10) & (t5e["status"] == "ranked")]
    A(f"Across the {len(t5c)} combinations ranked under both (county and state): top-10 overlap "
      f"{_fmt(int(t5c['overlap'].min()))}-{_fmt(int(t5c['overlap'].max()))}, Kendall tau-b "
      f"{_fmt(t5c['kendall_tau_b'].min())}-{_fmt(t5c['kendall_tau_b'].max())}; "
      f"{int((t5e['status'] != 'ranked').sum() // 2)} combinations are not ranked under evidence_weighted "
      "(performance UNKNOWN or partial; `test5_contrasts.csv`, contrast `evidence_weighted`).")
    A("")
    # ---- S5
    A("### 12.4 S5 with evidence: does the adapter swap now differ because of evidence?")
    A("")
    A(md_table(s5, ["weight_set", "measurement_id_b", "performance_status_b", "status", "n_compared",
                    "spearman_composite", "kendall_tau_b", "top10_overlap", "top25_overlap"], nd=3))
    A("")
    jp = _p(jc, condition_id=PKEY[0], geo_level=PKEY[2])
    A("Joint ranking over region x measurement pairs (primary condition set, county; "
      "`scoring_joint_ranking_composition.csv`, pairs in `scoring_joint_ranking_top25_primary.csv`):")
    A("")
    cols = ["weight_set", "top_n", "n_pairs_ranked", "n_distinct_regions"] + [f"n_{m}" for m in O.MEASUREMENTS]
    A(md_table(jp, cols))
    A("")
    j_eq = jp[(jp["weight_set"] == "equal") & (jp["top_n"] == 25)].iloc[0]
    j_ew = jp[(jp["weight_set"] == O.EVIDENCE_WEIGHT_SET) & (jp["top_n"] == 25)].iloc[0]
    differs = any(int(j_eq[f"n_{m}"]) != int(j_ew[f"n_{m}"]) for m in O.MEASUREMENTS)
    unr = [m for m in O.MEASUREMENTS if j_eq.get(f"ranked_{m}") and not j_ew.get(f"ranked_{m}")]
    A(f"Pre-specified reading: the ranking differs because of evidence if the joint top-25 bundle composition differs "
      f"or a bundle ranked under equal is unranked under evidence_weighted. Composition differs: {differs}; "
      f"unranked for lack of performance evidence: {', '.join(unr) or 'none'}. "
      "**Honest reading:** within a bundle, evidence cannot reorder counties (a within-bundle comparison under "
      "evidence_weighted exists only where both bundles have known performance); what changed is that bundles without "
      "a performance record for this condition set are no longer ranked as if they were equivalent. "
      + evidence_pair_text(jp))
    A("")
    # ---- yields
    A("### 12.5 Expected yield and false positives")
    A("")
    A("Count path (a defensible count exists: long COVID at state level, HPS % of all adults x ACS adults 18+). "
      "Planning estimates for a candidate pilot at the scored operating point, not predictions of diagnoses; "
      "performance is against healthy controls and overstates real-world performance:")
    A("")
    lw = lcy[lcy["measurement_id"] == PKEY[1]].sort_values("rank_evidence_weighted").head(10)
    A(md_table(lw, ["geo_name", "burden_value", "adults_18plus", "burden_count", "reach", "op_sensitivity",
                    "op_specificity", "expected_detectable_cases", "expected_yield_mc_p05", "expected_yield_mc_p95",
                    "expected_false_positives", "expected_ppv", "false_positives_per_detected_case",
                    "rank_evidence_weighted"], nd=3))
    A("")
    fpd = lcy["false_positives_per_detected_case"]
    A(f"Across the long-COVID state rows with known performance, every state has more expected false positives than "
      f"detectable cases (false positives per detected case {_fmt(float(fpd.min()), 2)}-{_fmt(float(fpd.max()), 2)}; "
      f"PPV {_fmt(float(lcy['expected_ppv'].min()), 2)}-{_fmt(float(lcy['expected_ppv'].max()), 2)}): at a prevalence "
      "of a few percent, a 90%-specificity test flags mostly people without the condition.")
    A("")
    A("Everywhere else (county rows, ME/CFS, condition sets) there is no defensible count at that resolution: the "
      "yield is an index (burden percentile x sensitivity x reach) and false positives are an upper bound "
      "((1 - specificity) x reached adults), shown per county in `metric_link_top10_primary.csv` and in the "
      "recommendation JSON.")
    A("")
    if len(nr):
        g = nr.groupby(["condition_id", "measurement_id", "geo_level"]).size()
        A(f"Rule 5.3: {int(nr['geo_id'].nunique())} ranking-eligible regions have no populated ZCTA whose internal point "
          f"lies in them, so reach (and the evidence-weighted composite) is UNKNOWN there: "
          + "; ".join(sorted(nr["geo_name"].unique())) + f" ({len(nr)} region x combination rows, "
          "`metric_link_unranked_no_reach.csv`). They keep their equal-weight ranks.")
        A("")
    A("Limitations of the link: two of the three tier-1 datasets are the same people (MUSCLE-ME and Appelman 2024); "
      "every scored own record is a case-control comparison against healthy controls; the MUSCLE-ME step result does "
      "not survive its own worst-case missing-control bound (AUROC 0.61); the sensitivity thresholds are chosen on the "
      "same controls; tier factors and weights are declared assumptions; reach assumes that an implementer-group "
      "facility within 50 km could run the measurement.")
    A("")
    return L


def _repo_relative(path) -> str:
    """A path relative to the checkout (whatever its directory is called); unchanged when outside it."""
    from pathlib import Path
    p = Path(path)
    try:
        return p.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def run(build: bool = True, workers: int = 6) -> None:
    t0 = time.time()
    if not table_exists(ML.TABLE):   # normally built by the pipeline step scoring_metric_link
        ML.build()
    if build or not table_exists(O.TABLE):
        O.build()
    print(f"[report] table ready ({time.time() - t0:.0f}s)", flush=True)
    run_test5()
    print(f"[report] test 5 ({time.time() - t0:.0f}s)", flush=True)
    run_test7()
    mc_decomposition()
    print(f"[report] test 7 ({time.time() - t0:.0f}s)", flush=True)
    run_sensitivity()
    run_metric_link()
    run_burden_rule()
    print(f"[report] sensitivity ({time.time() - t0:.0f}s)", flush=True)
    run_test6(workers=workers)
    print(f"[report] test 6 ({time.time() - t0:.0f}s)", flush=True)
    demo = REC.run()
    for name in ("demo_top10_regions", "demo_top10_candidate_sites", "demo_top10_research_evidence",
                 "demo_top10_regions_evidence_weighted"):
        T[name] = pd.read_csv(TABLES / f"{name}.csv")
    maps = [str(p) for p in sorted((O.PROCESSED.parent.parent / "results" / "maps" / "drafts").glob("*.png"))
            if p.name.startswith(("demo_", "test5_"))]
    figs = draft_figures()
    write_markdown(demo, figs, maps)
    print(f"[report] done ({time.time() - t0:.0f}s)", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Scoring: deployment opportunities, Tests 5-7, demo recommendation")
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args(argv)
    run(build=not a.no_build, workers=min(a.workers, 6))


if __name__ == "__main__":
    main()
