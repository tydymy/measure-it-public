"""Test 7 (ranking robustness across weight sets) and the pre-specified sensitivity analyses S1-S5.

Pre-specified in docs/ANALYSIS_PLAN_SCORING.md sections 6-7.

Test 7, per condition (or set) x measurement x level:
    * ranks under the 7 configured weight sets + `saturation_adjusted`, and under 1,000 flat Dirichlet(1,...,1)
      random weight vectors over the five components (observed burden);
    * pairwise Kendall tau-b between named sets; Kendall's W over named sets (with / without burden_only) and over
      random draws; per-region rank min/max/range/IQR (named) and 5th-95th percentile (random); fraction of named
      sets / random draws with the region in the top-10 / top-25;
    * unstable region: top-25 under >= 1 named set AND (top-25 in < 50% of named sets OR < 50% of random draws).
Sensitivity (primary query unless stated):
    S1 Monte Carlo with tau = 0 (pure state inheritance)          S2 long-COVID county burden = PLACES PHLTH proxy
    S3 saturation_adjusted weight set (inside Test 7)             S4 clinic_capacity without primary care
    S5 measurement swap (same pipeline, different measurement bundle)
    S1b inherited-burden uncertainty multiplier 1.0 vs the configured value (added 2026-09-23 with the multiplier;
        not part of the pre-specified plan)
    S6 specialist_only: relevant-provider density without primary care in the diagnostic desert AND clinic_capacity
        (declared 2026-09-24 before it was computed; plan deviation 8), with the two single swaps as decomposition
Burden completeness rule (plan deviation 7, 2026-09-24): before/after rankings of the old rule ('available': a set
burden from the members with a value) and the current one ('complete'), regenerated in every run.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from ..config import SEED
from . import opportunity as O
from .controls import jaccard, rank_corr, spearman, top_set

N_RANDOM = 1000
UNSTABLE_FRAC = 0.5
TOP = 25


def kendalls_w(R: np.ndarray) -> float:
    """Kendall's coefficient of concordance for m raters x n objects of untied ranks 1..n."""
    R = np.asarray(R, dtype=float)
    m, n = R.shape
    if m < 2 or n < 2:
        return np.nan
    Rs = R.sum(axis=0)
    S = float(((Rs - Rs.mean()) ** 2).sum())
    return 12.0 * S / (m ** 2 * (n ** 3 - n))


def _common_ranks(C: np.ndarray, eligible: np.ndarray) -> np.ndarray:
    """Re-rank each row (rater) among the columns that are eligible and non-NaN in every row."""
    ok = eligible & ~np.isnan(C).any(axis=0)
    return O.rank_rows(C[:, ok], np.ones(int(ok.sum()), bool))


def random_composites(cb: O.Combo, n: int = N_RANDOM, seed: int = SEED, chunk: int = 200) -> np.ndarray:
    rng = np.random.default_rng(seed + 7000)
    out = []
    done = 0
    while done < n:
        D = min(chunk, n - done)
        W = rng.dirichlet(np.ones(len(O.COMPONENTS)), size=D)
        out.append(O.evaluate(cb, weights=W, order=O.COMPONENTS)["composite"])
        done += D
    return np.vstack(out)


def test7_combo(cb: O.Combo, n_random: int = N_RANDOM, seed: int = SEED) -> tuple[dict, list[dict], pd.DataFrame]:
    ev = O.evaluate(cb)
    named = O.named_composites(cb, ev)
    names = list(named)
    elig = cb.eligible
    geo = cb.geo
    key = {"condition_id": cb.cond["condition_id"], "measurement_id": cb.meas["measurement_id"],
           "geo_level": cb.level}
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            ok = elig & ~np.isnan(named[a]) & ~np.isnan(named[b])
            tau = kendalltau(named[a][ok], named[b][ok]).statistic if ok.sum() > 2 else np.nan
            ta, tb = top_set(named[a], elig, TOP), top_set(named[b], elig, TOP)
            pairs.append({**key, "weight_set_a": a, "weight_set_b": b, "kendall_tau_b": tau, "n": int(ok.sum()),
                          "top25_overlap": len(set(ta) & set(tb)), "top25_jaccard": jaccard(ta, tb)})
    Cn = np.vstack([named[k] for k in names])
    defined = [k for k in names if (~np.isnan(named[k][elig])).any()]
    W_named = kendalls_w(_common_ranks(np.vstack([named[k] for k in defined]), elig))
    no_bo = [k for k in defined if k != "burden_only"]
    W_named_nb = kendalls_w(_common_ranks(np.vstack([named[k] for k in no_bo]), elig))
    Crand = random_composites(cb, n_random, seed)
    W_rand = kendalls_w(_common_ranks(Crand, elig))
    Rn = np.vstack([O.rank_rows(named[k], elig)[0] for k in names])
    Rr = O.rank_rows(Crand, elig)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        reg = pd.DataFrame({
            **{k: v for k, v in key.items()}, "geo_id": geo["geo_id"].to_numpy(), "geo_name": geo["geo_name"].to_numpy(),
            "state_abbr": geo["state_abbr"].to_numpy(), "eligible": elig,
            **{f"rank_{k}": Rn[i] for i, k in enumerate(names)},
            "n_named_sets_ranked": (~np.isnan(Rn)).sum(axis=0),
            "rank_named_min": np.nanmin(Rn, axis=0), "rank_named_max": np.nanmax(Rn, axis=0),
            "rank_named_iqr": np.nanpercentile(Rn, 75, axis=0) - np.nanpercentile(Rn, 25, axis=0),
            "frac_named_top10": np.nanmean(np.where(np.isnan(Rn), np.nan, Rn <= 10), axis=0),
            "frac_named_top25": np.nanmean(np.where(np.isnan(Rn), np.nan, Rn <= TOP), axis=0),
            "rank_random_p05": np.nanpercentile(Rr, 5, axis=0), "rank_random_p50": np.nanpercentile(Rr, 50, axis=0),
            "rank_random_p95": np.nanpercentile(Rr, 95, axis=0),
            "rank_random_iqr": np.nanpercentile(Rr, 75, axis=0) - np.nanpercentile(Rr, 25, axis=0),
            "frac_random_top10": np.nanmean(np.where(np.isnan(Rr), np.nan, Rr <= 10), axis=0),
            "frac_random_top25": np.nanmean(np.where(np.isnan(Rr), np.nan, Rr <= TOP), axis=0),
        })
    reg["rank_named_range"] = reg["rank_named_max"] - reg["rank_named_min"]
    reg["in_top25_any_named"] = (Rn <= TOP).any(axis=0)
    reg["in_top25_default"] = Rn[names.index("equal")] <= TOP
    reg["unstable"] = reg["in_top25_any_named"] & ((reg["frac_named_top25"] < UNSTABLE_FRAC) |
                                                  (reg["frac_random_top25"] < UNSTABLE_FRAC))
    reasons = []
    for pos, r in enumerate(reg.itertuples(index=False)):
        if not r.unstable:
            reasons.append(None)
            continue
        sets_in = [k for i, k in enumerate(names) if Rn[i][pos] <= TOP]
        why = []
        if r.frac_named_top25 < UNSTABLE_FRAC:
            why.append(f"top-25 in {r.frac_named_top25:.0%} of named weight sets ({', '.join(sets_in)})")
        if r.frac_random_top25 < UNSTABLE_FRAC:
            why.append(f"top-25 in {r.frac_random_top25:.0%} of 1,000 random Dirichlet weight draws")
        why.append(f"named-set rank range {r.rank_named_min:.0f}-{r.rank_named_max:.0f}")
        reasons.append("; ".join(why))
    reg["unstable_reason"] = reasons
    top_def = reg[reg["in_top25_default"]]
    summ = {**key, "n_eligible": int(elig.sum()), "weight_sets": "|".join(names),
            "weight_sets_defined": "|".join(defined),
            "kendalls_w_named": W_named, "kendalls_w_named_excl_burden_only": W_named_nb,
            "kendalls_w_random": W_rand, "n_random_draws": n_random,
            "pairwise_tau_median": float(np.nanmedian([p["kendall_tau_b"] for p in pairs])),
            "pairwise_tau_min": float(np.nanmin([p["kendall_tau_b"] for p in pairs])),
            "pairwise_tau_min_excl_burden_only": float(np.nanmin(
                [p["kendall_tau_b"] for p in pairs if "burden_only" not in (p["weight_set_a"], p["weight_set_b"])])),
            "n_regions_top25_any_named": int(reg["in_top25_any_named"].sum()),
            "n_unstable": int(reg["unstable"].sum()),
            "n_default_top25_unstable": int(top_def["unstable"].sum()),
            "n_default_top10_in_top10_all_named": int(((reg["rank_equal"] <= 10) & (reg["frac_named_top10"] == 1)).sum()),
            "median_named_range_default_top25": float(top_def["rank_named_range"].median()),
            "median_random_p95_default_top25": float(top_def["rank_random_p95"].median())}
    return summ, pairs, reg


# ------------------------------------------------------------------------------------------------------------------
# S1-S5
# ------------------------------------------------------------------------------------------------------------------

def s1_tau(cb: O.Combo, n_draws: int = O.MC_DRAWS) -> pd.DataFrame:
    """Monte Carlo rank intervals with and without the per-county deviation (tau) for inherited burden."""
    ev = O.evaluate(cb)
    rank = O.rank_rows(ev["composite"][0], cb.eligible)[0]
    R1 = O.monte_carlo(cb, n_draws, SEED)
    R0 = O.monte_carlo(cb, n_draws, SEED, tau_scale=0.0)
    a, b = O.summarize_ranks(R1, "rank_mc"), O.summarize_ranks(R0, "rank_mc")
    top = np.flatnonzero(rank <= TOP)
    top = top[np.argsort(rank[top])]
    return pd.DataFrame({"condition_id": cb.cond["condition_id"], "measurement_id": cb.meas["measurement_id"],
                         "geo_level": cb.level, "geo_id": cb.geo["geo_id"].to_numpy()[top],
                         "geo_name": cb.geo["geo_name"].to_numpy()[top], "rank_equal": rank[top],
                         "tau_pp": next((m.tau for m in cb.members if m.inherited.any()), np.nan),
                         "rank_mc_p05_with_tau": a["rank_mc_p05"][top], "rank_mc_p95_with_tau": a["rank_mc_p95"][top],
                         "p_top25_with_tau": a["p_top25_mc"][top],
                         "rank_mc_p05_tau0": b["rank_mc_p05"][top], "rank_mc_p95_tau0": b["rank_mc_p95"][top],
                         "p_top25_tau0": b["p_top25_mc"][top]})


def s1b_inherited_multiplier(cb: O.Combo, n_draws: int = O.MC_DRAWS, alt: float = 1.0) -> pd.DataFrame:
    """Monte Carlo rank intervals with the configured inherited-burden multiplier vs `alt` (1.0 = no inflation of
    the state estimate's error when it is carried to a county; the behaviour before 2026-09-23)."""
    ev = O.evaluate(cb)
    rank = O.rank_rows(ev["composite"][0], cb.eligible)[0]
    cfg = O.inherited_multiplier()
    R1 = O.monte_carlo(cb, n_draws, SEED)
    R0 = O.monte_carlo(cb, n_draws, SEED, inherited_multiplier=alt)
    a, b = O.summarize_ranks(R1, "rank_mc"), O.summarize_ranks(R0, "rank_mc")
    top = np.flatnonzero(rank <= TOP)
    top = top[np.argsort(rank[top])]
    return pd.DataFrame({"condition_id": cb.cond["condition_id"], "measurement_id": cb.meas["measurement_id"],
                         "geo_level": cb.level, "geo_id": cb.geo["geo_id"].to_numpy()[top],
                         "geo_name": cb.geo["geo_name"].to_numpy()[top], "rank_equal": rank[top],
                         "inherited_multiplier_configured": cfg, "inherited_multiplier_alt": alt,
                         "rank_mc_p05_configured": a["rank_mc_p05"][top], "rank_mc_p95_configured": a["rank_mc_p95"][top],
                         "p_top10_configured": a["p_top10_mc"][top], "p_top25_configured": a["p_top25_mc"][top],
                         "rank_mc_p05_alt": b["rank_mc_p05"][top], "rank_mc_p95_alt": b["rank_mc_p95"][top],
                         "p_top10_alt": b["p_top10_mc"][top], "p_top25_alt": b["p_top25_mc"][top]})


def compare(cb_a: O.Combo, cb_b: O.Combo, label: str, extra: dict | None = None) -> dict:
    ea, eb = O.evaluate(cb_a), O.evaluate(cb_b)
    ca, cbb = ea["composite"][0], eb["composite"][0]
    elig = cb_a.eligible & cb_b.eligible
    out = {"analysis": label, "condition_id_a": cb_a.cond["condition_id"], "measurement_id_a": cb_a.meas["measurement_id"],
           "condition_id_b": cb_b.cond["condition_id"], "measurement_id_b": cb_b.meas["measurement_id"],
           "geo_level": cb_a.level, "spearman_composite": spearman(ca, cbb, elig),
           "top10_overlap": len(set(top_set(ca, elig, 10)) & set(top_set(cbb, elig, 10))),
           "top25_overlap": len(set(top_set(ca, elig, TOP)) & set(top_set(cbb, elig, TOP))),
           "spearman_clinic_capacity": spearman(cb_a.cr["clinic_capacity_pct"], cb_b.cr["clinic_capacity_pct"], elig),
           "spearman_research_readiness": spearman(cb_a.cr["research_readiness_pct"],
                                                   cb_b.cr["research_readiness_pct"], elig),
           "spearman_burden_pct": spearman(ea["burden_pct"][0], eb["burden_pct"][0], elig),
           "n_compared": int(elig.sum()), "kendall_tau_b": rank_corr(ca, cbb, elig)[0]}
    out.update(extra or {})
    return out


# ------------------------------------------------------------------------------------------------------------------
# ranking agreement (S6, burden-rule before/after)
# ------------------------------------------------------------------------------------------------------------------

def ranking_agreement(ca: np.ndarray, cb_: np.ndarray, elig: np.ndarray, top_ns=(10, TOP)) -> dict:
    """Top-N overlap / Jaccard and Kendall tau-b / Spearman of two composites over the regions eligible in both."""
    tau, rho, n = rank_corr(np.asarray(ca, dtype=float), np.asarray(cb_, dtype=float), elig)
    out = {"n_compared": n, "kendall_tau_b": tau, "spearman_composite": rho}
    for N in top_ns:
        ta, tb = top_set(ca, elig, N), top_set(cb_, elig, N)
        out[f"top{N}_overlap"] = len(set(ta) & set(tb))
        out[f"top{N}_jaccard"] = jaccard(ta, tb)
    return out


S6_VARIANTS = {
    # declared analysis: primary care removed from relevant-provider density everywhere the composite uses it
    "specialist_only": {"provider_variant": "specialists_only", "exclude_groups": (O.SPECIALIST_EXCLUDED_GROUP,)},
    # decomposition (secondary): one swap at a time
    "desert_provider_term_only": {"provider_variant": "specialists_only"},
    "clinic_capacity_only": {"exclude_groups": (O.SPECIALIST_EXCLUDED_GROUP,)},
}
S6_LABEL = {
    "specialist_only": "S6 specialist_only: desert provider term + clinic_capacity without primary care (declared)",
    "desert_provider_term_only": "S6 decomposition: desert provider term without primary care only",
    "clinic_capacity_only": "S6 decomposition: clinic_capacity without primary care only (= S4)",
}


def s6_specialist_only(cid: str, mid: str, level: str, regions: bool = False) -> tuple[list[dict], pd.DataFrame | None]:
    """S6 vs the primary (equal-weight) ranking of the same combination; `regions` -> per-region ranks too."""
    a = O.make_combo(cid, mid, level)
    ea = O.evaluate(a)
    ca = ea["composite"][0]
    rank_a = O.rank_rows(ca, a.eligible)[0]
    rows, reg = [], None
    for name, kw in S6_VARIANTS.items():
        b = O.make_combo(cid, mid, level, **kw)
        eb = O.evaluate(b)
        cb_ = eb["composite"][0]
        elig = a.eligible & b.eligible
        rank_b = O.rank_rows(cb_, b.eligible)[0]
        ag = ranking_agreement(ca, cb_, elig)
        top10_b = set(top_set(cb_, elig, 10))
        rows.append({"analysis": S6_LABEL[name], "variant": name, "declared": name == "specialist_only",
                     "condition_id": cid, "measurement_id": mid, "geo_level": level, **ag,
                     "spearman_diagnostic_desert_pct": spearman(ea["diagnostic_desert_pct"][0],
                                                                eb["diagnostic_desert_pct"][0], elig),
                     "spearman_clinic_capacity_pct": spearman(a.cr["clinic_capacity_pct"], b.cr["clinic_capacity_pct"],
                                                              elig),
                     "entering_top10": "; ".join(a.geo["geo_name"].iloc[i] for i in
                                                 sorted(top10_b - set(top_set(ca, elig, 10)), key=lambda i: rank_b[i]))})
        if regions and name == "specialist_only":
            top = np.flatnonzero((rank_a <= TOP) | (rank_b <= TOP))
            top = top[np.lexsort((rank_b[top], rank_a[top]))]
            reg = pd.DataFrame({"condition_id": cid, "measurement_id": mid, "geo_level": level,
                                "geo_id": a.geo["geo_id"].to_numpy()[top], "geo_name": a.geo["geo_name"].to_numpy()[top],
                                "state_abbr": a.geo["state_abbr"].to_numpy()[top],
                                "rank_primary": rank_a[top], "rank_specialist_only": rank_b[top],
                                "diagnostic_desert_pct_primary": ea["diagnostic_desert_pct"][0][top],
                                "diagnostic_desert_pct_specialist_only": eb["diagnostic_desert_pct"][0][top],
                                "clinic_capacity_pct_primary": a.cr["clinic_capacity_pct"][top],
                                "clinic_capacity_pct_specialist_only": b.cr["clinic_capacity_pct"][top],
                                "composite_primary": ca[top], "composite_specialist_only": cb_[top]})
    return rows, reg


def s6_provider_qa() -> pd.DataFrame:
    """Relevant-provider counts of this module (from `providers`) vs the geography module and provider_density_county,
    and the primary-care share of relevant NPIs (what S6 removes)."""
    from ..store import read_table
    f = read_table("geo_condition_features", columns=["geo_id", "geo_level", "condition_id", "relevant_providers_n",
                                                      "relevant_specialists_n", "in_analysis_universe"])
    f = f[f["in_analysis_universe"]]
    rows = []
    for level in O.LEVELS:
        geo = O.spatial().levels[level].geo
        for cid in O.SINGLE_CONDITIONS:
            g = f[(f["condition_id"] == cid) & (f["geo_level"] == level)].set_index("geo_id").reindex(geo["geo_id"])
            n_all, _ = O.relevant_provider_counts(cid, level)
            n_sp, _ = O.relevant_provider_counts(cid, level, exclude_groups=(O.SPECIALIST_EXCLUDED_GROUP,))
            for what, here, col in (("relevant NPIs, all core groups", n_all, "relevant_providers_n"),
                                    ("relevant NPIs, core groups minus primary_care", n_sp, "relevant_specialists_n")):
                ref = g[col].to_numpy(float)
                ok = ~np.isnan(ref)
                rows.append({"check": f"{what}: this module (providers) vs geo_condition_features.{col}",
                             "condition_id": cid, "geo_level": level, "n_geographies": int(ok.sum()),
                             "n_identical": int((here[ok] == ref[ok]).sum()),
                             "max_abs_difference": float(np.abs(here[ok] - ref[ok]).max()) if ok.any() else np.nan,
                             "total_here": float(here[ok].sum()), "total_reference": float(ref[ok].sum())})
            rows.append({"check": "primary-care share of relevant NPIs (1 - specialists / all core groups), "
                                  "summed over the geographies",
                         "condition_id": cid, "geo_level": level, "n_geographies": len(geo),
                         "share_primary_care": float(1 - n_sp.sum() / n_all.sum()) if n_all.sum() else np.nan,
                         "total_relevant_npis": float(n_all.sum()), "total_specialist_npis": float(n_sp.sum())})
    pdc = read_table("provider_density_county", columns=["county_fips", "specialty_group", "n_individual_providers"])
    geo = O.spatial().levels["county"].geo
    for grp in (O.SPECIALIST_EXCLUDED_GROUP, "neurology", "cardiology"):
        ref = pdc[pdc["specialty_group"] == grp].set_index("county_fips")["n_individual_providers"].reindex(
            geo["geo_id"]).to_numpy(float)
        here, _ = O.provider_counts_for_groups([grp], "county")
        ok = ~np.isnan(ref)
        rows.append({"check": f"{grp} NPIs per county: this module (providers) vs provider_density_county",
                     "condition_id": None, "geo_level": "county", "n_geographies": int(ok.sum()),
                     "n_identical": int((here[ok] == ref[ok]).sum()),
                     "max_abs_difference": float(np.abs(here[ok] - ref[ok]).max()) if ok.any() else np.nan,
                     "total_here": float(here[ok].sum()), "total_reference": float(ref[ok].sum())})
    return pd.DataFrame(rows)


def set_burden_rule_before_after(cid: str, mid: str, level: str) -> tuple[dict, pd.DataFrame]:
    """Old rule ('available') vs current rule ('complete'): eligible regions, top-10/25, Kendall tau, top-10 lists."""
    new = O.make_combo(cid, mid, level)
    old = O.make_combo(cid, mid, level, burden_rule="available")
    cn, co = O.evaluate(new)["composite"][0], O.evaluate(old)["composite"][0]
    rn, ro = O.rank_rows(cn, new.eligible)[0], O.rank_rows(co, old.eligible)[0]
    bo = {"burden": 1.0}
    en, eo = O.evaluate(new, weights=bo, order=["burden"]), O.evaluate(old, weights=bo, order=["burden"])
    rbn, rbo = O.rank_rows(en["composite"][0], new.eligible)[0], O.rank_rows(eo["composite"][0], old.eligible)[0]
    inc = ~new.burden_complete
    common = new.eligible & old.eligible
    ag = ranking_agreement(co, cn, common)
    geo = new.geo
    summ = {"condition_id": cid, "measurement_id": mid, "geo_level": level,
            "n_incomplete": int(inc.sum()), "n_incomplete_population_eligible": int((inc & new.pop_eligible).sum()),
            "n_ranked_before": int(old.eligible.sum()), "n_ranked_after": int(new.eligible.sum()),
            "top10_overlap_before_after": len(set(np.flatnonzero(ro <= 10)) & set(np.flatnonzero(rn <= 10))),
            "top25_overlap_before_after": len(set(np.flatnonzero(ro <= TOP)) & set(np.flatnonzero(rn <= TOP))),
            "kendall_tau_b_common": ag["kendall_tau_b"], "spearman_common": ag["spearman_composite"],
            "n_common": ag["n_compared"],
            "incomplete_in_top10_before": "; ".join(f"{geo['geo_name'].iloc[i]} ({ro[i]:.0f})"
                                                    for i in np.flatnonzero(inc & (ro <= 10))),
            "incomplete_in_top25_before": "; ".join(f"{geo['geo_name'].iloc[i]} ({ro[i]:.0f})"
                                                    for i in np.flatnonzero(inc & (ro <= TOP))),
            "incomplete_in_burden_only_top25_before": "; ".join(f"{geo['geo_name'].iloc[i]} ({rbo[i]:.0f})"
                                                                for i in np.flatnonzero(inc & (rbo <= TOP)))}
    idx = np.flatnonzero((ro <= 10) | (rn <= 10))
    idx = idx[np.lexsort((np.nan_to_num(ro[idx], nan=1e9), np.nan_to_num(rn[idx], nan=1e9)))]
    status = []
    for i in idx:
        if ro[i] <= 10 and rn[i] <= 10:
            status.append("in top-10 before and after")
        elif ro[i] <= 10 and inc[i]:
            status.append("left: incomplete burden (excluded from the ranking)")
        elif ro[i] <= 10:
            status.append("left: moved down")
        else:
            status.append("entered")
    top = pd.DataFrame({"condition_id": cid, "measurement_id": mid, "geo_level": level,
                        "geo_id": geo["geo_id"].to_numpy()[idx], "geo_name": geo["geo_name"].to_numpy()[idx],
                        "state_abbr": geo["state_abbr"].to_numpy()[idx], "rank_before": ro[idx], "rank_after": rn[idx],
                        "rank_burden_only_before": rbo[idx], "rank_burden_only_after": rbn[idx],
                        "burden_incomplete": inc[idx], "status": status})
    return summ, top


# ------------------------------------------------------------------------------------------------------------------
# S5 with evidence (the metric -> translation link; docs/ANALYSIS_PLAN_METRIC_LINK.md section 6, added 2026-09-25)
# ------------------------------------------------------------------------------------------------------------------

def s5_evidence(cid: str = O.PRIMARY["condition"], level: str = "county") -> pd.DataFrame:
    """Wearable bundle vs each other bundle under `equal` and under `evidence_weighted`: ranking agreement where both
    rankings exist; 'not ranked' with the performance status where one does not (never compared as if equivalent)."""
    base = O.make_combo(cid, O.MEASUREMENTS[0], level)
    eb = O.evaluate(base)
    nb = O.named_composites(base, eb)
    rows = []
    for mid in O.MEASUREMENTS[1:]:
        other = O.make_combo(cid, mid, level)
        eo = O.evaluate(other)
        no = O.named_composites(other, eo)
        elig = base.eligible & other.eligible
        for ws in ("equal", O.EVIDENCE_WEIGHT_SET):
            a, b = nb[ws], no[ws]
            st_a = (base.perf or {}).get("performance_status")
            st_b = (other.perf or {}).get("performance_status")
            row = {"condition_id": cid, "geo_level": level, "weight_set": ws,
                   "measurement_id_a": base.meas["measurement_id"], "measurement_id_b": mid,
                   "performance_status_a": st_a, "performance_status_b": st_b,
                   "performance_tier_a": (base.perf or {}).get("quality_tier"),
                   "performance_tier_b": (other.perf or {}).get("quality_tier"),
                   "measurement_evidence_a": (base.perf or {}).get("measurement_evidence"),
                   "measurement_evidence_b": (other.perf or {}).get("measurement_evidence"),
                   "op_sensitivity_a": (base.perf or {}).get("op_sensitivity"),
                   "op_sensitivity_b": (other.perf or {}).get("op_sensitivity"),
                   "adapter_b": other.meas.get("adapter_name"), "adapter_status_b": other.meas.get("adapter_status")}
            ok_a, ok_b = (elig & ~np.isnan(a)).any(), (elig & ~np.isnan(b)).any()
            if ok_a and ok_b:
                row.update({"status": "both ranked", **ranking_agreement(a, b, elig)})
                row["spearman_expected_yield"] = spearman(np.asarray(eb["components"]["expected_yield"]).ravel(),
                                                          np.asarray(eo["components"]["expected_yield"]).ravel(), elig)
            else:
                miss = [m for m, ok in ((base.meas["measurement_id"], ok_a), (mid, ok_b)) if not ok]
                row.update({"status": "not ranked under " + ws + ": " + ", ".join(
                    f"{m} (measurement performance {st_b if m == mid else st_a})" for m in miss)})
            rows.append(row)
    return pd.DataFrame(rows)


def joint_composition(df: pd.DataFrame, top_ns=(10, TOP)) -> pd.DataFrame:
    """Bundle composition of the joint (region x measurement) top-N per condition x level under equal and
    evidence_weighted (rank_<ws>_joint of deployment_opportunities)."""
    rows = []
    for (cid, lev), g in df.groupby(["condition_id", "geo_level"]):
        for ws in ("equal", O.EVIDENCE_WEIGHT_SET):
            col = f"rank_{ws}_joint"
            if col not in g:
                continue
            n_pairs = int(g[col].notna().sum())
            for N in top_ns:
                top = g[g[col] <= N]
                row = {"condition_id": cid, "geo_level": lev, "weight_set": ws, "top_n": N, "n_pairs_ranked": n_pairs,
                       "n_distinct_regions": int(top["geo_id"].nunique())}
                for mid in O.MEASUREMENTS:
                    row[f"n_{mid}"] = int((top["measurement_id"] == mid).sum())
                    row[f"ranked_{mid}"] = bool(g.loc[g["measurement_id"] == mid, col].notna().any())
                rows.append(row)
    return pd.DataFrame(rows)
