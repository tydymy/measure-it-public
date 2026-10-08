"""Test 5 (does facility/research information change where we would deploy?) and Test 6 (negative controls).

Pre-specified in docs/ANALYSIS_PLAN_SCORING.md section 6.

Test 5 contrasts, each against the full `equal` ranking, among ranking-eligible regions:
    T5-1 burden_only           (SPEC's comparison)
    T5-2 need_only             burden + vulnerability (0.5 / 0.5): no facility or research information at all
    T5-3 need_plus_desert      burden + vulnerability + desert (1/3 each): isolates clinic_capacity + readiness
    T5-4 evidence_weighted     the metric-link weight set vs equal (docs/ANALYSIS_PLAN_METRIC_LINK.md; added
                               2026-09-25): adds measurement_evidence + expected_yield
Metrics: top-10/25 overlap, Jaccard, expected overlap N^2/M, tie-averaged overlap, Kendall tau-b, Spearman, and
regions entering / leaving the top-N with their component contributions as the reason.

Test 6 (scoring):
    (a) geography-label shuffles of burden (all members of a set with one permutation) and of vulnerability,
        national and within state
    (b) clinic-location shuffle: facility and individual-provider location tuples permuted (national, within
        state); clinic_capacity, research_readiness and technology_saturation recomputed for every region
    (c) phenotype-label permutation results restated from the wearable modules' tables
"""
from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

from ..config import SEED, TABLES
from . import opportunity as O

CONTRASTS = {
    "burden_only": {"burden": 1.0},
    "need_only": {"burden": 0.5, "vulnerability": 0.5},
    "need_plus_desert": {"burden": 1 / 3, "vulnerability": 1 / 3, "diagnostic_desert": 1 / 3},
}
CONTRAST_LABEL = {
    "burden_only": "T5-1 burden only vs full (SPEC comparison)",
    "need_only": "T5-2 burden + vulnerability (no facility/research information) vs full",
    "need_plus_desert": "T5-3 burden + vulnerability + desert vs full (adds clinic_capacity + research_readiness)",
}
TOP_NS = (10, 25)
MATERIAL_JACCARD = 0.5
N_TIE_BREAKS = 1000
N_LABEL_SHUFFLES = 200
N_CLINIC_NATIONAL = 200
N_CLINIC_WITHIN = 100
WORKERS = 6


# ------------------------------------------------------------------------------------------------------------------
# pure helpers
# ------------------------------------------------------------------------------------------------------------------

def top_set(score: np.ndarray, eligible: np.ndarray, n: int) -> np.ndarray:
    """Indices of the top-n eligible rows by score (ties by index = geo_id order); NaN never enters."""
    r = O.rank_rows(score, eligible)[0]
    return np.flatnonzero(r <= n)


def jaccard(a, b) -> float:
    a, b = set(map(int, a)), set(map(int, b))
    return len(a & b) / len(a | b) if (a | b) else np.nan


def tie_averaged_overlap(score: np.ndarray, eligible: np.ndarray, ref_top: np.ndarray, n: int,
                         n_rand: int = N_TIE_BREAKS, seed: int = SEED) -> tuple[float, int]:
    """Expected overlap of the top-n of `score` with ref_top when ties in `score` are broken at random.

    Returns (mean overlap, number of eligible rows tied at the top-n boundary value)."""
    s = np.asarray(score, dtype=float)
    idx = np.flatnonzero(eligible & ~np.isnan(s))
    if len(idx) <= n:
        return float(len(set(idx) & set(ref_top))), 0
    vals = s[idx]
    thr = np.sort(vals)[::-1][n - 1]
    above = idx[vals > thr]
    tied = idx[vals == thr]
    k = n - len(above)
    ref = set(map(int, ref_top))
    base = len(set(map(int, above)) & ref)
    if len(tied) <= k or k <= 0:
        return float(base + len(set(map(int, tied[:max(k, 0)])) & ref)), int(len(tied))
    # hypergeometric expectation of tied members of ref drawn into the k free slots
    t_in_ref = len(set(map(int, tied)) & ref)
    return float(base + k * t_in_ref / len(tied)), int(len(tied))


def rank_corr(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> tuple[float, float, int]:
    ok = mask & ~np.isnan(a) & ~np.isnan(b)
    if ok.sum() < 3:
        return np.nan, np.nan, int(ok.sum())
    tau = kendalltau(a[ok], b[ok]).statistic
    rho = spearmanr(a[ok], b[ok]).statistic
    return float(tau), float(rho), int(ok.sum())


def spearman(a, b, mask) -> float:
    ok = mask & ~np.isnan(a) & ~np.isnan(b)
    if ok.sum() < 3:
        return np.nan
    return float(spearmanr(a[ok], b[ok]).statistic)


def contributions(comps: dict, weights: dict, eligible: np.ndarray) -> dict[str, np.ndarray]:
    """w_k * (x_k - mean_k over eligible rows) per component (reason attribution for Test 5)."""
    out = {}
    for k, w in weights.items():
        x = np.asarray(comps[k], dtype=float).ravel()
        v = x[eligible]
        mu = np.nanmean(v) if (~np.isnan(v)).any() else np.nan  # all-NaN (level D burden) -> NaN contribution
        out[k] = w * (x - mu)
    return out


# ------------------------------------------------------------------------------------------------------------------
# Test 5
# ------------------------------------------------------------------------------------------------------------------

def combo_key(cb: O.Combo) -> dict:
    return {"condition_id": cb.cond["condition_id"], "measurement_id": cb.meas["measurement_id"],
            "geo_level": cb.level}


def test5_combo(cb: O.Combo, ev: dict | None = None) -> tuple[list[dict], list[dict]]:
    ev = ev or O.evaluate(cb)
    comps = {k: np.asarray(v, dtype=float).ravel() for k, v in ev["components"].items()}
    full = O.weighted_composite({k: comps[k] for k in O.COMPONENTS}, O.default_weights(), O.COMPONENTS)[0]
    elig = cb.eligible
    geo = cb.geo
    M_elig = int((elig & ~np.isnan(full)).sum())
    summary, moves = [], []
    contrib = contributions({k: comps[k] for k in O.COMPONENTS}, O.default_weights(), elig)
    full_rank = O.rank_rows(full, elig)[0]
    for cname, w in CONTRASTS.items():
        order = list(w)
        sc = O.weighted_composite({k: comps[k] for k in order}, w, order)[0]
        sc_rank = O.rank_rows(sc, elig)[0]
        tau, rho, n_corr = rank_corr(sc, full, elig)
        for N in TOP_NS:
            ft = top_set(full, elig, N)
            ct = top_set(sc, elig, N)
            ov = len(set(ft) & set(ct))
            tie_ov, n_tied = tie_averaged_overlap(sc, elig, ft, N)
            M_c = int((elig & ~np.isnan(sc)).sum())
            jac = jaccard(ft, ct)
            summary.append({**combo_key(cb), "contrast": cname, "contrast_label": CONTRAST_LABEL[cname],
                            "top_n": N, "n_eligible_full": M_elig, "n_eligible_contrast": M_c,
                            "overlap": ov, "jaccard": jac, "overlap_tie_averaged": tie_ov,
                            "jaccard_tie_averaged": tie_ov / (2 * N - tie_ov) if N else np.nan,
                            "n_tied_at_boundary_contrast": n_tied,
                            "expected_overlap_independent": N * N / M_elig if M_elig else np.nan,
                            "kendall_tau_b": tau, "spearman": rho, "n_rank_corr": n_corr,
                            "n_entering": N - ov, "n_leaving": N - ov,
                            "material_change_jaccard_lt_0_5": bool(jac < MATERIAL_JACCARD)})
            for kind, rows in (("entering_full_top", sorted(set(ft) - set(ct), key=lambda i: full_rank[i])),
                               ("leaving_contrast_top", sorted(set(ct) - set(ft), key=lambda i: sc_rank[i]))):
                for i in rows:
                    c = {k: float(contrib[k][i]) for k in O.COMPONENTS}
                    ranked = sorted([kv for kv in c.items() if kv[1] == kv[1]], key=lambda kv: kv[1],
                                    reverse=(kind == "entering_full_top"))
                    drivers = [f"{k} pct {comps[k][i]:.2f} ({v:+.3f})" for k, v in ranked[:2]]
                    if kind == "entering_full_top":
                        reason = (f"enters the full top-{N} (rank {full_rank[i]:.0f}; {cname} rank "
                                  f"{sc_rank[i]:.0f}" + (" (no burden: level D)" if np.isnan(sc[i]) else "") +
                                  "); largest positive contributions: " + "; ".join(drivers))
                    else:
                        reason = (f"leaves: {cname} rank {sc_rank[i]:.0f}, full rank {full_rank[i]:.0f}; "
                                  "smallest (most negative) contributions: " + "; ".join(drivers))
                    # (since 2026-09-24 a region whose set burden would rest on some members alone is incomplete and
                    # not ranked, so no ranked row needs the old '[set burden rests on the other member(s) alone]'
                    # note; incomplete regions are listed in results/tables/scoring_incomplete_burden.csv)
                    moves.append({**combo_key(cb), "contrast": cname, "top_n": N, "movement": kind,
                                  "geo_id": geo["geo_id"].iloc[i], "geo_name": geo["geo_name"].iloc[i],
                                  "state_abbr": geo["state_abbr"].iloc[i], "rank_full": full_rank[i],
                                  "rank_contrast": sc_rank[i],
                                  **{f"{k}_pct": comps[k][i] for k in O.COMPONENTS},
                                  **{f"contribution_{k}": c[k] for k in O.COMPONENTS}, "reason": reason})
    return summary, moves


EVIDENCE_LABEL = ("T5-4 evidence_weighted (adds measurement_evidence + expected_yield; docs/ANALYSIS_PLAN_METRIC_LINK.md) "
                  "vs equal")


def test5_evidence_combo(cb: O.Combo, ev: dict | None = None) -> tuple[list[dict], list[dict]]:
    """T5-4: the evidence_weighted ranking vs the equal-weight ranking of the same combination (plan section 6).

    Within one combination measurement_evidence is a constant, so any change comes from expected_yield and the weight
    change. A combination whose measurement performance is not known has no evidence_weighted ranking: one row per
    top-N says so (overlap etc. NaN) instead of comparing it as if it were ranked."""
    ev = ev or O.evaluate(cb)
    named = O.named_composites(cb, ev)
    ws = O.weight_sets().get(O.EVIDENCE_WEIGHT_SET)
    elig = cb.eligible
    geo = cb.geo
    perf = cb.perf or {}
    status = perf.get("performance_status")
    base = {**combo_key(cb), "contrast": O.EVIDENCE_WEIGHT_SET, "contrast_label": EVIDENCE_LABEL,
            "measurement_performance_status": status, "measurement_performance_tier": perf.get("quality_tier"),
            "measurement_performance_record_id": perf.get("selected_record_id")}
    full = named["equal"]
    ew = named.get(O.EVIDENCE_WEIGHT_SET)
    if ws is None or ew is None or not (elig & ~np.isnan(ew)).any():
        note = ("not ranked under evidence_weighted: measurement performance "
                + ("UNKNOWN" if status not in ("known", "partial") else status))
        return [{**base, "top_n": N, "status": note, "overlap": np.nan, "jaccard": np.nan} for N in TOP_NS], []
    comps = {k: np.asarray(v, dtype=float).ravel() for k, v in ev["components"].items()}
    contrib = contributions({k: comps[k] for k in ws}, ws, elig)
    full_rank = O.rank_rows(full, elig)[0]
    ew_rank = O.rank_rows(ew, elig)[0]
    tau, rho, n_corr = rank_corr(ew, full, elig)
    summary, moves = [], []
    M_elig = int((elig & ~np.isnan(full)).sum())
    for N in TOP_NS:
        ft, ct = top_set(full, elig, N), top_set(ew, elig, N)
        ov = len(set(ft) & set(ct))
        jac = jaccard(ft, ct)
        summary.append({**base, "top_n": N, "status": "ranked", "n_eligible_full": M_elig,
                        "n_eligible_contrast": int((elig & ~np.isnan(ew)).sum()), "overlap": ov, "jaccard": jac,
                        "expected_overlap_independent": N * N / M_elig if M_elig else np.nan,
                        "kendall_tau_b": tau, "spearman": rho, "n_rank_corr": n_corr, "n_entering": N - ov,
                        "n_leaving": N - ov, "material_change_jaccard_lt_0_5": bool(jac < MATERIAL_JACCARD)})
        for kind, rows in (("entering_evidence_weighted_top", sorted(set(ct) - set(ft), key=lambda i: ew_rank[i])),
                           ("leaving_equal_top", sorted(set(ft) - set(ct), key=lambda i: full_rank[i]))):
            for i in rows:
                c = {k: float(contrib[k][i]) for k in ws}
                ranked = sorted([kv for kv in c.items() if kv[1] == kv[1]], key=lambda kv: kv[1],
                                reverse=(kind == "entering_evidence_weighted_top"))
                drivers = [f"{k} {comps[k][i]:.2f} ({v:+.3f})" for k, v in ranked[:2]]
                reason = ((f"enters the evidence_weighted top-{N} (rank {ew_rank[i]:.0f}; equal rank "
                           f"{full_rank[i]:.0f}); largest positive contributions: ") if kind.startswith("entering")
                          else (f"leaves: equal rank {full_rank[i]:.0f}, evidence_weighted rank {ew_rank[i]:.0f}; "
                                "smallest (most negative) contributions: ")) + "; ".join(drivers)
                moves.append({**combo_key(cb), "contrast": O.EVIDENCE_WEIGHT_SET, "top_n": N, "movement": kind,
                              "geo_id": geo["geo_id"].iloc[i], "geo_name": geo["geo_name"].iloc[i],
                              "state_abbr": geo["state_abbr"].iloc[i], "rank_full": full_rank[i],
                              "rank_contrast": ew_rank[i], **{f"{k}_value": comps[k][i] for k in ws},
                              **{f"contribution_{k}": c[k] for k in ws}, "reason": reason})
    return summary, moves


# ------------------------------------------------------------------------------------------------------------------
# Test 6 (a): geography-label shuffles
# ------------------------------------------------------------------------------------------------------------------

def permutation(n: int, groups: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    """A random permutation of 0..n-1, national (groups=None) or within groups (perm[i] in the group of i)."""
    if groups is None:
        return rng.permutation(n)
    codes = pd.factorize(np.asarray(groups))[0]
    a = np.argsort(codes, kind="stable")              # positions, grouped
    b = np.lexsort((rng.random(n), codes))            # same group blocks, random order inside each block
    perm = np.empty(n, dtype=np.int64)
    perm[a] = b
    return perm


def permutation_among(mask: np.ndarray, groups: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    """A permutation of 0..n-1 that moves only the rows in `mask` among themselves (national, or within groups);
    rows outside the mask map to themselves. With an all-True mask it equals permutation(n, groups, rng) draw for
    draw. Used for the burden shuffle, so incomplete-burden regions stay unranked in every draw."""
    mask = np.asarray(mask, dtype=bool)
    n = len(mask)
    if mask.all():
        return permutation(n, groups, rng)
    idx = np.flatnonzero(mask)
    perm = np.arange(n)
    perm[idx] = idx[permutation(len(idx), None if groups is None else np.asarray(groups)[idx], rng)]
    return perm


def _stats(cb: O.Combo, comp: np.ndarray, obs: np.ndarray, obs_top: np.ndarray, real: dict) -> dict:
    elig = cb.eligible
    out = {"rho_with_observed_composite": spearman(comp, obs, elig),
           "top25_overlap_with_observed": len(set(top_set(comp, elig, 25)) & set(obs_top)),
           "top10_overlap_with_observed": len(set(top_set(comp, elig, 10)) & set(top_set(obs, elig, 10)))}
    for k, v in real.items():
        out[f"S_rho_with_real_{k}"] = spearman(comp, v, elig)
    return out


def label_shuffle(cb: O.Combo, what: str, scheme: str, n_draws: int = N_LABEL_SHUFFLES, seed: int = SEED
                  ) -> tuple[dict, list[dict]]:
    ev = O.evaluate(cb)
    obs = ev["composite"][0]
    elig = cb.eligible
    obs_top = top_set(obs, elig, 25)
    real = {"burden_pct": ev["burden_pct"][0], "vulnerability_pct": cb.vulnerability_pct}
    key = real["burden_pct"] if what == "burden" else real["vulnerability_pct"]
    rng = np.random.default_rng(seed + (11 if what == "burden" else 13) + (0 if scheme == "national" else 1000))
    groups = None if scheme == "national" else cb.geo["state_fips"].astype(str).to_numpy()
    draws = []
    complete = cb.burden_complete if cb.burden_complete is not None else np.ones(len(cb.geo), bool)
    for d in range(n_draws):
        if what == "burden":
            # burden values move only among complete-burden regions (deviation 7, 2026-09-24): the ranked set is
            # the same in every draw
            perm = permutation_among(complete, groups, rng)
            e2 = O.evaluate(cb, member_values=[m.value[perm] for m in cb.members])
        else:
            perm = permutation(len(cb.geo), groups, rng)
            e2 = O.evaluate(cb, vulnerability_pct=cb.vulnerability_pct[perm])
        draws.append({"draw": d + 1, **_stats(cb, e2["composite"][0], obs, obs_top, real)})
    dd = pd.DataFrame(draws)
    s_obs = spearman(obs, key, elig)
    col = f"S_rho_with_real_{'burden_pct' if what == 'burden' else 'vulnerability_pct'}"
    null = dd[col].to_numpy(float)
    ok = ~np.isnan(null)
    summ = {**combo_key(cb), "shuffled": what, "scheme": scheme, "n_draws": n_draws,
            "S_observed": s_obs, "S_null_mean": float(np.nanmean(null)) if ok.any() else np.nan,
            "S_null_p2_5": float(np.nanpercentile(null, 2.5)) if ok.any() else np.nan,
            "S_null_p97_5": float(np.nanpercentile(null, 97.5)) if ok.any() else np.nan,
            "p_one_sided": float((1 + np.sum(null[ok] >= s_obs)) / (1 + ok.sum())) if ok.any() else np.nan,
            "rho_with_observed_mean": float(dd["rho_with_observed_composite"].mean()),
            "rho_with_observed_p2_5": float(dd["rho_with_observed_composite"].quantile(0.025)),
            "rho_with_observed_p97_5": float(dd["rho_with_observed_composite"].quantile(0.975)),
            "top25_overlap_mean": float(dd["top25_overlap_with_observed"].mean()),
            "top25_overlap_p2_5": float(dd["top25_overlap_with_observed"].quantile(0.025)),
            "top25_overlap_p97_5": float(dd["top25_overlap_with_observed"].quantile(0.975)),
            "top10_overlap_mean": float(dd["top10_overlap_with_observed"].mean())}
    summ["loses_structure"] = bool(s_obs > summ["S_null_p97_5"]) if ok.any() else None
    summ["ranking_unchanged_all_draws"] = bool((dd["rho_with_observed_composite"] > 0.999999).all())
    for d in draws:
        d.update({**combo_key(cb), "shuffled": what, "scheme": scheme})
    return summ, draws


# ------------------------------------------------------------------------------------------------------------------
# Test 6 (b): clinic-location shuffle (recomputes clinic_capacity / research_readiness / saturation)
# ------------------------------------------------------------------------------------------------------------------

_W: dict = {}


def _clinic_init(cid: str, mid: str, level: str):
    cb = O.make_combo(cid, mid, level)
    ev = O.evaluate(cb)
    S = O.spatial()
    fac_state = S.loc["state"].to_numpy()[S.fac_loc]
    npi_state = S.loc["state"].to_numpy()[S.npi_loc]
    _W.update(cb=cb, obs=ev["composite"][0], S=S, fac_state=fac_state, npi_state=npi_state,
              real={"clinic_capacity_pct": cb.cr["clinic_capacity_pct"],
                    "research_readiness_pct": cb.cr["research_readiness_pct"],
                    "burden_pct": ev["burden_pct"][0]})


def _perm_groups(n: int, groups: np.ndarray | None, rng) -> np.ndarray:
    return permutation(n, groups, rng)


def _clinic_draw(args) -> dict:
    scheme, d = args
    cb, S = _W["cb"], _W["S"]
    rng = np.random.default_rng(SEED + (20001 if scheme == "national" else 40001) + d)
    fp = _perm_groups(len(S.fac_loc), None if scheme == "national" else _W["fac_state"], rng)
    npp = _perm_groups(len(S.npi_loc), None if scheme == "national" else _W["npi_state"], rng)
    fm = O.facility_measures(S, cb.level, tuple(cb.cond["members"]), cb.meas, fac_perm=fp, npi_perm=npp)
    cr = O.capacity_readiness(fm, cb.level)
    if cb.cr.get("clinic_capacity_basis") == "measurement_activity":
        # activity basis: permute the activity clinicians' location slots the same way (national / within state)
        from ..facilities import activity as A
        st = A.slot_states()
        ap = _perm_groups(len(st), None if scheme == "national" else st, rng)
        cap = A.capacity_pct_for(cb.meas["measurement_id"], cb.level, cb.geo["geo_id"].to_numpy(), perm=ap)
        cr["clinic_capacity_pct"] = O.percentile_rank(cap)
    e2 = O.evaluate(cb, capacity_pct=cr["clinic_capacity_pct"], readiness_pct=cr["research_readiness_pct"],
                    saturation_pct=cr["technology_saturation_pct"])
    obs = _W["obs"]
    st = _stats(cb, e2["composite"][0], obs, top_set(obs, cb.eligible, 25), _W["real"])
    st["rho_capacity_shuffled_vs_observed"] = spearman(cr["clinic_capacity_pct"], _W["real"]["clinic_capacity_pct"],
                                                       cb.eligible)
    st["rho_readiness_shuffled_vs_observed"] = spearman(cr["research_readiness_pct"],
                                                        _W["real"]["research_readiness_pct"], cb.eligible)
    st["pool_size_changed"] = int(np.sum(fm["pool_facilities_n"] != cb.fm["pool_facilities_n"]))
    return {"scheme": scheme, "draw": d + 1, **st}


def clinic_shuffle(cid: str, mid: str, level: str = "county", n_national: int = N_CLINIC_NATIONAL,
                   n_within: int = N_CLINIC_WITHIN, workers: int = WORKERS) -> tuple[list[dict], list[dict]]:
    _clinic_init(cid, mid, level)
    cb = _W["cb"]
    obs = _W["obs"]
    tasks = [("national", d) for d in range(n_national)] + [("within_state", d) for d in range(n_within)]
    t0 = time.time()
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers, initializer=_clinic_init, initargs=(cid, mid, level)) as ex:
            draws = list(ex.map(_clinic_draw, tasks, chunksize=4))
    else:
        draws = [_clinic_draw(t) for t in tasks]
    print(f"[scoring] clinic shuffle {cid} x {mid} x {level}: {len(draws)} draws ({time.time() - t0:.0f}s)",
          flush=True)
    dd = pd.DataFrame(draws)
    summ = []
    for scheme, g in dd.groupby("scheme"):
        for comp in ("clinic_capacity_pct", "research_readiness_pct"):
            s_obs = spearman(obs, _W["real"][comp], cb.eligible)
            null = g[f"S_rho_with_real_{comp}"].to_numpy(float)
            corr_col = ("rho_capacity_shuffled_vs_observed" if comp == "clinic_capacity_pct"
                        else "rho_readiness_shuffled_vs_observed")
            summ.append({"condition_id": cid, "measurement_id": mid, "geo_level": level, "shuffled":
                         "clinic_locations (facilities + individual providers)", "scheme": scheme,
                         "component": comp, "n_draws": len(g), "S_observed": s_obs,
                         "S_null_mean": float(np.nanmean(null)), "S_null_p2_5": float(np.nanpercentile(null, 2.5)),
                         "S_null_p97_5": float(np.nanpercentile(null, 97.5)),
                         "p_one_sided": float((1 + np.sum(null >= s_obs)) / (1 + len(null))),
                         "loses_structure": bool(s_obs > np.nanpercentile(null, 97.5)),
                         "component_rho_shuffled_vs_observed_mean": float(g[corr_col].mean()),
                         "component_rho_shuffled_vs_observed_p2_5": float(g[corr_col].quantile(0.025)),
                         "component_rho_shuffled_vs_observed_p97_5": float(g[corr_col].quantile(0.975)),
                         "rho_with_observed_mean": float(g["rho_with_observed_composite"].mean()),
                         "rho_with_observed_p2_5": float(g["rho_with_observed_composite"].quantile(0.025)),
                         "rho_with_observed_p97_5": float(g["rho_with_observed_composite"].quantile(0.975)),
                         "top25_overlap_mean": float(g["top25_overlap_with_observed"].mean()),
                         "top25_overlap_p2_5": float(g["top25_overlap_with_observed"].quantile(0.025)),
                         "top25_overlap_p97_5": float(g["top25_overlap_with_observed"].quantile(0.975)),
                         "top10_overlap_mean": float(g["top10_overlap_with_observed"].mean()),
                         "pool_size_invariance_violations": int(g["pool_size_changed"].sum())})
    for d in draws:
        d.update({"condition_id": cid, "measurement_id": mid, "geo_level": level})
    return summ, draws


# ------------------------------------------------------------------------------------------------------------------
# Test 6 (c): phenotype-label permutations (restated from the wearable modules' own tables)
# ------------------------------------------------------------------------------------------------------------------

def phenotype_label_nulls() -> pd.DataFrame:
    rows = []
    f = TABLES / "nhanes_permutation_null.csv"
    if f.exists():
        for r in pd.read_csv(f).to_dict("records"):
            rows.append({"source_file": f.name, "dataset": "nhanes", "target": r["target"],
                         "model_or_score": f"{r['model']} / {r['feature_set']}", "statistic": "AUROC (repeat 1)",
                         "control": "label permutation", "observed": r["observed_auroc_repeat1"],
                         "n_permutations": r["n_permutations"], "null_mean": r["null_mean"],
                         "null_upper": r["null_p95"], "null_upper_kind": "95th percentile",
                         "p": r["empirical_p"]})
    f = TABLES / "nhanes_shuffle_control.csv"
    if f.exists():
        for r in pd.read_csv(f).to_dict("records"):
            rows.append({"source_file": f.name, "dataset": "nhanes", "target": r["target"],
                         "model_or_score": f"{r['model']} / clinical + wearable",
                         "statistic": "delta AUROC (clinical+wearable minus clinical)",
                         "control": "wearable-row shuffle", "observed": r["observed_delta_repeat1"],
                         "n_permutations": r["n_shuffles"], "null_mean": r["shuffled_delta_mean"],
                         "null_upper": r["shuffled_delta_p95"], "null_upper_kind": "95th percentile",
                         "p": r["empirical_p"]})
    f = TABLES / "stanford_uwakwe_permutation_nulls.csv"
    if f.exists():
        for r in pd.read_csv(f).to_dict("records"):
            rows.append({"source_file": f.name, "dataset": "stanford_longcovid_uwakwe2025",
                         "target": "long_covid_self_report", "model_or_score": f"{r['model']} / {r['feature_set']}",
                         "statistic": r["statistic"], "control": r["null_type"], "observed": r["observed"],
                         "n_permutations": r["n_perm"], "null_mean": r["null_mean"], "null_upper": r["null_q975"],
                         "null_upper_kind": "97.5th percentile", "p": r["p_one_sided"]})
    f = TABLES / "stanford_acute_detection.csv"
    if f.exists():
        d = pd.read_csv(f)
        d = d[d["perm_null_mean"].notna()]
        for r in d.to_dict("records"):
            rows.append({"source_file": f.name, "dataset": "stanford_covid_mishra2020_alavi2022",
                         "target": f"acute infection {r['window']}", "model_or_score": r["score"],
                         "statistic": "AUROC", "control": "label permutation", "observed": r["auroc"],
                         "n_permutations": np.nan, "null_mean": r["perm_null_mean"],
                         "null_upper": r["perm_null_q975"], "null_upper_kind": "97.5th percentile",
                         "p": r["p_one_sided"]})
    out = pd.DataFrame(rows)
    if len(out):
        out["observed_above_null_p_le_0_05"] = out["p"] <= 0.05
        out["reading"] = np.where(out["observed_above_null_p_le_0_05"],
                                  "observed exceeds its permutation null (p <= 0.05): the control removes the signal",
                                  "NOT distinguishable from its null (p > 0.05)")
    return out
