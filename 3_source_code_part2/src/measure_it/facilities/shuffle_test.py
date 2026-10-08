"""Validation Test 6 (clinic locations): shuffle facility locations and watch the matched sets change.

Pre-specified in docs/ANALYSIS_PLAN_FACILITIES.md §5. A facility's location tuple (coordinates, county, state and
the county-level HPSA / RUCC flags that belong to the place) is permuted across all geocoded candidate facilities,
nationally (primary, 200 draws) or within state (secondary, 100 draws). Facility-intrinsic attributes (taxonomy
groups, FQHC status, trial history, NIH awards) stay with the facility. For a fixed, seeded sample of target
counties the clinic matcher (measure_it.facilities.matching: same ranking and round-robin code) is re-run on the
shuffled locations and the matched sets' characteristics are compared with the observed ones.

Because a permutation only reassigns location slots, each geography's pool is the same set of slots in every
draw; the pool is precomputed once and the facilities occupying those slots change. The number of facilities in
each pool is therefore invariant (checked).

Statistics. The empirical p compares the observed mean over the sampled counties with the permutation null
(conditional on that sample). Added at review: BH-FDR q-values within each (pair, null) family of 40 tests
(`q_bh_fdr_family`), and a county bootstrap (resampling within metro/nonmetro, 4,000 draws) of the per-county
observed-minus-null difference (`geo_boot_*`), which asks whether another county sample would show the same
difference. `robust_q05_and_geo_ci` requires both.

Outputs: results/tables/test6_clinic_shuffle.csv (summary), test6_clinic_shuffle_draws.csv,
test6_clinic_shuffle_geographies.csv, test6_clinic_shuffle_checks.csv; draft figure
results/figures/drafts/test6_clinic_shuffle.png.

Reproduce: `uv run python -m measure_it.facilities.shuffle_test [--draws 200 --state-draws 100 --workers 6]`
"""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from ..config import FIGURES, SEED, TABLES, load_config, utc_now_iso
from ..store import read_table
from . import matching as M
from .registry import haversine_km

PAIRS = [("long_covid", "wearable_autonomic_activity_monitoring", "primary"),
         ("me_cfs", "wearable_autonomic_activity_monitoring", "secondary"),
         ("pots", "autonomic_function_testing", "secondary")]
DEMO_COUNTIES = ["06073"]
N_PER_STRATUM = 100
METRICS = ["mean_condition_trials", "share_with_condition_trial", "share_condition_specialty",
           "share_with_nih_condition_project", "share_with_technology_experience", "share_fqhc_or_lookalike",
           "share_in_hpsa_county", "mean_distance_km", "n_eligible_in_pool", "km_to_nearest_condition_trial_site"]
TERRITORIES = {"60", "66", "69", "72", "78"}


def target_geographies() -> pd.DataFrame:
    """Seeded stratified sample of 2024 counties (50 states + DC, population >= min_population)."""
    minpop = load_config("scoring").get("ranking", {}).get("min_population", 10000)
    r = read_table("geo_context__rucc", columns=["geo_id", "state_fips", "state_abbr", "county_name",
                                                  "population_2020", "rucc_2023"])
    g = read_table("geographies", columns=["geo_id", "geo_level", "name", "lat", "lon", "ct_legacy"])
    g = g[(g["geo_level"] == "county") & ~g["ct_legacy"]]
    r = r.merge(g, on="geo_id", how="inner")
    r = r[~r["state_fips"].isin(TERRITORIES) & r["rucc_2023"].notna() & (r["population_2020"] >= minpop)]
    r = r.sort_values("geo_id").reset_index(drop=True)
    r["stratum"] = np.where(r["rucc_2023"] <= 3, "metro", "nonmetro")
    rng = np.random.default_rng(SEED)
    parts = []
    for s in ("metro", "nonmetro"):
        sub = r[r["stratum"] == s]
        parts.append(sub.iloc[np.sort(rng.choice(len(sub), size=min(N_PER_STRATUM, len(sub)), replace=False))])
    t = pd.concat(parts)
    extra = r[r["geo_id"].isin(DEMO_COUNTIES) & ~r["geo_id"].isin(t["geo_id"])]
    t = pd.concat([t, extra.assign(demo=True)]).reset_index(drop=True)
    t["demo"] = t["geo_id"].isin(DEMO_COUNTIES)
    t["n_eligible_universe"] = t["stratum"].map(r["stratum"].value_counts())
    return t


# --------------------------------------------------------------------------------------------------------------
# worker state
# --------------------------------------------------------------------------------------------------------------

_W: dict = {}


def _init(radius_km: float, max_sites: int, geos: pd.DataFrame):
    b = M.load_base()
    feats = {}
    for cid, mid, _ in PAIRS:
        meas = M.resolve_measurement(mid)
        feats[(cid, mid)] = M.features(cid, meas)
    n = len(b.fac)
    slots, dists = [], []
    geo_dicts = []
    for r in geos.itertuples(index=False):
        geo = {"geo_id": r.geo_id, "geo_level": "county", "lat": float(r.lat), "lon": float(r.lon),
               "name": r.name, "ct_legacy": False}
        d = haversine_km(geo["lat"], geo["lon"], b.lat, b.lon)
        s = np.flatnonzero((d <= radius_km) | (b.county == r.geo_id))
        slots.append(s)
        dists.append(d[s])
        geo_dicts.append(geo)
    _W.update(b=b, feats=feats, n=n, slots=slots, dists=dists, geos=geo_dicts, radius=radius_km,
              max_sites=max_sites, stratum=geos["stratum"].to_numpy(), geo_ids=geos["geo_id"].to_numpy())


def _run_perm(perm: np.ndarray | None) -> dict:
    """perm[j] = location slot given to facility j (None = observed). Returns metrics and selections."""
    b, n = _W["b"], _W["n"]
    if perm is None:
        perm = np.arange(n)
    inv = np.empty(n, dtype=np.int64)
    inv[perm] = np.arange(n)
    hpsa = b.hpsa[perm]
    nonmetro = b.nonmetro[perm]
    lat_p, lon_p = b.lat[perm], b.lon[perm]
    dist_buf = np.full(n, np.inf)
    out = {}
    for key, F in _W["feats"].items():
        access = F.intrinsic_access + hpsa.astype(int) + nonmetro.astype(int)
        trial_fac = np.flatnonzero(F.trials_cond > 0)
        rows, sels, pools = [], [], []
        for gi, geo in enumerate(_W["geos"]):
            s = _W["slots"][gi]
            fac = inv[s]
            dist_buf[fac] = _W["dists"][gi]
            orders = M.characteristic_orders(F, fac, dist_buf, access, True)
            sel, _ = M.round_robin(orders, _W["max_sites"])
            S = np.array(sel, dtype=int)
            if len(trial_fac):
                near = float(np.min(haversine_km(geo["lat"], geo["lon"], lat_p[trial_fac], lon_p[trial_fac])))
            else:
                near = np.nan
            if len(S):
                row = [F.trials_cond[S].mean(), (F.trials_cond[S] > 0).mean(), (F.n_cond_groups[S] > 0).mean(),
                       (F.nih_core[S] > 0).mean(), (F.tech[S] > 0).mean(), b.fqhc[S].mean(), hpsa[S].mean(),
                       dist_buf[S].mean(), float(F.eligible[fac].sum()), near]
            else:
                row = [np.nan] * 7 + [np.nan, float(F.eligible[fac].sum()), near]
            rows.append(row)
            sels.append(frozenset(b.fac["facility_id"].to_numpy()[S]) if len(S) else frozenset())
            pools.append(len(fac))
            dist_buf[fac] = np.inf
        out[key] = {"metrics": np.array(rows, dtype=float), "selected": sels, "pool_sizes": np.array(pools)}
    return out


def _perm_national(seed: int) -> np.ndarray:
    return np.random.default_rng(seed).permutation(_W["n"])


def _perm_within_state(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    st = _W["b"].state
    perm = np.arange(_W["n"])
    for s in np.unique(st):
        idx = np.flatnonzero(st == s)
        perm[idx] = idx[rng.permutation(len(idx))]
    return perm


def _draw(args) -> tuple:
    kind, d, obs_sel = args
    seed = SEED + (1 if kind == "national" else 100000) + d
    perm = _perm_national(seed) if kind == "national" else _perm_within_state(seed)
    res = _run_perm(perm)
    packed = {}
    for key, v in res.items():
        jac = [len(a & o) / len(a | o) if (a | o) else np.nan for a, o in zip(v["selected"], obs_sel[key])]
        packed[key] = {"metrics": v["metrics"], "jaccard": np.array(jac), "pool_sizes": v["pool_sizes"]}
    return kind, d, packed


def _worker_init(radius_km, max_sites, geos):
    _init(radius_km, max_sites, geos)


# --------------------------------------------------------------------------------------------------------------
# summaries
# --------------------------------------------------------------------------------------------------------------

def _stratum_means(metrics: np.ndarray, stratum: np.ndarray) -> dict:
    out = {"all": np.nanmean(metrics, axis=0)}
    for s in ("metro", "nonmetro"):
        out[s] = np.nanmean(metrics[stratum == s], axis=0)
    out["metro_minus_nonmetro"] = out["metro"] - out["nonmetro"]
    return out


N_GEO_BOOT = 4000
Q_ALPHA = 0.05


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values (NaN-safe). Pure function (unit-tested)."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def geo_bootstrap(diff: np.ndarray, stratum: np.ndarray, n_boot: int = N_GEO_BOOT, seed: int = SEED + 7) -> dict:
    """Uncertainty from WHICH counties were sampled (the permutation null conditions on the sample).

    diff = per-geography observed value minus that geography's mean over the null draws (geos x metrics).
    Counties are resampled with replacement within stratum (the sampling design); returns, per summary stratum,
    (point estimate, 2.5th, 97.5th percentile) arrays over metrics."""
    rng = np.random.default_rng(seed)
    idx = {s: np.flatnonzero(stratum == s) for s in ("metro", "nonmetro")}

    def stats(im, inn):
        m = np.nanmean(diff[im], axis=0)
        n = np.nanmean(diff[inn], axis=0)
        return {"all": np.nanmean(diff[np.concatenate([im, inn])], axis=0), "metro": m, "nonmetro": n,
                "metro_minus_nonmetro": m - n}

    est = stats(idx["metro"], idx["nonmetro"])
    boots = {k: [] for k in est}
    for _ in range(n_boot):
        b = stats(rng.choice(idx["metro"], len(idx["metro"])), rng.choice(idx["nonmetro"], len(idx["nonmetro"])))
        for k, v in b.items():
            boots[k].append(v)
    out = {}
    for k, v in est.items():
        arr = np.array(boots[k])
        out[k] = (v, np.nanpercentile(arr, 2.5, axis=0), np.nanpercentile(arr, 97.5, axis=0))
    return out


def summarize(obs: dict, draws: list, stratum: np.ndarray, geos: pd.DataFrame) -> tuple:
    summary, draw_rows, geo_rows = [], [], []
    for (cid, mid, role) in PAIRS:
        key = (cid, mid)
        om = _stratum_means(obs[key]["metrics"], stratum)
        for kind in ("national", "within_state"):
            ds = [p for k, d, p in draws if k == kind]
            if not ds:
                continue
            nm = [_stratum_means(p[key]["metrics"], stratum) for p in ds]
            stack = np.stack([p[key]["metrics"] for p in ds])  # draws x geos x metrics
            null_geo = np.nanmean(stack, axis=0)
            gb = geo_bootstrap(obs[key]["metrics"] - null_geo, stratum)
            jac = np.array([np.nanmean(p[key]["jaccard"]) for p in ds])
            jac_s = {s: np.array([np.nanmean(p[key]["jaccard"][stratum == s]) for p in ds]) for s in ("metro", "nonmetro")}
            for di, (p, m) in enumerate(zip(ds, nm)):
                for s, vals in m.items():
                    for mi, metric in enumerate(METRICS):
                        draw_rows.append({"condition_id": cid, "measurement_id": mid, "null": kind, "draw": di + 1,
                                          "stratum": s, "metric": metric, "value": vals[mi]})
                draw_rows.append({"condition_id": cid, "measurement_id": mid, "null": kind, "draw": di + 1,
                                  "stratum": "all", "metric": "jaccard_vs_observed", "value": jac[di]})
            for s in ("all", "metro", "nonmetro", "metro_minus_nonmetro"):
                for mi, metric in enumerate(METRICS):
                    null = np.array([m[s][mi] for m in nm])
                    o = om[s][mi]
                    nmean = float(np.nanmean(null))
                    p = (1 + np.sum(np.abs(null - nmean) >= abs(o - nmean))) / (1 + np.sum(~np.isnan(null)))
                    lo, hi = np.nanpercentile(null, [2.5, 97.5])
                    summary.append({"condition_id": cid, "measurement_id": mid, "pair_role": role, "null": kind,
                                    "n_draws": len(ds), "stratum": s, "metric": metric, "observed": o,
                                    "null_mean": nmean, "null_p2_5": lo, "null_p97_5": hi,
                                    "observed_minus_null": o - nmean, "p_two_sided_empirical": p,
                                    "outside_null_95": bool(o < lo or o > hi),
                                    "geo_boot_diff": float(gb[s][0][mi]),
                                    "geo_boot_ci_lo": float(gb[s][1][mi]), "geo_boot_ci_hi": float(gb[s][2][mi]),
                                    "n_geographies": int(np.sum(stratum == s)) if s in ("metro", "nonmetro") else
                                    int(len(stratum))})
                if s in ("all", "metro", "nonmetro"):
                    j = jac if s == "all" else jac_s[s]
                    summary.append({"condition_id": cid, "measurement_id": mid, "pair_role": role, "null": kind,
                                    "n_draws": len(ds), "stratum": s, "metric": "jaccard_vs_observed",
                                    "observed": 1.0, "null_mean": float(np.nanmean(j)),
                                    "null_p2_5": float(np.nanpercentile(j, 2.5)),
                                    "null_p97_5": float(np.nanpercentile(j, 97.5)),
                                    "observed_minus_null": 1.0 - float(np.nanmean(j)),
                                    "p_two_sided_empirical": np.nan, "outside_null_95": bool(np.nanpercentile(j, 97.5) < 1),
                                    "n_geographies": int(len(stratum)) if s == "all" else int(np.sum(stratum == s))})
            if kind == "national":
                for gi, g in enumerate(geos.itertuples(index=False)):
                    row = {"condition_id": cid, "measurement_id": mid, "geo_id": g.geo_id, "name": g.name,
                           "state_abbr": g.state_abbr, "stratum": g.stratum, "rucc_2023": g.rucc_2023,
                           "population_2020": g.population_2020, "demo": bool(g.demo),
                           "pool_size": int(obs[key]["pool_sizes"][gi]),
                           "jaccard_mean_national": float(np.nanmean([p[key]["jaccard"][gi] for p in ds]))}
                    for mi, metric in enumerate(METRICS):
                        row[f"obs__{metric}"] = obs[key]["metrics"][gi, mi]
                        row[f"null__{metric}"] = null_geo[gi, mi]
                    geo_rows.append(row)
    summ = pd.DataFrame(summary)
    # Multiplicity (added at review): 10 metrics x 4 strata = 40 p-values per (pair, null) family; BH-FDR within it.
    test = summ["metric"] != "jaccard_vs_observed"
    summ["q_bh_fdr_family"] = np.nan
    for _, idx in summ[test].groupby(["condition_id", "measurement_id", "null"]).groups.items():
        summ.loc[idx, "q_bh_fdr_family"] = bh_qvalues(summ.loc[idx, "p_two_sided_empirical"].to_numpy(float))
    summ["geo_boot_excludes_0"] = test & ((summ["geo_boot_ci_lo"] > 0) | (summ["geo_boot_ci_hi"] < 0))
    # A difference is called robust only if it survives BOTH the multiplicity adjustment (permutation null,
    # conditional on the sampled counties) AND resampling of the counties themselves.
    summ["robust_q05_and_geo_ci"] = test & (summ["q_bh_fdr_family"] <= Q_ALPHA) & summ["geo_boot_excludes_0"]
    return summ, pd.DataFrame(draw_rows), pd.DataFrame(geo_rows)


def draft_figure(summary: pd.DataFrame, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = summary[(summary["pair_role"] == "primary") & (summary["null"] == "national") &
                summary["stratum"].isin(["metro", "nonmetro"]) & (summary["metric"] != "jaccard_vs_observed")]
    metrics = [m for m in METRICS if m in set(s["metric"])]
    fig, axes = plt.subplots(2, 5, figsize=(17, 6.5))
    for ax, m in zip(axes.flat, metrics):
        d = s[s["metric"] == m].set_index("stratum").reindex(["metro", "nonmetro"])
        x = np.arange(len(d))
        ax.errorbar(x + 0.08, d["null_mean"], yerr=[d["null_mean"] - d["null_p2_5"], d["null_p97_5"] - d["null_mean"]],
                    fmt="o", color="#9aa0a6", capsize=4, label="shuffled locations (mean, 95% of 200 draws)")
        ax.plot(x - 0.08, d["observed"], "D", color="#1f5fa8", label="observed")
        ax.set_xticks(x, d.index)
        ax.set_title(m.replace("_", " "), fontsize=9)
        ax.grid(alpha=0.3)
    axes.flat[0].legend(fontsize=7, loc="best")
    fig.suptitle("Test 6 (draft): matched-set characteristics, observed vs national location shuffle\n"
                 "Long COVID x wearable autonomic/activity monitoring; 100 metro + 100 nonmetro sampled counties "
                 "(+ San Diego); radius 50 km, 15 sites", fontsize=10)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def run(draws: int = 200, state_draws: int = 100, workers: int = 6, radius_km: float | None = None,
        max_sites: int | None = None) -> dict:
    t0 = time.time()
    cfg = load_config("scoring").get("clinic_matching", {})
    radius_km = float(radius_km if radius_km is not None else cfg.get("default_radius_km", 50))
    max_sites = int(max_sites if max_sites is not None else cfg.get("max_sites_per_geography", 15))
    geos = target_geographies()
    _init(radius_km, max_sites, geos)
    stratum = _W["stratum"]
    obs = _run_perm(None)
    # check: the precomputed-slot path reproduces matching.select exactly on observed locations
    b = _W["b"]
    mism = 0
    for (cid, mid), F in _W["feats"].items():
        for gi, geo in enumerate(_W["geos"][:25]):
            r = M.select(F, b.lat, b.lon, b.county, b.state, b.hpsa, b.nonmetro, geo, radius_km, max_sites)
            ids = frozenset(b.fac["facility_id"].to_numpy()[np.array(r["selected"], dtype=int)]) if r["selected"] else frozenset()
            mism += int(ids != obs[(cid, mid)]["selected"][gi])
    obs_sel = {k: v["selected"] for k, v in obs.items()}
    tasks = [("national", d, obs_sel) for d in range(draws)] + [("within_state", d, obs_sel) for d in range(state_draws)]
    results = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init,
                             initargs=(radius_km, max_sites, geos)) as ex:
        for i, r in enumerate(ex.map(_draw, tasks, chunksize=4)):
            results.append(r)
            if (i + 1) % 25 == 0:
                print(f"{i + 1}/{len(tasks)} draws ({time.time() - t0:.0f}s)", flush=True)
    # invariance: pool sizes identical under every permutation
    inv_viol = 0
    for k, d, p in results:
        for key in p:
            inv_viol += int(np.sum(p[key]["pool_sizes"] != obs[key]["pool_sizes"]))
    summary, draw_df, geo_df = summarize(obs, results, stratum, geos)
    TABLES.mkdir(parents=True, exist_ok=True)
    summary.to_csv(TABLES / "test6_clinic_shuffle.csv", index=False)
    draw_df.to_csv(TABLES / "test6_clinic_shuffle_draws.csv", index=False)
    geo_df.to_csv(TABLES / "test6_clinic_shuffle_geographies.csv", index=False)
    checks = pd.DataFrame([
        {"check": "pool_size_invariance_violations", "value": inv_viol,
         "note": "number of (draw, pair, geography) pools whose size changed under permutation; must be 0"},
        {"check": "slot_path_vs_matching_select_mismatches_first25_geos", "value": mism,
         "note": "observed selections from the precomputed-slot path that differ from matching.select; must be 0"},
        {"check": "n_target_geographies", "value": len(geos), "note": json.dumps(geos["stratum"].value_counts().to_dict())},
        {"check": "n_national_draws", "value": draws, "note": f"seed base {SEED}+1"},
        {"check": "n_within_state_draws", "value": state_draws, "note": f"seed base {SEED}+100000"},
        {"check": "radius_km", "value": radius_km, "note": "scoring.yaml clinic_matching.default_radius_km"},
        {"check": "max_sites", "value": max_sites, "note": "scoring.yaml clinic_matching.max_sites_per_geography"},
        {"check": "n_candidate_facilities", "value": _W["n"], "note": "geocoded clinic-candidate facilities permuted"},
        {"check": "runtime_s", "value": round(time.time() - t0, 1), "note": utc_now_iso()},
    ])
    checks.to_csv(TABLES / "test6_clinic_shuffle_checks.csv", index=False)
    geos.to_csv(TABLES / "test6_clinic_shuffle_target_geographies.csv", index=False)
    draft_figure(summary, FIGURES / "drafts" / "test6_clinic_shuffle.png")
    print(checks.to_string(index=False))
    return {"summary_rows": len(summary), "invariance_violations": inv_viol, "select_mismatches": mism}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Test 6: clinic-location shuffle negative control")
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--state-draws", type=int, default=100)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args(argv)
    run(draws=a.draws, state_draws=a.state_draws, workers=min(a.workers, 6))


if __name__ == "__main__":
    main()
